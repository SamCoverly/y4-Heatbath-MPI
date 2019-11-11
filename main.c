#include <stdio.h>
#include <stdbool.h>
#include <mpi.h>
#include <stdlib.h>

#define BEGIN 1
#define END 2
#define EQUI 3
#define MASTER 0

int main(int argc, char *argv[]) {
    double Rinner,Router,Linner,Louter; //define inner and outer radius, inner and outer length
    int taskid,numthreads,numworkers,iterations,gridslice,rc,i,j;
    MPI_Status status;
    Rinner = 10;
    Router = 100;
    Linner = 10;
    Louter = 100;
    int Rmid = ((Rinner/Router) * 100); //calculates grid points at which inner conductor ends
    int Lmid = ((Linner/Louter) * 100);
    double psi0 = 100; //define potential of inner conductor
    double psiold; //stores old value of potential, used to determine whether equilibrium is met
    double psi[101][101] = {1}; //grid storing potential values
    double diff; //difference between old and new value of potential for a single grid point
    bool equilibrium[numthreads] = false;//whether or not equilibrium has been met, one value for each processor
    bool above = false;

    for (int i = 0; i < Lmid; i++){ //set potential values for inner conductor to psi0
        for (int k = 0; k < Rmid; k++){
            psi[i][k] = psi0;
        }
    }

    do { //do while loop continues to iterate until equilibrium is achieved
        iterations++;
        MPI_Init(&argc,&argv); //initialise MPI
        MPI_Comm_size(MPI_COMM_WORLD,&numthreads);
        MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
        gridslice = 100/numthreads; //determines and indexes number of rows to be sent to each thread
        if (taskid == MASTER){ //for master thread
            if (numthreads != 1 && numthreads != 2 && numthreads != 4 && numthreads != 5 && numthreads != 10) {
                printf("Please use either 1,2,4,5 or 10 threads \n");
                MPI_Abort(MPI_COMM_WORLD, rc); //aborts program if incompatible number of threads being used
                exit(1);
            }
            for(i = 1;i < numthreads; i++){
                ierr = MPI_Send(&psi[(i*gridslice)][0],((gridslice+1)*100),MPI_DOUBLE,i,BEGIN,MPI_COMM_WORLD); //send segment of master array to each thread
            }
        }
        if (taskid != MASTER){ //for worker thread
            ierr = MPI_Recv(&psi[(i*gridslice)][0],((gridslice+1)*100),MPI_DOUBLE,MASTER,BEGIN,MPI_COMM_WORLD,&status); //receive slice of grid
        }
            if((gridslice*taskid)<Lmid){ //iterates over different parts of grid, depending on whether grid slice contains part or all of inner conductor
                    for(int i = (gridslice*taskid); i < Lmid; i++){ //calculates grid values to the right of inner conductor
                        for (int k = Rmid; k < 100;k++){
                            psiold = psi[i][k]; //stores old value of potential
                            psi[i][k] = (((0.25)*(psi[i+1][k]+psi[i-1][k]+psi[i][k+1]+psi[i][k-1]))+((0.125/k)*(psi[i][k+1]-psi[i][k-1]))); //calculate grid value
                            diff = psi[i][k] - psiold; //determines difference between old and new value
                            if (diff <= -0.0001 || diff >= 0.0001){ //if difference is greater than some arbitrary precision, equilibrium = false
                                equilibrium[taskid] = false;
                            }
                        }
                    }
                    for(int i = Lmid; i < (gridslice*(taskid+1));i++){ //calculates r=0 grid values using special case
                        psiold = psi[i][0]; //stores old value of potential
                        psi[i][0] = (((0.667)*psi[i][1])+((0.167)*(psi[i+1][0]+psi[i-1][0]))); //calculate grid value
                        diff = psi[i][0] - psiold; //determines difference between old and new value
                        if (diff <= -0.0001 || diff >= 0.0001){ //if difference is greater than some arbitrary precision, equilibrium = false
                            equilibrium[taskid] = false;
                        }
                    }
                    for(int i = Lmid; i < (gridslice*(taskid+1)); i++){ //calculates grid values with L greater than Lmid
                    for (int k = 1; k < 100;k++){
                        psiold = psi[i][k];//stores old value of potential
                        psi[i][k] = ((0.25*(psi[i+1][k]+psi[i-1][k]+psi[i][k+1]+psi[i][k-1]))+((0.125/k)*(psi[i][k+1]-psi[i][k-1]))); //calculate grid value
                        diff = psi[i][k] - psiold; //determines difference between old and new value
                        if (diff <= -0.0001 || diff >= 0.0001){ //if difference is greater than some arbitrary precision, equilibrium = false
                            equilibrium = false;
                        }
                    }
                }
            }else{
                for(int i = (gridslice*taskid); i < (gridslice*(taskid+1));i++){ //calculates r=0 grid values using special case
                    psiold = psi[i][0]; //stores old value of potential
                    psi[i][0] = (((0.667)*psi[i][1])+((0.167)*(psi[i+1][0]+psi[i-1][0]))); //calculate grid value
                    diff = psi[i][0] - psiold; //determines difference between old and new value
                    if (diff <= -0.0001 || diff >= 0.0001){ //if difference is greater than some arbitrary precision, equilibrium = false
                        equilibrium[taskid] = false;
                    }
                }
                for(int i = (gridslice*taskid); i < (gridslice*(taskid+1)); i++){ //calculates remaining grid values
                    for (int k = 1; k < 100;k++){
                        psiold = psi[i][k]; //stores old value of potential
                        psi[i][k] = ((0.25*(psi[i+1][k]+psi[i-1][k]+psi[i][k+1]+psi[i][k-1]))+((0.125/k)*(psi[i][k+1]-psi[i][k-1]))); //calculate grid value
                        diff = psi[i][k] - psiold;  //determines difference between old and new value
                        if (diff <= -0.0001 || diff >= 0.0001){ //if difference is greater than some arbitrary precision, equilibrium = false
                            equilibrium = false;
                        }
                    }
                }
            }
            if (taskid == MASTER){ //master thread retrieves all calculated grid values from other threads, as well as equilibrium condition from each thread
                for(i = 1; i < numthreads; i++){
                    ierr = MPI_Recv(&psi[(i*gridslice)][0],((gridslice+1)*100),MPI_DOUBLE,i,END,MPI_COMM_WORLD,&status); //gather updated values into master array
                    ierr = MPI_Recv(&equilibrium[i],1,MPI_C_BOOL,i,EQUI,MPI_COMM_WORLD,&status);
                    if (equilibrium[i] == false){ //if equilibrium is not achieved for any single thread, master equililbrium value = false
                        equilibrium[MASTER] = false;
                    }
                }
                MPI_Finalize();
            }
            if (taskid != MASTER){ //worked threads send updated values back to master thread
                ierr = MPI_Send(&psi[i*gridslice][0],((gridslice+1)*100),MPI_DOUBLE,MASTER,END,MPI_COMM_WORLD);
                ierr = MPI_Send(&equilibrium[taskid],1,MPI_C_BOOL,MASTER,EQUI,MPI_COMM_WORLD);
                MPI_Finalize();
            }
             
    } while (equilibrium[MASTER] == false);
    printf("after %d iterations: \n",iterations); //prints number of iterations
    for (int i = 0; i < 101; i++){ //prints grid with values
        for (int j = 0; j < 101; j++){
            printf("%f ",psi[i][j]);
        }
        printf("\n");
    }
}