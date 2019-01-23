#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);
 
  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
 
  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
 
 
  if(world_rank == 0)
  {
    // int length = ;
    // int * arr = (int *)malloc(length*sizeof(int));
    // MPI_Send(arr,length,MPI_INT,destination,0,MPI_COMM_WORLD);


    int a [] = {42,43,44};
    int destination = 1;
    MPI_Send(&a,3,MPI_INT,destination,0,MPI_COMM_WORLD);
    printf("Successfully sent!\n");
  }
  else
  {
    int b[3];
    MPI_Status st;
    int source = 0;
    MPI_Recv(&b,3,MPI_INT,source,0,MPI_COMM_WORLD,&st);

    printf("Successfully received!\n");

    for(int i = 0; i<3; i++)
      printf("%d\n",b[i]); 


  }
 
  // Finalize the MPI environment.
  MPI_Finalize();
}
