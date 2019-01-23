#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
// you need to install openmpi or mpich
//to compile: mpicxx code7_mpi.cpp -o code7_mpi
//to run: mpirun -n "number of processes" ./code7_mpi

int main(int argc, char** argv) {
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);
 
  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
 
  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
 
  // Get the name of the processor
  // char processor_name[MPI_MAX_PROCESSOR_NAME];
  // int name_len;
  // MPI_Get_processor_name(processor_name, &name_len);
 
  if(world_rank == 0)
  {
  // Print off a hello world message
  // printf("Hello world from rank %d"
         // " out of %d processes\n",
         // world_rank, world_size);
    printf("Hola!\n");
  }
  else 
  {
    printf("Hello from rank %d"
         " out of %d processes\n",
         world_rank, world_size);
  }


 
  // Finalize the MPI environment.
  MPI_Finalize();

  return 0;
}
