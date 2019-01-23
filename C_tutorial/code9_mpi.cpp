#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char** argv) {
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);
 
  // Get the number of processes
  int numProcesses;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
 
  // Get the rank of the process
  int processId;
  MPI_Comm_rank(MPI_COMM_WORLD, &processId);
 
 
  int arr[3];
  int size = 3;
  for (int i = 0; i < 3; ++i)
  {
    arr[i] = processId*10 + i;
  }

  int * tempBuf = (int *)malloc(size*sizeof(int));

  int dest = (processId-1+numProcesses)%numProcesses;
  int src  = (processId+1)%numProcesses;


  printf("Process %d is sending %d, %d, %d to Process %d\n",processId, arr[0],arr[1],arr[2], dest);

  MPI_Status st;
  MPI_Sendrecv(arr, size, MPI_INT,
                  dest, 0, //send part
                  tempBuf, size, MPI_INT,
                  src, 0, //receive part
                  MPI_COMM_WORLD, &st);

  //copy incoming data
  memcpy(arr, tempBuf, size*sizeof(int));

  printf("Process %d has received %d, %d, %d from Process %d\n",processId, arr[0],arr[1],arr[2], src);

  // Finalize the MPI environment.
  MPI_Finalize();

  free(tempBuf);
}
