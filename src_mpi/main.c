#include <stdio.h>
#include <mpi.h>

int main(int argc, char **argv) {
  int process_rank, size_of_cluster;

  MPI_Init(&argc, &argv);

  return 0;
}
