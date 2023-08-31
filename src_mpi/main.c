#include <mpi/mpi.h>
#include <stdio.h>

int main(int argc, char **argv) {
  int process_rank, size_of_cluster;
  int distro_array[16] = {39, 72, 129, 42, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2};
  int scattered_data;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size_of_cluster);
  MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);

  MPI_Scatter(&distro_array, 1, MPI_INT, &scattered_data, 1, MPI_INT, 0,
              MPI_COMM_WORLD);

  printf("Process %d received %d\n", process_rank, scattered_data);

  MPI_Finalize();

  return 0;
}
