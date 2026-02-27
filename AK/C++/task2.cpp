#include <iostream>
#include <mpi.h>
#include <vector>

using namespace std;

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int n = 2; // размер локального куска
  vector<int> B(n, rank);
  vector<int> A(n * size);

  // Сбор данных "каждый каждому"
  MPI_Allgather(B.data(), n, MPI_INT, A.data(), n, MPI_INT, MPI_COMM_WORLD);

  if (rank == 0) {
    cout << "Результат Allgather в процессе 0: ";
    for (int val : A)
      cout << val << " ";
    cout << endl;
  }

  MPI_Finalize();
  return 0;
}
