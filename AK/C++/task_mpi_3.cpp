#include <iostream>
#include <mpi.h>
#include <vector>

using namespace std;

void mat_sq(const vector<double> &src, vector<double> &res, int n, int start,
            int end) {
  for (int i = start; i < end; ++i) {
    for (int j = 0; j < n; ++j) {
      double sum = 0;
      for (int k = 0; k < n; ++k)
        sum += src[i * n + k] * src[k * n + j];
      res[(i - start) * n + j] = sum;
    }
  }
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int N = 100;
  int rows = N / size;
  int start = rank * rows;
  int end = (rank == size - 1) ? N : (rank + 1) * rows;
  int local_n = (end - start) * N;

  vector<double> A(N * N, 1.0), B(N * N, 2.0), C(N * N, 3.0), D(N * N);
  vector<double> pA(local_n), pB(local_n), pC(local_n), pD(local_n);

  mat_sq(A, pA, N, start, end);
  mat_sq(B, pB, N, start, end);
  mat_sq(C, pC, N, start, end);

  for (int i = 0; i < local_n; i++)
    pD[i] = pA[i] + pB[i] + pC[i];

  vector<int> recvcounts(size), displs(size);
  if (rank == 0) {
    for (int i = 0; i < size; ++i) {
      int s = i * (N / size);
      int e = (i == size - 1) ? N : (i + 1) * (N / size);
      recvcounts[i] = (e - s) * N;
      displs[i] = s * N;
    }
  }

  MPI_Gatherv(pD.data(), local_n, MPI_DOUBLE, D.data(), recvcounts.data(),
              displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (rank == 0)
    cout << "D[0][0] = " << D[0] << " (ожидалось "
         << (1 * 1 * N + 2 * 2 * N + 3 * 3 * N) << ")" << endl;

  MPI_Finalize();
  return 0;
}
