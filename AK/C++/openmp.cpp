#include <iostream>
#include <omp.h>
#include <vector>

using namespace std;

const int N = 400;

void multiply(const vector<vector<double>> &M_in,
              vector<vector<double>> &M_out) {
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++) {
      double s = 0;
      for (int k = 0; k < N; k++)
        s += M_in[i][k] * M_in[k][j];
      M_out[i][j] = s;
    }
}

int main() {
  vector<vector<double>> A(N, vector<double>(N, 1.0)),
      B(N, vector<double>(N, 2.0)), C(N, vector<double>(N, 3.0)),
      D(N, vector<double>(N, 0.0)), A2(N, vector<double>(N)),
      B2(N, vector<double>(N)), C2(N, vector<double>(N));

  omp_set_num_threads(4);

  // Способ А: Parallel For
  double t = omp_get_wtime();
#pragma omp parallel
  {
#pragma omp for nowait
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++) {
        double s = 0;
        for (int k = 0; k < N; k++)
          s += A[i][k] * A[k][j];
        A2[i][j] = s;
      }
#pragma omp for nowait
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++) {
        double s = 0;
        for (int k = 0; k < N; k++)
          s += B[i][k] * B[k][j];
        B2[i][j] = s;
      }
#pragma omp for nowait
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++) {
        double s = 0;
        for (int k = 0; k < N; k++)
          s += C[i][k] * C[k][j];
        C2[i][j] = s;
      }
#pragma omp barrier
#pragma omp for
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++)
        D[i][j] = A2[i][j] + B2[i][j] + C2[i][j];
  }
  cout << "Time Parallel For: " << omp_get_wtime() - t << endl;

  // Способ Б: Manual (равномерное распределение по нитям)
  t = omp_get_wtime();
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    int n_thr = omp_get_num_threads();
    for (int i = tid; i < N; i += n_thr) { // Циклический шаг
      for (int j = 0; j < N; j++) {
        double sA = 0, sB = 0, sC = 0;
        for (int k = 0; k < N; k++) {
          sA += A[i][k] * A[k][j];
          sB += B[i][k] * B[k][j];
          sC += C[i][k] * C[k][j];
        }
        D[i][j] = sA + sB + sC;
      }
    }
  }
  cout << "Time Manual: " << omp_get_wtime() - t << endl;

  // Способ В: Sections
  t = omp_get_wtime();
#pragma omp parallel
  {
#pragma omp sections
    {
#pragma omp section
      multiply(A, A2);
#pragma omp section
      multiply(B, B2);
#pragma omp section
      multiply(C, C2);
    }
#pragma omp for
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++)
        D[i][j] = A2[i][j] + B2[i][j] + C2[i][j];
  }
  cout << "Time Sections: " << omp_get_wtime() - t << endl;

  return 0;
}

