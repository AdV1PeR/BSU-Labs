#include <iostream>
#include <mpi.h>
#include <vector>
using namespace std;

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n = 2;
    vector<int> B(n, rank);
    vector<int> A(n * size);

    // Каждый процесс рассылает B всем
    for (int i = 0; i < size; ++i) {
        if (i != rank) {
            MPI_Send(B.data(), n, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    }

    // Каждый процесс принимает B от всех
    for (int i = 0; i < size; ++i) {
        if (i == rank) {
            // копируем свой B
            for (int j = 0; j < n; ++j)
                A[i*n + j] = B[j];
        } else {
            MPI_Recv(&A[i*n], n, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    if (rank == 0) {
        cout << "Результат Send/Recv в процессе 0: ";
        for (int x : A) cout << x << " ";
        cout << endl;
    }

    MPI_Finalize();
    return 0;
}
