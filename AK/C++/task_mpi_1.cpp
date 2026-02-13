#include <iostream>
#include <mpi.h>
using namespace std;

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int i_node = 0;        // Центр "звезды"
  int j_node = size - 1; // Один из лучей

  // 1. Передача i -> j и обратно
  if (rank == i_node) {
    int send_val = 123, recv_val;
    MPI_Send(&send_val, 1, MPI_INT, j_node, 0, MPI_COMM_WORLD);
    MPI_Recv(&recv_val, 1, MPI_INT, j_node, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    cout << "[P" << rank << "] Обмен с P" << j_node
         << " завершен. Получено: " << recv_val << endl;
  } else if (rank == j_node) {
    int recv_val, send_val = 321;
    MPI_Recv(&recv_val, 1, MPI_INT, i_node, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    MPI_Send(&send_val, 1, MPI_INT, i_node, 0, MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  // 2. Рассылка различных данных от i остальным
  if (rank == i_node) {
    for (int k = 0; k < size; k++) {
      if (k != i_node) {
        int data = 1000 + k;
        MPI_Send(&data, 1, MPI_INT, k, 1, MPI_COMM_WORLD);
      }
    }
  } else {
    int received_data;
    MPI_Recv(&received_data, 1, MPI_INT, i_node, 1, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    cout << "[P" << rank << "] Получил уникальные данные: " << received_data
         << endl;
  }

  MPI_Finalize();
  return 0;
}
