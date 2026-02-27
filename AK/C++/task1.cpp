#include <iostream>
#include <mpi.h>
#include <string>

using namespace std;

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int i_node = 0; // ЦЕНТР (твой 0)
  int j_node = 2; // ЛУЧ (твой 2 - "Справа")

  // 1. Прямой обмен между центром (i) и конкретным лучом (j)
  if (rank == i_node) {
    int send_val = 123, recv_val;
    MPI_Send(&send_val, 1, MPI_INT, j_node, 0, MPI_COMM_WORLD);
    MPI_Recv(&recv_val, 1, MPI_INT, j_node, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    cout << "[P" << rank << " Центр] Обмен с P" << j_node << " (Справа) завершен." << endl;
  }
  else if (rank == j_node) {
    int recv_val, send_val = 321;
    MPI_Recv(&recv_val, 1, MPI_INT, i_node, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(&send_val, 1, MPI_INT, i_node, 0, MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  // 2. Рассылка разных данных от центра ВСЕМ лучам (1-Верх, 2-Право, 3-Низ, 4-Слева)
  if (rank == i_node) {
    for (int k = 1; k < size; k++) {
      int data = 1000 + k;
      MPI_Send(&data, 1, MPI_INT, k, 1, MPI_COMM_WORLD);
    }
  }
  else {
    int received_data;
    // ЛУЧ принимает данные ТОЛЬКО от центра (i_node)
    MPI_Recv(&received_data, 1, MPI_INT, i_node, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Массив названий для красоты вывода согласно твоей схеме
    string pos = "Луч";
    if (rank == 1) pos = "Верх";
    else if (rank == 2) pos = "Право";
    else if (rank == 3) pos = "Низ";
    else if (rank == 4) pos = "Слева";

    cout << "[P" << rank << " " << pos << "] Получил данные: " << received_data << endl;
  }

  MPI_Finalize();
  return 0;
}
