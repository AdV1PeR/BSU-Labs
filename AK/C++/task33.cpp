#include <iostream>
#include <mpi.h>
#include <vector>
#include <iomanip>

using namespace std;

// Функция для умножения части матрицы (строки от start до end) на всю матрицу
// Вычисляет кусок результирующей матрицы R = M * M
void multiply_chunk(const vector<double>& M, vector<double>& res_chunk, int n, int start_row, int end_row) {
    for (int i = start_row; i < end_row; ++i) {
        for (int j = 0; j < n; ++j) {
            double sum = 0;
            for (int k = 0; k < n; ++k) {
                sum += M[i * n + k] * M[k * n + j];
            }
            res_chunk[(i - start_row) * n + j] = sum;
        }
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int N = 600; // Размер матрицы (можно менять для тестов: 500, 800, 1000)
    
    // Векторы для хранения полных матриц (только у корня)
    vector<double> A, B, C, D;
    
    if (rank == 0) {
        A.resize(N * N, 1.0); // Заполнение тестовыми данными
        B.resize(N * N, 2.0);
        C.resize(N * N, 3.0);
        D.resize(N * N);
        cout << "Запуск на " << size << " процессах. Размер матрицы: " << N << "x" << N << endl;
    } else {
        // У остальных память под полные матрицы A, B, C для перемножения
        A.resize(N * N);
        B.resize(N * N);
        C.resize(N * N);
    }

    // Замеряем время начала вычислений
    double start_time = MPI_Wtime();

    // 1. Рассылаем исходные матрицы всем процессам (т.к. для A^2 нужны все данные)
    MPI_Bcast(A.data(), N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(B.data(), N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(C.data(), N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // 2. Определяем диапазон строк для каждого процесса
    int rows_per_proc = N / size;
    int start_row = rank * rows_per_proc;
    int end_row = (rank == size - 1) ? N : (rank + 1) * rows_per_proc;
    int num_rows = end_row - start_row;

    // Временные буферы для локальных вычислений
    vector<double> pA2(num_rows * N), pB2(num_rows * N), pC2(num_rows * N), pD(num_rows * N);

    // 3. Вычисляем части A^2, B^2, C^2
    multiply_chunk(A, pA2, N, start_row, end_row);
    multiply_chunk(B, pB2, N, start_row, end_row);
    multiply_chunk(C, pC2, N, start_row, end_row);

    // 4. Складываем результаты локально: D_part = A2 + B2 + C2
    for (int i = 0; i < num_rows * N; ++i) {
        pD[i] = pA2[i] + pB2[i] + pC2[i];
    }

    // 5. Подготовка к сбору данных (Gatherv)
    vector<int> recvcounts(size), displs(size);
    if (rank == 0) {
        for (int i = 0; i < size; ++i) {
            int s = i * rows_per_proc;
            int e = (i == size - 1) ? N : (i + 1) * rows_per_proc;
            recvcounts[i] = (e - s) * N;
            displs[i] = s * N;
        }
    }

    // Сборка итоговой матрицы D в корневой процесс
    MPI_Gatherv(pD.data(), num_rows * N, MPI_DOUBLE, D.data(), 
                recvcounts.data(), displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double end_time = MPI_Wtime();

    if (rank == 0) {
        double duration = end_time - start_time;
        cout << "Время выполнения: " << fixed << setprecision(4) << duration << " сек." << endl;
        // Проверка (первый элемент): 1^2*N + 2^2*N + 3^2*N = (1+4+9)*N = 14 * N
        cout << "Контрольное значение D[0]: " << D[0] << " (Ожидалось: " << 14 * N << ")" << endl;
    }

    MPI_Finalize();
    return 0;
}

