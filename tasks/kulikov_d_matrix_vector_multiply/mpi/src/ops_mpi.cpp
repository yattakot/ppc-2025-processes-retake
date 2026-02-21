#include "kulikov_d_matrix_vector_multiply/mpi/include/ops_mpi.hpp"

#include <mpi.h>

#include <numeric>
#include <vector>

#include "kulikov_d_matrix_vector_multiply/common/include/common.hpp"
#include "util/include/util.hpp"

namespace kulikov_d_matrix_vector_multiply {

KulikovDMatrixMultiplyMPI::KulikovDMatrixMultiplyMPI(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KulikovDMatrixMultiplyMPI::ValidationImpl() {
  const auto &input = GetInput();

  if (input.rows <= 0 || input.cols <= 0) {
    return false;
  }

  if (input.matrix.size() != static_cast<size_t>(input.rows * input.cols)) {
    return false;
  }

  if (input.vector.size() != static_cast<size_t>(input.cols)) {
    return false;
  }

  return GetOutput().empty();
}

bool KulikovDMatrixMultiplyMPI::PreProcessingImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    GetOutput().assign(GetInput().rows, 0);
  }

  return true;
}

bool KulikovDMatrixMultiplyMPI::RunImpl() {
  const auto &input = GetInput();

  int rank = 0, size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int rows = input.rows;
  int cols = input.cols;

  std::vector<int> vec = input.vector;
  MPI_Bcast(vec.data(), cols, MPI_INT, 0, MPI_COMM_WORLD);

  int base_rows = rows / size;
  int remainder = rows % size;
  int local_rows = base_rows + (rank < remainder ? 1 : 0);
  int start_row = rank * base_rows + std::min(rank, remainder);

  std::vector<int> local_result(local_rows, 0);
  for (int i = 0; i < local_rows; ++i) {
    int sum = 0;
    int global_row = start_row + i;
    for (int j = 0; j < cols; ++j) {
      sum += input.matrix[global_row * cols + j] * vec[j];
    }
    local_result[i] = sum;
  }

  std::vector<int> recv_counts(size, base_rows);
  for (int p = 0; p < remainder; ++p) {
    recv_counts[p] += 1;
  }

  std::vector<int> displs(size, 0);
  for (int i = 1; i < size; ++i) {
    displs[i] = displs[i - 1] + recv_counts[i - 1];
  }

  if (rank == 0) {
    GetOutput().assign(rows, 0);
  }

  MPI_Gatherv(local_result.data(), local_rows, MPI_INT, GetOutput().data(), recv_counts.data(), displs.data(), MPI_INT,
              0, MPI_COMM_WORLD);

  return true;
}

bool KulikovDMatrixMultiplyMPI::PostProcessingImpl() {
  return true;
}

}  // namespace kulikov_d_matrix_vector_multiply
