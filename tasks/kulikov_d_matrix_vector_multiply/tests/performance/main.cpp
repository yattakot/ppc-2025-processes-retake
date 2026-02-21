#include <gtest/gtest.h>
#include <mpi.h>

#include <random>
#include <string>
#include <vector>

#include "kulikov_d_matrix_vector_multiply/common/include/common.hpp"
#include "kulikov_d_matrix_vector_multiply/mpi/include/ops_mpi.hpp"
#include "kulikov_d_matrix_vector_multiply/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace kulikov_d_matrix_vector_multiply {

class KulikovMatrixMultiplyRunPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  InType input_data_;
  std::vector<int> expected_data_;

  void SetUp() override {
    const auto &param = GetParam();
    const std::string &test_name = std::get<1>(param);

    const bool is_task_run = test_name.find("task_run") != std::string::npos;

    int rows = 0;
    int cols = 0;

    if (is_task_run) {
      rows = 300;
      cols = 300;
    } else {
      rows = 1000;
      cols = 1000;
    }

    input_data_.rows = rows;
    input_data_.cols = cols;
    input_data_.matrix.resize(static_cast<size_t>(rows * cols));
    input_data_.vector.resize(static_cast<size_t>(cols));
    expected_data_.resize(static_cast<size_t>(rows));

    std::mt19937 gen(42);
    std::uniform_int_distribution<int> dist(-5, 5);

    for (int j = 0; j < cols; ++j) {
      input_data_.vector[j] = dist(gen);
    }

    for (int i = 0; i < rows; ++i) {
      int sum = 0;
      for (int j = 0; j < cols; ++j) {
        const int val = dist(gen);
        input_data_.matrix[i * cols + j] = val;
        sum += val * input_data_.vector[j];
      }
      expected_data_[i] = sum;
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank != 0) {
      return true;
    }

    if (output_data.size() != expected_data_.size()) {
      return false;
    }

    for (size_t i = 0; i < output_data.size(); ++i) {
      if (output_data[i] != expected_data_[i]) {
        return false;
      }
    }

    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(KulikovMatrixMultiplyRunPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, KulikovDMatrixMultiplyMPI, KulikovDMatrixMultiplySEQ>(
    PPC_SETTINGS_kulikov_d_matrix_vector_multiply);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

INSTANTIATE_TEST_SUITE_P(RunModeTests, KulikovMatrixMultiplyRunPerfTests, kGtestValues,
                         KulikovMatrixMultiplyRunPerfTests::CustomPerfTestName);

}  // namespace kulikov_d_matrix_vector_multiply
