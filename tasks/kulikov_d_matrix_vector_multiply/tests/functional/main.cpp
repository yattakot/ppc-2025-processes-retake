#include <gtest/gtest.h>
#include <mpi.h>

#include <stdexcept>
#include <tuple>
#include <vector>

#include "kulikov_d_matrix_vector_multiply/common/include/common.hpp"
#include "kulikov_d_matrix_vector_multiply/mpi/include/ops_mpi.hpp"
#include "kulikov_d_matrix_vector_multiply/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"

namespace kulikov_d_matrix_vector_multiply {

class KulikovMatrixMultiplyRunFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintFuncTestName(
      const testing::TestParamInfo<ppc::util::FuncTestParam<InType, OutType, TestType>> &info) {
    const auto &test_param = std::get<2>(info.param);
    const int case_id = std::get<0>(test_param);
    const std::string &base_name = std::get<1>(test_param);

    return std::to_string(case_id) + "_" + base_name + "_" + std::to_string(info.index);
  }

 protected:
  void SetUp() override {
    const auto &param_tuple = GetParam();
    const TestType &test_case = std::get<2>(param_tuple);  // третий элемент
    int case_id = std::get<0>(test_case);

    switch (case_id) {
      case 0:  // 1x1
        input_data_ = {1, 1, {5}, {3}};
        expected_ = {15};
        break;
      case 1:  // 1x4
        input_data_ = {1, 4, {1, 2, 3, 4}, {1, 1, 1, 1}};
        expected_ = {10};
        break;
      case 2:  // 4x1
        input_data_ = {4, 1, {1, 2, 3, 4}, {2}};
        expected_ = {2, 4, 6, 8};
        break;
      case 3:  // 3x3
        input_data_ = {3, 3, {1, 2, 3, 4, 5, 6, 7, 8, 9}, {1, 2, 3}};
        expected_ = {14, 32, 50};
        break;
      case 4:  // zeros
        input_data_ = {3, 3, std::vector<int>(9, 0), std::vector<int>(3, 0)};
        expected_ = {0, 0, 0};
        break;
      default:
        throw std::runtime_error("Unknown test case");
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank != 0) {
      return true;
    }
    return output_data == expected_;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_;
};

TEST_P(KulikovMatrixMultiplyRunFuncTests, RunFunctionalTest) {
  ExecuteTest(GetParam());
}

namespace {

const std::array<TestType, 5> kTestParams = {
    {std::make_tuple(0, std::string("Single")), std::make_tuple(1, std::string("SingleRow")),
     std::make_tuple(2, std::string("SingleCol")), std::make_tuple(3, std::string("Square")),
     std::make_tuple(4, std::string("Zeros"))}};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<KulikovDMatrixMultiplyMPI, InType>(
                                               kTestParams, PPC_SETTINGS_kulikov_d_matrix_vector_multiply),
                                           ppc::util::AddFuncTask<KulikovDMatrixMultiplySEQ, InType>(
                                               kTestParams, PPC_SETTINGS_kulikov_d_matrix_vector_multiply));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

INSTANTIATE_TEST_SUITE_P(MatrixVectorMultiply, KulikovMatrixMultiplyRunFuncTests, kGtestValues,
                         KulikovMatrixMultiplyRunFuncTests::PrintFuncTestName);

}  // namespace

}  // namespace kulikov_d_matrix_vector_multiply
