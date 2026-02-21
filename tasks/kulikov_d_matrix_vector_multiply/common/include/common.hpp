#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace kulikov_d_matrix_vector_multiply {

using ValueType = int;

struct MatrixType {
  int rows;
  int cols;
  std::vector<ValueType> matrix;
  std::vector<ValueType> vector;
};

using InType = MatrixType;
using OutType = std::vector<ValueType>;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace kulikov_d_matrix_vector_multiply
