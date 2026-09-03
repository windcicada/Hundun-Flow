// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/exec_execution.hpp"
#include "hundun/rt_mpi_context.hpp"

#include <cstddef>

namespace hundun::linear {

struct DotProductPair final {
  execution::VectorView<const double> left;
  execution::VectorView<const double> right;
};

class VectorOps final {
 public:
  explicit VectorOps(const execution::ExecutionContext& context);

  void fill(execution::VectorView<double> output, double value) const;
  void copy(execution::VectorView<const double> input,
            execution::VectorView<double> output) const;
  void scale(double alpha, execution::VectorView<double> values) const;
  void axpy(double alpha, execution::VectorView<const double> x,
            execution::VectorView<double> y) const;
  void linear_combination(
      double alpha, execution::VectorView<const double> x,
      double beta, execution::VectorView<const double> y,
      execution::VectorView<double> output) const;

  double norm(execution::VectorView<const double> values,
              const runtime::MpiContext& context) const;

  void dot_batch(const DotProductPair* pairs, std::size_t pair_count,
                 execution::VectorView<double> results,
                 const runtime::MpiContext& context) const;

 private:
  execution::BackendIdentity backend_identity_{};
  execution::ExecutionSpace space_{execution::ExecutionSpace::host};
};

}  // namespace hundun::linear
