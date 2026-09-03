// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/lin_system.hpp"

#include <cstdint>
#include <memory>

namespace hundun::linear {

#ifdef HUNDUN_LINEAR_ENABLE_TEST_ACCESS
namespace test {
class PreconditionerTestAccess;
}
#endif

class IdentityPreconditioner final : public Preconditioner {
 public:
  explicit IdentityPreconditioner(
      execution::ExecutionContext& context);
  ~IdentityPreconditioner() noexcept override;
  IdentityPreconditioner(const IdentityPreconditioner&) = delete;
  IdentityPreconditioner& operator=(const IdentityPreconditioner&) = delete;

  void update(const LinearOperator&, std::uint64_t revision) override;
  execution::ExecutionEvent apply(
      execution::VectorView<const double> r,
      execution::VectorView<double> z) const override;

 private:
  struct State;
  std::unique_ptr<State> state_;
};

class JacobiPreconditioner final : public Preconditioner {
 public:
  explicit JacobiPreconditioner(
      execution::ExecutionContext& context);
  ~JacobiPreconditioner() noexcept override;
  JacobiPreconditioner(const JacobiPreconditioner&) = delete;
  JacobiPreconditioner& operator=(const JacobiPreconditioner&) = delete;

  void update(const LinearOperator&, std::uint64_t revision) override;
  execution::ExecutionEvent apply(
      execution::VectorView<const double> r,
      execution::VectorView<double> z) const override;

 private:
  struct State;
  std::unique_ptr<State> state_;

#ifdef HUNDUN_LINEAR_ENABLE_TEST_ACCESS
  friend class test::PreconditionerTestAccess;
#endif
};

}  // namespace hundun::linear
