// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/rt_collective_status.hpp"

#include <cstdint>

namespace hundun::flow::detail {

enum class ReactingPressureDomain : std::uint8_t {
  open,
  closed,
  partially_closed
};

class ReactingPressureState final {
public:
  ReactingPressureState(double history_p0_pa, double committed_p0_pa);
  double history_p0_pa() const noexcept;
  double committed_p0_pa() const noexcept;
  double trial_p0_pa() const noexcept;
  bool attempt_active() const noexcept;
  void begin_attempt();
  void set_trial_p0_pa(double);
  bool commit(const runtime::CollectiveStatus &);
  void rollback();

private:
  double history_p0_pa_{};
  double committed_p0_pa_{};
  double trial_p0_pa_{};
  bool attempt_active_{};
};

struct ReactingPressureConstraintRequest final {
  ReactingPressureDomain domain{ReactingPressureDomain::open};
  double duration_s{};
  double local_pressure_capacity_m3_per_pa{};
  double local_predictor_expansion_m3{};
  double local_corrector_expansion_m3{};
  double local_boundary_volume_outflow_m3{};
};

struct ReactingPressureConstraintReport final {
  ReactingPressureDomain domain{ReactingPressureDomain::open};
  double previous_p0_pa{};
  double predictor_p0_pa{};
  double corrected_p0_pa{};
  double midpoint_dp0_dt_pa_per_s{};
  double integrated_expansion_m3{};
  double integrated_boundary_outflow_m3{};
  std::uint32_t global_reduction_count{};
};

ReactingPressureConstraintReport attempt_reacting_pressure_constraint(
    ReactingPressureState &, const runtime::MpiContext &,
    const ReactingPressureConstraintRequest &);

} // namespace hundun::flow::detail
