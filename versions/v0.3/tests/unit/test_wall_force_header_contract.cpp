// SPDX-License-Identifier: Apache-2.0

#include "hundun/ib_wall_force.hpp"

#include "hundun/flow_immersed.hpp"

#include <type_traits>

namespace {

using hundun::immersed::WallForceIntegrator;
using hundun::immersed::WallForceSample;

using Create =
    WallForceIntegrator (*)(const hundun::immersed::WallQuadraturePlan &,
                            const hundun::runtime::MpiContext &);
using IntegrateUnconstrained = WallForceSample (WallForceIntegrator::*)(
    const hundun::runtime::FieldView<const double> &,
    const hundun::runtime::FieldView<const double> &,
    const hundun::runtime::FieldView<const double> &,
    const hundun::runtime::FieldView<const double> &) const;

static_assert(std::is_final_v<hundun::immersed::ForceComponents>);
static_assert(std::is_final_v<hundun::immersed::MomentComponents>);
static_assert(std::is_final_v<WallForceSample>);
static_assert(std::is_final_v<WallForceIntegrator>);
static_assert(std::is_same_v<decltype(&WallForceIntegrator::create), Create>);
static_assert(std::is_same_v<decltype(&WallForceIntegrator::integrate),
                             IntegrateUnconstrained>);
static_assert(
    std::is_same_v<decltype(hundun::flow::ImmersedFlowStepAttemptReport{}.force),
                   std::optional<hundun::flow::ForceAttemptReport>>);

} // namespace

int main() {
  hundun::immersed::ForceComponents force{};
  hundun::immersed::MomentComponents moment{};
  WallForceSample sample{};
  return force.pressure_N.x == 0.0 && moment.pressure_N_m.x == 0.0 &&
                 sample.lowest_failing_rank == -1
             ? 0
             : 1;
}
