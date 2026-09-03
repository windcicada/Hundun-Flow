// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "flow_reacting_coupling_detail.hpp"

#include "hundun/rt_error.hpp"

#include <cmath>
#include <stdexcept>

namespace hundun::flow {
namespace {

void require_update_shape(const ReactingOperatorUpdate &update,
                          std::size_t cells, std::size_t species,
                          bool chemistry) {
  if (!update.ok)
    return;
  if ((!update.species_delta_kg_per_m3.empty() &&
       update.species_delta_kg_per_m3.size() != cells * species) ||
      (!update.enthalpy_delta_j_per_m3.empty() &&
       update.enthalpy_delta_j_per_m3.size() != cells))
    throw runtime::Error("reacting operator update shape is invalid");
  for (double value : update.species_delta_kg_per_m3)
    if (!std::isfinite(value))
      throw runtime::Error("reacting operator species update is non-finite");
  for (double value : update.enthalpy_delta_j_per_m3)
    if (!std::isfinite(value) || (chemistry && value != 0.0))
      throw runtime::Error("reacting operator enthalpy update is invalid");
  if (chemistry && !update.species_delta_kg_per_m3.empty()) {
    for (std::size_t cell = 0; cell < cells; ++cell) {
      double mass{};
      for (std::size_t k = 0; k < species; ++k)
        mass += update.species_delta_kg_per_m3[cell * species + k];
      if (std::abs(mass) > 1.0e-12)
        throw runtime::Error("reacting chemistry update does not conserve mass");
    }
  }
}

void apply_update(detail::ReactingSourceTransaction &transaction,
                  detail::ReactingLayerState &working,
                  const ReactingOperatorUpdate &update,
                  const detail::ReactingSourceIdentity &source,
                  std::size_t cells, std::size_t species) {
  if (!update.species_delta_kg_per_m3.empty()) {
    for (std::size_t cell = 0; cell < cells; ++cell) {
      for (std::size_t k = 0; k < species; ++k) {
        const std::size_t index = cell * species + k;
        const double delta = update.species_delta_kg_per_m3[index];
        transaction.add_species(source, cell, k, delta);
        working.rho_y_kg_per_m3[index] += delta;
        if (!std::isfinite(working.rho_y_kg_per_m3[index]) ||
            working.rho_y_kg_per_m3[index] < 0.0)
          throw runtime::Error("reacting operator produced invalid species state");
      }
    }
  }
  if (!update.enthalpy_delta_j_per_m3.empty()) {
    for (std::size_t cell = 0; cell < cells; ++cell) {
      const double delta = update.enthalpy_delta_j_per_m3[cell];
      transaction.add_enthalpy(source, cell, delta);
      working.rho_h_tc_j_per_m3[cell] += delta;
      if (!std::isfinite(working.rho_h_tc_j_per_m3[cell]))
        throw runtime::Error("reacting operator produced invalid enthalpy state");
    }
  }
}

} // namespace

ReactingStepReport attempt_open_reacting_step(
    detail::ReactingAttemptState &state, double start_time_s,
    double duration_s, const OpenReactingStepOperators &operators) {
  if (!std::isfinite(start_time_s) || !std::isfinite(duration_s) ||
      duration_s <= 0.0 || !operators.chemistry ||
      !operators.scalar_transport || !operators.pressure_corrector ||
      !operators.collective_validation)
    throw std::invalid_argument("open reacting step contract is invalid");

  ReactingStepReport report;
  report.p0_before_pa = state.committed().p0_pa;
  const std::size_t cells = state.cell_count();
  const std::size_t species = state.species_count();
  state.begin_attempt();
  detail::ReactingSourceTransaction transaction(state);
  detail::ReactingLayerState working = state.committed();

  const auto reject = [&](std::string stage) {
    report.failure_stage = std::move(stage);
    report.p0_after_pa = state.committed().p0_pa;
    transaction.rollback();
    return report;
  };
  const auto collectively_accept = [&](bool local_ok,
                                       const std::string &message) {
    const auto status = operators.collective_validation(local_ok, message);
    return local_ok && status.ok;
  };

  try {
    const double half = 0.5 * duration_s;
    auto chemistry1 = operators.chemistry(1U, start_time_s, half, working);
    ++report.chemistry_call_count;
    if (!collectively_accept(chemistry1.ok, chemistry1.message))
      return reject("chemistry-1");
    require_update_shape(chemistry1, cells, species, true);
    apply_update(transaction, working, chemistry1,
                 {"chemistry-half-1", detail::ReactingSourceKind::chemistry},
                 cells, species);

    auto transport =
        operators.scalar_transport(start_time_s + half, duration_s, working);
    ++report.scalar_transport_count;
    if (!collectively_accept(transport.ok, transport.message))
      return reject("scalar-transport");
    require_update_shape(transport, cells, species, false);
    apply_update(transaction, working, transport,
                 {"scalar-transport", detail::ReactingSourceKind::transport},
                 cells, species);

    std::string piso_message;
    const bool piso1 = operators.pressure_corrector(
        1U, start_time_s + half, 1U, working,
        transaction.species_delta_kg_per_m3(),
        transaction.enthalpy_delta_j_per_m3(), piso_message);
    ++report.pressure_corrector_count;
    if (!collectively_accept(piso1, piso_message))
      return reject("piso-1");

    auto chemistry2 =
        operators.chemistry(2U, start_time_s + half, half, working);
    ++report.chemistry_call_count;
    if (!collectively_accept(chemistry2.ok, chemistry2.message))
      return reject("chemistry-2");
    require_update_shape(chemistry2, cells, species, true);
    apply_update(transaction, working, chemistry2,
                 {"chemistry-half-2", detail::ReactingSourceKind::chemistry},
                 cells, species);
    report.post_chemistry2_epoch = 2U;

    piso_message.clear();
    const bool piso2 = operators.pressure_corrector(
        2U, start_time_s + duration_s, report.post_chemistry2_epoch, working,
        transaction.species_delta_kg_per_m3(),
        transaction.enthalpy_delta_j_per_m3(), piso_message);
    ++report.pressure_corrector_count;
    if (!collectively_accept(piso2, piso_message))
      return reject("piso-2");
    report.piso2_consumed_epoch = report.post_chemistry2_epoch;
    report.predictor_to_final_delta_flux_applied = true;

    report.integrated_species_delta_kg_per_m3 =
        transaction.species_delta_kg_per_m3();
    report.integrated_enthalpy_delta_j_per_m3 =
        transaction.enthalpy_delta_j_per_m3();
    if (working.p0_pa != report.p0_before_pa)
      return reject("p0-authority");
    const auto final_status = operators.collective_validation(true, {});
    if (!final_status.ok)
      return reject("final-validation");
    report.accepted = transaction.commit(final_status);
    report.p0_after_pa = state.committed().p0_pa;
    if (!report.accepted || report.p0_after_pa != report.p0_before_pa)
      throw runtime::Error("open reacting step violated p0 authority");
    return report;
  } catch (...) {
    if (state.attempt_active())
      transaction.rollback();
    throw;
  }
}

} // namespace hundun::flow
