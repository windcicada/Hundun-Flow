// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#pragma once
#include "models_kerosene_detail.hpp"
#include "cantera/base/Solution.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/numerics/FuncEval.h"
#include "cantera/numerics/Integrator.h"
#include "cantera/thermo/ThermoPhase.h"
#include <memory>
#include <string>
#include <vector>

namespace hundun::v04::chemistry::detail {
// One exclusive chemistry lane owns its fixed-material interval and CVODE
// state. Public coordinates remain PH/Y; the interval integrates kmol/kg.
class KeroseneMaterial final : public Cantera::FuncEval {
public:
  static constexpr const char* model = "kerosene_4_step_v1/frozen_material";
  explicit KeroseneMaterial(const Cantera::ThermoPhase& phase)
      : integrator_(Cantera::newIntegrator("CVODE")) {
    if (phase.type() != "ideal-gas" || phase.nSpecies() != 7 || phase.nElements() != 4)
      throw std::invalid_argument("kerosene model requires its seven ideal-gas species");
    const char* names[]{"H2", "H2O", "CO", "CO2", "O2", "N2", "C12H23"};
    const char* elements[]{"O", "H", "C", "N"};
    constexpr double atoms[4][7]{{0,1,1,2,2,0,0}, {2,2,0,0,0,0,23},
                                {0,0,1,1,0,0,12}, {0,0,0,0,0,2,0}};
    for (unsigned i = 0; i < 7; ++i) {
      indices_[i] = phase.speciesIndex(names[i]);
      if (indices_[i] == Cantera::npos)
        throw std::invalid_argument("kerosene model species identity mismatch");
      weights_[i] = phase.molecularWeight(indices_[i]);
      if (!(weights_[i] > 0) || !std::isfinite(weights_[i]))
        throw std::invalid_argument("kerosene model molecular weight mismatch");
      for (unsigned e = 0; e < 4; ++e) {
        const auto index = phase.elementIndex(elements[e]);
        if (index == Cantera::npos || phase.nAtoms(indices_[i], index) != atoms[e][i])
          throw std::invalid_argument("kerosene model elemental identity mismatch");
      }
    }
    integrator_->setMethod(Cantera::BDF_Method);
    integrator_->setLinearSolverType("DENSE");
    suppressErrors(true);
  }

  // Query preparation uses a local kernel, leaving an interval's frozen
  // coefficients intact. Outputs are kmol/(m3 s), in the asset's order.
  void molar_rates(const Cantera::ThermoPhase& phase, double* output) const {
    KeroseneFourStep kernel;
    const auto initial = prepare_kernel(phase, kernel);
    KeroseneFourStep::Tuple rates;
    const double density = phase.density();
    if (!kernel.evaluate(density, initial, rates))
      throw std::runtime_error("kerosene instantaneous source admission");
    for (unsigned i = 0; i < 7; ++i) output[indices_[i]] = density * rates[i];
  }

  void integrate(const Cantera::ThermoPhase& phase, double duration,
                 double relative, double absolute, int max_steps,
                 std::vector<double>& output) {
    steps_ = 0;
    if (output.size() != 7) throw std::invalid_argument("kerosene output capacity");
    initial_ = prepare_kernel(phase, kernel_);
    density_ = phase.density();
    std::array<double, 7> tolerances;
    for (unsigned i = 0; i < 7; ++i) tolerances[i] = absolute / weights_[i];
    integrator_->setTolerances(relative, 7, tolerances.data());
    integrator_->setMaxSteps(max_steps);
    if (initialized_) integrator_->reinitialize(0., *this);
    else {
      integrator_->initialize(0., *this);
      initialized_ = true;
    }
    try {
      integrator_->integrate(duration);
    } catch (...) {
      record_steps();
      throw;
    }
    record_steps();
    std::array<double, 7> candidate;
    for (unsigned i = 0; i < 7; ++i) {
      candidate[indices_[i]] = integrator_->solution(i) * weights_[i];
      if (!std::isfinite(candidate[indices_[i]]))
        throw std::runtime_error("kerosene interval finite-state admission");
    }
    std::copy(candidate.begin(), candidate.end(), output.begin());
  }

  long steps() const noexcept { return steps_; }
  std::size_t neq() const override { return 7; }
  void getState(double* y) override { std::copy(initial_.begin(), initial_.end(), y); }
  void eval(double, double* y, double* ydot, double*) override {
    KeroseneFourStep::Tuple trial, rates;
    std::copy_n(y, 7, trial.begin());
    if (!kernel_.evaluate(density_, trial, rates))
      throw Cantera::CanteraError("KeroseneMaterial::eval", "RHS admission failed");
    std::copy(rates.begin(), rates.end(), ydot);
  }

private:
  KeroseneFourStep::Tuple prepare_kernel(const Cantera::ThermoPhase& phase,
                                         KeroseneFourStep& kernel) const {
    std::array<double, 7> fractions, reference_gibbs;
    phase.getMassFractions(fractions.data());
    // Reference Gibbs values exclude the ideal-gas log(p/p_ref) contribution.
    // The four-step equilibrium expression supplies its own standard factor.
    phase.getGibbs_RT_ref(reference_gibbs.data());
    KeroseneFourStep::Tuple moles, gibbs;
    for (unsigned i = 0; i < 7; ++i) {
      moles[i] = fractions[indices_[i]] / weights_[i];
      gibbs[i] = reference_gibbs[indices_[i]];
    }
    if (!kernel.prepare(phase.temperature(), Cantera::GasConstant, weights_[6], moles, gibbs))
      throw std::runtime_error("kerosene coefficient admission");
    return moles;
  }
  void record_steps() noexcept {
    try { steps_ = std::max(0L, integrator_->solverStats()["steps"].asInt()); }
    catch (...) { steps_ = -1; }
  }
  std::array<std::size_t, 7> indices_{};
  std::array<double, 7> weights_{};
  KeroseneFourStep kernel_;
  KeroseneFourStep::Tuple initial_{};
  double density_{};
  long steps_{};
  bool initialized_{};
  std::unique_ptr<Cantera::Integrator> integrator_;
};

inline std::unique_ptr<KeroseneMaterial> kerosene_material(
    const std::shared_ptr<Cantera::Solution>& solution) {
  const auto& phase = *solution->thermo();
  const auto& input = phase.input();
  if (!input.hasKey("hundun-chemistry")) return {};
  const auto& spec = input.at("hundun-chemistry").as<Cantera::AnyMap>();
  for (const auto& entry : spec)
    if (entry.first != "model" && entry.first != "interval")
      throw std::invalid_argument("unknown kerosene chemistry model parameter");
  if (spec.at("model").asString() != "kerosene_4_step_v1" ||
      spec.at("interval").asString() != "frozen_material" ||
      solution->kinetics()->nReactions() != 0)
    throw std::invalid_argument("kerosene chemistry model identity mismatch");
  return std::make_unique<KeroseneMaterial>(phase);
}
} // namespace hundun::v04::chemistry::detail
