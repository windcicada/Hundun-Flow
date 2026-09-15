// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "cantera/base/Solution.h"
#include "cantera/kinetics/Arrhenius.h"
#include "cantera/kinetics/KineticsFactory.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/thermo/ThermoPhase.h"
#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

namespace hundun::v04::chemistry::detail {

// Every chemistry query and reactor RHS reads this lane's current concentrations.
// Storage is prepared with the kinetics evaluator, outside time integration.
struct OrderData final : Cantera::ReactionData {
  std::vector<double> concentrations;
  using Cantera::ReactionData::update;
  void resize(std::size_t species, std::size_t, std::size_t) override {
    concentrations.resize(species);
  }
  bool update(const Cantera::ThermoPhase &phase,
              const Cantera::Kinetics &) override {
    Cantera::ReactionData::update(phase.temperature());
    phase.getConcentrations(concentrations.data());
    return true;
  }
};

struct OrderTerm final {
  std::size_t species;
  double exponent;
};

class OrderRate final : public Cantera::ReactionRate {
public:
  OrderRate(const Cantera::ArrheniusRate &rate, std::vector<OrderTerm> terms,
            double negative_floor, double fractional_floor)
      : arrhenius_(rate), terms_(std::move(terms)),
        negative_floor_(negative_floor), fractional_floor_(fractional_floor) {
    m_valid = true;
    setCompositionDependence(true);
  }
  const std::string type() const override {
    return "hundun-concentration-c1-v1";
  }
  std::unique_ptr<Cantera::MultiRateBase> newMultiRate() const override {
    return std::make_unique<Cantera::MultiRate<OrderRate, OrderData>>();
  }
  double evalFromStruct(const OrderData &data) const {
    double rate = arrhenius_.evalRate(data.logT, data.recipT);
    for (const auto &term : terms_) {
      const double c = data.concentrations[term.species];
      const double exponent = term.exponent;
      double factor;
      if (exponent < 0) {
        factor = std::pow(std::max(c, negative_floor_), exponent);
      } else if (c < fractional_floor_) {
        // Match value and slope of c^exponent at the positive threshold.
        // Signed continuation supplies a restoring derivative during Newton
        // trials; accepted compositions still pass the common state checks.
        factor = c * std::pow(fractional_floor_, exponent - 1) *
            (2 - exponent + (exponent - 1) *
                                  std::max(c, 0.) / fractional_floor_);
      } else {
        factor = std::pow(c, exponent);
      }
      rate *= factor;
    }
    return rate;
  }

private:
  Cantera::ArrheniusRate arrhenius_;
  std::vector<OrderTerm> terms_;
  double negative_floor_, fractional_floor_;
};

// The opt-in law lives in the hashed mechanism asset. It applies to explicitly
// tagged irreversible Arrhenius reactions and keeps original SI rate constants.
// Selected concentration powers move into the rate evaluator; their external
// mass-action orders become zero, avoiding 0 * infinity and double application.
inline void install_order_rates(const std::shared_ptr<Cantera::Solution> &solution) {
  constexpr const char *key = "hundun-order-regularization";
  const auto original = solution->kinetics();
  bool enabled = false;
  for (std::size_t i = 0; i < original->nReactions(); ++i)
    enabled = enabled || original->reaction(i)->input.hasKey(key);
  if (!enabled)
    return;
  if (original->kineticsType() != "bulk" || original->nPhases() != 1 ||
      solution->thermo()->type() != "ideal-gas")
    throw std::invalid_argument("order regularization requires ideal-gas bulk kinetics");

  auto phase = solution->thermo()->input();
  phase["kinetics"] = "gas";
  phase["reactions"] = "none";
  auto kinetics = Cantera::newKinetics({solution->thermo()}, phase);
  for (std::size_t i = 0; i < original->nReactions(); ++i) {
    auto reaction = std::make_shared<Cantera::Reaction>(*original->reaction(i));
    if (reaction->input.hasKey(key)) {
      const auto &config = reaction->input.at(key).as<Cantera::AnyMap>();
      for (const auto &entry : config)
        if (entry.first != "model" && entry.first != "negative-floor" &&
            entry.first != "fractional-floor")
          throw std::invalid_argument("unknown order regularization parameter");
      if (config.at("model").asString() != "concentration_c1_v1")
        throw std::invalid_argument("unsupported order regularization model");
      const double negative = config.convert("negative-floor", "kmol/m^3");
      const double fractional = config.convert("fractional-floor", "kmol/m^3");
      if (!(negative > 0) || !std::isfinite(negative) ||
          !(fractional > 0) || !std::isfinite(fractional))
        throw std::invalid_argument("invalid order regularization concentration");
      const auto rate = std::dynamic_pointer_cast<Cantera::ArrheniusRate>(reaction->rate());
      if (!rate || reaction->rate()->type() != "Arrhenius" ||
          reaction->reversible || reaction->thirdBody())
        throw std::invalid_argument("order regularization requires irreversible Arrhenius rates");
      auto orders = reaction->reactants;
      for (const auto &entry : reaction->orders)
        orders[entry.first] = entry.second;
      std::vector<OrderTerm> terms;
      for (const auto &entry : orders) {
        if (!std::isfinite(entry.second))
          throw std::invalid_argument("invalid reaction order");
        if (entry.second < 0 || (entry.second > 0 && entry.second < 1)) {
          const auto species = solution->thermo()->speciesIndex(entry.first);
          if (species == Cantera::npos)
            throw std::invalid_argument("unknown regularized species");
          terms.push_back({species, entry.second});
          reaction->orders[entry.first] = 0;
        }
      }
      if (terms.empty())
        throw std::invalid_argument("order regularization requires a fractional or negative order");
      reaction->setRate(std::make_shared<OrderRate>(*rate, std::move(terms),
                                                   negative, fractional));
    }
    if (!kinetics->addReaction(reaction))
      throw std::invalid_argument("regularized kinetics reaction admission failed");
  }
  kinetics->checkDuplicates();
  solution->setKinetics(kinetics);
}

} // namespace hundun::v04::chemistry::detail
