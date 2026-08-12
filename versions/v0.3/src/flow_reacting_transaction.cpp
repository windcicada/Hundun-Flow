// SPDX-License-Identifier: Apache-2.0

#include "flow_reacting_transaction_detail.hpp"

#include "hundun/rt_error.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace hundun::flow::detail {
namespace {

void swap_layer(ReactingLayerState &left,
                ReactingLayerState &right) noexcept {
  left.rho_y_kg_per_m3.swap(right.rho_y_kg_per_m3);
  left.rho_h_tc_j_per_m3.swap(right.rho_h_tc_j_per_m3);
  std::swap(left.p0_pa, right.p0_pa);
}

void validate_layer(const ReactingLayerState &layer, std::size_t cell_count,
                    std::size_t species_count) {
  if (species_count == 0U ||
      cell_count > std::numeric_limits<std::size_t>::max() / species_count ||
      layer.rho_y_kg_per_m3.size() != cell_count * species_count ||
      layer.rho_h_tc_j_per_m3.size() != cell_count ||
      !std::isfinite(layer.p0_pa) || layer.p0_pa <= 0.0) {
    throw runtime::Error("invalid reacting attempt-state layer");
  }
  const auto finite = [](double value) { return std::isfinite(value); };
  if (!std::all_of(layer.rho_y_kg_per_m3.begin(),
                   layer.rho_y_kg_per_m3.end(), finite) ||
      !std::all_of(layer.rho_h_tc_j_per_m3.begin(),
                   layer.rho_h_tc_j_per_m3.end(), finite)) {
    throw runtime::Error("reacting attempt-state layer is non-finite");
  }
}

} // namespace

bool operator==(const ReactingLayerState &left,
                const ReactingLayerState &right) noexcept {
  return left.rho_y_kg_per_m3 == right.rho_y_kg_per_m3 &&
         left.rho_h_tc_j_per_m3 == right.rho_h_tc_j_per_m3 &&
         left.p0_pa == right.p0_pa;
}

bool operator!=(const ReactingLayerState &left,
                const ReactingLayerState &right) noexcept {
  return !(left == right);
}

ReactingAttemptState::ReactingAttemptState(std::size_t cell_count,
                                           std::size_t species_count,
                                           ReactingLayerState history,
                                           ReactingLayerState committed)
    : cell_count_(cell_count), species_count_(species_count),
      history_(std::move(history)), committed_(std::move(committed)),
      trial_(committed_) {
  if (cell_count_ == 0U) {
    throw runtime::Error("reacting attempt state requires cells");
  }
  validate_layer(history_, cell_count_, species_count_);
  validate_layer(committed_, cell_count_, species_count_);
}

std::size_t ReactingAttemptState::cell_count() const noexcept {
  return cell_count_;
}

std::size_t ReactingAttemptState::species_count() const noexcept {
  return species_count_;
}

const ReactingLayerState &ReactingAttemptState::history() const noexcept {
  return history_;
}

const ReactingLayerState &ReactingAttemptState::committed() const noexcept {
  return committed_;
}

bool ReactingAttemptState::attempt_active() const noexcept {
  return attempt_active_;
}

void ReactingAttemptState::begin_attempt() {
  if (attempt_active_) {
    throw runtime::Error("reacting attempt is already active");
  }
  trial_ = committed_;
  attempt_active_ = true;
}

void ReactingAttemptState::rollback_attempt() {
  if (!attempt_active_) {
    throw runtime::Error("reacting rollback requires an active attempt");
  }
  trial_ = committed_;
  attempt_active_ = false;
}

void ReactingAttemptState::publish_attempt() {
  ReactingLayerState next_history = committed_;
  ReactingLayerState next_committed = trial_;
  ReactingLayerState next_trial = trial_;
  swap_layer(history_, next_history);
  swap_layer(committed_, next_committed);
  swap_layer(trial_, next_trial);
  attempt_active_ = false;
}

ReactingSourceTransaction::ReactingSourceTransaction(
    ReactingAttemptState &state)
    : state_(&state),
      species_delta_(state.cell_count() * state.species_count(), 0.0),
      enthalpy_delta_(state.cell_count(), 0.0),
      chemistry_mass_delta_(state.cell_count(), 0.0) {
  if (!state.attempt_active()) {
    throw runtime::Error(
        "reacting source transaction requires an active attempt");
  }
}

void ReactingSourceTransaction::require_open() const {
  if (closed_ || state_ == nullptr || !state_->attempt_active()) {
    throw runtime::Error("reacting source transaction is closed");
  }
}

void ReactingSourceTransaction::add_species(
    const ReactingSourceIdentity &source, std::size_t cell,
    std::size_t species, double delta_kg_per_m3) {
  require_open();
  if (source.id.empty() || cell >= state_->cell_count() ||
      species >= state_->species_count() || !std::isfinite(delta_kg_per_m3)) {
    throw runtime::Error("invalid reacting species source delta");
  }
  const std::size_t index = cell * state_->species_count() + species;
  species_delta_[index] += delta_kg_per_m3;
  if (source.kind == ReactingSourceKind::chemistry) {
    chemistry_mass_delta_[cell] += delta_kg_per_m3;
  }
  records_.push_back({source, ReactingSourceQuantity::species_density, cell,
                      species, delta_kg_per_m3, "kg/m^3"});
}

void ReactingSourceTransaction::add_enthalpy(
    const ReactingSourceIdentity &source, std::size_t cell,
    double delta_j_per_m3) {
  require_open();
  if (source.id.empty() || source.kind == ReactingSourceKind::chemistry ||
      cell >= state_->cell_count() || !std::isfinite(delta_j_per_m3)) {
    throw runtime::Error("invalid reacting enthalpy source delta");
  }
  enthalpy_delta_[cell] += delta_j_per_m3;
  records_.push_back({source, ReactingSourceQuantity::enthalpy_density, cell,
                      0U, delta_j_per_m3, "J/m^3"});
}

const std::vector<double> &
ReactingSourceTransaction::species_delta_kg_per_m3() const noexcept {
  return species_delta_;
}

const std::vector<double> &
ReactingSourceTransaction::enthalpy_delta_j_per_m3() const noexcept {
  return enthalpy_delta_;
}

const std::vector<ReactingSourceRecord> &
ReactingSourceTransaction::records() const noexcept {
  return records_;
}

void ReactingSourceTransaction::rollback_noexcept() noexcept {
  if (state_ != nullptr && state_->attempt_active()) {
    std::copy(state_->committed_.rho_y_kg_per_m3.begin(),
              state_->committed_.rho_y_kg_per_m3.end(),
              state_->trial_.rho_y_kg_per_m3.begin());
    std::copy(state_->committed_.rho_h_tc_j_per_m3.begin(),
              state_->committed_.rho_h_tc_j_per_m3.end(),
              state_->trial_.rho_h_tc_j_per_m3.begin());
    state_->trial_.p0_pa = state_->committed_.p0_pa;
    state_->attempt_active_ = false;
  }
  closed_ = true;
}

bool ReactingSourceTransaction::commit(
    const runtime::CollectiveStatus &status) {
  require_open();
  if (!status.ok) {
    rollback_noexcept();
    return false;
  }
  try {
    for (const double chemistry_mass : chemistry_mass_delta_) {
      if (!std::isfinite(chemistry_mass) ||
          std::abs(chemistry_mass) > 1.0e-12) {
        throw runtime::Error("chemistry species source does not conserve mass");
      }
    }
    for (std::size_t index = 0; index < species_delta_.size(); ++index) {
      state_->trial_.rho_y_kg_per_m3[index] += species_delta_[index];
      if (!std::isfinite(state_->trial_.rho_y_kg_per_m3[index]) ||
          state_->trial_.rho_y_kg_per_m3[index] < 0.0) {
        throw runtime::Error("reacting species source produced invalid state");
      }
    }
    for (std::size_t cell = 0; cell < enthalpy_delta_.size(); ++cell) {
      state_->trial_.rho_h_tc_j_per_m3[cell] += enthalpy_delta_[cell];
      if (!std::isfinite(state_->trial_.rho_h_tc_j_per_m3[cell])) {
        throw runtime::Error("reacting enthalpy source produced invalid state");
      }
    }
    state_->publish_attempt();
    closed_ = true;
    return true;
  } catch (...) {
    rollback_noexcept();
    throw;
  }
}

void ReactingSourceTransaction::rollback() {
  require_open();
  state_->rollback_attempt();
  closed_ = true;
}

} // namespace hundun::flow::detail
