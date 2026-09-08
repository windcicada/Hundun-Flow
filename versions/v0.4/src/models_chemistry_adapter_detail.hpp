// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "hundun/v04_chemistry.hpp"
#include "hundun/v04_combustion.hpp"
#include "hundun/v04_portable.hpp"

namespace hundun::v04::chemistry::detail {

// Preparation validates both fingerprints and the explicit species/element
// mapping. The provider and backend must outlive this exclusive-lane adapter.
class BackendAdapter final : public combustion::ChemistryRepresentationAdapter,
                             public combustion::ChemicalRateQueryProvider {
public:
  BackendAdapter(ChemistryBackend &, portable::GasQueryProvider &,
                 combustion::ChemistryIdentity, portable::Revision);
  const combustion::ChemistryIdentity &identity() const noexcept override;
  combustion::ChemistryRepresentationKind
  representation() const noexcept override;
  combustion::ChemistryAdvanceReport
  integrate(const combustion::ChemistryAdvanceRequest &) noexcept override;
  combustion::ChemicalRateQueryReport
  query(const combustion::ChemicalRateQuery &) noexcept override;

private:
  ChemistryBackend &backend_;
  portable::GasQueryProvider &gas_;
  combustion::ChemistryIdentity identity_;
  portable::Revision revision_;
  std::vector<double> diffusion_, enthalpies_, rates_;
};

struct GasBatchReport {
  portable::Status status{portable::Status::invalid_input};
  std::size_t failure_index{};
  std::size_t count{};
  bool available{};
};
// Preparation allocates fixed capacities. query() allocates no HUNDUN memory;
// provider/third-party internal allocation is outside that guarantee.
// Borrowed result pointers remain valid until the next query or destruction.
class GasBatchWorkspace final {
public:
  GasBatchWorkspace(std::size_t states, std::size_t species);
  GasBatchWorkspace(const GasBatchWorkspace &) = delete;
  GasBatchWorkspace &operator=(const GasBatchWorkspace &) = delete;
  GasBatchReport query(portable::GasQueryProvider &, const portable::GasQuery *,
                       std::size_t, portable::Revision) noexcept;
  const portable::GasQueryOutput *results() const noexcept;

private:
  std::size_t species_{};
  bool available_{};
  std::vector<double> storage_;
  std::vector<portable::GasQueryOutput> outputs_;
};

// Independent synthetic A -> B constant-cp ideal-gas backend. Species have
// equal molecular weights/elements; k is first-order s^-1, no real fuel claim.
class AnalyticIsomerBackend final : public ChemistryBackend,
                                    public portable::GasQueryProvider,
                                    public portable::GasAdvanceProvider {
public:
  AnalyticIsomerBackend(double rate_s = 2.0, bool reversed_order = false,
                        double synthetic_cp_j_per_kg_k = 1000.0);
  const CompositionIdentity &composition() const noexcept override;
  const portable::GasIdentity &gas_identity() const noexcept override;
  const combustion::ChemistryIdentity &closure_identity() const noexcept;
  portable::Status query_gas(const portable::GasQuery &,
                             portable::GasQueryOutput &) noexcept override;
  ChemistryIntervalReport integrate(const ChemistryIntervalRequest &) override;
  portable::Status advance_gas(const portable::GasAdvanceQuery &,
                               portable::GasAdvanceOutput &) noexcept override;

private:
  CompositionIdentity composition_;
  portable::GasIdentity gas_;
  combustion::ChemistryIdentity closure_;
  double rate_{};
  double cp_{1000};
  std::size_t reactant_{}, product_{1};
};
} // namespace hundun::v04::chemistry::detail
