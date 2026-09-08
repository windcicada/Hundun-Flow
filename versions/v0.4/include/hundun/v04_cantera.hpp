// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/v04_chemistry.hpp"
#include "hundun/v04_combustion.hpp"
#include "hundun/v04_portable.hpp"

#include <filesystem>

#include <cstddef>
#include <memory>
#include <string_view>

namespace hundun::v04::chemistry {

// Backend-local values, not a product case schema. No flow/field/checkpoint
// state.
struct CanteraMechanismConfig final {
  std::filesystem::path file;
  std::string sha256;
  std::string phase;
};
struct ChemistrySolverConfig final {
  double relative_tolerance{};
  double absolute_tolerance{};
  int maximum_internal_steps{};
};
struct CanteraBackendConfig final {
  CanteraMechanismConfig mechanism;
  ChemistrySolverConfig chemistry;
  std::vector<std::string> species_names;
};

// Runtime/pool construction may throw. Each backend owns one exclusive lane;
// the pool must outlive its backends. No method publishes product state.
// Thermodynamic/transport errors throw; integrate also returns typed failures.
// Vector allocation can throw and must be contained by the eventual caller.
class CanteraBackend;
class CanteraWorkspacePool;

class CanteraBackendRuntime final {
public:
  explicit CanteraBackendRuntime(const CanteraBackendConfig &);
  ~CanteraBackendRuntime() noexcept;
  CanteraBackendRuntime(const CanteraBackendRuntime &) = delete;
  CanteraBackendRuntime &operator=(const CanteraBackendRuntime &) = delete;

  const CompositionIdentity &composition() const noexcept;
  std::string_view mechanism_sha256() const noexcept;
  std::string_view mechanism_phase() const noexcept;

private:
  bool matches(const CanteraBackendConfig &) const;
  struct Impl;
  std::unique_ptr<Impl> impl_;
  friend class CanteraWorkspacePool;
  friend class CanteraBackend;
  friend std::unique_ptr<CanteraBackend>
  make_cantera_backend(const CanteraBackendConfig &, CanteraWorkspacePool &);
};

class CanteraWorkspacePool final {
public:
  CanteraWorkspacePool(std::shared_ptr<const CanteraBackendRuntime>,
                       std::size_t workspace_count);
  ~CanteraWorkspacePool() noexcept;
  CanteraWorkspacePool(const CanteraWorkspacePool &) = delete;
  CanteraWorkspacePool &operator=(const CanteraWorkspacePool &) = delete;

  std::size_t workspace_count() const noexcept;
  bool workspaces_are_distinct() const noexcept;

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
  friend class CanteraBackend;
  friend std::unique_ptr<CanteraBackend>
  make_cantera_backend(const CanteraBackendConfig &, CanteraWorkspacePool &);
};

class CanteraBackend final : public ThermodynamicsService,
                             public TransportPropertyService,
                             public ChemistryBackend,
                             public portable::GasQueryProvider,
                             public portable::GasAdvanceProvider {
public:
  ~CanteraBackend() override;
  CanteraBackend(const CanteraBackend &) = delete;
  CanteraBackend &operator=(const CanteraBackend &) = delete;

  std::size_t lane_index() const noexcept;
  const portable::GasIdentity &gas_identity() const noexcept override;
  portable::Status query_gas(const portable::GasQuery &,
                             portable::GasQueryOutput &) noexcept override;
  const combustion::ChemistryIdentity &closure_identity() const noexcept;
  portable::Status advance_gas(const portable::GasAdvanceQuery &,
                               portable::GasAdvanceOutput &) noexcept override;
  const CompositionIdentity &composition() const noexcept override;
  ThermodynamicProperties evaluate(const ThermochemicalPoint &) const override;
  TransportProperties evaluate(const ThermochemicalPoint &,
                               const ThermodynamicProperties &) const override;
  ChemistryIntervalReport integrate(const ChemistryIntervalRequest &) override;

private:
  struct Impl;
  explicit CanteraBackend(std::unique_ptr<Impl>) noexcept;
  std::unique_ptr<Impl> impl_;
  friend std::unique_ptr<CanteraBackend>
  make_cantera_backend(const CanteraBackendConfig &, CanteraWorkspacePool &);
};

std::unique_ptr<CanteraBackend>
make_cantera_backend(const CanteraBackendConfig &, CanteraWorkspacePool &);

} // namespace hundun::v04::chemistry
