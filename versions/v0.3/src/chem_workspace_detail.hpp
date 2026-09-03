// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "chem_backend_detail.hpp"
#include "chem_thermodynamics_service_detail.hpp"

#include "hundun/cfg_resolved_case_v4.hpp"

#include <cstddef>
#include <memory>
#include <string_view>

namespace hundun::chemistry {

class CanteraBackend;
class CanteraWorkspacePool;

class CanteraBackendRuntime final {
public:
  explicit CanteraBackendRuntime(const config::ResolvedReactingCaseV4 &);
  ~CanteraBackendRuntime() noexcept;
  CanteraBackendRuntime(const CanteraBackendRuntime &) = delete;
  CanteraBackendRuntime &operator=(const CanteraBackendRuntime &) = delete;

  const CompositionIdentity &composition() const noexcept;
  std::string_view mechanism_sha256() const noexcept;
  std::string_view mechanism_phase() const noexcept;

private:
  bool matches(const config::ResolvedReactingCaseV4 &) const;
  struct Impl;
  std::unique_ptr<Impl> impl_;
  friend class CanteraWorkspacePool;
  friend std::unique_ptr<CanteraBackend>
  make_cantera_backend(const config::ResolvedReactingCaseV4 &,
                       CanteraWorkspacePool &);
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
  make_cantera_backend(const config::ResolvedReactingCaseV4 &,
                       CanteraWorkspacePool &);
};

class CanteraBackend final : public ThermodynamicsService,
                             public TransportPropertyService,
                             public ChemistryBackend {
public:
  ~CanteraBackend() override;
  CanteraBackend(const CanteraBackend &) = delete;
  CanteraBackend &operator=(const CanteraBackend &) = delete;

  std::size_t lane_index() const noexcept;
  const CompositionIdentity &composition() const noexcept override;
  ThermodynamicProperties
  evaluate(const ThermochemicalPoint &) const override;
  TransportProperties
  evaluate(const ThermochemicalPoint &,
           const ThermodynamicProperties &) const override;
  ChemistryIntervalReport
  integrate(const ChemistryIntervalRequest &) override;

private:
  struct Impl;
  explicit CanteraBackend(std::unique_ptr<Impl>) noexcept;
  std::unique_ptr<Impl> impl_;
  friend std::unique_ptr<CanteraBackend>
  make_cantera_backend(const config::ResolvedReactingCaseV4 &,
                       CanteraWorkspacePool &);
};

std::unique_ptr<CanteraBackend>
make_cantera_backend(const config::ResolvedReactingCaseV4 &,
                     CanteraWorkspacePool &);

} // namespace hundun::chemistry
