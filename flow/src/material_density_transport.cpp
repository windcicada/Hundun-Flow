// SPDX-License-Identifier: Apache-2.0

#include "hundun/flow/material_density_transport.hpp"

#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/flow/flow_state.hpp"
#include "hundun/mesh/mesh_geometry.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/field_access_plan.hpp"
#include "hundun/runtime/field_registry.hpp"
#include "hundun/runtime/field_storage.hpp"
#include "hundun/runtime/halo_exchange.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_operation_error.hpp"
#include "hundun/runtime/structured_decomposition.hpp"
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
#include "material_density_transport_test_access.hpp"
#endif

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <utility>

namespace hundun::flow {
namespace {

constexpr runtime::PhaseId kMaterialPhase = 1900U;
constexpr runtime::ActorId kMaterialActor = 1900U;
constexpr runtime::PhaseId kScratchPhase = 1901U;
constexpr runtime::ActorId kScratchActor = 1901U;

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
struct GateOverride final {
  std::optional<double> value;
  std::size_t field{};
  int rank{-1};
};
struct MaterialTestControl final {
  GateOverride density_residual;
  GateOverride transport_residual;
  GateOverride mass_conservation;
  GateOverride transport_conservation;
};
MaterialTestControl material_test_control;
#endif

runtime::FieldDescriptor
scratch_cell(std::string name, std::uint32_t components, int ghost_width) {
  return {std::move(name),
          "1",
          "material-density-transport",
          runtime::FunctionSpace::cell_average,
          runtime::ScalarType::float64,
          components,
          ghost_width,
          false,
          runtime::RestartPolicy::transient,
          runtime::OutputPolicy::never};
}

runtime::FieldDescriptor scratch_face(std::string name) {
  return {std::move(name),
          "1",
          "material-density-transport",
          runtime::FunctionSpace::face_value,
          runtime::ScalarType::float64,
          1U,
          0,
          false,
          runtime::RestartPolicy::transient,
          runtime::OutputPolicy::never};
}

bool same(runtime::Int3 left, runtime::Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

runtime::Int3 local_index(const mesh::MeshTopology &topology,
                          mesh::LocalCellId cell) {
  const auto global = topology.global_cell(cell);
  const auto box = topology.owned_global_box();
  return {global.x - box.begin.x, global.y - box.begin.y,
          global.z - box.begin.z};
}

template <class Function>
void for_each_owned_cell(const mesh::MeshTopology &topology,
                         Function &&function) {
  for (mesh::LocalCellId cell = 0; cell < topology.local_cell_count(); ++cell) {
    if (topology.cell_ownership(cell) == mesh::EntityOwnership::owned)
      function(cell, local_index(topology, cell));
  }
}

void require_conserved_cell(const runtime::FieldRegistry &registry,
                            runtime::FieldId field, std::string_view unit,
                            const char *name) {
  const auto &descriptor = registry.descriptor(field);
  if (descriptor.space != runtime::FunctionSpace::cell_average ||
      descriptor.scalar_type != runtime::ScalarType::float64 ||
      descriptor.components != 1U || descriptor.ghost_width < 2 ||
      !descriptor.conservative ||
      descriptor.restart != runtime::RestartPolicy::persistent ||
      descriptor.unit != unit) {
    throw runtime::Error(std::string("material transport invalid ") + name +
                         " descriptor");
  }
}

std::string scalar_field_id(std::size_t index) {
  std::ostringstream output;
  output.imbue(std::locale::classic());
  output << "rho_phi.s" << std::setfill('0') << std::setw(20) << index;
  return output.str();
}

std::string global_cell_layout_fingerprint(runtime::Int3 extent) {
  std::ostringstream output;
  output.imbue(std::locale::classic());
  output << "cell.f64.c1.g2plus.global." << extent.x << '.' << extent.y << '.'
         << extent.z;
  return output.str();
}

std::string owned_cell_layout_fingerprint(runtime::Box3 box) {
  std::ostringstream output;
  output.imbue(std::locale::classic());
  output << "cell.f64.c1.g2plus.owned." << box.begin.x << '.' << box.begin.y
         << '.' << box.begin.z << '.' << box.end.x << '.' << box.end.y << '.'
         << box.end.z;
  return output.str();
}

struct SynchronizedFailure final {
  MaterialTransportFailureReason reason{MaterialTransportFailureReason::none};
  int rank{-1};
};

class TransportComputationFailure final : public runtime::Error {
public:
  TransportComputationFailure(SynchronizedFailure failure,
                              const std::string &message)
      : runtime::Error(message), failure_(failure) {}

  SynchronizedFailure failure() const noexcept { return failure_; }

private:
  SynchronizedFailure failure_;
};

SynchronizedFailure
synchronize_failure(const runtime::MpiContext &mpi,
                    MaterialTransportFailureReason local_reason) {
  int failing = local_reason == MaterialTransportFailureReason::none
                    ? mpi.size()
                    : mpi.rank();
  int lowest = mpi.size();
  runtime::check_mpi_result(
      MPI_Allreduce(&failing, &lowest, 1, MPI_INT, MPI_MIN, mpi.comm()),
      "MPI_Allreduce(material failure rank)");
  if (lowest == mpi.size())
    return {};
  int reason = mpi.rank() == lowest ? static_cast<int>(local_reason) : 0;
  runtime::check_mpi_result(MPI_Bcast(&reason, 1, MPI_INT, lowest, mpi.comm()),
                            "MPI_Bcast(material failure reason)");
  return {static_cast<MaterialTransportFailureReason>(reason), lowest};
}

double relative_defect(double raw, double scale) noexcept {
  return std::abs(raw) / std::max(scale, std::numeric_limits<double>::min());
}

double conservation_denominator(double current_integral, double next_integral,
                                double history_integral,
                                double history_coefficient,
                                double boundary_scale,
                                double cancellation_current,
                                double cancellation_next,
                                double cancellation_history) noexcept {
  const double ordinary = std::max(
      {std::abs(current_integral), std::abs(next_integral),
       std::abs(history_coefficient * (history_integral - current_integral)),
       std::abs(boundary_scale), std::numeric_limits<double>::min()});
  const double cancellation_scale =
      std::max({cancellation_current, cancellation_next, cancellation_history});
  if (cancellation_scale > 0.0 &&
      ordinary <=
          64.0 * std::numeric_limits<double>::epsilon() * cancellation_scale)
    return cancellation_scale;
  return ordinary;
}

void hash_u64(std::uint64_t &hash, std::uint64_t value) noexcept {
  for (unsigned shift = 0; shift < 64U; shift += 8U) {
    hash ^= (value >> shift) & 0xffU;
    hash *= UINT64_C(1099511628211);
  }
}

std::uint64_t double_bits(double value) noexcept {
  std::uint64_t bits{};
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

class CompensatedSum final {
public:
  void add(double value) noexcept {
    const double adjusted = value - correction_;
    const double next = sum_ + adjusted;
    correction_ = (next - sum_) - adjusted;
    sum_ = next;
  }
  double value() const noexcept { return sum_; }

private:
  double sum_{};
  double correction_{};
};

} // namespace

struct MaterialFaceMassFlux::Impl final {
  Impl(const runtime::FieldRegistry &supplied_registry,
       const mesh::MeshTopology &supplied_topology,
       finite_volume::FaceMassFlux supplied_handle,
       runtime::FaceFieldView<const double> supplied_view,
       MaterialFluxProvenance supplied_provenance)
      : registry(&supplied_registry), topology(&supplied_topology),
        handle(std::move(supplied_handle)), view(std::move(supplied_view)),
        provenance(supplied_provenance) {}
  const runtime::FieldRegistry *registry;
  const mesh::MeshTopology *topology;
  finite_volume::FaceMassFlux handle;
  runtime::FaceFieldView<const double> view;
  MaterialFluxProvenance provenance;
};

MaterialFaceMassFlux MaterialFaceMassFlux::acquire(
    const runtime::FieldRegistry &registry,
    const runtime::FieldStorage &storage,
    const runtime::FieldAccessPlan &access_plan, runtime::PhaseId phase,
    runtime::ActorId actor, runtime::FieldId field,
    const mesh::MeshTopology &topology, MaterialFluxProvenance provenance) {
  auto handle = finite_volume::FaceMassFlux::acquire(
      registry, storage, access_plan, phase, actor, field, topology);
  auto view =
      storage.acquire_face_read<double>(access_plan, phase, actor, field);
  return MaterialFaceMassFlux(std::make_unique<Impl>(
      registry, topology, std::move(handle), std::move(view), provenance));
}
MaterialFaceMassFlux::MaterialFaceMassFlux(std::unique_ptr<Impl> impl) noexcept
    : impl_(std::move(impl)) {}
MaterialFaceMassFlux::~MaterialFaceMassFlux() noexcept = default;
MaterialFaceMassFlux::MaterialFaceMassFlux(MaterialFaceMassFlux &&) noexcept =
    default;
runtime::FieldId MaterialFaceMassFlux::field_id() const noexcept {
  return impl_ ? impl_->handle.field_id() : runtime::FieldId{};
}
std::size_t MaterialFaceMassFlux::face_count() const noexcept {
  return impl_ ? impl_->handle.face_count() : 0U;
}
MaterialFluxProvenance MaterialFaceMassFlux::provenance() const noexcept {
  return impl_ ? impl_->provenance : MaterialFluxProvenance::predictor;
}

MaterialTransportDisposition
MaterialDensityTransportReport::disposition() const noexcept {
  return disposition_;
}
MaterialTransportFailureReason
MaterialDensityTransportReport::reason() const noexcept {
  return reason_;
}
int MaterialDensityTransportReport::lowest_failing_rank() const noexcept {
  return lowest_failing_rank_;
}
const MomentumTimeStencil &
MaterialDensityTransportReport::stencil() const noexcept {
  return stencil_;
}
MaterialFluxProvenance
MaterialDensityTransportReport::flux_provenance() const noexcept {
  return flux_provenance_;
}
std::uint64_t
MaterialDensityTransportReport::attempt_identity() const noexcept {
  return attempt_identity_;
}
std::uint64_t
MaterialDensityTransportReport::finalization_identity() const noexcept {
  return finalization_identity_;
}
runtime::FieldId
MaterialDensityTransportReport::shared_face_mass_flux_field() const noexcept {
  return shared_face_mass_flux_field_;
}
bool MaterialDensityTransportReport::density_residual_available()
    const noexcept {
  return density_residual_available_;
}
double MaterialDensityTransportReport::density_normalized_l2() const noexcept {
  return density_normalized_l2_;
}
const std::vector<std::uint8_t> &
MaterialDensityTransportReport::transport_residual_availability()
    const noexcept {
  return transport_residual_available_;
}
const std::vector<double> &
MaterialDensityTransportReport::transport_normalized_l2() const noexcept {
  return transport_normalized_l2_;
}
bool MaterialDensityTransportReport::mass_conservation_available()
    const noexcept {
  return mass_conservation_available_;
}
double MaterialDensityTransportReport::mass_relative_conservation_defect()
    const noexcept {
  return mass_relative_conservation_defect_;
}
const std::vector<std::uint8_t> &
MaterialDensityTransportReport::transport_conservation_availability()
    const noexcept {
  return transport_conservation_available_;
}
const std::vector<double> &
MaterialDensityTransportReport::transport_relative_conservation_defect()
    const noexcept {
  return transport_relative_conservation_defect_;
}
bool MaterialDensityTransportReport::minimum_density_available()
    const noexcept {
  return minimum_density_available_;
}
double
MaterialDensityTransportReport::minimum_density_kg_per_m3() const noexcept {
  return minimum_density_kg_per_m3_;
}
mesh::GlobalCellId
MaterialDensityTransportReport::minimum_density_global_cell() const noexcept {
  return minimum_density_global_cell_;
}
int MaterialDensityTransportReport::minimum_density_rank() const noexcept {
  return minimum_density_rank_;
}

std::uint64_t MaterialDensityTransportReport::compute_seal() const noexcept {
  std::uint64_t hash = UINT64_C(14695981039346656037);
  hash_u64(hash, static_cast<std::uint64_t>(disposition_));
  hash_u64(hash, static_cast<std::uint64_t>(reason_));
  hash_u64(hash, static_cast<std::uint64_t>(lowest_failing_rank_));
  hash_u64(hash, static_cast<std::uint64_t>(stencil_.order));
  hash_u64(hash, double_bits(stencil_.dt_s));
  hash_u64(hash, double_bits(stencil_.previous_dt_s));
  hash_u64(hash, double_bits(stencil_.alpha0));
  hash_u64(hash, double_bits(stencil_.alpha1));
  hash_u64(hash, double_bits(stencil_.alpha2));
  hash_u64(hash, static_cast<std::uint64_t>(flux_provenance_));
  hash_u64(hash, attempt_identity_);
  hash_u64(hash, finalization_identity_);
  hash_u64(hash, static_cast<std::uint64_t>(shared_face_mass_flux_field_));
  hash_u64(hash, density_residual_available_ ? 1U : 0U);
  hash_u64(hash, double_bits(density_normalized_l2_));
  hash_u64(hash, transport_residual_available_.size());
  for (const auto value : transport_residual_available_)
    hash_u64(hash, value);
  hash_u64(hash, transport_normalized_l2_.size());
  for (const double value : transport_normalized_l2_)
    hash_u64(hash, double_bits(value));
  hash_u64(hash, mass_conservation_available_ ? 1U : 0U);
  hash_u64(hash, double_bits(mass_relative_conservation_defect_));
  hash_u64(hash, transport_conservation_available_.size());
  for (const auto value : transport_conservation_available_)
    hash_u64(hash, value);
  hash_u64(hash, transport_relative_conservation_defect_.size());
  for (const double value : transport_relative_conservation_defect_)
    hash_u64(hash, double_bits(value));
  hash_u64(hash, minimum_density_available_ ? 1U : 0U);
  hash_u64(hash, double_bits(minimum_density_kg_per_m3_));
  hash_u64(hash, minimum_density_global_cell_);
  hash_u64(hash, static_cast<std::uint64_t>(minimum_density_rank_));
  return hash == 0U ? 1U : hash;
}
void MaterialDensityTransportReport::seal() noexcept { seal_ = compute_seal(); }
bool MaterialDensityTransportReport::authenticated() const noexcept {
  return seal_ != 0U && seal_ == compute_seal();
}

struct MaterialDensityTransport::Impl final {
  Impl(const runtime::FieldRegistry &supplied_registry,
       const runtime::StructuredDecomposition &supplied_decomposition,
       const mesh::MeshTopology &supplied_topology,
       const mesh::MeshGeometry &supplied_geometry,
       const boundary::BoundaryRegistry &supplied_boundaries,
       const runtime::MpiContext &supplied_mpi,
       runtime::HaloExchange &supplied_halo, FlowFieldIds supplied_fields,
       MaterialDensityTransportSpec supplied_specification)
      : registry(&supplied_registry), decomposition(&supplied_decomposition),
        topology(&supplied_topology), geometry(&supplied_geometry),
        boundaries(&supplied_boundaries), mpi(&supplied_mpi),
        halo(&supplied_halo), fields(std::move(supplied_fields)),
        specification(std::move(supplied_specification)),
        access(supplied_registry),
        operators(finite_volume::CellCenteredFvmOperators::create(
            supplied_topology, supplied_geometry)) {
    const auto box = supplied_topology.owned_global_box();
    local_extent = {box.end.x - box.begin.x, box.end.y - box.begin.y,
                    box.end.z - box.begin.z};
    global_layout_fingerprint =
        global_cell_layout_fingerprint(supplied_topology.global_extent());
    owned_layout_fingerprint = owned_cell_layout_fingerprint(box);
    scratch_intensive_n =
        scratch_registry.declare_field(scratch_cell("material_q_n", 1U, 2));
    scratch_intensive_h = scratch_registry.declare_field(
        scratch_cell("material_q_history", 1U, 2));
    scratch_gradient_n =
        scratch_registry.declare_field(scratch_cell("material_grad_n", 3U, 2));
    scratch_gradient_h = scratch_registry.declare_field(
        scratch_cell("material_grad_history", 3U, 2));
    scratch_face_n =
        scratch_registry.declare_field(scratch_face("material_face_n"));
    scratch_face_h =
        scratch_registry.declare_field(scratch_face("material_face_history"));
    scratch_gamma =
        scratch_registry.declare_field(scratch_face("material_gamma"));
    scratch_mass_residual = scratch_registry.declare_field(
        scratch_cell("material_mass_residual", 1U, 0));
    scratch_residual_n = scratch_registry.declare_field(
        scratch_cell("material_residual_n", 1U, 0));
    scratch_residual_h = scratch_registry.declare_field(
        scratch_cell("material_residual_history", 1U, 0));
    scratch_registry.freeze();
    scratch_access =
        std::make_unique<runtime::FieldAccessPlan>(scratch_registry);
    for (runtime::FieldId field = 0;
         field < static_cast<runtime::FieldId>(scratch_registry.size());
         ++field)
      scratch_access->declare_access(kScratchPhase, kScratchActor, field,
                                     runtime::AccessMode::read_write);
    scratch_access->freeze();
    scratch_storage = std::make_unique<runtime::FieldStorage>(
        scratch_registry,
        runtime::FieldLayoutSet{local_extent,
                                supplied_topology.local_face_count()});

    access.declare_access(kMaterialPhase, kMaterialActor, fields.density,
                          runtime::AccessMode::read_write);
    access.declare_access(kMaterialPhase, kMaterialActor, fields.face_mass_flux,
                          runtime::AccessMode::read_write);
    access.declare_access(kMaterialPhase, kMaterialActor,
                          specification.enthalpy_density,
                          runtime::AccessMode::read_write);
    for (const auto field : specification.scalar_densities)
      access.declare_access(kMaterialPhase, kMaterialActor, field,
                            runtime::AccessMode::read_write);
    access.freeze();

    fingerprint_ids.push_back("rho");
    fingerprint_ids.push_back("rho_h");
    for (std::size_t index = 0; index < specification.scalar_densities.size();
         ++index)
      fingerprint_ids.push_back(scalar_field_id(index));
  }

  const runtime::FieldRegistry *registry;
  const runtime::StructuredDecomposition *decomposition;
  const mesh::MeshTopology *topology;
  const mesh::MeshGeometry *geometry;
  const boundary::BoundaryRegistry *boundaries;
  const runtime::MpiContext *mpi;
  runtime::HaloExchange *halo;
  FlowFieldIds fields;
  MaterialDensityTransportSpec specification;
  runtime::FieldAccessPlan access;
  finite_volume::CellCenteredFvmOperators operators;
  runtime::Int3 local_extent{};
  runtime::FieldRegistry scratch_registry;
  std::unique_ptr<runtime::FieldAccessPlan> scratch_access;
  std::unique_ptr<runtime::FieldStorage> scratch_storage;
  runtime::FieldId scratch_intensive_n{};
  runtime::FieldId scratch_intensive_h{};
  runtime::FieldId scratch_gradient_n{};
  runtime::FieldId scratch_gradient_h{};
  runtime::FieldId scratch_face_n{};
  runtime::FieldId scratch_face_h{};
  runtime::FieldId scratch_gamma{};
  runtime::FieldId scratch_mass_residual{};
  runtime::FieldId scratch_residual_n{};
  runtime::FieldId scratch_residual_h{};
  std::vector<std::string> fingerprint_ids;
  std::string global_layout_fingerprint;
  std::string owned_layout_fingerprint;
  mutable std::uint64_t finalization_identity{};
  mutable const FlowState *last_state{};
  mutable std::uint64_t last_attempt_identity{};
  mutable std::uint64_t last_finalization_identity{};
  mutable std::uint64_t last_report_seal{};
  mutable bool active{};
};

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
namespace test {
void MaterialDensityTransportTestAccess::reset() noexcept {
  material_test_control = {};
}
void MaterialDensityTransportTestAccess::set_density_residual(
    double value, int rank) noexcept {
  material_test_control.density_residual = {{value}, 0U, rank};
}
void MaterialDensityTransportTestAccess::set_transport_residual(
    std::size_t field, double value, int rank) noexcept {
  material_test_control.transport_residual = {{value}, field, rank};
}
void MaterialDensityTransportTestAccess::set_mass_conservation_defect(
    double value, int rank) noexcept {
  material_test_control.mass_conservation = {{value}, 0U, rank};
}
void MaterialDensityTransportTestAccess::set_transport_conservation_defect(
    std::size_t field, double value, int rank) noexcept {
  material_test_control.transport_conservation = {{value}, field, rank};
}
void MaterialDensityTransportTestAccess::force_finalization_identity_wrap(
    MaterialDensityTransport &transport) noexcept {
  transport.impl_->finalization_identity =
      std::numeric_limits<std::uint64_t>::max();
}
double MaterialDensityTransportTestAccess::conservation_denominator(
    const MaterialConservationScaleInput &input) noexcept {
  return ::hundun::flow::conservation_denominator(
      input.current_integral, input.next_integral, input.history_integral,
      input.history_coefficient, input.boundary_scale,
      input.cancellation_current, input.cancellation_next,
      input.cancellation_history);
}
void MaterialDensityTransportTestAccess::corrupt_report(
    MaterialDensityTransportReport &report,
    MaterialReportCorruption corruption) noexcept {
  switch (corruption) {
  case MaterialReportCorruption::scalar:
    ++report.finalization_identity_;
    break;
  case MaterialReportCorruption::vector_size:
    report.transport_normalized_l2_.push_back(0.0);
    break;
  case MaterialReportCorruption::vector_element:
    if (!report.transport_normalized_l2_.empty())
      report.transport_normalized_l2_.front() =
          std::nextafter(report.transport_normalized_l2_.front(), 1.0);
    break;
  case MaterialReportCorruption::availability:
    report.density_residual_available_ = !report.density_residual_available_;
    break;
  case MaterialReportCorruption::seal:
    report.seal_ ^= 1U;
    break;
  }
}
} // namespace test
#endif

namespace {

void validate_specification(const runtime::FieldRegistry &registry,
                            const MaterialDensityTransportSpec &specification,
                            const boundary::BoundaryRegistry &boundaries) {
  require_conserved_cell(registry, specification.enthalpy_density, "J/m3",
                         "enthalpy-density");
  if (specification.scalar_densities.size() !=
          specification.scalar_diffusivities_kg_per_m_s.size() ||
      specification.scalar_densities.size() != boundaries.scalar_count())
    throw runtime::Error("material transport scalar counts do not match");
  if (!std::isfinite(specification.enthalpy_diffusivity_kg_per_m_s) ||
      specification.enthalpy_diffusivity_kg_per_m_s < 0.0)
    throw runtime::Error("material transport enthalpy diffusivity is invalid");
  for (std::size_t index = 0; index < specification.scalar_densities.size();
       ++index) {
    require_conserved_cell(registry, specification.scalar_densities[index],
                           "kg/m3", "scalar-density");
    const double gamma = specification.scalar_diffusivities_kg_per_m_s[index];
    if (!std::isfinite(gamma) || gamma < 0.0)
      throw runtime::Error("material transport scalar diffusivity is invalid");
  }
}

template <class Implementation>
std::vector<runtime::FieldId> transport_fields(const Implementation &impl) {
  std::vector<runtime::FieldId> result{impl.specification.enthalpy_density};
  result.insert(result.end(), impl.specification.scalar_densities.begin(),
                impl.specification.scalar_densities.end());
  return result;
}

finite_volume::FiniteVolumeQuantity quantity(std::size_t field) {
  return field == 0U ? finite_volume::FiniteVolumeQuantity::enthalpy()
                     : finite_volume::FiniteVolumeQuantity::scalar(field - 1U);
}

template <class Implementation>
double diffusivity(const Implementation &impl, std::size_t field) {
  return field == 0U
             ? impl.specification.enthalpy_diffusivity_kg_per_m_s
             : impl.specification.scalar_diffusivities_kg_per_m_s[field - 1U];
}

struct TransportResidual final {
  std::vector<double> current;
  std::vector<double> history;
  double boundary_current{};
  double boundary_history{};
  double boundary_abs_current{};
  double boundary_abs_history{};
};

template <class Implementation>
TransportResidual compute_transport_residual(
    Implementation &impl, const finite_volume::FaceMassFlux &mass_flux,
    const FlowState &state, runtime::FieldId conserved_field,
    std::size_t field_index) {
  auto &scratch = *impl.scratch_storage;
  auto qn = scratch.template acquire_write<double>(*impl.scratch_access,
                                                   kScratchPhase, kScratchActor,
                                                   impl.scratch_intensive_n);
  auto qh = scratch.template acquire_write<double>(*impl.scratch_access,
                                                   kScratchPhase, kScratchActor,
                                                   impl.scratch_intensive_h);
  const auto rho_n =
      state.layer(FlowLayer::committed)
          .acquire_read<double>(impl.access, kMaterialPhase, kMaterialActor,
                                state.fields().density);
  const auto rho_h =
      state.layer(FlowLayer::history)
          .acquire_read<double>(impl.access, kMaterialPhase, kMaterialActor,
                                state.fields().density);
  const auto conserved_n =
      state.layer(FlowLayer::committed)
          .acquire_read<double>(impl.access, kMaterialPhase, kMaterialActor,
                                conserved_field);
  const auto conserved_h =
      state.layer(FlowLayer::history)
          .acquire_read<double>(impl.access, kMaterialPhase, kMaterialActor,
                                conserved_field);
  MaterialTransportFailureReason local_failure =
      MaterialTransportFailureReason::none;
  try {
    for_each_owned_cell(*impl.topology, [&](mesh::LocalCellId,
                                            runtime::Int3 index) {
      const double density_n = rho_n(index.x, index.y, index.z, 0);
      const double density_h = rho_h(index.x, index.y, index.z, 0);
      const double value_n = conserved_n(index.x, index.y, index.z, 0);
      const double value_h = conserved_h(index.x, index.y, index.z, 0);
      if (!(density_n > 0.0) || !(density_h > 0.0) ||
          !std::isfinite(density_n) || !std::isfinite(density_h) ||
          !std::isfinite(value_n) || !std::isfinite(value_h))
        throw runtime::Error("material transport intensive state is invalid");
      qn(index.x, index.y, index.z, 0) = value_n / density_n;
      qh(index.x, index.y, index.z, 0) = value_h / density_h;
    });
  } catch (const runtime::MpiOperationError &) {
    throw;
  } catch (const std::exception &) {
    local_failure = MaterialTransportFailureReason::non_finite_state;
  }
  auto synchronized = synchronize_failure(*impl.mpi, local_failure);
  if (synchronized.reason != MaterialTransportFailureReason::none)
    throw TransportComputationFailure(
        synchronized, "material transport intensive state is invalid");
  impl.halo->exchange(scratch, impl.scratch_intensive_n);
  impl.halo->exchange(scratch, impl.scratch_intensive_h);

  const auto qn_read = std::as_const(scratch).template acquire_read<double>(
      *impl.scratch_access, kScratchPhase, kScratchActor,
      impl.scratch_intensive_n);
  const auto qh_read = std::as_const(scratch).template acquire_read<double>(
      *impl.scratch_access, kScratchPhase, kScratchActor,
      impl.scratch_intensive_h);
  auto grad_n = scratch.template acquire_write<double>(
      *impl.scratch_access, kScratchPhase, kScratchActor,
      impl.scratch_gradient_n);
  auto grad_h = scratch.template acquire_write<double>(
      *impl.scratch_access, kScratchPhase, kScratchActor,
      impl.scratch_gradient_h);
  local_failure = MaterialTransportFailureReason::none;
  try {
    impl.operators.compute_gradient(finite_volume::GradientScheme::green_gauss,
                                    quantity(field_index), *impl.boundaries,
                                    qn_read, grad_n);
    impl.operators.compute_gradient(finite_volume::GradientScheme::green_gauss,
                                    quantity(field_index), *impl.boundaries,
                                    qh_read, grad_h);
  } catch (const runtime::MpiOperationError &) {
    throw;
  } catch (const std::exception &) {
    local_failure = MaterialTransportFailureReason::non_finite_state;
  }
  synchronized = synchronize_failure(*impl.mpi, local_failure);
  if (synchronized.reason != MaterialTransportFailureReason::none)
    throw TransportComputationFailure(synchronized,
                                      "material transport gradient is invalid");
  impl.halo->exchange(scratch, impl.scratch_gradient_n);
  impl.halo->exchange(scratch, impl.scratch_gradient_h);

  auto face_n = scratch.template acquire_face_write<double>(
      *impl.scratch_access, kScratchPhase, kScratchActor, impl.scratch_face_n);
  auto face_h = scratch.template acquire_face_write<double>(
      *impl.scratch_access, kScratchPhase, kScratchActor, impl.scratch_face_h);
  impl.operators.reconstruct_transport_faces(
      quantity(field_index), *impl.boundaries, mass_flux, qn_read, face_n);
  impl.operators.reconstruct_transport_faces(
      quantity(field_index), *impl.boundaries, mass_flux, qh_read, face_h);
  auto gamma = scratch.template acquire_face_write<double>(
      *impl.scratch_access, kScratchPhase, kScratchActor, impl.scratch_gamma);
  for (std::size_t face = 0; face < impl.topology->local_face_count(); ++face)
    gamma(face, 0) = diffusivity(impl, field_index);

  const auto grad_n_read = std::as_const(scratch).template acquire_read<double>(
      *impl.scratch_access, kScratchPhase, kScratchActor,
      impl.scratch_gradient_n);
  const auto grad_h_read = std::as_const(scratch).template acquire_read<double>(
      *impl.scratch_access, kScratchPhase, kScratchActor,
      impl.scratch_gradient_h);
  const auto face_n_read =
      std::as_const(scratch).template acquire_face_read<double>(
          *impl.scratch_access, kScratchPhase, kScratchActor,
          impl.scratch_face_n);
  const auto face_h_read =
      std::as_const(scratch).template acquire_face_read<double>(
          *impl.scratch_access, kScratchPhase, kScratchActor,
          impl.scratch_face_h);
  const auto gamma_read =
      std::as_const(scratch).template acquire_face_read<double>(
          *impl.scratch_access, kScratchPhase, kScratchActor,
          impl.scratch_gamma);

  auto residual_n = scratch.template acquire_write<double>(
      *impl.scratch_access, kScratchPhase, kScratchActor,
      impl.scratch_residual_n);
  auto residual_h = scratch.template acquire_write<double>(
      *impl.scratch_access, kScratchPhase, kScratchActor,
      impl.scratch_residual_h);
  for_each_owned_cell(*impl.topology,
                      [&](mesh::LocalCellId, runtime::Int3 index) {
                        residual_n(index.x, index.y, index.z, 0) = 0.0;
                        residual_h(index.x, index.y, index.z, 0) = 0.0;
                      });
  impl.operators.accumulate_convective_residual(mass_flux, face_n_read,
                                                residual_n);
  impl.operators.accumulate_scalar_diffusive_residual(
      quantity(field_index), *impl.boundaries, qn_read, grad_n_read, gamma_read,
      residual_n);
  impl.operators.accumulate_convective_residual(mass_flux, face_h_read,
                                                residual_h);
  impl.operators.accumulate_scalar_diffusive_residual(
      quantity(field_index), *impl.boundaries, qh_read, grad_h_read, gamma_read,
      residual_h);

  TransportResidual result;
  result.current.reserve(impl.topology->owned_cell_count());
  result.history.reserve(impl.topology->owned_cell_count());
  for_each_owned_cell(
      *impl.topology, [&](mesh::LocalCellId, runtime::Int3 index) {
        result.current.push_back(residual_n(index.x, index.y, index.z, 0));
        result.history.push_back(residual_h(index.x, index.y, index.z, 0));
      });
  std::vector<finite_volume::PhysicalBoundaryTransportContribution> current;
  std::vector<finite_volume::PhysicalBoundaryTransportContribution> history;
  impl.operators.physical_boundary_transport_contributions(
      quantity(field_index), *impl.boundaries, mass_flux, face_n_read, qn_read,
      grad_n_read, gamma_read, current);
  impl.operators.physical_boundary_transport_contributions(
      quantity(field_index), *impl.boundaries, mass_flux, face_h_read, qh_read,
      grad_h_read, gamma_read, history);
  CompensatedSum boundary_current;
  CompensatedSum boundary_history;
  CompensatedSum boundary_abs_current;
  CompensatedSum boundary_abs_history;
  for (const auto &entry : current) {
    const double term = entry.convective + entry.diffusive;
    boundary_current.add(term);
    boundary_abs_current.add(std::abs(entry.convective) +
                             std::abs(entry.diffusive));
  }
  for (const auto &entry : history) {
    const double term = entry.convective + entry.diffusive;
    boundary_history.add(term);
    boundary_abs_history.add(std::abs(entry.convective) +
                             std::abs(entry.diffusive));
  }
  result.boundary_current = boundary_current.value();
  result.boundary_history = boundary_history.value();
  result.boundary_abs_current = boundary_abs_current.value();
  result.boundary_abs_history = boundary_abs_history.value();
  return result;
}

} // namespace

MaterialDensityTransport MaterialDensityTransport::create(
    const runtime::FieldRegistry &registry,
    const runtime::StructuredDecomposition &decomposition,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const boundary::BoundaryRegistry &boundaries,
    const runtime::MpiContext &mpi, runtime::HaloExchange &halo,
    FlowFieldIds expected_fields, MaterialDensityTransportSpec specification) {
  validate_specification(registry, specification, boundaries);
  require_conserved_cell(registry, expected_fields.density, "kg/m3", "density");
  finite_volume::require_face_mass_flux_field(registry,
                                              expected_fields.face_mass_flux);
  if (expected_fields.transported_cell_fields.size() !=
          1U + specification.scalar_densities.size() ||
      expected_fields.transported_cell_fields.front() !=
          specification.enthalpy_density ||
      !std::equal(specification.scalar_densities.begin(),
                  specification.scalar_densities.end(),
                  expected_fields.transported_cell_fields.begin() + 1U))
    throw runtime::Error("material transport field order is invalid");
  std::vector<runtime::FieldId> conserved{expected_fields.density,
                                          specification.enthalpy_density};
  conserved.insert(conserved.end(), specification.scalar_densities.begin(),
                   specification.scalar_densities.end());
  if (std::set<runtime::FieldId>(conserved.begin(), conserved.end()).size() !=
      conserved.size())
    throw runtime::Error("material transport fields must be unique");
  geometry.require_compatible(topology);
  if (!halo.is_compatible_with(decomposition) || halo.ghost_width() < 2)
    throw runtime::Error("material transport Halo is incompatible");
  return MaterialDensityTransport(std::make_unique<Impl>(
      registry, decomposition, topology, geometry, boundaries, mpi, halo,
      std::move(expected_fields), std::move(specification)));
}
MaterialDensityTransport::MaterialDensityTransport(
    std::unique_ptr<Impl> impl) noexcept
    : impl_(std::move(impl)) {}
MaterialDensityTransport::~MaterialDensityTransport() noexcept = default;
MaterialDensityTransport::MaterialDensityTransport(
    MaterialDensityTransport &&) noexcept = default;

MaterialDensityTransportReport MaterialDensityTransport::finalize_trial(
    FlowState &state, const MaterialFaceMassFlux &final_mass_flux,
    const MomentumTimeStencil &stencil) const {
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  const MaterialTestControl overrides =
      std::exchange(material_test_control, MaterialTestControl{});
#endif
  MaterialDensityTransportReport report;
  const std::size_t transport_field_count =
      1U + impl_->specification.scalar_densities.size();
  report.transport_residual_available_.assign(transport_field_count, 0U);
  report.transport_normalized_l2_.assign(transport_field_count, 0.0);
  report.transport_conservation_available_.assign(transport_field_count, 0U);
  report.transport_relative_conservation_defect_.assign(transport_field_count,
                                                        0.0);
  report.stencil_ = stencil;
  report.flux_provenance_ = final_mass_flux.provenance();
  report.shared_face_mass_flux_field_ = final_mass_flux.field_id();
  if (impl_->active)
    throw runtime::Error("material transport call overlaps another call");
  if (impl_->finalization_identity == std::numeric_limits<std::uint64_t>::max())
    throw runtime::Error("material finalization identity would wrap");
  ++impl_->finalization_identity;
  report.finalization_identity_ = impl_->finalization_identity;
  report.attempt_identity_ = state.attempt_identity();
  const auto finish_report = [&]() {
    report.seal();
    impl_->last_state = &state;
    impl_->last_attempt_identity = report.attempt_identity_;
    impl_->last_finalization_identity = report.finalization_identity_;
    impl_->last_report_seal = report.seal_;
    return report;
  };
  struct Guard final {
    bool &active;
    ~Guard() { active = false; }
  } guard{impl_->active};
  impl_->active = true;

  MaterialTransportFailureReason local = MaterialTransportFailureReason::none;
  if (!state.attempt_active() || &state.solver_registry() != impl_->registry ||
      state.fields().density != impl_->fields.density ||
      state.fields().face_mass_flux != impl_->fields.face_mass_flux ||
      !final_mass_flux.impl_ ||
      (final_mass_flux.impl_ &&
       (final_mass_flux.impl_->registry != impl_->registry ||
        final_mass_flux.impl_->topology != impl_->topology)) ||
      final_mass_flux.provenance() != MaterialFluxProvenance::final_corrected ||
      final_mass_flux.field_id() != state.fields().face_mass_flux ||
      final_mass_flux.face_count() != impl_->topology->local_face_count() ||
      !same(state.layer(FlowLayer::trial).interior_extent(),
            impl_->local_extent) ||
      state.fields().transported_cell_fields.size() !=
          1U + impl_->specification.scalar_densities.size() ||
      state.fields().transported_cell_fields.front() !=
          impl_->specification.enthalpy_density ||
      !std::equal(impl_->specification.scalar_densities.begin(),
                  impl_->specification.scalar_densities.end(),
                  state.fields().transported_cell_fields.begin() + 1))
    local = MaterialTransportFailureReason::invalid_input;
  if (local == MaterialTransportFailureReason::none) {
    try {
      if (final_mass_flux.face_count() != 0U)
        static_cast<void>(final_mass_flux.impl_->view(0U, 0));
    } catch (const runtime::Error &) {
      local = MaterialTransportFailureReason::invalid_input;
    }
  }
  try {
    const auto expected = make_momentum_time_stencil(
        stencil.order, stencil.dt_s, stencil.previous_dt_s);
    if (expected.alpha0 != stencil.alpha0 ||
        expected.alpha1 != stencil.alpha1 || expected.alpha2 != stencil.alpha2)
      local = MaterialTransportFailureReason::invalid_input;
  } catch (const runtime::Error &) {
    local = MaterialTransportFailureReason::invalid_input;
  }
  auto failure = synchronize_failure(*impl_->mpi, local);
  if (failure.reason != MaterialTransportFailureReason::none) {
    report.disposition_ = MaterialTransportDisposition::non_retryable_failure;
    report.reason_ = failure.reason;
    report.lowest_failing_rank_ = failure.rank;
    return finish_report();
  }

  const auto fields = transport_fields(*impl_);
  const auto rho_n =
      state.layer(FlowLayer::committed)
          .acquire_read<double>(impl_->access, kMaterialPhase, kMaterialActor,
                                state.fields().density);
  const auto rho_h =
      state.layer(FlowLayer::history)
          .acquire_read<double>(impl_->access, kMaterialPhase, kMaterialActor,
                                state.fields().density);
  bool finite_state = true;
  bool positive_density = true;
  for_each_owned_cell(
      *impl_->topology, [&](mesh::LocalCellId, runtime::Int3 index) {
        const double current = rho_n(index.x, index.y, index.z, 0);
        const double history = rho_h(index.x, index.y, index.z, 0);
        finite_state =
            finite_state && std::isfinite(current) && std::isfinite(history);
        positive_density = positive_density && current > 0.0 && history > 0.0;
      });
  for (const auto field : fields) {
    const auto current =
        state.layer(FlowLayer::committed)
            .acquire_read<double>(impl_->access, kMaterialPhase, kMaterialActor,
                                  field);
    const auto history =
        state.layer(FlowLayer::history)
            .acquire_read<double>(impl_->access, kMaterialPhase, kMaterialActor,
                                  field);
    for_each_owned_cell(
        *impl_->topology, [&](mesh::LocalCellId, runtime::Int3 index) {
          finite_state = finite_state &&
                         std::isfinite(current(index.x, index.y, index.z, 0)) &&
                         std::isfinite(history(index.x, index.y, index.z, 0));
        });
  }
  local = !finite_state
              ? MaterialTransportFailureReason::non_finite_state
              : (!positive_density
                     ? MaterialTransportFailureReason::non_positive_density
                     : MaterialTransportFailureReason::none);
  failure = synchronize_failure(*impl_->mpi, local);
  if (failure.reason != MaterialTransportFailureReason::none) {
    report.disposition_ = MaterialTransportDisposition::recoverable_failure;
    report.reason_ = failure.reason;
    report.lowest_failing_rank_ = failure.rank;
    return finish_report();
  }

  auto &scratch = *impl_->scratch_storage;
  auto mass_residual = scratch.acquire_write<double>(
      *impl_->scratch_access, kScratchPhase, kScratchActor,
      impl_->scratch_mass_residual);
  for_each_owned_cell(*impl_->topology,
                      [&](mesh::LocalCellId, runtime::Int3 index) {
                        mass_residual(index.x, index.y, index.z, 0) = 0.0;
                      });

  std::vector<double> candidate_density;
  std::vector<std::vector<double>> candidates(fields.size());
  std::vector<TransportResidual> residuals;
  local = MaterialTransportFailureReason::none;
  try {
    impl_->operators.accumulate_mass_residual(final_mass_flux.impl_->handle,
                                              mass_residual);
    candidate_density.reserve(impl_->topology->owned_cell_count());
    for_each_owned_cell(*impl_->topology, [&](mesh::LocalCellId cell,
                                              runtime::Int3 index) {
      const double volume = impl_->geometry->cell_volume_m3(cell);
      const double candidate =
          -(stencil.alpha1 * rho_n(index.x, index.y, index.z, 0) +
            stencil.alpha2 * rho_h(index.x, index.y, index.z, 0) +
            stencil.dt_s * mass_residual(index.x, index.y, index.z, 0) /
                volume) /
          stencil.alpha0;
      if (!std::isfinite(candidate))
        throw runtime::Error("material density candidate is non-finite");
      if (!(candidate > 0.0))
        throw std::domain_error("material density candidate is non-positive");
      candidate_density.push_back(candidate);
    });
  } catch (const std::domain_error &) {
    local = MaterialTransportFailureReason::non_positive_density;
  } catch (const runtime::MpiOperationError &) {
    throw;
  } catch (const std::exception &) {
    local = MaterialTransportFailureReason::non_finite_state;
  }
  failure = synchronize_failure(*impl_->mpi, local);
  if (failure.reason != MaterialTransportFailureReason::none) {
    report.disposition_ = MaterialTransportDisposition::recoverable_failure;
    report.reason_ = failure.reason;
    report.lowest_failing_rank_ = failure.rank;
    return finish_report();
  }

  for (std::size_t field = 0; field < fields.size(); ++field) {
    local = MaterialTransportFailureReason::none;
    std::optional<TransportResidual> residual;
    std::vector<double> candidate_values;
    bool failure_already_synchronized = false;
    try {
      residual.emplace(compute_transport_residual(
          *impl_, final_mass_flux.impl_->handle, state, fields[field], field));
      const auto q_n = state.layer(FlowLayer::committed)
                           .acquire_read<double>(impl_->access, kMaterialPhase,
                                                 kMaterialActor, fields[field]);
      const auto q_h = state.layer(FlowLayer::history)
                           .acquire_read<double>(impl_->access, kMaterialPhase,
                                                 kMaterialActor, fields[field]);
      std::size_t owned = 0U;
      for_each_owned_cell(*impl_->topology, [&](mesh::LocalCellId cell,
                                                runtime::Int3 index) {
        const double spatial =
            stencil.order == MomentumTimeOrder::backward_euler
                ? residual->current[owned]
                : 2.0 * residual->current[owned] - residual->history[owned];
        const double candidate =
            -(stencil.alpha1 * q_n(index.x, index.y, index.z, 0) +
              stencil.alpha2 * q_h(index.x, index.y, index.z, 0) +
              stencil.dt_s * spatial / impl_->geometry->cell_volume_m3(cell)) /
            stencil.alpha0;
        if (!std::isfinite(candidate))
          throw runtime::Error("material transport candidate is non-finite");
        candidate_values.push_back(candidate);
        ++owned;
      });
    } catch (const TransportComputationFailure &error) {
      failure = error.failure();
      failure_already_synchronized = true;
    } catch (const runtime::MpiOperationError &) {
      throw;
    } catch (const std::exception &) {
      local = MaterialTransportFailureReason::non_finite_state;
    }
    if (!failure_already_synchronized)
      failure = synchronize_failure(*impl_->mpi, local);
    if (failure.reason != MaterialTransportFailureReason::none) {
      report.disposition_ = MaterialTransportDisposition::recoverable_failure;
      report.reason_ = failure.reason;
      report.lowest_failing_rank_ = failure.rank;
      return finish_report();
    }
    residuals.push_back(std::move(*residual));
    candidates[field] = std::move(candidate_values);
  }

  auto trial_rho = state.trial_layer().acquire_write<double>(
      impl_->access, kMaterialPhase, kMaterialActor, state.fields().density);
  std::size_t owned = 0U;
  for_each_owned_cell(
      *impl_->topology, [&](mesh::LocalCellId, runtime::Int3 index) {
        trial_rho(index.x, index.y, index.z, 0) = candidate_density[owned++];
      });
  for (std::size_t field = 0; field < fields.size(); ++field) {
    auto trial = state.trial_layer().acquire_write<double>(
        impl_->access, kMaterialPhase, kMaterialActor, fields[field]);
    owned = 0U;
    for_each_owned_cell(
        *impl_->topology, [&](mesh::LocalCellId, runtime::Int3 index) {
          trial(index.x, index.y, index.z, 0) = candidates[field][owned++];
        });
  }
  auto trial_flux = state.trial_layer().acquire_face_write<double>(
      impl_->access, kMaterialPhase, kMaterialActor,
      state.fields().face_mass_flux);
  for (std::size_t face = 0; face < impl_->topology->local_face_count(); ++face)
    trial_flux(face, 0) = final_mass_flux.impl_->view(face, 0);

  constexpr std::size_t mass_terms = 9U;
  constexpr std::size_t field_terms = 7U;
  std::vector<CompensatedSum> local_sums(mass_terms +
                                         field_terms * fields.size());
  double local_min = std::numeric_limits<double>::infinity();
  mesh::GlobalCellId local_min_id = 0U;
  owned = 0U;
  for_each_owned_cell(
      *impl_->topology, [&](mesh::LocalCellId cell, runtime::Int3 index) {
        const double volume = impl_->geometry->cell_volume_m3(cell);
        const double rn = rho_n(index.x, index.y, index.z, 0);
        const double rh = rho_h(index.x, index.y, index.z, 0);
        const double next = candidate_density[owned];
        const double residual =
            stencil.alpha0 * next + stencil.alpha1 * rn + stencil.alpha2 * rh +
            stencil.dt_s * mass_residual(index.x, index.y, index.z, 0) / volume;
        local_sums[0].add(residual * residual);
        local_sums[1].add(next * next);
        local_sums[2].add(volume * rn);
        local_sums[3].add(volume * rh);
        local_sums[4].add(volume * next);
        local_sums[6].add(std::abs(volume * rn));
        local_sums[7].add(std::abs(volume * rh));
        local_sums[8].add(std::abs(volume * next));
        if (next < local_min ||
            (next == local_min &&
             impl_->topology->global_cell_id(cell) < local_min_id)) {
          local_min = next;
          local_min_id = impl_->topology->global_cell_id(cell);
        }
        ++owned;
      });
  CompensatedSum boundary_mass;
  CompensatedSum boundary_mass_abs;
  for (mesh::LocalFaceId face = 0; face < impl_->topology->local_face_count();
       ++face) {
    const auto patch = impl_->topology->patch_id(face);
    if (!patch ||
        impl_->boundaries->patch(*patch).kind() ==
            boundary::BoundaryKind::periodic ||
        impl_->topology->cell_ownership(impl_->topology->owner(face)) !=
            mesh::EntityOwnership::owned)
      continue;
    const double value = final_mass_flux.impl_->view(face, 0);
    boundary_mass.add(value);
    boundary_mass_abs.add(std::abs(value));
  }
  local_sums[5].add(boundary_mass.value());
  for (std::size_t field = 0; field < fields.size(); ++field) {
    const auto qn = state.layer(FlowLayer::committed)
                        .acquire_read<double>(impl_->access, kMaterialPhase,
                                              kMaterialActor, fields[field]);
    const auto qh = state.layer(FlowLayer::history)
                        .acquire_read<double>(impl_->access, kMaterialPhase,
                                              kMaterialActor, fields[field]);
    owned = 0U;
    const std::size_t offset = mass_terms + field_terms * field;
    for_each_owned_cell(*impl_->topology, [&](mesh::LocalCellId cell,
                                              runtime::Int3 index) {
      const double volume = impl_->geometry->cell_volume_m3(cell);
      const double spatial = stencil.order == MomentumTimeOrder::backward_euler
                                 ? residuals[field].current[owned]
                                 : 2.0 * residuals[field].current[owned] -
                                       residuals[field].history[owned];
      const double equation =
          stencil.alpha0 * candidates[field][owned] +
          stencil.alpha1 * qn(index.x, index.y, index.z, 0) +
          stencil.alpha2 * qh(index.x, index.y, index.z, 0) +
          stencil.dt_s * spatial / volume;
      const double current = qn(index.x, index.y, index.z, 0);
      const double history = qh(index.x, index.y, index.z, 0);
      const double next = candidates[field][owned];
      local_sums[offset].add(equation * equation);
      local_sums[offset + 1U].add(next * next);
      local_sums[offset + 2U].add(volume * current);
      local_sums[offset + 3U].add(volume * history);
      local_sums[offset + 4U].add(std::abs(volume * current));
      local_sums[offset + 5U].add(std::abs(volume * history));
      local_sums[offset + 6U].add(std::abs(volume * next));
      ++owned;
    });
  }
  std::vector<double> sums(local_sums.size());
  std::transform(local_sums.begin(), local_sums.end(), sums.begin(),
                 [](const CompensatedSum &sum) { return sum.value(); });
  double boundary_mass_abs_value = boundary_mass_abs.value();
  local = std::all_of(sums.begin(), sums.end(),
                      [](double value) { return std::isfinite(value); }) &&
                  std::isfinite(boundary_mass_abs_value)
              ? MaterialTransportFailureReason::none
              : MaterialTransportFailureReason::non_finite_state;
  failure = synchronize_failure(*impl_->mpi, local);
  if (failure.reason != MaterialTransportFailureReason::none) {
    report.disposition_ = MaterialTransportDisposition::recoverable_failure;
    report.reason_ = failure.reason;
    report.lowest_failing_rank_ = failure.rank;
    return finish_report();
  }
  impl_->mpi->allreduce_fp64_in_place(sums.data(), sums.size(),
                                      runtime::Fp64ReductionOperation::sum);
  impl_->mpi->allreduce_fp64_in_place(&boundary_mass_abs_value, 1U,
                                      runtime::Fp64ReductionOperation::sum);
  report.density_normalized_l2_ = std::sqrt(
      sums[0] / std::max(sums[1], std::numeric_limits<double>::min()));
  report.density_residual_available_ = true;
  report.mass_relative_conservation_defect_ = relative_defect(
      sums[4] - sums[2] +
          (stencil.alpha2 / stencil.alpha0) * (sums[3] - sums[2]) +
          (stencil.dt_s / stencil.alpha0) * sums[5],
      conservation_denominator(
          sums[2], sums[4], sums[3], stencil.alpha2 / stencil.alpha0,
          (stencil.dt_s / stencil.alpha0) * boundary_mass_abs_value, sums[6],
          sums[8], stencil.order == MomentumTimeOrder::bdf2 ? sums[7] : 0.0));
  for (std::size_t field = 0; field < fields.size(); ++field) {
    const std::size_t offset = mass_terms + field_terms * field;
    report.transport_normalized_l2_[field] =
        std::sqrt(sums[offset] / std::max(sums[offset + 1U],
                                          std::numeric_limits<double>::min()));
    report.transport_residual_available_[field] = 1U;
    CompensatedSum next_integral_sum;
    owned = 0U;
    for_each_owned_cell(
        *impl_->topology, [&](mesh::LocalCellId cell, runtime::Int3) {
          next_integral_sum.add(impl_->geometry->cell_volume_m3(cell) *
                                candidates[field][owned++]);
        });
    double next_integral = next_integral_sum.value();
    impl_->mpi->allreduce_fp64_in_place(&next_integral, 1U,
                                        runtime::Fp64ReductionOperation::sum);
    const double effective_boundary =
        stencil.order == MomentumTimeOrder::backward_euler
            ? residuals[field].boundary_current
            : 2.0 * residuals[field].boundary_current -
                  residuals[field].boundary_history;
    double boundary_values[3]{effective_boundary,
                              residuals[field].boundary_abs_current,
                              residuals[field].boundary_abs_history};
    local = std::all_of(boundary_values, boundary_values + 3,
                        [](double value) { return std::isfinite(value); })
                ? MaterialTransportFailureReason::none
                : MaterialTransportFailureReason::non_finite_state;
    failure = synchronize_failure(*impl_->mpi, local);
    if (failure.reason != MaterialTransportFailureReason::none) {
      report.disposition_ = MaterialTransportDisposition::recoverable_failure;
      report.reason_ = failure.reason;
      report.lowest_failing_rank_ = failure.rank;
      return finish_report();
    }
    impl_->mpi->allreduce_fp64_in_place(boundary_values, 3U,
                                        runtime::Fp64ReductionOperation::sum);
    const double raw = next_integral - sums[offset + 2U] +
                       (stencil.alpha2 / stencil.alpha0) *
                           (sums[offset + 3U] - sums[offset + 2U]) +
                       (stencil.dt_s / stencil.alpha0) * boundary_values[0];
    const double boundary_scale =
        stencil.order == MomentumTimeOrder::backward_euler
            ? boundary_values[1]
            : 2.0 * boundary_values[1] + boundary_values[2];
    report.transport_relative_conservation_defect_[field] = relative_defect(
        raw, conservation_denominator(
                 sums[offset + 2U], next_integral, sums[offset + 3U],
                 stencil.alpha2 / stencil.alpha0,
                 (stencil.dt_s / stencil.alpha0) * boundary_scale,
                 sums[offset + 4U], sums[offset + 6U],
                 stencil.order == MomentumTimeOrder::bdf2 ? sums[offset + 5U]
                                                          : 0.0));
    report.transport_conservation_available_[field] = 1U;
  }
  report.mass_conservation_available_ = true;

  double global_min = local_min;
  runtime::check_mpi_result(MPI_Allreduce(MPI_IN_PLACE, &global_min, 1,
                                          MPI_DOUBLE, MPI_MIN,
                                          impl_->mpi->comm()),
                            "MPI_Allreduce(material minimum density)");
  int minimum_rank =
      local_min == global_min ? impl_->mpi->rank() : impl_->mpi->size();
  runtime::check_mpi_result(MPI_Allreduce(MPI_IN_PLACE, &minimum_rank, 1,
                                          MPI_INT, MPI_MIN, impl_->mpi->comm()),
                            "MPI_Allreduce(material minimum rank)");
  std::uint64_t minimum_id =
      impl_->mpi->rank() == minimum_rank ? local_min_id : 0U;
  runtime::check_mpi_result(
      MPI_Bcast(&minimum_id, 1, MPI_UINT64_T, minimum_rank, impl_->mpi->comm()),
      "MPI_Bcast(material minimum cell)");
  report.minimum_density_kg_per_m3_ = global_min;
  report.minimum_density_rank_ = minimum_rank;
  report.minimum_density_global_cell_ = minimum_id;
  report.minimum_density_available_ = true;

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  const auto applies = [&](const GateOverride &item) {
    return item.value.has_value() &&
           (item.rank < 0 || item.rank == impl_->mpi->rank());
  };
  if (applies(overrides.density_residual))
    report.density_normalized_l2_ = *overrides.density_residual.value;
  if (applies(overrides.transport_residual) &&
      overrides.transport_residual.field <
          report.transport_normalized_l2_.size())
    report.transport_normalized_l2_[overrides.transport_residual.field] =
        *overrides.transport_residual.value;
  if (applies(overrides.mass_conservation))
    report.mass_relative_conservation_defect_ =
        *overrides.mass_conservation.value;
  if (applies(overrides.transport_conservation) &&
      overrides.transport_conservation.field <
          report.transport_relative_conservation_defect_.size())
    report.transport_relative_conservation_defect_
        [overrides.transport_conservation.field] =
        *overrides.transport_conservation.value;
#endif

  local = MaterialTransportFailureReason::none;
  if (!std::isfinite(report.density_normalized_l2_) ||
      report.density_normalized_l2_ > 1.0e-10)
    local = MaterialTransportFailureReason::final_density_residual;
  if (local == MaterialTransportFailureReason::none) {
    for (double residual : report.transport_normalized_l2_) {
      if (!std::isfinite(residual) || residual > 1.0e-9) {
        local = MaterialTransportFailureReason::final_transport_residual;
        break;
      }
    }
  }
  if (local == MaterialTransportFailureReason::none &&
      (!std::isfinite(report.mass_relative_conservation_defect_) ||
       report.mass_relative_conservation_defect_ > 5.0e-11))
    local = MaterialTransportFailureReason::final_conservation_defect;
  if (local == MaterialTransportFailureReason::none) {
    for (double defect : report.transport_relative_conservation_defect_) {
      if (!std::isfinite(defect) || defect > 5.0e-11) {
        local = MaterialTransportFailureReason::final_conservation_defect;
        break;
      }
    }
  }
  failure = synchronize_failure(*impl_->mpi, local);
  if (failure.reason == MaterialTransportFailureReason::none) {
    report.disposition_ = MaterialTransportDisposition::finalized;
    report.reason_ = MaterialTransportFailureReason::none;
    report.lowest_failing_rank_ = -1;
  } else {
    report.disposition_ = MaterialTransportDisposition::recoverable_failure;
    report.reason_ = failure.reason;
    report.lowest_failing_rank_ = failure.rank;
  }
  return finish_report();
}

struct MaterialDensityDiagnosticSource::Impl final {
  Impl(const MaterialDensityTransport::Impl &supplied_transport,
       const FlowState &supplied_state,
       MaterialDensityTransportReport supplied_report)
      : transport(&supplied_transport), state(&supplied_state),
        report(std::move(supplied_report)) {}
  const MaterialDensityTransport::Impl *transport{};
  const FlowState *state{};
  MaterialDensityTransportReport report;
};

MaterialDensityDiagnosticSource MaterialDensityTransport::diagnostic_source(
    const FlowState &state,
    const MaterialDensityTransportReport &report) const {
  if (impl_->last_state != &state || !state.attempt_active() ||
      !report.authenticated() || impl_->last_report_seal != report.seal_ ||
      state.attempt_identity() != report.attempt_identity_ ||
      impl_->last_attempt_identity != report.attempt_identity_ ||
      impl_->last_finalization_identity != report.finalization_identity_ ||
      report.finalization_identity_ == 0U)
    throw runtime::Error("material diagnostic source identity is stale");
  auto source = std::make_unique<MaterialDensityDiagnosticSource::Impl>(
      *impl_, state, report);
  return MaterialDensityDiagnosticSource(std::move(source));
}

MaterialDensityDiagnosticSource::MaterialDensityDiagnosticSource(
    std::unique_ptr<Impl> impl) noexcept
    : impl_(std::move(impl)) {}
MaterialDensityDiagnosticSource::~MaterialDensityDiagnosticSource() noexcept =
    default;
MaterialDensityDiagnosticSource::MaterialDensityDiagnosticSource(
    MaterialDensityDiagnosticSource &&) noexcept = default;

namespace {
template <class SourceImplementation>
void validate_source(const SourceImplementation &source, bool attempt_active,
                     std::uint64_t attempt_identity, bool report_authenticated,
                     std::uint64_t report_seal,
                     std::uint64_t report_attempt_identity,
                     std::uint64_t report_finalization_identity) {
  if (source.transport->last_state != source.state || !attempt_active ||
      !report_authenticated ||
      source.transport->last_report_seal != report_seal ||
      attempt_identity != report_attempt_identity ||
      source.transport->last_attempt_identity != report_attempt_identity ||
      source.transport->last_finalization_identity !=
          report_finalization_identity)
    throw runtime::Error("material diagnostic source is stale");
}
} // namespace

#define HUNDUN_VALIDATE_MATERIAL_SOURCE()                                      \
  validate_source(*impl_, impl_->state->attempt_active(),                      \
                  impl_->state->attempt_identity(),                            \
                  impl_->report.authenticated(), impl_->report.seal_,          \
                  impl_->report.attempt_identity_,                             \
                  impl_->report.finalization_identity_)

std::size_t MaterialDensityDiagnosticSource::fingerprint_field_count() const {
  HUNDUN_VALIDATE_MATERIAL_SOURCE();
  return impl_->transport->fingerprint_ids.size();
}
std::string_view
MaterialDensityDiagnosticSource::fingerprint_field_id(std::size_t index) const {
  HUNDUN_VALIDATE_MATERIAL_SOURCE();
  if (index >= impl_->transport->fingerprint_ids.size())
    throw runtime::Error("material diagnostic field index is invalid");
  return impl_->transport->fingerprint_ids[index];
}
std::string_view
MaterialDensityDiagnosticSource::field_unit(std::size_t index) const {
  HUNDUN_VALIDATE_MATERIAL_SOURCE();
  if (index >= impl_->transport->fingerprint_ids.size())
    throw runtime::Error("material diagnostic field index is invalid");
  return index == 0U ? "kg/m3" : (index == 1U ? "J/m3" : "kg/m3");
}
std::string_view
MaterialDensityDiagnosticSource::owned_cell_layout_fingerprint() const {
  HUNDUN_VALIDATE_MATERIAL_SOURCE();
  return impl_->transport->owned_layout_fingerprint;
}
std::string_view
MaterialDensityDiagnosticSource::global_cell_layout_fingerprint() const {
  HUNDUN_VALIDATE_MATERIAL_SOURCE();
  return impl_->transport->global_layout_fingerprint;
}
std::size_t MaterialDensityDiagnosticSource::owned_cell_count() const {
  HUNDUN_VALIDATE_MATERIAL_SOURCE();
  return impl_->transport->topology->owned_cell_count();
}
mesh::GlobalCellId
MaterialDensityDiagnosticSource::global_cell_id(std::size_t local_cell) const {
  HUNDUN_VALIDATE_MATERIAL_SOURCE();
  if (local_cell >= impl_->transport->topology->owned_cell_count())
    throw runtime::Error("material diagnostic cell index is invalid");
  return impl_->transport->topology->global_cell_id(local_cell);
}
double
MaterialDensityDiagnosticSource::cell_volume_m3(std::size_t local_cell) const {
  HUNDUN_VALIDATE_MATERIAL_SOURCE();
  if (local_cell >= impl_->transport->topology->owned_cell_count())
    throw runtime::Error("material diagnostic cell index is invalid");
  return impl_->transport->geometry->cell_volume_m3(local_cell);
}
double
MaterialDensityDiagnosticSource::field_value(std::size_t field,
                                             std::size_t local_cell) const {
  HUNDUN_VALIDATE_MATERIAL_SOURCE();
  if (field >= impl_->transport->fingerprint_ids.size() ||
      local_cell >= impl_->transport->topology->owned_cell_count())
    throw runtime::Error("material diagnostic state index is invalid");
  const auto field_id =
      field == 0U
          ? impl_->state->fields().density
          : (field == 1U ? impl_->transport->specification.enthalpy_density
                         : impl_->transport->specification
                               .scalar_densities[field - 2U]);
  const auto view =
      impl_->state->layer(FlowLayer::trial)
          .acquire_read<double>(impl_->transport->access, kMaterialPhase,
                                kMaterialActor, field_id);
  const auto index = local_index(*impl_->transport->topology, local_cell);
  return view(index.x, index.y, index.z, 0);
}
const MaterialDensityTransportReport &
MaterialDensityDiagnosticSource::report() const {
  HUNDUN_VALIDATE_MATERIAL_SOURCE();
  return impl_->report;
}

#undef HUNDUN_VALIDATE_MATERIAL_SOURCE

} // namespace hundun::flow
