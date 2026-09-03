// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/fvm_immersed_reconstruction.hpp"

#include "hundun/bc_basic_boundary.hpp"
#include "hundun/ib_ghost_stencil_plan.hpp"
#include "hundun/ib_domain.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_collective_status.hpp"
#include "hundun/rt_error.hpp"
#include "hundun/rt_halo_exchange.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "fvm_immersed_boundary_authority_detail.hpp"
#include "ib_quadratic_reconstruction_detail.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
#include <mutex>
#include <unordered_map>
#endif

namespace hundun::finite_volume {
namespace {

using immersed::CellRegion;
using mesh::EntityOwnership;
using mesh::LocalCellId;
using mesh::LocalFaceId;
using runtime::Error;
using runtime::Int3;
using runtime::Real3;

constexpr std::uint64_t kInvalidLink =
    std::numeric_limits<std::uint64_t>::max();
constexpr std::size_t kInvalidFace = std::numeric_limits<std::size_t>::max();
constexpr std::uint64_t kFnvOffset = UINT64_C(14695981039346656037);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);

thread_local const std::vector<detail::ImmersedWallNormalGradient>
    *scoped_wall_normal_gradients = nullptr;

class ScopedWallNormalGradients final {
public:
  explicit ScopedWallNormalGradients(
      const std::vector<detail::ImmersedWallNormalGradient> &gradients) {
    if (scoped_wall_normal_gradients != nullptr)
      throw Error(
          "immersed reconstruction wall authority scope is already active");
    scoped_wall_normal_gradients = &gradients;
  }
  ~ScopedWallNormalGradients() { scoped_wall_normal_gradients = nullptr; }

  ScopedWallNormalGradients(const ScopedWallNormalGradients &) = delete;
  ScopedWallNormalGradients &
  operator=(const ScopedWallNormalGradients &) = delete;
};

#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
thread_local std::vector<std::uint8_t> execution_trace;
thread_local std::size_t inactive_attempts{};
thread_local bool reject_inactive_reads{};
thread_local bool force_momentum_fallback_for_test{};
std::mutex reconstruction_snapshots_mutex;
std::unordered_map<const void *, std::vector<std::uint64_t>>
    reconstruction_snapshots;

void register_reconstruction_snapshot(
    const void *identity, std::vector<std::uint64_t> snapshot) {
  const std::lock_guard<std::mutex> lock(reconstruction_snapshots_mutex);
  const auto [unused, inserted] =
      reconstruction_snapshots.emplace(identity, std::move(snapshot));
  static_cast<void>(unused);
  if (!inserted)
    throw Error("immersed reconstruction snapshot identity is already active");
}

void erase_reconstruction_snapshot(const void *identity) noexcept {
  try {
    const std::lock_guard<std::mutex> lock(reconstruction_snapshots_mutex);
    reconstruction_snapshots.erase(identity);
  } catch (...) {
    std::terminate();
  }
}

class ReconstructionSnapshotRollback final {
public:
  explicit ReconstructionSnapshotRollback(const void *identity) noexcept
      : identity_(identity) {}

  ~ReconstructionSnapshotRollback() noexcept {
    if (active_)
      erase_reconstruction_snapshot(identity_);
  }

  void arm() noexcept { active_ = true; }

  void rollback() noexcept {
    if (active_) {
      erase_reconstruction_snapshot(identity_);
      active_ = false;
    }
  }

  void release() noexcept { active_ = false; }

  ReconstructionSnapshotRollback(const ReconstructionSnapshotRollback &) =
      delete;
  ReconstructionSnapshotRollback &
  operator=(const ReconstructionSnapshotRollback &) = delete;

private:
  const void *identity_{};
  bool active_{};
};

std::vector<std::uint64_t>
lookup_reconstruction_snapshot(const void *identity) {
  const std::lock_guard<std::mutex> lock(reconstruction_snapshots_mutex);
  const auto found = reconstruction_snapshots.find(identity);
  return found == reconstruction_snapshots.end() ? std::vector<std::uint64_t>{}
                                                  : found->second;
}
#endif

#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
void trace_stage(std::uint8_t stage) {
  execution_trace.push_back(stage);
}
#endif

bool same(Int3 lhs, Int3 rhs) noexcept {
  return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

bool same(runtime::Box3 lhs, runtime::Box3 rhs) noexcept {
  return same(lhs.begin, rhs.begin) && same(lhs.end, rhs.end);
}

bool finite(Real3 value) noexcept {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

Real3 add(Real3 lhs, Real3 rhs) noexcept {
  return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

Real3 subtract(Real3 lhs, Real3 rhs) noexcept {
  return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

Real3 multiply(double scale, Real3 value) noexcept {
  return {scale * value.x, scale * value.y, scale * value.z};
}

double dot(Real3 lhs, Real3 rhs) noexcept {
  return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

void hash_u64(std::uint64_t &hash, std::uint64_t value) noexcept {
  for (unsigned shift = 0U; shift < 64U; shift += 8U) {
    hash ^= static_cast<std::uint8_t>(value >> shift);
    hash *= kFnvPrime;
  }
}

void check_mpi(int code, const char *operation) {
  if (code == MPI_SUCCESS)
    return;
  char text[MPI_MAX_ERROR_STRING]{};
  int length = 0;
  MPI_Error_string(code, text, &length);
  throw Error(std::string(operation) + ": " +
              std::string(text, static_cast<std::size_t>(length)));
}

[[noreturn]] void throw_collective(const runtime::CollectiveStatus &status) {
  throw Error(status.message + " (lowest failing rank " +
              std::to_string(status.failing_rank) + ")");
}

void require_collective(const runtime::MpiContext &mpi, bool local_ok,
                        const std::string &message) {
  const auto status = runtime::collective_status(mpi, local_ok, message);
  if (!status.ok)
    throw_collective(status);
}

struct GlobalLink final {
  std::uint64_t id{};
  mesh::GlobalCellId fluid{};
  mesh::GlobalCellId solid{};
};

struct GlobalLinkPair final {
  mesh::GlobalCellId first{};
  mesh::GlobalCellId second{};
  std::uint64_t id{};
};

std::vector<std::uint64_t>
allgather_u64(const std::vector<std::uint64_t> &local,
              const runtime::MpiContext &mpi) {
  const bool count_ok =
      local.size() <= static_cast<std::size_t>(std::numeric_limits<int>::max());
  require_collective(mpi, count_ok,
                     "immersed reconstruction gather count is unsupported");
  const int local_count = static_cast<int>(local.size());
  std::vector<int> counts(static_cast<std::size_t>(mpi.size()));
  check_mpi(MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT,
                          mpi.comm()),
            "MPI_Allgather immersed reconstruction counts");
  std::vector<int> displacements(counts.size(), 0);
  std::uint64_t total = 0U;
  bool extent_ok = true;
  for (std::size_t rank = 0U; rank < counts.size(); ++rank) {
    if (counts[rank] < 0 ||
        total > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
      extent_ok = false;
      break;
    }
    displacements[rank] = static_cast<int>(total);
    total += static_cast<std::uint64_t>(counts[rank]);
  }
  extent_ok = extent_ok && total <= static_cast<std::uint64_t>(
                                        std::numeric_limits<int>::max());
  require_collective(mpi, extent_ok,
                     "immersed reconstruction gather extent is unsupported");
  std::vector<std::uint64_t> gathered(static_cast<std::size_t>(total));
  check_mpi(MPI_Allgatherv(local.data(), local_count, MPI_UINT64_T,
                           gathered.data(), counts.data(), displacements.data(),
                           MPI_UINT64_T, mpi.comm()),
            "MPI_Allgatherv immersed reconstruction records");
  return gathered;
}

std::vector<GlobalLink> gather_links(const immersed::ImmersedDomain &domain,
                                     const runtime::MpiContext &mpi) {
  std::vector<std::uint64_t> encoded;
  encoded.reserve(domain.links().size() * 3U);
  for (const auto &link : domain.links()) {
    encoded.push_back(link.id);
    encoded.push_back(link.fluid_cell);
    encoded.push_back(link.solid_cell);
  }
  const auto gathered = allgather_u64(encoded, mpi);
  require_collective(mpi, gathered.size() % 3U == 0U,
                     "immersed reconstruction link record is malformed");
  std::vector<GlobalLink> result;
  result.reserve(gathered.size() / 3U);
  for (std::size_t i = 0U; i < gathered.size(); i += 3U)
    result.push_back({gathered[i], gathered[i + 1U], gathered[i + 2U]});
  std::sort(result.begin(), result.end(),
            [](const auto &a, const auto &b) { return a.id < b.id; });
  bool valid = !result.empty();
  for (std::size_t i = 0U; valid && i < result.size(); ++i)
    valid = result[i].id == i && result[i].fluid != result[i].solid;
  require_collective(mpi, valid,
                     "immersed reconstruction link inventory is invalid");
  return result;
}

std::vector<std::uint8_t>
gather_global_active(const immersed::ImmersedDomain &domain,
                     const mesh::MeshTopology &topology,
                     const runtime::MpiContext &mpi) {
  std::vector<std::uint64_t> local;
  local.reserve(domain.active_cells().owned_active_count());
  for (LocalCellId cell = 0U; cell < topology.owned_cell_count(); ++cell)
    if (domain.region(cell) == CellRegion::fluid)
      local.push_back(topology.global_cell_id(cell));
  const auto gathered = allgather_u64(local, mpi);
  std::vector<std::uint8_t> active(
      static_cast<std::size_t>(topology.global_cell_count()), 0U);
  bool valid = true;
  for (const auto id : gathered) {
    if (id >= active.size() || active[static_cast<std::size_t>(id)] != 0U) {
      valid = false;
      break;
    }
    active[static_cast<std::size_t>(id)] = 1U;
  }
  require_collective(
      mpi, valid, "immersed reconstruction active-cell inventory is invalid");
  return active;
}

void validate_quantity(FiniteVolumeQuantity quantity,
                       const boundary::BoundaryRegistry &boundaries,
                       bool allow_velocity) {
  switch (quantity.kind) {
  case FiniteVolumeQuantityKind::density:
  case FiniteVolumeQuantityKind::enthalpy:
    if (quantity.scalar_index != 0U)
      throw Error("non-scalar finite-volume quantity has a scalar index");
    return;
  case FiniteVolumeQuantityKind::scalar:
    if (quantity.scalar_index >= boundaries.scalar_count())
      throw Error("finite-volume scalar index is out of bounds");
    return;
  case FiniteVolumeQuantityKind::velocity:
    if (!allow_velocity || quantity.scalar_index != 0U)
      throw Error("finite-volume quantity is invalid for this operation");
    return;
  case FiniteVolumeQuantityKind::pressure:
    if (!allow_velocity || quantity.scalar_index != 0U)
      throw Error("finite-volume quantity is invalid for this operation");
    return;
  }
  throw Error("invalid finite-volume quantity kind");
}

double boundary_scalar(const boundary::BoundaryRegistry &boundaries,
                       FiniteVolumeQuantity quantity, std::uint32_t patch,
                       double interior, bool exterior) {
  boundary::ScalarBoundaryValues values{};
  switch (quantity.kind) {
  case FiniteVolumeQuantityKind::density:
    values = boundaries.evaluate_density(patch, interior);
    break;
  case FiniteVolumeQuantityKind::enthalpy:
    values = boundaries.evaluate_enthalpy(patch, interior);
    break;
  case FiniteVolumeQuantityKind::scalar:
    values = boundaries.evaluate_scalar(patch, quantity.scalar_index, interior);
    break;
  case FiniteVolumeQuantityKind::pressure:
    values = boundaries.evaluate_pressure(patch, interior);
    break;
  default:
    throw Error("scalar boundary evaluation received velocity");
  }
  return exterior ? values.exterior : values.face;
}

struct CellFaceRef final {
  LocalFaceId face{kInvalidFace};
  bool as_owner{};
};

struct CanonicalFaceLine final {
  LocalCellId p{};
  LocalCellId n{};
  Int3 pm1{};
  Int3 np1{};
  bool reversed{};
};

struct OperationGuard final {
  explicit OperationGuard(bool &active) : active_(&active) {
    if (active)
      throw Error("immersed reconstruction operation is already active");
    active = true;
  }
  ~OperationGuard() { *active_ = false; }
  bool *active_;
};

std::string exception_message() {
  try {
    throw;
  } catch (const Error &error) {
    return error.what();
  } catch (const std::exception &error) {
    return error.what();
  } catch (...) {
    return "immersed reconstruction operation failed";
  }
}

} // namespace

struct ImmersedReconstruction::Impl final {
  const mesh::MeshTopology *topology{};
  const mesh::MeshGeometry *geometry{};
  const boundary::BoundaryRegistry *boundaries{};
  const immersed::ImmersedDomain *domain{};
  const immersed::GhostStencilPlan *ghost_plan{};
  const runtime::StructuredDecomposition *decomposition{};
  const runtime::MpiContext *mpi{};
  runtime::HaloExchange *halo{};
  runtime::Box3 owned_box{};
  Int3 local_extent{};
  Int3 global_extent{};
  std::uint32_t halo_reach{};
  std::vector<GlobalLink> links;
  std::vector<GlobalLinkPair> link_pairs;
  std::vector<immersed::QuadraticReconstruction> reconstructions;
  std::vector<std::uint8_t> global_active;
  std::vector<std::uint8_t> local_active;
  std::vector<std::uint8_t> local_interface;
  std::vector<std::uint64_t> local_cell_link;
  std::vector<std::uint64_t> owned_row_link;
  std::vector<std::vector<std::uint64_t>> owned_wall_links;
  std::vector<LocalCellId> interior_rows;
  std::vector<LocalCellId> partition_rows;
  std::vector<LocalCellId> interface_rows;
  std::vector<LocalFaceId> interior_faces;
  std::vector<LocalFaceId> partition_faces;
  std::vector<LocalFaceId> interface_faces;
  std::vector<std::uint64_t> face_p_ghost_link;
  std::vector<std::uint64_t> face_n_ghost_link;
  std::vector<std::array<CellFaceRef, 6>> cell_faces;
  mutable std::vector<double> cell_scratch;
  mutable std::vector<double> face_scratch;
  mutable std::vector<double> face_p_ghost_symbols;
  mutable std::vector<double> face_n_ghost_symbols;
  mutable std::vector<double> wall_normal_gradient_scratch;
  mutable bool active{};
  std::uint64_t fingerprint{};

  Int3 global(LocalCellId cell) const { return topology->global_cell(cell); }

  Int3 local_index(LocalCellId cell) const { return field_index(global(cell)); }

  mesh::GlobalCellId global_id(Int3 point) const {
    return topology->global_cell_id(point);
  }

  Int3 logical(mesh::GlobalCellId cell) const {
    if (cell >= topology->global_cell_count())
      throw Error("immersed reconstruction global cell is out of range");
    const auto nx = static_cast<std::uint64_t>(global_extent.x);
    const auto ny = static_cast<std::uint64_t>(global_extent.y);
    const auto plane = nx * ny;
    return {static_cast<int>(cell % nx), static_cast<int>((cell / nx) % ny),
            static_cast<int>(cell / plane)};
  }

  Int3 wrap(Int3 point) const {
    const auto periodic = decomposition->periodic();
    int *coordinates[3]{&point.x, &point.y, &point.z};
    const int extents[3]{global_extent.x, global_extent.y, global_extent.z};
    for (std::size_t axis = 0U; axis < 3U; ++axis) {
      if (*coordinates[axis] < 0 || *coordinates[axis] >= extents[axis]) {
        if (!periodic[axis])
          throw Error("finite-volume stencil leaves a nonperiodic mesh");
        *coordinates[axis] %= extents[axis];
        if (*coordinates[axis] < 0)
          *coordinates[axis] += extents[axis];
      }
    }
    return point;
  }

  Int3 field_index(Int3 global_point) const {
    const auto periodic = decomposition->periodic();
    const int coordinates[3]{global_point.x, global_point.y, global_point.z};
    const int begins[3]{owned_box.begin.x, owned_box.begin.y,
                        owned_box.begin.z};
    const int local_sizes[3]{local_extent.x, local_extent.y, local_extent.z};
    const int global_sizes[3]{global_extent.x, global_extent.y,
                              global_extent.z};
    int result[3]{};
    const int reach = static_cast<int>(halo_reach);
    for (std::size_t axis = 0U; axis < 3U; ++axis) {
      result[axis] = coordinates[axis] - begins[axis];
      if (!periodic[axis])
        continue;
      const int lower_image = result[axis] - global_sizes[axis];
      const int upper_image = result[axis] + global_sizes[axis];
      const auto in_ghosted = [&](int value) {
        return value >= -reach && value < local_sizes[axis] + reach;
      };
      if (!in_ghosted(result[axis]) && in_ghosted(lower_image))
        result[axis] = lower_image;
      else if (!in_ghosted(result[axis]) && in_ghosted(upper_image))
        result[axis] = upper_image;
    }
    return {result[0], result[1], result[2]};
  }

  bool inside_or_periodic(Int3 point) const noexcept {
    const auto periodic = decomposition->periodic();
    const int coordinates[3]{point.x, point.y, point.z};
    const int extents[3]{global_extent.x, global_extent.y, global_extent.z};
    for (std::size_t axis = 0U; axis < 3U; ++axis)
      if ((coordinates[axis] < 0 || coordinates[axis] >= extents[axis]) &&
          !periodic[axis])
        return false;
    return true;
  }

  CanonicalFaceLine face_line(LocalFaceId face) const {
    const auto neighbour = topology->neighbour(face);
    if (!neighbour.has_value())
      throw Error("immersed reconstruction face line has no neighbour");
    const auto periodic = topology->periodic_pair(face);
    const bool reversed =
        periodic.has_value() && topology->global_face_id(face) > *periodic;
    const LocalCellId p = reversed ? *neighbour : topology->owner(face);
    const LocalCellId n = reversed ? topology->owner(face) : *neighbour;
    const auto p_global = topology->global_cell(p);
    const auto n_global = topology->global_cell(n);
    Int3 delta{n_global.x - p_global.x, n_global.y - p_global.y,
               n_global.z - p_global.z};
    if (periodic.has_value()) {
      if (std::abs(delta.x) > 1)
        delta.x = delta.x > 0 ? -1 : 1;
      if (std::abs(delta.y) > 1)
        delta.y = delta.y > 0 ? -1 : 1;
      if (std::abs(delta.z) > 1)
        delta.z = delta.z > 0 ? -1 : 1;
    }
    return {p,
            n,
            {p_global.x - delta.x, p_global.y - delta.y, p_global.z - delta.z},
            {n_global.x + delta.x, n_global.y + delta.y, n_global.z + delta.z},
            reversed};
  }

  double read_global(const runtime::FieldView<const double> &field, Int3 point,
                     int component) const {
    point = wrap(point);
    const auto id = global_id(point);
    if (id >= global_active.size() ||
        global_active[static_cast<std::size_t>(id)] == 0U) {
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
      ++inactive_attempts;
      if (reject_inactive_reads)
        throw Error("immersed reconstruction attempted an inactive cell read");
#endif
      throw Error("immersed reconstruction stencil contains an inactive cell");
    }
    const auto index = field_index(point);
    const double value = field(index.x, index.y, index.z, component);
    if (!std::isfinite(value))
      throw Error("immersed reconstruction active donor is non-finite");
    return value;
  }

  double read_local(const runtime::FieldView<const double> &field,
                    LocalCellId cell, int component) const {
    if (cell >= local_active.size() || local_active[cell] == 0U) {
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
      ++inactive_attempts;
      if (reject_inactive_reads)
        throw Error("immersed reconstruction attempted an inactive cell read");
#endif
      throw Error("immersed reconstruction attempted an inactive cell read");
    }
    const auto point = local_index(cell);
    const double value = field(point.x, point.y, point.z, component);
    if (!std::isfinite(value))
      throw Error("immersed reconstruction active cell is non-finite");
    return value;
  }

  double ghost_value(FiniteVolumeQuantity quantity, std::uint64_t link,
                     const runtime::FieldView<const double> &field,
                     int component) const {
    const std::vector<immersed::WeightedDonor> *donors = nullptr;
    switch (quantity.kind) {
    case FiniteVolumeQuantityKind::density:
      donors = &ghost_plan->density_extrapolation(link).donors;
      break;
    case FiniteVolumeQuantityKind::enthalpy:
    case FiniteVolumeQuantityKind::scalar:
    case FiniteVolumeQuantityKind::pressure:
      donors = &ghost_plan->zero_normal_constraint(link).donors;
      break;
    case FiniteVolumeQuantityKind::velocity:
      donors =
          &ghost_plan
               ->velocity_constraint(link, static_cast<std::size_t>(component))
               .donors;
      break;
    }
    if (donors == nullptr)
      throw Error("immersed reconstruction Ghost quantity is invalid");
    double result = 0.0;
    for (const auto &donor : *donors)
      result += donor.weight *
                read_global(field, logical(donor.global_cell), component);
    if (!std::isfinite(result))
      throw Error("immersed reconstruction Ghost symbol is non-finite");
    return result;
  }

  double stencil_value(FiniteVolumeQuantity quantity, std::uint64_t link,
                       const runtime::FieldView<const double> &field,
                       Int3 point, int component, double ghost_symbol) const {
    point = wrap(point);
    const auto id = global_id(point);
    if (id < global_active.size() &&
        global_active[static_cast<std::size_t>(id)] != 0U)
      return read_global(field, point, component);
    if (link == kInvalidLink)
      throw Error("immersed reconstruction inactive stencil has no link");
    if (!std::isfinite(ghost_symbol))
      throw Error("immersed reconstruction Ghost symbol was not updated");
    static_cast<void>(quantity);
    return ghost_symbol;
  }

  void update_face_ghost_symbols(FiniteVolumeQuantity quantity,
                                 const runtime::FieldView<const double> &field,
                                 std::uint32_t components) const {
    const double missing = std::numeric_limits<double>::quiet_NaN();
    std::fill(face_p_ghost_symbols.begin(), face_p_ghost_symbols.end(),
              missing);
    std::fill(face_n_ghost_symbols.begin(), face_n_ghost_symbols.end(),
              missing);
    for (const auto face : interface_faces) {
      for (std::uint32_t component = 0U; component < components; ++component) {
        const std::size_t index =
            face * 3U + static_cast<std::size_t>(component);
        if (face_p_ghost_link[face] != kInvalidLink)
          face_p_ghost_symbols[index] =
              ghost_value(quantity, face_p_ghost_link[face], field,
                          static_cast<int>(component));
        if (face_n_ghost_link[face] != kInvalidLink)
          face_n_ghost_symbols[index] =
              ghost_value(quantity, face_n_ghost_link[face], field,
                          static_cast<int>(component));
      }
    }
  }

  Real3 ordinary_gradient(GradientScheme scheme, FiniteVolumeQuantity quantity,
                          const runtime::FieldView<const double> &field,
                          LocalCellId cell, int component,
                          bool homogeneous_pressure_dirichlet = false) const {
    const double p = read_local(field, cell, component);
    Real3 result{};
    std::array<double, 6> matrix{};
    Real3 rhs{};
    bool constant = true;
    for (const auto ref : cell_faces[cell]) {
      if (ref.face == kInvalidFace)
        throw Error("immersed reconstruction cell face plan is incomplete");
      const auto owner = topology->owner(ref.face);
      const auto neighbour = topology->neighbour(ref.face);
      const auto patch = topology->patch_id(ref.face);
      Real3 area = geometry->face_area_vector_m2(
          ref.face,
          ref.as_owner ? mesh::FaceSide::owner : mesh::FaceSide::neighbour);
      double q = 0.0;
      double face_value = 0.0;
      Real3 displacement{};
      if (patch.has_value() && !topology->periodic_pair(ref.face).has_value()) {
        if (!ref.as_owner)
          throw Error("physical boundary is attached as neighbour");
        if (quantity.kind == FiniteVolumeQuantityKind::velocity) {
          const Real3 velocity{read_local(field, cell, 0),
                               read_local(field, cell, 1),
                               read_local(field, cell, 2)};
          const auto values =
              boundaries->evaluate_velocity(*patch, velocity, area);
          const double faces[3]{values.face.x, values.face.y, values.face.z};
          const double exteriors[3]{values.exterior.x, values.exterior.y,
                                    values.exterior.z};
          face_value = faces[component];
          q = exteriors[component];
        } else {
          const bool homogeneous_prescribed_pressure =
              homogeneous_pressure_dirichlet &&
              quantity.kind == FiniteVolumeQuantityKind::pressure &&
              boundaries->patch(*patch).pressure_rule() ==
                  boundary::PressureRule::prescribed_value;
          if (homogeneous_prescribed_pressure) {
            face_value = 0.0;
            q = -p;
          } else {
            face_value =
                boundary_scalar(*boundaries, quantity, *patch, p, false);
            q = boundary_scalar(*boundaries, quantity, *patch, p, true);
          }
        }
        displacement = multiply(2.0, subtract(geometry->face_center_m(ref.face),
                                              geometry->cell_center_m(cell)));
      } else {
        if (!neighbour.has_value())
          throw Error("internal reconstruction face has no neighbour");
        const LocalCellId other = ref.as_owner ? *neighbour : owner;
        q = read_local(field, other, component);
        const auto p_center = geometry->cell_center_m(cell);
        displacement = geometry->face_displacement_m(ref.face);
        if (!ref.as_owner)
          displacement = multiply(-1.0, displacement);
        const auto owner_to_face =
            subtract(geometry->face_center_m(ref.face), p_center);
        const double a = std::sqrt(dot(owner_to_face, owner_to_face));
        const auto face_to_other = subtract(displacement, owner_to_face);
        const double b = std::sqrt(dot(face_to_other, face_to_other));
        if (!(a + b > 0.0) || !std::isfinite(a + b))
          throw Error("immersed reconstruction face distance is invalid");
        face_value = (b / (a + b)) * p + (a / (a + b)) * q;
      }
      if (!std::isfinite(q) || !std::isfinite(face_value) ||
          !finite(displacement))
        throw Error("immersed reconstruction gradient stencil is non-finite");
      if (scheme == GradientScheme::green_gauss) {
        result = add(result, multiply(face_value, area));
        constant = constant && face_value == p;
      } else {
        const double d2 = dot(displacement, displacement);
        if (!(d2 > 0.0) || !std::isfinite(d2))
          throw Error("immersed reconstruction least-squares distance invalid");
        const double w = 1.0 / d2;
        matrix[0] += w * displacement.x * displacement.x;
        matrix[1] += w * displacement.x * displacement.y;
        matrix[2] += w * displacement.x * displacement.z;
        matrix[3] += w * displacement.y * displacement.y;
        matrix[4] += w * displacement.y * displacement.z;
        matrix[5] += w * displacement.z * displacement.z;
        rhs = add(rhs, multiply(w * (q - p), displacement));
      }
    }
    if (scheme == GradientScheme::green_gauss) {
      if (constant)
        return {};
      return multiply(1.0 / geometry->cell_volume_m3(cell), result);
    }
    if (scheme != GradientScheme::weighted_least_squares)
      throw Error("invalid immersed reconstruction gradient scheme");
    const double scale_x = std::sqrt(matrix[0]);
    const double scale_y = std::sqrt(matrix[3]);
    const double scale_z = std::sqrt(matrix[5]);
    if (!(scale_x > 0.0) || !(scale_y > 0.0) || !(scale_z > 0.0) ||
        !std::isfinite(scale_x) || !std::isfinite(scale_y) ||
        !std::isfinite(scale_z))
      throw Error("immersed reconstruction least-squares system is singular");
    const double c10 = matrix[1] / (scale_x * scale_y);
    const double c20 = matrix[2] / (scale_x * scale_z);
    const double c21 = matrix[4] / (scale_y * scale_z);
    const double p11 = 1.0 - c10 * c10;
    constexpr double pivot_floor =
        128.0 * std::numeric_limits<double>::epsilon();
    if (!(p11 > pivot_floor) || !std::isfinite(p11))
      throw Error("immersed reconstruction least-squares system is singular");
    const double l11 = std::sqrt(p11);
    const double l21 = (c21 - c20 * c10) / l11;
    const double p22 = 1.0 - c20 * c20 - l21 * l21;
    if (!(p22 > pivot_floor) || !std::isfinite(p22))
      throw Error("immersed reconstruction least-squares system is singular");
    const double l22 = std::sqrt(p22);
    const double scaled_rhs_x = rhs.x / scale_x;
    const double scaled_rhs_y = rhs.y / scale_y;
    const double scaled_rhs_z = rhs.z / scale_z;
    const double y0 = scaled_rhs_x;
    const double y1 = (scaled_rhs_y - c10 * y0) / l11;
    const double y2 = (scaled_rhs_z - c20 * y0 - l21 * y1) / l22;
    const double scaled_result_z = y2 / l22;
    const double scaled_result_y = (y1 - l21 * scaled_result_z) / l11;
    const double scaled_result_x =
        y0 - c10 * scaled_result_y - c20 * scaled_result_z;
    result.x = scaled_result_x / scale_x;
    result.y = scaled_result_y / scale_y;
    result.z = scaled_result_z / scale_z;
    if (!finite(result))
      throw Error("immersed reconstruction gradient is non-finite");
    return result;
  }

  Real3 gradient(GradientScheme scheme, FiniteVolumeQuantity quantity,
                 const runtime::FieldView<const double> &field,
                 LocalCellId cell, int component) const {
    if (quantity.kind == FiniteVolumeQuantityKind::velocity &&
        cell < owned_wall_links.size() && !owned_wall_links[cell].empty()) {
      Real3 result{};
      for (const auto wall_link : owned_wall_links[cell])
        result = add(result,
                     immersed::detail::gradient_with_origin_constraint(
                         reconstructions[static_cast<std::size_t>(wall_link)],
                         geometry->cell_center_m(cell), field,
                         static_cast<std::size_t>(component), 0.0));
      return multiply(1.0 / static_cast<double>(owned_wall_links[cell].size()),
                      result);
    }
    const auto link = owned_row_link[cell];
    if (link != kInvalidLink) {
      const auto &reconstruction =
          reconstructions[static_cast<std::size_t>(link)];
      if (quantity.kind == FiniteVolumeQuantityKind::velocity)
        return immersed::detail::gradient_with_origin_constraint(
            reconstruction, geometry->cell_center_m(cell), field,
            static_cast<std::size_t>(component), 0.0);
      return reconstruction.gradient(geometry->cell_center_m(cell), field,
                                     static_cast<std::size_t>(component));
    }
    return ordinary_gradient(scheme, quantity, field, cell, component);
  }

  Real3 pressure_gradient_with_wall_normal_constraint(
      GradientScheme scheme, const runtime::FieldView<const double> &field,
      LocalCellId cell,
      const std::vector<double> &normal_gradient_by_link) const {
    if (cell >= owned_wall_links.size())
      throw Error("immersed pressure reconstruction row is outside layout");
    const auto &wall_links = owned_wall_links[cell];
    if (wall_links.empty()) {
      const auto link = owned_row_link[cell];
      if (link != kInvalidLink) {
        if (link >= normal_gradient_by_link.size() ||
            !std::isfinite(normal_gradient_by_link[link]))
          throw Error(
              "immersed pressure reconstruction wall condition is missing");
        return immersed::detail::gradient_with_origin_normal_gradient(
            reconstructions[static_cast<std::size_t>(link)],
            geometry->cell_center_m(cell), field, 0U,
            normal_gradient_by_link[link]);
      }
      return ordinary_gradient(scheme, FiniteVolumeQuantity::pressure(), field,
                               cell, 0, true);
    }
    Real3 result{};
    for (const auto link : wall_links) {
      if (link >= normal_gradient_by_link.size() ||
          !std::isfinite(normal_gradient_by_link[link]))
        throw Error(
            "immersed pressure reconstruction wall condition is missing");
      result =
          add(result, immersed::detail::gradient_with_origin_normal_gradient(
                          reconstructions[static_cast<std::size_t>(link)],
                          geometry->cell_center_m(cell), field, 0U,
                          normal_gradient_by_link[link]));
    }
    return multiply(1.0 / static_cast<double>(wall_links.size()), result);
  }

  double transport_face(FiniteVolumeQuantity quantity,
                        const FaceMassFlux &mass_flux,
                        const runtime::FieldView<const double> &field,
                        LocalFaceId face) const {
    const auto owner = topology->owner(face);
    const auto neighbour = topology->neighbour(face);
    const auto patch = topology->patch_id(face);
    if (patch.has_value() && !topology->periodic_pair(face).has_value()) {
      if (local_active[owner] == 0U)
        return 0.0;
      const double p = read_local(field, owner, 0);
      return boundary_scalar(*boundaries, quantity, *patch, p, false);
    }
    if (!neighbour.has_value())
      throw Error("transport reconstruction face has no neighbour");
    const bool owner_active = local_active[owner] != 0U;
    const bool neighbour_active = local_active[*neighbour] != 0U;
    if (!owner_active || !neighbour_active)
      return 0.0;
    const auto line = face_line(face);
    const double qp = read_local(field, line.p, 0);
    const double qn = read_local(field, line.n, 0);
    const double qm =
        inside_or_periodic(line.pm1)
            ? stencil_value(quantity, face_p_ghost_link[face], field, line.pm1,
                            0, face_p_ghost_symbols[face * 3U])
            : qp;
    const double qplus =
        inside_or_periodic(line.np1)
            ? stencil_value(quantity, face_n_ghost_link[face], field, line.np1,
                            0, face_n_ghost_symbols[face * 3U])
            : qn;
    const double slope_p = monotonized_central(qp - qm, qn - qp);
    const double slope_n = monotonized_central(qn - qp, qplus - qn);
    const double flux = mass_flux.value_for_immersed_reconstruction(face);
    const double canonical_flux = line.reversed ? -flux : flux;
    return canonical_flux >= 0.0 ? qp + 0.5 * slope_p : qn - 0.5 * slope_n;
  }

  double momentum_face(const FaceMassFlux &mass_flux,
                       const runtime::FieldView<const double> &field,
                       LocalFaceId face, int component,
                       bool interface_fallback) const {
    const auto owner = topology->owner(face);
    const auto neighbour = topology->neighbour(face);
    const auto patch = topology->patch_id(face);
    if (patch.has_value() && !topology->periodic_pair(face).has_value()) {
      if (local_active[owner] == 0U)
        return 0.0;
      const Real3 p{read_local(field, owner, 0), read_local(field, owner, 1),
                    read_local(field, owner, 2)};
      const auto values = boundaries->evaluate_velocity(
          *patch, p,
          geometry->face_area_vector_m2(face, mesh::FaceSide::owner));
      return component == 0   ? values.face.x
             : component == 1 ? values.face.y
                              : values.face.z;
    }
    if (!neighbour.has_value())
      throw Error("momentum reconstruction face has no neighbour");
    if (local_active[owner] == 0U || local_active[*neighbour] == 0U)
      return 0.0;
    bool use_centered = !interface_fallback &&
                        geometry->face_skewness(face) <= 0.25 &&
                        geometry->face_non_orthogonality_degrees(face) <= 70.0;
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
    use_centered = use_centered && !force_momentum_fallback_for_test;
#endif
    if (use_centered)
      return 0.5 * read_local(field, owner, component) +
             0.5 * read_local(field, *neighbour, component);

    const auto line = face_line(face);
    const double qp = read_local(field, line.p, component);
    const double qn = read_local(field, line.n, component);
    const double qm =
        inside_or_periodic(line.pm1)
            ? stencil_value(
                  FiniteVolumeQuantity::velocity(), face_p_ghost_link[face],
                  field, line.pm1, component,
                  face_p_ghost_symbols[face * 3U +
                                       static_cast<std::size_t>(component)])
            : qp;
    const double qplus =
        inside_or_periodic(line.np1)
            ? stencil_value(
                  FiniteVolumeQuantity::velocity(), face_n_ghost_link[face],
                  field, line.np1, component,
                  face_n_ghost_symbols[face * 3U +
                                       static_cast<std::size_t>(component)])
            : qn;
    const double slope_p = monotonized_central(qp - qm, qn - qp);
    const double slope_n = monotonized_central(qn - qp, qplus - qn);
    const double flux = mass_flux.value_for_immersed_reconstruction(face);
    return (line.reversed ? -flux : flux) >= 0.0 ? qp + 0.5 * slope_p
                                                 : qn - 0.5 * slope_n;
  }
};

namespace {

void assign_link(std::vector<std::uint64_t> &links, LocalCellId cell,
                 std::uint64_t link) {
  if (cell >= links.size())
    throw Error("immersed reconstruction link cell is outside layout");
  links[cell] = std::min(links[cell], link);
}

std::uint64_t link_for_pair(const std::vector<GlobalLinkPair> &links,
                            mesh::GlobalCellId first,
                            mesh::GlobalCellId second) {
  const auto key = std::pair{std::min(first, second), std::max(first, second)};
  const auto found = std::lower_bound(
      links.begin(), links.end(), key,
      [](const GlobalLinkPair &candidate, const auto &target) {
        return std::pair{candidate.first, candidate.second} < target;
      });
  if (found != links.end() && std::pair{found->first, found->second} == key)
    return found->id;
  return kInvalidLink;
}

template <class ImplType, class View>
void validate_common_cell_view(const ImplType &impl, const View &view,
                               std::uint32_t components, bool require_halo) {
  if (!same(view.interior_extent(), impl.local_extent) ||
      view.components() != components ||
      (require_halo &&
       view.ghost_width() < static_cast<int>(impl.halo_reach))) {
    throw Error("immersed reconstruction cell field layout is invalid");
  }
}

template <class ImplType>
void validate_bindings_collective(const ImplType &impl,
                                  const ReconstructionFieldBinding &input,
                                  const ReconstructionFieldBinding &output,
                                  bool face_output,
                                  std::uint32_t input_components,
                                  std::uint32_t output_components) {
  bool local_ok = true;
  std::string message;
  try {
    if (&input.storage == &output.storage && input.field == output.field)
      throw Error("immersed reconstruction input and output fields alias");
    const auto in = input.storage.acquire_read<double>(
        input.access_plan, input.phase, input.actor, input.field);
    validate_common_cell_view(impl, in, input_components, true);
    if (face_output) {
      const auto out = output.storage.acquire_face_write<double>(
          output.access_plan, output.phase, output.actor, output.field);
      if (out.face_count() != impl.topology->local_face_count() ||
          out.components() != output_components)
        throw Error("immersed reconstruction face field layout is invalid");
    } else {
      const auto out = output.storage.acquire_write<double>(
          output.access_plan, output.phase, output.actor, output.field);
      validate_common_cell_view(impl, out, output_components, false);
    }
  } catch (...) {
    local_ok = false;
    message = exception_message();
  }
  const auto status = runtime::collective_status(*impl.mpi, local_ok, message);
  if (!status.ok)
    throw_collective(status);
}

template <class ImplType, class Validation>
void validate_operation_collective(const ImplType &impl,
                                   Validation &&validation) {
  bool local_ok = true;
  std::string message;
  try {
    validation();
  } catch (...) {
    local_ok = false;
    message = exception_message();
  }
  const auto status = runtime::collective_status(*impl.mpi, local_ok, message);
  if (!status.ok)
    throw_collective(status);
}

template <class ImplType, class Interior, class RemoteSymbols, class Partition,
          class Interface>
void execute_order(ImplType &impl, const ReconstructionFieldBinding &input,
                   const char *operation, Interior &&interior,
                   RemoteSymbols &&remote_symbols, Partition &&partition,
                   Interface &&interface) {
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
  trace_stage(0U);
#endif
  impl.halo->begin(input.storage, input.field);
  bool local_ok = true;
  std::string message;
  try {
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
    trace_stage(1U);
#endif
    interior();
  } catch (...) {
    local_ok = false;
    message = exception_message();
  }
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
  trace_stage(2U);
#endif
  try {
    impl.halo->wait(input.storage, input.field);
  } catch (...) {
    local_ok = false;
    if (message.empty())
      message = exception_message();
  }
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
  trace_stage(3U);
#endif
  if (local_ok) {
    try {
      remote_symbols();
    } catch (...) {
      local_ok = false;
      message = exception_message();
    }
  }
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
  trace_stage(4U);
#endif
  if (local_ok) {
    try {
      partition();
    } catch (...) {
      local_ok = false;
      message = exception_message();
    }
  }
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
  trace_stage(5U);
#endif
  if (local_ok) {
    try {
      interface();
    } catch (...) {
      local_ok = false;
      message = exception_message();
    }
  }
  if (!local_ok)
    message = std::string(operation) + ": " + message;
  const auto status = runtime::collective_status(*impl.mpi, local_ok, message);
  if (!status.ok)
    throw_collective(status);
}

} // namespace

ImmersedReconstruction ImmersedReconstruction::create(
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const boundary::BoundaryRegistry &boundaries,
    const immersed::ImmersedDomain &domain,
    const immersed::GhostStencilPlan &ghost_plan,
    const runtime::StructuredDecomposition &decomposition,
    const runtime::MpiContext &mpi, runtime::HaloExchange &halo) {
  bool local_ok = true;
  std::string message;
  try {
    if (mpi.comm() == MPI_COMM_NULL || decomposition.comm() == MPI_COMM_NULL)
      throw Error("immersed reconstruction MPI dependency is invalid");
    int communicator_relation = MPI_UNEQUAL;
    check_mpi(MPI_Comm_compare(mpi.comm(), decomposition.comm(),
                               &communicator_relation),
              "MPI_Comm_compare immersed reconstruction dependencies");
    if (communicator_relation != MPI_IDENT &&
        communicator_relation != MPI_CONGRUENT)
      throw Error("immersed reconstruction communicators are incompatible");
    if (!geometry.compatible(topology) ||
        !same(topology.global_extent(), decomposition.global_extent()) ||
        !same(topology.owned_global_box(), decomposition.owned_box()) ||
        !same(topology.owned_global_box(), geometry.owned_global_box()) ||
        !halo.is_compatible_with(decomposition))
      throw Error("immersed reconstruction dependencies are incompatible");
    if (ghost_plan.maximum_halo_reach() == 0U ||
        ghost_plan.maximum_halo_reach() >
            static_cast<std::uint32_t>(halo.ghost_width()))
      throw Error("immersed reconstruction Halo width is insufficient");
  } catch (...) {
    local_ok = false;
    message = exception_message();
  }
  const auto compatibility = runtime::collective_status(mpi, local_ok, message);
  if (!compatibility.ok)
    throw_collective(compatibility);

  auto impl = std::make_unique<Impl>();
  impl->topology = &topology;
  impl->geometry = &geometry;
  impl->boundaries = &boundaries;
  impl->domain = &domain;
  impl->ghost_plan = &ghost_plan;
  impl->decomposition = &decomposition;
  impl->mpi = &mpi;
  impl->halo = &halo;
  impl->owned_box = decomposition.owned_box();
  impl->local_extent = decomposition.local_extent();
  impl->global_extent = decomposition.global_extent();
  impl->halo_reach = ghost_plan.maximum_halo_reach();
  impl->links = gather_links(domain, mpi);
  impl->link_pairs.reserve(impl->links.size());
  for (const auto &link : impl->links)
    impl->link_pairs.push_back({std::min(link.fluid, link.solid),
                                std::max(link.fluid, link.solid), link.id});
  std::sort(impl->link_pairs.begin(), impl->link_pairs.end(),
            [](const auto &left, const auto &right) {
              return std::tuple{left.first, left.second, left.id} <
                     std::tuple{right.first, right.second, right.id};
            });
  if (std::adjacent_find(impl->link_pairs.begin(), impl->link_pairs.end(),
                         [](const auto &left, const auto &right) {
                           return left.first == right.first &&
                                  left.second == right.second;
                         }) != impl->link_pairs.end())
    throw Error("immersed reconstruction cell pair has duplicate links");
  impl->global_active = gather_global_active(domain, topology, mpi);
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
  std::vector<std::uint64_t> reconstruction_snapshot;
  const void *const reconstruction_identity =
      static_cast<const void *>(impl.get());
  ReconstructionSnapshotRollback reconstruction_snapshot_rollback(
      reconstruction_identity);
#endif
  local_ok = true;
  message.clear();
  try {
    impl->reconstructions.reserve(impl->links.size());
    for (const auto &link : impl->links)
      impl->reconstructions.push_back(ghost_plan.reconstruction(link.id));
    impl->wall_normal_gradient_scratch.resize(
        impl->reconstructions.size(), std::numeric_limits<double>::quiet_NaN());

    impl->local_active.resize(topology.local_cell_count(), 0U);
    impl->local_interface.resize(topology.local_cell_count(), 0U);
    impl->local_cell_link.resize(topology.local_cell_count(), kInvalidLink);
    impl->owned_row_link.resize(topology.owned_cell_count(), kInvalidLink);
    impl->owned_wall_links.resize(topology.owned_cell_count());
    for (LocalCellId cell = 0U; cell < topology.local_cell_count(); ++cell)
      impl->local_active[cell] =
          domain.region(cell) == CellRegion::fluid ? 1U : 0U;

    for (const auto &link : impl->links) {
      const auto &donors = ghost_plan.density_extrapolation(link.id).donors;
      const auto row_can_access_reconstruction = [&](mesh::GlobalCellId row) {
        const auto row_logical = impl->logical(row);
        const auto periodic = decomposition.periodic();
        const int extents[3]{impl->global_extent.x, impl->global_extent.y,
                             impl->global_extent.z};
        const int row_coordinates[3]{row_logical.x, row_logical.y,
                                     row_logical.z};
        for (const auto &donor : donors) {
          const auto donor_logical = impl->logical(donor.global_cell);
          const int donor_coordinates[3]{donor_logical.x, donor_logical.y,
                                         donor_logical.z};
          for (std::size_t axis = 0U; axis < 3U; ++axis) {
            int distance =
                std::abs(donor_coordinates[axis] - row_coordinates[axis]);
            if (periodic[axis])
              distance = std::min(distance, extents[axis] - distance);
            if (distance > static_cast<int>(impl->halo_reach))
              return false;
          }
        }
        return true;
      };
      if (const auto local = topology.find_local_cell(link.fluid);
          local.has_value() && impl->local_active[*local] != 0U) {
        impl->local_interface[*local] = 1U;
        if (!row_can_access_reconstruction(link.fluid))
          throw Error(
              "immersed reconstruction wall row exceeds its Halo reach");
        if (*local < topology.owned_cell_count())
          impl->owned_wall_links[*local].push_back(link.id);
        assign_link(impl->local_cell_link, *local, link.id);
      }
      for (const auto &donor : donors) {
        if (donor.global_cell >= impl->global_active.size() ||
            impl->global_active[static_cast<std::size_t>(donor.global_cell)] ==
                0U)
          throw Error("immersed reconstruction donor is not active");
        if (const auto local = topology.find_local_cell(donor.global_cell);
            local.has_value() && impl->local_active[*local] != 0U) {
          impl->local_interface[*local] = 1U;
          if (row_can_access_reconstruction(donor.global_cell))
            assign_link(impl->local_cell_link, *local, link.id);
        }
      }
    }
    for (LocalCellId cell = 0U; cell < topology.owned_cell_count(); ++cell) {
      impl->owned_row_link[cell] = impl->local_cell_link[cell];
    }
    for (auto &links : impl->owned_wall_links) {
      std::sort(links.begin(), links.end());
      if (std::adjacent_find(links.begin(), links.end()) != links.end())
        throw Error("immersed pressure reconstruction wall link is duplicated");
    }

    impl->cell_faces.resize(topology.owned_cell_count());
    for (auto &faces : impl->cell_faces)
      for (auto &face : faces)
        face = CellFaceRef{};
    std::vector<std::size_t> face_counts(topology.owned_cell_count(), 0U);
    std::vector<std::uint8_t> partition(topology.owned_cell_count(), 0U);
    for (LocalFaceId face = 0U; face < topology.local_face_count(); ++face) {
      const auto owner = topology.owner(face);
      const auto neighbour = topology.neighbour(face);
      const auto add_ref = [&](LocalCellId cell, bool as_owner) {
        if (cell >= topology.owned_cell_count())
          return;
        if (face_counts[cell] >= 6U)
          throw Error("immersed reconstruction cell has too many faces");
        impl->cell_faces[cell][face_counts[cell]++] = {face, as_owner};
      };
      add_ref(owner, true);
      const bool patch_face = topology.patch_id(face).has_value();
      if (neighbour.has_value() && !patch_face)
        add_ref(*neighbour, false);
      if (owner < topology.owned_cell_count() && neighbour.has_value() &&
          topology.cell_ownership(*neighbour) == EntityOwnership::ghost &&
          impl->local_active[owner] != 0U &&
          impl->local_active[*neighbour] != 0U)
        partition[owner] = 1U;
      if (!patch_face && neighbour.has_value() &&
          *neighbour < topology.owned_cell_count() &&
          topology.cell_ownership(owner) == EntityOwnership::ghost &&
          impl->local_active[owner] != 0U &&
          impl->local_active[*neighbour] != 0U)
        partition[*neighbour] = 1U;
    }
    for (LocalCellId cell = 0U; cell < topology.owned_cell_count(); ++cell) {
      if (impl->local_active[cell] == 0U)
        continue;
      if (face_counts[cell] != 6U)
        throw Error("immersed reconstruction active row has incomplete faces");
      if (impl->local_interface[cell] != 0U)
        impl->interface_rows.push_back(cell);
      else if (partition[cell] != 0U)
        impl->partition_rows.push_back(cell);
      else
        impl->interior_rows.push_back(cell);
    }

    impl->face_p_ghost_link.resize(topology.local_face_count(), kInvalidLink);
    impl->face_n_ghost_link.resize(topology.local_face_count(), kInvalidLink);
    const auto stencil_link = [&](LocalCellId fluid, Int3 point) {
      if (!impl->inside_or_periodic(point))
        return kInvalidLink;
      point = impl->wrap(point);
      const auto point_id = topology.global_cell_id(point);
      if (impl->global_active[static_cast<std::size_t>(point_id)] != 0U)
        return kInvalidLink;
      const auto link = link_for_pair(impl->link_pairs,
                                      topology.global_cell_id(fluid), point_id);
      if (link == kInvalidLink)
        throw Error("immersed reconstruction Ghost stencil link is missing");
      return link;
    };
    const auto active_stencil_is_remote = [&](Int3 point) {
      if (!impl->inside_or_periodic(point))
        return false;
      point = impl->wrap(point);
      const auto point_id = topology.global_cell_id(point);
      if (impl->global_active[static_cast<std::size_t>(point_id)] == 0U)
        return false;
      return point.x < impl->owned_box.begin.x ||
             point.x >= impl->owned_box.end.x ||
             point.y < impl->owned_box.begin.y ||
             point.y >= impl->owned_box.end.y ||
             point.z < impl->owned_box.begin.z ||
             point.z >= impl->owned_box.end.z;
    };
    for (LocalFaceId face = 0U; face < topology.local_face_count(); ++face) {
      const auto owner = topology.owner(face);
      const auto neighbour = topology.neighbour(face);
      bool extended_stencil_is_remote = false;
      if (neighbour.has_value() && impl->local_active[owner] != 0U &&
          impl->local_active[*neighbour] != 0U) {
        const auto line = impl->face_line(face);
        impl->face_p_ghost_link[face] = stencil_link(line.p, line.pm1);
        impl->face_n_ghost_link[face] = stencil_link(line.n, line.np1);
        extended_stencil_is_remote = active_stencil_is_remote(line.pm1) ||
                                     active_stencil_is_remote(line.np1);
      }
      const bool interface =
          impl->local_interface[owner] != 0U ||
          (neighbour.has_value() && impl->local_interface[*neighbour] != 0U) ||
          (neighbour.has_value() &&
           impl->local_active[owner] != impl->local_active[*neighbour]);
      if (interface) {
        impl->interface_faces.push_back(face);
      } else if ((topology.cell_ownership(owner) == EntityOwnership::ghost) ||
                 (neighbour.has_value() &&
                  topology.cell_ownership(*neighbour) ==
                      EntityOwnership::ghost) ||
                 extended_stencil_is_remote) {
        impl->partition_faces.push_back(face);
      } else {
        impl->interior_faces.push_back(face);
      }
    }

    impl->cell_scratch.resize(topology.owned_cell_count() * 9U, 0.0);
    impl->face_scratch.resize(topology.local_face_count() * 3U, 0.0);
    impl->face_p_ghost_symbols.resize(topology.local_face_count() * 3U, 0.0);
    impl->face_n_ghost_symbols.resize(topology.local_face_count() * 3U, 0.0);
    std::uint64_t fingerprint = kFnvOffset;
    hash_u64(fingerprint, UINT64_C(0x48554e4449524333));
    hash_u64(fingerprint, domain.classification_fingerprint());
    hash_u64(fingerprint, domain.active_boundaries().fingerprint());
    hash_u64(fingerprint, ghost_plan.fingerprint());
    hash_u64(fingerprint, static_cast<std::uint64_t>(impl->halo_reach));
    const auto grid = decomposition.process_grid();
    hash_u64(fingerprint, static_cast<std::uint64_t>(grid.x));
    hash_u64(fingerprint, static_cast<std::uint64_t>(grid.y));
    hash_u64(fingerprint, static_cast<std::uint64_t>(grid.z));
    for (const auto &link : impl->links) {
      hash_u64(fingerprint, link.id);
      hash_u64(fingerprint, link.fluid);
      hash_u64(fingerprint, link.solid);
    }
    impl->fingerprint = fingerprint;
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
    reconstruction_snapshot.reserve(6U + impl->interior_rows.size() +
                                    impl->partition_rows.size() +
                                    impl->interface_rows.size() +
                                    impl->interface_faces.size());
    reconstruction_snapshot.push_back(
        static_cast<std::uint64_t>(impl->interior_rows.size()));
    reconstruction_snapshot.push_back(
        static_cast<std::uint64_t>(impl->partition_rows.size()));
    reconstruction_snapshot.push_back(
        static_cast<std::uint64_t>(impl->interface_rows.size()));
    reconstruction_snapshot.push_back(
        static_cast<std::uint64_t>(impl->interface_faces.size()));
    reconstruction_snapshot.push_back(static_cast<std::uint64_t>(
        impl->domain->active_cells().owned_active_count()));
    reconstruction_snapshot.push_back(impl->fingerprint);
    const auto append_cells = [&](const std::vector<LocalCellId> &rows) {
      for (const auto cell : rows)
        reconstruction_snapshot.push_back(
            impl->topology->global_cell_id(cell));
    };
    append_cells(impl->interior_rows);
    append_cells(impl->partition_rows);
    append_cells(impl->interface_rows);
    for (const auto face : impl->interface_faces)
      reconstruction_snapshot.push_back(
          impl->topology->global_face_id(face));
    register_reconstruction_snapshot(reconstruction_identity,
                                     std::move(reconstruction_snapshot));
    reconstruction_snapshot_rollback.arm();
#endif
  } catch (...) {
    local_ok = false;
    message = exception_message();
  }
  const auto plan_status = runtime::collective_status(mpi, local_ok, message);
  if (!plan_status.ok) {
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
    reconstruction_snapshot_rollback.rollback();
#endif
    throw_collective(plan_status);
  }
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
  auto result = ImmersedReconstruction(std::move(impl));
  reconstruction_snapshot_rollback.release();
  return result;
#else
  return ImmersedReconstruction(std::move(impl));
#endif
}

ImmersedReconstruction::ImmersedReconstruction(
    std::unique_ptr<Impl> impl) noexcept
    : impl_(std::move(impl)) {}
ImmersedReconstruction::~ImmersedReconstruction() noexcept {
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
  if (impl_)
    erase_reconstruction_snapshot(static_cast<const void *>(impl_.get()));
#endif
}
ImmersedReconstruction::ImmersedReconstruction(
    ImmersedReconstruction &&) noexcept = default;

ImmersedReconstruction::Impl &ImmersedReconstruction::require_impl() const {
  if (!impl_)
    throw Error("immersed reconstruction has been moved from");
  return *impl_;
}

void ImmersedReconstruction::compute_gradient(
    GradientScheme scheme, FiniteVolumeQuantity quantity,
    const ReconstructionFieldBinding &cell_values,
    const ReconstructionFieldBinding &cell_gradients) const {
  auto &impl = require_impl();
  OperationGuard guard(impl.active);
  if (scoped_wall_normal_gradients != nullptr) {
    auto &normal_gradient_by_link = impl.wall_normal_gradient_scratch;
    std::fill(normal_gradient_by_link.begin(), normal_gradient_by_link.end(),
              std::numeric_limits<double>::quiet_NaN());
    validate_operation_collective(impl, [&] {
      if (quantity.kind != FiniteVolumeQuantityKind::pressure)
        throw Error(
            "immersed wall-normal constraints require pressure quantity");
      if (scheme != GradientScheme::green_gauss &&
          scheme != GradientScheme::weighted_least_squares)
        throw Error("invalid immersed reconstruction gradient scheme");
      for (const auto &condition : *scoped_wall_normal_gradients) {
        if (condition.link >= normal_gradient_by_link.size() ||
            !std::isfinite(condition.value))
          throw Error(
              "immersed pressure reconstruction wall condition is invalid");
        if (std::isfinite(normal_gradient_by_link[condition.link]))
          throw Error(
              "immersed pressure reconstruction wall condition is duplicated");
        normal_gradient_by_link[condition.link] = condition.value;
      }
      for (const auto cell : impl.interface_rows)
        for (const auto link : impl.owned_wall_links[cell])
          if (link >= normal_gradient_by_link.size() ||
              !std::isfinite(normal_gradient_by_link[link]))
            throw Error(
                "immersed pressure reconstruction wall condition is missing");
    });
    validate_bindings_collective(impl, cell_values, cell_gradients, false, 1U,
                                 3U);
    const auto input = cell_values.storage.acquire_read<double>(
        cell_values.access_plan, cell_values.phase, cell_values.actor,
        cell_values.field);
    const auto output = cell_gradients.storage.acquire_write<double>(
        cell_gradients.access_plan, cell_gradients.phase, cell_gradients.actor,
        cell_gradients.field);
    std::fill(impl.cell_scratch.begin(), impl.cell_scratch.end(), 0.0);
    const auto ordinary_rows = [&](const std::vector<LocalCellId> &selected) {
      for (const auto cell : selected) {
        const auto value = impl.ordinary_gradient(
            scheme, FiniteVolumeQuantity::pressure(), input, cell, 0, true);
        const std::size_t base = cell * 9U;
        impl.cell_scratch[base] = value.x;
        impl.cell_scratch[base + 1U] = value.y;
        impl.cell_scratch[base + 2U] = value.z;
      }
    };
    const auto interface_rows = [&] {
      for (const auto cell : impl.interface_rows) {
        const auto value = impl.pressure_gradient_with_wall_normal_constraint(
            scheme, input, cell, normal_gradient_by_link);
        const std::size_t base = cell * 9U;
        impl.cell_scratch[base] = value.x;
        impl.cell_scratch[base + 1U] = value.y;
        impl.cell_scratch[base + 2U] = value.z;
      }
    };
    execute_order(
        impl, cell_values, "immersed constrained pressure gradient",
        [&] { ordinary_rows(impl.interior_rows); }, [] {},
        [&] { ordinary_rows(impl.partition_rows); }, interface_rows);
    for (LocalCellId cell = 0U; cell < impl.topology->owned_cell_count();
         ++cell) {
      if (impl.local_active[cell] == 0U)
        continue;
      const auto index = impl.local_index(cell);
      for (std::uint32_t component = 0U; component < 3U; ++component)
        output(index.x, index.y, index.z, static_cast<int>(component)) =
            impl.cell_scratch[cell * 9U + component];
    }
    return;
  }
  validate_operation_collective(impl, [&] {
    validate_quantity(quantity, *impl.boundaries, true);
    if (scheme != GradientScheme::green_gauss &&
        scheme != GradientScheme::weighted_least_squares)
      throw Error("invalid immersed reconstruction gradient scheme");
  });
  const std::uint32_t components =
      quantity.kind == FiniteVolumeQuantityKind::velocity ? 3U : 1U;
  validate_bindings_collective(impl, cell_values, cell_gradients, false,
                               components, components * 3U);
  const auto input = cell_values.storage.acquire_read<double>(
      cell_values.access_plan, cell_values.phase, cell_values.actor,
      cell_values.field);
  const auto output = cell_gradients.storage.acquire_write<double>(
      cell_gradients.access_plan, cell_gradients.phase, cell_gradients.actor,
      cell_gradients.field);
  std::fill(impl.cell_scratch.begin(), impl.cell_scratch.end(), 0.0);
  const auto rows = [&](const std::vector<LocalCellId> &selected) {
    for (const auto cell : selected)
      for (std::uint32_t component = 0U; component < components; ++component) {
        const auto value = impl.gradient(scheme, quantity, input, cell,
                                         static_cast<int>(component));
        const std::size_t base =
            cell * 9U + static_cast<std::size_t>(component) * 3U;
        impl.cell_scratch[base] = value.x;
        impl.cell_scratch[base + 1U] = value.y;
        impl.cell_scratch[base + 2U] = value.z;
      }
  };
  execute_order(
      impl, cell_values, "immersed reconstruction gradient",
      [&] { rows(impl.interior_rows); }, [] {},
      [&] { rows(impl.partition_rows); }, [&] { rows(impl.interface_rows); });
  for (LocalCellId cell = 0U; cell < impl.topology->owned_cell_count();
       ++cell) {
    if (impl.local_active[cell] == 0U)
      continue;
    const auto index = impl.local_index(cell);
    for (std::uint32_t component = 0U; component < components * 3U; ++component)
      output(index.x, index.y, index.z, static_cast<int>(component)) =
          impl.cell_scratch[cell * 9U + component];
  }
}

void ImmersedReconstruction::reconstruct_transport_faces(
    FiniteVolumeQuantity quantity, const FaceMassFlux &mass_flux,
    const ReconstructionFieldBinding &cell_values,
    const ReconstructionFieldBinding &face_values) const {
  auto &impl = require_impl();
  OperationGuard guard(impl.active);
  validate_operation_collective(impl, [&] {
    validate_quantity(quantity, *impl.boundaries, false);
    mass_flux.require_immersed_reconstruction_compatible(*impl.topology);
  });
  validate_bindings_collective(impl, cell_values, face_values, true, 1U, 1U);
  const auto input = cell_values.storage.acquire_read<double>(
      cell_values.access_plan, cell_values.phase, cell_values.actor,
      cell_values.field);
  const auto output = face_values.storage.acquire_face_write<double>(
      face_values.access_plan, face_values.phase, face_values.actor,
      face_values.field);
  std::fill(impl.face_scratch.begin(), impl.face_scratch.end(), 0.0);
  const auto faces = [&](const std::vector<LocalFaceId> &selected) {
    for (const auto face : selected)
      impl.face_scratch[face * 3U] =
          impl.transport_face(quantity, mass_flux, input, face);
  };
  execute_order(
      impl, cell_values, "immersed transport face reconstruction",
      [&] { faces(impl.interior_faces); },
      [&] { impl.update_face_ghost_symbols(quantity, input, 1U); },
      [&] { faces(impl.partition_faces); },
      [&] { faces(impl.interface_faces); });
  for (LocalFaceId face = 0U; face < impl.topology->local_face_count(); ++face)
    output(face, 0) = impl.face_scratch[face * 3U];
}

void ImmersedReconstruction::reconstruct_momentum_faces(
    const FaceMassFlux &mass_flux, const ReconstructionFieldBinding &velocity,
    const ReconstructionFieldBinding &face_velocity) const {
  auto &impl = require_impl();
  OperationGuard guard(impl.active);
  validate_operation_collective(impl, [&] {
    mass_flux.require_immersed_reconstruction_compatible(*impl.topology);
  });
  validate_bindings_collective(impl, velocity, face_velocity, true, 3U, 3U);
  const auto input = velocity.storage.acquire_read<double>(
      velocity.access_plan, velocity.phase, velocity.actor, velocity.field);
  const auto output = face_velocity.storage.acquire_face_write<double>(
      face_velocity.access_plan, face_velocity.phase, face_velocity.actor,
      face_velocity.field);
  std::fill(impl.face_scratch.begin(), impl.face_scratch.end(), 0.0);
  const auto faces = [&](const std::vector<LocalFaceId> &selected,
                         bool interface_fallback) {
    for (const auto face : selected)
      for (int component = 0; component < 3; ++component)
        impl.face_scratch[face * 3U + static_cast<std::size_t>(component)] =
            impl.momentum_face(mass_flux, input, face, component,
                               interface_fallback);
  };
  execute_order(
      impl, velocity, "immersed momentum face reconstruction",
      [&] { faces(impl.interior_faces, false); },
      [&] {
        impl.update_face_ghost_symbols(FiniteVolumeQuantity::velocity(), input,
                                       3U);
      },
      [&] { faces(impl.partition_faces, false); },
      [&] { faces(impl.interface_faces, true); });
  for (LocalFaceId face = 0U; face < impl.topology->local_face_count(); ++face)
    for (int component = 0; component < 3; ++component)
      output(face, component) =
          impl.face_scratch[face * 3U + static_cast<std::size_t>(component)];
}

std::uint64_t ImmersedReconstruction::dependency_fingerprint() const noexcept {
  return impl_ ? impl_->fingerprint : 0U;
}

void ImmersedReconstruction::require_immersed_operator_compatible(
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const immersed::ImmersedDomain &domain,
    const immersed::GhostStencilPlan &ghost_plan) const {
  const auto &impl = require_impl();
  if (impl.topology != &topology || impl.geometry != &geometry ||
      impl.domain != &domain || impl.ghost_plan != &ghost_plan)
    throw Error(
        "immersed operator reconstruction dependencies are incompatible");
}

const runtime::MpiContext &
ImmersedReconstruction::mpi_for_immersed_operator() const {
  return *require_impl().mpi;
}

namespace detail {

void compute_pressure_gradient_with_wall_normal_constraints(
    const ImmersedReconstruction &provider, GradientScheme scheme,
    const ReconstructionFieldBinding &cell_values,
    const ReconstructionFieldBinding &cell_gradients,
    const std::vector<ImmersedWallNormalGradient> &conditions) {
  ScopedWallNormalGradients scope(conditions);
  provider.compute_gradient(scheme, FiniteVolumeQuantity::pressure(),
                            cell_values, cell_gradients);
}

} // namespace detail

#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
namespace detail {

std::vector<std::uint64_t>
immersed_reconstruction_snapshot_raw(const void *identity) {
  return lookup_reconstruction_snapshot(identity);
}

void immersed_reconstruction_reset_trace_raw() {
  execution_trace.clear();
  inactive_attempts = 0U;
}

std::vector<std::uint8_t> immersed_reconstruction_trace_raw() {
  return execution_trace;
}

std::size_t immersed_reconstruction_inactive_read_attempts_raw() noexcept {
  return inactive_attempts;
}

void immersed_reconstruction_fail_on_inactive_read_raw(
    bool enabled) noexcept {
  reject_inactive_reads = enabled;
}

void immersed_reconstruction_force_momentum_fallback_raw(
    bool enabled) noexcept {
  force_momentum_fallback_for_test = enabled;
}

} // namespace detail
#endif

} // namespace hundun::finite_volume
