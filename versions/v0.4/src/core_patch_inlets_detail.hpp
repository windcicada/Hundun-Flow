// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
// Cold-only imported geometry and inlet binding. See the GTMC port ledger.
#pragma once

#include "hundun/v04_ibm.hpp"
#include "hundun/v04_physics.hpp"

#include <algorithm>
#include <array>
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdint>
#include <fcntl.h>
#include <limits>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

namespace hundun::v04::detail {

// Reopen and verify the immutable input at product freeze, rather than trusting
// a file that may have changed since case validation. Paths remain direct and
// symlinks are rejected on both the directory and the referenced file.
inline Status read_case_binary(const std::filesystem::path& root,
    const std::filesystem::path& name, std::uint64_t bytes,
    PlanFingerprint expected, std::vector<std::uint8_t>& out) noexcept try {
  struct Descriptor {
    int value{-1};
    ~Descriptor() { if (value >= 0) ::close(value); }
  };
  if (name.empty() || name.has_parent_path() || expected == 0U ||
      bytes > std::numeric_limits<std::size_t>::max())
    return {StatusCode::invalid_plan, 15901U};
  Descriptor directory{::open(root.c_str(), O_RDONLY | O_DIRECTORY | O_NOFOLLOW | O_CLOEXEC)};
  if (directory.value < 0) return {StatusCode::invalid_case, 15902U};
  Descriptor file{::openat(directory.value, name.c_str(), O_RDONLY | O_NOFOLLOW | O_CLOEXEC)};
  struct stat metadata {};
  if (file.value < 0 || ::fstat(file.value, &metadata) != 0 ||
      !S_ISREG(metadata.st_mode) || metadata.st_size < 0 ||
      static_cast<std::uint64_t>(metadata.st_size) != bytes)
    return {StatusCode::invalid_case, 15902U};
  std::vector<std::uint8_t> data(static_cast<std::size_t>(bytes));
  std::size_t offset = 0U;
  while (offset < data.size()) {
    const ssize_t count = ::read(file.value, data.data() + offset,
                                std::min(data.size() - offset, std::size_t{1048576U}));
    if (count < 0 && errno == EINTR) continue;
    if (count <= 0) return {StatusCode::invalid_case, 15902U};
    offset += static_cast<std::size_t>(count);
  }
  std::uint8_t extra = 0U;
  if (::read(file.value, &extra, 1U) != 0)
    return {StatusCode::invalid_case, 15902U};
  std::uint64_t hash = UINT64_C(14695981039346656037);
  for (const auto byte : data) {
    hash ^= byte;
    hash *= UINT64_C(1099511628211);
  }
  if ((hash == 0U ? 1U : hash) != expected)
    return {StatusCode::invalid_case, 15903U};
  out = std::move(data);
  return {};
} catch (...) { return {StatusCode::invalid_plan, 15904U}; }

struct CompiledPatchInlets {
  std::vector<PatchInletSpec> patches;
  // -1 denotes zero inlet velocity outside all labelled supports on a face
  // which is explicitly partitioned by patches. Other faces retain the old BC.
  std::array<std::vector<std::int32_t>, 6U> face_patch;
  std::vector<double> global_area;
  std::vector<double> inlet_enthalpy;
  std::vector<double> inlet_density;
  std::vector<std::vector<double>> independent_species;
  std::vector<IbmInterfaceInletState> immersed_states;
  // Preallocated collective scratch: no allocation in the boundary hot loop.
  std::vector<double> local_capacity;
  std::vector<double> global_capacity;
};

inline std::size_t patch_flat(Int3 cells, Int3 p) noexcept {
  return (static_cast<std::size_t>(p.z) * cells.y + p.y) * cells.x + p.x;
}

inline Status compile_patch_inlets(const ValidatedModel& model,
    const std::filesystem::path& root, const CartesianGeometryPlan& geometry,
    MeshPatch local, const EBTopology* topology,
    const ThermodynamicsPlan& thermo, CompiledPatchInlets& out) noexcept try {
  if (!model.patch_inlets) return {};
  const auto& spec = *model.patch_inlets;
  const Int3 cells = geometry.global_cells();
  const std::uint64_t count = static_cast<std::uint64_t>(cells.x) * cells.y * cells.z;
  if (cells.x <= 0 || cells.y <= 0 || cells.z <= 0 ||
      spec.patches.empty() || spec.patches.size() > 256U)
    return {StatusCode::invalid_plan, 15905U};
  std::vector<std::uint8_t> raw;
  Status status = read_case_binary(root, spec.labels_file, count * 4U,
                                  spec.labels_fingerprint, raw);
  if (!status) return status;
  const auto label_at = [&](Int3 p) -> std::int64_t {
    const std::size_t offset = patch_flat(cells, p) * 4U;
    std::uint32_t value = 0U;
    for (unsigned byte = 0U; byte < 4U; ++byte)
      value |= static_cast<std::uint32_t>(raw[offset + byte]) << (8U * byte);
    return value <= INT32_MAX ? value : static_cast<std::int64_t>(value) - INT64_C(4294967296);
  };
  const auto patch_index = [&](std::int64_t label) -> std::size_t {
    for (std::size_t i = 0U; i < spec.patches.size(); ++i)
      if (spec.patches[i].label == label) return i;
    return spec.patches.size();
  };
  const auto component = [](Int3 p, unsigned axis) {
    return axis == 0U ? p.x : (axis == 1U ? p.y : p.z);
  };
  const auto owned = [&](Int3 p) {
    return p.x >= local.begin.x && p.x < local.begin.x + local.cells.x &&
           p.y >= local.begin.y && p.y < local.begin.y + local.cells.y &&
           p.z >= local.begin.z && p.z < local.begin.z + local.cells.z;
  };
  const auto area = [&](Int3 p, unsigned face) {
    if (face < 2U) return geometry.y().widths().data[p.y] * geometry.z().widths().data[p.z];
    if (face < 4U) return geometry.x().widths().data[p.x] * geometry.z().widths().data[p.z];
    return geometry.x().widths().data[p.x] * geometry.y().widths().data[p.y];
  };
  CompiledPatchInlets plan;
  plan.patches = spec.patches;
  plan.global_area.assign(spec.patches.size(), 0.0);
  plan.local_capacity.assign(spec.patches.size(), 0.0);
  plan.global_capacity.assign(spec.patches.size(), 0.0);
  plan.inlet_enthalpy.resize(spec.patches.size());
  plan.inlet_density.resize(spec.patches.size());
  plan.independent_species.resize(spec.patches.size());
  std::array<double, 6U> face_mass{};
  for (std::size_t i = 0U; i < spec.patches.size(); ++i) {
    const auto& p = spec.patches[i];
    if (p.face >= 6U || p.label <= 0 || !(p.boundary.mass_flow_rate > 0.0))
      return {StatusCode::invalid_plan, 15905U};
    if (!p.immersed) {
      const auto& parent = model.boundaries[p.face];
      // The present outer scalar BC has one target per face. Fail closed for
      // heterogeneous outer thermodynamic targets until that plan is extended.
      if (parent.temperature != p.boundary.temperature ||
          parent.thermal_kind != p.boundary.thermal_kind ||
          parent.scalars.size() != p.boundary.scalars.size())
        return {StatusCode::invalid_case, 15906U};
      for (const auto& scalar : p.boundary.scalars) {
        const auto found = std::find_if(parent.scalars.begin(), parent.scalars.end(),
            [&](const auto& value) { return value.stable_name == scalar.stable_name; });
        if (found == parent.scalars.end() || found->kind != scalar.kind ||
            found->value != scalar.value)
          return {StatusCode::invalid_case, 15906U};
      }
      face_mass[p.face] += p.boundary.mass_flow_rate;
      const std::size_t face_cells = p.face < 2U ?
          static_cast<std::size_t>(local.cells.y) * local.cells.z :
          (p.face < 4U ? static_cast<std::size_t>(local.cells.x) * local.cells.z :
                        static_cast<std::size_t>(local.cells.x) * local.cells.y);
      plan.face_patch[p.face].assign(face_cells, -1);
    } else if (topology == nullptr) {
      return {StatusCode::invalid_case, 15905U};
    }
    for (const auto& scalar : model.transported_scalars) {
      const auto found = std::find_if(p.boundary.scalars.begin(), p.boundary.scalars.end(),
          [&](const auto& value) { return value.stable_name == scalar.stable_name; });
      if (found == p.boundary.scalars.end() || found->kind != ScalarBoundaryKind::dirichlet ||
          (p.immersed && scalar.role != TransportedScalarRole::species))
        return {StatusCode::invalid_case, 15906U};
      if (scalar.role == TransportedScalarRole::species)
        plan.independent_species[i].push_back(found->value);
    }
    double cp = 0.0, gas = 0.0;
    const auto& y = plan.independent_species[i];
    status = thermo.mixture_enthalpy(p.boundary.temperature, {y.data(), y.size()},
                                    plan.inlet_enthalpy[i], cp, gas);
    if (!status) return status;
    plan.inlet_density[i] = p.boundary.pressure / (gas * p.boundary.temperature);
    if (!(plan.inlet_density[i] > 0.0) || !std::isfinite(plan.inlet_density[i]))
      return {StatusCode::invalid_case, 15906U};
  }
  for (unsigned face = 0U; face < 6U; ++face)
    if (!plan.face_patch[face].empty() &&
        std::abs(face_mass[face] - model.boundaries[face].mass_flow_rate) >
            64.0 * std::numeric_limits<double>::epsilon() * face_mass[face])
      return {StatusCode::invalid_case, 15907U};
  std::vector<std::uint64_t> owned_sources(spec.patches.size(), 0U), bound_sources(spec.patches.size(), 0U);
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 p{x, y, z};
        const auto label = label_at(p);
        if (label == 0) continue;
        const auto index = patch_index(label > 0 ? label : -label);
        if (index == spec.patches.size()) return {StatusCode::invalid_case, 15908U};
        const auto& inlet = spec.patches[index];
        if (label < 0) {
          if (!inlet.immersed || (owned(p) && topology->is_fluid_global(p)))
            return {StatusCode::invalid_case, 15908U};
          continue;
        }
        if (owned(p) && topology != nullptr && !topology->is_fluid_global(p))
          return {StatusCode::invalid_case, 15908U};
        const unsigned axis = inlet.face / 2U;
        if (inlet.immersed) {
          Int3 neighbor = p;
          const int shift = inlet.face % 2U == 0U ? -1 : 1;
          if (axis == 0U) neighbor.x += shift;
          else if (axis == 1U) neighbor.y += shift;
          else neighbor.z += shift;
          if (component(neighbor, axis) < 0 || component(neighbor, axis) >= component(cells, axis) ||
              label_at(neighbor) != -label ||
              (owned(p) && topology->is_fluid_global(neighbor)))
            return {StatusCode::invalid_case, 15908U};
          if (owned(p)) ++owned_sources[index];
        } else {
          if (component(p, axis) != (inlet.face % 2U == 0U ? 0 : component(cells, axis) - 1))
            return {StatusCode::invalid_case, 15908U};
          if (owned(p)) {
            const Int3 q{x - local.begin.x, y - local.begin.y, z - local.begin.z};
            const std::size_t offset = inlet.face < 2U ?
                static_cast<std::size_t>(q.z) * local.cells.y + q.y :
                (inlet.face < 4U ? static_cast<std::size_t>(q.z) * local.cells.x + q.x :
                                  static_cast<std::size_t>(q.y) * local.cells.x + q.x);
            plan.face_patch[inlet.face][offset] = static_cast<std::int32_t>(index);
          }
        }
        plan.global_area[index] += area(p, inlet.face);
      }
  for (const auto value : plan.global_area)
    if (!(value > 0.0) || !std::isfinite(value)) return {StatusCode::invalid_case, 15908U};
  if (topology != nullptr) {
    const auto links = topology->links();
    for (std::size_t l = 0U; l < links.size; ++l) {
      const auto& link = links.data[l];
      const auto index = patch_index(label_at(link.fluid_global_index));
      if (index == spec.patches.size() || !spec.patches[index].immersed) continue;
      const auto& inlet = spec.patches[index];
      if (static_cast<unsigned>(link.direction) != inlet.face) continue;
      ++bound_sources[index];
      const auto& y = plan.independent_species[index];
      const double inward = inlet.face < 2U ? inlet.boundary.direction.x :
          (inlet.face < 4U ? inlet.boundary.direction.y : inlet.boundary.direction.z);
      const double sign = inlet.face % 2U == 0U ? 1.0 : -1.0;
      const double speed = inlet.boundary.mass_flow_rate /
          (plan.inlet_density[index] * plan.global_area[index] * sign * inward);
      if (!(speed > 0.0) || !std::isfinite(speed)) return {StatusCode::invalid_case, 15908U};
      plan.immersed_states.push_back({link.global_link,
          sign * inlet.boundary.mass_flow_rate * link.cartesian_control_face_area / plan.global_area[index],
          {speed * inlet.boundary.direction.x, speed * inlet.boundary.direction.y,
           speed * inlet.boundary.direction.z}, plan.inlet_enthalpy[index], {y.data(), y.size()}});
    }
  }
  if (owned_sources != bound_sources) return {StatusCode::invalid_case, 15908U};
  out = std::move(plan);
  return {};
} catch (...) { return {StatusCode::invalid_plan, 15904U}; }

}  // namespace hundun::v04::detail
