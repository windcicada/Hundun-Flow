// SPDX-License-Identifier: Apache-2.0

#include "hundun/runtime/vtk_legacy.hpp"

#include "hundun/mesh/uniform_structured_mesh.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/field_registry.hpp"
#include "hundun/runtime/field_storage.hpp"
#include "hundun/runtime/types.hpp"

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <string>
#include <string_view>

namespace hundun::runtime {
namespace {

bool same(Int3 lhs, Int3 rhs) noexcept {
  return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

bool finite(Real3 value) noexcept {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

bool positive(Real3 value) noexcept {
  return value.x > 0.0 && value.y > 0.0 && value.z > 0.0;
}

std::uint64_t checked_cell_count(Int3 extent) {
  if (extent.x <= 0 || extent.y <= 0 || extent.z <= 0) {
    throw Error("VTK local cell extents must be positive");
  }
  const auto x = static_cast<std::uint64_t>(extent.x);
  const auto y = static_cast<std::uint64_t>(extent.y);
  const auto z = static_cast<std::uint64_t>(extent.z);
  if (x > std::numeric_limits<std::uint64_t>::max() / y ||
      x * y > std::numeric_limits<std::uint64_t>::max() / z) {
    throw Error("VTK cell count overflows uint64");
  }
  return x * y * z;
}

int checked_point_dimension(int cells) {
  if (cells <= 0 || cells == std::numeric_limits<int>::max()) {
    throw Error("VTK point dimension is outside the supported int range");
  }
  return cells + 1;
}

bool vtk_safe_name(std::string_view name) noexcept {
  if (name.empty()) {
    return false;
  }
  for (const char character : name) {
    const auto byte = static_cast<unsigned char>(character);
    if (byte <= 0x20U || byte == 0x7FU) {
      return false;
    }
  }
  return true;
}

std::string vtk_filename(std::int64_t step, int rank) {
  std::ostringstream output;
  output.imbue(std::locale::classic());
  output << "scalar.step" << std::setfill('0') << std::setw(8) << step
         << ".rank" << std::setw(6) << rank << ".vtk";
  if (!output.good()) {
    throw Error("unable to format VTK output filename");
  }
  return output.str();
}

void write_vtk_rank_impl(const std::filesystem::path &output_directory,
                         std::int64_t step, int rank,
                         const mesh::UniformStructuredMesh &mesh,
                         const FieldRegistry &registry,
                         const FieldStorage &storage, FieldId field) {
  if (step < 0 || rank < 0) {
    throw Error("VTK step and rank must be nonnegative");
  }
  if (output_directory.empty()) {
    throw Error("VTK output directory must not be empty");
  }
  if (!registry.frozen()) {
    throw Error("VTK output requires a frozen field registry");
  }
  const FieldDescriptor &descriptor = registry.descriptor(field);
  if (descriptor.space != FunctionSpace::cell_average ||
      descriptor.scalar_type != ScalarType::float64 ||
      descriptor.components != 1U ||
      (descriptor.output != OutputPolicy::selected &&
       descriptor.output != OutputPolicy::always)) {
    throw Error("Stage 1 VTK supports selected scalar Float64 cell-average "
                "fields only");
  }
  if (!vtk_safe_name(descriptor.name)) {
    throw Error("VTK field name must be nonempty and contain no ASCII space or "
                "control");
  }

  const Int3 extent = mesh.local_extent();
  if (!same(extent, storage.interior_extent())) {
    throw Error("VTK storage extent does not match the mesh local extent");
  }
  const auto values = storage.view<double>(field);
  if (!same(values.interior_extent(), extent) ||
      values.components() != descriptor.components ||
      values.ghost_width() != descriptor.ghost_width) {
    throw Error("VTK field storage layout does not match the registry");
  }
  const int point_x = checked_point_dimension(extent.x);
  const int point_y = checked_point_dimension(extent.y);
  const int point_z = checked_point_dimension(extent.z);
  const std::uint64_t cell_count = checked_cell_count(extent);
  const Real3 spacing = mesh.spacing_m();
  const Real3 first_center = mesh.cell_center(Int3{0, 0, 0});
  const Real3 origin{first_center.x - 0.5 * spacing.x,
                     first_center.y - 0.5 * spacing.y,
                     first_center.z - 0.5 * spacing.z};
  if (!finite(spacing) || !positive(spacing) || !finite(origin)) {
    throw Error("VTK origin and spacing must be finite and spacing positive");
  }

  std::filesystem::create_directories(output_directory);
  const std::filesystem::path final_path =
      output_directory / vtk_filename(step, rank);
  std::filesystem::path temporary_path = final_path;
  temporary_path += ".tmp";
  if (std::filesystem::exists(final_path) ||
      std::filesystem::exists(temporary_path)) {
    throw Error("VTK output refuses to overwrite an existing file");
  }

  std::ofstream output(temporary_path, std::ios::out | std::ios::trunc);
  if (!output.is_open()) {
    throw Error("unable to open VTK temporary file: " +
                temporary_path.string());
  }
  output.imbue(std::locale::classic());
  output << std::setprecision(17) << "# vtk DataFile Version 3.0\n"
         << "HUNDUN-FLOW Stage 1\n"
         << "ASCII\n"
         << "DATASET STRUCTURED_POINTS\n"
         << "DIMENSIONS " << point_x << ' ' << point_y << ' ' << point_z << '\n'
         << "ORIGIN " << origin.x << ' ' << origin.y << ' ' << origin.z << '\n'
         << "SPACING " << spacing.x << ' ' << spacing.y << ' ' << spacing.z
         << '\n'
         << "CELL_DATA " << cell_count << '\n'
         << "SCALARS " << descriptor.name << " double 1\n"
         << "LOOKUP_TABLE default\n";
  for (int k = 0; k < extent.z; ++k) {
    for (int j = 0; j < extent.y; ++j) {
      for (int i = 0; i < extent.x; ++i) {
        output << values(i, j, k, 0) << '\n';
      }
    }
  }
  output.flush();
  if (!output.good()) {
    throw Error("unable to write VTK temporary file: " +
                temporary_path.string());
  }
  output.close();
  if (output.fail()) {
    throw Error("unable to close VTK temporary file: " +
                temporary_path.string());
  }
  std::filesystem::rename(temporary_path, final_path);
}

} // namespace

void write_vtk_rank(const std::filesystem::path &output_directory,
                    std::int64_t step, int rank,
                    const mesh::UniformStructuredMesh &mesh,
                    const FieldRegistry &registry, const FieldStorage &storage,
                    FieldId field) {
  try {
    write_vtk_rank_impl(output_directory, step, rank, mesh, registry, storage,
                        field);
  } catch (const Error &) {
    throw;
  } catch (const std::exception &error) {
    throw Error(std::string("VTK output failed: ") + error.what());
  }
}

} // namespace hundun::runtime
