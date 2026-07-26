// SPDX-License-Identifier: Apache-2.0

#include "hundun/diagnostics/mesh_diagnostic_v2.hpp"

#include "hundun/runtime/collective_status.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/mpi_operation_error.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <locale>
#include <map>
#include <sstream>
#include <string>
#include <system_error>
#include <type_traits>

namespace hundun::diagnostics {
namespace {

constexpr std::array<char, 8> kMagic{'H', 'F', 'M', 'E', 'S', 'H', 'V', '2'};
constexpr std::uint32_t kVersion = 2U;
constexpr std::uint32_t kEndian = 0x01020304U;

[[noreturn]] void fail(std::string message) {
  throw runtime::Error("meshdiag v2: " + std::move(message));
}

std::uint64_t checked_add(std::uint64_t left, std::uint64_t right,
                          std::string_view subject) {
  if (left > std::numeric_limits<std::uint64_t>::max() - right)
    fail(std::string(subject) + " addition overflows u64");
  return left + right;
}

std::uint64_t checked_multiply(std::uint64_t left, std::uint64_t right,
                               std::string_view subject) {
  if (right != 0U &&
      left > std::numeric_limits<std::uint64_t>::max() / right)
    fail(std::string(subject) + " multiplication overflows u64");
  return left * right;
}

void require_directory_agreement(const runtime::MpiContext& mpi,
                                 const std::filesystem::path& directory) {
  std::string text;
  bool local_ok = true;
  std::string message;
  try {
    auto resolved =
        std::filesystem::absolute(directory).lexically_normal();
    if (resolved.filename().empty() &&
        resolved.parent_path() != resolved)
      resolved = resolved.parent_path();
    text = resolved.generic_string();
  } catch (const std::exception& error) {
    local_ok = false;
    message = error.what();
  }
  auto status = runtime::collective_status(
      mpi, local_ok,
      message.empty() ? "unable to resolve meshdiag output path" : message);
  if (!status.ok)
    fail(status.message);
  std::uint64_t root_size =
      mpi.rank() == 0 ? static_cast<std::uint64_t>(text.size()) : 0U;
  runtime::check_mpi_result(
      MPI_Bcast(&root_size, 1, MPI_UINT64_T, 0, mpi.comm()),
      "MPI_Bcast meshdiag directory size");
  local_ok =
      root_size <= static_cast<std::uint64_t>(
                       std::numeric_limits<std::size_t>::max()) &&
      root_size <= static_cast<std::uint64_t>(
                       std::numeric_limits<int>::max());
  std::string root;
  if (local_ok) {
    try {
      root.resize(static_cast<std::size_t>(root_size));
    } catch (...) {
      local_ok = false;
    }
  }
  status = runtime::collective_status(
      mpi, local_ok, "unable to compare meshdiag output paths");
  if (!status.ok)
    fail(status.message);
  if (mpi.rank() == 0)
    std::copy(text.begin(), text.end(), root.begin());
  runtime::check_mpi_result(
      MPI_Bcast(root.data(), static_cast<int>(root_size), MPI_BYTE, 0,
                mpi.comm()),
      "MPI_Bcast meshdiag directory bytes");
  status = runtime::collective_status(
      mpi, text == root, "meshdiag output directory differs across ranks");
  if (!status.ok)
    fail(status.message);
}

template <class Unsigned>
void append_unsigned(std::vector<std::byte>& out, Unsigned value) {
  static_assert(std::is_unsigned_v<Unsigned>);
  for (std::size_t byte = 0; byte < sizeof(Unsigned); ++byte)
    out.push_back(static_cast<std::byte>(
        (value >> (8U * static_cast<unsigned>(byte))) &
        static_cast<Unsigned>(0xffU)));
}

void append_i32(std::vector<std::byte>& out, int value) {
  static_assert(sizeof(int) >= sizeof(std::int32_t));
  if (value < std::numeric_limits<std::int32_t>::min() ||
      value > std::numeric_limits<std::int32_t>::max())
    fail("integer outside i32 domain");
  std::uint32_t bits{};
  const auto narrowed = static_cast<std::int32_t>(value);
  std::memcpy(&bits, &narrowed, sizeof(bits));
  append_unsigned(out, bits);
}

void append_f64(std::vector<std::byte>& out, double value) {
  std::uint64_t bits{};
  std::memcpy(&bits, &value, sizeof(bits));
  append_unsigned(out, bits);
}

std::uint64_t crc64(const std::byte* data, std::size_t size) noexcept {
  constexpr std::uint64_t polynomial = 0x42F0E1EBA9EA3693ULL;
  std::uint64_t result{};
  for (std::size_t index = 0; index < size; ++index) {
    result ^= static_cast<std::uint64_t>(
                  static_cast<unsigned char>(data[index]))
              << 56U;
    for (int bit = 0; bit < 8; ++bit)
      result = (result & (1ULL << 63U)) != 0U
                   ? (result << 1U) ^ polynomial
                   : result << 1U;
  }
  return result;
}

double norm(runtime::Real3 value) noexcept {
  return std::hypot(value.x, std::hypot(value.y, value.z));
}

std::uint64_t global_cell_id(runtime::Int3 extent,
                             runtime::Int3 cell) {
  if (cell.x < 0 || cell.y < 0 || cell.z < 0 || cell.x >= extent.x ||
      cell.y >= extent.y || cell.z >= extent.z)
    fail("logical cell is outside the mesh");
  const auto nx = static_cast<std::uint64_t>(extent.x);
  const auto ny = static_cast<std::uint64_t>(extent.y);
  const auto row = checked_add(
      checked_multiply(static_cast<std::uint64_t>(cell.z), ny,
                       "global cell ID"),
      static_cast<std::uint64_t>(cell.y), "global cell ID");
  return checked_add(checked_multiply(row, nx, "global cell ID"),
                     static_cast<std::uint64_t>(cell.x),
                     "global cell ID");
}

std::uint64_t global_face_id(runtime::Int3 extent, mesh::FaceAxis axis,
                             runtime::Int3 logical) {
  const auto nx = static_cast<std::uint64_t>(extent.x);
  const auto ny = static_cast<std::uint64_t>(extent.y);
  const auto nz = static_cast<std::uint64_t>(extent.z);
  const auto x_count = checked_multiply(
      checked_multiply(checked_add(nx, 1U, "x face extent"), ny,
                       "x face count"),
      nz, "x face count");
  const auto y_count = checked_multiply(
      checked_multiply(nx, checked_add(ny, 1U, "y face extent"),
                       "y face count"),
      nz, "y face count");
  const auto i = static_cast<std::uint64_t>(logical.x);
  const auto j = static_cast<std::uint64_t>(logical.y);
  const auto k = static_cast<std::uint64_t>(logical.z);
  if (axis == mesh::FaceAxis::x) {
    if (logical.x < 0 || logical.x > extent.x || logical.y < 0 ||
        logical.y >= extent.y || logical.z < 0 || logical.z >= extent.z)
      fail("x face logical coordinate is invalid");
    const auto row =
        checked_add(checked_multiply(k, ny, "x face ID"), j, "x face ID");
    return checked_add(
        checked_multiply(row, checked_add(nx, 1U, "x face extent"),
                         "x face ID"),
        i, "x face ID");
  }
  if (axis == mesh::FaceAxis::y) {
    if (logical.x < 0 || logical.x >= extent.x || logical.y < 0 ||
        logical.y > extent.y || logical.z < 0 || logical.z >= extent.z)
      fail("y face logical coordinate is invalid");
    const auto row = checked_add(
        checked_multiply(k, checked_add(ny, 1U, "y face extent"),
                         "y face ID"),
        j, "y face ID");
    return checked_add(
        x_count,
        checked_add(checked_multiply(row, nx, "y face ID"), i,
                    "y face ID"),
        "y face ID");
  }
  if (axis == mesh::FaceAxis::z) {
    if (logical.x < 0 || logical.x >= extent.x || logical.y < 0 ||
        logical.y >= extent.y || logical.z < 0 || logical.z > extent.z)
      fail("z face logical coordinate is invalid");
    const auto row =
        checked_add(checked_multiply(k, ny, "z face ID"), j, "z face ID");
    const auto local =
        checked_add(checked_multiply(row, nx, "z face ID"), i, "z face ID");
    return checked_add(checked_add(x_count, y_count, "z face offset"), local,
                       "z face ID");
  }
  fail("face axis is invalid");
}

bool contains(runtime::Box3 box, runtime::Int3 cell) noexcept {
  return cell.x >= box.begin.x && cell.x < box.end.x &&
         cell.y >= box.begin.y && cell.y < box.end.y &&
         cell.z >= box.begin.z && cell.z < box.end.z;
}

struct ExpectedFace final {
  std::uint64_t global_id{};
  mesh::FaceAxis axis{};
  runtime::Int3 logical{};
  std::uint64_t owner{};
  std::optional<std::uint64_t> neighbour;
  std::optional<std::uint32_t> patch;
  std::optional<std::uint64_t> periodic_pair;
};

ExpectedFace expected_face(runtime::Int3 extent, runtime::Box3 owned_box,
                           const std::array<bool, 3>& periodic,
                           mesh::FaceAxis axis, runtime::Int3 logical) {
  static_cast<void>(global_face_id(extent, axis, logical));
  runtime::Int3 lower = logical;
  runtime::Int3 upper = logical;
  int plane{};
  int cells{};
  std::size_t axis_index{};
  if (axis == mesh::FaceAxis::x) {
    --lower.x;
    plane = logical.x;
    cells = extent.x;
    axis_index = 0U;
  } else if (axis == mesh::FaceAxis::y) {
    --lower.y;
    plane = logical.y;
    cells = extent.y;
    axis_index = 1U;
  } else {
    --lower.z;
    plane = logical.z;
    cells = extent.z;
    axis_index = 2U;
  }
  runtime::Int3 owner =
      plane == 0 ? upper : lower;
  std::optional<runtime::Int3> neighbour;
  if (plane > 0 && plane < cells) {
    neighbour = upper;
  } else if (periodic[axis_index]) {
    neighbour = owner;
    if (axis == mesh::FaceAxis::x)
      neighbour->x = plane == 0 ? extent.x - 1 : 0;
    else if (axis == mesh::FaceAxis::y)
      neighbour->y = plane == 0 ? extent.y - 1 : 0;
    else
      neighbour->z = plane == 0 ? extent.z - 1 : 0;
  }
  std::optional<std::uint32_t> patch;
  if (plane == 0)
    patch = static_cast<std::uint32_t>(2U * axis_index);
  else if (plane == cells)
    patch = static_cast<std::uint32_t>(2U * axis_index + 1U);
  std::optional<std::uint64_t> pair;
  if (patch && periodic[axis_index]) {
    auto opposite = logical;
    if (axis == mesh::FaceAxis::x)
      opposite.x = plane == 0 ? extent.x : 0;
    else if (axis == mesh::FaceAxis::y)
      opposite.y = plane == 0 ? extent.y : 0;
    else
      opposite.z = plane == 0 ? extent.z : 0;
    pair = global_face_id(extent, axis, opposite);
  }
  ExpectedFace result{global_face_id(extent, axis, logical),
                      axis,
                      logical,
                      global_cell_id(extent, owner),
                      std::nullopt,
                      patch,
                      pair};
  if (neighbour)
    result.neighbour = global_cell_id(extent, *neighbour);
  if (!contains(owned_box, owner))
    result.global_id = std::numeric_limits<std::uint64_t>::max();
  return result;
}

bool patch_matches_logical_face(std::uint32_t patch,
                                const MeshDiagnosticFaceV2& face,
                                runtime::Int3 extent) noexcept {
  switch (patch) {
  case 0U:
    return face.axis == mesh::FaceAxis::x && face.logical.x == 0;
  case 1U:
    return face.axis == mesh::FaceAxis::x &&
           face.logical.x == extent.x;
  case 2U:
    return face.axis == mesh::FaceAxis::y && face.logical.y == 0;
  case 3U:
    return face.axis == mesh::FaceAxis::y &&
           face.logical.y == extent.y;
  case 4U:
    return face.axis == mesh::FaceAxis::z && face.logical.z == 0;
  case 5U:
    return face.axis == mesh::FaceAxis::z &&
           face.logical.z == extent.z;
  default:
    return false;
  }
}

bool finite(runtime::Real3 value) noexcept {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

bool negative_within_tolerance(runtime::Real3 owner,
                               runtime::Real3 neighbour) noexcept {
  const double scale =
      std::abs(owner.x) + std::abs(owner.y) + std::abs(owner.z) +
      std::abs(neighbour.x) + std::abs(neighbour.y) +
      std::abs(neighbour.z);
  const double error =
      norm({owner.x + neighbour.x, owner.y + neighbour.y,
            owner.z + neighbour.z});
  return error <= 256.0 * std::numeric_limits<double>::epsilon() *
                      std::max(1.0, scale);
}

void require_positive_extent(runtime::Int3 extent) {
  if (extent.x < 1 || extent.y < 1 || extent.z < 1)
    fail("global extent must be positive");
}

void require_valid_box(runtime::Box3 box, runtime::Int3 extent) {
  if (box.begin.x < 0 || box.begin.y < 0 || box.begin.z < 0 ||
      box.end.x < box.begin.x || box.end.y < box.begin.y ||
      box.end.z < box.begin.z || box.end.x > extent.x ||
      box.end.y > extent.y || box.end.z > extent.z)
    fail("owned box is outside global extent");
}

std::array<runtime::Int3, 4> face_vertices(mesh::FaceAxis axis,
                                          runtime::Int3 c) {
  if (axis == mesh::FaceAxis::x)
    return {{{c.x, c.y, c.z},
             {c.x, c.y + 1, c.z},
             {c.x, c.y + 1, c.z + 1},
             {c.x, c.y, c.z + 1}}};
  if (axis == mesh::FaceAxis::y)
    return {{{c.x, c.y, c.z},
             {c.x, c.y, c.z + 1},
             {c.x + 1, c.y, c.z + 1},
             {c.x + 1, c.y, c.z}}};
  return {{{c.x, c.y, c.z},
           {c.x + 1, c.y, c.z},
           {c.x + 1, c.y + 1, c.z},
           {c.x, c.y + 1, c.z}}};
}

void validate(const MeshDiagnosticV2& value) {
  if (value.rank < 0 || value.rank_count < 1 || value.rank >= value.rank_count)
    fail("rank metadata is invalid");
  require_positive_extent(value.global_extent);
  require_valid_box(value.owned_box, value.global_extent);
  if (!finite(value.origin_m) || !finite(value.length_m) ||
      value.length_m.x <= 0.0 || value.length_m.y <= 0.0 ||
      value.length_m.z <= 0.0)
    fail("mapping extents are invalid");
  if (value.mapping_kind != mesh::MappingKind::uniform_box &&
      value.mapping_kind != mesh::MappingKind::analytic_warped_box)
    fail("mapping kind is invalid");

  std::uint64_t previous{};
  bool first = true;
  std::map<std::uint64_t, runtime::Int3> vertices;
  for (const auto& vertex : value.vertices) {
    if ((!first && vertex.global_id <= previous) ||
        vertex.global_id != mesh_diagnostic_global_vertex_id(
                                value.global_extent, vertex.logical) ||
        !finite(vertex.position_m))
      fail("vertex table is not canonical");
    first = false;
    previous = vertex.global_id;
    vertices.emplace(vertex.global_id, vertex.logical);
  }

  double volume_sum{};
  double maximum_closure{};
  std::vector<std::uint64_t> expected_cells;
  for (int z = value.owned_box.begin.z; z < value.owned_box.end.z; ++z)
    for (int y = value.owned_box.begin.y; y < value.owned_box.end.y; ++y)
      for (int x = value.owned_box.begin.x; x < value.owned_box.end.x; ++x)
        expected_cells.push_back(
            global_cell_id(value.global_extent, {x, y, z}));
  if (value.cells.size() != expected_cells.size())
    fail("cell table does not match the owned box");
  first = true;
  previous = 0U;
  for (std::size_t index = 0; index < value.cells.size(); ++index) {
    const auto& cell = value.cells[index];
    if ((!first && cell.global_id <= previous) || !finite(cell.centre_m) ||
        !finite(cell.closure_m2) || !std::isfinite(cell.volume_m3) ||
        !std::isfinite(cell.minimum_jacobian_m3) || cell.volume_m3 <= 0.0 ||
        cell.minimum_jacobian_m3 <= 0.0 ||
        cell.global_id != expected_cells[index])
      fail("cell table is invalid");
    first = false;
    previous = cell.global_id;
    volume_sum += cell.volume_m3;
    maximum_closure = std::max(maximum_closure, norm(cell.closure_m2));
  }
  if (volume_sum != value.owned_volume_sum ||
      maximum_closure != value.maximum_cell_closure_norm)
    fail("cell conservation summary differs from records");

  std::uint64_t reciprocal_checks{};
  std::uint64_t reciprocal_failures{};
  std::array<bool, 3> periodic{};
  for (const auto& face : value.faces) {
    if (face.axis != mesh::FaceAxis::x &&
        face.axis != mesh::FaceAxis::y &&
        face.axis != mesh::FaceAxis::z)
      fail("face axis is invalid");
    if (face.periodic_pair)
      periodic[static_cast<std::size_t>(face.axis)] = true;
  }
  std::map<std::uint64_t, ExpectedFace> expected_faces;
  const auto add_expected = [&](mesh::FaceAxis axis,
                                runtime::Int3 logical) {
    const auto expected = expected_face(value.global_extent, value.owned_box,
                                        periodic, axis, logical);
    if (expected.global_id != std::numeric_limits<std::uint64_t>::max())
      expected_faces.emplace(expected.global_id, expected);
  };
  for (int z = value.owned_box.begin.z; z < value.owned_box.end.z; ++z)
    for (int y = value.owned_box.begin.y; y < value.owned_box.end.y; ++y)
      for (int x = value.owned_box.begin.x; x < value.owned_box.end.x; ++x) {
        add_expected(mesh::FaceAxis::x, {x, y, z});
        add_expected(mesh::FaceAxis::x, {x + 1, y, z});
        add_expected(mesh::FaceAxis::y, {x, y, z});
        add_expected(mesh::FaceAxis::y, {x, y + 1, z});
        add_expected(mesh::FaceAxis::z, {x, y, z});
        add_expected(mesh::FaceAxis::z, {x, y, z + 1});
      }
  if (value.faces.size() != expected_faces.size())
    fail("face table does not match the complete owned face set");

  std::map<std::uint64_t, runtime::Int3> referenced_vertices;
  for (int z = value.owned_box.begin.z; z <= value.owned_box.end.z; ++z)
    for (int y = value.owned_box.begin.y; y <= value.owned_box.end.y; ++y)
      for (int x = value.owned_box.begin.x; x <= value.owned_box.end.x; ++x) {
        const runtime::Int3 logical{x, y, z};
        referenced_vertices.emplace(
            mesh_diagnostic_global_vertex_id(value.global_extent, logical),
            logical);
      }
  first = true;
  previous = 0U;
  const std::uint64_t global_cells = checked_multiply(
      checked_multiply(static_cast<std::uint64_t>(value.global_extent.x),
                       static_cast<std::uint64_t>(value.global_extent.y),
                       "global cell count"),
      static_cast<std::uint64_t>(value.global_extent.z),
      "global cell count");
  for (const auto& face : value.faces) {
    if ((!first && face.global_id <= previous) ||
        (face.axis != mesh::FaceAxis::x &&
         face.axis != mesh::FaceAxis::y &&
         face.axis != mesh::FaceAxis::z) ||
        !finite(face.centre_m) || !finite(face.owner_area_vector_m2) ||
        !std::isfinite(face.area_m2) || face.area_m2 <= 0.0 ||
        !std::isfinite(face.skewness) || face.skewness < 0.0 ||
        !std::isfinite(face.non_orthogonality_degrees) ||
        face.non_orthogonality_degrees < 0.0 ||
        face.non_orthogonality_degrees > 180.0 ||
        face.neighbour_global_cell.has_value() !=
            face.neighbour_area_vector_m2.has_value() ||
        face.owner_global_cell >= global_cells ||
        (face.neighbour_global_cell.has_value() &&
         (*face.neighbour_global_cell >= global_cells ||
          *face.neighbour_global_cell == face.owner_global_cell)) ||
        (face.patch_id.has_value() && *face.patch_id > 5U) ||
        (face.patch_id.has_value() &&
         face.neighbour_global_cell.has_value() &&
         !face.periodic_pair.has_value()) ||
        (face.periodic_pair.has_value() && !face.patch_id.has_value()))
      fail("face table is invalid");
    const auto expected = expected_faces.find(face.global_id);
    if (expected == expected_faces.end() ||
        face.global_id !=
            global_face_id(value.global_extent, face.axis, face.logical) ||
        expected->second.axis != face.axis ||
        expected->second.logical.x != face.logical.x ||
        expected->second.logical.y != face.logical.y ||
        expected->second.logical.z != face.logical.z ||
        expected->second.owner != face.owner_global_cell ||
        expected->second.neighbour != face.neighbour_global_cell ||
        expected->second.patch != face.patch_id ||
        expected->second.periodic_pair != face.periodic_pair)
      fail("face logical topology is invalid");
    if (face.patch_id &&
        !patch_matches_logical_face(*face.patch_id, face,
                                    value.global_extent))
      fail("face patch relation is invalid");
    const double area_vector_norm = norm(face.owner_area_vector_m2);
    const double area_tolerance =
        256.0 * std::numeric_limits<double>::epsilon() *
        std::max({1.0, area_vector_norm, face.area_m2});
    if (std::abs(area_vector_norm - face.area_m2) > area_tolerance)
      fail("face area differs from area-vector norm");
    const auto logical_vertices = face_vertices(face.axis, face.logical);
    for (std::size_t corner = 0; corner < logical_vertices.size(); ++corner) {
      const auto expected_vertex = mesh_diagnostic_global_vertex_id(
          value.global_extent, logical_vertices[corner]);
      if (face.vertex_ids[corner] != expected_vertex ||
          vertices.find(expected_vertex) == vertices.end())
        fail("face vertex ordering is invalid");
      referenced_vertices.emplace(expected_vertex,
                                  logical_vertices[corner]);
    }
    if (face.neighbour_area_vector_m2) {
      ++reciprocal_checks;
      if (!finite(*face.neighbour_area_vector_m2) ||
          !negative_within_tolerance(face.owner_area_vector_m2,
                                     *face.neighbour_area_vector_m2))
        ++reciprocal_failures;
    }
    first = false;
    previous = face.global_id;
  }
  if (reciprocal_checks != value.reciprocal_check_count ||
      reciprocal_failures != value.reciprocal_failure_count)
    fail("face reciprocity summary differs from records");
  if (reciprocal_failures != 0U)
    fail("face reciprocity check failed");
  double maximum_area{};
  for (const auto& face : value.faces)
    maximum_area = std::max(maximum_area, face.area_m2);
  const double closure_limit =
      256.0 * std::numeric_limits<double>::epsilon() *
      std::max(1.0, 6.0 * maximum_area);
  if (value.maximum_cell_closure_norm > closure_limit)
    fail("cell closure exceeds the accepted geometry tolerance");
  const bool same_vertices =
      vertices.size() == referenced_vertices.size() &&
      std::equal(vertices.begin(), vertices.end(), referenced_vertices.begin(),
                 [](const auto& left, const auto& right) {
                   return left.first == right.first &&
                          left.second.x == right.second.x &&
                          left.second.y == right.second.y &&
                          left.second.z == right.second.z;
                 });
  if (!same_vertices)
    fail("vertex table is not the referenced canonical union");
}

class Reader final {
 public:
  explicit Reader(const std::vector<std::byte>& bytes) : bytes_(bytes) {}

  template <class Unsigned>
  Unsigned unsigned_value() {
    static_assert(std::is_unsigned_v<Unsigned>);
    require(sizeof(Unsigned));
    Unsigned result{};
    for (std::size_t byte = 0; byte < sizeof(Unsigned); ++byte)
      result |= static_cast<Unsigned>(
                    static_cast<unsigned char>(bytes_[position_ + byte]))
                << (8U * static_cast<unsigned>(byte));
    position_ += sizeof(Unsigned);
    return result;
  }

  std::uint8_t u8() { return unsigned_value<std::uint8_t>(); }
  std::uint32_t u32() { return unsigned_value<std::uint32_t>(); }
  std::uint64_t u64() { return unsigned_value<std::uint64_t>(); }
  int i32() {
    const auto raw = u32();
    std::int32_t value{};
    std::memcpy(&value, &raw, sizeof(value));
    return value;
  }
  double f64() {
    const auto raw = u64();
    double value{};
    std::memcpy(&value, &raw, sizeof(value));
    return value;
  }
  runtime::Int3 int3() { return {i32(), i32(), i32()}; }
  runtime::Real3 real3() { return {f64(), f64(), f64()}; }

  template <class T, class Read>
  std::optional<T> optional(Read&& read) {
    const auto tag = u8();
    if (tag == 0U)
      return std::nullopt;
    if (tag != 1U)
      fail("invalid option tag");
    return read();
  }

  std::size_t count(std::size_t minimum_record_bytes) {
    const auto raw = u64();
    if (raw > static_cast<std::uint64_t>(
                  std::numeric_limits<std::size_t>::max()) ||
        (minimum_record_bytes != 0U &&
         raw > static_cast<std::uint64_t>(
                   (bytes_.size() - position_) / minimum_record_bytes)))
      fail("record count exceeds remaining bytes");
    return static_cast<std::size_t>(raw);
  }

  std::size_t position() const noexcept { return position_; }

 private:
  void require(std::size_t count) {
    if (count > bytes_.size() - position_)
      fail("truncated input");
  }
  const std::vector<std::byte>& bytes_;
  std::size_t position_{};
};

std::filesystem::path mesh_path(const std::filesystem::path& directory,
                                int rank) {
  std::ostringstream name;
  name.imbue(std::locale::classic());
  name << "meshdiag.v2.rank-" << std::setw(6) << std::setfill('0') << rank
       << ".bin";
  return directory / name.str();
}

void write_bytes(const std::filesystem::path& path,
                 const std::vector<std::byte>& bytes) {
  const auto temporary = std::filesystem::path(path.string() + ".tmp");
  std::filesystem::create_directories(path.parent_path());
  {
    std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
    if (!stream)
      fail("unable to open temporary output");
    stream.write(reinterpret_cast<const char*>(bytes.data()),
                 static_cast<std::streamsize>(bytes.size()));
    stream.flush();
    if (!stream)
      fail("unable to write temporary output");
  }
  std::error_code error;
  std::filesystem::rename(temporary, path, error);
  if (error)
    fail("unable to publish output");
}

}  // namespace

std::uint64_t mesh_diagnostic_global_vertex_id(
    runtime::Int3 extent, runtime::Int3 vertex) {
  require_positive_extent(extent);
  if (vertex.x < 0 || vertex.y < 0 || vertex.z < 0 ||
      vertex.x > extent.x || vertex.y > extent.y || vertex.z > extent.z)
    fail("logical vertex is outside the mesh");
  const auto nx = static_cast<std::uint64_t>(extent.x) + 1U;
  const auto ny = static_cast<std::uint64_t>(extent.y) + 1U;
  const auto x = static_cast<std::uint64_t>(vertex.x);
  const auto y = static_cast<std::uint64_t>(vertex.y);
  const auto z = static_cast<std::uint64_t>(vertex.z);
  const auto row = checked_add(
      checked_multiply(z, ny, "global vertex ID"), y,
      "global vertex ID");
  return checked_add(checked_multiply(row, nx, "global vertex ID"), x,
                     "global vertex ID");
}

MeshDiagnosticV2 make_mesh_diagnostic_v2(
    int rank, int rank_count, const mesh::MeshTopology& topology,
    const mesh::MeshGeometry& geometry) {
  geometry.require_compatible(topology);
  MeshDiagnosticV2 result;
  result.rank = rank;
  result.rank_count = rank_count;
  result.global_extent = topology.global_extent();
  result.owned_box = topology.owned_global_box();
  result.mapping_kind = geometry.mapping_kind();
  result.origin_m = geometry.origin_m();
  result.length_m = geometry.length_m();

  std::map<std::uint64_t, runtime::Int3> vertices;
  for (mesh::LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    const auto c = topology.global_cell(cell);
    for (int dz = 0; dz < 2; ++dz)
      for (int dy = 0; dy < 2; ++dy)
        for (int dx = 0; dx < 2; ++dx) {
          const runtime::Int3 v{c.x + dx, c.y + dy, c.z + dz};
          vertices.emplace(mesh_diagnostic_global_vertex_id(
                               result.global_extent, v),
                           v);
        }
    const auto closure = geometry.cell_closure_m2(cell);
    result.cells.push_back(
        {topology.global_cell_id(cell), geometry.cell_center_m(cell),
         geometry.cell_volume_m3(cell),
         geometry.minimum_jacobian_determinant_m3(cell), closure});
    result.owned_volume_sum += result.cells.back().volume_m3;
    result.maximum_cell_closure_norm =
        std::max(result.maximum_cell_closure_norm, norm(closure));
  }
  std::sort(result.cells.begin(), result.cells.end(),
            [](const auto& left, const auto& right) {
              return left.global_id < right.global_id;
            });

  for (mesh::LocalFaceId face = 0; face < topology.local_face_count(); ++face) {
    if (topology.face_ownership(face) != mesh::EntityOwnership::owned)
      continue;
    const auto logical = topology.logical_face(face);
    MeshDiagnosticFaceV2 item;
    item.global_id = topology.global_face_id(face);
    item.axis = logical.axis;
    item.logical = logical.coordinate;
    item.owner_global_cell = topology.global_cell_id(topology.owner(face));
    if (const auto neighbour = topology.neighbour(face)) {
      item.neighbour_global_cell = topology.global_cell_id(*neighbour);
      item.neighbour_area_vector_m2 =
          geometry.face_area_vector_m2(face, mesh::FaceSide::neighbour);
      ++result.reciprocal_check_count;
      if (!negative_within_tolerance(
              geometry.face_area_vector_m2(face, mesh::FaceSide::owner),
              *item.neighbour_area_vector_m2))
        ++result.reciprocal_failure_count;
    }
    item.patch_id = topology.patch_id(face);
    item.periodic_pair = topology.periodic_pair(face);
    const auto logical_vertices = face_vertices(item.axis, item.logical);
    for (std::size_t corner = 0; corner < 4U; ++corner) {
      item.vertex_ids[corner] = mesh_diagnostic_global_vertex_id(
          result.global_extent, logical_vertices[corner]);
      vertices.emplace(item.vertex_ids[corner], logical_vertices[corner]);
    }
    item.centre_m = geometry.face_center_m(face);
    item.owner_area_vector_m2 =
        geometry.face_area_vector_m2(face, mesh::FaceSide::owner);
    item.area_m2 = geometry.face_area_m2(face);
    item.skewness = geometry.face_skewness(face);
    item.non_orthogonality_degrees =
        geometry.face_non_orthogonality_degrees(face);
    result.faces.push_back(std::move(item));
  }
  std::sort(result.faces.begin(), result.faces.end(),
            [](const auto& left, const auto& right) {
              return left.global_id < right.global_id;
            });
  for (const auto& [id, logical] : vertices)
    result.vertices.push_back(
        {id, logical, geometry.vertex_position_m(logical)});
  validate(result);
  return result;
}

std::vector<std::byte> encode_mesh_diagnostic_v2(
    const MeshDiagnosticV2& value) {
  validate(value);
  std::vector<std::byte> out;
  out.reserve(128U + value.vertices.size() * 44U +
              value.cells.size() * 72U + value.faces.size() * 180U);
  for (char byte : kMagic)
    out.push_back(static_cast<std::byte>(byte));
  append_unsigned(out, kVersion);
  append_unsigned(out, kEndian);
  append_i32(out, value.rank);
  append_i32(out, value.rank_count);
  for (int component : {value.global_extent.x, value.global_extent.y,
                        value.global_extent.z, value.owned_box.begin.x,
                        value.owned_box.begin.y, value.owned_box.begin.z,
                        value.owned_box.end.x, value.owned_box.end.y,
                        value.owned_box.end.z})
    append_i32(out, component);
  append_unsigned(out, static_cast<std::uint8_t>(value.mapping_kind));
  for (double component :
       {value.origin_m.x, value.origin_m.y, value.origin_m.z,
        value.length_m.x, value.length_m.y, value.length_m.z})
    append_f64(out, component);
  append_unsigned(out, static_cast<std::uint64_t>(value.vertices.size()));
  for (const auto& vertex : value.vertices) {
    append_unsigned(out, vertex.global_id);
    append_i32(out, vertex.logical.x);
    append_i32(out, vertex.logical.y);
    append_i32(out, vertex.logical.z);
    append_f64(out, vertex.position_m.x);
    append_f64(out, vertex.position_m.y);
    append_f64(out, vertex.position_m.z);
  }
  append_unsigned(out, static_cast<std::uint64_t>(value.cells.size()));
  for (const auto& cell : value.cells) {
    append_unsigned(out, cell.global_id);
    for (double component :
         {cell.centre_m.x, cell.centre_m.y, cell.centre_m.z, cell.volume_m3,
          cell.minimum_jacobian_m3, cell.closure_m2.x, cell.closure_m2.y,
          cell.closure_m2.z})
      append_f64(out, component);
  }
  append_unsigned(out, static_cast<std::uint64_t>(value.faces.size()));
  for (const auto& face : value.faces) {
    append_unsigned(out, face.global_id);
    append_unsigned(out, static_cast<std::uint8_t>(face.axis));
    append_i32(out, face.logical.x);
    append_i32(out, face.logical.y);
    append_i32(out, face.logical.z);
    append_unsigned(out, face.owner_global_cell);
    const auto append_optional = [&](const auto& optional,
                                     const auto& append_value) {
      append_unsigned(out,
                      static_cast<std::uint8_t>(optional ? 1U : 0U));
      if (optional)
        append_value(*optional);
    };
    append_optional(face.neighbour_global_cell,
                    [&](auto v) { append_unsigned(out, v); });
    append_optional(face.patch_id,
                    [&](auto v) { append_unsigned(out, v); });
    append_optional(face.periodic_pair,
                    [&](auto v) { append_unsigned(out, v); });
    for (auto vertex : face.vertex_ids)
      append_unsigned(out, vertex);
    for (double component :
         {face.centre_m.x, face.centre_m.y, face.centre_m.z,
          face.owner_area_vector_m2.x, face.owner_area_vector_m2.y,
          face.owner_area_vector_m2.z})
      append_f64(out, component);
    append_optional(face.neighbour_area_vector_m2, [&](runtime::Real3 v) {
      append_f64(out, v.x);
      append_f64(out, v.y);
      append_f64(out, v.z);
    });
    append_f64(out, face.area_m2);
    append_f64(out, face.skewness);
    append_f64(out, face.non_orthogonality_degrees);
  }
  append_f64(out, value.owned_volume_sum);
  append_f64(out, value.maximum_cell_closure_norm);
  append_unsigned(out, value.reciprocal_check_count);
  append_unsigned(out, value.reciprocal_failure_count);
  append_unsigned(out, crc64(out.data(), out.size()));
  return out;
}

MeshDiagnosticV2 decode_mesh_diagnostic_v2(
    const std::vector<std::byte>& bytes) {
  if (bytes.size() < 8U + 8U)
    fail("truncated input");
  const auto payload_size = bytes.size() - sizeof(std::uint64_t);
  std::uint64_t stored_crc{};
  for (std::size_t byte = 0; byte < sizeof(stored_crc); ++byte)
    stored_crc |=
        static_cast<std::uint64_t>(
            static_cast<unsigned char>(bytes[payload_size + byte]))
        << (8U * static_cast<unsigned>(byte));
  if (crc64(bytes.data(), payload_size) != stored_crc)
    fail("CRC-64/ECMA-182 mismatch");

  Reader reader(bytes);
  for (char expected : kMagic)
    if (reader.u8() != static_cast<std::uint8_t>(expected))
      fail("magic mismatch");
  if (reader.u32() != kVersion || reader.u32() != kEndian)
    fail("version or endian marker mismatch");
  MeshDiagnosticV2 result;
  result.rank = reader.i32();
  result.rank_count = reader.i32();
  result.global_extent = reader.int3();
  result.owned_box = {reader.int3(), reader.int3()};
  const auto mapping = reader.u8();
  if (mapping > 1U)
    fail("mapping kind is invalid");
  result.mapping_kind = static_cast<mesh::MappingKind>(mapping);
  result.origin_m = reader.real3();
  result.length_m = reader.real3();
  const auto vertex_count = reader.count(44U);
  result.vertices.reserve(vertex_count);
  for (std::size_t index = 0; index < vertex_count; ++index)
    result.vertices.push_back(
        {reader.u64(), reader.int3(), reader.real3()});
  const auto cell_count = reader.count(72U);
  result.cells.reserve(cell_count);
  for (std::size_t index = 0; index < cell_count; ++index)
    result.cells.push_back({reader.u64(), reader.real3(), reader.f64(),
                            reader.f64(), reader.real3()});
  const auto face_count = reader.count(117U);
  result.faces.reserve(face_count);
  for (std::size_t index = 0; index < face_count; ++index) {
    MeshDiagnosticFaceV2 face;
    face.global_id = reader.u64();
    const auto axis = reader.u8();
    if (axis > 2U)
      fail("face axis is invalid");
    face.axis = static_cast<mesh::FaceAxis>(axis);
    face.logical = reader.int3();
    face.owner_global_cell = reader.u64();
    face.neighbour_global_cell =
        reader.optional<std::uint64_t>([&] { return reader.u64(); });
    face.patch_id =
        reader.optional<std::uint32_t>([&] { return reader.u32(); });
    face.periodic_pair =
        reader.optional<std::uint64_t>([&] { return reader.u64(); });
    for (auto& vertex : face.vertex_ids)
      vertex = reader.u64();
    face.centre_m = reader.real3();
    face.owner_area_vector_m2 = reader.real3();
    face.neighbour_area_vector_m2 =
        reader.optional<runtime::Real3>([&] { return reader.real3(); });
    face.area_m2 = reader.f64();
    face.skewness = reader.f64();
    face.non_orthogonality_degrees = reader.f64();
    result.faces.push_back(std::move(face));
  }
  result.owned_volume_sum = reader.f64();
  result.maximum_cell_closure_norm = reader.f64();
  result.reciprocal_check_count = reader.u64();
  result.reciprocal_failure_count = reader.u64();
  if (reader.position() != payload_size)
    fail("trailing bytes or malformed record size");
  validate(result);
  return result;
}

void write_mesh_diagnostic_v2_file(const std::filesystem::path& path,
                                   const MeshDiagnosticV2& value) {
  write_bytes(path, encode_mesh_diagnostic_v2(value));
}

MeshDiagnosticV2 read_mesh_diagnostic_v2_file(
    const std::filesystem::path& path) {
  std::ifstream stream(path, std::ios::binary);
  if (!stream)
    fail("unable to open input");
  std::vector<char> raw((std::istreambuf_iterator<char>(stream)),
                        std::istreambuf_iterator<char>());
  if (stream.bad())
    fail("unable to read input");
  std::vector<std::byte> bytes(raw.size());
  std::transform(raw.begin(), raw.end(), bytes.begin(),
                 [](char value) { return static_cast<std::byte>(value); });
  return decode_mesh_diagnostic_v2(bytes);
}

void write_mesh_diagnostic_v2(
    const runtime::MpiContext& mpi, const std::filesystem::path& directory,
    const MeshDiagnosticV2& value) {
  bool local_ok = value.rank == mpi.rank() &&
                  value.rank_count == mpi.size();
  std::string message = local_ok ? "" : "meshdiag rank metadata mismatch";
  auto status = runtime::collective_status(mpi, local_ok, message);
  if (!status.ok)
    fail(status.message);
  require_directory_agreement(mpi, directory);
  const auto path = mesh_path(directory, mpi.rank());
  const auto temporary = std::filesystem::path(path.string() + ".tmp");
  try {
    const auto bytes = encode_mesh_diagnostic_v2(value);
    std::filesystem::create_directories(directory);
    std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
    if (!stream)
      fail("unable to open temporary output");
    stream.write(reinterpret_cast<const char*>(bytes.data()),
                 static_cast<std::streamsize>(bytes.size()));
    stream.flush();
    if (!stream)
      fail("unable to stage output");
    stream.close();
    local_ok = static_cast<bool>(stream);
    message = local_ok ? "" : "unable to close staged output";
  } catch (const std::exception& error) {
    local_ok = false;
    message = error.what();
  }
  status = runtime::collective_status(mpi, local_ok, message);
  if (!status.ok) {
    std::error_code ignored;
    std::filesystem::remove(temporary, ignored);
    fail(status.message);
  }
  std::error_code error;
  std::filesystem::rename(temporary, path, error);
  status = runtime::collective_status(
      mpi, !error, error ? "unable to publish meshdiag output" : "");
  if (!status.ok) {
    std::error_code ignored;
    std::filesystem::remove(temporary, ignored);
    fail(status.message);
  }
}

}  // namespace hundun::diagnostics
