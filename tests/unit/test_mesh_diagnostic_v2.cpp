// SPDX-License-Identifier: Apache-2.0

#include "hundun/diagnostics/mesh_diagnostic_v2.hpp"
#include "hundun/runtime/error.hpp"
#include "tests/support/test_main.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <optional>
#include <string_view>
#include <vector>

namespace {

std::uint64_t bits(double value) {
  std::uint64_t result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

bool exact(hundun::runtime::Int3 left, hundun::runtime::Int3 right) {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool exact(hundun::runtime::Real3 left, hundun::runtime::Real3 right) {
  return bits(left.x) == bits(right.x) &&
         bits(left.y) == bits(right.y) &&
         bits(left.z) == bits(right.z);
}

bool exact(hundun::runtime::Box3 left, hundun::runtime::Box3 right) {
  return exact(left.begin, right.begin) && exact(left.end, right.end);
}

bool exact(const std::optional<hundun::runtime::Real3>& left,
           const std::optional<hundun::runtime::Real3>& right) {
  return left.has_value() == right.has_value() &&
         (!left.has_value() || exact(*left, *right));
}

bool exact(const hundun::diagnostics::MeshDiagnosticV2& left,
           const hundun::diagnostics::MeshDiagnosticV2& right) {
  if (left.rank != right.rank || left.rank_count != right.rank_count ||
      !exact(left.global_extent, right.global_extent) ||
      !exact(left.owned_box, right.owned_box) ||
      left.mapping_kind != right.mapping_kind ||
      !exact(left.origin_m, right.origin_m) ||
      !exact(left.length_m, right.length_m) ||
      left.vertices.size() != right.vertices.size() ||
      left.cells.size() != right.cells.size() ||
      left.faces.size() != right.faces.size() ||
      bits(left.owned_volume_sum) != bits(right.owned_volume_sum) ||
      bits(left.maximum_cell_closure_norm) !=
          bits(right.maximum_cell_closure_norm) ||
      left.reciprocal_check_count != right.reciprocal_check_count ||
      left.reciprocal_failure_count != right.reciprocal_failure_count)
    return false;
  for (std::size_t i = 0; i < left.vertices.size(); ++i)
    if (left.vertices[i].global_id != right.vertices[i].global_id ||
        !exact(left.vertices[i].logical, right.vertices[i].logical) ||
        !exact(left.vertices[i].position_m, right.vertices[i].position_m))
      return false;
  for (std::size_t i = 0; i < left.cells.size(); ++i)
    if (left.cells[i].global_id != right.cells[i].global_id ||
        !exact(left.cells[i].centre_m, right.cells[i].centre_m) ||
        bits(left.cells[i].volume_m3) != bits(right.cells[i].volume_m3) ||
        bits(left.cells[i].minimum_jacobian_m3) !=
            bits(right.cells[i].minimum_jacobian_m3) ||
        !exact(left.cells[i].closure_m2, right.cells[i].closure_m2))
      return false;
  for (std::size_t i = 0; i < left.faces.size(); ++i) {
    const auto& a = left.faces[i];
    const auto& b = right.faces[i];
    if (a.global_id != b.global_id || a.axis != b.axis ||
        !exact(a.logical, b.logical) ||
        a.owner_global_cell != b.owner_global_cell ||
        a.neighbour_global_cell != b.neighbour_global_cell ||
        a.patch_id != b.patch_id || a.periodic_pair != b.periodic_pair ||
        a.vertex_ids != b.vertex_ids || !exact(a.centre_m, b.centre_m) ||
        !exact(a.owner_area_vector_m2, b.owner_area_vector_m2) ||
        !exact(a.neighbour_area_vector_m2, b.neighbour_area_vector_m2) ||
        bits(a.area_m2) != bits(b.area_m2) ||
        bits(a.skewness) != bits(b.skewness) ||
        bits(a.non_orthogonality_degrees) !=
            bits(b.non_orthogonality_degrees))
      return false;
  }
  return true;
}

hundun::diagnostics::MeshDiagnosticV2 fixture() {
  using namespace hundun;
  diagnostics::MeshDiagnosticV2 value;
  value.rank = 0;
  value.rank_count = 1;
  value.global_extent = {1, 1, 1};
  value.owned_box = {{0, 0, 0}, {1, 1, 1}};
  value.mapping_kind = mesh::MappingKind::uniform_box;
  value.origin_m = {0.0, 0.0, 0.0};
  value.length_m = {1.0, 1.0, 1.0};
  value.vertices = {
      {0, {0, 0, 0}, {0.0, 0.0, 0.0}},
      {1, {1, 0, 0}, {1.0, 0.0, 0.0}},
      {2, {0, 1, 0}, {0.0, 1.0, 0.0}},
      {3, {1, 1, 0}, {1.0, 1.0, 0.0}},
      {4, {0, 0, 1}, {0.0, 0.0, 1.0}},
      {5, {1, 0, 1}, {1.0, 0.0, 1.0}},
      {6, {0, 1, 1}, {0.0, 1.0, 1.0}},
      {7, {1, 1, 1}, {1.0, 1.0, 1.0}},
  };
  value.cells = {{0, {0.5, 0.5, 0.5}, 1.0, 1.0, {0.0, 0.0, 0.0}}};
  value.faces = {
      {0,
       mesh::FaceAxis::x,
       {0, 0, 0},
       0,
       std::nullopt,
       0U,
       std::nullopt,
       {0, 2, 6, 4},
       {0.0, 0.5, 0.5},
       {-1.0, 0.0, 0.0},
       std::nullopt,
       1.0,
       0.0,
       0.0},
      {1,
       mesh::FaceAxis::x,
       {1, 0, 0},
       0,
       std::nullopt,
       1U,
       std::nullopt,
       {1, 3, 7, 5},
       {1.0, 0.5, 0.5},
       {1.0, 0.0, 0.0},
       std::nullopt,
       1.0,
       0.0,
       0.0},
      {2,
       mesh::FaceAxis::y,
       {0, 0, 0},
       0,
       std::nullopt,
       2U,
       std::nullopt,
       {0, 4, 5, 1},
       {0.5, 0.0, 0.5},
       {0.0, -1.0, 0.0},
       std::nullopt,
       1.0,
       0.0,
       0.0},
      {3,
       mesh::FaceAxis::y,
       {0, 1, 0},
       0,
       std::nullopt,
       3U,
       std::nullopt,
       {2, 6, 7, 3},
       {0.5, 1.0, 0.5},
       {0.0, 1.0, 0.0},
       std::nullopt,
       1.0,
       0.0,
       0.0},
      {4,
       mesh::FaceAxis::z,
       {0, 0, 0},
       0,
       std::nullopt,
       4U,
       std::nullopt,
       {0, 1, 3, 2},
       {0.5, 0.5, 0.0},
       {0.0, 0.0, -1.0},
       std::nullopt,
       1.0,
       0.0,
       0.0},
      {5,
       mesh::FaceAxis::z,
       {0, 0, 1},
       0,
       std::nullopt,
       5U,
       std::nullopt,
       {4, 5, 7, 6},
       {0.5, 0.5, 1.0},
       {0.0, 0.0, 1.0},
       std::nullopt,
       1.0,
       0.0,
       0.0},
  };
  value.owned_volume_sum = 1.0;
  value.maximum_cell_closure_norm = 0.0;
  return value;
}

hundun::diagnostics::MeshDiagnosticV2 rescaled_three_cell_fixture() {
  using namespace hundun;
  diagnostics::MeshDiagnosticV2 value;
  value.rank = 0;
  value.rank_count = 1;
  value.global_extent = {3, 1, 1};
  value.owned_box = {{0, 0, 0}, {3, 1, 1}};
  value.mapping_kind = mesh::MappingKind::uniform_box;
  value.origin_m = {0.0, 0.0, 0.0};
  value.length_m = {3.0, 1.0, 1.0};
  for (int z = 0; z <= 1; ++z)
    for (int y = 0; y <= 1; ++y)
      for (int x = 0; x <= 3; ++x) {
        const auto id = static_cast<std::uint64_t>((z * 2 + y) * 4 + x);
        value.vertices.push_back(
            {id,
             {x, y, z},
             {static_cast<double>(x), static_cast<double>(y),
              static_cast<double>(z)}});
      }
  for (int x = 0; x < 3; ++x)
    value.cells.push_back({static_cast<std::uint64_t>(x),
                           {static_cast<double>(x) + 0.5, 0.5, 0.5},
                           1.0,
                           1.0,
                           {0.0, 0.0, 0.0}});
  value.owned_volume_sum = 3.0;

  const auto vertex_id = [](int x, int y, int z) {
    return static_cast<std::uint64_t>((z * 2 + y) * 4 + x);
  };
  const auto add_face = [&](std::uint64_t id, mesh::FaceAxis axis,
                            runtime::Int3 logical, std::uint64_t owner,
                            std::optional<std::uint64_t> neighbour,
                            std::optional<std::uint32_t> patch, double area) {
    diagnostics::MeshDiagnosticFaceV2 face;
    face.global_id = id;
    face.axis = axis;
    face.logical = logical;
    face.owner_global_cell = owner;
    face.neighbour_global_cell = neighbour;
    face.patch_id = patch;
    if (axis == mesh::FaceAxis::x) {
      face.vertex_ids = {vertex_id(logical.x, 0, 0), vertex_id(logical.x, 1, 0),
                         vertex_id(logical.x, 1, 1),
                         vertex_id(logical.x, 0, 1)};
      face.centre_m = {static_cast<double>(logical.x), 0.5, 0.5};
      face.owner_area_vector_m2 = {logical.x == 0 ? -area : area, 0.0, 0.0};
    } else if (axis == mesh::FaceAxis::y) {
      face.vertex_ids = {vertex_id(logical.x, logical.y, 0),
                         vertex_id(logical.x, logical.y, 1),
                         vertex_id(logical.x + 1, logical.y, 1),
                         vertex_id(logical.x + 1, logical.y, 0)};
      face.centre_m = {static_cast<double>(logical.x) + 0.5,
                       static_cast<double>(logical.y), 0.5};
      face.owner_area_vector_m2 = {0.0, logical.y == 0 ? -area : area, 0.0};
    } else {
      face.vertex_ids = {vertex_id(logical.x, 0, logical.z),
                         vertex_id(logical.x + 1, 0, logical.z),
                         vertex_id(logical.x + 1, 1, logical.z),
                         vertex_id(logical.x, 1, logical.z)};
      face.centre_m = {static_cast<double>(logical.x) + 0.5, 0.5,
                       static_cast<double>(logical.z)};
      face.owner_area_vector_m2 = {0.0, 0.0, logical.z == 0 ? -area : area};
    }
    if (neighbour)
      face.neighbour_area_vector_m2 = runtime::Real3{
          -face.owner_area_vector_m2.x, -face.owner_area_vector_m2.y,
          -face.owner_area_vector_m2.z};
    face.area_m2 = area;
    value.faces.push_back(std::move(face));
  };

  constexpr double small = 1.0e-6;
  add_face(0U, mesh::FaceAxis::x, {0, 0, 0}, 0U, std::nullopt, 0U, small);
  add_face(1U, mesh::FaceAxis::x, {1, 0, 0}, 0U, 1U, std::nullopt, small);
  add_face(2U, mesh::FaceAxis::x, {2, 0, 0}, 1U, 2U, std::nullopt, 1.0);
  add_face(3U, mesh::FaceAxis::x, {3, 0, 0}, 2U, std::nullopt, 1U, 1.0e6);
  for (int y = 0; y <= 1; ++y)
    for (int x = 0; x < 3; ++x)
      add_face(static_cast<std::uint64_t>(4 + y * 3 + x), mesh::FaceAxis::y,
               {x, y, 0}, static_cast<std::uint64_t>(x), std::nullopt,
               static_cast<std::uint32_t>(2 + y),
               x == 0 ? small : (x == 1 ? 1.0 : 1.0e6));
  for (int z = 0; z <= 1; ++z)
    for (int x = 0; x < 3; ++x)
      add_face(static_cast<std::uint64_t>(10 + z * 3 + x), mesh::FaceAxis::z,
               {x, 0, z}, static_cast<std::uint64_t>(x), std::nullopt,
               static_cast<std::uint32_t>(4 + z),
               x == 0 ? small : (x == 1 ? 1.0 : 1.0e6));
  value.reciprocal_check_count = 2U;
  return value;
}

hundun::diagnostics::MeshDiagnosticV2 partition_middle_cell_fixture() {
  auto value = rescaled_three_cell_fixture();
  value.rank = 1;
  value.rank_count = 3;
  value.owned_box = {{1, 0, 0}, {2, 1, 1}};
  value.cells = {value.cells[1]};
  value.owned_volume_sum = value.cells.front().volume_m3;
  value.maximum_cell_closure_norm = 0.0;
  value.faces.erase(
      std::remove_if(value.faces.begin(), value.faces.end(),
                     [](const auto& face) {
                       return face.global_id != 1U &&
                              face.global_id != 2U &&
                              face.global_id != 5U &&
                              face.global_id != 8U &&
                              face.global_id != 11U &&
                              face.global_id != 14U;
                     }),
      value.faces.end());
  value.vertices.erase(
      std::remove_if(value.vertices.begin(), value.vertices.end(),
                     [](const auto& vertex) {
                       return vertex.logical.x < 1 || vertex.logical.x > 2;
                     }),
      value.vertices.end());
  value.reciprocal_check_count = 2U;
  value.reciprocal_failure_count = 0U;
  return value;
}

template <class Function> bool rejects(Function &&operation) {
  try {
    operation();
  } catch (const hundun::runtime::Error&) {
    return true;
  }
  return false;
}

std::uint64_t crc64(const std::vector<std::byte>& bytes,
                    std::size_t size) noexcept {
  std::uint64_t crc = 0U;
  constexpr std::uint64_t polynomial = UINT64_C(0x42f0e1eba9ea3693);
  for (std::size_t index = 0; index < size; ++index) {
    crc ^= static_cast<std::uint64_t>(
               static_cast<unsigned char>(bytes[index]))
           << 56U;
    for (int bit = 0; bit < 8; ++bit)
      crc = (crc & UINT64_C(0x8000000000000000)) != 0U
                ? (crc << 1U) ^ polynomial
                : crc << 1U;
  }
  return crc;
}

void refresh_crc(std::vector<std::byte>& bytes) {
  const auto payload = bytes.size() - sizeof(std::uint64_t);
  const auto crc = crc64(bytes, payload);
  for (std::size_t byte = 0; byte < sizeof(crc); ++byte)
    bytes[payload + byte] =
        static_cast<std::byte>((crc >> (8U * byte)) & UINT64_C(0xff));
}

}  // namespace

int main() {
  return hundun::test::run([] {
    auto value = fixture();
    const auto bytes = hundun::diagnostics::encode_mesh_diagnostic_v2(value);
    HUNDUN_CHECK(bytes.size() > 8U);
    HUNDUN_CHECK(bytes[0] == std::byte{'H'});
    HUNDUN_CHECK(bytes[7] == std::byte{'2'});
    const auto decoded =
        hundun::diagnostics::decode_mesh_diagnostic_v2(bytes);
    HUNDUN_CHECK(exact(value, decoded));
    HUNDUN_CHECK(hundun::diagnostics::encode_mesh_diagnostic_v2(decoded) ==
                 bytes);

    auto changed = decoded;
    changed.cells.front().volume_m3 = 2.0;
    HUNDUN_CHECK(!exact(decoded, changed));
    auto nested_changed = decoded;
    nested_changed.faces.front().vertex_ids[0] =
        nested_changed.faces.front().vertex_ids[1];
    HUNDUN_CHECK(!exact(decoded, nested_changed));

    const std::string_view crc_reference = "123456789";
    std::vector<std::byte> crc_reference_bytes;
    for (const char byte : crc_reference)
      crc_reference_bytes.push_back(
          static_cast<std::byte>(static_cast<unsigned char>(byte)));
    HUNDUN_CHECK(crc64(crc_reference_bytes, crc_reference_bytes.size()) ==
                 UINT64_C(0x6c40df5f0b497347));

    auto corrupted = bytes;
    corrupted.back() ^= std::byte{1};
    HUNDUN_CHECK(rejects([&] {
      static_cast<void>(
          hundun::diagnostics::decode_mesh_diagnostic_v2(corrupted));
    }));
    for (std::size_t size = 0; size < bytes.size(); ++size) {
      std::vector<std::byte> truncated(
          bytes.begin(),
          bytes.begin() + static_cast<std::ptrdiff_t>(size));
      HUNDUN_CHECK(rejects([&] {
        static_cast<void>(
            hundun::diagnostics::decode_mesh_diagnostic_v2(truncated));
      }));
    }
    auto trailing = bytes;
    trailing.push_back(std::byte{0});
    HUNDUN_CHECK(rejects([&] {
      static_cast<void>(
          hundun::diagnostics::decode_mesh_diagnostic_v2(trailing));
    }));

    for (const auto offset : {0U, 8U, 12U, 16U, 565U, 586U, 587U, 592U}) {
      auto structurally_corrupt = bytes;
      structurally_corrupt[offset] ^= std::byte{0x7f};
      refresh_crc(structurally_corrupt);
      HUNDUN_CHECK(rejects([&] {
        static_cast<void>(hundun::diagnostics::decode_mesh_diagnostic_v2(
            structurally_corrupt));
      }));
    }

    auto duplicate_vertex = value;
    duplicate_vertex.vertices[1].global_id =
        duplicate_vertex.vertices[0].global_id;
    HUNDUN_CHECK(rejects([&] {
      static_cast<void>(
          hundun::diagnostics::encode_mesh_diagnostic_v2(duplicate_vertex));
    }));
    auto negative_volume = value;
    negative_volume.cells.front().volume_m3 = -1.0;
    HUNDUN_CHECK(rejects([&] {
      static_cast<void>(
          hundun::diagnostics::encode_mesh_diagnostic_v2(negative_volume));
    }));
    auto wrong_summary = value;
    wrong_summary.owned_volume_sum = 2.0;
    HUNDUN_CHECK(rejects([&] {
      static_cast<void>(
          hundun::diagnostics::encode_mesh_diagnostic_v2(wrong_summary));
    }));
    auto wrong_patch = value;
    wrong_patch.faces.front().patch_id = 1U;
    HUNDUN_CHECK(rejects([&] {
      static_cast<void>(
          hundun::diagnostics::encode_mesh_diagnostic_v2(wrong_patch));
    }));
    auto wrong_owner = value;
    wrong_owner.global_extent = {2, 1, 1};
    wrong_owner.owned_box = {{0, 0, 0}, {1, 1, 1}};
    wrong_owner.faces.front().owner_global_cell = 1U;
    HUNDUN_CHECK(rejects([&] {
      static_cast<void>(
          hundun::diagnostics::encode_mesh_diagnostic_v2(wrong_owner));
    }));
    auto wrong_neighbour = value;
    wrong_neighbour.faces.front().neighbour_global_cell = 0U;
    HUNDUN_CHECK(rejects([&] {
      static_cast<void>(
          hundun::diagnostics::encode_mesh_diagnostic_v2(wrong_neighbour));
    }));
    auto wrong_face_id = value;
    wrong_face_id.faces.front().global_id = 5U;
    HUNDUN_CHECK(rejects([&] {
      static_cast<void>(
          hundun::diagnostics::encode_mesh_diagnostic_v2(wrong_face_id));
    }));
    auto nonzero_reciprocal_failure = value;
    nonzero_reciprocal_failure.reciprocal_failure_count = 1U;
    HUNDUN_CHECK(rejects([&] {
      static_cast<void>(hundun::diagnostics::encode_mesh_diagnostic_v2(
          nonzero_reciprocal_failure));
    }));
    auto invalid_closure = value;
    invalid_closure.cells.front().closure_m2 = {1.0, 0.0, 0.0};
    invalid_closure.maximum_cell_closure_norm = 1.0;
    HUNDUN_CHECK(rejects([&] {
      static_cast<void>(
          hundun::diagnostics::encode_mesh_diagnostic_v2(invalid_closure));
    }));
    auto bad_small_cell = rescaled_three_cell_fixture();
    bad_small_cell.cells.front().closure_m2 = {1.0e-12, 0.0, 0.0};
    bad_small_cell.maximum_cell_closure_norm = 1.0e-12;
    HUNDUN_CHECK(rejects([&] {
      static_cast<void>(
          hundun::diagnostics::encode_mesh_diagnostic_v2(bad_small_cell));
    }));
    auto just_inside = rescaled_three_cell_fixture();
    const double per_cell_limit =
        256.0 * std::numeric_limits<double>::epsilon() * 6.0e-6;
    just_inside.cells.front().closure_m2 = {std::nextafter(per_cell_limit, 0.0),
                                            0.0, 0.0};
    just_inside.maximum_cell_closure_norm =
        just_inside.cells.front().closure_m2.x;
    HUNDUN_CHECK(!rejects([&] {
      static_cast<void>(
          hundun::diagnostics::encode_mesh_diagnostic_v2(just_inside));
    }));
    const auto partition_middle = partition_middle_cell_fixture();
    HUNDUN_CHECK(partition_middle.faces.size() == 6U);
    HUNDUN_CHECK(partition_middle.faces.front().owner_global_cell == 0U);
    HUNDUN_CHECK(partition_middle.faces.front().neighbour_global_cell ==
                 std::optional<std::uint64_t>{1U});
    HUNDUN_CHECK(!rejects([&] {
      const auto bytes =
          hundun::diagnostics::encode_mesh_diagnostic_v2(partition_middle);
      const auto decoded =
          hundun::diagnostics::decode_mesh_diagnostic_v2(bytes);
      HUNDUN_CHECK(decoded.faces.size() == 6U);
      HUNDUN_CHECK(exact(decoded, partition_middle));
    }));
    auto partition_mutation = partition_middle;
    const double complete_area_sum = 1.0e-6 + 5.0;
    const double old_owned_face_area_sum = 5.0;
    const double complete_limit =
        256.0 * std::numeric_limits<double>::epsilon() * complete_area_sum;
    const double old_limit =
        256.0 * std::numeric_limits<double>::epsilon() *
        old_owned_face_area_sum;
    HUNDUN_CHECK(complete_limit > old_limit);
    partition_mutation.cells.front().closure_m2 = {
        old_limit + 0.5 * (complete_limit - old_limit), 0.0, 0.0};
    partition_mutation.maximum_cell_closure_norm =
        partition_mutation.cells.front().closure_m2.x;
    HUNDUN_CHECK(!rejects([&] {
      static_cast<void>(
          hundun::diagnostics::encode_mesh_diagnostic_v2(partition_mutation));
    }));
    auto missing_face = value;
    missing_face.faces.pop_back();
    HUNDUN_CHECK(rejects([&] {
      static_cast<void>(
          hundun::diagnostics::encode_mesh_diagnostic_v2(missing_face));
    }));
    auto extra_face = value;
    extra_face.faces.push_back(extra_face.faces.back());
    extra_face.faces.back().global_id = 6U;
    HUNDUN_CHECK(rejects([&] {
      static_cast<void>(
          hundun::diagnostics::encode_mesh_diagnostic_v2(extra_face));
    }));
    HUNDUN_CHECK(rejects([&] {
      static_cast<void>(hundun::diagnostics::mesh_diagnostic_global_vertex_id(
          {std::numeric_limits<int>::max(),
           std::numeric_limits<int>::max(),
           std::numeric_limits<int>::max()},
          {std::numeric_limits<int>::max(),
           std::numeric_limits<int>::max(),
           std::numeric_limits<int>::max()}));
    }));

    const auto directory = std::filesystem::temp_directory_path() /
                           "hundun-task24-meshdiag-unit";
    std::filesystem::remove_all(directory);
    std::filesystem::create_directories(directory);
    const auto path = directory / "meshdiag.v2.rank-000000.bin";
    hundun::diagnostics::write_mesh_diagnostic_v2_file(path, value);
    HUNDUN_CHECK(exact(
        value, hundun::diagnostics::read_mesh_diagnostic_v2_file(path)));
    std::filesystem::remove_all(directory);
  });
}
