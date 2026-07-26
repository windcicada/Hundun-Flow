// SPDX-License-Identifier: Apache-2.0

#include "hundun/diagnostics/mesh_diagnostic_v2.hpp"
#include "hundun/runtime/error.hpp"
#include "tests/support/test_main.hpp"

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <optional>
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
  };
  value.owned_volume_sum = 1.0;
  value.maximum_cell_closure_norm = 0.0;
  return value;
}

template <class Function>
bool rejects(Function&& operation) {
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

    auto corrupted = bytes;
    corrupted.back() ^= std::byte{1};
    HUNDUN_CHECK(rejects([&] {
      static_cast<void>(
          hundun::diagnostics::decode_mesh_diagnostic_v2(corrupted));
    }));
    auto truncated = bytes;
    truncated.pop_back();
    HUNDUN_CHECK(rejects([&] {
      static_cast<void>(
          hundun::diagnostics::decode_mesh_diagnostic_v2(truncated));
    }));
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
    wrong_owner.faces.front().owner_global_cell = 1U;
    HUNDUN_CHECK(rejects([&] {
      static_cast<void>(
          hundun::diagnostics::encode_mesh_diagnostic_v2(wrong_owner));
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
