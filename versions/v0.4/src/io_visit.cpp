// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_io.hpp"

#include "io_output_detail.hpp"

#include <mpi.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <locale>
#include <new>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace hundun::v04 {
namespace {

std::string padded(std::uint64_t value, std::size_t width) {
  std::string number = std::to_string(value);
  return std::string(width - std::min(width, number.size()), '0') + number;
}

std::string xml_escape(std::string_view input) {
  std::string output;
  for (const char character : input) {
    switch (character) {
      case '&':
        output += "&amp;";
        break;
      case '<':
        output += "&lt;";
        break;
      case '>':
        output += "&gt;";
        break;
      case '"':
        output += "&quot;";
        break;
      case '\'':
        output += "&apos;";
        break;
      default:
        if (static_cast<unsigned char>(character) < 0x20U)
          throw std::invalid_argument("invalid XML field name");
        output += character;
    }
  }
  return output;
}

void append_u64(std::uint8_t*& cursor, std::uint64_t value) noexcept {
  for (unsigned shift = 0U; shift < 64U; shift += 8U)
    *cursor++ = static_cast<std::uint8_t>(value >> shift);
}

void append_real(std::uint8_t*& cursor, double value) noexcept {
  std::uint64_t bits = 0U;
  std::memcpy(&bits, &value, sizeof(bits));
  append_u64(cursor, bits);
}

std::size_t cell_count(Int3 cells) noexcept {
  return static_cast<std::size_t>(cells.x) *
         static_cast<std::size_t>(cells.y) *
         static_cast<std::size_t>(cells.z);
}

Status build_visit_file(const CommittedOutputSnapshot& snapshot,
                        std::size_t capacity,
                        std::vector<std::uint8_t>& out,
                        std::string& extension) {
  const bool uniform = snapshot.geometry->kind() == GeometryKind::uniform;
  extension = uniform ? ".vti" : ".vtr";
  std::size_t appended_bytes = 0U;
  std::vector<std::uint64_t> field_offsets(snapshot.fields.size);
  const std::size_t cells = cell_count(snapshot.patch.cells);
  const auto add_array = [&](std::size_t values, std::uint64_t& offset) {
    if (appended_bytes > capacity ||
        capacity - appended_bytes < sizeof(double) ||
        values > (capacity - appended_bytes - sizeof(double)) / sizeof(double))
      return false;
    offset = appended_bytes;
    appended_bytes += sizeof(std::uint64_t) + values * sizeof(double);
    return true;
  };
  for (std::size_t field_index = 0U; field_index < snapshot.fields.size;
       ++field_index) {
    const ConstFieldView field = snapshot.fields.data[field_index].values;
    // validate_output_snapshot has already checked cells * components.
    if (!add_array(cells * field.components, field_offsets[field_index]))
      return {StatusCode::invalid_plan, detail::kOutputCapacity};
  }
  std::array<std::uint64_t, 3U> coordinate_offsets{};
  if (!uniform) {
    const std::int32_t counts[3]{snapshot.patch.cells.x,
                                 snapshot.patch.cells.y,
                                 snapshot.patch.cells.z};
    for (std::size_t axis = 0U; axis < 3U; ++axis) {
      if (!add_array(static_cast<std::size_t>(counts[axis]) + 1U,
                     coordinate_offsets[axis]))
        return {StatusCode::invalid_plan, detail::kOutputCapacity};
    }
  }

  std::string prefix;
  {
    const Int3 begin = snapshot.patch.begin;
    const Int3 end{begin.x + snapshot.patch.cells.x,
                   begin.y + snapshot.patch.cells.y,
                   begin.z + snapshot.patch.cells.z};
    const Int3 global = snapshot.geometry->global_cells();
    std::ostringstream header;
    header.exceptions(std::ios::badbit | std::ios::failbit);
    header.imbue(std::locale::classic());
    header << std::setprecision(17);
    header << "<?xml version=\"1.0\"?>\n"
           << "<VTKFile type=\"" << (uniform ? "ImageData" : "RectilinearGrid")
           << "\" version=\"1.0\" byte_order=\"LittleEndian\" "
              "header_type=\"UInt64\">\n";
    if (uniform) {
      const AxisMetrics& x = snapshot.geometry->x();
      const AxisMetrics& y = snapshot.geometry->y();
      const AxisMetrics& z = snapshot.geometry->z();
      header << "<ImageData WholeExtent=\"0 " << global.x << " 0 " << global.y
             << " 0 " << global.z << "\" Origin=\""
             << snapshot.geometry->lower().x << ' '
             << snapshot.geometry->lower().y << ' '
             << snapshot.geometry->lower().z << "\" Spacing=\""
             << x.uniform_width() << ' ' << y.uniform_width() << ' '
             << z.uniform_width() << "\">\n";
    } else {
      header << "<RectilinearGrid WholeExtent=\"0 " << global.x << " 0 "
             << global.y << " 0 " << global.z << "\">\n";
    }
    header << "<Piece Extent=\"" << begin.x << ' ' << end.x << ' ' << begin.y
           << ' ' << end.y << ' ' << begin.z << ' ' << end.z << "\">\n"
           << "<CellData>\n";
    for (std::size_t index = 0U; index < snapshot.fields.size; ++index) {
      header << "<DataArray type=\"Float64\" Name=\""
             << xml_escape(snapshot.fields.data[index].stable_name)
             << "\" NumberOfComponents=\""
             << static_cast<unsigned>(
                    snapshot.fields.data[index].values.components)
             << "\" format=\"appended\" offset=\"" << field_offsets[index]
             << "\"/>\n";
    }
    header << "</CellData>\n";
    if (!uniform) {
      header << "<Coordinates>\n";
      for (std::size_t axis = 0U; axis < 3U; ++axis)
        header << "<DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                  "format=\"appended\" offset=\""
               << coordinate_offsets[axis] << "\"/>\n";
      header << "</Coordinates>\n";
    }
    header << "</Piece>\n"
           << (uniform ? "</ImageData>\n" : "</RectilinearGrid>\n")
           << "<AppendedData encoding=\"raw\">_";
    prefix = header.str();
  }  // Release the formatting stream before allocating the binary payload.
  constexpr std::string_view suffix{"</AppendedData>\n</VTKFile>\n"};
  if (prefix.size() > capacity || appended_bytes > capacity - prefix.size() ||
      suffix.size() > capacity - prefix.size() - appended_bytes)
    return {StatusCode::invalid_plan, detail::kOutputCapacity};
  const std::size_t total = prefix.size() + appended_bytes + suffix.size();
  // Account for the metadata still alive beside the single file buffer.
  if (prefix.capacity() > capacity - total ||
      field_offsets.capacity() >
          (capacity - total - prefix.capacity()) / sizeof(std::uint64_t))
    return {StatusCode::invalid_plan, detail::kOutputCapacity};
  out.resize(total);
  std::memcpy(out.data(), prefix.data(), prefix.size());
  std::uint8_t* cursor = out.data() + prefix.size();
  for (std::size_t index = 0U; index < snapshot.fields.size; ++index) {
    const ConstFieldView field = snapshot.fields.data[index].values;
    append_u64(cursor, cells * field.components * sizeof(double));
    for (std::int32_t z = 0; z < field.interior.z; ++z)
      for (std::int32_t y = 0; y < field.interior.y; ++y)
        for (std::int32_t x = 0; x < field.interior.x; ++x)
          for (std::uint8_t component = 0U; component < field.components;
               ++component) {
            const double value = field.unchecked({x, y, z}, component);
            if (!std::isfinite(value))
              return {StatusCode::numerical_failure, detail::kOutputInput};
            append_real(cursor, value);
          }
  }
  if (!uniform) {
    const AxisMetrics* axes[3]{&snapshot.geometry->x(), &snapshot.geometry->y(),
                               &snapshot.geometry->z()};
    const std::int32_t begins[3]{snapshot.patch.begin.x, snapshot.patch.begin.y,
                                 snapshot.patch.begin.z};
    const std::int32_t counts[3]{snapshot.patch.cells.x, snapshot.patch.cells.y,
                                 snapshot.patch.cells.z};
    for (std::size_t axis = 0U; axis < 3U; ++axis) {
      const std::size_t count = static_cast<std::size_t>(counts[axis]) + 1U;
      append_u64(cursor, count * sizeof(double));
      const Span<const double> faces = axes[axis]->faces();
      for (std::size_t index = 0U; index < count; ++index)
        append_real(cursor,
                    faces.data[static_cast<std::size_t>(begins[axis]) + index]);
    }
  }
  std::memcpy(cursor, suffix.data(), suffix.size());
  return {};
}

}  // namespace

Status VisitWriter::write(MPI_Comm communicator,
                          const std::filesystem::path& visit_directory,
                          const IoServicePlan& services,
                          const CommittedOutputSnapshot& snapshot) noexcept {
  int rank = 0;
  int size = 0;
  if (communicator == MPI_COMM_NULL ||
      MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0 ||
      rank < 0)
    return {StatusCode::invalid_plan, detail::kOutputInput};
  Status status = detail::validate_output_snapshot(
      services, RuntimeServiceKind::visit, snapshot);
  if (visit_directory.empty())
    status = {StatusCode::invalid_plan, detail::kOutputInput};
  status = detail::output_collective_status(communicator, status);
  if (!status) return status;
  status = detail::output_create_directory(communicator, rank, visit_directory);
  if (!status) return status;
  const RuntimeServiceCapacity* capacity =
      detail::output_service(services, RuntimeServiceKind::visit);
  std::vector<std::uint8_t> file;
  std::string extension;
  status = detail::output_collective_stage(communicator, [&]() -> Status {
    const Status built = build_visit_file(
        snapshot, capacity->maximum_staging_bytes_per_rank, file, extension);
    if (!built) return built;
    const std::string base = "step-" + padded(snapshot.step, 20U) + "-rank-" +
                             padded(static_cast<std::uint64_t>(rank), 8U) +
                             extension;
    return detail::output_write_file(visit_directory / base, file)
               ? Status{}
               : Status{StatusCode::io_failure, detail::kOutputFile};
  });
  if (!status) return status;
  std::vector<std::uint8_t>{}.swap(file);
  return detail::output_collective_stage(communicator, [&]() -> Status {
    if (rank == 0) {
      std::string index = "!NBLOCKS " + std::to_string(size) + "\n";
      for (int source = 0; source < size; ++source)
        index += "step-" + padded(snapshot.step, 20U) + "-rank-" +
                 padded(static_cast<std::uint64_t>(source), 8U) + extension +
                 "\n";
      if (!detail::output_write_file(
              visit_directory /
                  ("step-" + padded(snapshot.step, 20U) + ".visit"),
              index) ||
          !detail::output_sync_directory(visit_directory))
        return {StatusCode::io_failure, detail::kOutputFile};
    }
    return {};
  });
}

}  // namespace hundun::v04
