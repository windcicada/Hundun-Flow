// SPDX-License-Identifier: Apache-2.0

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

void append_u64(std::vector<std::uint8_t>& bytes,
                std::uint64_t value) {
  for (unsigned shift = 0U; shift < 64U; shift += 8U)
    bytes.push_back(static_cast<std::uint8_t>(value >> shift));
}

void append_real(std::vector<std::uint8_t>& bytes, double value) {
  std::uint64_t bits = 0U;
  std::memcpy(&bits, &value, sizeof(bits));
  append_u64(bytes, bits);
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
  std::vector<std::uint8_t> appended;
  std::vector<std::uint64_t> field_offsets(snapshot.fields.size);
  const std::size_t cells = cell_count(snapshot.patch.cells);
  for (std::size_t field_index = 0U; field_index < snapshot.fields.size;
       ++field_index) {
    const ConstFieldView field = snapshot.fields.data[field_index].values;
    field_offsets[field_index] = appended.size();
    const std::uint64_t payload =
        static_cast<std::uint64_t>(cells) * field.components * sizeof(double);
    append_u64(appended, payload);
    for (std::int32_t z = 0; z < field.interior.z; ++z)
      for (std::int32_t y = 0; y < field.interior.y; ++y)
        for (std::int32_t x = 0; x < field.interior.x; ++x)
          for (std::uint8_t component = 0U; component < field.components;
               ++component) {
            const double value = field.unchecked({x, y, z}, component);
            if (!std::isfinite(value))
              return {StatusCode::numerical_failure, detail::kOutputInput};
            append_real(appended, value);
          }
  }
  std::array<std::uint64_t, 3U> coordinate_offsets{};
  if (!uniform) {
    const AxisMetrics* axes[3]{&snapshot.geometry->x(),
                               &snapshot.geometry->y(),
                               &snapshot.geometry->z()};
    const std::int32_t begins[3]{snapshot.patch.begin.x,
                                 snapshot.patch.begin.y,
                                 snapshot.patch.begin.z};
    const std::int32_t counts[3]{snapshot.patch.cells.x,
                                 snapshot.patch.cells.y,
                                 snapshot.patch.cells.z};
    for (std::size_t axis = 0U; axis < 3U; ++axis) {
      coordinate_offsets[axis] = appended.size();
      append_u64(appended,
                 static_cast<std::uint64_t>(counts[axis] + 1) *
                     sizeof(double));
      const Span<const double> faces = axes[axis]->faces();
      for (std::int32_t index = 0; index <= counts[axis]; ++index)
        append_real(appended,
                    faces.data[static_cast<std::size_t>(begins[axis] + index)]);
    }
  }

  const Int3 begin = snapshot.patch.begin;
  const Int3 end{begin.x + snapshot.patch.cells.x,
                 begin.y + snapshot.patch.cells.y,
                 begin.z + snapshot.patch.cells.z};
  const Int3 global = snapshot.geometry->global_cells();
  std::ostringstream header;
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
           << static_cast<unsigned>(snapshot.fields.data[index].values.components)
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
  const std::string prefix = header.str();
  constexpr std::string_view suffix{"</AppendedData>\n</VTKFile>\n"};
  if (prefix.size() > capacity || appended.size() > capacity - prefix.size() ||
      suffix.size() > capacity - prefix.size() - appended.size())
    return {StatusCode::invalid_plan, detail::kOutputCapacity};
  out.reserve(prefix.size() + appended.size() + suffix.size());
  out.insert(out.end(), prefix.begin(), prefix.end());
  out.insert(out.end(), appended.begin(), appended.end());
  out.insert(out.end(), suffix.begin(), suffix.end());
  return {};
}

}  // namespace

Status VisitWriter::write(MPI_Comm communicator,
                          const std::filesystem::path& visit_directory,
                          const IoServicePlan& services,
                          const CommittedOutputSnapshot& snapshot) noexcept try {
  int rank = 0;
  int size = 0;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0 ||
      visit_directory.empty())
    return {StatusCode::invalid_plan, detail::kOutputInput};
  Status status = detail::validate_output_snapshot(
      services, RuntimeServiceKind::visit, snapshot);
  status = detail::output_collective_status(communicator, status);
  if (!status) return status;
  status = detail::output_create_directory(communicator, rank, visit_directory);
  if (!status) return status;
  const RuntimeServiceCapacity* capacity =
      detail::output_service(services, RuntimeServiceKind::visit);
  std::vector<std::uint8_t> file;
  std::string extension;
  status = build_visit_file(snapshot, capacity->maximum_staging_bytes_per_rank,
                            file, extension);
  const std::string base = "step-" + padded(snapshot.step, 20U) + "-rank-" +
                           padded(static_cast<std::uint64_t>(rank), 8U) +
                           extension;
  if (status && !detail::output_write_file(visit_directory / base, file))
    status = {StatusCode::io_failure, detail::kOutputFile};
  status = detail::output_collective_status(communicator, status);
  if (!status) return status;
  if (rank == 0) {
    std::string index = "!NBLOCKS " + std::to_string(size) + "\n";
    for (int source = 0; source < size; ++source)
      index += "step-" + padded(snapshot.step, 20U) + "-rank-" +
               padded(static_cast<std::uint64_t>(source), 8U) + extension +
               "\n";
    status = detail::output_write_file(
                 visit_directory /
                     ("step-" + padded(snapshot.step, 20U) + ".visit"),
                 index)
                 ? Status{}
                 : Status{StatusCode::io_failure, detail::kOutputFile};
    if (status && !detail::output_sync_directory(visit_directory))
      status = {StatusCode::io_failure, detail::kOutputFile};
  }
  return detail::output_collective_status(communicator, status);
} catch (const std::bad_alloc&) {
  return {StatusCode::allocation_failure, detail::kOutputCapacity};
} catch (...) {
  return {StatusCode::io_failure, detail::kOutputFile};
}

}  // namespace hundun::v04
