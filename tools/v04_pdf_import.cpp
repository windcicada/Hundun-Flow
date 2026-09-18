// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn
// Current-state PDF bridge. Missing temporal histories are explicitly rebuilt.
#include "hundun/v04_app.hpp"
#include "core_tcr_history_detail.hpp"
#include "esf_count_detail.hpp"
#include "hundun/v04_io.hpp"
#include "hundun/v04_mesh.hpp"
#include "hundun/v04_physics.hpp"
#include "yyjson.h"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <set>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

namespace {

using namespace hundun::v04;

enum : std::uint32_t {
  kCollective = 24101,
  kAllocation = 24102,
  kInputContract = 24103
};
bool checked_product(Int3 cells, std::size_t &out) noexcept {
  if (cells.x <= 0 || cells.y <= 0 || cells.z <= 0)
    return false;
  const auto x = static_cast<std::size_t>(cells.x);
  const auto y = static_cast<std::size_t>(cells.y);
  const auto z = static_cast<std::size_t>(cells.z);
  if (x > std::numeric_limits<std::size_t>::max() / y)
    return false;
  const std::size_t xy = x * y;
  if (xy > std::numeric_limits<std::size_t>::max() / z)
    return false;
  out = xy * z;
  return true;
}

bool same_int3(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool same_patch(MeshPatch left, MeshPatch right) noexcept {
  return same_int3(left.begin, right.begin) &&
         same_int3(left.cells, right.cells) &&
         same_int3(left.process_grid, right.process_grid) &&
         same_int3(left.process_coord, right.process_coord);
}

Status consensus(MPI_Comm communicator, Status local) noexcept {
  int rank = -1;
  int size = 0;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0)
    return {StatusCode::mpi_failure, kCollective};
  int failing = local ? size : rank;
  int lowest = size;
  if (MPI_Allreduce(&failing, &lowest, 1, MPI_INT, MPI_MIN, communicator) !=
      MPI_SUCCESS)
    return {StatusCode::mpi_failure, kCollective};
  if (lowest == size)
    return {};
  std::array<std::uint32_t, 2U> payload{};
  if (rank == lowest) {
    payload[0U] = static_cast<std::uint32_t>(local.code);
    payload[1U] = local.detail;
  }
  if (MPI_Bcast(payload.data(), static_cast<int>(payload.size()), MPI_UINT32_T,
                lowest, communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kCollective};
  return {static_cast<StatusCode>(payload[0U]), payload[1U]};
}

template <class Function>
Status local_stage(MPI_Comm communicator, Function &&function) noexcept {
  Status local;
  try {
    local = function();
  } catch (const std::bad_alloc &) {
    local = {StatusCode::allocation_failure, kAllocation};
  } catch (...) {
    local = {StatusCode::invalid_plan, kInputContract};
  }
  return consensus(communicator, local);
}

struct Header {
  unsigned version{};
  Int3 cells{};
  std::uint64_t step{};
  double time{}, dt{}, pressure{};
  std::size_t nf{}, ns{};
  std::vector<std::string> names;
};
Status invalid(unsigned detail = 24104) {
  return {StatusCode::invalid_case, detail};
}
Status read_header(const std::filesystem::path &root, Header &h) {
  std::ifstream in(root / "state.txt");
  std::string magic;
  in >> magic >> h.version >> h.cells.x >> h.cells.y >> h.cells.z >> h.step >>
      h.time >> h.dt >> h.pressure >> h.nf >> h.ns;
  std::size_t count = 0;
  if (!in || magic != "HUNDUN_PDF_TRANSFER" ||
      (h.version != 1 && h.version != 2) ||
      !checked_product(h.cells, count) || !esf::valid_field_count(h.nf) ||
      h.ns < 2 || h.ns > 64 || h.step == 0 || !std::isfinite(h.time) ||
      h.time < 0 || !std::isfinite(h.dt) || h.dt <= 0 ||
      !std::isfinite(h.pressure) || h.pressure <= 0)
    return invalid();
  h.names.resize(h.ns);
  for (auto &name : h.names)
    if (!(in >> name))
      return invalid();
  std::string extra;
  if (in >> extra)
    return invalid();
  return {};
}
struct Block {
  Int3 begin{}, cells{};
  std::vector<double> flow, reference_density, density, mean, cache;
  std::vector<std::vector<double>> pdf;
  std::vector<std::uint8_t> fluid;
  std::size_t offset(Int3 g) const {
    return (g.x - begin.x) +
           std::size_t(cells.x) *
               ((g.y - begin.y) + std::size_t(cells.y) * (g.z - begin.z));
  }
};

struct ImportedParcel {
  std::uint64_t id_high{}, id_low{}, material{}, owner{};
  std::array<double, 3> position{}, velocity{};
  double mass{}, diameter{}, multiplicity{}, temperature{};
};
struct ImportedInjector {
  std::uint64_t id{}, next_ordinal{};
  double residual_mass{};
};
struct SprayTransfer {
  std::vector<ImportedParcel> parcels;
  std::vector<ImportedInjector> injectors;
};

bool json_number(yyjson_val *value, double &out) {
  if (!value || !yyjson_is_num(value)) return false;
  out = yyjson_get_real(value);
  return std::isfinite(out);
}
bool json_u64(yyjson_val *value, std::uint64_t &out) {
  if (!value || !yyjson_is_uint(value)) return false;
  out = yyjson_get_uint(value);
  return true;
}
bool json_text(yyjson_val *value, std::string_view &out) {
  if (!value || !yyjson_is_str(value)) return false;
  out = {yyjson_get_str(value), yyjson_get_len(value)};
  return true;
}
bool json_vector3(yyjson_val *value, std::array<double, 3> &out) {
  if (!value || !yyjson_is_arr(value) || yyjson_arr_size(value) != 3)
    return false;
  for (std::size_t i = 0; i < out.size(); ++i)
    if (!json_number(yyjson_arr_get(value, i), out[i])) return false;
  return true;
}
std::uint64_t hex64(std::string_view text) {
  if (text.size() < 16) return 0;
  std::uint64_t value = 0;
  for (std::size_t i = 0; i < 16; ++i) {
    const char c = text[i];
    const unsigned digit = c >= '0' && c <= '9' ? unsigned(c - '0')
                           : c >= 'a' && c <= 'f' ? unsigned(c - 'a' + 10)
                                                  : 16U;
    if (digit > 15U) return 0;
    value = (value << 4U) | digit;
  }
  return value;
}
Status parse_json(const std::string &text, yyjson_doc *&out) {
  out = yyjson_read(text.data(), text.size(), 0);
  return out ? Status{} : invalid(24114);
}
Status read_text(const std::filesystem::path &path, std::string &out) {
  std::ifstream in(path, std::ios::binary);
  if (!in) return invalid(24114);
  in.seekg(0, std::ios::end);
  const auto n = in.tellg();
  if (n < 0 || n > std::streamoff(64U * 1024U * 1024U))
    return invalid(24114);
  out.resize(static_cast<std::size_t>(n));
  in.seekg(0);
  if (n && !in.read(out.data(), n)) return invalid(24114);
  return {};
}
Status read_spray_transfer(const std::filesystem::path &root,
                           const ValidatedModel &model, SprayTransfer &out) {
  if (!model.spray) return invalid(24114);
  std::string manifest_text;
  auto status = read_text(root / "restart.json", manifest_text);
  yyjson_doc *manifest_doc = nullptr;
  if (status) status = parse_json(manifest_text, manifest_doc);
  if (!status) return status;
  auto *manifest = yyjson_doc_get_root(manifest_doc);
  std::string_view format, file;
  std::uint64_t parcel_count = 0;
  auto *policy = yyjson_obj_get(manifest, "history_policy");
  const auto reset = [&](const char *key) {
    std::string_view value;
    return json_text(policy ? yyjson_obj_get(policy, key) : nullptr, value) &&
           value == "reset_zero";
  };
  bool valid = yyjson_is_obj(manifest) &&
               json_text(yyjson_obj_get(manifest, "format"), format) &&
               format == "hundun_spray_restart_import_v1" &&
               json_text(yyjson_obj_get(manifest, "parcels_file"), file) &&
               file == "parcels.jsonl" &&
               json_u64(yyjson_obj_get(manifest, "parcels"), parcel_count) &&
               parcel_count <= model.spray->maximum_local_parcels * UINT64_C(1048576) &&
               yyjson_is_obj(policy) && reset("age") &&
               reset("tab") && reset("breakup_ordinal") && reset("sgs");
  auto *injectors = yyjson_obj_get(manifest, "injectors");
  valid = valid && yyjson_is_arr(injectors) &&
          yyjson_arr_size(injectors) == model.spray->injectors.size();
  if (valid) {
    out.injectors.reserve(yyjson_arr_size(injectors));
    for (std::size_t i = 0; i < yyjson_arr_size(injectors); ++i) {
      auto *row = yyjson_arr_get(injectors, i);
      ImportedInjector value;
      valid = yyjson_is_obj(row) &&
              json_u64(yyjson_obj_get(row, "id"), value.id) && value.id != 0 &&
              json_u64(yyjson_obj_get(row, "next_ordinal"), value.next_ordinal) &&
              json_number(yyjson_obj_get(row, "residual_mass_kg"),
                          value.residual_mass) &&
              value.residual_mass >= 0;
      if (!valid) break;
      out.injectors.push_back(value);
    }
  }
  yyjson_doc_free(manifest_doc);
  if (!valid) return invalid(24114);
  std::sort(out.injectors.begin(), out.injectors.end(),
            [](const auto &a, const auto &b) { return a.id < b.id; });
  for (std::size_t i = 0; i < model.spray->injectors.size(); ++i) {
    const auto id = model.spray->injectors[i].id;
    const auto found = std::lower_bound(
        out.injectors.begin(), out.injectors.end(), id,
        [](const ImportedInjector &a, std::uint64_t b) { return a.id < b; });
    if (found == out.injectors.end() || found->id != id ||
        found->residual_mass >=
            model.spray->injectors[i].represented_mass_per_parcel_kg)
      return invalid(24114);
  }
  std::ifstream rows(root / file);
  if (!rows) return invalid(24114);
  std::set<std::pair<std::uint64_t, std::uint64_t>> ids;
  std::string line;
  while (std::getline(rows, line)) {
    if (line.empty() || line.size() > 4U * 1024U * 1024U)
      return invalid(24114);
    yyjson_doc *doc = nullptr;
    status = parse_json(line, doc);
    if (!status) return status;
    auto *row = yyjson_doc_get_root(doc);
    auto *source = yyjson_obj_get(row, "source");
    auto *target = yyjson_obj_get(row, "target_thermodynamics");
    ImportedParcel value;
    std::string_view sha, kind;
    std::uint64_t source_rank = 0, ordinal = 0;
    valid = yyjson_is_obj(row) && yyjson_is_obj(source) &&
            yyjson_is_obj(target) &&
            json_number(yyjson_obj_get(row, "droplet_mass_kg"), value.mass) &&
            json_vector3(yyjson_obj_get(source, "position_m"), value.position) &&
            json_vector3(yyjson_obj_get(source, "velocity_m_per_s"), value.velocity) &&
            json_number(yyjson_obj_get(source, "multiplicity"), value.multiplicity) &&
            json_text(yyjson_obj_get(source, "source_sha256"), sha) &&
            sha.size() == 64 &&
            json_u64(yyjson_obj_get(source, "source_rank"), source_rank) &&
            source_rank < UINT64_C(65536) &&
            json_u64(yyjson_obj_get(source, "source_ordinal"), ordinal) &&
            ordinal > 0 && ordinal < (UINT64_C(1) << 47) &&
            json_text(yyjson_obj_get(source, "kind"), kind) &&
            (kind == "retained" || kind == "daughter") &&
            json_number(yyjson_obj_get(target, "droplet_diameter_m"),
                        value.diameter) &&
            json_number(yyjson_obj_get(target, "temperature_k"),
                        value.temperature) &&
            json_u64(yyjson_obj_get(target, "liquid_material_fingerprint"),
                     value.material) &&
            value.mass > 0 && value.diameter > 0 && value.multiplicity > 0 &&
            value.temperature > 0 &&
            value.material == model.spray->liquid_fingerprint;
    value.id_high = hex64(sha);
    value.id_low = (source_rank << 48U) |
                   (std::uint64_t(kind == "daughter") << 47U) | ordinal;
    valid = valid && value.id_high != 0 && value.id_low != 0 &&
            ids.emplace(value.id_high, value.id_low).second;
    yyjson_doc_free(doc);
    if (!valid) return invalid(24114);
    out.parcels.push_back(value);
  }
  if (!rows.eof() || out.parcels.size() != parcel_count)
    return invalid(24114);
  return {};
}

int locate_axis(Span<const double> faces, double value) {
  if (!std::isfinite(value) || faces.size < 2 || value < faces.data[0] ||
      value >= faces.data[faces.size - 1])
    return -1;
  const auto *found = std::upper_bound(faces.data, faces.data + faces.size,
                                       value);
  return static_cast<int>(found - faces.data - 1);
}
bool patch_cell(MeshPatch patch, Int3 global, std::size_t &local) {
  const int x = global.x - patch.begin.x, y = global.y - patch.begin.y,
            z = global.z - patch.begin.z;
  if (x < 0 || y < 0 || z < 0 || x >= patch.cells.x || y >= patch.cells.y ||
      z >= patch.cells.z)
    return false;
  local = std::size_t(x) + std::size_t(patch.cells.x) *
                                (std::size_t(y) +
                                 std::size_t(patch.cells.y) * std::size_t(z));
  return true;
}
void put_u(std::uint8_t *&p, std::uint64_t value, unsigned bytes = 8) {
  for (unsigned i = 0; i < bytes; ++i)
    *p++ = std::uint8_t(value >> (8U * i));
}
void put_real(std::uint8_t *&p, double value) {
  std::uint64_t bits = 0;
  std::memcpy(&bits, &value, 8);
  put_u(p, bits);
}
// The transfer format uses little-endian binary64, cell-major components and
// x-fastest cells.
template <class T>
Status read_box(const std::filesystem::path &path, const Header &h,
                const Block &b, std::size_t components, std::vector<T> &out) {
  std::size_t global = 0, local = 0;
  if (!checked_product(h.cells, global) || !checked_product(b.cells, local) ||
      components == 0 ||
      global > std::numeric_limits<std::size_t>::max() / components / sizeof(T))
    return invalid();
  std::error_code error;
  if (std::filesystem::file_size(path, error) !=
          global * components * sizeof(T) ||
      error)
    return invalid(24105);
  out.resize(local * components);
  std::ifstream in(path, std::ios::binary);
  for (int z = 0; z < b.cells.z; ++z)
    for (int y = 0; y < b.cells.y; ++y) {
      std::size_t g =
          b.begin.x +
          std::size_t(h.cells.x) *
              (b.begin.y + y + std::size_t(h.cells.y) * (b.begin.z + z));
      std::size_t l = std::size_t(b.cells.x) * (y + std::size_t(b.cells.y) * z);
      in.seekg(g * components * sizeof(T));
      in.read(reinterpret_cast<char *>(out.data() + l * components),
              b.cells.x * components * sizeof(T));
      if (!in)
        return invalid(24105);
    }
  return {};
}
Status read_block(const std::filesystem::path &root, const Header &h,
                  MeshPatch patch, Block &b) {
  const std::uint16_t endian = 1;
  if (*reinterpret_cast<const std::uint8_t *>(&endian) != 1 ||
      sizeof(double) != 8 || !std::numeric_limits<double>::is_iec559)
    return invalid();
  b.begin = {std::max(0, patch.begin.x - 1), std::max(0, patch.begin.y - 1),
             std::max(0, patch.begin.z - 1)};
  b.cells = {std::min(h.cells.x, patch.begin.x + patch.cells.x + 1) - b.begin.x,
             std::min(h.cells.y, patch.begin.y + patch.cells.y + 1) - b.begin.y,
             std::min(h.cells.z, patch.begin.z + patch.cells.z + 1) -
                 b.begin.z};
  auto s = read_box(root / "flow.f64", h, b, 4, b.flow);
  if (s)
    s = read_box(root / "rho_ref.f64", h, b, 1, b.reference_density);
  if (s)
    s = read_box(root / "fluid.u8", h, b, 1, b.fluid);
  b.pdf.resize(h.nf);
  for (std::size_t f = 0; s && f < h.nf; ++f)
    s = read_box(root / ("pdf" + std::to_string(f) + ".f64"), h, b, h.ns + 1,
                 b.pdf[f]);
  return s;
}
Status reconstruct(const Header &h, const ThermodynamicsPlan &thermo,
                   const TransportPlan &transport,
                   const std::vector<std::size_t> &independent_source,
                   Block &b) {
  const std::size_t stride = h.ns + 1, count = b.fluid.size();
  if (independent_source.size() + 1U != h.ns)
    return invalid(24110);
  b.mean.assign(count * stride, 0);
  b.density.resize(count);
  b.cache.assign(count * 4, 0);
  std::vector<double> solid(h.ns, 0);
  std::vector<double> independent(h.ns - 1U, 0);
  double sh = 0, cp = 0, r = 0;
  auto status = Status{};
  if (std::find(b.fluid.begin(), b.fluid.end(), std::uint8_t{0}) !=
      b.fluid.end()) {
    const auto oxygen = std::find(h.names.begin(), h.names.end(), "O2");
    const auto nitrogen = std::find(h.names.begin(), h.names.end(), "N2");
    if (oxygen == h.names.end() || nitrogen == h.names.end())
      return invalid();
    // Match the native material's air definition; solid rows carry a valid
    // 295 K state.
    solid[oxygen - h.names.begin()] =
        .21 * 31.998 / (.21 * 31.998 + .79 * 28.014);
    solid[nitrogen - h.names.begin()] =
        1 - solid[oxygen - h.names.begin()];
    for (std::size_t species = 0; species < independent.size(); ++species)
      independent[species] = solid[independent_source[species]];
    status =
        thermo.mixture_enthalpy(295, {independent.data(), independent.size()},
                                sh, cp, r);
    if (!status)
      return status;
  }
  for (std::size_t i = 0; i < count; ++i) {
    if (b.fluid[i] > 1)
      return invalid();
    if (!b.fluid[i]) {
      for (std::size_t a = 0; a < 3; ++a)
        b.flow[4 * i + a] = 0;
      b.flow[4 * i + 3] = h.version == 2 ? 0.0 : h.pressure;
    }
    const double mechanical = b.flow[4 * i + 3];
    const double eos_pressure =
        h.version == 2 ? h.pressure : thermo.eos_pressure(mechanical);
    for (unsigned c = 0; c < 4; ++c)
      if (!std::isfinite(b.flow[4 * i + c]))
        return invalid(24106);
    if ((h.version == 1 && mechanical <= 0) ||
        !std::isfinite(eos_pressure) || eos_pressure <= 0 ||
        !std::isfinite(b.reference_density[i]) ||
        b.reference_density[i] <= 0)
      return invalid(24106);
    for (auto &pdf : b.pdf) {
      double *row = pdf.data() + i * stride;
      if (!b.fluid[i]) {
        std::copy(solid.begin(), solid.end(), row);
        row[h.ns] = sh;
      }
      double sum = 0;
      for (std::size_t s = 0; s < h.ns; ++s) {
        if (!std::isfinite(row[s]) || row[s] < 0 || row[s] > 1)
          return invalid(24107);
        sum += row[s];
      }
      if (std::abs(sum - 1) > 2e-12 || !std::isfinite(row[h.ns]))
        return invalid(24107);
      for (std::size_t c = 0; c < stride; ++c)
        b.mean[i * stride + c] += row[c] / h.nf;
    }
    const double *row = b.mean.data() + i * stride;
    for (std::size_t species = 0; species < independent.size(); ++species)
      independent[species] = row[independent_source[species]];
    ThermoState state;
    status = thermo.evaluate(eos_pressure, row[h.ns],
                             {independent.data(), independent.size()}, {},
                             state);
    if (!status)
      return status;
    b.density[i] = state.rho;
    MolecularTransportState molecular;
    status = transport.evaluate(
        state.temperature, {independent.data(), independent.size()}, molecular);
    if (!status)
      return status;
    double conductivity = 0, gamma = 0;
    if (transport.has_effective_enthalpy_transport()) {
      status = transport.effective_enthalpy_transport(
          molecular.viscosity, molecular.viscosity, state.cp, conductivity,
          gamma);
      if (!status)
        return status;
    } else {
      conductivity = molecular.conductivity;
      gamma = conductivity / state.cp;
      if (!std::isfinite(gamma) || gamma <= 0)
        return invalid(24113);
    }
    b.cache[4 * i] = gamma;
    b.cache[4 * i + 2] = molecular.viscosity;
  }
  return {};
}
// Face history is an explicit interpolation candidate. Native boundary/IBM
// admission follows.
double flux(const Header &h, const CartesianGeometryPlan &g, const Block &b,
            Int3 f, unsigned a) {
  const int index = a == 0 ? f.x : a == 1 ? f.y : f.z;
  const int limit = a == 0 ? h.cells.x : a == 1 ? h.cells.y : h.cells.z;
  Int3 l = f, r = f;
  if (a == 0)
    --l.x;
  else if (a == 1)
    --l.y;
  else
    --l.z;
  auto momentum = [&](Int3 c) {
    auto i = b.offset(c);
    return b.fluid[i] ? b.density[i] * b.flow[4 * i + a] : 0.;
  };
  double m = 0;
  if (index == 0)
    m = momentum(r);
  else if (index == limit)
    m = momentum(l);
  else {
    if (!b.fluid[b.offset(l)] || !b.fluid[b.offset(r)])
      return 0;
    const auto &metric = a == 0 ? g.x() : a == 1 ? g.y() : g.z();
    const double dl =
        metric.faces().data[index] - metric.centres().data[index - 1];
    const double dr = metric.centres().data[index] - metric.faces().data[index];
    m = (dr * momentum(l) + dl * momentum(r)) / (dl + dr);
  }
  return m * (a == 0   ? g.y().widths().data[f.y] * g.z().widths().data[f.z]
              : a == 1 ? g.x().widths().data[f.x] * g.z().widths().data[f.z]
                       : g.x().widths().data[f.x] * g.y().widths().data[f.y]);
}
Status fill(const Header &h, const RestartExpected &e,
            const CartesianGeometryPlan &g, const Block &b,
            const std::vector<std::size_t> &independent_source,
            RestartImage &image) {
  image.global_cells = e.global_cells;
  image.patch = e.target_patch;
  image.plan = e.plan;
  image.schema = e.schema;
  image.geometry = e.geometry;
  image.step = h.step;
  image.time = h.time;
  image.dt = h.dt;
  image.pressure_reference = h.pressure;
  image.controller_state = 1;
  image.source_format_version = 1;
  image.backward_euler_recovery = true;
  image.final_mass_flux_revision = 1;
  std::size_t scalar = 0, stochastic = 0, count = 0;
  if (!checked_product(image.patch.cells, count))
    return invalid();
  const auto c = image.patch.cells;
  for (std::size_t f = 0; f < e.fields.size; ++f) {
    auto d = e.fields.data[f];
    RestartImageField out;
    out.role = d.role;
    out.field = d.field;
    out.components = d.components;
    out.values.resize(count * d.components);
    for (int z = 0; z < c.z; ++z)
      for (int y = 0; y < c.y; ++y)
        for (int x = 0; x < c.x; ++x) {
          const std::size_t cell =
              x + std::size_t(c.x) * (y + std::size_t(c.y) * z);
          const auto i =
              b.offset({x + image.patch.begin.x, y + image.patch.begin.y,
                        z + image.patch.begin.z});
          auto *v = out.values.data() + cell * d.components;
          switch (d.role) {
          case RestartFieldRole::velocity:
            if (d.components != 3)
              return invalid();
            std::copy_n(b.flow.data() + 4 * i, 3, v);
            break;
          case RestartFieldRole::pressure_perturbation:
            if (d.components != 1)
              return invalid();
            *v = h.version == 2 ? b.flow[4 * i + 3]
                                : b.flow[4 * i + 3] - h.pressure;
            break;
          case RestartFieldRole::pressure_absolute:
            if (d.components != 1)
              return invalid();
            *v = h.version == 2 ? h.pressure + b.flow[4 * i + 3]
                                : b.flow[4 * i + 3];
            break;
          case RestartFieldRole::enthalpy:
            if (d.components != 1)
              return invalid();
            *v = b.mean[i * (h.ns + 1) + h.ns];
            break;
          case RestartFieldRole::independent_species:
            if (d.components != 1 || scalar >= independent_source.size())
              return invalid();
            *v = b.mean[i * (h.ns + 1) + independent_source[scalar]];
            break;
          case RestartFieldRole::stochastic_field:
            if (d.components != h.ns + 1 || stochastic >= h.nf)
              return invalid();
            std::copy_n(b.pdf[stochastic].data() + i * (h.ns + 1), h.ns + 1, v);
            break;
          case RestartFieldRole::stochastic_auxiliary:
            if (d.components != h.ns + 1) return invalid();
            // Source PDF checkpoint carries stochastic fields; reconstruct
            // auxiliary method state from the physical mean at import.
            std::copy_n(b.mean.data() + i * (h.ns + 1), h.ns + 1, v);
            break;
          case RestartFieldRole::stochastic_transport:
            if (d.components != 4)
              return invalid();
            std::copy_n(b.cache.data() + 4 * i, 4, v);
            break;
          default:
            return invalid(24108);
          }
        }
    if (d.role == RestartFieldRole::independent_species)
      ++scalar;
    if (d.role == RestartFieldRole::stochastic_field)
      ++stochastic;
    image.fields.push_back(std::move(out));
  }
  if (scalar != independent_source.size() || stochastic != h.nf)
    return invalid(24108);
  for (unsigned a = 0; a < 3; ++a) {
    auto ext = c;
    if (a == 0)
      ++ext.x;
    else if (a == 1)
      ++ext.y;
    else
      ++ext.z;
    if (!checked_product(ext, count))
      return invalid();
    image.final_mass_flux[a].resize(count);
    for (int z = 0; z < ext.z; ++z)
      for (int y = 0; y < ext.y; ++y)
        for (int x = 0; x < ext.x; ++x)
          image.final_mass_flux[a][x + std::size_t(ext.x) *
                                           (y + std::size_t(ext.y) * z)] =
              flux(h, g, b,
                   {x + image.patch.begin.x, y + image.patch.begin.y,
                    z + image.patch.begin.z},
                   a);
  }
  return {};
}
Status assemble_spray_history(MPI_Comm comm, const ValidatedModel &model,
                              const CartesianGeometryPlan &geometry,
                              const RestartExpected &expected,
                              SprayTransfer &transfer, RestartImage &image) {
  if (!model.spray || !model.reaction.esf ||
      expected.cell_record_identity == 0 || expected.cell_record_bytes != 0)
    return invalid(24115);
  const auto global = geometry.global_cells();
  const auto patch = image.patch;
  std::size_t cells = 0;
  if (!checked_product(patch.cells, cells)) return invalid(24115);
  const auto locate = [&](const std::array<double, 3> &point, Int3 &cell) {
    cell = {locate_axis(geometry.x().faces(), point[0]),
            locate_axis(geometry.y().faces(), point[1]),
            locate_axis(geometry.z().faces(), point[2])};
    return cell.x >= 0 && cell.y >= 0 && cell.z >= 0;
  };
  std::vector<std::size_t> local_parcels;
  local_parcels.reserve(transfer.parcels.size());
  std::vector<std::uint32_t> parcel_counts(cells), injector_counts(cells);
  for (std::size_t i = 0; i < transfer.parcels.size(); ++i) {
    Int3 cell;
    if (!locate(transfer.parcels[i].position, cell)) return invalid(24115);
    transfer.parcels[i].owner =
        std::uint64_t(cell.x) + std::uint64_t(global.x) *
                                    (std::uint64_t(cell.y) +
                                     std::uint64_t(global.y) * cell.z);
    std::size_t local = 0;
    if (patch_cell(patch, cell, local)) {
      if (parcel_counts[local] == UINT32_MAX) return invalid(24115);
      ++parcel_counts[local];
      local_parcels.push_back(i);
    }
  }
  std::uint64_t local_count = local_parcels.size(), global_count = 0;
  const int local_capacity_exceeded =
      local_parcels.size() > model.spray->maximum_local_parcels ? 1 : 0;
  int any_capacity_exceeded = 0;
  if (MPI_Allreduce(&local_count, &global_count, 1, MPI_UINT64_T, MPI_SUM,
                    comm) != MPI_SUCCESS ||
      MPI_Allreduce(&local_capacity_exceeded, &any_capacity_exceeded, 1,
                    MPI_INT, MPI_MAX, comm) != MPI_SUCCESS ||
      global_count != transfer.parcels.size())
    return invalid(24115);
  // Capacity is rank-local, so rejection must happen only after every rank
  // has participated in the ownership collectives above.
  if (any_capacity_exceeded)
    return {StatusCode::allocation_failure, 24115};

  struct LocalInjector {
    std::size_t cell{};
    const ImportedInjector *state{};
  };
  std::vector<LocalInjector> local_injectors;
  for (const auto &spec : model.spray->injectors) {
    const auto state = std::lower_bound(
        transfer.injectors.begin(), transfer.injectors.end(), spec.id,
        [](const ImportedInjector &a, std::uint64_t id) { return a.id < id; });
    if (state == transfer.injectors.end() || state->id != spec.id)
      return invalid(24115);
    Int3 cell;
    const std::array<double, 3> origin{spec.origin_m.x, spec.origin_m.y,
                                       spec.origin_m.z};
    if (!locate(origin, cell)) return invalid(24115);
    std::size_t local = 0;
    if (patch_cell(patch, cell, local)) {
      if (injector_counts[local] == UINT32_MAX) return invalid(24115);
      ++injector_counts[local];
      local_injectors.push_back({local, &*state});
    }
  }
  local_count = local_injectors.size();
  if (MPI_Allreduce(&local_count, &global_count, 1, MPI_UINT64_T, MPI_SUM,
                    comm) != MPI_SUCCESS ||
      global_count != model.spray->injectors.size())
    return invalid(24115);

  std::uint32_t tcr_width = 0;
  const auto &tcr = model.reaction.esf->tcr;
  if (tcr.mode != TcrMode::off) {
    if (tcr.model == TcrModel::cdphyso_dynamic_v1) {
      const auto ns = model.thermophysics.species.size();
      if (ns > (UINT32_MAX - 24U) / 40U)
        return invalid(24115);
      tcr_width = static_cast<std::uint32_t>(24U + 40U * ns + 24U);
    } else {
      tcr_width = 120U;
    }
  }
  std::vector<std::uint8_t> tcr_template;
  if (tcr_width) {
    detail::ProductTcrHistory seed;
    if (tcr.model == TcrModel::cdphyso_dynamic_v1)
      seed.configure_dynamic(1, 1, model.thermophysics.species.size(),
                             2.0 / model.reaction.mixing_c_z);
    else
      seed.configure(1, 1, tcr.initialization_sign);
    const auto records = seed.snapshot();
    if (records.record_bytes != tcr_width || records.values.size != tcr_width)
      return invalid(24115);
    tcr_template.assign(records.values.data,
                        records.values.data + records.values.size);
    auto *clock = tcr_template.data();
    put_u(clock, image.step);
    if (tcr.model == TcrModel::cdphyso_dynamic_v1) {
      clock = tcr_template.data() + 8U;
      put_u(clock, image.step % 4U);
    }
  }
  const std::uint32_t parcel_width =
      model.spray->breakup == SprayBreakupModel::stochastic_sgs ? 192U : 144U;
  const std::uint32_t version = parcel_width == 192U ? 2U : 1U;
  image.cell_record_identity = expected.cell_record_identity;
  image.cell_record_bytes = 0;
  image.cell_record_lengths.resize(cells);
  std::vector<std::size_t> offsets(cells + 1), cursors(cells);
  for (std::size_t cell = 0; cell < cells; ++cell) {
    const std::uint64_t length =
        (tcr_width || parcel_counts[cell] || injector_counts[cell]
             ? 24U + std::uint64_t(tcr_width)
             : 0U) +
        std::uint64_t(parcel_counts[cell]) * parcel_width +
        std::uint64_t(injector_counts[cell]) * 24U;
    if (length > UINT32_MAX || offsets[cell] > SIZE_MAX - length)
      return invalid(24115);
    image.cell_record_lengths[cell] = static_cast<std::uint32_t>(length);
    offsets[cell + 1] = offsets[cell] + static_cast<std::size_t>(length);
  }
  image.cell_records.resize(offsets.back());
  for (std::size_t cell = 0; cell < cells; ++cell) {
    if (!image.cell_record_lengths[cell]) continue;
    auto *p = image.cell_records.data() + offsets[cell];
    put_u(p, image.step);
    put_u(p, parcel_counts[cell], 4);
    put_u(p, injector_counts[cell], 4);
    put_u(p, tcr_width, 4);
    put_u(p, version, 4);
    if (tcr_width) {
      std::memcpy(p, tcr_template.data(), tcr_template.size());
      p += tcr_template.size();
    }
    cursors[cell] = static_cast<std::size_t>(p - image.cell_records.data());
  }
  for (const auto index : local_parcels) {
    const auto &value = transfer.parcels[index];
    const auto gx = int(value.owner % std::uint64_t(global.x));
    const auto gy = int((value.owner / std::uint64_t(global.x)) %
                        std::uint64_t(global.y));
    const auto gz = int(value.owner /
                        (std::uint64_t(global.x) * std::uint64_t(global.y)));
    std::size_t local = 0;
    if (!patch_cell(patch, {gx, gy, gz}, local)) return invalid(24115);
    auto *p = image.cell_records.data() + cursors[local];
    put_u(p, value.id_high);
    put_u(p, value.id_low);
    for (const double v : value.position) put_real(p, v);
    for (const double v : value.velocity) put_real(p, v);
    put_real(p, value.mass);
    put_real(p, value.diameter);
    put_real(p, value.multiplicity);
    put_real(p, value.temperature);
    put_u(p, value.material);
    put_u(p, value.owner);
    put_real(p, 0);
    put_real(p, 0);
    put_real(p, 0);
    put_u(p, 0);
    if (version == 2U) {
      // Migration deliberately resets all SGS exposure clocks, but a parcel
      // participating in stochastic-SGS breakup must still carry the active
      // history schema.  An all-zero wire image means "SGS absent" (version
      // zero) and is rejected by the first runtime advance.
      put_u(p, 1);
      for (unsigned lane = 1; lane < 6; ++lane) put_u(p, 0);
    }
    cursors[local] = static_cast<std::size_t>(p - image.cell_records.data());
  }
  for (const auto &value : local_injectors) {
    auto *p = image.cell_records.data() + cursors[value.cell];
    put_u(p, value.state->id);
    put_real(p, value.state->residual_mass);
    put_u(p, value.state->next_ordinal);
    cursors[value.cell] = static_cast<std::size_t>(p - image.cell_records.data());
  }
  for (std::size_t cell = 0; cell < cells; ++cell)
    if (image.cell_record_lengths[cell] && cursors[cell] != offsets[cell + 1])
      return invalid(24115);
  if (tcr_width) {
    detail::ProductTcrHistory probe;
    if (tcr.model == TcrModel::cdphyso_dynamic_v1)
      probe.configure_dynamic(1, 1, model.thermophysics.species.size(),
                              2.0 / model.reaction.mixing_c_z);
    else
      probe.configure(1, 1, tcr.initialization_sign);
    const auto checked = probe.stage_restore_records(
        {image.cell_records.data() + 24U, tcr_width}, image.step);
    if (!checked) return invalid(24117);
  }
  return {};
}

Status make_complete_migration(MPI_Comm comm, const ValidatedModel &model,
                               const CartesianGeometryPlan &geometry,
                               const RestartExpected &expected,
                               double closed_mass, SprayTransfer &transfer,
                               RestartImage &image) {
  if (!std::isfinite(closed_mass) || closed_mass <= 0 ||
      expected.method_history_signature == 0)
    return invalid(24116);
  std::size_t cells = 0;
  if (!checked_product(image.patch.cells, cells)) return invalid(24116);
  image.source_format_version = 5;
  image.backward_euler_recovery = false;
  image.previous_fields = image.fields;
  image.accepted_rate_fields.clear();
  image.previous_rate_fields.clear();
  for (std::size_t i = 0; i < expected.rate_fields.size; ++i) {
    const auto descriptor = expected.rate_fields.data[i];
    if (descriptor.components == 0 ||
        cells > SIZE_MAX / descriptor.components)
      return invalid(24116);
    RestartImageField field;
    field.role = descriptor.role;
    field.field = descriptor.field;
    field.components = descriptor.components;
    field.values.assign(cells * descriptor.components, 0.0);
    image.accepted_rate_fields.push_back(field);
    image.previous_rate_fields.push_back(std::move(field));
  }
  image.previous_mass_flux = image.final_mass_flux;
  image.previous_pressure_reference = image.pressure_reference;
  image.closed_mass_target = closed_mass;
  image.previous_mass_flux_revision = 16;
  image.final_mass_flux_revision = 17;
  image.method_history_signature = expected.method_history_signature;
  return assemble_spray_history(comm, model, geometry, expected, transfer,
                                image);
}
Status exact_physical_fields(const RestartImage &image,
                             const RestartSnapshot &snapshot) {
  if (image.fields.size() != snapshot.fields.size)
    return invalid(24109);
  for (std::size_t f = 0; f < image.fields.size(); ++f) {
    const auto &field = image.fields[f];
    const auto view = snapshot.fields.data[f].values;
    if (field.role == RestartFieldRole::stochastic_transport)
      continue;
    if (field.role != snapshot.fields.data[f].role ||
        field.components != view.components)
      return invalid(24109);
    std::size_t i = 0;
    for (int z = 0; z < view.interior.z; ++z)
      for (int y = 0; y < view.interior.y; ++y)
        for (int x = 0; x < view.interior.x; ++x)
          for (unsigned c = 0; c < view.components; ++c)
            if (field.values[i++] != view.unchecked({x, y, z}, c))
              return invalid(24109);
    if (i != field.values.size())
      return invalid(24109);
  }
  return {};
}
int run(const char *case_root, const char *transfer, const char *output,
        const char *spray_root, int rank) {
  const auto comm = MPI_COMM_WORLD;
  auto stage = [&](const char *name, Status s) {
    s = consensus(comm, s);
    if (rank == 0)
      std::cout << name << " status=" << unsigned(s.code) << '/' << s.detail
                << std::endl;
    return bool(s);
  };
  Header h;
  if (!stage("header",
             local_stage(comm, [&] { return read_header(transfer, h); })))
    return 3;
  ValidatedModel model;
  auto s = CaseCompiler::load_and_compile(comm, case_root, model);
  std::vector<std::size_t> independent_source;
  if (s && (!model.reaction.esf || model.reaction.esf->fields != h.nf ||
            model.transported_scalars.size() != h.ns - 1 ||
            model.thermophysics.species.size() != h.ns))
    s = invalid(24110);
  for (std::size_t i = 0; s && i < h.ns; ++i) {
    if (model.thermophysics.species[i].stable_name != h.names[i])
      s = invalid(24110);
  }
  for (const auto &scalar : model.transported_scalars) {
    if (!s)
      break;
    const auto found = std::find(h.names.begin(), h.names.end(),
                                 scalar.stable_name);
    if (scalar.role != TransportedScalarRole::species ||
        found == h.names.end()) {
      s = invalid(24110);
      break;
    }
    const auto source = static_cast<std::size_t>(found - h.names.begin());
    if (std::find(independent_source.begin(), independent_source.end(),
                  source) != independent_source.end()) {
      s = invalid(24110);
      break;
    }
    independent_source.push_back(source);
  }
  if (s && h.version == 2 &&
      (model.thermophysics.fixed_pressure_pa <= 0 ||
       model.thermophysics.fixed_pressure_pa != h.pressure))
    s = invalid(24112);
  if (s && (model.spray.has_value() != (spray_root != nullptr)))
    s = invalid(24114);
  if (!stage("case", s))
    return 4;
  SprayTransfer spray;
  if (spray_root &&
      !stage("spray_transfer", local_stage(comm, [&] {
               return read_spray_transfer(spray_root, model, spray);
             })))
    return 4;
  ThermodynamicsPlan thermo;
  TransportPlan transport_plan;
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  s = ThermodynamicsPlan::compile(
      model.thermophysics,
      {model.transported_scalars.data(), model.transported_scalars.size()},
      thermo);
  if (s)
    s = TransportPlan::compile(model.thermophysics, thermo, transport_plan);
  if (!stage("thermo", s))
    return 4;
  s = CartesianGeometryCompiler::compile(comm, model.mesh, {}, geometry, patch);
  if (s && !same_int3(h.cells, geometry.global_cells()))
    s = invalid();
  if (!stage("geometry", s))
    return 4;
  Block b;
  if (!stage("read", local_stage(comm, [&] {
               return read_block(transfer, h, patch, b);
             })))
    return 5;
  if (!stage("reconstruct", local_stage(comm, [&] {
               return reconstruct(h, thermo, transport_plan,
                                  independent_source, b);
             })))
    return 5;
  // Quantify canonical native mean EOS relative to the COAST startup PDF EOS.
  std::array<double, 6> audit{};
  std::array<double, 2> pressure_range{
      std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity()};
  for (int z = 0; z < patch.cells.z; ++z)
    for (int y = 0; y < patch.cells.y; ++y)
      for (int x = 0; x < patch.cells.x; ++x) {
        const Int3 global{patch.begin.x + x, patch.begin.y + y,
                          patch.begin.z + z};
        const auto i = b.offset(global);
        if (!b.fluid[i])
          continue;
        const double v = geometry.x().widths().data[global.x] *
                         geometry.y().widths().data[global.y] *
                         geometry.z().widths().data[global.z];
        const double ke = .5 * (b.flow[4 * i] * b.flow[4 * i] +
                                b.flow[4 * i + 1] * b.flow[4 * i + 1] +
                                b.flow[4 * i + 2] * b.flow[4 * i + 2]);
        const double eh = b.mean[i * (h.ns + 1) + h.ns] + ke;
        const double eos_pressure = h.version == 2
                                        ? h.pressure
                                        : thermo.eos_pressure(b.flow[4 * i + 3]);
        audit[0] += b.reference_density[i] * v;
        audit[1] += b.density[i] * v;
        audit[2] += (b.reference_density[i] * eh - eos_pressure) * v;
        audit[3] += (b.density[i] * eh - eos_pressure) * v;
        audit[4] +=
            std::abs(b.reference_density[i] * eh - eos_pressure) * v;
        audit[5] += 1;
        pressure_range[0] = std::min(pressure_range[0], b.flow[4 * i + 3]);
        pressure_range[1] = std::max(pressure_range[1], b.flow[4 * i + 3]);
      }
  MPI_Allreduce(MPI_IN_PLACE, audit.data(), 6, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, pressure_range.data(), 1, MPI_DOUBLE, MPI_MIN,
                comm);
  MPI_Allreduce(MPI_IN_PLACE, pressure_range.data() + 1, 1, MPI_DOUBLE,
                MPI_MAX, comm);
  if (audit[5] <= 0 || !std::isfinite(pressure_range[0]) ||
      !std::isfinite(pressure_range[1]))
    s = invalid(24111);
  if (rank == 0) {
    std::ofstream file(std::filesystem::path(transfer) / "native.json");
    file << std::setprecision(17)
         << "{\"transfer_format_version\":" << h.version
         << ",\"pressure_semantics\":\""
         << (h.version == 2 ? "mechanical_perturbation" : "absolute")
         << "\",\"eos_pressure_pa\":";
    if (h.version == 2)
      file << h.pressure;
    else
      file << "null";
    file << ",\"pressure_field_range_pa\":[" << pressure_range[0] << ','
         << pressure_range[1] << ']';
    if (h.version == 2)
      file << ",\"mechanical_pressure_range_pa\":[" << pressure_range[0]
           << ',' << pressure_range[1] << ']';
    file << ",\"reference_mass_kg\":" << audit[0]
         << ",\"native_mass_kg\":" << audit[1]
         << ",\"relative_mass_change\":" << (audit[1] - audit[0]) / audit[0]
         << ",\"reference_energy_J\":" << audit[2]
         << ",\"native_energy_J\":" << audit[3]
         << ",\"relative_energy_change\":" << (audit[3] - audit[2]) / audit[4]
         << ",\"fluid_cells\":" << audit[5]
         << ",\"scope\":\"fixed PDF h/Y and declared pressure semantics, "
            "canonical native mean EOS versus "
            "COAST startup harmonic PDF EOS\"}\n";
    if (!file)
      s = invalid(24111);
  }
  if (!stage("audit", s))
    return 5;
  CompiledCasePlan plan;
  ProductDriver driver;
  s = ProductCompiler::compile(comm, model, case_root, plan);
  if (s)
    s = ProductDriver::create(comm, std::move(plan), driver);
  if (!stage("product", s))
    return 6;
  RestartExpected expected;
  s = driver.restart_expected(
      expected, RestartStorageCompatibility::strict,
      spray_root ? RestartHistoryPolicy::rebuild_method_history
                 : RestartHistoryPolicy::require_compatible);
  if (s && !same_patch(expected.target_patch, patch))
    s = invalid();
  RestartImage image;
  if (s)
    s = local_stage(comm, [&] {
      return fill(h, expected, geometry, b, independent_source, image);
    });
  if (s && spray_root)
    s = local_stage(comm, [&] {
      return make_complete_migration(comm, model, geometry, expected, audit[1],
                                     spray, image);
    });
  if (!stage("image", s))
    return 6;
  b = Block{};
  s = driver.initialize_restart(image, RestartStorageCompatibility::strict,
                                RestartHistoryPolicy::rebuild_method_history);
  if (!stage("restore", s))
    return 7;
  RestartSnapshot snapshot;
  s = driver.committed_restart_snapshot(snapshot);
  if (s)
    s = exact_physical_fields(image, snapshot);
  if (!stage("physical_fields", s))
    return 7;
  // The gas-only bridge deliberately remains V1. Spray migration has already
  // staged a complete V5 ownership transaction and must retain its records.
  if (!spray_root) {
    snapshot.previous_fields = {};
    snapshot.accepted_rate_fields = {};
    snapshot.previous_rate_fields = {};
    snapshot.previous_mass_flux = {};
    snapshot.previous_pressure_reference = 0;
    snapshot.closed_mass_target = 0;
    snapshot.method_history_signature = 0;
  }
  s = RestartWriter::write(comm, output, snapshot, {1U, nullptr});
  if (!stage("write", s))
    return 8;
  // Read into the same product's expected partition after the synchronous
  // write.
  RestartImage back;
  RestartReadReport report;
  s = RestartReader::load(comm, output, expected, back, &report);
  const unsigned wanted_version = spray_root ? 5U : 1U;
  if (s && (back.source_format_version != wanted_version ||
            back.backward_euler_recovery == (spray_root != nullptr) ||
            back.step != h.step || back.time != h.time || back.dt != h.dt))
    s = invalid(24109);
  if (s && back.fields.size() != image.fields.size())
    s = invalid(24109);
  for (std::size_t f = 0; s && f < image.fields.size(); ++f) {
    if (image.fields[f].role != RestartFieldRole::stochastic_transport &&
        back.fields[f].values != image.fields[f].values)
      s = invalid(24109);
  }
  if (!stage("readback", s))
    return 9;
  driver = ProductDriver{};
  image = RestartImage{};
  s = ProductCompiler::compile(comm, model, case_root, plan);
  if (s)
    s = ProductDriver::create(comm, std::move(plan), driver);
  if (s)
    s = driver.initialize_restart(back);
  if (s)
    s = driver.committed_restart_snapshot(snapshot);
  if (s)
    s = exact_physical_fields(back, snapshot);
  if (!stage("native_restart", s))
    return 10;
  if (rank == 0)
    std::cout << "PDF_IMPORT_OK step=" << h.step << " fields=" << h.nf
              << " species=" << h.ns
              << " pressure="
              << (h.version == 2 ? "mechanical_perturbation" : "absolute")
              << " eos="
              << (h.version == 2 ? "fixed_thermodynamic"
                                 : "coupled_absolute")
              << " physical_readback=exact auxiliary=physical_mean_reconstruction history="
              << (spray_root ? "V5_migration spray_parcels=" +
                                     std::to_string(spray.parcels.size())
                             : "V1_recovery")
              << std::endl;
  return 0;
}
} // namespace
int main(int argc, char **argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
    return 2;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int result = 2;
  if (argc != 4 && argc != 5) {
    if (rank == 0)
      std::cerr << "usage: v04_pdf_import CASE TRANSFER RESTART [SPRAY_TRANSFER]\n";
  } else {
    try {
      result = run(argv[1], argv[2], argv[3], argc == 5 ? argv[4] : nullptr,
                   rank);
    } catch (const std::exception &e) {
      std::cerr << "rank=" << rank << " exception=" << e.what() << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 11);
    }
  }
  MPI_Finalize();
  return result;
}
