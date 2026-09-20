// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn
// Current-state PDF bridge. Missing temporal histories are explicitly rebuilt.
#include "hundun/v04_app.hpp"
#include "app_import_detail.hpp"
#include "core_tcr_dyn711_history_detail.hpp"
#include "core_tcr_dynamic_history_detail.hpp"
#include "esf_count_detail.hpp"
#include "hundun/v04_io.hpp"
#include "hundun/v04_mesh.hpp"
#include "hundun/v04_physics.hpp"

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
  Int3 cells{};
  std::uint64_t step{};
  double time{}, dt{}, pressure{};
  std::size_t nf{}, ns{};
  std::vector<std::string> names;
  std::vector<std::size_t> independent, passive_mapping;
  std::vector<std::string> passive_names;
  unsigned version{};
};
Status invalid(unsigned detail = 24104) {
  return {StatusCode::invalid_case, detail};
}
Status read_header(const std::filesystem::path &root, Header &h) {
  std::ifstream in(root / "state.txt");
  std::string magic;
  unsigned version = 0;
  in >> magic >> version >> h.cells.x >> h.cells.y >> h.cells.z >> h.step >>
      h.time >> h.dt >> h.pressure >> h.nf >> h.ns;
  std::size_t count = 0;
  if (!in || magic != "HUNDUN_PDF_TRANSFER" || (version != 1 && version != 2) ||
      !checked_product(h.cells, count) || !esf::valid_field_count(h.nf) ||
      h.ns < 2 || h.ns > 64 || h.step == 0 || !std::isfinite(h.time) ||
      h.time < 0 || !std::isfinite(h.dt) || h.dt <= 0 ||
      !std::isfinite(h.pressure) || h.pressure <= 0)
    return invalid();
  h.version=version;
  std::size_t passives=0;
  if(version==2 && (!(in>>passives) || passives>64)) return invalid();
  h.passive_names.resize(passives);
  h.names.resize(h.ns);
  for (auto &name : h.names)
    if (!(in >> name))
      return invalid();
  for(auto& name:h.passive_names) if(!(in>>name)) return invalid();
  for(std::size_t i=0;i<h.passive_names.size();++i) {
    if(std::find(h.names.begin(),h.names.end(),h.passive_names[i])!=h.names.end()) return invalid();
    if(std::find(h.passive_names.begin(),h.passive_names.begin()+i,h.passive_names[i])!=h.passive_names.begin()+i) return invalid();
  }
  std::string extra;
  if (in >> extra)
    return invalid();
  return {};
}
struct Block {
  Int3 begin{}, cells{};
  std::vector<double> flow, reference_density, density, mean, cache, passives;
  std::vector<std::vector<double>> pdf;
  std::vector<std::uint8_t> fluid;
  std::size_t offset(Int3 g) const {
    return (g.x - begin.x) +
           std::size_t(cells.x) *
               ((g.y - begin.y) + std::size_t(cells.y) * (g.z - begin.z));
  }
};
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
  if(s && !h.passive_names.empty())
    s=read_box(root / "passive.f64",h,b,h.passive_names.size(),b.passives);
  return s;
}
Status reconstruct(const Header &h, const ThermodynamicsPlan &thermo,
                   const TransportPlan &transport, Block &b) {
  const std::size_t stride = h.ns + 1, count = b.fluid.size();
  for(double value:b.passives) if(!std::isfinite(value)) return invalid(24107);
  b.mean.assign(count * stride, 0);
  b.density.resize(count);
  b.cache.assign(count * 4, 0);
  std::vector<double> solid(h.ns, 0);
  const auto oxygen = std::find(h.names.begin(), h.names.end(), "O2");
  const auto nitrogen = std::find(h.names.begin(), h.names.end(), "N2");
  if (oxygen == h.names.end() || nitrogen == h.names.end())
    return invalid();
  // Match the native material's air definition; solid rows carry a valid 295 K
  // state.
  solid[oxygen - h.names.begin()] =
      .21 * 31.998 / (.21 * 31.998 + .79 * 28.014);
  solid[nitrogen - h.names.begin()] = 1 - solid[oxygen - h.names.begin()];
  std::vector<double> independent(h.independent.size());
  for (std::size_t s = 0; s < independent.size(); ++s)
    independent[s] = solid[h.independent[s]];
  double sh = 0, cp = 0, r = 0;
  auto status =
      thermo.mixture_enthalpy(295, {independent.data(), independent.size()}, sh, cp, r);
  if (!status)
    return status;
  for (std::size_t i = 0; i < count; ++i) {
    if (b.fluid[i] > 1)
      return invalid();
    if (!b.fluid[i]) {
      for (std::size_t a = 0; a < 3; ++a)
        b.flow[4 * i + a] = 0;
      b.flow[4 * i + 3] = h.pressure;
    }
    const double p = b.flow[4 * i + 3];
    for (unsigned c = 0; c < 4; ++c)
      if (!std::isfinite(b.flow[4 * i + c]))
        return invalid(24106);
    if (thermo.eos_pressure(p) <= 0 || !std::isfinite(b.reference_density[i]) ||
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
    for (std::size_t s = 0; s < independent.size(); ++s)
      independent[s] = row[h.independent[s]];
    ThermoState state;
    status = thermo.evaluate(p, row[h.ns], {independent.data(), independent.size()}, {}, state);
    if (!status)
      return status;
    b.density[i] = state.rho;
    MolecularTransportState molecular;
    status = transport.evaluate(state.temperature, {independent.data(), independent.size()}, molecular);
    if (!status)
      return status;
    double conductivity = 0, gamma = 0;
    status = transport.effective_enthalpy_transport(
        molecular.viscosity, molecular.viscosity, state.cp, conductivity,
        gamma);
    if (!status)
      return status;
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
  std::size_t scalar = 0, passive = 0, stochastic = 0, count = 0;
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
            *v = b.flow[4 * i + 3] - h.pressure;
            break;
          case RestartFieldRole::pressure_absolute:
            if (d.components != 1)
              return invalid();
            *v = b.flow[4 * i + 3];
            break;
          case RestartFieldRole::enthalpy:
            if (d.components != 1)
              return invalid();
            *v = b.mean[i * (h.ns + 1) + h.ns];
            break;
          case RestartFieldRole::independent_species:
            if (d.components != 1 || scalar >= h.ns - 1)
              return invalid();
            *v = b.mean[i * (h.ns + 1) + h.independent[scalar]];
            break;
          case RestartFieldRole::transported_scalar:
            if(d.components!=1 || passive>=h.passive_mapping.size()) return invalid();
            *v=b.passives[i*h.passive_names.size()+h.passive_mapping[passive]];
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
    if (d.role == RestartFieldRole::transported_scalar) ++passive;
    if (d.role == RestartFieldRole::stochastic_field)
      ++stochastic;
    image.fields.push_back(std::move(out));
  }
  if (scalar != h.ns - 1 || stochastic != h.nf || passive!=h.passive_names.size())
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
        int rank, bool legacy_report, unsigned expected_version, bool initialize_model_history) {
  const auto comm = MPI_COMM_WORLD;
  auto stage = [&](const char *name, Status s) {
    s = consensus(comm, s);
    if (rank == 0)
      std::cout << name << " status=" << unsigned(s.code) << '/' << s.detail
                << std::endl;
    return bool(s);
  };
  // Admission is collective and precedes every output mutation.
  if (!stage("paths", local_stage(comm, [&]() -> Status {
        if (rank != 0) return {};
        namespace fs = std::filesystem;
        const auto destination = fs::weakly_canonical(output);
        const auto case_path = fs::weakly_canonical(case_root);
        const auto source = fs::weakly_canonical(transfer);
        auto contains = [](const fs::path& parent, const fs::path& child) {
          auto p=parent.begin(), c=child.begin();
          for (; p!=parent.end() && c!=child.end() && *p==*c; ++p,++c) {}
          return p==parent.end();
        };
        if (fs::exists(destination) || contains(destination,source) ||
            contains(source,destination) || contains(destination,case_path) ||
            contains(case_path,destination)) return invalid(24112);
        return {};
      }))) return 3;
  Header h;
  if (!stage("header",
             local_stage(comm, [&] {
               auto status=read_header(transfer,h);
               if(status && expected_version && h.version!=expected_version) return invalid();
               return status;
             })))
    return 3;
  ValidatedModel model;
  auto s = CaseCompiler::load_and_compile(comm, case_root, model);
  if (s && (!model.reaction.esf || model.reaction.esf->fields != h.nf ||
            model.transported_scalars.size() != h.ns - 1 + h.passive_names.size() ||
            model.thermophysics.species.size() != h.ns))
    s = invalid(24110);
  for (std::size_t i = 0; s && i < h.ns; ++i) {
    if (model.thermophysics.species[i].stable_name != h.names[i])
      s = invalid(24110);
  }
  s = local_stage(comm, [&]() -> Status {
    if (!s) return s;
    for (const auto& scalar : model.transported_scalars) {
      if(scalar.role==TransportedScalarRole::passive_scalar) {
        const auto found=std::find(h.passive_names.begin(),h.passive_names.end(),scalar.stable_name);
        if(found==h.passive_names.end()) return invalid(24110);
        h.passive_mapping.push_back(static_cast<std::size_t>(found-h.passive_names.begin()));
      } else {
        const auto found = std::find(h.names.begin(), h.names.end(), scalar.stable_name);
        if (scalar.role != TransportedScalarRole::species || found == h.names.end())
          return invalid(24110);
        h.independent.push_back(static_cast<std::size_t>(found-h.names.begin()));
      }
    }
    if(h.independent.size()!=h.ns-1 || h.passive_mapping.size()!=h.passive_names.size()) return invalid(24110);
    return {};
  });
  if(s && initialize_model_history &&
      (!model.reaction.esf || model.reaction.esf->tcr.mode==TcrMode::off ||
       !dynamic_tcr_model(model.reaction.esf->tcr.model))) s=invalid(24114);
  if (!stage("case", s))
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
               return reconstruct(h, thermo, transport_plan, b);
             })))
    return 5;
  // Quantify canonical native mean EOS relative to the REFERENCE startup PDF EOS.
  std::array<double, 6> audit{};
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
        const double eos_pressure = thermo.eos_pressure(b.flow[4 * i + 3]);
        audit[0] += b.reference_density[i] * v;
        audit[1] += b.density[i] * v;
        audit[2] += (b.reference_density[i] * eh - eos_pressure) * v;
        audit[3] += (b.density[i] * eh - eos_pressure) * v;
        audit[4] +=
            std::abs(b.reference_density[i] * eh - eos_pressure) * v;
        audit[5] += 1;
      }
  if (MPI_Allreduce(MPI_IN_PLACE, audit.data(), 6, MPI_DOUBLE, MPI_SUM, comm)!=MPI_SUCCESS)
    s={StatusCode::mpi_failure,kCollective};
  if (s && (!(audit[0]>0) || !(audit[5]>0) ||
            !std::all_of(audit.begin(),audit.end(),[](double v){return std::isfinite(v);})))
    s=invalid(24113);
  if (!stage("inventory",s)) return 5;
  std::ostringstream report_text;
  if (rank == 0) {
    auto& file=report_text;
    file << std::setprecision(17) << "{\"reference_mass_kg\":" << audit[0]
         << ",\"native_mass_kg\":" << audit[1]
         << ",\"relative_mass_change\":" << (audit[1] - audit[0]) / audit[0]
         << ",\"reference_energy_J\":" << audit[2]
         << ",\"native_energy_J\":" << audit[3]
         << ",\"relative_energy_change\":";
    if (audit[4]>0) file << (audit[3]-audit[2])/audit[4];
    else file << "null";
    file << ",\"energy_normalization_J\":" << audit[4]
         << ",\"absolute_energy_change_J\":" << audit[3]-audit[2]
         << ",\"fluid_cells\":" << audit[5]
         << ",\"pressure_model\":\""
         << (thermo.fixed_pressure_pa()>0 ? "fixed_thermodynamic" : "coupled_eos")
         << "\",\"mechanical_pressure_reference_pa\":" << h.pressure
         << ",\"thermodynamic_pressure_pa\":";
    if (thermo.fixed_pressure_pa()>0) file << thermo.fixed_pressure_pa();
    else file << "null";
    file << ",\"mechanical_pressure\":\"flow.f64 component 3 = reference + pi\""
         << ",\"energy_definition\":\"rho*(mean_h+kinetic)-p_eos\""
         << ",\"source_step\":" << h.step << ",\"source_time_s\":" << h.time
         << ",\"source_dt_s\":" << h.dt
         << ",\"auxiliary\":\"physical_mean_reconstruction\""
         << ",\"history\":\"" << (initialize_model_history ? "current_state_v6_model_initialization" : "current_state_v1_rebuild")
         << "\",\"model_history_policy\":\"" << (initialize_model_history ? "initialize_target_model" : "source_current_state_only")
         << "\",\"model_initialization\":";
    if(initialize_model_history) {
      const bool dyn=model.reaction.esf->tcr.model==TcrModel::dyn711_v1;
      file << "{\"model\":\"" << (dyn ? "dyn711_v1" : "cdphyso_dynamic_v1")
           << "\",\"flow_step\":" << h.step << ",\"initial_kappa\":" << (dyn ? .2 : 1.)
           << ",\"initial_mixing_coefficient\":" << 2/model.reaction.mixing_c_z
           << ",\"statistics_calls\":";
      if(dyn)file << 0;else file << "null";
      file << ",\"flow_step_mod4\":" << h.step%4 << '}';
    } else file << "null";
    file << ",\"species_order\":[";
    for (std::size_t i=0; i<h.names.size(); ++i)
      file << (i ? "," : "") << std::quoted(h.names[i]);
    file << "],\"independent_species_order\":[";
    for (std::size_t i=0; i<h.independent.size(); ++i)
      file << (i ? "," : "") << std::quoted(h.names[h.independent[i]]);
    file << "],\"transfer_version\":" << h.version << ",\"passive_scalar_order\":[";
    for(std::size_t i=0;i<h.passive_mapping.size();++i)
      file << (i ? "," : "") << std::quoted(h.passive_names[h.passive_mapping[i]]);
    file << "]"
         << ",\"scope\":\"fixed PDF p/h/Y, canonical native mean EOS versus "
            "REFERENCE startup harmonic PDF EOS\"}\n";
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
  s = driver.restart_expected(expected);
  if (s && !same_patch(expected.target_patch, patch))
    s = invalid();
  RestartImage image;
  if (s)
    s = local_stage(comm,
                    [&] { return fill(h, expected, geometry, b, image); });
  if(s && initialize_model_history) s=local_stage(comm,[&]() -> Status {
    std::size_t cells=0;
    if(!checked_product(image.patch.cells,cells) || !expected.cell_record_identity) return invalid(24114);
    detail::Dyn711History dyn711;
    detail::DynamicTcrHistory dynamic;
    RestartCellRecordsView records;
    const double initial_coefficient=2/model.reaction.mixing_c_z;
    // The target compiler owns the identity; these typed constructors supply
    // its documented fresh statistics at the imported global flow step.
    if(model.reaction.esf->tcr.model==TcrModel::dyn711_v1) {
      dyn711.configure(0,cells,h.ns,initial_coefficient,h.step);records=dyn711.snapshot();
    } else {
      dynamic.configure(0,cells,h.ns,initial_coefficient,h.step);records=dynamic.snapshot();
    }
    if(records.record_bytes!=expected.cell_record_bytes) return invalid(24114);
    image.cell_record_identity=expected.cell_record_identity;
    image.cell_record_bytes=records.record_bytes;
    image.cell_records.assign(records.values.data,records.values.data+records.values.size);
    image.source_format_version=6;
    return {};
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
  // V1 explicitly marks absent history; the next native step recovers with BE.
  snapshot.previous_fields = {};
  snapshot.accepted_rate_fields = {};
  snapshot.previous_rate_fields = {};
  snapshot.previous_mass_flux = {};
  snapshot.previous_pressure_reference = 0;
  snapshot.closed_mass_target = 0;
  snapshot.method_history_signature = 0;
  s = RestartWriter::write(comm, output, snapshot, {1U, nullptr});
  if (!stage("write", s))
    return 8;
  // Read into the same product's expected partition after the synchronous
  // write.
  RestartImage back;
  RestartReadReport report;
  s = RestartReader::load(comm, output, expected, back, &report);
  if (s && (back.source_format_version != (initialize_model_history ? 6U : 1U) || !back.backward_euler_recovery ||
            back.step != h.step || back.time != h.time || back.dt != h.dt))
    s = invalid(24109);
  if (s && back.fields.size() != image.fields.size())
    s = invalid(24109);
  if (s && initialize_model_history && back.cell_records!=image.cell_records)
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
  s=local_stage(comm,[&]() -> Status {
    if(rank!=0) return {};
    const auto report_path=legacy_report
        ? std::filesystem::path(transfer)/"native.json"
        : std::filesystem::path(output)/"import.json";
    std::ofstream file(report_path,std::ios::binary|std::ios::trunc);
    file << report_text.str();
    file.close();
    return file ? Status{} : invalid(24111);
  });
  if (!stage("report",s)) return 10;
  if (rank == 0)
    std::cout << "PDF_IMPORT_OK step=" << h.step << " fields=" << h.nf
              << " species=" << h.ns
              << " physical_readback=exact auxiliary=physical_mean_reconstruction history=" << (initialize_model_history ? "V6_model_initialization" : "V1_recovery") << std::endl;
  return 0;
}
} // namespace

namespace hundun::v04::detail {
int import_pdf_transfer(const char* case_root, const char* transfer,
                        const char* output, bool legacy_report, unsigned expected_version, bool initialize_model_history) {
  int rank=0;
  if (MPI_Comm_rank(MPI_COMM_WORLD,&rank)!=MPI_SUCCESS) return 2;
  try {
    return run(case_root,transfer,output,rank,legacy_report,expected_version,initialize_model_history);
  } catch (const std::exception& error) {
    std::cerr << "rank=" << rank << " exception=" << error.what() << std::endl;
    MPI_Abort(MPI_COMM_WORLD,11);
    return 11;
  }
}
} // namespace hundun::v04::detail
