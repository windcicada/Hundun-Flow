// SPDX-License-Identifier: Apache-2.0
#include "core_spray_advance_detail.hpp"
#include "models_chemistry_adapter_detail.hpp"
#include <climits>
#include <cstdlib>
#include <iostream>
#include <new>

namespace {
bool count_allocations{};
std::size_t allocations{};
} // namespace
void *operator new(std::size_t n) {
  if (count_allocations)
    ++allocations;
  if (auto p = std::malloc(n ? n : 1))
    return p;
  throw std::bad_alloc{};
}
void *operator new[](std::size_t n) { return ::operator new(n); }
void operator delete(void *p) noexcept { std::free(p); }
void operator delete[](void *p) noexcept { std::free(p); }
void operator delete(void *p, std::size_t) noexcept { std::free(p); }
void operator delete[](void *p, std::size_t) noexcept { std::free(p); }
using namespace hundun::v04;
namespace {
struct Gas final : portable::GasQueryProvider {
  chemistry::detail::AnalyticIsomerBackend backend;
  unsigned calls{}, fail_after{UINT_MAX};
  const portable::GasIdentity &gas_identity() const noexcept override {
    return backend.gas_identity();
  }
  portable::Status query_gas(const portable::GasQuery &q,
                             portable::GasQueryOutput &out) noexcept override {
    if (++calls > fail_after)
      return portable::Status::unavailable;
    return backend.query_gas(q, out);
  }
};
std::uint64_t id(Int3 p, Int3 n) {
  return std::uint64_t(p.x) +
         std::uint64_t(n.x) * (std::uint64_t(p.y) + std::uint64_t(n.y) * p.z);
}
bool run(int rank, int ranks, const char *path) {
  Gas gas;
  const auto loaded = spray::detail::load_liquid_asset(
      path, UINT64_C(6004043157121730787), gas.gas_identity());
  if (!loaded.available)
    return false;
  const auto &asset = loaded.asset;
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {0, 0, 0};
  mesh.upper = {1, 1, 1};
  mesh.has_exact_cells = true;
  mesh.exact_cells = {8, 6, 4};
  mesh.minimum_spacing = {1e-9, 1e-9, 1e-9};
  mesh.max_growth_ratio = 1;
  mesh.limits = {1000000, 1073741824};
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  if (!CartesianGeometryCompiler::compile(MPI_COMM_WORLD, mesh, {}, geometry,
                                          patch))
    return false;
  const portable::Revision revision{3, UINT64_C(9007199254740993), 1};
  const double duration = 1e-7;
  std::array<double, 2> y{.01, .99}, d{}, h{}, w{};
  portable::GasQuery q{revision,
                       gas.gas_identity().composition_fingerprint,
                       portable::GasStateCoordinates::pressure_enthalpy,
                       100000,
                       101000,
                       0,
                       y.data(),
                       2};
  portable::GasQueryOutput o{{}, d.data(), h.data(), w.data(), 2};
  if (gas.query_gas(q, o) != portable::Status::success)
    return false;
  std::array<std::vector<double>, 5> storage;
  std::array<FieldView, 5> fields{};
  for (std::size_t f = 0; f < 5; ++f) {
    const unsigned components = f == 2 ? 3 : 1;
    const std::size_t sy = patch.cells.x + 4, sz = sy * (patch.cells.y + 4),
                      stride = sz * (patch.cells.z + 4);
    const double value = f == 0   ? 0
                         : f == 1 ? 101000
                         : f == 3 ? .01
                         : f == 4 ? o.sample.density_kg_per_m3
                                  : 0;
    storage[f].assign(components * stride, value);
    auto &v = fields[f];
    v.base = storage[f].data() + 2 + 2 * sy + 2 * sz;
    v.interior = patch.cells;
    v.ghosts = {2, 2, 2};
    v.components = components;
    v.stride_y = sy;
    v.stride_z = sz;
    v.component_stride = stride;
    v.field = f + 1;
    v.revision = 9;
    v.storage_identity = 100 + f;
    v.revision_domain = 1000;
    if (f == 2)
      std::fill(storage[f].begin(), storage[f].begin() + stride, 1.);
  }
  const std::size_t independent = 0;
  detail::ProductParcelGas sampler;
  std::array<ConstFieldView, 1> species{as_const(fields[3])};
  if (!sampler.configure(geometry, patch, {true, true, true},
                         gas.gas_identity().composition_fingerprint,
                         {&independent, 1}, 1) ||
      !sampler.bind(revision, duration, 100000, as_const(fields[0]),
                    as_const(fields[1]), as_const(fields[2]),
                    {species.data(), species.size()}))
    return false;
  detail::ProductParcelGeometry events;
  if (!events.configure(geometry, {true, true, true}))
    return false;
  spray::detail::FilmEnvironmentBridge film(asset, gas, sampler, revision);
  detail::ProductParcelAdvance advance(sampler, events, asset, film);
  SpraySpec spec;
  spec.liquid_fingerprint = asset.content_fingerprint;
  spec.maximum_local_parcels = 32;
  spec.maximum_local_segments = 8192;
  spec.maximum_substep_s = duration;
  spec.minimum_substep_s = 1e-12;
  spec.relative_tolerance = 1e-6;
  auto status = advance.configure(MPI_COMM_WORLD, geometry, patch, spec,
                                  UINT64_C(67108864));
  if (!status) {
    std::cerr << "configure " << status.detail << '\n';
    return false;
  }
  const double mass = 800 * std::acos(-1.) / 6 * 1e-12;
  detail::ProductParcelAdvance::Parcel initial;
  auto &p = initial.parcel;
  p.id = {std::uint64_t(rank + 1), UINT64_C(9007199254740997)};
  p.position_m = {geometry.x().faces().data[patch.begin.x + patch.cells.x] -
                      1e-6,
                  geometry.y().centres().data[patch.begin.y],
                  geometry.z().centres().data[patch.begin.z]};
  p.velocity_m_per_s = {20, 0, 0};
  p.droplet_mass_kg = mass;
  p.droplet_diameter_m = 1e-4;
  p.multiplicity = 1;
  p.temperature_k = 298.15;
  p.liquid_material_fingerprint = asset.pack.material_fingerprint;
  p.owner_global_cell =
      id({patch.begin.x + patch.cells.x - 1, patch.begin.y, patch.begin.z},
         geometry.global_cells());
  initial.breakup_ordinal = UINT64_C(9007199254740999);
  spray::detail::DeterministicInjector injector;
  spray::detail::InjectorSpec injection;
  injection.seed = 71;
  injection.injector_id = 93;
  injection.origin_m = {.05, .05, .05};
  injection.axis = {1, 0, 0};
  injection.injection_speed_m_per_s = 2;
  injection.represented_mass_per_parcel_kg = 10 * mass;
  injection.mass_flow_rate_kg_per_s =
      injection.represented_mass_per_parcel_kg / duration;
  injection.droplet_mass_kg = mass;
  injection.droplet_diameter_m = 1e-4;
  injection.temperature_k = 298.15;
  injection.liquid_material_fingerprint = asset.pack.material_fingerprint;
  if (!injector.reserve(32) || !injector.configure(injection))
    return false;
  auto *injector_ptr = &injector;
  Span<spray::detail::DeterministicInjector *const> injectors{
      &injector_ptr, rank == 0 ? 1U : 0U};
  allocations = 0;
  count_allocations = true;
  status = advance.prepare({&initial, 1}, injectors, revision, duration,
                           as_const(fields[4]), as_const(fields[2]));
  count_allocations = false;
  bool ok = bool(status) && allocations == 0;
  if (!status)
    std::cerr << "prepare rank " << rank << " detail " << status.detail << '\n';
  const auto parcels = advance.parcels();
  const auto exchange = advance.exchange();
  double local[6]{}, global[6]{};
  local[0] = parcels.size;
  for (std::size_t i = 0; i < parcels.size; ++i) {
    const auto &v = parcels.data[i];
    const auto hliq =
        spray::detail::evaluate_liquid_enthalpy(asset, v.parcel.temperature_k);
    const double m = v.parcel.droplet_mass_kg * v.parcel.multiplicity;
    local[1] += m;
    local[2] += m * hliq.liquid_enthalpy_j_per_kg;
    for (double u : v.parcel.velocity_m_per_s)
      local[2] += .5 * m * u * u;
    ok &= v.parcel.droplet_mass_kg < mass && v.parcel.age_s == duration;
    spray::detail::ParcelLocation location;
    ok &= bool(sampler.locate(v.parcel.position_m, revision, location)) &&
          location.owner_rank == rank &&
          location.global_cell == v.parcel.owner_global_cell;
  }
  if (exchange.available) {
    for (std::size_t i = 0; i < exchange.cell_count; ++i) {
      const auto &v = exchange.cells[i];
      local[3] += v.gas.gas_mass_delta_kg;
      local[4] += v.gas.gas_thermochemical_enthalpy_delta_j +
                  v.gas.gas_kinetic_energy_delta_j;
      ok &= v.gas_species_mass_delta_kg[0] == v.gas.gas_mass_delta_kg &&
            v.gas_species_mass_delta_kg[1] == 0;
    }
    local[5] = exchange.external.mass_kg;
  } else
    ok = false;
  MPI_Allreduce(local, global, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  const double injected = 10 * mass,
               initial_energy =
                   ranks * mass * (-100000 + 200) + injected * (-100000 + 2);
  ok &= global[0] == ranks + 1 && global[3] > 0 &&
        std::abs(global[1] + global[3] - (ranks * mass + injected)) < 2e-22 &&
        std::abs(global[2] + global[4] - initial_energy) < 2e-15 &&
        std::abs(global[5] - injected) < 2e-23;
  if (rank == 0)
    ok &=
        injector.committed_state().next_ordinal == 0 && injector.trial_active();
  std::vector<detail::ProductParcelAdvance::Parcel> reference;
  if (parcels.size)
    reference.assign(parcels.data, parcels.data + parcels.size);
  if (injector.trial_active())
    (void)injector.rollback_trial();
  advance.discard();
  // A late PH-provider failure on only one rank must withdraw all candidates
  // and injector trials. The caller's accepted parcel remains untouched.
  gas.calls = 0;
  gas.fail_after = rank == 0 ? 100 : UINT_MAX;
  status = advance.prepare({&initial, 1}, injectors, revision, duration,
                           as_const(fields[4]), as_const(fields[2]));
  ok &= !status && advance.parcels().size == 0 &&
        !advance.exchange().available && !injector.trial_active() &&
        injector.committed_state().next_ordinal == 0;
  gas.calls = 0;
  gas.fail_after = UINT_MAX;
  status = advance.prepare({&initial, 1}, injectors, revision, duration,
                           as_const(fields[4]), as_const(fields[2]));
  ok &= bool(status) && advance.parcels().size == reference.size();
  if (advance.parcels().size == reference.size())
    for (std::size_t i = 0; i < reference.size(); ++i) {
      const auto &a = reference[i], &b = advance.parcels().data[i];
      ok &= a.parcel.id == b.parcel.id &&
            a.parcel.position_m == b.parcel.position_m &&
            a.parcel.droplet_mass_kg == b.parcel.droplet_mass_kg &&
            a.parcel.temperature_k == b.parcel.temperature_k &&
            a.breakup_ordinal == b.breakup_ordinal;
    }
  if (injector.trial_active())
    (void)injector.rollback_trial();
  // Actual TAB threshold -> child construction -> remaining-time A--S and
  // checked migration. The deformation/surface reservoir stays out of gas H.
  detail::ProductParcelAdvance tab_advance(sampler, events, asset, film);
  spec.tab_breakup = true;
  ok &= bool(tab_advance.configure(MPI_COMM_WORLD, geometry, patch, spec,
                                   UINT64_C(67108864)));
  auto triggered = initial;
  triggered.tab_deformation = 1;
  status = tab_advance.prepare({&triggered, 1}, {}, revision, duration,
                               as_const(fields[4]), as_const(fields[2]));
  ok &= bool(status);
  const auto breakup = tab_advance.breakup_energy();
  double tab_local[]{double(tab_advance.parcels().size),
                     double(breakup.event_count),
                     breakup.deformation_consumed_j, breakup.surface_increase_j,
                     breakup.bulk_kinetic_increase_j},
      tab_global[5]{};
  MPI_Allreduce(tab_local, tab_global, 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  ok &= tab_global[0] == 2 * ranks && tab_global[1] == ranks &&
        tab_global[2] > 0 &&
        std::abs(tab_global[3] + tab_global[4] - tab_global[2]) < 1e-15;
  if (!status || tab_global[0] != 2 * ranks)
    std::cerr << "TAB status " << status.detail << " children " << tab_global[0]
              << '\n';
  if (ranks > 1) {
    const double other_duration = rank == 0 ? duration * .5 : duration;
    ok &= !advance.prepare({&initial, 1}, injectors, revision, other_duration,
                           as_const(fields[4]), as_const(fields[2]));
  }
  if (!ok)
    std::cerr << "native spray advance failed rank " << rank << " masses "
              << global[1] << ' ' << global[3] << " count " << global[0]
              << '\n';
  return ok;
}
} // namespace
int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank = 0, ranks = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  int local = argc == 2 && run(rank, ranks, argv[1]), global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Finalize();
  return global ? 0 : 1;
}
