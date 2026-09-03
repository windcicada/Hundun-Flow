// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_app.hpp"

#include "../support/product_fixture.hpp"
#include "core_product_freeze_detail.hpp"
#include "solver_thermophysical_predictor_detail.hpp"

#include <mpi.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <string>
#include <unistd.h>
#include <vector>

namespace {

using namespace hundun::v04;

bool same_u64(std::uint64_t value, MPI_Comm communicator) {
  std::uint64_t minimum = 0U;
  std::uint64_t maximum = 0U;
  return MPI_Allreduce(&value, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                       communicator) == MPI_SUCCESS &&
         MPI_Allreduce(&value, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                       communicator) == MPI_SUCCESS &&
         minimum == maximum;
}

bool expect(bool condition, int rank, const char* description) {
  if (!condition) std::cerr << "rank " << rank << " FAIL: " << description << '\n';
  return condition;
}

bool collective(bool local, MPI_Comm communicator) {
  const int value = local ? 1 : 0;
  int result = 0;
  return MPI_Allreduce(&value, &result, 1, MPI_INT, MPI_MIN, communicator) ==
             MPI_SUCCESS &&
         result != 0;
}

std::uint64_t wire_bits(double value) {
  std::uint64_t bits = 0U;
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

struct CommittedProductBits {
  std::uint64_t time{};
  std::uint64_t dt{};
  std::uint64_t pressure_reference{};
  std::uint64_t step{};
  std::uint64_t controller_state{};
  std::vector<std::uint64_t> field_metadata;
  std::vector<std::uint64_t> field_values;
  std::array<std::vector<std::uint64_t>, 3U> flux_values;
  RevisionToken flux_revision{};
  FaceFluxCertificate flux_certificate{};
};

struct CommittedStateBits {
  std::uint64_t time{};
  std::uint64_t step{};
  std::vector<std::uint64_t> field_metadata;
  std::vector<std::uint64_t> field_values;
};

struct FreshCommittedNumericBits {
  std::uint64_t time{};
  std::uint64_t step{};
  std::vector<std::uint64_t> field_values;
  std::array<std::vector<std::uint64_t>, 3U> flux_values;
};

bool capture_fresh_committed_numeric_bits(ProductDriver& driver,
                                          FreshCommittedNumericBits& out) {
  CommittedOutputSnapshot snapshot;
  ConstFaceFluxView flux;
  if (!driver.committed_output_snapshot(snapshot) || !snapshot.committed ||
      !driver.committed_final_mass_flux_for_test(flux))
    return false;
  FreshCommittedNumericBits candidate;
  candidate.time = wire_bits(snapshot.time);
  candidate.step = snapshot.step;
  for (std::size_t index = 0U; index < snapshot.fields.size; ++index) {
    const ConstFieldView field = snapshot.fields.data[index].values;
    for (std::uint8_t component = 0U; component < field.components;
         ++component)
      for (std::int32_t z = 0; z < field.interior.z; ++z)
        for (std::int32_t y = 0; y < field.interior.y; ++y)
          for (std::int32_t x = 0; x < field.interior.x; ++x)
            candidate.field_values.push_back(
                wire_bits(field.unchecked({x, y, z}, component)));
  }
  const std::array<ConstFaceFieldView, 3U> faces{flux.x, flux.y, flux.z};
  for (std::size_t axis = 0U; axis < faces.size(); ++axis)
    for (std::int32_t z = 0; z < faces[axis].extents.z; ++z)
      for (std::int32_t y = 0; y < faces[axis].extents.y; ++y)
        for (std::int32_t x = 0; x < faces[axis].extents.x; ++x)
          candidate.flux_values[axis].push_back(
              wire_bits(faces[axis].unchecked({x, y, z})));
  out = std::move(candidate);
  return true;
}

bool same_fresh_committed_numeric_bits(
    const FreshCommittedNumericBits& left,
    const FreshCommittedNumericBits& right) {
  return left.time == right.time && left.step == right.step &&
         left.field_values == right.field_values &&
         left.flux_values == right.flux_values;
}

bool capture_committed_state_bits(ProductDriver& driver,
                                  CommittedStateBits& out) {
  CommittedOutputSnapshot snapshot;
  if (!driver.committed_output_snapshot(snapshot) || !snapshot.committed)
    return false;
  CommittedStateBits candidate;
  candidate.time = wire_bits(snapshot.time);
  candidate.step = snapshot.step;
  for (std::size_t index = 0U; index < snapshot.fields.size; ++index) {
    const SnapshotFieldView& field = snapshot.fields.data[index];
    candidate.field_metadata.push_back(field.values.field);
    candidate.field_metadata.push_back(field.values.revision);
    candidate.field_metadata.push_back(field.accepted_revision);
    candidate.field_metadata.push_back(field.values.storage_identity);
    candidate.field_metadata.push_back(field.values.revision_domain);
    for (std::uint8_t component = 0U;
         component < field.values.components; ++component)
      for (std::int32_t z = 0; z < field.values.interior.z; ++z)
        for (std::int32_t y = 0; y < field.values.interior.y; ++y)
          for (std::int32_t x = 0; x < field.values.interior.x; ++x)
            candidate.field_values.push_back(wire_bits(
                field.values.unchecked({x, y, z}, component)));
  }
  out = std::move(candidate);
  return true;
}

bool same_committed_state_bits(const CommittedStateBits& left,
                               const CommittedStateBits& right) {
  return left.time == right.time && left.step == right.step &&
         left.field_metadata == right.field_metadata &&
         left.field_values == right.field_values;
}

bool capture_committed_product_bits(ProductDriver& driver,
                                    CommittedProductBits& out) {
  RestartSnapshot snapshot;
  const Status status = driver.committed_restart_snapshot(snapshot);
  if (!status) return false;
  CommittedProductBits candidate;
  candidate.time = wire_bits(snapshot.time);
  candidate.dt = wire_bits(snapshot.dt);
  candidate.pressure_reference = wire_bits(snapshot.pressure_reference);
  candidate.step = snapshot.step;
  candidate.controller_state = snapshot.controller_state;
  for (std::size_t index = 0U; index < snapshot.fields.size; ++index) {
    const RestartFieldView& field = snapshot.fields.data[index];
    candidate.field_metadata.push_back(
        static_cast<std::uint64_t>(field.role));
    candidate.field_metadata.push_back(field.values.field);
    candidate.field_metadata.push_back(field.values.revision);
    candidate.field_metadata.push_back(field.values.storage_identity);
    candidate.field_metadata.push_back(field.values.revision_domain);
    for (std::uint8_t component = 0U;
         component < field.values.components; ++component)
      for (std::int32_t z = 0; z < field.values.interior.z; ++z)
        for (std::int32_t y = 0; y < field.values.interior.y; ++y)
          for (std::int32_t x = 0; x < field.values.interior.x; ++x)
            candidate.field_values.push_back(wire_bits(
                field.values.unchecked({x, y, z}, component)));
  }
  const std::array<ConstFaceFieldView, 3U> flux{{
      snapshot.final_mass_flux.x, snapshot.final_mass_flux.y,
      snapshot.final_mass_flux.z}};
  for (std::size_t axis = 0U; axis < flux.size(); ++axis)
    for (std::int32_t z = 0; z < flux[axis].extents.z; ++z)
      for (std::int32_t y = 0; y < flux[axis].extents.y; ++y)
        for (std::int32_t x = 0; x < flux[axis].extents.x; ++x)
          candidate.flux_values[axis].push_back(
              wire_bits(flux[axis].unchecked({x, y, z})));
  candidate.flux_revision = snapshot.final_mass_flux.revision;
  candidate.flux_certificate = snapshot.final_mass_flux.certificate;
  out = std::move(candidate);
  return true;
}

bool same_committed_product_bits(const CommittedProductBits& left,
                                 const CommittedProductBits& right) {
  return left.time == right.time && left.dt == right.dt &&
         left.pressure_reference == right.pressure_reference &&
         left.step == right.step &&
         left.controller_state == right.controller_state &&
         left.field_metadata == right.field_metadata &&
         left.field_values == right.field_values &&
         left.flux_values == right.flux_values &&
         left.flux_revision == right.flux_revision &&
         left.flux_certificate == right.flux_certificate;
}

ValidatedModel multispecies_open_model(bool immersed = false) {
  ValidatedModel model = test::product_model(
      immersed ? Int3{16, 16, 16} : Int3{8, 7, 6});
  model.time.initial_dt = immersed ? 2.5e-6 : 1.0e-8;
  model.time.minimum_dt = 1.0e-12;
  model.pressure_reference = PressureReferenceKind::boundary_absolute;
  const SpeciesThermophysicalSpec prototype =
      model.thermophysics.species.front();
  model.thermophysics.species.clear();
  for (std::size_t index = 0U; index < 3U; ++index) {
    SpeciesThermophysicalSpec species = prototype;
    species.stable_name = "sp" + std::to_string(index);
    species.molecular_weight = 24.0 + 8.0 * static_cast<double>(index);
    model.thermophysics.species.push_back(std::move(species));
  }
  model.transported_scalars = {
      {"sp0", TransportedScalarRole::species, 1.0, 1.0},
      {"sp1", TransportedScalarRole::species, 1.0, 1.0},
      {"tracer", TransportedScalarRole::passive_scalar, 1.0, 1.0}};
  for (BoundaryFaceSpec& face : model.boundaries) {
    face.flow_kind = BoundaryKind::symmetry;
    face.thermal_kind = BoundaryKind::none;
    face.pressure = 98000.0;
    face.temperature = 315.0;
    face.backflow_temperature = 315.0;
    face.scalars = {
        {"sp0", ScalarBoundaryKind::zero_gradient, 0.0,
         ScalarBoundaryKind::zero_gradient, 0.0},
        {"sp1", ScalarBoundaryKind::zero_gradient, 0.0,
         ScalarBoundaryKind::zero_gradient, 0.0},
        {"tracer", ScalarBoundaryKind::zero_gradient, 0.0,
         ScalarBoundaryKind::zero_gradient, 0.0}};
  }
  BoundaryFaceSpec& inlet = model.boundaries[0U];
  inlet.flow_kind = BoundaryKind::velocity_inlet;
  inlet.velocity = {1.0e-4, 0.0, 0.0};
  inlet.direction = {1.0, 0.0, 0.0};
  inlet.scalars = {
      {"sp0", ScalarBoundaryKind::dirichlet, 0.1001,
       ScalarBoundaryKind::dirichlet, 0.1001},
      {"sp1", ScalarBoundaryKind::dirichlet, 0.2001,
       ScalarBoundaryKind::dirichlet, 0.2001},
      {"tracer", ScalarBoundaryKind::dirichlet, 0.0501,
       ScalarBoundaryKind::dirichlet, 0.0501}};
  BoundaryFaceSpec& outlet = model.boundaries[1U];
  outlet.flow_kind = BoundaryKind::pressure_outlet;
  outlet.allow_backflow = true;
  outlet.backflow_velocity = {-1.0e-4, 0.0, 0.0};
  outlet.scalars = {
      {"sp0", ScalarBoundaryKind::zero_gradient, 0.0,
       ScalarBoundaryKind::dirichlet, 0.1001},
      {"sp1", ScalarBoundaryKind::zero_gradient, 0.0,
       ScalarBoundaryKind::dirichlet, 0.2001},
      {"tracer", ScalarBoundaryKind::zero_gradient, 0.0,
       ScalarBoundaryKind::dirichlet, 0.0501}};
  if (immersed) {
    model.mesh.lower = {-2.0, -2.0, -2.0};
    model.mesh.upper = {2.0, 2.0, 2.0};
    model.mesh.minimum_spacing = {0.25, 0.25, 0.25};
    model.immersed_boundary = ImmersedBoundarySpec{
        "cylinder_ascii.stl", ImmersedFluidSide::outside};
  }
  return model;
}

bool run_fresh_initialize_invalid_input_contract(int rank) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  const int failing_rank = size - 1;
  ValidatedModel model = multispecies_open_model();
  model.fingerprint = UINT64_C(0x18000cb11);
  const std::array<double, 3U> canonical_scalars{{0.19, 0.23, 0.031}};
  const auto canonical_initial = [&]() {
    DriverInitialState initial;
    initial.pressure_reference = 101325.0;
    initial.temperature = 289.0;
    initial.velocity = {0.0125, -0.0025, 0.00125};
    initial.transported_scalars =
        {canonical_scalars.data(), canonical_scalars.size()};
    initial.start_time = 0.0;
    return initial;
  };
  const auto create_driver = [&](ProductDriver& driver) {
    CompiledCasePlan plan;
    Status status =
        ProductCompiler::compile(MPI_COMM_WORLD, model, {}, plan);
    return status ? ProductDriver::create(MPI_COMM_WORLD, std::move(plan),
                                          driver)
                  : status;
  };

  ProductDriver reference;
  Status status = create_driver(reference);
  DriverInitialState reference_initial = canonical_initial();
  if (status) status = reference.initialize(reference_initial);
  FreshCommittedNumericBits reference_bits;
  const bool captured_reference =
      status && capture_fresh_committed_numeric_bits(reference, reference_bits);
  bool passed = status && captured_reference;

  for (std::size_t selected = 0U; selected < 2U; ++selected) {
    ProductDriver driver;
    Status setup = create_driver(driver);
    std::array<double, 3U> poisoned_scalars = canonical_scalars;
    if (rank == failing_rank)
      poisoned_scalars[1U] =
          selected == 0U ? std::numeric_limits<double>::quiet_NaN() : 1.25;
    DriverInitialState poisoned = canonical_initial();
    poisoned.transported_scalars =
        {poisoned_scalars.data(), poisoned_scalars.size()};
    const Status rejected = setup ? driver.initialize(poisoned) : setup;
    const std::uint64_t rejection_wire =
        (static_cast<std::uint64_t>(rejected.code) << 32U) | rejected.detail;
    const bool collective_rejection = same_u64(rejection_wire, MPI_COMM_WORLD);
    const bool untouched = !driver.initialized();

    DriverInitialState retry_initial = canonical_initial();
    const Status retry = setup ? driver.initialize(retry_initial) : setup;
    FreshCommittedNumericBits retry_bits;
    const bool captured_retry =
        retry && capture_fresh_committed_numeric_bits(driver, retry_bits);
    const bool local = setup && !rejected &&
                       rejected.code == StatusCode::numerical_failure &&
                       collective_rejection && untouched && retry &&
                       captured_retry && captured_reference &&
                       same_fresh_committed_numeric_bits(reference_bits,
                                                         retry_bits);
    if (!local)
      std::cerr << "rank " << rank << " fresh-init-invalid=" << selected
                << " setup/rejected/retry="
                << static_cast<unsigned>(setup.code) << '/' << setup.detail
                << ':' << static_cast<unsigned>(rejected.code) << '/'
                << rejected.detail << ':'
                << static_cast<unsigned>(retry.code) << '/' << retry.detail
                << " collective/untouched/captured/equal="
                << collective_rejection << '/' << untouched << '/'
                << captured_retry << '/'
                << (captured_retry && captured_reference &&
                    same_fresh_committed_numeric_bits(reference_bits,
                                                      retry_bits))
                << '\n';
    passed = passed && local;
  }
  if (size > 1) {
    for (std::size_t selected = 0U; selected < 5U; ++selected) {
      ProductDriver driver;
      Status setup = create_driver(driver);
      std::array<double, 3U> differing_scalars = canonical_scalars;
      DriverInitialState differing = canonical_initial();
      differing.transported_scalars =
          {differing_scalars.data(), differing_scalars.size()};
      if (rank == failing_rank) {
        if (selected == 0U)
          differing.pressure_reference += 17.0;
        else if (selected == 1U)
          differing.temperature += 3.0;
        else if (selected == 2U)
          differing.velocity.x += 0.004;
        else if (selected == 3U)
          differing_scalars[0U] += 0.013;
        else
          differing.start_time += 1.0e-7;
      }
      const Status rejected = setup ? driver.initialize(differing) : setup;
      const std::uint64_t rejection_wire =
          (static_cast<std::uint64_t>(rejected.code) << 32U) |
          rejected.detail;
      const bool collective_rejection =
          same_u64(rejection_wire, MPI_COMM_WORLD);
      const bool untouched = !driver.initialized();

      DriverInitialState retry_initial = canonical_initial();
      const Status retry = setup ? driver.initialize(retry_initial) : setup;
      FreshCommittedNumericBits retry_bits;
      const bool captured_retry =
          retry && capture_fresh_committed_numeric_bits(driver, retry_bits);
      const bool local = setup && !rejected &&
                         rejected.code == StatusCode::invalid_plan &&
                         collective_rejection && untouched && retry &&
                         captured_retry && captured_reference &&
                         same_fresh_committed_numeric_bits(reference_bits,
                                                           retry_bits);
      if (!local)
        std::cerr << "rank " << rank << " fresh-init-mismatch=" << selected
                  << " setup/rejected/retry="
                  << static_cast<unsigned>(setup.code) << '/' << setup.detail
                  << ':' << static_cast<unsigned>(rejected.code) << '/'
                  << rejected.detail << ':'
                  << static_cast<unsigned>(retry.code) << '/' << retry.detail
                  << " collective/untouched/captured/equal="
                  << collective_rejection << '/' << untouched << '/'
                  << captured_retry << '/'
                  << (captured_retry && captured_reference &&
                      same_fresh_committed_numeric_bits(reference_bits,
                                                        retry_bits))
                  << '\n';
      passed = passed && local;
    }
  }
  return passed;
}

bool run_restart_restore_allocation_contract(int rank) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  const int failing_rank = size - 1;
  ValidatedModel model = multispecies_open_model();
  model.fingerprint = UINT64_C(0x18000cb12);
  model.pressure_reference = PressureReferenceKind::closed_mass;
  for (BoundaryFaceSpec& face : model.boundaries) {
    face.flow_kind = BoundaryKind::symmetry;
    face.allow_backflow = false;
    for (ScalarBoundarySpec& scalar : face.scalars) {
      scalar.kind = ScalarBoundaryKind::zero_gradient;
      scalar.backflow_kind = ScalarBoundaryKind::zero_gradient;
    }
  }
  const std::array<double, 2U> composition{{0.10, 0.20}};
  ThermodynamicsPlan oracle;
  Status status = ThermodynamicsPlan::compile(
      model.thermophysics,
      {model.transported_scalars.data(), model.transported_scalars.size()},
      oracle);
  double enthalpy = 0.0;
  double heat_capacity = 0.0;
  double gas_constant = 0.0;
  if (status)
    status = oracle.mixture_enthalpy(
        315.0, {composition.data(), composition.size()}, enthalpy,
        heat_capacity, gas_constant);
  const auto create_driver = [&](ProductDriver& driver) {
    CompiledCasePlan plan;
    Status created =
        ProductCompiler::compile(MPI_COMM_WORLD, model, {}, plan);
    return created ? ProductDriver::create(
                         MPI_COMM_WORLD, std::move(plan), driver)
                   : created;
  };

  ProductDriver reference;
  if (status) status = create_driver(reference);
  RestartExpected expected;
  if (status) status = reference.restart_expected(expected);
  RestartImage image;
  if (status) {
    image.global_cells = expected.global_cells;
    image.patch = expected.target_patch;
    image.plan = expected.plan;
    image.schema = expected.schema;
    image.geometry = expected.geometry;
    image.time = 1.0e-8;
    image.dt = 1.0e-8;
    image.pressure_reference = 101325.0;
    image.step = 1U;
    image.controller_state = 1U;
    image.backward_euler_recovery = true;
    const std::size_t cells =
        static_cast<std::size_t>(image.patch.cells.x) *
        static_cast<std::size_t>(image.patch.cells.y) *
        static_cast<std::size_t>(image.patch.cells.z);
    std::size_t species = 0U;
    for (std::size_t index = 0U; index < expected.fields.size; ++index) {
      const RestartExpectedField descriptor = expected.fields.data[index];
      RestartImageField field;
      field.role = descriptor.role;
      field.field = descriptor.field;
      field.components = descriptor.components;
      field.values.assign(cells * descriptor.components, 0.0);
      if (descriptor.role == RestartFieldRole::velocity) {
        for (std::size_t cell = 0U; cell < cells; ++cell)
          field.values[cell * descriptor.components] = 1.0e-4;
      } else if (descriptor.role == RestartFieldRole::enthalpy) {
        std::fill(field.values.begin(), field.values.end(), enthalpy);
      } else if (descriptor.role ==
                 RestartFieldRole::independent_species) {
        const double value = composition[species++];
        std::fill(field.values.begin(), field.values.end(), value);
      } else if (descriptor.role == RestartFieldRole::transported_scalar) {
        std::fill(field.values.begin(), field.values.end(), 0.05);
      }
      image.fields.push_back(std::move(field));
    }
    const Int3 cells3 = image.patch.cells;
    image.final_mass_flux[0U].assign(
        static_cast<std::size_t>(cells3.x + 1) * cells3.y * cells3.z, 0.0);
    image.final_mass_flux[1U].assign(
        static_cast<std::size_t>(cells3.x) * (cells3.y + 1) * cells3.z,
        0.0);
    image.final_mass_flux[2U].assign(
        static_cast<std::size_t>(cells3.x) * cells3.y * (cells3.z + 1),
        0.0);
  }
  if (status) status = reference.initialize_restart(image);
  FreshCommittedNumericBits reference_bits;
  const bool captured_reference =
      status && capture_fresh_committed_numeric_bits(reference, reference_bits);
  bool passed = status && captured_reference;

  const std::array<detail::RestartRestoreAllocationPoint, 2U> points{{
      detail::RestartRestoreAllocationPoint::thermophysical_staging,
      detail::RestartRestoreAllocationPoint::scalar_seed,
  }};
  for (std::size_t selected = 0U; selected < points.size(); ++selected) {
    ProductDriver driver;
    const Status setup = create_driver(driver);
    detail::arm_restart_restore_allocation_failure_once_for_test(
        failing_rank, points[selected]);
    const Status rejected = setup ? driver.initialize_restart(image) : setup;
    const int lowest =
        detail::restart_restore_allocation_lowest_failing_rank_for_test();
    const std::uint64_t rejection_wire =
        (static_cast<std::uint64_t>(rejected.code) << 32U) | rejected.detail;
    const bool collective_rejection = same_u64(rejection_wire, MPI_COMM_WORLD);
    const bool collective_lowest =
        same_u64(static_cast<std::uint64_t>(lowest + 1), MPI_COMM_WORLD);
    const bool untouched = !driver.initialized();
    detail::clear_restart_restore_allocation_failure_for_test();

    const Status retry = setup ? driver.initialize_restart(image) : setup;
    FreshCommittedNumericBits retry_bits;
    const bool captured_retry =
        retry && capture_fresh_committed_numeric_bits(driver, retry_bits);
    const bool local =
        setup && !rejected &&
        rejected.code == StatusCode::allocation_failure &&
        collective_rejection && collective_lowest && lowest == failing_rank &&
        untouched && retry && captured_retry && captured_reference &&
        same_fresh_committed_numeric_bits(reference_bits, retry_bits);
    if (!local)
      std::cerr << "rank " << rank << " restart-allocation=" << selected
                << " setup/rejected/retry="
                << static_cast<unsigned>(setup.code) << '/' << setup.detail
                << ':' << static_cast<unsigned>(rejected.code) << '/'
                << rejected.detail << ':'
                << static_cast<unsigned>(retry.code) << '/' << retry.detail
                << " collective/lowest/untouched/captured/equal="
                << collective_rejection << '/' << lowest << '/' << untouched
                << '/' << captured_retry << '/'
                << (captured_retry && captured_reference &&
                    same_fresh_committed_numeric_bits(reference_bits,
                                                      retry_bits))
                << '\n';
    passed = passed && local;
  }
  detail::clear_restart_restore_allocation_failure_for_test();
  return passed;
}

bool valid_candidate_storage_lineage(const CompiledCasePlan& plan, int rank) {
  detail::PressureEnergyCandidateStorageDiagnostic diagnostic;
  bool valid = detail::pressure_energy_candidate_storage_diagnostic_for_test(
                   diagnostic) &&
               diagnostic.valid && diagnostic.plan == plan.fingerprint() &&
               diagnostic.lineage_fingerprint != 0U &&
               same_u64(diagnostic.lineage_fingerprint, MPI_COMM_WORLD) &&
               diagnostic.coupled_state_halo != 0U &&
               diagnostic.candidate_state_halo != 0U &&
               diagnostic.coupled_state_halo !=
                   diagnostic.candidate_state_halo &&
               diagnostic.candidate_finalizer_state_halo != 0U &&
               diagnostic.candidate_finalizer_state_halo !=
                   diagnostic.candidate_state_halo &&
               diagnostic.candidate_finalizer_state_halo !=
                   diagnostic.coupled_state_halo &&
               diagnostic.correction_halo != 0U &&
               diagnostic.candidate_correction_halo != 0U &&
               diagnostic.correction_halo !=
                   diagnostic.candidate_correction_halo &&
               diagnostic.workspace_flux_capacity.replicas == 6U &&
               diagnostic.workspace_flux_capacity.directional_blocks == 18U &&
               diagnostic.final_flux_capacity.replicas == 3U &&
               diagnostic.final_flux_capacity.directional_blocks == 9U &&
               diagnostic.execution_graph != 0U &&
               diagnostic.corrector_one_resource_contract
                       .merged_halo_messages > 0U &&
               diagnostic.corrector_one_resource_contract
                       .merged_halo_bytes > 0U &&
               diagnostic.corrector_two_resource_contract
                       .merged_halo_messages ==
                   diagnostic.corrector_one_resource_contract
                       .merged_halo_messages &&
               diagnostic.corrector_two_resource_contract
                       .merged_halo_bytes ==
                   diagnostic.corrector_one_resource_contract
                       .merged_halo_bytes;
  const auto disjoint = [](std::uintptr_t left_begin,
                           std::uintptr_t left_end,
                           std::uintptr_t right_begin,
                           std::uintptr_t right_end) {
    return left_begin != 0U && left_begin < left_end && right_begin != 0U &&
           right_begin < right_end &&
           (left_end <= right_begin || right_end <= left_begin);
  };
  const auto same_halo_capacity = [](const HaloPlanStats& left,
                                     const HaloPlanStats& right) {
    return left.transport_peer_count == right.transport_peer_count &&
           left.local_peer_count == right.local_peer_count &&
           left.persistent_request_count == right.persistent_request_count &&
           left.send_capacity_doubles == right.send_capacity_doubles &&
           left.receive_capacity_doubles == right.receive_capacity_doubles &&
           left.maximum_messages_per_exchange ==
               right.maximum_messages_per_exchange &&
           left.maximum_bytes_per_exchange ==
               right.maximum_bytes_per_exchange &&
           left.maximum_tag == right.maximum_tag;
  };
  valid &= same_halo_capacity(diagnostic.coupled_state_halo_plan,
                              diagnostic.candidate_state_halo_plan) &&
           same_halo_capacity(diagnostic.correction_halo_plan,
                              diagnostic.candidate_correction_halo_plan);
  for (std::size_t left = 0U; left < diagnostic.fields.size(); ++left) {
    const auto& field = diagnostic.fields[left];
    valid &= field.candidate_field != field.trial_field &&
             field.candidate_revision != field.trial_revision &&
             field.candidate_storage == field.trial_storage &&
             field.candidate_revision_domain == field.trial_revision_domain &&
             disjoint(field.candidate_begin, field.candidate_end,
                      field.trial_begin, field.trial_end);
    for (std::size_t right = left + 1U;
         right < diagnostic.fields.size(); ++right)
      valid &= disjoint(field.candidate_begin, field.candidate_end,
                        diagnostic.fields[right].candidate_begin,
                        diagnostic.fields[right].candidate_end);
  }
  for (std::size_t left = 0U; left < diagnostic.scratch_fields.size(); ++left) {
    const auto& field = diagnostic.scratch_fields[left];
    valid &= field.candidate_field != field.trial_field &&
             field.candidate_revision != field.trial_revision &&
             field.candidate_storage == field.trial_storage &&
             field.candidate_revision_domain == field.trial_revision_domain &&
             disjoint(field.candidate_begin, field.candidate_end,
                      field.trial_begin, field.trial_end);
    for (const auto& primitive : diagnostic.fields)
      valid &= disjoint(field.candidate_begin, field.candidate_end,
                        primitive.candidate_begin, primitive.candidate_end);
    for (std::size_t right = left + 1U;
         right < diagnostic.scratch_fields.size(); ++right)
      valid &= disjoint(field.candidate_begin, field.candidate_end,
                        diagnostic.scratch_fields[right].candidate_begin,
                        diagnostic.scratch_fields[right].candidate_end);
  }
  for (std::size_t left = 0U; left < diagnostic.flux_replicas.size(); ++left)
    for (std::size_t right = left + 1U;
         right < diagnostic.flux_replicas.size(); ++right)
      valid &= disjoint(
          diagnostic.flux_replicas[left].replica_begin,
          diagnostic.flux_replicas[left].replica_end,
          diagnostic.flux_replicas[right].replica_begin,
          diagnostic.flux_replicas[right].replica_end);
  return expect(valid, rank,
                "candidate state/material/halo/flux lineage is rank-consistent");
}

bool run_candidate_storage_lineage_only(int rank) {
  ValidatedModel model = test::product_model({17, 11, 7});
  model.fingerprint = UINT64_C(0x18000c501);
  CompiledCasePlan plan;
  Status status = ProductCompiler::compile(MPI_COMM_WORLD, model, {}, plan);
  bool passed = status && valid_candidate_storage_lineage(plan, rank) &&
                same_u64(plan.fingerprint(), MPI_COMM_WORLD);
  detail::PressureEnergyCandidateStorageDiagnostic diagnostic;
  if (status)
    passed &= detail::pressure_energy_candidate_storage_diagnostic_for_test(
        diagnostic);
  ProductDriver driver;
  if (status)
    status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
  DriverInitialState initial;
  initial.pressure_reference = 101325.0;
  initial.temperature = 300.0;
  if (status) status = driver.initialize(initial);
  DriverStepReport step;
  if (status)
    status = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, step);
  RestartSnapshot snapshot;
  if (status) status = driver.committed_restart_snapshot(snapshot);
  if (status) {
    const auto& candidate_flux = diagnostic.flux_replicas[2U];
    const std::array<ConstFaceFieldView, 3U> final_flux{{
        snapshot.final_mass_flux.x, snapshot.final_mass_flux.y,
        snapshot.final_mass_flux.z}};
    for (std::size_t axis = 0U; axis < final_flux.size(); ++axis)
      passed &= final_flux[axis].base != nullptr &&
                reinterpret_cast<std::uintptr_t>(final_flux[axis].base) !=
                    candidate_flux.bases[axis] &&
                final_flux[axis].storage_identity != candidate_flux.storage;
  }
  passed &= status && step.accepted;
  return collective(passed, MPI_COMM_WORLD);
}

bool run_localized_immersed_compile_only(int rank, bool inject_bind_failure) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 2 && size != 4) return true;

  // Place the immersed cube entirely in the low-x half of the canonical
  // decomposition. At least one rank therefore owns no local IBM stencil,
  // while every rank must bind the collective donor reach and enter the
  // candidate-boundary collectives in the same order.
  ValidatedModel model = test::product_model({32, 16, 16});
  model.fingerprint = UINT64_C(0x18000c502);
  model.mesh.lower = {-2.0, -2.0, -2.0};
  model.mesh.upper = {6.0, 2.0, 2.0};
  model.mesh.minimum_spacing = {0.25, 0.25, 0.25};
  model.mesh.limits.max_global_cells = UINT64_C(32) * 16U * 16U;
  model.pressure_reference = PressureReferenceKind::boundary_absolute;
  model.boundaries[0U].flow_kind = BoundaryKind::velocity_inlet;
  model.boundaries[0U].velocity = {1.0, 0.0, 0.0};
  model.boundaries[0U].direction = {1.0, 0.0, 0.0};
  model.boundaries[0U].temperature = 300.0;
  model.boundaries[1U].flow_kind = BoundaryKind::pressure_outlet;
  model.boundaries[1U].pressure = 101325.0;
  model.boundaries[1U].allow_backflow = false;
  model.boundaries[1U].backflow_temperature = 300.0;
  model.immersed_boundary = ImmersedBoundarySpec{
      "cylinder_ascii.stl", ImmersedFluidSide::outside};
  const std::filesystem::path data_root =
      std::filesystem::path{HUNDUN_V04_SOURCE_ROOT} / "tests" / "data";

  if (inject_bind_failure)
    detail::arm_product_piso_bind_failure_once_for_test(size - 1);
  CompiledCasePlan plan;
  const Status status =
      ProductCompiler::compile(MPI_COMM_WORLD, model, data_root, plan);
  detail::clear_product_piso_bind_failure_for_test();
  if (inject_bind_failure) {
    const std::uint64_t status_identity =
        (static_cast<std::uint64_t>(status.code) << 32U) | status.detail;
    const bool passed =
        status.code == StatusCode::invalid_plan && status.detail == 1503U &&
        plan.fingerprint() == 0U &&
        same_u64(status_identity, MPI_COMM_WORLD);
    if (!passed && rank == 0)
      std::cerr << "localized-immersed-bind-fail-close status="
                << static_cast<unsigned>(status.code) << '/' << status.detail
                << '\n';
    return collective(passed, MPI_COMM_WORLD);
  }
  detail::PressureEnergyCandidateStorageDiagnostic diagnostic;
  const bool observed =
      status && detail::pressure_energy_candidate_storage_diagnostic_for_test(
                    diagnostic);
  const int local_empty =
      observed && diagnostic.local_ibm_reconstruction_reach == 0U ? 1 : 0;
  const int local_active =
      observed && diagnostic.local_ibm_reconstruction_reach > 0U ? 1 : 0;
  int any_empty = 0;
  int any_active = 0;
  unsigned local_reach =
      observed ? diagnostic.local_ibm_reconstruction_reach : 0U;
  unsigned maximum_local_reach = 0U;
  MPI_Allreduce(&local_empty, &any_empty, 1, MPI_INT, MPI_MAX,
                MPI_COMM_WORLD);
  MPI_Allreduce(&local_active, &any_active, 1, MPI_INT, MPI_MAX,
                MPI_COMM_WORLD);
  MPI_Allreduce(&local_reach, &maximum_local_reach, 1, MPI_UNSIGNED, MPI_MAX,
                MPI_COMM_WORLD);
  const bool passed =
      status && observed && diagnostic.valid && plan.fingerprint() != 0U &&
      diagnostic.candidate_pressure_donor_plan != 0U && any_empty != 0 &&
      any_active != 0 && diagnostic.candidate_pressure_donor_reach > 0U &&
      diagnostic.candidate_pressure_donor_reach == maximum_local_reach &&
      same_u64(diagnostic.candidate_pressure_donor_reach, MPI_COMM_WORLD) &&
      same_u64(plan.fingerprint(), MPI_COMM_WORLD);
  if (!passed && rank == 0)
    std::cerr << "localized-immersed-compile status="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << " observed/empty/active=" << observed << '/' << any_empty
              << '/' << any_active << '\n';
  return collective(passed, MPI_COMM_WORLD);
}

bool valid_momentum_solve(
    const MomentumPredictorSolveReport& report,
    const std::array<std::uint32_t, 3U>& iterations,
    const std::array<std::uint64_t, 3U>& operator_applies,
    const std::array<std::uint64_t, 3U>& preconditioner_applies,
    const std::array<std::uint64_t, 3U>& reduction_calls) {
  if (report.solve_calls != 3U) return false;
  for (std::size_t component = 0U; component < report.components.size();
       ++component) {
    const LinearSolveResult& solve = report.components[component];
    const bool accepted =
        solve.termination == LinearTermination::converged ||
        solve.termination == LinearTermination::zero_rhs;
    if (!solve.status || !accepted ||
        solve.iterations != iterations[component] ||
        solve.operator_applies != operator_applies[component] ||
        solve.preconditioner_applies !=
            preconditioner_applies[component] ||
        solve.reduction_calls != reduction_calls[component] ||
        !std::isfinite(solve.initial_true_residual) ||
        !std::isfinite(solve.final_true_residual) ||
        solve.final_true_residual > solve.initial_true_residual ||
        solve.recycle_offered_directions != 0U ||
        solve.recycle_retained_directions != 0U ||
        solve.recycle_operator_applies != 0U ||
        solve.recycle_reduction_calls != 0U ||
        solve.recycle_projection_attempted ||
        solve.recycle_projection_accepted ||
        solve.recycle_cycle_corrections != 0U)
      return false;
  }
  return true;
}

bool run_temporal_method_fallback_product(int rank) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  // Keep one global physical/discrete fixture and change decomposition only.
  ValidatedModel model = test::product_model({8, 8, 8});
  model.fingerprint = UINT64_C(0x18000f001);
  model.time.initial_dt = 1.0e-3;
  model.time.minimum_dt = 1.0e-3;
  model.time.maximum_dt = 1.0e-3;
  model.time.maximum_growth = 1.0;
  model.time.maximum_retries = 1U;
  CompiledCasePlan plan;
  Status status = ProductCompiler::compile(MPI_COMM_WORLD, model, {}, plan);
  ProductDriver driver;
  if (status)
    status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
  DriverInitialState initial;
  initial.pressure_reference = 101325.0;
  initial.temperature = 300.0;
  if (status) status = driver.initialize(initial);
  DriverStepReport first;
  if (status)
    status = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, first);
  if (status)
    detail::arm_low_bdf_source_base_failure_once_for_test();
  DriverStepReport second;
  if (status)
    status = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, second);
  detail::clear_low_bdf_source_base_failure_for_test();

  const bool passed =
      status && first.accepted && first.attempts == 1U &&
      first.proposal.bdf.order == 1U && first.effective_bdf.order == 1U &&
      first.thermophysical_predictor_calls == 1U &&
      !first.temporal_method_fallback &&
      first.piso.pressure_solve_calls == 2U &&
      second.accepted && second.attempts == 1U &&
      second.proposal.bdf.order == 2U && second.effective_bdf.order == 1U &&
      second.thermophysical_predictor_calls == 2U &&
      second.temporal_method_fallback &&
      second.failure.code == StatusCode::ok &&
      second.piso.pressure_solve_calls == 2U && second.accepted_step == 2U &&
      std::abs(second.accepted_time - 2.0e-3) <= 1.0e-16;
  if (!passed)
    std::cerr << "rank " << rank << " temporal-fallback="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << " first=" << first.accepted << '/'
              << static_cast<unsigned>(first.proposal.bdf.order) << '/'
              << static_cast<unsigned>(first.effective_bdf.order) << '/'
              << static_cast<unsigned>(first.thermophysical_predictor_calls)
              << '/' << first.temporal_method_fallback << " second="
              << second.accepted << '/' << second.attempts << '/'
              << static_cast<unsigned>(second.proposal.bdf.order) << '/'
              << static_cast<unsigned>(second.effective_bdf.order) << '/'
              << static_cast<unsigned>(second.thermophysical_predictor_calls)
              << '/' << second.temporal_method_fallback << '/'
              << static_cast<unsigned>(second.failure.code) << '/'
              << second.failure.detail << '/'
              << static_cast<unsigned>(second.piso.pressure_solve_calls) << '/'
              << second.accepted_step << '/' << second.accepted_time << '\n';
  return passed;
}

bool run_implicit_enthalpy_product(int rank) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  // One fixed global problem is decomposed across 1/2/4 ranks.  Thirty-two
  // x cells retain at least eight local cells at the largest test size, so
  // the MG hierarchy and the physical endpoint are both rank invariant.
  ValidatedModel model = test::product_model({32, 8, 8});
  model.fingerprint = UINT64_C(0x18000b001);
  model.solver.pressure.krylov_restart = 64U;
  model.time.initial_dt = 0.006;
  model.time.minimum_dt = 0.006;
  model.time.maximum_dt = 0.006;
  model.time.maximum_growth = 1.0;
  model.time.maximum_retries = 1U;
  // Keep a finite physical conductivity, but make advection the isolated
  // endpoint mechanism after Restart reconstructs the checkpoint-consistent
  // thermal-diffusion history.  This mild divergence-free transport then
  // exercises the implicit-upwind positivity path without borrowing the
  // incompatible div(phi)=80 stress from the globalization rejection fixture.
  model.thermophysics.species[0U].conductivity = 1.625e-5;
  constexpr double streamwise_velocity = 3.75e-4;
  constexpr double molecular_weight = 28.96546;
  constexpr double gas_constant =
      kUniversalGasConstant / molecular_weight;
  constexpr double minimum_enthalpy = 3.5 * gas_constant * 200.0;
  constexpr double density = 101325.0 / (gas_constant * 200.0);
  constexpr double x_face_area = (1.0 / 8.0) * (0.5 / 8.0);
  // A tiny divergence-free periodic transport state exercises the implicit
  // enthalpy endpoint without borrowing the deliberately incompatible
  // div(phi)=80 restart used by the globalization rejection fixture.  The
  // two-joule perturbation sits just inside the compiled thermodynamic lower
  // bound, so the high-order endpoint takes the implicit-upwind route while
  // the final EOS/continuity/energy/mass/gauge audit remains stringent.
  CompiledCasePlan plan;
  Status status = ProductCompiler::compile(MPI_COMM_WORLD, model, {}, plan);
  ProductDriver driver;
  if (status)
    status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
  RestartExpected expected;
  if (status) status = driver.restart_expected(expected);

  RestartImage image;
  if (status) {
    image.global_cells = expected.global_cells;
    image.patch = expected.target_patch;
    image.plan = expected.plan;
    image.schema = expected.schema;
    image.geometry = expected.geometry;
    image.time = 0.006;
    image.dt = 0.006;
    image.pressure_reference = 101325.0;
    image.step = 1U;
    image.controller_state = 1U;
    image.backward_euler_recovery = true;
    const std::size_t cells =
        static_cast<std::size_t>(image.patch.cells.x) *
        static_cast<std::size_t>(image.patch.cells.y) *
        static_cast<std::size_t>(image.patch.cells.z);
    image.fields.reserve(expected.fields.size);
    for (std::size_t index = 0U; index < expected.fields.size; ++index) {
      const RestartExpectedField descriptor = expected.fields.data[index];
      RestartImageField field;
      field.role = descriptor.role;
      field.field = descriptor.field;
      field.components = descriptor.components;
      field.values.assign(cells * descriptor.components, 0.0);
      if (descriptor.role == RestartFieldRole::velocity) {
        for (std::size_t cell = 0U; cell < cells; ++cell)
          field.values[cell * 3U] = streamwise_velocity;
      } else if (descriptor.role == RestartFieldRole::enthalpy) {
        std::fill(field.values.begin(), field.values.end(),
                  minimum_enthalpy + 1.0);
        if (image.patch.begin.y == 0 && image.patch.begin.z == 0) {
          for (std::int32_t global_x : {0, 1}) {
            if (global_x < image.patch.begin.x ||
                global_x >= image.patch.begin.x + image.patch.cells.x)
              continue;
            const std::size_t local_x = static_cast<std::size_t>(
                global_x - image.patch.begin.x);
            field.values[local_x] =
                global_x == 0 ? minimum_enthalpy + 1.0e-10
                              : minimum_enthalpy + 2.0;
          }
        }
      }
      image.fields.push_back(std::move(field));
    }
    const Int3 local = image.patch.cells;
    image.final_mass_flux[0U].assign(
        static_cast<std::size_t>(local.x + 1) * local.y * local.z, 0.0);
    image.final_mass_flux[1U].assign(
        static_cast<std::size_t>(local.x) * (local.y + 1) * local.z, 0.0);
    image.final_mass_flux[2U].assign(
        static_cast<std::size_t>(local.x) * local.y * (local.z + 1), 0.0);
    std::fill(image.final_mass_flux[0U].begin(),
              image.final_mass_flux[0U].end(),
              density * streamwise_velocity * x_face_area);
  }
  if (status) status = driver.initialize_restart(image);
  DriverStepReport step;
  if (status)
    status = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, step);
  const LinearSolveResult& solve =
      step.thermophysical_predictor.enthalpy_endpoint;
  const bool passed =
      status && step.accepted && step.attempts == 1U &&
      step.thermophysical_predictor.low_state ==
          ThermophysicalLowStateKind::implicit_upwind &&
      step.thermophysical_predictor.enthalpy_solve_calls == 1U &&
      step.thermophysical_predictor.enthalpy_endpoint_alpha == 1.0 &&
      step.thermophysical_predictor.bdf_endpoint_alpha == 1.0 &&
      step.thermophysical_predictor.source_endpoint_alpha == 1.0 &&
      step.thermophysical_predictor.mass_flux_scale == 1.0 &&
      step.thermophysical_predictor.blocking_collectives == 13U &&
      step.thermophysical_predictor.theta == 0.0 && solve.status &&
      solve.termination == LinearTermination::converged &&
      solve.iterations == 0U && solve.norm_breakdown_restarts == 0U &&
      solve.reduction_calls == 4U && solve.operator_applies == 1U &&
      solve.preconditioner_applies == 0U &&
      std::isfinite(solve.initial_true_residual) &&
      std::isfinite(solve.final_true_residual) &&
      wire_bits(solve.final_true_residual) ==
          wire_bits(solve.initial_true_residual) &&
      step.piso.eos_residual <= model.solver.terminal.eos &&
      step.piso.continuity_residual <= model.solver.terminal.continuity &&
      step.piso.energy_residual <= model.solver.terminal.continuity &&
      step.piso.closed_mass_residual <= model.solver.terminal.closed_mass &&
      step.piso.gauge_residual <= model.solver.terminal.gauge;
  if (rank == 0) {
    std::cout << std::setprecision(17)
              << "implicit-enthalpy-positive size=" << size
              << " h-min/base/peak=" << minimum_enthalpy + 1.0e-10 << '/'
              << minimum_enthalpy + 1.0 << '/'
              << minimum_enthalpy + 2.0
              << " dt/U/phi=" << model.time.initial_dt << '/'
              << streamwise_velocity << '/'
              << density * streamwise_velocity * x_face_area
              << " solve=" << solve.initial_true_residual << '/'
              << solve.final_true_residual << " terminal="
              << step.piso.eos_residual << '/'
              << step.piso.continuity_residual << '/'
              << step.piso.energy_residual << '/'
              << step.piso.closed_mass_residual << '/'
              << step.piso.gauge_residual << '\n';
  }
  if (!passed)
    std::cerr << "rank " << rank << " implicit-enthalpy-product="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << " failed-stage=" << step.failed_stage << " failure="
              << static_cast<unsigned>(step.failure.code) << '/'
              << step.failure.detail
              << " accepted=" << step.accepted
              << " attempts=" << step.attempts
              << " low="
              << static_cast<unsigned>(
                     step.thermophysical_predictor.low_state)
              << " solves="
              << static_cast<unsigned>(
                     step.thermophysical_predictor.enthalpy_solve_calls)
              << " alpha="
              << step.thermophysical_predictor.enthalpy_endpoint_alpha
              << " theta=" << step.thermophysical_predictor.theta
              << " iterations=" << solve.iterations
              << " reductions=" << solve.reduction_calls
              << " operators=" << solve.operator_applies
              << " preconditioners=" << solve.preconditioner_applies
              << " solve_status="
              << static_cast<unsigned>(solve.status.code)
              << " termination="
              << static_cast<unsigned>(solve.termination)
              << " initial=" << solve.initial_true_residual
              << " residual=" << solve.final_true_residual
              << " piso=" << step.piso.eos_residual << '/'
              << step.piso.continuity_residual << '/'
              << step.piso.energy_residual << '/'
              << step.piso.closed_mass_residual << '/'
              << step.piso.gauge_residual << '\n';
  return passed;
}

bool run_pressure_energy_candidate_globalization_red(int rank) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  // Keep one global physical fixture while changing only the decomposition.
  ValidatedModel model = test::product_model({8, 8, 8});
  model.fingerprint = UINT64_C(0x18000b101);
  model.solver.pressure.krylov_restart = 64U;
  // Keep the RED at the first target time.  A rejected full correction cannot
  // be hidden by the normal retry chain because half dt is below the compiled
  // minimum.
  model.time.initial_dt = 0.006;
  model.time.minimum_dt = 0.006;
  model.time.maximum_dt = 0.006;
  model.time.maximum_growth = 1.0;
  model.time.maximum_retries = 1U;
  CompiledCasePlan plan;
  Status status = ProductCompiler::compile(MPI_COMM_WORLD, model, {}, plan);
  ProductDriver driver;
  if (status)
    status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
  RestartExpected expected;
  if (status) status = driver.restart_expected(expected);

  RestartImage image;
  if (status) {
    image.global_cells = expected.global_cells;
    image.patch = expected.target_patch;
    image.plan = expected.plan;
    image.schema = expected.schema;
    image.geometry = expected.geometry;
    image.time = 0.006;
    image.dt = 0.006;
    image.pressure_reference = 101325.0;
    image.step = 1U;
    image.controller_state = 1U;
    image.backward_euler_recovery = true;
    const std::size_t cells =
        static_cast<std::size_t>(image.patch.cells.x) *
        static_cast<std::size_t>(image.patch.cells.y) *
        static_cast<std::size_t>(image.patch.cells.z);
    image.fields.reserve(expected.fields.size);
    for (std::size_t index = 0U; index < expected.fields.size; ++index) {
      const RestartExpectedField descriptor = expected.fields.data[index];
      RestartImageField field;
      field.role = descriptor.role;
      field.field = descriptor.field;
      field.components = descriptor.components;
      field.values.assign(cells * descriptor.components, 0.0);
      if (descriptor.role == RestartFieldRole::enthalpy) {
        std::fill(field.values.begin(), field.values.end(), 300000.0);
        // Anchor the stress to two global cells.  Repeating local offsets
        // 0/1 on every rank changes the physical initial state with the MPI
        // decomposition and cannot certify collective invariance.
        if (image.patch.begin.y == 0 && image.patch.begin.z == 0) {
          for (std::int32_t global_x : {0, 1}) {
            if (global_x < image.patch.begin.x ||
                global_x >= image.patch.begin.x + image.patch.cells.x)
              continue;
            const std::size_t local_x = static_cast<std::size_t>(
                global_x - image.patch.begin.x);
            field.values[local_x] =
                global_x == 0 ? 220000.0 : 380000.0;
          }
        }
      }
      image.fields.push_back(std::move(field));
    }
    const Int3 local = image.patch.cells;
    image.final_mass_flux[0U].assign(
        static_cast<std::size_t>(local.x + 1) * local.y * local.z, 0.0);
    image.final_mass_flux[1U].assign(
        static_cast<std::size_t>(local.x) * (local.y + 1) * local.z, 0.0);
    image.final_mass_flux[2U].assign(
        static_cast<std::size_t>(local.x) * local.y * (local.z + 1), 0.0);
    constexpr double divergence = 80.0;
    constexpr double cell_volume = 1.0 / 512.0;
    if (rank == 0)
      image.final_mass_flux[0U][1U] = divergence * cell_volume;
  }
  if (status) status = driver.initialize_restart(image);
  CommittedProductBits committed_before;
  const bool captured_before =
      status && capture_committed_product_bits(driver, committed_before);
  if (status)
    detail::arm_pressure_energy_candidate_globalization_once_for_test();
  DriverStepReport step;
  if (status)
    status = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, step);
  detail::PressureEnergyCandidateGlobalizationDiagnostic diagnostic;
  const bool observed =
      detail::pressure_energy_candidate_globalization_diagnostic_for_test(
          diagnostic);
  detail::clear_pressure_energy_candidate_globalization_for_test();
  CommittedProductBits committed_after;
  const bool captured_after =
      capture_committed_product_bits(driver, committed_after);
  // The rejected fatal proposal is retired by advancing the controller ticket.
  // That generation is not committed flow state, so prove its exact one-step
  // transition separately before comparing all physical state bit-for-bit.
  const bool fatal_controller_recovery =
      captured_before && captured_after &&
      committed_before.controller_state !=
          std::numeric_limits<std::uint64_t>::max() &&
      committed_after.controller_state ==
          committed_before.controller_state + 1U;
  if (captured_before && captured_after)
    committed_before.controller_state = committed_after.controller_state;
  const bool committed_rollback_exact =
      captured_before && captured_after && fatal_controller_recovery &&
      same_committed_product_bits(committed_before, committed_after);

  const std::size_t reported_sample_count = std::min<std::size_t>(
      diagnostic.sample_count,
      detail::kPressureEnergyCandidateGlobalizationSampleCapacity);
  bool smaller_admissible = false;
  bool admissible_strict_merit_decrease = false;
  bool alpha_sequence = reported_sample_count == diagnostic.sample_count;
  const double baseline_merit =
      std::hypot(diagnostic.baseline_normalized_continuity,
                 diagnostic.baseline_normalized_energy);
  for (std::size_t index = 0U; index < reported_sample_count; ++index) {
    alpha_sequence &= diagnostic.samples[index].alpha ==
                      std::ldexp(1.0, -static_cast<int>(index));
    if (index != 0U) {
      smaller_admissible |= diagnostic.samples[index].admissible;
      if (diagnostic.samples[index].admissible) {
        const double merit = std::hypot(
            diagnostic.samples[index].normalized_continuity,
            diagnostic.samples[index].normalized_energy);
        admissible_strict_merit_decrease |= merit < baseline_merit;
      }
    }
  }
  const auto wire_double = [](double value) {
    std::uint64_t wire = 0U;
    static_assert(sizeof(wire) == sizeof(value));
    std::memcpy(&wire, &value, sizeof(wire));
    return wire;
  };
  const auto mix = [](std::uint64_t hash, std::uint64_t value) {
    hash ^= value;
    hash *= UINT64_C(1099511628211);
    return hash;
  };
  std::uint64_t diagnostic_fingerprint = UINT64_C(1469598103934665603);
  diagnostic_fingerprint = mix(diagnostic_fingerprint, diagnostic.valid);
  diagnostic_fingerprint =
      mix(diagnostic_fingerprint, diagnostic.production_candidate_loop);
  diagnostic_fingerprint =
      mix(diagnostic_fingerprint, diagnostic.baseline_commit);
  diagnostic_fingerprint =
      mix(diagnostic_fingerprint, diagnostic.selection_valid);
  diagnostic_fingerprint =
      mix(diagnostic_fingerprint, diagnostic.replay_valid);
  diagnostic_fingerprint = mix(diagnostic_fingerprint, diagnostic.committed);
  diagnostic_fingerprint =
      mix(diagnostic_fingerprint, diagnostic.attempted_step);
  diagnostic_fingerprint = mix(diagnostic_fingerprint, diagnostic.generation);
  diagnostic_fingerprint = mix(diagnostic_fingerprint, diagnostic.attempt);
  diagnostic_fingerprint = mix(diagnostic_fingerprint, diagnostic.corrector);
  diagnostic_fingerprint =
      mix(diagnostic_fingerprint, diagnostic.sample_count);
  diagnostic_fingerprint =
      mix(diagnostic_fingerprint, diagnostic.first_admissible_sample);
  diagnostic_fingerprint =
      mix(diagnostic_fingerprint,
          wire_double(diagnostic.maximum_absolute_pressure_correction));
  diagnostic_fingerprint =
      mix(diagnostic_fingerprint,
          wire_double(diagnostic.maximum_absolute_enthalpy_correction));
  diagnostic_fingerprint = mix(
      diagnostic_fingerprint,
      wire_double(diagnostic.alpha_zero_maximum_velocity_difference));
  diagnostic_fingerprint = mix(
      diagnostic_fingerprint,
      wire_double(diagnostic.alpha_zero_maximum_density_difference));
  diagnostic_fingerprint = mix(
      diagnostic_fingerprint,
      wire_double(diagnostic.alpha_zero_maximum_temperature_difference));
  diagnostic_fingerprint = mix(
      diagnostic_fingerprint,
      wire_double(diagnostic.alpha_zero_maximum_pressure_difference));
  diagnostic_fingerprint = mix(
      diagnostic_fingerprint,
      wire_double(diagnostic.alpha_zero_maximum_pressure_relative_error));
  diagnostic_fingerprint = mix(
      diagnostic_fingerprint,
      wire_double(diagnostic.alpha_zero_maximum_density_relative_error));
  diagnostic_fingerprint = mix(
      diagnostic_fingerprint,
      wire_double(diagnostic.alpha_zero_maximum_temperature_relative_error));
  diagnostic_fingerprint = mix(
      diagnostic_fingerprint,
      diagnostic.alpha_zero_maximum_density_ulp_difference);
  diagnostic_fingerprint = mix(
      diagnostic_fingerprint,
      diagnostic.alpha_zero_maximum_temperature_ulp_difference);
  diagnostic_fingerprint = mix(
      diagnostic_fingerprint,
      wire_double(diagnostic.alpha_zero_material_oracle_error));
  diagnostic_fingerprint = mix(
      diagnostic_fingerprint,
      wire_double(diagnostic.alpha_zero_gradient_oracle_error));
  diagnostic_fingerprint = mix(
      diagnostic_fingerprint,
      wire_double(
          diagnostic.alpha_zero_effective_viscosity_oracle_error));
  diagnostic_fingerprint = mix(
      diagnostic_fingerprint,
      diagnostic.alpha_zero_oracle_numeric_lineage);
  diagnostic_fingerprint = mix(
      diagnostic_fingerprint, diagnostic.alpha_zero_velocity_bitwise_equal);
  diagnostic_fingerprint = mix(
      diagnostic_fingerprint, diagnostic.alpha_zero_density_bitwise_equal);
  diagnostic_fingerprint = mix(
      diagnostic_fingerprint,
      diagnostic.alpha_zero_temperature_bitwise_equal);
  diagnostic_fingerprint = mix(
      diagnostic_fingerprint,
      wire_double(diagnostic.baseline_normalized_continuity));
  diagnostic_fingerprint = mix(
      diagnostic_fingerprint,
      wire_double(diagnostic.baseline_normalized_energy));
  diagnostic_fingerprint = mix(
      diagnostic_fingerprint,
      wire_double(diagnostic.selected_normalized_continuity));
  diagnostic_fingerprint = mix(
      diagnostic_fingerprint,
      wire_double(diagnostic.selected_normalized_energy));
  diagnostic_fingerprint =
      mix(diagnostic_fingerprint, wire_double(diagnostic.selected_alpha));
  diagnostic_fingerprint = mix(
      diagnostic_fingerprint,
      wire_double(diagnostic.corrector_one_baseline_normalized_continuity));
  diagnostic_fingerprint = mix(
      diagnostic_fingerprint,
      wire_double(diagnostic.corrector_one_baseline_normalized_energy));
  diagnostic_fingerprint = mix(
      diagnostic_fingerprint,
      wire_double(diagnostic.corrector_one_selected_normalized_continuity));
  diagnostic_fingerprint = mix(
      diagnostic_fingerprint,
      wire_double(diagnostic.corrector_one_selected_normalized_energy));
  diagnostic_fingerprint = mix(
      diagnostic_fingerprint,
      wire_double(diagnostic.corrector_one_selected_alpha));
  diagnostic_fingerprint = mix(
      diagnostic_fingerprint,
      wire_double(
          diagnostic.corrector_one_linear_predicted_normalized_continuity));
  diagnostic_fingerprint = mix(
      diagnostic_fingerprint,
      wire_double(diagnostic.corrector_one_linear_predicted_normalized_energy));
  diagnostic_fingerprint = mix(
      diagnostic_fingerprint,
      wire_double(
          diagnostic.corrector_two_linear_predicted_normalized_continuity));
  diagnostic_fingerprint = mix(
      diagnostic_fingerprint,
      wire_double(diagnostic.corrector_two_linear_predicted_normalized_energy));
  diagnostic_fingerprint =
      mix(diagnostic_fingerprint, diagnostic.correction_direction);
  diagnostic_fingerprint =
      mix(diagnostic_fingerprint, diagnostic.selected_state_provenance);
  diagnostic_fingerprint =
      mix(diagnostic_fingerprint, diagnostic.selected_mass_flux_provenance);
  for (std::size_t index = 0U; index < reported_sample_count; ++index) {
    const auto& sample = diagnostic.samples[index];
    diagnostic_fingerprint =
        mix(diagnostic_fingerprint, wire_double(sample.alpha));
    diagnostic_fingerprint = mix(diagnostic_fingerprint, sample.admissible);
    diagnostic_fingerprint =
        mix(diagnostic_fingerprint, sample.first_failing_global_cell);
    diagnostic_fingerprint = mix(
        diagnostic_fingerprint,
        static_cast<std::uint8_t>(sample.first_failure_reason));
    diagnostic_fingerprint = mix(
        diagnostic_fingerprint, wire_double(sample.minimum_absolute_pressure));
    diagnostic_fingerprint = mix(
        diagnostic_fingerprint, wire_double(sample.minimum_temperature));
    diagnostic_fingerprint =
        mix(diagnostic_fingerprint, wire_double(sample.minimum_density));
    diagnostic_fingerprint = mix(
        diagnostic_fingerprint,
        wire_double(sample.normalized_continuity));
    diagnostic_fingerprint =
        mix(diagnostic_fingerprint, wire_double(sample.normalized_energy));
    diagnostic_fingerprint =
        mix(diagnostic_fingerprint, sample.state_and_flux_finite);
  }
  const bool diagnostic_identical =
      same_u64(diagnostic_fingerprint, MPI_COMM_WORLD);
  const detail::PressureEnergyCandidateAlphaDiagnostic& full =
      diagnostic.samples[0U];
  const bool first_admissible_valid =
      diagnostic.first_admissible_sample < reported_sample_count;
  bool first_admissible_is_first = first_admissible_valid;
  for (std::size_t index = 0U;
       index < diagnostic.first_admissible_sample &&
       index < reported_sample_count;
       ++index)
    first_admissible_is_first &= !diagnostic.samples[index].admissible;
  const bool passed =
      !status && status.code == StatusCode::rejected_step && !step.accepted &&
      status.detail == 5792U && step.failed_stage == 53U &&
      step.attempts == 1U && observed && diagnostic.valid &&
      diagnostic.production_candidate_loop && !diagnostic.selection_valid &&
      !diagnostic.replay_valid && !diagnostic.committed &&
      diagnostic_identical && committed_rollback_exact &&
      diagnostic.attempted_step == 2U && diagnostic.generation != 0U &&
      diagnostic.attempt == 0U && diagnostic.corrector == 2U &&
      diagnostic.sample_count ==
          detail::kPressureEnergyCandidateGlobalizationSampleCapacity &&
      alpha_sequence && full.alpha == 1.0 && !full.admissible &&
      !full.state_and_flux_finite &&
      full.first_failing_global_cell ==
          std::numeric_limits<std::uint64_t>::max() &&
      full.first_failure_reason ==
          detail::PressureEnergyCandidateFailureReason::
              production_candidate_evaluation &&
      std::isfinite(diagnostic.maximum_absolute_pressure_correction) &&
      diagnostic.maximum_absolute_pressure_correction > 0.0 &&
      std::isfinite(diagnostic.maximum_absolute_enthalpy_correction) &&
      diagnostic.maximum_absolute_enthalpy_correction > 0.0 &&
      smaller_admissible &&
      diagnostic.first_admissible_sample > 0U &&
      first_admissible_valid &&
      first_admissible_is_first &&
      diagnostic.samples[diagnostic.first_admissible_sample].admissible &&
      !admissible_strict_merit_decrease &&
      diagnostic.selected_alpha == 0.0 &&
      diagnostic.selected_normalized_continuity == 0.0 &&
      diagnostic.selected_normalized_energy == 0.0 &&
      diagnostic.selected_state_provenance == 0U &&
      diagnostic.selected_mass_flux_provenance == 0U &&
      diagnostic.corrector_one_selected_alpha > 0.0 &&
      diagnostic.corrector_one_selected_alpha < 1.0 &&
      diagnostic.corrector_one_selected_normalized_energy <
          diagnostic.corrector_one_baseline_normalized_energy &&
      std::hypot(diagnostic.corrector_one_selected_normalized_continuity,
                 diagnostic.corrector_one_selected_normalized_energy) <
          std::hypot(diagnostic.corrector_one_baseline_normalized_continuity,
                     diagnostic.corrector_one_baseline_normalized_energy) &&
      diagnostic.corrector_one_linear_predicted_normalized_continuity <
          1.0e-8 &&
      diagnostic.corrector_one_linear_predicted_normalized_energy < 1.0e-7 &&
      diagnostic.corrector_two_linear_predicted_normalized_continuity <
          1.0e-8 &&
      diagnostic.corrector_two_linear_predicted_normalized_energy < 1.0e-7 &&
      diagnostic.alpha_zero_velocity_bitwise_equal &&
      diagnostic.alpha_zero_density_bitwise_equal &&
      diagnostic.alpha_zero_temperature_bitwise_equal &&
      diagnostic.alpha_zero_maximum_velocity_difference == 0.0 &&
      diagnostic.alpha_zero_maximum_density_ulp_difference >= 1U &&
      diagnostic.alpha_zero_maximum_density_ulp_difference <= 2U &&
      (size != 1 ||
       diagnostic.alpha_zero_maximum_density_ulp_difference == 1U) &&
      diagnostic.alpha_zero_maximum_temperature_ulp_difference == 0U &&
      diagnostic.alpha_zero_oracle_numeric_lineage != 0U;
  if (rank == 0 && observed) {
    const double first_alpha =
        first_admissible_valid
            ? diagnostic.samples[diagnostic.first_admissible_sample].alpha
            : std::numeric_limits<double>::quiet_NaN();
    std::cout << "candidate-globalization-red attempt/corrector="
              << diagnostic.attempt << '/'
              << static_cast<unsigned>(diagnostic.corrector)
              << " dp/dh="
              << diagnostic.maximum_absolute_pressure_correction << '/'
              << diagnostic.maximum_absolute_enthalpy_correction
              << " alpha0-dU/drho/dT="
              << diagnostic.alpha_zero_maximum_velocity_difference << '/'
              << diagnostic.alpha_zero_maximum_density_difference << '/'
              << diagnostic.alpha_zero_maximum_temperature_difference
              << " bits="
              << diagnostic.alpha_zero_velocity_bitwise_equal << '/'
              << diagnostic.alpha_zero_density_bitwise_equal << '/'
              << diagnostic.alpha_zero_temperature_bitwise_equal
              << " oracle-rel-rho/T="
              << diagnostic.alpha_zero_maximum_density_relative_error << '/'
              << diagnostic.alpha_zero_maximum_temperature_relative_error
              << " oracle-ulp-rho/T="
              << diagnostic.alpha_zero_maximum_density_ulp_difference << '/'
              << diagnostic.alpha_zero_maximum_temperature_ulp_difference
              << " oracle-p="
              << diagnostic.alpha_zero_maximum_pressure_difference << '/'
              << diagnostic.alpha_zero_maximum_pressure_relative_error
              << " oracle-material/grad/mu_eff="
              << diagnostic.alpha_zero_material_oracle_error << '/'
              << diagnostic.alpha_zero_gradient_oracle_error << '/'
              << diagnostic.alpha_zero_effective_viscosity_oracle_error
              << " oracle-lineage="
              << diagnostic.alpha_zero_oracle_numeric_lineage
              << " RC/RE=" << diagnostic.baseline_normalized_continuity
              << '/' << diagnostic.baseline_normalized_energy << " -> "
              << diagnostic.selected_normalized_continuity << '/'
              << diagnostic.selected_normalized_energy
              << " C1="
              << diagnostic.corrector_one_baseline_normalized_continuity
              << '/'
              << diagnostic.corrector_one_baseline_normalized_energy << " -> "
              << diagnostic.corrector_one_selected_normalized_continuity
              << '/'
              << diagnostic.corrector_one_selected_normalized_energy
              << '@' << diagnostic.corrector_one_selected_alpha
              << " linear-C1/C2="
              << diagnostic.corrector_one_linear_predicted_normalized_continuity
              << '/'
              << diagnostic.corrector_one_linear_predicted_normalized_energy
              << ','
              << diagnostic.corrector_two_linear_predicted_normalized_continuity
              << '/'
              << diagnostic.corrector_two_linear_predicted_normalized_energy
              << " alpha=" << diagnostic.selected_alpha
              << " alpha1-failure=" << full.first_failing_global_cell << '/'
              << static_cast<unsigned>(full.first_failure_reason)
              << " min-p/T/rho=" << full.minimum_absolute_pressure << '/'
              << full.minimum_temperature << '/' << full.minimum_density
              << " first-admissible=" << first_alpha << '\n';
  }
  if (!passed) {
    std::cerr << "rank " << rank << " candidate-globalization="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << " failed-stage=" << step.failed_stage << " failure="
              << static_cast<unsigned>(step.failure.code) << '/'
              << step.failure.detail
              << " accepted=" << step.accepted
              << " attempts=" << step.attempts
              << " observed=" << observed << '/' << diagnostic.valid
              << " rollback=" << captured_before << '/' << captured_after
              << '/' << committed_rollback_exact
              << " loop/select/replay/commit="
              << diagnostic.production_candidate_loop << '/'
              << diagnostic.selection_valid << '/'
              << diagnostic.replay_valid << '/' << diagnostic.committed
              << " attempt/corrector=" << diagnostic.attempt << '/'
              << static_cast<unsigned>(diagnostic.corrector)
              << " samples/first="
              << static_cast<unsigned>(diagnostic.sample_count) << '/'
              << static_cast<unsigned>(diagnostic.first_admissible_sample)
              << " dp/dh="
              << diagnostic.maximum_absolute_pressure_correction << '/'
              << diagnostic.maximum_absolute_enthalpy_correction
              << " alpha0-dU/drho/dT="
              << diagnostic.alpha_zero_maximum_velocity_difference << '/'
              << diagnostic.alpha_zero_maximum_density_difference << '/'
              << diagnostic.alpha_zero_maximum_temperature_difference
              << " bits="
              << diagnostic.alpha_zero_velocity_bitwise_equal << '/'
              << diagnostic.alpha_zero_density_bitwise_equal << '/'
              << diagnostic.alpha_zero_temperature_bitwise_equal
              << " oracle-rel-rho/T="
              << diagnostic.alpha_zero_maximum_density_relative_error << '/'
              << diagnostic.alpha_zero_maximum_temperature_relative_error
              << " oracle-ulp-rho/T="
              << diagnostic.alpha_zero_maximum_density_ulp_difference << '/'
              << diagnostic.alpha_zero_maximum_temperature_ulp_difference
              << " oracle-p="
              << diagnostic.alpha_zero_maximum_pressure_difference << '/'
              << diagnostic.alpha_zero_maximum_pressure_relative_error
              << " oracle-material/grad/mu_eff="
              << diagnostic.alpha_zero_material_oracle_error << '/'
              << diagnostic.alpha_zero_gradient_oracle_error << '/'
              << diagnostic.alpha_zero_effective_viscosity_oracle_error
              << " oracle-lineage="
              << diagnostic.alpha_zero_oracle_numeric_lineage
              << " RC/RE=" << diagnostic.baseline_normalized_continuity
              << '/' << diagnostic.baseline_normalized_energy << " -> "
              << diagnostic.selected_normalized_continuity << '/'
              << diagnostic.selected_normalized_energy
              << " C1="
              << diagnostic.corrector_one_baseline_normalized_continuity
              << '/'
              << diagnostic.corrector_one_baseline_normalized_energy << " -> "
              << diagnostic.corrector_one_selected_normalized_continuity
              << '/'
              << diagnostic.corrector_one_selected_normalized_energy
              << '@' << diagnostic.corrector_one_selected_alpha
              << " linear-C1/C2="
              << diagnostic.corrector_one_linear_predicted_normalized_continuity
              << '/'
              << diagnostic.corrector_one_linear_predicted_normalized_energy
              << ','
              << diagnostic.corrector_two_linear_predicted_normalized_continuity
              << '/'
              << diagnostic.corrector_two_linear_predicted_normalized_energy
              << " alpha=" << diagnostic.selected_alpha
              << " full=" << full.admissible << '/'
              << full.first_failing_global_cell << '/'
              << static_cast<unsigned>(full.first_failure_reason) << '/'
              << full.minimum_absolute_pressure << '/'
              << full.minimum_temperature << '/'
              << full.minimum_density << " terminal="
              << step.piso.eos_residual << '/'
              << step.piso.continuity_residual << '/'
              << step.piso.energy_residual << '/'
              << step.piso.closed_mass_residual << '/'
              << step.piso.gauge_residual << '\n';
    if (rank == 0) {
      for (std::size_t index = 0U; index < reported_sample_count; ++index) {
        const auto& sample = diagnostic.samples[index];
        std::cerr << "  ladder[" << index << "] alpha=" << sample.alpha
                  << " admissible=" << sample.admissible
                  << " finite=" << sample.state_and_flux_finite
                  << " RC/RE=" << sample.normalized_continuity << '/'
                  << sample.normalized_energy << '\n';
      }
    }
  }
  return passed;
}

bool run_local_donor_product(int rank) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  ValidatedModel model = test::product_model({8 * size, 8, 8});
  model.fingerprint = UINT64_C(0x18000d001);
  // Exercise the local-donor paired flux under a fixed absolute-pressure
  // authority.  The predictor density and pressure/EOS density intentionally
  // have distinct revisions while retaining one storage lineage.
  model.pressure_reference = PressureReferenceKind::boundary_absolute;
  for (BoundaryFaceSpec& face : model.boundaries) {
    face.flow_kind = BoundaryKind::symmetry;
    face.thermal_kind = BoundaryKind::none;
    face.pressure = 101325.0;
    face.temperature = 300.0;
  }
  model.boundaries[1U].flow_kind = BoundaryKind::pressure_outlet;
  model.time.initial_dt = 0.006;
  model.time.minimum_dt = 0.006;
  model.time.maximum_dt = 0.006;
  model.time.maximum_growth = 1.0;
  model.time.maximum_retries = 1U;
  CompiledCasePlan plan;
  Status status = ProductCompiler::compile(MPI_COMM_WORLD, model, {}, plan);
  ProductDriver driver;
  if (status)
    status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
  RestartExpected expected;
  if (status) status = driver.restart_expected(expected);

  RestartImage image;
  if (status) {
    image.global_cells = expected.global_cells;
    image.patch = expected.target_patch;
    image.plan = expected.plan;
    image.schema = expected.schema;
    image.geometry = expected.geometry;
    image.time = 0.006;
    image.dt = 0.006;
    image.pressure_reference = 101325.0;
    image.step = 1U;
    image.controller_state = 1U;
    image.backward_euler_recovery = true;
    const std::size_t cells =
        static_cast<std::size_t>(image.patch.cells.x) *
        static_cast<std::size_t>(image.patch.cells.y) *
        static_cast<std::size_t>(image.patch.cells.z);
    image.fields.reserve(expected.fields.size);
    for (std::size_t index = 0U; index < expected.fields.size; ++index) {
      const RestartExpectedField descriptor = expected.fields.data[index];
      RestartImageField field;
      field.role = descriptor.role;
      field.field = descriptor.field;
      field.components = descriptor.components;
      field.values.assign(cells * descriptor.components, 0.0);
      if (descriptor.role == RestartFieldRole::enthalpy) {
        std::fill(field.values.begin(), field.values.end(), 300000.0);
      }
      image.fields.push_back(std::move(field));
    }
    const Int3 local = image.patch.cells;
    image.final_mass_flux[0U].assign(
        static_cast<std::size_t>(local.x + 1) * local.y * local.z, 0.0);
    image.final_mass_flux[1U].assign(
        static_cast<std::size_t>(local.x) * (local.y + 1) * local.z, 0.0);
    image.final_mass_flux[2U].assign(
        static_cast<std::size_t>(local.x) * local.y * (local.z + 1), 0.0);
    constexpr double divergence = 200.0;
    const double cell_volume = 1.0 / (512.0 * static_cast<double>(size));
    image.final_mass_flux[0U][1U] = divergence * cell_volume;
  }
  if (status) status = driver.initialize_restart(image);
  CommittedProductBits committed_before;
  const bool captured_before =
      status && capture_committed_product_bits(driver, committed_before);
  DriverStepReport step;
  if (status)
    status = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, step);
  CommittedProductBits committed_after;
  const bool captured_after =
      capture_committed_product_bits(driver, committed_after);
  // A fatal proposal is intentionally retired and replaced by a distinct BE
  // recovery generation.  That ticket is controller state, not committed
  // flow state.  Compare the latter bit-for-bit while independently proving
  // the one-step generation transition instead of hiding it behind the old
  // V1 Restart image's constant-zero controller placeholder.
  const bool fatal_controller_recovery =
      captured_before && captured_after &&
      committed_before.controller_state !=
          std::numeric_limits<std::uint64_t>::max() &&
      committed_after.controller_state ==
          committed_before.controller_state + 1U;
  if (captured_before && captured_after)
    committed_before.controller_state = committed_after.controller_state;
  const bool rollback_exact =
      captured_before && captured_after && fatal_controller_recovery &&
      same_committed_product_bits(committed_before, committed_after);
  std::uint64_t factor_bits = 0U;
  static_assert(sizeof(factor_bits) == sizeof(
                    step.thermophysical_predictor.mass_flux_scale));
  std::memcpy(&factor_bits,
              &step.thermophysical_predictor.mass_flux_scale,
              sizeof(factor_bits));
  const bool factor_rank_identical =
      same_u64(factor_bits, MPI_COMM_WORLD);
  const bool terminal_physics_rejected =
      step.piso.eos_residual > model.solver.terminal.eos ||
      step.piso.continuity_residual > model.solver.terminal.continuity ||
      step.piso.energy_residual > model.solver.terminal.continuity ||
      step.piso.closed_mass_residual >
          model.solver.terminal.closed_mass ||
      step.piso.gauge_residual > model.solver.terminal.gauge;
  const auto converged = [](const LinearSolveResult& solve) noexcept {
    return solve.status &&
           solve.termination == LinearTermination::converged;
  };
  // This family keeps eight x-cells per rank, so increasing the rank count
  // also refines the physical problem.  The finer cases can prove the forced
  // history flux incompatible at the provisional CFL authority, during C1
  // candidate globalization, or at the final coupled audit.  All three are
  // explicit fail-closed authorities, and none may publish state or flux.
  const MomentumAdvectiveCflCertificate& provisional_cfl =
      step.momentum_predictor_limiter.advective_cfl;
  const bool provisional_cfl_rejected =
      status.detail == 10211U && step.failed_stage == 30U &&
      step.piso.pressure_solve_calls == 0U && provisional_cfl.valid() &&
      provisional_cfl.failure_witness.valid &&
      provisional_cfl.out_max > provisional_cfl.limit;
  const bool candidate_rejected =
      status.detail == 5792U && step.failed_stage == 44U &&
      step.piso.pressure_solve_calls == 1U &&
      converged(step.piso.pressure[0U]);
  const bool terminal_rejected =
      status.detail == 1505U && step.failed_stage == 60U &&
      step.piso.pressure_solve_calls == 2U &&
      converged(step.piso.pressure[0U]) &&
      converged(step.piso.pressure[1U]) && terminal_physics_rejected;
  const bool passed =
      !status && status.code == StatusCode::rejected_step &&
      !step.accepted && step.attempts == 1U && rollback_exact &&
      (provisional_cfl_rejected || candidate_rejected || terminal_rejected) &&
      step.thermophysical_predictor.low_state ==
          ThermophysicalLowStateKind::bdf_local_donor_flux &&
      step.thermophysical_predictor.enthalpy_solve_calls == 0U &&
      step.thermophysical_predictor.enthalpy_endpoint_alpha == 1.0 &&
      step.thermophysical_predictor.bdf_endpoint_alpha == 1.0 &&
      step.thermophysical_predictor.source_endpoint_alpha == 1.0 &&
      std::isfinite(step.thermophysical_predictor.mass_flux_scale) &&
      step.thermophysical_predictor.mass_flux_scale >= 0.0 &&
      step.thermophysical_predictor.mass_flux_scale < 1.0 &&
      step.thermophysical_predictor.theta >= 0.0 &&
      step.thermophysical_predictor.theta < 1.0 &&
      step.thermophysical_predictor.low_order_substeps == 1U &&
      step.thermophysical_predictor.low_order_halo_exchanges == 1U &&
      step.thermophysical_predictor.low_order_transport_passes == 1U &&
      step.thermophysical_predictor.blocking_collectives > 1U &&
      factor_rank_identical;
  if (rank == 0) {
    std::cout << std::setprecision(17)
              << "local-donor-candidate-red size=" << size << " status="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << " low/scale/theta="
              << static_cast<unsigned>(
                     step.thermophysical_predictor.low_state)
              << '/' << step.thermophysical_predictor.mass_flux_scale << '/'
              << step.thermophysical_predictor.theta << " terminal="
              << step.piso.eos_residual << '/'
              << step.piso.continuity_residual << '/'
              << step.piso.energy_residual << '/'
              << step.piso.closed_mass_residual << '/'
              << step.piso.gauge_residual
              << " rollback=" << rollback_exact << '\n';
  }
  if (!passed)
    std::cerr << "rank " << rank << " local-donor-product="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << " accepted=" << step.accepted
              << " attempts=" << step.attempts
              << " stage=" << step.failed_stage << " low="
              << static_cast<unsigned>(
                     step.thermophysical_predictor.low_state)
              << " bdf_alpha="
              << step.thermophysical_predictor.bdf_endpoint_alpha
              << " mass_scale="
              << step.thermophysical_predictor.mass_flux_scale
              << " theta=" << step.thermophysical_predictor.theta
              << " collectives="
              << step.thermophysical_predictor.blocking_collectives
              << " pressure_calls="
              << static_cast<unsigned>(step.piso.pressure_solve_calls)
              << " C2="
              << static_cast<unsigned>(step.piso.pressure[1U].status.code)
              << '/' << step.piso.pressure[1U].status.detail << '/'
              << static_cast<unsigned>(step.piso.pressure[1U].termination)
              << " audit="
              << step.piso.pressure[1U].final_convergence_metric << '/'
              << step.piso.pressure[1U].convergence_limit << '/'
              << step.piso.pressure[1U].convergence_application_scale << '/'
              << step.piso.pressure[1U].convergence_unscaled_metric << '/'
              << step.piso.pressure[1U].convergence_maximum_depletion << '/'
              << step.piso.pressure[1U].convergence_operator_parity_error
              << " terminal=" << step.piso.eos_residual << '/'
              << step.piso.continuity_residual << '/'
              << step.piso.energy_residual << '/'
              << step.piso.closed_mass_residual << '/'
              << step.piso.gauge_residual
              << " rollback=" << rollback_exact
              << '\n';
  return passed;
}

bool run_mass_flow_product(int rank) {
  constexpr double kBaselineDt = 5.0e-4;
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  ValidatedModel model = test::product_model({17, 11, 7});
  // Retry recovery is covered by the dedicated driver test.  This product
  // baseline must remain a one-transaction BE -> BDF2 advance after the
  // BDF-history-aware pressure-energy normalization removed the former false
  // terminal-energy rejection.
  model.time.initial_dt = kBaselineDt;
  model.pressure_reference = PressureReferenceKind::boundary_absolute;
  for (BoundaryFaceSpec& face : model.boundaries) {
    face.flow_kind = BoundaryKind::symmetry;
    face.thermal_kind = BoundaryKind::none;
    face.pressure = 98000.0;
    face.temperature = 315.0;
  }
  model.boundaries[0U].flow_kind = BoundaryKind::mass_flow_inlet;
  model.boundaries[0U].direction = {1.0, 0.0, 0.0};
  model.boundaries[0U].velocity = {1.0, 0.0, 0.0};
  model.boundaries[0U].mass_flow_rate = 0.25;
  model.boundaries[1U].flow_kind = BoundaryKind::pressure_outlet;
  CompiledCasePlan plan;
  Status status = ProductCompiler::compile(MPI_COMM_WORLD, model, {}, plan);
  PlanSummary summary;
  if (status) summary = plan.summary();
  ProductDriver driver;
  if (status)
    status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
  DriverInitialState initial;
  initial.pressure_reference = 98000.0;
  initial.temperature = 315.0;
  if (status) status = driver.initialize(initial);
  DriverStepReport first;
  if (status)
    status = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, first);
  DriverStepReport second;
  if (status)
    status = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, second);
  RestartSnapshot restart;
  if (status) status = driver.committed_restart_snapshot(restart);
  double local_inlet = 0.0;
  if (status && restart.patch.begin.x == 0) {
    for (std::int32_t z = 0; z < restart.patch.cells.z; ++z)
      for (std::int32_t y = 0; y < restart.patch.cells.y; ++y)
        local_inlet += restart.final_mass_flux.x.unchecked({0, y, z});
  }
  double global_inlet = 0.0;
  if (MPI_Allreduce(&local_inlet, &global_inlet, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD) != MPI_SUCCESS)
    status = {StatusCode::mpi_failure, 1U};
  if (!status)
    std::cerr << "rank " << rank << " mass-flow="
              << static_cast<unsigned>(status.code) << "/" << status.detail
              << " first=" << first.accepted << '/' << first.attempts
              << '/' << first.failed_stage << '/'
              << static_cast<unsigned>(first.failure.code) << ':'
              << first.failure.detail << " residuals="
              << first.piso.eos_residual << '/'
              << first.piso.continuity_residual << '/'
              << first.piso.energy_residual << '/'
              << first.piso.closed_mass_residual << '/'
              << first.piso.gauge_residual << " second="
              << second.accepted << '/' << second.attempts << '/'
              << second.failed_stage << '/'
              << static_cast<unsigned>(second.failure.code) << ':'
              << second.failure.detail << " residuals="
              << second.piso.eos_residual << '/'
              << second.piso.continuity_residual << '/'
              << second.piso.energy_residual << '/'
              << second.piso.closed_mass_residual << '/'
              << second.piso.gauge_residual << '\n';
  if (status && std::abs(global_inlet - 0.25) > 1.0e-10)
    std::cerr << std::setprecision(17) << "rank " << rank
              << " inlet=" << global_inlet
              << " error=" << global_inlet - 0.25 << '\n';

  const LinearSolveResult& mass_c1 = first.piso.pressure[0U];
  const LinearSolveResult& mass_c2 = first.piso.pressure[1U];
  const auto accepted_solve = [](const LinearSolveResult& solve) {
    return solve.status &&
           (solve.termination == LinearTermination::converged ||
            solve.termination == LinearTermination::zero_rhs) &&
           std::isfinite(solve.initial_true_residual) &&
           std::isfinite(solve.final_true_residual) &&
           solve.final_true_residual <= solve.initial_true_residual;
  };
  const auto capture_role = [&](const LinearSolveResult& solve) {
    const std::uint64_t cycles = solve.recycle_cycle_corrections;
    return accepted_solve(solve) && solve.convergence_audits == 0U &&
           solve.convergence_rejections == 0U &&
           solve.recycle_offered_directions == 0U &&
           !solve.recycle_projection_attempted &&
           solve.recycle_retained_directions == 0U &&
           solve.recycle_operator_applies == 0U &&
           solve.recycle_reduction_calls == 0U &&
           !solve.recycle_projection_accepted && cycles > 0U &&
           solve.recycle_capture_cycle_attempts == cycles &&
           solve.recycle_capture_vector_passes == 2U * cycles &&
           solve.recycle_capture_reduction_calls == cycles &&
           solve.recycle_capture_blocking_operations == 2U * cycles;
  };
  const auto projection_role = [&](const LinearSolveResult& solve,
                                   std::uint64_t offered) {
    return accepted_solve(solve) && solve.convergence_audits == 0U &&
           solve.convergence_rejections == 0U && offered > 0U &&
           solve.recycle_offered_directions == offered &&
           solve.recycle_projection_attempted &&
           solve.recycle_retained_directions == offered &&
           solve.recycle_operator_applies == offered + 1U &&
           solve.recycle_reduction_calls > 0U &&
           solve.recycle_projection_accepted &&
           solve.recycle_cycle_corrections == 0U &&
           solve.recycle_capture_vector_passes == 0U &&
           solve.recycle_capture_cycle_attempts == 0U &&
           solve.recycle_capture_reduction_calls == 0U &&
           solve.recycle_capture_blocking_operations == 0U &&
           std::isfinite(solve.recycle_projected_true_residual) &&
           solve.recycle_projected_true_residual <
               solve.initial_true_residual;
  };
  const bool mass_role_matrix =
      first.piso.pressure_solve_calls == 2U && capture_role(mass_c1) &&
      projection_role(mass_c2, mass_c1.recycle_cycle_corrections);

  const LinearSolveResult& bdf_c1 = second.piso.pressure[0U];
  const LinearSolveResult& bdf_c2 = second.piso.pressure[1U];
  const bool bdf_role_matrix =
      second.piso.pressure_solve_calls == 2U && capture_role(bdf_c1) &&
      projection_role(bdf_c2, bdf_c1.recycle_cycle_corrections);
  // Bind semantic resource relationships instead of a brittle count
  // snapshot.  The BDF2 transaction includes momentum work beyond its two
  // pressure solves, while open flow remains free of IBM traffic.
  const bool bdf_resource_relationship =
      second.resources.structured_exchanges > 0U &&
      second.resources.reduction_collectives > 0U &&
      second.resources.reduction_logical_bytes > 0U &&
      second.resources.linear_iterations >
          bdf_c1.iterations + bdf_c2.iterations &&
      second.resources.preconditioner_applications >=
          bdf_c1.preconditioner_applies + bdf_c2.preconditioner_applies &&
      second.resources.exact_numeric_refills >= 2U &&
      second.resources.hierarchy_rebuilds == 0U &&
      second.resources.ibm_exchanges == 0U &&
      second.resources.ibm_messages == 0U &&
      second.resources.ibm_bytes == 0U &&
      (size == 1
           ? second.resources.structured_messages == 0U &&
                 second.resources.structured_bytes == 0U &&
                 second.resources.mg_blocking_collectives == 0U &&
                 second.resources.mg_collective_logical_bytes == 0U
           : second.resources.structured_messages > 0U &&
                 second.resources.structured_bytes > 0U &&
                 second.resources.mg_blocking_collectives > 0U &&
                 second.resources.mg_collective_logical_bytes > 0U);
  const auto terminal_physics = [&](const DriverStepReport& step) {
    return std::isfinite(step.piso.eos_residual) &&
           step.piso.eos_residual <= summary.terminal_eos_tolerance &&
           std::isfinite(step.piso.continuity_residual) &&
           step.piso.continuity_residual <=
               summary.terminal_continuity_tolerance &&
           std::isfinite(step.piso.energy_residual) &&
           step.piso.energy_residual <=
               summary.terminal_continuity_tolerance &&
           std::isfinite(step.piso.closed_mass_residual) &&
           step.piso.closed_mass_residual <=
               summary.terminal_closed_mass_tolerance &&
           std::isfinite(step.piso.gauge_residual) &&
           step.piso.gauge_residual <= summary.terminal_gauge_tolerance;
  };
  const bool physical_terminal_gates =
      terminal_physics(first) && terminal_physics(second);
  const auto committed_cfl = [&](const DriverStepReport& step) {
    constexpr double slack =
        64.0 * std::numeric_limits<double>::epsilon();
    const CommittedConvectiveCflCertificate& certificate =
        step.piso.committed_convective_cfl;
    return step.piso.final_flux_revision != 0U && certificate.valid() &&
           certificate.final_flux == step.piso.final_flux_revision &&
           certificate.dt == step.proposal.dt &&
           certificate.out_max ==
               step.piso.committed_convective_cfl_out_max &&
           certificate.absolute_max ==
               step.piso.committed_convective_cfl_abs_max &&
           certificate.limit == step.piso.committed_convective_cfl_limit &&
           certificate.out_winner.out == certificate.out_max &&
           certificate.absolute_winner.absolute ==
               certificate.absolute_max &&
           !certificate.failure_witness.valid &&
           std::isfinite(step.piso.committed_convective_cfl_out_max) &&
           step.piso.committed_convective_cfl_out_max > 0.0 &&
           std::isfinite(step.piso.committed_convective_cfl_abs_max) &&
           step.piso.committed_convective_cfl_abs_max > 0.0 &&
           step.piso.committed_convective_cfl_limit ==
               model.time.convective_cfl &&
           step.piso.committed_convective_cfl_out_max <=
               model.time.convective_cfl * (1.0 + slack) &&
           same_u64(wire_bits(step.piso.committed_convective_cfl_out_max),
                    MPI_COMM_WORLD) &&
           same_u64(wire_bits(step.piso.committed_convective_cfl_abs_max),
                    MPI_COMM_WORLD) &&
           same_u64(wire_bits(step.piso.committed_convective_cfl_limit),
                    MPI_COMM_WORLD) &&
           same_u64(certificate.density, MPI_COMM_WORLD) &&
           same_u64(certificate.final_flux, MPI_COMM_WORLD) &&
           same_u64(certificate.activity_collective, MPI_COMM_WORLD) &&
           same_u64(static_cast<std::uint64_t>(
                        certificate.out_winner.global_cell.x),
                    MPI_COMM_WORLD) &&
           same_u64(static_cast<std::uint64_t>(
                        certificate.out_winner.global_cell.y),
                    MPI_COMM_WORLD) &&
           same_u64(static_cast<std::uint64_t>(
                        certificate.out_winner.global_cell.z),
                    MPI_COMM_WORLD) &&
           same_u64(static_cast<std::uint64_t>(certificate.out_winner.rank),
                    MPI_COMM_WORLD) &&
           same_u64(static_cast<std::uint64_t>(
                        certificate.absolute_winner.rank),
                    MPI_COMM_WORLD) &&
           same_u64(wire_bits(certificate.out_winner.density_volume),
                    MPI_COMM_WORLD) &&
           same_u64(wire_bits(
                        certificate.absolute_winner.density_volume),
                    MPI_COMM_WORLD);
  };
  const bool committed_cfl_role =
      committed_cfl(first) && committed_cfl(second);
  const auto advective_cfl = [&](const DriverStepReport& step) {
    constexpr double slack =
        64.0 * std::numeric_limits<double>::epsilon();
    const MomentumAdvectiveCflCertificate& cfl =
        step.momentum_predictor_limiter.advective_cfl;
    return cfl.valid() && cfl.dt == step.proposal.dt &&
           cfl.face_flux != step.piso.final_flux_revision &&
           cfl.out_max >= 0.0 && cfl.absolute_max >= 0.0 &&
           cfl.limit == model.time.convective_cfl &&
           cfl.out_max <= cfl.limit * (1.0 + slack) &&
           !cfl.failure_witness.valid &&
           same_u64(cfl.plan, MPI_COMM_WORLD) &&
           same_u64(cfl.time, MPI_COMM_WORLD) &&
           same_u64(cfl.density, MPI_COMM_WORLD) &&
           same_u64(cfl.face_flux, MPI_COMM_WORLD) &&
           same_u64(wire_bits(cfl.out_max), MPI_COMM_WORLD) &&
           same_u64(wire_bits(cfl.absolute_max), MPI_COMM_WORLD);
  };
  const bool advective_cfl_role =
      advective_cfl(first) && advective_cfl(second);
  const auto limiter_metrics = [&](const DriverStepReport& step) {
    const MomentumPredictorLimiterReport& limiter =
        step.momentum_predictor_limiter;
    const bool applicable = limiter.correction_metrics_applicable;
    return (applicable
                ? limiter.active_correction_faces > 0U &&
                      limiter.activations <=
                          limiter.active_correction_faces &&
                      limiter.minimum_face_alpha >= 0.0 &&
                      limiter.minimum_face_alpha <= 1.0 &&
                      limiter.limited_face_fraction ==
                          static_cast<double>(limiter.activations) /
                              limiter.active_correction_faces
                : limiter.active_correction_faces == 0U &&
                      limiter.activations == 0U &&
                      limiter.minimum_face_alpha == 0.0 &&
                      limiter.limited_face_fraction == 0.0 &&
                      !limiter.limited) &&
           same_u64(limiter.active_correction_faces, MPI_COMM_WORLD) &&
           same_u64(limiter.activations, MPI_COMM_WORLD) &&
           same_u64(wire_bits(limiter.minimum_face_alpha), MPI_COMM_WORLD) &&
           same_u64(wire_bits(limiter.limited_face_fraction),
                    MPI_COMM_WORLD);
  };
  const bool limiter_role = limiter_metrics(first) && limiter_metrics(second) &&
      !first.momentum_predictor_limiter.limited &&
      first.momentum_predictor_limiter.theta == 1.0 &&
      second.momentum_predictor_limiter.limited &&
      second.momentum_predictor_limiter.theta > 0.0 &&
      second.momentum_predictor_limiter.theta < 1.0;
  const bool momentum_role =
      valid_momentum_solve(first.momentum_predictor_solve,
                           {1U, 0U, 0U}, {3U, 0U, 0U},
                           {1U, 0U, 0U}, {7U, 1U, 1U}) &&
      valid_momentum_solve(second.momentum_predictor_solve,
                           {2U, 0U, 0U}, {4U, 1U, 1U},
                           {2U, 0U, 0U}, {8U, 4U, 4U});
  const bool final_flux_role =
      restart.final_mass_flux.certificate.valid() &&
      restart.final_mass_flux.revision == second.piso.final_flux_revision;
  const bool inlet_role = std::abs(global_inlet - 0.25) <= 1.0e-10;
  const bool accepted_bdf_role =
      status && first.accepted && first.attempts == 1U &&
      first.failure.code == StatusCode::ok && first.failed_stage == 0U &&
      first.proposal.origin == StepOrigin::fresh_start &&
      first.proposal.dt == kBaselineDt &&
      first.proposal.bdf.order == 1U && first.effective_bdf.order == 1U &&
      second.accepted && second.attempts == 1U &&
      second.failure.code == StatusCode::ok && second.failed_stage == 0U &&
      second.proposal.origin == StepOrigin::accepted &&
      second.proposal.attempt == 0U && second.proposal.dt > 0.0 &&
      second.proposal.bdf.order == 2U && second.effective_bdf.order == 2U &&
      std::abs(second.accepted_time -
               (first.accepted_time + second.proposal.dt)) <=
          16.0 * std::numeric_limits<double>::epsilon() *
              std::max(1.0, std::abs(second.accepted_time));
  if (!accepted_bdf_role || !mass_role_matrix || !bdf_role_matrix ||
                 !bdf_resource_relationship || !physical_terminal_gates ||
                 !committed_cfl_role || !advective_cfl_role ||
                 !limiter_role || !momentum_role || !final_flux_role ||
                 !inlet_role) {
    std::cerr << std::setprecision(17) << "rank " << rank
              << " BE C1 iter/op/pre/red/capture="
              << mass_c1.iterations << '/' << mass_c1.operator_applies << '/'
              << mass_c1.preconditioner_applies << '/'
              << mass_c1.reduction_calls << '/'
              << mass_c1.recycle_cycle_corrections
              << " C2 iter/op/pre/red/projection=" << mass_c2.iterations
              << '/' << mass_c2.operator_applies << '/'
              << mass_c2.preconditioner_applies << '/'
              << mass_c2.reduction_calls << '/'
              << mass_c2.recycle_offered_directions << '/'
              << mass_c2.recycle_retained_directions << '/'
              << mass_c2.recycle_operator_applies << '/'
              << mass_c2.recycle_reduction_calls << '\n'
              << "rank " << rank
              << " BDF2 C1 init/iter/op/pre/red/capture="
              << bdf_c1.initial_true_residual << '/' << bdf_c1.iterations
              << '/' << bdf_c1.operator_applies << '/'
              << bdf_c1.preconditioner_applies << '/'
              << bdf_c1.reduction_calls << '/'
              << bdf_c1.recycle_cycle_corrections
              << " C2 iter/op/pre/red/projection=" << bdf_c2.iterations
              << '/' << bdf_c2.operator_applies << '/'
              << bdf_c2.preconditioner_applies << '/'
              << bdf_c2.reduction_calls << '/'
              << bdf_c2.recycle_offered_directions << '/'
              << bdf_c2.recycle_retained_directions << '/'
              << bdf_c2.recycle_operator_applies << '/'
              << bdf_c2.recycle_reduction_calls
              << " resources=" << second.resources.structured_exchanges
              << '/' << second.resources.structured_messages << '/'
              << second.resources.structured_bytes << '/'
              << second.resources.reduction_collectives << '/'
              << second.resources.reduction_logical_bytes << '/'
              << second.resources.mg_blocking_collectives << '/'
              << second.resources.mg_collective_logical_bytes << '/'
              << second.resources.linear_iterations << '/'
              << second.resources.preconditioner_applications
              << " checks=" << mass_role_matrix << '/'
              << bdf_role_matrix << '/' << bdf_resource_relationship
              << " audits=" << mass_c2.convergence_audits << '/'
              << mass_c2.convergence_rejections << '/'
              << bdf_c2.convergence_audits << '/'
              << bdf_c2.convergence_rejections
              << " terminal=" << first.piso.eos_residual << '/'
              << first.piso.continuity_residual << '/'
              << first.piso.energy_residual << '/'
              << first.piso.closed_mass_residual << '/'
              << first.piso.gauge_residual << ';'
              << second.piso.eos_residual << '/'
              << second.piso.continuity_residual << '/'
              << second.piso.energy_residual << '/'
              << second.piso.closed_mass_residual << '/'
              << second.piso.gauge_residual
              << " cfl=" << first.piso.committed_convective_cfl_out_max
              << '/' << first.piso.committed_convective_cfl_abs_max << ';'
              << second.piso.committed_convective_cfl_out_max << '/'
              << second.piso.committed_convective_cfl_abs_max << '/'
              << second.piso.committed_convective_cfl_limit
              << " cert=" << first.piso.committed_convective_cfl.valid()
              << ':' << first.piso.committed_convective_cfl.plan << ':'
              << first.piso.committed_convective_cfl.correction_state
              << ':' << first.piso.committed_convective_cfl.density << ':'
              << first.piso.committed_convective_cfl.final_flux << ':'
              << first.piso.committed_convective_cfl.out_winner.rank << ':'
              << first.piso.committed_convective_cfl.absolute_winner.rank
              << ';' << second.piso.committed_convective_cfl.valid() << ':'
              << second.piso.committed_convective_cfl.plan << ':'
              << second.piso.committed_convective_cfl.correction_state << ':'
              << second.piso.committed_convective_cfl.density << ':'
              << second.piso.committed_convective_cfl.final_flux << ':'
              << second.piso.committed_convective_cfl.out_winner.rank << ':'
              << second.piso.committed_convective_cfl.absolute_winner.rank
              << " advective="
              << first.momentum_predictor_limiter.advective_cfl.out_max
              << '/'
              << first.momentum_predictor_limiter.advective_cfl.absolute_max
              << ';'
              << second.momentum_predictor_limiter.advective_cfl.out_max
              << '/'
              << second.momentum_predictor_limiter.advective_cfl.absolute_max
              << " final=" << accepted_bdf_role << '/' << limiter_role << '/'
              << momentum_role << '/' << final_flux_role << '/'
              << committed_cfl_role << '/'
              << advective_cfl_role << '/'
              << inlet_role << " base="
              << static_cast<unsigned>(status.code) << ':' << status.detail
              << '/' << first.accepted << ':' << first.attempts << ':'
              << first.failed_stage << ':'
              << static_cast<unsigned>(first.failure.code) << ':'
              << first.failure.detail << ':'
              << static_cast<unsigned>(first.proposal.origin) << ':'
              << first.proposal.attempt << ':' << first.proposal.dt << '/'
              << second.accepted << ':'
              << second.attempts << ':'
              << second.failed_stage << ':'
              << static_cast<unsigned>(second.failure.code) << ':'
              << second.failure.detail << ':'
              << static_cast<unsigned>(second.proposal.origin) << ':'
              << second.proposal.attempt << ':' << second.proposal.dt << ':'
              << static_cast<unsigned>(second.proposal.bdf.order) << ':'
              << static_cast<unsigned>(second.effective_bdf.order) << '\n';
    std::cerr << "rank " << rank << " momentum theta="
              << first.momentum_predictor_limiter.theta << '/'
              << second.momentum_predictor_limiter.theta;
    for (const DriverStepReport* step : {&first, &second}) {
      std::cerr << " momentum_solve["
                << static_cast<unsigned>(step->momentum_predictor_solve.solve_calls)
                << ']';
      for (const LinearSolveResult& solve :
           step->momentum_predictor_solve.components) {
        std::cerr << ' ' << solve.iterations << '/'
                  << solve.operator_applies << '/'
                  << solve.preconditioner_applies << '/'
                  << solve.reduction_calls << '/'
                  << solve.initial_true_residual << '/'
                  << solve.final_true_residual;
      }
    }
    std::cerr << '\n';
  }
  return accepted_bdf_role && mass_role_matrix &&
         bdf_role_matrix && bdf_resource_relationship &&
         physical_terminal_gates && committed_cfl_role &&
         advective_cfl_role &&
         limiter_role && momentum_role && final_flux_role && inlet_role;
}

bool run_multispecies_open_product(int rank) {
  const std::filesystem::path data_root =
      std::filesystem::path{HUNDUN_V04_SOURCE_ROOT} / "tests" / "data";
  ValidatedModel model = multispecies_open_model();
  ThermodynamicsPlan oracle_thermodynamics;
  Status status = ThermodynamicsPlan::compile(
      model.thermophysics,
      {model.transported_scalars.data(), model.transported_scalars.size()},
      oracle_thermodynamics);
  CompiledCasePlan plan;
  if (status)
    status = ProductCompiler::compile(MPI_COMM_WORLD, model, data_root, plan);
  PlanSummary summary;
  if (status) summary = plan.summary();
  ProductDriver driver;
  if (status)
    status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
  const std::array<double, 3U> initial_scalars{{0.10, 0.20, 0.05}};
  DriverInitialState initial;
  initial.pressure_reference = 98000.0;
  initial.temperature = 315.0;
  initial.velocity = {1.0e-4, 0.0, 0.0};
  initial.transported_scalars =
      {initial_scalars.data(), initial_scalars.size()};
  if (status) status = driver.initialize(initial);
  if (status)
    detail::arm_pressure_energy_candidate_globalization_once_for_test();
  DriverStepReport step;
  if (status)
    status = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, step);
  detail::PressureEnergyCandidateGlobalizationDiagnostic diagnostic;
  const bool observed =
      detail::pressure_energy_candidate_globalization_diagnostic_for_test(
          diagnostic);
  detail::clear_pressure_energy_candidate_globalization_for_test();
  RestartSnapshot restart;
  if (status) status = driver.committed_restart_snapshot(restart);

  double local_inlet = 0.0;
  double local_outlet = 0.0;
  double local_expected_inlet = 0.0;
  std::uint64_t local_inlet_faces = 0U;
  std::uint64_t local_outlet_faces = 0U;
  std::uint64_t local_outlet_backflow_faces = 0U;
  const RestartFieldView* pressure_field = nullptr;
  if (status) {
    for (std::size_t index = 0U; index < restart.fields.size; ++index)
      if (restart.fields.data[index].role ==
          RestartFieldRole::pressure_perturbation) {
        pressure_field = restart.fields.data + index;
        break;
      }
    if (pressure_field == nullptr)
      status = {StatusCode::invalid_plan, 1U};
  }
  const std::array<double, 2U> inlet_composition{{0.1001, 0.2001}};
  double inlet_enthalpy = 0.0;
  double inlet_cp = 0.0;
  double inlet_gas_constant = 0.0;
  if (status)
    status = oracle_thermodynamics.mixture_enthalpy(
        315.0, {inlet_composition.data(), inlet_composition.size()},
        inlet_enthalpy, inlet_cp, inlet_gas_constant);
  const double inlet_area =
      (model.mesh.upper.y - model.mesh.lower.y) /
      static_cast<double>(restart.global_cells.y) *
      (model.mesh.upper.z - model.mesh.lower.z) /
      static_cast<double>(restart.global_cells.z);
  if (status && restart.patch.begin.x == 0) {
    for (std::int32_t z = 0; z < restart.patch.cells.z; ++z)
      for (std::int32_t y = 0; y < restart.patch.cells.y; ++y) {
        local_inlet += restart.final_mass_flux.x.unchecked({0, y, z});
        ThermoState inlet_state;
        status = oracle_thermodynamics.evaluate(
            restart.pressure_reference +
                pressure_field->values.unchecked({0, y, z}, 0U),
            inlet_enthalpy,
            {inlet_composition.data(), inlet_composition.size()},
            {1.0e-4, 0.0, 0.0}, inlet_state, 315.0);
        if (!status) break;
        local_expected_inlet += inlet_state.rho * 1.0e-4 * inlet_area;
        ++local_inlet_faces;
      }
  }
  if (status &&
      restart.patch.begin.x + restart.patch.cells.x ==
          restart.global_cells.x) {
    for (std::int32_t z = 0; z < restart.patch.cells.z; ++z)
      for (std::int32_t y = 0; y < restart.patch.cells.y; ++y) {
        const double value = restart.final_mass_flux.x.unchecked(
            {restart.patch.cells.x, y, z});
        local_outlet += value;
        if (value < 0.0) ++local_outlet_backflow_faces;
        ++local_outlet_faces;
      }
  }
  std::array<double, 3U> local_flux{{local_inlet, local_outlet,
                                     local_expected_inlet}};
  std::array<double, 3U> global_flux{};
  std::array<std::uint64_t, 3U> local_faces{{
      local_inlet_faces, local_outlet_faces,
      local_outlet_backflow_faces}};
  std::array<std::uint64_t, 3U> global_faces{};
  if (MPI_Allreduce(local_flux.data(), global_flux.data(), 3, MPI_DOUBLE,
                    MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Allreduce(local_faces.data(), global_faces.data(), 3, MPI_UINT64_T,
                    MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS)
    status = {StatusCode::mpi_failure, 1U};

  std::uint64_t local_species_fields = 0U;
  std::uint64_t local_passive_fields = 0U;
  bool local_scalars_finite_and_bounded = true;
  bool local_scalar_changed = false;
  if (status) {
    for (std::size_t index = 0U; index < restart.fields.size; ++index) {
      const RestartFieldView& field = restart.fields.data[index];
      std::size_t scalar_index = initial_scalars.size();
      if (field.role == RestartFieldRole::independent_species) {
        scalar_index = static_cast<std::size_t>(local_species_fields++);
      } else if (field.role == RestartFieldRole::transported_scalar) {
        scalar_index = 2U + static_cast<std::size_t>(local_passive_fields++);
      } else {
        continue;
      }
      if (scalar_index >= initial_scalars.size()) {
        local_scalars_finite_and_bounded = false;
        continue;
      }
      for (std::int32_t z = 0; z < field.values.interior.z; ++z)
        for (std::int32_t y = 0; y < field.values.interior.y; ++y)
          for (std::int32_t x = 0; x < field.values.interior.x; ++x) {
            const double value = field.values.unchecked({x, y, z}, 0U);
            local_scalars_finite_and_bounded &=
                std::isfinite(value) && value >= 0.0 && value <= 1.0;
            local_scalar_changed |=
                wire_bits(value) != wire_bits(initial_scalars[scalar_index]);
          }
    }
  }
  const bool scalar_roles =
      local_species_fields == 2U && local_passive_fields == 1U;
  const bool terminal =
      status && std::isfinite(step.piso.eos_residual) &&
      step.piso.eos_residual <= summary.terminal_eos_tolerance &&
      std::isfinite(step.piso.continuity_residual) &&
      step.piso.continuity_residual <=
          summary.terminal_continuity_tolerance &&
      std::isfinite(step.piso.energy_residual) &&
      step.piso.energy_residual <=
          summary.terminal_continuity_tolerance &&
      std::isfinite(step.piso.gauge_residual) &&
      step.piso.gauge_residual <= summary.terminal_gauge_tolerance;
  const bool final_boundary =
      observed && diagnostic.valid && diagnostic.production_candidate_loop &&
      diagnostic.replay_valid && diagnostic.committed &&
      diagnostic.final_boundary_flux_certified &&
      diagnostic.final_boundary_canonical_lineage != 0U &&
      diagnostic.final_physical_flux_provenance != 0U;
  const bool boundary_flux =
      global_faces[0U] > 0U && global_faces[1U] > 0U &&
      global_faces[2U] == 0U && global_flux[0U] > 0.0 &&
      global_flux[1U] > 0.0 &&
      std::abs(global_flux[0U] - global_flux[2U]) <=
          128.0 * std::numeric_limits<double>::epsilon() *
              std::max({1.0, std::abs(global_flux[0U]),
                        std::abs(global_flux[2U])});
  const bool passed =
      status && step.accepted && terminal && final_boundary && boundary_flux &&
      scalar_roles && local_scalars_finite_and_bounded &&
      local_scalar_changed && restart.final_mass_flux.certificate.valid();
  if (!passed && rank == 0) {
    std::cerr << std::setprecision(17)
              << "multispecies-open status="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << " accepted/attempts/stage=" << step.accepted << '/'
              << step.attempts << '/' << step.failed_stage
              << " failure=" << static_cast<unsigned>(step.failure.code)
              << '/' << step.failure.detail << " gates="
              << step.piso.eos_residual << '/'
              << step.piso.continuity_residual << '/'
              << step.piso.energy_residual << '/'
              << step.piso.gauge_residual << " final=" << observed << '/'
              << diagnostic.production_candidate_loop << '/'
              << diagnostic.replay_valid << '/' << diagnostic.committed << '/'
              << diagnostic.final_boundary_flux_certified << " phi="
              << global_flux[0U] << '/' << global_flux[1U] << "/expected="
              << global_flux[2U] << " faces=" << global_faces[0U] << '/'
              << global_faces[1U] << "/backflow=" << global_faces[2U]
              << " scalars="
              << local_species_fields << '/' << local_passive_fields << '/'
              << local_scalars_finite_and_bounded << '/'
              << local_scalar_changed << " predictor="
              << step.thermophysical_predictor.failure.valid << '/'
              << static_cast<unsigned>(
                     step.thermophysical_predictor.failure.reason)
              << '/'
              << static_cast<unsigned>(
                     step.thermophysical_predictor.failure.field)
              << '/'
              << static_cast<unsigned>(
                     step.thermophysical_predictor.failure.constraint)
              << " cell="
              << step.thermophysical_predictor.failure.global_index.x << ','
              << step.thermophysical_predictor.failure.global_index.y << ','
              << step.thermophysical_predictor.failure.global_index.z
              << " scalar-mask/value/div="
              << step.thermophysical_predictor.failure.scalar_mask << '/'
              << step.thermophysical_predictor.failure.observed_value << '/'
              << step.thermophysical_predictor.failure.divergence << '\n';
  }
  return passed;
}

enum class FreshOpenFluxCase : std::uint8_t {
  velocity_inlet,
  pressure_outlet_backflow,
  mass_flow_inlet,
};

bool run_fresh_open_boundary_eos_flux_case(int rank,
                                           FreshOpenFluxCase selected) {
  const std::filesystem::path data_root =
      std::filesystem::path{HUNDUN_V04_SOURCE_ROOT} / "tests" / "data";
  ValidatedModel model = multispecies_open_model();
  model.fingerprint = selected == FreshOpenFluxCase::velocity_inlet
                          ? UINT64_C(0x18000cb01)
                      : selected ==
                                FreshOpenFluxCase::pressure_outlet_backflow
                          ? UINT64_C(0x18000cb02)
                          : UINT64_C(0x18000cb03);

  BoundaryFaceSpec& inlet = model.boundaries[0U];
  inlet.temperature = 421.0;
  inlet.velocity = {2.5e-2, 0.0, 0.0};
  inlet.scalars[0U].value = 0.4103;
  inlet.scalars[1U].value = 0.1307;
  if (selected == FreshOpenFluxCase::mass_flow_inlet) {
    inlet.flow_kind = BoundaryKind::mass_flow_inlet;
    inlet.direction = {1.0, 0.0, 0.0};
    inlet.mass_flow_rate = 0.037;
  }

  BoundaryFaceSpec& outlet = model.boundaries[1U];
  outlet.pressure = 87321.0;
  outlet.allow_backflow = true;
  outlet.backflow_temperature = 367.0;
  outlet.backflow_velocity = {-1.75e-2, 0.0, 0.0};
  outlet.scalars[0U].backflow_kind = ScalarBoundaryKind::dirichlet;
  outlet.scalars[0U].backflow_value = 0.1209;
  outlet.scalars[1U].backflow_kind = ScalarBoundaryKind::dirichlet;
  outlet.scalars[1U].backflow_value = 0.3111;

  ThermodynamicsPlan oracle;
  const char* phase = "oracle-compile";
  Status status = ThermodynamicsPlan::compile(
      model.thermophysics,
      {model.transported_scalars.data(), model.transported_scalars.size()},
      oracle);
  CompiledCasePlan plan;
  if (status) {
    phase = "product-compile";
    status = ProductCompiler::compile(MPI_COMM_WORLD, model, data_root, plan);
  }
  ProductDriver driver;
  if (status) {
    phase = "create";
    status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
  }

  const std::array<double, 3U> initial_scalars{{0.1901, 0.2302, 0.031}};
  DriverInitialState initial;
  initial.pressure_reference = 101325.0;
  initial.temperature = 289.0;
  initial.velocity = selected == FreshOpenFluxCase::pressure_outlet_backflow
                         ? Real3{-1.75e-2, 0.0, 0.0}
                         : Real3{2.5e-2, 0.0, 0.0};
  initial.transported_scalars =
      {initial_scalars.data(), initial_scalars.size()};
  if (status) {
    phase = "initialize";
    status = driver.initialize(initial);
  }

  CommittedOutputSnapshot snapshot;
  ConstFaceFluxView committed_flux;
  if (status) {
    phase = "snapshot";
    status = driver.committed_output_snapshot(snapshot);
  }
  if (status) {
    phase = "final-flux";
    status = driver.committed_final_mass_flux_for_test(committed_flux);
  }
  const SnapshotFieldView* pressure = nullptr;
  if (status) {
    for (std::size_t field = 0U; field < snapshot.fields.size; ++field)
      if (snapshot.fields.data[field].stable_name == "pi") {
        pressure = snapshot.fields.data + field;
        break;
      }
    if (pressure == nullptr)
      status = {StatusCode::invalid_plan, 1U};
  }

  const bool selected_outlet =
      selected == FreshOpenFluxCase::pressure_outlet_backflow;
  const BoundaryFaceSpec& expected_boundary = selected_outlet ? outlet : inlet;
  const std::array<double, 2U> expected_composition{{
      selected_outlet ? outlet.scalars[0U].backflow_value
                      : inlet.scalars[0U].value,
      selected_outlet ? outlet.scalars[1U].backflow_value
                      : inlet.scalars[1U].value,
  }};
  const double expected_temperature = selected_outlet
                                          ? outlet.backflow_temperature
                                          : inlet.temperature;
  const Real3 expected_velocity =
      selected_outlet ? outlet.backflow_velocity : inlet.velocity;
  double expected_enthalpy = 0.0;
  double expected_cp = 0.0;
  double expected_gas = 0.0;
  if (status)
    status = oracle.mixture_enthalpy(
        expected_temperature,
        {expected_composition.data(), expected_composition.size()},
        expected_enthalpy, expected_cp, expected_gas);

  const Int3 global_cells =
      snapshot.geometry == nullptr ? Int3{} : snapshot.geometry->global_cells();
  const double face_area =
      global_cells.y > 0 && global_cells.z > 0
          ? (model.mesh.upper.y - model.mesh.lower.y) /
                static_cast<double>(global_cells.y) *
                (model.mesh.upper.z - model.mesh.lower.z) /
                static_cast<double>(global_cells.z)
          : 0.0;
  double local_maximum_relative_error = 0.0;
  double local_mass_flow = 0.0;
  std::uint64_t local_face_count = 0U;
  std::uint64_t local_wrong_sign = 0U;
  const bool owns_selected_face = selected_outlet
                                      ? snapshot.patch.begin.x +
                                                snapshot.patch.cells.x ==
                                            global_cells.x
                                      : snapshot.patch.begin.x == 0;
  if (status && owns_selected_face) {
    const std::int32_t face_x =
        selected_outlet ? snapshot.patch.cells.x : 0;
    const std::int32_t owner_x =
        selected_outlet ? snapshot.patch.cells.x - 1 : 0;
    for (std::int32_t z = 0; z < snapshot.patch.cells.z && status; ++z)
      for (std::int32_t y = 0; y < snapshot.patch.cells.y; ++y) {
        const double actual = committed_flux.x.unchecked({face_x, y, z});
        local_mass_flow += actual;
        if (selected == FreshOpenFluxCase::mass_flow_inlet) {
          if (!(actual > 0.0)) ++local_wrong_sign;
          ++local_face_count;
          continue;
        }
        const double absolute_pressure =
            selected_outlet
                ? expected_boundary.pressure
                : initial.pressure_reference +
                      pressure->values.unchecked({owner_x, y, z}, 0U);
        ThermoState expected;
        status = oracle.evaluate(
            absolute_pressure, expected_enthalpy,
            {expected_composition.data(), expected_composition.size()},
            expected_velocity, expected, expected_temperature);
        if (!status) break;
        const double oracle_flux =
            expected.rho * expected_velocity.x * face_area;
        const double scale =
            std::max({1.0, std::abs(oracle_flux), std::abs(actual)});
        local_maximum_relative_error = std::max(
            local_maximum_relative_error, std::abs(actual - oracle_flux) / scale);
        if (selected_outlet && !(actual < 0.0))
          ++local_wrong_sign;
        ++local_face_count;
      }
  }
  double global_maximum_relative_error = 0.0;
  double global_mass_flow = 0.0;
  std::array<std::uint64_t, 2U> local_counts{{local_face_count,
                                               local_wrong_sign}};
  std::array<std::uint64_t, 2U> global_counts{};
  if (MPI_Allreduce(&local_maximum_relative_error,
                    &global_maximum_relative_error, 1, MPI_DOUBLE, MPI_MAX,
                    MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Allreduce(&local_mass_flow, &global_mass_flow, 1, MPI_DOUBLE,
                    MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Allreduce(local_counts.data(), global_counts.data(), 2,
                    MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS)
    status = {StatusCode::mpi_failure, 1U};

  const bool flux_matches =
      selected == FreshOpenFluxCase::mass_flow_inlet
          ? std::abs(global_mass_flow - inlet.mass_flow_rate) <=
                128.0 * std::numeric_limits<double>::epsilon() *
                    std::max(1.0, std::abs(inlet.mass_flow_rate))
          : global_maximum_relative_error <=
                128.0 * std::numeric_limits<double>::epsilon();
  const bool passed =
      status && snapshot.step == 0U && committed_flux.certificate.valid() &&
      global_counts[0U] > 0U && global_counts[1U] == 0U &&
      flux_matches;
  if (!passed && rank == 0)
    std::cerr << std::setprecision(17) << "fresh-open-eos phase=" << phase
              << " case="
              << static_cast<unsigned>(selected) << " status="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << " step/cert=" << snapshot.step << '/'
              << committed_flux.certificate.valid()
              << " faces/wrong-sign=" << global_counts[0U] << '/'
              << global_counts[1U] << " relative-error="
              << global_maximum_relative_error << " mass-flow="
              << global_mass_flow << '/' << inlet.mass_flow_rate << '\n';
  return passed;
}

bool run_fresh_open_boundary_eos_flux(int rank) {
  const bool inlet = run_fresh_open_boundary_eos_flux_case(
      rank, FreshOpenFluxCase::velocity_inlet);
  const bool outlet = run_fresh_open_boundary_eos_flux_case(
      rank, FreshOpenFluxCase::pressure_outlet_backflow);
  const bool mass_flow = run_fresh_open_boundary_eos_flux_case(
      rank, FreshOpenFluxCase::mass_flow_inlet);
  return inlet && outlet && mass_flow;
}

bool run_unsupported_candidate_boundary_product(int rank) {
  enum class UnsupportedCase : std::uint8_t {
    static_inlet,
    total_inlet,
    nscbc_outlet,
  };
  const std::array<UnsupportedCase, 3U> cases{
      UnsupportedCase::static_inlet, UnsupportedCase::total_inlet,
      UnsupportedCase::nscbc_outlet};
  bool passed = true;
  for (const UnsupportedCase selected : cases) {
    ValidatedModel model = test::product_model({8, 7, 6});
    model.time.initial_dt = 1.0e-5;
    model.pressure_reference = PressureReferenceKind::boundary_absolute;
    for (BoundaryFaceSpec& face : model.boundaries) {
      face.flow_kind = BoundaryKind::symmetry;
      face.thermal_kind = BoundaryKind::none;
      face.pressure = 98000.0;
      face.temperature = 315.0;
    }
    BoundaryFaceSpec& inlet = model.boundaries[0U];
    inlet.direction = {1.0, 0.0, 0.0};
    inlet.velocity = {0.1, 0.0, 0.0};
    inlet.pressure = 98100.0;
    inlet.temperature = 315.0;
    BoundaryFaceSpec& outlet = model.boundaries[1U];
    outlet.flow_kind = BoundaryKind::pressure_outlet;
    if (selected == UnsupportedCase::static_inlet) {
      inlet.flow_kind = BoundaryKind::static_state_inlet;
    } else if (selected == UnsupportedCase::total_inlet) {
      inlet.flow_kind = BoundaryKind::total_state_inlet;
      inlet.total_pressure = 98200.0;
      inlet.total_temperature = 315.0;
    } else {
      inlet.flow_kind = BoundaryKind::velocity_inlet;
      outlet.flow_kind = BoundaryKind::nscbc_outlet;
      outlet.pressure = 98000.0;
      outlet.temperature = 315.0;
      outlet.relaxation = 0.1;
      outlet.mach_limit = 0.9;
      outlet.allow_backflow = false;
    }

    // The static/total models are rejected by the boundary compiler before
    // a distributed runtime exists.  Validate that compile-time rejection
    // independently on every rank, then use the distributed communicator for
    // the NSCBC model that reaches the ProductDriver candidate-scope gate.
    const MPI_Comm case_communicator =
        selected == UnsupportedCase::nscbc_outlet ? MPI_COMM_WORLD
                                                   : MPI_COMM_SELF;
    CompiledCasePlan plan;
    Status setup =
        ProductCompiler::compile(case_communicator, model, {}, plan);
    ProductDriver driver;
    if (setup)
      setup = ProductDriver::create(case_communicator, std::move(plan),
                                    driver);
    DriverInitialState initial;
    initial.pressure_reference = 98000.0;
    initial.temperature = 315.0;
    initial.velocity = {0.1, 0.0, 0.0};
    if (setup) setup = driver.initialize(initial);
    CommittedStateBits before;
    const bool captured_before =
        setup && capture_committed_state_bits(driver, before);
    DriverStepReport step;
    Status advance_status = setup;
    if (advance_status)
      advance_status =
          driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, step);
    CommittedStateBits after;
    const bool captured_after =
        setup && capture_committed_state_bits(driver, after);
    const Status observed = setup ? advance_status : setup;
    const bool rank_consensus =
        same_u64((static_cast<std::uint64_t>(observed.code) << 32U) |
                     observed.detail,
                 MPI_COMM_WORLD) &&
        (!setup ||
         same_u64((static_cast<std::uint64_t>(step.failed_stage) << 32U) |
                      step.failure.detail,
                  MPI_COMM_WORLD));
    const bool compile_rejected =
        !setup && setup.code == StatusCode::invalid_plan && rank_consensus;
    const bool attempt_rejected =
        setup && !advance_status &&
        advance_status.code == StatusCode::invalid_plan && !step.accepted &&
        step.attempts == 1U && step.failure.code == StatusCode::invalid_plan &&
        step.failed_stage < 46U && captured_before && captured_after &&
        same_committed_state_bits(before, after) && rank_consensus;
    const bool rejected_before_legacy_publication =
        compile_rejected || attempt_rejected;
    if (!rejected_before_legacy_publication) {
      std::cerr << "rank " << rank << " unsupported-boundary="
                << static_cast<unsigned>(selected) << " setup="
                << static_cast<unsigned>(setup.code) << '/' << setup.detail
                << " advance="
                << static_cast<unsigned>(advance_status.code) << '/'
                << advance_status.detail << " accepted/attempts/stage="
                << step.accepted << '/' << step.attempts << '/'
                << step.failed_stage << " failure="
                << static_cast<unsigned>(step.failure.code) << '/'
                << step.failure.detail << " rollback=" << captured_before
                << '/' << captured_after << '/'
                << (captured_before && captured_after &&
                    same_committed_state_bits(before, after))
                << '\n';
    }
    passed = passed && rejected_before_legacy_publication;
  }
  return passed;
}

bool run_warm_start_lifecycle_product(int rank) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  ValidatedModel model = test::product_model({8, 7, 6});
  model.pressure_reference = PressureReferenceKind::boundary_absolute;
  model.time.initial_dt = 1.0e-5;
  for (BoundaryFaceSpec& face : model.boundaries) {
    face.flow_kind = BoundaryKind::symmetry;
    face.thermal_kind = BoundaryKind::none;
    face.pressure = 98000.0;
    face.temperature = 315.0;
  }
  model.boundaries[0U].flow_kind = BoundaryKind::mass_flow_inlet;
  model.boundaries[0U].direction = {1.0, 0.0, 0.0};
  model.boundaries[0U].velocity = {0.01, 0.0, 0.0};
  model.boundaries[0U].mass_flow_rate = 0.01;
  model.boundaries[1U].flow_kind = BoundaryKind::pressure_outlet;
  // This fixture deliberately starts with a negative streamwise velocity to
  // exercise the warm/cold Krylov seed paths.  Declare the resulting outlet
  // inflow explicitly so Fresh initialization remains a physically complete
  // boundary state under the final-flux contract.
  model.boundaries[1U].allow_backflow = true;
  model.boundaries[1U].backflow_velocity = {-2.0, 0.0, 0.0};
  model.boundaries[1U].backflow_temperature = 315.0;
  CompiledCasePlan plan;
  Status status = ProductCompiler::compile(MPI_COMM_WORLD, model, {}, plan);
  ProductDriver driver;
  if (status)
    status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
  DriverInitialState initial;
  initial.pressure_reference = 98000.0;
  initial.temperature = 315.0;
  initial.velocity = {-2.0, 0.0, 0.0};
  if (status) status = driver.initialize(initial);
  DriverStepReport first;
  if (status)
    status = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, first);
  DriverStepReport second;
  if (status)
    status = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, second);

  const LinearSolveResult& c1 = second.piso.pressure[0U];
  const LinearSolveResult& c2 = second.piso.pressure[1U];
  detail::PressureCorrectionWarmStartDiagnostic warm_start;
  const bool observed_warm_start =
      detail::pressure_correction_warm_start_diagnostic_for_test(warm_start);
  CompiledCasePlan cold_plan;
  Status cold_status =
      ProductCompiler::compile(MPI_COMM_WORLD, model, {}, cold_plan);
  ProductDriver cold_driver;
  if (cold_status)
    cold_status = ProductDriver::create(MPI_COMM_WORLD, std::move(cold_plan),
                                        cold_driver);
  if (cold_status) cold_status = cold_driver.initialize(initial);
  DriverStepReport cold_first;
  if (cold_status)
    cold_status = cold_driver.advance({1.0, 1.0, 1.0, 1.0, 1.0},
                                      cold_first);
  DriverStepReport cold_second;
  if (cold_status) {
    detail::suppress_pressure_correction_warm_start_once_for_test();
    cold_status = cold_driver.advance({1.0, 1.0, 1.0, 1.0, 1.0},
                                      cold_second);
  }
  detail::PressureCorrectionWarmStartDiagnostic cold_warm_start;
  const bool observed_cold_warm_start =
      detail::pressure_correction_warm_start_diagnostic_for_test(
          cold_warm_start);
  const LinearSolveResult& cold_c1 = cold_second.piso.pressure[0U];
  const LinearSolveResult& cold_c2 = cold_second.piso.pressure[1U];
  const auto converged_solve = [](const LinearSolveResult& solve) {
    return solve.status &&
           solve.termination == LinearTermination::converged &&
           std::isfinite(solve.initial_true_residual) &&
           std::isfinite(solve.final_true_residual) &&
           solve.final_true_residual <= solve.initial_true_residual &&
           solve.convergence_rejections == 0U;
  };
  const auto coherent_capture = [](const LinearSolveResult& solve) {
    return solve.recycle_capture_vector_passes ==
               2U * solve.recycle_capture_cycle_attempts &&
           solve.recycle_capture_reduction_calls ==
               solve.recycle_capture_cycle_attempts &&
           solve.recycle_capture_blocking_operations ==
               2U * solve.recycle_capture_cycle_attempts;
  };
  const auto terminal_physics = [&](const DriverStepReport& step) {
    return step.piso.eos_residual <= model.solver.terminal.eos &&
           step.piso.continuity_residual <=
               model.solver.terminal.continuity &&
           step.piso.energy_residual <= model.solver.terminal.continuity &&
           step.piso.closed_mass_residual <=
               model.solver.terminal.closed_mass &&
           step.piso.gauge_residual <= model.solver.terminal.gauge;
  };
  const auto momentum_solve_semantic =
      [](const MomentumPredictorSolveReport& report) {
        bool valid = report.solve_calls == 3U;
        for (const LinearSolveResult& solve : report.components)
          valid = valid && solve.status &&
                  (solve.termination == LinearTermination::converged ||
                   solve.termination == LinearTermination::zero_rhs) &&
                  std::isfinite(solve.initial_true_residual) &&
                  std::isfinite(solve.final_true_residual) &&
                  solve.final_true_residual <= solve.initial_true_residual;
        return valid;
      };
  const bool warm_seed_consumed =
      observed_warm_start && warm_start.valid &&
      warm_start.origin == StepOrigin::accepted && warm_start.attempt == 0U &&
      warm_start.authority_available && warm_start.used;
  const bool independent_cold_comparison =
      cold_status && cold_first.accepted && cold_second.accepted &&
      cold_first.attempts == 1U && cold_second.attempts == 1U &&
      cold_second.proposal.bdf.order == second.proposal.bdf.order &&
      observed_cold_warm_start && cold_warm_start.valid &&
      cold_warm_start.origin == StepOrigin::accepted &&
      cold_warm_start.attempt == 0U &&
      cold_warm_start.authority_available && !cold_warm_start.used &&
      converged_solve(cold_c1) && converged_solve(cold_c2) &&
      coherent_capture(cold_c1) && coherent_capture(cold_c2) &&
      terminal_physics(cold_first) && terminal_physics(cold_second);
  const bool c1_semantic =
      converged_solve(c1) && coherent_capture(c1) &&
      c1.recycle_cycle_corrections > 0U &&
      c1.recycle_capture_cycle_attempts > 0U;
  const bool c2_semantic =
      converged_solve(c2) && coherent_capture(c2) &&
      c2.recycle_offered_directions > 0U &&
      c2.recycle_retained_directions > 0U &&
      c2.recycle_retained_directions <= c2.recycle_offered_directions &&
      c2.recycle_projection_attempted && c2.recycle_projection_accepted;
  const bool resources_semantic =
      second.resources.structured_exchanges > 0U &&
      second.resources.reduction_collectives > 0U &&
      second.resources.reduction_logical_bytes > 0U &&
      second.resources.linear_iterations >= c1.iterations + c2.iterations &&
      second.resources.preconditioner_applications >=
          c1.preconditioner_applies + c2.preconditioner_applies &&
      second.resources.ibm_exchanges == 0U &&
      second.resources.ibm_messages == 0U &&
      second.resources.ibm_bytes == 0U;
  const bool passed =
      status && first.accepted && first.attempts == 1U && second.accepted &&
      second.attempts == 1U && second.proposal.bdf.order == 2U &&
      second.piso.pressure_solve_calls == 2U && warm_seed_consumed &&
      independent_cold_comparison && c1_semantic && c2_semantic &&
      resources_semantic &&
      terminal_physics(first) && terminal_physics(second) &&
      !first.momentum_predictor_limiter.limited &&
      first.momentum_predictor_limiter.theta == 1.0 &&
      second.momentum_predictor_limiter.limited &&
      second.momentum_predictor_limiter.activations > 0U &&
      second.momentum_predictor_limiter.theta > 0.0 &&
      second.momentum_predictor_limiter.theta < 1.0 &&
      momentum_solve_semantic(first.momentum_predictor_solve) &&
      momentum_solve_semantic(second.momentum_predictor_solve);
  if (!passed) {
    std::cerr << std::setprecision(17) << "rank " << rank
              << " warm-start lifecycle status="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << " accepted/attempts=" << first.accepted << '/'
              << first.attempts << ' ' << second.accepted << '/'
              << second.attempts << " warm=" << observed_warm_start << '/'
              << warm_start.authority_available << '/' << warm_start.used
              << " cold=" << static_cast<unsigned>(cold_status.code) << '/'
              << observed_cold_warm_start << '/'
              << cold_warm_start.authority_available << '/'
              << cold_warm_start.used << '/' << cold_c1.initial_true_residual
              << " C1 init/final/iter/op/pre/red="
              << c1.initial_true_residual << '/' << c1.final_true_residual << '/'
              << c1.iterations << '/' << c1.operator_applies << '/'
              << c1.preconditioner_applies << '/' << c1.reduction_calls
              << " C2 iter/op/pre/red=" << c2.iterations << '/'
              << c2.operator_applies << '/' << c2.preconditioner_applies
              << '/' << c2.reduction_calls << " resources="
              << second.resources.structured_exchanges << '/'
              << second.resources.structured_messages << '/'
              << second.resources.structured_bytes << '/'
              << second.resources.reduction_collectives << '/'
              << second.resources.reduction_logical_bytes << '/'
              << second.resources.mg_blocking_collectives << '/'
              << second.resources.mg_collective_logical_bytes << '/'
              << second.resources.linear_iterations << '/'
              << second.resources.preconditioner_applications
              << " terminal2=" << second.piso.eos_residual << '/'
              << second.piso.continuity_residual << '/'
              << second.piso.energy_residual << '/'
              << second.piso.closed_mass_residual << '/'
              << second.piso.gauge_residual
              << " momentum theta="
              << first.momentum_predictor_limiter.theta << '/'
              << second.momentum_predictor_limiter.theta;
    for (const DriverStepReport* step : {&first, &second}) {
      std::cerr << " momentum_solve["
                << static_cast<unsigned>(step->momentum_predictor_solve.solve_calls)
                << ']';
      for (const LinearSolveResult& solve :
           step->momentum_predictor_solve.components) {
        std::cerr << ' ' << solve.iterations << '/'
                  << solve.operator_applies << '/'
                  << solve.preconditioner_applies << '/'
                  << solve.reduction_calls << '/'
                  << solve.initial_true_residual << '/'
                  << solve.final_true_residual;
      }
    }
    std::cerr << '\n';
  }
  return passed;
}

bool run_immersed_product(int rank) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  ValidatedModel model = test::product_model({16, 16, 16});
  model.mesh.lower = {-2.0, -2.0, -2.0};
  model.mesh.upper = {2.0, 2.0, 2.0};
  model.mesh.minimum_spacing = {0.25, 0.25, 0.25};
  model.time.initial_dt = 2.5e-4;
  model.immersed_boundary = ImmersedBoundarySpec{
      "cylinder_ascii.stl", ImmersedFluidSide::outside};
  const std::filesystem::path data_root =
      std::filesystem::path{HUNDUN_V04_SOURCE_ROOT} / "tests" / "data";
  CompiledCasePlan plan;
  Status status =
      ProductCompiler::compile(MPI_COMM_WORLD, model, data_root, plan);
  PlanSummary summary;
  if (status) summary = plan.summary();
  const bool reconstruction_audit =
      summary.immersed && summary.ibm_boundary_reconstruction.valid &&
      summary.ibm_surface_reconstruction.valid &&
      summary.ibm_boundary_reconstruction.policy ==
          IbmReconstructionPolicy::strict_quadratic &&
      summary.ibm_surface_reconstruction.policy ==
          IbmReconstructionPolicy::strict_quadratic &&
      summary.ibm_boundary_reconstruction.group_count > 0U &&
      summary.ibm_surface_reconstruction.group_count > 0U &&
      summary.ibm_boundary_reconstruction.group_count ==
          summary.ibm_boundary_reconstruction.quadratic_groups &&
      summary.ibm_surface_reconstruction.group_count ==
          summary.ibm_surface_reconstruction.quadratic_groups &&
      summary.ibm_boundary_reconstruction.linear_groups == 0U &&
      summary.ibm_surface_reconstruction.linear_groups == 0U;
  detail::PressureEnergyCandidateStorageDiagnostic storage;
  const bool storage_observed =
      status && detail::pressure_energy_candidate_storage_diagnostic_for_test(
                    storage);
  const bool candidate_donor_lineage =
      storage_observed && storage.valid &&
      storage.candidate_pressure_donor_plan != 0U &&
      storage.candidate_velocity_donor_plan != 0U &&
      storage.candidate_rate_donor_plan != 0U &&
      storage.candidate_pressure_donor_plan !=
          storage.candidate_velocity_donor_plan &&
      storage.candidate_pressure_donor_plan !=
          storage.candidate_rate_donor_plan &&
      storage.candidate_velocity_donor_plan !=
          storage.candidate_rate_donor_plan;
  ProductDriver driver;
  if (status)
    status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
  DriverInitialState initial;
  initial.pressure_reference = 98000.0;
  initial.temperature = 315.0;
  if (status) status = driver.initialize(initial);
  if (status)
    detail::arm_pressure_energy_candidate_globalization_once_for_test();
  DriverStepReport step;
  if (status)
    status = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, step);
  detail::PressureEnergyCandidateGlobalizationDiagnostic diagnostic;
  const bool candidate_observed =
      detail::pressure_energy_candidate_globalization_diagnostic_for_test(
          diagnostic);
  detail::clear_pressure_energy_candidate_globalization_for_test();
  SurfaceForce force;
  FinalForceCertificate certificate;
  if (status) status = driver.committed_surface_force(force, certificate);
  std::uint64_t global_ibm_messages = 0U;
  if (status &&
      MPI_Allreduce(&step.resources.ibm_messages, &global_ibm_messages, 1,
                    MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS)
    status = {StatusCode::mpi_failure, 2U};
  std::uint64_t global_ibm_interface[3U]{};
  const std::uint64_t local_ibm_interface[3U]{
      diagnostic.local_ibm_interface_links,
      diagnostic.local_ibm_interface_nonzero_count,
      diagnostic.local_ibm_interface_negative_zero_count};
  if (status &&
      MPI_Allreduce(local_ibm_interface, global_ibm_interface, 3,
                    MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS)
    status = {StatusCode::mpi_failure, 2U};
  if (!status)
    std::cerr << "rank " << rank << " immersed="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << " stage=" << step.failed_stage
              << " attempts=" << step.attempts << '\n';
  const LinearSolveResult& ibm_c1 = step.piso.pressure[0U];
  const LinearSolveResult& ibm_c2 = step.piso.pressure[1U];
  const bool ibm_role_matrix =
      step.piso.pressure_solve_calls == 2U && ibm_c1.status &&
      ibm_c1.termination == LinearTermination::converged &&
      ibm_c1.iterations == 4U && ibm_c1.operator_applies == 6U &&
      ibm_c1.preconditioner_applies == 4U &&
      ibm_c1.reduction_calls == 14U &&
      ibm_c1.recycle_offered_directions == 0U &&
      !ibm_c1.recycle_projection_attempted &&
      ibm_c1.recycle_retained_directions == 0U &&
      ibm_c1.recycle_operator_applies == 0U &&
      ibm_c1.recycle_reduction_calls == 0U &&
      !ibm_c1.recycle_projection_accepted &&
      ibm_c1.recycle_cycle_corrections == 1U &&
      ibm_c1.recycle_capture_cycle_attempts == 1U &&
      ibm_c1.recycle_capture_vector_passes == 2U &&
      ibm_c1.recycle_capture_reduction_calls == 1U &&
      ibm_c1.recycle_capture_blocking_operations == 2U && ibm_c2.status &&
      ibm_c2.termination == LinearTermination::converged &&
      ibm_c2.iterations == 0U && ibm_c2.operator_applies == 1U &&
      ibm_c2.preconditioner_applies == 0U &&
      ibm_c2.reduction_calls == 4U &&
      ibm_c2.recycle_offered_directions == 1U &&
      ibm_c2.recycle_cycle_corrections == 0U &&
      ibm_c2.recycle_capture_vector_passes == 0U &&
      ibm_c2.recycle_capture_cycle_attempts == 0U &&
      ibm_c2.recycle_capture_reduction_calls == 0U &&
      ibm_c2.recycle_capture_blocking_operations == 0U &&
      (ibm_c2.recycle_projection_attempted
           ? ibm_c2.recycle_retained_directions <=
                     ibm_c2.recycle_offered_directions &&
                 ibm_c2.recycle_operator_applies > 0U &&
                 ibm_c2.recycle_reduction_calls > 0U &&
                 ibm_c2.recycle_projection_accepted &&
                 std::isfinite(ibm_c2.recycle_projected_true_residual) &&
                 ibm_c2.recycle_projected_true_residual <
                     ibm_c2.initial_true_residual
           : ibm_c2.recycle_retained_directions == 0U &&
                 ibm_c2.recycle_operator_applies == 0U &&
                 ibm_c2.recycle_reduction_calls == 0U &&
                 !ibm_c2.recycle_projection_accepted);
  constexpr double kAlphaZeroOracleTolerance = 64.0 *
      std::numeric_limits<double>::epsilon();
  const bool candidate_relations =
      candidate_observed && diagnostic.valid &&
      diagnostic.production_candidate_loop && diagnostic.replay_valid &&
      diagnostic.committed && diagnostic.corrector == 2U &&
      diagnostic.alpha_zero_gradient_oracle_error <=
          kAlphaZeroOracleTolerance &&
      diagnostic.alpha_zero_energy_residual_oracle_error <=
          kAlphaZeroOracleTolerance &&
      global_ibm_interface[0U] > 0U &&
      global_ibm_interface[1U] == 0U &&
      global_ibm_interface[2U] == 0U &&
      diagnostic.candidate_runtime_halo_messages <=
          diagnostic.candidate_sealed_halo_messages &&
      diagnostic.candidate_runtime_halo_bytes <=
          diagnostic.candidate_sealed_halo_bytes &&
      diagnostic.candidate_sealed_halo_messages ==
          storage.corrector_two_resource_contract.merged_halo_messages &&
      diagnostic.candidate_sealed_halo_bytes ==
          storage.corrector_two_resource_contract.merged_halo_bytes;
  const bool resource_relations =
      step.resources.structured_exchanges > 0U &&
      step.resources.ibm_exchanges >= 12U &&
      step.resources.reduction_collectives > 0U &&
      step.resources.predictor_blocking_collectives == 1U &&
      step.resources.linear_iterations ==
          ibm_c1.iterations + ibm_c2.iterations &&
      step.resources.exact_numeric_refills == 1U &&
      step.resources.hierarchy_rebuilds == 1U &&
      step.resources.preconditioner_applications ==
          ibm_c1.preconditioner_applies + ibm_c2.preconditioner_applies &&
      (size == 1 || step.resources.ibm_messages > 0U);
  const bool passed =
      status && step.accepted && step.attempts == 1U &&
      step.failure.code == StatusCode::ok && ibm_role_matrix &&
      certificate.valid() && reconstruction_audit && candidate_donor_lineage &&
      candidate_relations && resource_relations &&
      (size == 1 ? global_ibm_messages == 0U
                 : global_ibm_messages > 0U) &&
      std::isfinite(force.total.x) && std::isfinite(force.total.y) &&
      std::isfinite(force.total.z);
  if (!passed && rank == 0) {
    std::cerr << std::setprecision(17) << "rank " << rank
              << " immersed status/accepted/attempts="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << '/' << step.accepted << '/' << step.attempts
              << " C1 term/init/final/iter/op/pre/red/capture="
              << static_cast<unsigned>(ibm_c1.termination) << '/'
              << ibm_c1.initial_true_residual << '/'
              << ibm_c1.final_true_residual << '/' << ibm_c1.iterations << '/'
              << ibm_c1.operator_applies << '/'
              << ibm_c1.preconditioner_applies << '/'
              << ibm_c1.reduction_calls << '/'
              << ibm_c1.recycle_cycle_corrections
              << " C2 term/init/final/iter/op/pre/red/projection="
              << static_cast<unsigned>(ibm_c2.termination) << '/'
              << ibm_c2.initial_true_residual << '/'
              << ibm_c2.final_true_residual << '/' << ibm_c2.iterations << '/'
              << ibm_c2.operator_applies << '/'
              << ibm_c2.preconditioner_applies << '/'
              << ibm_c2.reduction_calls << '/'
              << ibm_c2.recycle_offered_directions << '/'
              << ibm_c2.recycle_retained_directions << '/'
              << ibm_c2.recycle_operator_applies << '/'
              << ibm_c2.recycle_reduction_calls << '/'
              << ibm_c2.recycle_projection_attempted << '/'
              << ibm_c2.recycle_projection_accepted
              << " resources=" << step.resources.structured_exchanges << '/'
              << step.resources.structured_messages << '/'
              << step.resources.structured_bytes << '/'
              << step.resources.ibm_exchanges << '/'
              << step.resources.ibm_messages << '/'
              << step.resources.ibm_bytes << '/'
              << step.resources.reduction_collectives << '/'
              << step.resources.predictor_blocking_collectives << '/'
              << step.resources.reduction_logical_bytes << '/'
              << step.resources.reduction_tree_messages << '/'
              << step.resources.mg_blocking_collectives << '/'
              << step.resources.mg_collective_logical_bytes << '/'
              << step.resources.linear_iterations << '/'
              << step.resources.exact_numeric_refills << '/'
              << step.resources.hierarchy_rebuilds << '/'
              << step.resources.preconditioner_applications
              << " certificate/messages=" << certificate.valid() << '/'
              << global_ibm_messages << " force=" << force.total.x << ','
              << force.total.y << ',' << force.total.z
              << " candidate=" << candidate_observed << '/'
              << diagnostic.production_candidate_loop << '/'
              << diagnostic.replay_valid << '/' << diagnostic.committed
              << " oracle=" << diagnostic.alpha_zero_gradient_oracle_error
              << '/'
              << diagnostic.alpha_zero_energy_residual_oracle_error
              << " ibm_phi=" << global_ibm_interface[0U] << '/'
              << global_ibm_interface[1U] << '/'
              << global_ibm_interface[2U]
              << " seal=" << diagnostic.candidate_runtime_halo_messages
              << '/' << diagnostic.candidate_sealed_halo_messages << ':'
              << diagnostic.candidate_runtime_halo_bytes << '/'
              << diagnostic.candidate_sealed_halo_bytes << '\n';
  }
  return passed;
}

bool run_rank_change_product(int rank) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 4) return true;
  std::string root;
  if (rank == 0) {
    root = (std::filesystem::temp_directory_path() /
            ("hundun-v04-driver-rank-change-" +
             std::to_string(::getpid())))
               .string();
    std::error_code error;
    std::filesystem::remove_all(root, error);
    if (!error)
      ApplicationService::initialize_case_directory(
          std::filesystem::path{root} / "case");
  }
  std::uint64_t length = root.size();
  MPI_Bcast(&length, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
  if (rank != 0) root.resize(length);
  MPI_Bcast(root.data(), static_cast<int>(length), MPI_CHAR, 0,
            MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  bool passed = true;
  std::filesystem::path restart;
  const std::array<int, 4U> ranks{{1, 2, 4, 1}};
  for (std::size_t phase = 0U; phase < ranks.size(); ++phase) {
    MPI_Comm active = MPI_COMM_NULL;
    MPI_Comm_split(MPI_COMM_WORLD, rank < ranks[phase] ? 0 : MPI_UNDEFINED,
                   rank, &active);
    bool local = true;
    const std::filesystem::path output =
        std::filesystem::path{root} / ("run-" + std::to_string(phase));
    if (active != MPI_COMM_NULL) {
      ApplicationRunOptions options;
      options.case_root = std::filesystem::path{root} / "case";
      options.run_directory = output;
      options.source_root = HUNDUN_V04_SOURCE_ROOT;
      options.restart_directory = restart;
      options.steps = 2U;
      options.output_interval = 1U;
      options.restart_interval = 2U;
      ApplicationRunReport report;
      const Status status = ApplicationService::run(active, options, report);
      local = status && report.accepted_steps == 2U * (phase + 1U) &&
              report.product != 0U;
      if (!local)
        std::cerr << "rank " << rank << " restart phase=" << phase
                  << " status=" << static_cast<unsigned>(status.code) << '/'
                  << status.detail << " accepted=" << report.accepted_steps
                  << " stage=" << report.failed_stage
                  << " attempts=" << report.attempts
                  << '\n';
      MPI_Comm_free(&active);
    }
    const int local_integer = local ? 1 : 0;
    int global_integer = 0;
    MPI_Allreduce(&local_integer, &global_integer, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    passed &= global_integer != 0;
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0 && phase == 1U) {
      std::ifstream evidence(output / "evidence.jsonl", std::ios::binary);
      const std::string text{std::istreambuf_iterator<char>(evidence),
                             std::istreambuf_iterator<char>()};
      passed &= text.find("\"step\":3") != std::string::npos &&
                text.find("\"bdf_order\":2") != std::string::npos &&
                text.find("\"step\":4") != std::string::npos &&
                text.find("\"bdf_order\":2") != std::string::npos &&
                text.find("\"restart_recovery\":true") ==
                    std::string::npos &&
                text.find("\"restart_recovery\":false") !=
                    std::string::npos;
    }
    passed = collective(passed, MPI_COMM_WORLD);
    restart = output / "Restart";
  }
  if (rank == 0) {
    std::error_code error;
    std::filesystem::remove_all(root, error);
    passed &= !error;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  return collective(passed, MPI_COMM_WORLD);
}

bool run_simple_coupling_compile(int rank) {
  ValidatedModel model = test::product_model({8, 8, 8});
  model.solver.coupling = CouplingKind::simple;
  model.fingerprint = UINT64_C(0x18000a11);
  CompiledCasePlan plan;
  const Status status =
      ProductCompiler::compile(MPI_COMM_WORLD, model, {}, plan);
  const PlanSummary summary = status ? plan.summary() : PlanSummary{};
  ProductDriver driver;
  Status runtime = status;
  if (runtime)
    runtime = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
  if (runtime) runtime = driver.initialize({});
  DriverStepReport step;
  if (runtime)
    runtime = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, step);
  return expect(runtime && summary.coupling == CouplingKind::simple &&
                    summary.pressure_correctors == 2U && summary.sealed &&
                    step.accepted && step.piso.pressure_solve_calls == 2U &&
                    step.momentum_predictor_solve.predictor_passes == 2U &&
                    step.piso.pressure[1U].termination ==
                        LinearTermination::zero_rhs &&
                    step.piso.pressure[1U].iterations == 0U &&
                    step.piso.pressure_energy_refinement_solve_calls == 0U,
                rank,
                "SIMPLE performs two momentum-pressure outer iterations");
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (argc == 2 &&
      std::strcmp(argv[1], "--simple-compile-only") == 0) {
    const bool passed =
        collective(run_simple_coupling_compile(rank), MPI_COMM_WORLD);
    MPI_Finalize();
    return passed ? 0 : 1;
  }
  if (argc == 2 &&
      std::strcmp(argv[1], "--candidate-storage-only") == 0) {
    const bool passed = run_candidate_storage_lineage_only(rank);
    MPI_Finalize();
    return passed ? 0 : 1;
  }
  if (argc == 2 &&
      std::strcmp(argv[1], "--localized-immersed-compile-only") == 0) {
    const bool passed = run_localized_immersed_compile_only(rank, false);
    MPI_Finalize();
    return passed ? 0 : 1;
  }
  if (argc == 2 &&
      std::strcmp(argv[1],
                  "--localized-immersed-bind-fail-close-only") == 0) {
    const bool passed = run_localized_immersed_compile_only(rank, true);
    MPI_Finalize();
    return passed ? 0 : 1;
  }
  if (argc == 2 &&
      std::strcmp(argv[1], "--fresh-initialize-contract-only") == 0) {
    const bool passed = collective(
        run_fresh_initialize_invalid_input_contract(rank), MPI_COMM_WORLD);
    MPI_Finalize();
    return passed ? 0 : 1;
  }
  if (argc == 2 &&
      std::strcmp(argv[1], "--restart-allocation-contract-only") == 0) {
    const bool passed = collective(
        run_restart_restore_allocation_contract(rank), MPI_COMM_WORLD);
    MPI_Finalize();
    return passed ? 0 : 1;
  }
  if (argc == 2 &&
      std::strcmp(argv[1], "--temporal-fallback-only") == 0) {
    const bool passed = collective(
        run_temporal_method_fallback_product(rank), MPI_COMM_WORLD);
    MPI_Finalize();
    return passed ? 0 : 1;
  }
  if (argc == 2 &&
      std::strcmp(argv[1], "--mass-flow-only") == 0) {
    const bool passed =
        collective(run_mass_flow_product(rank), MPI_COMM_WORLD);
    MPI_Finalize();
    return passed ? 0 : 1;
  }
  if (argc == 2 &&
      std::strcmp(argv[1], "--multispecies-open-only") == 0) {
    const bool passed =
        collective(run_multispecies_open_product(rank), MPI_COMM_WORLD);
    MPI_Finalize();
    return passed ? 0 : 1;
  }
  if (argc == 2 &&
      std::strcmp(argv[1], "--fresh-open-eos-flux-only") == 0) {
    const bool passed = collective(
        run_fresh_open_boundary_eos_flux(rank), MPI_COMM_WORLD);
    MPI_Finalize();
    return passed ? 0 : 1;
  }
  if (argc == 2 &&
      std::strcmp(argv[1], "--unsupported-boundary-only") == 0) {
    const bool passed = collective(
        run_unsupported_candidate_boundary_product(rank), MPI_COMM_WORLD);
    MPI_Finalize();
    return passed ? 0 : 1;
  }
  if (argc == 2 &&
      std::strcmp(argv[1], "--implicit-enthalpy-only") == 0) {
    const bool passed =
        collective(run_implicit_enthalpy_product(rank), MPI_COMM_WORLD);
    MPI_Finalize();
    return passed ? 0 : 1;
  }
  if (argc == 2 &&
      std::strcmp(argv[1], "--candidate-globalization-red-only") == 0) {
    const bool passed = collective(
        run_pressure_energy_candidate_globalization_red(rank),
        MPI_COMM_WORLD);
    MPI_Finalize();
    return passed ? 0 : 1;
  }
  if (argc == 2 &&
      std::strcmp(argv[1], "--local-donor-only") == 0) {
    const bool passed =
        collective(run_local_donor_product(rank), MPI_COMM_WORLD);
    MPI_Finalize();
    return passed ? 0 : 1;
  }
  if (argc == 2 &&
      std::strcmp(argv[1], "--warm-start-only") == 0) {
    const bool passed = collective(
        run_warm_start_lifecycle_product(rank), MPI_COMM_WORLD);
    MPI_Finalize();
    return passed ? 0 : 1;
  }
  if (argc == 2 &&
      std::strcmp(argv[1], "--immersed-only") == 0) {
    const bool passed =
        collective(run_immersed_product(rank), MPI_COMM_WORLD);
    MPI_Finalize();
    return passed ? 0 : 1;
  }
  if (argc == 2 &&
      std::strcmp(argv[1], "--rank-change-only") == 0) {
    const bool passed =
        collective(run_rank_change_product(rank), MPI_COMM_WORLD);
    MPI_Finalize();
    return passed ? 0 : 1;
  }
  CompiledCasePlan plan;
  const Status status = ProductCompiler::compile(
      MPI_COMM_WORLD, test::product_model({17, 11, 7}), {}, plan);
  bool passed = true;
  if (!status) {
    std::cerr << "rank " << rank << " status="
              << static_cast<unsigned>(status.code)
              << " detail=" << status.detail << '\n';
  }
  passed &= expect(static_cast<bool>(status), rank,
                   "nondivisible product freezes collectively");
  if (status) {
    passed &= expect(same_u64(plan.fingerprint(), MPI_COMM_WORLD), rank,
                     "semantic product fingerprint is rank invariant");
    passed &= valid_candidate_storage_lineage(plan, rank);
    const Span<const ProductFreezePhase> order = plan.freeze_order();
    passed &= expect(order.size == test::kFreezeOrder.size(), rank,
                     "every rank records all phases");
    for (std::size_t index = 0U;
         index < order.size && index < test::kFreezeOrder.size(); ++index) {
      passed &= expect(order.data[index] == test::kFreezeOrder[index], rank,
                       "freeze order is rank invariant");
    }
    const PlanSummary summary = plan.summary();
    const std::uint64_t local_cells =
        static_cast<std::uint64_t>(summary.local_cells.x) *
        static_cast<std::uint64_t>(summary.local_cells.y) *
        static_cast<std::uint64_t>(summary.local_cells.z);
    std::uint64_t global_cells = 0U;
    MPI_Allreduce(&local_cells, &global_cells, 1, MPI_UINT64_T, MPI_SUM,
                  MPI_COMM_WORLD);
    passed &= expect(global_cells == UINT64_C(17) * 11U * 7U, rank,
                     "nondivisible local patches cover the global domain once");
    passed &= expect(summary.coupling == CouplingKind::piso &&
                         summary.pressure_correctors == 2U && summary.sealed &&
                         !summary.exact_numeric_certified &&
                         !summary.preconditioner_setup_certified,
                     rank, "all ranks seal the same lifecycle state");
    const FrozenExecutionGraph* graph = plan.execution_graph();
    const FrozenStage* predictor_stage =
        graph == nullptr ? nullptr : graph->stage(10U);
    passed &= expect(predictor_stage != nullptr &&
                         predictor_stage->resources.linear_iterations == 128U,
                     rank,
                     "predictor graph budgets the frozen enthalpy endpoint");
    passed &= expect(plan.state_storage_address() != 0U &&
                         plan.krylov_storage_address() != 0U &&
                         plan.mg_storage_address() != 0U,
                     rank, "rank-local storages are bound");
    ProductDriver driver;
    const char* runtime_phase = "create";
    Status runtime =
        ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
    DriverInitialState initial;
    initial.pressure_reference = 98000.0;
    initial.temperature = 315.0;
    if (runtime) {
      runtime_phase = "initialize";
      runtime = driver.initialize(initial);
    }
    DriverStepReport first;
    if (runtime) {
      runtime_phase = "first";
      runtime = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, first);
    }
    DriverStepReport second;
    if (runtime) {
      runtime_phase = "second";
      runtime = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, second);
    }
    const auto coherent_coupled_role = [](const LinearSolveResult& result) {
      // Exact operation counts are decomposition- and seed-path details.  The
      // product contract is that the coupled solve converges without an audit
      // rejection and that any recycled/captured direction has a complete,
      // internally consistent accounting trail.
      return result.status &&
             result.termination == LinearTermination::converged &&
             result.convergence_rejections == 0U &&
             result.recycle_retained_directions <=
                 result.recycle_offered_directions &&
             (!result.recycle_projection_accepted ||
              result.recycle_projection_attempted) &&
             result.recycle_capture_vector_passes ==
                 2U * result.recycle_capture_cycle_attempts &&
             result.recycle_capture_reduction_calls ==
                 result.recycle_capture_cycle_attempts &&
             result.recycle_capture_blocking_operations ==
                 2U * result.recycle_capture_cycle_attempts &&
             std::isfinite(result.initial_true_residual) &&
             std::isfinite(result.final_true_residual);
    };
    RestartSnapshot committed;
    if (runtime) runtime = driver.committed_restart_snapshot(committed);
    passed &= expect(runtime && first.accepted && second.accepted &&
                         first.attempts == 1U && second.attempts == 1U &&
                         first.proposal.bdf.order == 1U &&
                         first.effective_bdf.order == 1U &&
                         first.thermophysical_predictor_calls == 1U &&
                         !first.temporal_method_fallback &&
                         second.proposal.bdf.order == 2U &&
                         second.effective_bdf.order == 2U &&
                         second.thermophysical_predictor_calls == 1U &&
                         !second.temporal_method_fallback &&
                         first.piso.pressure_solve_calls == 2U &&
                         second.piso.pressure_solve_calls == 2U &&
                         first.thermophysical_predictor.enthalpy_solve_calls ==
                             0U &&
                         second.thermophysical_predictor.enthalpy_solve_calls ==
                             0U &&
                         first.thermophysical_predictor
                                 .enthalpy_endpoint_alpha == 1.0 &&
                         second.thermophysical_predictor
                                 .enthalpy_endpoint_alpha == 1.0 &&
                         first.thermophysical_predictor
                                 .bdf_endpoint_alpha == 1.0 &&
                         second.thermophysical_predictor
                                 .bdf_endpoint_alpha == 1.0 &&
                         first.thermophysical_predictor
                                 .source_endpoint_alpha == 1.0 &&
                         second.thermophysical_predictor
                                 .source_endpoint_alpha == 1.0 &&
                         first.failure.code == StatusCode::ok &&
                         second.failure.code == StatusCode::ok &&
                         coherent_coupled_role(first.piso.pressure[0U]) &&
                         coherent_coupled_role(first.piso.pressure[1U]) &&
                         coherent_coupled_role(second.piso.pressure[0U]) &&
                         coherent_coupled_role(second.piso.pressure[1U]) &&
                         second.accepted_step == 2U &&
                         committed.final_mass_flux.certificate.valid() &&
                         committed.final_mass_flux.revision ==
                             second.piso.final_flux_revision,
                     rank,
                     "nondivisible product advances BE then BDF2 collectively");
    if (!runtime)
      std::cerr << "rank " << rank << " runtime="
                << static_cast<unsigned>(runtime.code) << "/"
                << runtime.detail << " stage=" << second.failed_stage
                << " phase=" << runtime_phase
                << '\n';
  }
  passed &= expect(run_mass_flow_product(rank), rank,
                   "mass-flow resolver advances collectively");
  passed &= expect(
      run_unsupported_candidate_boundary_product(rank), rank,
      "static, total-state, and NSCBC boundary-absolute cases reject before legacy pressure-only publication");
  passed &= expect(run_implicit_enthalpy_product(rank), rank,
                   "product reports the forced implicit enthalpy endpoint");
  passed &= expect(run_local_donor_product(rank), rank,
                   "incompatible forced local-donor flux is rejected by a "
                   "coupled candidate/terminal gate, rolls back committed "
                   "flow state exactly, and retires the fatal controller "
                   "ticket");
  passed &= expect(
      run_temporal_method_fallback_product(rank), rank,
      "BDF2 low-base failure falls back to BE within one attempt");
  passed &= expect(
      run_warm_start_lifecycle_product(rank), rank,
      "accepted C2 seed is available and consumed on the next C1 while an "
      "independent same-state cold path suppresses it coherently");
  passed &= expect(run_immersed_product(rank), rank,
                   "immersed product advances with compact remote donors");
  passed &= expect(run_rank_change_product(rank), rank,
                   "product restart advances across 1->2->4->1 ranks");
  passed &= expect(run_simple_coupling_compile(rank), rank,
                   "SIMPLE coupling freezes collectively");
  passed = collective(passed, MPI_COMM_WORLD);
  MPI_Finalize();
  return passed ? 0 : 1;
}
