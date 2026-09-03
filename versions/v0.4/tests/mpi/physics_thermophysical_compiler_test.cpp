// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_physics.hpp"

#include "physics_input_detail.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string_view>

namespace {

using namespace hundun::v04;

bool expect(bool condition, int rank, std::string_view description) {
  if (!condition) {
    std::cerr << "rank " << rank << " FAIL: " << description << '\n';
  }
  return condition;
}

bool all_true(bool local) {
  const int input = local ? 1 : 0;
  int output = 0;
  return MPI_Allreduce(&input, &output, 1, MPI_INT, MPI_MIN,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         output != 0;
}

bool identical_u64(std::uint64_t value) {
  std::uint64_t minimum = 0U;
  std::uint64_t maximum = 0U;
  return MPI_Allreduce(&value, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         MPI_Allreduce(&value, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         minimum == maximum;
}

std::uint64_t packed_status(Status status) {
  return (static_cast<std::uint64_t>(status.code) << 32U) |
         static_cast<std::uint64_t>(status.detail);
}

SpeciesThermophysicalSpec constant_species(std::string_view name,
                                           double molecular_weight,
                                           double dimensionless_cp,
                                           double viscosity,
                                           double conductivity) {
  SpeciesThermophysicalSpec species;
  species.stable_name.assign(name.data(), name.size());
  species.molecular_weight = molecular_weight;
  species.temperature_switch = 1000.0;
  species.nasa7_low[0U] = dimensionless_cp;
  species.nasa7_high[0U] = dimensionless_cp;
  species.transport_law = TransportLaw::constant;
  species.viscosity_reference = viscosity;
  species.conductivity = conductivity;
  return species;
}

ValidatedModel valid_model() {
  ValidatedModel model;
  model.fingerprint = 0x54a8c3e77612019bULL;
  model.thermophysics.data_file = "thermophysics.d";
  model.thermophysics.minimum_temperature = 200.0;
  model.thermophysics.maximum_temperature = 3000.0;
  model.thermophysics.temperature_relative_tolerance = 1.0e-12;
  model.thermophysics.maximum_temperature_iterations = 80U;
  model.thermophysics.closed_mass_relative_tolerance = 1.0e-12;
  model.thermophysics.maximum_closed_mass_iterations = 32U;
  model.thermophysics.maximum_closed_mass_relative_step = 0.25;
  model.thermophysics.species.push_back(
      constant_species("A", 20.0, 3.5, 1.0e-5, 0.020));
  model.thermophysics.species.push_back(
      constant_species("B", 40.0, 4.0, 2.0e-5, 0.030));
  model.transported_scalars.push_back(
      TransportedScalarSpec{"A", TransportedScalarRole::species});
  model.transported_scalars.push_back(
      TransportedScalarSpec{"tracer",
                            TransportedScalarRole::passive_scalar});
  model.transported_scalars.push_back(
      TransportedScalarSpec{"dye", TransportedScalarRole::passive_scalar});
  return model;
}

struct PublishedOutputs {
  ThermophysicalSpec spec;
  ThermodynamicsPlan thermodynamics;
  TransportPlan transport;
};

bool published_unchanged(const PublishedOutputs& outputs,
                         PlanFingerprint spec_fingerprint,
                         PlanFingerprint thermodynamics_fingerprint,
                         PlanFingerprint transport_fingerprint) {
  return detail::thermophysical_spec_fingerprint(outputs.spec) ==
             spec_fingerprint &&
         outputs.thermodynamics.fingerprint() ==
             thermodynamics_fingerprint &&
         outputs.transport.fingerprint() == transport_fingerprint &&
         outputs.spec.data_file == "thermophysics.d" &&
         outputs.spec.species.size() == 2U &&
         outputs.spec.species[0U].stable_name == "A" &&
         outputs.spec.maximum_temperature_iterations == 80U;
}

enum class Mutation : std::uint8_t {
  species,
  nasa,
  controls,
  path,
  invalid_spec,
  bad_case_fingerprint,
  transported_mapping,
  passive_name,
  passive_order,
  passive_count
};

void mutate_model(ValidatedModel& model, Mutation mutation, int rank,
                  int size, int mutation_rank) {
  const bool selected =
      mutation == Mutation::transported_mapping
          ? (rank == mutation_rank || (size >= 4 && rank == size - 1))
          : rank == mutation_rank;
  if (!selected) {
    return;
  }
  switch (mutation) {
    case Mutation::species:
      model.thermophysics.species[0U].molecular_weight =
          size == 1 ? 0.0 : 20.125;
      break;
    case Mutation::nasa:
      if (size == 1) {
        model.thermophysics.species[0U].nasa7_low[0U] = 0.5;
      } else {
        model.thermophysics.species[0U].nasa7_low[0U] = 3.75;
        model.thermophysics.species[0U].nasa7_high[0U] = 3.75;
      }
      break;
    case Mutation::controls:
      model.thermophysics.maximum_temperature_iterations =
          size == 1 ? 0U : 81U;
      break;
    case Mutation::path:
      model.thermophysics.data_file =
          size == 1 ? "" : "rank-local-thermophysics.d";
      break;
    case Mutation::invalid_spec:
      model.thermophysics.minimum_temperature = -1.0;
      break;
    case Mutation::bad_case_fingerprint:
      model.fingerprint = 0U;
      break;
    case Mutation::transported_mapping:
      model.transported_scalars[0U].stable_name =
          rank == mutation_rank ? "missing_species" : "also_missing";
      break;
    case Mutation::passive_name:
      model.transported_scalars[1U].stable_name = "renamed_tracer";
      break;
    case Mutation::passive_order:
      std::swap(model.transported_scalars[1U],
                model.transported_scalars[2U]);
      break;
    case Mutation::passive_count:
      model.transported_scalars.push_back(
          TransportedScalarSpec{"extra_passive",
                                TransportedScalarRole::passive_scalar});
      break;
  }
}

double dimensionless_h(const std::array<double, 7U>& coefficients,
                       double temperature) {
  const double t2 = temperature * temperature;
  const double t3 = t2 * temperature;
  const double t4 = t3 * temperature;
  return coefficients[0U] * temperature + coefficients[1U] * t2 * 0.5 +
         coefficients[2U] * t3 / 3.0 +
         coefficients[3U] * t4 * 0.25 +
         coefficients[4U] * t4 * temperature * 0.2 + coefficients[5U];
}

bool test_canonicalized_success(int rank) {
  ValidatedModel model = valid_model();
  constexpr double accepted_a6_mismatch = 1.0e-6;
  model.thermophysics.species[0U].nasa7_high[5U] =
      accepted_a6_mismatch;
  const double original_high_a6 =
      model.thermophysics.species[0U].nasa7_high[5U];

  PublishedOutputs outputs;
  ThermophysicalCompileDiagnostics diagnostics;
  diagnostics.lowest_failing_rank = 91;
  const Status status = ThermophysicalCompiler::load_and_compile(
      MPI_COMM_WORLD, model, outputs.spec, outputs.thermodynamics,
      outputs.transport, &diagnostics);

  bool passed = expect(static_cast<bool>(status), rank,
                       "within-tolerance high-a6 mismatch compiles");
  passed &= expect(identical_u64(packed_status(status)), rank,
                   "canonicalized success status is identical");
  passed &= expect(model.thermophysics.species[0U].nasa7_high[5U] ==
                           original_high_a6 &&
                       outputs.spec.species[0U].nasa7_high[5U] == 0.0,
                   rank,
                   "publication canonicalizes high a6 without mutating input");
  const auto& published = outputs.spec.species[0U];
  const double switch_temperature = published.temperature_switch;
  const double low_h =
      dimensionless_h(published.nasa7_low, switch_temperature);
  const double high_h =
      dimensionless_h(published.nasa7_high, switch_temperature);
  passed &= expect(low_h == high_h && std::isfinite(low_h) &&
                       outputs.thermodynamics.kernel() ==
                           ThermodynamicsKernel::constant_cp &&
                       diagnostics.lowest_failing_rank == -1,
                   rank,
                   "canonical spec has exactly continuous switch enthalpy");
  passed &= expect(
      identical_u64(detail::thermophysical_spec_fingerprint(outputs.spec)) &&
          identical_u64(outputs.thermodynamics.fingerprint()) &&
          identical_u64(outputs.transport.fingerprint()),
      rank, "canonical spec and plans have collective identities");
  return all_true(passed);
}

bool test_success(int rank, PublishedOutputs& outputs) {
  const ValidatedModel model = valid_model();
  ThermophysicalCompileDiagnostics diagnostics;
  diagnostics.lowest_failing_rank = 91;
  const Status status = ThermophysicalCompiler::load_and_compile(
      MPI_COMM_WORLD, model, outputs.spec, outputs.thermodynamics,
      outputs.transport, &diagnostics);

  const PlanFingerprint spec_fingerprint =
      detail::thermophysical_spec_fingerprint(outputs.spec);
  bool passed = expect(static_cast<bool>(status), rank,
                       "valid typed model compiles without I/O");
  passed &= expect(identical_u64(packed_status(status)), rank,
                   "success status is identical on every rank");
  passed &= expect(spec_fingerprint != 0U &&
                       identical_u64(spec_fingerprint),
                   rank, "published thermophysical spec is identical");
  passed &= expect(outputs.thermodynamics.fingerprint() != 0U &&
                       identical_u64(outputs.thermodynamics.fingerprint()),
                   rank, "published thermodynamics plan is identical");
  passed &= expect(outputs.transport.fingerprint() != 0U &&
                       identical_u64(outputs.transport.fingerprint()),
                   rank, "published transport plan is identical");
  passed &= expect(outputs.spec.species.size() == 2U &&
                       outputs.spec.species[0U].stable_name == "A" &&
                       outputs.spec.species[1U].stable_name == "B" &&
                       outputs.thermodynamics.species_count() == 2U &&
                       outputs.thermodynamics.independent_species_count() ==
                           1U &&
                       outputs.thermodynamics.dependent_species_index() ==
                           1U &&
                       outputs.thermodynamics.kernel() ==
                           ThermodynamicsKernel::constant_cp &&
                       outputs.transport.kernel() == TransportKernel::constant,
                   rank, "typed N-1 species mapping publishes expected plans");
  passed &= expect(diagnostics.lowest_failing_rank == -1, rank,
                   "successful compile reports no failing rank");
  return all_true(passed);
}

bool test_failure(Mutation mutation, std::string_view description, int rank,
                  int size, PublishedOutputs& outputs) {
  const PlanFingerprint retained_spec =
      detail::thermophysical_spec_fingerprint(outputs.spec);
  const PlanFingerprint retained_thermodynamics =
      outputs.thermodynamics.fingerprint();
  const PlanFingerprint retained_transport = outputs.transport.fingerprint();
  ValidatedModel model = valid_model();
  const int mutation_rank = size > 1 ? 1 : 0;
  mutate_model(model, mutation, rank, size, mutation_rank);
  const PlanFingerprint unchanged_case_fingerprint =
      0x54a8c3e77612019bULL;

  ThermophysicalCompileDiagnostics diagnostics;
  diagnostics.lowest_failing_rank = 91;
  const Status status = ThermophysicalCompiler::load_and_compile(
      MPI_COMM_WORLD, model, outputs.spec, outputs.thermodynamics,
      outputs.transport, &diagnostics);

  const int expected_lowest = mutation_rank;
  bool passed = expect(!status && status.code == StatusCode::invalid_plan,
                       rank, description);
  passed &= expect(identical_u64(packed_status(status)), rank,
                   "collective failure status is bitwise identical");
  passed &= expect(diagnostics.lowest_failing_rank == expected_lowest &&
                       identical_u64(static_cast<std::uint64_t>(
                           diagnostics.lowest_failing_rank + 1)),
                   rank, "diagnostics report the lowest failing rank");
  passed &= expect(published_unchanged(
                       outputs, retained_spec, retained_thermodynamics,
                       retained_transport),
                   rank, "collective rejection preserves all prior outputs");
  if (mutation != Mutation::bad_case_fingerprint) {
    passed &= expect(model.fingerprint == unchanged_case_fingerprint, rank,
                     "typed mutation does not rely on case-fingerprint change");
  }
  return all_true(passed);
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  PublishedOutputs outputs;
  bool passed = test_success(rank, outputs);
  passed &= test_canonicalized_success(rank);
  passed &= test_failure(Mutation::species,
                         "rank-local species mutation rejects collectively",
                         rank, size, outputs);
  passed &= test_failure(Mutation::nasa,
                         "rank-local NASA mutation rejects collectively",
                         rank, size, outputs);
  passed &= test_failure(Mutation::controls,
                         "rank-local control mutation rejects collectively",
                         rank, size, outputs);
  passed &= test_failure(Mutation::path,
                         "rank-local data path mutation rejects collectively",
                         rank, size, outputs);
  passed &= test_failure(Mutation::invalid_spec,
                         "invalid local spec rejects collectively", rank, size,
                         outputs);
  passed &= test_failure(Mutation::bad_case_fingerprint,
                         "zero case fingerprint rejects collectively", rank,
                         size, outputs);
  passed &= test_failure(
      Mutation::transported_mapping,
      "rank-local transported-species mapping rejects collectively", rank,
      size, outputs);
  if (size > 1) {
    passed &= test_failure(
        Mutation::passive_name,
        "rank-local legal passive-scalar name rejects collectively", rank,
        size, outputs);
    passed &= test_failure(
        Mutation::passive_order,
        "rank-local legal passive-scalar order rejects collectively", rank,
        size, outputs);
    passed &= test_failure(
        Mutation::passive_count,
        "rank-local legal passive-scalar count rejects collectively", rank,
        size, outputs);
  }

  passed = all_true(passed);
  if (rank == 0 && passed) {
    std::cout << "v0.4 thermophysical compiler MPI tests passed\n";
  }
  MPI_Finalize();
  return passed ? 0 : 1;
}
