// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "../support/turbulence_fixture.hpp"

#include <mpi.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <limits>
#include <string_view>

namespace {

using namespace hundun::v04;
using namespace hundun::v04::test;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

TurbulenceOwnedField gradient_field(const TurbulenceFixture& fixture,
                                    FieldId field, RevisionToken revision,
                                    StorageIdentity storage) {
  TurbulenceOwnedField result =
      make_turbulence_field(field, fixture.patch.cells, revision, storage);
  const std::size_t scalar_size = result.storage.size();
  result.storage.resize(9U * scalar_size, 0.0);
  result.view.base = result.storage.data();
  result.view.components = 9U;
  result.view.component_stride = scalar_size;
  for (std::size_t cell = 0U; cell < fixture.gradients.size(); ++cell) {
    const std::int32_t x = static_cast<std::int32_t>(
        cell % static_cast<std::size_t>(fixture.patch.cells.x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(fixture.patch.cells.x);
    const std::int32_t y = static_cast<std::int32_t>(
        yz % static_cast<std::size_t>(fixture.patch.cells.y));
    const std::int32_t z = static_cast<std::int32_t>(
        yz / static_cast<std::size_t>(fixture.patch.cells.y));
    for (std::uint8_t component = 0U; component < 9U; ++component) {
      result.view.unchecked({x, y, z}, component) =
          fixture.gradients[cell].value[component];
    }
  }
  return result;
}

bool same_bytes(const std::vector<double>& left,
                const std::vector<double>& right) {
  return left.size() == right.size() &&
         std::memcmp(left.data(), right.data(),
                     left.size() * sizeof(double)) == 0;
}

bool test_vreman_invariants() {
  bool passed = true;
  VelocityGradient shear;
  shear.value[1U] = 3.0;
  double viscosity = -1.0;
  passed &= expect(vreman_kinematic_viscosity(
                       shear, {0.1, 0.2, 0.3}, 0.07, viscosity) &&
                       viscosity == 0.0,
                   "Vreman is exactly zero for a rank-one shear gradient");

  VelocityGradient diagonal;
  diagonal.value[0U] = 1.0;
  diagonal.value[4U] = -2.0;
  diagonal.value[8U] = 0.5;
  const Real3 widths{0.1, 0.2, 0.3};
  const double beta_0 = widths.x * widths.x;
  const double beta_1 = widths.y * widths.y * 4.0;
  const double beta_2 = widths.z * widths.z * 0.25;
  const double invariant =
      beta_0 * beta_1 + beta_0 * beta_2 + beta_1 * beta_2;
  const double expected = 0.07 * std::sqrt(invariant / 5.25);
  passed &= expect(vreman_kinematic_viscosity(diagonal, widths, 0.07,
                                              viscosity) &&
                       std::abs(viscosity - expected) < 1.0e-15,
                   "Vreman anisotropic invariant matches independent oracle");
  return passed;
}

bool test_static_model_binding() {
  TurbulencePlanSpec spec;
  TurbulenceFixture fixture;
  bool passed = expect(fixture.initialize(spec),
                       "default turbulence is Vreman plus wall function");
  passed &= expect(fixture.plan.kind() ==
                           TurbulenceKind::vreman_wall_function &&
                       fixture.plan.subgrid_kind() == SubgridKind::vreman &&
                       fixture.plan.wall_treatment() ==
                           WallTreatmentKind::equilibrium_wall_function,
                   "default plan freezes one subgrid and wall branch cold");
  for (VelocityGradient& gradient : fixture.gradients) {
    gradient.value[0U] = 0.7;
    gradient.value[4U] = -0.4;
    gradient.value[8U] = -0.3;
  }
  TurbulenceCertificate certificate;
  passed &= expect(fixture.plan.update(fixture.input(), fixture.effective.view,
                                       certificate) &&
                       certificate.valid(),
                   "Vreman publishes mu_eff through its sole authority");
  for (double value : fixture.effective.storage) {
    passed &= expect(std::isfinite(value) && value > 1.8e-5,
                     "nonzero Vreman invariant raises mu_eff");
  }
  const std::vector<double> span_result = fixture.effective.storage;
  TurbulenceOwnedField gradient_field = make_turbulence_field(
      3U, fixture.patch.cells, 105U, 205U);
  const std::size_t scalar_size = gradient_field.storage.size();
  gradient_field.storage.resize(9U * scalar_size, 0.0);
  gradient_field.view.base = gradient_field.storage.data();
  gradient_field.view.components = 9U;
  gradient_field.view.component_stride = scalar_size;
  for (std::size_t cell = 0U; cell < fixture.gradients.size(); ++cell) {
    const std::int32_t x = static_cast<std::int32_t>(
        cell % static_cast<std::size_t>(fixture.patch.cells.x));
    const std::size_t yz = cell /
                           static_cast<std::size_t>(fixture.patch.cells.x);
    const std::int32_t y = static_cast<std::int32_t>(
        yz % static_cast<std::size_t>(fixture.patch.cells.y));
    const std::int32_t z = static_cast<std::int32_t>(
        yz / static_cast<std::size_t>(fixture.patch.cells.y));
    for (std::uint8_t component = 0U; component < 9U; ++component) {
      gradient_field.view.unchecked({x, y, z}, component) =
          fixture.gradients[cell].value[component];
    }
  }
  TurbulenceUpdateInput field_input;
  field_input.density = as_const(fixture.density.view);
  field_input.molecular_viscosity = as_const(fixture.molecular.view);
  field_input.gradient_revision = gradient_field.view.revision;
  field_input.velocity_gradient_field = as_const(gradient_field.view);
  passed &= expect(
      fixture.plan.update(field_input, fixture.effective.view, certificate) &&
          fixture.effective.storage == span_result,
      "SoA gradient field path matches the legacy AoS oracle without copy");

  TurbulencePlanSpec none_spec;
  none_spec.kind = TurbulenceKind::none;
  TurbulenceFixture none;
  passed &= expect(none.initialize(none_spec) &&
                       none.plan.update(none.input(), none.effective.view,
                                        certificate) &&
                       none.effective.storage == none.molecular.storage,
                   "verification none path copies molecular viscosity exactly");

  CartesianGeometryPlan geometry;
  MeshPatch patch;
  ContributionRegistry contributions;
  const std::array<FieldId, 3U> fields{0U, 1U, 2U};
  TurbulencePlan rejected;
  TurbulencePlanSpec invalid = spec;
  invalid.kind = static_cast<TurbulenceKind>(255U);
  passed &= expect(CartesianGeometryCompiler::compile(
                       MPI_COMM_SELF, turbulence_mesh(GeometryKind::uniform),
                       {}, geometry, patch) &&
                       contributions.configure({fields.data(), fields.size()}) &&
                       TurbulencePlan::compile(MPI_COMM_SELF, invalid, geometry,
                                               patch, 2U, 41U,
                                               contributions, rejected)
                               .code == StatusCode::invalid_plan &&
                       rejected.fingerprint() == 0U,
                   "every unapproved turbulence combination rejects cold");
  return passed;
}

bool test_candidate_effective_viscosity_is_stateless_and_model_exact() {
  std::array<TurbulencePlanSpec, 3U> specs{};
  specs[0].kind = TurbulenceKind::none;
  specs[1].kind = TurbulenceKind::wale;
  specs[2].kind = TurbulenceKind::vreman_wall_function;

  bool passed = true;
  for (std::size_t model = 0U; model < specs.size(); ++model) {
    TurbulenceFixture fixture;
    passed &= expect(fixture.initialize(specs[model]),
                     "candidate turbulence model compiles");
    for (std::size_t cell = 0U; cell < fixture.gradients.size(); ++cell) {
      VelocityGradient& gradient = fixture.gradients[cell];
      gradient.value[0U] = 0.7 + 0.001 * static_cast<double>(cell);
      gradient.value[1U] = 0.2;
      gradient.value[4U] = -0.4;
      gradient.value[8U] = -0.3;
    }

    TurbulenceCertificate live_certificate;
    passed &= expect(static_cast<bool>(fixture.plan.update(
                         fixture.input(), fixture.effective.view,
                         live_certificate)),
                     "production turbulence update succeeds");
    const std::vector<double> live_before = fixture.effective.storage;
    const std::uint64_t count_before = fixture.plan.update_count();

    TurbulenceOwnedField gradient = gradient_field(
        fixture, static_cast<FieldId>(10U + model),
        static_cast<RevisionToken>(310U + model),
        static_cast<StorageIdentity>(410U + model));
    TurbulenceOwnedField candidate = make_turbulence_field(
        static_cast<FieldId>(20U + model), fixture.patch.cells,
        static_cast<RevisionToken>(320U + model),
        static_cast<StorageIdentity>(420U + model));
    fill(candidate, -7.0);
    const TurbulenceCandidateInput input{
        as_const(fixture.density.view), as_const(fixture.molecular.view),
        as_const(gradient.view), gradient.view.revision};
    TurbulenceCandidateCertificate candidate_certificate;
    passed &= expect(
        fixture.plan.evaluate_candidate_effective_viscosity(
            input, candidate.view, candidate_certificate) &&
            candidate_certificate.valid() &&
            candidate_certificate.matches(
                fixture.plan.fingerprint(), input, as_const(candidate.view)) &&
            candidate_certificate.no_state_mutation &&
            candidate_certificate.allocation_free &&
            candidate_certificate.effective_viscosity.field !=
                fixture.effective.view.field,
        "candidate certificate binds a non-authority output and exact views");
    passed &= expect(candidate.storage == live_before,
                     "none, WALE, and Vreman candidate values are bitwise exact");
    passed &= expect(fixture.effective.storage == live_before &&
                         fixture.plan.update_count() == count_before,
                     "candidate evaluation leaves live values and count unchanged");

    TurbulenceCertificate cached_certificate;
    passed &= expect(
        fixture.plan.update(fixture.input(), fixture.effective.view,
                            cached_certificate) &&
            fixture.plan.update_count() == count_before &&
            cached_certificate.state == live_certificate.state,
        "candidate evaluation leaves the production active cache intact");
  }
  return passed;
}

bool test_candidate_effective_viscosity_rejects_without_writes_and_recovers() {
  TurbulenceFixture fixture;
  bool passed = expect(fixture.initialize(TurbulencePlanSpec{}),
                       "candidate failure fixture compiles");
  for (VelocityGradient& gradient : fixture.gradients) {
    gradient.value[0U] = 0.7;
    gradient.value[1U] = 0.2;
    gradient.value[4U] = -0.4;
    gradient.value[8U] = -0.3;
  }
  TurbulenceCertificate live_certificate;
  passed &= expect(static_cast<bool>(fixture.plan.update(
                       fixture.input(), fixture.effective.view,
                       live_certificate)),
                   "live update establishes active cache before failures");
  const std::vector<double> live_before = fixture.effective.storage;
  const std::uint64_t count_before = fixture.plan.update_count();

  TurbulenceOwnedField gradient =
      gradient_field(fixture, 31U, 331U, 431U);
  TurbulenceOwnedField direct =
      make_turbulence_field(32U, fixture.patch.cells, 332U, 432U);
  const TurbulenceCandidateInput clean_input{
      as_const(fixture.density.view), as_const(fixture.molecular.view),
      as_const(gradient.view), gradient.view.revision};
  TurbulenceCandidateCertificate certificate;
  passed &= expect(static_cast<bool>(
                       fixture.plan.evaluate_candidate_effective_viscosity(
                           clean_input, direct.view, certificate)),
                   "direct clean candidate evaluation succeeds");
  const std::vector<double> direct_result = direct.storage;
  ConstFieldView changed_replica = as_const(direct.view);
  ++changed_replica.replica;
  ConstFieldView changed_base = as_const(direct.view);
  ++changed_base.base;
  ConstFieldView changed_storage = as_const(direct.view);
  ++changed_storage.storage_identity;
  passed &= expect(!certificate.matches(fixture.plan.fingerprint(), clean_input,
                                        changed_replica) &&
                       !certificate.matches(fixture.plan.fingerprint(),
                                            clean_input, changed_base) &&
                       !certificate.matches(fixture.plan.fingerprint(),
                                            clean_input, changed_storage),
                   "candidate certificate rejects swapped replica/base/storage");

  TurbulenceOwnedField candidate =
      make_turbulence_field(33U, fixture.patch.cells, 333U, 433U);
  fill(candidate, -7.25);
  const std::vector<double> sentinel = candidate.storage;
  auto rejects_unchanged = [&](const TurbulenceCandidateInput& input,
                               FieldView output,
                               StatusCode expected) {
    certificate = {};
    const Status status = fixture.plan.evaluate_candidate_effective_viscosity(
        input, output, certificate);
    return status.code == expected && !certificate.valid() &&
           same_bytes(candidate.storage, sentinel);
  };

  FieldView authority_output = candidate.view;
  authority_output.field = fixture.effective.view.field;
  passed &= expect(rejects_unchanged(clean_input, authority_output,
                                     StatusCode::invalid_plan),
                   "candidate cannot occupy the production authority FieldId");
  FieldView density_alias = fixture.density.view;
  density_alias.field = candidate.view.field;
  density_alias.revision = candidate.view.revision;
  passed &= expect(rejects_unchanged(clean_input, density_alias,
                                     StatusCode::invalid_plan),
                   "candidate output cannot alias density");
  FieldView molecular_alias = fixture.molecular.view;
  molecular_alias.field = candidate.view.field;
  molecular_alias.revision = candidate.view.revision;
  passed &= expect(rejects_unchanged(clean_input, molecular_alias,
                                     StatusCode::invalid_plan),
                   "candidate output cannot alias molecular viscosity");
  FieldView gradient_alias = gradient.view;
  gradient_alias.components = 1U;
  gradient_alias.field = candidate.view.field;
  gradient_alias.revision = candidate.view.revision;
  passed &= expect(rejects_unchanged(clean_input, gradient_alias,
                                     StatusCode::invalid_plan),
                   "candidate output cannot alias velocity gradient");

  TurbulenceCandidateInput stale = clean_input;
  ++stale.gradient_revision;
  passed &= expect(rejects_unchanged(stale, candidate.view,
                                     StatusCode::invalid_plan),
                   "stale gradient revision rejects before output writes");
  TurbulenceCandidateInput malformed = clean_input;
  malformed.velocity_gradient.components = 8U;
  passed &= expect(rejects_unchanged(malformed, candidate.view,
                                     StatusCode::invalid_plan),
                   "malformed gradient shape rejects before output writes");
  TurbulenceCandidateInput foreign_storage = clean_input;
  foreign_storage.velocity_gradient.storage_identity = 0U;
  passed &= expect(rejects_unchanged(foreign_storage, candidate.view,
                                     StatusCode::invalid_plan),
                   "unbound gradient storage rejects before output writes");
  FieldView foreign_shape = candidate.view;
  ++foreign_shape.interior.x;
  passed &= expect(rejects_unchanged(clean_input, foreign_shape,
                                     StatusCode::invalid_plan),
                   "foreign output shape rejects before output writes");
  FieldView unbound_output = candidate.view;
  unbound_output.storage_identity = 0U;
  passed &= expect(rejects_unchanged(clean_input, unbound_output,
                                     StatusCode::invalid_plan),
                   "unbound output storage rejects before output writes");

  const double saved_density = fixture.density.storage.back();
  fixture.density.storage.back() =
      std::numeric_limits<double>::quiet_NaN();
  passed &= expect(rejects_unchanged(clean_input, candidate.view,
                                     StatusCode::numerical_failure),
                   "non-finite density on a late cell causes bitwise zero writes");
  fixture.density.storage.back() = saved_density;
  const double saved_gradient = gradient.storage.back();
  gradient.storage.back() = std::numeric_limits<double>::infinity();
  passed &= expect(rejects_unchanged(clean_input, candidate.view,
                                     StatusCode::numerical_failure),
                   "non-finite gradient on a late cell causes bitwise zero writes");
  gradient.storage.back() = saved_gradient;

  passed &= expect(static_cast<bool>(
                       fixture.plan.evaluate_candidate_effective_viscosity(
                           clean_input, candidate.view, certificate)) &&
                       same_bytes(candidate.storage, direct_result),
                   "clean evaluation after poison equals direct clean result");
  TurbulenceCertificate cached_certificate;
  passed &= expect(static_cast<bool>(fixture.plan.update(
                       fixture.input(), fixture.effective.view,
                       cached_certificate)) &&
                       fixture.plan.update_count() == count_before &&
                       same_bytes(fixture.effective.storage, live_before),
                   "candidate failures leave live production cache untouched");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  bool passed = test_vreman_invariants();
  passed &= test_static_model_binding();
  passed &= test_candidate_effective_viscosity_is_stateless_and_model_exact();
  passed &=
      test_candidate_effective_viscosity_rejects_without_writes_and_recovers();
  MPI_Finalize();
  return passed ? 0 : 1;
}
