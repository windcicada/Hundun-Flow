// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_app.hpp"
#include <algorithm>
#include <iomanip>
#include <ostream>

namespace hundun::v04 {

Status write_step_completion_failure(std::ostream& stream,
                                     const StepCompletionReport& completion) noexcept {
  try {
    stream << std::setprecision(17) << "step_failure_completion_v1";
    const auto status = [&](const char* key, Status value) {
      stream << ' ' << key << '=' << static_cast<unsigned>(value.code)
             << '/' << value.detail;
    };
    status("first", completion.first_failure.failure);
    stream << " first_attempt=" << completion.first_failure.attempt
           << " first_stage=" << completion.first_failure.stage
           << " first_dt=" << completion.first_failure.dt;
    status("last", completion.last_failure.failure);
    stream << " last_attempt=" << completion.last_failure.attempt
           << " last_stage=" << completion.last_failure.stage
           << " last_dt=" << completion.last_failure.dt;
    status("stop", completion.stop_reason);
    status("outcome", completion.outcome);
    stream << '\n';
    return stream ? Status{} : Status{StatusCode::io_failure, 10507U};
  } catch (...) {
    return {StatusCode::io_failure, 10507U};
  }
}

Status write_numerical_failure(std::ostream& stream,
    const NumericalFailureContext& failure) noexcept {
  try {
  if (!failure.valid) return {};
  stream << std::setprecision(17)
            << "numerical_failure_context_v1"
            << " field=" << static_cast<unsigned>(failure.field)
            << " field_id=" << failure.field_id
            << " status=" << static_cast<unsigned>(failure.failure.code)
            << '/' << failure.failure.detail
            << " stage=" << failure.stage
            << " global_cell=" << failure.global_cell
            << " global_index=" << failure.global_index.x << ','
            << failure.global_index.y << ',' << failure.global_index.z
            << " rank=" << failure.rank
            << " attempted_step=" << failure.attempted_step
            << " generation=" << failure.generation
            << " time=" << failure.time
            << " target_time=" << failure.target_time
            << " dt_before=" << failure.dt_before
            << " method_before="
            << static_cast<unsigned>(failure.method_before.order) << ','
            << failure.method_before.a0 << ',' << failure.method_before.a1
            << ',' << failure.method_before.a2
            << " origin_before="
            << static_cast<unsigned>(failure.origin_before)
            << " retry_proposed=" << (failure.retry_proposed ? 1 : 0)
            << " dt_after=" << failure.dt_after
            << " method_after="
            << static_cast<unsigned>(failure.method_after.order) << ','
            << failure.method_after.a0 << ',' << failure.method_after.a1
            << ',' << failure.method_after.a2
            << " origin_after=" << static_cast<unsigned>(failure.origin_after)
            << " failed_value=" << failure.failed_value
            << " allowed=" << failure.allowed_minimum << ','
            << failure.allowed_maximum
            << " p_abs=" << failure.pressure_absolute
            << " T_before=" << failure.temperature_before
            << " T_estimate=" << failure.temperature_estimate
            << " rho=" << failure.density_before << ','
            << failure.density_previous << ',' << failure.density_predicted
            << " h=" << failure.enthalpy_before << ','
            << failure.enthalpy_previous
            << " cp_before=" << failure.cp_before
            << " Y=";
  for (std::size_t index = 0U;
       index < std::min(failure.mass_fraction_count, failure.mass_fractions.size()); ++index) {
    if (index != 0U) stream << ',';
    stream << failure.mass_fractions[index];
  }
  stream << " Y_truncated="
            << (failure.mass_fractions_truncated ||
                failure.mass_fraction_count > failure.mass_fractions.size() ? 1 : 0)
            << " revisions=" << failure.enthalpy_accepted_revision << ','
            << failure.enthalpy_previous_revision << ','
            << failure.nonadvective_accepted_revision << ','
            << failure.nonadvective_previous_revision << ','
            << failure.mass_flux_accepted_revision << ','
            << failure.mass_flux_previous_revision << ','
            << failure.temperature_accepted_revision << ','
            << failure.conductivity_revision << ','
            << failure.velocity_gradient_revision << ','
            << failure.effective_viscosity_revision
            << " boundary_stencil="
            << (failure.physical_boundary_stencil ? 1 : 0) << ','
            << (failure.mpi_boundary_stencil ? 1 : 0) << ','
            << (failure.immersed_interface_cell ? 1 : 0)
            << " mass_div=" << failure.mass_divergence_accepted << ','
            << failure.mass_divergence_previous
            << " advection=" << failure.advection_accepted << ','
            << failure.advection_previous
            << " nonadvective=" << failure.nonadvective_accepted << ','
            << failure.nonadvective_previous
            << " predictor_terms=" << failure.conservative_history_value
            << ',' << failure.accepted_advection_delta << ','
            << failure.accepted_nonadvective_delta << ','
            << failure.previous_advection_delta << ','
            << failure.previous_nonadvective_delta << ','
            << failure.reconstructed_value
            << " rate_terms=" << failure.diffusion_accepted << ','
            << failure.pressure_work_accepted << ','
            << failure.viscous_dissipation_accepted << ','
            << failure.explicit_source_accepted << ','
            << failure.implicit_sink_accepted
            << " rate_complete="
            << (failure.rate_breakdown_complete ? 1 : 0)
            << " first_bad="
            << static_cast<unsigned>(failure.first_bad_contributor)
            << " tvd_envelope="
            << (failure.face_envelope_checked ? 1 : 0) << ','
            << (failure.face_envelope_valid ? 1 : 0) << ','
            << failure.maximum_face_envelope_violation << ','
            << failure.selected_face_value << ','
            << failure.selected_donor_minimum << ','
            << failure.selected_donor_maximum
            << " be_counterfactual=" << failure.counterfactual_be_density
            << ',' << failure.counterfactual_be_enthalpy << ','
            << (failure.counterfactual_be_admissible ? 1 : 0)
            << " runtime_region=" << static_cast<unsigned>(failure.runtime_region)
            << " inversion=" << static_cast<unsigned>(failure.inversion.outcome)
            << " inversion_accepted=" << failure.inversion.accepted
            << " inversion_iterations=" << failure.inversion.iterations
            << " inversion_h=" << failure.inversion.input_enthalpy
            << " inversion_h_bounds=" << failure.inversion.minimum_enthalpy
            << ',' << failure.inversion.maximum_enthalpy
            << " inversion_T_bracket=" << failure.inversion.lower_temperature
            << ',' << failure.inversion.upper_temperature
            << " inversion_residual=" << failure.inversion.residual
            << " inversion_roundoff_bound=" << failure.inversion.roundoff_bound
            << '\n';

    return stream ? Status{} : Status{StatusCode::io_failure, 10507U};
  } catch (...) {
    return {StatusCode::io_failure, 10507U};
  }
}

}  // namespace hundun::v04
