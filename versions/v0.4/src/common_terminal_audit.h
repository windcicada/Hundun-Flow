/* SPDX-License-Identifier: Apache-2.0 */
#ifndef HF_COAST_COMMON_TERMINAL_AUDIT_H
#define HF_COAST_COMMON_TERMINAL_AUDIT_H

#ifdef __cplusplus
extern "C" {
#endif

int hf_coast_common_terminal_cell_v1(
    double rho, double rho_eos, double rho_accepted, double rho_previous,
    double volume, double bdf_a0, double bdf_a1, double bdf_a2,
    double flux_x_minus, double flux_x_plus, double flux_y_minus,
    double flux_y_plus, double flux_z_minus, double flux_z_plus,
    double pressure_perturbation, int closed_mass,
    double drho_dp_at_fixed_h_y, double *eos_residual,
    double *continuity_residual, double *mass_contribution,
    double *volume_contribution, double *absolute_pressure_perturbation,
    double *compressibility_pressure_moment,
    double *compressibility_weight);

int hf_coast_common_terminal_finalize_v1(
    int closed_mass, double global_mass, double global_volume,
    double closed_mass_target, double global_compressibility_pressure_moment,
    double global_compressibility_weight,
    double global_max_absolute_pressure_perturbation,
    double boundary_closure_residual, double *closed_mass_residual,
    double *gauge_residual);

int hf_coast_common_terminal_outlet_v1(double absolute_pressure,
                                       double target_absolute_pressure,
                                       double *gauge_residual);

int hf_coast_common_terminal_accept_v1(
    double eos_residual, double continuity_residual,
    double closed_mass_residual, double gauge_residual, double eos_tolerance,
    double continuity_tolerance, double closed_mass_tolerance,
    double gauge_tolerance);

#ifdef __cplusplus
}
#endif

#endif
