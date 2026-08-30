/* SPDX-License-Identifier: Apache-2.0 */
#include "common_terminal_audit.h"

#include <float.h>
#include <math.h>
#include <stddef.h>

static int finite_nonnegative(double value) {
  return isfinite(value) && value >= 0.0;
}

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
    double *compressibility_weight) {
  double unsteady;
  double flux_sum;
  double scale;
  double eos;
  double continuity;
  if (eos_residual == NULL || continuity_residual == NULL ||
      mass_contribution == NULL || volume_contribution == NULL ||
      absolute_pressure_perturbation == NULL ||
      compressibility_pressure_moment == NULL ||
      compressibility_weight == NULL || !(rho > 0.0) ||
      !(rho_eos > 0.0) || !(volume > 0.0) || !(bdf_a0 > 0.0) ||
      !isfinite(rho_accepted) || !isfinite(rho_previous) ||
      !isfinite(bdf_a1) || !isfinite(bdf_a2) ||
      !isfinite(flux_x_minus) || !isfinite(flux_x_plus) ||
      !isfinite(flux_y_minus) || !isfinite(flux_y_plus) ||
      !isfinite(flux_z_minus) || !isfinite(flux_z_plus) ||
      !isfinite(pressure_perturbation) ||
      (closed_mass && !(drho_dp_at_fixed_h_y > 0.0)) ||
      (closed_mass && !isfinite(drho_dp_at_fixed_h_y))) {
    return 1;
  }
  unsteady = volume *
             (bdf_a0 * rho + bdf_a1 * rho_accepted +
              bdf_a2 * rho_previous);
  flux_sum = (flux_x_plus - flux_x_minus) +
             (flux_y_plus - flux_y_minus) +
             (flux_z_plus - flux_z_minus);
  scale = fabs(volume * bdf_a0 * rho) +
          fabs(volume * bdf_a1 * rho_accepted) +
          fabs(volume * bdf_a2 * rho_previous) + fabs(flux_x_minus) +
          fabs(flux_x_plus) + fabs(flux_y_minus) + fabs(flux_y_plus) +
          fabs(flux_z_minus) + fabs(flux_z_plus);
  if (scale < DBL_MIN) scale = DBL_MIN;
  eos = fabs(rho - rho_eos) / fmax(1.0, fabs(rho_eos));
  continuity = fabs(unsteady + flux_sum) / scale;
  if (!finite_nonnegative(eos) || !finite_nonnegative(continuity)) return 1;

  *eos_residual = eos;
  *continuity_residual = continuity;
  *mass_contribution = volume * rho;
  *volume_contribution = volume;
  *absolute_pressure_perturbation = fabs(pressure_perturbation);
  *compressibility_pressure_moment =
      closed_mass ? volume * drho_dp_at_fixed_h_y * pressure_perturbation
                  : 0.0;
  *compressibility_weight =
      closed_mass ? volume * drho_dp_at_fixed_h_y : 0.0;
  return 0;
}

int hf_coast_common_terminal_finalize_v1(
    int closed_mass, double global_mass, double global_volume,
    double closed_mass_target, double global_compressibility_pressure_moment,
    double global_compressibility_weight,
    double global_max_absolute_pressure_perturbation,
    double boundary_closure_residual, double *closed_mass_residual,
    double *gauge_residual) {
  double mass;
  double gauge;
  if (closed_mass_residual == NULL || gauge_residual == NULL ||
      !isfinite(global_mass) || !(global_volume > 0.0) ||
      !isfinite(global_max_absolute_pressure_perturbation) ||
      global_max_absolute_pressure_perturbation < 0.0) {
    return 1;
  }
  if (closed_mass) {
    if (!(closed_mass_target > 0.0) ||
        !(global_compressibility_weight > 0.0) ||
        !isfinite(global_compressibility_pressure_moment)) {
      return 1;
    }
    mass = fabs(global_mass - closed_mass_target) / closed_mass_target;
    gauge = fabs(global_compressibility_pressure_moment /
                 global_compressibility_weight) /
            fmax(1.0, global_max_absolute_pressure_perturbation);
  } else {
    if (!finite_nonnegative(boundary_closure_residual)) return 1;
    mass = 0.0;
    gauge = boundary_closure_residual;
  }
  if (!finite_nonnegative(mass) || !finite_nonnegative(gauge)) return 1;
  *closed_mass_residual = mass;
  *gauge_residual = gauge;
  return 0;
}

int hf_coast_common_terminal_outlet_v1(double absolute_pressure,
                                       double target_absolute_pressure,
                                       double *gauge_residual) {
  double residual;
  if (gauge_residual == NULL || !isfinite(absolute_pressure) ||
      !isfinite(target_absolute_pressure)) {
    return 1;
  }
  residual = fabs(absolute_pressure - target_absolute_pressure) /
             fmax(1.0, fabs(target_absolute_pressure));
  if (!finite_nonnegative(residual)) return 1;
  *gauge_residual = residual;
  return 0;
}

int hf_coast_common_terminal_accept_v1(
    double eos_residual, double continuity_residual,
    double closed_mass_residual, double gauge_residual, double eos_tolerance,
    double continuity_tolerance, double closed_mass_tolerance,
    double gauge_tolerance) {
  if (!finite_nonnegative(eos_residual) ||
      !finite_nonnegative(continuity_residual) ||
      !finite_nonnegative(closed_mass_residual) ||
      !finite_nonnegative(gauge_residual) || !(eos_tolerance > 0.0) ||
      !(continuity_tolerance > 0.0) || !(closed_mass_tolerance > 0.0) ||
      !(gauge_tolerance > 0.0) || !isfinite(eos_tolerance) ||
      !isfinite(continuity_tolerance) ||
      !isfinite(closed_mass_tolerance) || !isfinite(gauge_tolerance)) {
    return -1;
  }
  return eos_residual <= eos_tolerance &&
                 continuity_residual <= continuity_tolerance &&
                 closed_mass_residual <= closed_mass_tolerance &&
                 gauge_residual <= gauge_tolerance
             ? 1
             : 0;
}
