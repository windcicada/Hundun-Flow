# SPDX-License-Identifier: Apache-2.0
"""Finite-iteration flow accounting; strict physical defects remain reported.

Sum(R_mass) = storage + boundary - phase.
Sum(R_energy) = total storage + advective/work fluxes - exchange.
Frozen PDF accounting is derived in the internal V1.1.0 implementation record.
The envelope scales with actual FP64 storage/rates, never a fitted defect.
"""
import copy
import math
import sys

EPS = sys.float_info.epsilon


def check_composition(payload):
    for group in ('species', 'elements'):
        for row in payload['composition_balance'][group]:
            terms = [row[k] for k in ('temporal_rate', 'transport_outflow', 'pressure_outflow',
                                     'noise_source', 'mixing_source', 'chemistry_source')]
            defect = math.fsum(terms[:3] + [-v for v in terms[3:]])
            explanation = [row[k] for k in ('advective_conversion', 'density_update',
                                            'transport_equation_residual')]
            scale = math.fsum(map(abs, terms + explanation))
            bound = row['storage_roundoff_bound'] + 32*EPS*scale
            unexplained = defect-math.fsum(explanation)
            assert all(math.isfinite(v) for v in terms+explanation+[bound, unexplained])
            assert abs(defect-row['defect']) <= bound, row
            assert abs(unexplained-row['unexplained_defect']) <= bound, row
            assert abs(unexplained) <= bound, row


def check_flow(payload):
    p = payload
    assert p['dt'] > 0 and p['mass_kg'] > 0
    mass = math.fsum([p['mass_bdf_rate_kg_s'], p['mass_outflow_kg_s'],
                      -p['phase_mass_input_kg_s']])
    energy_terms = [p['total_energy_bdf_rate_W'], p['enthalpy_outflow_W'],
                    p['kinetic_energy_outflow_W'], p['mechanical_pressure_work_outflow_W'],
                    p['statistical_enthalpy_outflow_W'], -p['statistical_enthalpy_source_W'],
                    -p['conductive_heat_input_W'], -p['viscous_work_input_W'],
                    -p['species_enthalpy_diffusion_input_W'], -p['phase_energy_input_W']]
    energy = math.fsum(energy_terms)
    mass_bound = 128*EPS*(2*p['mass_kg']/p['dt'] + abs(p['mass_bdf_rate_kg_s']) +
                          abs(p['mass_outflow_kg_s']) + abs(p['phase_mass_input_kg_s']))
    energy_bound = 128*EPS*(2*(abs(p['internal_energy_J'])+abs(p['kinetic_energy_J']))/p['dt'] +
                            math.fsum(map(abs, energy_terms)))
    for actual, report, residual, recorded, bound in (
        (mass, p['mass_balance_defect_kg_s'], p['continuity_equation_sum_kg_s'],
         p['mass_budget_identity_kg_s'], mass_bound),
        (energy, p['total_energy_balance_defect_W'], p['energy_equation_sum_W'],
         p['energy_budget_identity_W'], energy_bound)):
        assert all(math.isfinite(v) for v in (actual, report, residual, recorded, bound))
        assert abs(actual-report) <= bound
        assert abs(actual-residual) <= bound, (actual, residual, bound)
        assert abs((report-residual)-recorded) <= bound
    # Retain the established flow mass-inventory accuracy floor independently.
    assert abs(mass*p['dt'])/p['mass_kg'] < 1e-6
    check_composition(p)


def check_omission_counterexamples(payload):
    for change in (-1, 1):
        bad = copy.deepcopy(payload)
        row = bad['composition_balance']['species'][0]
        row['chemistry_source'] += change*max(1., abs(row['chemistry_source']))*.01
        try:
            check_composition(bad)
        except AssertionError:
            pass
        else:
            raise AssertionError('a missing or duplicated chemical source escaped accounting')
