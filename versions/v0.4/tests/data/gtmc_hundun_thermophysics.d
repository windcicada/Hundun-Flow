# HUNDUN development-state transfer for COAST GTMC cold CH4/air.
# Molecular weights and NASA7 coefficients are copied from the audited COAST
# therm.d snapshot. HUNDUN intentionally retains Ru=8314.46261815324.
HUNDUN_THERMOPHYSICS_V1
temperature_bounds 200 6000
temperature_inversion 1e-12 64
closed_mass_newton 1e-12 32 0.2
species_count 3
species CH4
molecular_weight 16.04308
temperature_switch 1000
nasa7_low 5.14825732 -0.013700241 4.93749414e-05 -4.91952339e-08 1.70097299e-11 -10245.3222 -4.63322726
nasa7_high 1.911786 0.0096026796 -3.38387841e-06 5.3879724e-10 -3.19306807e-14 -10099.2136 8.48241861
transport_sutherland 1.0945738056963247e-05 295 167.42379305670073 0.72
end_species
species O2
molecular_weight 31.9988
temperature_switch 1000
nasa7_low 3.78245636 -0.00299673416 9.84730201e-06 -9.68129509e-09 3.24372837e-12 -1063.94356 3.65767573
nasa7_high 3.66096065 0.000656365811 -1.41149627e-07 2.05797935e-11 -1.29913436e-15 -1215.97718 3.41536279
transport_sutherland 2.018407740068663e-05 295 136.02421948058128 0.72
end_species
species N2
molecular_weight 28.0134
temperature_switch 1000
nasa7_low 3.53100528 -0.000123660988 -5.02999433e-07 2.43530612e-09 -1.40881235e-12 -1046.97628 2.96747038
nasa7_high 2.95257637 0.0013969004 -4.92631603e-07 7.86010195e-11 -4.60755204e-15 -923.948688 5.87188762
transport_sutherland 1.7625721726908212e-05 295 111.96444353041214 0.72
end_species
end
