C-----------------------------------------------------------------------
C   COAST1: settle to 28000, then fresh time averages through 38000
C-----------------------------------------------------------------------
false                           /use if number blocks > number of cores
C-----------------------------------------------------------------------
C   INPUT/OUTPUT DIRECTORIES:
C-----------------------------------------------------------------------
/Decomp/
/Restart/
/Restart/
C-------------------------------------------------------------------
C   COMMAND PARAMETERS
C   Restart may be remapped automatically when mpirun -np differs from active mesh.
true                            /form
true true /restart, write to restart
false                            /restart_regrid_replace_mesh
6000 /lstep
500 /stepsave
1000 /stepplot
C-----------------------------------------------------------------------
C   LES CONTROL - LARGE EDDY SIMULATION
C-----------------------------------------------------------------------
vreman                          /dyn_stress_piomelli,dyn_stress_lilly smagorinsky,dyn_energy
false                           /dyn_restart_stress
0.12                            /smagorinsky constant Cs
false false 0 /turbstat,turbread,statistics_start_step
false                           /reserved legacy flag; keep false
C-----------------------------------------------------------------------
C... Some miscellaneous variables
C-----------------------------------------------------------------------
2                               /itnum
10000000                        /integration time
10000000                        /nskip
C-----------------------------------------------------------------------
C... CFL variables
C-----------------------------------------------------------------------
0.2 0.3 0.25                    /cflmin,cflmax,cfl_ok
C-----------------------------------------------------------------------
C... Wall clock execution time in hours
C-----------------------------------------------------------------------
200                             /TEXEC
C-----------------------------------------------------------------------
C   PHYSICS CONTROL - PHYSICAL PROPERTIES
C-----------------------------------------------------------------------
1.18 0.663 1.0E-5 1.0e+05       /rhoair, rhojet, visco, pressure
C-----------------------------------------------------------------------
C   Variable
C-----------------------------------------------------------------------
  50    1.0e-04                 /u-vel 1 icycl(),tol()
  50    1.0e-04                 /v-vel 2 icycl(),tol()
  50    1.0e-04                 /w-vel 3 icycl(),tol()
  500   1.0e-04                 /dp    4 icycl(),tol()
  20    1.0e-03                 /f     5 icycl(),tol()
  20    1.0e-04                 /h     6 icycl(),tol()
C-----------------------------------------------------------------------
C... Combustion variables
C-----------------------------------------------------------------------
methane                         /fuel
jl4 /reaction mechanism
! fuel/mechanism is loaded here; burn and PDF switches are set below
0  0                            /CO2_content,H2_content - Volume fraction
!
true true vode /burn,chemkin format,solver
C-----------------------------------------------------------------------
C... PDF equation variables
C-----------------------------------------------------------------------
true /solve for pdf
4 /number of stochastic fields
true true /read_pdf,write_pdf
true                            /stochastic noise reduction
  10    1.0e-04                 /icycl(),tol()
C-----------------------------------------------------------------------
C   Species output: mass fraction or mole fraction
C-----------------------------------------------------------------------
mass_fraction                   /mole_fraction; mass_fraction
C-----------------------------------------------------------------------
C   Global heat release rate
C-----------------------------------------------------------------------
false /read global heat release rate: true/false
C-----------------------------------------------------------------------
C   Radiative Heat Transfer
C-----------------------------------------------------------------------
false                           /Include radiate heat transfer: true/false
C-----------------------------------------------------------------------
C   Ignition
C-----------------------------------------------------------------------
false                           /Ignite mixture: true/false
C-----------------------------------------------------------------------
C   Phase Averaging
C-----------------------------------------------------------------------
false false                     /phase averaging; read restart file: true/false
400  0.2872750                  /frequency from PSD; phase averaging start time
C-----------------------------------------------------------------------
C... Compressible flow
C-----------------------------------------------------------------------
true                            /compressible flow
C-----------------------------------------------------------------------
C... Effective boundary
C-----------------------------------------------------------------------
0.28 0.005266                   /Sigma, Characteristic length
! Focused migration fixture: source core plus representative extension records.
IEM manual 2.0 /micro_mixing_model, cd_mode, cd_value
&time_control
 dt_cap_enabled = .true.
 dt_cap = 9.9d-6
/
chemistry_backend 'name/with!punctuation' ! preserve quoted text
