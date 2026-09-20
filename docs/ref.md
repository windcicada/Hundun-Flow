# Internal reference history

### COAST functional reference

- Read-only location: `/home/wyf/code_dev/Coast_software`
- Observation date: 2026-08-13
- Reference scope: replacement capabilities and observed lifecycle only, including local
  absolute-pressure EOS, pressure density derivative, sound speed/Mach, NSCBC, SIMPLE-style
  transient iterations, Vreman/wall-function paths, flat `.d` inputs, STL scan, Restart,
  Visit, screen, and ICCG.
- Boundary: user-specified read-only reference; no general redistribution or license grant is
  presumed. Task 1 copies no COAST source. Old COAST Fortran may not be copied or translated,
  and no other COAST C++ reuse is authorized. Its SIMPLE loop and lack of periodic product
  capability are reference limitations, not HUNDUN v0.4 algorithm requirements.
- Sole exception: for Task 5, the user authorizes porting or adapting the `imb_mesh_y.cpp`
  scanning mathematics/method without provenance boilerplate. This narrow authorization does
  not extend to any other COAST C++ or to old COAST Fortran.
