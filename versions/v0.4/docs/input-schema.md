<!-- SPDX-License-Identifier: Apache-2.0 -->

# HUNDUN-FLOW v0.4 case-input schema

Status: Task 3 cold-path schema, 2026-08-13

## Authority and directory layout

`case.json` is the only configuration authority. A case is a flat directory:

```text
case-root/
  case.json
  mesh-focus.d
  inlet-profile.d
  body.stl
```

Referenced `.d` and STL files must be direct children of the case root. Absolute paths,
parent traversal, nested paths, missing files, non-regular files, wrong suffixes, duplicate
references (including hard-link aliases), and symbolic links are rejected. Rank 0 opens every
file relative to one case-root directory descriptor, verifies the opened descriptor with
`fstat`, and reads or hashes through that same descriptor. A `.d`
file contains bulk arrays or tables only; it cannot override JSON configuration.

Rank 0 alone canonicalizes the case root, opens `case.json`, validates the schema, and reads
referenced bytes. It broadcasts a bounded typed model. Other ranks neither parse JSON nor
touch case files. `CaseCompiler::load_and_compile` requires an initialized MPI runtime and a
valid communicator.

## Complete Task 3 schema

Every object is closed: unknown, duplicate, missing, or incorrectly typed keys are errors.
The sole optional block is `turbulence`; omitting it selects
`vreman_wall_function`. The schema is deliberately small until later capability tasks add
their own registered input.

```json
{
  "schema_version": 1,
  "units": "SI",
  "mesh": {
    "kind": "uniform",
    "data_files": [],
    "stl_file": null
  },
  "flow": {
    "model": "single_phase_low_mach_compressible",
    "pressure_closure": "local_absolute_pressure_drho_dp",
    "reacting": false
  },
  "solver": {
    "coupling": "PISO",
    "pressure_correctors": 2
  },
  "turbulence": {
    "model": "vreman_wall_function"
  },
  "time": {
    "control": "adaptive_flow"
  }
}
```

Accepted values are:

| Key | Accepted value |
| --- | --- |
| `schema_version` | unsigned integer `1` |
| `units` | `SI` |
| `mesh.kind` | `uniform`, `tensor_stretched` |
| `mesh.data_files[]` | unique direct-root filenames ending in `.d` |
| `mesh.stl_file` | `null` or one unique direct-root filename ending in `.stl` |
| `flow.model` | `single_phase_low_mach_compressible` |
| `flow.pressure_closure` | `local_absolute_pressure_drho_dp` |
| `flow.reacting` | `false` |
| `solver.coupling` | `PISO` |
| `solver.pressure_correctors` | unsigned integer `2` |
| `turbulence.model` | optional block; default `vreman_wall_function`; explicit `vreman_wall_function`, `wale`, `none` |
| `time.control` | `fixed`, `adaptive_flow`, `adaptive_acoustic` |

`vreman_wall_function` is the production default chosen by generated examples and case
templates. Verification cases may explicitly select `none`; retained WALE cases explicitly
select `wale`. When the optional turbulence block is absent, the compiler resolves and hashes
the versioned `vreman_wall_function` default as part of the typed model.

The following selections are rejected in v0.4: constant density, body-fitted or multiblock
geometry, AMR/nonmatching refinement, SIMPLE, PIMPLE, a pressure-corrector count other than
two, reacting flow, shock physics, and turbulence models outside the three values above.

## Cold-path limits

Limits are checked before an unbounded allocation or collective transfer:

| Resource | Limit |
| --- | ---: |
| `case.json` bytes | 1 MiB |
| JSON nesting depth | 32 |
| referenced files | 256 |
| direct-root relative filename | 255 bytes |
| each referenced `.d` or STL file | 64 MiB |
| canonical typed MPI wire payload | 128 KiB |

Referenced files are hashed through a fixed-size streaming buffer; their contents are not
placed in the model wire. Exceeding a limit returns a deterministic rejected input status on
all ranks.

## Canonical model and fingerprint

The broadcast model contains only typed enum values, normalized relative filenames, an
optional normalized STL filename, and a nonzero 64-bit plan fingerprint. The fingerprint is
computed in declared order from:

1. a versioned semantic-contract identifier covering schema version, SI units, low-Mach flow,
   local-absolute-pressure EOS closure, nonreacting physics, PISO coupling, and two pressure
   correctors;
2. geometry, turbulence, and time-control enum values;
3. normalized reference type and direct-root name; and
4. the exact bytes and byte count of every referenced `.d` and STL file.

Changing referenced data therefore changes the fingerprint without broadcasting the data.
JSON formatting and object-key order do not affect it because only the validated typed model
is hashed. Every rank validates the bounded wire before publishing its `ValidatedModel`.

## Field registration boundary

`FieldRegistry::declare_field` assigns deterministic IDs in registration order and rejects
empty names, duplicate names, unsupported component/ghost counts, ID exhaustion, and mutation
after freeze. `freeze_for_test()` exists only to test snapshot mechanics on a synthetic
registry. Task 3 does not freeze the production `FieldSchema`; IBM, turbulence, equation,
solver, workspace, and diagnostic capabilities must register first. Task 18 performs the one
production freeze.

No case name, geometry name, measurement station, or benchmark-specific field is compiled
into product source. New runtime cases use this input contract and capability registration,
not source edits.
