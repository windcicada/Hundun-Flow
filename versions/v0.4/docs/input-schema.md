<!-- SPDX-License-Identifier: Apache-2.0 -->

# HUNDUN-FLOW v0.4 case-input schema

Status: Task 5 closed cold-path schema, 2026-08-13

## Authority and directory layout

`case.json` is the only configuration authority. Referenced `.d` tables and STL geometry
are bytes-only inputs; they cannot override JSON. Every referenced file must be a unique,
regular direct child of the case root. Absolute paths, parent traversal, nested paths,
incorrect suffixes, missing files, symbolic links, and hard-link aliases are rejected.

Rank 0 alone opens the case-root descriptor, `case.json`, `.d`, and STL files. It parses and
normalizes the closed schema, hashes the typed mesh and referenced bytes, and broadcasts a
bounded wire-v2 typed model. Other ranks never inspect the filesystem or parse JSON. Failed
compilation leaves the caller's output model unchanged on every rank.

## Complete closed schema

Every object below is closed. Unknown, duplicate, missing, or incorrectly typed keys are
errors. Only `turbulence` is optional; omission selects `vreman_wall_function`. JSON
`schema_version` remains `1`; wire version `2` is an internal typed-model format.

```json
{
  "schema_version": 1,
  "units": "SI",
  "mesh": {
    "kind": "tensor_stretched",
    "domain": {
      "lower": [0.0, 0.0, 0.0],
      "upper": [2.0, 1.0, 1.0]
    },
    "exact_cells": [64, 32, 32],
    "base_spacing": [0.05, 0.05, 0.05],
    "minimum_spacing": [0.01, 0.01, 0.01],
    "max_growth_ratio": 1.2,
    "focus_regions": [
      {
        "lower": [0.2, 0.2, 0.2],
        "upper": [0.8, 0.8, 0.8],
        "target_spacing": [0.02, 0.02, 0.02]
      }
    ],
    "limits": {
      "max_global_cells": 1000000,
      "max_memory_bytes_per_rank": 1073741824
    },
    "data_files": ["profile.d"],
    "stl_file": "body.stl"
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
  "turbulence": {"model": "vreman_wall_function"},
  "time": {"control": "adaptive_flow"}
}
```

All mesh numbers must be finite IEEE-754 values. Each domain upper coordinate must be greater
than its lower coordinate. Every spacing is positive, `max_growth_ratio >= 1`, every cell
count is a positive signed-32-bit representable integer, and both limits are positive.
When `exact_cells` is present, its checked product cannot exceed `max_global_cells`.

### Uniform Cartesian mesh

For `kind: "uniform"`:

- `exact_cells` is a required three-integer vector;
- `base_spacing` is `null`;
- `focus_regions` is empty;
- `max_growth_ratio` is exactly `1.0`; and
- each `minimum_spacing` component is no greater than the corresponding domain span divided
  by the exact cell count.

### Tensor-stretched Cartesian mesh

For `kind: "tensor_stretched"`, `base_spacing` is a required positive vector and each
component is at least `minimum_spacing`. `exact_cells` may be `null` or a three-integer vector;
when present it fixes the global counts while base spacing and focus regions define the target
distribution. Each focus region is a closed object with `lower`, `upper`, and
`target_spacing` three-vectors. Target spacing must be component-wise within
`[minimum_spacing, base_spacing]`.

A focus box must overlap the domain with positive extent in all three axes. Partly external
boxes are clipped to the domain; wholly external or merely touching boxes are rejected.
After clipping, regions are sorted lexicographically by lower bound, upper bound, and target
spacing, then exact duplicates are removed. Negative zero is normalized to positive zero.
Consequently JSON focus order, duplicate entries, equivalent clipped bounds, and `-0.0` do
not change the typed model or fingerprint.

## Remaining accepted values

| Key | Accepted value |
| --- | --- |
| `units` | `SI` |
| `mesh.data_files[]` | unique direct-root names ending in `.d` |
| `mesh.stl_file` | `null` or a unique direct-root name ending in `.stl` |
| `flow.model` | `single_phase_low_mach_compressible` |
| `flow.pressure_closure` | `local_absolute_pressure_drho_dp` |
| `flow.reacting` | `false` |
| `solver.coupling` | `PISO` |
| `solver.pressure_correctors` | unsigned integer `2` |
| `turbulence.model` | optional block; `vreman_wall_function`, `wale`, or `none` |
| `time.control` | `fixed`, `adaptive_flow`, or `adaptive_acoustic` |

Body-fitted, multiblock, AMR/nonmatching refinement, constant-density flow, SIMPLE/PIMPLE,
reacting physics, other corrector counts, and unlisted turbulence or time controls are
rejected.

## Cold-path bounds and fingerprint

| Resource | Limit |
| --- | ---: |
| `case.json` | 1 MiB |
| JSON nesting depth | 32 |
| focus regions before normalization | 1024 |
| referenced files | 256 |
| direct-root filename | 255 bytes |
| each referenced file | 64 MiB |
| typed MPI wire payload | 128 KiB |

The deterministic nonzero fingerprint hashes a versioned semantic contract followed by the
complete canonical typed mesh: kind; domain bounds; exact/base presence flags and values;
minimum spacing; growth ratio; normalized focus count and every focus bound/target; both
limits; turbulence and time enums; normalized reference kinds/names; and exact referenced
file bytes plus byte counts. Floating-point values are hashed by their canonical IEEE-754
binary64 bit patterns. JSON whitespace, object-key order, and normalization-equivalent focus
input therefore do not affect identity, while any typed mesh or referenced-byte change does.

`FieldRegistry::freeze_for_test()` remains a synthetic-test mechanism only. No production
case compile freezes the field schema; all capabilities must register before Task 18 performs
the single production freeze. Product code contains no source-case names or benchmark-specific
field declarations.
