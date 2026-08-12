# Stage 1 runtime architecture

## Dependency direction and application flow

The Stage 1 dependency direction is:

```text
applications -> solver / mesh / sdk -> runtime
applications -> config
```

`hundun` parses the command line before MPI initialization. Rank zero loads and
strictly validates JSON, and the application broadcasts the resulting typed
configuration rather than JSON text. Validation-only modes exit before field
allocation. A numerical run then constructs the decomposition, concrete
uniform mesh, frozen field registry, storage, exchange plan, Halo instance,
and passive-scalar solver in dependency-safe lifetime order.

## MPI context and decomposition

`MpiEnvironment` requests `MPI_THREAD_FUNNELED` and owns initialization only
when necessary. `MpiContext` collectively duplicates a user-supplied
intracommunicator, sets `MPI_ERRORS_RETURN` on the owned duplicate, caches
relative rank, size, and thread level, and centralizes project MPI errors.
Communicator owners are destroyed before the environment. Cached values remain
local after finalization, while operations that call MPI reject an inactive
runtime.

The structured decomposition accepts a complete explicit process grid or
chooses a feasible grid using the global cell shape and periodicity. Owned
boxes are half-open. Local cells use `[0,n)` while lower and upper ghosts use
`[-g,0)` and `[n,n+g)`. Global cell IDs are independent of the process grid.

## Fields and Halo exchange

The registry owns typed descriptors and is frozen before storage allocation.
Stage 1 allocates cell averages. A `FieldView` is a checked, borrowed view over
one storage allocation; callers do not cache it across an operation, a time
step, or a storage lifetime/change boundary.

An immutable `ExchangePlan` describes all 26 face, edge, and corner regions at
an arbitrary Ghost width. Each `HaloExchange` owns a duplicate communicator,
reusable buffers, and request storage, so instances isolate their direction
tags by communicator context. Creation checks `MPI_TAG_UB`; payload and request
counts are split within MPI-3 `int` limits. Synchronous `exchange()` wraps the
single-in-flight `begin()`/`wait()` interface. Stage 1 does not fuse fields or
run a background progress thread.

## Numerical and file I/O contracts

The passive scalar uses conservative MUSCL reconstruction with an MC limiter
and SSPRK2 time integration. The application declares exactly two Float64
cell-average fields: persistent selected field `scalar`, and transient hidden
field `stage`, both with Ghost width two.

Restart v1 is a collective same-partition protocol. Every step directory has
`manifest.v1.bin`, one rank file per relative rank, and `COMPLETED`. The
manifest records the global topology, process grid, owned boxes, persistent
schema, byte sizes, and CRC-64/ECMA-182 values. Rank files and the manifest use
temporary-name/rename writes; rank zero writes the marker last. Reading first
validates the marker, manifest, every record, and the local rank checksum, then
commits staged values only after all ranks succeed. Version 1 supports only the
same rank count, process grid, owned partition, and persistent-field schema. It
does not repartition and does not claim power-loss durability.

Primitive legacy VTK output contains one selected Float64 cell-average scalar
per rank. It is a Stage 1 verification artifact, not a general parallel I/O
framework.

## Public interface limits

The C plugin v1 interface performs metadata discovery and compatibility
negotiation only. It is not a model lifecycle or model callback ABI. No MPI
handle, field view, allocator, model handle, or callback crosses it.

`UniformStructuredMesh` is a concrete Stage 1 baseline, not a general mesh
interface. The current `FieldView` has no liveness epoch, acquisition-time
read/write capability, or public unchecked kernel view. These limitations are
explicit; Stage 1 does not add placeholder abstractions for later physics.

## Mandatory inputs to a future Stage 2 plan

Stage 2 requires a new approved plan. Its required planning order is:

1. measurement and performance contracts;
2. field epoch, permission, checked-view, and kernel-view safety;
3. separated mesh topology, geometry, and boundary patches;
4. backend-neutral linear algebra;
5. variable-density operators;
6. PISO and Rhie-Chow coupling;
7. common correctness and performance gates.

That future plan must define project-owned, replaceable `LinearOperator`,
`Preconditioner`, and `LinearSolver` contracts without vendor types. It must
evaluate CPU reference, CPU optimized, and device backend categories; explicit
execution-space ownership, residency, transfer, lifetime, and error semantics;
matrix-free and batched operations; asynchronous Halo with capability-detected
device-aware MPI and a correct host-staging fallback; and mixed precision
followed by project-owned double-precision residual correction. Every backend
must pass the same residual, conservation, convergence,
decomposition-invariance, and failure contracts.

Future measurements include strong and weak scaling, allocated bytes per cell,
Halo payload and bandwidth, collective kind/count/bytes, I/O throughput, and
linear-solver iterations and residuals. Stage 1 sets no portable wall-clock
threshold and implements none of these future interfaces or backends.
