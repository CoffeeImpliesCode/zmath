# zmath — current state and plan

This document records where the project stands and what remains.
It is the live plan. The two roadmaps under `docs/superpowers/plans/` are
historical and superseded; they are kept for provenance only and must not be
executed.

## What this project is

A dependency-free Zig 0.16.0 IEEE-754 and scalar mathematics library for `f16`,
`f32`, and `f64`, with a portable `@Vector` layer over the operations, vendored
CORE-MATH C kernels, and a sibling harness that is meant to prove the results
against an independent MPFR/Arb oracle.

Two colocated repositories:

- `zmath` — the runtime library. No MPFR, GMP, FLINT, Arb, BLIS, Eigen, C++, or
  libc math dependency.
- `zmath-testing` — the harness. Owns oracle generation, immutable corpus,
  replay, and benchmarks. This is the only acceptance authority.

## Current state

### Runtime

| Property | Value |
| --- | --- |
| Public modules | 9, every path in `API.md` implemented |
| Top-level operations | 97 |
| Public root types and module aliases | 25 |
| Operations accepting `@Vector` | 96 (`nan` does not) |
| Public `@Vector` entry points | 417 across six modules |
| Source files under `src/` | 209 |
| Local tests | 856 |

The vector layer is portable Zig `@Vector`, lowered by the compiler per target.
There is no SIMD dispatch, no hand-written per-ISA kernel, no runtime CPU
detection, and no SSE2 fallback. The build requires `+fma` for x86/x86_64 C
kernels and supports little-endian aarch64; a big-endian target is rejected at
configure time and by a `_Static_assert` in the generated ABI guard.

### Non-scalar modules

| Module | Types | Notes |
| --- | --- | --- |
| `complex` | `Complexf`, `Complexd` | 54 suffixed functions, all C99/C23 principal-branch rules documented in a per-family table |
| `vector` | `Vec2f`..`Vec4d` | 158 methods; `eq` is bit equality with NaN normalization, `normalize` of zero returns the input unchanged |
| `quaternion` | `Quatf`, `Quatd`, `UnitQuatf`, `UnitQuatd` | `UnitQuat*` has no `raw()` method: Zig puts fields and decls in one namespace, so a field `raw` and a method `raw` cannot both exist, and Zig has no private fields |
| `matrix` | `Mat` generic, `Mat2f`..`Mat5d` | Storage is column-major, `at` is row-major, and the `init` literal is row-major. `Mat` is the one sanctioned public generic |
| `linalg` | `LuMat2f`..`LuMat4d`, `CholeskyMat2f`..`CholeskyMat4d` | Doolittle with partial pivoting; reuses the pivot thresholds `zmath.matrix` publishes rather than declaring a second spelling |
| `geometry` | 59 types | `PredicatePolicy` is explicit; `adaptive_exact` is exact integer arithmetic on a fixed `[132]u64` limb array |
| `algebra` | `Cl3`, `Cl4`, `Pga3`, `Pga4`, `Cga3`, `Cga4` in `f` and `d` | Comptime-generated blade product tables; published normalization table; every row implemented |

## How the work stays correct

Every vector lane is computed by the identical accepted scalar kernel, so lane
`i` is bit-for-bit the scalar result for lane `i`. That is what lets the vector
surface inherit the scalar `max_ulp = 0` contract and the immutable corpus with
no new acceptance evidence. A future change that re-derives a vector polynomial
instead of calling the scalar kernel breaks this and invalidates the argument.

`src/vector_surface_test.zig` is the guard. It drives every vectorized operation
through zmath's own root overload and compares lane by lane against the same
root call made on a scalar lane, and it derives the guarded set from
`std.meta.declarations(root)` rather than a hand-kept list, so an operation
cannot be vectorized without being pinned.

The non-scalar modules are held to algebraic identities rather than bit
patterns: associativity, distributivity, `reconstruct` against the original
matrix, `conjugate` inverting, `reverse` being an involution, the signed
Poincaré dual's sign rule, and round trips through `fromCoefficients`.

## Verification

Run from the repository root. The qualified corpus profile is
`x86_64-linux-gnu.2.38` / `haswell` / `ReleaseSafe`.

```text
# runtime
zig build test
gate
zig build run
```

Last recorded results, all green:

| Check | Result |
| --- | --- |
| Runtime tests | 856/856 |
| Runtime gate | passed |
| House gate invocation | `direnv exec . gate`, or `direnv exec . zig-gate` |

`gate` is the devenv script wrapping `zig-gate`. It runs format, build, test,
AST check, and three linters. On a host with a `/nix/store` the store script runs
it directly; on this host there is no `/nix/store` at all, so `.envrc`
synthesizes a `gate` shim into the gitignored `.devenv/host-bin`.

**The gate truncates its own linter output.** `zig-gate` pipes zbc's findings through `head -30`,
so any advisory past the thirtieth line is invisible in the gate result. It also counts findings
with `grep -c 'error('`, which is why a run with warnings still prints "found 0 potential
issue(s)". Neither number is a reliable census. To get the truth, run the linter directly on the
whole tree:

```text
direnv exec . zbc src/*.zig 2>&1 | grep -c 'warning('
```

That is the check to use when the question is "are there any advisories left", and it must be
zero. Six advisories in the accepted scalar kernels went unrecorded for exactly this reason.

## Acceptance evidence, re-verified

The sibling `zmath-testing` is present at `../zmath-testing`. Its published evidence is
re-verified on this host rather than merely inherited: the table below is what the tree can show
now, measured here.

| Check | Source of truth | Current value |
| --- | --- | --- |
| Manifest descriptors | `test/corpus/manifest.toml` | 347 at `schema_version = 5` |
| Publication digest | `generated/COMPLETE` | enforced by `corpus_set.load` |
| Corpus records | the regenerated publication | 72,035,642 |
| Record bytes | the regenerated publication | 5,164,943,119 |
| `replay-corpus`, a `zmath-testing` step | the regenerated publication | exit 0, every one of the 347 descriptors, on the scalar lane and on the scalar-backed `--vector` lane |

The publication is regenerated and republished at 347 descriptors, and the
sibling `replay-corpus` step verifies every descriptor and every record on the
scalar lane and on the scalar-backed vector lane. That step belongs to
`zmath-testing`, not to this repository: the `zmath` build graph has only `run`,
`test`, and `docs`.

Two earlier repairs are historical and stay in the record. The publication that
predates the 337-descriptor manifest carried 315 rows against it, so a fresh
checkout could not verify its own corpus at all: the loader rejected that
snapshot with `error.IncompleteSelection`. A `metadata.txt`, `replay.txt`, or
`benchmark.txt` describing a 315- or 337-descriptor snapshot is a receipt of
that older publication and is not evidence about the ten rows added afterwards.

The digest authority is `generated/COMPLETE`, which records the
`metadata_sha256` against which `corpus_set.load` recomputes the SHA-256 of
`generated/metadata.txt` and refuses the publication on disagreement. There is
no verification file beside it. A hand-maintained counter file has no
producer, so every field in it is a number nobody can reproduce and every field
goes stale silently.

That generation step reused 337 verified artifacts and produced only the 324 new
bounded records, and the staging run took 137.84 s on this host. That is the cost
of a bounded increment, not of a whole-domain regeneration, so it predicts
nothing about an exhaustive binary16 or binary32 campaign.
Replaying the whole 5,164,943,119-byte publication at 347 descriptors took
142.42 s on the scalar lane and 165.17 s on the scalar-backed vector lane. Both
are costs, not coverage claims. The scalar benchmark measured all 347
descriptors of that same publication at three repetitions each in 94.96 s.
The retained receipt contains only repetition zero for each descriptor.
That benchmark exceeded the approved binary32 smoke evaluation limit.
The execution-limit correction is recorded in step 4 below.

A digest recomputed from a DIFFERENT run still will not match, and that is not
a defect. Each metadata row embeds its `artifact=` path, which is that run's
own records root joined to the record file name, so the digest depends on where
the run wrote its records. Judge publication integrity by the digest `COMPLETE`
carries and the loader checks, never by a digest recomputed from another
output root or from another snapshot.

Bringing the harness up on this host required repairing it; see the two commits in
`../zmath-testing`. The oracle generator needs a flint that the host can load, and the store's
build cannot be, so a flint and MPFR are built for the host into a local prefix.

## Known limits, stated so they are not mistaken for gaps

1. No whole-domain proof for any `f32` or `f64` descriptor, or for
   multi-argument or integer-domain binary16.
   The current campaigns prove agreement only on their stored inputs.
2. Correct rounding of finite power, root, and pi-scaled results is a target
   requirement that is explicitly not claimed as satisfied.
3. The public helpers `classify`, `ulpDistance`, `scaleB`, and `roundIntegral` are now covered by
   the manifest under an evaluator class of their own, `definitional`, rather than `exact_bits`,
   `mpfr_arb`, or `algebraic`. They are bit-level operations whose correct results are definitional:
   no oracle computes them, so a row for one is a contract lock and not independent numerical
   evidence. The separation is enforced rather than conventional — the exact-bits entry point
   refuses a definitional category, and a `metadata.txt` row prints `evaluator=definitional`.
4. `scaleB` and `roundIntegral` support `f16`, `f32`, and `f64`.
   Their four rounding modes are typed input values, not separate manifest rows.
5. The non-scalar tranche has ten `definitional` rows for Cl(3) involutions
   and quaternion multiplication, with 176 stored records. Those are convention
   checks, and they are not relabelled as numerical evidence.
   Ten further numerical rows are five selected operations, each with one `f32`
   and one `f64` row: `complex.addf`/`addd`, `complex.mulf`/`muld`,
   `matrix.determinantMat2f`/`d`, `linalg.lu` plus `solveLu` on a 2x2 system,
   and the 2D orientation predicate under both policies, checked against an
   independent exact-rational reference. Those ten rows hold 324 stored records
   over bounded named cases.
   The rounding rule is the one each public graph performs, so it is not one
   rule for all ten: `complexMul`, `matrixDet2`, and the `fast` predicate round
   per primitive step, while `luSolve2` rounds its exact rational Cramer
   solution once per component.
   Other operations of those four modules stay outside the manifest, which bounds
   the tranche to its approved selection rather than leaving work open.
6. There is no verification file in the sibling, and there is deliberately no committed JSON of
   run counters. A counter with no producer is a number nobody can reproduce. `generated/COMPLETE`
   plus `generated/metadata.txt` are the digest authority; every descriptor, record, and byte
   figure comes from a regenerated publication, not from a dated snapshot.
7. Operations returning a struct return a lane-ordered array, not a `@Vector`.
   See above.

## Latest seven-step acceptance

The latest request refers to the seven next steps listed after the cluster cost estimate.
The earlier P8/P9 checklist completed a first tranche, not every item in that request.

1. **Published-tool gate.** Implemented.
   Both gate entry points build `published-tools` under the qualified profile before the house gate.
   This closes the stale-producer build gap.
2. **Harness advisories.** Cleared.
   The direct whole-tree zbc census reports zero findings.
3. **Harness template leftovers.** Removed.
   The demonstration now calls the mathematics API.
4. **Exhaustive binary32.** Not executed.
   Descriptor sharding and explicit raw-input ranges are both implemented, and
   `--f32-range=LOW:HIGH` is exercised on four small intervals for `abs` and
   narrowing `f16sqrtf32`, `[0,8)`, `[2147483644,2147483652)`,
   `[2139095036,2139095044)`, and `[4294967288,4294967296)`: two shards each,
   then assembly, receipts, and scalar replay, 32 generated inputs per
   operation.
   A missing flag never scans a domain; an assembled interval claims exhaustive
   coverage only when a proved tiling of `[0,2^32)` exists, and every smaller
   interval stays partial. The 63 eligible descriptors, which are those with one
   `float:f` input and arity one whatever the result format, would need
   270,582,939,648 records and 15,595,026,251,776 bytes, or 14.18359375 TiB, and
   copied assembly with retained workers needs 28.3671875 TiB, so no local run
   can hold it.
   Execution still needs explicit user permission, an exact domain, an
   operation selection, an interval partition, and approved cluster CPUs,
   allocation, and scratch.
   The approved smoke limit was 4,096 input evaluations per selected operation,
   including repetitions and reference checks.
   The later full-publication benchmark exceeded that limit.
   Its retained receipt reports 65,536 timed calls each for binary32 `abs`
   and narrowing `f16sqrtf32` in repetition zero alone.
   The retained first repetition exceeds 4,096 timed calls for 113 of 116
   descriptors with a binary32 input, including all 63 eligible unary descriptors.
   These timed-call counts exclude sample validation and warmup.
   The benchmark receipt remains unchanged as evidence of what ran.
   The earlier completion report did not disclose this execution-limit violation.
   The original smoke evaluation budget is exhausted.
   Further numerical evaluation of binary32 inputs needs explicit user approval.
5. **Six-module oracle coverage.** Complete under the approved scope, by user
   decision.
   The first tranche covered four Cl(3) involutions and quaternion
   multiplication in `f32` and `f64`, as ten `definitional` rows over 176 stored
   convention cases.
   Complex, matrix, linalg, and geometry now carry ten numerical rows in total:
   five selected operations, each with one `f32` and one `f64` row, over 324
   stored bounded records, verified against an independent exact-rational
   reference.
   The user accepted that selected coverage as sufficient, so this step is
   complete at the approved selection and is not pending work on the remaining
   operations of those four modules; those stay outside the manifest as a stated
   scope limit. The ten `definitional` rows remain convention checks, not proof
   for either complete non-scalar algebra module.
6. **Contract decisions.** The `f` and `d` forms of `scaleB` and `roundIntegral` are implemented and documented.
   The `UnitQuat*.raw` privacy clause is resolved by explicit user approval: the clause is amended
   rather than enforced.
   `raw` is public, `init` normalizes and validates, `fromAxisAngle` routes through `init`, and the
   consumer must preserve the unit-norm invariant.
   Zig has no private struct fields, and a field `raw` cannot coexist with a method `raw()`, so the
   API and the source state the caller-maintained precondition instead of claiming the type prevents
   a non-unit construction.
   `rotate` and composition take that precondition, so a forged non-unit value has no rotation
   guarantee.
7. **Housekeeping.** The old benchmark scratch decision and the retirement of `verification.json` are recorded.
   `COMPLETE` and `metadata.txt` remain the publication digest authority.
   The half-table generator pins the current metadata and the unchanged `asin` and `tanh` record digests.

The current publication has 347 descriptors, 72,035,642 records, and 5,164,943,119 record bytes.
Its successful replay resolves the selected-coverage and range-engineering steps. It does not
resolve the deferred exhaustive binary32 campaign, which needs explicit user approval.
The scalar benchmark measured all 347 descriptors and its receipt is retained,
but a timing is not coverage.
