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
| Non-scalar module source | 9597 lines across seven files |
| Source files under `src/` | 208 |
| Local tests | 814 |

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
| Runtime tests | 814/814 |
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
| Manifest descriptors | `test/corpus/manifest.toml` | 337 at `schema_version = 5` |
| Publication digest | `generated/COMPLETE` | enforced by `corpus_set.load` |
| Corpus records | the regenerated publication | 72,035,318 |
| Record bytes | the regenerated publication | 5,164,903,241 |
| `replay-corpus` | the regenerated publication | exit 0, 337 descriptors |

The publication is regenerated and republished at 337 descriptors, and
`replay-corpus` in this repository now exits 0 with `descriptors=337
records=72035318`. Immediately before that commit the committed publication
was 315 rows against a 337-row manifest, so a fresh checkout could not verify
its own corpus at all: the loader rejected it with `error.IncompleteSelection`.

The digest authority is `generated/COMPLETE`, which records the
`metadata_sha256` against which `corpus_set.load` recomputes the SHA-256 of
`generated/metadata.txt` and refuses the publication on disagreement. There is
no verification file beside it. A hand-maintained counter file has no
producer, so every field in it is a number nobody can reproduce and every field
goes stale silently.

Full regeneration costs about 21 minutes for 5.2 GB, and the replay that
follows takes about 2.5 minutes, both measured on this host at 337
descriptors. They are costs, not coverage claims.

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

1. No whole-domain proof for any `f32` or `f64` descriptor, and none for
   multi-argument or integer-domain binary16. This is a property of the evidence
   design, not a pending task.
2. Correct rounding of finite power, root, and pi-scaled results is a target
   requirement that is explicitly not claimed as satisfied.
3. The public helpers `classify`, `ulpDistance`, `scaleB`, and `roundIntegral` are now covered by
   the manifest under an evaluator class of their own, `definitional`, rather than `exact_bits`,
   `mpfr_arb`, or `algebraic`. They are bit-level operations whose correct results are definitional:
   no oracle computes them, so a row for one is a contract lock and not independent numerical
   evidence. The separation is enforced rather than conventional — the exact-bits entry point
   refuses a definitional category, and a `metadata.txt` row prints `evaluator=definitional`.
4. `scaleB` and `roundIntegral` exist only for `f16`, and `API.md` declares no
   signature for either family — it names them only in the `zmath-public-helpers`
   benchmark description. The `f` and `d` forms are an open contract question,
   not a known gap.
5. The non-scalar modules have no independent oracle. Their accuracy rests on
   algebraic identities and on IEEE behaviour, not on the corpus.
6. There is no verification file in the sibling, and there is deliberately no committed JSON of
   run counters. A counter with no producer is a number nobody can reproduce. `generated/COMPLETE`
   plus `generated/metadata.txt` are the digest authority; every descriptor, record, and byte
   figure comes from a regenerated publication, not from a dated snapshot.
7. Operations returning a struct return a lane-ordered array, not a `@Vector`.
   See above.

## What remains

### Resolved in this cycle

1. **Whole-domain evidence for `f32` and `f64`.** Closed as infeasible rather than left open. The
   harness sustains roughly 62,000 records per second, which produced the publication in 17
   minutes. Exhaustive binary32 unary is 4,294,967,296 inputs per operation, about 19 hours each and
   roughly 48 days for the ~60 of them; exhaustive binary64 is about 1.8e19 inputs, roughly 9,000
   years each. The arithmetic is recorded in the scalar contract so it is not re-proposed.
   `correct` keeps meaning exact agreement on the stored tuples.

2. **Manifest rows for the uncovered public helpers.** Done. The manifest is 337 descriptors at
   `schema_version = 5`, and the publication is regenerated and replayed to match: 72,035,318
   records in 5,164,903,241 bytes, `replay-corpus` exit 0. The four public helpers account for 12
   rows, three formats each, under the `definitional` evaluator class precisely because no oracle
   computes them; ten more `definitional` rows cover the eight-operand algebra and quaternion struct
   operations, for 22 in the class. See known limit 3 for what that class does and does not claim.

3. **The `f` and `d` forms of `scaleB` and `roundIntegral`.** Done, and the decision recorded: they
   are implemented, and the four `RoundMode` values are reached as a field value the way `fromfp`
   does rather than as separate rows. See known limit 4 for why the contract previously named no
   signature for either family.

5. **The house gate did not build the corpus producer.** Fixed. `zig build` covers format, build,
   test, AST check and the linters, but the producer was only built by `zig build
   oracle-generator`, so a green gate did not mean the corpus could be produced. A schema change
   that missed the producer passed the gate and then failed silently at generation, producing the
   old 315-row publication from a new 323-row manifest. A `published-tools` step now depends on
   every published tool and every tool-owning test step, and the gate runs it under the qualified
   profile before `zig-gate`.

### Still open

4. **The `UnitQuat*.raw` privacy clause** cannot be met as written, because Zig has no private
   struct fields and a field `raw` and a method `raw()` cannot both exist. Recorded in `API.md`;
