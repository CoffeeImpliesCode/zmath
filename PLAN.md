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

## The harness is not present

**`zmath-testing` does not exist on this machine and is not recoverable from
here.** Verified absent: no `zmath-testing` directory under the home tree at any
depth, no `manifest.toml`, `corpus_runner.zig`, or `oracle_generator.zig`
anywhere, and no repository for it under the only authenticated GitHub account
(`CoffeeImpliesCode`), either directly or by global search.

Consequence: the recorded acceptance evidence in the sibling project — 315
descriptors, 63,160,084 records, and the `verification.json` digests — is
neither confirmed nor contradicted here. It is simply uncheckable. Treat it as
historical, not as current evidence.

To restore it, clone the sibling into the documented location and regenerate the
gitignored corpus:

```text
git clone <zmath-testing> ../zmath-testing
cd ../zmath-testing && zig build generate-corpus \
    -Dtarget=x86_64-linux-gnu.2.38 -Dcpu=haswell -Doptimize=ReleaseSafe
```

## Known limits, stated so they are not mistaken for gaps

1. No whole-domain proof for any `f32` or `f64` descriptor, and none for
   multi-argument or integer-domain binary16. This is a property of the evidence
   design, not a pending task.
2. Correct rounding of finite power, root, and pi-scaled results is a target
   requirement that is explicitly not claimed as satisfied.
3. The public helpers `classify`, `ulpDistance`, the four `scaleB` modes, the
   four `roundIntegral` modes, and integer `add` have no manifest row. `add` was
   in that list and has been removed as Zig-init template text, so the list is
   now `classify`, `ulpDistance`, and the eight rounding modes.
4. `scaleB` and `roundIntegral` exist only for `f16`, and `API.md` declares no
   signature for either family — it names them only in the `zmath-public-helpers`
   benchmark description. The `f` and `d` forms are an open contract question,
   not a known gap.
5. The non-scalar modules have no independent oracle. Their accuracy rests on
   algebraic identities and on IEEE behaviour, not on the corpus.
6. `generated/verification.json` in the sibling is a dated 2026-09-09 snapshot
   with no checked-in producer. Its counters are not current.
7. Operations returning a struct return a lane-ordered array, not a `@Vector`.
   See above.

## What remains

### Blocked on the missing harness

1. **Whole-domain evidence for `f32` and `f64`.** The largest remaining evidence
   gap. It is a budget question, not a design question: the campaigns exist in
   the sibling project, they are simply enormous, and the sibling is not here.
2. **Manifest rows for the uncovered public helpers.** A small, well-bounded
   piece of work that would close a named coverage hole. The manifest is a
   sibling file, so this cannot be started from this repository.

Both are blocked on the same missing prerequisite, and neither has a reachable
half: the deliverable for each is a harness artifact, not a runtime change.

### Open on this repository

3. **Commit the seven executed items.** The work is verified in the working copy
   and uncommitted.
4. **Decide whether `f` and `d` forms of `scaleB` and `roundIntegral` belong in
   the contract.** See known limit 4. That is a contract decision, and the
   contract is the thing that is under-specified.
5. **The `UnitQuat*.raw` privacy clause cannot be met as written.** Recorded in
   `API.md`; no further work unless the contract is amended.
