# zmath — current state and plan

This document records where the project stands at this point and what remains.
It is the live plan. The two roadmaps under `docs/superpowers/plans/` are
historical and superseded; they are kept for provenance only and must not be
executed.

## What this project is

A dependency-free Zig 0.16.0 IEEE-754 and scalar mathematics library for `f16`,
`f32`, and `f64`, with a portable `@Vector` layer over the operations, vendored
CORE-MATH C kernels, and a sibling harness that proves the results against an
independent MPFR/Arb oracle.

Two colocated repositories:

- `zmath` — the runtime library. No MPFR, GMP, FLINT, Arb, BLIS, Eigen, C++, or
  libc math dependency.
- `zmath-testing` — the harness. Owns oracle generation, immutable corpus,
  replay, and benchmarks. This is the only acceptance authority.

## Current state

### Runtime

| Property | Value |
| --- | --- |
| Top-level operations | 99 |
| Operations accepting `@Vector` | 96 (`nan`, `printAnotherMessage`, `add` do not) |
| Vector dispatch lines in `root.zig` | 102 |
| Public vector entry points | 417 across six modules |
| Source files under `src/` | 194 |
| Commits | 729 |
| Tracked files | 238 |

The vector layer is portable Zig `@Vector`, lowered by the compiler per target.
There is no SIMD dispatch, no hand-written per-ISA kernel, no runtime CPU
detection, and no SSE2 fallback. The build still requires `+fma` for x86/x86_64
C kernels and supports little-endian aarch64.

### Vector coverage

`src/vector.zig` holds the 55 elementwise float-to-float operations and
re-exports five sibling modules:

| Module | Operations | Result shape |
| --- | --- | --- |
| `core` | shared helpers | — |
| `struct_family` | 13 result-struct operations | lane-ordered `[N]Result` |
| `intarg` | 11 operations with a non-float extra argument | `@Vector`, or `[N]Result` for the power family |
| `narrow` | 6 C23 narrowing operations plus ternary `fma` | `@Vector(N, f32)` from `@Vector(N, f64)` |
| `ieee_family` | 16 IEEE, extrema, integer, and boolean operations | `@Vector` of float, `i32`, `i64`, or `bool` |

**The one shape that is not a `@Vector`:** Zig 0.16 rejects a struct as a `@Vector`
element type (`expected integer, float, bool, or pointer for the vector element
type`). Every operation whose scalar kernel returns a result struct therefore
returns a lane-ordered array. This buys no throughput, because the kernel still
runs once per lane; the value is a uniform lane-parallel API. Do not read any
vector result as a correctness or performance claim.

### Acceptance

All 315 manifest descriptors are `status = "correct"` with `max_ulp = 0`. The
evaluator split is 150 `exact_bits`, 126 `mpfr_arb`, 39 `algebraic`.

Corpus coverage is not whole-domain. Only the 61 exhaustive binary16 unary
campaigns cover their entire domain. Every `f32` and `f64` row is exact agreement
on its stored selected tuples, and multi-argument or integer-domain binary16 rows
are selected too. A `correct` row is a statement about the stored tuples, never
a theorem about the operation.

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
cannot be vectorized without being pinned. It also pins three known traps: the
`op + type + v` naming collision, the binary16 kernel trap where root `sin` binds
`f16` to `scalar.sinh` and not `scalar.sinhh`, and the exclusion of
struct-returning operations from the float-vector surface.

## Verification

Run from the repository root. The qualified profile is
`x86_64-linux-gnu.2.38` / `haswell` / `ReleaseSafe`.

```text
# runtime
zig build test
gate

# harness, from the sibling directory
zig build test
gate
zig build replay-corpus   -Dtarget=x86_64-linux-gnu.2.38 -Dcpu=haswell -Doptimize=ReleaseSafe
zig build replay-vector   -Dtarget=x86_64-linux-gnu.2.38 -Dcpu=haswell -Doptimize=ReleaseSafe
```

Last recorded results, all green:

| Check | Result |
| --- | --- |
| Runtime gate | passed, all seven checks |
| Runtime tests | 604/604 |
| Harness gate | passed, all seven checks |
| Harness tests | 36/36 |
| Scalar corpus replay | 315 descriptors, 63,160,084 records |
| Vector corpus replay | 165 descriptors, 3,179,250 records |

`gate` is the devenv script wrapping `zig-gate`. It runs format, build, test,
AST check, and three linters. The harness gate sets `GATE_EXCLUDE` so gitignored
benchmark artifacts are not linted as source.

## Known limits, stated so they are not mistaken for gaps

1. No whole-domain proof for any `f32` or `f64` descriptor, and none for
   multi-argument or integer-domain binary16. This is a property of the evidence
   design, not a pending task.
2. Correct rounding of finite power, root, and pi-scaled results is a target
   requirement that is explicitly not claimed as satisfied.
3. The public helpers `classify`, `ulpDistance`, the four `scaleB` modes, the
   four `roundIntegral` modes, and integer `add` have no manifest row.
4. `generated/verification.json` is a dated 2026-09-09 snapshot with no
   checked-in producer. Its counters are not current. It is deliberately left
   unedited; the corpus README says so.
5. The eight-gigabyte corpus under `test/corpus/generated/` is gitignored and
   must be regenerated with `generate-corpus` on a new machine.
6. Operations returning a struct return a lane-ordered array, not a `@Vector`.
   See above.

## What remains

Nothing is in flight. The working copy is clean at the commit that added the
new-family surface guard.

Decisions that are genuinely open, in the order they are worth taking:

1. **Whole-domain evidence for `f32` and `f64`.** The largest remaining
   evidence gap. It is a budget question, not a design question: the campaigns
   exist, they are simply enormous.
2. **Manifest rows for the uncovered public helpers.** A small, well-bounded
   piece of work that would close a named coverage hole.
3. **The planned non-scalar modules.** `complex`, the linear-algebra `Vec*` value
   types, `quaternion`, `matrix`, `linalg`, `geometry`, and `algebra` are
   specified in `API.md` and entirely unimplemented. The implemented
   `zmath.vector` module is the `@Vector` layer and is not those types.
4. **Remove the Zig-init template leftovers.** `printAnotherMessage` and `add`
   are project-template boilerplate in a float-only library, currently
   documented as exports in `API.md`. Removing them is a public API change.
