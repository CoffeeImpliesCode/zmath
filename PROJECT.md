# zmath

`zmath` is a Zig 0.16.0 mathematics module for `f16`, `f32`, and `f64`.

## Implemented

All nine module paths named by the [API contract](API.md) exist and are exported from the root.

| Module | Source | Surface |
| --- | --- | --- |
| `zmath.ieee` | `src/ieee.zig` | The C23 IEEE primitives |
| `zmath.scalar` | `src/scalar.zig` | The accepted scalar kernels |
| `zmath.vector` | `src/vector.zig` plus six siblings and `src/vector_value.zig` | 417 `@Vector` entry points and the six `Vec*` value types |
| `zmath.complex` | `src/complex.zig` | `Complexf`, `Complexd`, and 54 suffixed functions |
| `zmath.quaternion` | `src/quaternion.zig` | `Quatf`, `Quatd`, `UnitQuatf`, `UnitQuatd` |
| `zmath.matrix` | `src/matrix.zig` | The `Mat` generic, eight size aliases, 19 module functions |
| `zmath.linalg` | `src/linalg.zig` | Six LU and six Cholesky factorization types, 10 module functions |
| `zmath.geometry` | `src/geometry.zig` | 2D/3D/4D point, vector, line, ray, segment, flat, ball, AABB families for `f` and `d`, 59 types |
| `zmath.algebra` | `src/algebra.zig` | `Cl3`, `Cl4`, `Pga3`, `Pga4`, `Cga3`, `Cga4` in `f` and `d`, with a published normalization table |

- Implementation and tests are flat under `src/`: `src/root.zig`, `src/ieee.zig`, `src/scalar*.zig`, `src/vector*.zig`, the six non-scalar modules, `src/*_test.zig`, and C kernels in `src/coremath_*.[ch]`.
- The root exports 97 top-level `pub fn` operations and 25 public types and module aliases.
- Each `@Vector` lane applies the identical accepted scalar kernel, so a vector result is bit-for-bit equal to the accepted scalar result for that lane and the `max_ulp = 0` contract is preserved. Zig forbids a struct as a `@Vector` element type, so the thirteen result-struct operations and the `pown`, `compoundn`, and `rootn` power family return a lane-ordered `[N]Result` array rather than a `@Vector`; that is a language limitation, and the array shape buys no throughput over a scalar loop.
- `zmath.vector` carries two unrelated layers that share the module path: the portable `@Vector` elementwise surface, and the named-component `Vec2f` through `Vec4d` value types from `src/vector_value.zig`. The elementwise operations use `v` or `vLanes` suffixes; no `Vec*` type name does.
- `Mat` is the one public structural generic. Every other public type in the library is concretely named, with private comptime generics behind it so the `f` and `d` forms cannot drift.
- `zmath.geometry` takes an explicit `PredicatePolicy`. The `adaptive_exact` path evaluates the predicate in exact integer arithmetic on a fixed `[132]u64` limb array, with no epsilon and no fallback to the `fast` answer.
- The runtime has no external oracle or math-library dependency. The vector layer uses portable Zig `@Vector` and is lowered by the compiler per target; it has no hand-written per-ISA kernels, no runtime dispatch, and no SSE2 fallback.
- `build.zig` makes the `zmath` module available to consumers and installs the `zig-out/bin/zmath` demonstration executable, not a compiled library artifact.
- C kernels require IEEE binary32/binary64 storage, a known little- or big-endian byte order, and `_BitInt(128)` or native `__int128`. The binary16 C ABI is bit-encoded `uint16_t`, with optional native conversion and direct portable RNE packing. C and Zig share a compile-time native-FMA policy and use bounded integer binary32/binary64 fusion where hardware fusion is absent. The x86/x86_64 `+fma` floor remains.

## Consumer interface

Add `b.dependency("zmath", .{ .target = target, .optimize = optimize }).module("zmath")`
to the consumer's imports. See [Use as a Zig dependency](API.md#use-as-a-zig-dependency) for a complete source example.
The build module carries the private C kernels and link settings.
The C headers do not define a supported standalone C product.

`zig build docs` generates the public module's documentation.
Serve the installed `docs` directory through HTTP.
The docs step stages regular Zig sources in the local build cache, including files whose
filesystem directory entries have an unknown type. This avoids incomplete source archives.
The fetched package includes `API.md`, this file, and the scalar corpus contract.
Local plans and setup files stay outside the package.

## API contract changes

Checked geometry distinguishes invalid supplied values from failed arithmetic.
Non-finite arguments and receiver fields return `error.NonFiniteInput`.
Non-finite results or required construction intermediates return `error.NonFiniteResult`.
Exhaustive error switches and narrower caller error sets must account for the new tag.
Identical rays now return `.coincident`, not their two-sided supporting `.line`.
Ball/ray intersections now return the closed common `segment`, not just two boundary points.
Subnormal directions normalize through a scaled vector; ordinary normalization retains division.
The intersection union has no new variants.

Numeric structs remain writable value types.
Constructors establish geometric and unit-quaternion invariants. Callers preserve those invariants after direct field writes.
No allocator, handle registry, teardown protocol, or runtime dispatch was added.

Owning declarations document status residues, sentinels, nullable failures, matrix input shapes, and fixed pivot thresholds.
The test-only matrix `entry_types` catalog and `acceptsEntryType` probe are removed.
Unsupported matrix entries fail at `Mat`; use the documented four entry types.
Root narrowing operations remain scalar-only. The portable-vector narrowing route remains `zmath.vector.narrow`.
These contract repairs do not change corpus acceptance or establish new accuracy or performance evidence.

## Big-endian lane publication

`zmath.vector` publishes logical lane order rather than the target's vector storage representation.
The publication rule lives in `src/vector_core.zig`, the `Vec` bridges in `src/vector_value.zig`
map components to and from lanes by logical lane index, and `src/vector_value_test.zig` and
`src/vector_surface_test.zig` are the receipt.

- The `Vec` bridges assign and read named components one lane at a time by logical lane index,
  without an intermediate array. That is the shape the parent commit already had: the experimental
  array-plus-`@bitCast` form of `toVector` was measured, changed nothing on `armeb-linux`, and was
  reverted, so `src/vector_value.zig` carries no change and no repair is claimed for it. No public
  method became an inline function, so the function-pointer form of every public operation is
  unchanged.

- On LLVM AArch64 big-endian only, indexed boolean result stores stage their verdicts in
  `@Vector(n, u8)` byte lanes, and the mask is then formed by one of three routes: a byte-aligned
  unsigned integer with one bit per lane outside `Debug`, a per-group lane reversal in `Debug` at
  power-of-two widths above one lane, and a comparison of the staged bytes at every other width.
  Every other target and result type keeps direct vector storage.
- `toVector publishes lane zero at the lowest address, on every target` reads the produced lanes
  the way a consumer does, by reinterpreting the vector's bytes rather than the returned
  register, and the boolean family is swept one-hot over runtime and comptime widths.

The publication holds everywhere it is executed, with two recorded limits. Native `Debug` and
native `ReleaseSafe` are both 956 of 956. The executed full-matrix `cross-test` result is `Debug`:
19 pass, 1 fail, 33 build_only across 53 targets, and the executed targets `loongarch64-linux`,
`aarch64_be-linux`, `thumbeb-linux`, `x86_64-linux`, `mips-linux`, `s390x-linux` and `riscv64-linux`
each report all 956 tests passed. A `ReleaseSafe` `cross-test` also ran, over `armeb-linux` and
`aarch64_be-linux`, with 2 pass, 0 fail and 0 build_only and both target logs reading all 956 tests
passed. That run is the first receipt for the integer publication route, which needs an `aarch64_be`
target and a non-`Debug` mode at once, and it fixes the register defect below as `Debug`-only. None
of this is whole-domain or oracle evidence.

The one failure is `armeb-linux` in `Debug` at 950 of 956. It is not a boolean-publication
failure: it is the pre-existing returned-register defect that
[`armeb-linux` Debug register-read defect](#armeb-linux-debug-register-read-defect) records, and
the same six tests fail at the parent commit. The other limit is that the three byte powers of
two above sixteen, and width 14, are swept by no test,
[Boolean mask publication at excluded widths](#boolean-mask-publication-at-excluded-widths)
records why, so no executed sweep covers every width the route keys on.

Two standalone-probe results disagree with that and are recorded as unresolved, not dismissed: a
`loongarch64-linux` `Debug` probe reports width 16 wrong on 7 of 16 marked lanes with runtime
operands and clean under `ReleaseSafe` in the same program, and an `aarch64_be-linux` `Debug`
probe reports every power-of-two width wrong for both ordering kernels with runtime operands.
Constant folding does not explain either one: the `loongarch64-linux` width 16 probe was re-run
with its operands built from a runtime loop index and reproduced the same failure. Both probes
differ from the project build in shape rather than in operands, since they omit the project's C
kernels and the `-ffreestanding -fno-math-errno -ffp-contract=off` cflags and compile a different
program. What is established is that boolean-mask publication on those big-endian targets is
sensitive to the shape of the surrounding program and not only to the width, so the sweep here is
a detector in the shipped configuration and not a guarantee for an arbitrarily shaped consumer
program.

## Boolean mask publication at excluded widths

Four `bool` result widths are swept by no test, for two distinct reasons, and each is left out of
a sweep rather than asserted weakly. Three of them are excluded because the compiler publishes a
wrong mask at exactly those widths on one executed target: they are 32, 64 and 128, the byte
powers of two above sixteen. The fourth is 14, excluded on a separate receipt, on a second target
and in a second optimization mode, recorded below.

`loongarch64-linux` in `Debug` is the measured target. The one-hot sweep marks one lane at a time so the marked lane holds `-0.0` against `+0.0` while every other lane holds `+1.0` against `+0.0`, and it requires the public `zmath.totalOrder` and `zmath.totalOrderMag` masks over `@Vector(n, T)` to answer the scalar kernel lane for lane and to hold exactly one true lane, on the marked index. At widths 32, 64 and 128 that target publishes a spurious true lane 0, and at no other swept width, so no marked lane other than lane 0 is one-hot there. Every other swept width is clean there.

The defect is upstream and predates this working tree: the identical failure appears when the same
program is built against a `git archive HEAD src` export, so nothing in this tree introduces it.
Reading a lane through `@shuffle` and `@reduce(.Or, ...)` instead of `mask[i]` gives the same wrong
verdicts, so the published mask value is wrong and not the consumer's read of it. Marking
`mapBinaryR` inline in a scratch copy changes nothing. Four candidate publication routes were
measured on that target: direct per-lane stores into `@Vector(n, bool)` and per-lane byte
comparisons are correct in isolation and fail once the same code sits inside a larger program, while
a `@bitCast` byte-packed integer publication and its `@shuffle`-narrowed form fail in both `Debug`
and `ReleaseSafe` with an exact full reversal. No candidate route is reliable at these widths on
that target, so the exclusion records a compiler limitation rather than a repair left out of
`src/vector_core.zig`.

`mask_widths`, `mask_boundary_widths` and the comptime width list in `src/vector_surface_test.zig`
omit 32, 64 and 128, `mask_widths` also omits 14, and every assertion that still runs is unchanged.
No publication route lost its last width, so no width is substituted: 8, 16 and 24 still fill a
whole lane group of staged bytes, 2, 4, 8 and 16 are still the power-of-two widths the `Debug` group
reversal acts inside, 1, 3, 5, 7, 9, 15 and 33 still publish by comparing staged lanes with 33
above one lane group and not byte aligned, and every swept width still reaches the direct store
route on each target that keeps direct vector storage.

Width 14 is excluded on a second, separate receipt, and it is a different failure from the three
above. The `cross-test` run for `aarch64_be-linux` in `ReleaseSafe` reports 955 of 956, and the
single failure is `vector_surface_test.test.the boolean family operations preserve one-hot masks
across result widths`, reporting "expected true, found false". A standalone probe built against
this library for the same target and the same mode, iterating the shipped `mask_widths` list,
reports `width=14 bad=6`, that is 6 of the 14 marked lanes wrong at width 14, and every other
width in that list is clean on the same target and in the same mode: 1, 2, 3, 4, 5, 7, 8, 9, 12,
13, 15, 16, 17, 20, 24, 28, 29, 30, 31, 33, 63, 65 and 127.

Width 14 is not a multiple of eight, so it reaches the staged-byte comparison route rather than
the integer route, which accepts only widths that are a multiple of eight, and the widths 8, 16
and 24 that do exercise the integer route are clean on that target and in that mode. The same
width 14 is clean on `aarch64_be-linux` in `Debug` and in isolated single-width programs, so this
is sensitive to the shape of the surrounding program in the same way the `loongarch64-linux`
receipts above are, and it is not a stable width-specific defect. The excluded set is therefore
not a stable boundary. It records what the shipped program shape measures now, and further widths
may have to be excluded as the program shape changes.

## `armeb-linux` Debug register-read defect

`armeb-linux` in `Debug` is expected to report 950 of 956, and the six failures are the
pre-existing `toVector` and `fromVector` lane-order tests in `src/vector_value_test.zig`. The same
six tests fail at the parent commit, which reports 948 of 954, so this tree neither introduces nor
repairs the defect and the change is not a regression.

The cause is LLVM's ARM big-endian backend, not the bridge. The IR is correct at both ends, and
the backend computes the big-endian element reversal of the returned `@Vector` into a scratch
register, `d16`, that is never read, then returns the unreversed register, `d0`, while the caller
applies one extra reversal. The net effect is one lane swap per forwarded result, so a test that
reads the returned register sees the wrong lanes.

Reading the published bytes instead of the returned register is correct, which is why
`toVector publishes lane zero at the lowest address, on every target` passes on `armeb-linux`
while the register-read tests fail there.

Three source forms of `toVector` were measured and all three leave those six tests failing at the
same count: the original per-lane `@Vector` assignment, a private `inline` variant, and an
array-plus-`@bitCast` variant. The array variant was reverted, so `src/vector_value.zig` keeps the
original `@Vector` form and is byte-equivalent to the parent commit there. No source change inside
`toVectorOf` avoids the defect, so it is recorded as an upstream compiler limitation rather than
left as an open repair in this tree.

## Verification status

- Native `Debug` and native `ReleaseSafe` both report 956 of 956.
- The executed full-matrix `cross-test` result is `Debug`: 19 pass, 1 fail, 33 build_only across the
  53 Zig 0.16.0 Tier 1/2 targets, run `20261004-153923-46001-17397` under the project profile with
  `--jobs 16 --build-jobs 3` in 489.79 s. The single failure is `armeb-linux` in `Debug` at 950 of
  956. A `build_only` target is a compile, not a test result, so the 33 are not passes.
- A `ReleaseSafe` `cross-test` ran over two targets with
  `--targets armeb-linux,aarch64_be-linux --jobs 2 --build-jobs 28 -- -Doptimize=ReleaseSafe`:
  2 pass, 0 fail, 0 build_only, and both target logs read all 956 tests passed. The byte-lane route
  therefore does have a cross-target optimized receipt.
- The authoritative receipt is the repository's own test binary, installed with
  `zig build test-build -Dtarget=<target> -Dcpu=baseline` and run under the `cross-test` qemu
  runner. A standalone probe is not a receipt for this tree.
- None of the above is whole-domain accuracy evidence, and no corpus or oracle result is claimed
  here; see [Numerical evidence](#numerical-evidence).
- The configured gate is not green. It also fails on frozen lint debt that predates this change:
  `ziglint` 265, `zlint` 269, `zbc` 145 and `zwanzig` 375 findings over the same 213 files. No
  green gate is claimed for this tree.
- Those figures are the campaign's frozen measured source snapshot, not the parent commit, so the
  two sets are not a parent-commit comparison.
- The gate measured on the current tree is fmt, build, test and ast-check OK, `ziglint` 263,
  `zlint` 263, `zbc` 146 and `zwanzig` 353 in 8 rule classes, `cross-test` 19 pass, 1 fail and 33
  build_only on the known `armeb-linux` `Debug` case, and exit 1.
- Per file, `zbc` moved `src/vector_core.zig` 13 to 6, `src/geometry.zig` 28 to 22 and
  `src/vector_test.zig` 2 to 0 against `src/vector_value_test.zig` 2 to 8, `src/linalg.zig` 19 to
  23, `src/vector_value.zig` 24 to 28 and `src/geometry_test.zig` 0 to 2, so the byte-staging
  change removed seven findings and the net rise of one belongs to the wider uncommitted change set
  rather than to the lane-publication files, and `src/vector_value.zig` is byte-equivalent to the
  parent commit.

## Not covered by evidence

The non-scalar modules have no general accuracy proof against an independent oracle.
The corpus has ten `definitional` rows for four Cl(3) involutions and quaternion multiplication in `f32` and `f64`.
These rows check 176 stored convention cases and do not establish either complete module's numerical accuracy.
Ten further numerical rows are five selected operations, each with one `f32` and one `f64` row: `complex` addition and multiplication, `matrix.determinantMat2`, `linalg.lu` with `solveLu` on a 2x2 system, and the 2D orientation predicate under both policies.
Those ten rows check 324 stored bounded cases against an independent exact-rational reference.
Each row follows the rounding of the graph it checks, so the rule is per primitive step for `complexMul`, `matrixDet2`, and the `fast` predicate, and one round per component for the exact rational `luSolve2` solution.
This is the selected sample the user accepted: the other operations of complex, matrix, linalg, and geometry stay outside the manifest as a scope limit.
Their unit tests do not substitute for oracle evidence.

## Numerical evidence

The [scalar contract](spec/core-math-2026-08-22.md) defines the corpus schema and target profile. The [checked-in manifest](../zmath-testing/test/corpus/manifest.toml) in sibling `zmath-testing` is the acceptance authority. Its oracle generation, corpus replay, and benchmarks are separate from this runtime. A source export or selected corpus does not prove whole-domain accuracy.

## Commands

From this repository, use `zig build` to install the demonstration and `zig build run` to run it. Use `zig build test` for local tests. This build graph has no `bench` step. Run qualified oracle and benchmark steps from `zmath-testing` as documented in [API.md: Tests and benchmarks](API.md#tests-and-benchmarks).

Use `cross-test` to attempt the complete 53-platform Zig 0.16.0 Tier 1/2
matrix inherited by `cross-test.profile.sh`; append `-- -Doptimize=ReleaseSafe`
for optimized correctness checks. `cross-test --list` shows detected runners.
Missing foreign-OS execution environments remain build-only, and
`--require-runners` makes incomplete execution coverage fail.

The only executed `cross-test` result in this tree is `Debug`: 19 pass, 1 fail and 33 build_only
across those 53 targets. Native `ReleaseSafe` is 956 of 956. See
[Verification status](#verification-status).

`zig build test-build` installs the same test executable without executing it.
These checks run no corpus or profiling workload and do not replace the gate.
