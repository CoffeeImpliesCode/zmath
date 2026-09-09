# Historical zmath Correctness Roadmap (superseded)

> **HISTORICAL / SUPERSEDED — not a live plan.** This August 2026 proposal records old names, paths, and commands. Its unchecked tasks remain unchecked; do not execute them or infer completion. The old uppercase function suffixes, `q` plan, and local `zig build bench`/`oracle` commands do not describe the current API or build graph.
>
> Use the current [API contract](../../../API.md), [schema-3 scalar contract](../../../spec/core-math-2026-08-22.md), and sibling [corpus manifest](../../../../zmath-testing/test/corpus/manifest.toml) for names, commands, and acceptance status. `zmath-testing`, not `zmath`, owns oracle and benchmark steps.

**Historical goal:** Build a dependency-free, fast, and correct f16/f32/f64 scalar foundation with exhaustive corner-case coverage and independent MPFR/Arb validation in `zmath-testing`.

**Proposed architecture:** The runtime library uses format-specific IEEE bit primitives and scalar kernels. `zmath-testing` owns external oracles, immutable corpora, fuzz replay, and benchmarks. Every exported numerical symbol gets a pinned manifest entry, red test, portable contract test, and oracle comparison before release.

**Historical tech stack:** Zig 0.16.0, Jujutsu, devenv/Nix, MPFR, GMP, FLINT/Arb.

## Proposed Global Constraints (historical)

- Runtime zmath links no MPFR, GMP, FLINT, Arb, BLIS, Eigen, C++, or libc math library.
- Public type suffixes are `h` for `f16`, `f` for `f32`, `d` for `f64`, and `q` for `f128`.
- Public function suffixes are `H`, `F`, `D`, and `Q`.
- `q` is reserved until f16/f32/f64 oracle coverage and range-reduction infrastructure are complete.
- bfloat16 has no public API or test fixture.
- No unsuffixed numeric aliases exist.
- Every production behavior starts with a failing test.
- Tests run with explicit output-bit, class, rounding, and NaN rules.
- Commit only green Jujutsu changes.

---

## Proposed tasks (not for execution)

### Task 1: Migrate public numeric naming

**Files:**
- Modify: `API.md`
- Modify: `src/ieee.zig`
- Modify: `src/scalar.zig`
- Modify: `src/root.zig`
- Create: `src/api_suffix_test.zig`

**Interfaces:**

```zig
pub fn classifyH(value: f16) Class;
pub fn classifyF(value: f32) Class;
pub fn classifyD(value: f64) Class;

pub fn absH(value: f16) f16;
pub fn absF(value: f32) f32;
pub fn absD(value: f64) f64;
```

- [ ] **Step 1: Write the failing public suffix test**

Create `src/api_suffix_test.zig` with compile-time references to every current H/F/D IEEE and scalar export. Do not reference Q before f128 exists.

```zig
test "numeric exports use H F and D suffixes" {
    _ = ieee.classifyH;
    _ = ieee.classifyF;
    _ = ieee.classifyD;
    _ = scalar.absH;
    _ = scalar.absF;
    _ = scalar.absD;
}
```

- [ ] **Step 2: Run the suffix test red**

Run:

```bash
zig test src/api_suffix_test.zig
```

Expected: failure because the f16 functions are absent.

- [ ] **Step 3: Rename and add format-specific exports**

Keep existing F/D behavior. Add H implementations only after their red tests exist. Remove any `F32` or `F64` public name. Do not add compatibility aliases.

- [ ] **Step 4: Update API.md**

Define lowercase type suffixes and uppercase function suffixes. State that Q names are reserved, not exported.

- [ ] **Step 5: Verify green**

Run:

```bash
zig test src/api_suffix_test.zig
zig build test
```

- [ ] **Step 6: Commit**

```bash
jj commit -m "refactor: standardize numeric suffixes"
```

### Task 2: Complete IEEE bit primitives for H/F/D

**Files:**
- Modify: `src/ieee.zig`
- Modify: `src/ieee_test.zig`
- Modify: `src/ieee_edge_test.zig`
- Create: `src/ieee_exhaustive_h_test.zig`
- Create: `src/ieee_transition_test.zig`

**Interfaces:**

```zig
pub fn toBitsH(value: f16) u16;
pub fn toBitsF(value: f32) u32;
pub fn toBitsD(value: f64) u64;
pub fn fromBitsH(bits: u16) f16;
pub fn fromBitsF(bits: u32) f32;
pub fn fromBitsD(bits: u64) f64;
pub fn nextUpH(value: f16) f16;
pub fn nextDownH(value: f16) f16;
pub fn scaleBH(mode: RoundMode, value: f16, exponent: i32) f16;
pub fn roundIntegralH(mode: RoundMode, value: f16) f16;
```

- [ ] **Step 1: Write H format red tests**

Test every H class and edge encoding: signed zero, minimum and maximum subnormal, minimum normal, one, each power-of-two transition, maximum finite, infinities, quiet NaN, and signaling NaN.

- [ ] **Step 2: Run H red tests**

Run:

```bash
zig test src/ieee_test.zig
```

Expected: missing H symbols.

- [ ] **Step 3: Implement H bit primitives**

Use `u16` sign, exponent, fraction, and quiet masks. Do not convert through f32.

- [ ] **Step 4: Add exhaustive H encoding test**

Loop over all `0x0000` through `0xffff` encodings. For each encoding, assert class, bit round trip, copy-sign behavior, successor/predecessor policy, and total-order consistency.

- [ ] **Step 5: Add F/D transition windows**

Test all encodings in fixed windows around zero, subnormal-normal boundaries, every exponent transition, maximum finite, infinities, and NaN payload boundaries. Shard F/D windows deterministically by sign and exponent.

- [ ] **Step 6: Verify green**

Run:

```bash
zig test src/ieee_test.zig
zig test src/ieee_edge_test.zig
zig test src/ieee_exhaustive_h_test.zig
zig test src/ieee_transition_test.zig
```

- [ ] **Step 7: Commit**

```bash
jj commit -m "feat: complete IEEE H F D primitives"
```

### Task 3: Freeze scalar error measurement

**Files:**
- Modify: `src/scalar.zig`
- Modify: `src/scalar_test.zig`
- Modify: `src/scalar_boundary_test.zig`
- Create: `src/scalar_property_test.zig`

**Interfaces:**

```zig
pub fn ulpDistanceH(a: f16, b: f16) u16;
pub fn ulpDistanceF(a: f32, b: f32) u32;
pub fn ulpDistanceD(a: f64, b: f64) u64;
```

- [ ] **Step 1: Write red H ULP tests**

Test equal values, opposite signed zero, adjacent values, minimum negative to minimum positive subnormal, infinity span, and NaN sentinel behavior.

- [ ] **Step 2: Run H ULP test red**

Run:

```bash
zig test src/scalar_test.zig
```

Expected: missing H symbols.

- [ ] **Step 3: Implement H absolute value and ULP distance**

Use H ordered-bit ranking. Define the same signed-zero and NaN sentinel policy as F/D.

- [ ] **Step 4: Add properties**

Assert symmetry, identity, adjacency, monotonic rank behavior, and no panic for all H bit patterns plus F/D transition corpora.

- [ ] **Step 5: Verify green**

Run:

```bash
zig test src/scalar_test.zig
zig test src/scalar_boundary_test.zig
zig test src/scalar_property_test.zig
```

- [ ] **Step 6: Commit**

```bash
jj commit -m "feat: add H scalar error primitives"
```

### Task 4: Build the pinned scalar manifest

**Files:**
- Create: `spec/core-math-<commit>.md`
- Create: `zmath-testing/test/corpus/manifest.toml`
- Create: `zmath-testing/test/corpus/edge_cases.tsv`
- Create: `zmath-testing/test/corpus/cases.tsv`
- Create: `zmath-testing/test/corpus/README.md`

**Interfaces:**

Each TSV row has:

```text
op	format	input_bits	output_bits	output_class	max_ulp
```

Format names are `h`, `f`, and `d`. Bit strings are lowercase hexadecimal. Q rows are rejected until Q exists.

- [ ] **Step 1: Add the manifest parser red test**

Create `zmath-testing/test/oracle/manifest_test.zig`. Require a pinned CORE-MATH commit, schema version, MPFR precision, Arb precision, and one function-format row per exported scalar function.

- [ ] **Step 2: Run manifest parser red**

Run:

```bash
zig build oracle
```

Expected: missing parser and corpus loader.

- [ ] **Step 3: Add immutable manifest and edge corpus**

Pin the exact upstream commit. Include H/F/D rows for `exp`, `log`, and `sin`. Mark a row `unsupported` until its implementation and oracle pass.

- [ ] **Step 4: Implement parser and validation**

Reject unknown formats, uppercase hexadecimal, duplicate rows, missing special-value rules, unpinned commits, and Q rows before Q export.

- [ ] **Step 5: Verify green**

Run:

```bash
zig build oracle
```

- [ ] **Step 6: Commit**

```bash
jj commit -m "test: add pinned scalar corpus manifest"
```

### Task 5: Build MPFR and Arb differential infrastructure

**Files:**
- Modify: `zmath-testing/build.zig`
- Create: `zmath-testing/test/oracle/oracle_bridge.h`
- Create: `zmath-testing/test/oracle/oracle_bridge.c`
- Create: `zmath-testing/test/oracle/oracle_corpus.zig`
- Create: `zmath-testing/test/support/ulp.zig`

**Interfaces:**

```c
int zmath_testing_mpfr_reference(
    const char *op,
    unsigned format,
    const uint64_t *input_bits,
    size_t arity,
    uint64_t *output_bits,
    unsigned *output_class
);
```

`format` is `1` for H, `2` for F, and `3` for D. The bridge rejects every other value.

- [ ] **Step 1: Write bridge ABI red tests**

Test malformed operation names, invalid formats, wrong arity, exact H bit decoding, MPFR result encoding, and Arb uncertainty failure.

- [ ] **Step 2: Run bridge tests red**

Run:

```bash
zig build oracle
```

Expected: missing bridge symbols.

- [ ] **Step 3: Implement MPFR target-format rounding**

Decode input bits exactly. Evaluate with pinned MPFR precision and requested rounding. Return canonical target-format bits and class.

- [ ] **Step 4: Implement Arb enclosure validation**

Fail closed when the Arb enclosure cannot prove the rounded target result. Do not treat raw point containment as proof.

- [ ] **Step 5: Add test tiers**

Add `-Dtest-tier=portable|oracle|all`. Default `zig build test` remains portable and does not link oracle libraries. `zig build oracle` selects oracle tests. `zig build all-tests` runs both groups.

- [ ] **Step 6: Verify green**

Run:

```bash
zig build test
zig build oracle
zig build all-tests
```

- [ ] **Step 7: Commit**

```bash
jj commit -m "test: add MPFR Arb scalar oracle bridge"
```

### Task 6: Implement scalar functions one family at a time

**Files:**
- Create: `src/scalar/exp.zig`
- Create: `src/scalar/log.zig`
- Create: `src/scalar/trig.zig`
- Modify: `src/scalar.zig`
- Modify: `src/root.zig`
- Create: `src/scalar_transcendental_test.zig`
- Create: `src/scalar_fuzz_test.zig`

**Interfaces:**

```zig
pub fn expH(mode: ieee.RoundMode, value: f16) f16;
pub fn expF(mode: ieee.RoundMode, value: f32) f32;
pub fn expD(mode: ieee.RoundMode, value: f64) f64;

pub fn logH(mode: ieee.RoundMode, value: f16) f16;
pub fn logF(mode: ieee.RoundMode, value: f32) f32;
pub fn logD(mode: ieee.RoundMode, value: f64) f64;

pub fn sinH(mode: ieee.RoundMode, value: f16) f16;
pub fn sinF(mode: ieee.RoundMode, value: f32) f32;
pub fn sinD(mode: ieee.RoundMode, value: f64) f64;
```

- [ ] **Step 1: Add Exp red tests**

For H, test every encoding against the pinned corpus. For F/D, test special-value tables, hard cases, corpus vectors, and manifest-supported modes.

- [ ] **Step 2: Run Exp tests red**

Run:

```bash
zig test src/scalar_transcendental_test.zig
zig build oracle -- --op exp
```

Expected: missing `expH/F/D`.

- [ ] **Step 3: Implement Exp minimal green path**

Port only the pinned dependency-free algorithm, constants, range reduction, and fallback required by the manifest row. Do not call host libm.

- [ ] **Step 4: Verify Exp**

Run portable, H exhaustive, MPFR/Arb, hard-case, property, and fuzz-replay tests.

- [ ] **Step 5: Commit Exp**

```bash
jj commit -m "feat: add manifest-backed exp kernels"
```

- [ ] **Step 6: Repeat the same red-green sequence for Log**

Use the `log` manifest rows, domain and pole table, hard cases, and corpus. Commit `feat: add manifest-backed log kernels`.

- [ ] **Step 7: Repeat the same red-green sequence for Sin**

Use the `sin` manifest rows, oddness/signed-zero table, range-reduction hard cases, and corpus. Commit `feat: add manifest-backed sin kernels`.

### Task 7: Add reproducible benchmarks

**Files:**
- Modify: `zmath-testing/build.zig`
- Create: `zmath-testing/bench/targets.zig`
- Create: `zmath-testing/bench/bench_cases.zig`
- Create: `zmath-testing/bench/target_info.zig`
- Modify: `zmath-testing/bench/public_api.zig`

**Interfaces:**

```text
zig build bench -Dbench-target=ieee-nextafter-f
zig build run-bench -Dbench-target=scalar-exp-h
```

- [ ] **Step 1: Add benchmark target parser red test**

Reject unknown targets. Require deterministic input corpus, operation count, and checksum.

- [ ] **Step 2: Run target parser red**

Run:

```bash
zig build bench
```

Expected: missing target table.

- [ ] **Step 3: Implement target table and metadata**

Print target name, corpus ID, iteration count, checksum, optimize mode, target triple, and sorted resolved CPU features. Do not probe the host CPU at runtime.

- [ ] **Step 4: Add benchmark targets**

Add IEEE next-after H/F/D and every released scalar function-format pair. Benchmarks do not link MPFR or Arb.

- [ ] **Step 5: Verify green**

Run:

```bash
zig build bench -Dbench-target=ieee-nextafter-f
zig build run-bench -Dbench-target=ieee-nextafter-f
```

- [ ] **Step 6: Commit**

```bash
jj commit -m "bench: add deterministic scalar targets"
```

## Proposed Coverage Matrix (historical)

| Tier | Scope | Pass condition |
| --- | --- | --- |
| Fast | H/F/D class and edge tables | Exact expected bits and classes |
| Exhaustive H | Every 16-bit encoding for each unary IEEE and released scalar H function | Every input passes exact manifest policy |
| F transition | Boundary windows and deterministic shards | Every configured shard passes |
| D transition | Boundary windows, exponents, payload extrema, and hard corpus | Every configured shard passes |
| Oracle | MPFR exact target rounding plus Arb enclosure | Exact bits or manifest ULP bound |
| Property | Algebraic and ordering invariants | No invariant violation |
| Fuzz | Seeded category-weighted bit/pair generators | No panic, hang, or invariant violation |
| Benchmark | Fixed corpus and deterministic metadata | Build succeeds and checksum is fixed |

## Plan Review (historical)

- Proposed implementation order at the time: IEEE H/F/D completion, manifest, oracle bridge, then Exp, Log, and Sin H/F/D.
- Corner cases: every IEEE class, payload boundary, transition, domain edge, rounding mode, hard case, property, and fuzz-replay path has a named tier.
- Infrastructure: the oracle bridge, corpus parser, tier selection, and deterministic benchmark metadata remain in `zmath-testing`.
- Dependency-free runtime: every external numerical dependency stays outside zmath.
- Proposed naming at the time: H/F/D/Q function suffixes and h/f/d/q type suffixes before scalar transcendental exports. This is not the current API.
