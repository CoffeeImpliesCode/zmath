# Scalar Acceptance Roadmap Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `subagent-driven-development` or `executing-plans` to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Convert zmath's broad C23 scalar implementation from reflection-tested API surface into a reproducible, manifest-backed, target-qualified numerical release foundation.

**Architecture:** `zmath` remains a dependency-free runtime. The sibling `zmath-testing` project owns manifest parsing, immutable corpora, MPFR/Arb production and proof, and benchmark artifacts. The implementation sequence starts by reconciling the contract, then creates a strict schema-v2 corpus pipeline, then accepts one function-format pair at a time. No SIMD, complex, vector, matrix, geometry, or algebra work begins until the scalar acceptance pipeline is reproducible.

**Tech Stack:** Zig 0.16.0, Jujutsu, Nix/devenv, MPFR 512-bit, Arb 512-bit via FLINT, vendored CORE-MATH C kernels.

## Global Constraints

- Runtime `zmath` must not link MPFR, GMP, FLINT, Arb, BLIS, Eigen, C++, libc math, `fenv`, or `errno` behavior.
- Public floating formats are only `f16`, `f32`, and `f64`; `q`, `f128`, and bfloat16 remain rejected.
- Public scalar suffixes are `h`, `f`, and `d`; generic root functions reject mismatched or unsupported floating types at compile time.
- Every acceptance transition requires exact special-value behavior, caller round-to-nearest/ties-to-even with gradual underflow, MPFR `MPFR_RNDN` at exactly 512 bits, and a unique matching Arb rounding proof at exactly 512 bits.
- Every corpus and benchmark input set is repository-relative, immutable, deterministic, and has an explicit selection rule and count.
- `unsupported` means no acceptance claim. Export presence, local unit tests, source provenance, or a passing build never changes row status.
- The accepted manifest has 315 sorted schema-v2 rows. Exactly six pairs are accepted: `exp/h`, `exp/f`, `exp/d`, `log/h`, `log/f`, and `sin/h`. The remaining 309 rows are `unsupported`.
- The accepted corpus profile is exactly target `x86_64-linux-gnu.2.38`, CPU `haswell`, and mode `ReleaseSafe`. Do not claim support for another profile.
- Implementation, manifest acceptance, corpus evidence, and whole-domain proof are separate facts. F and D selected corpora are not exhaustive whole-domain proof.
- Schema-v2 `bounded` rows permit any non-negative `max_ulp`, including zero. `correct` rows require `max_ulp = 0`.

## Current integration state

The source, strict schema-v2 manifest, corpus loader, deterministic replay lane,
and benchmark metadata for the six accepted pairs are integrated.
This roadmap now records re-verification of that work.
It does not invent historical red runs or retroactive failure evidence.
Verification and exact-power resource tasks remain open until measured evidence is
recorded.

The broad C23 campaign and evaluator branch remains preserved at the
`preserve-c23-campaigns` bookmark.
Only strict manifest and schema rejection regressions move into this roadmap.
The obsolete approximate `atan2` split and its invalid audit constants are superseded by the canonical implementation and migrated regression cases.
The current `lgamma` contract governs zero and finite negative-integer poles.

---

### Task 1: Reconcile the scalar acceptance contract

**Files:**
- Modify: `spec/core-math-2026-08-22.md:486-622`
- Modify: `API.md:641-658`
- Modify: `~/projects/zmath-testing/test/corpus/manifest.toml`
- Test: `~/projects/zmath-testing/test/integration/infrastructure.zig`

**Decision:** Keep schema version 2 and the canonical 315 sorted typed rows. The manifest accepts exactly six pairs: `exp/h`, `exp/f`, `exp/d`, `log/h`, `log/f`, and `sin/h`. The remaining 309 rows stay `unsupported` until their own implementation and evidence exist.

**Interfaces:**
- Root manifest keys are exactly: `schema_version`, `core_math_repository`, `core_math_revision`, `mpfr_precision_bits`, `arb_precision_bits`, `rounding`, `formats`, `reserved_formats`, `rejected_formats`, `operation_spelling`, `columns`, `hexadecimal`, `uncertainty_policy`, and `special_policy`.
- Every function row contains exactly: `op`, `provenance`, `status`, `arity`, `inputs`, `result`, `oracle`, `special_policy`, and `max_ulp`.
- Define the allowed `operation_spelling` value and the declared provenance vocabulary in the spec before parser implementation. The parser must reject values outside those sets.

- [x] **Step 1: Write contract regression tests**
  - Integrated strict rejection tests for the old schema-v1 layout, missing `operation_spelling`, unknown provenance, and invalid `correct` row relations.
- [ ] **Step 2: Reverify the infrastructure test**
  - Run `zig build test` in `~/projects/zmath-testing`.
  - This roadmap records no historical red result. Reverification remains open.
- [x] **Step 3: Make the specification unambiguous**
  - Reconcile the current six accepted rows, the 309 unsupported rows, exact `operation_spelling`, provenance, target profile, and evidence boundaries.
- [x] **Step 4: Use the checked-in sorted schema-v2 manifest**
  - Keep all 315 descriptors in deterministic order.
  - Keep the six accepted rows at `status = "correct"`, `oracle = "mpfr_arb"`, and `max_ulp = 0`.
  - Keep the remaining 309 rows at `status = "unsupported"`, `oracle = "none"`, and `max_ulp = 0`.
- [ ] **Step 5: Verify the contract gate**
  - Run `zig build test` in `~/projects/zmath-testing`.
### Task 2: Re-verify strict schema-v2 parsing and artifact loading


**Files:**
- Modify: `~/projects/zmath-testing/test/support/corpus.zig`
- Modify: `~/projects/zmath-testing/test/integration/infrastructure.zig`
- Modify: `~/projects/zmath-testing/build.zig`
- Create: `~/projects/zmath-testing/test/corpus/exp_h.tsv`
- Create: `~/projects/zmath-testing/test/corpus/exp_f.tsv`

**Interfaces:**
```zig
pub const FunctionStatus = enum { unsupported, correct, bounded };
pub const CorpusRow = struct {
    op: Operation,
    format: Format,
    input_bits: u64,
    output_bits: u64,
    output_class: OutputClass,
    max_ulp: u64,
};
pub const TargetMode = struct {
    target: []const u8,
    cpu: []const u8,
    optimize: []const u8,
};
pub fn parseEmbeddedManifest() ManifestError!Manifest;
pub fn loadEmbeddedCorpus(comptime path: []const u8) CorpusError![]const CorpusRow;
```

- [x] **Step 1: Write strict parser tests**
  - Test missing, unknown, duplicate, out-of-order, and duplicate function keys.
  - Test every descriptor form used by the 315 rows and reject unsupported formats, malformed widths, non-lowercase hexadecimal, duplicate corpus tuples, and bounded rows without a declared non-negative bound.
- [ ] **Step 2: Reverify the infrastructure test**
  - Run `zig build test` in `~/projects/zmath-testing`.
  - This roadmap records no historical red result. Reverification remains open.
- [x] **Step 3: Replace fixed schema-v1 capacities and operation tables**
  - Parse the manifest into bounded fixed arrays sized for the checked-in 315 rows.
  - Keep parser work allocation-free or use explicit test-only allocation with bounded input size.
- [x] **Step 4: Embed immutable artifacts through the build graph**
  - Extend the `addWriteFiles` corpus package to copy every referenced TSV artifact.
  - Make missing artifacts a build failure.
- [ ] **Step 5: Verify parser and artifact integrity**
  - Run `zig build test` in `~/projects/zmath-testing`.
  - Confirm that the manifest and every embedded corpus parse, have deterministic counts, and bind to exactly one manifest row.

### Task 3: Re-verify differential execution

**Files:**
- Modify: `~/projects/zmath-testing/test/oracle/corpus_runner.zig`
- Modify: `~/projects/zmath-testing/test/oracle/oracle_bridge_test.zig`
- Modify: `~/projects/zmath-testing/test/oracle/oracle_bridge.c`
- Modify: `~/projects/zmath-testing/build.zig`
- Modify: `~/projects/zmath-testing/test/corpus/README.md`

**Interfaces:**
```zig
pub fn runCorpusRow(row: corpus.CorpusRow) DifferentialError!void;
pub fn runAcceptedCorpus(
    manifest: corpus.Manifest,
    target: TargetMode,
) DifferentialError!void;
```

- [x] **Step 1: Write deterministic-runner tests**
  - Integrated checks require runtime bits, stored corpus bits, MPFR candidate bits, and Arb proof bits to agree for every finite row.
  - Integrated checks require special-value rows to use stored exact target bits and reject an Arb ambiguity, MPFR disagreement, wrong target/mode, or wrong ULP limit.
- [ ] **Step 2: Reverify the oracle lane**
  - Run `zig build oracle -Dtarget=x86_64-linux-gnu.2.38 -Dcpu=haswell -Doptimize=ReleaseSafe` in `~/projects/zmath-testing`.
  - This roadmap records no historical red result. Reverification remains open.
- [x] **Step 3: Remove random acceptance sampling**
  - Replace random acceptance inputs with checked-in corpus rows.
  - Keep fuzzing as a separate explicitly seeded exploratory lane. It must not alter acceptance status.
- [x] **Step 4: Extend bridge tests before bridge families**
  - Cover H/F/D decoding, target-format conversion, NaN and signed-zero handling, MPFR/Arb disagreement, and unique-rounding failure for unary Exp/Log/Sin.
  - Do not add C23 binary, ternary, or integer-parameter bridge operations in this task.
- [x] **Step 5: Record corpus provenance**
  - Document the qualified producer and replay commands, pinned inputs, output class policy, tuple selection, counts, and target profile in `test/corpus/README.md`.
- [ ] **Step 6: Reverify reproducibility**
  - Run the qualified oracle command twice with target `x86_64-linux-gnu.2.38`, CPU `haswell`, and mode `ReleaseSafe`.
  - Confirm identical row count, selected corpus, and pass/fail result.

### Task 4: Re-verify the six accepted unary rows

**Files:**
- Modify: `~/projects/zmath-testing/test/corpus/manifest.toml`
- Modify: `~/projects/zmath-testing/test/oracle/corpus_runner.zig`
- Create: `~/projects/zmath-testing/test/corpus/log_h.tsv`
- Create: `~/projects/zmath-testing/test/corpus/sin_h.tsv`
- Modify: `~/projects/zmath-testing/test/corpus/README.md`

**Interfaces:**
- Each accepted row is `exp/h`, `exp/f`, `exp/d`, `log/h`, `log/f`, or `sin/h`.
- Each accepted row uses `status = "correct"`, `oracle = "mpfr_arb"`, `special_policy = "zmath_c23_v1"`, and `max_ulp = 0`.

- [x] **Step 1: Generate immutable corpus artifacts**
  - Enumerate all 65,536 H inputs for `exp/h`, `log/h`, and `sin/h`.
  - Add deterministic F/D `exp` hard-case rows and the selected `log/f` rows with exact bits.
  - Treat F/D artifacts as evidence for stored inputs, not exhaustive whole-domain proof.
- [ ] **Step 2: Reverify each accepted corpus**
  - Run the qualified MPFR and Arb replay for each accepted pair.
  - Keep any disputed row unsupported until its stored row, MPFR result, and Arb enclosure agree.
  - This roadmap records no historical red result.
- [x] **Step 3: Fix runtime only when an oracle row proves a defect**
  - Touch only the exact kernel or wrapper owning a checked-in mismatch. Preserve pure-value behavior and caller-RNE requirements.
- [x] **Step 4: Keep one manifest status per accepted proof**
  - The current manifest accepts all six pairs independently. Do not change adjacent statuses when a future pair changes.
- [ ] **Step 5: Verify the accepted corpus matrix**
  - Run `zig build test` in `zmath`.
  - Run `zig build oracle -Dtarget=x86_64-linux-gnu.2.38 -Dcpu=haswell -Doptimize=ReleaseSafe` in `zmath-testing`.
  - Confirm exact zero-ULP evidence for all six accepted rows and `unsupported` status for the remaining 309 rows.

### Task 5: Re-verify the accepted `log/f` pair

**Files:**
- Modify: `~/projects/zmath-testing/test/corpus/manifest.toml`
- Create: `~/projects/zmath-testing/test/corpus/log_f.tsv`
- Modify: `~/projects/zmath-testing/test/oracle/corpus_runner.zig`
- Modify: `~/projects/zmath-testing/test/oracle/oracle_bridge_test.zig`
- Modify: `src/scalar_log_test.zig`

**Interfaces:**
```zig
// Existing runtime path; do not add an alias.
pub fn logf(value: f32) f32;
```

- [x] **Step 1: Add the immutable `log/f` corpus**
  - Include signed zero, subnormal-normal boundaries, `-1`, the domain interval below zero, infinities, signaling and quiet NaNs, adjacent values around one, rounding midpoints, and existing CORE-MATH fallback vectors.
- [ ] **Step 2: Reverify the `log/f` oracle lane**
  - Run `zig build oracle -Dtarget=x86_64-linux-gnu.2.38 -Dcpu=haswell -Doptimize=ReleaseSafe` in `~/projects/zmath-testing`.
  - This roadmap records no historical red result. Reverification remains open.
- [x] **Step 3: Add `log/f` corpus dispatch**
  - Route only `Operation.log` plus `Format.f` to `zmath.scalar.logf`.
  - Compare runtime, stored, MPFR, and Arb results. Reject any mismatch.
- [x] **Step 4: Repair only oracle-proven runtime defects**
  - Modify `src/scalar.zig` only if a checked-in failing corpus row demonstrates a mismatch.
  - Add the exact regression input to `src/scalar_log_test.zig` before changing the kernel.
- [x] **Step 5: Keep exactly one accepted row**
  - The current manifest sets only `log/f` in this pair's transition to `correct`, `mpfr_arb`, and zero ULP.
- [ ] **Step 6: Verify no adjacent status changed**
  - Run `zig build test` in `zmath`.
  - Run `zig build oracle -Dtarget=x86_64-linux-gnu.2.38 -Dcpu=haswell -Doptimize=ReleaseSafe` in `zmath-testing`.
  - Confirm `log/f` is accepted. Confirm `log/d`, `sin/f`, and `sin/d` stay unsupported.

### Task 6: Re-verify target qualification and scalar benchmarks

**Files:**
- Modify: `build.zig`
- Modify: `~/projects/zmath-testing/build.zig`
- Create: `~/projects/zmath-testing/bench/targets.zig`
- Create: `~/projects/zmath-testing/bench/target_info.zig`
- Modify: `~/projects/zmath-testing/bench/public_api.zig`

**Interfaces:**
```text
zig build bench -Dbench-target=scalar-log-f -Dtarget=x86_64-linux-gnu.2.38 -Dcpu=haswell -Doptimize=ReleaseSafe
zig build run-bench -Dbench-target=scalar-log-f -Dtarget=x86_64-linux-gnu.2.38 -Dcpu=haswell -Doptimize=ReleaseSafe
```

- [x] **Step 1: Write target-policy tests**
  - Reject targets outside the explicit C-kernel ABI allowlist.
  - Reject a benchmark target without a corpus ID, checked-in corpus path, operation count, checksum, target triple, sorted feature set, and optimize mode.
- [ ] **Step 2: Reverify benchmark target parsing**
  - Run `zig build bench -Dbench-target=scalar-log-f -Dtarget=x86_64-linux-gnu.2.38 -Dcpu=haswell -Doptimize=ReleaseSafe` in `~/projects/zmath-testing`.
  - This roadmap records no historical red result. Reverification remains open.
- [x] **Step 3: Make target restrictions explicit**
  - Enforce exactly target `x86_64-linux-gnu.2.38`, CPU `haswell`, and mode `ReleaseSafe` for the accepted corpus profile.
  - Do not claim support for another target, CPU, or optimize mode.
- [x] **Step 4: Add deterministic `scalar-log-f` benchmarking**
  - Use the immutable `log_f.tsv` corpus, fixed iteration count, checksum, and no MPFR/Arb linkage.
  - Print corpus ID, input count, target triple, resolved CPU features, optimize mode, and checksum.
- [ ] **Step 5: Measure exact-power resource bounds separately**
  - Add ReleaseSafe and ReleaseFast tests around minimum-subnormal inputs and degree ±64 for `pown`, `compoundn`, and `rootn`.
  - Record stack and latency evidence before accepting any power family.
  - Do not claim a resource bound or benchmark result before measurement.
- [ ] **Step 6: Verify target-qualified release lanes**
  - Run local tests, the qualified oracle command, and the qualified scalar benchmark.
  - Confirm runtime remains dependency-free. Confirm external numerical libraries link only into oracle artifacts.

## Deferred Scope

- SIMD dispatch, AVX2/AVX-512/NEON kernels, and native CPU benchmarking remain deferred.
- Complex numbers, vectors, quaternions, matrices, linear algebra, geometry, and algebra modules remain deferred.
- Oracle bridge support for binary, ternary, integer-parameter, and C-backed C23 families remains deferred.
- Whole-domain proof for F and D scalar families remains deferred. Their selected corpora are not exhaustive.

The broad C23 campaign and evaluator branch stays at the
`preserve-c23-campaigns` bookmark.
Only strict manifest and schema rejection regressions migrated into this lane.
Do not migrate descriptor runtime or broad campaign changes.
The deferred work starts only after the accepted scalar lanes have re-verification
evidence and the exact-power resource envelope has measured evidence.

## Plan Self-Review

- **Spec coverage:** schema v2, all 315 row descriptors, six accepted pairs, 309 unsupported rows, corpus artifacts, MPFR/Arb proof, target/mode metadata, benchmark metadata, and exact-power resource limits map to Tasks 1-6.
- **Evidence boundaries:** implementation, manifest acceptance, stored corpus evidence, and whole-domain proof remain separate. F and D samples do not establish exhaustive proof.
- **Dependency order:** no future status transition precedes contract reconciliation, parser migration, immutable corpus loading, deterministic oracle execution, or target qualification.
- **Scope:** the plan preserves the dependency-free runtime boundary and excludes premature SIMD, broad C23 campaign migration, and planned non-scalar modules.

