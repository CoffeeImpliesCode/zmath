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
- The current working copy is user-owned and already contains an uncommitted C23 integration. Execute each task in an isolated Jujutsu change after that integration is checkpointed.

---

### Task 1: Reconcile the scalar acceptance contract

**Files:**
- Modify: `spec/core-math-2026-08-22.md:486-622`
- Modify: `API.md:641-658`
- Modify: `~/projects/zmath-testing/test/corpus/manifest.toml`
- Test: `~/projects/zmath-testing/test/integration/infrastructure.zig`

**Decision:** Use the conservative baseline: all 315 schema-v2 rows begin `unsupported` until their evidence is generated and reviewed. Do not preserve the existing five `correct` rows without the required immutable corpus, allowlisted target/mode matrix, and benchmark record.

**Interfaces:**
- Root manifest keys are exactly: `schema_version`, `core_math_repository`, `core_math_revision`, `mpfr_precision_bits`, `arb_precision_bits`, `rounding`, `formats`, `reserved_formats`, `rejected_formats`, `operation_spelling`, `columns`, `hexadecimal`, `uncertainty_policy`, and `special_policy`.
- Every function row contains exactly: `op`, `provenance`, `status`, `arity`, `inputs`, `result`, `oracle`, `special_policy`, and `max_ulp`.
- Define the allowed `operation_spelling` value and the declared provenance vocabulary in the spec before parser implementation. The parser must reject values outside those sets.

- [ ] **Step 1: Write contract regression tests**
  - Add tests that reject the old schema-v1 root layout, a missing `operation_spelling`, an unknown provenance item, and a `correct` row without `oracle = "mpfr_arb"` and `max_ulp = 0`.
- [ ] **Step 2: Run the new infrastructure test red**
  - Run `zig build test` in `~/projects/zmath-testing`.
  - Expected: parser or fixture failure because the checked-in manifest is schema v1.
- [ ] **Step 3: Make the specification unambiguous**
  - Replace the contradictory statements that both mark five rows correct and require every exp row to remain unsupported.
  - Define exact `operation_spelling`, provenance, and target/mode allowlist semantics.
- [ ] **Step 4: Replace the checked-in manifest with sorted schema-v2 rows**
  - Emit all 315 descriptors in deterministic order.
  - Keep every status `unsupported`, every oracle `none`, and every `max_ulp` zero until later evidence tasks change one row.
- [ ] **Step 5: Verify the contract gate**
  - Run `zig build test` in `~/projects/zmath-testing`.
  - Expected: the parser accepts the strict schema-v2 manifest and rejects malformed variants.

### Task 2: Implement strict schema-v2 parsing and artifact loading

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

- [ ] **Step 1: Write strict parser tests**
  - Test missing, unknown, duplicate, out-of-order, and duplicate function keys.
  - Test every descriptor form used by the 315 rows and reject unsupported formats, malformed widths, non-lowercase hexadecimal, duplicate corpus tuples, and bounded rows without a nonzero declared bound.
- [ ] **Step 2: Run infrastructure tests red**
  - Run `zig build test` in `~/projects/zmath-testing`.
  - Expected: absent schema-v2 parser types and embedded corpus loader.
- [ ] **Step 3: Replace fixed schema-v1 capacities and operation tables**
  - Parse the manifest into bounded fixed arrays sized for the checked-in 315 rows.
  - Keep parser work allocation-free or use explicit test-only allocation with bounded input size.
- [ ] **Step 4: Embed immutable artifacts through the build graph**
  - Extend the `addWriteFiles` corpus package to copy every referenced TSV artifact.
  - Make missing artifacts a build failure.
- [ ] **Step 5: Verify parser and artifact integrity**
  - Run `zig build test` in `~/projects/zmath-testing`.
  - Expected: manifest and every embedded corpus parse, have deterministic counts, and bind to exactly one manifest row.

### Task 3: Make differential execution reproducible

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

- [ ] **Step 1: Write deterministic-runner tests**
  - Require runtime bits, stored corpus bits, MPFR candidate bits, and Arb proof bits to agree for every finite row.
  - Require special-value rows to use stored exact target bits and reject an Arb ambiguity, MPFR disagreement, wrong target/mode, or wrong ULP limit.
- [ ] **Step 2: Run the oracle test red**
  - Run `zig build oracle` in `~/projects/zmath-testing`.
  - Expected: no TSV consumer and no deterministic target/mode record.
- [ ] **Step 3: Remove random acceptance sampling**
  - Replace `std.testing.random_seed` inputs with checked-in corpus rows.
  - Keep fuzzing as a separate explicitly seeded exploratory lane; it must not alter acceptance status.
- [ ] **Step 4: Extend bridge tests before bridge families**
  - Prove H/F/D decoding, target-format conversion, NaN and signed-zero handling, MPFR/Arb disagreement, and unique-rounding failure for unary Exp/Log/Sin.
  - Do not add C23 binary/ternary bridge operations in this task.
- [ ] **Step 5: Record corpus provenance**
  - Document producer command, pinned inputs, output class policy, tuple selection, counts, and allowlisted target/mode combinations in `test/corpus/README.md`.
- [ ] **Step 6: Verify reproducibility**
  - Run `zig build oracle` twice with the same explicit supported target and optimize mode.
  - Expected: identical row count, selected corpus, and pass/fail result.

### Task 4: Re-establish the initial five accepted unary rows

**Files:**
- Modify: `~/projects/zmath-testing/test/corpus/manifest.toml`
- Modify: `~/projects/zmath-testing/test/oracle/corpus_runner.zig`
- Create: `~/projects/zmath-testing/test/corpus/log_h.tsv`
- Create: `~/projects/zmath-testing/test/corpus/sin_h.tsv`
- Modify: `~/projects/zmath-testing/test/corpus/README.md`

**Interfaces:**
- Re-accept only `exp/h`, `exp/f`, `exp/d`, `log/h`, and `sin/h` after Task 3 proves their checked-in artifacts.
- Each row uses `status = "correct"`, `oracle = "mpfr_arb"`, `special_policy = "zmath_c23_v1"`, and `max_ulp = 0`.

- [ ] **Step 1: Generate immutable corpus artifacts**
  - Enumerate all 65,536 H inputs for `exp/h`, `log/h`, and `sin/h`.
  - Add F/D Exp special-value, underflow/overflow, transition, midpoint, adjacent, and hard-case rows with exact bits.
- [ ] **Step 2: Run each candidate corpus against MPFR and Arb red**
  - Keep every manifest row unsupported while any stored row, MPFR result, or Arb enclosure disagrees.
- [ ] **Step 3: Fix runtime only when an oracle row proves a defect**
  - Touch the exact kernel or wrapper owning the mismatched result; preserve pure-value behavior and caller-RNE requirements.
- [ ] **Step 4: Change one manifest row per green proof**
  - Transition `exp/h`, `exp/f`, `exp/d`, `log/h`, and `sin/h` independently. Do not batch statuses before each operation-format corpus passes.
- [ ] **Step 5: Verify the accepted corpus matrix**
  - Run `zig build test` in `zmath`, then `zig build oracle` in `zmath-testing` for every allowlisted target/mode.
  - Expected: every accepted row has exact zero ULP evidence and every remaining row remains unsupported.

### Task 5: Accept `log/f` as the first new pair

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

- [ ] **Step 1: Add a red immutable `log/f` corpus test**
  - Include signed zero, subnormal-normal boundaries, `-1`, the domain interval below zero, infinities, signaling and quiet NaNs, adjacent values around one, rounding midpoints, and existing CORE-MATH fallback vectors.
- [ ] **Step 2: Run the log/f oracle lane red**
  - Run `zig build oracle` in `~/projects/zmath-testing`.
  - Expected: the manifest remains unsupported or the runner rejects `log/f` dispatch.
- [ ] **Step 3: Add `log/f` corpus dispatch**
  - Route only `Operation.log` plus `Format.f` to `zmath.scalar.logf`.
  - Compare runtime, stored, MPFR, and Arb results; reject any mismatch.
- [ ] **Step 4: Repair only oracle-proven runtime defects**
  - Modify `src/scalar.zig` only if a checked-in failing corpus row demonstrates a mismatch.
  - Add the exact regression input to `src/scalar_log_test.zig` before changing the kernel.
- [ ] **Step 5: Transition exactly one row**
  - Set only `log/f` to `correct`, `mpfr_arb`, and zero ULP after the full deterministic corpus passes.
- [ ] **Step 6: Verify no adjacent status changed**
  - Run `zig build test` in `zmath` and `zig build oracle` in `zmath-testing`.
  - Expected: `log/f` is accepted; `log/d`, `sin/f`, and `sin/d` stay unsupported.

### Task 6: Add target qualification and scalar benchmarks

**Files:**
- Modify: `build.zig`
- Modify: `~/projects/zmath-testing/build.zig`
- Create: `~/projects/zmath-testing/bench/targets.zig`
- Create: `~/projects/zmath-testing/bench/target_info.zig`
- Modify: `~/projects/zmath-testing/bench/public_api.zig`

**Interfaces:**
```text
zig build bench -Dbench-target=scalar-log-f
zig build run-bench -Dbench-target=scalar-log-f
```

- [ ] **Step 1: Write target-policy tests**
  - Reject targets outside the explicit C-kernel ABI allowlist.
  - Reject a benchmark target without a corpus ID, checked-in corpus path, operation count, checksum, target triple, sorted feature set, and optimize mode.
- [ ] **Step 2: Run benchmark target parsing red**
  - Run `zig build bench` in `~/projects/zmath-testing`.
  - Expected: no scalar target table exists.
- [ ] **Step 3: Make target restrictions explicit**
  - Document and enforce the supported FMA-capable x86/x86_64 and aarch64 configurations.
  - Add compile/ABI probes for `_Float16`, 128-bit integer support, endianness assumptions, and every linked C symbol before declaring an allowlisted target.
- [ ] **Step 4: Add deterministic `scalar-log-f` benchmarking**
  - Use the immutable `log_f.tsv` corpus, fixed iteration count, checksum, and no MPFR/Arb linkage.
  - Print corpus ID, input count, target triple, resolved CPU features, optimize mode, and checksum.
- [ ] **Step 5: Measure exact-power resource bounds separately**
  - Add ReleaseSafe and ReleaseFast tests around minimum-subnormal inputs and degree ±64 for `pown`, `compoundn`, and `rootn`.
  - Record stack and latency evidence before accepting any power family.
- [ ] **Step 6: Verify target-qualified release lanes**
  - Run local tests, oracle tests, and the scalar benchmark for each allowlisted target/mode.
  - Expected: runtime remains dependency-free; external numerical libraries link only into oracle artifacts.

## Deferred Scope

- SIMD dispatch, AVX2/AVX-512/NEON kernels, and native CPU benchmarking.
- Complex numbers, vectors, quaternions, matrices, linear algebra, geometry, and algebra modules.
- Oracle bridge support for binary, ternary, integer-parameter, and C-backed C23 families.

These begin only after Tasks 1-6 establish deterministic scalar acceptance, target qualification, and a measured exact-power resource envelope.

## Plan Self-Review

- **Spec coverage:** contract contradictions, schema v2, all 315 row descriptors, corpus artifacts, MPFR/Arb proof, target/mode evidence, benchmark metadata, and exact-power resource limits map to Tasks 1-6.
- **Dependency order:** no status transition precedes contract reconciliation, parser migration, immutable corpus loading, deterministic oracle execution, or target qualification.
- **Scope:** the plan preserves the current dependency-free runtime boundary and excludes premature SIMD and planned non-scalar modules.
