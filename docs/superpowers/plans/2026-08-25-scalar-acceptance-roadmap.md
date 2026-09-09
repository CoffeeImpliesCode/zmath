# Historical Scalar Acceptance Roadmap (superseded)

> **HISTORICAL / SUPERSEDED — not a live plan.** This August 2026 schema-2 snapshot had six accepted rows. It records past tasks, commands, counts, and verification. Keep its recorded evidence and checkbox states as history. Do not execute its steps or use its row statuses as current acceptance policy.
>
> The live sources are the [API contract](../../../API.md), [schema-3 scalar contract](../../../spec/core-math-2026-08-22.md), and sibling [corpus manifest](../../../../zmath-testing/test/corpus/manifest.toml). Use the API contract for current `zmath-testing` steps; `zmath` has no `bench` step.

**Historical goal:** Convert zmath's broad C23 scalar implementation from reflection-tested API surface into a reproducible, manifest-backed, target-qualified numerical release foundation.

**Historical architecture:** `zmath` remains a dependency-free runtime. The sibling `zmath-testing` project owns manifest parsing, immutable corpora, MPFR/Arb production and proof, and benchmark artifacts. This plan proposed schema v2 and six incremental acceptances before non-scalar work.

**Historical tech stack:** Zig 0.16.0, Jujutsu, Nix/devenv, MPFR 512-bit, Arb 512-bit via FLINT, vendored CORE-MATH C kernels.

## Global Constraints (historical snapshot)

- Runtime `zmath` must not link MPFR, GMP, FLINT, Arb, BLIS, Eigen, C++, libc math, `fenv`, or `errno` behavior.
- Public floating formats are only `f16`, `f32`, and `f64`; `q`, `f128`, and bfloat16 remain rejected.
- Public scalar suffixes are `h`, `f`, and `d`; generic root functions reject mismatched or unsupported floating types at compile time.
- Every acceptance transition requires exact special-value behavior, caller round-to-nearest/ties-to-even with gradual underflow, MPFR `MPFR_RNDN` at exactly 512 bits, and a unique matching Arb rounding proof at exactly 512 bits.
- Every corpus and benchmark input set is repository-relative, immutable, deterministic, and has an explicit selection rule and count.
- `unsupported` means no acceptance claim. Export presence, local unit tests, source provenance, or a passing build never changes row status.
- The then-accepted manifest had 315 sorted schema-v2 rows. Six pairs were accepted: `exp/h`, `exp/f`, `exp/d`, `log/h`, `log/f`, and `sin/h`. The other 309 rows were `unsupported` at that time.
- That corpus profile was target `x86_64-linux-gnu.2.38`, CPU `haswell`, and mode `ReleaseSafe`. These historical statuses do not qualify any other profile.
- Implementation, manifest acceptance, corpus evidence, and whole-domain proof are separate facts. F and D selected corpora are not exhaustive whole-domain proof.
- Under schema v2, `bounded` rows permitted any non-negative `max_ulp`, including zero. `correct` rows required `max_ulp = 0`.

## Recorded integration state (August 2026)

At that time, the source, strict schema-v2 manifest, corpus loader, deterministic replay lane,
and benchmark metadata for the six accepted pairs were recorded as integrated.
This roadmap recorded their re-verification, not the current manifest status.
It did not invent historical red runs or retroactive failure evidence.
All six tasks had execution evidence recorded below.
Exact-power measurements were observations, not full-domain resource bounds.

The broad C23 campaign and evaluator branch were then recorded at the
`preserve-c23-campaigns` bookmark. This does not describe the branch's current state.
Only strict manifest and schema rejection regressions entered this historical lane.
The obsolete approximate `atan2` split and its invalid audit constants had been superseded by the canonical implementation and migrated regression cases.
The then-current `lgamma` contract governed zero and finite negative-integer poles.

---

### Task 1: Reconcile the scalar acceptance contract

**Files:**
- Modify: `spec/core-math-2026-08-22.md:486-622`
- Modify: `API.md:641-658`
- Modify: `~/projects/zmath-testing/test/corpus/manifest.toml`
- Test: `~/projects/zmath-testing/test/integration/infrastructure.zig`

**Decision at the time:** Keep schema version 2 and the 315 sorted typed rows. The manifest then accepted exactly six pairs: `exp/h`, `exp/f`, `exp/d`, `log/h`, `log/f`, and `sin/h`. The other 309 rows were `unsupported`.

**Interfaces:**
- Root manifest keys are exactly: `schema_version`, `core_math_repository`, `core_math_revision`, `mpfr_precision_bits`, `arb_precision_bits`, `rounding`, `formats`, `reserved_formats`, `rejected_formats`, `operation_spelling`, `columns`, `hexadecimal`, `uncertainty_policy`, and `special_policy`.
- Every function row contains exactly: `op`, `provenance`, `status`, `arity`, `inputs`, `result`, `oracle`, `special_policy`, and `max_ulp`.
- Define the allowed `operation_spelling` value and the declared provenance vocabulary in the spec before parser implementation. The parser must reject values outside those sets.

- [x] **Step 1: Write contract regression tests**
  - Integrated strict rejection tests for the old schema-v1 layout, missing `operation_spelling`, unknown provenance, and invalid `correct` row relations.
- [x] **Step 2: Reverify the infrastructure test**
  - Run `zig build test` in `~/projects/zmath-testing`.
  - Reverification passed. No historical red result is inferred.
- [x] **Step 3: Make the specification unambiguous**
  - Reconcile the then-current six accepted rows, the 309 unsupported rows, exact `operation_spelling`, provenance, target profile, and evidence boundaries.
- [x] **Step 4: Use the checked-in sorted schema-v2 manifest**
  - Keep all 315 descriptors in deterministic order.
  - Keep the six accepted rows at `status = "correct"`, `oracle = "mpfr_arb"`, and `max_ulp = 0`.
  - Keep the remaining 309 rows at `status = "unsupported"`, `oracle = "none"`, and `max_ulp = 0`.
- [x] **Step 5: Verify the contract gate**
  - Run `zig build test` in `~/projects/zmath-testing`.

### Task 2: Re-verify strict schema-v2 parsing and artifact loading


**Files:**
- Modify: `~/projects/zmath-testing/test/support/corpus.zig`
- Modify: `~/projects/zmath-testing/test/integration/infrastructure.zig`
- Modify: `~/projects/zmath-testing/build.zig`
- Create: `~/projects/zmath-testing/test/corpus/exp_h.tsv`
- Create: `~/projects/zmath-testing/test/corpus/exp_f.tsv`

**Interfaces:**
- `corpus.TargetMode` carries `target`, `cpu`, and `mode`.
- `corpus.parseEmbeddedManifest()` returns a validated manifest.
- `corpus.loadEmbeddedArtifact()` validates the manifest binding, source, metadata, and active target mode.
- `Artifact.rows()` returns a bounded row iterator. Replay does not allocate a row array.

- [x] **Step 1: Write strict parser tests**
  - Test missing, unknown, duplicate, out-of-order, and duplicate function keys.
  - Test every descriptor form used by the 315 rows and reject unsupported formats, malformed widths, non-lowercase hexadecimal, duplicate corpus tuples, and bounded rows without a declared non-negative bound.
- [x] **Step 2: Reverify the infrastructure test**
  - Run `zig build test` in `~/projects/zmath-testing`.
  - Reverification passed. No historical red result is inferred.
- [x] **Step 3: Replace fixed schema-v1 capacities and operation tables**
  - Parse the manifest into bounded fixed arrays sized for the checked-in 315 rows.
  - Keep parser work allocation-free or use explicit test-only allocation with bounded input size.
- [x] **Step 4: Embed immutable artifacts through the build graph**
  - Extend the `addWriteFiles` corpus package to copy every referenced TSV artifact.
  - Make missing artifacts a build failure.
- [x] **Step 5: Verify parser and artifact integrity**
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
- `corpus_runner.runAccepted(active_target_mode)` replays all embedded accepted artifacts.
- Each row must match runtime bits, stored bits, the MPFR candidate, and the Arb proof or exact special-value rule.

- [x] **Step 1: Write deterministic-runner tests**
  - Require a bridge proof for every row, including finite-input overflow. Only exact nonfinite-input and domain-special rules bypass interval evaluation.
  - Integrated checks require special-value rows to use stored exact target bits and reject an Arb ambiguity, MPFR disagreement, wrong target/mode, or wrong ULP limit.
- [x] **Step 2: Reverify the oracle lane**
  - Run `zig build oracle -Dtarget=x86_64-linux-gnu.2.38 -Dcpu=haswell -Doptimize=ReleaseSafe` in `~/projects/zmath-testing`.
  - The new finite-overflow regression failed before the bridge repair; the repaired oracle lane passed 25/25.
- [x] **Step 3: Remove random acceptance sampling**
  - Replace random acceptance inputs with checked-in corpus rows.
  - Keep fuzzing as a separate explicitly seeded exploratory lane. It must not alter acceptance status.
- [x] **Step 4: Extend bridge tests before bridge families**
  - Cover H/F/D decoding, target-format conversion, NaN and signed-zero handling, MPFR/Arb disagreement, and unique-rounding failure for unary Exp/Log/Sin.
  - Do not add C23 binary, ternary, or integer-parameter bridge operations in this task.
- [x] **Step 5: Record corpus provenance**
  - Document the qualified producer and replay commands, pinned inputs, output class policy, tuple selection, counts, and target profile in `test/corpus/README.md`.
- [x] **Step 6: Reverify reproducibility**
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
- The six accepted rows in this snapshot were `exp/h`, `exp/f`, `exp/d`, `log/h`, `log/f`, and `sin/h`.
- Each used `status = "correct"`, `oracle = "mpfr_arb"`, `special_policy = "zmath_c23_v1"`, and `max_ulp = 0`.

- [x] **Step 1: Generate immutable corpus artifacts**
  - Enumerate all 65,536 H inputs for `exp/h`, `log/h`, and `sin/h`.
  - Add deterministic F/D `exp` hard-case rows and the selected `log/f` rows with exact bits.
  - Treat F/D artifacts as evidence for stored inputs, not exhaustive whole-domain proof.
- [x] **Step 2: Reverify each accepted corpus**
  - Run the qualified MPFR and Arb replay for each accepted pair.
  - Keep any disputed row unsupported until its stored row, MPFR result, and Arb enclosure agree.
  - This roadmap records no historical red result.
- [x] **Step 3: Fix runtime only when an oracle row proves a defect**
  - Touch only the exact kernel or wrapper owning a checked-in mismatch. Preserve pure-value behavior and caller-RNE requirements.
- [x] **Step 4: Keep one manifest status per accepted proof**
  - The then-current manifest accepted all six pairs independently. The plan did not change adjacent statuses when a pair changed.
- [x] **Step 5: Verify the accepted corpus matrix**
  - Run `zig build test` in `zmath`.
  - Run `zig build oracle -Dtarget=x86_64-linux-gnu.2.38 -Dcpu=haswell -Doptimize=ReleaseSafe` in `zmath-testing`.
  - At the time, confirm zero-ULP evidence for the six then-accepted rows and `unsupported` status for the other 309 rows.

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
- [x] **Step 2: Reverify the `log/f` oracle lane**
  - Run `zig build oracle -Dtarget=x86_64-linux-gnu.2.38 -Dcpu=haswell -Doptimize=ReleaseSafe` in `~/projects/zmath-testing`.
  - Reverification passed. No historical red result is inferred.
- [x] **Step 3: Add `log/f` corpus dispatch**
  - Route only `Operation.log` plus `Format.f` to `zmath.scalar.logf`.
  - Compare runtime, stored, MPFR, and Arb results. Reject any mismatch.
- [x] **Step 4: Repair only oracle-proven runtime defects**
  - Modify `src/scalar.zig` only if a checked-in failing corpus row demonstrates a mismatch.
  - Add the exact regression input to `src/scalar_log_test.zig` before changing the kernel.
- [x] **Step 5: Keep exactly one accepted row**
  - The then-current manifest set `log/f` alone to `correct`, `mpfr_arb`, and zero ULP in this pair's transition.
- [x] **Step 6: Verify no adjacent status changed**
  - Run `zig build test` in `zmath`.
  - Run `zig build oracle -Dtarget=x86_64-linux-gnu.2.38 -Dcpu=haswell -Doptimize=ReleaseSafe` in `zmath-testing`.
  - At the time, confirm `log/f` was accepted and `log/d`, `sin/f`, and `sin/d` were unsupported.

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
- [x] **Step 2: Reverify benchmark target parsing**
  - Run `zig build bench -Dbench-target=scalar-log-f -Dtarget=x86_64-linux-gnu.2.38 -Dcpu=haswell -Doptimize=ReleaseSafe` in `~/projects/zmath-testing`.
  - All six targets built and ran; unknown targets and altered profiles were rejected.
- [x] **Step 3: Make target restrictions explicit**
  - Enforce exactly target `x86_64-linux-gnu.2.38`, CPU `haswell`, and mode `ReleaseSafe` for the accepted corpus profile.
  - Do not claim support for another target, CPU, or optimize mode.
- [x] **Step 4: Add deterministic scalar benchmarking**
  - Cover `scalar-exp-h`, `scalar-exp-f`, `scalar-exp-d`, `scalar-log-h`, `scalar-log-f`, and `scalar-sin-h`.
  - Validate each immutable corpus and SHA-256 before timing. Use a fixed operation count and per-call input barriers.
  - Print corpus ID, input count, target triple, resolved CPU features, optimize mode, elapsed time, and checksum.
- [x] **Step 5: Measure exact-power resource use separately**
  - Add ReleaseSafe and ReleaseFast tests around minimum-subnormal inputs and degree ±64 for `pown`, `compoundn`, and `rootn`.
  - Record observed stack and native batch latency before accepting any power family.
  - State measurement scope. A sampled peak is not a full-domain resource bound.
- [x] **Step 6: Verify target-qualified release lanes**
  - Run local tests, the qualified oracle command, and the qualified scalar benchmark.
  - Confirm runtime remains dependency-free. Confirm external numerical libraries link only into oracle artifacts.

## Recorded verification (August 2026 snapshot, not current)

Toolchain: Zig 0.16.0. Accepted corpus profile:
`x86_64-linux-gnu.2.38`, exact populated `haswell` features, `ReleaseSafe`.

| Check | Result |
| --- | --- |
| Library tests, native Debug | 520/520 |
| Library tests, qualified target/CPU, ReleaseSafe | 520/520 |
| Library tests, same target/CPU, ReleaseFast | 520/520; not an additional corpus profile |
| Harness tests, native Debug | 42/42 |
| Harness `test oracle oracle-fuzz`, qualified profile | 68/68: 42 ordinary, 25 oracle, 1 exploratory smoke |
| Repeated repaired oracle replay | 25/25 on both runs |
| Six-corpus regeneration with full Arb proof | All six SHA-256 values unchanged; 196,883 rows |
| Manifest after regeneration | 315 rows: six `correct`, 309 `unsupported` |
| Six scalar benchmark targets | Built and executed with fixed counts, hashes, and checksums |
| Bounded Exp fuzz, seed `0x7a6d6174`, requested budget 100 | 112 runs; no crash |
| Producer wrong libc, altered Haswell features, or ReleaseFast | Rejected before artifact replacement |
| Library without FMA; unknown benchmark; unqualified benchmark/oracle | Rejected |
| Formatting and library build | Passed |

Recorded red-to-green evidence from August 2026:

- Missing embedded-global metadata was accepted before its guard was restored.
- The finite-input overflow regression failed with `UNCERTAIN`; the repaired
  bridge proved the exact overflow rounding interval at 512 bits.
- A forced `fclose` failure exposed a second close attempt. The repaired
  producer reported `FileWriteFailed` and left no temporary or final artifact.
- The exploratory target first failed to link, then exposed incompatible C/Zig
  coverage tables. Its Exp-only Zig campaign ran with a separately linked,
  uninstrumented oracle bridge and made no C-kernel coverage claim.

Power measurements covered 36 subnormal/degree-±64 scenarios in each release mode.
Outputs and checksums agreed. Massif observed whole-program stack peaks of
362,832 bytes in ReleaseSafe and 362,976 bytes in ReleaseFast, with zero heap
bytes. Native four-call `compoundn` batches reached 14.678 seconds and 3.233
seconds respectively on a shared host. These are samples, not latency bounds
or a controlled comparison between modes. Power rows were unsupported in this snapshot.

The detailed commands, benchmark samples, counts, and measurement caveats were
recorded in `zmath-testing/test/corpus/README.md` at the time. The eight measured
runtime executables had no dynamic section and no MPFR/GMP/Arb/FLINT symbols.

## Deferred Scope (historical snapshot, not current)

- SIMD dispatch, AVX2/AVX-512/NEON kernels, and native CPU benchmarking were deferred in this plan.
- Complex numbers, vectors, quaternions, matrices, linear algebra, geometry, and algebra modules were deferred in this plan.
- Oracle bridge support for binary, ternary, integer-parameter, and C-backed C23 families was deferred at the time. This is not a claim about the current harness.
- Whole-domain proof for F and D scalar families was deferred. Their selected corpora were not exhaustive.

At the time, the broad C23 campaign and evaluator branch were recorded at the
`preserve-c23-campaigns` bookmark.
Only strict manifest and schema rejection regressions entered this lane.
This plan excluded descriptor runtime and broad campaign changes.
It conditioned deferred work on re-verification of the six then-accepted scalar lanes
and measured evidence for exact-power resource use. See the current contract and manifest above instead.

## Plan Self-Review (historical)

- **Spec coverage:** schema v2, all 315 row descriptors, six accepted pairs, 309 unsupported rows, corpus artifacts, MPFR/Arb proof, target/mode metadata, benchmark metadata, and exact-power resource limits mapped to Tasks 1-6 at the time.
- **Evidence boundaries:** implementation, manifest acceptance, stored corpus evidence, and whole-domain proof remain separate. F and D samples do not establish exhaustive proof.
- **Dependency order:** the plan placed contract reconciliation, parser migration, immutable corpus loading, deterministic oracle execution, and target qualification before later status transitions.
- **Scope:** this plan preserved the dependency-free runtime boundary and excluded SIMD, broad C23 campaign migration, and non-scalar modules.

