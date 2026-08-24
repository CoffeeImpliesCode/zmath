# zmath

**Stack**: Zig 0.16.0 · SIMD / Vector math · `jj` (colocated git)  
**Status**: Active — Foundational high-performance mathematical primitives

---

## 1. Overview & Architecture

`zmath` provides deterministic, high-throughput scalar and SIMD mathematical functions across double, single, and half-precision floating point representations.

```
├── src/
│   ├── root.zig          # Root module entry point
│   ├── scalar/           # Scalar elementary transcendental functions
│   ├── simd/             # Vectorized AVX2 / AVX-512 / NEON kernels
│   └── tests/            # High-precision accuracy & boundary test suites
├── build.zig             # Build definition & benchmark targets
└── build.zig.zon         # Package metadata
```

---

## 2. Invariants & Mathematical Contracts

1. **Numerical Accuracy & Edge Cases**:
   - Every function must handle IEEE-754 special values (NaN, $+0$, $-0$, $+\infty$, $-\infty$, subnormals) with strict, documented behavior.
   - Max Error (ULP) bounds must be verified against high-precision reference oracles (`mpfr`/`quadmath`).
2. **SIMD & Fallbacks**:
   - Vectorized routines must maintain identical numerical results to scalar fallback paths.
   - Explicit CPU target feature gating (e.g. AVX2/FMA vs baseline SSE2).
3. **TigerStyle & Explicitness**:
   - Zero dynamic allocations inside mathematical evaluation loops. All scratch memory must be statically bounded or caller-supplied.

---

## 3. Commands & Verification

```bash
zig build                      # Build library
zig build test                 # Run unit and accuracy tests
zig build bench                # Run throughput benchmark (baseline)
zig build bench -Dcpu=native   # Benchmark with native vector extensions
zig fmt --check .              # Verify formatting
gate                           # Run house gate script
jj st                          # VCS status (NEVER raw git)
```

## NixOS, Development Environment & Shared Notes

- This host uses NixOS. Never install system or development libraries imperatively.
- Declare dependencies, headers, runtime libraries, and tools in `devenv.nix`.
- Use the `devenv.nix` environment through direnv or devenv.
- Run `direnv allow` after an environment change.
- Use `devenv shell` for a one-off shell.
- Store durable project architecture, design, protocol, and experiment notes under `docs/`.
- Store cross-project or private agent handoff notes in `~/.agents/agent-notes/`. Never commit these notes.
