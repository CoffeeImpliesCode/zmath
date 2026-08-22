# zmath
Zig 0.16 project (house scaffold). Foundational mathematical library

## Stack
- **Language**: Zig 0.16 (module API, `std.Io`, Juicy Main)
- **VCS**: jj — DO NOT USE GIT COMMANDS
- **Dev env**: devenv + direnv (`.envrc` → `use devenv`)
- **Linting**: zlint (`zlint.json` rules) + `zig fmt --check` + `ziglint` + `zbc`
- **Toolchain**: `zig` from anyzig
- **Lint/analyze tools**: prefer the locally built versions (`~/.local/bin/zlint`, `~/tools/ziglint`, `~/tools/zbc`, `~/tools/zcov`) over devenv's nixpkgs ones — they may be patched for this setup.

## Structure

## Commands

```bash
direnv allow        # enter devenv shell (first time / after devenv.nix edits)
zig build           # Debug build → zig-out/bin/ztex
zig build test      # unit + integration + regression tests (all modules)
zig build bench     # ReleaseFast benchmark (parse+layout+flatten throughput)
zig build bench -Dcpu=native   # unlock AVX2/512 (default = baseline SSE)
zig-cov test        # line coverage (installed zig-cov; ~49% baseline)
zig fmt --check .   # formatting gate
zlint               # lint (prefer ~/.local/bin/zlint; devenv nixpkgs is fallback)
gate                # all gates: fmt --check, build, test, ast-check, ziglint, zlint (local ~/.local/bin); zbc by default (advisory; strict via ZBC_STRICT=1); optional `--cov`; exits 1 on issues
hyperfine 'zig build bench'  # stable cross-run timing
jj st               # working-copy status (never git)
```

## Conventions

- Tiger style: explicit allocators, no hidden allocations, `defer`/`errdefer`
  cleanup, named `gpa`/`arena` allocator params
- Zig 0.16 idioms only: `std.Io` threaded through I/O, unmanaged containers
  with explicit allocator args, `.empty` decl literals, `@import("root")`
- Error sets explicit; no `anyerror` leaks; no `_ =` discards
- `zig fmt` before finishing any change; lines ≤ 100 columns
- Unit tests live next to code in `src/`; run via `zig build test`
- No backwards compatibility: delete dead code, change directly

## Version control workflow (jj)

- Start every task from a CLEAN working state: run `jj st` first and confirm
  no uncommitted changes (nothing staged, nothing modified). If dirty, commit
  or stash the pending work BEFORE starting new work.
- One logical change per commit: after each fix or feature, immediately
  `jj commit -m "..."` with a focused message describing the change. NEVER
  accumulate uncommitted changes across multiple tasks.
- Commit only green work (gates pass: fmt, build, test, zlint). Never commit
  a red tree.
- Never use git commands.

## Development

- Use tools, always use `edit`/`write` for file changes — NEVER patch files
  with inline python/perl/sed/awk scripts, heredocs, or `echo >` redirection
- Read with `read`, search with `grep`/`glob`; shell is for builds/tests and
  one-shot fact queries only
- NEVER use web_* tools (broken) — for web search use the local server at
  `localhost:4444`; prefer reading primary sources directly
- Delegate aggressively: many small subagents over few big ones; give full
  instructions; assume subagents are dumb
- Parallelize: N independent files → N workers
