{ pkgs, ... }:

rec {
  packages = with pkgs; [
    # bash
    # pkg-config

    hyperfine        # CLI timing of benchmarks
    valgrind         # memory / cachegrind analysis
    perf-tools       # linux perf profiler

  ];

  env = {
    LD_LIBRARY_PATH = "${pkgs.lib.makeLibraryPath packages}";
    INCLUDE_DIRECTORIES = "${pkgs.lib.makeIncludePath packages}";
  };

  # `zig-gate` is a host tool in ~/.local/bin, not a nix package, so it is not
  # in `packages`; the exported environment reaches it through PATH instead.
  # Where /nix/store holds this devshell, the store script already runs it.
  scripts.gate.exec = "zig-gate";
}
