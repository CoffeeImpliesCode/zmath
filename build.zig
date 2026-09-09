const std = @import("std");

// Although this function looks imperative, it does not perform the build
// directly and instead it mutates the build graph (`b`) that will be then
// executed by an external runner. The functions in `std.Build` implement a DSL
// for defining build steps and express dependencies between them, allowing the
// build runner to parallelize the build automatically (and the cache system to
// know when a step doesn't need to be re-run).
fn requireCKernelTarget(target: std.Build.ResolvedTarget) void {
    switch (target.result.cpu.arch) {
        .x86, .x86_64 => {
            if (!target.result.cpu.has(.x86, .fma)) {
                @panic(
                    "zmath C kernels require target feature +fma on x86/x86_64",
                );
            }
        },
        .aarch64 => {},
        else => @panic(
            "zmath C kernels support only x86/x86_64 with +fma and little-endian aarch64",
        ),
    }
}

pub fn build(b: *std.Build) void {
    // The C kernels have target-specific ABI and FMA assumptions. Reject an
    // unsupported target before constructing artifacts so every build step
    // uses the same deterministic policy.
    const target = b.standardTargetOptions(.{});
    requireCKernelTarget(target);
    // Standard optimization options allow the person running `zig build` to select
    // between Debug, ReleaseSafe, ReleaseFast, and ReleaseSmall. Here we do not
    // set a preferred release mode, allowing the user to decide how to optimize.
    const optimize = b.standardOptimizeOption(.{});
    // It's also possible to define more custom flags to toggle optional features
    // of this build script using `b.option()`. All defined flags (including
    // target and optimize options) will be listed when running `zig build --help`
    // in this directory.

    // This creates a module, which represents a collection of source files alongside
    // some compilation options, such as optimization mode and linked system libraries.
    // Zig modules are the preferred way of making Zig code available to consumers.
    // addModule defines a module that we intend to make available for importing
    // to our consumers. We must give it a name because a Zig package can expose
    // multiple modules and consumers will need to be able to specify which
    // module they want to access.
    const mod = b.addModule("zmath", .{
        // The root source file is the "entry point" of this module. Users of
        // this module will only be able to access public declarations contained
        // in this file, which means that if you have declarations that you
        // intend to expose to consumers that were defined in other files part
        // of this module, you will have to make sure to re-export them from
        // the root file.
        .root_source_file = b.path("src/root.zig"),
        // Later on we'll use this module as the root module of a test executable
        // which requires us to specify a target.
        .target = target,
        .optimize = optimize,
    });
    const core_math_erfc_flags = &.{
        "-ffreestanding",
        "-fno-math-errno",
        "-ffp-contract=off",
    };
    const c_kernel_abi_guard = b.addWriteFiles().add("c_kernel_abi_guard.c",
        \\#include <stdint.h>
        \\#if !defined(__FLT16_MANT_DIG__)
        \\#error "zmath C kernels require compiler support for _Float16"
        \\#else
        \\_Static_assert(sizeof(_Float16) == 2, "zmath C kernels require 16-bit _Float16");
        \\_Static_assert(__FLT16_MANT_DIG__ == 11, "zmath C kernels require IEEE-754 binary16 precision");
        \\_Static_assert(__FLT16_MAX_EXP__ == 16, "zmath C kernels require IEEE-754 binary16 exponent range");
        \\#endif
        \\
        \\#if !defined(__SIZEOF_INT128__)
        \\#error "zmath C kernels require compiler support for 128-bit integers"
        \\#else
        \\_Static_assert(sizeof(__int128) == 16, "zmath C kernels require signed 128-bit integers");
        \\_Static_assert(sizeof(unsigned __int128) == 16, "zmath C kernels require unsigned 128-bit integers");
        \\#endif
        \\
        \\#if !defined(__BYTE_ORDER__) || !defined(__ORDER_LITTLE_ENDIAN__)
        \\#error "zmath C kernels require compiler byte-order macros"
        \\#else
        \\_Static_assert(__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__, "zmath C kernels require little-endian targets");
        \\#endif
        \\
        \\extern _Float16 cr_atan2f16(_Float16, _Float16);
        \\extern float zmath_core_atan2f_finite(float, float);
        \\extern int zmath_core_atan2d_checked(double, double, double *);
        \\extern int zmath_core_atan2d_forced_ambiguity(double, double, double *);
        \\extern void zmath_core_atan2d_test_subtract_two_limbs(uint64_t *, uint64_t *, uint64_t *);
        \\extern _Float16 cr_cospif16(_Float16);
        \\extern float cr_cospif(float);
        \\extern double cr_cospi(double);
        \\extern _Float16 cr_erfcf16(_Float16);
        \\extern float cr_erfcf(float);
        \\extern double zmath_core_erfcd_negative(double);
        \\extern double zmath_core_erfcd_positive(double);
        \\extern _Float16 zmath_core_lgammah(_Float16, int *);
        \\extern float zmath_core_lgammaf(float, int *);
        \\extern double zmath_core_lgammad(double, int *);
        \\extern _Float16 cr_log10p1f16(_Float16);
        \\extern float zmath_core_log10p1f_finite(float);
        \\extern double zmath_core_log10p1_finite(double);
        \\extern _Float16 cr_log2p1f16(_Float16);
        \\extern float zmath_core_log2p1f_finite(float);
        \\extern double zmath_core_log2p1_finite(double);
        \\extern _Float16 zmath_core_asinpih(_Float16);
        \\extern float zmath_core_asinpif(float);
        \\extern double zmath_core_asinpid_finite_bits(double, uint64_t);
        \\extern _Float16 zmath_core_acospih(_Float16);
        \\extern float zmath_core_acospif(float);
        \\extern double zmath_core_acospid(double);
        \\extern _Float16 zmath_core_atanpih(_Float16);
        \\extern float zmath_core_atanpif(float);
        \\extern double zmath_core_atanpid(double);
        \\extern _Float16 zmath_core_atan2pih(_Float16, _Float16);
        \\extern float zmath_core_atan2pif(float, float);
        \\extern double zmath_core_atan2pid_finite_bits(double, double, uint64_t, uint64_t);
        \\extern int zmath_core_atan2pid_forced_ambiguity(double, double, double *);
        \\extern _Float16 cr_powf16(_Float16, _Float16);
        \\extern float cr_powfFiniteBits(uint32_t, float);
        \\extern double cr_powFinite(double, double);
        \\extern _Float16 cr_tgammaf16(_Float16);
        \\extern float cr_tgammaf(float);
        \\extern double cr_tgamma_finite(double);
        \\
        \\typedef void (*zmath_c_kernel_symbol)(void);
        \\__attribute__((used, retain))
        \\static zmath_c_kernel_symbol const zmath_c_kernel_symbols[] = {
        \\    (zmath_c_kernel_symbol)cr_atan2f16,
        \\    (zmath_c_kernel_symbol)zmath_core_atan2f_finite,
        \\    (zmath_c_kernel_symbol)zmath_core_atan2d_checked,
        \\    (zmath_c_kernel_symbol)zmath_core_atan2d_forced_ambiguity,
        \\    (zmath_c_kernel_symbol)zmath_core_atan2d_test_subtract_two_limbs,
        \\    (zmath_c_kernel_symbol)cr_cospif16,
        \\    (zmath_c_kernel_symbol)cr_cospif,
        \\    (zmath_c_kernel_symbol)cr_cospi,
        \\    (zmath_c_kernel_symbol)cr_erfcf16,
        \\    (zmath_c_kernel_symbol)cr_erfcf,
        \\    (zmath_c_kernel_symbol)zmath_core_erfcd_negative,
        \\    (zmath_c_kernel_symbol)zmath_core_erfcd_positive,
        \\    (zmath_c_kernel_symbol)zmath_core_lgammah,
        \\    (zmath_c_kernel_symbol)zmath_core_lgammaf,
        \\    (zmath_c_kernel_symbol)zmath_core_lgammad,
        \\    (zmath_c_kernel_symbol)cr_log10p1f16,
        \\    (zmath_c_kernel_symbol)zmath_core_log10p1f_finite,
        \\    (zmath_c_kernel_symbol)zmath_core_log10p1_finite,
        \\    (zmath_c_kernel_symbol)cr_log2p1f16,
        \\    (zmath_c_kernel_symbol)zmath_core_log2p1f_finite,
        \\    (zmath_c_kernel_symbol)zmath_core_log2p1_finite,
        \\    (zmath_c_kernel_symbol)zmath_core_asinpih,
        \\    (zmath_c_kernel_symbol)zmath_core_asinpif,
        \\    (zmath_c_kernel_symbol)zmath_core_asinpid_finite_bits,
        \\    (zmath_c_kernel_symbol)zmath_core_acospih,
        \\    (zmath_c_kernel_symbol)zmath_core_acospif,
        \\    (zmath_c_kernel_symbol)zmath_core_acospid,
        \\    (zmath_c_kernel_symbol)zmath_core_atanpih,
        \\    (zmath_c_kernel_symbol)zmath_core_atanpif,
        \\    (zmath_c_kernel_symbol)zmath_core_atanpid,
        \\    (zmath_c_kernel_symbol)zmath_core_atan2pih,
        \\    (zmath_c_kernel_symbol)zmath_core_atan2pif,
        \\    (zmath_c_kernel_symbol)zmath_core_atan2pid_finite_bits,
        \\    (zmath_c_kernel_symbol)zmath_core_atan2pid_forced_ambiguity,
        \\    (zmath_c_kernel_symbol)cr_powf16,
        \\    (zmath_c_kernel_symbol)cr_powfFiniteBits,
        \\    (zmath_c_kernel_symbol)cr_powFinite,
        \\    (zmath_c_kernel_symbol)cr_tgammaf16,
        \\    (zmath_c_kernel_symbol)cr_tgammaf,
        \\    (zmath_c_kernel_symbol)cr_tgamma_finite,
        \\};
        \\
    );
    mod.addCSourceFile(.{
        .file = c_kernel_abi_guard,
        .flags = core_math_erfc_flags,
    });
    mod.addCSourceFile(.{
        .file = b.path("src/coremath_erfch.c"),
        .flags = core_math_erfc_flags,
    });
    mod.addCSourceFile(.{
        .file = b.path("src/coremath_erfcf.c"),
        .flags = core_math_erfc_flags,
    });
    mod.addCSourceFile(.{
        .file = b.path("src/coremath_erfcd.c"),
        .flags = core_math_erfc_flags,
    });
    mod.addCSourceFile(.{
        .file = b.path("src/coremath_tgammah.c"),
        .flags = core_math_erfc_flags,
    });
    mod.addCSourceFile(.{
        .file = b.path("src/coremath_tgammaf.c"),
        .flags = core_math_erfc_flags,
    });
    mod.addCSourceFile(.{
        .file = b.path("src/coremath_tgammad.c"),
        .flags = core_math_erfc_flags,
    });
    mod.addCSourceFile(.{
        .file = b.path("src/coremath_cospih.c"),
        .flags = core_math_erfc_flags,
    });
    mod.addCSourceFile(.{
        .file = b.path("src/coremath_cospif.c"),
        .flags = core_math_erfc_flags,
    });
    mod.addCSourceFile(.{
        .file = b.path("src/coremath_cospid.c"),
        .flags = core_math_erfc_flags,
    });
    mod.addCSourceFile(.{
        .file = b.path("src/coremath_log2p1h.c"),
        .flags = core_math_erfc_flags,
    });
    mod.addCSourceFile(.{
        .file = b.path("src/coremath_log2p1f.c"),
        .flags = core_math_erfc_flags,
    });
    mod.addCSourceFile(.{
        .file = b.path("src/coremath_log2p1d.c"),
        .flags = core_math_erfc_flags,
    });
    mod.addCSourceFile(.{
        .file = b.path("src/coremath_log10p1h.c"),
        .flags = core_math_erfc_flags,
    });
    mod.addCSourceFile(.{
        .file = b.path("src/coremath_log10p1f.c"),
        .flags = core_math_erfc_flags,
    });
    mod.addCSourceFile(.{
        .file = b.path("src/coremath_log10p1d.c"),
        .flags = core_math_erfc_flags,
    });
    mod.addCSourceFile(.{
        .file = b.path("src/coremath_powh.c"),
        .flags = core_math_erfc_flags,
    });
    mod.addCSourceFile(.{
        .file = b.path("src/coremath_powf.c"),
        .flags = core_math_erfc_flags,
    });
    mod.addCSourceFile(.{
        .file = b.path("src/coremath_powd.c"),
        .flags = core_math_erfc_flags,
    });
    mod.addCSourceFile(.{
        .file = b.path("src/coremath_atan2h.c"),
        .flags = core_math_erfc_flags,
    });
    mod.addCSourceFile(.{
        .file = b.path("src/coremath_atan2f.c"),
        .flags = core_math_erfc_flags,
    });
    mod.addCSourceFile(.{
        .file = b.path("src/coremath_atan2d.c"),
        .flags = core_math_erfc_flags,
    });
    mod.addCSourceFile(.{
        .file = b.path("src/coremath_lgammah.c"),
        .flags = core_math_erfc_flags,
    });
    mod.addCSourceFile(.{
        .file = b.path("src/coremath_lgammaf.c"),
        .flags = core_math_erfc_flags,
    });
    mod.addCSourceFile(.{
        .file = b.path("src/coremath_lgammad.c"),
        .flags = core_math_erfc_flags,
    });
    mod.addCSourceFile(.{
        .file = b.path("src/coremath_pi_inverseh.c"),
        .flags = core_math_erfc_flags,
    });
    mod.addCSourceFile(.{
        .file = b.path("src/coremath_pi_inversef.c"),
        .flags = core_math_erfc_flags,
    });
    mod.addCSourceFile(.{
        .file = b.path("src/coremath_pi_inversed.c"),
        .flags = core_math_erfc_flags,
    });

    // Here we define an executable. An executable needs to have a root module
    // which needs to expose a `main` function. While we could add a main function
    // to the module defined above, it's sometimes preferable to split business
    // logic and the CLI into two separate modules.
    //
    // If your goal is to create a Zig library for others to use, consider if
    // it might benefit from also exposing a CLI tool. A parser library for a
    // data serialization format could also bundle a CLI syntax checker, for example.
    //
    // If instead your goal is to create an executable, consider if users might
    // be interested in also being able to embed the core functionality of your
    // program in their own executable in order to avoid the overhead involved in
    // subprocessing your CLI tool.
    //
    // If neither case applies to you, feel free to delete the declaration you
    // don't need and to put everything under a single module.
    const exe = b.addExecutable(.{
        .name = "zmath",
        .root_module = b.createModule(.{
            // b.createModule defines a new module just like b.addModule but,
            // unlike b.addModule, it does not expose the module to consumers of
            // this package, which is why in this case we don't have to give it a name.
            .root_source_file = b.path("src/main.zig"),
            // Target and optimization levels must be explicitly wired in when
            // defining an executable or library (in the root module), and you
            // can also hardcode a specific target for an executable or library
            // definition if desireable (e.g. firmware for embedded devices).
            .target = target,
            .optimize = optimize,
            // List of modules available for import in source files part of the
            // root module.
            .imports = &.{
                // Here "zmath" is the name you will use in your source code to
                // import this module (e.g. `@import("zmath")`). The name is
                // repeated because you are allowed to rename your imports, which
                // can be extremely useful in case of collisions (which can happen
                // importing modules from different packages).
                .{ .name = "zmath", .module = mod },
            },
        }),
    });

    // This declares intent for the executable to be installed into the
    // install prefix when running `zig build` (i.e. when executing the default
    // step). By default the install prefix is `zig-out/` but can be overridden
    // by passing `--prefix` or `-p`.
    b.installArtifact(exe);

    // This creates a top level step. Top level steps have a name and can be
    // invoked by name when running `zig build` (e.g. `zig build run`).
    // This will evaluate the `run` step rather than the default step.
    // For a top level step to actually do something, it must depend on other
    // steps (e.g. a Run step, as we will see in a moment).
    const run_step = b.step("run", "Run the app");

    // This creates a RunArtifact step in the build graph. A RunArtifact step
    // invokes an executable compiled by Zig. Steps will only be executed by the
    // runner if invoked directly by the user (in the case of top level steps)
    // or if another step depends on it, so it's up to you to define when and
    // how this Run step will be executed. In our case we want to run it when
    // the user runs `zig build run`, so we create a dependency link.
    const run_cmd = b.addRunArtifact(exe);
    run_step.dependOn(&run_cmd.step);

    // By making the run step depend on the default step, it will be run from the
    // installation directory rather than directly from within the cache directory.
    run_cmd.step.dependOn(b.getInstallStep());

    // Zig 0.17 development builds temporarily removed Build.args. The coverage
    // lane does not invoke `run`; retain argument forwarding where the field exists.
    if (comptime @hasField(std.Build, "args")) {
        if (@field(b, "args")) |args| {
            run_cmd.addArgs(args);
        }
    }

    // Creates an executable that will run `test` blocks from the provided module.
    // Here `mod` needs to define a target, which is why earlier we made sure to
    // set the releative field.
    const mod_tests = b.addTest(.{
        .root_module = mod,
    });

    // A run step that will run the test executable.
    const run_mod_tests = b.addRunArtifact(mod_tests);

    // The demonstration has no tests. Run the library's numerical suites.
    const test_step = b.step("test", "Run tests");
    test_step.dependOn(&run_mod_tests.step);
    const install_docs = b.addInstallDirectory(.{
        .source_dir = exe.getEmittedDocs(),
        .install_dir = .prefix,
        .install_subdir = "docs",
    });
    const docs_step = b.step("docs", "Generate documentation");
    docs_step.dependOn(&install_docs.step);

    // Just like flags, top level steps are also listed in the `--help` menu.
    //
    // The Zig build system is entirely implemented in userland, which means
    // that it cannot hook into private compiler APIs. All compilation work
    // orchestrated by the build system will result in other Zig compiler
    // subcommands being invoked with the right flags defined. You can observe
    // these invocations when one fails (or you pass a flag to increase
    // verbosity) to validate assumptions and diagnose problems.
    //
    // Lastly, the Zig build system is relatively simple and self-contained,
    // and reading its source code will allow you to master it.
}
