const std = @import("std");

// One compile-time policy for Zig intrinsics and the C kernel flags. FMA is
// part of the base floating-point ISA on RISC-V, LoongArch, PowerPC and s390x;
// absence of an x86-specific predefine does not imply a software operation.
pub fn hasNativeFma32(cpu: std.Target.Cpu) bool {
    return switch (cpu.arch.family()) {
        .x86 => cpu.has(.x86, .fma),
        // 32-bit ARM reaches no fused instruction, so it owns none. On
        // `-mcpu=baseline`, which resolves to `.fpu neon` and therefore to
        // Zig's vfp4, vfp4d16, vfp4sp and vfp4d16sp at once, the recorded
        // lowering turns `__builtin_fmaf` into a tail call to the C library,
        // `b fmaf`, in every 32-bit ARM entry. A reduced-register variant
        // that carries the instruction the toolchain never selects is not a
        // native fused target, so soft-float and VFPv4 cores alike keep the
        // exact integer dyadic engine.
        .arm => false,
        .aarch64 => cpu.has(.aarch64, .fp_armv8),
        .hexagon => cpu.has(.hexagon, .v5),
        .loongarch => cpu.has(.loongarch, .f),
        .powerpc => cpu.has(.powerpc, .hard_float) and !cpu.has(.powerpc, .spe),
        .riscv => cpu.has(.riscv, .f),
        .s390x => !cpu.has(.s390x, .soft_float),
        else => false,
    };
}

/// Binary64 needs double-precision fused instructions: Hexagon v5 fuses only
/// single-precision products, and the other listed families require their
/// double-precision floating-point features. 32-bit ARM requires none of
/// them, because its recorded lowering is a library call for both widths.
pub fn hasNativeFma64(cpu: std.Target.Cpu) bool {
    return switch (cpu.arch.family()) {
        .x86 => cpu.has(.x86, .fma),
        // The same four entries end in `b fma`, so binary64 on 32-bit ARM is
        // a library call too and takes the exact path.
        .arm => false,
        .aarch64 => cpu.has(.aarch64, .fp_armv8),
        .loongarch => cpu.has(.loongarch, .d),
        .powerpc => cpu.has(.powerpc, .hard_float) and !cpu.has(.powerpc, .spe),
        .riscv => cpu.has(.riscv, .d),
        .s390x => !cpu.has(.s390x, .soft_float),
        else => false,
    };
}
