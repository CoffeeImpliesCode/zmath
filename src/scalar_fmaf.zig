const fma_impl = @import("scalar_fma.zig");

// Correctly rounded fused multiply-add for binary32 values.
// The exceptional cases are settled here so that the result never depends on
// which operand a target's fused instruction chooses to propagate; the single
// rounding itself is the shared helper's job.
// Uses @mulAdd with IEEE round-to-nearest, ties-to-even where the target CPU has
// a fused instruction, and exact integer dyadic arithmetic everywhere else.
// No host libm, no fenv mutation, no dynamic allocation.

/// Computes the correctly rounded fused multiply-add of three binary32 values.
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn fmaf(a: f32, b: f32, c: f32) f32 {
    @setFloatMode(.strict);
    if (comptime fma_impl.has_native_fma) {
        const a_bits: u32 = @bitCast(a);
        const b_bits: u32 = @bitCast(b);
        const c_bits: u32 = @bitCast(c);
        if (fma_impl.exceptionalResult(a_bits, b_bits, c_bits)) |result| return @bitCast(result);
    }

    // Otherwise, compute a*b + c with one rounding.
    return fma_impl.fma(a, b, c);
}
