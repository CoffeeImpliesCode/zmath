const fma_impl = @import("scalar_fma.zig");

// Correctly rounded fused multiply-add for binary64 values.
// The exceptional cases are settled here so that the result never depends on
// which operand a target's fused instruction chooses to propagate. The single
// rounding itself, and the whole software path of a target without a binary64
// fused instruction, is the shared helper's job; this file must not screen a
// case that the helper screens again.
// Uses the target's fused instruction where it has a binary64 one, and exact
// integer dyadic arithmetic everywhere else.
// No host libm, no fenv mutation, no dynamic allocation.

/// Computes the correctly rounded fused multiply-add of three binary64 values.
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn fmad(a: f64, b: f64, c: f64) f64 {
    @setFloatMode(.strict);
    if (comptime fma_impl.has_native_fma64) {
        const a_bits: u64 = @bitCast(a);
        const b_bits: u64 = @bitCast(b);
        const c_bits: u64 = @bitCast(c);
        if (fma_impl.exceptionalResult64(a_bits, b_bits, c_bits)) |result| return @bitCast(result);
    }

    // Otherwise, compute a*b + c with one rounding.
    return fma_impl.fma64(a, b, c);
}
