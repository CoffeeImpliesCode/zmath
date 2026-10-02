// Binary64 division with exact subnormal results on software-float targets.
//
// A single IEEE division rounds the exact quotient once, and that rounding is
// correct for every normal or infinite result on every target, so this module
// keeps it. Only the subnormal range needs a second path: soft-float targets
// lower the binary64 divide to compiler-rt `__divdf3`, whose subnormal branch
// returns zero ("Flush denormals to zero"), so a subnormal quotient comes back
// there as a signed zero. The exact dyadic ratio is then packed once by the
// shared narrowing engine of `scalar_narrow_common.zig`, the same engine and
// the same single packer that `scalar_narrow.zig` already uses for a narrowing
// divide, so both paths agree bit for bit wherever the divide already works.
//
// No host libm, no fenv mutation, no dynamic allocation, no runtime dispatch.

const std = @import("std");
const builtin = @import("builtin");
const fma_target = @import("scalar_fma_target.zig");
const narrow = @import("scalar_narrow_common.zig");

const format = narrow.binary64_format;

/// Smallest binary64 magnitude that is not subnormal. Zero and every subnormal
/// magnitude compare below it, an infinite one compares above.
const minimum_normal_magnitude: u64 = @as(u64, 1) << format.fraction_bits;

/// Divide with native exceptional-value behavior and exact finite underflow.
///
/// Targets with native binary64 FMA also have a conforming hardware divide.
/// Other targets retain their divide for normal, infinite and NaN results;
/// only a finite zero or subnormal quotient needs the exact ratio packer.
pub fn div64(a: f64, b: f64) f64 {
    @setFloatMode(.strict);
    if (comptime fma_target.hasNativeFma64(builtin.cpu)) return a / b;

    const quotient: u64 = @bitCast(a / b);
    const sign_mask = narrow.Format.signMask(format);
    if ((quotient & ~sign_mask) >= minimum_normal_magnitude) return @bitCast(quotient);
    const divisor: u64 = @bitCast(b);
    if ((divisor & ~sign_mask) >= 0x7ff0_0000_0000_0000) return @bitCast(quotient);
    return @bitCast(exactQuotientBits(@as(u64, @bitCast(a)), divisor));
}

/// C ABI entry point for `zmath_core_finite_div64_bits`, which the binary64
/// atan2 quadrant calls once its quotient lands in the subnormal range. It
/// always takes the exact path, so the C caller owns the range decision.
export fn zmath_core_finite_div64_bits(a_bits: u64, b_bits: u64) callconv(.c) u64 {
    return exactQuotientBits(a_bits, b_bits);
}

/// Correctly rounded quotient of two finite binary64 patterns whose divisor is
/// nonzero. No intermediate is a float, so the result does not depend on the
/// ambient floating-point mode.
fn exactQuotientBits(a_bits: u64, b_bits: u64) u64 {
    const numerator = narrow.decode(format, a_bits);
    const denominator = narrow.decode(format, b_bits);
    std.debug.assert(numerator.class == .finite or numerator.class == .zero);
    std.debug.assert(denominator.class == .finite);

    const top = narrow.finiteDyadic(format, numerator);
    const bottom = narrow.finiteDyadic(format, denominator);
    const negative = numerator.negative != denominator.negative;
    // A zero numerator is an exact zero quotient, and a zero ratio leaves the
    // packer nothing to round.
    if (top.magnitude == 0) return narrow.zero(format, negative);
    return narrow.packRatio(
        format,
        negative,
        top.magnitude,
        bottom.magnitude,
        top.exponent - bottom.exponent,
    );
}