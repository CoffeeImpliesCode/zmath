// Exact fused multiply-add for binary32 and binary64.
//
// `@mulAdd` is one instruction on a CPU that has one and a library call on
// every other target. The recorded portability evidence shows that library
// call is not correctly rounded for binary32 or for binary64: an exact
// midpoint plus a tiny perturbation loses the perturbation and returns the
// tie. Every target without a fused instruction therefore evaluates the
// product, the sum and the single final rounding in the integer dyadic engine
// of `scalar_narrow_common.zig`, which keeps the discarded alignment bits as a
// sticky bit and the borrow of an opposite-sign addition. Both widths share
// that engine and differ only in the format description, so a binary64
// product is the same exact product as a binary32 one over a wider significand,
// and the 128-bit headroom of the shared engine already covers it.
//
// No host libm, no fenv mutation, no dynamic allocation, no runtime dispatch.

const builtin = @import("builtin");
const narrow = @import("scalar_narrow_common.zig");
const fma_target = @import("scalar_fma_target.zig");

const format = narrow.binary32_format;
const wide_format = narrow.binary64_format;

/// Uses the same target policy that configures `coremath_fma.h` in build.zig.
pub const has_native_fma = fma_target.hasNativeFma32(builtin.target.cpu);

/// Binary64 fuses in a different set of instructions than binary32, so it
/// carries its own target policy instead of sharing the binary32 one.
pub const has_native_fma64 = fma_target.hasNativeFma64(builtin.target.cpu);

/// The results IEEE 754 settles before any rounding happens: a NaN operand,
/// taken in argument order, and the two invalid products. `null` leaves an
/// exact multiplication, an exact addition and a single rounding.
fn exceptionalBits(
    comptime width: narrow.Format,
    a_bits: u64,
    b_bits: u64,
    c_bits: u64,
) ?width.Bits() {
    const sign_mask = narrow.Format.signMask(width);
    const infinity = narrow.Format.infinityBits(width);
    const magnitude_mask = narrow.Format.fractionMask(width) | infinity;
    const quiet_bit = narrow.Format.quietBit(width);

    const a_magnitude = a_bits & magnitude_mask;
    const b_magnitude = b_bits & magnitude_mask;
    const c_magnitude = c_bits & magnitude_mask;

    if (a_magnitude > infinity) return narrow.fromBits(width, a_bits | quiet_bit);
    if (b_magnitude > infinity) return narrow.fromBits(width, b_bits | quiet_bit);
    if (c_magnitude > infinity) return narrow.fromBits(width, c_bits | quiet_bit);

    // An infinite factor times a zero factor makes the product invalid.
    if ((a_magnitude == infinity and b_magnitude == 0) or
        (a_magnitude == 0 and b_magnitude == infinity))
    {
        return narrow.fromBits(width, infinity | quiet_bit);
    }

    // An infinite product added to an infinity of the other sign is invalid.
    const product_infinite = (a_magnitude == infinity and b_magnitude != 0) or
        (b_magnitude == infinity and a_magnitude != 0);
    if (product_infinite and c_magnitude == infinity and
        ((a_bits ^ b_bits) & sign_mask != 0) != (c_bits & sign_mask != 0))
    {
        return narrow.fromBits(width, infinity | quiet_bit);
    }
    return null;
}

/// The binary32 case of `exceptionalBits`. It stays public because the binary32
/// public kernel settles these cases on a target that has a fused instruction
/// too, and the shared engine settles them on every other target.
pub fn exceptionalResult(a_bits: u32, b_bits: u32, c_bits: u32) ?u32 {
    return exceptionalBits(format, a_bits, b_bits, c_bits);
}

/// The binary64 case of `exceptionalBits`.
pub fn exceptionalResult64(a_bits: u64, b_bits: u64, c_bits: u64) ?u64 {
    return exceptionalBits(wide_format, a_bits, b_bits, c_bits);
}

/// Exact software fused multiply-add, computed on the integer significands of
/// the format alone. No intermediate is a float, so the result does not depend
/// on the ambient floating-point mode.
fn softwareFmaBits(
    comptime width: narrow.Format,
    a_bits: u64,
    b_bits: u64,
    c_bits: u64,
) width.Bits() {
    if (exceptionalBits(width, a_bits, b_bits, c_bits)) |result| return result;

    const sign_mask = narrow.Format.signMask(width);
    const infinity = narrow.Format.infinityBits(width);
    const magnitude_mask = narrow.Format.fractionMask(width) | infinity;

    // The screen rejected every invalid case, so a surviving infinity is a real
    // term of the sum and its sign settles the result.
    const c_magnitude = c_bits & magnitude_mask;
    if (c_magnitude == infinity) return narrow.fromBits(width, c_bits);
    const a_magnitude = a_bits & magnitude_mask;
    const b_magnitude = b_bits & magnitude_mask;
    if (a_magnitude == infinity or b_magnitude == infinity) {
        return narrow.fromBits(width, ((a_bits ^ b_bits) & sign_mask) | infinity);
    }

    // A zero product adds nothing to a nonzero addend, and two zeros of
    // opposite signs add up to +0 under round-to-nearest-even.
    if (a_magnitude == 0 or b_magnitude == 0) {
        if (c_magnitude != 0) return narrow.fromBits(width, c_bits);
        const product_negative = (a_bits ^ b_bits) & sign_mask != 0;
        return narrow.fromBits(width, if (product_negative and (c_bits & sign_mask != 0)) sign_mask else 0);
    }

    const product = narrow.multiplyDyadic(
        narrow.finiteDyadic(width, narrow.decode(width, a_bits)),
        narrow.finiteDyadic(width, narrow.decode(width, b_bits)),
    );
    var sum = narrow.addDyadic(
        product,
        narrow.finiteDyadic(width, narrow.decode(width, c_bits)),
    );
    // Two nonzero terms of opposite signs that cancel exactly are +0, whichever
    // of them the shared engine happened to keep.
    if (sum.magnitude == 0) sum.negative = false;
    return narrow.packDyadic(width, sum);
}

/// Fused multiply-add of three binary32 values with a single rounding. A CPU
/// with a fused instruction keeps it, and every other target computes the same
/// single rounding exactly.
pub fn fma(a: f32, b: f32, c: f32) f32 {
    @setFloatMode(.strict);
    if (comptime has_native_fma) {
        return @mulAdd(f32, a, b, c);
    }
    return softwareFma32(a, b, c);
}

/// The binary64 case of `fma`.
pub fn fma64(a: f64, b: f64, c: f64) f64 {
    @setFloatMode(.strict);
    if (comptime has_native_fma64) {
        return @mulAdd(f64, a, b, c);
    }
    return softwareFma64(a, b, c);
}

fn softwareFma32(a: f32, b: f32, c: f32) f32 {
    return @bitCast(softwareFmaBits(
        format,
        @as(u32, @bitCast(a)),
        @as(u32, @bitCast(b)),
        @as(u32, @bitCast(c)),
    ));
}

fn softwareFma64(a: f64, b: f64, c: f64) f64 {
    return @bitCast(softwareFmaBits(
        wide_format,
        @as(u64, @bitCast(a)),
        @as(u64, @bitCast(b)),
        @as(u64, @bitCast(c)),
    ));
}

/// C ABI entry point for `zmath_fmaf` in `coremath_fma.h`, which calls it on
/// every target without a fused instruction. The exported symbol always runs
/// the exact software path so that both paths stay independently checkable.
export fn zmath_core_fmaf_soft(a: f32, b: f32, c: f32) callconv(.c) f32 {
    return softwareFma32(a, b, c);
}

/// C ABI entry point for `zmath_fma` in `coremath_fma.h`, the binary64 half of
/// the same seam.
export fn zmath_core_fma_soft(a: f64, b: f64, c: f64) callconv(.c) f64 {
    return softwareFma64(a, b, c);
}
