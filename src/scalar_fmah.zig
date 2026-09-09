const std = @import("std");
const narrow = @import("scalar_narrow_common.zig");

const sign_mask_h: u16 = 0x8000;
const magnitude_mask_h: u16 = 0x7fff;
const infinity_bits_h: u16 = 0x7c00;
const quiet_nan_bit_h: u16 = 0x0200;
const canonical_nan_bits_h: u16 = 0x7e00;

// Correctly rounded fused multiply-add for binary16 values.
// Finite inputs use exact integer dyadics and one final round-to-nearest,
// ties-to-even packing step. No host libm, fenv mutation, or allocation.

/// Computes the correctly rounded fused multiply-add of three binary16 values.
pub fn fmah(a: f16, b: f16, c: f16) f16 {
    @setFloatMode(.strict);
    const a_bits: u16 = @bitCast(a);
    const b_bits: u16 = @bitCast(b);
    const c_bits: u16 = @bitCast(c);
    const a_magnitude: u16 = a_bits & magnitude_mask_h;
    const b_magnitude: u16 = b_bits & magnitude_mask_h;
    const c_magnitude: u16 = c_bits & magnitude_mask_h;

    // Check for NaN in a, b, c (first NaN in argument order)
    if (a_magnitude > infinity_bits_h) {
        return @bitCast(a_bits | quiet_nan_bit_h);
    }
    if (b_magnitude > infinity_bits_h) {
        return @bitCast(b_bits | quiet_nan_bit_h);
    }
    if (c_magnitude > infinity_bits_h) {
        return @bitCast(c_bits | quiet_nan_bit_h);
    }

    // Check for (inf * 0) or (0 * inf) which produces a NaN in the product
    if ((a_magnitude == infinity_bits_h and b_magnitude == 0) or
        (a_magnitude == 0 and b_magnitude == infinity_bits_h))
    {
        return @bitCast(canonical_nan_bits_h);
    }

    // Check for product infinite and c infinite with opposite sign
    const a_sign: bool = (a_bits & sign_mask_h) != 0;
    const b_sign: bool = (b_bits & sign_mask_h) != 0;
    const product_sign: bool = a_sign != b_sign;
    const product_infinite: bool = (a_magnitude == infinity_bits_h and b_magnitude != 0) or
        (b_magnitude == infinity_bits_h and a_magnitude != 0);
    const c_sign: bool = (c_bits & sign_mask_h) != 0;
    if (product_infinite and
        (c_magnitude == infinity_bits_h) and
        (product_sign != c_sign))
    {
        return @bitCast(canonical_nan_bits_h);
    }

    if (product_infinite) {
        const result_bits = infinity_bits_h | (if (product_sign) sign_mask_h else 0);
        return @bitCast(result_bits);
    }
    if (c_magnitude == infinity_bits_h) return c;

    // Each finite operand is M * 2^E with M <= 2^11 - 1 and -24 <= E <= 5.
    // The product magnitude is below 2^22 with -48 <= E <= 10. When the
    // terms are aligned at the smaller exponent, the largest shift is
    // 5 - (-48) = 53, so the larger addend occupies at most 64 bits and
    // their exact sum at most 65 bits. addDyadic cannot discard a bit here;
    // packDyadic performs the only RNE rounding.
    const first = narrow.finiteDyadic(narrow.H, narrow.decode(narrow.H, a_bits));
    const second = narrow.finiteDyadic(narrow.H, narrow.decode(narrow.H, b_bits));
    const addend = narrow.finiteDyadic(narrow.H, narrow.decode(narrow.H, c_bits));
    const product = narrow.multiplyDyadic(first, second);
    var result = narrow.addDyadic(product, addend);
    if (result.magnitude == 0) {
        // Exact cancellation produces +0. Only two negative zero terms produce -0.
        result.negative = product_sign and c_sign;
    }
    return @bitCast(narrow.packDyadic(narrow.H, result));
}

test "binary16 fused multiply-add avoids double rounding at the normal boundary" {
    var inputs = [_]f16{
        @bitCast(@as(u16, 0x03ff)),
        @bitCast(@as(u16, 0x1001)),
    };
    std.mem.doNotOptimizeAway(&inputs);
    const result: u16 = @bitCast(fmah(inputs[0], inputs[1], inputs[0]));
    try std.testing.expectEqual(@as(u16, 0x03ff), result);
}

test "binary16 fused multiply-add rounds an exact underflow tie to even" {
    var inputs = [_]f16{
        @bitCast(@as(u16, 0x0001)),
        @bitCast(@as(u16, 0x3800)),
        @bitCast(@as(u16, 0x0000)),
    };
    std.mem.doNotOptimizeAway(&inputs);
    // 2^-24 * 2^-1 is halfway between zero and the minimum subnormal; zero is even.
    const result: u16 = @bitCast(fmah(inputs[0], inputs[1], inputs[2]));
    try std.testing.expectEqual(@as(u16, 0x0000), result);
}

test "binary16 fused multiply-add rounds an exact overflow tie to infinity" {
    var inputs = [_]f16{
        @bitCast(@as(u16, 0x7bff)),
        @bitCast(@as(u16, 0x3c00)),
        @bitCast(@as(u16, 0x4c00)),
    };
    std.mem.doNotOptimizeAway(&inputs);
    // 65504 + 16 is halfway to the next binade; the even carry overflows to infinity.
    const result: u16 = @bitCast(fmah(inputs[0], inputs[1], inputs[2]));
    try std.testing.expectEqual(@as(u16, 0x7c00), result);
}

test "binary16 fused multiply-add preserves exact zero signs" {
    var cancellation = [_]f16{
        @bitCast(@as(u16, 0x3c00)),
        @bitCast(@as(u16, 0x3c00)),
        @bitCast(@as(u16, 0xbc00)),
    };
    var negative_zeros = [_]f16{
        @bitCast(@as(u16, 0x8000)),
        @bitCast(@as(u16, 0x3c00)),
        @bitCast(@as(u16, 0x8000)),
    };
    std.mem.doNotOptimizeAway(&cancellation);
    std.mem.doNotOptimizeAway(&negative_zeros);
    const cancellation_result: u16 = @bitCast(fmah(
        cancellation[0],
        cancellation[1],
        cancellation[2],
    ));
    const negative_zero_result: u16 = @bitCast(fmah(
        negative_zeros[0],
        negative_zeros[1],
        negative_zeros[2],
    ));
    // Opposite signed zero terms cancel to +0; two negative zero terms produce -0.
    try std.testing.expectEqual(@as(u16, 0x0000), cancellation_result);
    try std.testing.expectEqual(@as(u16, 0x8000), negative_zero_result);
}
