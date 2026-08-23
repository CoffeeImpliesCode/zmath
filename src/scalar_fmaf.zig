const sign_mask_f: u32 = 0x8000_0000;
const magnitude_mask_f: u32 = 0x7fff_ffff;
const infinity_bits_f: u32 = 0x7f80_0000;
const quiet_nan_bit_f: u32 = 0x0040_0000;
const canonical_nan_bits_f: u32 = 0x7fc0_0000;

// Correctly rounded fused multiply-add for binary32 values.
// Uses @mulAdd with IEEE round-to-nearest, ties-to-even.
// No host libm, no fenv mutation, no dynamic allocation.

/// Computes the correctly rounded fused multiply-add of three binary32 values.
// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn fmaf(a: f32, b: f32, c: f32) f32 {
    @setFloatMode(.strict);
    const a_bits: u32 = @bitCast(a);
    const b_bits: u32 = @bitCast(b);
    const c_bits: u32 = @bitCast(c);
    const a_magnitude: u32 = a_bits & magnitude_mask_f;
    const b_magnitude: u32 = b_bits & magnitude_mask_f;
    const c_magnitude: u32 = c_bits & magnitude_mask_f;

    // Check for NaN in a, b, c (first NaN in argument order)
    if (a_magnitude > infinity_bits_f) {
        return @bitCast(a_bits | quiet_nan_bit_f);
    }
    if (b_magnitude > infinity_bits_f) {
        return @bitCast(b_bits | quiet_nan_bit_f);
    }
    if (c_magnitude > infinity_bits_f) {
        return @bitCast(c_bits | quiet_nan_bit_f);
    }

    // Check for (inf * 0) or (0 * inf) which produces a NaN in the product
    if ((a_magnitude == infinity_bits_f and b_magnitude == 0) or
        (a_magnitude == 0 and b_magnitude == infinity_bits_f))
    {
        return @bitCast(canonical_nan_bits_f);
    }

    // Check for product infinite and c infinite with opposite sign
    const a_sign: bool = (a_bits & sign_mask_f) != 0;
    const b_sign: bool = (b_bits & sign_mask_f) != 0;
    const product_sign: bool = a_sign != b_sign;
    const product_infinite: bool = (a_magnitude == infinity_bits_f and b_magnitude != 0) or
        (b_magnitude == infinity_bits_f and a_magnitude != 0);
    const c_sign: bool = (c_bits & sign_mask_f) != 0;
    if (product_infinite and
        (c_magnitude == infinity_bits_f) and
        (product_sign != c_sign))
    {
        return @bitCast(canonical_nan_bits_f);
    }

    // Otherwise, compute a*b + c with one rounding
    return @mulAdd(f32, a, b, c);
}
