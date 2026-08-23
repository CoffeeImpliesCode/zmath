const sign_mask_h: u16 = 0x8000;
const magnitude_mask_h: u16 = 0x7fff;
const infinity_bits_h: u16 = 0x7c00;
const quiet_nan_bit_h: u16 = 0x0200;
const canonical_nan_bits_h: u16 = 0x7e00;

// Correctly rounded fused multiply-add for binary16 values.
// Uses @mulAdd with IEEE round-to-nearest, ties-to-even.
// No host libm, no fenv mutation, no dynamic allocation.

/// Computes the correctly rounded fused multiply-add of three binary16 values.
// Callers must keep IEEE round-to-nearest, ties-to-even active.
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

    // Otherwise, compute a*b + c with one rounding
    return @mulAdd(f16, a, b, c);
}
