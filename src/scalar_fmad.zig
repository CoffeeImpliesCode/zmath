// Correctly rounded fused multiply-add for binary64 values.
// Uses @mulAdd with IEEE round-to-nearest, ties-to-even.
// No host libm, no fenv mutation, no dynamic allocation.

const sign_mask_d: u64 = 0x8000_0000_0000_0000;
const magnitude_mask_d: u64 = 0x7fff_ffff_ffff_ffff;
const infinity_bits_d: u64 = 0x7ff0_0000_0000_0000;
const quiet_nan_bit_d: u64 = 0x0008_0000_0000_0000;
const canonical_nan_bits_d: u64 = 0x7ff8_0000_0000_0000;

/// Computes the correctly rounded fused multiply-add of three binary64 values.
// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn fmad(a: f64, b: f64, c: f64) f64 {
    @setFloatMode(.strict);
    const a_bits: u64 = @bitCast(a);
    const b_bits: u64 = @bitCast(b);
    const c_bits: u64 = @bitCast(c);
    const a_magnitude: u64 = a_bits & magnitude_mask_d;
    const b_magnitude: u64 = b_bits & magnitude_mask_d;
    const c_magnitude: u64 = c_bits & magnitude_mask_d;

    // Check for NaN in a, b, c (first NaN in argument order)
    if (a_magnitude > infinity_bits_d) {
        return @bitCast(a_bits | quiet_nan_bit_d);
    }
    if (b_magnitude > infinity_bits_d) {
        return @bitCast(b_bits | quiet_nan_bit_d);
    }
    if (c_magnitude > infinity_bits_d) {
        return @bitCast(c_bits | quiet_nan_bit_d);
    }
    // Check for (inf * 0) or (0 * inf) which produces a NaN in the product
    if ((a_magnitude == infinity_bits_d and b_magnitude == 0) or
        (a_magnitude == 0 and b_magnitude == infinity_bits_d))
    {
        return @bitCast(canonical_nan_bits_d);
    }

    // Check for product infinite and c infinite with opposite sign
    const a_sign: bool = (a_bits & sign_mask_d) != 0;
    const b_sign: bool = (b_bits & sign_mask_d) != 0;
    const product_sign: bool = a_sign != b_sign;
    const product_infinite: bool = (a_magnitude == infinity_bits_d and b_magnitude != 0) or
        (b_magnitude == infinity_bits_d and a_magnitude != 0);
    const c_sign: bool = (c_bits & sign_mask_d) != 0;
    if (product_infinite and
        (c_magnitude == infinity_bits_d) and
        (product_sign != c_sign))
    {
        return @bitCast(canonical_nan_bits_d);
    }

    // Otherwise, compute a*b + c with one rounding
    return @mulAdd(f64, a, b, c);
}
