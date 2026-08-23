const sign_mask_h: u16 = 0x8000;
const magnitude_mask_h: u16 = 0x7fff;
const infinity_bits_h: u16 = 0x7c00;
const quiet_nan_bit_h: u16 = 0x0200;
const canonical_nan_bits_h: u16 = 0x7e00;

const sign_mask_f: u32 = 0x8000_0000;
const magnitude_mask_f: u32 = 0x7fff_ffff;
const infinity_bits_f: u32 = 0x7f80_0000;
const quiet_nan_bit_f: u32 = 0x0040_0000;
const canonical_nan_bits_f: u32 = 0x7fc0_0000;

const sign_mask_d: u64 = 0x8000_0000_0000_0000;
const magnitude_mask_d: u64 = 0x7fff_ffff_ffff_ffff;
const infinity_bits_d: u64 = 0x7ff0_0000_0000_0000;
const quiet_nan_bit_d: u64 = 0x0008_0000_0000_0000;
const canonical_nan_bits_d: u64 = 0x7ff8_0000_0000_0000;

// Correctly rounded square root for binary16 values.
//
// Uses @sqrt with IEEE round-to-nearest, ties-to-even.
// No host libm, no fenv mutation, no dynamic allocation.

/// Computes the correctly rounded square root of a binary16 value.
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn sqrth(value: f16) f16 {
    @setFloatMode(.strict);
    const bits: u16 = @bitCast(value);
    const magnitude = bits & magnitude_mask_h;

    if (magnitude == 0) {
        // +/-0 preserves sign
        return value;
    }
    if (magnitude > infinity_bits_h) {
        // NaN: quiet and preserve payload
        return @bitCast(bits | quiet_nan_bit_h);
    }
    if ((bits & sign_mask_h) != 0) {
        // negative finite or -Inf
        return @bitCast(canonical_nan_bits_h);
    }
    // positive +Inf or normal positive finite: @sqrt is exact for these
    return @sqrt(value);
}

/// Computes the correctly rounded square root of a binary32 value.
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn sqrtf(value: f32) f32 {
    @setFloatMode(.strict);
    const bits: u32 = @bitCast(value);
    const magnitude = bits & magnitude_mask_f;

    if (magnitude == 0) {
        // +/-0 preserves sign
        return value;
    }
    if (magnitude > infinity_bits_f) {
        // NaN: quiet and preserve payload
        return @bitCast(bits | quiet_nan_bit_f);
    }
    if ((bits & sign_mask_f) != 0) {
        // negative finite or -Inf
        return @bitCast(canonical_nan_bits_f);
    }
    // positive +Inf or normal positive finite: @sqrt is exact for these
    return @sqrt(value);
}

/// Computes the correctly rounded square root of a binary64 value.
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn sqrtd(value: f64) f64 {
    @setFloatMode(.strict);
    const bits: u64 = @bitCast(value);
    const magnitude = bits & magnitude_mask_d;

    if (magnitude == 0) {
        // +/-0 preserves sign
        return value;
    }
    if (magnitude > infinity_bits_d) {
        // NaN: quiet and preserve payload
        return @bitCast(bits | quiet_nan_bit_d);
    }
    if ((bits & sign_mask_d) != 0) {
        // negative finite or -Inf
        return @bitCast(canonical_nan_bits_d);
    }
    // positive +Inf or normal positive finite: @sqrt is exact for these
    return @sqrt(value);
}
