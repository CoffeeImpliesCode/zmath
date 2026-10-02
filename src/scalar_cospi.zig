// Correctly rounded cosine of pi times the input for binary16, binary32, and binary64.
//
// The finite kernels are derived directly from CORE-MATH at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa and are used under the MIT License.
// Exceptional results are constructed from target-format bits so NaN payloads,
// signs, and infinity behavior do not depend on host library conventions.

extern fn zmath_core_cospih_bits(value: u16) callconv(.c) u16;
extern fn cr_cospif(value: f32) callconv(.c) f32;
extern fn cr_cospi(value: f64) callconv(.c) f64;

const magnitude_mask_h: u16 = 0x7fff;
const magnitude_mask_f: u32 = 0x7fff_ffff;
const magnitude_mask_d: u64 = 0x7fff_ffff_ffff_ffff;
const infinity_bits_h: u16 = 0x7c00;
const infinity_bits_f: u32 = 0x7f80_0000;
const infinity_bits_d: u64 = 0x7ff0_0000_0000_0000;
const quiet_nan_bit_h: u16 = 0x0200;
const quiet_nan_bit_f: u32 = 0x0040_0000;
const quiet_nan_bit_d: u64 = 0x0008_0000_0000_0000;
const canonical_nan_bits_h: u16 = 0x7e00;
const canonical_nan_bits_f: u32 = 0x7fc0_0000;
const canonical_nan_bits_d: u64 = 0x7ff8_0000_0000_0000;

/// Correctly rounded cos(pi * value) for a binary16 value.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn cospih(value: f16) f16 {
    @setFloatMode(.strict);
    const bits: u16 = @bitCast(value);
    const magnitude: u16 = bits & magnitude_mask_h;
    if (magnitude > infinity_bits_h) return @bitCast(bits | quiet_nan_bit_h);
    if (magnitude == infinity_bits_h) return @bitCast(canonical_nan_bits_h);
    return @bitCast(zmath_core_cospih_bits(bits));
}

/// Correctly rounded cos(pi * value) for a binary32 value.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn cospif(value: f32) f32 {
    @setFloatMode(.strict);
    const bits: u32 = @bitCast(value);
    const magnitude: u32 = bits & magnitude_mask_f;
    if (magnitude > infinity_bits_f) return @bitCast(bits | quiet_nan_bit_f);
    if (magnitude == infinity_bits_f) return @bitCast(canonical_nan_bits_f);
    return cr_cospif(value);
}

/// Correctly rounded cos(pi * value) for a binary64 value.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn cospid(value: f64) f64 {
    @setFloatMode(.strict);
    const bits: u64 = @bitCast(value);
    const magnitude: u64 = bits & magnitude_mask_d;
    if (magnitude > infinity_bits_d) return @bitCast(bits | quiet_nan_bit_d);
    if (magnitude == infinity_bits_d) return @bitCast(canonical_nan_bits_d);
    return cr_cospi(value);
}
