// Correctly rounded complementary error function for binary16, binary32, and binary64.
//
// The kernels are derived directly from CORE-MATH at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa and are used under the MIT License.
// They are compiled freestanding and do not call host libm, fenv, or allocation APIs.

extern fn cr_erfcf16(value: f16) callconv(.c) f16;
extern fn cr_erfcf(value: f32) callconv(.c) f32;
extern fn cr_erfc(value: f64) callconv(.c) f64;

const magnitude_mask_h: u16 = 0x7fff;
const magnitude_mask_f: u32 = 0x7fff_ffff;
const magnitude_mask_d: u64 = 0x7fff_ffff_ffff_ffff;
const infinity_bits_h: u16 = 0x7c00;
const infinity_bits_f: u32 = 0x7f80_0000;
const infinity_bits_d: u64 = 0x7ff0_0000_0000_0000;
const quiet_nan_bit_h: u16 = 0x0200;
const quiet_nan_bit_f: u32 = 0x0040_0000;
const quiet_nan_bit_d: u64 = 0x0008_0000_0000_0000;

/// Correctly rounded complementary error function for binary16 values.
pub fn erfch(value: f16) f16 {
    @setFloatMode(.strict);
    const bits: u16 = @bitCast(value);
    if ((bits & magnitude_mask_h) > infinity_bits_h) return @bitCast(bits | quiet_nan_bit_h);
    return cr_erfcf16(value);
}

/// Correctly rounded complementary error function for binary32 values.
pub fn erfcf(value: f32) f32 {
    @setFloatMode(.strict);
    const bits: u32 = @bitCast(value);
    if ((bits & magnitude_mask_f) > infinity_bits_f) return @bitCast(bits | quiet_nan_bit_f);
    return cr_erfcf(value);
}

/// Correctly rounded complementary error function for binary64 values.
pub fn erfcd(value: f64) f64 {
    @setFloatMode(.strict);
    const bits: u64 = @bitCast(value);
    if ((bits & magnitude_mask_d) > infinity_bits_d) return @bitCast(bits | quiet_nan_bit_d);
    return cr_erfc(value);
}
