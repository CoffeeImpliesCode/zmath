// Correctly rounded complementary error function for binary16, binary32, and binary64.
//
// The kernels are derived directly from CORE-MATH at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa and are used under the MIT License.
// They are compiled freestanding and do not call host libm, fenv, or allocation APIs.

extern fn zmath_core_erfch_bits(value: u16) callconv(.c) u16;
extern fn cr_erfcf(value: f32) callconv(.c) f32;
extern fn zmath_core_erfcd_negative(value: f64) callconv(.c) f64;
extern fn zmath_core_erfcd_positive(value: f64) callconv(.c) f64;

const magnitude_mask_h: u16 = 0x7fff;
const magnitude_mask_f: u32 = 0x7fff_ffff;
const magnitude_mask_d: u64 = 0x7fff_ffff_ffff_ffff;
const infinity_bits_h: u16 = 0x7c00;
const infinity_bits_f: u32 = 0x7f80_0000;
const infinity_bits_d: u64 = 0x7ff0_0000_0000_0000;
const quiet_nan_bit_h: u16 = 0x0200;
const quiet_nan_bit_f: u32 = 0x0040_0000;
const quiet_nan_bit_d: u64 = 0x0008_0000_0000_0000;
const erfc_neg_tiny_d: u64 = 0x3c9c_5bf8_91b4_ef6a;
const erfc_pos_tiny_d: u64 = 0x3c8c_5bf8_91b4_ef6a;
const erfc_neg_saturation_d: u64 = 0x4017_744f_8f74_e94b;
const erfc_pos_saturation_d: u64 = 0x403b_39dc_41e4_8bfd;

/// Correctly rounded complementary error function for binary16 values.
pub fn erfch(value: f16) f16 {
    @setFloatMode(.strict);
    const bits: u16 = @bitCast(value);
    if ((bits & magnitude_mask_h) > infinity_bits_h) return @bitCast(bits | quiet_nan_bit_h);
    return @bitCast(zmath_core_erfch_bits(bits));
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
    const magnitude = bits & magnitude_mask_d;
    if ((bits & 0x8000_0000_0000_0000) != 0) {
        if (magnitude >= erfc_neg_saturation_d) {
            if (magnitude > infinity_bits_d) return @bitCast(bits | quiet_nan_bit_d);
            return 2.0;
        }
        if (magnitude <= erfc_neg_tiny_d) return 1.0;
        return zmath_core_erfcd_negative(value);
    }
    if (magnitude >= erfc_pos_saturation_d) {
        if (magnitude > infinity_bits_d) return @bitCast(bits | quiet_nan_bit_d);
        return 0.0;
    }
    if (magnitude <= erfc_pos_tiny_d) return 1.0;
    return zmath_core_erfcd_positive(value);
}
