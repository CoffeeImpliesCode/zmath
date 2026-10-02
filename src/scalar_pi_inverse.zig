// C23 inverse trigonometric functions measured in half revolutions.
//
// Exceptional values are handled here because the direct kernels expose
// operation-specific special cases. Finite interior values are delegated to
// the pinned, correctly-rounded CORE-MATH kernels without host math calls.

const common = @import("scalar_c23_common.zig");
const c23 = @import("scalar_c23_types.zig");
const tiny = @import("scalar_tiny.zig");

extern fn zmath_core_asinpih_bits(value: u16) callconv(.c) u16;
extern fn zmath_core_asinpif(value: f32) callconv(.c) f32;
extern fn zmath_core_asinpid_finite_bits(value: f64, bits: u64) callconv(.c) f64;
extern fn zmath_core_acospih_bits(value: u16) callconv(.c) u16;
extern fn zmath_core_acospif(value: f32) callconv(.c) f32;
extern fn zmath_core_acospid(value: f64) callconv(.c) f64;
extern fn zmath_core_atanpih_bits(value: u16) callconv(.c) u16;
extern fn zmath_core_atanpif(value: f32) callconv(.c) f32;
extern fn zmath_core_atanpid(value: f64) callconv(.c) f64;
extern fn zmath_core_atan2pih_bits(y: u16, x: u16) callconv(.c) u16;
extern fn zmath_core_atan2pif(y: f32, x: f32) callconv(.c) f32;
extern fn zmath_core_atan2pid_finite_bits(
    y: f64,
    x: f64,
    y_bits: u64,
    x_bits: u64,
) callconv(.c) f64;
extern fn zmath_core_atan2pid_forced_ambiguity(
    y: f64,
    x: f64,
    value: *f64,
) callconv(.c) c_int;

const sign_mask_h: u16 = 0x8000;
const magnitude_mask_h: u16 = 0x7fff;
const one_bits_h: u16 = 0x3c00;
const quarter_bits_h: u16 = 0x3400;
const half_bits_h: u16 = 0x3800;
const three_quarter_bits_h: u16 = 0x3a00;
const one_result_bits_h: u16 = 0x3c00;

const sign_mask_f: u32 = 0x8000_0000;
const magnitude_mask_f: u32 = 0x7fff_ffff;
const one_bits_f: u32 = 0x3f80_0000;
const quarter_bits_f: u32 = 0x3e80_0000;
const half_bits_f: u32 = 0x3f00_0000;
const three_quarter_bits_f: u32 = 0x3f40_0000;
const one_result_bits_f: u32 = 0x3f80_0000;

const sign_mask_d: u64 = 0x8000_0000_0000_0000;
const magnitude_mask_d: u64 = 0x7fff_ffff_ffff_ffff;
const infinity_bits_d: u64 = 0x7ff0_0000_0000_0000;
const one_bits_d: u64 = 0x3ff0_0000_0000_0000;
const quarter_bits_d: u64 = 0x3fd0_0000_0000_0000;
const half_bits_d: u64 = 0x3fe0_0000_0000_0000;
const three_quarter_bits_d: u64 = 0x3fe8_0000_0000_0000;
const one_result_bits_d: u64 = 0x3ff0_0000_0000_0000;

fn signedH(magnitude_bits: u16, sign_source_bits: u16) f16 {
    return @bitCast(magnitude_bits | (sign_source_bits & sign_mask_h));
}

fn signedF(magnitude_bits: u32, sign_source_bits: u32) f32 {
    return @bitCast(magnitude_bits | (sign_source_bits & sign_mask_f));
}

fn signedD(magnitude_bits: u64, sign_source_bits: u64) f64 {
    return @bitCast(magnitude_bits | (sign_source_bits & sign_mask_d));
}

fn atan2ResolvedH(value: f16) c23.Atan2Result(f16) {
    return .{ .value = value, .resolved = true };
}

fn atan2ResolvedF(value: f32) c23.Atan2Result(f32) {
    return .{ .value = value, .resolved = true };
}

fn atan2ResolvedD(value: f64) c23.Atan2Result(f64) {
    return .{ .value = value, .resolved = true };
}

fn atan2UnresolvedD() c23.Atan2Result(f64) {
    return .{
        .value = @bitCast(@as(u64, 0x7ff8_0000_0000_0000)),
        .resolved = false,
    };
}

/// Return asin(x) / pi for binary16.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn asinpih(value: f16) f16 {
    @setFloatMode(.strict);
    const bits = common.H.toBits(value);
    const magnitude = bits & magnitude_mask_h;
    if (common.H.isNaNBits(bits)) return common.H.fromBits(common.H.quietNaNBits(bits));
    if (common.H.isInfinityBits(bits)) return common.H.canonicalNaN();
    if (magnitude > one_bits_h) return common.H.canonicalNaN();
    if (common.H.isZeroBits(bits)) return value;
    if (magnitude == one_bits_h) return signedH(half_bits_h, bits);
    return @bitCast(zmath_core_asinpih_bits(bits));
}

/// Return asin(x) / pi for binary32.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn asinpif(value: f32) f32 {
    @setFloatMode(.strict);
    const bits = common.F.toBits(value);
    const magnitude = bits & magnitude_mask_f;
    if (common.F.isNaNBits(bits)) return common.F.fromBits(common.F.quietNaNBits(bits));
    if (magnitude <= one_bits_f) return zmath_core_asinpif(value);
    return common.F.canonicalNaN();
}

/// Return asin(x) / pi for binary64.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn asinpid(value: f64) f64 {
    @setFloatMode(.strict);
    const bits = common.D.toBits(value);
    const magnitude = bits & magnitude_mask_d;
    if (magnitude >= one_bits_d) {
        if (magnitude == one_bits_d) return signedD(half_bits_d, bits);
        if (common.D.isNaNBits(bits)) {
            return common.D.fromBits(common.D.quietNaNBits(bits));
        }
        return common.D.canonicalNaN();
    }
    if (magnitude < tiny.magnitude_limit) {
        // For this range, |asin(x) - x| / pi < |x|^3 < |x| * 2^-2042.
        // This is smaller than the extra |x| * 2^-128 remainder accepted by scale().
        return tiny.scale(.inv_pi, value);
    }
    return zmath_core_asinpid_finite_bits(value, bits);
}

/// Return acos(x) / pi for binary16.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn acospih(value: f16) f16 {
    @setFloatMode(.strict);
    const bits = common.H.toBits(value);
    const magnitude = bits & magnitude_mask_h;
    if (common.H.isNaNBits(bits)) return common.H.fromBits(common.H.quietNaNBits(bits));
    if (common.H.isInfinityBits(bits)) return common.H.canonicalNaN();
    if (magnitude > one_bits_h) return common.H.canonicalNaN();
    if (common.H.isZeroBits(bits)) return @bitCast(half_bits_h);
    if (magnitude == one_bits_h) {
        if ((bits & sign_mask_h) == 0) return @bitCast(@as(u16, 0x0000));
        return @bitCast(one_result_bits_h);
    }
    return @bitCast(zmath_core_acospih_bits(bits));
}

/// Return acos(x) / pi for binary32.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn acospif(value: f32) f32 {
    @setFloatMode(.strict);
    const bits = common.F.toBits(value);
    const magnitude = bits & magnitude_mask_f;
    if (common.F.isNaNBits(bits)) return common.F.fromBits(common.F.quietNaNBits(bits));
    if (common.F.isInfinityBits(bits)) return common.F.canonicalNaN();
    if (magnitude > one_bits_f) return common.F.canonicalNaN();
    if (common.F.isZeroBits(bits)) return @bitCast(half_bits_f);
    if (magnitude == one_bits_f) {
        if ((bits & sign_mask_f) == 0) return @bitCast(@as(u32, 0x0000_0000));
        return @bitCast(one_result_bits_f);
    }
    return zmath_core_acospif(value);
}

/// Return acos(x) / pi for binary64.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn acospid(value: f64) f64 {
    @setFloatMode(.strict);
    const bits = common.D.toBits(value);
    const magnitude = bits & magnitude_mask_d;
    // For |x| < 0x1.921fbp-54, acos(x)/pi rounds to 0.5, including both zeros.
    if (magnitude < 0x3c99_21fb_0000_0000) return @bitCast(half_bits_d);
    if (common.D.isNaNBits(bits)) return common.D.fromBits(common.D.quietNaNBits(bits));
    if (common.D.isInfinityBits(bits)) return common.D.canonicalNaN();
    if (magnitude > one_bits_d) return common.D.canonicalNaN();
    if (magnitude == one_bits_d) {
        if ((bits & sign_mask_d) == 0) return @bitCast(@as(u64, 0x0000_0000_0000_0000));
        return @bitCast(one_result_bits_d);
    }
    return zmath_core_acospid(value);
}

/// Return atan(x) / pi for binary16.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn atanpih(value: f16) f16 {
    @setFloatMode(.strict);
    const bits = common.H.toBits(value);
    const magnitude = bits & magnitude_mask_h;
    if (common.H.isNaNBits(bits)) return common.H.fromBits(common.H.quietNaNBits(bits));
    if (common.H.isZeroBits(bits)) return value;
    if (common.H.isInfinityBits(bits)) return signedH(half_bits_h, bits);
    if (magnitude == one_bits_h) return signedH(quarter_bits_h, bits);
    return @bitCast(zmath_core_atanpih_bits(bits));
}

/// Return atan(x) / pi for binary32.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn atanpif(value: f32) f32 {
    @setFloatMode(.strict);
    const bits = common.F.toBits(value);
    if (common.F.isNaNBits(bits)) return common.F.fromBits(common.F.quietNaNBits(bits));
    return zmath_core_atanpif(value);
}

/// Return atan(x) / pi for binary64.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn atanpid(value: f64) f64 {
    @setFloatMode(.strict);
    const bits = common.D.toBits(value);
    const magnitude = bits & magnitude_mask_d;
    if (common.D.isNaNBits(bits)) return common.D.fromBits(common.D.quietNaNBits(bits));
    if (magnitude < tiny.magnitude_limit) {
        // For this range, |atan(x) - x| / pi < |x|^3 / 3 < |x| * 2^-2043.
        // This is smaller than the extra |x| * 2^-128 remainder accepted by scale().
        return tiny.scale(.inv_pi, value);
    }
    return zmath_core_atanpid(value);
}

/// Return atan2(y, x) / pi for binary16.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn atan2pih(y: f16, x: f16) c23.Atan2Result(f16) {
    @setFloatMode(.strict);
    if (common.H.firstNaN2(y, x)) |nan| return atan2ResolvedH(nan);

    const y_bits = common.H.toBits(y);
    const x_bits = common.H.toBits(x);
    const y_magnitude = y_bits & magnitude_mask_h;
    const x_magnitude = x_bits & magnitude_mask_h;

    if (common.H.isInfinityBits(y_bits)) {
        if (common.H.isInfinityBits(x_bits)) {
            const magnitude = if ((x_bits & sign_mask_h) == 0)
                quarter_bits_h
            else
                three_quarter_bits_h;
            return atan2ResolvedH(signedH(magnitude, y_bits));
        }
        return atan2ResolvedH(signedH(half_bits_h, y_bits));
    }
    if (common.H.isInfinityBits(x_bits)) {
        if ((x_bits & sign_mask_h) == 0) return atan2ResolvedH(signedH(0, y_bits));
        return atan2ResolvedH(signedH(one_result_bits_h, y_bits));
    }
    if (common.H.isZeroBits(y_bits)) {
        if (common.H.isZeroBits(x_bits) or (x_bits & sign_mask_h) != 0) {
            const magnitude = if ((x_bits & sign_mask_h) == 0) 0 else one_result_bits_h;
            return atan2ResolvedH(signedH(magnitude, y_bits));
        }
        return atan2ResolvedH(y);
    }
    if (common.H.isZeroBits(x_bits)) return atan2ResolvedH(signedH(half_bits_h, y_bits));
    if (y_magnitude == x_magnitude) {
        const magnitude = if ((x_bits & sign_mask_h) == 0)
            quarter_bits_h
        else
            three_quarter_bits_h;
        return atan2ResolvedH(signedH(magnitude, y_bits));
    }
    return atan2ResolvedH(@bitCast(zmath_core_atan2pih_bits(y_bits, x_bits)));
}

/// Return atan2(y, x) / pi for binary32.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn atan2pif(y: f32, x: f32) c23.Atan2Result(f32) {
    @setFloatMode(.strict);
    if (common.F.firstNaN2(y, x)) |nan| return atan2ResolvedF(nan);

    const y_bits = common.F.toBits(y);
    const x_bits = common.F.toBits(x);
    const y_magnitude = y_bits & magnitude_mask_f;
    const x_magnitude = x_bits & magnitude_mask_f;

    if (common.F.isInfinityBits(y_bits)) {
        if (common.F.isInfinityBits(x_bits)) {
            const magnitude = if ((x_bits & sign_mask_f) == 0)
                quarter_bits_f
            else
                three_quarter_bits_f;
            return atan2ResolvedF(signedF(magnitude, y_bits));
        }
        return atan2ResolvedF(signedF(half_bits_f, y_bits));
    }
    if (common.F.isInfinityBits(x_bits)) {
        if ((x_bits & sign_mask_f) == 0) return atan2ResolvedF(signedF(0, y_bits));
        return atan2ResolvedF(signedF(one_result_bits_f, y_bits));
    }
    if (common.F.isZeroBits(y_bits)) {
        if (common.F.isZeroBits(x_bits) or (x_bits & sign_mask_f) != 0) {
            const magnitude = if ((x_bits & sign_mask_f) == 0) 0 else one_result_bits_f;
            return atan2ResolvedF(signedF(magnitude, y_bits));
        }
        return atan2ResolvedF(y);
    }
    if (common.F.isZeroBits(x_bits)) return atan2ResolvedF(signedF(half_bits_f, y_bits));
    if (y_magnitude == x_magnitude) {
        const magnitude = if ((x_bits & sign_mask_f) == 0)
            quarter_bits_f
        else
            three_quarter_bits_f;
        return atan2ResolvedF(signedF(magnitude, y_bits));
    }
    return atan2ResolvedF(zmath_core_atan2pif(y, x));
}

/// Return atan2(y, x) / pi for binary64.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn atan2pid(y: f64, x: f64) c23.Atan2Result(f64) {
    @setFloatMode(.strict);
    const y_bits = common.D.toBits(y);
    const x_bits = common.D.toBits(x);
    const y_magnitude = y_bits & magnitude_mask_d;
    const x_magnitude = x_bits & magnitude_mask_d;

    if (y_magnitude >= infinity_bits_d or x_magnitude >= infinity_bits_d) {
        if (y_magnitude > infinity_bits_d) {
            return atan2ResolvedD(common.D.fromBits(common.D.quietNaNBits(y_bits)));
        }
        if (x_magnitude > infinity_bits_d) {
            return atan2ResolvedD(common.D.fromBits(common.D.quietNaNBits(x_bits)));
        }
        if (y_magnitude == infinity_bits_d) {
            if (x_magnitude == infinity_bits_d) {
                const magnitude = if ((x_bits & sign_mask_d) == 0)
                    quarter_bits_d
                else
                    three_quarter_bits_d;
                return atan2ResolvedD(signedD(magnitude, y_bits));
            }
            return atan2ResolvedD(signedD(half_bits_d, y_bits));
        }
        if ((x_bits & sign_mask_d) == 0) return atan2ResolvedD(signedD(0, y_bits));
        return atan2ResolvedD(signedD(one_result_bits_d, y_bits));
    }
    if (y_magnitude == 0 or x_magnitude == 0) {
        if (y_magnitude == 0) {
            const magnitude: u64 = if ((x_bits & sign_mask_d) == 0)
                0
            else
                one_result_bits_d;
            return atan2ResolvedD(signedD(magnitude, y_bits));
        }
        return atan2ResolvedD(signedD(half_bits_d, y_bits));
    }
    if (y_magnitude == x_magnitude) {
        const magnitude = if ((x_bits & sign_mask_d) == 0)
            quarter_bits_d
        else
            three_quarter_bits_d;
        return atan2ResolvedD(signedD(magnitude, y_bits));
    }

    // The finite C ABI reserves canonical qNaN for an unresolved result.
    const value = zmath_core_atan2pid_finite_bits(y, x, y_bits, x_bits);
    return .{
        .value = value,
        .resolved = !common.D.isNaNBits(common.D.toBits(value)),
    };
}

pub const testing = struct {
    pub fn forcedAmbiguity(y: f64, x: f64) c23.Atan2Result(f64) {
        // SAFETY: The C call writes `value` before returning nonzero; zero exits before the read.
        var value: f64 = undefined;
        if (zmath_core_atan2pid_forced_ambiguity(y, x, &value) == 0) return atan2UnresolvedD();
        return atan2ResolvedD(value);
    }
};
