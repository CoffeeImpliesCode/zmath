// Correctly rounded atan2 for binary16, binary32, and binary64.
//
// The finite kernels are ported from CORE-MATH revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa:
//   src/binary16/atan2/atan2f16.c (Copyright (c) 2022-2025 Alexei Sibidanov and Paul Zimmermann).
//   src/binary32/atan2/atan2f.c (Copyright (c) 2022-2025 Alexei Sibidanov and Paul Zimmermann).
//   src/binary64/atan2/atan2.c (Copyright (c) 2024-2025 Paul Zimmermann and Alexei Sibidanov).
//
// The wrappers handle C23 special-value precedence and exact quadrant constants
// from target-format bits. Finite, nonzero arguments go to the direct kernels.

const c23 = @import("scalar_c23_types.zig");

extern fn cr_atan2f16(y: f16, x: f16) callconv(.c) f16;
extern fn cr_atan2f(y: f32, x: f32) callconv(.c) f32;
extern fn zmath_core_atan2d_checked(y: f64, x: f64, out: *f64) callconv(.c) c_int;
extern fn zmath_core_atan2d_forced_ambiguity(y: f64, x: f64, out: *f64) callconv(.c) c_int;

const sign_mask_h: u16 = 0x8000;
const magnitude_mask_h: u16 = 0x7fff;
const infinity_bits_h: u16 = 0x7c00;
const quiet_nan_bit_h: u16 = 0x0200;
const pi_bits_h: u16 = 0x4248;
const pi_over_two_bits_h: u16 = 0x3e48;
const pi_over_four_bits_h: u16 = 0x3a48;
const three_pi_over_four_bits_h: u16 = 0x40b6;

const sign_mask_f: u32 = 0x8000_0000;
const magnitude_mask_f: u32 = 0x7fff_ffff;
const infinity_bits_f: u32 = 0x7f80_0000;
const quiet_nan_bit_f: u32 = 0x0040_0000;
const pi_bits_f: u32 = 0x4049_0fdb;
const pi_over_two_bits_f: u32 = 0x3fc9_0fdb;
const pi_over_four_bits_f: u32 = 0x3f49_0fdb;
const three_pi_over_four_bits_f: u32 = 0x4016_cbe4;

const sign_mask_d: u64 = 0x8000_0000_0000_0000;
const magnitude_mask_d: u64 = 0x7fff_ffff_ffff_ffff;
const infinity_bits_d: u64 = 0x7ff0_0000_0000_0000;
const quiet_nan_bit_d: u64 = 0x0008_0000_0000_0000;
const pi_bits_d: u64 = 0x4009_21fb_5444_2d18;
const pi_over_two_bits_d: u64 = 0x3ff9_21fb_5444_2d18;
const pi_over_four_bits_d: u64 = 0x3fe9_21fb_5444_2d18;
const three_pi_over_four_bits_d: u64 = 0x4002_d97c_7f33_21d2;

fn signedH(magnitude_bits: u16, sign_bits: u16) f16 {
    return @bitCast(magnitude_bits | (sign_bits & sign_mask_h));
}

fn signedF(magnitude_bits: u32, sign_bits: u32) f32 {
    return @bitCast(magnitude_bits | (sign_bits & sign_mask_f));
}

fn signedD(magnitude_bits: u64, sign_bits: u64) f64 {
    return @bitCast(magnitude_bits | (sign_bits & sign_mask_d));
}

fn unresolvedD() c23.Atan2Result(f64) {
    return .{
        .value = @bitCast(@as(u64, 0x7ff8_0000_0000_0000)),
        .resolved = false,
    };
}

/// Test-only access to the checked ambiguity hook.
pub const testing = struct {
    pub fn forcedAmbiguity(y: f64, x: f64) c23.Atan2Result(f64) {
        // SAFETY: The C call writes `value` before returning nonzero; zero exits before the read.
        var value: f64 = undefined;
        if (zmath_core_atan2d_forced_ambiguity(y, x, &value) == 0)
            return unresolvedD();
        return .{ .value = value, .resolved = true };
    }
};

/// Return an `Atan2Result(f16)` for correctly rounded atan2(y, x).
///
/// `value` contains the result and `resolved` is true for every binary16
/// result, including the special-value branches.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn atan2h(y: f16, x: f16) c23.Atan2Result(f16) {
    @setFloatMode(.strict);

    const y_bits: u16 = @bitCast(y);
    const x_bits: u16 = @bitCast(x);
    const y_magnitude = y_bits & magnitude_mask_h;
    const x_magnitude = x_bits & magnitude_mask_h;

    if (y_magnitude > infinity_bits_h) return .{ .value = @bitCast(y_bits | quiet_nan_bit_h), .resolved = true };
    if (x_magnitude > infinity_bits_h) return .{ .value = @bitCast(x_bits | quiet_nan_bit_h), .resolved = true };

    if (y_magnitude == infinity_bits_h) {
        if (x_magnitude == infinity_bits_h) {
            return .{
                .value = if ((x_bits & sign_mask_h) != 0)
                    signedH(three_pi_over_four_bits_h, y_bits)
                else
                    signedH(pi_over_four_bits_h, y_bits),
                .resolved = true,
            };
        }
        return .{ .value = signedH(pi_over_two_bits_h, y_bits), .resolved = true };
    }

    if (x_magnitude == infinity_bits_h) {
        return .{
            .value = if ((x_bits & sign_mask_h) != 0)
                signedH(pi_bits_h, y_bits)
            else
                @bitCast(y_bits & sign_mask_h),
            .resolved = true,
        };
    }

    if (y_magnitude == 0) {
        if (x_magnitude == 0) {
            return .{
                .value = if ((x_bits & sign_mask_h) != 0)
                    signedH(pi_bits_h, y_bits)
                else
                    y,
                .resolved = true,
            };
        }
        return .{
            .value = if ((x_bits & sign_mask_h) != 0)
                signedH(pi_bits_h, y_bits)
            else
                y,
            .resolved = true,
        };
    }

    if (x_magnitude == 0) return .{ .value = signedH(pi_over_two_bits_h, y_bits), .resolved = true };
    return .{ .value = cr_atan2f16(y, x), .resolved = true };
}

/// Return an `Atan2Result(f32)` for correctly rounded atan2(y, x).
///
/// `value` contains the result and `resolved` is true for every binary32
/// result, including the special-value branches.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn atan2f(y: f32, x: f32) c23.Atan2Result(f32) {
    @setFloatMode(.strict);

    const y_bits: u32 = @bitCast(y);
    const x_bits: u32 = @bitCast(x);
    const y_magnitude = y_bits & magnitude_mask_f;
    const x_magnitude = x_bits & magnitude_mask_f;

    if (y_magnitude > infinity_bits_f) return .{ .value = @bitCast(y_bits | quiet_nan_bit_f), .resolved = true };
    if (x_magnitude > infinity_bits_f) return .{ .value = @bitCast(x_bits | quiet_nan_bit_f), .resolved = true };

    if (y_magnitude == infinity_bits_f) {
        if (x_magnitude == infinity_bits_f) {
            return .{
                .value = if ((x_bits & sign_mask_f) != 0)
                    signedF(three_pi_over_four_bits_f, y_bits)
                else
                    signedF(pi_over_four_bits_f, y_bits),
                .resolved = true,
            };
        }
        return .{ .value = signedF(pi_over_two_bits_f, y_bits), .resolved = true };
    }

    if (x_magnitude == infinity_bits_f) {
        return .{
            .value = if ((x_bits & sign_mask_f) != 0)
                signedF(pi_bits_f, y_bits)
            else
                @bitCast(y_bits & sign_mask_f),
            .resolved = true,
        };
    }

    if (y_magnitude == 0) {
        if (x_magnitude == 0) {
            return .{
                .value = if ((x_bits & sign_mask_f) != 0)
                    signedF(pi_bits_f, y_bits)
                else
                    y,
                .resolved = true,
            };
        }
        return .{
            .value = if ((x_bits & sign_mask_f) != 0)
                signedF(pi_bits_f, y_bits)
            else
                y,
            .resolved = true,
        };
    }

    if (x_magnitude == 0) return .{ .value = signedF(pi_over_two_bits_f, y_bits), .resolved = true };
    return .{ .value = cr_atan2f(y, x), .resolved = true };
}

/// Return an `Atan2Result(f64)` for correctly rounded atan2(y, x).
///
/// `value` contains the result when `resolved` is true. If the binary64
/// kernel's ambiguity predicate cannot prove one rounded bit, `value` is the
/// canonical positive quiet NaN and `resolved` is false.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn atan2d(y: f64, x: f64) c23.Atan2Result(f64) {
    @setFloatMode(.strict);

    const y_bits: u64 = @bitCast(y);
    const x_bits: u64 = @bitCast(x);
    const y_magnitude = y_bits & magnitude_mask_d;
    const x_magnitude = x_bits & magnitude_mask_d;

    if (y_magnitude > infinity_bits_d) return .{ .value = @bitCast(y_bits | quiet_nan_bit_d), .resolved = true };
    if (x_magnitude > infinity_bits_d) return .{ .value = @bitCast(x_bits | quiet_nan_bit_d), .resolved = true };

    if (y_magnitude == infinity_bits_d) {
        if (x_magnitude == infinity_bits_d) {
            return .{
                .value = if ((x_bits & sign_mask_d) != 0)
                    signedD(three_pi_over_four_bits_d, y_bits)
                else
                    signedD(pi_over_four_bits_d, y_bits),
                .resolved = true,
            };
        }
        return .{ .value = signedD(pi_over_two_bits_d, y_bits), .resolved = true };
    }

    if (x_magnitude == infinity_bits_d) {
        return .{
            .value = if ((x_bits & sign_mask_d) != 0)
                signedD(pi_bits_d, y_bits)
            else
                @bitCast(y_bits & sign_mask_d),
            .resolved = true,
        };
    }

    if (y_magnitude == 0) {
        if (x_magnitude == 0) {
            return .{
                .value = if ((x_bits & sign_mask_d) != 0)
                    signedD(pi_bits_d, y_bits)
                else
                    y,
                .resolved = true,
            };
        }
        return .{
            .value = if ((x_bits & sign_mask_d) != 0)
                signedD(pi_bits_d, y_bits)
            else
                y,
            .resolved = true,
        };
    }

    if (x_magnitude == 0) return .{ .value = signedD(pi_over_two_bits_d, y_bits), .resolved = true };

    // SAFETY: The C call writes `value` before returning nonzero; zero exits before the read.
    var value: f64 = undefined;
    if (zmath_core_atan2d_checked(y, x, &value) == 0)
        return unresolvedD();
    return .{ .value = value, .resolved = true };
}
