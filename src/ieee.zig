const std = @import("std");
const c23_core = @import("ieee_c23_core.zig");
const c23_optional = @import("ieee_c23_optional.zig");
const c23_types = @import("scalar_c23_types.zig");

pub const CanonicalizeResult = c23_types.CanonicalizeResult;
pub const PayloadResult = c23_types.PayloadResult;

pub const canonicalizeh = c23_core.canonicalizeh;
pub const canonicalizef = c23_core.canonicalizef;
pub const canonicalized = c23_core.canonicalized;
pub const nanh = c23_core.nanh;
pub const nanf = c23_core.nanf;
pub const nand = c23_core.nand;
pub const nextTowardh = c23_core.nextTowardh;
pub const nextTowardf = c23_core.nextTowardf;
pub const nextTowardd = c23_core.nextTowardd;
pub const nextUph = c23_core.nextUph;
pub const nextUpf = c23_core.nextUpf;
pub const nextUpd = c23_core.nextUpd;
pub const nextDownh = c23_core.nextDownh;
pub const nextDownf = c23_core.nextDownf;
pub const nextDownd = c23_core.nextDownd;

pub const totalOrderMagh = c23_optional.totalOrderMagh;
pub const totalOrderMagf = c23_optional.totalOrderMagf;
pub const totalOrderMagd = c23_optional.totalOrderMagd;
pub const getPayloadh = c23_optional.getPayloadh;
pub const getPayloadf = c23_optional.getPayloadf;
pub const getPayloadd = c23_optional.getPayloadd;
pub const setPayloadh = c23_optional.setPayloadh;
pub const setPayloadf = c23_optional.setPayloadf;
pub const setPayloadd = c23_optional.setPayloadd;
pub const setPayloadSigh = c23_optional.setPayloadSigh;
pub const setPayloadSigf = c23_optional.setPayloadSigf;
pub const setPayloadSigd = c23_optional.setPayloadSigd;

pub const RoundMode = enum {
    nearest_even,
    toward_zero,
    toward_positive,
    toward_negative,
};

pub const Class = enum {
    zero,
    subnormal,
    normal,
    infinity,
    quiet_nan,
    signaling_nan,
};

const sign_mask_h: u16 = 0x8000;
const exponent_mask_h: u16 = 0x7c00;
const fraction_mask_h: u16 = 0x03ff;
const quiet_mask_h: u16 = 0x0200;

const sign_mask_f: u32 = 0x8000_0000;
const exponent_mask_f: u32 = 0x7f80_0000;
const fraction_mask_f: u32 = 0x007f_ffff;
const quiet_mask_f: u32 = 0x0040_0000;

const sign_mask_d: u64 = 0x8000_0000_0000_0000;
const exponent_mask_d: u64 = 0x7ff0_0000_0000_0000;
const fraction_mask_d: u64 = 0x000f_ffff_ffff_ffff;
const quiet_mask_d: u64 = 0x0008_0000_0000_0000;

// The encoding parameters the two mode-parameterised cores derive every bound
// from. Each entry names the masks declared above, so a field width is written
// once per file, and `bias` is the one number no mask recovers: the exponent
// field is `fraction_bits + 3` bits wide with its top value reserved, which
// makes the largest finite biased exponent twice the bias.
const HFormat = struct {
    const Float = f16;
    const Bits = u16;
    const sign_mask = sign_mask_h;
    const exponent_mask = exponent_mask_h;
    const fraction_mask = fraction_mask_h;
    const quiet_nan_bit = quiet_mask_h;
    const bias: i32 = 15;

    fn toBits(value: Float) Bits {
        return @bitCast(value);
    }

    fn fromBits(bits: Bits) Float {
        return @bitCast(bits);
    }

    fn isNaNBits(bits: Bits) bool {
        return bits & ~sign_mask > exponent_mask;
    }

    fn quietNaN(bits: Bits) Float {
        return @bitCast(bits | quiet_nan_bit);
    }
};

const FFormat = struct {
    const Float = f32;
    const Bits = u32;
    const sign_mask = sign_mask_f;
    const exponent_mask = exponent_mask_f;
    const fraction_mask = fraction_mask_f;
    const quiet_nan_bit = quiet_mask_f;
    const bias: i32 = 127;

    fn toBits(value: Float) Bits {
        return @bitCast(value);
    }

    fn fromBits(bits: Bits) Float {
        return @bitCast(bits);
    }

    fn isNaNBits(bits: Bits) bool {
        return bits & ~sign_mask > exponent_mask;
    }

    fn quietNaN(bits: Bits) Float {
        return @bitCast(bits | quiet_nan_bit);
    }
};

const DFormat = struct {
    const Float = f64;
    const Bits = u64;
    const sign_mask = sign_mask_d;
    const exponent_mask = exponent_mask_d;
    const fraction_mask = fraction_mask_d;
    const quiet_nan_bit = quiet_mask_d;
    const bias: i32 = 1023;

    fn toBits(value: Float) Bits {
        return @bitCast(value);
    }

    fn fromBits(bits: Bits) Float {
        return @bitCast(bits);
    }

    fn isNaNBits(bits: Bits) bool {
        return bits & ~sign_mask > exponent_mask;
    }

    fn quietNaN(bits: Bits) Float {
        return @bitCast(bits | quiet_nan_bit);
    }
};

pub fn toBitsh(value: f16) u16 {
    return @bitCast(value);
}

pub fn fromBitsh(bits: u16) f16 {
    return @bitCast(bits);
}

pub fn classifyh(value: f16) Class {
    const bits: u16 = @bitCast(value);
    const exponent = bits & exponent_mask_h;
    const fraction = bits & fraction_mask_h;

    if (exponent == 0) {
        if (fraction == 0) return .zero;
        return .subnormal;
    }
    if (exponent != exponent_mask_h) return .normal;
    if (fraction == 0) return .infinity;
    if (fraction & quiet_mask_h != 0) return .quiet_nan;
    return .signaling_nan;
}

pub fn classifyf(value: f32) Class {
    const bits: u32 = @bitCast(value);
    const exponent = bits & exponent_mask_f;
    const fraction = bits & fraction_mask_f;

    if (exponent == 0) {
        if (fraction == 0) return .zero;
        return .subnormal;
    }
    if (exponent != exponent_mask_f) return .normal;
    if (fraction == 0) return .infinity;
    if (fraction & quiet_mask_f != 0) return .quiet_nan;
    return .signaling_nan;
}

pub fn classifyd(value: f64) Class {
    const bits: u64 = @bitCast(value);
    const exponent = bits & exponent_mask_d;
    const fraction = bits & fraction_mask_d;

    if (exponent == 0) {
        if (fraction == 0) return .zero;
        return .subnormal;
    }
    if (exponent != exponent_mask_d) return .normal;
    if (fraction == 0) return .infinity;
    if (fraction & quiet_mask_d != 0) return .quiet_nan;
    return .signaling_nan;
}

pub fn copySignh(value: f16, sign: f16) f16 {
    const value_bits: u16 = @bitCast(value);
    const sign_bits: u16 = @bitCast(sign);
    const result_bits = (value_bits & ~sign_mask_h) | (sign_bits & sign_mask_h);
    return @bitCast(result_bits);
}

pub fn copySignf(value: f32, sign: f32) f32 {
    const value_bits: u32 = @bitCast(value);
    const sign_bits: u32 = @bitCast(sign);
    const result_bits = (value_bits & ~sign_mask_f) | (sign_bits & sign_mask_f);
    return @bitCast(result_bits);
}

pub fn copySignd(value: f64, sign: f64) f64 {
    const value_bits: u64 = @bitCast(value);
    const sign_bits: u64 = @bitCast(sign);
    const result_bits = (value_bits & ~sign_mask_d) | (sign_bits & sign_mask_d);
    return @bitCast(result_bits);
}

pub fn nextAfterh(from: f16, to: f16) f16 {
    const from_bits: u16 = @bitCast(from);
    const to_bits: u16 = @bitCast(to);

    if (isNaNh(from_bits)) return quietNaNh(from_bits);
    if (isNaNh(to_bits)) return quietNaNh(to_bits);
    if (from_bits == to_bits) return to;

    const from_magnitude = from_bits & ~sign_mask_h;
    const to_magnitude = to_bits & ~sign_mask_h;
    if (from_magnitude == 0 and to_magnitude == 0) return to;
    if (from_magnitude == 0) {
        return @bitCast((to_bits & sign_mask_h) | 1);
    }

    const toward_greater = totalOrderKeyH(from_bits) < totalOrderKeyH(to_bits);
    const from_negative = from_bits & sign_mask_h != 0;
    const result_bits = if (toward_greater == from_negative)
        from_bits - 1
    else
        from_bits + 1;
    return @bitCast(result_bits);
}

pub fn nextAfterf(from: f32, to: f32) f32 {
    const from_bits: u32 = @bitCast(from);
    const to_bits: u32 = @bitCast(to);

    if (isNaNf(from_bits)) return quietNaNf(from_bits);
    if (isNaNf(to_bits)) return quietNaNf(to_bits);
    if (from_bits == to_bits) return to;

    const from_magnitude = from_bits & ~sign_mask_f;
    const to_magnitude = to_bits & ~sign_mask_f;
    if (from_magnitude == 0 and to_magnitude == 0) return to;
    if (from_magnitude == 0) {
        return @bitCast((to_bits & sign_mask_f) | 1);
    }

    const toward_greater = totalOrderKeyF(from_bits) < totalOrderKeyF(to_bits);
    const from_negative = from_bits & sign_mask_f != 0;
    const result_bits = if (toward_greater == from_negative)
        from_bits - 1
    else
        from_bits + 1;
    return @bitCast(result_bits);
}

pub fn nextAfterd(from: f64, to: f64) f64 {
    const from_bits: u64 = @bitCast(from);
    const to_bits: u64 = @bitCast(to);

    if (isNaNd(from_bits)) return quietNaNd(from_bits);
    if (isNaNd(to_bits)) return quietNaNd(to_bits);
    if (from_bits == to_bits) return to;

    const from_magnitude = from_bits & ~sign_mask_d;
    const to_magnitude = to_bits & ~sign_mask_d;
    if (from_magnitude == 0 and to_magnitude == 0) return to;
    if (from_magnitude == 0) {
        return @bitCast((to_bits & sign_mask_d) | 1);
    }

    const toward_greater = totalOrderKeyD(from_bits) < totalOrderKeyD(to_bits);
    const from_negative = from_bits & sign_mask_d != 0;
    const result_bits = if (toward_greater == from_negative)
        from_bits - 1
    else
        from_bits + 1;
    return @bitCast(result_bits);
}

pub fn scaleBh(mode: RoundMode, value: f16, exponent: i32) f16 {
    return scaleB(HFormat, mode, value, exponent);
}

pub fn scaleBf(mode: RoundMode, value: f32, exponent: i32) f32 {
    return scaleB(FFormat, mode, value, exponent);
}

pub fn scaleBd(mode: RoundMode, value: f64, exponent: i32) f64 {
    return scaleB(DFormat, mode, value, exponent);
}

pub fn roundIntegralh(mode: RoundMode, value: f16) f16 {
    return roundIntegral(HFormat, mode, value);
}

pub fn roundIntegralf(mode: RoundMode, value: f32) f32 {
    return roundIntegral(FFormat, mode, value);
}

pub fn roundIntegrald(mode: RoundMode, value: f64) f64 {
    return roundIntegral(DFormat, mode, value);
}

pub fn totalOrderh(left: f16, right: f16) bool {
    const left_bits: u16 = @bitCast(left);
    const right_bits: u16 = @bitCast(right);
    return totalOrderKeyH(left_bits) <= totalOrderKeyH(right_bits);
}

pub fn totalOrderf(left: f32, right: f32) bool {
    const left_bits: u32 = @bitCast(left);
    const right_bits: u32 = @bitCast(right);
    return totalOrderKeyF(left_bits) <= totalOrderKeyF(right_bits);
}

pub fn totalOrderd(left: f64, right: f64) bool {
    const left_bits: u64 = @bitCast(left);
    const right_bits: u64 = @bitCast(right);
    return totalOrderKeyD(left_bits) <= totalOrderKeyD(right_bits);
}

fn isNaNh(bits: u16) bool {
    return HFormat.isNaNBits(bits);
}

fn isNaNf(bits: u32) bool {
    return FFormat.isNaNBits(bits);
}

fn isNaNd(bits: u64) bool {
    return DFormat.isNaNBits(bits);
}

fn quietNaNh(bits: u16) f16 {
    return HFormat.quietNaN(bits);
}

fn quietNaNf(bits: u32) f32 {
    return FFormat.quietNaN(bits);
}

fn quietNaNd(bits: u64) f64 {
    return DFormat.quietNaN(bits);
}

// The saturation magnitude of an overflow, for one mode and one sign.
// Round-to-nearest and the directed modes all stop at the format boundary:
// infinity when the mode lets the magnitude grow, the largest finite value
// when it does not. The exponent mask is the whole exponent field, so one
// unit below it is the largest finite magnitude.
fn overflow(comptime Format: type, mode: RoundMode, negative: bool) Format.Float {
    const maximum_finite = Format.exponent_mask - 1;
    const magnitude: Format.Bits = switch (mode) {
        .nearest_even => Format.exponent_mask,
        .toward_zero => maximum_finite,
        .toward_positive => if (negative) maximum_finite else Format.exponent_mask,
        .toward_negative => if (negative) Format.exponent_mask else maximum_finite,
    };
    return Format.fromBits((if (negative) Format.sign_mask else 0) | magnitude);
}

// The single rounding a right shift needs. The significand is in
// [1, 2^(fraction_bits+1) - 1] and the shift is the whole exponent drop,
// which no i32 keeps small.
fn roundShiftRight(
    comptime Format: type,
    mode: RoundMode,
    negative: bool,
    significand: Format.Bits,
    shift: i64,
) Format.Bits {
    const Bits = Format.Bits;
    const Shift = std.math.Log2Int(Bits);
    // The fraction is the run of low bits under the exponent field, so its
    // width is the trailing-zero count of the exponent mask: 10, 23, 52.
    const fraction_bits: Shift = @intCast(@ctz(Format.exponent_mask));

    // A drop of more than fraction_bits + 1 leaves less than half of the
    // smallest subnormal quantum, so only the mode's own direction reaches
    // one quantum and the other two reach none.
    if (shift > @as(i64, fraction_bits) + 1) {
        return switch (mode) {
            .nearest_even, .toward_zero => 0,
            .toward_positive => if (negative) 0 else 1,
            .toward_negative => if (negative) 1 else 0,
        };
    }

    // The return above bounds the shift to [1, fraction_bits + 1]. That upper
    // bound is 11, 24, or 53 against shift-type maxima 15, 31, and 63, so both
    // casts hold and no shift count reaches the bit width of its own type.
    const shift_amount: Shift = @intCast(shift);
    const quotient = significand >> shift_amount;
    const discarded_mask = (@as(Bits, 1) << shift_amount) - 1;
    const discarded = significand & discarded_mask;
    const half_shift: Shift = @intCast(shift - 1);
    const half = @as(Bits, 1) << half_shift;
    const increment = switch (mode) {
        .nearest_even => discarded > half or
            (discarded == half and quotient & 1 != 0),
        .toward_zero => false,
        .toward_positive => !negative and discarded != 0,
        .toward_negative => negative and discarded != 0,
    };
    // The quotient is at most (2^(fraction_bits+1) - 1) >> 1, which is
    // 2^fraction_bits - 1, so the carry stops at the hidden bit: a caller
    // that reads it into a zero exponent field gets the minimum normal.
    return if (increment) quotient + 1 else quotient;
}

// The shared `scaleB` core: an exact power-of-two scale with one rounding, in
// the mode the caller names. The exponent carried here is the biased exponent
// of the width being compiled, so every bound below is a statement about that
// width and not about binary16.
fn scaleB(
    comptime Format: type,
    mode: RoundMode,
    value: Format.Float,
    exponent: i32,
) Format.Float {
    const Bits = Format.Bits;
    const Shift = std.math.Log2Int(Bits);
    // The fraction is the run of low bits under the exponent field, so its
    // width is the trailing-zero count of the exponent mask: 10, 23, 52. Each
    // is below the bit width it shifts inside, so all three fit the shift type.
    const fraction_bits: Shift = @intCast(@ctz(Format.exponent_mask));
    const hidden_bit: Bits = @as(Bits, 1) << fraction_bits;
    // The largest finite biased exponent. The exponent field is
    // `fraction_bits + 3` bits wide and its top value is reserved, which makes
    // twice the bias 30, 254, or 2046 the last field a finite result can hold.
    const maximum_biased: i64 = @as(i64, Format.bias) * 2;

    const bits = Format.toBits(value);
    const sign = bits & Format.sign_mask;
    const magnitude = bits & ~Format.sign_mask;

    if (Format.isNaNBits(bits)) return Format.quietNaN(bits);
    if (magnitude == 0 or magnitude == Format.exponent_mask) return value;

    var significand: Bits = (magnitude & Format.fraction_mask) | hidden_bit;
    // NaN and infinity returned above, so a finite nonzero input carries a
    // biased exponent in [1, maximum_biased], and that is at most 2046.
    var exponent_field: i64 = @intCast((magnitude & Format.exponent_mask) >> fraction_bits);
    if (exponent_field == 0) {
        significand = magnitude & Format.fraction_mask;
        exponent_field = 1;
        // Zero returned above, so the fraction is in 1..2^fraction_bits - 1. The
        // CLZ difference is 1..fraction_bits and normalizes without losing
        // bits. Subtracting it preserves the value and leaves exponent_field in
        // [1 - fraction_bits, 0].
        const normalization_shift: Shift = @intCast(@clz(significand) - @clz(hidden_bit));
        significand <<= normalization_shift;
        exponent_field -= normalization_shift;
    }

    // The finite exponent field is in [1 - fraction_bits, maximum_biased], so
    // adding any i32 stays inside i64.
    const scaled_exponent = exponent_field + @as(i64, exponent);
    if (scaled_exponent > maximum_biased) {
        return overflow(Format, mode, sign != 0);
    }
    if (scaled_exponent >= 1) {
        // This branch bounds scaled_exponent to [1, maximum_biased], at most
        // 2046, so the cast holds. The shift cannot reach the sign bit either:
        // maximum_biased << fraction_bits is 0x7800, 0x7f00_0000, and
        // 0x7fe0_0000_0000_0000, all inside the format.
        const result_exponent: Bits = @intCast(scaled_exponent);
        const result_bits = sign |
            (result_exponent << fraction_bits) |
            (significand & Format.fraction_mask);
        return Format.fromBits(result_bits);
    }

    // The exact result is subnormal or zero, and one right shift is the only
    // rounding it can need. The shift is at most 1 + 2^31 + fraction_bits,
    // which stays far inside i64.
    const shift = 1 - scaled_exponent;
    const rounded_significand = roundShiftRight(Format, mode, sign != 0, significand, shift);
    return Format.fromBits(sign | rounded_significand);
}

// The shared `roundIntegral` core: one rounding to an integer, in the mode
// the caller names.
fn roundIntegral(
    comptime Format: type,
    mode: RoundMode,
    value: Format.Float,
) Format.Float {
    const Bits = Format.Bits;
    const Shift = std.math.Log2Int(Bits);
    // The fraction is the run of low bits under the exponent field, so its
    // width is the trailing-zero count of the exponent mask: 10, 23, 52.
    const fraction_bits: Shift = @intCast(@ctz(Format.exponent_mask));
    // One has a zero fraction at the unbiased exponent zero, so its field is
    // the bias shifted up by the fraction width: 0x3c00, 0x3f80_0000, and
    // 0x3ff0_0000_0000_0000.
    const one: Bits = @as(Bits, Format.bias) << fraction_bits;

    const bits = Format.toBits(value);
    const sign = bits & Format.sign_mask;
    const magnitude = bits & ~Format.sign_mask;

    if (Format.isNaNBits(bits)) return Format.quietNaN(bits);
    if (magnitude == 0 or magnitude == Format.exponent_mask) return value;

    // NaN and infinity returned above, so the field is in [0, 2 * bias], which
    // is at most 2046 and so fits i32 with the bias subtracted from it.
    const exponent_field: i32 = @intCast((magnitude & Format.exponent_mask) >> fraction_bits);
    const unbiased_exponent = exponent_field - Format.bias;
    // A value at or above one ulp of the bias exponent is a multiple of one
    // already, so it needs no rounding.
    if (unbiased_exponent >= @as(i32, fraction_bits)) return value;

    if (unbiased_exponent < 0) {
        const rounded_magnitude: Bits = switch (mode) {
            .nearest_even => if (unbiased_exponent == -1 and
                magnitude & Format.fraction_mask != 0)
                one
            else
                0,
            .toward_zero => 0,
            .toward_positive => if (sign == 0) one else 0,
            .toward_negative => if (sign != 0) one else 0,
        };
        return Format.fromBits(sign | rounded_magnitude);
    }

    // The two returns above leave unbiased_exponent in [0, fraction_bits - 1],
    // so the discarded width is in [1, fraction_bits]: 1..10, 1..23, or 1..52,
    // each inside its own shift type, whose maxima are 15, 31, and 63.
    const fractional_bits: Shift = @intCast(@as(i32, fraction_bits) - unbiased_exponent);
    const discarded_mask = (@as(Bits, 1) << fractional_bits) - 1;
    const discarded = magnitude & discarded_mask;
    if (discarded == 0) return value;

    const truncated = magnitude & ~discarded_mask;
    const unit = discarded_mask + 1;
    const half = unit >> 1;
    const increment = switch (mode) {
        .nearest_even => discarded > half or
            (discarded == half and (truncated >> fractional_bits) & 1 != 0),
        .toward_zero => false,
        .toward_positive => sign == 0,
        .toward_negative => sign != 0,
    };
    // The carry is one whole quantum and the truncated magnitude has a zero
    // fraction, so an incremented magnitude stays inside its own exponent
    // field and never reaches the reserved top one.
    const rounded_magnitude = truncated + if (increment) unit else 0;
    return Format.fromBits(sign | rounded_magnitude);
}

fn totalOrderKeyH(bits: u16) u16 {
    if (bits & sign_mask_h != 0) return ~bits;
    return bits | sign_mask_h;
}

fn totalOrderKeyF(bits: u32) u32 {
    if (bits & sign_mask_f != 0) return ~bits;
    return bits | sign_mask_f;
}

fn totalOrderKeyD(bits: u64) u64 {
    if (bits & sign_mask_d != 0) return ~bits;
    return bits | sign_mask_d;
}
