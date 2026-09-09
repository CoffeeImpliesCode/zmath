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
    const bits = toBitsh(value);
    const sign = bits & sign_mask_h;
    const magnitude = bits & ~sign_mask_h;

    if (isNaNh(bits)) return quietNaNh(bits);
    if (magnitude == 0 or magnitude == exponent_mask_h) return value;

    var significand = (magnitude & fraction_mask_h) | 0x0400;
    var exponent_field: i64 = @intCast((magnitude & exponent_mask_h) >> 10);
    if (exponent_field == 0) {
        significand = magnitude & fraction_mask_h;
        exponent_field = 1;
        // Zero returned above, so the fraction is in 1..0x03ff. The CLZ
        // difference is 1..10 and normalizes without losing bits. Subtracting
        // it preserves the value and leaves exponent_field in [-9, 0].
        const normalization_shift: u4 =
            @intCast(@clz(significand) - @clz(@as(u16, 0x0400)));
        significand <<= normalization_shift;
        exponent_field -= normalization_shift;
    }

    // The finite exponent field is in [-9, 30], so adding any i32 fits i64.
    const scaled_exponent = exponent_field + @as(i64, exponent);
    if (scaled_exponent > 30) {
        return overflowH(mode, sign != 0);
    }
    if (scaled_exponent >= 1) {
        const result_exponent: u16 = @intCast(scaled_exponent);
        const result_bits = sign |
            (result_exponent << 10) |
            (significand & fraction_mask_h);
        return fromBitsh(result_bits);
    }

    const shift = 1 - scaled_exponent;
    const rounded_significand = roundShiftRightH(
        mode,
        sign != 0,
        significand,
        shift,
    );
    return fromBitsh(sign | rounded_significand);
}

pub fn roundIntegralh(mode: RoundMode, value: f16) f16 {
    const bits = toBitsh(value);
    const sign = bits & sign_mask_h;
    const magnitude = bits & ~sign_mask_h;

    if (isNaNh(bits)) return quietNaNh(bits);
    if (magnitude == 0 or magnitude == exponent_mask_h) return value;

    const exponent_field: i32 = @intCast((magnitude & exponent_mask_h) >> 10);
    const unbiased_exponent = exponent_field - 15;
    if (unbiased_exponent >= 10) return value;

    if (unbiased_exponent < 0) {
        const rounded_magnitude: u16 = switch (mode) {
            .nearest_even => if (unbiased_exponent == -1 and
                magnitude & fraction_mask_h != 0)
                0x3c00
            else
                0,
            .toward_zero => 0,
            .toward_positive => if (sign == 0) 0x3c00 else 0,
            .toward_negative => if (sign != 0) 0x3c00 else 0,
        };
        return fromBitsh(sign | rounded_magnitude);
    }

    const fractional_bits: u4 = @intCast(10 - unbiased_exponent);
    const discarded_mask = (@as(u16, 1) << fractional_bits) - 1;
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
    const rounded_magnitude = truncated + if (increment) unit else 0;
    return fromBitsh(sign | rounded_magnitude);
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
    return bits & ~sign_mask_h > exponent_mask_h;
}

fn isNaNf(bits: u32) bool {
    return bits & ~sign_mask_f > exponent_mask_f;
}

fn isNaNd(bits: u64) bool {
    return bits & ~sign_mask_d > exponent_mask_d;
}

fn quietNaNh(bits: u16) f16 {
    return @bitCast(bits | quiet_mask_h);
}

fn overflowH(mode: RoundMode, negative: bool) f16 {
    const magnitude: u16 = switch (mode) {
        .nearest_even => exponent_mask_h,
        .toward_zero => exponent_mask_h - 1,
        .toward_positive => if (negative) exponent_mask_h - 1 else exponent_mask_h,
        .toward_negative => if (negative) exponent_mask_h else exponent_mask_h - 1,
    };
    return fromBitsh((if (negative) sign_mask_h else 0) | magnitude);
}

fn roundShiftRightH(
    mode: RoundMode,
    negative: bool,
    significand: u16,
    shift: i64,
) u16 {
    if (shift > 11) {
        return switch (mode) {
            .nearest_even, .toward_zero => 0,
            .toward_positive => if (negative) 0 else 1,
            .toward_negative => if (negative) 1 else 0,
        };
    }

    const shift_amount: u4 = @intCast(shift);
    const quotient = significand >> shift_amount;
    const discarded_mask = (@as(u16, 1) << shift_amount) - 1;
    const discarded = significand & discarded_mask;
    const half_shift: u4 = @intCast(shift - 1);
    const half = @as(u16, 1) << half_shift;
    const increment = switch (mode) {
        .nearest_even => discarded > half or
            (discarded == half and quotient & 1 != 0),
        .toward_zero => false,
        .toward_positive => !negative and discarded != 0,
        .toward_negative => negative and discarded != 0,
    };
    return if (increment) quotient + 1 else quotient;
}

fn quietNaNf(bits: u32) f32 {
    return @bitCast(bits | quiet_mask_f);
}

fn quietNaNd(bits: u64) f64 {
    return @bitCast(bits | quiet_mask_d);
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
