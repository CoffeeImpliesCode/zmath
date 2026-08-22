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

const sign_mask_f: u32 = 0x8000_0000;
const exponent_mask_f: u32 = 0x7f80_0000;
const fraction_mask_f: u32 = 0x007f_ffff;
const quiet_mask_f: u32 = 0x0040_0000;

const sign_mask_d: u64 = 0x8000_0000_0000_0000;
const exponent_mask_d: u64 = 0x7ff0_0000_0000_0000;
const fraction_mask_d: u64 = 0x000f_ffff_ffff_ffff;
const quiet_mask_d: u64 = 0x0008_0000_0000_0000;

pub fn classifyF(value: f32) Class {
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

pub fn classifyD(value: f64) Class {
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

pub fn copySignF(value: f32, sign: f32) f32 {
    const value_bits: u32 = @bitCast(value);
    const sign_bits: u32 = @bitCast(sign);
    const result_bits = (value_bits & ~sign_mask_f) | (sign_bits & sign_mask_f);
    return @bitCast(result_bits);
}

pub fn copySignD(value: f64, sign: f64) f64 {
    const value_bits: u64 = @bitCast(value);
    const sign_bits: u64 = @bitCast(sign);
    const result_bits = (value_bits & ~sign_mask_d) | (sign_bits & sign_mask_d);
    return @bitCast(result_bits);
}

pub fn nextAfterF(from: f32, to: f32) f32 {
    const from_bits: u32 = @bitCast(from);
    const to_bits: u32 = @bitCast(to);

    if (isNaNF(from_bits)) return quietNaNF(from_bits);
    if (isNaNF(to_bits)) return quietNaNF(to_bits);
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

pub fn nextAfterD(from: f64, to: f64) f64 {
    const from_bits: u64 = @bitCast(from);
    const to_bits: u64 = @bitCast(to);

    if (isNaND(from_bits)) return quietNaND(from_bits);
    if (isNaND(to_bits)) return quietNaND(to_bits);
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

pub fn totalOrderF(left: f32, right: f32) bool {
    const left_bits: u32 = @bitCast(left);
    const right_bits: u32 = @bitCast(right);
    return totalOrderKeyF(left_bits) <= totalOrderKeyF(right_bits);
}

pub fn totalOrderD(left: f64, right: f64) bool {
    const left_bits: u64 = @bitCast(left);
    const right_bits: u64 = @bitCast(right);
    return totalOrderKeyD(left_bits) <= totalOrderKeyD(right_bits);
}

fn isNaNF(bits: u32) bool {
    return bits & ~sign_mask_f > exponent_mask_f;
}

fn isNaND(bits: u64) bool {
    return bits & ~sign_mask_d > exponent_mask_d;
}

fn quietNaNF(bits: u32) f32 {
    return @bitCast(bits | quiet_mask_f);
}

fn quietNaND(bits: u64) f64 {
    return @bitCast(bits | quiet_mask_d);
}

fn totalOrderKeyF(bits: u32) u32 {
    if (bits & sign_mask_f != 0) return ~bits;
    return bits | sign_mask_f;
}

fn totalOrderKeyD(bits: u64) u64 {
    if (bits & sign_mask_d != 0) return ~bits;
    return bits | sign_mask_d;
}
