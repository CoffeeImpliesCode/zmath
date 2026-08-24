const common = @import("scalar_c23_common.zig");
const types = @import("scalar_c23_types.zig");

const double_sign_mask: u64 = 0x8000_0000_0000_0000;
const double_payload_mask: u64 = 0x0007_ffff_ffff_ffff;

fn canonicalize(
    comptime Format: type,
    value: Format.Float,
) types.CanonicalizeResult(Format.Float) {
    const bits = Format.toBits(value);
    return .{
        .value = if (Format.isNaNBits(bits))
            Format.fromBits(Format.quietNaNBits(bits))
        else
            value,
        .ok = true,
    };
}

fn parseNanPayload(comptime Bits: type, comptime payload_mask: Bits, tag: []const u8) ?Bits {
    if (tag.len == 0) return null;

    var payload: Bits = 0;
    for (tag) |character| {
        const digit: Bits = switch (character) {
            '0'...'9' => character - '0',
            'A'...'F' => character - 'A' + 10,
            'a'...'f' => character - 'a' + 10,
            else => return null,
        };
        if (payload > payload_mask >> 4) return null;
        payload = (payload << 4) | digit;
        if (payload > payload_mask) return null;
    }
    return payload;
}

fn nan(comptime Format: type, comptime payload_mask: Format.Bits, tag: []const u8) Format.Float {
    const payload = parseNanPayload(Format.Bits, payload_mask, tag) orelse
        return Format.canonicalNaN();
    const canonical_bits = Format.toBits(Format.canonicalNaN());
    return Format.fromBits(canonical_bits | payload);
}

fn nextUp(
    comptime Format: type,
    comptime sign_mask: Format.Bits,
    value: Format.Float,
) Format.Float {
    const bits = Format.toBits(value);

    if (Format.isNaNBits(bits)) return Format.fromBits(Format.quietNaNBits(bits));
    if (Format.isZeroBits(bits)) return Format.fromBits(1);
    if (Format.isPositiveInfinityBits(bits)) return value;
    if (bits & sign_mask != 0) return Format.fromBits(bits - 1);
    return Format.fromBits(bits + 1);
}

fn nextDown(
    comptime Format: type,
    comptime sign_mask: Format.Bits,
    value: Format.Float,
) Format.Float {
    const bits = Format.toBits(value);

    if (Format.isNaNBits(bits)) return Format.fromBits(Format.quietNaNBits(bits));
    if (Format.isZeroBits(bits)) return Format.fromBits(sign_mask | 1);
    if (Format.isNegativeInfinityBits(bits)) return value;
    if (bits & sign_mask != 0) return Format.fromBits(bits + 1);
    return Format.fromBits(bits - 1);
}

/// When a `nextToward` target is wider than H or F, preserve the target sign
/// and high-order representable payload bits, then force the target qNaN bit.
fn narrowTowardNaN(
    comptime Format: type,
    comptime sign_mask: Format.Bits,
    comptime payload_mask: Format.Bits,
    comptime payload_shift: u6,
    toward_bits: u64,
) Format.Float {
    const sign: Format.Bits = if (toward_bits & double_sign_mask != 0) sign_mask else 0;
    const payload: Format.Bits = @truncate((toward_bits & double_payload_mask) >> payload_shift);
    const canonical_bits = Format.toBits(Format.canonicalNaN());
    return Format.fromBits(sign | canonical_bits | (payload & payload_mask));
}

fn nextToward(
    comptime Format: type,
    comptime sign_mask: Format.Bits,
    comptime payload_mask: Format.Bits,
    comptime payload_shift: u6,
    value: Format.Float,
    toward: f64,
) Format.Float {
    const value_bits = Format.toBits(value);
    if (Format.isNaNBits(value_bits)) {
        return Format.fromBits(Format.quietNaNBits(value_bits));
    }

    const toward_bits = common.D.toBits(toward);
    if (common.D.isNaNBits(toward_bits)) {
        return narrowTowardNaN(Format, sign_mask, payload_mask, payload_shift, toward_bits);
    }

    const value_wide: f64 = value;
    if (value_wide == toward) {
        if (value_wide == 0.0) {
            const sign = if (toward_bits & double_sign_mask != 0) sign_mask else 0;
            return Format.fromBits(sign);
        }
        return value;
    }
    if (value_wide < toward) return nextUp(Format, sign_mask, value);
    return nextDown(Format, sign_mask, value);
}

pub fn canonicalizeh(value: f16) types.CanonicalizeResult(f16) {
    return canonicalize(common.H, value);
}

pub fn canonicalizef(value: f32) types.CanonicalizeResult(f32) {
    return canonicalize(common.F, value);
}

pub fn canonicalized(value: f64) types.CanonicalizeResult(f64) {
    return canonicalize(common.D, value);
}

pub fn nanh(tag: []const u8) f16 {
    return nan(common.H, 0x01ff, tag);
}

pub fn nanf(tag: []const u8) f32 {
    return nan(common.F, 0x003f_ffff, tag);
}

pub fn nand(tag: []const u8) f64 {
    return nan(common.D, double_payload_mask, tag);
}

pub fn nextTowardh(value: f16, toward: f64) f16 {
    return nextToward(common.H, 0x8000, 0x01ff, 42, value, toward);
}

pub fn nextTowardf(value: f32, toward: f64) f32 {
    return nextToward(common.F, 0x8000_0000, 0x003f_ffff, 29, value, toward);
}

pub fn nextTowardd(value: f64, toward: f64) f64 {
    return nextToward(common.D, double_sign_mask, double_payload_mask, 0, value, toward);
}

pub fn nextUph(value: f16) f16 {
    return nextUp(common.H, 0x8000, value);
}

pub fn nextUpf(value: f32) f32 {
    return nextUp(common.F, 0x8000_0000, value);
}

pub fn nextUpd(value: f64) f64 {
    return nextUp(common.D, double_sign_mask, value);
}

pub fn nextDownh(value: f16) f16 {
    return nextDown(common.H, 0x8000, value);
}

pub fn nextDownf(value: f32) f32 {
    return nextDown(common.F, 0x8000_0000, value);
}

pub fn nextDownd(value: f64) f64 {
    return nextDown(common.D, double_sign_mask, value);
}
