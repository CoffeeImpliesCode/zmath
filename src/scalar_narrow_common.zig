const std = @import("std");

/// The narrowing engine keeps every finite operand as an integer significand
/// multiplied by an explicit power of two. A hardware square root seeds the
/// integer root, but exact integer comparisons still decide every sqrt result.
/// Internal alignment may compress bits below the fixed 128-bit headroom
/// into one sticky bit; only the final destination packers round.
pub const Format = struct {
    bits: u8,
    fraction_bits: u8,
    exponent_bits: u8,
    bias: i32,

    pub fn Bits(comptime format: Format) type {
        return switch (format.bits) {
            16 => u16,
            32 => u32,
            64 => u64,
            else => @compileError("unsupported narrowing format"),
        };
    }

    pub fn signMask(comptime format: Format) u64 {
        return @as(u64, 1) << (format.bits - 1);
    }

    pub fn fractionMask(comptime format: Format) u64 {
        return (@as(u64, 1) << format.fraction_bits) - 1;
    }

    pub fn exponentMask(comptime format: Format) u64 {
        return ((@as(u64, 1) << format.exponent_bits) - 1) << format.fraction_bits;
    }

    pub fn infinityBits(comptime format: Format) u64 {
        return exponentMask(format);
    }

    pub fn quietBit(comptime format: Format) u64 {
        return @as(u64, 1) << (format.fraction_bits - 1);
    }

    pub fn maxExponent(comptime format: Format) i32 {
        return ((@as(i32, 1) << format.exponent_bits) - 2) - format.bias;
    }

    pub fn minNormalExponent(comptime format: Format) i32 {
        return 1 - format.bias;
    }
};

pub const binary16_format: Format = .{ .bits = 16, .fraction_bits = 10, .exponent_bits = 5, .bias = 15 };
pub const binary32_format: Format = .{ .bits = 32, .fraction_bits = 23, .exponent_bits = 8, .bias = 127 };
pub const binary64_format: Format = .{ .bits = 64, .fraction_bits = 52, .exponent_bits = 11, .bias = 1023 };

pub const Class = enum {
    zero,
    finite,
    infinity,
    nan,
};

pub const Decoded = struct {
    class: Class,
    negative: bool,
    significand: u64,
    exponent: i32,
    bits: u64,
};

pub const Dyadic = struct {
    negative: bool,
    magnitude: u128,
    exponent: i32,
    sticky: bool,
};

pub fn decode(comptime format: Format, bits: u64) Decoded {
    const sign = bits & Format.signMask(format) != 0;
    const exponent_field = (bits & Format.exponentMask(format)) >> format.fraction_bits;
    const fraction = bits & Format.fractionMask(format);
    const all_ones = (@as(u64, 1) << format.exponent_bits) - 1;

    if (exponent_field == 0) {
        if (fraction == 0) {
            return .{ .class = .zero, .negative = sign, .significand = 0, .exponent = 0, .bits = bits };
        }
        return .{
            .class = .finite,
            .negative = sign,
            .significand = fraction,
            .exponent = 1 - format.bias - @as(i32, format.fraction_bits),
            .bits = bits,
        };
    }
    if (exponent_field == all_ones) {
        return .{
            .class = if (fraction == 0) .infinity else .nan,
            .negative = sign,
            .significand = fraction,
            .exponent = 0,
            .bits = bits,
        };
    }
    return .{
        .class = .finite,
        .negative = sign,
        .significand = (@as(u64, 1) << format.fraction_bits) | fraction,
        .exponent = @as(i32, @intCast(exponent_field)) - format.bias - @as(i32, format.fraction_bits),
        .bits = bits,
    };
}

pub fn toBits(comptime format: Format, value: anytype) u64 {
    return @as(format.Bits(), @bitCast(value));
}

pub fn fromBits(comptime format: Format, bits: anytype) format.Bits() {
    return @bitCast(@as(format.Bits(), @intCast(bits)));
}

pub fn invalid(comptime format: Format) format.Bits() {
    return fromBits(format, Format.infinityBits(format) | Format.quietBit(format));
}

/// Preserve the sign and the high representable payload bits before quieting.
pub fn convertedNaN(
    comptime source: Format,
    comptime destination: Format,
    bits: u64,
) destination.Bits() {
    const source_fraction = bits & Format.fractionMask(source);
    const dropped = source.fraction_bits - destination.fraction_bits;
    const payload = (source_fraction >> dropped) | Format.quietBit(destination);
    const sign = bits & Format.signMask(source) != 0;
    const result = Format.infinityBits(destination) | payload |
        (if (sign) Format.signMask(destination) else 0);
    return fromBits(destination, result);
}

pub fn zero(comptime format: Format, negative: bool) format.Bits() {
    return fromBits(format, if (negative) Format.signMask(format) else 0);
}

fn bitLength(value: u128) u32 {
    if (value == 0) return 0;
    return 128 - @as(u32, @intCast(@clz(value)));
}

fn topExponent(value: Dyadic) i32 {
    if (value.magnitude == 0) return value.exponent;
    return value.exponent + @as(i32, @intCast(bitLength(value.magnitude) - 1));
}

fn lowBitsNonzero(value: u128, shift: u32) bool {
    if (shift == 0) return false;
    if (shift >= 128) return value != 0;
    const mask = (@as(u128, 1) << @as(std.math.Log2Int(u128), @intCast(shift))) - 1;
    return value & mask != 0;
}

fn roundRight(value: u128, shift: u32, sticky: bool) struct { value: u128, up: bool } {
    if (shift == 0) return .{ .value = value, .up = false };
    if (shift > 128) return .{ .value = 0, .up = false };
    if (shift == 128) {
        // The truncated value is zero, so an exact half rounds down. A
        // nonzero tail makes that half strictly above the tie.
        const half = @as(u128, 1) << 127;
        return .{
            .value = 0,
            .up = value > half or (value == half and sticky),
        };
    }

    const shift_type: std.math.Log2Int(u128) = @intCast(shift);
    const truncated = value >> shift_type;
    const remainder_mask = (@as(u128, 1) << shift_type) - 1;
    const remainder = value & remainder_mask;
    const half = @as(u128, 1) << @as(std.math.Log2Int(u128), @intCast(shift - 1));
    const above = remainder > half;
    const tie = remainder == half;
    return .{
        .value = truncated,
        .up = above or (tie and (sticky or (truncated & 1) != 0)),
    };
}

fn packRounded(
    comptime destination: Format,
    negative: bool,
    rounded: u128,
    exponent: i32,
    subnormal: bool,
) destination.Bits() {
    const sign_field = if (negative) Format.signMask(destination) else 0;
    const fraction_mask = Format.fractionMask(destination);
    const hidden = @as(u128, 1) << destination.fraction_bits;

    if (rounded == 0) return fromBits(destination, sign_field);
    if (subnormal) {
        if (rounded < hidden) {
            return fromBits(destination, sign_field | @as(u64, @intCast(rounded)));
        }
        return packRounded(destination, negative, rounded, Format.minNormalExponent(destination), false);
    }

    var significand = rounded;
    var result_exponent = exponent;
    if (significand >= (hidden << 1)) {
        significand >>= 1;
        result_exponent += 1;
    }
    if (result_exponent > Format.maxExponent(destination)) {
        return fromBits(destination, sign_field | Format.infinityBits(destination));
    }
    if (result_exponent < Format.minNormalExponent(destination)) {
        return packRounded(destination, negative, significand, result_exponent, true);
    }
    std.debug.assert(significand >= hidden);
    std.debug.assert(significand < (hidden << 1));
    const exponent_field = @as(u64, @intCast(result_exponent + destination.bias)) << destination.fraction_bits;
    return fromBits(destination, sign_field | exponent_field | (@as(u64, @intCast(significand)) & fraction_mask));
}

/// Round an exact dyadic magnitude once to the destination format.
pub fn packDyadic(comptime destination: Format, value: Dyadic) destination.Bits() {
    if (value.magnitude == 0) return zero(destination, value.negative);

    const top = topExponent(value);
    const minimum = Format.minNormalExponent(destination);
    if (top < minimum) {
        const scale = minimum - @as(i32, destination.fraction_bits);
        var rounded: u128 = 0;
        if (value.exponent >= scale) {
            const shift: u32 = @intCast(value.exponent - scale);
            std.debug.assert(shift < 128);
            const shifted = value.magnitude << @as(std.math.Log2Int(u128), @intCast(shift));
            rounded = shifted;
        } else {
            const shift: u32 = @intCast(scale - value.exponent);
            const result = roundRight(value.magnitude, shift, value.sticky);
            rounded = result.value + @as(u128, @intFromBool(result.up));
        }
        return packRounded(destination, value.negative, rounded, minimum, true);
    }

    const shift_signed = top - @as(i32, destination.fraction_bits) - value.exponent;
    // SAFETY: Both shift branches assign `rounded` before `packRounded` reads it.
    var rounded: u128 = undefined;
    if (shift_signed > 0) {
        const shift: u32 = @intCast(shift_signed);
        const result = roundRight(value.magnitude, shift, value.sticky);
        rounded = result.value + @as(u128, @intFromBool(result.up));
    } else {
        // `shift_signed` is `bitLength(magnitude) - 1 - fraction_bits` here, so
        // it lies in `[-52, 0]` and its magnitude is the left shift. `@abs` takes
        // that magnitude without the `-x` overflow at `minInt(i32)`, and the
        // assert below is the width bound that keeps the shift inside a `u128`.
        const shift: u32 = @abs(shift_signed);
        std.debug.assert(shift < 128);
        rounded = value.magnitude << @as(std.math.Log2Int(u128), @intCast(shift));
    }
    return packRounded(destination, value.negative, rounded, top, false);
}

fn exactDyadic(comptime source: Format, value: Decoded) Dyadic {
    _ = source;
    return .{
        .negative = value.negative,
        .magnitude = value.significand,
        .exponent = value.exponent,
        .sticky = false,
    };
}
fn alignDyadic(value: Dyadic, exponent: i32) struct { magnitude: u128, sticky: bool } {
    if (value.exponent >= exponent) {
        const shift: u32 = @intCast(value.exponent - exponent);
        std.debug.assert(bitLength(value.magnitude) + shift <= 128);
        return .{
            .magnitude = value.magnitude << @as(std.math.Log2Int(u128), @intCast(shift)),
            .sticky = value.sticky,
        };
    }
    const shift: u32 = @intCast(exponent - value.exponent);
    if (shift >= 128) return .{ .magnitude = 0, .sticky = value.magnitude != 0 or value.sticky };
    const rounded = value.magnitude >> @as(std.math.Log2Int(u128), @intCast(shift));
    return .{
        .magnitude = rounded,
        .sticky = value.sticky or lowBitsNonzero(value.magnitude, shift),
    };
}

/// Keep 128 high-order bits. Every supported destination needs at most 53
/// significant bits, so discarded alignment bits affect only the sticky bit.
pub fn addDyadic(first: Dyadic, second: Dyadic) Dyadic {
    if (first.magnitude == 0) return second;
    if (second.magnitude == 0) return first;

    const first_top = topExponent(first);
    const second_top = topExponent(second);
    const high = if (first_top >= second_top) first else second;
    const low = if (first_top >= second_top) second else first;
    const minimum = @min(high.exponent, low.exponent);
    const exponent = @max(minimum, topExponent(high) - 127);
    const high_aligned = alignDyadic(high, exponent);
    const low_aligned = alignDyadic(low, exponent);
    if (high.negative == low.negative) {
        const sticky = high_aligned.sticky or low_aligned.sticky;
        if (high_aligned.magnitude > ~@as(u128, 0) - low_aligned.magnitude) {
            const wrapped = high_aligned.magnitude +% low_aligned.magnitude;
            return .{
                .negative = high.negative,
                .magnitude = (wrapped >> 1) | (@as(u128, 1) << 127),
                .exponent = exponent + 1,
                .sticky = sticky or (wrapped & 1) != 0,
            };
        }
        return .{
            .negative = high.negative,
            .magnitude = high_aligned.magnitude + low_aligned.magnitude,
            .exponent = exponent,
            .sticky = sticky,
        };
    }
    // The high operand is left-aligned by construction. For an opposite-sign
    // result, only the low operand can have discarded bits. Those bits are
    // subtracted, so delta - epsilon is represented as (delta - 1) plus the
    // positive complement of epsilon. Keeping that borrow is required for
    // ties at the destination boundary and is shared by every FMA path.
    std.debug.assert(!high_aligned.sticky);
    if (high_aligned.magnitude >= low_aligned.magnitude) {
        var magnitude = high_aligned.magnitude - low_aligned.magnitude;
        var result_sticky = false;
        if (low_aligned.sticky) {
            std.debug.assert(magnitude != 0);
            magnitude -= 1;
            result_sticky = true;
        }
        return .{
            .negative = high.negative,
            .magnitude = magnitude,
            .exponent = exponent,
            .sticky = result_sticky,
        };
    }
    std.debug.assert(!low_aligned.sticky);
    return .{
        .negative = low.negative,
        .magnitude = low_aligned.magnitude - high_aligned.magnitude,
        .exponent = exponent,
        .sticky = false,
    };
}

pub fn multiplyDyadic(first: Dyadic, second: Dyadic) Dyadic {
    return .{
        .negative = first.negative != second.negative,
        .magnitude = first.magnitude * second.magnitude,
        .exponent = first.exponent + second.exponent,
        .sticky = first.sticky or second.sticky,
    };
}

fn ratioTop(numerator: u128, denominator: u128) i32 {
    const numerator_bits = bitLength(numerator);
    const denominator_bits = bitLength(denominator);
    var candidate: i32 = @as(i32, @intCast(numerator_bits)) - @as(i32, @intCast(denominator_bits));
    if (candidate >= 0) {
        const shift: u32 = @intCast(candidate);
        const shifted = denominator << @as(std.math.Log2Int(u128), @intCast(shift));
        if (numerator < shifted) candidate -= 1;
    } else {
        const shift: u32 = @intCast(-candidate);
        const shifted = numerator << @as(std.math.Log2Int(u128), @intCast(shift));
        if (shifted < denominator) candidate -= 1;
    }
    return candidate;
}

fn roundRatioAtScale(numerator: u128, denominator: u128, shift: i32) struct { value: u128, up: bool } {
    if (shift >= 0) {
        const shift_count: u32 = @intCast(shift);
        std.debug.assert(bitLength(numerator) + shift_count <= 128);
        const scaled = numerator << @as(std.math.Log2Int(u128), @intCast(shift_count));
        const quotient = scaled / denominator;
        const remainder = scaled % denominator;
        const doubled = remainder * 2;
        return .{ .value = quotient, .up = doubled > denominator or (doubled == denominator and (quotient & 1) != 0) };
    }

    const shift_count: u32 = @abs(shift);
    if (bitLength(denominator) + shift_count > 128) {
        return .{ .value = 0, .up = false };
    }
    const scaled_denominator = denominator << @as(std.math.Log2Int(u128), @intCast(shift_count));
    const quotient = numerator / scaled_denominator;
    const remainder = numerator % scaled_denominator;
    const half = scaled_denominator >> 1;
    const above = remainder > half;
    const tie = (scaled_denominator & 1) == 0 and remainder == half;
    return .{ .value = quotient, .up = above or (tie and (quotient & 1) != 0) };
}

/// Round an exact quotient once. The exponent is the power of two outside the
/// integer numerator/denominator ratio.
pub fn packRatio(
    comptime destination: Format,
    negative: bool,
    numerator: u128,
    denominator: u128,
    exponent: i32,
) destination.Bits() {
    std.debug.assert(numerator != 0);
    std.debug.assert(denominator != 0);
    const top = ratioTop(numerator, denominator);
    const result_top = exponent + top;
    const minimum = Format.minNormalExponent(destination);
    if (result_top < minimum) {
        const scale = minimum - @as(i32, destination.fraction_bits);
        const rounded_result = roundRatioAtScale(numerator, denominator, exponent - scale);
        const rounded = rounded_result.value + @as(u128, @intFromBool(rounded_result.up));
        return packRounded(destination, negative, rounded, minimum, true);
    }
    const scale = result_top - @as(i32, destination.fraction_bits);
    const rounded_result = roundRatioAtScale(numerator, denominator, exponent - scale);
    const rounded = rounded_result.value + @as(u128, @intFromBool(rounded_result.up));
    return packRounded(destination, negative, rounded, result_top, false);
}

fn compareShifted(left: u128, shift: u32, right: u128) i8 {
    if (left == 0) return if (right == 0) 0 else -1;
    if (right == 0) return 1;
    const left_width = bitLength(left) + shift;
    const right_width = bitLength(right);
    if (left_width < right_width) return -1;
    if (left_width > right_width) return 1;
    if (left_width > 128) return 1;
    const shifted = left << @as(std.math.Log2Int(u128), @intCast(shift));
    return if (shifted < right) -1 else if (shifted > right) 1 else 0;
}

fn integerSqrt(value: u64) u64 {
    @setFloatMode(.strict);
    // Conversion is exact below 2^53. Above it, the input error is at most
    // 2^10 and both roots are at least 2^26, giving a root error <= 2^-17.
    // Rounding the f64 square root adds at most 2^-21 (root <= 2^32).
    // Thus the truncated seed differs from floor(sqrt(value)) by at most one.
    const seed: u64 = @intFromFloat(@sqrt(@as(f64, @floatFromInt(value))));
    const square = @as(u128, seed) * seed;
    if (square > value) return seed - 1;
    if (@as(u128, value) - square > (@as(u128, seed) << 1)) return seed + 1;
    return seed;
}

const SqrtRoundResult = struct {
    value: u128,
    up: bool,
};

fn sqrtScaledPositive(base: u64, shift_half: u32) SqrtRoundResult {
    const root = integerSqrt(base);
    const square = @as(u128, root) * root;
    const remainder = @as(u128, base) - square;
    const scaled_root = @as(u128, root) << @as(std.math.Log2Int(u128), @intCast(shift_half));
    if (remainder == 0) return .{ .value = scaled_root, .up = false };

    const left_shift = shift_half * 2 + 2;
    const right = @as(u128, root) << @as(std.math.Log2Int(u128), @intCast(shift_half + 2));
    const compare = if (bitLength(remainder) + left_shift > 128)
        @as(i8, 1)
    else
        @as(
            i8,
            if ((remainder << @as(std.math.Log2Int(u128), @intCast(left_shift))) > right + 1)
                1
            else if ((remainder << @as(std.math.Log2Int(u128), @intCast(left_shift))) < right + 1)
                -1
            else
                0,
        );
    return .{ .value = scaled_root, .up = compare > 0 or (compare == 0 and (scaled_root & 1) != 0) };
}

fn sqrtScaledNegative(base: u64, denominator_shift: u32) SqrtRoundResult {
    std.debug.assert((denominator_shift & 1) == 0);
    const root = integerSqrt(base);
    // The caller uses denominator_shift = 2*k. Shifting the exact integer
    // root by k is already the exact floor of the scaled square root.
    const quotient: u128 = if (denominator_shift >= 64)
        0
    else
        root >> @as(std.math.Log2Int(u64), @intCast(denominator_shift / 2));
    const doubled = (quotient * 2) + 1;
    const half_compare = compareShifted(doubled * doubled, denominator_shift, @as(u128, base) * 4);
    return .{ .value = quotient, .up = half_compare < 0 or (half_compare == 0 and (quotient & 1) != 0) };
}

/// Round the exact square root of `significand * 2^exponent` once.
pub fn packSqrt(comptime destination: Format, significand: u64, exponent: i32) destination.Bits() {
    std.debug.assert(significand != 0);
    const source_top = exponent + @as(i32, @intCast(bitLength(significand) - 1));
    const result_top = @divFloor(source_top, 2);
    const minimum = Format.minNormalExponent(destination);
    const scale = if (result_top < minimum)
        minimum - @as(i32, destination.fraction_bits)
    else
        result_top - @as(i32, destination.fraction_bits);
    const square_shift = exponent - 2 * scale;

    // SAFETY: Both square-shift branches assign `rounded_result` before its fields are read.
    var rounded_result: SqrtRoundResult = undefined;
    if (square_shift >= 0) {
        const shift: u32 = @intCast(square_shift);
        const odd = shift & 1 != 0;
        const base: u64 = if (odd) significand * 2 else significand;
        rounded_result = sqrtScaledPositive(base, (if (odd) shift - 1 else shift) / 2);
    } else {
        // `square_shift < 0` here, so its magnitude is the denominator shift and
        // `@abs` takes it without the `-x` overflow at `minInt(i32)`. `exponent`
        // is a decoded binary exponent, so it is at least `-1074`: the `minimum`
        // branch above then gives `exponent - 2 * scale >= -1074 + 2 * (14 + 10)`,
        // and the `result_top` branch gives `exponent - 2 * scale <= 1 + 2 * 52`
        // because `exponent - 2 * result_top` is `-bit` or `1 - bit`. So the
        // magnitude is at most 1026, far inside `u32`, and `sqrtScaledNegative`
        // already folds any shift of 64 or more to zero.
        const shift: u32 = @abs(square_shift);
        const odd = shift & 1 != 0;
        const base: u64 = if (odd) significand * 2 else significand;
        const denominator_shift = if (odd) shift + 1 else shift;
        rounded_result = sqrtScaledNegative(base, denominator_shift);
    }
    const rounded = rounded_result.value + @as(u128, @intFromBool(rounded_result.up));
    return packRounded(
        destination,
        false,
        rounded,
        if (result_top < minimum) minimum else result_top,
        result_top < minimum,
    );
}

pub fn finiteDyadic(comptime source: Format, value: Decoded) Dyadic {
    std.debug.assert(value.class == .finite or value.class == .zero);
    return exactDyadic(source, value);
}

test "integer square root is exact at range ends and adjacent squares" {
    const precision_root: u64 = (@as(u64, 1) << 27) - 1;
    const precision_square = precision_root * precision_root;
    const maximum_root: u64 = (@as(u64, 1) << 32) - 1;
    const maximum_square = maximum_root * maximum_root;
    const cases = [_]struct { value: u64, root: u64 }{
        .{ .value = 0, .root = 0 },
        .{ .value = 1, .root = 1 },
        .{ .value = 3, .root = 1 },
        .{ .value = 4, .root = 2 },
        .{ .value = 8, .root = 2 },
        .{ .value = 9, .root = 3 },
        .{ .value = precision_square - 1, .root = precision_root - 1 },
        .{ .value = precision_square, .root = precision_root },
        .{ .value = precision_square + 1, .root = precision_root },
        .{ .value = maximum_square - 1, .root = maximum_root - 1 },
        .{ .value = maximum_square, .root = maximum_root },
        .{ .value = maximum_square + 1, .root = maximum_root },
        .{ .value = ~@as(u64, 0), .root = maximum_root },
    };
    for (cases) |case| {
        try std.testing.expectEqual(case.root, integerSqrt(case.value));
    }
}
