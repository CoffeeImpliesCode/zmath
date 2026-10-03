// C23 integer-power, restricted-power, compound-power, and integer-root families.
//
// Exceptional results are assembled from target-format bits.  General finite
// pown, rootn, and compoundn values in the accepted degree range are delegated
// to scalar_power_exact, which retains input dyadics as fixed-limb integers
// and rounds only once.  Binary64 single-operation identities use strict IEEE
// arithmetic directly, with a binary64 reciprocal routed through
// scalar_div.div64 so the subnormal range survives a software-float divide
// that drops it.

const common = @import("scalar_c23_common.zig");
const c23 = @import("scalar_c23_types.zig");
const direct_pow = @import("scalar_pow.zig");
const exact = @import("scalar_power_exact.zig");
const scalar_div = @import("scalar_div.zig");

fn signBits(comptime Format: type) Format.Bits {
    const negative_zero: Format.Float = -0.0;
    return Format.toBits(negative_zero);
}

fn infinityBits(comptime Format: type) Format.Bits {
    if (Format == common.H) return 0x7c00;
    if (Format == common.F) return 0x7f80_0000;
    if (Format == common.D) return 0x7ff0_0000_0000_0000;
    @compileError("unsupported C23 scalar format");
}

fn oneBits(comptime Format: type) Format.Bits {
    if (Format == common.H) return 0x3c00;
    if (Format == common.F) return 0x3f80_0000;
    if (Format == common.D) return 0x3ff0_0000_0000_0000;
    @compileError("unsupported C23 scalar format");
}

fn tinyCompoundLimitBits(comptime Format: type) Format.Bits {
    // Exactly 2^(-p-8), where p is the destination precision.
    if (Format == common.H) return 0x0020;
    if (Format == common.F) return 0x2f80_0000;
    if (Format == common.D) return 0x3c20_0000_0000_0000;
    @compileError("unsupported C23 scalar format");
}

fn signedZero(comptime Format: type, negative: bool) Format.Float {
    const bits: Format.Bits = if (negative) signBits(Format) else 0;
    return Format.fromBits(bits);
}

fn signedInfinity(comptime Format: type, negative: bool) Format.Float {
    const bits: Format.Bits = infinityBits(Format) |
        if (negative) signBits(Format) else 0;
    return Format.fromBits(bits);
}

fn absoluteValue(comptime Format: type, value: Format.Float) Format.Float {
    return Format.fromBits(Format.toBits(value) & ~signBits(Format));
}

fn copySign(comptime Format: type, magnitude: Format.Float, sign_source: Format.Bits) Format.Float {
    const magnitude_bits = Format.toBits(magnitude) & ~signBits(Format);
    return Format.fromBits(magnitude_bits | (sign_source & signBits(Format)));
}

fn quietNaN(comptime Format: type, value: Format.Float) Format.Float {
    return Format.fromBits(Format.quietNaNBits(Format.toBits(value)));
}

fn isOddExponent(exponent: i64) bool {
    return (@as(u64, @bitCast(exponent)) & 1) != 0;
}

fn inBound(exponent: i64) bool {
    return exponent >= -exact.max_degree and exponent <= exact.max_degree;
}

fn success(comptime Format: type, value: Format.Float) c23.PowerResult(Format.Float) {
    return .{ .value = value, .ok = true };
}

fn rejected(comptime Format: type) c23.PowerResult(Format.Float) {
    return .{ .value = Format.canonicalNaN(), .ok = false };
}

fn directPow(comptime Format: type, base: Format.Float, exponent: Format.Float) Format.Float {
    if (Format == common.H) return direct_pow.powh(base, exponent);
    if (Format == common.F) return direct_pow.powf(base, exponent);
    if (Format == common.D) return direct_pow.powd(base, exponent);
    @compileError("unsupported C23 scalar format");
}

fn pownCore(
    comptime Format: type,
    base: Format.Float,
    exponent: i64,
) c23.PowerResult(Format.Float) {
    @setFloatMode(.strict);
    if (!inBound(exponent)) return rejected(Format);

    const base_bits = Format.toBits(base);
    const category = Format.classifyBits(base_bits);
    if (category == .signaling_nan) return success(Format, quietNaN(Format, base));
    if (exponent == 0) return success(Format, 1.0);
    if (category == .quiet_nan) return success(Format, quietNaN(Format, base));

    const base_negative = (base_bits & signBits(Format)) != 0;
    const negative_result = base_negative and isOddExponent(exponent);
    if (category == .zero) {
        if (exponent > 0) return success(Format, signedZero(Format, negative_result));
        return success(Format, signedInfinity(Format, negative_result));
    }
    if (category == .infinity) {
        if (exponent > 0) return success(Format, signedInfinity(Format, negative_result));
        return success(Format, signedZero(Format, negative_result));
    }

    if (exponent == 1) return success(Format, base);
    if ((base_bits & ~signBits(Format)) == oneBits(Format)) {
        return success(
            Format,
            copySign(Format, 1.0, if (negative_result) signBits(Format) else 0),
        );
    }
    // IEEE multiplication rounds the exact product once in the destination format.
    if (exponent == 2) return success(Format, base * base);
    // IEEE division rounds the exact reciprocal once; div64 rebuilds the
    if (Format == common.D and exponent == -1) return success(Format, scalar_div.div64(1.0, base));

    const magnitude_result = exact.pownFinite(Format, absoluteValue(Format, base), exponent);
    return success(
        Format,
        copySign(Format, magnitude_result, if (negative_result) signBits(Format) else 0),
    );
}

fn powrCore(comptime Format: type, base: Format.Float, exponent: Format.Float) Format.Float {
    @setFloatMode(.strict);
    const base_bits = Format.toBits(base);
    const exponent_bits = Format.toBits(exponent);
    const base_is_nan = Format.isNaNBits(base_bits);
    const exponent_is_nan = Format.isNaNBits(exponent_bits);

    // powr is the exp(y * log(x)) operation.  Its domain excludes negative
    // bases, 0^0, 1^infinity, and infinity^0 before ordinary NaN handling.
    if (!base_is_nan) {
        if (base < 0.0) return Format.canonicalNaN();
        if (Format.isZeroBits(base_bits) and Format.isZeroBits(exponent_bits)) {
            return Format.canonicalNaN();
        }
        if (base == 1.0 and Format.isInfinityBits(exponent_bits)) {
            return Format.canonicalNaN();
        }
        if (Format.isInfinityBits(base_bits) and Format.isZeroBits(exponent_bits)) {
            return Format.canonicalNaN();
        }
    }
    if (base_is_nan) return quietNaN(Format, base);
    if (exponent_is_nan) return quietNaN(Format, exponent);

    return directPow(Format, absoluteValue(Format, base), exponent);
}

fn rootnCore(
    comptime Format: type,
    value: Format.Float,
    exponent: i64,
) c23.PowerResult(Format.Float) {
    @setFloatMode(.strict);
    if (!inBound(exponent)) return rejected(Format);

    const value_bits = Format.toBits(value);
    const category = Format.classifyBits(value_bits);
    if (exponent == 0) return success(Format, Format.canonicalNaN());
    if (category == .signaling_nan or category == .quiet_nan) {
        return success(Format, quietNaN(Format, value));
    }

    const value_negative = (value_bits & signBits(Format)) != 0;
    const odd_exponent = isOddExponent(exponent);
    // Negative zero is not less than zero, so it remains valid for even roots.
    if (value_negative and category != .zero and !odd_exponent) {
        return success(Format, Format.canonicalNaN());
    }
    if (category == .infinity) {
        if (exponent > 0) return success(Format, copySign(Format, signedInfinity(Format, false), value_bits));
        return success(Format, signedZero(Format, value_negative and odd_exponent));
    }
    if (category == .zero) {
        if (exponent > 0) return success(Format, signedZero(Format, value_negative and odd_exponent));
        return success(Format, signedInfinity(Format, value_negative and odd_exponent));
    }

    if (exponent == 1) return success(Format, value);
    if ((value_bits & ~signBits(Format)) == oneBits(Format)) return success(Format, value);
    if (Format == common.D) {
        // Each identity rounds once; degree -2 still needs exact comparisons.
        if (exponent == -1) return success(Format, scalar_div.div64(1.0, value));
        if (exponent == 2) return success(Format, @sqrt(value));
    }
    const magnitude_result = exact.rootFinite(Format, absoluteValue(Format, value), exponent);
    return success(Format, copySign(Format, magnitude_result, value_bits));
}

fn compoundnCore(
    comptime Format: type,
    value: Format.Float,
    exponent: i64,
) c23.PowerResult(Format.Float) {
    @setFloatMode(.strict);
    if (!inBound(exponent)) return rejected(Format);

    const value_bits = Format.toBits(value);
    const category = Format.classifyBits(value_bits);
    const value_is_nan = Format.isNaNBits(value_bits);

    if (category == .signaling_nan) return success(Format, quietNaN(Format, value));

    const magnitude_bits = value_bits & ~signBits(Format);
    const value_negative = (value_bits & signBits(Format)) != 0;
    if (!value_is_nan and value_negative and magnitude_bits > oneBits(Format)) {
        return success(Format, Format.canonicalNaN());
    }
    if (exponent == 0) return success(Format, 1.0);
    if (category == .quiet_nan) return success(Format, quietNaN(Format, value));

    if (value_negative and magnitude_bits == oneBits(Format)) {
        if (exponent > 0) return success(Format, signedZero(Format, false));
        return success(Format, signedInfinity(Format, false));
    }
    if (category == .infinity) {
        if (exponent > 0) return success(Format, signedInfinity(Format, false));
        return success(Format, signedZero(Format, false));
    }

    // IEEE addition rounds the exact degree-one compound result once.
    if (Format == common.D and exponent == 1) return success(Format, 1.0 + value);

    // Let t=2^(-p-8) and |x|<=t. For |n|<=64, both signs and
    // both exponent directions obey |(1+x)^n-1| <= (1-t)^-64-1.
    // Bernoulli gives (1-t)^64>=1-64t, hence the error is at most
    // 64t/(1-64t)<128t=2^(-p-1). This is below the narrower
    // midpoint immediately below one; the above-one midpoint is 2^-p.
    if (magnitude_bits <= tinyCompoundLimitBits(Format)) return success(Format, 1.0);

    return success(Format, exact.compoundFinite(Format, value, exponent));
}

/// Computes C23 pown for binary16 values.
/// Requires round-to-nearest, ties-to-even. Exponents outside [-64, 64]
/// return `.ok = false` with canonical NaN; all in-range cases set `.ok = true`,
/// including overflow, underflow, and quieted NaNs. Exponent zero returns one
/// except for a signaling-NaN base. Negative results require an odd exponent;
/// a negative exponent takes the reciprocal, including signed zero poles.
pub fn pownh(base: f16, exponent: i64) c23.PowerResult(f16) {
    @setFloatMode(.strict);
    return pownCore(common.H, base, exponent);
}

/// Computes C23 pown for binary32 values.
/// Requires round-to-nearest, ties-to-even. Exponents outside [-64, 64]
/// return `.ok = false` with canonical NaN; all in-range cases set `.ok = true`,
/// including overflow, underflow, and quieted NaNs. Exponent zero returns one
/// except for a signaling-NaN base. Negative results require an odd exponent;
/// a negative exponent takes the reciprocal, including signed zero poles.
pub fn pownf(base: f32, exponent: i64) c23.PowerResult(f32) {
    @setFloatMode(.strict);
    return pownCore(common.F, base, exponent);
}

/// Computes C23 pown for binary64 values.
/// Requires round-to-nearest, ties-to-even. Exponents outside [-64, 64]
/// return `.ok = false` with canonical NaN; all in-range cases set `.ok = true`,
/// including overflow, underflow, and quieted NaNs. Exponent zero returns one
/// except for a signaling-NaN base. Negative results require an odd exponent;
/// a negative exponent takes the reciprocal, including signed zero poles.
pub fn pownd(base: f64, exponent: i64) c23.PowerResult(f64) {
    @setFloatMode(.strict);
    return pownCore(common.D, base, exponent);
}

/// Computes C23 powr for binary16 values.
/// Requires round-to-nearest, ties-to-even. Negative bases, zero to zero,
/// one to an infinite exponent, and infinity to zero produce canonical NaN.
/// These domain cases precede NaN propagation; other NaNs are quieted, with
/// the base taking precedence. Negative zero is treated as positive zero.
pub fn powrh(base: f16, exponent: f16) f16 {
    @setFloatMode(.strict);
    return powrCore(common.H, base, exponent);
}

/// Computes C23 powr for binary32 values.
/// Requires round-to-nearest, ties-to-even. Negative bases, zero to zero,
/// one to an infinite exponent, and infinity to zero produce canonical NaN.
/// These domain cases precede NaN propagation; other NaNs are quieted, with
/// the base taking precedence. Negative zero is treated as positive zero.
pub fn powrf(base: f32, exponent: f32) f32 {
    @setFloatMode(.strict);
    return powrCore(common.F, base, exponent);
}

/// Computes C23 powr for binary64 values.
/// Requires round-to-nearest, ties-to-even. Negative bases, zero to zero,
/// one to an infinite exponent, and infinity to zero produce canonical NaN.
/// These domain cases precede NaN propagation; other NaNs are quieted, with
/// the base taking precedence. Negative zero is treated as positive zero.
pub fn powrd(base: f64, exponent: f64) f64 {
    @setFloatMode(.strict);
    return powrCore(common.D, base, exponent);
}

/// Computes C23 compoundn, (1 + x)^n, for binary16 values.
/// Requires round-to-nearest, ties-to-even. Exponents outside [-64, 64]
/// return `.ok = false` with canonical NaN. In-range domain failures (`value < -1`)
/// return canonical NaN with `.ok = true`; quieted NaNs and infinities also
/// keep `.ok = true`. Exponent zero returns one for `value >= -1` and quiet NaN,
/// but a signaling NaN is quieted first. At -1, negative exponents yield +infinity.
pub fn compoundnh(value: f16, exponent: i64) c23.PowerResult(f16) {
    @setFloatMode(.strict);
    return compoundnCore(common.H, value, exponent);
}

/// Computes C23 compoundn, (1 + x)^n, for binary32 values.
/// Requires round-to-nearest, ties-to-even. Exponents outside [-64, 64]
/// return `.ok = false` with canonical NaN. In-range domain failures (`value < -1`)
/// return canonical NaN with `.ok = true`; quieted NaNs and infinities also
/// keep `.ok = true`. Exponent zero returns one for `value >= -1` and quiet NaN,
/// but a signaling NaN is quieted first. At -1, negative exponents yield +infinity.
pub fn compoundnf(value: f32, exponent: i64) c23.PowerResult(f32) {
    @setFloatMode(.strict);
    return compoundnCore(common.F, value, exponent);
}

/// Computes C23 compoundn, (1 + x)^n, for binary64 values.
/// Requires round-to-nearest, ties-to-even. Exponents outside [-64, 64]
/// return `.ok = false` with canonical NaN. In-range domain failures (`value < -1`)
/// return canonical NaN with `.ok = true`; quieted NaNs and infinities also
/// keep `.ok = true`. Exponent zero returns one for `value >= -1` and quiet NaN,
/// but a signaling NaN is quieted first. At -1, negative exponents yield +infinity.
pub fn compoundnd(value: f64, exponent: i64) c23.PowerResult(f64) {
    @setFloatMode(.strict);
    return compoundnCore(common.D, value, exponent);
}

/// Computes C23 rootn for binary16 values.
/// Requires round-to-nearest, ties-to-even. Degrees outside [-64, 64] return
/// `.ok = false` with canonical NaN. Degree zero and a negative nonzero value
/// with an even degree produce canonical NaN with `.ok = true`; propagated
/// quieted NaNs, overflow, and underflow also keep `.ok = true`.
/// A negative degree takes the reciprocal; signed zero keeps its sign only
/// for an odd degree.
pub fn rootnh(value: f16, exponent: i64) c23.PowerResult(f16) {
    @setFloatMode(.strict);
    return rootnCore(common.H, value, exponent);
}

/// Computes C23 rootn for binary32 values.
/// Requires round-to-nearest, ties-to-even. Degrees outside [-64, 64] return
/// `.ok = false` with canonical NaN. Degree zero and a negative nonzero value
/// with an even degree produce canonical NaN with `.ok = true`; propagated
/// quieted NaNs, overflow, and underflow also keep `.ok = true`.
/// A negative degree takes the reciprocal; signed zero keeps its sign only
/// for an odd degree.
pub fn rootnf(value: f32, exponent: i64) c23.PowerResult(f32) {
    @setFloatMode(.strict);
    return rootnCore(common.F, value, exponent);
}

/// Computes C23 rootn for binary64 values.
/// Requires round-to-nearest, ties-to-even. Degrees outside [-64, 64] return
/// `.ok = false` with canonical NaN. Degree zero and a negative nonzero value
/// with an even degree produce canonical NaN with `.ok = true`; propagated
/// quieted NaNs, overflow, and underflow also keep `.ok = true`.
/// A negative degree takes the reciprocal; signed zero keeps its sign only
/// for an odd degree.
pub fn rootnd(value: f64, exponent: i64) c23.PowerResult(f64) {
    @setFloatMode(.strict);
    return rootnCore(common.D, value, exponent);
}
