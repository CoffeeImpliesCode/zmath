//! Narrowing operations: the exact source-format result, rounded once to the
//! destination format.
//!
//! Every kernel below forms the exact mathematical result over the source
//! operands and rounds it once, round-to-nearest ties-to-even, to the
//! destination format. None of them is a source-format operation followed by a
//! narrowing conversion: `fadd(1 + 2^-25, -1)` is `2^-25`, where a binary32 add
//! of the two narrowed operands is `+0`. The fused forms round the whole of
//! `x * y + z` once, exactly as the two-operand forms round their single
//! expression.
//!
//! The result is a bare float and the only channel. Nothing is reported
//! through `errno`, through a floating-point exception flag, or through any
//! retained status, and no kernel reads or changes the process floating-point
//! environment: the rounding direction is fixed here, not selected by the
//! ambient mode.
//!
//! A NaN operand is converted to the destination width, keeping its sign and
//! the top payload bits that fit in the destination and setting the quiet bit,
//! so a signaling NaN is quieted rather than discarded. When more than one
//! operand is NaN, the leftmost one is the one converted.
//!
//! An operation that is invalid for the destination gives the positive
//! canonical quiet NaN of that format, `0x7e00` for binary16 and `0x7fc00000`
//! for binary32, and each declaration below names the input shape that is
//! invalid for its own operation. A result that rounds above the largest
//! finite destination value is an infinity of the same sign instead. No square
//! root does that, and a square root is never negative.
//!
//! An exact zero result takes the IEC sign. A sum that cancels is `+0` unless
//! both addends carry the same sign, which keeps it, and in a subtraction the
//! second addend is `-y`. A product and a quotient take the sign of the
//! product, a fused multiply-add keeps a sign only when the product and the
//! addend share it, and a square root of a signed zero keeps the input sign.
const common = @import("scalar_narrow_common.zig");

fn firstNaN2(
    comptime source: common.Format,
    comptime destination: common.Format,
    first_bits: u64,
    second_bits: u64,
) ?destination.Bits() {
    const first = common.decode(source, first_bits);
    if (first.class == .nan) return common.convertedNaN(source, destination, first_bits);
    const second = common.decode(source, second_bits);
    if (second.class == .nan) return common.convertedNaN(source, destination, second_bits);
    return null;
}

fn firstNaN3(
    comptime source: common.Format,
    comptime destination: common.Format,
    first_bits: u64,
    second_bits: u64,
    third_bits: u64,
) ?destination.Bits() {
    if (firstNaN2(source, destination, first_bits, second_bits)) |result| return result;
    const third = common.decode(source, third_bits);
    if (third.class == .nan) return common.convertedNaN(source, destination, third_bits);
    return null;
}

fn addBits(
    comptime source: common.Format,
    comptime destination: common.Format,
    first_bits: u64,
    second_bits: u64,
    subtract: bool,
) destination.Bits() {
    if (firstNaN2(source, destination, first_bits, second_bits)) |result| return result;
    const first = common.decode(source, first_bits);
    var second = common.decode(source, second_bits);
    if (subtract) second.negative = !second.negative;

    if (first.class == .infinity) {
        if (second.class == .infinity and first.negative != second.negative) return common.invalid(destination);
        return common.fromBits(destination, common.Format.infinityBits(destination) |
            (if (first.negative) common.Format.signMask(destination) else 0));
    }
    if (second.class == .infinity) {
        return common.fromBits(destination, common.Format.infinityBits(destination) |
            (if (second.negative) common.Format.signMask(destination) else 0));
    }

    var result = common.addDyadic(common.finiteDyadic(source, first), common.finiteDyadic(source, second));
    if (result.magnitude == 0) {
        // Round-to-nearest gives +0 for cancellation. Equal signed zeros keep
        // their sign, which is the IEC result for exact zero inputs.
        result.negative = first.negative and second.negative;
    }
    return common.packDyadic(destination, result);
}

fn multiplyBits(
    comptime source: common.Format,
    comptime destination: common.Format,
    first_bits: u64,
    second_bits: u64,
) destination.Bits() {
    if (firstNaN2(source, destination, first_bits, second_bits)) |result| return result;
    const first = common.decode(source, first_bits);
    const second = common.decode(source, second_bits);
    const negative = first.negative != second.negative;

    if ((first.class == .infinity and second.class == .zero) or
        (first.class == .zero and second.class == .infinity)) return common.invalid(destination);
    if (first.class == .infinity or second.class == .infinity) {
        return common.fromBits(destination, common.Format.infinityBits(destination) |
            (if (negative) common.Format.signMask(destination) else 0));
    }
    if (first.class == .zero or second.class == .zero) return common.zero(destination, negative);

    return common.packDyadic(destination, common.multiplyDyadic(
        common.finiteDyadic(source, first),
        common.finiteDyadic(source, second),
    ));
}

fn divideBits(
    comptime source: common.Format,
    comptime destination: common.Format,
    numerator_bits: u64,
    denominator_bits: u64,
) destination.Bits() {
    if (firstNaN2(source, destination, numerator_bits, denominator_bits)) |result| return result;
    const numerator = common.decode(source, numerator_bits);
    const denominator = common.decode(source, denominator_bits);
    const negative = numerator.negative != denominator.negative;

    if ((numerator.class == .zero and denominator.class == .zero) or
        (numerator.class == .infinity and denominator.class == .infinity)) return common.invalid(destination);
    if (numerator.class == .infinity) {
        return common.fromBits(destination, common.Format.infinityBits(destination) |
            (if (negative) common.Format.signMask(destination) else 0));
    }
    if (denominator.class == .infinity) return common.zero(destination, negative);
    if (denominator.class == .zero) {
        return common.fromBits(destination, common.Format.infinityBits(destination) |
            (if (negative) common.Format.signMask(destination) else 0));
    }
    if (numerator.class == .zero) return common.zero(destination, negative);

    return common.packRatio(
        destination,
        negative,
        numerator.significand,
        denominator.significand,
        numerator.exponent - denominator.exponent,
    );
}

fn fmaBits(
    comptime source: common.Format,
    comptime destination: common.Format,
    first_bits: u64,
    second_bits: u64,
    addend_bits: u64,
) destination.Bits() {
    if (firstNaN3(source, destination, first_bits, second_bits, addend_bits)) |result| return result;
    const first = common.decode(source, first_bits);
    const second = common.decode(source, second_bits);
    const addend = common.decode(source, addend_bits);
    const product_negative = first.negative != second.negative;

    if ((first.class == .infinity and second.class == .zero) or
        (first.class == .zero and second.class == .infinity)) return common.invalid(destination);
    if (first.class == .infinity or second.class == .infinity) {
        if (addend.class == .infinity and product_negative != addend.negative) return common.invalid(destination);
        return common.fromBits(destination, common.Format.infinityBits(destination) |
            (if (product_negative) common.Format.signMask(destination) else 0));
    }
    if (addend.class == .infinity) {
        return common.fromBits(destination, common.Format.infinityBits(destination) |
            (if (addend.negative) common.Format.signMask(destination) else 0));
    }

    const product = common.multiplyDyadic(
        common.finiteDyadic(source, first),
        common.finiteDyadic(source, second),
    );
    var result = common.addDyadic(product, common.finiteDyadic(source, addend));
    if (result.magnitude == 0) {
        result.negative = product_negative and addend.negative;
    }
    return common.packDyadic(destination, result);
}

fn sqrtBits(
    comptime source: common.Format,
    comptime destination: common.Format,
    input_bits: u64,
) destination.Bits() {
    const input = common.decode(source, input_bits);
    if (input.class == .nan) return common.convertedNaN(source, destination, input_bits);
    if (input.class == .infinity) {
        if (input.negative) return common.invalid(destination);
        return common.fromBits(destination, common.Format.infinityBits(destination));
    }
    if (input.class == .zero) return common.zero(destination, input.negative);
    if (input.negative) return common.invalid(destination);
    return common.packSqrt(destination, input.significand, input.exponent);
}

fn sourceBits(comptime source: common.Format, value: anytype) u64 {
    return common.toBits(source, value);
}

/// The binary64 sum of `x` and `y`, rounded once to binary32 rather than
/// computed on the two narrowed operands. Opposite-signed infinities are the
/// invalid case.
pub fn fadd(x: f64, y: f64) f32 {
    return @bitCast(common.fromBits(
        common.binary32_format,
        addBits(
            common.binary64_format,
            common.binary32_format,
            sourceBits(common.binary64_format, x),
            sourceBits(common.binary64_format, y),
            false,
        ),
    ));
}

/// The binary64 difference of `x` and `y`, rounded once to binary32 rather
/// than computed on the two narrowed operands. An infinity minus itself is
/// the invalid case.
pub fn fsub(x: f64, y: f64) f32 {
    return @bitCast(common.fromBits(
        common.binary32_format,
        addBits(
            common.binary64_format,
            common.binary32_format,
            sourceBits(common.binary64_format, x),
            sourceBits(common.binary64_format, y),
            true,
        ),
    ));
}

/// The binary64 product of `x` and `y`, rounded once to binary32 rather than
/// computed on the two narrowed operands. An infinite operand together with a
/// zero operand is the invalid case.
pub fn fmul(x: f64, y: f64) f32 {
    return @bitCast(common.fromBits(
        common.binary32_format,
        multiplyBits(
            common.binary64_format,
            common.binary32_format,
            sourceBits(common.binary64_format, x),
            sourceBits(common.binary64_format, y),
        ),
    ));
}

/// The binary64 quotient of `x` by `y`, rounded once to binary32 rather than
/// computed on the two narrowed operands. `0/0` and `inf/inf` are the invalid
/// cases.
pub fn fdiv(x: f64, y: f64) f32 {
    return @bitCast(common.fromBits(
        common.binary32_format,
        divideBits(
            common.binary64_format,
            common.binary32_format,
            sourceBits(common.binary64_format, x),
            sourceBits(common.binary64_format, y),
        ),
    ));
}

/// The binary64 value of `x * y + z`, rounded once to binary32 over the whole
/// expression, so the product is not rounded before the addend joins it. An
/// infinite operand together with a zero operand, or an infinite product
/// opposed in sign by an infinite addend, is the invalid case.
pub fn ffma(x: f64, y: f64, z: f64) f32 {
    return @bitCast(common.fromBits(
        common.binary32_format,
        fmaBits(
            common.binary64_format,
            common.binary32_format,
            sourceBits(common.binary64_format, x),
            sourceBits(common.binary64_format, y),
            sourceBits(common.binary64_format, z),
        ),
    ));
}

/// The binary64 square root of `x`, rounded once to binary32 rather than taken
/// of a binary32 rounding of `x`. A negative `x` is the invalid case, and a
/// negative zero is not.
pub fn fsqrt(x: f64) f32 {
    return @bitCast(common.fromBits(
        common.binary32_format,
        sqrtBits(common.binary64_format, common.binary32_format, sourceBits(common.binary64_format, x)),
    ));
}

/// The binary32 sum of `x` and `y`, rounded once to binary16. Opposite-signed
/// infinities are the invalid case.
pub fn f16addf32(x: f32, y: f32) f16 {
    return @bitCast(common.fromBits(
        common.binary16_format,
        addBits(
            common.binary32_format,
            common.binary16_format,
            sourceBits(common.binary32_format, x),
            sourceBits(common.binary32_format, y),
            false,
        ),
    ));
}

/// The binary32 difference of `x` and `y`, rounded once to binary16. An
/// infinity minus itself is the invalid case.
pub fn f16subf32(x: f32, y: f32) f16 {
    return @bitCast(common.fromBits(
        common.binary16_format,
        addBits(
            common.binary32_format,
            common.binary16_format,
            sourceBits(common.binary32_format, x),
            sourceBits(common.binary32_format, y),
            true,
        ),
    ));
}

/// The binary32 product of `x` and `y`, rounded once to binary16. An infinite
/// operand together with a zero operand is the invalid case.
pub fn f16mulf32(x: f32, y: f32) f16 {
    return @bitCast(common.fromBits(
        common.binary16_format,
        multiplyBits(
            common.binary32_format,
            common.binary16_format,
            sourceBits(common.binary32_format, x),
            sourceBits(common.binary32_format, y),
        ),
    ));
}

/// The binary32 quotient of `x` by `y`, rounded once to binary16. `0/0` and
/// `inf/inf` are the invalid cases.
pub fn f16divf32(x: f32, y: f32) f16 {
    return @bitCast(common.fromBits(
        common.binary16_format,
        divideBits(
            common.binary32_format,
            common.binary16_format,
            sourceBits(common.binary32_format, x),
            sourceBits(common.binary32_format, y),
        ),
    ));
}

/// The binary32 value of `x * y + z`, rounded once to binary16 over the whole
/// expression, so the product is not rounded before the addend joins it. An
/// infinite operand together with a zero operand, or an infinite product
/// opposed in sign by an infinite addend, is the invalid case.
pub fn f16fmaf32(x: f32, y: f32, z: f32) f16 {
    return @bitCast(common.fromBits(
        common.binary16_format,
        fmaBits(
            common.binary32_format,
            common.binary16_format,
            sourceBits(common.binary32_format, x),
            sourceBits(common.binary32_format, y),
            sourceBits(common.binary32_format, z),
        ),
    ));
}

/// The binary32 square root of `x`, rounded once to binary16. A negative `x`
/// is the invalid case, and a negative zero is not.
pub fn f16sqrtf32(x: f32) f16 {
    return @bitCast(common.fromBits(
        common.binary16_format,
        sqrtBits(common.binary32_format, common.binary16_format, sourceBits(common.binary32_format, x)),
    ));
}

/// The binary64 sum of `x` and `y`, rounded once to binary16. Opposite-signed
/// infinities are the invalid case.
pub fn f16addf64(x: f64, y: f64) f16 {
    return @bitCast(common.fromBits(
        common.binary16_format,
        addBits(
            common.binary64_format,
            common.binary16_format,
            sourceBits(common.binary64_format, x),
            sourceBits(common.binary64_format, y),
            false,
        ),
    ));
}

/// The binary64 difference of `x` and `y`, rounded once to binary16. An
/// infinity minus itself is the invalid case.
pub fn f16subf64(x: f64, y: f64) f16 {
    return @bitCast(common.fromBits(
        common.binary16_format,
        addBits(
            common.binary64_format,
            common.binary16_format,
            sourceBits(common.binary64_format, x),
            sourceBits(common.binary64_format, y),
            true,
        ),
    ));
}

/// The binary64 product of `x` and `y`, rounded once to binary16. An infinite
/// operand together with a zero operand is the invalid case.
pub fn f16mulf64(x: f64, y: f64) f16 {
    return @bitCast(common.fromBits(
        common.binary16_format,
        multiplyBits(
            common.binary64_format,
            common.binary16_format,
            sourceBits(common.binary64_format, x),
            sourceBits(common.binary64_format, y),
        ),
    ));
}

/// The binary64 quotient of `x` by `y`, rounded once to binary16. `0/0` and
/// `inf/inf` are the invalid cases.
pub fn f16divf64(x: f64, y: f64) f16 {
    return @bitCast(common.fromBits(
        common.binary16_format,
        divideBits(
            common.binary64_format,
            common.binary16_format,
            sourceBits(common.binary64_format, x),
            sourceBits(common.binary64_format, y),
        ),
    ));
}

/// The binary64 value of `x * y + z`, rounded once to binary16 over the whole
/// expression, so the product is not rounded before the addend joins it. An
/// infinite operand together with a zero operand, or an infinite product
/// opposed in sign by an infinite addend, is the invalid case.
pub fn f16fmaf64(x: f64, y: f64, z: f64) f16 {
    return @bitCast(common.fromBits(
        common.binary16_format,
        fmaBits(
            common.binary64_format,
            common.binary16_format,
            sourceBits(common.binary64_format, x),
            sourceBits(common.binary64_format, y),
            sourceBits(common.binary64_format, z),
        ),
    ));
}

/// The binary64 square root of `x`, rounded once to binary16. A negative `x`
/// is the invalid case, and a negative zero is not.
pub fn f16sqrtf64(x: f64) f16 {
    return @bitCast(common.fromBits(
        common.binary16_format,
        sqrtBits(common.binary64_format, common.binary16_format, sourceBits(common.binary64_format, x)),
    ));
}

/// The binary64 sum of `x` and `y`, rounded once to binary32. Opposite-signed
/// infinities are the invalid case.
pub fn f32addf64(x: f64, y: f64) f32 {
    return @bitCast(common.fromBits(
        common.binary32_format,
        addBits(
            common.binary64_format,
            common.binary32_format,
            sourceBits(common.binary64_format, x),
            sourceBits(common.binary64_format, y),
            false,
        ),
    ));
}

/// The binary64 difference of `x` and `y`, rounded once to binary32. An
/// infinity minus itself is the invalid case.
pub fn f32subf64(x: f64, y: f64) f32 {
    return @bitCast(common.fromBits(
        common.binary32_format,
        addBits(
            common.binary64_format,
            common.binary32_format,
            sourceBits(common.binary64_format, x),
            sourceBits(common.binary64_format, y),
            true,
        ),
    ));
}

/// The binary64 product of `x` and `y`, rounded once to binary32. An infinite
/// operand together with a zero operand is the invalid case.
pub fn f32mulf64(x: f64, y: f64) f32 {
    return @bitCast(common.fromBits(
        common.binary32_format,
        multiplyBits(
            common.binary64_format,
            common.binary32_format,
            sourceBits(common.binary64_format, x),
            sourceBits(common.binary64_format, y),
        ),
    ));
}

/// The binary64 quotient of `x` by `y`, rounded once to binary32. `0/0` and
/// `inf/inf` are the invalid cases.
pub fn f32divf64(x: f64, y: f64) f32 {
    return @bitCast(common.fromBits(
        common.binary32_format,
        divideBits(
            common.binary64_format,
            common.binary32_format,
            sourceBits(common.binary64_format, x),
            sourceBits(common.binary64_format, y),
        ),
    ));
}

/// The binary64 value of `x * y + z`, rounded once to binary32 over the whole
/// expression, so the product is not rounded before the addend joins it. An
/// infinite operand together with a zero operand, or an infinite product
/// opposed in sign by an infinite addend, is the invalid case.
pub fn f32fmaf64(x: f64, y: f64, z: f64) f32 {
    return @bitCast(common.fromBits(
        common.binary32_format,
        fmaBits(
            common.binary64_format,
            common.binary32_format,
            sourceBits(common.binary64_format, x),
            sourceBits(common.binary64_format, y),
            sourceBits(common.binary64_format, z),
        ),
    ));
}

/// The binary64 square root of `x`, rounded once to binary32. A negative `x`
/// is the invalid case, and a negative zero is not.
pub fn f32sqrtf64(x: f64) f32 {
    return @bitCast(common.fromBits(
        common.binary32_format,
        sqrtBits(common.binary64_format, common.binary32_format, sourceBits(common.binary64_format, x)),
    ));
}
