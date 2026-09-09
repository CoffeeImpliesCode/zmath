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

pub fn fsqrt(x: f64) f32 {
    return @bitCast(common.fromBits(
        common.binary32_format,
        sqrtBits(common.binary64_format, common.binary32_format, sourceBits(common.binary64_format, x)),
    ));
}

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

pub fn f16sqrtf32(x: f32) f16 {
    return @bitCast(common.fromBits(
        common.binary16_format,
        sqrtBits(common.binary32_format, common.binary16_format, sourceBits(common.binary32_format, x)),
    ));
}

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

pub fn f16sqrtf64(x: f64) f16 {
    return @bitCast(common.fromBits(
        common.binary16_format,
        sqrtBits(common.binary64_format, common.binary16_format, sourceBits(common.binary64_format, x)),
    ));
}

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

pub fn f32sqrtf64(x: f64) f32 {
    return @bitCast(common.fromBits(
        common.binary32_format,
        sqrtBits(common.binary64_format, common.binary32_format, sourceBits(common.binary64_format, x)),
    ));
}
