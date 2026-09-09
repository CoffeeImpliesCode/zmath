const ieee = @import("ieee.zig");

fn expectBits(expected: u16, actual: f16) !void {
    if (ieee.toBitsh(actual) != expected) return error.TestExpectedEqual;
}

test "H bit conversions preserve representative IEEE encodings" {
    const encodings = [_]u16{
        0x0000, // Positive zero.
        0x8000, // Negative zero.
        0x0001, // Minimum positive subnormal.
        0x83ff, // Maximum negative subnormal.
        0x0400, // Minimum positive normal.
        0xbc00, // Negative one.
        0x7bff, // Maximum positive finite.
        0xfbff, // Maximum negative finite.
        0x7c00, // Positive infinity.
        0xfc00, // Negative infinity.
        0x7e55, // Positive quiet NaN with payload.
        0xfd55, // Negative signaling NaN with payload.
    };

    for (encodings) |bits| {
        if (ieee.toBitsh(ieee.fromBitsh(bits)) != bits) return error.TestExpectedEqual;
    }
}

test "nextUph and nextDownh step across both signed zeros" {
    const positive_zero = ieee.fromBitsh(0x0000);
    const negative_zero = ieee.fromBitsh(0x8000);

    try expectBits(0x0001, ieee.nextUph(positive_zero));
    try expectBits(0x0001, ieee.nextUph(negative_zero));
    try expectBits(0x8001, ieee.nextDownh(positive_zero));
    try expectBits(0x8001, ieee.nextDownh(negative_zero));
}

test "nextUph and nextDownh handle infinity boundaries" {
    const positive_infinity = ieee.fromBitsh(0x7c00);
    const negative_infinity = ieee.fromBitsh(0xfc00);

    try expectBits(0x7c00, ieee.nextUph(positive_infinity));
    try expectBits(0x7bff, ieee.nextDownh(positive_infinity));
    try expectBits(0xfbff, ieee.nextUph(negative_infinity));
    try expectBits(0xfc00, ieee.nextDownh(negative_infinity));
}

test "nextUph and nextDownh quiet signaling NaNs" {
    const positive_signaling_nan = ieee.fromBitsh(0x7d55);
    const negative_signaling_nan = ieee.fromBitsh(0xfd55);

    try expectBits(0x7f55, ieee.nextUph(positive_signaling_nan));
    try expectBits(0x7f55, ieee.nextDownh(positive_signaling_nan));
    try expectBits(0xff55, ieee.nextUph(negative_signaling_nan));
    try expectBits(0xff55, ieee.nextDownh(negative_signaling_nan));
}

test "scaleBh rounds positive and negative overflow in every mode" {
    const positive_maximum = ieee.fromBitsh(0x7bff);
    const negative_maximum = ieee.fromBitsh(0xfbff);
    const cases = [_]struct {
        mode: ieee.RoundMode,
        positive_bits: u16,
        negative_bits: u16,
    }{
        .{ .mode = .nearest_even, .positive_bits = 0x7c00, .negative_bits = 0xfc00 },
        .{ .mode = .toward_zero, .positive_bits = 0x7bff, .negative_bits = 0xfbff },
        .{ .mode = .toward_positive, .positive_bits = 0x7c00, .negative_bits = 0xfbff },
        .{ .mode = .toward_negative, .positive_bits = 0x7bff, .negative_bits = 0xfc00 },
    };

    for (cases) |case| {
        try expectBits(case.positive_bits, ieee.scaleBh(case.mode, positive_maximum, 1));
        try expectBits(case.negative_bits, ieee.scaleBh(case.mode, negative_maximum, 1));
    }
}

test "scaleBh rounds positive and negative underflow in every mode" {
    const positive_minimum = ieee.fromBitsh(0x0001);
    const negative_minimum = ieee.fromBitsh(0x8001);
    const cases = [_]struct {
        mode: ieee.RoundMode,
        positive_bits: u16,
        negative_bits: u16,
    }{
        .{ .mode = .nearest_even, .positive_bits = 0x0000, .negative_bits = 0x8000 },
        .{ .mode = .toward_zero, .positive_bits = 0x0000, .negative_bits = 0x8000 },
        .{ .mode = .toward_positive, .positive_bits = 0x0001, .negative_bits = 0x8000 },
        .{ .mode = .toward_negative, .positive_bits = 0x0000, .negative_bits = 0x8001 },
    };

    for (cases) |case| {
        try expectBits(case.positive_bits, ieee.scaleBh(case.mode, positive_minimum, -1));
        try expectBits(case.negative_bits, ieee.scaleBh(case.mode, negative_minimum, -1));
    }
}

test "scaleBh promotes nonzero subnormals exactly in every mode" {
    const modes = [_]ieee.RoundMode{
        .nearest_even,
        .toward_zero,
        .toward_positive,
        .toward_negative,
    };
    for (modes) |mode| {
        try expectBits(0x0400, ieee.scaleBh(mode, ieee.fromBitsh(0x0001), 10));
        try expectBits(0x8400, ieee.scaleBh(mode, ieee.fromBitsh(0x8001), 10));
        try expectBits(0x0400, ieee.scaleBh(mode, ieee.fromBitsh(0x0200), 1));
        try expectBits(0x07fe, ieee.scaleBh(mode, ieee.fromBitsh(0x03ff), 1));
        try expectBits(0x03ff, ieee.scaleBh(mode, ieee.fromBitsh(0x03ff), 0));
    }
}

test "scaleBh rounds subnormal ties by parity in every mode" {
    const cases = [_]struct {
        mode: ieee.RoundMode,
        even_positive_bits: u16,
        even_negative_bits: u16,
        odd_positive_bits: u16,
        odd_negative_bits: u16,
    }{
        .{
            .mode = .nearest_even,
            .even_positive_bits = 0x0200,
            .even_negative_bits = 0x8200,
            .odd_positive_bits = 0x0202,
            .odd_negative_bits = 0x8202,
        },
        .{
            .mode = .toward_zero,
            .even_positive_bits = 0x0200,
            .even_negative_bits = 0x8200,
            .odd_positive_bits = 0x0201,
            .odd_negative_bits = 0x8201,
        },
        .{
            .mode = .toward_positive,
            .even_positive_bits = 0x0201,
            .even_negative_bits = 0x8200,
            .odd_positive_bits = 0x0202,
            .odd_negative_bits = 0x8201,
        },
        .{
            .mode = .toward_negative,
            .even_positive_bits = 0x0200,
            .even_negative_bits = 0x8201,
            .odd_positive_bits = 0x0201,
            .odd_negative_bits = 0x8202,
        },
    };

    for (cases) |case| {
        try expectBits(
            case.even_positive_bits,
            ieee.scaleBh(case.mode, ieee.fromBitsh(0x0401), -1),
        );
        try expectBits(
            case.even_negative_bits,
            ieee.scaleBh(case.mode, ieee.fromBitsh(0x8401), -1),
        );
        try expectBits(
            case.odd_positive_bits,
            ieee.scaleBh(case.mode, ieee.fromBitsh(0x0403), -1),
        );
        try expectBits(
            case.odd_negative_bits,
            ieee.scaleBh(case.mode, ieee.fromBitsh(0x8403), -1),
        );
    }
}

test "scaleBh rounds a subnormal carry into the minimum normal in every mode" {
    const cases = [_]struct {
        mode: ieee.RoundMode,
        positive_bits: u16,
        negative_bits: u16,
    }{
        .{ .mode = .nearest_even, .positive_bits = 0x0400, .negative_bits = 0x8400 },
        .{ .mode = .toward_zero, .positive_bits = 0x03ff, .negative_bits = 0x83ff },
        .{ .mode = .toward_positive, .positive_bits = 0x0400, .negative_bits = 0x83ff },
        .{ .mode = .toward_negative, .positive_bits = 0x03ff, .negative_bits = 0x8400 },
    };

    for (cases) |case| {
        try expectBits(case.positive_bits, ieee.scaleBh(case.mode, ieee.fromBitsh(0x07ff), -1));
        try expectBits(case.negative_bits, ieee.scaleBh(case.mode, ieee.fromBitsh(0x87ff), -1));
    }
}

test "scaleBh keeps special values and saturates extreme exponents in every mode" {
    const cases = [_]struct {
        mode: ieee.RoundMode,
        positive_overflow_bits: u16,
        negative_overflow_bits: u16,
        positive_underflow_bits: u16,
        negative_underflow_bits: u16,
    }{
        .{
            .mode = .nearest_even,
            .positive_overflow_bits = 0x7c00,
            .negative_overflow_bits = 0xfc00,
            .positive_underflow_bits = 0x0000,
            .negative_underflow_bits = 0x8000,
        },
        .{
            .mode = .toward_zero,
            .positive_overflow_bits = 0x7bff,
            .negative_overflow_bits = 0xfbff,
            .positive_underflow_bits = 0x0000,
            .negative_underflow_bits = 0x8000,
        },
        .{
            .mode = .toward_positive,
            .positive_overflow_bits = 0x7c00,
            .negative_overflow_bits = 0xfbff,
            .positive_underflow_bits = 0x0001,
            .negative_underflow_bits = 0x8000,
        },
        .{
            .mode = .toward_negative,
            .positive_overflow_bits = 0x7bff,
            .negative_overflow_bits = 0xfc00,
            .positive_underflow_bits = 0x0000,
            .negative_underflow_bits = 0x8001,
        },
    };

    for (cases) |case| {
        try expectBits(
            case.positive_overflow_bits,
            ieee.scaleBh(case.mode, ieee.fromBitsh(0x0001), 2_147_483_647),
        );
        try expectBits(
            case.negative_overflow_bits,
            ieee.scaleBh(case.mode, ieee.fromBitsh(0x8001), 2_147_483_647),
        );
        try expectBits(
            case.positive_underflow_bits,
            ieee.scaleBh(case.mode, ieee.fromBitsh(0x0001), -2_147_483_648),
        );
        try expectBits(
            case.negative_underflow_bits,
            ieee.scaleBh(case.mode, ieee.fromBitsh(0x8001), -2_147_483_648),
        );
        try expectBits(
            0x0000,
            ieee.scaleBh(case.mode, ieee.fromBitsh(0x0000), -2_147_483_648),
        );
        try expectBits(
            0x8000,
            ieee.scaleBh(case.mode, ieee.fromBitsh(0x8000), 2_147_483_647),
        );
        try expectBits(
            0x7c00,
            ieee.scaleBh(case.mode, ieee.fromBitsh(0x7c00), -2_147_483_648),
        );
        try expectBits(
            0xfc00,
            ieee.scaleBh(case.mode, ieee.fromBitsh(0xfc00), 2_147_483_647),
        );
        try expectBits(
            0x7e55,
            ieee.scaleBh(case.mode, ieee.fromBitsh(0x7e55), 2_147_483_647),
        );
        try expectBits(
            0xff55,
            ieee.scaleBh(case.mode, ieee.fromBitsh(0xfd55), -2_147_483_648),
        );
    }
}

test "roundIntegralh rounds positive and negative ties in every mode" {
    const positive_tie = ieee.fromBitsh(0x4100);
    const negative_tie = ieee.fromBitsh(0xc100);
    const cases = [_]struct {
        mode: ieee.RoundMode,
        positive_bits: u16,
        negative_bits: u16,
    }{
        .{ .mode = .nearest_even, .positive_bits = 0x4000, .negative_bits = 0xc000 },
        .{ .mode = .toward_zero, .positive_bits = 0x4000, .negative_bits = 0xc000 },
        .{ .mode = .toward_positive, .positive_bits = 0x4200, .negative_bits = 0xc000 },
        .{ .mode = .toward_negative, .positive_bits = 0x4000, .negative_bits = 0xc200 },
    };

    for (cases) |case| {
        try expectBits(case.positive_bits, ieee.roundIntegralh(case.mode, positive_tie));
        try expectBits(case.negative_bits, ieee.roundIntegralh(case.mode, negative_tie));
    }
}
