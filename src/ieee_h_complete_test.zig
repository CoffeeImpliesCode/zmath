const ieee = @import("ieee.zig");

fn expectBits(expected: u16, actual: f16) !void {
    if (ieee.toBitsH(actual) != expected) return error.TestExpectedEqual;
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
        if (ieee.toBitsH(ieee.fromBitsH(bits)) != bits) return error.TestExpectedEqual;
    }
}

test "nextUpH and nextDownH step across both signed zeros" {
    const positive_zero = ieee.fromBitsH(0x0000);
    const negative_zero = ieee.fromBitsH(0x8000);

    try expectBits(0x0001, ieee.nextUpH(positive_zero));
    try expectBits(0x0001, ieee.nextUpH(negative_zero));
    try expectBits(0x8001, ieee.nextDownH(positive_zero));
    try expectBits(0x8001, ieee.nextDownH(negative_zero));
}

test "nextUpH and nextDownH handle infinity boundaries" {
    const positive_infinity = ieee.fromBitsH(0x7c00);
    const negative_infinity = ieee.fromBitsH(0xfc00);

    try expectBits(0x7c00, ieee.nextUpH(positive_infinity));
    try expectBits(0x7bff, ieee.nextDownH(positive_infinity));
    try expectBits(0xfbff, ieee.nextUpH(negative_infinity));
    try expectBits(0xfc00, ieee.nextDownH(negative_infinity));
}

test "nextUpH and nextDownH quiet signaling NaNs" {
    const positive_signaling_nan = ieee.fromBitsH(0x7d55);
    const negative_signaling_nan = ieee.fromBitsH(0xfd55);

    try expectBits(0x7f55, ieee.nextUpH(positive_signaling_nan));
    try expectBits(0x7f55, ieee.nextDownH(positive_signaling_nan));
    try expectBits(0xff55, ieee.nextUpH(negative_signaling_nan));
    try expectBits(0xff55, ieee.nextDownH(negative_signaling_nan));
}

test "scaleBH rounds positive and negative overflow in every mode" {
    const positive_maximum = ieee.fromBitsH(0x7bff);
    const negative_maximum = ieee.fromBitsH(0xfbff);
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
        try expectBits(case.positive_bits, ieee.scaleBH(case.mode, positive_maximum, 1));
        try expectBits(case.negative_bits, ieee.scaleBH(case.mode, negative_maximum, 1));
    }
}

test "scaleBH rounds positive and negative underflow in every mode" {
    const positive_minimum = ieee.fromBitsH(0x0001);
    const negative_minimum = ieee.fromBitsH(0x8001);
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
        try expectBits(case.positive_bits, ieee.scaleBH(case.mode, positive_minimum, -1));
        try expectBits(case.negative_bits, ieee.scaleBH(case.mode, negative_minimum, -1));
    }
}

test "roundIntegralH rounds positive and negative ties in every mode" {
    const positive_tie = ieee.fromBitsH(0x4100);
    const negative_tie = ieee.fromBitsH(0xc100);
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
        try expectBits(case.positive_bits, ieee.roundIntegralH(case.mode, positive_tie));
        try expectBits(case.negative_bits, ieee.roundIntegralH(case.mode, negative_tie));
    }
}
