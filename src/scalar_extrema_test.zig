const std = @import("std");
const common = @import("scalar_c23_common.zig");
const extrema = @import("scalar_extrema.zig");

fn Operations(comptime T: type) type {
    const Binary = *const fn (T, T) T;
    return struct {
        fdim: Binary,
        fmax: Binary,
        fmin: Binary,
        fmaximum: Binary,
        fminimum: Binary,
        fmaximum_mag: Binary,
        fminimum_mag: Binary,
        fmaximum_num: Binary,
        fminimum_num: Binary,
        fmaximum_mag_num: Binary,
        fminimum_mag_num: Binary,
    };
}

fn NaNPatterns(comptime Format: type) type {
    return struct {
        quiet_a: Format.Bits,
        quiet_b: Format.Bits,
        signaling_a: Format.Bits,
        signaling_a_quieted: Format.Bits,
        signaling_b: Format.Bits,
        signaling_b_quieted: Format.Bits,
    };
}

const operations_h: Operations(f16) = .{
    .fdim = &extrema.fdimh,
    .fmax = &extrema.fmaxh,
    .fmin = &extrema.fminh,
    .fmaximum = &extrema.fmaximumh,
    .fminimum = &extrema.fminimumh,
    .fmaximum_mag = &extrema.fmaximumMagh,
    .fminimum_mag = &extrema.fminimumMagh,
    .fmaximum_num = &extrema.fmaximumNumh,
    .fminimum_num = &extrema.fminimumNumh,
    .fmaximum_mag_num = &extrema.fmaximumMagNumh,
    .fminimum_mag_num = &extrema.fminimumMagNumh,
};

const operations_f: Operations(f32) = .{
    .fdim = &extrema.fdimf,
    .fmax = &extrema.fmaxf,
    .fmin = &extrema.fminf,
    .fmaximum = &extrema.fmaximumf,
    .fminimum = &extrema.fminimumf,
    .fmaximum_mag = &extrema.fmaximumMagf,
    .fminimum_mag = &extrema.fminimumMagf,
    .fmaximum_num = &extrema.fmaximumNumf,
    .fminimum_num = &extrema.fminimumNumf,
    .fmaximum_mag_num = &extrema.fmaximumMagNumf,
    .fminimum_mag_num = &extrema.fminimumMagNumf,
};

const operations_d: Operations(f64) = .{
    .fdim = &extrema.fdimd,
    .fmax = &extrema.fmaxd,
    .fmin = &extrema.fmind,
    .fmaximum = &extrema.fmaximumd,
    .fminimum = &extrema.fminimumd,
    .fmaximum_mag = &extrema.fmaximumMagd,
    .fminimum_mag = &extrema.fminimumMagd,
    .fmaximum_num = &extrema.fmaximumNumd,
    .fminimum_num = &extrema.fminimumNumd,
    .fmaximum_mag_num = &extrema.fmaximumMagNumd,
    .fminimum_mag_num = &extrema.fminimumMagNumd,
};

fn expectBits(
    comptime Format: type,
    expected: Format.Bits,
    actual: Format.Float,
) !void {
    try std.testing.expectEqual(expected, Format.toBits(actual));
}

fn expectOrderedInfinityAndZeroBehavior(
    comptime Format: type,
    operations: Operations(Format.Float),
) !void {
    const T = Format.Float;
    const Bits = Format.Bits;
    const Binary = *const fn (T, T) T;

    const positive_zero: T = 0.0;
    const negative_zero: T = -0.0;
    const three: T = 3.0;
    const five: T = 5.0;
    const positive_infinity = std.math.inf(T);
    const negative_infinity = -std.math.inf(T);

    try expectBits(Format, Format.toBits(@as(T, 2.0)), operations.fdim(five, three));
    try expectBits(Format, Format.toBits(positive_zero), operations.fdim(three, five));
    try expectBits(
        Format,
        Format.toBits(positive_zero),
        operations.fdim(positive_infinity, positive_infinity),
    );
    try expectBits(
        Format,
        Format.toBits(positive_zero),
        operations.fdim(negative_infinity, negative_infinity),
    );
    try expectBits(
        Format,
        Format.toBits(positive_infinity),
        operations.fdim(positive_infinity, three),
    );
    try expectBits(
        Format,
        Format.toBits(positive_infinity),
        operations.fdim(three, negative_infinity),
    );
    try expectBits(
        Format,
        Format.toBits(positive_zero),
        operations.fdim(negative_infinity, three),
    );

    const maximum_operations = [_]Binary{
        operations.fmax,
        operations.fmaximum,
        operations.fmaximum_num,
    };
    const minimum_operations = [_]Binary{
        operations.fmin,
        operations.fminimum,
        operations.fminimum_num,
    };
    for (maximum_operations) |operation| {
        try expectBits(Format, Format.toBits(five), operation(three, five));
        try expectBits(Format, Format.toBits(five), operation(five, three));
        try expectBits(
            Format,
            Format.toBits(positive_infinity),
            operation(negative_infinity, positive_infinity),
        );
        try expectBits(
            Format,
            Format.toBits(positive_infinity),
            operation(positive_infinity, negative_infinity),
        );
    }
    for (minimum_operations) |operation| {
        try expectBits(Format, Format.toBits(three), operation(three, five));
        try expectBits(Format, Format.toBits(three), operation(five, three));
        try expectBits(
            Format,
            Format.toBits(negative_infinity),
            operation(negative_infinity, positive_infinity),
        );
        try expectBits(
            Format,
            Format.toBits(negative_infinity),
            operation(positive_infinity, negative_infinity),
        );
    }

    const magnitude_maximum_operations = [_]Binary{
        operations.fmaximum_mag,
        operations.fmaximum_mag_num,
    };
    const magnitude_minimum_operations = [_]Binary{
        operations.fminimum_mag,
        operations.fminimum_mag_num,
    };
    for (magnitude_maximum_operations) |operation| {
        try expectBits(Format, Format.toBits(-five), operation(-five, three));
        try expectBits(Format, Format.toBits(-five), operation(three, -five));
        try expectBits(
            Format,
            Format.toBits(negative_infinity),
            operation(negative_infinity, five),
        );
        try expectBits(
            Format,
            Format.toBits(negative_infinity),
            operation(five, negative_infinity),
        );
    }
    for (magnitude_minimum_operations) |operation| {
        try expectBits(Format, Format.toBits(three), operation(-five, three));
        try expectBits(Format, Format.toBits(three), operation(three, -five));
        try expectBits(Format, Format.toBits(five), operation(negative_infinity, five));
        try expectBits(Format, Format.toBits(five), operation(five, negative_infinity));
    }

    const ZeroCase = struct {
        first: T,
        second: T,
        maximum_bits: Bits,
        minimum_bits: Bits,
    };
    const zero_cases = [_]ZeroCase{
        .{
            .first = positive_zero,
            .second = positive_zero,
            .maximum_bits = Format.toBits(positive_zero),
            .minimum_bits = Format.toBits(positive_zero),
        },
        .{
            .first = positive_zero,
            .second = negative_zero,
            .maximum_bits = Format.toBits(positive_zero),
            .minimum_bits = Format.toBits(negative_zero),
        },
        .{
            .first = negative_zero,
            .second = positive_zero,
            .maximum_bits = Format.toBits(positive_zero),
            .minimum_bits = Format.toBits(negative_zero),
        },
        .{
            .first = negative_zero,
            .second = negative_zero,
            .maximum_bits = Format.toBits(negative_zero),
            .minimum_bits = Format.toBits(negative_zero),
        },
    };
    const all_maximum_operations = [_]Binary{
        operations.fmax,
        operations.fmaximum,
        operations.fmaximum_mag,
        operations.fmaximum_num,
        operations.fmaximum_mag_num,
    };
    const all_minimum_operations = [_]Binary{
        operations.fmin,
        operations.fminimum,
        operations.fminimum_mag,
        operations.fminimum_num,
        operations.fminimum_mag_num,
    };
    for (zero_cases) |zero_case| {
        try expectBits(
            Format,
            Format.toBits(positive_zero),
            operations.fdim(zero_case.first, zero_case.second),
        );
        for (all_maximum_operations) |operation| {
            try expectBits(
                Format,
                zero_case.maximum_bits,
                operation(zero_case.first, zero_case.second),
            );
        }
        for (all_minimum_operations) |operation| {
            try expectBits(
                Format,
                zero_case.minimum_bits,
                operation(zero_case.first, zero_case.second),
            );
        }
    }
}

fn expectNaNBehavior(
    comptime Format: type,
    operations: Operations(Format.Float),
    patterns: NaNPatterns(Format),
) !void {
    const T = Format.Float;
    const Bits = Format.Bits;
    const Binary = *const fn (T, T) T;
    const NanValue = struct {
        value: T,
        quiet_bits: Bits,
    };

    const nan_values = [_]NanValue{
        .{
            .value = Format.fromBits(patterns.quiet_a),
            .quiet_bits = patterns.quiet_a,
        },
        .{
            .value = Format.fromBits(patterns.quiet_b),
            .quiet_bits = patterns.quiet_b,
        },
        .{
            .value = Format.fromBits(patterns.signaling_a),
            .quiet_bits = patterns.signaling_a_quieted,
        },
        .{
            .value = Format.fromBits(patterns.signaling_b),
            .quiet_bits = patterns.signaling_b_quieted,
        },
    };
    const numeric: T = -3.25;
    const numeric_bits = Format.toBits(numeric);

    const propagating_operations = [_]Binary{
        operations.fdim,
        operations.fmaximum,
        operations.fminimum,
        operations.fmaximum_mag,
        operations.fminimum_mag,
    };
    const number_operations = [_]Binary{
        operations.fmaximum_num,
        operations.fminimum_num,
        operations.fmaximum_mag_num,
        operations.fminimum_mag_num,
    };
    const legacy_operations = [_]Binary{
        operations.fmax,
        operations.fmin,
    };

    for (nan_values) |nan_value| {
        for (propagating_operations) |operation| {
            try expectBits(Format, nan_value.quiet_bits, operation(nan_value.value, numeric));
            try expectBits(Format, nan_value.quiet_bits, operation(numeric, nan_value.value));
        }
        for (number_operations) |operation| {
            try expectBits(Format, numeric_bits, operation(nan_value.value, numeric));
            try expectBits(Format, numeric_bits, operation(numeric, nan_value.value));
        }
        for (legacy_operations) |operation| {
            try expectBits(Format, numeric_bits, operation(nan_value.value, numeric));
            try expectBits(Format, numeric_bits, operation(numeric, nan_value.value));
        }
    }

    for (nan_values) |first| {
        for (nan_values) |second| {
            for (propagating_operations) |operation| {
                try expectBits(Format, first.quiet_bits, operation(first.value, second.value));
            }
            for (number_operations) |operation| {
                try expectBits(Format, first.quiet_bits, operation(first.value, second.value));
            }
            for (legacy_operations) |operation| {
                try expectBits(Format, first.quiet_bits, operation(first.value, second.value));
            }
        }
    }
}

fn expectEqualMagnitudeBehavior(
    comptime Format: type,
    operations: Operations(Format.Float),
) !void {
    const T = Format.Float;
    const Binary = *const fn (T, T) T;
    const three: T = 3.0;
    const positive_infinity = std.math.inf(T);
    const negative_infinity = -std.math.inf(T);
    const maximum_operations = [_]Binary{
        operations.fmaximum_mag,
        operations.fmaximum_mag_num,
    };
    const minimum_operations = [_]Binary{
        operations.fminimum_mag,
        operations.fminimum_mag_num,
    };

    for (maximum_operations) |operation| {
        try expectBits(Format, Format.toBits(three), operation(three, -three));
        try expectBits(Format, Format.toBits(three), operation(-three, three));
        try expectBits(
            Format,
            Format.toBits(positive_infinity),
            operation(positive_infinity, negative_infinity),
        );
        try expectBits(
            Format,
            Format.toBits(positive_infinity),
            operation(negative_infinity, positive_infinity),
        );
    }
    for (minimum_operations) |operation| {
        try expectBits(Format, Format.toBits(-three), operation(three, -three));
        try expectBits(Format, Format.toBits(-three), operation(-three, three));
        try expectBits(
            Format,
            Format.toBits(negative_infinity),
            operation(positive_infinity, negative_infinity),
        );
        try expectBits(
            Format,
            Format.toBits(negative_infinity),
            operation(negative_infinity, positive_infinity),
        );
    }
}

test "extrema expose every H F and D signature" {
    comptime {
        std.debug.assert(@TypeOf(extrema.fdimh) == fn (f16, f16) f16);
        std.debug.assert(@TypeOf(extrema.fdimf) == fn (f32, f32) f32);
        std.debug.assert(@TypeOf(extrema.fdimd) == fn (f64, f64) f64);
        std.debug.assert(@TypeOf(extrema.fmaxh) == fn (f16, f16) f16);
        std.debug.assert(@TypeOf(extrema.fmaxf) == fn (f32, f32) f32);
        std.debug.assert(@TypeOf(extrema.fmaxd) == fn (f64, f64) f64);
        std.debug.assert(@TypeOf(extrema.fminh) == fn (f16, f16) f16);
        std.debug.assert(@TypeOf(extrema.fminf) == fn (f32, f32) f32);
        std.debug.assert(@TypeOf(extrema.fmind) == fn (f64, f64) f64);
        std.debug.assert(@TypeOf(extrema.fmaximumh) == fn (f16, f16) f16);
        std.debug.assert(@TypeOf(extrema.fmaximumf) == fn (f32, f32) f32);
        std.debug.assert(@TypeOf(extrema.fmaximumd) == fn (f64, f64) f64);
        std.debug.assert(@TypeOf(extrema.fminimumh) == fn (f16, f16) f16);
        std.debug.assert(@TypeOf(extrema.fminimumf) == fn (f32, f32) f32);
        std.debug.assert(@TypeOf(extrema.fminimumd) == fn (f64, f64) f64);
        std.debug.assert(@TypeOf(extrema.fmaximumMagh) == fn (f16, f16) f16);
        std.debug.assert(@TypeOf(extrema.fmaximumMagf) == fn (f32, f32) f32);
        std.debug.assert(@TypeOf(extrema.fmaximumMagd) == fn (f64, f64) f64);
        std.debug.assert(@TypeOf(extrema.fminimumMagh) == fn (f16, f16) f16);
        std.debug.assert(@TypeOf(extrema.fminimumMagf) == fn (f32, f32) f32);
        std.debug.assert(@TypeOf(extrema.fminimumMagd) == fn (f64, f64) f64);
        std.debug.assert(@TypeOf(extrema.fmaximumNumh) == fn (f16, f16) f16);
        std.debug.assert(@TypeOf(extrema.fmaximumNumf) == fn (f32, f32) f32);
        std.debug.assert(@TypeOf(extrema.fmaximumNumd) == fn (f64, f64) f64);
        std.debug.assert(@TypeOf(extrema.fminimumNumh) == fn (f16, f16) f16);
        std.debug.assert(@TypeOf(extrema.fminimumNumf) == fn (f32, f32) f32);
        std.debug.assert(@TypeOf(extrema.fminimumNumd) == fn (f64, f64) f64);
        std.debug.assert(@TypeOf(extrema.fmaximumMagNumh) == fn (f16, f16) f16);
        std.debug.assert(@TypeOf(extrema.fmaximumMagNumf) == fn (f32, f32) f32);
        std.debug.assert(@TypeOf(extrema.fmaximumMagNumd) == fn (f64, f64) f64);
        std.debug.assert(@TypeOf(extrema.fminimumMagNumh) == fn (f16, f16) f16);
        std.debug.assert(@TypeOf(extrema.fminimumMagNumf) == fn (f32, f32) f32);
        std.debug.assert(@TypeOf(extrema.fminimumMagNumd) == fn (f64, f64) f64);
    }
}

test "f16 extrema cover ordered values infinities and every signed-zero tie" {
    try expectOrderedInfinityAndZeroBehavior(common.H, operations_h);
}

test "f32 extrema cover ordered values infinities and every signed-zero tie" {
    try expectOrderedInfinityAndZeroBehavior(common.F, operations_f);
}

test "f64 extrema cover ordered values infinities and every signed-zero tie" {
    try expectOrderedInfinityAndZeroBehavior(common.D, operations_d);
}

test "f16 extrema preserve first NaN bits and prefer numeric operands" {
    try expectNaNBehavior(common.H, operations_h, .{
        .quiet_a = 0xfe55,
        .quiet_b = 0x7e2a,
        .signaling_a = 0xfd33,
        .signaling_a_quieted = 0xff33,
        .signaling_b = 0x7c44,
        .signaling_b_quieted = 0x7e44,
    });
}

test "f32 extrema preserve first NaN bits and prefer numeric operands" {
    try expectNaNBehavior(common.F, operations_f, .{
        .quiet_a = 0xffc1_2345,
        .quiet_b = 0x7fc5_4321,
        .signaling_a = 0xff81_1357,
        .signaling_a_quieted = 0xffc1_1357,
        .signaling_b = 0x7f82_2468,
        .signaling_b_quieted = 0x7fc2_2468,
    });
}

test "f64 extrema preserve first NaN bits and prefer numeric operands" {
    try expectNaNBehavior(common.D, operations_d, .{
        .quiet_a = 0xfff8_1234_5678_9abc,
        .quiet_b = 0x7ffc_0fed_cba9_8765,
        .signaling_a = 0xfff0_1357_9bdf_2468,
        .signaling_a_quieted = 0xfff8_1357_9bdf_2468,
        .signaling_b = 0x7ff0_2468_ace0_1357,
        .signaling_b_quieted = 0x7ff8_2468_ace0_1357,
    });
}

test "f16 magnitude extrema resolve equal magnitudes through maximum and minimum" {
    try expectEqualMagnitudeBehavior(common.H, operations_h);
}

test "f32 magnitude extrema resolve equal magnitudes through maximum and minimum" {
    try expectEqualMagnitudeBehavior(common.F, operations_f);
}

test "f64 magnitude extrema resolve equal magnitudes through maximum and minimum" {
    try expectEqualMagnitudeBehavior(common.D, operations_d);
}
