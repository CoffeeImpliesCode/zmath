const std = @import("std");
const common = @import("scalar_c23_common.zig");
const remainder = @import("scalar_remainder.zig");
const types = @import("scalar_c23_types.zig");
// These cases assert exact quotient selection, payload policy, and signed-zero preservation.

fn Operations(comptime T: type) type {
    return struct {
        fmod: *const fn (T, T) T,
        remainder: *const fn (T, T) T,
        remquo: *const fn (T, T) types.RemquoResult(T),
    };
}

fn SpecialBits(comptime Bits: type) type {
    return struct {
        positive_infinity: Bits,
        negative_infinity: Bits,
        quiet_nan_x: Bits,
        signaling_nan_x: Bits,
        quieted_nan_x: Bits,
        quiet_nan_y: Bits,
        signaling_nan_y: Bits,
        quieted_nan_y: Bits,
    };
}

fn expectBits(
    comptime Ieee: type,
    expected: Ieee.Bits,
    actual: Ieee.Float,
) !void {
    try std.testing.expectEqual(expected, Ieee.toBits(actual));
}

fn expectFiniteCases(
    comptime T: type,
    comptime Ieee: type,
    operations: Operations(T),
    huge_power: T,
) !void {
    const negative_zero_bits = Ieee.toBits(@as(T, -0.0));

    try expectBits(Ieee, 0, operations.fmod(6.0, 3.0));
    try expectBits(Ieee, negative_zero_bits, operations.fmod(-6.0, -3.0));
    try expectBits(Ieee, 0, operations.remainder(6.0, -3.0));
    try expectBits(Ieee, negative_zero_bits, operations.remainder(-6.0, 3.0));

    const exact_positive = operations.remquo(6.0, -3.0);
    try expectBits(Ieee, 0, exact_positive.value);
    try std.testing.expectEqual(@as(i32, -2), exact_positive.quotient);

    const exact_negative = operations.remquo(-6.0, -3.0);
    try expectBits(Ieee, negative_zero_bits, exact_negative.value);
    try std.testing.expectEqual(@as(i32, 2), exact_negative.quotient);

    try expectBits(Ieee, negative_zero_bits, operations.fmod(-0.0, 3.0));
    try expectBits(Ieee, negative_zero_bits, operations.remainder(-0.0, -3.0));
    const zero_quotient = operations.remquo(-0.0, -3.0);
    try expectBits(Ieee, negative_zero_bits, zero_quotient.value);
    try std.testing.expectEqual(@as(i32, 0), zero_quotient.quotient);

    try expectBits(Ieee, 0, operations.fmod(6.0, 3.0));
    try expectBits(Ieee, 0, operations.fmod(6.0, -3.0));
    try expectBits(Ieee, negative_zero_bits, operations.fmod(-6.0, 3.0));
    try expectBits(Ieee, negative_zero_bits, operations.fmod(-6.0, -3.0));
    try std.testing.expectEqual(@as(T, 1.0), operations.fmod(7.0, 2.0));
    try std.testing.expectEqual(@as(T, -1.0), operations.fmod(-7.0, 2.0));
    try std.testing.expectEqual(@as(T, 1.0), operations.fmod(7.0, -2.0));
    try std.testing.expectEqual(@as(T, -1.0), operations.fmod(-7.0, -2.0));
    try std.testing.expectEqual(@as(T, 1.0), operations.remainder(5.0, 2.0));
    try std.testing.expectEqual(@as(T, -1.0), operations.remainder(7.0, 2.0));
    try std.testing.expectEqual(@as(T, -1.0), operations.remainder(-5.0, 2.0));
    try std.testing.expectEqual(@as(T, 1.0), operations.remainder(-7.0, 2.0));

    const tie_even_zero = operations.remquo(2.0, 4.0);
    try std.testing.expectEqual(@as(T, 2.0), tie_even_zero.value);
    try std.testing.expectEqual(@as(i32, 0), tie_even_zero.quotient);

    const above_half = operations.remquo(3.0, 4.0);
    try std.testing.expectEqual(@as(T, -1.0), above_half.value);
    try std.testing.expectEqual(@as(i32, 1), above_half.quotient);

    const above_half_negative = operations.remquo(-3.0, 4.0);
    try std.testing.expectEqual(@as(T, 1.0), above_half_negative.value);
    try std.testing.expectEqual(@as(i32, -1), above_half_negative.quotient);

    const signs_pp = operations.remquo(7.0, 2.0);
    try std.testing.expectEqual(@as(T, -1.0), signs_pp.value);
    try std.testing.expectEqual(@as(i32, 4), signs_pp.quotient);

    const signs_np = operations.remquo(-7.0, 2.0);
    try std.testing.expectEqual(@as(T, 1.0), signs_np.value);
    try std.testing.expectEqual(@as(i32, -4), signs_np.quotient);

    const signs_pn = operations.remquo(7.0, -2.0);
    try std.testing.expectEqual(@as(T, -1.0), signs_pn.value);
    try std.testing.expectEqual(@as(i32, -4), signs_pn.quotient);

    const signs_nn = operations.remquo(-7.0, -2.0);
    try std.testing.expectEqual(@as(T, 1.0), signs_nn.value);
    try std.testing.expectEqual(@as(i32, 4), signs_nn.quotient);

    const huge_positive = operations.remquo(huge_power, 3.0);
    try std.testing.expectEqual(@as(T, 1.0), huge_positive.value);
    try std.testing.expect(huge_positive.quotient > 0);
    try std.testing.expectEqual(@as(i32, 5), huge_positive.quotient & 7);

    const huge_negative = operations.remquo(-huge_power, 3.0);
    try std.testing.expectEqual(@as(T, -1.0), huge_negative.value);
    try std.testing.expect(huge_negative.quotient < 0);
    try std.testing.expectEqual(@as(i32, 5), (-huge_negative.quotient) & 7);
}

fn expectSpecialCases(
    comptime T: type,
    comptime Ieee: type,
    operations: Operations(T),
    bits: SpecialBits(Ieee.Bits),
) !void {
    const positive_infinity = Ieee.fromBits(bits.positive_infinity);
    const negative_infinity = Ieee.fromBits(bits.negative_infinity);
    const quiet_nan_x = Ieee.fromBits(bits.quiet_nan_x);
    const signaling_nan_x = Ieee.fromBits(bits.signaling_nan_x);
    const quiet_nan_y = Ieee.fromBits(bits.quiet_nan_y);
    const signaling_nan_y = Ieee.fromBits(bits.signaling_nan_y);
    const positive_zero: T = 0.0;
    const negative_zero: T = -0.0;
    const canonical_nan_bits = Ieee.toBits(Ieee.canonicalNaN());

    try expectBits(Ieee, bits.quiet_nan_x, operations.fmod(quiet_nan_x, signaling_nan_y));
    try expectBits(Ieee, bits.quiet_nan_x, operations.fmod(quiet_nan_x, positive_zero));
    try expectBits(
        Ieee,
        bits.quieted_nan_x,
        operations.remainder(signaling_nan_x, quiet_nan_y),
    );

    const second_signaling = operations.remquo(1.0, signaling_nan_y);
    try expectBits(Ieee, bits.quieted_nan_y, second_signaling.value);
    try std.testing.expectEqual(@as(i32, 0), second_signaling.quotient);

    const nan_precedes_invalid = operations.remquo(positive_infinity, quiet_nan_y);
    try expectBits(Ieee, bits.quiet_nan_y, nan_precedes_invalid.value);
    try std.testing.expectEqual(@as(i32, 0), nan_precedes_invalid.quotient);

    try expectBits(Ieee, canonical_nan_bits, operations.fmod(positive_infinity, 2.0));
    try expectBits(Ieee, canonical_nan_bits, operations.remainder(negative_infinity, 2.0));
    try expectBits(Ieee, canonical_nan_bits, operations.fmod(2.0, positive_zero));
    try expectBits(Ieee, canonical_nan_bits, operations.remainder(2.0, negative_zero));

    const invalid_infinity = operations.remquo(positive_infinity, 2.0);
    try expectBits(Ieee, canonical_nan_bits, invalid_infinity.value);
    try std.testing.expectEqual(@as(i32, 0), invalid_infinity.quotient);

    const invalid_zero = operations.remquo(-2.0, negative_zero);
    try expectBits(Ieee, canonical_nan_bits, invalid_zero.value);
    try std.testing.expectEqual(@as(i32, 0), invalid_zero.quotient);

    const finite: T = -5.0;
    const finite_bits = Ieee.toBits(finite);
    try expectBits(Ieee, finite_bits, operations.fmod(finite, positive_infinity));
    try expectBits(Ieee, finite_bits, operations.remainder(finite, negative_infinity));

    const infinite_divisor = operations.remquo(finite, negative_infinity);
    try expectBits(Ieee, finite_bits, infinite_divisor.value);
    try std.testing.expectEqual(@as(i32, 0), infinite_divisor.quotient);
}

const operations_h: Operations(f16) = .{
    .fmod = &remainder.fmodh,
    .remainder = &remainder.remainderh,
    .remquo = &remainder.remquoh,
};

const operations_f: Operations(f32) = .{
    .fmod = &remainder.fmodf,
    .remainder = &remainder.remainderf,
    .remquo = &remainder.remquof,
};

const operations_d: Operations(f64) = .{
    .fmod = &remainder.fmodd,
    .remainder = &remainder.remainderd,
    .remquo = &remainder.remquod,
};

test "f16 fmod remainder and remquo use exact quotient semantics" {
    try expectFiniteCases(f16, common.H, operations_h, common.H.fromBits(0x7400));
}

test "f32 fmod remainder and remquo use exact quotient semantics" {
    try expectFiniteCases(f32, common.F, operations_f, common.F.fromBits(0x7180_0000));
}

test "f64 fmod remainder and remquo use exact quotient semantics" {
    try expectFiniteCases(f64, common.D, operations_d, common.D.fromBits(0x7e70_0000_0000_0000));
}

test "f16 remainder operations implement C23 special values" {
    try expectSpecialCases(f16, common.H, operations_h, .{
        .positive_infinity = 0x7c00,
        .negative_infinity = 0xfc00,
        .quiet_nan_x = 0xfe55,
        .signaling_nan_x = 0x7d55,
        .quieted_nan_x = 0x7f55,
        .quiet_nan_y = 0x7e33,
        .signaling_nan_y = 0xfd2a,
        .quieted_nan_y = 0xff2a,
    });
}

test "f32 remainder operations implement C23 special values" {
    try expectSpecialCases(f32, common.F, operations_f, .{
        .positive_infinity = 0x7f80_0000,
        .negative_infinity = 0xff80_0000,
        .quiet_nan_x = 0xffc0_0123,
        .signaling_nan_x = 0x7f80_0456,
        .quieted_nan_x = 0x7fc0_0456,
        .quiet_nan_y = 0x7fc0_0abc,
        .signaling_nan_y = 0xff80_0789,
        .quieted_nan_y = 0xffc0_0789,
    });
}

test "f64 remainder operations implement C23 special values" {
    try expectSpecialCases(f64, common.D, operations_d, .{
        .positive_infinity = 0x7ff0_0000_0000_0000,
        .negative_infinity = 0xfff0_0000_0000_0000,
        .quiet_nan_x = 0xfff8_0000_0000_0123,
        .signaling_nan_x = 0x7ff0_0000_0000_0456,
        .quieted_nan_x = 0x7ff8_0000_0000_0456,
        .quiet_nan_y = 0x7ff8_0000_0000_0abc,
        .signaling_nan_y = 0xfff0_0000_0000_0789,
        .quieted_nan_y = 0xfff8_0000_0000_0789,
    });
}

test "remainder public signatures use shared RemquoResult" {
    comptime {
        std.debug.assert(@TypeOf(remainder.fmodh) == fn (f16, f16) f16);
        std.debug.assert(@TypeOf(remainder.fmodf) == fn (f32, f32) f32);
        std.debug.assert(@TypeOf(remainder.fmodd) == fn (f64, f64) f64);
        std.debug.assert(@TypeOf(remainder.remainderh) == fn (f16, f16) f16);
        std.debug.assert(@TypeOf(remainder.remainderf) == fn (f32, f32) f32);
        std.debug.assert(@TypeOf(remainder.remainderd) == fn (f64, f64) f64);
        std.debug.assert(
            @TypeOf(remainder.remquoh) == fn (f16, f16) types.RemquoResult(f16),
        );
        std.debug.assert(
            @TypeOf(remainder.remquof) == fn (f32, f32) types.RemquoResult(f32),
        );
        std.debug.assert(
            @TypeOf(remainder.remquod) == fn (f64, f64) types.RemquoResult(f64),
        );
    }
}
