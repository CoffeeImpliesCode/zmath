const std = @import("std");
const common = @import("scalar_c23_common.zig");

fn expectHBits(expected: u16, value: f16) !void {
    try std.testing.expectEqual(expected, common.H.toBits(value));
}

fn expectFBits(expected: u32, value: f32) !void {
    try std.testing.expectEqual(expected, common.F.toBits(value));
}

fn expectDBits(expected: u64, value: f64) !void {
    try std.testing.expectEqual(expected, common.D.toBits(value));
}

test "f16 preserves NaN payloads and quiets signaling NaNs" {
    try std.testing.expectEqual(common.Class.quiet_nan, common.H.classifyBits(0x7e55));
    try std.testing.expectEqual(common.Class.quiet_nan, common.H.classifyBits(0xfe55));
    try std.testing.expectEqual(common.Class.signaling_nan, common.H.classifyBits(0x7d55));
    try std.testing.expectEqual(common.Class.signaling_nan, common.H.classifyBits(0xfd55));

    try std.testing.expectEqual(@as(u16, 0x7e55), common.H.quietNaNBits(0x7e55));
    try std.testing.expectEqual(@as(u16, 0xfe55), common.H.quietNaNBits(0xfe55));
    try std.testing.expectEqual(@as(u16, 0x7f55), common.H.quietNaNBits(0x7d55));
    try std.testing.expectEqual(@as(u16, 0xff55), common.H.quietNaNBits(0xfd55));
    try expectHBits(0x7e00, common.H.canonicalNaN());

    const first = common.H.firstNaN3(
        common.H.fromBits(0x3c00),
        common.H.fromBits(0xfd55),
        common.H.fromBits(0x7e55),
    );
    try std.testing.expect(first != null);
    try expectHBits(0xff55, first.?);

    const second = common.H.firstNaN2(
        common.H.fromBits(0x3c00),
        common.H.fromBits(0x7d55),
    );
    try std.testing.expect(second != null);
    try expectHBits(0x7f55, second.?);
    try std.testing.expect(common.H.firstNaN2(
        common.H.fromBits(0x3c00),
        common.H.fromBits(0x4000),
    ) == null);
}

test "f16 recognizes signed zero and infinity exactly" {
    try std.testing.expect(common.H.isZeroBits(0x0000));
    try std.testing.expect(common.H.isZeroBits(0x8000));
    try std.testing.expect(common.H.isPositiveZeroBits(0x0000));
    try std.testing.expect(common.H.isNegativeZeroBits(0x8000));
    try std.testing.expect(!common.H.isPositiveZeroBits(0x8000));
    try std.testing.expect(!common.H.isNegativeZeroBits(0x0000));

    try std.testing.expect(common.H.isInfinityBits(0x7c00));
    try std.testing.expect(common.H.isInfinityBits(0xfc00));
    try std.testing.expect(common.H.isPositiveInfinityBits(0x7c00));
    try std.testing.expect(common.H.isNegativeInfinityBits(0xfc00));
    try std.testing.expect(!common.H.isPositiveInfinityBits(0xfc00));
    try std.testing.expect(!common.H.isNegativeInfinityBits(0x7c00));
}

test "f16 finite ordering keys span every encoding boundary" {
    const bits = [_]u16{
        0xfbff,
        0x8400,
        0x83ff,
        0x8000,
        0x0000,
        0x0001,
        0x03ff,
        0x0400,
        0x7bff,
    };
    const expected_keys = [_]u16{
        0x0400,
        0x7bff,
        0x7c00,
        0x7fff,
        0x8000,
        0x8001,
        0x83ff,
        0x8400,
        0xfbff,
    };
    for (bits, expected_keys) |value_bits, expected| {
        try std.testing.expectEqual(expected, common.H.finiteOrderKey(value_bits));
    }
    for (expected_keys[0 .. expected_keys.len - 1], expected_keys[1..]) |left, right| {
        try std.testing.expect(left < right);
    }
}

test "f32 preserves NaN payloads and quiets signaling NaNs" {
    try std.testing.expectEqual(common.Class.quiet_nan, common.F.classifyBits(0x7fc0_0123));
    try std.testing.expectEqual(common.Class.quiet_nan, common.F.classifyBits(0xffc0_0123));
    try std.testing.expectEqual(common.Class.signaling_nan, common.F.classifyBits(0x7f80_0123));
    try std.testing.expectEqual(common.Class.signaling_nan, common.F.classifyBits(0xff80_0123));

    try std.testing.expectEqual(@as(u32, 0x7fc0_0123), common.F.quietNaNBits(0x7fc0_0123));
    try std.testing.expectEqual(@as(u32, 0xffc0_0123), common.F.quietNaNBits(0xffc0_0123));
    try std.testing.expectEqual(@as(u32, 0x7fc0_0123), common.F.quietNaNBits(0x7f80_0123));
    try std.testing.expectEqual(@as(u32, 0xffc0_0123), common.F.quietNaNBits(0xff80_0123));
    try expectFBits(0x7fc0_0000, common.F.canonicalNaN());

    const first = common.F.firstNaN3(
        common.F.fromBits(0x3f80_0000),
        common.F.fromBits(0xff80_0123),
        common.F.fromBits(0x7fc0_0456),
    );
    try std.testing.expect(first != null);
    try expectFBits(0xffc0_0123, first.?);

    const second = common.F.firstNaN2(
        common.F.fromBits(0x3f80_0000),
        common.F.fromBits(0x7f80_0456),
    );
    try std.testing.expect(second != null);
    try expectFBits(0x7fc0_0456, second.?);
    try std.testing.expect(common.F.firstNaN2(
        common.F.fromBits(0x3f80_0000),
        common.F.fromBits(0x4000_0000),
    ) == null);
}

test "f32 recognizes signed zero and infinity exactly" {
    try std.testing.expect(common.F.isZeroBits(0x0000_0000));
    try std.testing.expect(common.F.isZeroBits(0x8000_0000));
    try std.testing.expect(common.F.isPositiveZeroBits(0x0000_0000));
    try std.testing.expect(common.F.isNegativeZeroBits(0x8000_0000));
    try std.testing.expect(!common.F.isPositiveZeroBits(0x8000_0000));
    try std.testing.expect(!common.F.isNegativeZeroBits(0x0000_0000));

    try std.testing.expect(common.F.isInfinityBits(0x7f80_0000));
    try std.testing.expect(common.F.isInfinityBits(0xff80_0000));
    try std.testing.expect(common.F.isPositiveInfinityBits(0x7f80_0000));
    try std.testing.expect(common.F.isNegativeInfinityBits(0xff80_0000));
    try std.testing.expect(!common.F.isPositiveInfinityBits(0xff80_0000));
    try std.testing.expect(!common.F.isNegativeInfinityBits(0x7f80_0000));
}

test "f32 finite ordering keys span every encoding boundary" {
    const bits = [_]u32{
        0xff7f_ffff,
        0x8080_0000,
        0x807f_ffff,
        0x8000_0000,
        0x0000_0000,
        0x0000_0001,
        0x007f_ffff,
        0x0080_0000,
        0x7f7f_ffff,
    };
    const expected_keys = [_]u32{
        0x0080_0000,
        0x7f7f_ffff,
        0x7f80_0000,
        0x7fff_ffff,
        0x8000_0000,
        0x8000_0001,
        0x807f_ffff,
        0x8080_0000,
        0xff7f_ffff,
    };
    for (bits, expected_keys) |value_bits, expected| {
        try std.testing.expectEqual(expected, common.F.finiteOrderKey(value_bits));
    }
    for (expected_keys[0 .. expected_keys.len - 1], expected_keys[1..]) |left, right| {
        try std.testing.expect(left < right);
    }
}

test "f64 preserves NaN payloads and quiets signaling NaNs" {
    try std.testing.expectEqual(
        common.Class.quiet_nan,
        common.D.classifyBits(0x7ff8_0000_0000_0123),
    );
    try std.testing.expectEqual(
        common.Class.quiet_nan,
        common.D.classifyBits(0xfff8_0000_0000_0123),
    );
    try std.testing.expectEqual(
        common.Class.signaling_nan,
        common.D.classifyBits(0x7ff0_0000_0000_0123),
    );
    try std.testing.expectEqual(
        common.Class.signaling_nan,
        common.D.classifyBits(0xfff0_0000_0000_0123),
    );

    try std.testing.expectEqual(
        @as(u64, 0x7ff8_0000_0000_0123),
        common.D.quietNaNBits(0x7ff8_0000_0000_0123),
    );
    try std.testing.expectEqual(
        @as(u64, 0xfff8_0000_0000_0123),
        common.D.quietNaNBits(0xfff8_0000_0000_0123),
    );
    try std.testing.expectEqual(
        @as(u64, 0x7ff8_0000_0000_0123),
        common.D.quietNaNBits(0x7ff0_0000_0000_0123),
    );
    try std.testing.expectEqual(
        @as(u64, 0xfff8_0000_0000_0123),
        common.D.quietNaNBits(0xfff0_0000_0000_0123),
    );
    try expectDBits(0x7ff8_0000_0000_0000, common.D.canonicalNaN());

    const first = common.D.firstNaN3(
        common.D.fromBits(0x3ff0_0000_0000_0000),
        common.D.fromBits(0xfff0_0000_0000_0123),
        common.D.fromBits(0x7ff8_0000_0000_0456),
    );
    try std.testing.expect(first != null);
    try expectDBits(0xfff8_0000_0000_0123, first.?);

    const second = common.D.firstNaN2(
        common.D.fromBits(0x3ff0_0000_0000_0000),
        common.D.fromBits(0x7ff0_0000_0000_0456),
    );
    try std.testing.expect(second != null);
    try expectDBits(0x7ff8_0000_0000_0456, second.?);
    try std.testing.expect(common.D.firstNaN2(
        common.D.fromBits(0x3ff0_0000_0000_0000),
        common.D.fromBits(0x4000_0000_0000_0000),
    ) == null);
}

test "f64 recognizes signed zero and infinity exactly" {
    try std.testing.expect(common.D.isZeroBits(0x0000_0000_0000_0000));
    try std.testing.expect(common.D.isZeroBits(0x8000_0000_0000_0000));
    try std.testing.expect(common.D.isPositiveZeroBits(0x0000_0000_0000_0000));
    try std.testing.expect(common.D.isNegativeZeroBits(0x8000_0000_0000_0000));
    try std.testing.expect(!common.D.isPositiveZeroBits(0x8000_0000_0000_0000));
    try std.testing.expect(!common.D.isNegativeZeroBits(0x0000_0000_0000_0000));

    try std.testing.expect(common.D.isInfinityBits(0x7ff0_0000_0000_0000));
    try std.testing.expect(common.D.isInfinityBits(0xfff0_0000_0000_0000));
    try std.testing.expect(common.D.isPositiveInfinityBits(0x7ff0_0000_0000_0000));
    try std.testing.expect(common.D.isNegativeInfinityBits(0xfff0_0000_0000_0000));
    try std.testing.expect(!common.D.isPositiveInfinityBits(0xfff0_0000_0000_0000));
    try std.testing.expect(!common.D.isNegativeInfinityBits(0x7ff0_0000_0000_0000));
}

test "f64 finite ordering keys span every encoding boundary" {
    const bits = [_]u64{
        0xffef_ffff_ffff_ffff,
        0x8010_0000_0000_0000,
        0x800f_ffff_ffff_ffff,
        0x8000_0000_0000_0000,
        0x0000_0000_0000_0000,
        0x0000_0000_0000_0001,
        0x000f_ffff_ffff_ffff,
        0x0010_0000_0000_0000,
        0x7fef_ffff_ffff_ffff,
    };
    const expected_keys = [_]u64{
        0x0010_0000_0000_0000,
        0x7fef_ffff_ffff_ffff,
        0x7ff0_0000_0000_0000,
        0x7fff_ffff_ffff_ffff,
        0x8000_0000_0000_0000,
        0x8000_0000_0000_0001,
        0x800f_ffff_ffff_ffff,
        0x8010_0000_0000_0000,
        0xffef_ffff_ffff_ffff,
    };
    for (bits, expected_keys) |value_bits, expected| {
        try std.testing.expectEqual(expected, common.D.finiteOrderKey(value_bits));
    }
    for (expected_keys[0 .. expected_keys.len - 1], expected_keys[1..]) |left, right| {
        try std.testing.expect(left < right);
    }
}
test "f16 exact boundary classes and NaN payload edges" {
    const cases = [_]struct {
        bits: u16,
        category: common.Class,
        is_nan: bool,
    }{
        .{ .bits = 0x0000, .category = .zero, .is_nan = false },
        .{ .bits = 0x8000, .category = .zero, .is_nan = false },
        .{ .bits = 0x0001, .category = .subnormal, .is_nan = false },
        .{ .bits = 0x03ff, .category = .subnormal, .is_nan = false },
        .{ .bits = 0x0400, .category = .normal, .is_nan = false },
        .{ .bits = 0x7bff, .category = .normal, .is_nan = false },
        .{ .bits = 0x7c00, .category = .infinity, .is_nan = false },
        .{ .bits = 0xfc00, .category = .infinity, .is_nan = false },
        .{ .bits = 0x7e00, .category = .quiet_nan, .is_nan = true },
        .{ .bits = 0xfe00, .category = .quiet_nan, .is_nan = true },
        .{ .bits = 0x7c01, .category = .signaling_nan, .is_nan = true },
        .{ .bits = 0xfc01, .category = .signaling_nan, .is_nan = true },
    };

    for (cases) |case| {
        const value = common.H.fromBits(case.bits);
        try std.testing.expectEqual(case.bits, common.H.toBits(value));
        try std.testing.expectEqual(case.category, common.H.classifyBits(case.bits));
        try std.testing.expectEqual(case.is_nan, common.H.isNaNBits(case.bits));
    }

    const payloads = [_]struct {
        signaling: u16,
        quiet: u16,
    }{
        .{ .signaling = 0x7c01, .quiet = 0x7e01 },
        .{ .signaling = 0x7d55, .quiet = 0x7f55 },
        .{ .signaling = 0x7dff, .quiet = 0x7fff },
        .{ .signaling = 0xfd55, .quiet = 0xff55 },
        .{ .signaling = 0xfdff, .quiet = 0xffff },
    };
    for (payloads) |payload| {
        try std.testing.expectEqual(
            common.Class.signaling_nan,
            common.H.classifyBits(payload.signaling),
        );
        try std.testing.expectEqual(
            common.Class.quiet_nan,
            common.H.classifyBits(payload.quiet),
        );
        try std.testing.expectEqual(payload.quiet, common.H.quietNaNBits(payload.signaling));
        try std.testing.expectEqual(payload.quiet, common.H.quietNaNBits(payload.quiet));
        try std.testing.expectEqual(payload.quiet, common.H.toBits(common.H.fromBits(payload.quiet)));
    }
    try expectHBits(0x7e00, common.H.canonicalNaN());

    const normal = common.H.fromBits(0x3c00);
    const signaling = common.H.fromBits(0x7d55);
    const quiet = common.H.fromBits(0x7e55);
    const first = common.H.firstNaN2(signaling, quiet);
    try std.testing.expect(first != null);
    try expectHBits(0x7f55, first.?);
    const second = common.H.firstNaN2(normal, signaling);
    try std.testing.expect(second != null);
    try expectHBits(0x7f55, second.?);
    try std.testing.expect(common.H.firstNaN2(normal, normal) == null);

    const first_three = common.H.firstNaN3(signaling, quiet, normal);
    try std.testing.expect(first_three != null);
    try expectHBits(0x7f55, first_three.?);
    const second_three = common.H.firstNaN3(normal, signaling, quiet);
    try std.testing.expect(second_three != null);
    try expectHBits(0x7f55, second_three.?);
    const third_three = common.H.firstNaN3(normal, normal, signaling);
    try std.testing.expect(third_three != null);
    try expectHBits(0x7f55, third_three.?);
    try std.testing.expect(common.H.firstNaN3(normal, normal, normal) == null);
}

test "f16 exact finite order key endpoints" {
    const cases = [_]struct {
        bits: u16,
        key: u16,
    }{
        .{ .bits = 0xfc00, .key = 0x03ff },
        .{ .bits = 0xfbff, .key = 0x0400 },
        .{ .bits = 0x8400, .key = 0x7bff },
        .{ .bits = 0x83ff, .key = 0x7c00 },
        .{ .bits = 0x8000, .key = 0x7fff },
        .{ .bits = 0x0000, .key = 0x8000 },
        .{ .bits = 0x0001, .key = 0x8001 },
        .{ .bits = 0x03ff, .key = 0x83ff },
        .{ .bits = 0x0400, .key = 0x8400 },
        .{ .bits = 0x7bff, .key = 0xfbff },
        .{ .bits = 0x7c00, .key = 0xfc00 },
    };
    for (cases) |case| {
        try std.testing.expectEqual(case.key, common.H.finiteOrderKey(case.bits));
    }
}

test "f32 exact boundary classes and NaN payload edges" {
    const cases = [_]struct {
        bits: u32,
        category: common.Class,
        is_nan: bool,
    }{
        .{ .bits = 0x0000_0000, .category = .zero, .is_nan = false },
        .{ .bits = 0x8000_0000, .category = .zero, .is_nan = false },
        .{ .bits = 0x0000_0001, .category = .subnormal, .is_nan = false },
        .{ .bits = 0x007f_ffff, .category = .subnormal, .is_nan = false },
        .{ .bits = 0x0080_0000, .category = .normal, .is_nan = false },
        .{ .bits = 0x7f7f_ffff, .category = .normal, .is_nan = false },
        .{ .bits = 0x7f80_0000, .category = .infinity, .is_nan = false },
        .{ .bits = 0xff80_0000, .category = .infinity, .is_nan = false },
        .{ .bits = 0x7fc0_0000, .category = .quiet_nan, .is_nan = true },
        .{ .bits = 0xffc0_0000, .category = .quiet_nan, .is_nan = true },
        .{ .bits = 0x7f80_0001, .category = .signaling_nan, .is_nan = true },
        .{ .bits = 0xff80_0001, .category = .signaling_nan, .is_nan = true },
    };

    for (cases) |case| {
        const value = common.F.fromBits(case.bits);
        try std.testing.expectEqual(case.bits, common.F.toBits(value));
        try std.testing.expectEqual(case.category, common.F.classifyBits(case.bits));
        try std.testing.expectEqual(case.is_nan, common.F.isNaNBits(case.bits));
    }

    const payloads = [_]struct {
        signaling: u32,
        quiet: u32,
    }{
        .{ .signaling = 0x7f80_0001, .quiet = 0x7fc0_0001 },
        .{ .signaling = 0x7f80_0123, .quiet = 0x7fc0_0123 },
        .{ .signaling = 0x7fbf_ffff, .quiet = 0x7fff_ffff },
        .{ .signaling = 0xff80_0123, .quiet = 0xffc0_0123 },
        .{ .signaling = 0xffbf_ffff, .quiet = 0xffff_ffff },
    };
    for (payloads) |payload| {
        try std.testing.expectEqual(
            common.Class.signaling_nan,
            common.F.classifyBits(payload.signaling),
        );
        try std.testing.expectEqual(
            common.Class.quiet_nan,
            common.F.classifyBits(payload.quiet),
        );
        try std.testing.expectEqual(payload.quiet, common.F.quietNaNBits(payload.signaling));
        try std.testing.expectEqual(payload.quiet, common.F.quietNaNBits(payload.quiet));
        try std.testing.expectEqual(payload.quiet, common.F.toBits(common.F.fromBits(payload.quiet)));
    }
    try expectFBits(0x7fc0_0000, common.F.canonicalNaN());

    const normal = common.F.fromBits(0x3f80_0000);
    const signaling = common.F.fromBits(0x7f80_0123);
    const quiet = common.F.fromBits(0x7fc0_0456);
    const first = common.F.firstNaN2(signaling, quiet);
    try std.testing.expect(first != null);
    try expectFBits(0x7fc0_0123, first.?);
    const second = common.F.firstNaN2(normal, signaling);
    try std.testing.expect(second != null);
    try expectFBits(0x7fc0_0123, second.?);
    try std.testing.expect(common.F.firstNaN2(normal, normal) == null);

    const first_three = common.F.firstNaN3(signaling, quiet, normal);
    try std.testing.expect(first_three != null);
    try expectFBits(0x7fc0_0123, first_three.?);
    const second_three = common.F.firstNaN3(normal, signaling, quiet);
    try std.testing.expect(second_three != null);
    try expectFBits(0x7fc0_0123, second_three.?);
    const third_three = common.F.firstNaN3(normal, normal, signaling);
    try std.testing.expect(third_three != null);
    try expectFBits(0x7fc0_0123, third_three.?);
    try std.testing.expect(common.F.firstNaN3(normal, normal, normal) == null);
}

test "f32 exact finite order key endpoints" {
    const cases = [_]struct {
        bits: u32,
        key: u32,
    }{
        .{ .bits = 0xff80_0000, .key = 0x007f_ffff },
        .{ .bits = 0xff7f_ffff, .key = 0x0080_0000 },
        .{ .bits = 0x8080_0000, .key = 0x7f7f_ffff },
        .{ .bits = 0x807f_ffff, .key = 0x7f80_0000 },
        .{ .bits = 0x8000_0000, .key = 0x7fff_ffff },
        .{ .bits = 0x0000_0000, .key = 0x8000_0000 },
        .{ .bits = 0x0000_0001, .key = 0x8000_0001 },
        .{ .bits = 0x007f_ffff, .key = 0x807f_ffff },
        .{ .bits = 0x0080_0000, .key = 0x8080_0000 },
        .{ .bits = 0x7f7f_ffff, .key = 0xff7f_ffff },
        .{ .bits = 0x7f80_0000, .key = 0xff80_0000 },
    };
    for (cases) |case| {
        try std.testing.expectEqual(case.key, common.F.finiteOrderKey(case.bits));
    }
}

test "f64 exact boundary classes and NaN payload edges" {
    const cases = [_]struct {
        bits: u64,
        category: common.Class,
        is_nan: bool,
    }{
        .{ .bits = 0x0000_0000_0000_0000, .category = .zero, .is_nan = false },
        .{ .bits = 0x8000_0000_0000_0000, .category = .zero, .is_nan = false },
        .{ .bits = 0x0000_0000_0000_0001, .category = .subnormal, .is_nan = false },
        .{ .bits = 0x000f_ffff_ffff_ffff, .category = .subnormal, .is_nan = false },
        .{ .bits = 0x0010_0000_0000_0000, .category = .normal, .is_nan = false },
        .{ .bits = 0x7fef_ffff_ffff_ffff, .category = .normal, .is_nan = false },
        .{ .bits = 0x7ff0_0000_0000_0000, .category = .infinity, .is_nan = false },
        .{ .bits = 0xfff0_0000_0000_0000, .category = .infinity, .is_nan = false },
        .{ .bits = 0x7ff8_0000_0000_0000, .category = .quiet_nan, .is_nan = true },
        .{ .bits = 0xfff8_0000_0000_0000, .category = .quiet_nan, .is_nan = true },
        .{ .bits = 0x7ff0_0000_0000_0001, .category = .signaling_nan, .is_nan = true },
        .{ .bits = 0xfff0_0000_0000_0001, .category = .signaling_nan, .is_nan = true },
    };

    for (cases) |case| {
        const value = common.D.fromBits(case.bits);
        try std.testing.expectEqual(case.bits, common.D.toBits(value));
        try std.testing.expectEqual(case.category, common.D.classifyBits(case.bits));
        try std.testing.expectEqual(case.is_nan, common.D.isNaNBits(case.bits));
    }

    const payloads = [_]struct {
        signaling: u64,
        quiet: u64,
    }{
        .{ .signaling = 0x7ff0_0000_0000_0001, .quiet = 0x7ff8_0000_0000_0001 },
        .{ .signaling = 0x7ff0_0000_0000_0123, .quiet = 0x7ff8_0000_0000_0123 },
        .{ .signaling = 0x7ff7_ffff_ffff_ffff, .quiet = 0x7fff_ffff_ffff_ffff },
        .{ .signaling = 0xfff0_0000_0000_0123, .quiet = 0xfff8_0000_0000_0123 },
        .{ .signaling = 0xfff7_ffff_ffff_ffff, .quiet = 0xffff_ffff_ffff_ffff },
    };
    for (payloads) |payload| {
        try std.testing.expectEqual(
            common.Class.signaling_nan,
            common.D.classifyBits(payload.signaling),
        );
        try std.testing.expectEqual(
            common.Class.quiet_nan,
            common.D.classifyBits(payload.quiet),
        );
        try std.testing.expectEqual(payload.quiet, common.D.quietNaNBits(payload.signaling));
        try std.testing.expectEqual(payload.quiet, common.D.quietNaNBits(payload.quiet));
        try std.testing.expectEqual(payload.quiet, common.D.toBits(common.D.fromBits(payload.quiet)));
    }
    try expectDBits(0x7ff8_0000_0000_0000, common.D.canonicalNaN());

    const normal = common.D.fromBits(0x3ff0_0000_0000_0000);
    const signaling = common.D.fromBits(0x7ff0_0000_0000_0123);
    const quiet = common.D.fromBits(0x7ff8_0000_0000_0456);
    const first = common.D.firstNaN2(signaling, quiet);
    try std.testing.expect(first != null);
    try expectDBits(0x7ff8_0000_0000_0123, first.?);
    const second = common.D.firstNaN2(normal, signaling);
    try std.testing.expect(second != null);
    try expectDBits(0x7ff8_0000_0000_0123, second.?);
    try std.testing.expect(common.D.firstNaN2(normal, normal) == null);

    const first_three = common.D.firstNaN3(signaling, quiet, normal);
    try std.testing.expect(first_three != null);
    try expectDBits(0x7ff8_0000_0000_0123, first_three.?);
    const second_three = common.D.firstNaN3(normal, signaling, quiet);
    try std.testing.expect(second_three != null);
    try expectDBits(0x7ff8_0000_0000_0123, second_three.?);
    const third_three = common.D.firstNaN3(normal, normal, signaling);
    try std.testing.expect(third_three != null);
    try expectDBits(0x7ff8_0000_0000_0123, third_three.?);
    try std.testing.expect(common.D.firstNaN3(normal, normal, normal) == null);
}

test "f64 exact finite order key endpoints" {
    const cases = [_]struct {
        bits: u64,
        key: u64,
    }{
        .{ .bits = 0xfff0_0000_0000_0000, .key = 0x000f_ffff_ffff_ffff },
        .{ .bits = 0xffef_ffff_ffff_ffff, .key = 0x0010_0000_0000_0000 },
        .{ .bits = 0x8010_0000_0000_0000, .key = 0x7fef_ffff_ffff_ffff },
        .{ .bits = 0x800f_ffff_ffff_ffff, .key = 0x7ff0_0000_0000_0000 },
        .{ .bits = 0x8000_0000_0000_0000, .key = 0x7fff_ffff_ffff_ffff },
        .{ .bits = 0x0000_0000_0000_0000, .key = 0x8000_0000_0000_0000 },
        .{ .bits = 0x0000_0000_0000_0001, .key = 0x8000_0000_0000_0001 },
        .{ .bits = 0x000f_ffff_ffff_ffff, .key = 0x800f_ffff_ffff_ffff },
        .{ .bits = 0x0010_0000_0000_0000, .key = 0x8010_0000_0000_0000 },
        .{ .bits = 0x7fef_ffff_ffff_ffff, .key = 0xffef_ffff_ffff_ffff },
        .{ .bits = 0x7ff0_0000_0000_0000, .key = 0xfff0_0000_0000_0000 },
    };
    for (cases) |case| {
        try std.testing.expectEqual(case.key, common.D.finiteOrderKey(case.bits));
    }
}
