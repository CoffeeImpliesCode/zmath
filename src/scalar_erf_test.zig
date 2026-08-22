const std = @import("std");
const zmath = @import("root.zig");

fn expectBits(comptime T: type, comptime U: type, expected: U, actual: T) !void {
    try std.testing.expectEqual(expected, @as(U, @bitCast(actual)));
}

test "Erf scalar families expose nearest-even signatures" {
    comptime {
        std.debug.assert(@TypeOf(zmath.scalar.erfh) == fn (f16) f16);
        std.debug.assert(@TypeOf(zmath.scalar.erff) == fn (f32) f32);
        std.debug.assert(@TypeOf(zmath.scalar.erfd) == fn (f64) f64);
    }
}

test "Erf scalar families preserve edge values" {
    try expectBits(f16, u16, 0x0000, zmath.scalar.erfh(@as(f16, 0.0)));
    try expectBits(f16, u16, 0x8000, zmath.scalar.erfh(@as(f16, -0.0)));
    try expectBits(f16, u16, 0x3c00, zmath.scalar.erfh(std.math.inf(f16)));
    try expectBits(f16, u16, 0xbc00, zmath.scalar.erfh(-std.math.inf(f16)));

    try expectBits(f32, u32, 0x0000_0000, zmath.scalar.erff(@as(f32, 0.0)));
    try expectBits(f32, u32, 0x8000_0000, zmath.scalar.erff(@as(f32, -0.0)));
    try expectBits(f32, u32, 0x3f80_0000, zmath.scalar.erff(std.math.inf(f32)));
    try expectBits(f32, u32, 0xbf80_0000, zmath.scalar.erff(-std.math.inf(f32)));

    try expectBits(f64, u64, 0x0000_0000_0000_0000, zmath.scalar.erfd(@as(f64, 0.0)));
    try expectBits(f64, u64, 0x8000_0000_0000_0000, zmath.scalar.erfd(@as(f64, -0.0)));
    try expectBits(f64, u64, 0x3ff0_0000_0000_0000, zmath.scalar.erfd(std.math.inf(f64)));
    try expectBits(f64, u64, 0xbff0_0000_0000_0000, zmath.scalar.erfd(-std.math.inf(f64)));
}

test "Erf scalar families quiet NaNs while preserving sign and payload" {
    try expectBits(f16, u16, 0x7f55, zmath.scalar.erfh(@bitCast(@as(u16, 0x7d55))));
    try expectBits(f16, u16, 0xff55, zmath.scalar.erfh(@bitCast(@as(u16, 0xfd55))));
    try expectBits(f16, u16, 0x7e55, zmath.scalar.erfh(@bitCast(@as(u16, 0x7e55))));
    try expectBits(f16, u16, 0xfe55, zmath.scalar.erfh(@bitCast(@as(u16, 0xfe55))));

    try expectBits(
        f32,
        u32,
        0x7fc1_2345,
        zmath.scalar.erff(@bitCast(@as(u32, 0x7f81_2345))),
    );
    try expectBits(
        f32,
        u32,
        0xffc1_2345,
        zmath.scalar.erff(@bitCast(@as(u32, 0xff81_2345))),
    );
    try expectBits(
        f32,
        u32,
        0x7fc1_2345,
        zmath.scalar.erff(@bitCast(@as(u32, 0x7fc1_2345))),
    );
    try expectBits(
        f32,
        u32,
        0xffc1_2345,
        zmath.scalar.erff(@bitCast(@as(u32, 0xffc1_2345))),
    );
    try expectBits(
        f64,
        u64,
        0x7ff8_1234_5678_9abc,
        zmath.scalar.erfd(@bitCast(@as(u64, 0x7ff0_1234_5678_9abc))),
    );
    try expectBits(
        f64,
        u64,
        0xfff8_1234_5678_9abc,
        zmath.scalar.erfd(@bitCast(@as(u64, 0xfff0_1234_5678_9abc))),
    );
    try expectBits(
        f64,
        u64,
        0x7ff8_1234_5678_9abc,
        zmath.scalar.erfd(@bitCast(@as(u64, 0x7ff8_1234_5678_9abc))),
    );
    try expectBits(
        f64,
        u64,
        0xfff8_1234_5678_9abc,
        zmath.scalar.erfd(@bitCast(@as(u64, 0xfff8_1234_5678_9abc))),
    );
}

test "Erf scalar families round representative values to nearest even" {
    try expectBits(f16, u16, 0x382a, zmath.scalar.erfh(@as(f16, 0.5)));
    try expectBits(f16, u16, 0x3abe, zmath.scalar.erfh(@as(f16, 1.0)));
    try expectBits(f16, u16, 0x3bf6, zmath.scalar.erfh(@as(f16, 2.0)));

    try expectBits(f32, u32, 0x3f05_3f7b, zmath.scalar.erff(@as(f32, 0.5)));
    try expectBits(f32, u32, 0x3f57_bb3d, zmath.scalar.erff(@as(f32, 1.0)));
    try expectBits(f32, u32, 0x3f7e_cd71, zmath.scalar.erff(@as(f32, 2.0)));

    try expectBits(f64, u64, 0x3fe0_a7ef_5c18_edd2, zmath.scalar.erfd(@as(f64, 0.5)));
    try expectBits(f64, u64, 0x3fea_f767_a741_088b, zmath.scalar.erfd(@as(f64, 1.0)));
    try expectBits(f64, u64, 0x3fef_d9ae_1427_95e3, zmath.scalar.erfd(@as(f64, 2.0)));
}

test "Erf scalar families resolve pinned CORE-MATH witnesses" {
    try expectBits(f32, u32, 0x3f7f_ffff, zmath.scalar.erff(@bitCast(@as(u32, 0x407a_d444))));
    try expectBits(f32, u32, 0x3f80_0000, zmath.scalar.erff(@bitCast(@as(u32, 0x407a_d445))));

    try expectBits(
        f64,
        u64,
        0x3c22_4cb3_732e_544b,
        zmath.scalar.erfd(@bitCast(@as(u64, 0x3c20_37b5_48d9_d7a6))),
    );
    try expectBits(
        f64,
        u64,
        0x3c2b_606e_e962_afaa,
        zmath.scalar.erfd(@bitCast(@as(u64, 0x3c28_430f_efa1_5cea))),
    );
    try expectBits(
        f64,
        u64,
        0x3cec_b772_7228_d5e1,
        zmath.scalar.erfd(@bitCast(@as(u64, 0x3ce9_730c_e1d0_f154))),
    );
}

test "Erf scalar families are odd" {
    const h: f16 = 0.75;
    const f: f32 = 0.75;
    const d: f64 = 0.75;
    try std.testing.expectEqual(-zmath.scalar.erfh(h), zmath.scalar.erfh(-h));
    try std.testing.expectEqual(-zmath.scalar.erff(f), zmath.scalar.erff(-f));
    try std.testing.expectEqual(-zmath.scalar.erfd(d), zmath.scalar.erfd(-d));
}

test "generic Erf dispatches without changing type" {
    const h = zmath.erf(@as(f16, 0.5));
    const f = zmath.erf(@as(f32, 0.5));
    const d = zmath.erf(@as(f64, 0.5));
    comptime std.debug.assert(@TypeOf(h) == f16);
    comptime std.debug.assert(@TypeOf(f) == f32);
    comptime std.debug.assert(@TypeOf(d) == f64);
    try std.testing.expectEqual(zmath.scalar.erfh(@as(f16, 0.5)), h);
    try std.testing.expectEqual(zmath.scalar.erff(@as(f32, 0.5)), f);
    try std.testing.expectEqual(zmath.scalar.erfd(@as(f64, 0.5)), d);
}
