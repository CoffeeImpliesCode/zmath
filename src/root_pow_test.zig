const std = @import("std");
const zmath = @import("root.zig");

fn bitsH(value: f16) u16 {
    return @bitCast(value);
}

fn bitsF(value: f32) u32 {
    return @bitCast(value);
}

fn bitsD(value: f64) u64 {
    return @bitCast(value);
}

test "root exports nearest-even Pow H F and D" {
    try std.testing.expectEqual(@as(u16, 0x4800), bitsH(zmath.pow(@as(f16, 2.0), @as(f16, 3.0))));
    try std.testing.expectEqual(@as(u32, 0x4100_0000), bitsF(zmath.pow(@as(f32, 2.0), @as(f32, 3.0))));
    try std.testing.expectEqual(
        @as(u64, 0x4020_0000_0000_0000),
        bitsD(zmath.pow(@as(f64, 2.0), @as(f64, 3.0))),
    );
}

test "root Pow preserves the semantic binary type" {
    const half_result: f16 = zmath.pow(@as(f16, 4.0), @as(f16, 0.5));
    const float_result: f32 = zmath.pow(@as(f32, 4.0), @as(f32, 0.5));
    const double_result: f64 = zmath.pow(@as(f64, 4.0), @as(f64, 0.5));
    try std.testing.expectEqual(@as(u16, 0x4000), bitsH(half_result));
    try std.testing.expectEqual(@as(u32, 0x4000_0000), bitsF(float_result));
    try std.testing.expectEqual(@as(u64, 0x4000_0000_0000_0000), bitsD(double_result));
}
