const std = @import("std");
const scalar = @import("scalar.zig");

test "atan2 special cases h" {
    const cases = [_][3]u16{
        .{ 0x7c01, 0x0000, 0x7e01 },
        .{ 0x7c00, 0x7c00, 0x40b6 },
        .{ 0xfc00, 0x7c00, 0xc0b6 },
        .{ 0x7c00, 0xfc00, 0x40b6 },
        .{ 0xfc00, 0xfc00, 0xc0b6 },
        .{ 0x3c00, 0x7c00, 0x0000 },
        .{ 0x3c00, 0xfc00, 0x4248 },
        .{ 0xbc00, 0xfc00, 0xc248 },
        .{ 0x0000, 0x0000, 0x0000 },
        .{ 0x8000, 0x0000, 0x8000 },
    };
    for (cases) |c| {
        const y: f16 = @bitCast(c[0]);
        const x: f16 = @bitCast(c[1]);
        const got: u16 = @bitCast(scalar.atan2h(y, x));
        try std.testing.expectEqual(c[2], got);
    }
}

test "atan2 special cases f" {
    const cases = [_][3]u32{
        .{ 0x7f800001, 0x00000000, 0x7fc00001 },
        .{ 0x7f800000, 0x7f800000, 0x4048f5c3 },
        .{ 0xff800000, 0x7f800000, 0xc048f5c3 },
        .{ 0x7f800000, 0xff800000, 0x4048f5c3 },
        .{ 0xff800000, 0xff800000, 0xc048f5c3 },
        .{ 0x40400000, 0x7f800000, 0x00000000 },
        .{ 0x40400000, 0xff800000, 0x4248f5c3 },
        .{ 0xc0400000, 0xff800000, 0xc248f5c3 },
        .{ 0x00000000, 0x00000000, 0x00000000 },
        .{ 0x80000000, 0x00000000, 0x80000000 },
    };
    for (cases) |c| {
        const y: f32 = @bitCast(c[0]);
        const x: f32 = @bitCast(c[1]);
        const got: u32 = @bitCast(scalar.atan2f(y, x));
        try std.testing.expectEqual(c[2], got);
    }
}

test "atan2 special cases d" {
    const cases = [_][3]u64{
        .{ 0x7ff0000000000001, 0x0000000000000000, 0x7ff8000000000001 },
        .{ 0x7ff0000000000000, 0x7ff0000000000000, 0x4002d97c7f3321d2 },
        .{ 0xfff0000000000000, 0x7ff0000000000000, 0xc002d97c7f3321d2 },
        .{ 0x7ff0000000000000, 0xfff0000000000000, 0x4002d97c7f3321d2 },
        .{ 0xfff0000000000000, 0xfff0000000000000, 0xc002d97c7f3321d2 },
        .{ 0x4012e00000000000, 0x7ff0000000000000, 0x0000000000000000 },
        .{ 0x4012e00000000000, 0xfff0000000000000, 0x4212e00000000000 },
        .{ 0xc012e00000000000, 0xfff0000000000000, 0xc212e00000000000 },
        .{ 0x0000000000000000, 0x0000000000000000, 0x0000000000000000 },
        .{ 0x8000000000000000, 0x0000000000000000, 0x8000000000000000 },
    };
    for (cases) |c| {
        const y: f64 = @bitCast(c[0]);
        const x: f64 = @bitCast(c[1]);
        const got: u64 = @bitCast(scalar.atan2d(y, x));
        try std.testing.expectEqual(c[2], got);
    }
}