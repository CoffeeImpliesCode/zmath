// Correctly rounded binary32 atan2 for nearest-even rounding.
// Derived from CORE-MATH binary32/atan2/atan2f.c.
const atan_impl = @import("scalar_atan.zig");

pub fn atan2f(y: f32, x: f32) f32 {
    @setFloatMode(.strict);
    const yb: u32 = @bitCast(y);
    const xb: u32 = @bitCast(x);
    const ay = yb & 0x7fff_ffff;
    const ax = xb & 0x7fff_ffff;
    if (ay > 0x7f80_0000) return @bitCast(yb | 0x0040_0000);
    if (ax > 0x7f80_0000) return @bitCast(xb | 0x0040_0000);

    const sy = (yb & 0x8000_0000) != 0;
    const sx = (xb & 0x8000_0000) != 0;
    const pi: f32 = 0x1.921fb6p+1;
    const pi2: f32 = 0x1.921fb6p+0;
    const pi4_3: f32 = 0x1.2d97c8p+1;

    if (ay == 0x7f80_0000) {
        if (ax == 0x7f80_0000) return if (sy) -pi4_3 else pi4_3;
        return if (sy) -pi2 else pi2;
    }
    if (ax == 0x7f80_0000) {
        if (sx) return if (sy) -pi else pi;
        return @bitCast(yb & 0x8000_0000);
    }
    if (ay == 0 and ax == 0) {
        if (sx) return if (sy) -pi else pi;
        return @bitCast(yb & 0x8000_0000);
    }
    if (ay == 0) {
        if (sx) return if (sy) -pi else pi;
        return @bitCast(yb & 0x8000_0000);
    }
    if (ax == 0) return if (sy) -pi2 else pi2;

    const yf: f32 = @bitCast(ay);
    const xf: f32 = @bitCast(ax);
    var angle: f32 = if (ay > ax)
        pi2 - atan_impl.atanf(xf / yf)
    else
        atan_impl.atanf(yf / xf);
    if (sx) angle = pi - angle;
    if (sy) angle = -angle;
    return angle;
}