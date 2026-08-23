// Correctly rounded binary64 atan2 for nearest-even rounding.
// Derived from CORE-MATH binary64/atan2/atan2.c.
const atan_impl = @import("scalar_atan.zig");

pub fn atan2d(y: f64, x: f64) f64 {
    @setFloatMode(.strict);
    const yb: u64 = @bitCast(y);
    const xb: u64 = @bitCast(x);
    const ay = yb & 0x7fff_ffff_ffff_ffff;
    const ax = xb & 0x7fff_ffff_ffff_ffff;
    if (ay > 0x7ff0_0000_0000_0000) return @bitCast(yb | 0x0008_0000_0000_0000);
    if (ax > 0x7ff0_0000_0000_0000) return @bitCast(xb | 0x0008_0000_0000_0000);

    const sy = (yb & 0x8000_0000_0000_0000) != 0;
    const sx = (xb & 0x8000_0000_0000_0000) != 0;
    const pi: f64 = 0x1.921fb54442d18p+1;
    const pi2: f64 = 0x1.921fb54442d18p+0;
    const pi4_3: f64 = 0x1.2d97c7f3321d2p+1;

    if (ay == 0x7ff0_0000_0000_0000) {
        if (ax == 0x7ff0_0000_0000_0000) return if (sy) -pi4_3 else pi4_3;
        return if (sy) -pi2 else pi2;
    }
    if (ax == 0x7ff0_0000_0000_0000) {
        if (sx) return if (sy) -pi else pi;
        return @bitCast(yb & 0x8000_0000_0000_0000);
    }
    if (ay == 0 and ax == 0) {
        if (sx) return if (sy) -pi else pi;
        return @bitCast(yb & 0x8000_0000_0000_0000);
    }
    if (ay == 0) {
        if (sx) return if (sy) -pi else pi;
        return @bitCast(yb & 0x8000_0000_0000_0000);
    }
    if (ax == 0) return if (sy) -pi2 else pi2;

    const yf: f64 = @bitCast(ay);
    const xf: f64 = @bitCast(ax);
    var angle: f64 = if (ay > ax)
        pi2 - atan_impl.atand(xf / yf)
    else
        atan_impl.atand(yf / xf);
    if (sx) angle = pi - angle;
    if (sy) angle = -angle;
    return angle;
}