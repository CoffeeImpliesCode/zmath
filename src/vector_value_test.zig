const std = @import("std");
const vv = @import("vector_value.zig");

const Vec2f = vv.Vec2f;
const Vec3f = vv.Vec3f;
const Vec4f = vv.Vec4f;
const Vec2d = vv.Vec2d;
const Vec3d = vv.Vec3d;
const Vec4d = vv.Vec4d;

/// The positive canonical quiet NaN of each format, from
/// `spec/core-math-2026-08-22.md:769-801`.
const qnan_f: u32 = 0x7fc0_0000;
const qnan_d: u64 = 0x7ff8_0000_0000_0000;

/// The smallest positive subnormal of each format, 2^-149 and 2^-1074.
const subnormal_f: f32 = 0x1p-149;
const subnormal_d: f64 = 0x1p-1074;

/// The smallest positive normal of each format, 2^-126 and 2^-1022.
const normal_f: f32 = 0x1p-126;
const normal_d: f64 = 0x1p-1022;

/// The reference norms the tests compare against.
const root_two: f64 = @sqrt(2.0);
const root_three: f64 = @sqrt(3.0);
const root_two_f: f32 = @sqrt(2.0);

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// The number of components in `V`.
fn widthOf(comptime V: type) usize {
    return @typeInfo(V).@"struct".fields.len;
}

/// The unsigned integer of the same width as `T`.
fn intOf(comptime T: type) type {
    return if (T == f32) u32 else u64;
}

/// The exact bits of a component value.
fn bits(comptime T: type, v: T) intOf(T) {
    return @bitCast(v);
}

/// The exact component bits of `v`. A float comparison cannot tell `+0` from
/// `-0` and cannot compare a NaN at all, so every test that cares about those
/// goes through here.
fn bitsOf(comptime V: type, v: V) [widthOf(V)]intOf(@typeInfo(V).@"struct".fields[0].type) {
    const T = @typeInfo(V).@"struct".fields[0].type;
    const width = comptime widthOf(V);
    var out: [width]intOf(T) = undefined;
    const array = v.toArray();
    inline for (0..width) |index| {
        out[index] = @bitCast(array[index]);
    }
    return out;
}

/// The lanes of a `@Vector` as an array, for a test that compares them.
///
/// This is the coercion from a `@Vector` to an array of the same length and
/// element type, so lane `i` becomes element `i`. It must not be a `@bitCast`:
/// that reinterprets the bytes of the vector, an array has a well-defined byte
/// layout but a vector has none, and on a big-endian target the pun returned
/// the lanes reversed.
fn lanesOf(comptime N: usize, comptime T: type, v: @Vector(N, T)) [N]T {
    const lanes: [N]T = v;
    return lanes;
}

/// The value one ulp above `v`, for a test that needs two distinct neighbours
/// of the same number.
fn nextUp(comptime T: type, v: T) T {
    return @bitCast(@as(intOf(T), @bitCast(v)) + 1);
}

/// True when the product of `a` and `b` overflows to infinity. Taking the
/// operands as parameters keeps the multiplication out of comptime folding, so
/// the test reports what the hardware does at run time.
fn productIsInfinite(comptime T: type, a: T, b: T) bool {
    return std.math.isInf(a * b);
}

/// True when the product of `a` and `b` flushes to zero.
fn productIsZero(comptime T: type, a: T, b: T) bool {
    return a * b == 0.0;
}

/// A normalized vector has a length of one up to the roundings in the division
/// and in the `hypot` fold. Eight eps covers those with room to spare, and is
/// still far tighter than the agreement the fold could offer.
fn expectUnit(comptime V: type, v: V) !void {
    const T = @typeInfo(V).@"struct".fields[0].type;
    try std.testing.expect(@abs(v.normalize().length() - 1.0) <= 8 * std.math.floatEps(T));
}

/// The `Vec*` and `@Vector` bridges move bits, so the round trip is exact even
/// for the values a float comparison cannot see.
fn expectVectorRoundTrip(comptime V: type, v: V) !void {
    try std.testing.expectEqual(bitsOf(V, v), bitsOf(V, V.fromVector(v.toVector())));
}

/// True when `info` mentions `anytype` in its return type or in a parameter,
/// directly or inside a function or optional type.
fn mentionsAnytype(comptime info: std.builtin.Type.Fn) bool {
    if (info.return_type == null) return true;
    inline for (info.params) |param| {
        if (param.is_generic) return true;
        const param_type = param.type orelse return true;
        switch (@typeInfo(param_type)) {
            .@"fn" => |nested| if (mentionsAnytype(nested)) return true,
            .optional => |nested| {
                const child = nested.child orelse continue;
                switch (@typeInfo(child)) {
                    .@"fn" => |inner| if (mentionsAnytype(inner)) return true,
                    else => {},
                }
            },
            else => {},
        }
    }
    return false;
}

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

test "init, fromArray, and toArray agree in x y z w order" {
    try std.testing.expectEqual(Vec2f.init(1, 2), Vec2f.fromArray([2]f32{ 1, 2 }));
    try std.testing.expectEqual(Vec3f.init(1, 2, 3), Vec3f.fromArray([3]f32{ 1, 2, 3 }));
    try std.testing.expectEqual(Vec4f.init(1, 2, 3, 4), Vec4f.fromArray([4]f32{ 1, 2, 3, 4 }));
    try std.testing.expectEqual(Vec2d.init(1, 2), Vec2d.fromArray([2]f64{ 1, 2 }));
    try std.testing.expectEqual(Vec3d.init(1, 2, 3), Vec3d.fromArray([3]f64{ 1, 2, 3 }));
    try std.testing.expectEqual(Vec4d.init(1, 2, 3, 4), Vec4d.fromArray([4]f64{ 1, 2, 3, 4 }));

    try std.testing.expectEqual([2]f32{ 1, 2 }, Vec2f.init(1, 2).toArray());
    try std.testing.expectEqual([3]f32{ 1, 2, 3 }, Vec3f.init(1, 2, 3).toArray());
    try std.testing.expectEqual([4]f32{ 1, 2, 3, 4 }, Vec4f.init(1, 2, 3, 4).toArray());
    try std.testing.expectEqual([2]f64{ 1, 2 }, Vec2d.init(1, 2).toArray());
    try std.testing.expectEqual([3]f64{ 1, 2, 3 }, Vec3d.init(1, 2, 3).toArray());
    try std.testing.expectEqual([4]f64{ 1, 2, 3, 4 }, Vec4d.init(1, 2, 3, 4).toArray());
}

test "the fourth component is w, the homogeneous component" {
    try std.testing.expectEqual(@as(f32, 4), Vec4f.init(1, 2, 3, 4).w);
    try std.testing.expectEqual(@as(f64, 4), Vec4d.init(1, 2, 3, 4).w);
    try std.testing.expectEqual(@as(f32, 4), Vec4f.fromArray([4]f32{ 1, 2, 3, 4 }).w);
    // w is the fourth lane, and 4.0 is 0x4080_0000 in binary32 and
    // 0x4010_0000_0000_0000 in binary64.
    try std.testing.expectEqual(@as(u32, 0x4080_0000), bitsOf(Vec4f, Vec4f.init(1, 2, 3, 4))[3]);
    try std.testing.expectEqual(
        @as(u64, 0x4010_0000_0000_0000),
        bitsOf(Vec4d, Vec4d.init(1, 2, 3, 4))[3],
    );
}

test "zero is the additive identity" {
    try std.testing.expect(Vec2d.init(1, -2).add(Vec2d.zero()).eq(Vec2d.init(1, -2)));
    try std.testing.expect(Vec2f.init(1, -2).add(Vec2f.zero()).eq(Vec2f.init(1, -2)));
    try std.testing.expect(Vec3d.init(1, -2, 3).add(Vec3d.zero()).eq(Vec3d.init(1, -2, 3)));
    try std.testing.expect(Vec3f.init(1, -2, 3).add(Vec3f.zero()).eq(Vec3f.init(1, -2, 3)));
    try std.testing.expect(Vec4d.init(1, -2, 3, -4).add(Vec4d.zero()).eq(Vec4d.init(1, -2, 3, -4)));
    try std.testing.expect(Vec4f.init(1, -2, 3, -4).add(Vec4f.zero()).eq(Vec4f.init(1, -2, 3, -4)));

    // The identity works in reverse too, and zero really is every component.
    try std.testing.expect(
        Vec4d.zero().add(Vec4d.init(1, -2, 3, -4)).eq(Vec4d.init(1, -2, 3, -4)),
    );
    try std.testing.expectEqual([4]u64{ 0, 0, 0, 0 }, bitsOf(Vec4d, Vec4d.zero()));
    try std.testing.expectEqual([2]u32{ 0, 0 }, bitsOf(Vec2f, Vec2f.zero()));
}

test "one is the scaleVec identity and not a unit vector" {
    try std.testing.expect(Vec2d.one().scaleVec(Vec2d.init(3, 4)).eq(Vec2d.init(3, 4)));
    try std.testing.expect(Vec2f.one().scaleVec(Vec2f.init(3, 4)).eq(Vec2f.init(3, 4)));
    try std.testing.expect(Vec3d.one().scaleVec(Vec3d.init(3, 4, 5)).eq(Vec3d.init(3, 4, 5)));
    try std.testing.expect(Vec3f.one().scaleVec(Vec3f.init(3, 4, 5)).eq(Vec3f.init(3, 4, 5)));
    try std.testing.expect(Vec4d.one().scaleVec(Vec4d.init(3, 4, 5, 6)).eq(Vec4d.init(3, 4, 5, 6)));
    try std.testing.expect(Vec4f.one().scaleVec(Vec4f.init(9, 8, 7, 6)).eq(Vec4f.init(9, 8, 7, 6)));

    try std.testing.expect(Vec2d.zero().isFinite());
    try std.testing.expect(Vec2d.one().isFinite());
    try std.testing.expect(Vec3f.zero().isFinite());

    // The length of the all-ones vector is sqrt(N), so it is 1 only for N = 1.
    try std.testing.expectEqual(@as(f64, 2), Vec2d.one().lengthSquared());
    try std.testing.expectEqual(@as(f64, 3), Vec3d.one().lengthSquared());
    try std.testing.expectEqual(@as(f64, 4), Vec4d.one().lengthSquared());
    try std.testing.expectEqual(@as(f32, 3), Vec3f.one().lengthSquared());

    try std.testing.expect(@abs(Vec2d.one().length() - root_two) <= 2 * std.math.floatEps(f64));
    try std.testing.expect(@abs(Vec3d.one().length() - root_three) <= 4 * std.math.floatEps(f64));
    try std.testing.expect(@abs(Vec4d.one().length() - 2.0) <= 4 * std.math.floatEps(f64));
}

// ---------------------------------------------------------------------------
// Arithmetic
// ---------------------------------------------------------------------------

test "add, sub, neg, scale, scaleVec, and divScalar are componentwise" {
    try std.testing.expectEqual(Vec3d.init(5, -1, 9), Vec3d.init(1, 2, 3).add(Vec3d.init(4, -3, 6)));
    try std.testing.expectEqual(Vec3d.init(-3, -1, -3), Vec3d.init(1, 2, 3).sub(Vec3d.init(4, 3, 6)));
    try std.testing.expectEqual(Vec3d.init(-1, -2, -3), Vec3d.init(1, 2, 3).neg());
    try std.testing.expectEqual(Vec3d.init(2, 4, 6), Vec3d.init(1, 2, 3).scale(2));
    try std.testing.expectEqual(Vec3d.init(4, 6, 6), Vec3d.init(1, 2, 3).scaleVec(Vec3d.init(4, 3, 2)));
    try std.testing.expectEqual(Vec3d.init(0.5, 1, 1.5), Vec3d.init(1, 2, 3).divScalar(2));

    try std.testing.expectEqual(Vec2f.init(0.5, 1.5), Vec2f.init(1, 3).divScalar(2));
    try std.testing.expectEqual(Vec4f.init(-1, -2, -3, -4), Vec4f.init(1, 2, 3, 4).neg());
    try std.testing.expectEqual(
        Vec4d.init(1, 4, 9, 16),
        Vec4d.init(1, 2, 3, 4).scaleVec(Vec4d.init(1, 2, 3, 4)),
    );
    try std.testing.expectEqual(Vec2d.init(6, 8), Vec2d.init(3, 4).scale(2));
    // Subtracting is adding the negation, component for component.
    try std.testing.expect(
        Vec4d.init(1, -2, 3, -4).sub(Vec4d.init(4, 3, 6, -1)).eq(
            Vec4d.init(1, -2, 3, -4).add(Vec4d.init(4, 3, 6, -1).neg()),
        ),
    );
}

test "negation and abs clear the sign of a negative zero" {
    try std.testing.expectEqual([2]u32{ 0, 0xbf80_0000 }, bitsOf(Vec2f, Vec2f.init(-0.0, 1.0).neg()));
    try std.testing.expectEqual(
        [2]u64{ 0, 0xbff0_0000_0000_0000 },
        bitsOf(Vec2d, Vec2d.init(-0.0, 1.0).neg()),
    );
    try std.testing.expectEqual(
        [3]u32{ 0x4040_0000, 0x4040_0000, 0 },
        bitsOf(Vec3f, Vec3f.init(-3, 3, -0.0).abs()),
    );
    try std.testing.expectEqual(
        [4]u64{ 0x4008_0000_0000_0000, 0x4008_0000_0000_0000, 0, 0 },
        bitsOf(Vec4d, Vec4d.init(3, -3, -0.0, 0.0).abs()),
    );
    try std.testing.expect(std.math.isNan(Vec2d.init(std.math.nan(f64), 0.0).abs().x));
    try std.testing.expect(std.math.isNan(Vec3f.init(-std.math.nan(f32), 0.0, 0.0).abs().x));
}

test "an invalid non-NaN input gives the positive canonical quiet NaN" {
    const inf_f: f32 = std.math.inf(f32);
    const inf_d: f64 = std.math.inf(f64);
    const one_f: u32 = 0x3f80_0000;
    const one_d: u64 = 0x3ff0_0000_0000_0000;

    // 0 / 0 is invalid, so every component of the zero vector divided by zero
    // is the positive canonical quiet NaN. The host unit alone would give the
    // negative indefinite here, so this pins the substitution.
    try std.testing.expectEqual(
        [2]u32{ qnan_f, qnan_f },
        bitsOf(Vec2f, Vec2f.zero().divScalar(0.0)),
    );
    try std.testing.expectEqual(
        [3]u64{ qnan_d, qnan_d, qnan_d },
        bitsOf(Vec3d, Vec3d.zero().divScalar(0.0)),
    );
    // infinity * 0 is invalid, and the other components stay as they were.
    // 1.0 is 0x3f80_0000 in binary32.
    try std.testing.expectEqual(
        [4]u32{ qnan_f, 0, 0, 0 },
        bitsOf(Vec4f, Vec4f.init(inf_f, 0, 0, 0).scale(0.0)),
    );
    try std.testing.expectEqual(
        [4]u32{ 0, 0, 0, 0 },
        bitsOf(Vec4f, Vec4f.init(1, 0, 0, 0).scale(0.0)),
    );
    // infinity - infinity and infinity + (-infinity) are invalid.
    try std.testing.expectEqual(
        [2]u32{ qnan_f, one_f },
        bitsOf(Vec2f, Vec2f.init(inf_f, 1.0).sub(Vec2f.init(inf_f, 0.0))),
    );
    try std.testing.expectEqual(
        [2]u64{ qnan_d, one_d },
        bitsOf(Vec2d, Vec2d.init(inf_d, 1.0).add(Vec2d.init(-inf_d, 0.0))),
    );
    // 0 * infinity is invalid from either side, and a componentwise product
    // does not hide it.
    try std.testing.expectEqual(
        [2]u32{ qnan_f, 0 },
        bitsOf(Vec2f, Vec2f.init(0.0, 0.0).scaleVec(Vec2f.init(inf_f, 1.0))),
    );
    // infinity / infinity is invalid.
    try std.testing.expectEqual(
        [2]u64{ qnan_d, 0 },
        bitsOf(Vec2d, Vec2d.init(inf_d, 0.0).divScalar(inf_d)),
    );
    // A dot whose terms are infinite and opposite gives the same NaN, and the
    // NaN then absorbs the terms that follow.
    try std.testing.expectEqual(
        @as(u64, qnan_d),
        bits(f64, Vec3d.init(inf_d, 1.0, 1.0).dot(Vec3d.init(1.0, -inf_d, 1.0))),
    );
    // A cross product whose difference is infinite minus infinite gives the same
    // NaN, in that component only.
    const crossed = Vec3d.init(1.0, 1.0, 1.0).cross(Vec3d.init(1.0, inf_d, inf_d));
    try std.testing.expectEqual(@as(u64, qnan_d), bits(f64, crossed.x));
    // An invalid product is caught before the difference ever sees it:
    // 0 * infinity is invalid, so the first component is that same NaN.
    const crossed2 = Vec3d.init(1.0, 0.0, 1.0).cross(Vec3d.init(1.0, 1.0, inf_d));
    try std.testing.expectEqual(@as(u64, qnan_d), bits(f64, crossed2.x));
}

test "a NaN input propagates with its quiet bit set" {
    // A signaling NaN becomes quiet and keeps its sign and payload.
    const signaling: f64 = @bitCast(@as(u64, 0x7ff0_0000_0000_0001));
    try std.testing.expectEqual(
        [2]u64{ 0x7ff8_0000_0000_0001, 0 },
        bitsOf(Vec2d, Vec2d.init(signaling, 0.0).add(Vec2d.zero())),
    );
    // A negative quiet NaN keeps its sign and payload.
    const negative: f64 = @bitCast(@as(u64, 0xfff8_0000_0000_0000));
    try std.testing.expectEqual(
        [2]u64{ 0xfff8_0000_0000_0000, 0 },
        bitsOf(Vec2d, Vec2d.init(negative, 0.0).add(Vec2d.zero())),
    );
    const signaling32: f32 = @bitCast(@as(u32, 0x7f80_0001));
    try std.testing.expectEqual(
        [3]u32{ 0x7fc0_0001, 0, 0 },
        bitsOf(Vec3f, Vec3f.init(signaling32, 0.0, 0.0).add(Vec3f.zero())),
    );
    // One NaN component does not spoil the others. 1.0 + 1.0 is 2.0, which is
    // 0x4000_0000_0000_0000.
    try std.testing.expectEqual(
        [2]u64{ qnan_d, 0x4000_0000_0000_0000 },
        bitsOf(Vec2d, Vec2d.init(std.math.nan(f64), 1.0).add(Vec2d.init(1.0, 1.0))),
    );
}

// ---------------------------------------------------------------------------
// dot
// ---------------------------------------------------------------------------

test "dot sums products left to right" {
    // (1, 2, 3) . (4, 5, 6) = 4 + 10 + 18 = 32. Flipping the sign of one
    // component of the right operand subtracts twice that product, which gives
    // all four sign combinations of the hand value.
    try std.testing.expectEqual(@as(f64, 32), Vec3d.init(1, 2, 3).dot(Vec3d.init(4, 5, 6)));
    try std.testing.expectEqual(@as(f64, 24), Vec3d.init(1, 2, 3).dot(Vec3d.init(-4, 5, 6)));
    try std.testing.expectEqual(@as(f64, 12), Vec3d.init(1, 2, 3).dot(Vec3d.init(4, -5, 6)));
    try std.testing.expectEqual(@as(f64, 4), Vec3d.init(1, 2, 3).dot(Vec3d.init(-4, -5, 6)));

    try std.testing.expectEqual(@as(f32, 32), Vec3f.init(1, 2, 3).dot(Vec3f.init(4, 5, 6)));
    try std.testing.expectEqual(@as(f32, 24), Vec3f.init(1, 2, 3).dot(Vec3f.init(-4, 5, 6)));
    try std.testing.expectEqual(@as(f32, 12), Vec3f.init(1, 2, 3).dot(Vec3f.init(4, -5, 6)));
    try std.testing.expectEqual(@as(f32, 4), Vec3f.init(1, 2, 3).dot(Vec3f.init(-4, -5, 6)));

    // The other four types share the same kernel and the same order.
    try std.testing.expectEqual(@as(f64, 19), Vec2d.init(1, 2).dot(Vec2d.init(3, 8)));
    try std.testing.expectEqual(@as(f64, -13), Vec2d.init(1, 2).dot(Vec2d.init(3, -8)));
    try std.testing.expectEqual(@as(f32, 19), Vec2f.init(1, 2).dot(Vec2f.init(3, 8)));
    // 1 + 4 + 9 + 16 = 30.
    try std.testing.expectEqual(@as(f64, 30), Vec4d.init(1, 2, 3, 4).dot(Vec4d.init(1, 2, 3, 4)));
    try std.testing.expectEqual(@as(f32, 30), Vec4f.init(1, 2, 3, 4).dot(Vec4f.init(1, 2, 3, 4)));
    // The fourth term pairs with the fourth component, not the first.
    try std.testing.expectEqual(@as(f64, 1), Vec4d.init(1, 2, 3, 4).dot(Vec4d.init(1, 0, 0, 0)));
    try std.testing.expectEqual(@as(f64, 4), Vec4d.init(1, 2, 3, 4).dot(Vec4d.init(0, 0, 0, 1)));
}

test "dot keeps the summation order, and the order changes the answer" {
    // One unit of the last place at 1e30 is about 1.4e4, so 1e30 + 1 rounds
    // back to 1e30. The left-to-right sum of (1e30, 1, -1e30) against
    // (1, 1, 1) is therefore 1e30 - 1e30 = 0. A compiler that grouped the
    // first and last terms first would answer 0 + 1 = 1, so 0 rules it out.
    try std.testing.expectEqual(@as(f64, 0.0), Vec3d.init(1e30, 1, -1e30).dot(Vec3d.one()));
    try std.testing.expectEqual(@as(f32, 0.0), Vec3f.init(1e20, 1, -1e20).dot(Vec3f.one()));
    try std.testing.expectEqual(@as(f64, 0.0), Vec4d.init(1e30, 1, -1e30, 0).dot(Vec4d.one()));
    try std.testing.expectEqual(@as(f32, 0.0), Vec4f.init(1e20, 1, -1e20, 0).dot(Vec4f.one()));

    // The declaration order, not the reverse: in (1e30, -1e30, 1) the first two
    // terms cancel before the third is added, so the answer is 1.
    try std.testing.expectEqual(@as(f64, 1.0), Vec3d.init(1e30, -1e30, 1).dot(Vec3d.one()));
    try std.testing.expectEqual(@as(f32, 1.0), Vec3f.init(1e20, -1e20, 1).dot(Vec3f.one()));
}

test "dot is symmetric and negates with an operand" {
    const a = Vec3d.init(1e-8, 2, -3);
    const b = Vec3d.init(-4, 5, 6);
    // Each product is the same in either order, and the same terms added in
    // the same order give the same bits.
    try std.testing.expectEqual(bits(f64, a.dot(b)), bits(f64, b.dot(a)));
    // Negating one operand negates every product, and negating a sum is
    // exact, so the result tracks bit for bit.
    try std.testing.expectEqual(bits(f64, a.dot(b.neg())), bits(f64, a.neg().dot(b)));
    try std.testing.expectEqual(bits(f64, a.neg().dot(b.neg())), bits(f64, a.dot(b)));
}

test "lengthSquared is the dot product with itself" {
    try std.testing.expectEqual(@as(f64, 14), Vec3d.init(1, 2, 3).lengthSquared());
    try std.testing.expectEqual(@as(f64, 5), Vec2d.init(1, 2).lengthSquared());
    try std.testing.expectEqual(@as(f32, 30), Vec4f.init(1, 2, 3, 4).lengthSquared());
    try std.testing.expectEqual(
        Vec3d.init(1, 2, 3).dot(Vec3d.init(1, 2, 3)),
        Vec3d.init(1, 2, 3).lengthSquared(),
    );
    try std.testing.expectEqual(
        Vec4f.init(1, 2, 3, 4).dot(Vec4f.init(1, 2, 3, 4)),
        Vec4f.init(1, 2, 3, 4).lengthSquared(),
    );
    // A zero vector has a positive zero squared length, not a negative one.
    try std.testing.expectEqual(@as(u32, 0), bits(f32, Vec2f.init(-0.0, -0.0).lengthSquared()));
}

// ---------------------------------------------------------------------------
// cross
// ---------------------------------------------------------------------------

test "cross is the right-handed cross product" {
    // (1, 2, 3) x (4, 5, 6) = (2*6 - 3*5, 3*4 - 1*6, 1*5 - 2*4) = (-3, 6, -3).
    try std.testing.expectEqual(
        Vec3d.init(-3, 6, -3),
        Vec3d.init(1, 2, 3).cross(Vec3d.init(4, 5, 6)),
    );
    try std.testing.expectEqual(
        Vec3f.init(-3, 6, -3),
        Vec3f.init(1, 2, 3).cross(Vec3f.init(4, 5, 6)),
    );
    // x cross y is z, so the basis is right-handed, and the other two follow.
    try std.testing.expectEqual(
        Vec3d.init(0, 0, 1),
        Vec3d.init(1, 0, 0).cross(Vec3d.init(0, 1, 0)),
    );
    try std.testing.expectEqual(
        Vec3d.init(1, 0, 0),
        Vec3d.init(0, 1, 0).cross(Vec3d.init(0, 0, 1)),
    );
    try std.testing.expectEqual(
        Vec3d.init(0, 1, 0),
        Vec3d.init(0, 0, 1).cross(Vec3d.init(1, 0, 0)),
    );
    try std.testing.expectEqual(
        Vec3f.init(0, 0, 1),
        Vec3f.init(1, 0, 0).cross(Vec3f.init(0, 1, 0)),
    );
    // A vector crossed with itself is a positive zero in every component.
    try std.testing.expectEqual(
        [3]u64{ 0, 0, 0 },
        bitsOf(Vec3d, Vec3d.init(1, 2, 3).cross(Vec3d.init(1, 2, 3))),
    );
}

test "cross is anti-commutative" {
    const pairs = [_][2]Vec3d{
        .{ Vec3d.init(1, 2, 3), Vec3d.init(4, 5, 6) },
        .{ Vec3d.init(-1, 0.5, 7), Vec3d.init(2, -3, 0.25) },
        .{ Vec3d.init(1e-8, 2e-8, 3e-8), Vec3d.init(-4e-8, 5e-8, -6e-8) },
    };
    for (pairs) |pair| {
        const a = pair[0];
        const b = pair[1];
        // Negation is exact, so the reversal is a bit-for-bit match.
        try std.testing.expectEqual(bitsOf(Vec3d, b.cross(a)), bitsOf(Vec3d, a.cross(b).neg()));
    }
}

test "cross is orthogonal to both inputs" {
    // Integer components cancel exactly.
    const a = Vec3d.init(1, 2, 3);
    const b = Vec3d.init(4, 5, 6);
    const c = a.cross(b);
    try std.testing.expectEqual(@as(f64, 0.0), c.dot(a));
    try std.testing.expectEqual(@as(f64, 0.0), c.dot(b));

    // Scaled components only cancel up to rounding. Each term of c . a has a
    // magnitude of about |a|^2 * |b|, and six roundings at 8 eps bound the
    // residual of a three-term sum.
    const s = Vec3d.init(1e-8, 2e-8, 3e-8);
    const t = Vec3d.init(-4e-8, 5e-8, -6e-8);
    const d = s.cross(t);
    const bound = 8 * std.math.floatEps(f64) * s.length() * s.length() * t.length();
    try std.testing.expect(@abs(d.dot(s)) <= bound);
    try std.testing.expect(@abs(d.dot(t)) <= bound);
}

test "cross exists only in three dimensions" {
    try std.testing.expect(@hasDecl(Vec3f, "cross"));
    try std.testing.expect(@hasDecl(Vec3d, "cross"));
    try std.testing.expect(!@hasDecl(Vec2f, "cross"));
    try std.testing.expect(!@hasDecl(Vec2d, "cross"));
    try std.testing.expect(!@hasDecl(Vec4f, "cross"));
    try std.testing.expect(!@hasDecl(Vec4d, "cross"));
}

// ---------------------------------------------------------------------------
// length, distance, and normalize
// ---------------------------------------------------------------------------

test "length does not overflow where a sum of squares does" {
    // 1e200 squared is past the top of binary64, so `sqrt(x*x + y*y)` loses the
    // answer to overflow while the `hypot` fold keeps it.
    try std.testing.expect(productIsInfinite(f64, 1e200, 1e200));

    const two = Vec2d.init(1e200, 1e200).length();
    try std.testing.expect(std.math.isFinite(two));
    try std.testing.expectApproxEqRel(root_two * 1e200, two, 4 * std.math.floatEps(f64));
    try std.testing.expectApproxEqRel(
        root_three * 1e200,
        Vec3d.init(1e200, 1e200, 1e200).length(),
        8 * std.math.floatEps(f64),
    );
    try std.testing.expectApproxEqRel(
        @as(f64, 2e200),
        Vec4d.init(1e200, 1e200, 1e200, 1e200).length(),
        8 * std.math.floatEps(f64),
    );

    try std.testing.expect(productIsInfinite(f32, 1e30, 1e30));
    try std.testing.expectApproxEqRel(
        @as(f32, 2e30),
        Vec4f.init(1e30, 1e30, 1e30, 1e30).length(),
        8 * std.math.floatEps(f32),
    );
    try std.testing.expectApproxEqRel(
        root_two_f * 1e30,
        Vec2f.init(1e30, 1e30).length(),
        4 * std.math.floatEps(f32),
    );
}

test "length does not flush a subnormal vector to zero" {
    // The accepted kernel returns one ulp for hypot of one ulp against one ulp,
    // which `scalar_hypot_test.zig` pins for both formats. The fold therefore
    // keeps a norm that a sum of squares loses, because one ulp squared is
    // below the bottom of the range.
    try std.testing.expectEqual(
        @as(f64, subnormal_d),
        Vec2d.init(subnormal_d, subnormal_d).length(),
    );
    try std.testing.expectEqual(
        @as(f32, subnormal_f),
        Vec2f.init(subnormal_f, subnormal_f).length(),
    );
    try std.testing.expectEqual(
        @as(f64, subnormal_d),
        Vec3d.init(subnormal_d, 0, 0).length(),
    );
    try std.testing.expectEqual(
        @as(f32, subnormal_f),
        Vec4f.init(subnormal_f, 0, 0, 0).length(),
    );
    try std.testing.expect(productIsZero(f64, subnormal_d, subnormal_d));

    // The smallest normal has a norm that is still a normal.
    try std.testing.expectApproxEqRel(
        root_two * normal_d,
        Vec2d.init(normal_d, normal_d).length(),
        4 * std.math.floatEps(f64),
    );
    try std.testing.expectApproxEqRel(
        root_two_f * normal_f,
        Vec2f.init(normal_f, normal_f).length(),
        4 * std.math.floatEps(f32),
    );
}

test "length ignores the sign of every component" {
    try std.testing.expectEqual(
        Vec3d.init(3, 4, 0).length(),
        Vec3d.init(-3, -4, -0.0).length(),
    );
    try std.testing.expectEqual(@as(f32, 5.0), Vec2f.init(-3, -4).length());
    // A zero norm is a positive zero whatever the signs of its components.
    try std.testing.expectEqual(@as(u64, 0), bits(f64, Vec2d.init(-0.0, -0.0).length()));
}

test "the first component enters the length fold as its exact magnitude" {
    // A vector whose other components are zero has the magnitude of its first
    // component, bit for bit: the fold cannot round it and the sign is gone.
    // That is what pins a subnormal, a normal, an overflow-scale value, and
    // both signs of each of them.
    const magnitudes_f = [_]u32{
        0x0000_0001, 0x007f_ffff, 0x0080_0000, 0x3f80_0000, 0x7f7f_ffff, 0x7f80_0000,
        0x8000_0001, 0x807f_ffff, 0x8080_0000, 0xbf80_0000, 0xff7f_ffff, 0xff80_0000,
    };
    for (magnitudes_f) |magnitude| {
        const x: f32 = @bitCast(magnitude);
        const cleared: u32 = magnitude & 0x7fff_ffff;
        try std.testing.expectEqual(cleared, bits(f32, Vec2f.init(x, 0.0).length()));
        try std.testing.expectEqual(cleared, bits(f32, Vec4f.init(x, 0.0, 0.0, 0.0).length()));
    }
    const magnitudes_d = [_]u64{
        0x0000_0000_0000_0001, 0x000f_ffff_ffff_ffff, 0x0010_0000_0000_0000,
        0x3ff0_0000_0000_0000, 0x7fef_ffff_ffff_ffff, 0x7ff0_0000_0000_0000,
        0x8000_0000_0000_0001, 0x800f_ffff_ffff_ffff, 0x8010_0000_0000_0000,
        0xbff0_0000_0000_0000, 0xffef_ffff_ffff_ffff, 0xfff0_0000_0000_0000,
    };
    for (magnitudes_d) |magnitude| {
        const x: f64 = @bitCast(magnitude);
        const cleared: u64 = magnitude & 0x7fff_ffff_ffff_ffff;
        try std.testing.expectEqual(cleared, bits(f64, Vec2d.init(x, 0.0).length()));
        try std.testing.expectEqual(cleared, bits(f64, Vec3d.init(x, 0.0, 0.0).length()));
    }
}

test "length quiets a signaling NaN and keeps its sign and payload" {
    // A signaling NaN in the first component comes back with its quiet bit set
    // and its sign and payload intact, so a caller reading the bits of a
    // non-finite norm sees the operand it passed in rather than a canonical
    // NaN. A signaling NaN must not survive into the result.
    const signaling32: f32 = @bitCast(@as(u32, 0x7f80_0001));
    try std.testing.expectEqual(
        @as(u32, 0x7fc0_0001),
        bits(f32, Vec3f.init(signaling32, 1.0, 1.0).length()),
    );
    const negative32: f32 = @bitCast(@as(u32, 0xff80_0001));
    try std.testing.expectEqual(
        @as(u32, 0xffc0_0001),
        bits(f32, Vec2f.init(negative32, 1.0).length()),
    );
    const signaling: f64 = @bitCast(@as(u64, 0x7ff0_0000_0000_0001));
    try std.testing.expectEqual(
        @as(u64, 0x7ff8_0000_0000_0001),
        bits(f64, Vec4d.init(signaling, 1.0, 1.0, 1.0).length()),
    );
    const negative: f64 = @bitCast(@as(u64, 0xfff0_0000_0000_0001));
    try std.testing.expectEqual(
        @as(u64, 0xfff8_0000_0000_0001),
        bits(f64, Vec2d.init(negative, 1.0).length()),
    );
    // A quiet NaN keeps its own sign, which the magnitude rule above cannot
    // express: a NaN is not `|x|`.
    const quiet: f64 = @bitCast(@as(u64, 0xfff8_0000_0000_0000));
    try std.testing.expectEqual(
        @as(u64, 0xfff8_0000_0000_0000),
        bits(f64, Vec2d.init(quiet, 0.0).length()),
    );
}

test "length follows the kernel's precedence between an infinity and a NaN" {
    // The accepted `hypot` kernel lets a signaling NaN beat an infinity but
    // loses to one, and a quiet NaN loses to an infinity in either position.
    // The length is that kernel applied repeatedly, so which rule applies
    // depends on the order of the components.
    const signaling32: f32 = @bitCast(@as(u32, 0x7f80_0001));
    const quiet32: f32 = @bitCast(@as(u32, 0x7fc0_0000));
    // A signaling NaN in the first position is already quiet when the infinity
    // arrives, and a quiet NaN loses to an infinity.
    try std.testing.expectEqual(
        @as(u32, 0x7f80_0000),
        bits(f32, Vec2f.init(signaling32, std.math.inf(f32)).length()),
    );
    // An infinity first and a signaling NaN second: the later signaling NaN
    // wins, and it arrives quieted.
    try std.testing.expectEqual(
        @as(u32, 0x7fc0_0001),
        bits(f32, Vec2f.init(std.math.inf(f32), signaling32).length()),
    );
    // An infinity first and a quiet NaN second: the infinity wins.
    try std.testing.expectEqual(
        @as(u32, 0x7f80_0000),
        bits(f32, Vec2f.init(std.math.inf(f32), quiet32).length()),
    );
    const signaling: f64 = @bitCast(@as(u64, 0x7ff0_0000_0000_0001));
    const quiet: f64 = @bitCast(@as(u64, 0x7ff8_0000_0000_0000));
    try std.testing.expectEqual(
        @as(u64, 0x7ff0_0000_0000_0000),
        bits(f64, Vec2d.init(signaling, std.math.inf(f64)).length()),
    );
    try std.testing.expectEqual(
        @as(u64, 0x7ff8_0000_0000_0001),
        bits(f64, Vec2d.init(std.math.inf(f64), signaling).length()),
    );
    try std.testing.expectEqual(
        @as(u64, 0x7ff0_0000_0000_0000),
        bits(f64, Vec2d.init(std.math.inf(f64), quiet).length()),
    );
    // Three components fold twice, so a NaN in the first position is already
    // the accumulator and the infinity only has to survive one more step.
    try std.testing.expectEqual(
        @as(u64, 0x7ff0_0000_0000_0000),
        bits(f64, Vec3d.init(signaling, std.math.inf(f64), 1.0).length()),
    );
}

test "a vector of signed zeros has a positive zero length" {
    try std.testing.expectEqual(@as(u32, 0), bits(f32, Vec2f.init(-0.0, 0.0).length()));
    try std.testing.expectEqual(@as(u32, 0), bits(f32, Vec3f.init(-0.0, -0.0, -0.0).length()));
    try std.testing.expectEqual(@as(u64, 0), bits(f64, Vec4d.init(-0.0, -0.0, -0.0, -0.0).length()));
    // A signed zero in the first position does not decide the sign of a norm
    // that the other components make positive.
    try std.testing.expectEqual(@as(u32, 0x3f80_0000), bits(f32, Vec3f.init(-0.0, 1.0, 0.0).length()));
    try std.testing.expectEqual(@as(u64, 0x3ff0_0000_0000_0000), bits(f64, Vec4d.init(-0.0, -0.0, 1.0, -0.0).length()));
}

test "distance inherits the quieting of the first component" {
    // The distance is the length of the rounded difference, and the difference
    // propagates the NaN with its sign cleared, so both directions of the
    // distance answer the same quieted NaN.
    const signaling: f64 = @bitCast(@as(u64, 0x7ff0_0000_0000_0001));
    const quieted: u64 = 0x7ff8_0000_0000_0001;
    try std.testing.expectEqual(
        quieted,
        bits(f64, Vec2d.init(signaling, 1.0).distance(Vec2d.zero())),
    );
    try std.testing.expectEqual(
        quieted,
        bits(f64, Vec2d.zero().distance(Vec2d.init(signaling, 1.0))),
    );
    const signaling32: f32 = @bitCast(@as(u32, 0x7f80_0001));
    try std.testing.expectEqual(
        @as(u32, 0x7fc0_0001),
        bits(f32, Vec3f.init(signaling32, 1.0, 1.0).distance(Vec3f.zero())),
    );
    // The normalizer divides by the length, so a quieted first component is
    // the divisor and the result is not finite.
    const unit = Vec2f.init(signaling32, 0.0).normalize();
    try std.testing.expect(!unit.isFinite());
    try std.testing.expect(std.math.isNan(unit.x));
    try std.testing.expect(std.math.isNan(unit.y));
}

test "normalize returns a unit vector" {
    try expectUnit(Vec2d, Vec2d.init(3, 4));
    try expectUnit(Vec3d, Vec3d.init(1, 2, 3));
    try expectUnit(Vec4d, Vec4d.init(0.1, 0.2, 0.3, 0.4));
    try expectUnit(Vec2f, Vec2f.init(-3, 4));
    try expectUnit(Vec3f, Vec3f.init(1, 2, 3));
    try expectUnit(Vec4f, Vec4f.init(0.1, 0.2, 0.3, 0.4));
    // The smallest normal is nowhere near zero, so it is still a direction.
    try expectUnit(Vec2d, Vec2d.init(normal_d, normal_d));
    try expectUnit(Vec2f, Vec2f.init(normal_f, normal_f));

    // The direction survives: a scalar multiple of a vector normalizes to the
    // same direction, so the two unit vectors are parallel.
    const parallel = Vec3d.init(1, 2, 3).normalize().cross(Vec3d.init(2, 4, 6).normalize());
    const bound = 8 * std.math.floatEps(f64);
    try std.testing.expect(@abs(parallel.x) <= bound);
    try std.testing.expect(@abs(parallel.y) <= bound);
    try std.testing.expect(@abs(parallel.z) <= bound);
}

test "normalize of zero returns the input unchanged" {
    try std.testing.expectEqual([2]u64{ 0, 0 }, bitsOf(Vec2d, Vec2d.zero().normalize()));
    try std.testing.expectEqual([2]u32{ 0, 0 }, bitsOf(Vec2f, Vec2f.zero().normalize()));
    try std.testing.expectEqual(
        [3]u32{ 0, 0x8000_0000, 0 },
        bitsOf(Vec3f, Vec3f.init(0.0, -0.0, 0.0).normalize()),
    );
    // The sign of a signed zero survives, which `a / 0` would not.
    try std.testing.expectEqual(
        [2]u64{ 0x8000_0000_0000_0000, 0x8000_0000_0000_0000 },
        bitsOf(Vec2d, Vec2d.init(-0.0, -0.0).normalize()),
    );
    try std.testing.expectEqual(
        [4]u64{ 0x8000_0000_0000_0000, 0, 0x8000_0000_0000_0000, 0 },
        bitsOf(Vec4d, Vec4d.init(-0.0, 0.0, -0.0, 0.0).normalize()),
    );
    try std.testing.expectEqual(@as(f64, 0.0), Vec4d.zero().normalize().length());
}

test "normalize does not overflow where a naive form would" {
    try std.testing.expect(productIsInfinite(f64, 1e300, 1e300));

    const unit = Vec2d.init(1e300, 1e300).normalize();
    try expectUnit(Vec2d, Vec2d.init(1e300, 1e300));
    try std.testing.expectEqual(unit.x, unit.y);
    try std.testing.expectApproxEqRel(1.0 / root_two, unit.x, 8 * std.math.floatEps(f64));

    const wide = Vec4f.init(1e30, 1e30, 1e30, 1e30).normalize();
    try expectUnit(Vec4f, Vec4f.init(1e30, 1e30, 1e30, 1e30));
    try std.testing.expectApproxEqRel(@as(f32, 0.5), wide.x, 8 * std.math.floatEps(f32));
    try std.testing.expectApproxEqRel(@as(f32, 0.5), wide.w, 8 * std.math.floatEps(f32));
}

test "normalize of a subnormal norm scales by a value with no relative precision" {
    // hypot of one ulp against one ulp is one ulp, the nearest representable
    // value to sqrt(2) times one ulp, so dividing by it gives the unit
    // components. The result is a direction but not a unit vector, and its
    // length is sqrt(2). The norm had no significant bits left to divide by,
    // so this is documented rather than fixed.
    const n = Vec2d.init(subnormal_d, subnormal_d).normalize();
    try std.testing.expect(n.eq(Vec2d.init(1.0, 1.0)));
    try std.testing.expectApproxEqRel(root_two, n.length(), 8 * std.math.floatEps(f64));

    // A norm that is exactly a component does give a unit result.
    const m = Vec2d.init(subnormal_d, 0.0).normalize();
    try std.testing.expectEqual(Vec2d.init(1.0, 0.0), m);
    try std.testing.expectEqual(@as(f64, 1.0), m.length());
}

test "normalize of an infinite component is undefined and gives the canonical NaN" {
    // The length of a vector with an infinite component is infinite, and
    // `infinity / infinity` is the invalid operation that produces the
    // positive canonical quiet NaN. A caller who needs a defined result for a
    // non-finite input has to check `isFinite` first.
    const v = Vec2d.init(std.math.inf(f64), 0.0);
    try std.testing.expect(std.math.isInf(v.length()));
    try std.testing.expect(!v.isFinite());
    try std.testing.expectEqual([2]u64{ qnan_d, 0 }, bitsOf(Vec2d, v.normalize()));
}

test "distance is the length of the difference" {
    // 3, 4, 5 by hand, in every arity that has one.
    try std.testing.expectEqual(@as(f64, 5.0), Vec2d.zero().distance(Vec2d.init(3, 4)));
    try std.testing.expectEqual(@as(f32, 5.0), Vec2f.zero().distance(Vec2f.init(3, 4)));
    try std.testing.expectEqual(@as(f64, 25.0), Vec2d.zero().distanceSquared(Vec2d.init(3, 4)));
    try std.testing.expectEqual(@as(f64, 5.0), Vec3d.init(0, 3, 4).distance(Vec3d.zero()));
    try std.testing.expectEqual(@as(f64, 5.0), Vec4d.init(0, 0, 3, 4).distance(Vec4d.zero()));
    try std.testing.expectEqual(@as(f32, 5.0), Vec3f.init(3, 0, 4).distance(Vec3f.zero()));

    const pairs = [_][2]Vec3d{
        .{ Vec3d.init(1, 2, 3), Vec3d.init(4, 5, 6) },
        .{ Vec3d.init(-1e200, 1e-200, 3), Vec3d.init(4e200, 5, -6) },
        .{ Vec3d.init(subnormal_d, 0, 0), Vec3d.init(0, subnormal_d, 0) },
    };
    for (pairs) |pair| {
        const a = pair[0];
        const b = pair[1];
        try std.testing.expectEqual(a.sub(b).length(), a.distance(b));
        try std.testing.expectEqual(a.sub(b).lengthSquared(), a.distanceSquared(b));
        // The difference is the exact negation either way round, so the
        // distance is symmetric bit for bit.
        try std.testing.expectEqual(bits(f64, a.distance(b)), bits(f64, b.distance(a)));
        try std.testing.expectEqual(
            bits(f64, a.distanceSquared(b)),
            bits(f64, b.distanceSquared(a)),
        );
        // A point has no distance from itself, and that distance is a positive
        // zero rather than a signed one.
        try std.testing.expectEqual(@as(u64, 0), bits(f64, a.distance(a)));
        try std.testing.expectEqual(@as(u64, 0), bits(f64, a.distanceSquared(a)));
    }
}

// ---------------------------------------------------------------------------
// lerp, min, max, clamp
// ---------------------------------------------------------------------------

test "lerp returns the exact endpoints" {
    // `a + (b - a) * 1` is not `b`: the difference 1 - 1e20 rounds to -1e20
    // before it is added back, so the general form would answer (0, 1e20).
    const a = Vec2d.init(1e20, 1.0);
    const b = Vec2d.init(1.0, 1e20);
    try std.testing.expectEqual(a, a.lerp(b, 0.0));
    try std.testing.expectEqual(b, a.lerp(b, 1.0));

    // The sign of a signed zero endpoint survives.
    try std.testing.expectEqual(
        [2]u64{ 0x8000_0000_0000_0000, 0 },
        bitsOf(Vec2d, Vec2d.init(-0.0, 0.0).lerp(b, 0.0)),
    );
    try std.testing.expectEqual(
        [2]u64{ 0, 0x8000_0000_0000_0000 },
        bitsOf(Vec2d, a.lerp(Vec2d.init(0.0, -0.0), 1.0)),
    );
    try std.testing.expectEqual(
        Vec2f.init(-0.0, 1.0),
        Vec2f.init(-0.0, 1.0).lerp(Vec2f.init(5.0, 5.0), 0.0),
    );
    try std.testing.expectEqual(
        Vec3f.init(1e20, 1.0, 2.0),
        Vec3f.init(1e20, 1.0, 2.0).lerp(Vec3f.init(1.0, 1e20, 3.0), 0.0),
    );
    try std.testing.expectEqual(
        Vec3d.init(1.0, 1e20, 3.0),
        Vec3d.init(1e20, 1.0, 2.0).lerp(Vec3d.init(1.0, 1e20, 3.0), 1.0),
    );
    try std.testing.expectEqual(
        Vec4d.init(1.0, 2.0, 3.0, 7.0),
        Vec4d.init(1.0, 2.0, 3.0, 7.0).lerp(Vec4d.init(5.0, 6.0, 7.0, 7.0), 0.0),
    );
    try std.testing.expectEqual(
        Vec4f.init(9.0, 8.0, 7.0, 6.0),
        Vec4f.init(1.0, 2.0, 3.0, 4.0).lerp(Vec4f.init(9.0, 8.0, 7.0, 6.0), 1.0),
    );
}

test "lerp reaches the points in between" {
    try std.testing.expectEqual(Vec3d.init(1, 2, 3), Vec3d.zero().lerp(Vec3d.init(2, 4, 6), 0.5));
    try std.testing.expectEqual(Vec2d.init(-1, 1), Vec2d.init(-3, -1).lerp(Vec2d.init(1, 3), 0.5));
    // From (0, 0) to (4, 8), a quarter of the way is (1, 2).
    try std.testing.expectEqual(Vec2f.init(1, 2), Vec2f.zero().lerp(Vec2f.init(4, 8), 0.25));
    try std.testing.expectEqual(
        Vec4d.init(0.5, 1.5, 2.5, 3.5),
        Vec4d.zero().lerp(Vec4d.init(1, 3, 5, 7), 0.5),
    );
    try std.testing.expectEqual(@as(f32, 1.0), Vec2f.zero().lerp(Vec2f.init(4, 0), 0.25).x);
    try std.testing.expectEqual(@as(f32, 3.0), Vec2f.zero().lerp(Vec2f.init(4, 0), 0.75).x);
    // The midpoint of a 3, 4, 5 triangle is 2.5 from each end.
    const midpoint = Vec2d.zero().lerp(Vec2d.init(3, 4), 0.5);
    try std.testing.expect(@abs(midpoint.length() - 2.5) <= 4 * std.math.floatEps(f64));
}

test "min and max are componentwise" {
    try std.testing.expectEqual(Vec2d.init(2, -1), Vec2d.init(3, -1).min(Vec2d.init(2, 2)));
    try std.testing.expectEqual(Vec2d.init(3, 2), Vec2d.init(3, -1).max(Vec2d.init(2, 2)));
    try std.testing.expectEqual(
        Vec3f.init(-1, 0, 2),
        Vec3f.init(0, 5, 7).min(Vec3f.init(-1, 0, 2)),
    );
    try std.testing.expectEqual(
        Vec3f.init(0, 5, 7),
        Vec3f.init(0, 5, 7).max(Vec3f.init(-1, 0, 2)),
    );
    try std.testing.expectEqual(Vec4d.init(0, 0, 0, 0), Vec4d.one().min(Vec4d.zero()));
    try std.testing.expectEqual(Vec4f.init(1, 1, 1, 1), Vec4f.one().max(Vec4f.zero()));
    // A larger magnitude wins whatever the sign. 1.0 is
    // 0x3ff0_0000_0000_0000 and -1.0 is 0xbff0_0000_0000_0000.
    try std.testing.expectEqual(
        @as(u64, 0xbff0_0000_0000_0000),
        bitsOf(Vec2d, Vec2d.init(0.0, 1.0).min(Vec2d.init(-0.0, -1.0)))[1],
    );
    try std.testing.expectEqual(
        @as(u64, 0x3ff0_0000_0000_0000),
        bitsOf(Vec2d, Vec2d.init(0.0, 1.0).max(Vec2d.init(-0.0, -1.0)))[1],
    );
}

test "min prefers the negative zero and max the positive one" {
    // fmin picks -0 and fmax picks +0, in either order of the arguments.
    // 1.0 is 0x3ff0_0000_0000_0000 in binary64 and 0x3f80_0000 in binary32.
    try std.testing.expectEqual(
        [2]u64{ 0x8000_0000_0000_0000, 0x3ff0_0000_0000_0000 },
        bitsOf(Vec2d, Vec2d.init(0.0, 1.0).min(Vec2d.init(-0.0, 1.0))),
    );
    try std.testing.expectEqual(
        [2]u64{ 0x8000_0000_0000_0000, 0x3ff0_0000_0000_0000 },
        bitsOf(Vec2d, Vec2d.init(-0.0, 1.0).min(Vec2d.init(0.0, 1.0))),
    );
    try std.testing.expectEqual(
        [2]u32{ 0, 0x3f80_0000 },
        bitsOf(Vec2f, Vec2f.init(0.0, 1.0).max(Vec2f.init(-0.0, 1.0))),
    );
    try std.testing.expectEqual(
        [2]u32{ 0, 0x3f80_0000 },
        bitsOf(Vec2f, Vec2f.init(-0.0, 1.0).max(Vec2f.init(0.0, 1.0))),
    );
    // fmin keeps each component's own zero when the two are the same value, so
    // a negative zero survives here, while fmax makes them all positive.
    try std.testing.expectEqual(
        [4]u32{ 0x8000_0000, 0, 0, 0x8000_0000 },
        bitsOf(Vec4f, Vec4f.init(-0.0, 0.0, 0.0, -0.0).min(Vec4f.zero())),
    );
    try std.testing.expectEqual(
        [4]u64{ 0, 0, 0, 0 },
        bitsOf(Vec4d, Vec4d.init(-0.0, 0.0, 0.0, -0.0).max(Vec4d.zero())),
    );
}

test "min and max are selection, so a NaN component loses" {
    const nan = std.math.nan(f64);
    // One NaN operand gives the other operand.
    try std.testing.expectEqual(
        Vec2d.init(2.0, 1.0),
        Vec2d.init(nan, 1.0).min(Vec2d.init(2.0, 3.0)),
    );
    try std.testing.expectEqual(
        Vec2d.init(2.0, 3.0),
        Vec2d.init(1.0, 3.0).max(Vec2d.init(2.0, nan)),
    );
    // Two NaN operands give a NaN, quieted, with the first sign.
    try std.testing.expectEqual(
        [2]u64{ qnan_d, 0x3ff0_0000_0000_0000 },
        bitsOf(Vec2d, Vec2d.init(nan, 1.0).min(Vec2d.init(nan, 1.0))),
    );
    try std.testing.expectEqual(
        [3]u32{ 0, qnan_f, 0 },
        bitsOf(Vec3f, Vec3f.init(0, nan, 0).max(Vec3f.init(0, nan, 0))),
    );
}

test "clamp is min after max, so the bounds decide the result" {
    try std.testing.expectEqual(
        Vec3d.init(1, 0, 0.5),
        Vec3d.init(5, -5, 0.5).clamp(Vec3d.zero(), Vec3d.one()),
    );
    // Bounds in the wrong order give the upper bound, because clamp is
    // min(max(v, lower), upper).
    try std.testing.expectEqual(
        Vec2d.init(1, 1),
        Vec2d.init(5, -5).clamp(Vec2d.init(2, 2), Vec2d.init(1, 1)),
    );
    // A value already inside the bounds is untouched.
    try std.testing.expectEqual(
        Vec4f.init(2, 3, 4, 5),
        Vec4f.init(2, 3, 4, 5).clamp(Vec4f.zero(), Vec4f.one().scale(6.0)),
    );
    try std.testing.expectEqual(
        Vec4d.init(0, 0, 0, 0),
        Vec4d.init(-1, -2, -3, -4).clamp(Vec4d.zero(), Vec4d.one()),
    );
    try std.testing.expectEqual(
        Vec4d.init(6, 6, 6, 6),
        Vec4d.init(7, 8, 9, 10).clamp(Vec4d.zero(), Vec4d.one().scale(6.0)),
    );
    try std.testing.expectEqual(
        Vec2f.init(0, 1),
        Vec2f.init(-1, 1).clamp(Vec2f.zero(), Vec2f.one()),
    );
}

test "clamp lets the bounds decide the signed zero" {
    try std.testing.expectEqual(
        [2]u64{ 0, 0x3ff0_0000_0000_0000 },
        bitsOf(Vec2d, Vec2d.init(0.0, 1.0).clamp(Vec2d.init(-0.0, 0.0), Vec2d.init(1.0, 1.0))),
    );
    // With a negative zero lower bound, a negative zero in the value survives
    // both the max and the min step.
    try std.testing.expectEqual(
        [2]u64{ 0x8000_0000_0000_0000, 0x3ff0_0000_0000_0000 },
        bitsOf(Vec2d, Vec2d.init(-0.0, 1.0).clamp(Vec2d.init(-0.0, 0.0), Vec2d.init(1.0, 1.0))),
    );
    // With a positive zero lower bound, fmax prefers the bound, so a negative
    // zero in the value becomes positive.
    try std.testing.expectEqual(
        [3]u32{ 0, 0, 0 },
        bitsOf(Vec3f, Vec3f.init(0.0, -0.0, -0.0).clamp(Vec3f.zero(), Vec3f.one())),
    );
    // A NaN component becomes the lower bound, and then the upper bound.
    try std.testing.expectEqual(
        [2]u64{ 0x3ff0_0000_0000_0000, 0 },
        bitsOf(Vec2d, Vec2d.init(std.math.nan(f64), 0.0).clamp(
            Vec2d.init(1.0, 0.0),
            Vec2d.init(3.0, 1.0),
        )),
    );
}

// ---------------------------------------------------------------------------
// eq and isFinite
// ---------------------------------------------------------------------------

test "eq is bit equality, so it is reflexive" {
    try std.testing.expect(Vec2d.init(1, -2).eq(Vec2d.init(1, -2)));
    try std.testing.expect(Vec2f.init(1, -2).eq(Vec2f.init(1, -2)));
    try std.testing.expect(Vec3d.init(1, -2, 3).eq(Vec3d.init(1, -2, 3)));
    try std.testing.expect(Vec3f.init(1, -2, 3).eq(Vec3f.init(1, -2, 3)));
    try std.testing.expect(Vec4d.init(1, -2, 3, -4).eq(Vec4d.init(1, -2, 3, -4)));
    try std.testing.expect(Vec4f.init(1, -2, 3, -4).eq(Vec4f.init(1, -2, 3, -4)));
    // Including on the values a float comparison cannot handle.
    try std.testing.expect(Vec2d.init(-0.0, subnormal_d).eq(Vec2d.init(-0.0, subnormal_d)));
    try std.testing.expect(
        Vec4d.init(0.0, -0.0, 0, -0.0).eq(Vec4d.init(0.0, -0.0, 0, -0.0)),
    );
    try std.testing.expect(
        Vec2d.init(std.math.nan(f64), 1.0).eq(Vec2d.init(std.math.nan(f64), 1.0)),
    );
}

test "eq separates the two zeros and one ulp" {
    // One ulp apart is a different component.
    try std.testing.expect(!Vec3d.init(1, 2, 3).eq(Vec3d.init(1, 2, nextUp(f64, 3.0))));
    try std.testing.expect(!Vec2f.init(1, 2).eq(Vec2f.init(1, nextUp(f32, 2.0))));
    try std.testing.expect(!Vec4d.init(1, 2, 3, 4).eq(Vec4d.init(1, 2, 3, nextUp(f64, 4.0))));
    // The two zeros are different components, which `==` cannot see.
    try std.testing.expect(!Vec2d.init(0.0, 1.0).eq(Vec2d.init(-0.0, 1.0)));
    try std.testing.expect(!Vec4f.init(1, 2, 3, -0.0).eq(Vec4f.init(1, 2, 3, 0.0)));
    try std.testing.expect(!Vec3f.init(-0.0, 1, 2).eq(Vec3f.init(0.0, 1, 2)));
    // A NaN is not equal to a number.
    try std.testing.expect(!Vec2d.init(std.math.nan(f64), 1.0).eq(Vec2d.init(0.0, 1.0)));
}

test "eq compares a NaN with a NaN, whatever the sign or the payload" {
    // 0x7ff8_0000_0000_0000 is the positive canonical quiet NaN,
    // 0xfff8_0000_0000_0001 a negative quiet NaN with a payload, and
    // 0x7ff0_0000_0000_0001 a signaling NaN, from
    // `spec/core-math-2026-08-22.md:769-801`.
    const canonical: f64 = @bitCast(@as(u64, qnan_d));
    const negative_payload: f64 = @bitCast(@as(u64, 0xfff8_0000_0000_0001));
    const signaling: f64 = @bitCast(@as(u64, 0x7ff0_0000_0000_0001));

    try std.testing.expect(
        Vec2d.init(canonical, 1.0).eq(Vec2d.init(negative_payload, 1.0)),
    );
    try std.testing.expect(Vec2d.init(canonical, 1.0).eq(Vec2d.init(signaling, 1.0)));
    try std.testing.expect(
        Vec4d.init(signaling, 0, 0, 0).eq(Vec4d.init(negative_payload, 0, 0, 0)),
    );

    const quiet32: f32 = @bitCast(@as(u32, qnan_f));
    const signaling32: f32 = @bitCast(@as(u32, 0x7f80_0001));
    try std.testing.expect(Vec3f.init(quiet32, 1, 2).eq(Vec3f.init(signaling32, 1, 2)));
}

test "eq is an equivalence relation across NaN payloads" {
    const a: f64 = @bitCast(@as(u64, 0x7ff8_0000_0000_0000));
    const b: f64 = @bitCast(@as(u64, 0xfff8_0000_0000_0001));
    const c: f64 = @bitCast(@as(u64, 0x7ff0_0000_0000_0001));
    const first = Vec4d.init(a, b, c, -0.0);
    const second = Vec4d.init(c, a, b, -0.0);
    const third = Vec4d.init(b, c, a, -0.0);
    try std.testing.expect(first.eq(second) and second.eq(first));
    try std.testing.expect(first.eq(third) and second.eq(third));
}

test "isFinite needs every component to be finite" {
    const inf_f: f32 = std.math.inf(f32);
    const inf_d: f64 = std.math.inf(f64);
    const nan_f: f32 = std.math.nan(f32);
    const nan_d: f64 = std.math.nan(f64);

    try std.testing.expect(Vec2d.init(1, 2).isFinite());
    try std.testing.expect(Vec2f.init(1, 2).isFinite());
    try std.testing.expect(Vec3d.init(subnormal_d, normal_d, 0).isFinite());
    try std.testing.expect(Vec4f.zero().isFinite());

    try std.testing.expect(!Vec2d.init(inf_d, 0).isFinite());
    try std.testing.expect(!Vec2d.init(0, -inf_d).isFinite());
    try std.testing.expect(!Vec2f.init(inf_f, 0).isFinite());
    try std.testing.expect(!Vec2f.init(0, nan_f).isFinite());
    try std.testing.expect(!Vec3d.init(1, 1, nan_d).isFinite());
    try std.testing.expect(!Vec4d.init(1, 1, 1, -inf_d).isFinite());
    try std.testing.expect(!Vec4f.init(1, 1, 1, nan_f).isFinite());
}

// ---------------------------------------------------------------------------
// Bridges to the @Vector layer, and the shape of the surface
// ---------------------------------------------------------------------------

test "toVector and fromVector move the components in lane order" {
    try std.testing.expectEqual([2]f32{ 1, 2 }, lanesOf(2, f32, Vec2f.init(1, 2).toVector()));
    try std.testing.expectEqual([2]f64{ 5, 6 }, lanesOf(2, f64, Vec2d.init(5, 6).toVector()));
    try std.testing.expectEqual([3]f32{ 7, 8, 9 }, lanesOf(3, f32, Vec3f.init(7, 8, 9).toVector()));
    try std.testing.expectEqual([3]f64{ 1, 2, 3 }, lanesOf(3, f64, Vec3d.init(1, 2, 3).toVector()));
    try std.testing.expectEqual(
        [4]f32{ 1, 2, 3, 4 },
        lanesOf(4, f32, Vec4f.init(1, 2, 3, 4).toVector()),
    );
    try std.testing.expectEqual(
        [4]f64{ 1, 2, 3, 4 },
        lanesOf(4, f64, Vec4d.init(1, 2, 3, 4).toVector()),
    );

    try std.testing.expectEqual(Vec2f.init(1, 2), Vec2f.fromVector(@Vector(2, f32){ 1, 2 }));
    try std.testing.expectEqual(Vec3d.init(1, 2, 3), Vec3d.fromVector(@Vector(3, f64){ 1, 2, 3 }));
    try std.testing.expectEqual(
        Vec4f.init(1, 2, 3, 4),
        Vec4f.fromVector(@Vector(4, f32){ 1, 2, 3, 4 }),
    );

    // The lanes are the components, in declaration order.
    const v = Vec4d.init(0.1, 0.2, 0.3, 0.4);
    const lanes = v.toVector();
    try std.testing.expectEqual(v.x, lanes[0]);
    try std.testing.expectEqual(v.y, lanes[1]);
    try std.testing.expectEqual(v.z, lanes[2]);
    try std.testing.expectEqual(v.w, lanes[3]);
}

test "toVector and fromVector round-trip the bits a float comparison cannot see" {
    const quiet32: f32 = @bitCast(@as(u32, qnan_f));
    try expectVectorRoundTrip(Vec2f, Vec2f.init(-0.0, subnormal_f));
    try expectVectorRoundTrip(Vec3f, Vec3f.init(quiet32, -0.0, subnormal_f));
    try expectVectorRoundTrip(
        Vec4f,
        Vec4f.init(subnormal_f, 0.1, -0.0, std.math.inf(f32)),
    );
    try expectVectorRoundTrip(Vec2d, Vec2d.init(-0.0, subnormal_d));
    try expectVectorRoundTrip(
        Vec3d,
        Vec3d.init(std.math.nan(f64), -0.0, subnormal_d),
    );
    try expectVectorRoundTrip(
        Vec4d,
        Vec4d.init(subnormal_d, -subnormal_d, 0.1, -0.0),
    );
}

test "the public surface is concrete and free of v-suffixed names" {
    const types = [_]type{ Vec2f, Vec3f, Vec4f, Vec2d, Vec3d, Vec4d };
    inline for (types) |V| {
        inline for (@typeInfo(V).@"struct".decls) |decl| {
            const name = decl.name;
            // Every `@Vector` entry point in `vector.zig` ends in `v`, so a
            // name that does is the only way the two layers could be confused.
            if (name.len != 0 and name[name.len - 1] == 'v') {
                std.debug.print("zmath Vec* name ends in v: {s}\n", .{name});
                return error.NameEndsInV;
            }
            switch (@typeInfo(@TypeOf(@field(V, name)))) {
                .@"fn" => |info| {
                    if (mentionsAnytype(info)) {
                        std.debug.print("zmath Vec* method takes anytype: {s}\n", .{name});
                        return error.AnytypeInSurface;
                    }
                },
                else => {},
            }
        }
    }
}
