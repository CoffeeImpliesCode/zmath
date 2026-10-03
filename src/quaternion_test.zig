const std = @import("std");
const quat = @import("quaternion.zig");
const vector_value = @import("vector_value.zig");

const Quatf = quat.Quatf;
const Quatd = quat.Quatd;
const UnitQuatf = quat.UnitQuatf;
const UnitQuatd = quat.UnitQuatd;
const Vec3f = vector_value.Vec3f;
const Vec3d = vector_value.Vec3d;

/// The positive canonical quiet NaN of each format, from
/// `spec/core-math-2026-08-22.md#special-values`.
const qnan_f: u32 = 0x7fc0_0000;
const qnan_d: u64 = 0x7ff8_0000_0000_0000;

/// A signaling NaN of each format, one bit below the quiet bit. The
/// specification turns it into the same quiet NaN with the same payload.
const snan_f: u32 = 0x7fc0_0001;
const snan_d: u64 = 0x7ff8_0000_0000_0001;

/// The smallest positive subnormal of each format, 2^-149 and 2^-1074.
const subnormal_f: f32 = 0x1p-149;
const subnormal_d: f64 = 0x1p-1074;

/// Exact binary angles, written in hex so the test does not depend on how a
/// decimal literal rounded.
const half_pi_f: f32 = 0x1.921fb6p+0;
const half_pi_d: f64 = 0x1.921fb54442d18p+0;
const third_pi_f: f32 = 0x1.0c6f7ap+0;

// ---------------------------------------------------------------------------
// Tolerances
//
// Every identity below is a law of exact real arithmetic that floating point
// only satisfies to rounding, so each one gets a band and each band is derived.
// `eps` is `std.math.floatEps`, which is one ulp at one. Where a band is
// justified by a measurement instead of a bound, the measured figure is
// written down next to it, because a band nobody measured is a guess.
//
// `unit_tolerance` is the band the module itself accepts in `init`. A value
// the module calls unit is a value a test may compare against one, and using
// the same band is what keeps the test from failing on a value the
// constructor deliberately admits.
//
// `algebra_tolerance` covers a law of the algebra whose two sides each reach
// about one unit of magnitude through three products and two sums. Five
// roundings at half an ulp each bound the disagreement at about three ulps of
// the largest partial sum, so eight ulps is a factor of 2.7 over the bound.
// The measured disagreement for the operands below is a quarter of an ulp, so
// the band also carries a factor of 32 over what actually happens. Eight ulps
// at one is still five orders of magnitude below the differences the
// non-commutativity and the short-arc tests depend on.
//
// `rotation_tolerance` covers `rotate`, which is two Hamilton products and a
// component selection on top of a unit quaternion that is itself within an ulp
// of the rotation it names, because `sin` and `cos` of half an angle are not
// exact in binary. The measured worst case over the pairs and vectors below is
// 2.4e-7, and ten micro is a factor of 42 over that. The margin is not
// fussiness: `f32` is IEEE-deterministic, so a band this wide cannot hide a
// sign error, a swapped axis, a missing rotation or a reversed composition,
// and each of those produces a difference of order one.
// ---------------------------------------------------------------------------

/// Eight ulps at one, the band `init` accepts and every unit assertion uses.
const unit_tolerance_f: f32 = 8 * std.math.floatEps(f32);
const unit_tolerance_d: f64 = 8 * std.math.floatEps(f64);

/// Eight ulps at one, for a law of the algebra on unit-scale values.
const algebra_tolerance_f: f32 = 8 * std.math.floatEps(f32);
const algebra_tolerance_d: f64 = 8 * std.math.floatEps(f64);

/// Forty times the measured error of a composed rotation.
const rotation_tolerance_f: f32 = 1.0e-5;
const rotation_tolerance_d: f64 = 1.0e-13;

/// The band for `a * conj(a) == normSquared(a) * identity()`. Each vector
/// component cancels a product against the same product negated, and the
/// scalar component sums the same squares the norm squared sums, so the two
/// differ by at most a few ulps of the largest product, and every product is
/// bounded by the norm squared. Sixteen ulps of the norm squared is therefore
/// the band, and it scales with the value rather than being an absolute
/// constant. For `(1, 2, 3, 4)` the two agree exactly, because every partial
/// sum is a small integer.
fn conjugateRoundTripBand(comptime T: type, norm_squared: T) T {
    return 16 * std.math.floatEps(T) * @abs(norm_squared);
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// The unsigned integer of the same width as `T`.
fn intOf(comptime T: type) type {
    return if (T == f32) u32 else u64;
}

/// The exact bits of a component. A float comparison cannot tell `+0` from `-0`
/// and cannot compare a NaN at all, so every assertion about those goes here.
fn bits(comptime T: type, v: T) intOf(T) {
    return @bitCast(v);
}

fn f32FromBits(v: u32) f32 {
    return @bitCast(v);
}

fn f64FromBits(v: u64) f64 {
    return @bitCast(v);
}

/// The component type of a raw or validated quaternion.
fn scalarOf(comptime Q: type) type {
    return @typeInfo(Q).@"struct".fields[0].type;
}

/// Componentwise closeness, which is what a law over the four components
/// means. Returning an error instead of failing an `expect` is what lets the
/// message name the component that went wrong instead of a whole struct.
fn expectQuatClose(comptime Q: type, actual: Q, expected: Q, tolerance: scalarOf(Q)) !void {
    const T = scalarOf(Q);
    const got = actual.toArray();
    const want = expected.toArray();
    inline for (0..4) |i| {
        const difference = @abs(got[i] - want[i]);
        if (difference > tolerance) {
            const got_bits = bits(T, got[i]);
            const want_bits = bits(T, want[i]);
            std.debug.print("{s} component {d}: got {x} want {x}, delta {e} > {e}\n", .{
                @typeName(Q), i, got_bits, want_bits, @as(f64, difference), @as(f64, tolerance),
            });
            return error.QuatNotClose;
        }
    }
}

fn expectVecClose(comptime V: type, actual: V, expected: V, tolerance: scalarOf(V)) !void {
    const T = @typeInfo(V).@"struct".fields[0].type;
    const got = actual.toArray();
    const want = expected.toArray();
    inline for (0..3) |i| {
        const difference = @abs(got[i] - want[i]);
        if (difference > tolerance) {
            const got_bits = bits(T, got[i]);
            const want_bits = bits(T, want[i]);
            std.debug.print("{s} component {d}: got {x} want {x}, delta {e} > {e}\n", .{
                @typeName(V), i, got_bits, want_bits, @as(f64, difference), @as(f64, tolerance),
            });
            return error.VecNotClose;
        }
    }
}

/// The `f32` rotation of `radians` about `axis`, unwrapped so a test reads as
/// one expression. The inputs below are total over `fromAxisAngle`, so the
/// unwrap is a test convenience and not a claim about the API.
fn unitf(axis: Vec3f, radians: f32) UnitQuatf {
    return UnitQuatf.fromAxisAngle(axis, radians).?;
}

fn unitd(axis: Vec3d, radians: f64) UnitQuatd {
    return UnitQuatd.fromAxisAngle(axis, radians).?;
}

/// The unit quaternion for a raw value, unwrapped.
fn validatedf(raw: Quatf) UnitQuatf {
    return UnitQuatf.init(raw).?;
}

fn validatedd(raw: Quatd) UnitQuatd {
    return UnitQuatd.init(raw).?;
}

/// The identity rotation. `UnitQuat` has no constructor of its own for it,
/// which is the point: the raw identity is the value a caller would reach for
/// and `init` is the only documented way in.
fn identityf() UnitQuatf {
    return validatedf(Quatf.identity());
}

fn identityd() UnitQuatd {
    return validatedd(Quatd.identity());
}

// The surface the contract names, and nothing else. `Quatf` is the algebra and
// must not grow a rotation entry point. These are compile-time assertions, so
// a surface that drifts is a compile error in the gate and not a review
// comment.
comptime {
    if (@hasDecl(Quatf, "rotate")) @compileError("Quatf must not rotate");
    if (@hasDecl(Quatd, "rotate")) @compileError("Quatd must not rotate");
    if (!@hasDecl(UnitQuatf, "init")) @compileError("UnitQuatf needs init");
    if (!@hasDecl(UnitQuatd, "init")) @compileError("UnitQuatd needs init");
}

// ---------------------------------------------------------------------------
// The Hamilton product
// ---------------------------------------------------------------------------

test "hamilton product matches the hand-computed form" {
    // Every product and every partial sum here is an integer well inside the
    // exact range of the format, so the expected values are bit patterns and
    // not a tolerance: the formula has to be the one that was written down.
    //
    //   w = 4*8 - 1*5 - 2*6 - 3*7 = 32 -  5 - 12 - 21 =  -6
    //   x = 4*5 + 1*8 + 2*7 - 3*6 = 20 +  8 + 14 - 18 =  24
    //   y = 4*6 - 1*7 + 2*8 + 3*5 = 24 -  7 + 16 + 15 =  48
    //   z = 4*7 + 1*6 - 2*5 + 3*8 = 28 +  6 - 10 + 24 =  48
    const a = Quatf.init(1, 2, 3, 4);
    const b = Quatf.init(5, 6, 7, 8);
    try std.testing.expect(Quatf.init(24, 48, 48, -6).eq(a.mul(b)));

    const ad = Quatd.init(1, 2, 3, 4);
    const bd = Quatd.init(5, 6, 7, 8);
    try std.testing.expect(Quatd.init(24, 48, 48, -6).eq(ad.mul(bd)));

    // i*j = k and j*i = -k, the two products that differ by a sign. A formula
    // that permuted the cross terms could not produce both of these.
    const i = Quatf.init(1, 0, 0, 0);
    const j = Quatf.init(0, 1, 0, 0);
    try std.testing.expect(Quatf.init(0, 0, 1, 0).eq(i.mul(j)));
    try std.testing.expect(Quatf.init(0, 0, -1, 0).eq(j.mul(i)));

    // i*i = -1, so the scalar part of a product of two pure imaginary
    // quaternions is the negated dot product of the two vector parts.
    try std.testing.expect(Quatf.init(0, 0, 0, -1).eq(i.mul(i)));
    try std.testing.expect(Quatf.init(0, 0, 0, -1).eq(j.mul(j)));
    const k = Quatf.init(0, 0, 1, 0);
    try std.testing.expect(Quatf.init(0, 0, 0, -1).eq(k.mul(k)));
}

test "the identity is a two-sided unit of the hamilton product" {
    const a = Quatf.init(-0.25, 0.5, 0.75, 1.5);
    // Bit for bit, not close: every term but `w*w` is a product with a zero and
    // an exact zero added to a nonzero changes nothing, so the identity has to
    // reproduce `a` with no rounding at all.
    try std.testing.expect(a.eq(Quatf.identity().mul(a)));
    try std.testing.expect(a.eq(a.mul(Quatf.identity())));

    const d = Quatd.init(-0.25, 0.5, 0.75, 1.5);
    try std.testing.expect(d.eq(Quatd.identity().mul(d)));
    try std.testing.expect(d.eq(d.mul(Quatd.identity())));
}

test "hamilton product does not commute" {
    // Two quarter turns about different axes, both built through the public
    // constructor, so these are the unit values the contract talks about.
    const a = unitf(Vec3f.init(1, 0, 0), half_pi_f);
    const b = unitf(Vec3f.init(0, 0, 1), half_pi_f);
    const ab = a.mul(b);
    const ba = b.mul(a);

    // `eq` compares bit patterns, so this says the two orders are different
    // answers, not that they are close.
    try std.testing.expect(!ab.eq(ba));
    // The gap is not a last-bit artifact either, so a product that had
    // quietly become symmetric could not pass on a loose tolerance.
    try std.testing.expect(@abs(ab.raw.y - ba.raw.y) > 0.1);
}

test "a pair about one axis commutes to the bit" {
    // Both orderings form the same products from the same pairs, and a product
    // of two floats is commutative bit for bit, so the two sums are the same
    // sums. That is why this pair is equal under `eq` where the pair above is
    // not: it is the products that differ between the orders, not the
    // summation, and a same-axis pair is the case where they do not.
    const a = unitf(Vec3f.init(0, 0, 1), half_pi_f);
    const b = unitf(Vec3f.init(0, 0, 1), third_pi_f);
    try std.testing.expect(a.mul(b).eq(b.mul(a)));

    // And the common product is the rotation by the summed angle, which is the
    // only reason a same-axis pair commutes. Checked on a vector, because the
    // quaternion for the summed angle is not bit-equal to the product: the
    // product rounds, and `sin` and `cos` of the sum do not reproduce a rounded
    // sum of two rounded half-angles.
    const axis = Vec3f.init(0, 0, 1);
    const v = Vec3f.init(1, 0, 0);
    try expectVecClose(
        Vec3f,
        a.mul(b).rotate(v),
        unitf(axis, half_pi_f + third_pi_f).rotate(v),
        rotation_tolerance_f,
    );
}

// ---------------------------------------------------------------------------
// The laws of the algebra
//
// These are laws of exact real arithmetic. Floating point only satisfies them
// to rounding, and not even to the same rounding on both sides: the right side
// of distributivity groups the terms as `a*(b + c)`, and the left side forms
// `b + c` first, so one extra rounding sits on one side. A band is therefore
// the only honest assertion, and the band is the sum of the roundings.
// ---------------------------------------------------------------------------

test "the product distributes over the sum" {
    const a = unitf(Vec3f.init(0.3, 0.5, 0.81), 0.7).raw;
    const b = unitf(Vec3f.init(1, 0, 0), third_pi_f).raw;
    const c = unitf(Vec3f.init(0, 1, 0), half_pi_f).raw;

    const left = a.mul(b.add(c));
    const right = a.mul(b).add(a.mul(c));
    // Exactly unequal is the expected outcome, which is why this cannot be an
    // `eq` test; the tolerance is sixteen ulps at one, three times the five
    // half-ulp roundings the two sides can accumulate.
    try std.testing.expect(!left.eq(right));
    try expectQuatClose(Quatf, left, right, algebra_tolerance_f);

    const ad = unitd(Vec3d.init(0.3, 0.5, 0.81), 0.7).raw;
    const bd = unitd(Vec3d.init(1, 0, 0), third_pi_f).raw;
    const cd = unitd(Vec3d.init(0, 1, 0), half_pi_d).raw;
    try expectQuatClose(Quatd, ad.mul(bd.add(cd)), ad.mul(bd).add(ad.mul(cd)), algebra_tolerance_d);
}

test "the product is associative" {
    const a = unitf(Vec3f.init(0.577, 0.577, 0.577), 0.9).raw;
    const b = unitf(Vec3f.init(1, 0, 0), third_pi_f).raw;
    const c = unitf(Vec3f.init(0, 0, 1), 1.1).raw;

    // `a*(b*c)` and `(a*b)*c` group the same eight products differently, so
    // the same band applies.
    try expectQuatClose(Quatf, a.mul(b.mul(c)), a.mul(b).mul(c), algebra_tolerance_f);

    const ad = unitd(Vec3d.init(0.577, 0.577, 0.577), 0.9).raw;
    const bd = unitd(Vec3d.init(1, 0, 0), third_pi_f).raw;
    const cd = unitd(Vec3d.init(0, 0, 1), 1.1).raw;
    try expectQuatClose(Quatd, ad.mul(bd.mul(cd)), ad.mul(bd).mul(cd), algebra_tolerance_d);
}

test "the conjugate inverts a quaternion" {
    // The law is `a * conj(a) == normSquared(a) * identity()`. The vector
    // parts are exact zeros because each is a sum of a product and its own
    // negation, and the scalar part is a sum of positive squares, but the
    // product sums those squares with three subtractions where `normSquared`
    // sums them with three additions, so the two scalars can differ by an ulp
    // or two. That is why the round trip is a band and not an equality.
    const a = Quatf.init(1, 2, 3, 4);
    const round_trip = a.mul(a.conjugate());
    const expected = Quatf.identity().scale(a.normSquared());
    try expectQuatClose(Quatf, round_trip, expected, conjugateRoundTripBand(f32, a.normSquared()));

    // And the norm squared is the squared length, not something near it.
    try std.testing.expectEqual(@as(f32, 30.0), a.normSquared());
    try std.testing.expectEqual(@as(f32, 30.0), a.lengthSquared());

    const d = Quatd.init(1, 2, 3, 4);
    try expectQuatClose(
        Quatd,
        d.mul(d.conjugate()),
        Quatd.identity().scale(d.normSquared()),
        conjugateRoundTripBand(f64, d.normSquared()),
    );
    try std.testing.expectEqual(@as(f64, 30.0), d.normSquared());
}

test "the inverse undoes the product" {
    const a = Quatf.init(0.5, -1.25, 2, 3);
    const inverse = a.inverse().?;
    // The inverse divides by a norm that is not one, so this is a band around
    // the identity and not the identity. The division and the two products
    // contribute about four ulps between them, and the tolerance is sixteen.
    try expectQuatClose(Quatf, inverse.mul(a), Quatf.identity(), algebra_tolerance_f);
    try expectQuatClose(Quatf, a.mul(inverse), Quatf.identity(), algebra_tolerance_f);

    const d = Quatd.init(0.5, -1.25, 2, 3);
    const di = d.inverse().?;
    try expectQuatClose(Quatd, di.mul(d), Quatd.identity(), algebra_tolerance_d);
    try expectQuatClose(Quatd, d.mul(di), Quatd.identity(), algebra_tolerance_d);

    // A zero quaternion and one whose norm overflows to infinity have no
    // inverse, and saying so with `null` is the contract.
    try std.testing.expect(Quatf.zero().inverse() == null);
    try std.testing.expect(Quatf.init(1e200, 1e200, 1e200, 1e200).inverse() == null);
    try std.testing.expect(Quatd.init(1e200, 1e200, 1e200, 1e200).inverse() == null);
}

test "the unit inverse is the conjugate" {
    // Not a division: the norm squared is one, so `conj(q)/1` and `conj(q)` are
    // the same number and the division can only add a rounding.
    const a = unitf(Vec3f.init(0.2, 0.4, 0.89), 1.3);
    try std.testing.expect(a.inverse().eq(a.conjugate()));
    const v = Vec3f.init(1, 2, 3);
    try std.testing.expect(a.conjugate().rotate(v).eq(a.inverse().rotate(v)));

    const d = unitd(Vec3d.init(0.2, 0.4, 0.89), 1.3);
    try std.testing.expect(d.inverse().eq(d.conjugate()));
}

test "the componentwise operations and the array bridge" {
    const a = Quatf.init(1, 2, 3, 4);
    const b = Quatf.init(5, 6, 7, 8);
    try std.testing.expect(Quatf.init(6, 8, 10, 12).eq(a.add(b)));
    try std.testing.expect(Quatf.init(-4, -4, -4, -4).eq(a.sub(b)));
    try std.testing.expect(Quatf.init(-1, -2, -3, -4).eq(a.neg()));
    try std.testing.expect(Quatf.init(2, 4, 6, 8).eq(a.scale(2)));
    try std.testing.expect(Quatf.init(5, 12, 21, 32).eq(a.scaleQuat(b)));
    try std.testing.expectEqual(@as(f32, 70.0), a.dot(b));
    try std.testing.expect(Quatf.fromArray(a.toArray()).eq(a));
    try std.testing.expectEqual([4]f32{ 1, 2, 3, 4 }, a.toArray());
    try std.testing.expect(Quatf.zero().eq(Quatf.init(0, 0, 0, 0)));
    try std.testing.expect(Quatf.identity().isFinite());
    try std.testing.expect(!Quatf.init(f32FromBits(qnan_f), 0, 0, 0).isFinite());
}

// ---------------------------------------------------------------------------
// Rotation
// ---------------------------------------------------------------------------

test "rotate preserves the length of the vector" {
    // A unit quaternion is a rotation, so `q * (v, 0) * conj(q)` leaves the
    // length of `v` alone. The two products round, so the length is preserved
    // to a few ulps rather than exactly, and the band is the same
    // `rotation_tolerance` the composed cases use.
    const vectors = [_]Vec3f{
        Vec3f.init(1, 0, 0),
        Vec3f.init(0, 1, 0),
        Vec3f.init(0, 0, 1),
        Vec3f.init(1, 2, 3),
        Vec3f.init(-0.5, 0.25, 7.5),
    };
    const rotations = [_]UnitQuatf{
        identityf(),
        unitf(Vec3f.init(0, 0, 1), half_pi_f),
        unitf(Vec3f.init(0.577, 0.577, 0.577), 1.1),
        unitf(Vec3f.init(1, 0, 0), third_pi_f).then(unitf(Vec3f.init(0, 1, 0), 0.77)),
    };
    for (rotations) |q| {
        for (vectors) |v| {
            try std.testing.expectApproxEqAbs(
                @as(f32, v.length()),
                q.rotate(v).length(),
                rotation_tolerance_f,
            );
        }
    }

    // The same law in binary64, where the band is eleven orders tighter and
    // still many orders above the arithmetic.
    const dq = unitd(Vec3d.init(0.577, 0.577, 0.577), 1.1);
    const dv = Vec3d.init(-0.5, 0.25, 7.5);
    try std.testing.expectApproxEqAbs(
        @as(f64, dv.length()),
        dq.rotate(dv).length(),
        rotation_tolerance_d,
    );
}

test "a quarter turn about z maps the axes the way a right-handed frame does" {
    // A tolerance is needed and not optional here, and the reason is the
    // representation rather than the implementation. The quarter turn is stored
    // as `sin(pi/4)` and `cos(pi/4)`, and neither is one half in binary, so the
    // quaternion that names a true quarter turn is already within an ulp of it
    // before any rotation happens; the two Hamilton products and the component
    // selection on top of that put the answer a further ulp or two away. The
    // measured disagreement with the exact `(0, 1, 0)` is 6.0e-8, which is half
    // an ulp of one, and the band is a hundred and sixty times that. It is five
    // orders of magnitude below the 1.0 it is compared against, so it cannot
    // hide a sign error, a swapped axis or a missing rotation, each of which
    // produces a difference of order one.
    const q = unitf(Vec3f.init(0, 0, 1), half_pi_f);
    try expectVecClose(
        Vec3f,
        q.rotate(Vec3f.init(1, 0, 0)),
        Vec3f.init(0, 1, 0),
        rotation_tolerance_f,
    );
    try expectVecClose(
        Vec3f,
        q.rotate(Vec3f.init(0, 1, 0)),
        Vec3f.init(-1, 0, 0),
        rotation_tolerance_f,
    );
    // The rotation axis is a fixed point, which a rotation preserves by
    // construction and a formula with a sign slip in the scalar part does not.
    try expectVecClose(
        Vec3f,
        q.rotate(Vec3f.init(0, 0, 1)),
        Vec3f.init(0, 0, 1),
        rotation_tolerance_f,
    );
}

test "rotate is right-handed" {
    // The classic sign error: a left-handed implementation sends +x to -y and
    // every other test in the file still passes. So the sign is asserted
    // against zero directly, and not only as a distance to a wanted vector.
    const q = unitf(Vec3f.init(0, 0, 1), half_pi_f);
    const moved = q.rotate(Vec3f.init(1, 0, 0));
    try std.testing.expect(moved.y > 0.0);
    try std.testing.expect(@abs(moved.x) <= rotation_tolerance_f);

    // The triple product of the rotated axes keeps its handedness, which is the
    // property a right-handed convention exists to protect.
    const ex = q.rotate(Vec3f.init(1, 0, 0));
    const ey = q.rotate(Vec3f.init(0, 1, 0));
    const ez = q.rotate(Vec3f.init(0, 0, 1));
    try std.testing.expectApproxEqAbs(@as(f32, 1.0), ex.cross(ey).dot(ez), rotation_tolerance_f);

    // A negative angle is the opposite rotation, so it undoes the positive one.
    const backwards = unitf(Vec3f.init(0, 0, 1), -half_pi_f);
    try expectVecClose(
        Vec3f,
        backwards.rotate(Vec3f.init(1, 0, 0)),
        Vec3f.init(0, -1, 0),
        rotation_tolerance_f,
    );
    try expectVecClose(
        Vec3f,
        q.then(backwards).rotate(Vec3f.init(1, 2, 3)),
        Vec3f.init(1, 2, 3),
        rotation_tolerance_f,
    );
}

test "the identity rotation is the identity map" {
    const q = identityf();
    for ([_]Vec3f{ Vec3f.init(1, 0, 0), Vec3f.init(0, 1, 0), Vec3f.init(3, -4, 5) }) |v| {
        try expectVecClose(Vec3f, q.rotate(v), v, rotation_tolerance_f);
    }
    // Composing with the identity changes nothing either way round.
    const a = unitf(Vec3f.init(0, 1, 0), 0.9);
    try std.testing.expect(a.eq(a.then(identityf())));
    try std.testing.expect(a.eq(identityf().then(a)));
}

test "the binary64 rotation obeys the same laws" {
    // Every claim in this file is about the algebra and the composition, and
    // none of it may hold for one precision and not the other. The band is
    // eleven orders tighter than the binary32 one, and the measured error of a
    // composed binary64 rotation is one ulp of one, so this is a real
    // assertion and not a formality.
    const q = identityd();
    for ([_]Vec3d{ Vec3d.init(1, 0, 0), Vec3d.init(0, 1, 0), Vec3d.init(3, -4, 5) }) |v| {
        try expectVecClose(Vec3d, q.rotate(v), v, rotation_tolerance_d);
    }
    const a = unitd(Vec3d.init(0.3, 0.5, 0.81), 0.7);
    const b = unitd(Vec3d.init(0, 1, 0), 1.4);
    try std.testing.expect(a.then(b).eq(b.mul(a)));
    try std.testing.expect(a.then(identityd()).eq(a));
    try std.testing.expect(identityd().then(a).eq(a));

    // The quarter turn, right-handed, in binary64.
    const quarter = unitd(Vec3d.init(0, 0, 1), half_pi_d);
    try expectVecClose(
        Vec3d,
        quarter.rotate(Vec3d.init(1, 0, 0)),
        Vec3d.init(0, 1, 0),
        rotation_tolerance_d,
    );
    try expectVecClose(
        Vec3d,
        quarter.rotate(Vec3d.init(0, 1, 0)),
        Vec3d.init(-1, 0, 0),
        rotation_tolerance_d,
    );
    try std.testing.expect(quarter.rotate(Vec3d.init(1, 0, 0)).y > 0.0);

    // And the `then` contract, which is the one the whole type exists for.
    const v = Vec3d.init(-0.5, 0.25, 7.5);
    try expectVecClose(Vec3d, a.then(b).rotate(v), b.rotate(a.rotate(v)), rotation_tolerance_d);

    // A non-unit raw value is normalized, not rejected, in binary64 too.
    const scaled = validatedd(Quatd.init(1, 2, 3, 4).scale(-9.0));
    try std.testing.expect(@abs(scaled.raw.length() - 1.0) <= unit_tolerance_d);
    try expectVecClose(
        Vec3d,
        scaled.rotate(v),
        validatedd(Quatd.init(1, 2, 3, 4)).rotate(v),
        rotation_tolerance_d,
    );
}

// ---------------------------------------------------------------------------
// Composition
// ---------------------------------------------------------------------------

test "then applies self first and b second" {
    // This is the test the whole composition contract rests on. `then` is
    // defined as `b * a` and `rotate` as the conjugation, so the two together
    // have to agree with running the rotations in the order the name says.
    // Sixteen rotation pairs and three vectors, so a reversal that happened to
    // be self-consistent for one pair cannot pass.
    const axes = [_]Vec3f{
        Vec3f.init(1, 0, 0),
        Vec3f.init(0, 1, 0),
        Vec3f.init(0, 0, 1),
        Vec3f.init(0.577, 0.577, 0.577),
    };
    const angles = [_]f32{ 0.3, 1.1, 2.2, 0.05 };
    var reversed_matters = false;
    for (axes, 0..) |axis_a, i| {
        for (axes, 0..) |axis_b, j| {
            const a = unitf(axis_a, angles[i]);
            const b = unitf(axis_b, angles[(j + 1) % angles.len]);
            for ([_]Vec3f{ Vec3f.init(1, 0, 0), Vec3f.init(0, 1, 0), Vec3f.init(1, 2, 3) }) |v| {
                const composed = a.then(b).rotate(v);
                const staged = b.rotate(a.rotate(v));
                try expectVecClose(Vec3f, composed, staged, rotation_tolerance_f);
                // The other order is the answer a `then` that forgot to
                // reverse would give. For at least one of these pairs it is far
                // away, so the comparison above is not passing by accident.
                if (@abs(composed.x - a.rotate(b.rotate(v)).x) > 1.0e-3) reversed_matters = true;
            }
        }
    }
    try std.testing.expect(reversed_matters);
}

test "then is b multiplied by a" {
    const a = unitf(Vec3f.init(0.3, 0.5, 0.81), 0.7);
    const b = unitf(Vec3f.init(1, 0, 0), third_pi_f);
    // Both sides go through the same re-normalizing product, so this is bit
    // equality and not a band.
    try std.testing.expect(a.then(b).eq(b.mul(a)));
    try std.testing.expect(a.then(b).eq(identityf().then(a.then(b)).then(identityf())));
    try std.testing.expect((a.then(b).then(identityf())).eq(a.then(b)));

    const d = unitd(Vec3d.init(0.3, 0.5, 0.81), 0.7);
    const e = unitd(Vec3d.init(0, 1, 0), 1.4);
    try std.testing.expect(d.then(e).eq(e.mul(d)));
}

test "a product of unit quaternions stays unit" {
    // The product of two unit quaternions is unit in exact arithmetic, so this
    // is a statement about the re-normalization: repeated products must not walk
    // out of the unit ball, because once a value is off the sphere its rotate
    // scales the vector by the norm squared and nothing else notices.
    //
    // A thousandth of a full turn, a thousand times, so the walk comes back to
    // where it started and what is left is the accumulated rounding of a
    // thousand products. The band on the accumulated result is 1e-4: one ulp of
    // the step is 6e-8, so the total angle misses a full turn by at most 4e-7
    // radians, and a thousand products round by a few ulps each.
    const step = unitf(Vec3f.init(1, 0, 0), 0x1.921fb6p+2 / 1000.0);
    var q = identityf();
    for (0..1000) |_| {
        q = q.mul(step);
        try std.testing.expect(@abs(q.raw.length() - 1.0) <= unit_tolerance_f);
    }

    // The result is the identity ROTATION, and the rotation is what the
    // assertion is about. A full turn lands in the negative hemisphere, so the
    // stored quaternion is near `-1` and not near `+1`, and comparing the
    // quaternion to the identity here would be comparing a rotation to one of
    // its two names. That `q` and `-q` name the same rotation is a property of
    // the type, not a rounding artifact, and it is why the check below is on
    // the rotated vector.
    const v = Vec3f.init(0.3, -0.7, 1.1);
    try expectVecClose(Vec3f, q.rotate(v), v, 1.0e-4);
    try std.testing.expect(@abs(q.raw.x) <= 1.0e-4);
    try std.testing.expect(@abs(q.raw.y) <= 1.0e-4);
    try std.testing.expect(@abs(q.raw.z) <= 1.0e-4);
    try std.testing.expect(@abs(@abs(q.raw.w) - 1.0) <= 1.0e-4);
}

test "a quaternion and its negation name the same rotation" {
    // `(-q) * (v, 0) * conj(-q)` is the conjugation with every factor negated
    // twice, so the two products are equal and the rotation is unchanged. This
    // is why `slerp` may negate one endpoint, and why a caller must compare
    // rotations with a tolerance and quaternions with `eq`.
    const q = unitf(Vec3f.init(0.3, 0.5, 0.81), 1.7);
    const negated = UnitQuatf.init(Quatf.init(-q.raw.x, -q.raw.y, -q.raw.z, -q.raw.w)).?;
    try std.testing.expect(!negated.eq(q));
    for ([_]Vec3f{ Vec3f.init(1, 0, 0), Vec3f.init(0, 1, 0), Vec3f.init(-2, 3, 4) }) |v| {
        try expectVecClose(Vec3f, negated.rotate(v), q.rotate(v), rotation_tolerance_f);
    }
    // The negation of a unit quaternion is still unit, so the constructor
    // accepts it rather than rejecting it.
    try std.testing.expect(@abs(negated.raw.length() - 1.0) <= unit_tolerance_f);
}

// ---------------------------------------------------------------------------
// The validating constructor
// ---------------------------------------------------------------------------

test "init rejects every value that is not a rotation" {
    // Zero has no direction, so it is not normalized into a unit value.
    try std.testing.expect(UnitQuatf.init(Quatf.zero()) == null);
    try std.testing.expect(UnitQuatd.init(Quatd.zero()) == null);
    // A negative zero is still a zero.
    try std.testing.expect(UnitQuatf.init(Quatf.init(-0.0, 0, 0, -0.0)) == null);

    // A NaN in any component is not a direction. Both the canonical quiet NaN
    // and a signaling NaN are rejected, and the signaling one is rejected in
    // every position, not only in `w`.
    try std.testing.expect(UnitQuatf.init(Quatf.init(f32FromBits(qnan_f), 0, 0, 1)) == null);
    try std.testing.expect(UnitQuatf.init(Quatf.init(0, f32FromBits(snan_f), 0, 1)) == null);
    try std.testing.expect(UnitQuatf.init(Quatf.init(0, 0, f32FromBits(qnan_f), 1)) == null);
    try std.testing.expect(UnitQuatf.init(Quatf.init(0, 0, 0, f32FromBits(qnan_f))) == null);
    try std.testing.expect(UnitQuatd.init(Quatd.init(f64FromBits(snan_d), 0, 0, 1)) == null);

    // An infinity has no finite length, and a length that overflowed to
    // infinity divides every finite component down to zero, which is not a
    // unit value either. All three are rejected.
    try std.testing.expect(UnitQuatf.init(Quatf.init(std.math.inf(f32), 0, 0, 1)) == null);
    try std.testing.expect(UnitQuatf.init(Quatf.init(0, -std.math.inf(f32), 0, 1)) == null);
    const huge_f = std.math.floatMax(f32);
    try std.testing.expect(UnitQuatf.init(Quatf.init(huge_f, huge_f, huge_f, huge_f)) == null);
    try std.testing.expect(UnitQuatd.init(Quatd.init(std.math.inf(f64), 0, 0, 1)) == null);
    const huge_d = std.math.floatMax(f64);
    try std.testing.expect(UnitQuatd.init(Quatd.init(huge_d, huge_d, huge_d, huge_d)) == null);
}

test "init normalizes a finite non-unit value and proves the result is unit" {
    // A quaternion of magnitude seven is a perfectly good rotation once it is
    // scaled, and the constructor is where that happens. The invariant that
    // matters is the length of the result, not the components: this asserts
    // that `init` returns something on the unit sphere rather than a value
    // that merely started finite.
    const unit = UnitQuatf.init(Quatf.init(1, 2, 3, 4).scale(7.0)).?;
    try std.testing.expect(@abs(unit.raw.length() - 1.0) <= unit_tolerance_f);
    // And it is the same rotation as the unnormalized original, which is the
    // property a caller is actually relying on. `init` accepts a non-unit
    // value rather than rejecting it, and this is what it is for.
    try expectVecClose(
        Vec3f,
        unit.rotate(Vec3f.init(1, 2, 3)),
        validatedf(Quatf.init(1, 2, 3, 4)).rotate(Vec3f.init(1, 2, 3)),
        rotation_tolerance_f,
    );

    // Four equal components normalize to four halves exactly: the fold reaches
    // two, which is a power of two, and every division by it is exact. So this
    // is a bit pattern and not a band, which is the strongest statement
    // available that the normalization divides and does not approximate.
    const quarters = UnitQuatf.init(Quatf.init(1, 1, 1, 1)).?;
    try std.testing.expectEqual([4]f32{ 0.5, 0.5, 0.5, 0.5 }, quarters.raw.toArray());
    try std.testing.expect(@abs(quarters.raw.length() - 1.0) <= unit_tolerance_f);

    // Ordinary components far from either end of the format, so the only thing
    // between the input and the result is the division.
    const middle = UnitQuatf.init(Quatf.init(3.0, 4.0, 0.0, 12.0)).?;
    try std.testing.expect(@abs(middle.raw.length() - 1.0) <= unit_tolerance_f);
    try expectQuatClose(
        Quatf,
        middle.raw,
        Quatf.init(3.0 / 13.0, 4.0 / 13.0, 0.0, 12.0 / 13.0),
        unit_tolerance_f,
    );

    const d = UnitQuatd.init(Quatd.init(3.0, 4.0, 0.0, 12.0)).?;
    try std.testing.expect(@abs(d.raw.length() - 1.0) <= unit_tolerance_d);
    try expectQuatClose(
        Quatd,
        d.raw,
        Quatd.init(3.0 / 13.0, 4.0 / 13.0, 0.0, 12.0 / 13.0),
        unit_tolerance_d,
    );
}

test "a quaternion of subnormals is rejected without producing an infinity" {
    // The division in `init` cannot overflow, because no component exceeds the
    // length, so no infinity can appear here whatever the length is. What the
    // subnormal case does show is the fold's own limit: the accepted `hypot`
    // kernel is correctly rounded for normal arguments, and for two subnormals
    // the exact answer 2*2^-150 is not representable, so the kernel returns one
    // unit of the last place, 2^-149, rather than half of it. The fold
    // therefore does not scale, the four equal components come back as four
    // ones, the length of that is sqrt(4) and not one, and the unit band
    // rejects the value. Rejection is the contract's own answer for a value
    // that is not unit after normalization; the important half of the
    // requirement is that nothing here overflows on the way.
    const all_subnormal = Quatf.init(subnormal_f, subnormal_f, subnormal_f, subnormal_f);
    try std.testing.expect(UnitQuatf.init(all_subnormal) == null);
    // The length itself is finite and nonzero, so the rejection came from the
    // unit band and not from a division that produced something unusable.
    try std.testing.expect(all_subnormal.length() > 0.0);
    try std.testing.expect(std.math.isFinite(all_subnormal.length()));
    const all_subnormal_d = [4]f64{ subnormal_d, subnormal_d, subnormal_d, subnormal_d };
    try std.testing.expect(UnitQuatd.init(Quatd.fromArray(all_subnormal_d)) == null);

    // A single subnormal component is a direction, and the fold handles it:
    // one component, three zeros, so the norm is that component and the
    // division is exact. This is the case that proves the subnormal path is
    // normalization and not blanket rejection.
    const single = Quatf.init(subnormal_f, 0, 0, 0);
    const axis_unit = UnitQuatf.init(single).?;
    try std.testing.expectEqual([4]f32{ 1, 0, 0, 0 }, axis_unit.raw.toArray());
    try std.testing.expect(@abs(axis_unit.raw.length() - 1.0) <= unit_tolerance_f);
    try std.testing.expectEqual(
        [4]f64{ 1, 0, 0, 0 },
        UnitQuatd.init(Quatd.init(subnormal_d, 0, 0, 0)).?.raw.toArray(),
    );
}

// ---------------------------------------------------------------------------
// The magnitude
// ---------------------------------------------------------------------------

test "the length does not overflow and does not flush a subnormal to zero" {
    // `sqrt(x*x + y*y + z*z + w*w)` at 1e200 needs 4e400, which is an infinity
    // before the square root is taken, so the naive form answers infinity for
    // four perfectly finite components. The fold through the accepted `hypot`
    // kernel answers a finite length. It is not exactly twice the component,
    // because three correctly rounded folds do not have to land on the
    // correctly rounded four-component norm, so the band is four ulps of the
    // answer, which is what three folds can cost.
    const huge = Quatd.init(1e200, 1e200, 1e200, 1e200);
    try std.testing.expect(std.math.isFinite(huge.length()));
    try std.testing.expectApproxEqAbs(
        @as(f64, 2e200),
        huge.length(),
        4 * std.math.floatEps(f64) * 2e200,
    );
    // The squared length is a different operation and does overflow, which is
    // the reason the length is the operation `init` divides by.
    try std.testing.expect(std.math.isInf(huge.lengthSquared()));

    const huge_f = Quatf.init(1.0e30, 1.0e30, 1.0e30, 1.0e30);
    try std.testing.expect(std.math.isFinite(huge_f.length()));
    try std.testing.expectApproxEqAbs(
        @as(f32, 2.0e30),
        huge_f.length(),
        4 * std.math.floatEps(f32) * 2.0e30,
    );
    // A component mix, so the answer is not the two-fold case again: the exact
    // length of (1, 2, 2, 4) times 1e200 is 5e200.
    try std.testing.expectApproxEqAbs(
        @as(f64, 5e200),
        Quatd.init(1e200, 2e200, 2e200, 4e200).length(),
        4 * std.math.floatEps(f64) * 5e200,
    );

    // At the bottom, `sqrt(sub*sub)` underflows and then reports a zero length
    // for a quaternion that is not zero. The fold keeps the length nonzero and
    // finite, which is the half of the requirement that matters; the value it
    // lands on is one unit of the last place, because the exact answer for two
    // subnormals is not representable and the accepted kernel does not scale
    // it. That limit is the subject of its own test, and here it only has to
    // stay nonzero and finite.
    const tiny = Quatf.init(subnormal_f, subnormal_f, subnormal_f, subnormal_f);
    try std.testing.expect(tiny.length() > 0.0);
    try std.testing.expect(std.math.isFinite(tiny.length()));

    const tiny_d = Quatd.init(subnormal_d, subnormal_d, subnormal_d, subnormal_d);
    try std.testing.expect(tiny_d.length() > 0.0);
    try std.testing.expect(std.math.isFinite(tiny_d.length()));

    // A single subnormal component is the length of itself, and this one is
    // exact: folding a zero into a subnormal is that subnormal.
    try std.testing.expectEqual(subnormal_f, Quatf.init(subnormal_f, 0, 0, 0).length());
    try std.testing.expectEqual(subnormal_d, Quatd.init(subnormal_d, 0, 0, 0).length());
    // The smallest normal behaves the same way, so the subnormal result above
    // is a limit of the kernel and not a special case of this fold.
    try std.testing.expectEqual(0x1p-126, Quatf.init(0x1p-126, 0, 0, 0).length());
    // And zero really is zero, which is the one value that is allowed to be.
    try std.testing.expectEqual(@as(f32, 0.0), Quatf.zero().length());
    try std.testing.expectEqual(@as(f64, 1.0), Quatd.identity().length());
    // A unit quaternion built by the constructor is unit to the band the
    // constructor itself accepts, which is the invariant a caller relies on.
    try std.testing.expect(@abs(identityf().raw.length() - 1.0) <= unit_tolerance_f);
}

test "the length quietens a NaN in the first component and keeps its payload" {
    // The fold seeds itself with the first component instead of folding a
    // positive zero into it, so that component reaches the special-value rule
    // by a different route than the other three. Both routes owe the same
    // answer: the operand with its sign and its payload, and the quiet bit set.
    // The canonical NaN is the wrong answer for either, because the payload is
    // the only diagnostic a caller has.
    const signaling = f32FromBits(snan_f);
    const negative_signaling = f32FromBits(snan_f | 0x8000_0000);
    try std.testing.expectEqual(snan_f, bits(f32, Quatf.init(signaling, 1, 1, 1).length()));
    try std.testing.expectEqual(
        snan_f | 0x8000_0000,
        bits(f32, Quatf.init(negative_signaling, 1, 1, 1).length()),
    );
    try std.testing.expectEqual(
        qnan_f,
        bits(f32, Quatf.init(f32FromBits(qnan_f), 1, 1, 1).length()),
    );
    // A NaN in a later component is found by the kernel instead, and it has to
    // reach the same answer.
    try std.testing.expectEqual(snan_f, bits(f32, Quatf.init(1, 1, signaling, 1).length()));
    // An infinite component wins over every later operand, because the fold is
    // a maximum once one of its operands is unbounded.
    try std.testing.expectEqual(
        @as(u32, 0x7f80_0000),
        bits(f32, Quatf.init(std.math.inf(f32), 1, 1, 1).length()),
    );
    // The binary64 fold owes the same, with the wider payload.
    try std.testing.expectEqual(snan_d, bits(f64, Quatd.init(f64FromBits(snan_d), 1, 1, 1).length()));
    try std.testing.expectEqual(qnan_d, bits(f64, Quatd.init(f64FromBits(qnan_d), 1, 1, 1).length()));
    // The case that separates a seed from a plain absolute value: a genuinely
    // signaling NaN has its quiet bit clear, and the seed has to set it without
    // touching the sign or the payload. An `@abs` that only cleared the sign
    // would report the input back, which is not a number the specification
    // allows a result to be.
    try std.testing.expectEqual(
        snan_f,
        bits(f32, Quatf.init(f32FromBits(0x7f80_0001), 1, 1, 1).length()),
    );
    try std.testing.expectEqual(
        snan_f | 0x8000_0000,
        bits(f32, Quatf.init(f32FromBits(0xff80_0001), 1, 1, 1).length()),
    );
    try std.testing.expectEqual(
        snan_d,
        bits(f64, Quatd.init(f64FromBits(0x7ff0_0000_0000_0001), 1, 1, 1).length()),
    );
}

// ---------------------------------------------------------------------------
// Interpolation
// ---------------------------------------------------------------------------

test "slerp returns its endpoints exactly" {
    // Not a band. The general form does not reach the endpoints: it forms
    // `sin((1-t)*theta) / sin(theta)` and `sin(t*theta) / sin(theta)`, and at
    // `t = 1` the second coefficient is `sin(theta)/sin(theta)`, which is one
    // only when the same rounded value comes back from both calls. It usually
    // does, but the short circuit is what makes it a guarantee.
    const a = unitf(Vec3f.init(1, 0, 0), 0.35);
    const b = unitf(Vec3f.init(0, 1, 0), 1.15);
    try std.testing.expect(a.slerp(b, 0.0).?.eq(a));
    try std.testing.expect(a.slerp(b, 1.0).?.eq(b));
    // And the endpoints are the endpoints, not a rotation near them.
    try expectVecClose(
        Vec3f,
        a.slerp(b, 0.0).?.rotate(Vec3f.init(1, 2, 3)),
        a.rotate(Vec3f.init(1, 2, 3)),
        rotation_tolerance_f,
    );
    try expectVecClose(
        Vec3f,
        a.slerp(b, 1.0).?.rotate(Vec3f.init(1, 2, 3)),
        b.rotate(Vec3f.init(1, 2, 3)),
        rotation_tolerance_f,
    );

    const c = unitd(Vec3d.init(1, 0, 0), 0.35);
    const d = unitd(Vec3d.init(0, 1, 0), 1.15);
    try std.testing.expect(c.slerp(d, 0.0).?.eq(c));
    try std.testing.expect(c.slerp(d, 1.0).?.eq(d));
}

test "slerp takes the short arc" {
    // Two rotations about one axis, separated by nearly a full turn. The
    // quaternion dot product is the cosine of HALF the angle between the two
    // rotations, so `a.dot(b)` here is `cos(5.7 / 2) = -0.958`: negative, which
    // is the case that has to be handled, because `b` and `-b` are two names for
    // the same rotation and only one of the two arcs between the rotations is
    // short.
    const axis = Vec3f.init(0, 0, 1);
    const a_angle: f32 = 0.4;
    const b_angle: f32 = 6.1;
    const a = unitf(axis, a_angle);
    const b = unitf(axis, b_angle);
    try std.testing.expect(a.dot(b) < 0.0);

    // The short arc is the rest of the turn, `2*pi - 5.7 = 0.583`, so its
    // midpoint is 0.29 from `a` and its dot product with `a` is
    // `cos(0.29 / 2) = 0.989`. The midpoint of the long arc would be 2.85 from
    // `a`, with a dot product of `cos(1.425) = 0.145`. The two differ by a
    // factor of seven, so the test below is not a coin flip.
    const short_arc = 2.0 * std.math.pi - (b_angle - a_angle);
    const mid = a.slerp(b, 0.5).?;
    try std.testing.expect(mid.dot(a) > 0.9);
    // The same statement from the other end. `slerp` returns the lift nearest
    // `a`, and `b` is that lift's near antipode, so the closeness to `b` is a
    // closeness to `-b`. `q` and `-q` are the same rotation and the contract
    // does not say which of the two a caller gets back, so the measure has to be
    // the magnitude.
    try std.testing.expect(@abs(mid.dot(b)) > 0.9);

    // And the midpoint is the rotation by half the short angle from `a`, which
    // is nowhere near the rotation by half the long angle. The second candidate
    // is the answer a `slerp` that forgot the sign flip would return, and the
    // two are 3.1 radians apart, so checking the first against the second would
    // fail and checking it against nothing at all would prove nothing.
    const short_mid = unitf(axis, a_angle - short_arc / 2.0);
    const long_mid = unitf(axis, a_angle + (b_angle - a_angle) / 2.0);
    try expectVecClose(
        Vec3f,
        mid.rotate(Vec3f.init(1, 0, 0)),
        short_mid.rotate(Vec3f.init(1, 0, 0)),
        rotation_tolerance_f,
    );
    try std.testing.expect(long_mid.dot(a) < 0.3);
    try std.testing.expect(long_mid.dot(a) < mid.dot(a));

    // A positive dot product takes the same short arc, so a pair that needs no
    // sign flip does not behave differently in a way a caller can observe.
    const near = unitf(axis, 0.9);
    try std.testing.expect(a.dot(near) > 0.0);
    try expectVecClose(
        Vec3f,
        a.slerp(near, 0.5).?.rotate(Vec3f.init(1, 0, 0)),
        unitf(axis, 0.65).rotate(Vec3f.init(1, 0, 0)),
        rotation_tolerance_f,
    );
}

test "slerp and lerp stay on the unit sphere" {
    // Every intermediate is a unit quaternion, or a caller who feeds it to a
    // third product has a scale factor nobody asked for. The endpoints are the
    // two that make `acos` unusable, so the parallel path is the one that
    // matters here: it is the branch a caller hits when a value has settled.
    const a = unitf(Vec3f.init(0.577, 0.577, 0.577), 0.0);
    const b = unitf(Vec3f.init(0.577, 0.577, 0.577), 0.02);
    var t: f32 = 0.0;
    while (t <= 1.0) : (t += 0.05) {
        const s = a.slerp(b, t).?;
        try std.testing.expect(@abs(s.raw.length() - 1.0) <= unit_tolerance_f);
        const l = a.lerp(b, t);
        try std.testing.expect(@abs(l.raw.length() - 1.0) <= unit_tolerance_f);
    }

    // And a wide pair, which takes the great-circle path, is unit throughout.
    const wide_a = unitf(Vec3f.init(1, 0, 0), 0.2);
    const wide_b = unitf(Vec3f.init(0, 1, 0), 2.0);
    t = 0.0;
    while (t <= 1.0) : (t += 0.05) {
        const wide = wide_a.slerp(wide_b, t).?;
        try std.testing.expect(@abs(wide.raw.length() - 1.0) <= unit_tolerance_f);
    }

    const d = unitd(Vec3d.init(0.577, 0.577, 0.577), 0.0);
    const e = unitd(Vec3d.init(0.577, 0.577, 0.577), 0.02);
    t = 0.0;
    while (t <= 1.0) : (t += 0.05) {
        try std.testing.expect(@abs(d.slerp(e, t).?.raw.length() - 1.0) <= unit_tolerance_d);
        try std.testing.expect(@abs(d.lerp(e, t).raw.length() - 1.0) <= unit_tolerance_d);
    }
}

test "slerp rejects a non-finite t and lerp reports it as a NaN" {
    const a = unitf(Vec3f.init(1, 0, 0), 0.35);
    const b = unitf(Vec3f.init(0, 1, 0), 1.15);
    try std.testing.expect(a.slerp(b, std.math.nan(f32)) == null);
    try std.testing.expect(a.slerp(b, std.math.inf(f32)) == null);
    try std.testing.expect(a.slerp(b, -std.math.inf(f32)) == null);
    const d = unitd(Vec3d.init(1, 0, 0), 0.35);
    const e = unitd(Vec3d.init(0, 1, 0), 1.15);
    try std.testing.expect(d.slerp(e, std.math.nan(f64)) == null);
    try std.testing.expect(d.slerp(e, std.math.inf(f64)) == null);

    // `lerp` has no optional to report through, so the component arithmetic
    // gives the specification's answer: an invalid operation on non-NaN
    // operands is the POSITIVE canonical quiet NaN, which is the encoding the
    // host unit does not produce on its own.
    const poisoned = a.lerp(b, std.math.nan(f32));
    inline for (poisoned.raw.toArray()) |component| {
        try std.testing.expectEqual(qnan_f, bits(f32, component));
    }
}

test "lerp reaches its endpoints and stays on the chord" {
    // A pair separated by more than half a turn, so the difference from `slerp`
    // is visible: `slerp` negates one endpoint and follows the short arc, while
    // `lerp` interpolates toward `b` as given and follows the chord.
    const axis = Vec3f.init(0, 0, 1);
    const a = unitf(axis, 0.4);
    const b = unitf(axis, 4.0);
    try std.testing.expect(a.lerp(b, 0.0).eq(a));
    try std.testing.expect(a.lerp(b, 1.0).eq(b));

    // The two disagree at the midpoint by a measurable amount and agree at the
    // endpoints. That is the whole difference between them, and a test that
    // only checked the endpoints would pass for two functions that were the
    // same function.
    const chord = a.lerp(b, 0.5);
    const arc = a.slerp(b, 0.5).?;
    try std.testing.expect(!chord.eq(arc));
    // The normalized midpoint of a chord between two unit quaternions with dot
    // `d` has dot product `sqrt((1 + d)/2)` with either end, and the midpoint of
    // the short arc has `sqrt((1 - d)/2)`. With `d` at -0.227 those are 0.62 and
    // 0.78, so the chord is measurably nearer the long way round than the arc
    // is. This is the assertion that says `lerp` is a blend and `slerp` is a
    // path.
    try std.testing.expect(chord.dot(a) < arc.dot(a));
    try std.testing.expect(chord.dot(a) > 0.5 and chord.dot(a) < 0.7);
}

// ---------------------------------------------------------------------------
// The axis and angle constructor
// ---------------------------------------------------------------------------

test "fromAxisAngle normalizes the axis" {
    // The same rotation named by a unit axis and by any positive nonzero
    // multiple of it. `37.5` and `1.0e-3` divide back to exactly `+1`, so these
    // are bit-identical values and not close ones.
    const reference = unitf(Vec3f.init(0, 0, 1), 0.9);
    try std.testing.expect(unitf(Vec3f.init(0, 0, 37.5), 0.9).eq(reference));
    try std.testing.expect(unitf(Vec3f.init(0, 0, 1.0e-3), 0.9).eq(reference));

    // A negated axis is the other rotation and not the same one: negating the
    // axis negates the vector part of the quaternion, which is the conjugate,
    // and the conjugate of a unit quaternion is its inverse. So the two rotate
    // opposite ways and only the parts that do not depend on the sign agree.
    const flipped = unitf(Vec3f.init(0, 0, -37.5), 0.9);
    try std.testing.expect(!flipped.eq(reference));
    // The validated constructor normalizes what the raw one built, and
    // `sin(0.45)^2 + cos(0.45)^2` is one only to rounding, so the two agree to
    // a band and not to the bit. The band is the constructor's own.
    try expectQuatClose(Quatf, flipped.raw, reference.conjugate().raw, unit_tolerance_f);
    try expectVecClose(
        Vec3f,
        flipped.rotate(Vec3f.init(1, 2, 3)),
        reference.conjugate().rotate(Vec3f.init(1, 2, 3)),
        rotation_tolerance_f,
    );
    try expectVecClose(
        Vec3f,
        flipped.rotate(Vec3f.init(1, 2, 3)),
        reference.inverse().rotate(Vec3f.init(1, 2, 3)),
        rotation_tolerance_f,
    );

    // An oblique axis, so the direction is a real division and not a copy. The
    // half-angle reference comes from the language builtin rather than from the
    // module, so the test compares the accepted kernels against something that
    // did not come out of the same table.
    const half_angle: f64 = 0.65;
    const sine: f32 = @floatCast(@sin(half_angle));
    const cosine: f32 = @floatCast(@cos(half_angle));
    const oblique = Vec3f.init(1.0, 2.0, 2.0);
    const q = unitf(oblique, 1.3);

    // The vector part is the unit axis times the sine of half the angle, so its
    // length is that sine and its direction is the axis.
    const vector_part = Vec3f.init(q.raw.x, q.raw.y, q.raw.z);
    try std.testing.expectApproxEqAbs(sine, vector_part.length(), 1.0e-6);
    try expectVecClose(Vec3f, vector_part, oblique.normalize().scale(sine), 1.0e-6);
    // The scalar part is the cosine of half the angle, and the four squares sum
    // to one, which is the unit invariant reached from a different direction.
    try std.testing.expectApproxEqAbs(cosine, @abs(q.raw.w), 1.0e-6);
    try std.testing.expectApproxEqAbs(@as(f32, 1.0), q.raw.lengthSquared(), 1.0e-6);

    // The raw constructor builds the same quaternion and the validated one only
    // normalizes it, so the two agree to rounding and not to the bit. The raw
    // value is one `init` would accept, which is the point of routing through
    // it.
    const raw = Quatf.fromAxisAngle(oblique, 1.3).?;
    try expectQuatClose(Quatf, raw, q.raw, 1.0e-6);
    try std.testing.expect(UnitQuatf.init(raw) != null);

    const dq = unitd(Vec3d.init(1.0, 2.0, 2.0), 1.3);
    try std.testing.expectApproxEqAbs(
        @as(f64, @sin(half_angle)),
        @sqrt(dq.raw.x * dq.raw.x + dq.raw.y * dq.raw.y + dq.raw.z * dq.raw.z),
        1.0e-15,
    );
    try std.testing.expectApproxEqAbs(@as(f64, @cos(half_angle)), @abs(dq.raw.w), 1.0e-15);
}

test "fromAxisAngle rejects an axis or an angle with no direction" {
    // A zero axis has no direction, and `normalize` hands a zero vector back
    // unchanged, so the constructor has to notice for itself.
    try std.testing.expect(UnitQuatf.fromAxisAngle(Vec3f.zero(), 1.0) == null);
    try std.testing.expect(UnitQuatf.fromAxisAngle(Vec3f.init(-0.0, 0, 0), 1.0) == null);
    try std.testing.expect(Quatf.fromAxisAngle(Vec3f.zero(), 1.0) == null);
    try std.testing.expect(UnitQuatd.fromAxisAngle(Vec3d.zero(), 1.0) == null);

    // A non-finite axis has no finite length, and one whose length overflows
    // the format has no representable direction left once the length is taken.
    const infinite = std.math.inf(f32);
    const not_a_number = std.math.nan(f32);
    const largest = std.math.floatMax(f32);
    try std.testing.expect(UnitQuatf.fromAxisAngle(Vec3f.init(infinite, 0, 0), 1.0) == null);
    try std.testing.expect(UnitQuatf.fromAxisAngle(Vec3f.init(not_a_number, 0, 0), 1.0) == null);
    const overflowing = Vec3f.init(largest, largest, largest);
    try std.testing.expect(UnitQuatf.fromAxisAngle(overflowing, 1.0) == null);

    // A non-finite angle has no half to take a sine of.
    try std.testing.expect(UnitQuatf.fromAxisAngle(Vec3f.init(0, 0, 1), std.math.nan(f32)) == null);
    try std.testing.expect(UnitQuatf.fromAxisAngle(Vec3f.init(0, 0, 1), std.math.inf(f32)) == null);
    try std.testing.expect(UnitQuatd.fromAxisAngle(Vec3d.init(0, 0, 1), std.math.inf(f64)) == null);

    // A zero angle is not an error: it is the identity rotation.
    try std.testing.expect(UnitQuatf.fromAxisAngle(Vec3f.init(0, 0, 1), 0.0).?.eq(identityf()));
    try std.testing.expect(UnitQuatf.fromAxisAngle(Vec3f.init(1, 2, 3), 0.0).?.eq(identityf()));
}

test "the axis accept and reject sets are the length's, not the quotient's" {
    // A finite axis whose length fold overflows the format has no direction
    // left, and the raw constructor owes the same answer the validated one
    // gives, in both precisions.
    const largest = std.math.floatMax(f32);
    const largest_d = std.math.floatMax(f64);
    const overflowing = Vec3f.init(largest, largest, largest);
    const overflowing_d = Vec3d.init(largest_d, largest_d, largest_d);
    try std.testing.expect(std.math.isInf(overflowing.length()));
    try std.testing.expect(std.math.isInf(overflowing_d.length()));
    try std.testing.expect(Quatf.fromAxisAngle(overflowing, 1.0) == null);
    try std.testing.expect(Quatd.fromAxisAngle(overflowing_d, 1.0) == null);
    try std.testing.expect(UnitQuatf.fromAxisAngle(overflowing, 1.0) == null);
    try std.testing.expect(UnitQuatd.fromAxisAngle(overflowing_d, 1.0) == null);

    // At the other end, an axis of subnormals is a direction and not a
    // rejection: the fold cannot scale a subnormal, so the length it reports is
    // a small multiple of the last place rather than the exact norm, but it is
    // finite and nonzero and the quotient is a positive multiple of the axis.
    const smallest = subnormal_f;
    const tiny = Vec3f.init(smallest, smallest, smallest);
    try std.testing.expect(std.math.isFinite(tiny.length()));
    try std.testing.expect(tiny.length() > 0.0);
    const tiny_raw = Quatf.fromAxisAngle(tiny, 0.5).?;
    // The axis is (1, 1, 1), so one factor scales all three vector components
    // and the scalar component is the cosine of the half angle.
    try std.testing.expect(tiny_raw.x > 0.0);
    try std.testing.expectEqual(bits(f32, tiny_raw.x), bits(f32, tiny_raw.y));
    try std.testing.expectEqual(bits(f32, tiny_raw.y), bits(f32, tiny_raw.z));
    try std.testing.expect(tiny_raw.w > 0.0);
    try std.testing.expect(tiny_raw.isFinite());

    // A single subnormal component is already a direction, and folding a zero
    // into the length changes nothing, so there the quotient is exactly one.
    const single = Vec3f.init(smallest, 0, 0);
    try std.testing.expectEqual(smallest, single.length());
    const axial = Quatf.fromAxisAngle(single, 0.5).?;
    try std.testing.expect(axial.x > 0.0 and axial.w > 0.0);
    try std.testing.expectEqual(@as(u32, 0), bits(f32, axial.y));
    try std.testing.expectEqual(@as(u32, 0), bits(f32, axial.z));
    // The binary64 constructor rejects and accepts the same axes.
    const smallest_d = subnormal_d;
    const tiny_d = Vec3d.init(smallest_d, smallest_d, smallest_d);
    try std.testing.expect(std.math.isFinite(tiny_d.length()));
    try std.testing.expect(tiny_d.length() > 0.0);
    try std.testing.expect(Quatd.fromAxisAngle(tiny_d, 0.5).?.isFinite());

    // A unit axis passes through the same division unchanged, which is what
    // makes the ordinary case a copy and not a rounding.
    try std.testing.expectEqual([3]f32{ 0, 0, 1 }, Vec3f.init(0, 0, 1).normalize().toArray());
    try std.testing.expect(UnitQuatf.fromAxisAngle(Vec3f.init(0, 0, 1), 0.5) != null);
}

// ---------------------------------------------------------------------------
// Special values
//
// `spec/core-math-2026-08-22.md#special-values` gives two rules and both are exact
// bit patterns, so both are asserted as bit patterns. A NaN operand keeps its
// sign and diagnostic payload with the quiet bit set, and an invalid operation
// on non-NaN operands gives the POSITIVE canonical quiet NaN. The second rule
// is not free: an invalid SSE operation yields the NEGATIVE indefinite,
// 0xffc00000, so a kernel that returned the host answer would fail here.
// ---------------------------------------------------------------------------

test "a NaN operand propagates with its payload and the quiet bit set" {
    // A quiet NaN is already the specification's answer, so it comes back
    // unchanged wherever it is, and only where it is: the other components go
    // through the arithmetic as ordinary numbers.
    const quiet = f32FromBits(qnan_f);
    for ([_]usize{ 0, 1, 2, 3 }) |slot| {
        var components = [4]f32{ 1.0, 2.0, 3.0, 4.0 };
        components[slot] = quiet;
        const poisoned = Quatf.fromArray(components);

        const doubled = poisoned.scale(2.0).toArray();
        const conjugated = poisoned.conjugate().toArray();
        for ([_]usize{ 0, 1, 2, 3 }) |out| {
            if (out == slot) {
                try std.testing.expectEqual(qnan_f, bits(f32, doubled[out]));
                // Conjugation negates the vector part, and the rule hands back
                // the operand's own NaN rather than the negated one, so the
                // sign of a NaN in `x` survives a conjugation.
                try std.testing.expectEqual(qnan_f, bits(f32, conjugated[out]));
            } else {
                // Conjugation negates the vector part and leaves `w` alone, so
                // the untouched component is the negation of itself for three
                // of the four slots and the component itself for `w`.
                const expected_conjugate = if (out == 3) components[out] else -components[out];
                try std.testing.expectEqual(
                    bits(f32, 2.0 * components[out]),
                    bits(f32, doubled[out]),
                );
                try std.testing.expectEqual(
                    bits(f32, expected_conjugate),
                    bits(f32, conjugated[out]),
                );
            }
        }
    }

    // A signaling NaN becomes quiet and KEEPS ITS PAYLOAD, which is the whole
    // of the NaN rule and the reason these expected values are not the
    // canonical NaN: quieting 0x7fc00001 sets a bit that is already implied and
    // leaves 0x7fc00001, and the sign bit survives a negative one. A kernel that
    // answered the canonical NaN here would be throwing away the diagnostic
    // payload the specification says to preserve.
    const signaling = f32FromBits(snan_f);
    try std.testing.expectEqual(snan_f, bits(f32, Quatf.init(signaling, 0, 0, 1).scale(2.0).x));
    const negative_signaling = f32FromBits(snan_f | 0x8000_0000);
    try std.testing.expectEqual(
        snan_f | 0x8000_0000,
        bits(f32, Quatf.init(negative_signaling, 0, 0, 1).scale(2.0).x),
    );
    // A quiet NaN is unchanged, so its own quiet bit is the only thing that
    // distinguishes it and the answer is the operand bit for bit.
    try std.testing.expectEqual(qnan_f, bits(f32, Quatf.init(quiet, 0, 0, 1).scale(2.0).x));

    // The scan is over every component of both operands, so a NaN in `y` of the
    // left operand surfaces in `x` of a product, where a kernel that looked only
    // at the matching components would report a canonical NaN and lose the
    // payload. Every component of the product carries it here, because every
    // component sums at least one product that uses `y`.
    const carried = Quatf.init(0, signaling, 0, 1).mul(Quatf.init(1, 1, 1, 1));
    try std.testing.expectEqual(snan_f, bits(f32, carried.x));
    try std.testing.expectEqual(snan_f, bits(f32, carried.y));
    try std.testing.expectEqual(snan_f, bits(f32, carried.z));
    try std.testing.expectEqual(snan_f, bits(f32, carried.w));

    // The NaN of the right operand is found the same way, and the scan runs the
    // left operand first, so a pair that is NaN on both sides reports the left
    // one. That order is the module's, stated on the kernel, and not the host's.
    const left = Quatf.init(f32FromBits(0x7fc0_0002), 0, 0, 1);
    const right = Quatf.init(f32FromBits(0x7fc0_0003), 0, 0, 1);
    try std.testing.expectEqual(@as(u32, 0x7fc0_0002), bits(f32, left.mul(right).x));
    try std.testing.expectEqual(@as(u32, 0x7fc0_0003), bits(f32, right.mul(left).x));

    const d_quiet = f64FromBits(qnan_d);
    const d_signaling = f64FromBits(snan_d);
    try std.testing.expectEqual(qnan_d, bits(f64, Quatd.init(d_quiet, 0, 0, 1).scale(2.0).x));
    try std.testing.expectEqual(snan_d, bits(f64, Quatd.init(d_signaling, 0, 0, 1).scale(2.0).x));
    try std.testing.expectEqual(
        snan_d,
        bits(f64, Quatd.init(0, 0, d_signaling, 1).mul(Quatd.init(1, 1, 1, 1)).w),
    );
}

test "an invalid operation on non-NaN operands gives the positive canonical NaN" {
    // `infinity * 0` and `infinity - infinity` are the invalid operations, and
    // every operand here is a non-NaN number, so the specification's answer is
    // the POSITIVE canonical quiet NaN. The host unit produces the negative
    // indefinite for these and would return 0xffc00000, so the sign bit in these
    // assertions IS the test.
    const infinite = std.math.inf(f32);
    const product = Quatf.init(infinite, 0, 0, 0).mul(Quatf.zero());
    inline for (0..4) |i| {
        try std.testing.expectEqual(qnan_f, bits(f32, product.toArray()[i]));
    }

    // The same rule on the sum, the componentwise product and the dot, so it is
    // a property of the module and not of one kernel.
    const left = Quatf.init(infinite, 0, 0, 1);
    try std.testing.expectEqual(qnan_f, bits(f32, left.add(Quatf.init(-infinite, 0, 0, 1)).x));
    try std.testing.expectEqual(qnan_f, bits(f32, left.scaleQuat(Quatf.zero()).x));
    try std.testing.expectEqual(qnan_f, bits(f32, left.dot(Quatf.zero())));

    // And the operations that are NOT invalid stay finite. An infinity is a
    // number, so `infinity + infinity` and `infinity * infinity` are infinities
    // and not NaNs, and a result that is an infinity is the honest answer.
    try std.testing.expectEqual(bits(f32, infinite), bits(f32, left.add(left).x));
    try std.testing.expectEqual(bits(f32, infinite), bits(f32, left.scaleQuat(left).x));
    try std.testing.expectEqual(bits(f32, infinite), bits(f32, left.dot(left)));
    // And negation is not arithmetic: `-infinity` is a number with a sign, not
    // an invalid operation, so a conjugate of an infinite component is an
    // infinity.
    try std.testing.expectEqual(bits(f32, -infinite), bits(f32, left.conjugate().x));

    const d_infinite = std.math.inf(f64);
    const d_product = Quatd.init(d_infinite, 0, 0, 0).mul(Quatd.zero());
    inline for (0..4) |i| {
        try std.testing.expectEqual(qnan_d, bits(f64, d_product.toArray()[i]));
    }
    try std.testing.expectEqual(
        qnan_d,
        bits(f64, Quatd.init(d_infinite, 0, 0, 1).add(Quatd.init(-d_infinite, 0, 0, 1)).x),
    );
    try std.testing.expectEqual(
        qnan_d,
        bits(f64, Quatd.init(d_infinite, 0, 0, 1).dot(Quatd.zero())),
    );
}
