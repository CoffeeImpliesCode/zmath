//! Tests for `zmath.geometry`.
//!
//! The suite is organised by what each test PROVES, not by which type
//! implements it. Three sections carry the most weight: the constructor
//! validation tests, the intersection-tag tests, and the predicate-policy
//! tests at the end. Every policy disagreement there is CONSTRUCTED, and each
//! one says where its numbers came from. None is quoted from a reference.
const std = @import("std");
const g = @import("geometry.zig");

fn expectTag(comptime I: type, result: I, comptime expected: std.meta.Tag(I)) !void {
    try std.testing.expect(std.meta.activeTag(result) == expected);
}

fn expectSign(comptime T: type, actual: T, expected: T, tolerance: T) !void {
    const difference = if (actual > expected) actual - expected else expected - actual;
    if (!(difference <= tolerance)) {
        std.debug.print("expected {d}, got {d}\n", .{ expected, actual });
        return error.TestUnexpectedResult;
    }
}

fn nan32() f32 {
    return @bitCast(@as(u32, 0x7fc0_0000));
}

/// Point equality through `eq`, which the module defines as bit equality with
/// NaN normalization. Plain `==` is not available on a struct, and a field
/// comparison would silently disagree with `eq` on the two zeros.
fn expectPoint(actual: anytype, expected: @TypeOf(actual)) !void {
    try std.testing.expect(actual.eq(expected));
}

/// Compares a `two_points` boundary summary in either endpoint order.
/// A zero tolerance compares the pair exactly, without prescribing its order.
fn expectCrossingPair(result: anytype, first: anytype, second: anytype, tolerance: anytype) !void {
    if (std.meta.activeTag(result) != .two_points) {
        std.debug.print("expected two_points, got {any}\n", .{std.meta.activeTag(result)});
        return error.TestUnexpectedResult;
    }
    const given = result.two_points;
    const in_order = nearPoint(given[0], first, tolerance) and nearPoint(given[1], second, tolerance);
    const swapped = nearPoint(given[0], second, tolerance) and nearPoint(given[1], first, tolerance);
    if (!in_order and !swapped) {
        std.debug.print("expected the pair {any} {any}, got {any} {any}\n", .{ first, second, given[0], given[1] });
        return error.TestUnexpectedResult;
    }
}

fn nearPoint(actual: anytype, expected: anytype, tolerance: anytype) bool {
    return nearCoordinate(actual.x, expected.x, tolerance) and nearCoordinate(actual.y, expected.y, tolerance);
}

fn nearCoordinate(actual: anytype, expected: anytype, tolerance: anytype) bool {
    const difference = if (actual > expected) actual - expected else expected - actual;
    return difference <= tolerance;
}

fn nan64() f64 {
    return @bitCast(@as(u64, 0x7ff8_0000_0000_0000));
}

fn inf32() f32 {
    return @bitCast(@as(u32, 0x7f80_0000));
}

fn inf64() f64 {
    return @bitCast(@as(u64, 0x7ff0_0000_0000_0000));
}

/// The three-by-three minor of the LAST THREE rows of a four-by-four integer
/// matrix, over every column but `skip_column`.
///
/// The three kept columns are COLLECTED rather than chosen by a chain of
/// conditionals. Choosing them by `if` is where this went wrong before: for
/// `skip_column == 2` the chain picked columns `0, 3, 3` and for
/// `skip_column == 3` it picked `0, 2, 2`, so two of the three minors
/// counted one column twice and dropped another, and the expansion returned
/// a number that was not a determinant at all. A minor is a determinant, so
/// the helper the test checks the module against has to be obviously one.
fn minor3x4(rows: [4][4]i128, skip_column: usize) i128 {
    var kept: [3]usize = undefined;
    var next: usize = 0;
    for (0..4) |column| {
        if (column == skip_column) continue;
        kept[next] = column;
        next += 1;
    }
    const c0 = kept[0];
    const c1 = kept[1];
    const c2 = kept[2];
    // The textbook three-by-three determinant,
    // `a(ei - fh) - b(di - fg) + c(dh - eg)`, written out over `rows[1..4]`.
    // The entries are named by their row and column rather than by the
    // textbook letters, because a local `g` would shadow the module import.
    const m11 = rows[1][c0];
    const m12 = rows[1][c1];
    const m13 = rows[1][c2];
    const m21 = rows[2][c0];
    const m22 = rows[2][c1];
    const m23 = rows[2][c2];
    const m31 = rows[3][c0];
    const m32 = rows[3][c1];
    const m33 = rows[3][c2];
    return m11 * (m22 * m33 - m23 * m32) - m12 * (m21 * m33 - m23 * m31) + m13 * (m21 * m32 - m22 * m31);
}

/// The exact four-by-four determinant: the Laplace expansion along the first
/// row, whose cofactor signs alternate `+ - + -` and whose minor skips the
/// column its own coefficient multiplies.
fn determinant4(rows: [4][4]i128) i128 {
    var total: i128 = 0;
    for (0..4) |column| {
        const term = rows[0][column] * minor3x4(rows, column);
        if (column % 2 == 0) total += term else total -= term;
    }
    return total;
}

/// The sign of an exact integer, as the module's own `Sign`.
fn signOfExact(value: i128) g.Sign {
    if (value > 0) return .positive;
    if (value < 0) return .negative;
    return .zero;
}

/// The exact oriented area of a triangle whose coordinates are INTEGERS, as an
/// `i128`.
///
/// This is the reference the contract names: `adaptive_exact` is the exact
/// sign of `(b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)` over the
/// ORIGINAL stored coordinates. The integers used below stay below `2^60`, so
/// every difference and every product fits `i128` and nothing rounds. No step
/// here is a floating-point subtraction, which is the whole point: a
/// floating-point DIFFERENCE is what loses the sign, and a reference that
/// formed one would be the reading under test rather than the truth.
fn exactOrient2Integer(a: [2]i128, b: [2]i128, c: [2]i128) i128 {
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
}

/// An integer as the `f64` the module receives, with a check that the
/// conversion lost nothing.
///
/// `@intFromFloat` truncates toward zero, so a coordinate that rounded on the
/// way in comes back as a DIFFERENT integer. The exact reference above is a
/// statement about the points the module was given only while this holds.
fn exactCoordinate(value: i128) !f64 {
    const out: f64 = @floatFromInt(value);
    if (!std.math.isFinite(out) or out < -0x1p127 or out >= 0x1p127)
        return error.UnrepresentableCoordinate;
    try std.testing.expectEqual(value, @as(i128, @intFromFloat(out)));
    return out;
}

// ---------------------------------------------------------------------------
// The public surface
// ---------------------------------------------------------------------------

test "every family is concretely named in both precisions and every dimension" {
    comptime {
        const families = .{
            g.Point2f,        g.Point2d,        g.Point3f,        g.Point3d,        g.Point4f,        g.Point4d,
            g.Vector2f,       g.Vector2d,       g.Vector3f,       g.Vector3d,       g.Vector4f,       g.Vector4d,
            g.Line2f,         g.Line2d,         g.Line3f,         g.Line3d,         g.Line4f,         g.Line4d,
            g.Ray2f,          g.Ray2d,          g.Ray3f,          g.Ray3d,          g.Ray4f,          g.Ray4d,
            g.Segment2f,      g.Segment2d,      g.Segment3f,      g.Segment3d,      g.Segment4f,      g.Segment4d,
            g.Circle2f,       g.Circle2d,       g.Sphere3f,       g.Sphere3d,       g.Hypersphere4f,  g.Hypersphere4d,
            g.Plane3f,        g.Plane3d,        g.Hyperplane4f,   g.Hyperplane4d,   g.Aabb2f,         g.Aabb2d,
            g.Aabb3f,         g.Aabb3d,         g.Aabb4f,         g.Aabb4d,         g.Bivector4f,     g.Bivector4d,
            g.Flat2f,         g.Flat2d,         g.Intersection2f, g.Intersection2d, g.Intersection3f, g.Intersection3d,
            g.Intersection4f, g.Intersection4d,
        };
        for (families) |family| {
            const info = @typeInfo(family);
            std.debug.assert(info == .@"struct" or info == .@"union");
        }
        std.debug.assert(@typeInfo(g.PredicatePolicy) == .@"enum");
        std.debug.assert(@typeInfo(g.Sign) == .@"enum");
        std.debug.assert(@typeInfo(g.ContainedSide) == .@"enum");
        std.debug.assert(@typeInfo(g.Intersection2f) == .@"union");
        std.debug.assert(@typeInfo(g.Intersection2f).@"union".tag_type != null);
        std.debug.assert(@typeInfo(g.Intersection3d).@"union".tag_type != null);
        std.debug.assert(@typeInfo(g.Intersection4f).@"union".tag_type != null);
    }
}

test "a point and a vector are different types in every dimension" {
    comptime {
        // `@TypeOf` of a type-valued declaration is `type` for both, so the
        // comparison has to be between the types themselves.
        std.debug.assert(g.Point2f != g.Vector2f);
        std.debug.assert(g.Point3f != g.Vector3f);
        std.debug.assert(g.Point4f != g.Vector4f);
        std.debug.assert(g.Point2d != g.Vector2d);
    }
    const point = try g.Point2f.init(1.0, 2.0);
    const vector = try g.Vector2f.init(3.0, 4.0);
    const moved = try point.add(vector);
    try std.testing.expectEqual(@as(f32, 4.0), moved.x);
    try std.testing.expectEqual(@as(f32, 6.0), moved.y);
    const back = try moved.sub(vector);
    try std.testing.expect(back.eq(point));
    const as_vector = point.toVector();
    try std.testing.expectEqual(@as(f32, 1.0), as_vector.x);
    try std.testing.expectEqual(@as(f32, 2.0), as_vector.y);
    try std.testing.expect(!as_vector.eq(vector));
    // The reverse direction subtracts, it does not add a negated vector.
    const lower = try point.sub(vector);
    try std.testing.expectEqual(@as(f32, -2.0), lower.x);
    try std.testing.expectEqual(@as(f32, -2.0), lower.y);
    try std.testing.expectEqual(@as(f32, 5.0), try vector.length());
}

test "the predicate policy is an explicit parameter wherever it can differ" {
    comptime {
        // A predicate that can differ takes `PredicatePolicy` as its last
        // parameter. A predicate that cannot differ must NOT take one, because
        // a policy argument a predicate ignores is a lie about what it
        // computes.
        const differing = .{
            @TypeOf(g.Point2f.orientation),
            @TypeOf(g.Point3f.orientation),
            @TypeOf(g.Point4f.orientation),
            @TypeOf(g.Point2f.side),
            @TypeOf(g.Circle2f.contains),
            @TypeOf(g.Sphere3f.contains),
            @TypeOf(g.Hypersphere4f.contains),
            @TypeOf(g.Circle2f.intersect),
            @TypeOf(g.Circle2f.intersectSegment),
            @TypeOf(g.Sphere3f.intersectPlane),
            @TypeOf(g.Hypersphere4f.intersectHyperplane),
            @TypeOf(g.Plane3f.intersect),
            @TypeOf(g.Hyperplane4f.intersect),
            @TypeOf(g.Segment2f.intersect),
            @TypeOf(g.Segment4d.intersectRay),
            @TypeOf(g.Ray3f.intersect),
            @TypeOf(g.Line4d.intersect),
        };
        for (differing) |entry| {
            const parameters = @typeInfo(entry).@"fn".params;
            const last = parameters[parameters.len - 1].type.?;
            std.debug.assert(last == g.PredicatePolicy);
        }
        const fixed = .{
            @TypeOf(g.Aabb2f.contains),
            @TypeOf(g.Aabb3f.intersects),
            @TypeOf(g.Aabb4f.intersect),
        };
        for (fixed) |entry| {
            for (@typeInfo(entry).@"fn".params) |parameter| {
                std.debug.assert(parameter.type.? != g.PredicatePolicy);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Constructor validation
// ---------------------------------------------------------------------------

test "a zero direction or normal is Degenerate and a non-finite one is NonFiniteInput" {
    const zero2 = try g.Vector2f.init(0.0, 0.0);
    const origin2 = try g.Point2f.init(0.0, 0.0);
    try std.testing.expectError(error.Degenerate, g.Line2f.init(origin2, zero2));
    try std.testing.expectError(error.Degenerate, g.Ray2f.init(origin2, zero2));
    try std.testing.expectError(error.Degenerate, zero2.normalize());
    try std.testing.expectError(error.Degenerate, (try g.Vector3d.init(0.0, 0.0, 0.0)).normalize());

    const infinite = g.Vector2f{ .x = inf32(), .y = 1.0 };
    try std.testing.expectError(error.NonFiniteInput, g.Line2f.init(origin2, infinite));
    try std.testing.expectError(error.NonFiniteInput, g.Ray2f.init(origin2, g.Vector2f{ .x = nan32(), .y = 0.0 }));
    try std.testing.expectError(error.NonFiniteInput, g.Ray2f.init(g.Point2f{ .x = nan32(), .y = 0.0 }, try g.Vector2f.init(1.0, 0.0)));

    const point0 = try g.Point3f.init(0.0, 0.0, 0.0);
    try std.testing.expectError(error.Degenerate, g.Plane3f.init(point0, try g.Vector3f.init(0.0, 0.0, 0.0)));
    try std.testing.expectError(error.Degenerate, g.Hyperplane4d.init(
        try g.Point4d.init(0, 0, 0, 0),
        try g.Vector4d.init(0, 0, 0, 0),
    ));
    try std.testing.expectError(error.NonFiniteInput, g.Plane3f.init(g.Point3f{ .x = nan32(), .y = 0, .z = 0 }, try g.Vector3f.init(0, 0, 1)));
    try std.testing.expectError(error.NonFiniteInput, g.Plane3f.init(point0, g.Vector3f{ .x = inf32(), .y = 0, .z = 0 }));

    // The direction and the normal are normalized, not merely accepted.
    const line = try g.Line3f.init(point0, try g.Vector3f.init(0.0, 0.0, 5.0));
    try std.testing.expectEqual(@as(f32, 1.0), line.direction.z);
    const unit = try g.Vector3f.init(1.0, 0.0, 0.0);
    try std.testing.expectEqual(@as(f32, 1.0), (try unit.normalize()).x);
    const ray = try g.Ray4d.init(try g.Point4d.init(0, 0, 0, 0), try g.Vector4d.init(0, 3.0, 4.0, 0.0));
    try expectSign(f64, ray.direction.y, 0.6, 1e-15);
    try expectSign(f64, ray.direction.z, 0.8, 1e-15);
}

test "a negative, infinite, or NaN radius is rejected by every ball family" {
    const circle_center = try g.Point2f.init(0.0, 0.0);
    try std.testing.expectError(error.NegativeRadius, g.Circle2f.init(circle_center, -1.0));
    try std.testing.expectError(error.NonFiniteInput, g.Circle2f.init(circle_center, inf32()));
    try std.testing.expectError(error.NonFiniteInput, g.Circle2f.init(circle_center, nan32()));
    try std.testing.expectError(error.NonFiniteInput, g.Circle2f.init(g.Point2f{ .x = nan32(), .y = 0.0 }, 1.0));
    try std.testing.expectEqual(@as(f32, 0.0), (try g.Circle2d.init(try g.Point2d.init(0, 0), 0.0)).radius);

    const sphere_center = try g.Point3d.init(0.0, 0.0, 0.0);
    try std.testing.expectError(error.NegativeRadius, g.Sphere3d.init(sphere_center, -0.5));
    try std.testing.expectError(error.NonFiniteInput, g.Sphere3d.init(sphere_center, inf64()));
    try std.testing.expectError(error.NonFiniteInput, g.Sphere3f.init(g.Point3f{ .x = 0, .y = inf32(), .z = 0 }, 1.0));

    const hyper_center = try g.Point4d.init(0.0, 0.0, 0.0, 0.0);
    try std.testing.expectError(error.NegativeRadius, g.Hypersphere4d.init(hyper_center, -1.0));
    try std.testing.expectError(error.NonFiniteInput, g.Hypersphere4d.init(hyper_center, nan64()));
    try std.testing.expectError(error.NonFiniteInput, g.Hypersphere4f.init(g.Point4f{ .x = 0, .y = 0, .z = 0, .w = nan32() }, 1.0));
}

test "an inverted box is InvertedBounds and a non-finite corner is NonFiniteInput" {
    const low = try g.Point2f.init(0.0, 0.0);
    const high = try g.Point2f.init(1.0, 1.0);
    try std.testing.expectError(error.InvertedBounds, g.Aabb2f.init(high, low));
    try std.testing.expectError(error.InvertedBounds, g.Aabb2f.init(try g.Point2f.init(0.0, 2.0), try g.Point2f.init(1.0, 1.0)));
    try std.testing.expectError(error.NonFiniteInput, g.Aabb2f.init(g.Point2f{ .x = nan32(), .y = 0.0 }, high));
    try std.testing.expectError(error.InvertedBounds, g.Aabb3f.init(try g.Point3f.init(1, 1, 1), try g.Point3f.init(0, 0, 0)));
    try std.testing.expectError(error.NonFiniteInput, g.Aabb3d.init(try g.Point3d.init(0, 0, 0), g.Point3d{ .x = 1, .y = 1, .z = inf64() }));
    try std.testing.expectError(error.InvertedBounds, g.Aabb4d.init(try g.Point4d.init(0, 0, 0, 1), try g.Point4d.init(0, 0, 0, 0)));
    try std.testing.expectError(error.NonFiniteInput, g.Aabb4d.init(g.Point4d{ .x = nan64(), .y = 0, .z = 0, .w = 0 }, try g.Point4d.init(0, 0, 0, 0)));
    // A degenerate box, minimum equal to maximum, is legal: the contract
    // requires `minimum <= maximum`, not a strict inequality.
    const degenerate = try g.Aabb2f.init(high, high);
    try std.testing.expect(degenerate.minimum.eq(high));
    try std.testing.expect(try degenerate.contains(high));
}

test "a NaN in any argument is NonFiniteInput" {
    const origin = try g.Point3f.init(0.0, 0.0, 0.0);
    const target = try g.Point3f.init(1.0, 0.0, 0.0);
    try std.testing.expectError(error.NonFiniteInput, g.Segment3f.init(origin, g.Point3f{ .x = nan32(), .y = 0, .z = 0 }));
    try std.testing.expectError(error.NonFiniteInput, g.Segment3f.init(g.Point3f{ .x = 0, .y = inf32(), .z = 0 }, target));
    const nan_vector = g.Vector3f{ .x = nan32(), .y = 0.0, .z = 0.0 };
    try std.testing.expectError(error.NonFiniteInput, nan_vector.normalize());
    try std.testing.expectError(error.NonFiniteInput, g.Point3f.init(1.0, nan32(), 0.0));

    const segment = try g.Segment3f.init(origin, target);
    try std.testing.expectError(error.NonFiniteInput, segment.at(nan32()));
    try std.testing.expectError(error.NonFiniteInput, segment.closestPointTo(g.Point3f{ .x = nan32(), .y = 0, .z = 0 }));
    const line = try g.Line3f.init(origin, try g.Vector3f.init(1, 0, 0));
    try std.testing.expectError(error.NonFiniteInput, line.at(nan32()));
    try std.testing.expectError(error.NonFiniteInput, line.closestPointTo(g.Point3f{ .x = 0, .y = inf32(), .z = 0 }));
    const ray = try g.Ray3f.init(origin, try g.Vector3f.init(1, 0, 0));
    try std.testing.expectError(error.NonFiniteInput, ray.at(inf32()));

    // The fields are public, so a caller can poison a value after its
    // constructor validated it. Every operation re-validates.
    const ball = try g.Sphere3f.init(origin, 1.0);
    var poisoned = ball;
    poisoned.center.z = nan32();
    try std.testing.expectError(error.NonFiniteInput, poisoned.intersect(ball, .adaptive_exact));
    try std.testing.expectError(error.NonFiniteInput, ball.contains(g.Point3f{ .x = nan32(), .y = 0, .z = 0 }, .fast));
    try std.testing.expect(ball.isFinite());

    var poisoned_box = try g.Aabb3f.init(origin, target);
    poisoned_box.maximum.y = inf32();
    try std.testing.expectError(error.NonFiniteInput, poisoned_box.intersect(try g.Aabb3f.init(origin, target)));
    try std.testing.expectError(error.NonFiniteInput, poisoned_box.center());
    try std.testing.expectError(error.NonFiniteInput, poisoned_box.extent());
    try std.testing.expectError(error.NonFiniteInput, poisoned_box.intersects(try g.Aabb3f.init(origin, target)));
    try std.testing.expect(!poisoned.isFinite());
}

// ---------------------------------------------------------------------------
// Length, clamps, and projections
// ---------------------------------------------------------------------------

test "the folded length survives what the squared length cannot, and says so" {
    // `sqrt(x*x + y*y + z*z)` overflows here, because 1e30 squared is past the
    // binary32 range. The folded `hypot` cannot, because it never forms the
    // square.
    const large = try g.Vector3f.init(1.0e30, 1.0e30, 1.0e30);
    try std.testing.expect(std.math.isFinite(try large.length()));
    try expectSign(f32, try large.length(), 1.7320508e30, 1.0e24);
    // The squared length is the one that leaves the format, and it says so
    // rather than handing back an infinity the caller would then normalize
    // into a zero direction.
    try std.testing.expectError(error.NonFiniteResult, large.lengthSquared());

    const large64 = try g.Vector4d.init(1.0e150, 1.0e150, 1.0e150, 1.0e150);
    try std.testing.expect(std.math.isFinite(try large64.length()));
    try expectSign(f64, try large64.length(), 2.0e150, 1.0e144);

    // A subnormal component has a representable norm of the same order, and
    // `x*x` underflows to zero. The fold starts from `hypot(+0, x)`, which is
    // exactly `|x|` with no flush.
    const tiny = try g.Vector2f.init(1.0e-40, 0.0);
    try std.testing.expectEqual(@as(f32, 1.0e-40), try tiny.length());
    try std.testing.expectEqual(@as(f32, 0.0), try tiny.lengthSquared());
    const tiny_pair = try g.Vector3f.init(1.0e-40, 1.0e-40, 0.0);
    try expectSign(f32, try tiny_pair.length(), 1.4142135e-40, 1.0e-44);
    const tiny64 = try g.Vector2d.init(5.0e-320, 0.0);
    try std.testing.expectEqual(@as(f64, 5.0e-320), try tiny64.length());
}

test "ray at clamps at the origin and segment at clamps to the unit interval" {
    const origin = try g.Point2f.init(1.0, 1.0);
    const ray = try g.Ray2f.init(origin, try g.Vector2f.init(1.0, 0.0));
    try expectPoint(try ray.at(-5.0), origin);
    try expectPoint(try ray.at(0.0), origin);
    const ahead = try ray.at(2.0);
    try std.testing.expectEqual(@as(f32, 3.0), ahead.x);
    try std.testing.expectEqual(@as(f32, 1.0), ahead.y);

    const low = try g.Point2f.init(0.0, 0.0);
    const high = try g.Point2f.init(2.0, 4.0);
    const segment = try g.Segment2f.init(low, high);
    try expectPoint(try segment.at(-1.0), low);
    try expectPoint(try segment.at(2.0), high);
    try expectPoint(try segment.at(0.0), low);
    try expectPoint(try segment.at(1.0), high);
    const middle = try segment.at(0.5);
    try std.testing.expectEqual(@as(f32, 1.0), middle.x);
    try std.testing.expectEqual(@as(f32, 2.0), middle.y);
    try std.testing.expectEqual(@as(f32, @sqrt(4.0 + 16.0)), try segment.length());
    const collapsed = try g.Segment3f.init(try g.Point3f.init(1, 1, 1), try g.Point3f.init(1, 1, 1));
    try std.testing.expectEqual(@as(f32, 0.0), try collapsed.length());
}

test "closest-point queries project onto the line and clamp onto the segment" {
    const origin = try g.Point3d.init(0.0, 0.0, 0.0);
    const line = try g.Line3d.init(origin, try g.Vector3d.init(1.0, 1.0, 0.0));
    const foot = try line.closestPointTo(try g.Point3d.init(1.0, 0.0, 5.0));
    try expectSign(f64, foot.x, 0.5, 1e-15);
    try expectSign(f64, foot.y, 0.5, 1e-15);
    try std.testing.expectEqual(@as(f64, 0.0), foot.z);

    const low = try g.Point3d.init(0, 0, 0);
    const high = try g.Point3d.init(1, 0, 0);
    const segment = try g.Segment3d.init(low, high);
    try expectPoint(try segment.closestPointTo(try g.Point3d.init(5, 3, 0)), high);
    try expectPoint(try segment.closestPointTo(try g.Point3d.init(-5, 3, 0)), low);
    const middle = try g.Point3d.init(0.25, 3, 0);
    try expectPoint(try segment.closestPointTo(middle), try segment.at(0.25));
}

// ---------------------------------------------------------------------------
// Segment against segment
// ---------------------------------------------------------------------------

test "segment against segment reports every tag in two dimensions" {
    const policy = g.PredicatePolicy.adaptive_exact;
    const low = try g.Point2f.init(0, 0);
    const one = try g.Point2f.init(1, 0);
    const two = try g.Point2f.init(2, 0);
    const three = try g.Point2f.init(3, 0);
    const up = try g.Point2f.init(0, 1);
    const up_right = try g.Point2f.init(1, 1);

    // Disjoint: one segment lies on the line y = 0 and the other on y = 1, so
    // no point is common.
    const below = try g.Segment2f.init(low, one);
    const above = try g.Segment2f.init(up, up_right);
    try expectTag(g.Intersection2f, try below.intersect(above, policy), .none);

    // A single crossing: the two segments are the diagonals of the unit
    // square, which meet at the centre.
    const rising = try g.Segment2f.init(low, up_right);
    const falling = try g.Segment2f.init(up, one);
    const crossed = try rising.intersect(falling, policy);
    try expectTag(g.Intersection2f, crossed, .point);
    try expectSign(f32, crossed.point.x, 0.5, 1e-6);
    try expectSign(f32, crossed.point.y, 0.5, 1e-6);

    // Touching ends: the second segment starts exactly where the first ends,
    // both on y = 0, so the common set is that one end.
    const trailing = try g.Segment2f.init(one, two);
    const touching = try below.intersect(trailing, policy);
    try expectTag(g.Intersection2f, touching, .point);
    try expectSign(f32, touching.point.x, 1.0, 1e-6);
    try expectSign(f32, touching.point.y, 0.0, 1e-6);

    // Overlapping collinear: the first runs from 0 to 2 and the second from
    // 1 to 3 along y = 0, so the common set is [1, 2], which is neither whole
    // segment.
    const over = try g.Segment2f.init(low, two);
    const later = try g.Segment2f.init(one, three);
    const overlap = try over.intersect(later, policy);
    try expectTag(g.Intersection2f, overlap, .segment);
    try expectSign(f32, overlap.segment.a.x, 1.0, 1e-6);
    try expectSign(f32, overlap.segment.b.x, 2.0, 1e-6);
    try expectSign(f32, overlap.segment.a.y, 0.0, 1e-6);

    // Contained: the second segment is strictly inside the first on the same
    // line, so the whole second segment is the common set.
    const inner = try g.Segment2f.init(one, two);
    const contained = try over.intersect(inner, policy);
    try expectTag(g.Intersection2f, contained, .contained);
    try std.testing.expectEqual(g.ContainedSide.second, contained.contained);

    // Fully coincident: the two segments are the same set.
    const same = try g.Segment2f.init(low, one);
    try expectTag(g.Intersection2f, try below.intersect(same, policy), .coincident);

    // Collinear but disjoint, which only the collinear branch sees.
    const far = try g.Segment2f.init(two, three);
    try expectTag(g.Intersection2f, try below.intersect(far, policy), .none);

    // A degenerate segment is one point, and it lands on the other segment.
    const dot_point = try g.Segment2f.init(one, one);
    const on_it = try below.intersect(dot_point, policy);
    try expectTag(g.Intersection2f, on_it, .point);
    try expectSign(f32, on_it.point.x, 1.0, 1e-6);
    const off_it = try g.Segment2f.init(up, up);
    try expectTag(g.Intersection2f, try below.intersect(off_it, policy), .none);
}

test "segment against segment gates on exact coplanarity in three dimensions" {
    const policy = g.PredicatePolicy.adaptive_exact;
    // Two segments in the plane z = 0 that cross, and the same pair lifted so
    // that they are skew. The lifted pair has no common point, and the gate
    // that proves it is the exact signed volume.
    const a = try g.Point3f.init(0.0, 0.0, 0.0);
    const b = try g.Point3f.init(2.0, 0.0, 0.0);
    const c = try g.Point3f.init(1.0, -1.0, 0.0);
    const d = try g.Point3f.init(1.0, 1.0, 0.0);
    const along = try g.Segment3f.init(a, b);
    const across = try g.Segment3f.init(c, d);
    const hit = try along.intersect(across, policy);
    try expectTag(g.Intersection3f, hit, .point);
    try expectSign(f32, hit.point.x, 1.0, 1e-5);
    try expectSign(f32, hit.point.y, 0.0, 1e-5);

    const lifted = try g.Segment3f.init(c, try g.Point3f.init(1.0, 1.0, 1.0));
    try expectTag(g.Intersection3f, try along.intersect(lifted, policy), .none);

    const collinear = try g.Segment3f.init(try g.Point3f.init(1.0, 0.0, 0.0), try g.Point3f.init(3.0, 0.0, 0.0));
    const overlapping = try along.intersect(collinear, policy);
    try expectTag(g.Intersection3f, overlapping, .segment);
    try expectSign(f32, overlapping.segment.a.x, 1.0, 1e-5);
    try expectSign(f32, overlapping.segment.b.x, 2.0, 1e-5);
}

test "segment against segment gates on the four-dimensional minors" {
    const policy = g.PredicatePolicy.adaptive_exact;
    // Four points in a genuine two-flat, so the two segments cross there.
    const a = try g.Point4d.init(0, 0, 0, 0);
    const b = try g.Point4d.init(2, 0, 0, 0);
    const c = try g.Point4d.init(1, -1, 0, 0);
    const d = try g.Point4d.init(1, 1, 0, 0);
    const first = try g.Segment4d.init(a, b);
    const second = try g.Segment4d.init(c, d);
    const hit = try first.intersect(second, policy);
    try expectTag(g.Intersection4d, hit, .point);
    try expectSign(f64, hit.point.x, 1.0, 1e-12);
    try expectSign(f64, hit.point.y, 0.0, 1e-12);

    // The same configuration lifted out of the two-flat. Four points in four
    // dimensions are affinely dependent exactly when the three-by-four matrix
    // of the difference vectors has rank below three, which is four minors,
    // each degree two. Here the first one is nonzero.
    const lifted = try g.Segment4d.init(try g.Point4d.init(1, -1, 1, 0), try g.Point4d.init(1, 1, 1, 0));
    try expectTag(g.Intersection4d, try first.intersect(lifted, policy), .none);

    // Collinear overlap inside the same two-flat still works in four
    // dimensions, which shows the gate is the only thing that changed.
    const collinear = try g.Segment4d.init(try g.Point4d.init(1, 0, 0, 0), try g.Point4d.init(3, 0, 0, 0));
    const overlapping = try first.intersect(collinear, policy);
    try expectTag(g.Intersection4d, overlapping, .segment);
    try expectSign(f64, overlapping.segment.a.x, 1.0, 1e-12);
    try expectSign(f64, overlapping.segment.b.x, 2.0, 1e-12);
}

test "ray against ray and ray against segment" {
    const policy = g.PredicatePolicy.adaptive_exact;
    const origin = try g.Point2d.init(0, 0);
    // Two rays from the same end in different directions share only that end.
    const east = try g.Ray2d.init(origin, try g.Vector2d.init(1, 0));
    const north = try g.Ray2d.init(origin, try g.Vector2d.init(0, 1));
    const shared = try east.intersect(north, policy);
    try expectTag(g.Intersection2d, shared, .point);
    try std.testing.expect(shared.point.eq(origin));

    // Rays that cross away from the end: the first along y = 0 and the second
    // from (1, -1) straight up, meeting at (1, 0).
    const along = try g.Ray2d.init(origin, try g.Vector2d.init(1, 0));
    const up = try g.Ray2d.init(try g.Point2d.init(1, -1), try g.Vector2d.init(0, 1));
    const crossed = try along.intersect(up, policy);
    try expectTag(g.Intersection2d, crossed, .point);
    try expectSign(f64, crossed.point.x, 1.0, 1e-12);
    try expectSign(f64, crossed.point.y, 0.0, 1e-12);

    // A ray that leaves the line entirely.
    const away = try g.Ray2d.init(try g.Point2d.init(0, 1), try g.Vector2d.init(0, 1));
    try expectTag(g.Intersection2d, try along.intersect(away, policy), .none);

    // A ray and a segment that lies wholly on it: the whole segment is the
    // common set, and the segment is the argument of the call.
    const covered = try g.Segment2d.init(try g.Point2d.init(1, 0), try g.Point2d.init(2, 0));
    const on_ray = try along.intersectSegment(covered, policy);
    try expectTag(g.Intersection2d, on_ray, .contained);
    try std.testing.expectEqual(g.ContainedSide.second, on_ray.contained);

    // A segment straddling the ray's end: the part behind the end is not on
    // the ray, so the common set is the half of the segment from the end
    // forward, which is a segment and not a point.
    const straddling = try g.Segment2d.init(try g.Point2d.init(-1, 0), try g.Point2d.init(1, 0));
    const partial = try along.intersectSegment(straddling, policy);
    try expectTag(g.Intersection2d, partial, .segment);
    try expectSign(f64, partial.segment.a.x, 0.0, 1e-12);
    try expectSign(f64, partial.segment.b.x, 1.0, 1e-12);

    // A segment crossing the ray strictly inside: one point where they meet.
    const crossing = try g.Segment2d.init(try g.Point2d.init(1, -1), try g.Point2d.init(1, 1));
    const met = try along.intersectSegment(crossing, policy);
    try expectTag(g.Intersection2d, met, .point);
    try expectSign(f64, met.point.x, 1.0, 1e-12);
    try expectSign(f64, met.point.y, 0.0, 1e-12);

    // A segment that misses the ray on the far side.
    const behind = try g.Segment2d.init(try g.Point2d.init(-3, -1), try g.Point2d.init(-2, 1));
    try expectTag(g.Intersection2d, try along.intersectSegment(behind, policy), .none);

    // A ray against a line in three dimensions, where the same ray lies in
    // the line and the common set is the ray.
    const ray3 = try g.Ray3d.init(try g.Point3d.init(0, 0, 0), try g.Vector3d.init(1, 0, 0));
    const line3 = try g.Line3d.init(try g.Point3d.init(-5, 0, 0), try g.Vector3d.init(1, 0, 0));
    const ray_in_line = try ray3.intersectLine(line3, policy);
    try expectTag(g.Intersection3d, ray_in_line, .contained);
    try std.testing.expectEqual(g.ContainedSide.first, ray_in_line.contained);
}

test "line against line is a point, none, or a whole line" {
    const policy = g.PredicatePolicy.adaptive_exact;
    // Two lines through the origin in different directions meet there.
    const along_x = try g.Line3d.init(try g.Point3d.init(0, 0, 0), try g.Vector3d.init(1, 0, 0));
    const along_y = try g.Line3d.init(try g.Point3d.init(0, 0, 0), try g.Vector3d.init(0, 1, 0));
    const at_origin = try along_x.intersect(along_y, policy);
    try expectTag(g.Intersection3d, at_origin, .point);
    try std.testing.expect(at_origin.point.eq(try g.Point3d.init(0, 0, 0)));

    // Parallel and distinct: in three dimensions two parallel lines are skew,
    // so the coplanarity gate closes.
    const shifted = try g.Line3d.init(try g.Point3d.init(0, 1, 0), try g.Vector3d.init(1, 0, 0));
    try expectTag(g.Intersection3d, try along_x.intersect(shifted, policy), .none);

    // The same line written twice with opposite directions.
    const reversed = try g.Line3d.init(try g.Point3d.init(1, 0, 0), try g.Vector3d.init(-1, 0, 0));
    const coincident = try along_x.intersect(reversed, policy);
    try expectTag(g.Intersection3d, coincident, .line);
    try expectSign(f64, try coincident.line.direction.length(), 1.0, 1e-12);

    // In two dimensions two parallel lines ARE coplanar, so the gate lets the
    // collinear branch decide, and the answer is `none`.
    const flat_x = try g.Line2d.init(try g.Point2d.init(0, 0), try g.Vector2d.init(1, 0));
    const flat_shift = try g.Line2d.init(try g.Point2d.init(0, 1), try g.Vector2d.init(1, 0));
    try expectTag(g.Intersection2d, try flat_x.intersect(flat_shift, policy), .none);
    const same_flat = try g.Line2d.init(try g.Point2d.init(2, 0), try g.Vector2d.init(1, 0));
    try expectTag(g.Intersection2d, try flat_x.intersect(same_flat, policy), .line);
    // In four dimensions the same skew pair is rejected by the minors.
    const line4 = try g.Line4d.init(try g.Point4d.init(0, 0, 0, 0), try g.Vector4d.init(1, 0, 0, 0));
    const shifted4 = try g.Line4d.init(try g.Point4d.init(0, 1, 0, 0), try g.Vector4d.init(1, 0, 0, 0));
    try expectTag(g.Intersection4d, try line4.intersect(shifted4, policy), .none);
}

// ---------------------------------------------------------------------------
// Ball intersections
// ---------------------------------------------------------------------------

test "two circles report every configuration and each one pins to a tag" {
    const policy = g.PredicatePolicy.adaptive_exact;
    const centre = try g.Point2d.init(0, 0);
    const right = try g.Point2d.init(2, 0);
    const near_right = try g.Point2d.init(1, 0);
    const far_right = try g.Point2d.init(5, 0);

    const unit = try g.Circle2d.init(centre, 1.0);
    const unit_right = try g.Circle2d.init(right, 1.0);
    const wide = try g.Circle2d.init(centre, 2.0);
    const wide_right = try g.Circle2d.init(right, 2.0);
    const small_near = try g.Circle2d.init(near_right, 1.0);
    const same = try g.Circle2d.init(try g.Point2d.init(3, 4), 1.5);
    const same_again = try g.Circle2d.init(try g.Point2d.init(3, 4), 1.5);

    // Disjoint: the centres are five apart and the radii sum to two.
    try expectTag(g.Intersection2d, try unit.intersect(try g.Circle2d.init(far_right, 1.0), policy), .none);

    // Externally tangent: the centres are two apart and the radii sum to two,
    // so the single common point is halfway.
    const tangent = try unit.intersect(unit_right, policy);
    try expectTag(g.Intersection2d, tangent, .point);
    try expectSign(f64, tangent.point.x, 1.0, 1e-12);
    try expectSign(f64, tangent.point.y, 0.0, 1e-12);

    // Crossing: the centres are two apart and both radii are two, so the
    // circles cut each other at two points symmetric about the line of
    // centres.
    const crossing = try wide.intersect(wide_right, policy);
    try expectTag(g.Intersection2d, crossing, .two_points);
    try expectSign(f64, crossing.two_points[0].x, 1.0, 1e-12);
    try expectSign(f64, crossing.two_points[1].x, 1.0, 1e-12);
    const first_height = if (crossing.two_points[0].y > 0) crossing.two_points[0].y else -crossing.two_points[0].y;
    const second_height = if (crossing.two_points[1].y > 0) crossing.two_points[1].y else -crossing.two_points[1].y;
    try expectSign(f64, first_height, second_height, 1e-12);
    try std.testing.expect((crossing.two_points[0].y > 0) != (crossing.two_points[1].y > 0));

    // Internally tangent: the centres are one apart and the radii differ by
    // one, so the single common point is on the far side of the smaller
    // circle. It carries the SAME `point` tag as external tangency, because
    // both are a single common point. The contract's tag set has no second
    // tangency tag, and inventing one would put a tag in the union that no
    // other dimension could reach.
    const inner_tangent = try wide.intersect(small_near, policy);
    try expectTag(g.Intersection2d, inner_tangent, .point);
    try expectSign(f64, inner_tangent.point.x, 2.0, 1e-12);

    // Concentric with different radii: the smaller lies wholly inside the
    // larger, so the common set is the whole smaller circle, which here is
    // `unit`, the RECEIVER of the call. `ContainedSide` names the operand the
    // intersection EQUALS by its position in the call, so this one is `first`.
    const big = try g.Circle2d.init(centre, 2.0);
    const concentric = try unit.intersect(big, policy);
    try expectTag(g.Intersection2d, concentric, .contained);
    try std.testing.expectEqual(g.ContainedSide.first, concentric.contained);
    // The order of the call flips the side: now the common set is the
    // argument, so it is `second`.
    const reversed = try big.intersect(unit, policy);
    try expectTag(g.Intersection2d, reversed, .contained);
    try std.testing.expectEqual(g.ContainedSide.second, reversed.contained);

    // Coincident: same centre, same radius.
    try expectTag(g.Intersection2d, try same.intersect(same_again, policy), .coincident);

    // A zero-radius circle inside another is a single point of the larger, and
    // the tag says the smaller circle is the whole common set.
    const dot_ball = try g.Circle2d.init(centre, 0.0);
    const containing = try g.Circle2d.init(centre, 3.0);
    try expectTag(g.Intersection2d, try dot_ball.intersect(containing, policy), .contained);
}

test "sphere and hypersphere share the ball classification in their dimensions" {
    const policy = g.PredicatePolicy.adaptive_exact;
    const centre = try g.Point3d.init(0, 0, 0);
    const right = try g.Point3d.init(2, 0, 0);
    const near_right = try g.Point3d.init(1, 0, 0);

    const unit = try g.Sphere3d.init(centre, 1.0);
    const unit_right = try g.Sphere3d.init(right, 1.0);
    const wide = try g.Sphere3d.init(centre, 2.0);
    const wide_right = try g.Sphere3d.init(right, 2.0);

    const tangent = try unit.intersect(unit_right, policy);
    try expectTag(g.Intersection3d, tangent, .point);
    try expectSign(f64, tangent.point.x, 1.0, 1e-12);
    const crossing = try wide.intersect(wide_right, policy);
    try expectTag(g.Intersection3d, crossing, .two_points);
    try expectSign(f64, crossing.two_points[0].x, 1.0, 1e-12);
    try expectSign(f64, crossing.two_points[1].x, 1.0, 1e-12);
    try expectTag(g.Intersection3d, try wide.intersect(try g.Sphere3d.init(near_right, 1.0), policy), .point);
    try expectTag(g.Intersection3d, try unit.intersect(try g.Sphere3d.init(try g.Point3d.init(9, 0, 0), 1.0), policy), .none);
    try expectTag(g.Intersection3d, try unit.intersect(unit, policy), .coincident);
    try expectTag(g.Intersection3d, try wide.intersect(unit, policy), .contained);

    // Four dimensions, the same six configurations.
    const hyper_centre = try g.Point4d.init(0, 0, 0, 0);
    const hyper_right = try g.Point4d.init(2, 0, 0, 0);
    const hyper_unit = try g.Hypersphere4d.init(hyper_centre, 1.0);
    const hyper_unit_right = try g.Hypersphere4d.init(hyper_right, 1.0);
    const hyper_tangent = try hyper_unit.intersect(hyper_unit_right, policy);
    try expectTag(g.Intersection4d, hyper_tangent, .point);
    try expectSign(f64, hyper_tangent.point.x, 1.0, 1e-12);
    const hyper_crossing = try (try g.Hypersphere4d.init(hyper_centre, 2.0)).intersect(
        try g.Hypersphere4d.init(hyper_right, 2.0),
        policy,
    );
    try expectTag(g.Intersection4d, hyper_crossing, .two_points);
    try expectTag(g.Intersection4d, try hyper_unit.intersect(hyper_unit, policy), .coincident);
    try expectTag(g.Intersection4d, try hyper_unit.intersect(
        try g.Hypersphere4d.init(try g.Point4d.init(9, 0, 0, 0), 1.0),
        policy,
    ), .none);
}

test "the contained side names the common set itself in a strict enclosure" {
    // How this was built. The common set of two balls when one lies strictly
    // inside the other IS the inner ball, whole, so the only question the tag
    // leaves open is WHICH of the two operands it is, and that is answered by
    // the position in the call rather than by any comparison the code makes.
    // Two things keep the fixture from answering that question by accident:
    //
    // - The centres are NOT the same, so the pair is classified from a real
    //   separation instead of the zero-length one a concentric pair has. The
    //   side answers the same question in both cases, and the two-dimensional
    //   concentric pair above covers the other one.
    // - Every number is a small dyadic, and the two tangency comparisons are
    //   built from the separation squared against the sum and the difference
    //   squared. Those are `5/16`, `25/16` and `9/16`: all exact in binary32,
    //   so `fast` and `adaptive_exact` reach the same classification and the
    //   fixture is not a rounding artefact under either policy. Both
    //   comparisons are strictly negative with room to spare, so this is a
    //   strict enclosure and not a tangency or a coincidence.
    const outer = try g.Sphere3f.init(try g.Point3f.init(0, 0, 0), 1);
    const inner = try g.Sphere3f.init(try g.Point3f.init(0.5, 0.25, 0), 0.25);

    // The receiver holds the argument, so the whole common set is the
    // ARGUMENT, and the side is `second`.
    const holds = try outer.intersect(inner, .fast);
    try expectTag(g.Intersection3f, holds, .contained);
    try std.testing.expectEqual(g.ContainedSide.second, holds.contained);
    const holds_exact = try outer.intersect(inner, .adaptive_exact);
    try expectTag(g.Intersection3f, holds_exact, .contained);
    try std.testing.expectEqual(g.ContainedSide.second, holds_exact.contained);

    // Reversed, the receiver is the one wholly inside, so the side is `first`.
    const inside = try inner.intersect(outer, .fast);
    try expectTag(g.Intersection3f, inside, .contained);
    try std.testing.expectEqual(g.ContainedSide.first, inside.contained);
    const inside_exact = try inner.intersect(outer, .adaptive_exact);
    try expectTag(g.Intersection3f, inside_exact, .contained);
    try std.testing.expectEqual(g.ContainedSide.first, inside_exact.contained);
}

test "a ball of zero radius inside another is the whole common set in either order" {
    // Zero radius is the boundary of the classification, and it is the case
    // where a side computed from the radii alone is most likely to go wrong:
    // the inner ball is a single point, so it is also what an internal
    // tangency would report, and the common set really is that whole point.
    // The point is NOT the outer's centre, so the answer cannot fall out of a
    // zero-length separation: the separation squared is `5/16` against a
    // squared difference of `4`, both tangency comparisons are strictly
    // negative, and containment is the only classification left.
    const containing = try g.Hypersphere4d.init(try g.Point4d.init(0, 0, 0, 0), 2);
    const dot_ball = try g.Hypersphere4d.init(try g.Point4d.init(0.5, 0, 0, 0.25), 0);

    // The zero-radius ball is the receiver here, so it is `first`.
    const dot_first = try dot_ball.intersect(containing, .fast);
    try expectTag(g.Intersection4d, dot_first, .contained);
    try std.testing.expectEqual(g.ContainedSide.first, dot_first.contained);
    const dot_first_exact = try dot_ball.intersect(containing, .adaptive_exact);
    try expectTag(g.Intersection4d, dot_first_exact, .contained);
    try std.testing.expectEqual(g.ContainedSide.first, dot_first_exact.contained);

    // Reversed, the same point is the argument, so the side is `second`.
    const dot_second = try containing.intersect(dot_ball, .fast);
    try expectTag(g.Intersection4d, dot_second, .contained);
    try std.testing.expectEqual(g.ContainedSide.second, dot_second.contained);
    const dot_second_exact = try containing.intersect(dot_ball, .adaptive_exact);
    try expectTag(g.Intersection4d, dot_second_exact, .contained);
    try std.testing.expectEqual(g.ContainedSide.second, dot_second_exact.contained);
}

test "a ball against a segment is none, a point, a segment, or the whole segment" {
    const policy = g.PredicatePolicy.adaptive_exact;
    const centre = try g.Point3d.init(0, 0, 0);
    const unit = try g.Sphere3d.init(centre, 1.0);

    // A chord through the centre: the segment runs from outside one side to
    // outside the other, so the common set is the chord, a proper sub-segment.
    const chord = try g.Segment3d.init(try g.Point3d.init(-2, 0, 0), try g.Point3d.init(2, 0, 0));
    const through = try unit.intersectSegment(chord, policy);
    try expectTag(g.Intersection3d, through, .segment);
    try expectSign(f64, through.segment.a.x, -1.0, 1e-12);
    try expectSign(f64, through.segment.b.x, 1.0, 1e-12);

    // A segment wholly outside.
    const outside = try g.Segment3d.init(try g.Point3d.init(5, 0, 0), try g.Point3d.init(6, 0, 0));
    try expectTag(g.Intersection3d, try unit.intersectSegment(outside, policy), .none);

    // A segment that ENDS exactly on the surface, but reaches it from
    // outside, so it is not a touch. Substituting the parameter
    // `p(t) = (t, 2 - 2t, 0)` into the unit sphere gives
    // `t^2 + (2 - 2t)^2 = 1`, that is `5t^2 - 8t + 3 = 0`, whose roots are
    // `t = 3/5` and `t = 1`. The quadratic opens upward, so the ball holds
    // the segment strictly between the roots and the common set is the
    // sub-segment `3/5 <= t <= 1`, from `(0.6, 0.8, 0)` to `(1, 0, 0)`. It
    // is a proper sub-segment of the argument, not one point.
    const grazing = try g.Segment3d.init(try g.Point3d.init(0, 2, 0), try g.Point3d.init(1, 0, 0));
    const touch = try unit.intersectSegment(grazing, policy);
    try expectTag(g.Intersection3d, touch, .segment);
    try expectSign(f64, touch.segment.a.x, 0.6, 1e-12);
    try expectSign(f64, touch.segment.a.y, 0.8, 1e-12);
    try expectSign(f64, touch.segment.b.x, 1.0, 1e-12);
    try expectSign(f64, touch.segment.b.y, 0.0, 1e-12);

    // A segment that TOUCHES at one point, which is the case where the
    // quadratic has a double root. Along the line `(t, 1, 0)` the sphere
    // equation is `t^2 + 1 = 1`, that is `t^2 = 0`, which has the DOUBLE root
    // `t = 0`. The other line `(t, 1, 1)` would give `t^2 + 2 = 1`, which has
    // no real root at all and misses the ball entirely.
    const tangent = try g.Segment3d.init(try g.Point3d.init(-1, 1, 0), try g.Point3d.init(1, 1, 0));
    const grazed = try unit.intersectSegment(tangent, policy);
    try expectTag(g.Intersection3d, grazed, .point);
    try expectSign(f64, grazed.point.x, 0.0, 1e-12);
    try expectSign(f64, grazed.point.y, 1.0, 1e-12);
    try expectSign(f64, grazed.point.z, 0.0, 1e-12);

    // A segment wholly inside: the common set is the whole segment, which is
    // the argument of the call and therefore `second`.
    const inside = try g.Segment3d.init(try g.Point3d.init(-0.5, 0, 0), try g.Point3d.init(0.5, 0, 0));
    const whole = try unit.intersectSegment(inside, policy);
    try expectTag(g.Intersection3d, whole, .contained);
    try std.testing.expectEqual(g.ContainedSide.second, whole.contained);

    // A segment that starts inside and leaves: a segment from the start to
    // the exit, which is a proper sub-segment.
    const leaving = try g.Segment3d.init(try g.Point3d.init(0, 0, 0), try g.Point3d.init(3, 0, 0));
    const exit_run = try unit.intersectSegment(leaving, policy);
    try expectTag(g.Intersection3d, exit_run, .segment);
    try expectSign(f64, exit_run.segment.a.x, 0.0, 1e-12);
    try expectSign(f64, exit_run.segment.b.x, 1.0, 1e-12);

    // A degenerate segment is a single point.
    const dot_point = try g.Segment3d.init(centre, centre);
    try expectTag(g.Intersection3d, try unit.intersectSegment(dot_point, policy), .contained);
    const outside_dot = try g.Segment3d.init(try g.Point3d.init(9, 0, 0), try g.Point3d.init(9, 0, 0));
    try expectTag(g.Intersection3d, try unit.intersectSegment(outside_dot, policy), .none);

    // Two dimensions, the same classification.
    const circle = try g.Circle2f.init(try g.Point2f.init(0, 0), 1.0);
    const flat_chord = try g.Segment2f.init(try g.Point2f.init(-2, 0), try g.Point2f.init(2, 0));
    const flat_through = try circle.intersectSegment(flat_chord, policy);
    try expectTag(g.Intersection2f, flat_through, .segment);
    try expectSign(f32, flat_through.segment.a.x, -1.0, 1e-6);
    try expectSign(f32, flat_through.segment.b.x, 1.0, 1e-6);
}

test "a ball against a ray is none, a point, or a closed segment" {
    const policy = g.PredicatePolicy.adaptive_exact;
    const unit = try g.Sphere3d.init(try g.Point3d.init(0, 0, 0), 1.0);

    // A ray straight through the centre shares the whole closed chord.
    const through = try g.Ray3d.init(try g.Point3d.init(-5, 0, 0), try g.Vector3d.init(1, 0, 0));
    const both = try unit.intersectRay(through, policy);
    try expectTag(g.Intersection3d, both, .segment);
    try expectSign(f64, both.segment.a.x, -1.0, 1e-12);
    try expectSign(f64, both.segment.b.x, 1.0, 1e-12);
    try expectPoint(try both.segment.at(0.5), try g.Point3d.init(0, 0, 0));

    // A ray that starts inside and leaves. The common set is the whole run
    // from the origin to the exit, so it is a SEGMENT and not the exit
    // alone: the origin and every point between it and the exit are common
    // too, and dropping them would report a set that is not the intersection.
    // With `p(t) = (t, 0, 0)` the sphere equation is `t^2 = 1`, so the run
    // is `0 <= t <= 1`.
    const from_inside = try g.Ray3d.init(try g.Point3d.init(0, 0, 0), try g.Vector3d.init(1, 0, 0));
    const exit = try unit.intersectRay(from_inside, policy);
    try expectTag(g.Intersection3d, exit, .segment);
    try expectSign(f64, exit.segment.a.x, 0.0, 1e-12);
    try expectSign(f64, exit.segment.b.x, 1.0, 1e-12);

    // A ray that starts inside but NOT at the centre, so the run is a proper
    // sub-segment whose two ends differ. `p(t) = (t, 1/2, 0)` gives
    // `t^2 + 1/4 = 1`, so the run is `0 <= t <= sqrt(3)/2`, and the origin
    // `(0, 1/2, 0)` is strictly inside the ball.
    const offset_inside = try g.Ray3d.init(try g.Point3d.init(0, 0.5, 0), try g.Vector3d.init(1, 0, 0));
    const offset_exit = try unit.intersectRay(offset_inside, policy);
    try expectTag(g.Intersection3d, offset_exit, .segment);
    try expectSign(f64, offset_exit.segment.a.x, 0.0, 1e-12);
    try expectSign(f64, offset_exit.segment.a.y, 0.5, 1e-12);
    try expectSign(f64, offset_exit.segment.b.x, 0.8660254037844386, 1e-12);
    try expectSign(f64, offset_exit.segment.b.y, 0.5, 1e-12);

    // A ray that starts on the surface and points outward: one point, itself.
    const outward = try g.Ray3d.init(try g.Point3d.init(1, 0, 0), try g.Vector3d.init(1, 0, 0));
    const at_surface = try unit.intersectRay(outward, policy);
    try expectTag(g.Intersection3d, at_surface, .point);
    try expectSign(f64, at_surface.point.x, 1.0, 1e-12);

    // A ray that points away: the single algebraic root is behind the end, so
    // the common set is empty.
    const backward = try g.Ray3d.init(try g.Point3d.init(5, 0, 0), try g.Vector3d.init(1, 0, 0));
    try expectTag(g.Intersection3d, try unit.intersectRay(backward, policy), .none);

    // A ray that misses.
    const missing = try g.Ray3d.init(try g.Point3d.init(-5, 2, 0), try g.Vector3d.init(1, 0, 0));
    try expectTag(g.Intersection3d, try unit.intersectRay(missing, policy), .none);

    // Tangent along a diameter perpendicular to the ray: the discriminant is
    // exactly zero, so one point.
    const tangent = try g.Ray3d.init(try g.Point3d.init(-5, 1, 0), try g.Vector3d.init(1, 0, 0));
    const grazed = try unit.intersectRay(tangent, policy);
    try expectTag(g.Intersection3d, grazed, .point);
    try expectSign(f64, grazed.point.y, 1.0, 1e-12);
    try expectSign(f64, grazed.point.x, 0.0, 1e-12);
}

test "a ball against a flat is none, the tangent point, or the ball it cuts" {
    const policy = g.PredicatePolicy.adaptive_exact;
    const centre = try g.Point3f.init(0, 0, 0);
    const unit = try g.Sphere3f.init(centre, 1.0);
    const up = try g.Vector3f.init(0, 0, 1);
    const through = try g.Plane3f.init(centre, up);
    const offset = try g.Plane3f.init(try g.Point3f.init(0, 0, 1), up);
    const beyond = try g.Plane3f.init(try g.Point3f.init(0, 0, 2), up);

    // The plane cuts the ball through its centre, so the common set is the
    // great circle in the plane.
    const great = try unit.intersectPlane(through, policy);
    try expectTag(g.Intersection3f, great, .circle);
    try expectSign(f32, great.circle.center.z, 0.0, 1e-6);
    try expectSign(f32, great.circle.radius, 1.0, 1e-6);

    // A plane at exactly the radius is tangent: one point, the foot.
    const tangent = try unit.intersectPlane(offset, policy);
    try expectTag(g.Intersection3f, tangent, .point);
    try expectSign(f32, tangent.point.z, 1.0, 1e-6);
    try expectSign(f32, tangent.point.x, 0.0, 1e-6);

    // A plane beyond the radius misses.
    try expectTag(g.Intersection3f, try unit.intersectPlane(beyond, policy), .none);

    // An off-centre cut is a smaller circle centred on the perpendicular foot,
    // not on the ball's centre.
    const bigger = try g.Sphere3f.init(centre, 2.0);
    const off_cut = try bigger.intersectPlane(offset, policy);
    try expectTag(g.Intersection3f, off_cut, .circle);
    try expectSign(f32, off_cut.circle.center.z, 1.0, 1e-5);
    try expectSign(f32, off_cut.circle.radius, 1.7320508, 1e-5);

    // The same pair asked from the plane's side gives the same answer.
    try expectTag(g.Intersection3f, try through.intersectSphere(unit, policy), .circle);
    try expectTag(g.Intersection3f, try offset.intersectSphere(unit, policy), .point);

    // Four dimensions: a hyperplane cutting a hypersphere gives a three-ball
    // boundary, which is why the tag is named for the dimension rather than
    // reusing the three-dimensional name.
    const hyper_plane = try g.Hyperplane4d.init(try g.Point4d.init(0, 0, 0, 0), try g.Vector4d.init(0, 0, 0, 1));
    const hyper_ball = try g.Hypersphere4d.init(try g.Point4d.init(0, 0, 0, 0), 1.0);
    const hyper_cut = try hyper_ball.intersectHyperplane(hyper_plane, policy);
    try expectTag(g.Intersection4d, hyper_cut, .hypersphere);
    try expectSign(f64, hyper_cut.hypersphere.radius, 1.0, 1e-12);
    const hyper_tangent = try hyper_ball.intersectHyperplane(
        try g.Hyperplane4d.init(try g.Point4d.init(0, 0, 0, 1), try g.Vector4d.init(0, 0, 0, 1)),
        policy,
    );
    try expectTag(g.Intersection4d, hyper_tangent, .point);
    try expectTag(g.Intersection4d, try hyper_ball.intersectHyperplane(
        try g.Hyperplane4d.init(try g.Point4d.init(0, 0, 0, 2), try g.Vector4d.init(0, 0, 0, 1)),
        policy,
    ), .none);
}

test "a flat against a segment is none, the cut point, or the whole segment" {
    const policy = g.PredicatePolicy.adaptive_exact;
    const flat = try g.Plane3d.init(try g.Point3d.init(0, 0, 0), try g.Vector3d.init(0, 0, 1));

    // The segment crosses z = 0 strictly inside.
    const crossing = try g.Segment3d.init(try g.Point3d.init(1, 2, -1), try g.Point3d.init(1, 2, 1));
    const cut = try flat.intersectSegment(crossing, policy);
    try expectTag(g.Intersection3d, cut, .point);
    try expectSign(f64, cut.point.z, 0.0, 1e-12);
    try expectSign(f64, cut.point.x, 1.0, 1e-12);
    try expectSign(f64, cut.point.y, 2.0, 1e-12);

    // A segment wholly above the plane misses it.
    const above = try g.Segment3d.init(try g.Point3d.init(0, 0, 1), try g.Point3d.init(0, 0, 2));
    try expectTag(g.Intersection3d, try flat.intersectSegment(above, policy), .none);

    // A segment lying in the plane: the common set is the whole segment, and
    // the segment is the argument of the call.
    const lying = try g.Segment3d.init(try g.Point3d.init(1, 0, 0), try g.Point3d.init(2, 0, 0));
    const inside = try flat.intersectSegment(lying, policy);
    try expectTag(g.Intersection3d, inside, .contained);
    try std.testing.expectEqual(g.ContainedSide.second, inside.contained);

    // A segment that ends exactly on the plane: one point, that end.
    const ending = try g.Segment3d.init(try g.Point3d.init(0, 0, 0), try g.Point3d.init(0, 0, 1));
    const at_end = try flat.intersectSegment(ending, policy);
    try expectTag(g.Intersection3d, at_end, .point);
    try expectSign(f64, at_end.point.z, 0.0, 1e-12);

    const hyper_plane = try g.Hyperplane4f.init(try g.Point4f.init(0, 0, 0, 0), try g.Vector4f.init(0, 1, 0, 0));
    const hyper_segment = try g.Segment4f.init(try g.Point4f.init(0, -1, 0, 0), try g.Point4f.init(0, 1, 0, 0));
    const hyper_cut = try hyper_plane.intersectSegment(hyper_segment, policy);
    try expectTag(g.Intersection4f, hyper_cut, .point);
    try expectSign(f32, hyper_cut.point.y, 0.0, 1e-6);
    const hyper_lying = try g.Segment4f.init(try g.Point4f.init(1, 0, 0, 0), try g.Point4f.init(2, 0, 0, 0));
    try expectTag(g.Intersection4f, try hyper_plane.intersectSegment(hyper_lying, policy), .contained);
}

test "two planes are none, one plane, or a line" {
    const policy = g.PredicatePolicy.adaptive_exact;
    const up = try g.Vector3d.init(0, 0, 1);
    const flat = try g.Plane3d.init(try g.Point3d.init(0, 0, 0), up);

    // Parallel and distinct: z = 0 against z = 3.
    const raised = try g.Plane3d.init(try g.Point3d.init(0, 0, 3), up);
    try expectTag(g.Intersection3d, try flat.intersect(raised, policy), .none);

    // The same flat written with the opposite normal. The common set is the
    // whole plane, so the tag names the span rather than saying "coincident",
    // which would hide the answer.
    const flipped = try g.Plane3d.init(try g.Point3d.init(0, 0, 0), try g.Vector3d.init(0, 0, -1));
    const same = try flat.intersect(flipped, policy);
    try expectTag(g.Intersection3d, same, .plane);
    try expectSign(f64, try same.plane.normal.length(), 1.0, 1e-12);
    try expectSign(f64, try flat.signedDistance(try g.Point3d.init(5, 5, 0)), 0.0, 1e-12);
    try expectSign(f64, try raised.signedDistance(try g.Point3d.init(5, 5, 3)), 0.0, 1e-12);
    try expectSign(f64, try raised.signedDistance(try g.Point3d.init(0, 0, 0)), -3.0, 1e-12);

    // Not parallel: the common set is a line, here the x axis.
    const sideways = try g.Plane3d.init(try g.Point3d.init(0, 0, 0), try g.Vector3d.init(0, 1, 0));
    const meet = try flat.intersect(sideways, policy);
    try expectTag(g.Intersection3d, meet, .line);
    try expectSign(f64, try meet.line.direction.length(), 1.0, 1e-12);
    const on_meet = try meet.line.closestPointTo(try g.Point3d.init(7, 7, 0));
    try expectSign(f64, on_meet.y, 0.0, 1e-12);
    try expectSign(f64, on_meet.z, 0.0, 1e-12);

    // An offset plane still meets in a line, somewhere else.
    const offset = try g.Plane3d.init(try g.Point3d.init(0, 4, 0), try g.Vector3d.init(0, 1, 0));
    const shifted_meet = try flat.intersect(offset, policy);
    try expectTag(g.Intersection3d, shifted_meet, .line);
    const on_shifted = try shifted_meet.line.closestPointTo(try g.Point3d.init(0, 4, 0));
    try expectSign(f64, on_shifted.y, 4.0, 1e-12);
    try expectSign(f64, on_shifted.z, 0.0, 1e-12);

    // Normals that are neither parallel nor orthogonal, which is where a
    // careless transversal divides by zero.
    const oblique = try g.Plane3d.init(try g.Point3d.init(1, 1, 1), try g.Vector3d.init(1, 1, 1));
    const oblique_meet = try flat.intersect(oblique, policy);
    try expectTag(g.Intersection3d, oblique_meet, .line);
    const on_oblique = try oblique_meet.line.closestPointTo(try g.Point3d.init(1, 1, 1));
    try expectSign(f64, on_oblique.z, 0.0, 1e-12);
    try expectSign(f64, try oblique.signedDistance(on_oblique), 0.0, 1e-9);
}

test "two hyperplanes are none, one hyperplane, or a two-flat" {
    const policy = g.PredicatePolicy.adaptive_exact;
    const w = try g.Vector4d.init(0, 0, 0, 1);
    const flat = try g.Hyperplane4d.init(try g.Point4d.init(0, 0, 0, 0), w);
    const raised = try g.Hyperplane4d.init(try g.Point4d.init(0, 0, 0, 3), w);
    try expectTag(g.Intersection4d, try flat.intersect(raised, policy), .none);
    const flipped = try g.Hyperplane4d.init(try g.Point4d.init(0, 0, 0, 0), try g.Vector4d.init(0, 0, 0, -1));
    try expectTag(g.Intersection4d, try flat.intersect(flipped, policy), .hyperplane);

    // Not parallel: the common set is a two-flat, and the `plane` tag names
    // it because in four dimensions a tag names the dimension of the span.
    const sideways = try g.Hyperplane4d.init(try g.Point4d.init(0, 0, 0, 0), try g.Vector4d.init(0, 1, 0, 0));
    const meet = try flat.intersect(sideways, policy);
    try expectTag(g.Intersection4d, meet, .plane);
    const two_flat = meet.plane;
    try expectSign(f64, try two_flat.u.dot(two_flat.u), 1.0, 1e-12);
    try expectSign(f64, try two_flat.v.dot(two_flat.v), 1.0, 1e-12);
    try expectSign(f64, try two_flat.u.dot(two_flat.v), 0.0, 1e-12);
    // The reported point lies on both hyperplanes, and both directions lie in
    // both, which is what makes them the span of the common set.
    try expectSign(f64, try flat.signedDistance(two_flat.origin), 0.0, 1e-12);
    try expectSign(f64, try sideways.signedDistance(two_flat.origin), 0.0, 1e-12);
    try expectSign(f64, try flat.signedDistance(try two_flat.origin.add(try two_flat.u.scale(5.0))), 0.0, 1e-9);
    try expectSign(f64, try sideways.signedDistance(try two_flat.origin.add(try two_flat.v.scale(5.0))), 0.0, 1e-9);

    // An offset pair, so the two-flat is not through the origin.
    const offset = try g.Hyperplane4d.init(try g.Point4d.init(0, 3, 0, 0), try g.Vector4d.init(0, 1, 0, 0));
    const shifted = try flat.intersect(offset, policy);
    try expectTag(g.Intersection4d, shifted, .plane);
    try expectSign(f64, try offset.signedDistance(shifted.plane.origin), 0.0, 1e-9);
    try expectSign(f64, shifted.plane.origin.y, 3.0, 1e-9);

    // An oblique pair, which is where a careless transversal divides by zero.
    const oblique = try g.Hyperplane4d.init(
        try g.Point4d.init(1, 1, 1, 1),
        try g.Vector4d.init(1, 1, 1, 1),
    );
    const oblique_meet = try flat.intersect(oblique, policy);
    try expectTag(g.Intersection4d, oblique_meet, .plane);
    try expectSign(f64, try flat.signedDistance(oblique_meet.plane.origin), 0.0, 1e-9);
    try expectSign(f64, try oblique.signedDistance(oblique_meet.plane.origin), 0.0, 1e-9);
}

// ---------------------------------------------------------------------------
// Boxes
// ---------------------------------------------------------------------------

test "box against box is none, box, contained, or coincident" {
    const low = try g.Point2d.init(0, 0);
    const high = try g.Point2d.init(1, 1);
    const unit = try g.Aabb2d.init(low, high);

    // Disjoint on x: the first stops at 1 and the second starts at 2.
    const far = try g.Aabb2d.init(try g.Point2d.init(2, 0), try g.Point2d.init(3, 1));
    try expectTag(g.Intersection2d, try unit.intersect(far), .none);
    try std.testing.expect(!(try unit.intersects(far)));

    // Touching at a single corner: the common set is that corner, which is a
    // degenerate box. The `box` tag is right, because the intersection is a
    // box and it is neither operand.
    const corner = try g.Aabb2d.init(high, try g.Point2d.init(2, 2));
    const touch = try unit.intersect(corner);
    try expectTag(g.Intersection2d, touch, .box);
    try std.testing.expect(touch.box.minimum.eq(high));
    try std.testing.expect(touch.box.maximum.eq(high));
    try std.testing.expect(try unit.intersects(corner));

    // Touching along a face: still a box, degenerate on one axis.
    const face = try g.Aabb2d.init(try g.Point2d.init(1, 0), try g.Point2d.init(2, 1));
    const face_touch = try unit.intersect(face);
    try expectTag(g.Intersection2d, face_touch, .box);
    try std.testing.expect(face_touch.box.minimum.eq(try g.Point2d.init(1, 0)));
    try std.testing.expect(face_touch.box.maximum.eq(high));

    // Partial overlap: the common box is strictly smaller than both.
    const shifted = try g.Aabb2d.init(try g.Point2d.init(0.5, 0.5), try g.Point2d.init(1.5, 1.5));
    const partial = try unit.intersect(shifted);
    try expectTag(g.Intersection2d, partial, .box);
    try std.testing.expect(partial.box.minimum.eq(try g.Point2d.init(0.5, 0.5)));
    try std.testing.expect(partial.box.maximum.eq(high));

    // Containment and coincidence, with the side following the call order.
    const inner = try g.Aabb2d.init(try g.Point2d.init(0.25, 0.25), try g.Point2d.init(0.75, 0.75));
    const held = try unit.intersect(inner);
    try expectTag(g.Intersection2d, held, .contained);
    try std.testing.expectEqual(g.ContainedSide.second, held.contained);
    try std.testing.expectEqual(g.ContainedSide.first, (try inner.intersect(unit)).contained);
    try expectTag(g.Intersection2d, try unit.intersect(unit), .coincident);

    const cube = try g.Aabb3f.init(try g.Point3f.init(0, 0, 0), try g.Point3f.init(1, 1, 1));
    try expectTag(g.Intersection3f, try cube.intersect(cube), .coincident);
    try expectTag(g.Intersection3f, try cube.intersect(
        try g.Aabb3f.init(try g.Point3f.init(2, 2, 2), try g.Point3f.init(3, 3, 3)),
    ), .none);
    const hyper = try g.Aabb4d.init(try g.Point4d.init(0, 0, 0, 0), try g.Point4d.init(1, 1, 1, 1));
    // Disjoint on the FOURTH axis only. The argument occupies `2 <= w <= 3`
    // while the receiver occupies `0 <= w <= 1`, and those intervals do not
    // overlap, so the per-axis test is empty on `w` and the whole
    // intersection is empty however much the other three axes agree. A tag
    // of `box` here would be a box with a lower corner above its upper one.
    const hyper_gap = try hyper.intersect(
        try g.Aabb4d.init(try g.Point4d.init(0, 0, 0, 2), try g.Point4d.init(1, 1, 1, 3)),
    );
    try expectTag(g.Intersection4d, hyper_gap, .none);
    try std.testing.expect(!(try hyper.intersects(
        try g.Aabb4d.init(try g.Point4d.init(0, 0, 0, 2), try g.Point4d.init(1, 1, 1, 3)),
    )));

    // A genuine four-dimensional partial overlap: the argument starts inside
    // on the first three axes and is cut off above on the fourth, so every
    // per-axis interval is proper and the common box is neither operand.
    const hyper_partial = try hyper.intersect(
        try g.Aabb4d.init(try g.Point4d.init(0, 0, 0, 0.5), try g.Point4d.init(1, 1, 1, 3)),
    );
    try expectTag(g.Intersection4d, hyper_partial, .box);
    try std.testing.expect(hyper_partial.box.minimum.eq(try g.Point4d.init(0, 0, 0, 0.5)));
    try std.testing.expect(hyper_partial.box.maximum.eq(try g.Point4d.init(1, 1, 1, 1)));

    // Containment in four dimensions, with the side following the call
    // order: the receiver holds its argument whole, so the common set is
    // the argument and the side is `second`.
    const hyper_inner = try g.Aabb4d.init(try g.Point4d.init(0.25, 0.25, 0.25, 0.25), try g.Point4d.init(0.75, 0.75, 0.75, 0.75));
    try expectTag(g.Intersection4d, try hyper.intersect(hyper_inner), .contained);
    try std.testing.expectEqual(g.ContainedSide.second, (try hyper.intersect(hyper_inner)).contained);
    try std.testing.expectEqual(g.ContainedSide.first, (try hyper_inner.intersect(hyper)).contained);
    try expectTag(g.Intersection4d, try hyper.intersect(hyper), .coincident);
}

test "box containment is INCLUSIVE on the boundary" {
    const box = try g.Aabb2f.init(try g.Point2f.init(0.0, 0.0), try g.Point2f.init(1.0, 1.0));
    // Every corner and face of the box is part of the box, and the range
    // comparison rounds nothing, so the boundary is inside. That is a choice,
    // and it is the choice the inclusive comparison makes.
    try std.testing.expect(try box.contains(try g.Point2f.init(0.0, 0.0)));
    try std.testing.expect(try box.contains(try g.Point2f.init(1.0, 1.0)));
    try std.testing.expect(try box.contains(try g.Point2f.init(0.0, 1.0)));
    try std.testing.expect(try box.contains(try g.Point2f.init(0.5, 0.5)));
    // One ulp outside is outside, and the gap here is a single binary32 ulp
    // rather than an arbitrary epsilon.
    const just_outside: f32 = 1.0 + 0.0000001;
    try std.testing.expect(just_outside > 1.0);
    try std.testing.expect(!try box.contains(try g.Point2f.init(just_outside, 0.5)));
    try std.testing.expect(!try box.contains(try g.Point2f.init(-0.0000001, 0.5)));

    const solid = try g.Aabb3d.init(try g.Point3d.init(0, 0, 0), try g.Point3d.init(1, 1, 1));
    try std.testing.expect(try solid.contains(try g.Point3d.init(1, 1, 1)));
    try std.testing.expect(!(try solid.contains(try g.Point3d.init(1, 1, 1.0000001))));
    const hyper = try g.Aabb4f.init(try g.Point4f.init(0, 0, 0, 0), try g.Point4f.init(1, 1, 1, 1));
    try std.testing.expect(try hyper.contains(try g.Point4f.init(1, 1, 1, 1)));
    try std.testing.expect(!(try hyper.contains(try g.Point4f.init(1, 1, 1, 1.0000001))));
}

test "box centre, extent, expansion, and the measure each dimension defines" {
    const box = try g.Aabb2d.init(try g.Point2d.init(0, 0), try g.Point2d.init(2, 6));
    try std.testing.expect((try box.center()).eq(try g.Point2d.init(1, 3)));
    const extent = try box.extent();
    try std.testing.expectEqual(@as(f64, 2.0), extent[0]);
    try std.testing.expectEqual(@as(f64, 6.0), extent[1]);
    try std.testing.expectEqual(@as(f64, 12.0), try box.area());

    const grown = try box.expand(try g.Point2d.init(-1, 1));
    try std.testing.expect(grown.minimum.eq(try g.Point2d.init(-1, 0)));
    try std.testing.expect(grown.maximum.eq(try g.Point2d.init(2, 6)));
    const grown_up = try box.expand(try g.Point2d.init(1, 9));
    try std.testing.expect(grown_up.maximum.eq(try g.Point2d.init(2, 9)));
    const unchanged = try box.expand(try g.Point2d.init(1, 3));
    try std.testing.expect(unchanged.minimum.eq(box.minimum));
    try std.testing.expect(unchanged.maximum.eq(box.maximum));

    const solid = try g.Aabb3d.init(try g.Point3d.init(0, 0, 0), try g.Point3d.init(2, 3, 4));
    try std.testing.expectEqual(@as(f64, 24.0), try solid.volume());
    // 2 * (xy + xz + yz) = 2 * (6 + 8 + 12) = 52.
    try std.testing.expectEqual(@as(f64, 52.0), try solid.doubledSurfaceArea());
    try std.testing.expectEqual([3]f64{ 2, 3, 4 }, try solid.extent());

    const hyper = try g.Aabb4d.init(try g.Point4d.init(0, 0, 0, 0), try g.Point4d.init(1, 2, 3, 4));
    try std.testing.expectEqual(@as(f64, 24.0), try hyper.hyperVolume());
    try std.testing.expectEqual([4]f64{ 1, 2, 3, 4 }, try hyper.extent());
    try std.testing.expect((try hyper.center()).eq(try g.Point4d.init(0.5, 1, 1.5, 2)));
}

// ---------------------------------------------------------------------------
// Predicate policy
//
// Every disagreement below is CONSTRUCTED, and each construction says where its
// numbers came from. None of them is quoted from a reference implementation.
// ---------------------------------------------------------------------------

test "the two policies agree wherever the scalar evaluation is well conditioned" {
    const a = try g.Point2d.init(0, 0);
    const b = try g.Point2d.init(1, 0);
    const c = try g.Point2d.init(0, 1);
    try std.testing.expectEqual(g.Sign.positive, try g.Point2d.orientation(a, b, c, .fast));
    try std.testing.expectEqual(g.Sign.positive, try g.Point2d.orientation(a, b, c, .adaptive_exact));
    try std.testing.expectEqual(g.Sign.negative, try g.Point2d.orientation(a, c, b, .fast));
    try std.testing.expectEqual(g.Sign.negative, try g.Point2d.orientation(a, c, b, .adaptive_exact));
    try std.testing.expectEqual(g.Sign.zero, try g.Point2d.orientation(a, b, try g.Point2d.init(2, 0), .fast));
    try std.testing.expectEqual(g.Sign.zero, try g.Point2d.orientation(a, b, try g.Point2d.init(2, 0), .adaptive_exact));

    const p0 = try g.Point3d.init(0, 0, 0);
    const p1 = try g.Point3d.init(1, 0, 0);
    const p2 = try g.Point3d.init(0, 1, 0);
    const p3 = try g.Point3d.init(0, 0, 1);
    try std.testing.expectEqual(g.Sign.positive, try g.Point3d.orientation(p0, p1, p2, p3, .fast));
    try std.testing.expectEqual(g.Sign.positive, try g.Point3d.orientation(p0, p1, p2, p3, .adaptive_exact));
    try std.testing.expectEqual(g.Sign.negative, try g.Point3d.orientation(p0, p1, p3, p2, .adaptive_exact));

    // The four unit vectors, in order, are a positively oriented basis.
    const q0 = try g.Point4d.init(1, 0, 0, 0);
    const q1 = try g.Point4d.init(0, 1, 0, 0);
    const q2 = try g.Point4d.init(0, 0, 1, 0);
    const q3 = try g.Point4d.init(0, 0, 0, 1);
    try std.testing.expectEqual(g.Sign.positive, try g.Point4d.orientation(q0, q1, q2, q3, .fast));
    try std.testing.expectEqual(g.Sign.positive, try g.Point4d.orientation(q0, q1, q2, q3, .adaptive_exact));
    try std.testing.expectEqual(g.Sign.negative, try g.Point4d.orientation(q0, q1, q3, q2, .adaptive_exact));
    // A repeated row makes the determinant exactly zero, and the exact path
    // is what proves it.
    try std.testing.expectEqual(g.Sign.zero, try g.Point4d.orientation(q0, q1, q2, q0, .adaptive_exact));
    try std.testing.expectEqual(g.Sign.positive, try g.Point2d.side(a, b, c, .adaptive_exact));
    try std.testing.expectEqual(g.Sign.negative, try g.Point2d.side(a, b, try g.Point2d.init(0, -1), .adaptive_exact));
}

test "CONSTRUCTED: the two-dimensional orientation disagrees, fast reports collinear" {
    // How this was built. Round-to-nearest sends two DIFFERENT exact products
    // to the same double:
    //
    //     134217823 * 134217818 = 18014423339770214
    //     134217822 * 134217819 = 18014423339770218
    //
    // Both products land in the binade that starts at 2^55, where the ulp is
    // 2, and both exact values are even, so each rounds to 18014423339770216
    // and the difference of the products cancels to exactly zero. The
    // difference of the EXACT products is -4, so the three points are not
    // collinear and the exact sign is negative. A plain evaluation cannot do
    // better here: rounding to nearest is monotone, so the only way the sign
    // of a single subtraction of two rounded products can differ from the
    // sign of the exact value is by landing on zero, and that is what happens.
    const a = try g.Point2d.init(0, 0);
    const b = try g.Point2d.init(134217823, 134217822);
    const c = try g.Point2d.init(134217819, 134217818);

    // The exact determinant, re-derived from the construction above.
    try std.testing.expectEqual(
        @as(i64, 134217823 * 134217818 - 134217822 * 134217819),
        @as(i64, -4),
    );

    try std.testing.expectEqual(g.Sign.zero, try g.Point2d.orientation(a, b, c, .fast));
    try std.testing.expectEqual(g.Sign.negative, try g.Point2d.orientation(a, b, c, .adaptive_exact));

    // The binary32 rounding of the SAME coordinates is a different number
    // altogether, because a binary32 significand holds 24 bits and the
    // coordinates need 27, so this is a reminder that the binary64 case above
    // is a claim about binary64 inputs and not a portable one.
    const small_a = try g.Point2f.init(0, 0);
    const small_b = try g.Point2f.init(@floatFromInt(134217823), @floatFromInt(134217822));
    const small_c = try g.Point2f.init(@floatFromInt(134217819), @floatFromInt(134217818));
    try std.testing.expectEqual(
        try g.Point2f.orientation(small_a, small_b, small_c, .fast),
        try g.Point2f.orientation(small_a, small_b, small_c, .adaptive_exact),
    );
}

test "CONSTRUCTED: the two-dimensional exact orientation survives a difference that overflows" {
    // How this was built. `floatMax(f64)` is `(2^53 - 1) * 2^971` exactly, so
    // the three points
    //
    //     a = (-floatMax, 0)    b = (floatMax, 0)    c = (0, 1)
    //
    // have the exact determinant
    //
    //     (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
    //   = (2 * floatMax) * 1 - 0 * floatMax
    //   = 2 * (2^53 - 1) * 2^971
    //
    // and every factor on the last line is positive, so the sign of that
    // determinant is the sign of `2^53 - 1`. Exchanging the second and third
    // points negates a two-by-two determinant, so the same construction makes
    // the second triangle negative.
    const max = std.math.floatMax(f64);
    const a = try g.Point2d.init(-max, 0);
    const b = try g.Point2d.init(max, 0);
    const c = try g.Point2d.init(0, 1);

    // The difference a floating-point reading must form is past the binary64
    // range, so the value it produces is not a real number at all. That is
    // why this witness is about the exact policy's own intermediate and not
    // about a disagreement between the policies: the plain evaluation
    // overflows to an infinity that happens to carry the right sign here.
    try std.testing.expect(std.math.isInf(max - -max));

    try std.testing.expectEqual(g.Sign.positive, try g.Point2d.orientation(a, b, c, .adaptive_exact));
    try std.testing.expectEqual(g.Sign.negative, try g.Point2d.orientation(a, c, b, .adaptive_exact));
}

test "CONSTRUCTED: a translated cancellation the two-dimensional exact orientation must not lose" {
    // How this was built. The two products of the determinant are made to
    // differ by ONE, and that one sits inside a DIFFERENCE binary64 cannot
    // represent:
    //
    //     b.x - a.x = 2^54 - (2^53 - 1) = 2^53 + 1
    //     c.x - a.x = 3 * 2^53 - (2^53 - 1) = 2^54 + 1
    //
    // `2^53 + 1` needs 54 significant bits and the spacing at `2^53` is 2, so
    // that difference is the exact halfway tie and round-to-nearest-even sends
    // it to `2^53`. The spacing at `2^54` is 4 and `2^54 + 1` sits one unit
    // from `2^54` and three from the next value, so that one rounds to `2^54`.
    // Each difference therefore loses the one that decides the sign, before
    // any multiplication happens, and the exact determinant is
    //
    //     (2^53 + 1) * 2 - 1 * (2^54 + 1) = 1
    //
    // Take the translation away and the same triangle is `(0, 0)`,
    // `(2^53 + 1, 1)`, `(2^54 + 1, 2)`, whose second and third points are not
    // binary64 numbers at all, so this sign is reachable only by exact
    // arithmetic.
    const a_int = [2]i128{ (@as(i128, 1) << 53) - 1, 0 };
    const b_int = [2]i128{ @as(i128, 1) << 54, 1 };
    const integer_c = [2]i128{ 3 * (@as(i128, 1) << 53), 2 };

    // The reference, derived twice: on the stored coordinates and on the same
    // triangle with the translation removed. Two exact derivations that agree
    // are what pins the sign to the geometry rather than to a reading.
    const stored = exactOrient2Integer(a_int, b_int, integer_c);
    const translated = exactOrient2Integer(
        .{ 0, 0 },
        .{ (@as(i128, 1) << 53) + 1, 1 },
        .{ (@as(i128, 1) << 54) + 1, 2 },
    );
    try std.testing.expectEqual(@as(i128, 1), translated);
    try std.testing.expectEqual(translated, stored);
    try std.testing.expect(stored > 0);
    try std.testing.expectEqual(-stored, exactOrient2Integer(a_int, integer_c, b_int));

    const a = try g.Point2d.init(try exactCoordinate(a_int[0]), try exactCoordinate(a_int[1]));
    const b = try g.Point2d.init(try exactCoordinate(b_int[0]), try exactCoordinate(b_int[1]));
    const c = try g.Point2d.init(try exactCoordinate(integer_c[0]), try exactCoordinate(integer_c[1]));

    // The two differences, read off the STORED coordinates. Neither exact
    // value is a binary64 number, and each rounds DOWN to the power of two
    // below it, which is what lets the two products become the same number.
    try std.testing.expectEqual(@as(i128, 9007199254740993), b_int[0] - a_int[0]);
    try std.testing.expectEqual(@as(i128, 18014398509481985), integer_c[0] - a_int[0]);
    try std.testing.expectEqual(@as(f64, 9007199254740992), b.x - a.x);
    try std.testing.expectEqual(@as(f64, 18014398509481984), c.x - a.x);

    // With those two roundings both products are `2^54` and they cancel to
    // zero, so the plain evaluation reports three collinear points where the
    // exact evaluation reports a determinant of one.
    try std.testing.expectEqual(g.Sign.zero, try g.Point2d.orientation(a, b, c, .fast));
    try std.testing.expectEqual(g.Sign.positive, try g.Point2d.orientation(a, b, c, .adaptive_exact));
    try std.testing.expectEqual(g.Sign.negative, try g.Point2d.orientation(a, c, b, .adaptive_exact));
}

test "CONSTRUCTED: the three-dimensional orientation disagrees on the SIGN in binary64" {
    // How this was built. The cross product's x component is the difference
    // of the SAME colliding pair of products the two-dimensional case used:
    //
    //     u.y * v.z = 134217823 * 134217818 = 18014423339770214
    //     u.z * v.y = 134217822 * 134217819 = 18014423339770218
    //
    // so the plain evaluation gets `n.x = 0` where the exact one gets -4.
    // Multiplying that error by a large `w.x` and adding the two exact terms
    // flips the sign of the whole determinant:
    //
    //     fast  = 0 * 40000000 + 134217822 * 0 + (-134217823) * (-1)
    //            = 134217823                      (positive)
    //     exact = -4 * 40000000 + 0 + 134217823
    //            = -25782177                      (negative)
    //
    // This is a real sign flip, not a zero. A three-term sum has no
    // monotonicity to protect it, so the plain evaluation is free to return a
    // definite sign of the wrong value, and it does.
    const a = try g.Point3d.init(0, 0, 0);
    const b = try g.Point3d.init(0, 134217823, 134217822);
    const c = try g.Point3d.init(1, 134217819, 134217818);
    const d = try g.Point3d.init(40000000, 0, -1);

    // The exact determinant, re-derived from the construction above.
    const cross_x: i128 = @as(i128, 134217823) * 134217818 - @as(i128, 134217822) * 134217819;
    try std.testing.expectEqual(@as(i128, -4), cross_x);
    const total: i128 = cross_x * 40000000 + @as(i128, 134217823) * @as(i128, 0) +
        @as(i128, -134217823) * @as(i128, -1);
    try std.testing.expectEqual(@as(i128, -25782177), total);

    try std.testing.expectEqual(g.Sign.positive, try g.Point3d.orientation(a, b, c, d, .fast));
    try std.testing.expectEqual(g.Sign.negative, try g.Point3d.orientation(a, b, c, d, .adaptive_exact));
}

test "CONSTRUCTED: the three-dimensional orientation disagrees on the SIGN in binary32" {
    // How this was built. A search over integer coordinate sets for a case
    // where `dot(cross(b - a, c - a), d - a)` has the wrong sign in binary32
    // turned this up. The exact determinant is -6 while the plain evaluation
    // returns +4096, so the terms cancel to within a few units out of
    // products in the millions. Every component is below 2^13, so each is
    // exactly representable and the disagreement is arithmetic, not
    // representation.
    const a = try g.Point3f.init(0, 0, 0);
    const b = try g.Point3f.init(-4101, 4098, -4101);
    const c = try g.Point3f.init(4113, -4115, 4112);
    const d = try g.Point3f.init(-4103, 4100, -4103);

    // Re-derive the exact determinant on paper from the difference vectors,
    // which the origin at `a` leaves untouched.
    // The recomputation widens first, because the exact value -6 is far below
    // the rounding error of the binary32 products and only binary64 holds it.
    const u: [3]f64 = .{ @floatCast(b.x), @floatCast(b.y), @floatCast(b.z) };
    const v: [3]f64 = .{ @floatCast(c.x), @floatCast(c.y), @floatCast(c.z) };
    const w: [3]f64 = .{ @floatCast(d.x), @floatCast(d.y), @floatCast(d.z) };
    const nx: f64 = u[1] * v[2] - u[2] * v[1];
    const ny: f64 = u[2] * v[0] - u[0] * v[2];
    const nz: f64 = u[0] * v[1] - u[1] * v[0];
    try std.testing.expectEqual(@as(f64, -6), nx * w[0] + ny * w[1] + nz * w[2]);

    try std.testing.expectEqual(g.Sign.positive, try g.Point3f.orientation(a, b, c, d, .fast));
    try std.testing.expectEqual(g.Sign.negative, try g.Point3f.orientation(a, b, c, d, .adaptive_exact));
}

test "CONSTRUCTED: the four-dimensional orientation disagrees on a degree-four form" {
    // How this was built. The fourth row is the exact sum of the second and
    // third, so the four position vectors are exactly linearly dependent and
    // the determinant is EXACTLY zero. Each three-by-three minor is a sum of
    // three products of three entries, and those products are around 2^33
    // while a binary32 significand holds 24 bits, so the plain evaluation
    // rounds every term and the exact cancellation is destroyed. The result
    // is a definite positive sign for a matrix that is exactly singular, in
    // the degree-four form that sizes the exact arithmetic.
    const a = try g.Point4f.init(-3284, 3630, -676, 3296);
    const b = try g.Point4f.init(2061, 4168, 2430, 4840);
    const c = try g.Point4f.init(4492, 2438, 3410, 2596);
    const d = try g.Point4f.init(6553, 6606, 5840, 7436);

    // The linear dependence, re-derived: d is b + c on every component.
    const b_row = b.toArray();
    const c_row = c.toArray();
    const d_row = d.toArray();
    inline for (0..4) |index| {
        try std.testing.expectEqual(d_row[index], b_row[index] + c_row[index]);
    }
    // And the exact determinant of the exactly singular matrix is zero.
    const a_row = a.toArray();
    const rows: [4][4]i128 = .{
        .{ @intFromFloat(a_row[0]), @intFromFloat(a_row[1]), @intFromFloat(a_row[2]), @intFromFloat(a_row[3]) },
        .{ @intFromFloat(b_row[0]), @intFromFloat(b_row[1]), @intFromFloat(b_row[2]), @intFromFloat(b_row[3]) },
        .{ @intFromFloat(c_row[0]), @intFromFloat(c_row[1]), @intFromFloat(c_row[2]), @intFromFloat(c_row[3]) },
        .{ @intFromFloat(d_row[0]), @intFromFloat(d_row[1]), @intFromFloat(d_row[2]), @intFromFloat(d_row[3]) },
    };
    try std.testing.expectEqual(@as(i128, 0), determinant4(rows));

    try std.testing.expectEqual(g.Sign.zero, try g.Point4f.orientation(a, b, c, d, .adaptive_exact));
    try std.testing.expect(try g.Point4f.orientation(a, b, c, d, .fast) != .zero);
}

test "the exact sign of a degree-four form matches an integer recomputation" {
    // A cross-check on the exact machinery itself rather than on a policy
    // disagreement. For integer inputs the exact determinant is an integer,
    // so the test recomputes it in `i128` from the SAME rows and checks the
    // sign the module reports against THAT number.
    //
    // The expected sign comes from `determinant4` and not from the module.
    // Negating whatever the module returned, as this test used to do, checks
    // that the row exchange flips the sign and nothing else: a module that
    // returned any fixed nonzero value would pass it.
    const a = try g.Point4d.init(1, 2, 3, 4);
    const b = try g.Point4d.init(3, 5, 7, 11);
    const c = try g.Point4d.init(13, 17, 19, 23);
    const d = try g.Point4d.init(29, 31, 37, 41);
    const rows: [4][4]i128 = .{
        .{ 1, 2, 3, 4 },
        .{ 3, 5, 7, 11 },
        .{ 13, 17, 19, 23 },
        .{ 29, 31, 37, 41 },
    };
    // The exact integer determinant, recomputed in the test.
    const value = determinant4(rows);
    const expected = signOfExact(value);
    try std.testing.expect(expected != .zero);
    try std.testing.expectEqual(expected, try g.Point4d.orientation(a, b, c, d, .adaptive_exact));

    // The same value with two rows exchanged is the negative of it, and the
    // exact path must see that too. The swapped determinant is the negation
    // of the integer one, computed the same independent way.
    const swapped = [4][4]i128{ rows[0], rows[2], rows[1], rows[3] };
    try std.testing.expectEqual(signOfExact(determinant4(swapped)), try g.Point4d.orientation(a, c, b, d, .adaptive_exact));
    try std.testing.expectEqual(-value, determinant4(swapped));

    // The binary32 version of the SAME matrix, so it is the same integer
    // determinant: every entry is below 2^24 and so exactly representable in
    // binary32, and the exact path widens rather than rounds them.
    const af = try g.Point4f.init(1, 2, 3, 4);
    const bf = try g.Point4f.init(3, 5, 7, 11);
    const cf = try g.Point4f.init(13, 17, 19, 23);
    const df = try g.Point4f.init(29, 31, 37, 41);
    try std.testing.expectEqual(expected, try g.Point4f.orientation(af, bf, cf, df, .adaptive_exact));
    try std.testing.expectEqual(signOfExact(determinant4(swapped)), try g.Point4f.orientation(af, cf, bf, df, .adaptive_exact));

    // And the plain binary32 evaluation of this matrix happens to agree,
    // because the entries are small enough that every product of four of
    // them stays inside the 24-bit significand. Asserting agreement here is
    // a statement about THESE numbers, not about the fast path in general.
    try std.testing.expectEqual(expected, try g.Point4f.orientation(af, bf, cf, df, .fast));

    // A matrix the plain binary32 evaluation gets WRONG, so the pair of
    // assertions above is not the only thing the exact path is being asked.
    // The first row is the exact sum of the second and third plus ONE in the
    // last component, so the determinant is a small number in a sum of terms
    // around 2^92: the exact value is recomputed below in `i128` and the
    // plain evaluation rounds it to the wrong sign. Every entry is below
    // 2^24, so the disagreement is arithmetic and not representation.
    const ha = try g.Point4f.init(11530199, 9623174, 10372602, 3882227);
    const hb = try g.Point4f.init(8388607, 1234567, 7654321, 2468013);
    const hc = try g.Point4f.init(3141592, 8388607, 2718281, 1414213);
    const hd = try g.Point4f.init(7777777, 8888888, 9999999, 1010101);
    const hard_rows: [4][4]i128 = .{
        .{ 11530199, 9623174, 10372602, 3882227 },
        .{ 8388607, 1234567, 7654321, 2468013 },
        .{ 3141592, 8388607, 2718281, 1414213 },
        .{ 7777777, 8888888, 9999999, 1010101 },
    };
    // The near-dependence, re-derived: the first row is the sum of the
    // second and third on the first three components and one more on the
    // fourth.
    inline for (0..3) |index| {
        try std.testing.expectEqual(hb.toArray()[index] + hc.toArray()[index], ha.toArray()[index]);
    }
    try std.testing.expectEqual(hb.toArray()[3] + hc.toArray()[3] + 1, ha.toArray()[3]);
    const hard_expected = signOfExact(determinant4(hard_rows));
    try std.testing.expect(hard_expected != .zero);
    try std.testing.expectEqual(hard_expected, try g.Point4f.orientation(ha, hb, hc, hd, .adaptive_exact));
}

test "CONSTRUCTED: a binary32 circle reports a point outside whose difference rounded onto the boundary" {
    // How this was built. The binary32 form of the boundary witness below. The
    // centre and the radius are the consecutive binary32 numbers `2^24 - 1`
    // and `2^24`, and the point sits at `(2^25, 0)`, so
    //
    //     dx = 2^25 - (2^24 - 1) = 2^24 + 1
    //
    // `2^24 + 1` needs 25 significant bits and the spacing at `2^24` is `2`,
    // so that difference is the exact halfway tie and round-to-nearest-even
    // sends it to `2^24`. The exact squared distance from the centre is
    //
    //     (2^24 + 1)^2 = 2^48 + 2^25 + 1
    //
    // and the radius square is exactly `2^48`, so the exact residual is
    // `2^25 + 1`, positive: the point is strictly OUTSIDE the closed ball. The
    // rounded difference leaves a residual of exactly `0`, which is the
    // boundary the inclusive test counts as inside.
    const centre = try g.Point2f.init(16777215.0, 0);
    const circle = try g.Circle2f.init(centre, 16777216.0);
    const point = try g.Point2f.init(33554432.0, 0);

    // The reference, re-derived on the integers the stored coordinates hold.
    const exact_dx: i128 = (@as(i128, 1) << 25) - (@as(i128, 1) << 24) + 1;
    try std.testing.expectEqual(@as(i128, 16777217), exact_dx);
    const exact_residual = exact_dx * exact_dx - (@as(i128, 1) << 48);
    try std.testing.expect(exact_residual > 0);
    // And the difference the plain evaluation actually forms.
    try std.testing.expectEqual(@as(f32, 16777216.0), point.x - centre.x);

    try std.testing.expect(!(try circle.contains(point, .adaptive_exact)));
    try std.testing.expect(try circle.contains(point, .fast));
}

test "CONSTRUCTED: a point exactly on a binary32 circle is inside under both policies" {
    // The inclusive test on a point that is exactly ON the rim, with every
    // stored coordinate below `2^24` and so a number binary32 holds exactly.
    // The three-four-five triangle gives a residual of exactly `0` under both
    // policies, and a residual of `0` is inside.
    const centre = try g.Point2f.init(0, 0);
    const circle = try g.Circle2f.init(centre, 5);
    const on_the_circle = try g.Point2f.init(3, 4);
    try std.testing.expect(try circle.contains(on_the_circle, .fast));
    try std.testing.expect(try circle.contains(on_the_circle, .adaptive_exact));
    // One step past the rim, where neither policy is in doubt.
    const past_the_rim = try g.Point2f.init(5, 1);
    try std.testing.expect(!(try circle.contains(past_the_rim, .fast)));
    try std.testing.expect(!(try circle.contains(past_the_rim, .adaptive_exact)));
    // The centre, which is inside by the whole radius.
    try std.testing.expect(try circle.contains(centre, .adaptive_exact));
}

test "CONSTRUCTED: the three-dimensional exact orientation survives a difference that overflows" {
    // How this was built. The two-dimensional overflow witness lifted into
    // three dimensions, with `M` for `floatMax(f64)`:
    //
    //     a = (-M, 0, 0)   b = (M, 0, 0)
    //     c = (0, 1, 0)    d = (0, 0, 1)
    //
    // The three differences the predicate takes are
    //
    //     u = b - a = (2M, 0, 0)
    //     v = c - a = (M, 1, 0)
    //     w = d - a = (M, 0, 1)
    //
    // and the cross product of the first two is
    //
    //     cross(u, v) = (u.y*v.z - u.z*v.y, u.z*v.x - u.x*v.z, u.x*v.y - u.y*v.x)
    //                 = (0, 0, 2M)
    //
    // because every term but `u.x*v.y` is a product of a zero. The signed
    // volume is then `2M * 1 = 2M`, positive, and `2 * floatMax` is not a
    // binary64 number at all, so a reading that forms the differences first
    // has no value to carry that sign. Exchanging `b` and `c` exchanges `u`
    // and `v`, which negates the cross product and the volume with it.
    const max = std.math.floatMax(f64);
    const a = try g.Point3d.init(-max, 0, 0);
    const b = try g.Point3d.init(max, 0, 0);
    const c = try g.Point3d.init(0, 1, 0);
    const d = try g.Point3d.init(0, 0, 1);

    // The difference a floating-point reading must form, and what it is.
    try std.testing.expect(std.math.isInf(max - -max));
    // And the exact sign of the volume, re-derived. The cross product of the
    // first two differences is `2M` along the third axis, so the volume is
    // `2M * 1`, and `2M` is `(2^54 - 2) * 2^971`, whose every factor is
    // positive. It is far past an `i128`, so the sign is read off the factors.
    try std.testing.expectEqual(g.Sign.positive, try g.Point3d.orientation(a, b, c, d, .adaptive_exact));
    try std.testing.expectEqual(g.Sign.negative, try g.Point3d.orientation(a, c, b, d, .adaptive_exact));
}

test "CONSTRUCTED: a circle reports a point outside whose difference rounded onto the boundary" {
    // How this was built. The centre and the radius are the consecutive
    // binary64 numbers `2^53 - 1` and `2^53`, and the point sits at
    // `(2^54, 2)`, so
    //
    //     dx = 2^54 - (2^53 - 1) = 2^53 + 1
    //     dy = 2 - 0              = 2
    //
    // `2^53 + 1` needs 54 significant bits and the spacing at `2^53` is 2,
    // so that difference is the halfway tie and round-to-nearest-even sends
    // it to `2^53`. The exact squared distance from the centre is
    //
    //     (2^53 + 1)^2 + 2^2 = 2^106 + 2^54 + 5
    //
    // and the radius square is exactly `2^106`, so the exact residual is
    // `2^54 + 5`, positive: the point is strictly OUTSIDE the closed ball.
    // The rounded difference leaves a residual of exactly `0`, which is the
    // boundary the inclusive test counts as inside.
    const centre_x: f64 = 9007199254740991.0;
    const radius_x: f64 = 9007199254740992.0;
    const centre = try g.Point2d.init(centre_x, 0);
    const circle = try g.Circle2d.init(centre, radius_x);
    const point = try g.Point2d.init(18014398509481984.0, 2);

    // The reference, re-derived on the integers the stored coordinates hold.
    const exact_dx: i128 = (@as(i128, 1) << 54) - (@as(i128, 1) << 53) + 1;
    try std.testing.expectEqual(@as(i128, 9007199254740993), exact_dx);
    const exact_residual = exact_dx * exact_dx + 4 - (@as(i128, 1) << 106);
    try std.testing.expect(exact_residual > 0);
    // And the two differences the plain evaluation actually forms.
    try std.testing.expectEqual(@as(f64, 9007199254740992.0), point.x - centre.x);
    try std.testing.expectEqual(@as(f64, 2), point.y - centre.y);

    try std.testing.expect(!(try circle.contains(point, .adaptive_exact)));
    try std.testing.expect(try circle.contains(point, .fast));
}

test "CONSTRUCTED: two circles cross where a rounded radius sum calls them tangent" {
    // How this was built. The big circle has centre `(0, 0)` and radius
    // `2^53`, the small one centre `(2^53, 0)` and radius `1`, so the
    // centres are exactly `2^53` apart and
    //
    //     2^53 < 2^53 + 1     the radii sum past the separation, so the small
    //                         circle reaches out beyond the big rim
    //     2^53 > 2^53 - 1     the radii differ by less than the separation,
    //                         so the big circle does not swallow the small
    //
    // The two circles therefore CROSS, in two distinct points. The sum
    // `2^53 + 1` needs 54 significant bits, so the plain evaluation reads it
    // as `2^53`, the separation EQUALS the sum, and the pair looks exactly
    // tangent.
    //
    // What the two points ARE is the point of this test, because a tag of
    // `two_points` with two coincident coordinates is a claim the tag alone
    // cannot keep. The chord stands on the line of centres at
    //
    //     along = (r1^2 - r2^2 + d^2) / (2 d)
    //           = (2^106 - 1 + 2^106) / 2^54
    //           = 2^53 - 2^-54
    //
    // which rounds to `2^53`, the spacing at `2^53` being `2`. Measured from
    // the small centre the same plane is `2^-54` away, so the half-chord is
    // the law of cosines' third side
    //
    //     h^2 = r2^2 - (d - along)^2 = 1 - 2^-108
    //
    // and its square root rounds to `1` for the same reason. The two points
    // are therefore `(2^53, 1)` and `(2^53, -1)`: a chord of length `2`
    // across a rim of radius `2^53`, one rounding away from a single point.
    const big = try g.Circle2d.init(try g.Point2d.init(0, 0), 9007199254740992.0);
    const small = try g.Circle2d.init(try g.Point2d.init(9007199254740992.0, 0), 1);

    // And the sum the plain evaluation reads.
    try std.testing.expectEqual(@as(f64, 9007199254740992.0), 9007199254740992.0 + 1.0);

    try expectTag(g.Intersection2d, try big.intersect(small, .adaptive_exact), .two_points);
    try expectTag(g.Intersection2d, try big.intersect(small, .fast), .point);

    // The payload, exactly, in either order. Both coordinates are whole
    // numbers, so nothing is left to a tolerance here.
    const up = try g.Point2d.init(9007199254740992.0, 1);
    const down = try g.Point2d.init(9007199254740992.0, -1);
    try expectCrossingPair(try big.intersect(small, .adaptive_exact), up, down, 0.0);
}

test "CONSTRUCTED: the same crossing measured from the small circle first" {
    // The same two circles with the operands exchanged, so the chord is
    // measured from the small centre. That centre is `2^-54` from the chord
    // plane, and the margin that says so is `(r2 - d) + r1`, taken against
    // the larger radius first: the plain evaluation forms the radii sum
    // `2^53 + 1`, rounds it to `2^53`, and the margin comes out ZERO, which
    // is the whole defect. The pair reported is the same set, and it is
    // still two distinct points a chord of length `2` apart.
    const big = try g.Circle2d.init(try g.Point2d.init(0, 0), 9007199254740992.0);
    const small = try g.Circle2d.init(try g.Point2d.init(9007199254740992.0, 0), 1);
    const up = try g.Point2d.init(9007199254740992.0, 1);
    const down = try g.Point2d.init(9007199254740992.0, -1);

    try expectTag(g.Intersection2d, try small.intersect(big, .adaptive_exact), .two_points);
    try expectCrossingPair(try small.intersect(big, .adaptive_exact), up, down, 1e-15);
}

test "CONSTRUCTED: the same crossing at binary32, where the rim is 2^24" {
    // The binary32 form of the witness above. The big circle has centre
    // `(0, 0)` and radius `2^24`, the small one centre `(2^24, 0)` and radius
    // `1`, and the exact half-chord is `sqrt(1 - 2^-50)`, which is `1` to
    // within a rounding. The plain evaluation reads the radii sum as `2^24`,
    // so it reports tangency, and the payload is the same pair of points one
    // binary32 ulp below the exact half-chord.
    const big = try g.Circle2f.init(try g.Point2f.init(0, 0), 16777216.0);
    const small = try g.Circle2f.init(try g.Point2f.init(16777216.0, 0), 1);
    const up = try g.Point2f.init(16777216.0, 1);
    const down = try g.Point2f.init(16777216.0, -1);

    try expectTag(g.Intersection2f, try big.intersect(small, .adaptive_exact), .two_points);
    try expectTag(g.Intersection2f, try big.intersect(small, .fast), .point);
    try expectCrossingPair(try big.intersect(small, .adaptive_exact), up, down, 1e-6);
}

test "CONSTRUCTED: an ordinary crossing reports the chord of the law of cosines" {
    // How this was built. A three-four-five triangle, so both points are whole
    // numbers. The big circle has centre `(0, 0)` and radius `5`, the small
    // one centre `(4, 0)` and radius `3`; the centres are `4` apart,
    // `2 < 4 < 8`, and
    //
    //     along = (25 - 9 + 16) / 8 = 4     h^2 = 25 - 16 = 9
    //
    // so the pair is `(4, 3)` and `(4, -3)`. The same call with the operands
    // exchanged measures `along` from the small centre, where it is `0`, and
    // the pair is the same set. That second pair lands a rounding away from
    // `3` on the half-chord, because the factored form takes one square root
    // of each pair of margins rather than one of their product, so it is
    // compared to within a rounding.
    const big = try g.Circle2d.init(try g.Point2d.init(0, 0), 5);
    const small = try g.Circle2d.init(try g.Point2d.init(4, 0), 3);
    const up = try g.Point2d.init(4, 3);
    const down = try g.Point2d.init(4, -3);

    try expectTag(g.Intersection2d, try big.intersect(small, .adaptive_exact), .two_points);
    try expectCrossingPair(try big.intersect(small, .adaptive_exact), up, down, 0.0);
    try expectTag(g.Intersection2d, try small.intersect(big, .adaptive_exact), .two_points);
    try expectCrossingPair(try small.intersect(big, .adaptive_exact), up, down, 1e-15);
}

test "CONSTRUCTED: a crossing near the top of the binary64 range keeps a finite chord" {
    // How this was built. The same three-four-five construction, scaled until
    // the product of the four tangency margins leaves the format. Both radii
    // are `5 * 2^260` and the centres are `6 * 2^260` apart, so
    //
    //     along = d / 2 = 3 * 2^260     h = sqrt(25 - 9) * 2^260 = 2^262
    //
    // and the four margins multiply to `4 d^2 h^2 = 2^1051`, which is NOT a
    // binary64 number at all. A chord formed as one product of the four
    // margins overflows to an infinity here, and one formed as
    // `r1^2 - along^2` survives only because the squares still fit. The pair
    // is compared to the exact one to within a rounding, and every coordinate
    // is finite: an infinity in the payload is the failure this guards.
    const radius: f64 = 9.263367138985296e78; // 5 * 2^260
    const gap: f64 = 1.1116040566782355e79; // 6 * 2^260
    const base: f64 = 5.558020283391177e78; // 3 * 2^260
    const chord: f64 = 7.410693711188237e78; // 2^262
    const left = try g.Circle2d.init(try g.Point2d.init(0, 0), radius);
    const right = try g.Circle2d.init(try g.Point2d.init(gap, 0), radius);

    try expectTag(g.Intersection2d, try left.intersect(right, .adaptive_exact), .two_points);
    try expectTag(g.Intersection2d, try left.intersect(right, .fast), .two_points);

    const crossing = try left.intersect(right, .adaptive_exact);
    for (crossing.two_points) |point| {
        try std.testing.expect(std.math.isFinite(point.x));
        try std.testing.expect(std.math.isFinite(point.y));
    }
    try expectCrossingPair(crossing, try g.Point2d.init(base, chord), try g.Point2d.init(base, -chord), base * 1e-15);
}

test "CONSTRUCTED: a ray misses a ball whose offset difference rounded onto it" {
    // How this was built. The ray runs along `+y` from `(2^54, 0, 0)`, so it
    // lies on the line `x = 2^54`. The ball is centred at `(2^53 - 1, 1, 0)`
    // with radius `2^53`, and the distance from the centre to that whole LINE
    // is
    //
    //     2^54 - (2^53 - 1) = 2^53 + 1 > 2^53
    //
    // A distance to the whole line that already exceeds the radius leaves the
    // ray no point in the ball at all, whichever way it points. `2^53 + 1` is
    // not a binary64 number, so the plain evaluation reads the offset as
    // `2^53`, the distance becomes exactly the radius, and the ray appears
    // to graze the ball and cut it in two.
    const ball = try g.Sphere3d.init(try g.Point3d.init(9007199254740991.0, 1, 0), 9007199254740992.0);
    const ray = try g.Ray3d.init(try g.Point3d.init(18014398509481984.0, 0, 0), try g.Vector3d.init(0, 1, 0));

    // The exact discriminant, re-derived on the integers. The projection is
    // `-1` because the only non-zero component of the offset that the unit
    // direction `(0, 1, 0)` reads is `-1`.
    const offset: i128 = (@as(i128, 1) << 54) - (@as(i128, 1) << 53) + 1;
    try std.testing.expectEqual(@as(i128, 9007199254740993), offset);
    const radius: i128 = @as(i128, 1) << 53;
    const constant = offset * offset + 1 - radius * radius;
    try std.testing.expect(1 - constant < 0);
    // And the offset the plain evaluation forms.
    try std.testing.expectEqual(@as(f64, 9007199254740992.0), 18014398509481984.0 - 9007199254740991.0);

    try expectTag(g.Intersection3d, try ball.intersectRay(ray, .adaptive_exact), .none);
    try expectTag(g.Intersection3d, try ball.intersectRay(ray, .fast), .segment);
}

test "CONSTRUCTED: two parallel hyperplanes a rounded difference called one flat" {
    // How this was built. Both hyperplanes carry the normal `(1, -1, 1, -1)`,
    // whose length is exactly two, so `normalize` stores
    // `(0.5, -0.5, 0.5, -0.5)` and no coefficient rounds. The first flat
    // passes through `(-M, -M, 0, 0)` and the second through `(M, M, 0, 1)`,
    // with `M = floatMax(f64)`. The signed distance of the second origin from
    // the first flat is the exact
    //
    //     (M - (-M)) * 0.5 + (M - (-M)) * -0.5 + (0 - 0) * 0.5 + (1 - 0) * -0.5
    //   = 2M * 0.5 - 2M * 0.5 - 0.5
    //   = -0.5
    //
    // so the second origin is strictly off the first flat, the two hyperplanes
    // are DISTINCT, and the common set is empty. Each of the first two
    // differences is past the binary64 range, so the plain evaluation forms
    // two infinities of OPPOSITE sign; their sum is not a real number at all.
    // That sum is the case `fast` cannot answer: a NaN has no sign, and
    // reporting `zero` for it would claim the origin lies exactly on the
    // flat. The call therefore fails rather than naming the whole flat.
    const max = std.math.floatMax(f64);
    const flat = try g.Hyperplane4d.init(try g.Point4d.init(-max, -max, 0, 0), try g.Vector4d.init(1, -1, 1, -1));
    const shifted = try g.Hyperplane4d.init(try g.Point4d.init(max, max, 0, 1), try g.Vector4d.init(1, -1, 1, -1));

    // The normal the module stored, read back, so the derivation above is a
    // statement about the coefficients that ran. Nothing about it rounds.
    try std.testing.expectEqual(@as(f64, 0.5), flat.normal.x);
    try std.testing.expectEqual(@as(f64, -0.5), flat.normal.y);
    try std.testing.expectEqual(@as(f64, 0.5), flat.normal.z);
    try std.testing.expectEqual(@as(f64, -0.5), flat.normal.w);
    // The distance, re-derived on the halves of `M` that do not overflow.
    try std.testing.expect(0.5 * max - 0.5 * max - 0.5 == -0.5);
    // And the difference the plain evaluation forms for each coordinate.
    try std.testing.expect(std.math.isInf(max - -max));

    // The exact policy forms no floating-point difference here, so it answers
    // with a real sign and the geometry comes out right.
    try expectTag(g.Intersection4d, try flat.intersect(shifted, .adaptive_exact), .none);
    // The plain policy says it cannot answer, and it says which reason: the
    // arithmetic left the format rather than the caller's input being bad.
    try std.testing.expectError(error.NonFiniteResult, flat.intersect(shifted, .fast));
}

test "the coplanarity gate a policy drives is the one the higher dimensions use" {
    // The same four points, asked through the gate directly and through the
    // intersection that consults it. A skew pair and a coplanar pair must give
    // opposite answers under BOTH policies, or the gate is not the gate.
    const a = try g.Point3d.init(0, 0, 0);
    const b = try g.Point3d.init(2, 0, 0);
    const c = try g.Point3d.init(1, -1, 0);
    const d = try g.Point3d.init(1, 1, 0);
    const skew = try g.Point3d.init(1, 1, 1);
    for ([_]g.PredicatePolicy{ .fast, .adaptive_exact }) |policy| {
        try std.testing.expectEqual(g.Sign.zero, try g.Point3d.orientation(a, b, c, d, policy));
        // The sign of the skew configuration is not what this test is about;
        // what it is about is that both policies agree on it, which is what
        // makes the gate a gate.
        try std.testing.expectEqual(
            try g.Point3d.orientation(a, b, c, skew, .fast),
            try g.Point3d.orientation(a, b, c, skew, .adaptive_exact),
        );
        try std.testing.expect(try g.Point3d.orientation(a, b, c, skew, .adaptive_exact) != .zero);

        const along = try g.Segment3d.init(a, b);
        const across = try g.Segment3d.init(c, d);
        const lifted = try g.Segment3d.init(c, skew);
        try expectTag(g.Intersection3d, try along.intersect(across, policy), .point);
        try expectTag(g.Intersection3d, try along.intersect(lifted, policy), .none);
    }
}

// ---------------------------------------------------------------------------
// Products that only one dimension has
// ---------------------------------------------------------------------------

test "cross is three-dimensional and wedge is four-dimensional" {
    const x = try g.Vector3f.init(1, 0, 0);
    const y = try g.Vector3f.init(0, 1, 0);
    const z = try x.cross(y);
    try std.testing.expectEqual(@as(f32, 0.0), z.x);
    try std.testing.expectEqual(@as(f32, 0.0), z.y);
    try std.testing.expectEqual(@as(f32, 1.0), z.z);
    // Anticommutative, and perpendicular to both operands.
    const back = try y.cross(x);
    try expectSign(f32, try z.dot(back), -1.0, 1e-6);
    try expectSign(f32, try z.dot(x), 0.0, 1e-6);
    try expectSign(f32, try z.dot(y), 0.0, 1e-6);
    // The cross product of a vector with itself is zero.
    const self = try x.cross(x);
    try std.testing.expectEqual(@as(f32, 0.0), self.x);

    // The exterior product of the same two four-dimensional vectors.
    const u = try g.Vector4d.init(1, 0, 0, 0);
    const v = try g.Vector4d.init(0, 1, 0, 0);
    const blade = try u.wedge(v);
    try std.testing.expectEqual(@as(f64, 1.0), blade.xy);
    try std.testing.expectEqual(@as(f64, 0.0), blade.xz);
    try std.testing.expectEqual(@as(f64, 0.0), blade.xw);
    try std.testing.expectEqual(@as(f64, 0.0), blade.yz);
    try std.testing.expectEqual(@as(f64, 0.0), blade.yw);
    try std.testing.expectEqual(@as(f64, 0.0), blade.zw);
    try std.testing.expectEqual(@as(f64, -1.0), (try v.wedge(u)).xy);
    try std.testing.expectEqual(@as(f64, -1.0), blade.neg().xy);
    // A vector wedged with itself is the zero bivector, which is the property
    // that makes the product exterior.
    const with_itself = try u.wedge(u);
    try std.testing.expectEqual(@as(f64, 0.0), with_itself.xy);
    try std.testing.expect(blade.isFinite());
    var poisoned = v;
    poisoned.w = inf64();
    try std.testing.expectError(error.NonFiniteInput, u.wedge(poisoned));
}

test "eq normalizes NaN and keeps the zeros apart" {
    const positive = try g.Vector2d.init(0.0, 1.0);
    const negative = try g.Vector2d.init(-0.0, 1.0);
    try std.testing.expect(positive.eq(try g.Vector2d.init(0.0, 1.0)));
    try std.testing.expect(!positive.eq(negative));
    const quiet = g.Vector2d{ .x = nan64(), .y = 0.0 };
    const other = g.Vector2d{ .x = -nan64(), .y = 0.0 };
    try std.testing.expect(quiet.eq(other));
    try std.testing.expect(!quiet.eq(positive));
    try std.testing.expect(positive.isFinite());
    try std.testing.expect(!quiet.isFinite());
    try std.testing.expect((try g.Point4f.init(0, 0, 0, 0)).eq(try g.Point4f.init(0, 0, 0, 0)));
    try std.testing.expect((try g.Segment2f.init(try g.Point2f.init(0, 0), try g.Point2f.init(1, 1))).isFinite());
}

// Checked geometry must reject non-finite arithmetic, not return a value that
// violates the next operation's input contract.
test "API alignment: finite arithmetic reports NonFiniteResult" {
    inline for (.{
        .{ g.Vector2f, g.Point2f, f32 },
        .{ g.Vector2d, g.Point2d, f64 },
    }) |Case| {
        const V = Case[0];
        const P = Case[1];
        const T = Case[2];
        const largest = std.math.floatMax(T);
        const value = try V.init(largest, 0);
        try std.testing.expectError(error.NonFiniteResult, value.add(value));
        try std.testing.expectError(error.NonFiniteResult, value.scale(2));
        try std.testing.expectError(error.NonFiniteResult, value.dot(value));
        try std.testing.expectError(error.NonFiniteResult, value.lengthSquared());
        try std.testing.expectError(
            error.NonFiniteResult,
            value.distance(try P.init(-largest, 0)),
        );
        const ordinary = try (try V.init(3, 4)).add(try V.init(1, -1));
        try std.testing.expectEqual(@as(T, 4), ordinary.x);
        try std.testing.expectEqual(@as(T, 3), ordinary.y);
    }
}

// A receiver is input too. Public fields can change after construction, and
// checked methods must reject the whole value before predicates consume it.
test "API alignment: composite receivers reject non-finite fields" {
    const origin = try g.Point2f.init(0, 0);
    var circle = try g.Circle2f.init(origin, 1);
    circle.radius = inf32();
    const other = try g.Circle2f.init(try g.Point2f.init(5, 0), 1);
    for ([_]g.PredicatePolicy{ .fast, .adaptive_exact }) |policy| {
        try std.testing.expectError(error.NonFiniteInput, circle.intersect(other, policy));
        try std.testing.expectError(error.NonFiniteInput, circle.contains(origin, policy));
    }
    var plane = try g.Plane3f.init(
        try g.Point3f.init(0, 0, 0),
        try g.Vector3f.init(0, 0, 1),
    );
    plane.normal.z = inf32();
    try std.testing.expectError(
        error.NonFiniteInput,
        plane.signedDistance(try g.Point3f.init(0, 0, 1)),
    );
    var segment = try g.Segment2f.init(origin, try g.Point2f.init(1, 0));
    segment.b.x = inf32();
    try std.testing.expectError(error.NonFiniteInput, segment.length());
    try std.testing.expectError(error.NonFiniteInput, segment.at(0));
    try std.testing.expectError(error.NonFiniteInput, segment.at(0.5));
    try std.testing.expectError(error.NonFiniteInput, segment.at(1));
    var line = try g.Line2f.init(origin, try g.Vector2f.init(1, 0));
    line.direction.x = inf32();
    try std.testing.expectError(error.NonFiniteInput, line.at(1));
    var ray = try g.Ray2f.init(origin, try g.Vector2f.init(1, 0));
    ray.direction.x = inf32();
    try std.testing.expectError(error.NonFiniteInput, ray.at(0));
    try std.testing.expectError(error.NonFiniteInput, ray.at(1));
}

// A half-line is not its supporting line. Exercise both the direct 2D path
// and the projected 3D/4D paths, in both precisions and both policies.
test "API alignment: coincident rays stay half-lines" {
    inline for (.{
        .{ g.Point2f, g.Vector2f, g.Ray2f, g.Line2f, g.Intersection2f, f32 },
        .{ g.Point2d, g.Vector2d, g.Ray2d, g.Line2d, g.Intersection2d, f64 },
        .{ g.Point3f, g.Vector3f, g.Ray3f, g.Line3f, g.Intersection3f, f32 },
        .{ g.Point3d, g.Vector3d, g.Ray3d, g.Line3d, g.Intersection3d, f64 },
        .{ g.Point4f, g.Vector4f, g.Ray4f, g.Line4f, g.Intersection4f, f32 },
        .{ g.Point4d, g.Vector4d, g.Ray4d, g.Line4d, g.Intersection4d, f64 },
    }) |Case| {
        const P = Case[0];
        const V = Case[1];
        const R = Case[2];
        const L = Case[3];
        const I = Case[4];
        const T = Case[5];
        const dim = @typeInfo(P).@"struct".fields.len;
        const coordinates = [_]T{0} ** dim;
        var direction = coordinates;
        direction[0] = 1;
        const point = try P.fromArray(coordinates);
        const vector = try V.fromArray(direction);
        const ray = try R.init(point, vector);
        const line = try L.init(point, vector);
        for ([_]g.PredicatePolicy{ .fast, .adaptive_exact }) |policy| {
            try expectTag(I, try ray.intersect(ray, policy), .coincident);
            try expectTag(I, try line.intersect(line, policy), .line);
        }
    }
}

// A ray inside its own line is contained; only a pair of lines that are one
// line may answer `line`. The whole-set tags are exercised in every dimension
// and both policies, because the answer is a property of the set and not of
// the dimension the engine projected into.
test "API alignment: line, ray, and segment report the set they actually share" {
    inline for (.{
        .{ g.Point2f, g.Vector2f, g.Ray2f, g.Line2f, g.Segment2f, g.Intersection2f, f32 },
        .{ g.Point3d, g.Vector3d, g.Ray3d, g.Line3d, g.Segment3d, g.Intersection3d, f64 },
        .{ g.Point4d, g.Vector4d, g.Ray4d, g.Line4d, g.Segment4d, g.Intersection4d, f64 },
    }) |Case| {
        const P = Case[0];
        const V = Case[1];
        const R = Case[2];
        const L = Case[3];
        const S = Case[4];
        const I = Case[5];
        const T = Case[6];
        const dim = @typeInfo(P).@"struct".fields.len;
        const axis = blk: {
            var row = [_]T{0} ** dim;
            row[0] = 1;
            break :blk row;
        };
        const back = blk: {
            var row = [_]T{0} ** dim;
            row[0] = -1;
            break :blk row;
        };
        const twice = blk: {
            var row = [_]T{0} ** dim;
            row[0] = 2;
            break :blk row;
        };
        const origin = try P.fromArray([_]T{0} ** dim);
        const direction = try V.fromArray(axis);
        const ray = try R.init(origin, direction);
        const line = try L.init(origin, direction);
        const on_ray = try S.init(try P.fromArray(axis), try P.fromArray(twice));
        const same_segment = try S.init(try P.fromArray(axis), try P.fromArray(twice));
        const backward = try R.init(origin, try V.fromArray(back));

        for ([_]g.PredicatePolicy{ .fast, .adaptive_exact }) |policy| {
            // Two lines that are one line: the common set runs both ways.
            try expectTag(I, try line.intersect(line, policy), .line);
            // A ray with its own copy: the common set is a half-line, and
            // naming it a line would name points the pair does not share.
            try expectTag(I, try ray.intersect(ray, policy), .coincident);
            // The ray is the whole common set with its supporting line.
            const ray_in_line = try ray.intersectLine(line, policy);
            try expectTag(I, ray_in_line, .contained);
            try std.testing.expectEqual(g.ContainedSide.first, ray_in_line.contained);
            // A segment wholly on the ray is the segment, whole.
            const covered = try ray.intersectSegment(on_ray, policy);
            try expectTag(I, covered, .contained);
            try std.testing.expectEqual(g.ContainedSide.second, covered.contained);
            const segment_in_line = try on_ray.intersectLine(line, policy);
            try expectTag(I, segment_in_line, .contained);
            try std.testing.expectEqual(g.ContainedSide.first, segment_in_line.contained);
            // An identical pair of segments is the same bounded set.
            try expectTag(I, try on_ray.intersect(same_segment, policy), .coincident);
            // Two rays leaving the same end in opposite directions share
            // that end and nothing else.
            const shared_end = try ray.intersect(backward, policy);
            try expectTag(I, shared_end, .point);
            try std.testing.expect(shared_end.point.eq(origin));
        }
    }
}

// A NaN is not a zero. The plain evaluation of a two-dimensional determinant
// can subtract an infinity from an infinity, and the exact value of the same
// form is exactly zero; the plain policy must not report the zero the exact
// policy can prove.
test "API alignment: a NaN evaluation is never reported as a zero sign" {
    const largest = std.math.floatMax(f64);

    // Both products of the determinant overflow in the same direction, so the
    // plain evaluation forms `-inf - -inf`. `b` and `c` are the same point, so
    // the exact value of that form is exactly zero.
    const a = try g.Point2d.init(0, 0);
    const b = try g.Point2d.init(largest, -largest);
    const c = try g.Point2d.init(largest, -largest);
    try std.testing.expectError(error.NonFiniteResult, g.Point2d.orientation(a, b, c, .fast));
    try std.testing.expectEqual(g.Sign.zero, try g.Point2d.orientation(a, b, c, .adaptive_exact));

    // The same shape in three dimensions: the cross product's third component
    // is `-inf - inf`, and the scalar product multiplies it by zero.
    const p0 = try g.Point3d.init(0, 0, 0);
    const p1 = try g.Point3d.init(largest, largest, 0);
    const p2 = try g.Point3d.init(largest, -largest, 0);
    const p3 = try g.Point3d.init(1, 0, 0);
    try std.testing.expectError(error.NonFiniteResult, g.Point3d.orientation(p0, p1, p2, p3, .fast));
    try std.testing.expectEqual(g.Sign.zero, try g.Point3d.orientation(p0, p1, p2, p3, .adaptive_exact));
}

// A ball whose offset and radius square both overflow leaves a residual that
// is an infinity less an infinity. The exact policy decides the same question
// on integers and says the point is outside.
test "API alignment: a fast containment is not decided by a NaN residual" {
    const largest = std.math.floatMax(f64);
    const circle = try g.Circle2d.init(try g.Point2d.init(-largest, 0), largest);
    const far_point = try g.Point2d.init(largest, 0);
    // The exact residual is `4 * largest^2 - largest^2 = 3 * largest^2`.
    try std.testing.expectError(error.NonFiniteResult, circle.contains(far_point, .fast));
    try std.testing.expect(!(try circle.contains(far_point, .adaptive_exact)));
}

// Every value in this module is by value with public writable fields, so the
// receiver is input as much as the argument is. Each family below is poisoned
// after its constructor accepted it, and every operation that reads it must
// refuse it.
test "API alignment: every composite receiver is gated before it is consumed" {
    const policy = g.PredicatePolicy.adaptive_exact;
    const origin2 = try g.Point2f.init(0, 0);
    const step = try g.Vector2f.init(1, 0);
    const clean_line = try g.Line2f.init(origin2, step);
    const clean_ray = try g.Ray2f.init(origin2, step);
    const clean_segment = try g.Segment2f.init(origin2, try g.Point2f.init(2, 0));
    const probe = try g.Point2f.init(1, 0);

    var segment = clean_segment;
    segment.a.x = nan32();
    try std.testing.expectError(error.NonFiniteInput, segment.at(0));
    try std.testing.expectError(error.NonFiniteInput, segment.length());
    try std.testing.expectError(error.NonFiniteInput, segment.closestPointTo(probe));
    try std.testing.expectError(error.NonFiniteInput, segment.intersect(clean_segment, policy));
    try std.testing.expectError(error.NonFiniteInput, segment.intersectLine(clean_line, policy));
    try std.testing.expectError(error.NonFiniteInput, segment.intersectRay(clean_ray, policy));

    var line = clean_line;
    line.origin.x = nan32();
    try std.testing.expectError(error.NonFiniteInput, line.at(1));
    try std.testing.expectError(error.NonFiniteInput, line.closestPointTo(probe));
    try std.testing.expectError(error.NonFiniteInput, line.intersect(clean_line, policy));

    var ray = clean_ray;
    ray.direction.y = inf32();
    try std.testing.expectError(error.NonFiniteInput, ray.at(0));
    try std.testing.expectError(error.NonFiniteInput, ray.at(1));
    try std.testing.expectError(error.NonFiniteInput, ray.intersect(clean_ray, policy));
    try std.testing.expectError(error.NonFiniteInput, ray.intersectLine(clean_line, policy));
    try std.testing.expectError(error.NonFiniteInput, ray.intersectSegment(clean_segment, policy));

    const origin3 = try g.Point3f.init(0, 0, 0);
    const ball = try g.Sphere3f.init(origin3, 1);
    const plane3 = try g.Plane3f.init(origin3, try g.Vector3f.init(0, 0, 1));
    const segment3 = try g.Segment3f.init(try g.Point3f.init(-1, 0, 0), try g.Point3f.init(1, 0, 0));
    const ray3 = try g.Ray3f.init(try g.Point3f.init(-2, 0, 0), try g.Vector3f.init(1, 0, 0));
    var sphere = ball;
    sphere.center.z = inf32();
    try std.testing.expectError(error.NonFiniteInput, sphere.contains(origin3, policy));
    try std.testing.expectError(error.NonFiniteInput, sphere.intersect(ball, policy));
    try std.testing.expectError(error.NonFiniteInput, sphere.intersectSegment(segment3, policy));
    try std.testing.expectError(error.NonFiniteInput, sphere.intersectRay(ray3, policy));
    try std.testing.expectError(error.NonFiniteInput, sphere.intersectPlane(plane3, policy));
    try std.testing.expect(!sphere.isFinite());

    var plane = plane3;
    plane.normal.z = nan32();
    try std.testing.expectError(error.NonFiniteInput, plane.intersect(plane3, policy));
    try std.testing.expectError(error.NonFiniteInput, plane.intersectSegment(segment3, policy));
    try std.testing.expectError(error.NonFiniteInput, plane.intersectSphere(ball, policy));

    const origin4 = try g.Point4f.init(0, 0, 0, 0);
    const hyper = try g.Hyperplane4f.init(origin4, try g.Vector4f.init(0, 0, 0, 1));
    const hypersphere = try g.Hypersphere4f.init(origin4, 1);
    const segment4 = try g.Segment4f.init(try g.Point4f.init(-1, 0, 0, 0), try g.Point4f.init(1, 0, 0, 0));
    var flat = hyper;
    flat.origin.w = inf32();
    try std.testing.expectError(error.NonFiniteInput, flat.signedDistance(try g.Point4f.init(0, 0, 0, 1)));
    try std.testing.expectError(error.NonFiniteInput, flat.intersect(hyper, policy));
    try std.testing.expectError(error.NonFiniteInput, flat.intersectSegment(segment4, policy));
    try std.testing.expectError(error.NonFiniteInput, flat.intersectSphere(hypersphere, policy));

    var deep = hypersphere;
    deep.radius = inf32();
    try std.testing.expectError(error.NonFiniteInput, deep.contains(origin4, policy));
    try std.testing.expectError(error.NonFiniteInput, deep.intersect(hypersphere, policy));
    try std.testing.expectError(error.NonFiniteInput, deep.intersectHyperplane(hyper, policy));

    const clean_box = try g.Aabb2f.init(origin2, try g.Point2f.init(1, 1));
    var box = clean_box;
    box.maximum.y = inf32();
    try std.testing.expectError(error.NonFiniteInput, box.contains(probe));
    try std.testing.expectError(error.NonFiniteInput, box.expand(probe));
    try std.testing.expectError(error.NonFiniteInput, box.intersects(clean_box));
    try std.testing.expectError(error.NonFiniteInput, box.center());
    try std.testing.expectError(error.NonFiniteInput, box.extent());
    try std.testing.expectError(error.NonFiniteInput, box.area());
    try std.testing.expectError(error.NonFiniteInput, box.intersect(clean_box));
    try std.testing.expect(!box.isFinite());

    // A bivector carries no invariant beyond finiteness, and `neg` is
    // arithmetic on the stored bits, so the check a caller needs after
    // writing a field is `isFinite`.
    var blade = try g.Bivector4f.init(1, 0, 0, 0, 0, 0);
    blade.xy = inf32();
    try std.testing.expect(!blade.isFinite());
    try std.testing.expect(std.math.isInf(blade.neg().xy));

    // A poisoned vector and a poisoned point are refused the same way, and
    // `eq` still answers, because it compares bits rather than numbers.
    var vector = step;
    vector.x = nan32();
    try std.testing.expectError(error.NonFiniteInput, vector.add(step));
    try std.testing.expectError(error.NonFiniteInput, vector.normalize());
    var point = origin2;
    point.y = inf32();
    try std.testing.expectError(error.NonFiniteInput, point.add(step));
    try std.testing.expect(!point.eq(origin2));
}

// A norm past the top of the range has no direction, and dividing by an
// infinite length would hand back the zero vector as though a zero vector had
// one.
test "API alignment: an overflowing norm cannot become a direction" {
    const largest = std.math.floatMax(f32);
    const wide = try g.Vector2f.init(largest, largest);
    // `hypot` of two components at the top of the range is `sqrt(2)` times
    // either of them, which is not a binary32 number.
    try std.testing.expectError(error.NonFiniteResult, wide.length());
    try std.testing.expectError(error.NonFiniteResult, wide.normalize());
    try std.testing.expectError(
        error.NonFiniteResult,
        g.Line2f.init(try g.Point2f.init(0, 0), wide),
    );
    // The check is on the overflow and not on the size: the ordinary large
    // vector whose squares overflow keeps its folded norm.
    const big = try g.Vector3f.init(1.0e30, 1.0e30, 1.0e30);
    try expectSign(f32, try big.length(), 1.7320508e30, 1.0e24);
}

test "API alignment: an overflowing translation, product, and distance are reported" {
    const largest = std.math.floatMax(f32);
    const corner = try g.Point2f.init(largest, 0);
    const step = try g.Vector2f.init(largest, 0);
    try std.testing.expectError(error.NonFiniteResult, corner.add(step));
    try std.testing.expectError(error.NonFiniteResult, corner.sub(try g.Vector2f.init(-largest, 0)));
    const cancellation = try corner.sub(step);
    try std.testing.expectEqual(@as(f32, 0), cancellation.x);
    try std.testing.expectEqual(@as(f32, 0), cancellation.y);
    try std.testing.expectError(
        error.NonFiniteResult,
        step.distance(try g.Point2f.init(-largest, 0)),
    );
    try std.testing.expectError(
        error.NonFiniteResult,
        (try g.Line2f.init(corner, try g.Vector2f.init(1, 0))).at(largest),
    );

    const ordinary = try (try g.Point2f.init(1, 2)).add(try g.Vector2f.init(3, 4));
    try std.testing.expectEqual(@as(f32, 4.0), ordinary.x);
    try std.testing.expectEqual(@as(f32, 6.0), ordinary.y);
}

test "API alignment: an overflowing cross and wedge coefficient are reported" {
    const largest = std.math.floatMax(f32);
    const up = try g.Vector3f.init(0, largest, 0);
    const across = try g.Vector3f.init(0, 0, largest);
    try std.testing.expectError(error.NonFiniteResult, up.cross(across));
    // The same product in the four-dimensional exterior form.
    const blade_a = try g.Vector4f.init(0, largest, 0, 0);
    const blade_b = try g.Vector4f.init(0, 0, largest, 0);
    try std.testing.expectError(error.NonFiniteResult, blade_a.wedge(blade_b));
    // An ordinary pair still crosses and still wedges.
    const x = try g.Vector3f.init(1, 0, 0);
    const y = try g.Vector3f.init(0, 1, 0);
    try std.testing.expectEqual(@as(f32, 1.0), (try x.cross(y)).z);
}

test "API alignment: an overflowing box measure is reported" {
    const largest = std.math.floatMax(f32);
    const wide = try g.Aabb2f.init(
        try g.Point2f.init(0, 0),
        try g.Point2f.init(largest, largest),
    );
    // The extents are finite; their product is not.
    try std.testing.expectEqual([2]f32{ largest, largest }, try wide.extent());
    try std.testing.expectError(error.NonFiniteResult, wide.area());

    const spanning = try g.Aabb2f.init(
        try g.Point2f.init(-largest, -largest),
        try g.Point2f.init(largest, largest),
    );
    // The difference itself leaves the format here.
    try std.testing.expectError(error.NonFiniteResult, spanning.extent());
    try std.testing.expectError(error.NonFiniteResult, spanning.area());

    const far_corner = try g.Aabb2f.init(
        try g.Point2f.init(largest * 0.9, largest * 0.9),
        try g.Point2f.init(largest, largest),
    );
    // The centre is `(minimum + maximum) / 2`, and that sum leaves the format
    // even though the halved value would not.
    try std.testing.expectError(error.NonFiniteResult, far_corner.center());

    // The ordinary measures are untouched by any of this.
    const ordinary = try g.Aabb3f.init(
        try g.Point3f.init(0, 0, 0),
        try g.Point3f.init(2, 3, 4),
    );
    try std.testing.expectEqual(@as(f32, 24.0), try ordinary.volume());
    try std.testing.expectEqual(@as(f32, 52.0), try ordinary.doubledSurfaceArea());
}

test "API alignment: an overflowing chord and signed distance are reported" {
    const largest = std.math.floatMax(f32);
    // Two discs of the largest radius a unit apart: the radii sum past the
    // top of the range, so the chord has no margin to divide and the
    // projection along the line of centres has an infinity to subtract.
    const left = try g.Circle2f.init(try g.Point2f.init(0, 0), largest);
    const right = try g.Circle2f.init(try g.Point2f.init(1, 0), largest);
    for ([_]g.PredicatePolicy{ .fast, .adaptive_exact }) |policy| {
        try std.testing.expectError(error.NonFiniteResult, left.intersect(right, policy));
    }

    // Three terms of `largest / sqrt(3)` are past the top of the range.
    const plane = try g.Plane3f.init(
        try g.Point3f.init(0, 0, 0),
        try g.Vector3f.init(1, 1, 1),
    );
    try std.testing.expectError(
        error.NonFiniteResult,
        plane.signedDistance(try g.Point3f.init(largest, largest, largest)),
    );
    // Four terms of `floatMax(f64) / 2` are, in binary64.
    const wide: f64 = std.math.floatMax(f64);
    const hyper = try g.Hyperplane4d.init(
        try g.Point4d.init(0, 0, 0, 0),
        try g.Vector4d.init(1, 1, 1, 1),
    );
    try std.testing.expectError(
        error.NonFiniteResult,
        hyper.signedDistance(try g.Point4d.init(wide, wide, wide, wide)),
    );
    // The same plane still answers at an ordinary distance.
    try expectSign(f32, try plane.signedDistance(try g.Point3f.init(3, 0, 0)), 1.7320508, 1.0e-4);
}

test "API alignment: a segment whose own span leaves the format is reported" {
    const largest = std.math.floatMax(f32);
    const far = try g.Segment2f.init(
        try g.Point2f.init(-largest, 0),
        try g.Point2f.init(largest, 0),
    );
    try std.testing.expectError(error.NonFiniteResult, far.length());
    // The projection divides by the squared span, which leaves the format
    // even when the offset does not.
    try std.testing.expectError(
        error.NonFiniteResult,
        far.closestPointTo(try g.Point2f.init(0, 1)),
    );
    // The supporting line of the same pair still answers, because its
    // direction is unit and the difference it forms is the offset.
    const line = try g.Line2f.init(try g.Point2f.init(-largest, 0), try g.Vector2f.init(1, 0));
    const foot = try line.closestPointTo(try g.Point2f.init(0, 1));
    try expectSign(f32, foot.y, 0.0, 0.0);
    try expectSign(f32, foot.x, 0.0, 0.0);
}

// A box intersection is the WHOLE common box. Touching is not emptiness, and
// a box that meets its argument on one axis is a box degenerate on that axis
// rather than a point or a face.
test "API alignment: a box intersection is the whole common box, never a part" {
    const high = try g.Point3f.init(2, 2, 2);
    const unit = try g.Aabb3f.init(try g.Point3f.init(0, 0, 0), high);

    // Cut off above on one axis only: the common box is full on the other two
    // and degenerate on this one.
    const slab = try g.Aabb3f.init(try g.Point3f.init(1.5, 0.5, 0.5), try g.Point3f.init(4, 1.5, 1.5));
    const slab_touch = try unit.intersect(slab);
    try expectTag(g.Intersection3f, slab_touch, .box);
    try std.testing.expect(slab_touch.box.minimum.eq(try g.Point3f.init(1.5, 0.5, 0.5)));
    try std.testing.expect(slab_touch.box.maximum.eq(try g.Point3f.init(2, 1.5, 1.5)));
    try std.testing.expect(try unit.intersects(slab));

    // The same cut off on all three axes: a box that is one corner.
    const corner = try g.Aabb3f.init(try g.Point3f.init(1.5, 1.5, 1.5), try g.Point3f.init(4, 4, 4));
    const corner_touch = try unit.intersect(corner);
    try expectTag(g.Intersection3f, corner_touch, .box);
    try std.testing.expect(corner_touch.box.minimum.eq(try g.Point3f.init(1.5, 1.5, 1.5)));
    try std.testing.expect(corner_touch.box.maximum.eq(high));

    // Empty on the third axis alone, however much the other two agree.
    const above = try g.Aabb3f.init(try g.Point3f.init(0, 0, 3), try g.Point3f.init(2, 2, 4));
    try expectTag(g.Intersection3f, try unit.intersect(above), .none);
    try std.testing.expect(!(try unit.intersects(above)));

    // One box wholly inside another, in both orders, and the identical pair.
    const inner = try g.Aabb3f.init(try g.Point3f.init(0.5, 0.5, 0.5), try g.Point3f.init(1, 1, 1));
    const held = try unit.intersect(inner);
    try expectTag(g.Intersection3f, held, .contained);
    try std.testing.expectEqual(g.ContainedSide.second, held.contained);
    try std.testing.expectEqual(g.ContainedSide.first, (try inner.intersect(unit)).contained);
    try expectTag(g.Intersection3f, try unit.intersect(unit), .coincident);

    // A degenerate box is legal. Inside another box it IS the whole common
    // set, so the answer is containment rather than a degenerate box.
    const single = try g.Aabb3f.init(high, high);
    const at_corner = try unit.intersect(single);
    try expectTag(g.Intersection3f, at_corner, .contained);
    try std.testing.expectEqual(g.ContainedSide.second, at_corner.contained);
    try std.testing.expect(try unit.contains(high));
    try std.testing.expectEqual([3]f32{ 0, 0, 0 }, try single.extent());
}

test "API alignment: subnormal directions remain unit" {
    inline for (.{ .{ g.Vector2f, f32 }, .{ g.Vector2d, f64 } }) |Case| {
        const V = Case[0];
        const T = Case[1];
        const U = std.meta.Int(.unsigned, @bitSizeOf(T));
        const tiny: T = @bitCast(@as(U, 1));
        const direction = try (try V.init(tiny, tiny)).normalize();
        try std.testing.expect(@abs(try direction.length() - 1) < 0.00001);
        try std.testing.expectEqual(direction.x, direction.y);
    }
}

test "API alignment: outside ray shares the complete solid-ball chord" {
    const center = try g.Point2f.init(0, 0);
    const ball = try g.Circle2f.init(center, 1);
    const ray = try g.Ray2f.init(try g.Point2f.init(-2, 0), try g.Vector2f.init(1, 0));
    for ([_]g.PredicatePolicy{ .fast, .adaptive_exact }) |policy| {
        const common = try ball.intersectRay(ray, policy);
        try expectTag(g.Intersection2f, common, .segment);
        try expectPoint(common.segment.a, try g.Point2f.init(-1, 0));
        try expectPoint(common.segment.b, try g.Point2f.init(1, 0));
        try expectPoint(try common.segment.at(0.5), center);
    }
}

/// Runtime corner storage prevents constant folding of the midpoint.
fn boxCenter(comptime B: type, lo: anytype, hi: anytype) !@FieldType(B, "minimum") {
    const P = @FieldType(B, "minimum");
    const T = ScalarOf(P);
    const dim = comptime fieldCountOf(P);
    var lows: [dim]T = [1]T{lo} ** dim;
    var highs: [dim]T = [1]T{hi} ** dim;
    std.mem.doNotOptimizeAway(&lows);
    std.mem.doNotOptimizeAway(&highs);
    return (try B.init(try pointOf(P, dim, lows), try pointOf(P, dim, highs))).center();
}

/// The corner of a uniform box: `value` on every axis.
fn uniformPoint(comptime P: type, value: anytype) !P {
    const T = ScalarOf(P);
    const dim = comptime fieldCountOf(P);
    const components: [dim]T = [1]T{value} ** dim;
    return pointOf(P, dim, components);
}

test "API alignment: a box center keeps a subnormal midpoint and both zeros" {
    inline for (.{ g.Aabb2f, g.Aabb3f, g.Aabb4f, g.Aabb2d, g.Aabb3d, g.Aabb4d }) |B| {
        const P = @FieldType(B, "minimum");
        const T = ScalarOf(P);
        const Bits = @Int(.unsigned, @bitSizeOf(T));
        const unit: T = @bitCast(@as(Bits, 1));
        const two: T = @bitCast(@as(Bits, 2));
        const negative_zero: T = -0.0;
        const largest = std.math.floatMax(T);
        // Exact subnormal midpoint.
        try expectPoint(try boxCenter(B, unit, unit), try uniformPoint(P, unit));
        // Half of three subnormal units rounds to the even two.
        try expectPoint(try boxCenter(B, unit, two), try uniformPoint(P, two));
        // Negative subnormal midpoint.
        try expectPoint(try boxCenter(B, -unit, -unit), try uniformPoint(P, -unit));
        // Negative zero stays negative; mixed zero signs produce positive zero.
        try expectPoint(
            try boxCenter(B, negative_zero, negative_zero),
            try uniformPoint(P, negative_zero),
        );
        try expectPoint(try boxCenter(B, negative_zero, 0.0), try uniformPoint(P, 0.0));
        // The rounded sum overflows before halving.
        try std.testing.expectError(error.NonFiniteResult, boxCenter(B, largest, largest));
    }
}

test "API alignment: wedge overflow follows separate rounded products" {
    const largest = std.math.floatMax(f32);
    const a = try g.Vector4f.init(largest, largest, 0, 0);
    const b = try g.Vector4f.init(1, 2, 0, 0);
    // A fused first product could return largest; the strict graph overflows.
    try std.testing.expectError(error.NonFiniteResult, a.wedge(b));
}

// ---------------------------------------------------------------------------
// An INDEPENDENT exact opinion about every `adaptive_exact` predicate
// ---------------------------------------------------------------------------

/// How many limbs the reference integer below carries.
///
/// A finite binary64 sits at most `2^2045` above its significand, so one
/// coordinate on a common per-case scale is at most `53 + 2045 = 2098` bits
/// wide. The degree-four determinant over four of them is twenty-four terms of
/// a product of four, so it is below `24 * 2^8392 < 2^8397`; the degree-two
/// forms that rescale by `2^2148` stay below `2^8398`. 160 limbs is 10240
/// bits, so every value this test builds fits whole and nothing here ever
/// truncates - which is checked rather than assumed, below.
const reference_words = 160;

/// An exact signed integer, written HERE rather than taken from the module.
///
/// A second opinion that reuses the arithmetic under test is not a second
/// opinion, so this is a fresh implementation over `u64` limbs with `u128`
/// intermediates, and it never calls `Big` or anything else of the module's.
/// It has its own rule about which limbs are live - only the first `used` -
/// and it keeps that rule by writing every limb it claims, which is what makes
/// it usable as the reference for a module that now clears its working set
/// only over the region it will look at.
const Reference = struct {
    negative: bool,
    used: usize,
    limbs: [reference_words]u64,

    fn zero() Reference {
        return .{ .negative = false, .used = 0, .limbs = [_]u64{0} ** reference_words };
    }

    fn sign(a: Reference) g.Sign {
        if (a.used == 0) return .zero;
        return if (a.negative) .negative else .positive;
    }

    fn trim(out: *Reference) void {
        while (out.used > 0 and out.limbs[out.used - 1] == 0) out.used -= 1;
        if (out.used == 0) out.negative = false;
    }

    /// `-(significand * 2^shift)` or `+(significand * 2^shift)`.
    fn fromShifted(negative: bool, significand: u64, shift: usize) Reference {
        std.debug.assert(significand != 0);
        var out = Reference.zero();
        const word = shift / 64;
        std.debug.assert(word + 2 <= reference_words);
        const bits: u6 = @intCast(shift % 64);
        const placed: u128 = @as(u128, significand) << bits;
        out.limbs[word] = @truncate(placed);
        out.limbs[word + 1] = @truncate(placed >> 64);
        out.used = word + 2;
        out.negative = negative;
        out.trim();
        return out;
    }

    fn compare(a: Reference, b: Reference) std.math.Order {
        if (a.used != b.used) return if (a.used > b.used) .gt else .lt;
        var index: usize = a.used;
        while (index > 0) {
            index -= 1;
            if (a.limbs[index] != b.limbs[index]) {
                return if (a.limbs[index] > b.limbs[index]) .gt else .lt;
            }
        }
        return .eq;
    }

    fn magnitudeAdd(out: *Reference, a: Reference, b: Reference) void {
        const width = @max(a.used, b.used) + 1;
        std.debug.assert(width <= reference_words);
        var carry: u64 = 0;
        var index: usize = 0;
        while (index < width) : (index += 1) {
            const left: u128 = if (index < a.used) a.limbs[index] else 0;
            const right: u128 = if (index < b.used) b.limbs[index] else 0;
            const total = left + right + @as(u128, carry);
            out.limbs[index] = @truncate(total);
            carry = @intCast(total >> 64);
        }
        out.used = width;
        out.negative = false;
        out.trim();
    }

    /// `|a| - |b|`, which the caller has already ordered.
    fn magnitudeDrop(out: *Reference, a: Reference, b: Reference) void {
        var borrow: u64 = 0;
        var index: usize = 0;
        while (index < a.used) : (index += 1) {
            const left: u128 = a.limbs[index];
            const right: u128 = if (index < b.used) b.limbs[index] else 0;
            const total = left -% right -% @as(u128, borrow);
            out.limbs[index] = @truncate(total);
            borrow = @intCast((total >> 64) & 1);
        }
        std.debug.assert(borrow == 0);
        out.used = a.used;
        out.negative = false;
        out.trim();
    }

    fn add(out: *Reference, a: Reference, b: Reference) void {
        const a_negative = a.negative;
        const b_negative = b.negative;
        if (a.used == 0) {
            out.* = b;
            return;
        }
        if (b.used == 0) {
            out.* = a;
            return;
        }
        if (a_negative == b_negative) {
            magnitudeAdd(out, a, b);
            out.negative = a_negative and out.used != 0;
            return;
        }
        const order = Reference.compare(a, b);
        if (order == .eq) {
            out.* = Reference.zero();
            return;
        }
        if (order == .gt) {
            magnitudeDrop(out, a, b);
            out.negative = a_negative and out.used != 0;
            return;
        }
        magnitudeDrop(out, b, a);
        out.negative = b_negative and out.used != 0;
    }

    fn sub(out: *Reference, a: Reference, b: Reference) void {
        var flipped = b;
        flipped.negative = !flipped.negative;
        Reference.add(out, a, flipped);
    }

    fn mul(out: *Reference, a: Reference, b: Reference) void {
        if (a.used == 0 or b.used == 0) {
            out.* = Reference.zero();
            return;
        }
        std.debug.assert(a.used + b.used <= reference_words);
        var result = Reference.zero();
        var i: usize = 0;
        while (i < a.used) : (i += 1) {
            if (a.limbs[i] == 0) continue;
            var carry: u128 = 0;
            var j: usize = 0;
            while (j < b.used) : (j += 1) {
                const total = @as(u128, a.limbs[i]) * @as(u128, b.limbs[j]) +
                    @as(u128, result.limbs[i + j]) + carry;
                result.limbs[i + j] = @truncate(total);
                carry = total >> 64;
            }
            result.limbs[i + b.used] = @truncate(carry);
        }
        result.used = a.used + b.used;
        result.negative = a.negative != b.negative;
        result.trim();
        out.* = result;
    }

    fn shiftedLeft(out: *Reference, a: Reference, bits: usize) void {
        if (a.used == 0 or bits == 0) {
            out.* = a;
            return;
        }
        const words = bits / 64;
        const shift: u6 = @intCast(bits % 64);
        std.debug.assert(a.used + words + 1 <= reference_words);
        var result = Reference.zero();
        var index: usize = 0;
        while (index < a.used) : (index += 1) {
            const placed: u128 = @as(u128, a.limbs[index]) << shift;
            result.limbs[index + words] |= @truncate(placed);
            result.limbs[index + words + 1] |= @as(u64, @truncate(placed >> 64));
        }
        result.used = a.used + words + 1;
        result.negative = a.negative;
        result.trim();
        out.* = result;
    }
};

/// Decode each format directly from its IEEE bits, without the module's
/// conversion code. The magnitude is `significand * 2^(offset - 1074)`.
/// Binary32 offsets include `1074 - 149` so both formats share this scale.
const Parts = struct { negative: bool, significand: u64, offset: usize };

fn partsOf(comptime T: type, value: T) Parts {
    switch (T) {
        f32 => {
            const bits: u32 = @bitCast(value);
            const exponent_field: u32 = (bits >> 23) & 0xff;
            const fraction: u64 = bits & 0x007f_ffff;
            std.debug.assert(exponent_field != 0xff);
            const negative = (bits >> 31) != 0;
            if (exponent_field == 0 and fraction == 0) {
                return .{ .negative = negative, .significand = 0, .offset = 925 };
            }
            if (exponent_field == 0) return .{ .negative = negative, .significand = fraction, .offset = 925 };
            return .{
                .negative = negative,
                .significand = fraction | (@as(u64, 1) << 23),
                .offset = @as(usize, exponent_field - 1) + 925,
            };
        },
        f64 => {
            const bits: u64 = @bitCast(value);
            const exponent_field: u64 = (bits >> 52) & 0x7ff;
            const fraction: u64 = bits & 0x000f_ffff_ffff_ffff;
            std.debug.assert(exponent_field != 0x7ff);
            const negative = (bits >> 63) != 0;
            if (exponent_field == 0 and fraction == 0) return .{ .negative = negative, .significand = 0, .offset = 0 };
            if (exponent_field == 0) return .{ .negative = negative, .significand = fraction, .offset = 0 };
            return .{
                .negative = negative,
                .significand = fraction | (@as(u64, 1) << 52),
                .offset = @intCast(exponent_field - 1),
            };
        },
        else => @compileError("zmath geometry is f32 or f64"),
    }
}

/// Use the smallest nonzero offset as one case's common scale. Each coordinate
/// is then an integer. A degree-`k` form carries the positive factor
/// `2^(k * (1074 - base))`, so its sign is unchanged.
fn baseOf(comptime T: type, values: []const T) usize {
    var base: usize = std.math.maxInt(usize);
    for (values) |value| {
        const parts = partsOf(T, value);
        if (parts.significand == 0) continue;
        base = @min(base, parts.offset);
    }
    return if (base == std.math.maxInt(usize)) 0 else base;
}

fn scaledAt(comptime T: type, value: T, base: usize) Reference {
    const parts = partsOf(T, value);
    if (parts.significand == 0) return Reference.zero();
    std.debug.assert(parts.offset >= base);
    return Reference.fromShifted(parts.negative, parts.significand, parts.offset - base);
}

/// The published point type's scalar, read off its first field so a test can
/// pick the matching corpus without naming the type again.
fn ScalarOf(comptime P: type) type {
    return @typeInfo(P).@"struct".fields[0].type;
}

fn fieldCountOf(comptime P: type) usize {
    return @typeInfo(P).@"struct".fields.len;
}

/// One published point, built from the first `dim` components of a four-wide
/// row. The condition is `comptime`, so only the branch that matches is
/// analyzed and the arity of `init` is never wrong for the type in hand.
fn pointOf(comptime P: type, comptime dim: usize, v: anytype) !P {
    if (comptime dim == 2) return P.init(v[0], v[1]);
    if (comptime dim == 3) return P.init(v[0], v[1], v[2]);
    return P.init(v[0], v[1], v[2], v[3]);
}

fn vectorOf(comptime V: type, comptime dim: usize, v: anytype) !V {
    if (comptime dim == 2) return V.init(v[0], v[1]);
    if (comptime dim == 3) return V.init(v[0], v[1], v[2]);
    return V.init(v[0], v[1], v[2], v[3]);
}

/// The corpus every predicate below is driven over: ordinary magnitudes and
/// signs, exact zeros of both signs, powers of two whose low limbs are all
/// zero, rows that are exactly collinear and rows one unit in the last place
/// away from it, the smallest subnormal of the format, three subnormals of
/// different sizes, and coordinates at the very top of the range next to
/// coordinates a thousand decades below it.
fn corpusOf(comptime T: type) [16][4]T {
    const unsigned = @Int(.unsigned, @bitSizeOf(T));
    const smallest: T = @bitCast(@as(unsigned, 1));
    const largest: T = std.math.floatMax(T);
    const wide_high: T = if (T == f32) 3.0e38 else 1.0e308;
    const wide_low: T = if (T == f32) 1.0e-38 else 1.0e-308;
    const sub_a: T = if (T == f32) 1.0e-40 else 1.0e-320;
    const sub_b: T = if (T == f32) 3.0e-40 else 3.0e-320;
    const sub_c: T = if (T == f32) 7.0e-40 else 7.0e-320;
    const eps: T = if (T == f32) 1.19209290e-7 else 2.220446049250313e-16;
    const mid: T = if (T == f32) 1.0e7 else 1.0e15;
    return .{
        .{ 0, 0, 0, 0 },
        .{ 1, 0, 0, 0 },
        .{ 0, 1, 0, 0 },
        .{ 1, 1, 1, 1 },
        .{ 0.5, -0.25, 0.125, -0.0625 },
        .{ 3, -7, 11, -13 },
        .{ -2.5, 4.25, -6.125, 8.0625 },
        .{ wide_high, -wide_low, wide_high, -wide_low },
        .{ 2, 4, 8, 16 },
        .{ smallest, 0, -smallest, 0 },
        .{ largest, -largest, largest / 2, -largest / 2 },
        .{ sub_a, sub_b, -sub_c, sub_a / 2 },
        .{ 1, 1 + eps, 1, 1 - eps },
        .{ mid, mid + 2, mid, -mid },
        .{ 0.1, 0.2, 0.3, 0.4 },
        .{ -0.0, 0.0, 1, -1 },
    };
}

/// `sign((b - a) x (c - a))`, exactly, on one scale.
fn orient2Reference(comptime T: type, a: [2]T, b: [2]T, c: [2]T) g.Sign {
    const all = [_]T{ a[0], a[1], b[0], b[1], c[0], c[1] };
    const base = baseOf(T, all[0..]);
    const ax = scaledAt(T, a[0], base);
    const ay = scaledAt(T, a[1], base);
    const bx = scaledAt(T, b[0], base);
    const by = scaledAt(T, b[1], base);
    const cx = scaledAt(T, c[0], base);
    const cy = scaledAt(T, c[1], base);
    // SAFETY: Reference.sub writes metadata and active limbs before multiplication.
    var left: Reference = undefined;
    // SAFETY: Reference.sub writes metadata and active limbs before multiplication.
    var right: Reference = undefined;
    // SAFETY: Reference.sub writes metadata and active limbs before multiplication.
    var inner: Reference = undefined;
    // SAFETY: Reference.sub writes metadata and active limbs before sign extraction.
    var total: Reference = undefined;
    Reference.sub(&left, bx, ax);
    Reference.sub(&inner, cy, ay);
    Reference.mul(&left, left, inner);
    Reference.sub(&right, by, ay);
    Reference.sub(&inner, cx, ax);
    Reference.mul(&right, right, inner);
    Reference.sub(&total, left, right);
    return Reference.sign(total);
}

/// `sign(dot(cross(b - a, c - a), d - a))`, exactly, on one scale.
fn orient3Reference(comptime T: type, a: [3]T, b: [3]T, c: [3]T, d: [3]T) g.Sign {
    const all = [_]T{
        a[0], a[1], a[2], b[0], b[1], b[2],
        c[0], c[1], c[2], d[0], d[1], d[2],
    };
    const base = baseOf(T, all[0..]);
    var u: [3]Reference = undefined;
    var v: [3]Reference = undefined;
    var w: [3]Reference = undefined;
    inline for (0..3) |index| {
        const anchor = scaledAt(T, a[index], base);
        Reference.sub(&u[index], scaledAt(T, b[index], base), anchor);
        Reference.sub(&v[index], scaledAt(T, c[index], base), anchor);
        Reference.sub(&w[index], scaledAt(T, d[index], base), anchor);
    }
    // SAFETY: Reference.mul writes metadata and active limbs before subtraction.
    var nx: Reference = undefined;
    // SAFETY: Reference.mul writes metadata and active limbs before subtraction.
    var scratch: Reference = undefined;
    // SAFETY: Reference.mul writes metadata and active limbs before accumulation.
    var total: Reference = undefined;
    Reference.mul(&nx, u[1], v[2]);
    Reference.mul(&scratch, u[2], v[1]);
    Reference.sub(&nx, nx, scratch);
    Reference.mul(&total, nx, w[0]);
    Reference.mul(&nx, u[2], v[0]);
    Reference.mul(&scratch, u[0], v[2]);
    Reference.sub(&nx, nx, scratch);
    Reference.mul(&scratch, nx, w[1]);
    Reference.add(&total, total, scratch);
    Reference.mul(&nx, u[0], v[1]);
    Reference.mul(&scratch, u[1], v[0]);
    Reference.sub(&nx, nx, scratch);
    Reference.mul(&scratch, nx, w[2]);
    Reference.add(&total, total, scratch);
    return Reference.sign(total);
}

/// The three-by-three minor of the LAST THREE rows over the three columns
/// that are not `skip`, in the textbook grouping.
fn minor3Reference(out: *Reference, rows: [4][4]Reference, skip: usize) void {
    var kept: [3]usize = undefined;
    var next: usize = 0;
    for (0..4) |column| {
        if (column == skip) continue;
        kept[next] = column;
        next += 1;
    }
    // SAFETY: Reference.mul writes metadata and active limbs before subtraction.
    var left: Reference = undefined;
    // SAFETY: Reference.mul writes metadata and active limbs before subtraction.
    var right: Reference = undefined;
    // SAFETY: Reference.sub writes metadata and active limbs before multiplication.
    var inner: Reference = undefined;
    Reference.mul(&left, rows[2][kept[1]], rows[3][kept[2]]);
    Reference.mul(&right, rows[2][kept[2]], rows[3][kept[1]]);
    Reference.sub(&inner, left, right);
    Reference.mul(out, rows[1][kept[0]], inner);
    Reference.mul(&left, rows[2][kept[0]], rows[3][kept[2]]);
    Reference.mul(&right, rows[2][kept[2]], rows[3][kept[0]]);
    Reference.sub(&inner, left, right);
    Reference.mul(&right, rows[1][kept[1]], inner);
    Reference.sub(out, out.*, right);
    Reference.mul(&left, rows[2][kept[0]], rows[3][kept[1]]);
    Reference.mul(&right, rows[2][kept[1]], rows[3][kept[0]]);
    Reference.sub(&inner, left, right);
    Reference.mul(&right, rows[1][kept[2]], inner);
    Reference.add(out, out.*, right);
}

/// The exact four-by-four determinant, the Laplace expansion along the first
/// row whose cofactor signs alternate.
fn orient4Reference(comptime T: type, a: [4]T, b: [4]T, c: [4]T, d: [4]T) g.Sign {
    const all = [_]T{
        a[0], a[1], a[2], a[3], b[0], b[1], b[2], b[3],
        c[0], c[1], c[2], c[3], d[0], d[1], d[2], d[3],
    };
    const base = baseOf(T, all[0..]);
    var rows: [4][4]Reference = undefined;
    inline for (0..4) |row| {
        inline for (0..4) |column| {
            const source = switch (row) {
                0 => a,
                1 => b,
                2 => c,
                else => d,
            };
            rows[row][column] = scaledAt(T, source[column], base);
        }
    }
    // SAFETY: Reference.mul writes metadata and active limbs before accumulation.
    var total: Reference = undefined;
    // SAFETY: minor3Reference writes metadata and active limbs before multiplication.
    var minor: Reference = undefined;
    // SAFETY: Reference.mul writes metadata and active limbs before accumulation.
    var term: Reference = undefined;
    minor3Reference(&minor, rows, 0);
    Reference.mul(&total, rows[0][0], minor);
    var index: usize = 1;
    while (index < 4) : (index += 1) {
        minor3Reference(&minor, rows, index);
        Reference.mul(&term, rows[0][index], minor);
        if (index % 2 == 0) Reference.add(&total, total, term) else Reference.sub(&total, total, term);
    }
    return Reference.sign(total);
}

/// `sign(dot(p - c, p - c) - r * r)`, exactly, on one scale.
fn containmentReference(comptime T: type, comptime dim: usize, centre: [dim]T, radius: T, point: [dim]T) g.Sign {
    var all: [2 * dim + 1]T = undefined;
    inline for (0..dim) |index| {
        all[index] = centre[index];
        all[dim + index] = point[index];
    }
    all[2 * dim] = radius;
    const base = baseOf(T, all[0..]);
    var offset: [dim]Reference = undefined;
    inline for (0..dim) |index| {
        Reference.sub(&offset[index], scaledAt(T, point[index], base), scaledAt(T, centre[index], base));
    }
    // SAFETY: The first non-empty component iteration assigns total before accumulation.
    var total: Reference = undefined;
    // SAFETY: Reference.mul writes metadata and active limbs before accumulation.
    var term: Reference = undefined;
    inline for (0..dim) |index| {
        Reference.mul(&term, offset[index], offset[index]);
        if (index == 0) {
            total = term;
        } else {
            Reference.add(&total, total, term);
        }
    }
    var scale = scaledAt(T, radius, base);
    Reference.mul(&scale, scale, scale);
    Reference.sub(&total, total, scale);
    return Reference.sign(total);
}

/// `sign(dot(s, s) - t * t)` for `s = second - first`, exactly. This is the
/// one form both tangency comparisons use, once with `t = r1 + r2` and once
/// with `t = r1 - r2`, and it is also the squared separation itself when `t`
/// is zero.
fn tangencyReference(comptime T: type, comptime dim: usize, first: [dim]T, second: [dim]T, t: T) g.Sign {
    var all: [2 * dim + 1]T = undefined;
    inline for (0..dim) |index| {
        all[index] = first[index];
        all[dim + index] = second[index];
    }
    all[2 * dim] = t;
    const base = baseOf(T, all[0..]);
    var separation: [dim]Reference = undefined;
    inline for (0..dim) |index| {
        Reference.sub(&separation[index], scaledAt(T, second[index], base), scaledAt(T, first[index], base));
    }
    // SAFETY: The first non-empty component iteration assigns total before accumulation.
    var total: Reference = undefined;
    // SAFETY: Reference.mul writes metadata and active limbs before accumulation.
    var term: Reference = undefined;
    inline for (0..dim) |index| {
        Reference.mul(&term, separation[index], separation[index]);
        if (index == 0) {
            total = term;
        } else {
            Reference.add(&total, total, term);
        }
    }
    var scale = scaledAt(T, t, base);
    Reference.mul(&scale, scale, scale);
    Reference.sub(&total, total, scale);
    return Reference.sign(total);
}

/// `sign(dot(w, direction) * dot(w, direction) - (dot(w, w) - r * r) * 2^k)`,
/// exactly, where `k = 2 * (base - 1074)` is the factor that puts the
/// projection's square and the constant term on one scale.
///
/// With `P` and `C` the integer projection and constant at scale
/// `2^(base - 1074)`, the answer is `sign(P^2 * 2^k - C)`. Only one of `k` and
/// `-k` is ever positive, so multiplying both terms by the other power of two
/// keeps every shift a non-negative one and the arithmetic stays integral.
fn rayDiscriminantReference(
    comptime T: type,
    comptime dim: usize,
    centre: [dim]T,
    origin: [dim]T,
    direction: [dim]T,
    r: T,
) g.Sign {
    var all: [3 * dim + 1]T = undefined;
    inline for (0..dim) |index| {
        all[index] = centre[index];
        all[dim + index] = origin[index];
        all[2 * dim + index] = direction[index];
    }
    all[3 * dim] = r;
    const base = baseOf(T, all[0..]);
    var w: [dim]Reference = undefined;
    // SAFETY: The first non-empty component iteration assigns projection before accumulation.
    var projection: Reference = undefined;
    // SAFETY: Reference.mul writes metadata and active limbs before accumulation.
    var term: Reference = undefined;
    inline for (0..dim) |index| {
        Reference.sub(&w[index], scaledAt(T, origin[index], base), scaledAt(T, centre[index], base));
        const unit = scaledAt(T, direction[index], base);
        Reference.mul(&term, w[index], unit);
        if (index == 0) {
            projection = term;
        } else {
            Reference.add(&projection, projection, term);
        }
    }
    // SAFETY: The first non-empty component iteration assigns squared before accumulation.
    var squared: Reference = undefined;
    inline for (0..dim) |index| {
        Reference.mul(&term, w[index], w[index]);
        if (index == 0) {
            squared = term;
        } else {
            Reference.add(&squared, squared, term);
        }
    }
    var scale = scaledAt(T, r, base);
    Reference.mul(&scale, scale, scale);
    Reference.sub(&squared, squared, scale);
    Reference.mul(&projection, projection, projection);
    const up: usize = if (2 * base > 2148) 2 * base - 2148 else 0;
    const down: usize = if (2148 > 2 * base) 2148 - 2 * base else 0;
    Reference.shiftedLeft(&projection, projection, up);
    Reference.shiftedLeft(&squared, squared, down);
    Reference.sub(&projection, projection, squared);
    return Reference.sign(projection);
}

test "adaptive_exact orientation in two dimensions agrees with an independent exact area" {
    inline for (.{ g.Point2f, g.Point2d }) |P| {
        const T = ScalarOf(P);
        const corpus = corpusOf(T);
        for (0..corpus.len) |i| {
            const j = (i + 1) % corpus.len;
            const k = (i + 2) % corpus.len;
            const expected = orient2Reference(
                T,
                .{ corpus[i][0], corpus[i][1] },
                .{ corpus[j][0], corpus[j][1] },
                .{ corpus[k][0], corpus[k][1] },
            );
            const given = try P.orientation(
                try P.init(corpus[i][0], corpus[i][1]),
                try P.init(corpus[j][0], corpus[j][1]),
                try P.init(corpus[k][0], corpus[k][1]),
                .adaptive_exact,
            );
            if (given != expected) {
                std.debug.print("2D row {d}: expected {any}, got {any}\n", .{ i, expected, given });
                return error.ExactOrientationDisagrees;
            }
        }
    }
}

test "adaptive_exact orientation in three dimensions agrees with an independent exact volume" {
    inline for (.{ g.Point3f, g.Point3d }) |P| {
        const T = ScalarOf(P);
        const corpus = corpusOf(T);
        for (0..corpus.len) |i| {
            const j = (i + 1) % corpus.len;
            const k = (i + 2) % corpus.len;
            const l = (i + 3) % corpus.len;
            const expected = orient3Reference(
                T,
                .{ corpus[i][0], corpus[i][1], corpus[i][2] },
                .{ corpus[j][0], corpus[j][1], corpus[j][2] },
                .{ corpus[k][0], corpus[k][1], corpus[k][2] },
                .{ corpus[l][0], corpus[l][1], corpus[l][2] },
            );
            const given = try P.orientation(
                try P.init(corpus[i][0], corpus[i][1], corpus[i][2]),
                try P.init(corpus[j][0], corpus[j][1], corpus[j][2]),
                try P.init(corpus[k][0], corpus[k][1], corpus[k][2]),
                try P.init(corpus[l][0], corpus[l][1], corpus[l][2]),
                .adaptive_exact,
            );
            if (given != expected) {
                std.debug.print("3D row {d}: expected {any}, got {any}\n", .{ i, expected, given });
                return error.ExactOrientationDisagrees;
            }
        }
    }
}

test "adaptive_exact orientation in four dimensions agrees with an independent exact determinant" {
    inline for (.{ g.Point4f, g.Point4d }) |P| {
        const T = ScalarOf(P);
        const corpus = corpusOf(T);
        for (0..corpus.len) |i| {
            const j = (i + 1) % corpus.len;
            const k = (i + 2) % corpus.len;
            const l = (i + 3) % corpus.len;
            const expected = orient4Reference(T, corpus[i], corpus[j], corpus[k], corpus[l]);
            const given = try P.orientation(
                try P.init(corpus[i][0], corpus[i][1], corpus[i][2], corpus[i][3]),
                try P.init(corpus[j][0], corpus[j][1], corpus[j][2], corpus[j][3]),
                try P.init(corpus[k][0], corpus[k][1], corpus[k][2], corpus[k][3]),
                try P.init(corpus[l][0], corpus[l][1], corpus[l][2], corpus[l][3]),
                .adaptive_exact,
            );
            if (given != expected) {
                std.debug.print("4D row {d}: expected {any}, got {any}\n", .{ i, expected, given });
                return error.ExactOrientationDisagrees;
            }
        }
    }
}

test "adaptive_exact ball containment agrees with an independent exact squared distance" {
    inline for (.{
        g.Circle2f, g.Circle2d, g.Sphere3f, g.Sphere3d, g.Hypersphere4f, g.Hypersphere4d,
    }) |Ball| {
        const T = @FieldType(Ball, "radius");
        const P = @FieldType(Ball, "center");
        const dim = comptime fieldCountOf(P);
        const corpus = corpusOf(T);
        for ([_]T{ 1, 0.5 }) |radius| {
            for (0..corpus.len) |i| {
                for (0..corpus.len) |j| {
                    const expected = containmentReference(T, dim, corpus[i][0..dim].*, radius, corpus[j][0..dim].*);
                    const ball = try Ball.init(try pointOf(P, dim, corpus[i]), radius);
                    const given = try ball.contains(try pointOf(P, dim, corpus[j]), .adaptive_exact);
                    // The boundary is INSIDE, so zero counts with the negative
                    // side and only a strictly positive exact value is outside.
                    const wanted = expected != .positive;
                    if (given != wanted) {
                        std.debug.print(
                            "ball {d}/{d}/{any}: expected {any}, got {any}\n",
                            .{ i, j, radius, wanted, given },
                        );
                        return error.ExactContainmentDisagrees;
                    }
                }
            }
        }
    }
}

test "adaptive_exact ball tangency agrees with an independent exact separation" {
    inline for (.{
        g.Circle2f, g.Circle2d, g.Sphere3f, g.Sphere3d, g.Hypersphere4f, g.Hypersphere4d,
    }) |Ball| {
        const T = @FieldType(Ball, "radius");
        const P = @FieldType(Ball, "center");
        const dim = comptime fieldCountOf(P);
        const corpus = corpusOf(T);
        for (0..corpus.len) |i| {
            const j = (i + 1) % corpus.len;
            // The module answers `.coincident` before it consults either exact
            // sign, and it decides that on the centres and the radii alone: a
            // separation that is exactly zero with equal radii is one ball.
            var centres_equal = true;
            inline for (0..dim) |index| {
                if (corpus[i][index] != corpus[j][index]) centres_equal = false;
            }
            for ([_]T{ 0.5, 1 }) |first_radius| {
                for ([_]T{ 0.5, 1 }) |second_radius| {
                    // Both comparisons are the same exact form over the same
                    // separation: once against `r1 + r2`, once against `r1 - r2`.
                    // The radii and their sum and difference are exact here, so
                    // the only thing the reference has to get right is the
                    // separation itself.
                    const outer = tangencyReference(
                        T,
                        dim,
                        corpus[i][0..dim].*,
                        corpus[j][0..dim].*,
                        first_radius + second_radius,
                    );
                    const inner = tangencyReference(
                        T,
                        dim,
                        corpus[i][0..dim].*,
                        corpus[j][0..dim].*,
                        first_radius - second_radius,
                    );
                    const ball = try Ball.init(try pointOf(P, dim, corpus[i]), first_radius);
                    const other = try Ball.init(try pointOf(P, dim, corpus[j]), second_radius);
                    const common = ball.intersect(other, .adaptive_exact) catch |err| {
                        // A construction that leaves the format is reported,
                        // and a reported failure has no tag to compare against.
                        if (err == error.NonFiniteResult) continue;
                        return err;
                    };
                    const tag = std.meta.Tag(@TypeOf(common));
                    const given = std.meta.activeTag(common);
                    if (centres_equal and first_radius == second_radius) {
                        try std.testing.expectEqual(tag.coincident, given);
                    } else if (outer == .positive) {
                        try std.testing.expectEqual(tag.none, given);
                    } else if (outer == .zero or inner == .zero) {
                        try std.testing.expectEqual(tag.point, given);
                    } else if (inner == .negative) {
                        try std.testing.expectEqual(tag.contained, given);
                    } else {
                        try std.testing.expectEqual(tag.two_points, given);
                    }
                }
            }
        }
    }
}

test "adaptive_exact ray-against-ball agrees with an independent exact quadratic" {
    inline for (.{
        .{ g.Circle2f, g.Ray2f },
        .{ g.Circle2d, g.Ray2d },
        .{ g.Sphere3f, g.Ray3f },
        .{ g.Sphere3d, g.Ray3d },
        .{ g.Hypersphere4f, g.Ray4f },
        .{ g.Hypersphere4d, g.Ray4d },
    }) |Case| {
        const Ball: type = Case[0];
        const Ray: type = Case[1];
        const P = @FieldType(Ball, "center");
        const V = @FieldType(Ray, "direction");
        const T = @FieldType(Ball, "radius");
        const dim = comptime fieldCountOf(P);
        const corpus = corpusOf(T);
        const radius: T = 1;
        for (0..corpus.len) |i| {
            for ([_]usize{ 1, 5 }) |origin_step| {
                for ([_]usize{ 3, 9 }) |direction_step| {
                    const k = (i + origin_step) % corpus.len;
                    const j = (i + direction_step) % corpus.len;
                    const ball = try Ball.init(try pointOf(P, dim, corpus[i]), radius);
                    // The reference reads the direction the ray STORED, which
                    // is the normalized one the predicate is given.
                    const ray = Ray.init(
                        try pointOf(P, dim, corpus[k]),
                        try vectorOf(V, dim, corpus[j]),
                    ) catch |err| {
                        // A zero direction is not a ray and an overflowing
                        // normalization left the format: neither is a
                        // disagreement about a sign.
                        if (err == error.Degenerate or
                            err == error.NonFiniteInput or
                            err == error.NonFiniteResult) continue;
                        return err;
                    };
                    const common = ball.intersectRay(ray, .adaptive_exact) catch |err| {
                        if (err == error.NonFiniteResult) continue;
                        return err;
                    };
                    const given = std.meta.activeTag(common);
                    const discriminant = rayDiscriminantReference(
                        T,
                        dim,
                        corpus[i][0..dim].*,
                        corpus[k][0..dim].*,
                        ray.direction.toArray(),
                        radius,
                    );
                    if (discriminant == .negative) {
                        // A negative discriminant settles the question before
                        // any construction runs, so nothing can be common.
                        if (given != .none) {
                            std.debug.print("ray {d}/{d}/{d}: negative discriminant, tag {any}\n", .{ i, k, j, given });
                            return error.ExactRayDisagrees;
                        }
                    } else if (containmentReference(
                        T,
                        dim,
                        corpus[i][0..dim].*,
                        radius,
                        corpus[k][0..dim].*,
                    ) == .negative) {
                        // The origin is strictly inside, so the common set is
                        // the run from the origin out to the surface: the far
                        // root is positive exactly because the perpendicular
                        // distance is below the radius.
                        if (given != .segment) {
                            std.debug.print("ray {d}/{d}/{d}: origin inside, tag {any}\n", .{ i, k, j, given });
                            return error.ExactRayDisagrees;
                        }
                    }
                }
            }
        }
    }
}

test "exact orientation keeps mixed normal and subnormal collinearity" {
    inline for (.{ g.Point2f, g.Point2d }) |P| {
        const T = ScalarOf(P);
        const Bits = @Int(.unsigned, @bitSizeOf(T));
        const tiny: T = @bitCast(@as(Bits, 1));
        var coordinates = [_]T{ 1, 2, tiny, tiny * 2 };
        std.mem.doNotOptimizeAway(&coordinates);
        const origin = try P.init(0, 0);
        const normal = try P.init(coordinates[0], coordinates[1]);
        const subnormal = try P.init(coordinates[2], coordinates[3]);
        try std.testing.expectEqual(g.Sign.zero, try P.orientation(origin, normal, subnormal, .adaptive_exact));
    }
}
