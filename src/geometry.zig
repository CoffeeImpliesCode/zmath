//! Concrete `f32`/`f64` families in two, three, and four dimensions.
//!
//! The "Geometry" section of `API.md` is the whole contract. This file spells
//! it out.
//!
//! ## Naming
//!
//! Every public type carries its dimension and its precision, in that order:
//! `Point2f`, `Vector3d`, `Line4f`, `Aabb2d`, `Hypersphere4f`. The `f`/`d`
//! suffix is the project's convention (`scalar.hypotf`, `scalar.hypotd`), and
//! the digit is the number of components. A two-component `f32` point is
//! `Point2f`; the same object in binary64 is `Point2d`. Nothing in the public
//! surface is generic.
//!
//! ## Using it
//!
//! The shortest sequence that reaches every part of the contract. It is a
//! consumer fragment: it names the package module, the concrete types, and
//! the errors, and it compiles as written in any Zig source file that can
//! `@import("zmath")`.
//!
//! ```zig
//! const std = @import("std");
//! const geometry = @import("zmath").geometry;
//!
//! fn cut(centre: geometry.Point3f, radius: f32, through: geometry.Point3f) !void {
//!     // Constructors reject non-finite fields and a negative radius.
//!     // Plane3f.init also normalizes its nonzero normal and can return
//!     // error.Degenerate or error.NonFiniteResult.
//!     const ball = try geometry.Sphere3f.init(centre, radius);
//!     const plane = try geometry.Plane3f.init(
//!         through,
//!         try geometry.Vector3f.init(0, 0, 1),
//!     );
//!
//!     // `error.NonFiniteInput` for a non-finite field anywhere in either
//!     // receiver or in the coordinate, and `error.NonFiniteResult` when
//!     // finite values make the arithmetic leave the format. The policy is
//!     // the caller's to pick: `.adaptive_exact` answers on integers and
//!     // `.fast` reports the plain evaluation.
//!     const common = try ball.intersectPlane(plane, .adaptive_exact);
//!     switch (common) {
//!         .point => |touch| std.debug.print("tangent at {any}\n", .{touch}),
//!         .circle => |section| {
//!             // `section` is a centre and a radius read INSIDE `plane`. It
//!             // carries no direction, so `plane` is what places this curve
//!             // in space.
//!             std.debug.print("section at {any} radius {d}\n", .{ section.center, section.radius });
//!         },
//!         .none => std.debug.print("the plane misses the ball\n", .{}),
//!         else => std.debug.print("an answer this pair does not produce\n", .{}),
//!     }
//! }
//! ```
//!
//! ## Finiteness
//!
//! Checked arithmetic validates every supplied field, including the whole
//! receiver, and returns `error.NonFiniteInput` for an infinity or NaN.
//! Fields are writable; a constructor's check does not validate later writes.
//! `isFinite`, `eq`, `toArray`, and the infallible `Bivector4.neg` do not
//! perform that validation. They can inspect or preserve non-finite fields.
//!
//! Checked numeric operations return `error.NonFiniteResult` for a non-finite
//! computed value or required construction intermediate. They never publish
//! an infinity, NaN, or zero direction produced by dividing by an infinite
//! norm. Sign predicates have the separate policy below.
//!
//! The two tags therefore answer different questions.
//! `error.NonFiniteInput` says the caller's own number was already not a
//! number. `error.NonFiniteResult` says the caller's numbers were all ordinary
//! and the arithmetic itself ran off the end of the format, which is the case
//! a caller scales geometry for and cannot predict from the values it holds.
//!
//! A sign is a different kind of result. `orientation` and every predicate
//! behind an intersection answer with a `Sign`, never with a float, so an
//! arithmetic that leaves the format has nothing to leak in the first place.
//! `adaptive_exact` reads the exact value and forms no floating-point
//! intermediate at all, so it answers for any finite inputs. `fast` reports
//! the sign of the plain scalar evaluation exactly as this module defines it,
//! and an INFINITY keeps the sign it has. A NaN is the one case it cannot
//! answer: every comparison with a NaN is false, so reporting `zero` for it
//! would claim the value is exactly zero, and `error.NonFiniteResult` says
//! instead that the plain evaluation left the format. The two policies
//! therefore differ in failure as well as in sign, and a caller that sees the
//! tag from `fast` can retry the same inputs under `adaptive_exact`.
//!
//! Each constructor adds the one failure mode its own contract names:
//!
//! - `Line`, `Ray`, `Plane`, `Hyperplane`, `Vector.normalize`:
//!   `error{ NonFiniteInput, NonFiniteResult, Degenerate }`. The direction or
//!   normal is normalized, and normalizing is arithmetic.
//! - `Circle`, `Sphere`, `Hypersphere`: `error{ NonFiniteInput,
//!   NegativeRadius }`. They store a centre and a radius and compute nothing.
//! - `Aabb`: `error{ NonFiniteInput, InvertedBounds }`. They compare two
//!   corners and compute nothing either.
//!
//! ## A point and a vector are different types
//!
//! `Point2f` and `Vector2f` are distinct structs. `Point2f.add` takes a
//! `Vector2f`; `Point2f.toVector` converts. `Line2f.init` takes a point and a
//! vector because a point and a direction are what a line is made of, and
//! neither type is accepted where the other is named.
//!
//! ## Predicate policy
//!
//! Every predicate that can DIFFER between policies takes the policy as an
//! argument, and nothing picks one for the caller. `fast` reports the sign of
//! the plain scalar evaluation, in the operation order this module documents,
//! with contraction disabled, and returns `error.NonFiniteResult` when that
//! evaluation is a NaN. `adaptive_exact` reports the exact sign for finite
//! `f32` and `f64` inputs. There is no epsilon anywhere: `adaptive_exact` is
//! exact, not tighter, and it never falls back to the `fast` answer, nor the
//! other way round.
//!
//! Predicates that CANNOT differ take no policy, because a policy argument a
//! predicate ignores is a lie. An `Aabb` contains an axis-aligned range
//! comparison, which rounds nothing, so `Aabb.contains`, `Aabb.intersects`,
//! and `Aabb.intersect` take no policy. Every polynomial predicate -
//! orientation, sphere containment, sphere to sphere, sphere to flat, and
//! every affine intersection that reduces to them - takes one.
//!
//! `adaptive_exact` is exact for the values handed to it. The three- and
//! four-dimensional affine intersections first decide coplanarity with an
//! exact predicate on the ORIGINAL coordinates, then project into a plane and
//! run the exact two-dimensional predicate on the PROJECTED coordinates. The
//! coplanarity gate, which is the decision a plain computation gets wrong, is
//! therefore exact on the inputs as given; the in-plane decisions are exact
//! for the frame the projection produced. The intersection PARAMETERS are
//! `f32`/`f64` values and the range tests on them are ordinary comparisons,
//! so a parameter or a projected coordinate that left the format is
//! `error.NonFiniteResult` rather than a comparison against a bound that
//! means nothing.
//!
//! ## Exact arithmetic
//!
//! Every finite `f64` is `m * 2^e` with `|m| < 2^53` and
//! `-1074 <= e <= 971`, so `x * 2^1074` is an integer below `2^2098`. An `f32`
//! widens to `f64` exactly, so one representation covers both precisions.
//!
//! A difference is NEVER rounded on the way in. Two stored coordinates are
//! converted and subtracted inside the exact arithmetic, because a
//! floating-point subtraction is where the exact policy used to lose the
//! sign, and on a pair of large opposite coordinates it overflows to an
//! infinity that is not a real number at all. Both operands carry the same
//! POSITIVE `2^1074`, so their difference carries it too and needs no
//! rescaling. A difference of two finite values is below `2 * floatMax <
//! 2^1025`, so the scaled difference is below `2^2099`: one binary digit more
//! than a stored coordinate, which is the headroom the limb bound below
//! carries.
//!
//! Every term of a predicate has the same degree, so every term of its integer
//! value carries the same power of two, and that power is positive. The sign
//! of the integer is the sign of the real value, and nothing rounded.
//!
//! Storage is a fixed `[132]u64` magnitude plus a sign, on the stack. 132
//! limbs is 8448 bits. The worst predicate here is the ray-against-ball
//! discriminant, and the ball-against-flat distance beside it. A projection
//! of exact differences onto exact coordinates is below
//! `4 * 2^1025 * 2^1024 = 2^2051`, so on the `2^2148` scale its square is
//! below `2^8398`; the constant term it is compared against is below
//! `2^4201` on the `2^2148` scale and is scaled up by that same positive
//! `2^2148`, leaving the whole value below `2^8399`. The degree-four
//! four-by-four determinant, whose 24 terms are products of four scaled
//! coordinates, is below `24 * 2^8392 < 2^8397`, and every degree-three and
//! degree-two form is smaller again. That leaves 49 bits of headroom in the
//! 8448 the storage holds. Each predicate converts every coordinate it reads
//! exactly once, into a fixed array of these, so the scratch one needs is a
//! compile-time constant number of them. Nothing allocates, and no limb
//! count depends on a value.
//!
//! ## A limb at or above `used` is never read
//!
//! Only the first `used` limbs of a `Big` carry its value, and nothing in this
//! file reads above that; the audit that establishes it is on the `Big`
//! declaration. A predicate therefore keeps ONE working set of exact values and
//! has every operation write into it, rather than having each operation build
//! a fresh 132-limb array and hand 1064 bytes of it back by value. Two rules
//! follow from the invariant, and each is stated where it bites: an operation
//! that READS its destination before writing it clears exactly the region its
//! closing scan reads, and an operation that writes every index it scans needs
//! no clear at all. Either way no operation above the one it was handed is ever
//! read, so a reused working set cannot carry a stale limb into a result.
//!
//! ## Intersection unions
//!
//! One union per dimension and precision, `Intersection2f` through
//! `Intersection4d`. Each is a tagged union; none of them uses `null` to
//! combine cases, and `none` is an ordinary tag like any other. The tags name
//! the DIMENSION of the intersection's affine span - zero `point`, one `line`,
//! two `plane`, three `hyperplane` - with the bounded and solid forms named
//! separately, plus three tags the contract's nine do not cover:
//!
//! - `box`: an AABB. The intersection of two boxes is a box, and none of the
//!   nine names it. A box is a set of points like any other; the tag is here
//!   because the SHAPE of the answer is not one the nine name.
//! - `circle` in three dimensions and `hypersphere` in four: a flat that cuts
//!   a ball gives back the cut, which is a curve rather than a single set of
//!   isolated points. The payload is a centre and a radius read INSIDE the
//!   flat, so it does not carry the flat itself; a caller that needs the cut
//!   as a set in space keeps the flat it passed in.
//!
//! `coincident` means the two operands are the same set and that set is not an
//! infinite flat, so it covers identical segments, rays, spheres, and boxes.
//! `contained` means the intersection is exactly one operand and the operands
//! are not the same set. A coincident pair of lines, planes, or hyperplanes
//! returns the tag that names its own span instead, because "the intersection
//! is this line" is the more useful answer than "the two are the same". A RAY
//! is not a line, and an identical pair of rays is a coincident pair of
//! half-lines: `line` is reserved for a common set that runs both ways.
//!
//! Every tag names the WHOLE common set, never a part of it, with two
//! exceptions the union cannot express and the declaration on the union says
//! so. The first settles the cases where a part would have been convenient:
//!
//! - A ray whose origin is inside a ball meets it in the whole run from the
//!   origin to the exit, so the tag is `segment`, not the exit alone. The
//!   origin is a common point, so reporting only the exit would report a set
//!   that is not the intersection. The run is a `point` when the origin is
//!   already on the surface and the ray points away, because then it has
//!   zero length. A ray whose origin is outside also reports the closed
//!   `segment` from entry to exit, or `point` when those endpoints coincide.
//! - A segment that merely ENDS on a ball or a flat is not a touch. A
//!   segment that reaches the surface from outside crosses into the solid and
//!   its common set is a proper sub-segment whose far end is that surface
//!   point. A true touch is a DOUBLE root, and only then is the common set
//!   the single tangent point.
//!
//! The second is what a pair of solid balls produces. Two closed balls that
//! cross share a lens, and two_points names where their BOUNDARIES cross
//! instead; a closed ball cut by a flat shares the disc its circle bounds, and
//! the `circle` and `hypersphere` tags name the bounding curve. Both are
//! boundary answers, stated here and on the union declarations rather than
//! hidden, because naming them any other way would be a claim the union does
//! not keep.
//!
//! A ball is the CLOSED solid ball, so its boundary belongs to it and a
//! segment ending exactly on the surface is inside up to that point.
//!
//! ## Deliberate omissions
//!
//! - `Vector.cross` exists in three dimensions only. In four dimensions there
//!   is no cross product, because the orthogonal complement of the plane two
//!   vectors span is a whole plane rather than a line. `Vector4.wedge`
//!   provides the exterior product `u ^ v`, which is the product that
//!   dimension actually has.
//! - `Point2.side` exists in two dimensions only. A line divides the plane
//!   into two sides; a line in three or four dimensions does not.
//! - `Ray` has no `closestPointTo`. The contract names `at` for a ray and
//!   `closestPointTo` for a segment, and that asymmetry is kept rather than
//!   smoothed over.
//!
//! `ContainedSide` names the operand the intersection EQUALS, and it names it
//! by its position in the call, so the receiver of a method is `first` and
//! its argument is `second`. A box that holds its argument has the ARGUMENT
//! as its whole common set, so it reports `contained = .second`, exactly as
//! a segment lying inside the ball or the flat it is tested against does. It
//! is `.first` only when the receiver is the one wholly inside, which for
//! `unit.intersect(inner)` means the reversed call `inner.intersect(unit)`.
const std = @import("std");
const scalar = @import("scalar.zig");

// ---------------------------------------------------------------------------
// Policy, sign, and the shared scalar kernels
// ---------------------------------------------------------------------------

/// Selects the sign evaluation for a predicate.
///
/// `fast` uses the rounded scalar polynomial. Infinity retains its sign;
/// NaN returns `error.NonFiniteResult`. `adaptive_exact` evaluates the
/// polynomial exactly on finite inputs. It can resolve a sign when `fast`
/// rounds to zero or fails. Floating constructions can still fail under
/// either policy.
pub const PredicatePolicy = enum { fast, adaptive_exact };

/// The sign of a real value.
pub const Sign = enum { negative, zero, positive };

/// Which operand an intersection equals, when the intersection is one whole
/// operand and the two operands are not the same set.
pub const ContainedSide = enum { first, second };

fn hypotOf(comptime T: type, a: T, b: T) T {
    return switch (T) {
        f32 => scalar.hypotf(a, b),
        f64 => scalar.hypotd(a, b),
        else => @compileError("zmath geometry is f32 or f64"),
    };
}

fn sqrtOf(comptime T: type, a: T) T {
    return switch (T) {
        f32 => scalar.sqrtf(a),
        f64 => scalar.sqrtd(a),
        else => @compileError("zmath geometry is f32 or f64"),
    };
}

fn minOf(comptime T: type, a: T, b: T) T {
    return switch (T) {
        f32 => scalar.fminf(a, b),
        f64 => scalar.fmind(a, b),
        else => @compileError("zmath geometry is f32 or f64"),
    };
}

fn maxOf(comptime T: type, a: T, b: T) T {
    return switch (T) {
        f32 => scalar.fmaxf(a, b),
        f64 => scalar.fmaxd(a, b),
        else => @compileError("zmath geometry is f32 or f64"),
    };
}

fn absOf(comptime T: type, x: T) T {
    return switch (T) {
        f32 => scalar.absf(x),
        f64 => scalar.absd(x),
        else => @compileError("zmath geometry is f32 or f64"),
    };
}

fn isFiniteOf(comptime T: type, x: T) bool {
    return switch (T) {
        f32 => std.math.isFinite(x),
        f64 => std.math.isFinite(x),
        else => @compileError("zmath geometry is f32 or f64"),
    };
}

/// The sign of a finite plain evaluation, and of an infinite one, which keeps
/// the sign it has. A NaN has no sign, so `fast` cannot answer from it and the
/// caller is told so instead of being handed a `zero` that reads as "exactly
/// zero". IEEE makes that the natural split: `inf > 0` and `inf < 0` are
/// ordinary comparisons that decide, and every comparison with a NaN is false,
/// so a `signOfScalar` that answered `zero` for a NaN would be reporting a
/// value the arithmetic never produced.
fn fastSign(comptime T: type, value: T) error{NonFiniteResult}!Sign {
    if (std.math.isNan(value)) return error.NonFiniteResult;
    if (value > 0) return .positive;
    if (value < 0) return .negative;
    return .zero;
}

/// The gate a computed value passes before it is returned or fed to a
/// construction. Supplied finite values can leave the format inside this
/// module, and a geometry call reports that instead of handing the caller an
/// infinity, a NaN, or a value silently normalized to zero because an
/// intermediate overflowed.
fn requireFiniteResult(comptime T: type, value: T) error{NonFiniteResult}!T {
    if (!isFiniteOf(T, value)) return error.NonFiniteResult;
    return value;
}

/// `x` widened to `f64`, which is exact for every finite `f32`.
fn widen(comptime T: type, x: T) f64 {
    return @floatCast(x);
}

fn clampRange(comptime T: type, value: T, lower: T, upper: T) T {
    if (value < lower) return lower;
    if (value > upper) return upper;
    return value;
}

fn clampUnitOf(comptime T: type, t: T) T {
    return clampRange(T, t, 0, 1);
}

/// The open end of the parameter range of a line or a ray. It is a range
/// bound only, never an input: nothing here validates it, because the
/// comparisons that use it accept an infinite bound by construction. The
/// products that reach it are likewise left alone: an infinite end of a
/// range is the geometry, not an overflow.
fn infinityOf(comptime T: type) T {
    return std.math.inf(T);
}

// ---------------------------------------------------------------------------
// Field-name plumbing
//
// Every value type is a struct of named components, so the shared algorithms
// walk the fields by name. The names and the order are the ones the value
// layer in `vector_value.zig` already uses, so the two layers cannot drift.
// ---------------------------------------------------------------------------

const component_names = [4][]const u8{ "x", "y", "z", "w" };

fn fieldCount(comptime V: type) usize {
    return @typeInfo(V).@"struct".field_names.len;
}

fn fieldName(comptime index: usize) []const u8 {
    return component_names[index];
}

fn componentType(comptime V: type) type {
    return @typeInfo(V).@"struct".field_types[0];
}

fn toArrayOf(comptime V: type, a: V) [fieldCount(V)]componentType(V) {
    const width = comptime fieldCount(V);
    var out: [width]componentType(V) = undefined;
    inline for (0..width) |index| {
        out[index] = @field(a, fieldName(index));
    }
    return out;
}

fn fromArrayOf(comptime V: type, a: [fieldCount(V)]componentType(V)) V {
    const width = comptime fieldCount(V);
    var out: V = undefined;
    inline for (0..width) |index| {
        @field(out, fieldName(index)) = a[index];
    }
    return out;
}

/// The input gate. It walks the value's OWN fields, so it covers a point, a
/// vector, and every composite - a line, a ray, a segment, a ball, a flat, a
/// box - without the caller naming the field's type at the call site. A
/// composite's fields are not all named after an axis, so the axis-named walk
/// that the point and vector types use cannot be shared with them.
///
/// It covers EVERY field, not the ones the calling operation reads, because the
/// fields are public and writable: an operation that validated half a
/// receiver would let a poisoned field reach the next one.
fn requireFiniteField(value: anytype) error{NonFiniteInput}!void {
    const T = @TypeOf(value);
    if (T == f32 or T == f64) {
        if (!std.math.isFinite(value)) return error.NonFiniteInput;
        return;
    }
    const width = comptime @typeInfo(T).@"struct".field_names.len;
    inline for (0..width) |index| {
        const name = comptime @typeInfo(T).@"struct".field_names[index];
        try requireFiniteField(@field(value, name));
    }
}

/// The result gate over a row of components, which is the shape every shared
/// algorithm carries between its steps.
fn requireFiniteRow(comptime dim: usize, comptime T: type, a: [dim]T) error{NonFiniteResult}![dim]T {
    inline for (0..dim) |index| {
        if (!isFiniteOf(T, a[index])) return error.NonFiniteResult;
    }
    return a;
}

/// The same gate over a point or a vector, reached by the axis names every
/// such value already carries.
fn requireFiniteComponents(comptime V: type, a: V) error{NonFiniteResult}!V {
    return fromArrayOf(V, try requireFiniteRow(fieldCount(V), componentType(V), toArrayOf(V, a)));
}

/// The input gate for one scalar.
fn requireFiniteScalar(comptime T: type, a: T) error{NonFiniteInput}!void {
    if (!isFiniteOf(T, a)) return error.NonFiniteInput;
}

fn isFiniteQuiet(comptime V: type, a: V) bool {
    const width = comptime fieldCount(V);
    var finite = true;
    inline for (0..width) |index| {
        finite = finite and isFiniteOf(componentType(V), @field(a, fieldName(index)));
    }
    return finite;
}

/// Bit equality with the NaN normalization the value layer uses: two NaNs
/// match, `+0` and `-0` do not.
fn eqOf(comptime V: type, a: V, b: V) bool {
    const T = componentType(V);
    const width = comptime fieldCount(V);
    var same = true;
    inline for (0..width) |index| {
        const left = @field(a, fieldName(index));
        const right = @field(b, fieldName(index));
        const left_nan = std.math.isNan(left);
        if (left_nan != std.math.isNan(right)) {
            same = false;
        } else if (left_nan) {
            same = same;
        } else {
            same = same and switch (T) {
                f32 => @as(u32, @bitCast(left)) == @as(u32, @bitCast(right)),
                f64 => @as(u64, @bitCast(left)) == @as(u64, @bitCast(right)),
                else => @compileError("zmath geometry is f32 or f64"),
            };
        }
    }
    return same;
}

// ---------------------------------------------------------------------------
// Exact sign arithmetic
// ---------------------------------------------------------------------------

/// Limbs of exact magnitude. The module comment proves 132 is enough for the
/// worst case, which is the ray-against-ball discriminant and not the
/// degree-four form.
const exact_limbs = 132;

/// The bits that lift a degree-two exact value onto the scale of the SQUARE
/// of a degree-one one. `fromFloat` scales every stored coordinate and every
/// exact difference by `2^1074`, so a projection of two of them carries
/// `2^2148` and its square carries `2^4296`, while a product of two `2^1074`
/// values carries `2^2148`. `Big` holds integers rather than rationals, so
/// such a comparison scales the SMALLER side up by this much instead of
/// taking the larger side down, and the sign is the same either way.
const projection_scale: u32 = 2148;

/// An exact signed integer, used only to decide a sign. A degree-`k`
/// polynomial's value here is the real value times `2^(1074 * k)`, a positive
/// factor, so the sign of the integer is the sign of the real value.
/// ## A limb at or above `used` is never read
///
/// `limbs` is a fixed 132-limb array, but only its first `used` limbs carry
/// the value, and the array is NOT required to be zero anywhere else. That is
/// what lets each producer below clear only the region its own later reads
/// touch instead of the whole array, and it is AUDITED, not assumed: every
/// `.limbs[` read in this file is inside `0..used` of the value it reads.
///
/// The audit, read by read:
///
/// - `mul` reads `a.limbs[i]` under `i < a.used` and `b.limbs[j]` under
///   `j < @min(b.used, exact_limbs - i)`, so `j < b.used`.
/// - `shiftedLeft` reads `a.limbs[index]` under `index < a.used`.
/// - `magnitudeCompare` reads `a.limbs[index]` and `b.limbs[index]` under
///   `index < a.used`, and it returns early unless `a.used == b.used`, so
///   `index < b.used` as well.
/// - `magnitudeSum` reads `a.limbs[index]` and `b.limbs[index]` under the
///   guards `index < a.used` and `index < b.used`.
/// - `magnitudeDifference` reads `a.limbs[index]` under `index < a.used` and
///   `b.limbs[index]` under the guard `index < b.used`.
///
/// There is no other read. `topUsed` is the one exception and it is not one:
/// it reads `0..limit` of a producer's OWN destination, and the producer is
/// the code that established `limit` as the top of the region it filled. Every
/// other consumer of a `Big` sees it only through `0..used`, which the scan
/// itself bounds: a downward scan that stops at `limit` returns an index whose
/// limbs above it, inside `0..limit`, are all zero, and nothing above `limit`
/// is read at all. A partially initialised array therefore cannot escape into
/// a result.
const Big = struct {
    negative: bool,
    used: u8,
    limbs: [exact_limbs]u64,

    /// The value zero, written into a destination that may still hold a
    /// previous answer. Only the two fields a reader looks at are stored:
    /// `sign` reports zero from `used`, `addInto` and `subInto` short-circuit
    /// on `used`, and `mulInto` and `shiftedLeftInto` return before they read
    /// a limb, so a magnitude under a zero `used` is never read and clearing
    /// it would be a write nobody reads.
    fn setZero(out: *Big) void {
        out.negative = false;
        out.used = 0;
    }

    fn sign(a: *const Big) Sign {
        if (a.used == 0) return .zero;
        return if (a.negative) .negative else .positive;
    }

    /// The exact integer `value * 2^1074`, written straight into `out`. The
    /// caller has already rejected a non-finite input, which is the only thing
    /// that reaches here.
    ///
    /// Clear the low prefix before placing the significand. Consumers read
    /// every limb below `used`, including the implicit zeros below `word`.
    /// The pair stores then define the remaining active limbs. The tail above
    /// `used` is not part of the value and needs no initialization.
    fn fromFloat(out: *Big, value: f64) void {
        const bits: u64 = @bitCast(value);
        const exponent_field: u64 = (bits >> 52) & 0x7ff;
        const fraction = bits & 0x000f_ffff_ffff_ffff;
        if (exponent_field == 0x7ff) return setZero(out);
        if (exponent_field == 0 and fraction == 0) return setZero(out);
        // A normal value is `significand * 2^(E - 1075)` and a subnormal is
        // `fraction * 2^-1074`; scaling by 2^1074 makes the shift a
        // non-negative limb offset.
        const significand = fraction | (if (exponent_field == 0) @as(u64, 0) else (@as(u64, 1) << 52));
        const offset: u32 = if (exponent_field == 0) 0 else @intCast(exponent_field - 1);
        const word: usize = offset / 64;
        const shift: u6 = @intCast(offset % 64);
        // A 128-bit shift writes the two limbs at once, which avoids a
        // shift-by-64 that the six-bit shift count cannot even name.
        const placed: u128 = @as(u128, significand) << shift;
        @memset(out.limbs[0..word], 0);
        out.limbs[word] = @truncate(placed);
        out.limbs[word + 1] = @truncate(placed >> 64);
        out.negative = (bits >> 63) != 0;
        // `placed` is not zero, so at least one of the two limbs just
        // written is nonzero, and nothing above them was touched: the
        // count is the index of the upper one plus one. The upper limb
        // is nonzero exactly when the shifted significand crosses the
        // 64-bit boundary.
        out.used = @intCast(word + 1 + @intFromBool(out.limbs[word + 1] != 0));
    }

    /// Copy only initialized magnitude limbs. Identical operands need no copy.
    /// The unused tail remains outside every arithmetic read.
    fn copyInto(out: *Big, a: *const Big) void {
        if (out == a) return;
        @memcpy(out.limbs[0..a.used], a.limbs[0..a.used]);
        out.used = a.used;
        out.negative = a.negative;
    }

    /// `-b`, written into `out`. Only the sign field changes, so the
    /// destination may be `b` itself. This is the one place a negation is
    /// written out; `subInto` below folds the flip into the sign it passes
    /// down rather than building a value to negate.
    fn negateInto(out: *Big, b: *const Big) void {
        copyInto(out, b);
        out.negative = !b.negative;
    }

    /// `a + b`, written into `out`. `out` MAY BE `a` OR `b`: every path below
    /// reads an index before it writes that same index and never comes back
    /// to it, and the two scalar fields are read before any store reaches
    /// them.
    fn addInto(out: *Big, a: *const Big, b: *const Big) void {
        if (a.used == 0) {
            copyInto(out, b);
            return;
        }
        if (b.used == 0) {
            copyInto(out, a);
            return;
        }
        if (a.negative == b.negative) return magnitudeSum(out, a, b, a.negative);
        const order = magnitudeCompare(a, b);
        if (order == .eq) return setZero(out);
        if (order == .gt) return magnitudeDifference(out, a, b, a.negative);
        return magnitudeDifference(out, b, a, b.negative);
    }

    /// `a - b`, written into `out`. This is `a + (-b)`, and a negation flips
    /// one field rather than building a value, so the branches are chosen from
    /// `b`'s sign directly instead of from a negated copy. `out` MAY BE `a` OR
    /// `b`, for the reason `addInto` gives.
    fn subInto(out: *Big, a: *const Big, b: *const Big) void {
        if (a.used == 0) return negateInto(out, b);
        if (b.used == 0) {
            copyInto(out, a);
            return;
        }
        if (a.negative != b.negative) return magnitudeSum(out, a, b, a.negative);
        const order = magnitudeCompare(a, b);
        if (order == .eq) return setZero(out);
        if (order == .gt) return magnitudeDifference(out, a, b, a.negative);
        return magnitudeDifference(out, b, a, !b.negative);
    }

    /// Schoolbook, truncated at the array, written into `out`. The bounds in
    /// the module comment prove the product fits, so the truncation is a bound
    /// check and not a loss.
    ///
    /// `out` MUST NOT BE `a` OR `b`: the loop reads `a.limbs[i]` at an index
    /// that climbs while it writes `out.limbs[i + j]` above it.
    ///
    /// It clears exactly the region it goes on to read. The loop READS
    /// `out.limbs[i + j]` before writing it, and `scan` below - the same
    /// `usedBound` expression the `used` line uses - is above every index the
    /// loop touches: `i < a.used` and `j < @min(b.used, exact_limbs - i)` give
    /// `i + j < i + row <= min(b.used + i, exact_limbs) <= scan`, and the
    /// closing add is at `i + row`, which the same line covers. Nothing above
    /// `scan` is read, so nothing above `scan` is cleared, and the cost tracks
    /// the width of the data rather than the width of the storage.
    ///
    /// A row whose LEFT limb is zero is a no-op, and skipping it is what makes
    /// an ordinary coordinate cheap to multiply: `fromFloat` writes its pair at
    /// limb `(exponent_field - 1) / 64` for a normal value and at limb zero for
    /// a subnormal, so in both cases every limb below the written pair is zero
    /// and only the top one or two limbs carry anything. The proof that
    /// the skipped row changed nothing is that `carry` enters the row at zero,
    /// so each inner total is `0 * b.limbs[j] + limbs[i + j] + 0`, which
    /// writes back what it read and leaves the carry at zero, and the closing
    /// add is then `limbs[i + row] + 0`. Every write in the row is the value it
    /// already held. The RIGHT operand gets no such treatment: there a zero
    /// limb still has to absorb the running carry, which is not zero.
    fn mulInto(out: *Big, a: *const Big, b: *const Big) void {
        if (a.used == 0 or b.used == 0) return setZero(out);
        const scan = usedBound(@as(usize, a.used) + @as(usize, b.used));
        @memset(out.limbs[0..scan], 0);
        var i: usize = 0;
        while (i < a.used and i < exact_limbs) : (i += 1) {
            if (a.limbs[i] == 0) continue;
            const row = @min(b.used, exact_limbs - i);
            var carry: u128 = 0;
            var j: usize = 0;
            while (j < row) : (j += 1) {
                const total: u128 = @as(u128, a.limbs[i]) * @as(u128, b.limbs[j]) +
                    @as(u128, out.limbs[i + j]) + carry;
                out.limbs[i + j] = @truncate(total);
                carry = total >> 64;
            }
            if (i + row < exact_limbs) {
                const total: u128 = @as(u128, out.limbs[i + row]) + carry;
                out.limbs[i + row] = @truncate(total);
            }
        }
        out.negative = a.negative != b.negative;
        // The row `i` writes only `i + j` for `j < row` and then `i + row`,
        // and `i + row` is at most `min(b.used + i, exact_limbs)`, so
        // nothing above `a.used + b.used` was touched. The sum is widened
        // before it is added, because two full-width operands overflow a
        // byte.
        out.used = topUsed(&out.limbs, scan);
    }

    /// `a` scaled up by `2^bits`, the sign-preserving rescale that puts two
    /// exact products on one scale. Whole limbs move as whole limbs, so this
    /// is a limb pass rather than a bit loop, and the bounds in the module
    /// comment prove the result fits.
    ///
    /// `out` MUST NOT BE `a`: the `|=` reads its target, and the clear below
    /// would take the operand's own limbs first. The clear is exactly the
    /// region that is read: the source limb `index` lands at `index + words`
    /// and `index + words + 1`, both below `a.used + words + 1`, the loop
    /// stops as soon as the first of those leaves the array, and `scan` below
    /// is the same `usedBound` expression the `used` line uses. Nothing at or
    /// above `scan` is read by the `|=` or by the scan, so nothing at or
    /// above it is cleared.
    fn shiftedLeftInto(out: *Big, a: *const Big, bits: u32) void {
        if (a.used == 0 or bits == 0) {
            copyInto(out, a);
            return;
        }
        const words: usize = @intCast(bits / 64);
        const shift: u6 = @intCast(bits % 64);
        const scan = usedBound(@as(usize, a.used) + words + 1);
        @memset(out.limbs[0..scan], 0);
        var index: usize = 0;
        while (index < a.used) : (index += 1) {
            const target = index + words;
            if (target >= exact_limbs) break;
            // A 128-bit shift places both halves at once, and the high half
            // of one limb is the low half of the next, so it accumulates.
            const placed: u128 = @as(u128, a.limbs[index]) << shift;
            out.limbs[target] |= @truncate(placed);
            if (target + 1 < exact_limbs) out.limbs[target + 1] |= @as(u64, @truncate(placed >> 64));
        }
        out.negative = a.negative;
        // The source limb `index` lands at `index + words` and at
        // `index + words + 1`, and the loop stops as soon as the first of
        // those leaves the array, so nothing above `a.used + words + 1` was
        // touched.
        out.used = topUsed(&out.limbs, scan);
    }
};

/// The number of limbs in use, found by scanning DOWN from `limit`.
///
/// `limit` is the bound the producer that calls this one established, and the
/// scan reads `0..limit` and nothing above it. That is the whole contract: the
/// destination is a REUSED value, so the array is not zero to begin with, and
/// each producer fills exactly `0..limit` - by assignment where it writes
/// every index, or by an explicit clear where it reads a destination before
/// writing it - before calling this. A scan that starts at the bound then
/// finds exactly the index a scan of the whole array would, because every limb
/// between that index and the bound is zero and every limb above the bound is
/// never read.
///
/// The bound is what makes the cost track the WIDTH OF THE DATA rather than
/// the width of the storage: a coordinate converted from a stored double
/// occupies at most a third of the array, and a two-by-two determinant pays
/// nothing like a full pass per operation.
fn topUsed(limbs: *const [exact_limbs]u64, limit: u8) u8 {
    var index: usize = limit;
    while (index > 0) {
        index -= 1;
        if (limbs[index] != 0) return @intCast(index + 1);
    }
    return 0;
}

/// The number of limbs a producer can have written, clamped to the storage.
fn usedBound(limit: usize) u8 {
    return @intCast(@min(limit, exact_limbs));
}

fn magnitudeCompare(a: *const Big, b: *const Big) std.math.Order {
    if (a.used != b.used) return if (a.used > b.used) .gt else .lt;
    var index: usize = a.used;
    while (index > 0) {
        index -= 1;
        const left = a.limbs[index];
        const right = b.limbs[index];
        if (left != right) return if (left > right) .gt else .lt;
    }
    return .eq;
}

/// `|a| + |b|` with the sign `negative`, written into `out`. `out` MAY BE `a`
/// OR `b`: each index is read before it is written and the index only climbs.
///
/// The bound is `usedBound(width + 1)`, the same expression the `used` below
/// is taken from, and the loop ASSIGNS every index in `0..width`. The only
/// index inside the bound the loop does not write is `width` itself, which the
/// carry line stores - as `0` when there is no carry, which is exactly what a
/// cleared array held - so the region the scan reads is fully defined without
/// a clear and no clear is paid for.
fn magnitudeSum(out: *Big, a: *const Big, b: *const Big, negative: bool) void {
    const width = @max(a.used, b.used);
    var carry: u64 = 0;
    var index: usize = 0;
    while (index < width) : (index += 1) {
        const left = if (index < a.used) a.limbs[index] else 0;
        const right = if (index < b.used) b.limbs[index] else 0;
        const total: u128 = @as(u128, left) + @as(u128, right) + @as(u128, carry);
        out.limbs[index] = @truncate(total);
        carry = @intCast(total >> 64);
    }
    // The loop writes `0..width` and the carry adds at most the one limb above
    // them. That store is unconditional so `width` is written either way: a
    // zero carry leaves the same zero a cleared array would have held, and the
    // scan below reads that index whatever it is.
    if (width < exact_limbs) out.limbs[width] = carry;
    out.negative = negative;
    out.used = topUsed(&out.limbs, usedBound(width + 1));
}

/// `|a| - |b|` in magnitude, written into `out`. The caller has already
/// ordered the operands, so the final borrow always resolves inside the loop
/// and the result fits. `out` MAY BE `a` OR `b`, for the reason
/// `magnitudeSum` gives.
///
/// The arithmetic is WRAPPING, not checked: an intermediate
/// `left - right - borrow` goes below zero exactly when it is negative, and a
/// checked unsigned subtraction would turn that ordinary borrow into a panic.
/// A wrapped result has every bit of `total` above the limb set, so the bit
/// just below them is the borrow out: it is ONE exactly when the intermediate
/// underflowed.
///
/// The bound is `a.used`, which is what the scan below is given. The loop
/// ASSIGNS every index in `0..a.used`, so every index the scan reads has been
/// written and nothing at or above `a.used` is cleared.
fn magnitudeDifference(out: *Big, a: *const Big, b: *const Big, negative: bool) void {
    var borrow: u64 = 0;
    var index: usize = 0;
    while (index < a.used) : (index += 1) {
        const left = a.limbs[index];
        const right = if (index < b.used) b.limbs[index] else 0;
        const total: u128 = @as(u128, left) -% @as(u128, right) -% @as(u128, borrow);
        out.limbs[index] = @truncate(total);
        borrow = @intCast((total >> 64) & 1);
    }
    // The loop writes `0..a.used` and nothing else, and a difference of
    // magnitudes can only shrink the high limbs.
    out.negative = negative;
    out.used = topUsed(&out.limbs, a.used);
}

/// `[dim]Big` of the exact values of a row of stored coordinates, built ONCE
/// and then shared.
///
/// Every predicate below takes its inputs as cached anchors and forms every
/// factor and every difference from them, which is what bounds the work: a
/// question that reads `k` stored coordinates converts `k` of them, exactly
/// once, however many products and differences it goes on to form. The
/// conversion itself writes two limbs straight into the destination, so a row
/// costs one array and two stores per coordinate rather than a cleared array
/// and a copy per coordinate.
///
/// A stored coordinate and an exact difference of two of them carry the SAME
/// `2^1074` scale, so one can stand in for the other wherever a predicate
/// asks for a factor.
fn exactRow(comptime dim: usize, comptime T: type, a: [dim]T) [dim]Big {
    var out: [dim]Big = undefined;
    exactRowInto(T, dim, &out, a);
    return out;
}

/// The same row written into a caller-owned destination, so a caller that
/// already has the array does not pay for a second one.
fn exactRowInto(comptime T: type, comptime dim: usize, out: *[dim]Big, a: [dim]T) void {
    inline for (0..dim) |index| {
        Big.fromFloat(&out[index], widen(T, a[index]));
    }
}

/// `[dim]Big` of the exact differences `to - from`, componentwise, written
/// into a caller-owned destination. Both rows are read through a pointer: a
/// row is `dim` whole limb arrays, and copying two of them into the call is
/// what a predicate used to spend most of its copies on.
///
/// The subtraction happens INSIDE the exact arithmetic, never in floating
/// point. A floating-point subtraction is where the exact policy used to lose
/// the sign, and on a pair of large opposite coordinates it overflows to an
/// infinity that is not a real number at all. Both operands carry the same
/// POSITIVE `2^1074`, so their difference carries it too and needs no
/// rescaling. A difference of two finite values is below `2 * floatMax <
/// 2^1025`, so the scaled result is below `2^2099`: one binary digit more
/// than a stored coordinate.
fn exactDeltasInto(comptime dim: usize, out: *[dim]Big, from: *const [dim]Big, to: *const [dim]Big) void {
    inline for (0..dim) |index| {
        Big.subInto(&out[index], &to[index], &from[index]);
    }
}

/// Both tangency combinations of one pair of radii, over ONE conversion of
/// each radius. Both combinations are exact values of their own, because a
/// rounded sum or a rounded difference is the whole difference between two
/// balls that cross and two that touch, and both come off the same pair of
/// anchors.
fn exactRadiusPairInto(out: *ExactRadiusPair, comptime T: type, a: T, b: T) void {
    // SAFETY: Big.fromFloat writes metadata and active limbs before use.
    var first: Big = undefined;
    // SAFETY: Big.fromFloat writes metadata and active limbs before use.
    var second: Big = undefined;
    Big.fromFloat(&first, widen(T, a));
    Big.fromFloat(&second, widen(T, b));
    Big.addInto(&out.sum, &first, &second);
    Big.subInto(&out.difference, &first, &second);
}

/// The exact value of `dot(x, y)` over exact components, degree two, written
/// into a caller-owned destination. Both rows are read through a pointer and
/// the products land in ONE slot the loop reuses, so the row costs no copy
/// and the products cost no array each.
fn exactDotInto(comptime dim: usize, out: *Big, x: *const [dim]Big, y: *const [dim]Big) void {
    // SAFETY: Big.mulInto writes metadata and active limbs before each addition.
    var product: Big = undefined;
    Big.setZero(out);
    inline for (0..dim) |index| {
        Big.mulInto(&product, &x[index], &y[index]);
        Big.addInto(out, out, &product);
    }
}

/// The exact value of `dot(w, w)` over exact components, degree two, written
/// into a caller-owned destination. Naming it once keeps the squared norm
/// from converting or copying the whole row a second time.
fn exactNormInto(comptime dim: usize, out: *Big, w: *const [dim]Big) void {
    // SAFETY: Big.mulInto writes metadata and active limbs before each addition.
    var product: Big = undefined;
    Big.setZero(out);
    inline for (0..dim) |index| {
        Big.mulInto(&product, &w[index], &w[index]);
        Big.addInto(out, out, &product);
    }
}

/// The exact sign of `p * q - r * s`, degree two, over four exact components
/// that already carry the `2^1074` scale. Every one of them is read through a
/// pointer, so a determinant never copies a whole limb array per coefficient.
fn exactDet2(p: *const Big, q: *const Big, r: *const Big, s: *const Big) Sign {
    // SAFETY: Big.mulInto writes metadata and active limbs before subtraction.
    var left: Big = undefined;
    // SAFETY: Big.mulInto writes metadata and active limbs before subtraction.
    var right: Big = undefined;
    Big.mulInto(&left, p, q);
    Big.mulInto(&right, r, s);
    Big.subInto(&left, &left, &right);
    return Big.sign(&left);
}

/// The exact sign of `dot(w, w) - k * k`, degree two, over an exact difference
/// row and an exact radius. Both carry the `2^1074` scale, so both squares
/// carry `2^2148` and nothing needs rescaling.
fn squaredNormMinusSquareExact(comptime dim: usize, w: *const [dim]Big, k: *const Big) Sign {
    // SAFETY: exactNormInto writes metadata and active limbs before subtraction.
    var total: Big = undefined;
    // SAFETY: Big.mulInto writes metadata and active limbs before subtraction.
    var square: Big = undefined;
    exactNormInto(dim, &total, w);
    Big.mulInto(&square, k, k);
    Big.subInto(&total, &total, &square);
    return Big.sign(&total);
}

/// The exact sign of the two-component oriented area `(b - a) x (c - a)`,
/// degree two, over the ORIGINAL stored coordinates. The four differences are
/// exact values, so a pair of coordinates whose difference is past the
/// binary64 range still decides its own sign. `a` is converted once and both
/// difference rows read that one anchor.
fn exactOrient2(comptime T: type, a: [2]T, b: [2]T, c: [2]T) Sign {
    const anchor = exactRow(2, T, a);
    const b_row = exactRow(2, T, b);
    const c_row = exactRow(2, T, c);
    var ab: [2]Big = undefined;
    var ac: [2]Big = undefined;
    exactDeltasInto(2, &ab, &anchor, &b_row);
    exactDeltasInto(2, &ac, &anchor, &c_row);
    return exactDet2(&ab[0], &ac[1], &ab[1], &ac[0]);
}

/// The exact sign of `dot(cross(b - a, c - a), d - a)`, degree three, over
/// the ORIGINAL stored coordinates, with `a` converted once and the three
/// difference rows read off that one anchor. The expression matches the fast
/// path one level down, so the two policies differ only in the arithmetic and
/// never in the algebra.
///
/// The three cross components are each a difference of two products, formed
/// one at a time through the SAME pair of slots, with the running total as the
/// third. The sequence of products, differences and additions is the one the
/// fast path reads, and no product of one component is ever live while another
/// is formed.
fn exactOrient3(comptime T: type, a: [3]T, b: [3]T, c: [3]T, d: [3]T) Sign {
    const anchor = exactRow(3, T, a);
    const b_row = exactRow(3, T, b);
    const c_row = exactRow(3, T, c);
    const d_row = exactRow(3, T, d);
    var u: [3]Big = undefined;
    var v: [3]Big = undefined;
    var w: [3]Big = undefined;
    exactDeltasInto(3, &u, &anchor, &b_row);
    exactDeltasInto(3, &v, &anchor, &c_row);
    exactDeltasInto(3, &w, &anchor, &d_row);
    // SAFETY: Big.mulInto writes metadata and active limbs before subtraction.
    var left: Big = undefined;
    // SAFETY: Big.mulInto writes metadata and active limbs before subtraction.
    var right: Big = undefined;
    // SAFETY: Big.mulInto writes metadata and active limbs before accumulation.
    var total: Big = undefined;
    Big.mulInto(&left, &u[1], &v[2]);
    Big.mulInto(&right, &u[2], &v[1]);
    Big.subInto(&left, &left, &right);
    Big.mulInto(&total, &left, &w[0]);
    Big.mulInto(&left, &u[2], &v[0]);
    Big.mulInto(&right, &u[0], &v[2]);
    Big.subInto(&left, &left, &right);
    Big.mulInto(&right, &left, &w[1]);
    Big.addInto(&total, &total, &right);
    Big.mulInto(&left, &u[0], &v[1]);
    Big.mulInto(&right, &u[1], &v[0]);
    Big.subInto(&left, &left, &right);
    Big.mulInto(&right, &left, &w[2]);
    Big.addInto(&total, &total, &right);
    return Big.sign(&total);
}

/// The four points as the sixteen exact values the determinant reads: one
/// conversion each, written straight into a caller-owned matrix, and every
/// minor the expansion asks about reads these cells.
///
/// `exactMinor3Columns` already builds this shape for the four-dimensional
/// coplanarity gate, and the degree-four determinant gets it for the same
/// reason: `fromFloat` was a conversion and a tail scan apiece, so reading a
/// cell out of the stored matrix re-converted the same coordinate once per
/// minor per column.
fn exactCells4Into(comptime T: type, out: *[4][4]Big, rows: [4][4]T) void {
    inline for (0..4) |row| {
        exactRowInto(T, 4, &out[row], rows[row]);
    }
}

/// The exact three-by-three minor of a four-by-four matrix, written into
/// `out`: the determinant of the three selected rows against the three
/// selected columns, expanded along the third column into three terms of two
/// products each. Each cell comes back as a POINTER, so reading one never
/// copies a limb array, and the two working slots below are reused across all
/// three terms. The products, differences and additions are the ones the plain
/// `minor3FastOf` evaluation reads, in the same order.
fn exactMinor3Into(
    out: *Big,
    cells: *const [4][4]Big,
    r0: usize,
    r1: usize,
    r2: usize,
    c0: usize,
    c1: usize,
    c2: usize,
) void {
    // SAFETY: Big.mulInto writes metadata and active limbs before subtraction.
    var left: Big = undefined;
    // SAFETY: Big.mulInto writes metadata and active limbs before subtraction.
    var right: Big = undefined;
    Big.mulInto(&left, cell4(cells, r1, c0), cell4(cells, r2, c1));
    Big.mulInto(&right, cell4(cells, r1, c1), cell4(cells, r2, c0));
    Big.subInto(&left, &left, &right);
    Big.mulInto(out, cell4(cells, r0, c2), &left);
    Big.mulInto(&left, cell4(cells, r0, c0), cell4(cells, r2, c1));
    Big.mulInto(&right, cell4(cells, r0, c1), cell4(cells, r2, c0));
    Big.subInto(&left, &left, &right);
    Big.mulInto(&right, cell4(cells, r1, c2), &left);
    Big.subInto(out, out, &right);
    Big.mulInto(&left, cell4(cells, r0, c0), cell4(cells, r1, c1));
    Big.mulInto(&right, cell4(cells, r0, c1), cell4(cells, r1, c0));
    Big.subInto(&left, &left, &right);
    Big.mulInto(&right, cell4(cells, r2, c2), &left);
    Big.addInto(out, out, &right);
}

fn cell4(cells: *const [4][4]Big, row: usize, column: usize) *const Big {
    return &cells[row][column];
}

/// The exact sign of the four-by-four determinant whose rows are the four
/// points, degree four. The expansion runs along the first row, so the
/// cofactor signs alternate and each minor skips the column its own
/// coefficient multiplies.
fn exactOrient4(comptime T: type, a: [4]T, b: [4]T, c: [4]T, d: [4]T) Sign {
    const rows = [4][4]T{ a, b, c, d };
    var cells: [4][4]Big = undefined;
    exactCells4Into(T, &cells, rows);
    // SAFETY: exactMinor3Into writes metadata and active limbs before multiplication.
    var minor: Big = undefined;
    // SAFETY: Big.mulInto writes metadata and active limbs before subtraction.
    var product: Big = undefined;
    // SAFETY: Big.mulInto writes metadata and active limbs before accumulation.
    var total: Big = undefined;
    exactMinor3Into(&minor, &cells, 1, 2, 3, 1, 2, 3);
    Big.mulInto(&total, &minor, cell4(&cells, 0, 0));
    exactMinor3Into(&minor, &cells, 1, 2, 3, 0, 2, 3);
    Big.mulInto(&product, &minor, cell4(&cells, 0, 1));
    Big.subInto(&total, &total, &product);
    exactMinor3Into(&minor, &cells, 1, 2, 3, 0, 1, 3);
    Big.mulInto(&product, &minor, cell4(&cells, 0, 2));
    Big.addInto(&total, &total, &product);
    exactMinor3Into(&minor, &cells, 1, 2, 3, 0, 1, 2);
    Big.mulInto(&product, &minor, cell4(&cells, 0, 3));
    Big.subInto(&total, &total, &product);
    return Big.sign(&total);
}

/// The plain evaluation of the same four-by-four determinant, with the same
/// expansion along the first row and the same grouping inside each minor, so
/// that the two policies differ only in the arithmetic and never in the
/// algebra.
fn fastOrient4(comptime T: type, a: [4]T, b: [4]T, c: [4]T, d: [4]T) T {
    @setFloatMode(.strict);
    const rows = [4][4]T{ a, b, c, d };
    return minor3FastOf(T, rows, 1, 2, 3, 1, 2, 3) * rows[0][0] -
        minor3FastOf(T, rows, 1, 2, 3, 0, 2, 3) * rows[0][1] +
        minor3FastOf(T, rows, 1, 2, 3, 0, 1, 3) * rows[0][2] -
        minor3FastOf(T, rows, 1, 2, 3, 0, 1, 2) * rows[0][3];
}

fn minor3FastOf(comptime T: type, rows: [4][4]T, r0: usize, r1: usize, r2: usize, c0: usize, c1: usize, c2: usize) T {
    @setFloatMode(.strict);
    return rows[r0][c2] * (rows[r1][c0] * rows[r2][c1] - rows[r1][c1] * rows[r2][c0]) -
        rows[r1][c2] * (rows[r0][c0] * rows[r2][c1] - rows[r0][c1] * rows[r2][c0]) +
        rows[r2][c2] * (rows[r0][c0] * rows[r1][c1] - rows[r0][c1] * rows[r1][c0]);
}

/// The sign of the degree-four form, from either policy.
/// `error.NonFiniteResult` is the `fast` policy saying its own evaluation came
/// out a NaN; `adaptive_exact` answers for any finite coordinates.
fn orient4Sign(comptime T: type, a: [4]T, b: [4]T, c: [4]T, d: [4]T, policy: PredicatePolicy) error{NonFiniteResult}!Sign {
    return switch (policy) {
        .fast => fastSign(T, fastOrient4(T, a, b, c, d)),
        .adaptive_exact => exactOrient4(T, a, b, c, d),
    };
}

/// The three COLUMNS of the four-dimensional coplanarity matrix, as exact
/// values: `a_finish - a_start`, `b_finish - b_start`, and `b_start -
/// a_start`. Each of the four endpoints is converted ONCE and every minor the
/// gate asks about reads these columns, so the gate costs sixteen conversions
/// instead of one per minor per column.
fn exactMinor3Columns(comptime T: type, a_start: [4]T, a_finish: [4]T, b_start: [4]T, b_finish: [4]T) [3][4]Big {
    const first_start = exactRow(4, T, a_start);
    const second_start = exactRow(4, T, b_start);
    const first_finish = exactRow(4, T, a_finish);
    const second_finish = exactRow(4, T, b_finish);
    var out: [3][4]Big = undefined;
    exactDeltasInto(4, &out[0], &first_start, &first_finish);
    exactDeltasInto(4, &out[1], &second_start, &second_finish);
    exactDeltasInto(4, &out[2], &first_start, &second_start);
    return out;
}

/// One entry of the three-by-four exact matrix, read through a pointer so a
/// single cell never copies the columns.
fn cellOf(comptime dim: usize, columns: *const [3][dim]Big, row: usize, column: usize) *const Big {
    return &columns.*[column][row];
}

/// The three-by-three determinant of three exact COLUMNS over the three
/// selected rows, expanded along the THIRD column into three terms of two
/// products each. The third column is the one the caller selected rows
/// against, so this is the cofactor expansion the four-by-four determinant
/// needs, and the exact and fast paths use the same grouping.
fn minor3ExactInto(out: *Big, comptime dim: usize, columns: *const [3][dim]Big, rows: [3]usize) void {
    // SAFETY: Big.mulInto writes metadata and active limbs before subtraction.
    var left: Big = undefined;
    // SAFETY: Big.mulInto writes metadata and active limbs before subtraction.
    var right: Big = undefined;
    Big.mulInto(&left, cellOf(dim, columns, rows[1], 0), cellOf(dim, columns, rows[2], 1));
    Big.mulInto(&right, cellOf(dim, columns, rows[1], 1), cellOf(dim, columns, rows[2], 0));
    Big.subInto(&left, &left, &right);
    Big.mulInto(out, cellOf(dim, columns, rows[0], 2), &left);
    Big.mulInto(&left, cellOf(dim, columns, rows[0], 1), cellOf(dim, columns, rows[2], 0));
    Big.mulInto(&right, cellOf(dim, columns, rows[0], 0), cellOf(dim, columns, rows[2], 1));
    Big.subInto(&left, &left, &right);
    Big.mulInto(&right, cellOf(dim, columns, rows[1], 2), &left);
    Big.addInto(out, out, &right);
    Big.mulInto(&left, cellOf(dim, columns, rows[0], 0), cellOf(dim, columns, rows[1], 1));
    Big.mulInto(&right, cellOf(dim, columns, rows[0], 1), cellOf(dim, columns, rows[1], 0));
    Big.subInto(&left, &left, &right);
    Big.mulInto(&right, cellOf(dim, columns, rows[2], 2), &left);
    Big.addInto(out, out, &right);
}

fn fastMinor3x4(comptime T: type, u: [4]T, v: [4]T, w: [4]T, rows: [3]usize) T {
    @setFloatMode(.strict);
    const columns = [3][4]T{ u, v, w };
    const m: [3][3]T = .{
        .{ columns[0][rows[0]], columns[1][rows[0]], columns[2][rows[0]] },
        .{ columns[0][rows[1]], columns[1][rows[1]], columns[2][rows[1]] },
        .{ columns[0][rows[2]], columns[1][rows[2]], columns[2][rows[2]] },
    };
    return m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]) +
        m[1][2] * (m[0][1] * m[2][0] - m[0][0] * m[2][1]) +
        m[2][2] * (m[0][0] * m[1][1] - m[0][1] * m[1][0]);
}

/// The two sign answers a segment question needs, one per end. A segment is
/// measured against ONE flat or ONE ball, so the shared centre is converted
/// once and both ends read that one anchor.
const EndSigns = struct { start: Sign, end: Sign };

/// The two tangency comparisons a pair of balls turns on, named so that the
/// sign of each is unambiguous at the call.
const TangencySigns = struct { outer: Sign, inner: Sign };

/// The two combinations of one pair of radii that a tangency comparison turns
/// on: the SUM and the DIFFERENCE.
const ExactRadiusPair = struct { sum: Big, difference: Big };

/// The sign of `dot(a - o, n)`, which is the signed distance of `a` from the
/// flat through `o` with unit normal `n`, degree one over the ORIGINAL stored
/// coordinates. `error.NonFiniteResult` is the `fast` policy saying its own
/// evaluation came out a NaN; on coordinates at opposite ends of the range
/// the differences overflow to infinities of opposite sign and their sum is
/// exactly that case.
fn flatSideSign(comptime dim: usize, comptime T: type, a: [dim]T, o: [dim]T, n: [dim]T, policy: PredicatePolicy) error{NonFiniteResult}!Sign {
    return switch (policy) {
        .fast => fastSign(T, flatDistanceFast(dim, T, a, o, n)),
        .adaptive_exact => blk: {
            const flat = exactRow(dim, T, o);
            const probe = exactRow(dim, T, a);
            const normal = exactRow(dim, T, n);
            var deltas: [dim]Big = undefined;
            // SAFETY: exactDotInto writes metadata and active limbs before sign extraction.
            var total: Big = undefined;
            exactDeltasInto(dim, &deltas, &flat, &probe);
            exactDotInto(dim, &total, &deltas, &normal);
            break :blk Big.sign(&total);
        },
    };
}

/// The two ends of a segment measured against ONE flat, over the ORIGINAL
/// stored coordinates. The flat point and the normal are converted once and
/// both ends read them. `error.NonFiniteResult` is the `fast` policy saying
/// its own evaluation came out a NaN.
fn segmentSideSigns(
    comptime dim: usize,
    comptime T: type,
    o: [dim]T,
    n: [dim]T,
    a: [dim]T,
    b: [dim]T,
    policy: PredicatePolicy,
) error{NonFiniteResult}!EndSigns {
    return switch (policy) {
        .fast => .{
            .start = try fastSign(T, flatDistanceFast(dim, T, a, o, n)),
            .end = try fastSign(T, flatDistanceFast(dim, T, b, o, n)),
        },
        .adaptive_exact => blk: {
            const flat = exactRow(dim, T, o);
            const normal = exactRow(dim, T, n);
            var deltas: [dim]Big = undefined;
            // SAFETY: exactDotInto writes metadata and active limbs before sign extraction.
            var total: Big = undefined;
            const first = exactRow(dim, T, a);
            exactDeltasInto(dim, &deltas, &flat, &first);
            exactDotInto(dim, &total, &deltas, &normal);
            const start = Big.sign(&total);
            const second = exactRow(dim, T, b);
            exactDeltasInto(dim, &deltas, &flat, &second);
            exactDotInto(dim, &total, &deltas, &normal);
            break :blk .{ .start = start, .end = Big.sign(&total) };
        },
    };
}

/// The signed distance of `a` from the flat through `o` with unit normal `n`,
/// the plain scalar evaluation whose sign `fast` reports. The terms are
/// summed left to right in declaration order.
fn flatDistanceFast(comptime dim: usize, comptime T: type, a: [dim]T, o: [dim]T, n: [dim]T) T {
    @setFloatMode(.strict);
    var sum: T = 0;
    inline for (0..dim) |index| {
        sum = sum + (a[index] - o[index]) * n[index];
    }
    return sum;
}

fn dotOfArray(comptime dim: usize, comptime T: type, a: [dim]T, b: [dim]T) T {
    @setFloatMode(.strict);
    var sum: T = 0;
    inline for (0..dim) |index| {
        sum = sum + a[index] * b[index];
    }
    return sum;
}

// ---------------------------------------------------------------------------
// Vector
// ---------------------------------------------------------------------------

/// The raw scalar kernels every arithmetic helper composes. Each sets strict
/// float mode, which is what stops the compiler from reassociating a sum or
/// contracting a product into the addition that follows it. The `fast`
/// predicate policy is DEFINED as the sign these produce, so a contraction
/// would silently change what `fast` means.
fn addRaw(comptime T: type, a: T, b: T) T {
    @setFloatMode(.strict);
    return a + b;
}

fn subRaw(comptime T: type, a: T, b: T) T {
    @setFloatMode(.strict);
    return a - b;
}

fn negRaw(comptime T: type, a: T) T {
    @setFloatMode(.strict);
    return -a;
}

fn scaleRaw(comptime T: type, a: T, b: T) T {
    @setFloatMode(.strict);
    return a * b;
}

fn divRaw(comptime T: type, a: T, b: T) T {
    @setFloatMode(.strict);
    return a / b;
}

/// Which componentwise operation a vector helper applies. A comptime tag
/// rather than a function value, because the raw kernels are generic over the
/// scalar type and a generic function does not coerce to a non-generic
/// function pointer.
const BinaryOp = enum { add, sub };
const UnaryOp = enum { neg };

fn applyBinary(comptime op: BinaryOp, comptime T: type, a: T, b: T) T {
    return switch (op) {
        .add => addRaw(T, a, b),
        .sub => subRaw(T, a, b),
    };
}

fn applyUnary(comptime op: UnaryOp, comptime T: type, a: T) T {
    return switch (op) {
        .neg => negRaw(T, a),
    };
}

/// The componentwise operation over two whole values: both inputs are gated,
/// and so is the result, because a sum or a difference of two finite values
/// can leave the format.
fn vectorBinary(comptime V: type, a: V, b: V, comptime op: BinaryOp) error{ NonFiniteInput, NonFiniteResult }!V {
    try requireFiniteField(a);
    try requireFiniteField(b);
    return vectorBinaryChecked(V, a, b, op);
}

/// The same operation over two values whose components an earlier gate already
/// covered. It checks only the RESULT, which is the half no input check can
/// have reached.
fn vectorBinaryChecked(comptime V: type, a: V, b: V, comptime op: BinaryOp) error{NonFiniteResult}!V {
    const width = comptime fieldCount(V);
    var out: V = undefined;
    inline for (0..width) |index| {
        const name = comptime fieldName(index);
        @field(out, name) = applyBinary(op, componentType(V), @field(a, name), @field(b, name));
    }
    return requireFiniteComponents(V, out);
}

/// The negation of a gated value. Negating a finite value cannot leave the
/// format, so this reports input failure only, and that is why its result set
/// is narrower than the binary helpers above.
fn vectorNegated(comptime V: type, a: V) error{NonFiniteInput}!V {
    try requireFiniteField(a);
    const width = comptime fieldCount(V);
    var out: V = undefined;
    inline for (0..width) |index| {
        @field(out, fieldName(index)) = applyUnary(.neg, componentType(V), @field(a, fieldName(index)));
    }
    return out;
}

fn vectorScaled(comptime V: type, a: V, factor: componentType(V)) error{ NonFiniteInput, NonFiniteResult }!V {
    try requireFiniteField(a);
    try requireFiniteScalar(componentType(V), factor);
    return vectorScaledChecked(V, a, factor);
}

/// The scaling over a value and a factor an earlier gate already covered, with
/// the result gated.
fn vectorScaledChecked(comptime V: type, a: V, factor: componentType(V)) error{NonFiniteResult}!V {
    const width = comptime fieldCount(V);
    var out: V = undefined;
    inline for (0..width) |index| {
        const name = comptime fieldName(index);
        @field(out, name) = scaleRaw(componentType(V), @field(a, name), factor);
    }
    return requireFiniteComponents(V, out);
}

/// The scalar product, summed left to right in declaration order. The order is
/// the contract: `fast` is defined to report the sign of exactly this
/// evaluation, so a caller that groups the terms differently gets a different
/// `fast` sign. That is why the order is written down and not left to the
/// optimizer. The sum can overflow, so the value is gated here rather than
/// where it is read.
fn vectorDot(comptime V: type, a: V, b: V) error{NonFiniteResult}!componentType(V) {
    @setFloatMode(.strict);
    const T = componentType(V);
    const width = comptime fieldCount(V);
    var sum: T = 0;
    inline for (0..width) |index| {
        const name = comptime fieldName(index);
        sum = sum + (@field(a, name) * @field(b, name));
    }
    return requireFiniteResult(T, sum);
}

fn vectorLengthSquared(comptime V: type, a: V) error{NonFiniteResult}!componentType(V) {
    return vectorDot(V, a, a);
}

/// The length, folded through the accepted `hypot` kernel from a positive
/// zero. The fold cannot flush a representable norm to zero, which
/// `sqrt(x * x + y * y)` does in both directions. Nesting a two-argument
/// kernel is not correctly rounded for three or four components; it is the
/// accepted kernel applied again, never a re-derived polynomial.
///
/// It can still leave the format: `hypot` of two components at the top of the
/// range is `sqrt(2)` times either of them, and that is not representable. So
/// the fold no longer claims the length cannot overflow, and a norm that runs
/// off the end is `error.NonFiniteResult` rather than an infinity the caller
/// would then normalize into a zero direction.
fn vectorLength(comptime V: type, a: V) error{NonFiniteResult}!componentType(V) {
    const T = componentType(V);
    var accumulator: T = 0;
    inline for (toArrayOf(V, a)) |component| {
        accumulator = hypotOf(T, accumulator, component);
    }
    return requireFiniteResult(T, accumulator);
}

/// The unit vector along `a`. A zero vector has no direction, an infinite
/// component has none either, and a norm past the top of the range cannot
/// produce one either: dividing by an infinite length would hand back the zero
/// vector as though a zero vector had a direction.
fn vectorNormalize(comptime V: type, a: V) error{ NonFiniteInput, NonFiniteResult, Degenerate }!V {
    try requireFiniteField(a);
    const length = try vectorLength(V, a);
    if (length == 0) return error.Degenerate;
    const T = componentType(V);
    var out: V = undefined;
    const width = comptime fieldCount(V);
    if (length < std.math.floatMin(T)) {
        // A subnormal norm can round away the component ratios. Normalize a
        // scaled vector instead; its largest component has magnitude one.
        var largest: T = 0;
        inline for (0..width) |index| {
            largest = maxOf(T, largest, absOf(T, @field(a, fieldName(index))));
        }
        inline for (0..width) |index| {
            @field(out, fieldName(index)) = divRaw(T, @field(a, fieldName(index)), largest);
        }
        const scaled_length = try vectorLength(V, out);
        inline for (0..width) |index| {
            @field(out, fieldName(index)) = divRaw(T, @field(out, fieldName(index)), scaled_length);
        }
    } else {
        inline for (0..width) |index| {
            @field(out, fieldName(index)) = divRaw(T, @field(a, fieldName(index)), length);
        }
    }
    return requireFiniteComponents(V, out);
}

fn vectorDistanceSquared(comptime V: type, a: V, b: V) error{NonFiniteResult}!componentType(V) {
    const difference = try vectorBinaryChecked(V, a, b, .sub);
    return vectorLengthSquared(V, difference);
}

fn vectorDistance(comptime V: type, a: V, b: V) error{NonFiniteResult}!componentType(V) {
    return vectorLength(V, try vectorBinaryChecked(V, a, b, .sub));
}

fn cross3Raw(comptime T: type, a: [3]T, b: [3]T) [3]T {
    @setFloatMode(.strict);
    return .{
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    };
}

/// A two-component `f32` vector.
pub const Vector2f = Vector2(f32);
/// A two-component `f64` vector.
pub const Vector2d = Vector2(f64);
/// A three-component `f32` vector.
pub const Vector3f = Vector3(f32);
/// A three-component `f64` vector.
pub const Vector3d = Vector3(f64);
/// A four-component `f32` vector.
pub const Vector4f = Vector4(f32);
/// A four-component `f64` vector.
pub const Vector4d = Vector4(f64);

fn VectorOf(comptime dim: usize, comptime T: type) type {
    return switch (dim) {
        2 => Vector2(T),
        3 => Vector3(T),
        4 => Vector4(T),
        else => @compileError("zmath geometry dimensions are 2, 3, and 4"),
    };
}

fn Vector2(comptime T: type) type {
    return struct {
        /// The first component. Public and writable; `init` established that it
        /// is finite and every operation below checks it again.
        x: T,
        /// The second component, under the same obligation as `x`.
        y: T,

        const Self = @This();

        /// Builds a vector from named components. `error.NonFiniteInput` when
        /// either component is an infinity or a NaN.
        pub fn init(x: T, y: T) error{NonFiniteInput}!Self {
            const out = Self{ .x = x, .y = y };
            try requireFiniteField(out);
            return out;
        }

        /// Builds a vector from `[x, y]`, under the same gate as `init`.
        pub fn fromArray(a: [2]T) error{NonFiniteInput}!Self {
            return init(a[0], a[1]);
        }

        /// The components as `[x, y]`, copied out.
        pub fn toArray(self: Self) [2]T {
            return toArrayOf(Self, self);
        }

        /// The componentwise sum. `error.NonFiniteResult` when a component of
        /// the sum is not finite, which two finite components reach together
        /// at the top of the range.
        pub fn add(self: Self, other: Self) error{ NonFiniteInput, NonFiniteResult }!Self {
            return vectorBinary(Self, self, other, .add);
        }

        /// The componentwise difference. `error.NonFiniteResult` when a
        /// component of the difference is not finite.
        pub fn sub(self: Self, other: Self) error{ NonFiniteInput, NonFiniteResult }!Self {
            return vectorBinary(Self, self, other, .sub);
        }

        /// The componentwise negation. Negating a finite component cannot
        /// leave the format, so this reports input failure only.
        pub fn neg(self: Self) error{NonFiniteInput}!Self {
            return vectorNegated(Self, self);
        }

        /// Every component times `factor`. `error.NonFiniteResult` when a
        /// product is not finite.
        pub fn scale(self: Self, factor: T) error{ NonFiniteInput, NonFiniteResult }!Self {
            return vectorScaled(Self, self, factor);
        }

        /// The scalar product, summed left to right in component order.
        /// `error.NonFiniteResult` when the sum is not finite.
        pub fn dot(self: Self, other: Self) error{ NonFiniteInput, NonFiniteResult }!T {
            try requireFiniteField(self);
            try requireFiniteField(other);
            return vectorDot(Self, self, other);
        }

        /// The squared length, which is the dot product with itself, and which
        /// overflows on a vector large enough that `length` does not.
        /// `error.NonFiniteResult` when the sum is not finite.
        pub fn lengthSquared(self: Self) error{ NonFiniteInput, NonFiniteResult }!T {
            try requireFiniteField(self);
            return vectorLengthSquared(Self, self);
        }

        /// The length, folded through the accepted `hypot` kernel, so a large
        /// vector keeps a finite norm and a subnormal one keeps a representable
        /// one. `error.NonFiniteResult` when the norm itself is past the top
        /// of the range.
        pub fn length(self: Self) error{ NonFiniteInput, NonFiniteResult }!T {
            try requireFiniteField(self);
            return vectorLength(Self, self);
        }

        /// The unit vector along `self`. `error.Degenerate` for a zero vector,
        /// and `error.NonFiniteResult` for a norm past the top of the range.
        pub fn normalize(self: Self) error{ NonFiniteInput, NonFiniteResult, Degenerate }!Self {
            return vectorNormalize(Self, self);
        }

        /// The squared distance to `point`, from the rounded difference.
        /// `error.NonFiniteResult` when the difference or the sum leaves the
        /// format.
        pub fn distanceSquared(self: Self, point: Point2(T)) error{ NonFiniteInput, NonFiniteResult }!T {
            try requireFiniteField(self);
            try requireFiniteField(point);
            return vectorDistanceSquared(Self, self, toVectorOfPoint(Point2(T), point));
        }

        /// The distance to `point`, the length of the rounded difference.
        /// `error.NonFiniteResult` when the difference leaves the format.
        pub fn distance(self: Self, point: Point2(T)) error{ NonFiniteInput, NonFiniteResult }!T {
            try requireFiniteField(self);
            try requireFiniteField(point);
            return vectorDistance(Self, self, toVectorOfPoint(Point2(T), point));
        }

        /// True when both components are finite. This is the check a caller
        /// needs after writing the fields directly.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(Self, self);
        }

        /// Bit equality, with every NaN equal to every other NaN and `+0`
        /// distinct from `-0`. It compares stored bits and never fails.
        pub fn eq(self: Self, other: Self) bool {
            return eqOf(Self, self, other);
        }
    };
}

fn Vector3(comptime T: type) type {
    return struct {
        /// The first component. Public and writable; `init` established that it
        /// is finite and every operation below checks it again.
        x: T,
        /// The second component, under the same obligation as `x`.
        y: T,
        /// The third component, under the same obligation as `x`.
        z: T,

        const Self = @This();

        /// Builds a vector from named components. `error.NonFiniteInput` when
        /// a component is an infinity or a NaN.
        pub fn init(x: T, y: T, z: T) error{NonFiniteInput}!Self {
            const out = Self{ .x = x, .y = y, .z = z };
            try requireFiniteField(out);
            return out;
        }

        /// Builds a vector from `[x, y, z]`, under the same gate as `init`.
        pub fn fromArray(a: [3]T) error{NonFiniteInput}!Self {
            return init(a[0], a[1], a[2]);
        }

        /// The components as `[x, y, z]`, copied out.
        pub fn toArray(self: Self) [3]T {
            return toArrayOf(Self, self);
        }

        /// The componentwise sum. `error.NonFiniteResult` when a component of
        /// the sum is not finite.
        pub fn add(self: Self, other: Self) error{ NonFiniteInput, NonFiniteResult }!Self {
            return vectorBinary(Self, self, other, .add);
        }

        /// The componentwise difference. `error.NonFiniteResult` when a
        /// component of the difference is not finite.
        pub fn sub(self: Self, other: Self) error{ NonFiniteInput, NonFiniteResult }!Self {
            return vectorBinary(Self, self, other, .sub);
        }

        /// The componentwise negation. Negating a finite component cannot
        /// leave the format, so this reports input failure only.
        pub fn neg(self: Self) error{NonFiniteInput}!Self {
            return vectorNegated(Self, self);
        }

        /// Every component times `factor`. `error.NonFiniteResult` when a
        /// product is not finite.
        pub fn scale(self: Self, factor: T) error{ NonFiniteInput, NonFiniteResult }!Self {
            return vectorScaled(Self, self, factor);
        }

        /// The scalar product, summed left to right in component order.
        /// `error.NonFiniteResult` when the sum is not finite.
        pub fn dot(self: Self, other: Self) error{ NonFiniteInput, NonFiniteResult }!T {
            try requireFiniteField(self);
            try requireFiniteField(other);
            return vectorDot(Self, self, other);
        }

        /// The squared length, which is the dot product with itself, and which
        /// overflows on a vector large enough that `length` does not.
        /// `error.NonFiniteResult` when the sum is not finite.
        pub fn lengthSquared(self: Self) error{ NonFiniteInput, NonFiniteResult }!T {
            try requireFiniteField(self);
            return vectorLengthSquared(Self, self);
        }

        /// The length, folded through the accepted `hypot` kernel, so a large
        /// vector keeps a finite norm and a subnormal one keeps a representable
        /// one. `error.NonFiniteResult` when the norm itself is past the top
        /// of the range.
        pub fn length(self: Self) error{ NonFiniteInput, NonFiniteResult }!T {
            try requireFiniteField(self);
            return vectorLength(Self, self);
        }

        /// The unit vector along `self`. `error.Degenerate` for a zero vector,
        /// and `error.NonFiniteResult` for a norm past the top of the range.
        pub fn normalize(self: Self) error{ NonFiniteInput, NonFiniteResult, Degenerate }!Self {
            return vectorNormalize(Self, self);
        }

        /// The right-handed cross product. It exists in three dimensions and
        /// nowhere else, because only there does the orthogonal complement of
        /// the plane two vectors span come out as a single vector.
        /// `error.NonFiniteResult` when a product or a difference of products
        /// is not finite.
        pub fn cross(self: Self, other: Self) error{ NonFiniteInput, NonFiniteResult }!Self {
            try requireFiniteField(self);
            try requireFiniteField(other);
            return requireFiniteComponents(
                Self,
                fromArrayOf(Self, cross3Raw(T, self.toArray(), other.toArray())),
            );
        }

        /// The squared distance to `point`, from the rounded difference.
        /// `error.NonFiniteResult` when the difference or the sum leaves the
        /// format.
        pub fn distanceSquared(self: Self, point: Point3(T)) error{ NonFiniteInput, NonFiniteResult }!T {
            try requireFiniteField(self);
            try requireFiniteField(point);
            return vectorDistanceSquared(Self, self, toVectorOfPoint(Point3(T), point));
        }

        /// The distance to `point`, the length of the rounded difference.
        /// `error.NonFiniteResult` when the difference leaves the format.
        pub fn distance(self: Self, point: Point3(T)) error{ NonFiniteInput, NonFiniteResult }!T {
            try requireFiniteField(self);
            try requireFiniteField(point);
            return vectorDistance(Self, self, toVectorOfPoint(Point3(T), point));
        }

        /// True when every component is finite. This is the check a caller
        /// needs after writing the fields directly.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(Self, self);
        }

        /// Bit equality, with every NaN equal to every other NaN and `+0`
        /// distinct from `-0`. It compares stored bits and never fails.
        pub fn eq(self: Self, other: Self) bool {
            return eqOf(Self, self, other);
        }
    };
}

fn Vector4(comptime T: type) type {
    return struct {
        /// The first component. Public and writable; `init` established that it
        /// is finite and every operation below checks it again.
        x: T,
        /// The second component, under the same obligation as `x`.
        y: T,
        /// The third component, under the same obligation as `x`.
        z: T,
        /// The fourth component, under the same obligation as `x`.
        w: T,

        const Self = @This();

        /// Builds a vector from named components. `error.NonFiniteInput` when
        /// a component is an infinity or a NaN.
        pub fn init(x: T, y: T, z: T, w: T) error{NonFiniteInput}!Self {
            const out = Self{ .x = x, .y = y, .z = z, .w = w };
            try requireFiniteField(out);
            return out;
        }

        /// Builds a vector from `[x, y, z, w]`, under the same gate as `init`.
        pub fn fromArray(a: [4]T) error{NonFiniteInput}!Self {
            return init(a[0], a[1], a[2], a[3]);
        }

        /// The components as `[x, y, z, w]`, copied out.
        pub fn toArray(self: Self) [4]T {
            return toArrayOf(Self, self);
        }

        /// The componentwise sum. `error.NonFiniteResult` when a component of
        /// the sum is not finite.
        pub fn add(self: Self, other: Self) error{ NonFiniteInput, NonFiniteResult }!Self {
            return vectorBinary(Self, self, other, .add);
        }

        /// The componentwise difference. `error.NonFiniteResult` when a
        /// component of the difference is not finite.
        pub fn sub(self: Self, other: Self) error{ NonFiniteInput, NonFiniteResult }!Self {
            return vectorBinary(Self, self, other, .sub);
        }

        /// The componentwise negation. Negating a finite component cannot
        /// leave the format, so this reports input failure only.
        pub fn neg(self: Self) error{NonFiniteInput}!Self {
            return vectorNegated(Self, self);
        }

        /// Every component times `factor`. `error.NonFiniteResult` when a
        /// product is not finite.
        pub fn scale(self: Self, factor: T) error{ NonFiniteInput, NonFiniteResult }!Self {
            return vectorScaled(Self, self, factor);
        }

        /// The scalar product, summed left to right in component order.
        /// `error.NonFiniteResult` when the sum is not finite.
        pub fn dot(self: Self, other: Self) error{ NonFiniteInput, NonFiniteResult }!T {
            try requireFiniteField(self);
            try requireFiniteField(other);
            return vectorDot(Self, self, other);
        }

        /// The squared length, which is the dot product with itself, and which
        /// overflows on a vector large enough that `length` does not.
        /// `error.NonFiniteResult` when the sum is not finite.
        pub fn lengthSquared(self: Self) error{ NonFiniteInput, NonFiniteResult }!T {
            try requireFiniteField(self);
            return vectorLengthSquared(Self, self);
        }

        /// The length, folded through the accepted `hypot` kernel, so a large
        /// vector keeps a finite norm and a subnormal one keeps a representable
        /// one. `error.NonFiniteResult` when the norm itself is past the top
        /// of the range.
        pub fn length(self: Self) error{ NonFiniteInput, NonFiniteResult }!T {
            try requireFiniteField(self);
            return vectorLength(Self, self);
        }

        /// The unit vector along `self`. `error.Degenerate` for a zero vector,
        /// and `error.NonFiniteResult` for a norm past the top of the range.
        pub fn normalize(self: Self) error{ NonFiniteInput, NonFiniteResult, Degenerate }!Self {
            return vectorNormalize(Self, self);
        }

        /// The exterior product `self ^ other`, which is the product four
        /// dimensions has. There is no cross product here: the orthogonal
        /// complement of the plane two vectors span is a whole plane, not a
        /// line, so no single vector completes the product.
        /// `error.NonFiniteResult` when a component of the product is not
        /// finite, which is a computed value rather than a supplied one.
        pub fn wedge(self: Self, other: Self) error{ NonFiniteInput, NonFiniteResult }!Bivector4(T) {
            @setFloatMode(.strict);
            try requireFiniteField(self);
            try requireFiniteField(other);
            const a = self.toArray();
            const b = other.toArray();
            const out = Bivector4(T){
                .xy = a[0] * b[1] - a[1] * b[0],
                .xz = a[0] * b[2] - a[2] * b[0],
                .xw = a[0] * b[3] - a[3] * b[0],
                .yz = a[1] * b[2] - a[2] * b[1],
                .yw = a[1] * b[3] - a[3] * b[1],
                .zw = a[2] * b[3] - a[3] * b[2],
            };
            try requireFiniteBivectorResult(T, out);
            return out;
        }

        /// The squared distance to `point`, from the rounded difference.
        /// `error.NonFiniteResult` when the difference or the sum leaves the
        /// format.
        pub fn distanceSquared(self: Self, point: Point4(T)) error{ NonFiniteInput, NonFiniteResult }!T {
            try requireFiniteField(self);
            try requireFiniteField(point);
            return vectorDistanceSquared(Self, self, toVectorOfPoint(Point4(T), point));
        }

        /// The distance to `point`, the length of the rounded difference.
        /// `error.NonFiniteResult` when the difference leaves the format.
        pub fn distance(self: Self, point: Point4(T)) error{ NonFiniteInput, NonFiniteResult }!T {
            try requireFiniteField(self);
            try requireFiniteField(point);
            return vectorDistance(Self, self, toVectorOfPoint(Point4(T), point));
        }

        /// True when every component is finite. This is the check a caller
        /// needs after writing the fields directly.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(Self, self);
        }

        /// Bit equality, with every NaN equal to every other NaN and `+0`
        /// distinct from `-0`. It compares stored bits and never fails.
        pub fn eq(self: Self, other: Self) bool {
            return eqOf(Self, self, other);
        }
    };
}

/// The exterior product of two four-dimensional vectors, in the six
/// independent components the order `xy, xz, xw, yz, yw, zw` gives.
///
/// A plain value of six `f32`. The fields are public and writable; `init`
/// established that they are finite, and `Vector4d.wedge` is the operation
/// that produces one. Nothing here is normalized: a bivector carries no unit
/// invariant, so a caller may write any six finite coefficients.
pub const Bivector4f = Bivector4(f32);
/// The exterior product of two four-dimensional vectors, in `f64`, under the
/// same field obligation as `Bivector4f`.
pub const Bivector4d = Bivector4(f64);

/// The six-component input gate for a bivector. A bivector names its
/// components after the axis pairs they belong to rather than by an axis, so
/// it cannot share the axis-named helper the point and vector types use.
fn requireFiniteBivector(comptime T: type, value: Bivector4(T)) error{NonFiniteInput}!void {
    if (!isFiniteOf(T, value.xy) or !isFiniteOf(T, value.xz) or
        !isFiniteOf(T, value.xw) or !isFiniteOf(T, value.yz) or
        !isFiniteOf(T, value.yw) or !isFiniteOf(T, value.zw))
    {
        return error.NonFiniteInput;
    }
}

/// The six-component RESULT gate for a bivector the module computed. A wedge
/// whose coefficients overflowed is not bad input, so it is not reported as
/// `error.NonFiniteInput`.
fn requireFiniteBivectorResult(comptime T: type, value: Bivector4(T)) error{NonFiniteResult}!void {
    if (!isFiniteOf(T, value.xy) or !isFiniteOf(T, value.xz) or
        !isFiniteOf(T, value.xw) or !isFiniteOf(T, value.yz) or
        !isFiniteOf(T, value.yw) or !isFiniteOf(T, value.zw))
    {
        return error.NonFiniteResult;
    }
}

fn Bivector4(comptime T: type) type {
    return struct {
        /// The `x y` component. Public and writable; nothing normalizes it.
        xy: T,
        /// The `x z` component.
        xz: T,
        /// The `x w` component.
        xw: T,
        /// The `y z` component.
        yz: T,
        /// The `y w` component.
        yw: T,
        /// The `z w` component.
        zw: T,

        const Self = @This();

        /// Builds a bivector from its six components. `error.NonFiniteInput`
        /// when a component is an infinity or a NaN.
        pub fn init(xy: T, xz: T, xw: T, yz: T, yw: T, zw: T) error{NonFiniteInput}!Self {
            const out = Self{ .xy = xy, .xz = xz, .xw = xw, .yz = yz, .yw = yw, .zw = zw };
            try requireFiniteBivector(T, out);
            return out;
        }

        /// The componentwise negation. Negating a finite component cannot
        /// leave the format, so this is infallible and returns a value that
        /// may hold a NaN the caller wrote into a field.
        pub fn neg(self: Self) Self {
            return .{
                .xy = negRaw(T, self.xy),
                .xz = negRaw(T, self.xz),
                .xw = negRaw(T, self.xw),
                .yz = negRaw(T, self.yz),
                .yw = negRaw(T, self.yw),
                .zw = negRaw(T, self.zw),
            };
        }

        /// True when every component is finite. The bivector's six components
        /// are named after the pairs they belong to rather than by an axis, so
        /// it does not share the axis-named helper the point and vector types
        /// use. This is the check a caller needs after writing the fields
        /// directly.
        pub fn isFinite(self: Self) bool {
            return isFiniteOf(T, self.xy) and isFiniteOf(T, self.xz) and
                isFiniteOf(T, self.xw) and isFiniteOf(T, self.yz) and
                isFiniteOf(T, self.yw) and isFiniteOf(T, self.zw);
        }
    };
}

// ---------------------------------------------------------------------------
// Point
// ---------------------------------------------------------------------------

fn toVectorOfPoint(comptime P: type, point: P) VectorOf(fieldCount(P), componentType(P)) {
    const V = VectorOf(fieldCount(P), componentType(P));
    return fromArrayOf(V, toArrayOf(P, point));
}

fn pointTranslate(comptime P: type, point: P, offset: VectorOf(fieldCount(P), componentType(P))) error{ NonFiniteInput, NonFiniteResult }!P {
    try requireFiniteField(point);
    try requireFiniteField(offset);
    return pointTranslateChecked(P, point, offset);
}

/// The translation over a point and an offset an earlier gate already covered,
/// with the result gated. Every construction in this module ends here, so this
/// is where a value that left the format is stopped before it is returned.
fn pointTranslateChecked(comptime P: type, point: P, offset: VectorOf(fieldCount(P), componentType(P))) error{NonFiniteResult}!P {
    const T = componentType(P);
    const components = toArrayOf(P, point);
    const deltas = toArrayOf(@TypeOf(offset), offset);
    var out: [fieldCount(P)]T = undefined;
    inline for (0..comptime fieldCount(P)) |index| {
        out[index] = addRaw(T, components[index], deltas[index]);
    }
    return fromArrayOf(P, try requireFiniteRow(fieldCount(P), T, out));
}

fn pointAtParameter(
    comptime P: type,
    origin: P,
    direction: VectorOf(fieldCount(P), componentType(P)),
    t: componentType(P),
) error{NonFiniteResult}!P {
    return pointTranslateChecked(P, origin, try vectorScaledChecked(@TypeOf(direction), direction, t));
}

/// A two-dimensional `f32` point.
///
/// A plain value of two `f32`. The fields are public and writable; `init`
/// established that they are finite and every operation below checks them
/// again. A point and a vector are different types, so a translation names the
/// direction as a `Vector2f` and never as a point.
pub const Point2f = Point2(f32);
/// A two-dimensional `f64` point, under the same field obligation as
/// `Point2f`.
pub const Point2d = Point2(f64);
/// A three-dimensional `f32` point, under the same field obligation as
/// `Point2f`.
pub const Point3f = Point3(f32);
/// A three-dimensional `f64` point, under the same field obligation as
/// `Point2f`.
pub const Point3d = Point3(f64);
/// A four-dimensional `f32` point, under the same field obligation as
/// `Point2f`.
pub const Point4f = Point4(f32);
/// A four-dimensional `f64` point, under the same field obligation as
/// `Point2f`.
pub const Point4d = Point4(f64);

fn PointOf(comptime dim: usize, comptime T: type) type {
    return switch (dim) {
        2 => Point2(T),
        3 => Point3(T),
        4 => Point4(T),
        else => @compileError("zmath geometry dimensions are 2, 3, and 4"),
    };
}

fn Point2(comptime T: type) type {
    return struct {
        /// The first component. Public and writable; `init` established that it
        /// is finite and every operation below checks it again.
        x: T,
        /// The second component, under the same obligation as `x`.
        y: T,

        const Self = @This();

        /// Builds a point from named components. `error.NonFiniteInput` when
        /// a component is an infinity or a NaN.
        pub fn init(x: T, y: T) error{NonFiniteInput}!Self {
            const out = Self{ .x = x, .y = y };
            try requireFiniteField(out);
            return out;
        }

        /// Builds a point from `[x, y]`, under the same gate as `init`.
        pub fn fromArray(a: [2]T) error{NonFiniteInput}!Self {
            return init(a[0], a[1]);
        }

        /// The components as `[x, y]`, copied out.
        pub fn toArray(self: Self) [2]T {
            return toArrayOf(Self, self);
        }

        /// The point `offset` from `self`, where `offset` is a direction and
        /// not a position. `error.NonFiniteResult` when a component of the
        /// translation is not finite.
        pub fn add(self: Self, offset: Vector2(T)) error{ NonFiniteInput, NonFiniteResult }!Self {
            return pointTranslate(Self, self, offset);
        }

        /// The point `offset` back from `self`. The offset is negated first,
        /// after it is gated, so a caller cannot reach the negation with a
        /// value the constructor would have refused.
        /// `error.NonFiniteResult` when a component of the translation is not
        /// finite.
        pub fn sub(self: Self, offset: Vector2(T)) error{ NonFiniteInput, NonFiniteResult }!Self {
            return pointTranslate(Self, self, try vectorNegated(Vector2(T), offset));
        }

        /// The same coordinates read as a vector. The read never fails, so a
        /// point a caller poisoned still converts; every operation that
        /// consumes the result gates it again.
        pub fn toVector(self: Self) Vector2(T) {
            return toVectorOfPoint(Self, self);
        }

        /// The sign of the oriented area of the triangle `a`, `b`, `c`: it is
        /// `negative`, `zero`, or `positive` as `c` lies to the right of, on,
        /// or to the left of the directed line from `a` to `b`.
        ///
        /// `adaptive_exact` forms no floating-point intermediate and answers
        /// for any finite inputs. `fast` answers from the plain evaluation,
        /// including when that evaluation overflows to an infinity, which
        /// carries a sign. `error.NonFiniteResult` names the one case `fast`
        /// cannot answer: an evaluation that is a NaN, which has no sign to
        /// report.
        pub fn orientation(a: Self, b: Self, c: Self, policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!Sign {
            try requireFiniteField(a);
            try requireFiniteField(b);
            try requireFiniteField(c);
            const av = a.toArray();
            const bv = b.toArray();
            const cv = c.toArray();
            return switch (policy) {
                .fast => fastSign(T, (bv[0] - av[0]) * (cv[1] - av[1]) -
                    (bv[1] - av[1]) * (cv[0] - av[0])),
                .adaptive_exact => exactOrient2(T, av, bv, cv),
            };
        }

        /// Which side of the directed line from `self` to `other` the `point`
        /// lies on. A two-dimensional line has two sides; a line in three or
        /// four dimensions does not, so this exists here only. The errors and
        /// the policy mean what they mean in `orientation`.
        pub fn side(self: Self, other: Self, point: Self, policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!Sign {
            return orientation(self, other, point, policy);
        }

        /// True when both components are finite. This is the check a caller
        /// needs after writing the fields directly.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(Self, self);
        }

        /// Bit equality, with every NaN equal to every other NaN and `+0`
        /// distinct from `-0`. It compares stored bits and never fails.
        pub fn eq(self: Self, other: Self) bool {
            return eqOf(Self, self, other);
        }
    };
}

fn Point3(comptime T: type) type {
    return struct {
        /// The first component. Public and writable; `init` established that it
        /// is finite and every operation below checks it again.
        x: T,
        /// The second component, under the same obligation as `x`.
        y: T,
        /// The third component, under the same obligation as `x`.
        z: T,

        const Self = @This();

        /// Builds a point from named components. `error.NonFiniteInput` when
        /// a component is an infinity or a NaN.
        pub fn init(x: T, y: T, z: T) error{NonFiniteInput}!Self {
            const out = Self{ .x = x, .y = y, .z = z };
            try requireFiniteField(out);
            return out;
        }

        /// Builds a point from `[x, y, z]`, under the same gate as `init`.
        pub fn fromArray(a: [3]T) error{NonFiniteInput}!Self {
            return init(a[0], a[1], a[2]);
        }

        /// The components as `[x, y, z]`, copied out.
        pub fn toArray(self: Self) [3]T {
            return toArrayOf(Self, self);
        }

        /// The point `offset` from `self`, where `offset` is a direction and
        /// not a position. `error.NonFiniteResult` when a component of the
        /// translation is not finite.
        pub fn add(self: Self, offset: Vector3(T)) error{ NonFiniteInput, NonFiniteResult }!Self {
            return pointTranslate(Self, self, offset);
        }

        /// The point `offset` back from `self`. The offset is negated first,
        /// after it is gated, so a caller cannot reach the negation with a
        /// value the constructor would have refused.
        /// `error.NonFiniteResult` when a component of the translation is not
        /// finite.
        pub fn sub(self: Self, offset: Vector3(T)) error{ NonFiniteInput, NonFiniteResult }!Self {
            return pointTranslate(Self, self, try vectorNegated(Vector3(T), offset));
        }

        /// The same coordinates read as a vector. The read never fails, so a
        /// point a caller poisoned still converts; every operation that
        /// consumes the result gates it again.
        pub fn toVector(self: Self) Vector3(T) {
            return toVectorOfPoint(Self, self);
        }

        /// The sign of the signed volume of the tetrahedron `a`, `b`, `c`,
        /// `d`, which is `dot(cross(b - a, c - a), d - a)`. It is zero
        /// exactly when the four points are coplanar, and that is the decision
        /// a plain determinant evaluation gets wrong when the terms cancel.
        ///
        /// The two policies differ in the errors they can return, exactly as
        /// `Point2d.orientation` documents: `adaptive_exact` needs no
        /// floating-point intermediate, and `fast` reports
        /// `error.NonFiniteResult` only when its evaluation is a NaN.
        pub fn orientation(a: Self, b: Self, c: Self, d: Self, policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!Sign {
            try requireFiniteField(a);
            try requireFiniteField(b);
            try requireFiniteField(c);
            try requireFiniteField(d);
            const av = a.toArray();
            const bv = b.toArray();
            const cv = c.toArray();
            const dv = d.toArray();
            return switch (policy) {
                .fast => fastSign(T, crossDotFast(T, av, bv, cv, dv)),
                .adaptive_exact => exactOrient3(T, av, bv, cv, dv),
            };
        }

        /// True when every component is finite. This is the check a caller
        /// needs after writing the fields directly.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(Self, self);
        }

        /// Bit equality, with every NaN equal to every other NaN and `+0`
        /// distinct from `-0`. It compares stored bits and never fails.
        pub fn eq(self: Self, other: Self) bool {
            return eqOf(Self, self, other);
        }
    };
}

/// `dot(cross(b - a, c - a), d - a)`, summed left to right. The grouping is
/// the definition of `fast` for this predicate.
fn crossDotFast(comptime T: type, a: [3]T, b: [3]T, c: [3]T, d: [3]T) T {
    @setFloatMode(.strict);
    const n = cross3Raw(T, .{
        b[0] - a[0],
        b[1] - a[1],
        b[2] - a[2],
    }, .{
        c[0] - a[0],
        c[1] - a[1],
        c[2] - a[2],
    });
    const w = [3]T{ d[0] - a[0], d[1] - a[1], d[2] - a[2] };
    return n[0] * w[0] + n[1] * w[1] + n[2] * w[2];
}

fn Point4(comptime T: type) type {
    return struct {
        /// The first component. Public and writable; `init` established that it
        /// is finite and every operation below checks it again.
        x: T,
        /// The second component, under the same obligation as `x`.
        y: T,
        /// The third component, under the same obligation as `x`.
        z: T,
        /// The fourth component, under the same obligation as `x`.
        w: T,

        const Self = @This();

        /// Builds a point from named components. `error.NonFiniteInput` when
        /// a component is an infinity or a NaN.
        pub fn init(x: T, y: T, z: T, w: T) error{NonFiniteInput}!Self {
            const out = Self{ .x = x, .y = y, .z = z, .w = w };
            try requireFiniteField(out);
            return out;
        }

        /// Builds a point from `[x, y, z, w]`, under the same gate as `init`.
        pub fn fromArray(a: [4]T) error{NonFiniteInput}!Self {
            return init(a[0], a[1], a[2], a[3]);
        }

        /// The components as `[x, y, z, w]`, copied out.
        pub fn toArray(self: Self) [4]T {
            return toArrayOf(Self, self);
        }

        /// The point `offset` from `self`, where `offset` is a direction and
        /// not a position. `error.NonFiniteResult` when a component of the
        /// translation is not finite.
        pub fn add(self: Self, offset: Vector4(T)) error{ NonFiniteInput, NonFiniteResult }!Self {
            return pointTranslate(Self, self, offset);
        }

        /// The point `offset` back from `self`. The offset is negated first,
        /// after it is gated, so a caller cannot reach the negation with a
        /// value the constructor would have refused.
        /// `error.NonFiniteResult` when a component of the translation is not
        /// finite.
        pub fn sub(self: Self, offset: Vector4(T)) error{ NonFiniteInput, NonFiniteResult }!Self {
            return pointTranslate(Self, self, try vectorNegated(Vector4(T), offset));
        }

        /// The same coordinates read as a vector. The read never fails, so a
        /// point a caller poisoned still converts; every operation that
        /// consumes the result gates it again.
        pub fn toVector(self: Self) Vector4(T) {
            return toVectorOfPoint(Self, self);
        }

        /// The sign of the four-by-four determinant whose rows are the four
        /// points. It is zero exactly when the four position vectors are
        /// linearly dependent, and its magnitude is the oriented
        /// four-dimensional volume they span from the origin. This is the
        /// degree-four predicate that sizes the exact arithmetic: its value
        /// needs every limb of the array.
        ///
        /// It is not the coplanarity test. Four four-dimensional points are
        /// coplanar through a different form, a rank test on three difference
        /// vectors, and the affine intersections use that one privately.
        ///
        /// The two policies differ in the errors they can return, exactly as
        /// `Point2d.orientation` documents: `adaptive_exact` needs no
        /// floating-point intermediate, and `fast` reports
        /// `error.NonFiniteResult` only when its evaluation is a NaN.
        pub fn orientation(a: Self, b: Self, c: Self, d: Self, policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!Sign {
            try requireFiniteField(a);
            try requireFiniteField(b);
            try requireFiniteField(c);
            try requireFiniteField(d);
            return orient4Sign(T, a.toArray(), b.toArray(), c.toArray(), d.toArray(), policy);
        }

        /// True when every component is finite. This is the check a caller
        /// needs after writing the fields directly.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(Self, self);
        }

        /// Bit equality, with every NaN equal to every other NaN and `+0`
        /// distinct from `-0`. It compares stored bits and never fails.
        pub fn eq(self: Self, other: Self) bool {
            return eqOf(Self, self, other);
        }
    };
}

// ---------------------------------------------------------------------------
// Line, ray, segment
// ---------------------------------------------------------------------------

fn lineBuild(comptime P: type, point: P, direction: VectorOf(fieldCount(P), componentType(P))) error{ NonFiniteInput, NonFiniteResult, Degenerate }!struct { origin: P, direction: VectorOf(fieldCount(P), componentType(P)) } {
    try requireFiniteField(point);
    const unit = try vectorNormalize(@TypeOf(direction), direction);
    return .{ .origin = point, .direction = unit };
}

/// The orthogonal projection of `point` onto the line through `origin` along
/// `direction`. All three values are gated: the receiver's fields are public
/// and writable, so the projection is exactly where a poisoned direction would
/// otherwise become a poisoned parameter.
fn projectOntoLine(comptime P: type, origin: P, direction: VectorOf(fieldCount(P), componentType(P)), point: P) error{ NonFiniteInput, NonFiniteResult }!P {
    try requireFiniteField(origin);
    try requireFiniteField(direction);
    try requireFiniteField(point);
    const V = @TypeOf(direction);
    const difference = try vectorBinaryChecked(V, toVectorOfPoint(P, point), toVectorOfPoint(P, origin), .sub);
    const t = try vectorDot(V, difference, direction);
    return pointTranslateChecked(P, origin, try vectorScaledChecked(V, direction, t));
}

/// A two-dimensional `f32` line.
///
/// A plain value of a point and a direction. `init` established that `origin`
/// is finite and that `direction` is a UNIT vector, and that unit invariant
/// is the caller's to keep: the fields are public, and every operation below
/// gates the finiteness of both without re-checking the length. Writing a
/// non-unit `direction` is legal and changes what the line means, so a caller
/// that rewrites the field must keep the direction unit itself.
pub const Line2f = LineOf(2, f32);
/// A two-dimensional `f64` line, under the same field obligation as `Line2f`.
pub const Line2d = LineOf(2, f64);
/// A three-dimensional `f32` line, under the same field obligation as `Line2f`.
pub const Line3f = LineOf(3, f32);
/// A three-dimensional `f64` line, under the same field obligation as `Line2f`.
pub const Line3d = LineOf(3, f64);
/// A four-dimensional `f32` line, under the same field obligation as `Line2f`.
pub const Line4f = LineOf(4, f32);
/// A four-dimensional `f64` line, under the same field obligation as `Line2f`.
pub const Line4d = LineOf(4, f64);

fn LineOf(comptime dim: usize, comptime T: type) type {
    const P = PointOf(dim, T);
    const V = VectorOf(dim, T);
    return struct {
        /// A point on the line. Public and writable; `init` established that it
        /// is finite.
        origin: P,
        /// The unit direction. Public and writable; the caller keeps it unit.
        direction: V,

        const Self = @This();

        /// Builds the line through `point` along `direction`, which is
        /// normalized. `error.Degenerate` for a zero direction,
        /// `error.NonFiniteInput` for a non-finite point or direction, and
        /// `error.NonFiniteResult` for a direction whose norm is past the top
        /// of the range.
        pub fn init(point: P, direction: V) error{ NonFiniteInput, NonFiniteResult, Degenerate }!Self {
            const built = try lineBuild(P, point, direction);
            return .{ .origin = built.origin, .direction = built.direction };
        }

        /// The point at `t`. There is no clamp: a line runs in both
        /// directions, so a negative `t` is a point behind `origin` and not an
        /// error. `error.NonFiniteResult` when `t * direction` or the
        /// translation leaves the format.
        pub fn at(self: Self, t: T) error{ NonFiniteInput, NonFiniteResult }!P {
            try requireFiniteField(self);
            try requireFiniteScalar(T, t);
            return pointAtParameter(P, self.origin, self.direction, t);
        }

        /// The point of the line closest to `point`, which is the orthogonal
        /// projection. `error.NonFiniteResult` when the offset, the parameter,
        /// or the translation leaves the format.
        pub fn closestPointTo(self: Self, point: P) error{ NonFiniteInput, NonFiniteResult }!P {
            return projectOntoLine(P, self.origin, self.direction, point);
        }

        /// The set of points common to two lines. Two lines in two or more
        /// dimensions share at most one point unless they are the same line,
        /// so the answer is `none`, `point`, or `line`. `line` is the whole
        /// two-sided line, and it appears only for two lines that are one
        /// line; a line against a ray or a segment reports `.contained`
        /// instead. `error.NonFiniteResult` when a derived point, parameter, or
        /// constructed result leaves the format.
        pub fn intersect(self: Self, other: Self, policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!IntersectionOf(dim, T) {
            try requireFiniteField(self);
            try requireFiniteField(other);
            return affineIntersectDispatch(dim, T, self.origin.toArray(), try pointAtArrayChecked(dim, T, self.origin.toArray(), self.direction.toArray(), 1), -infinityOf(T), infinityOf(T), other.origin.toArray(), try pointAtArrayChecked(dim, T, other.origin.toArray(), other.direction.toArray(), 1), -infinityOf(T), infinityOf(T), policy);
        }

        /// True when the origin and the direction are finite. It says nothing
        /// about the direction being unit, which is the caller's invariant to
        /// keep rather than one this module re-checks.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(P, self.origin) and isFiniteQuiet(V, self.direction);
        }
    };
}

/// A two-dimensional `f32` ray.
///
/// A plain value of an end and a direction. `init` established that `origin`
/// is finite and that `direction` is a UNIT vector; as with a line, the unit
/// invariant is the caller's to keep and the fields are public. What
/// separates a ray from a line is the parameter range `[0, +infinity)`, and
/// that is what decides whether two identical rays report `coincident` or a
/// line reports `line`.
pub const Ray2f = RayOf(2, f32);
/// A two-dimensional `f64` ray, under the same field obligation as `Ray2f`.
pub const Ray2d = RayOf(2, f64);
/// A three-dimensional `f32` ray, under the same field obligation as `Ray2f`.
pub const Ray3f = RayOf(3, f32);
/// A three-dimensional `f64` ray, under the same field obligation as `Ray2f`.
pub const Ray3d = RayOf(3, f64);
/// A four-dimensional `f32` ray, under the same field obligation as `Ray2f`.
pub const Ray4f = RayOf(4, f32);
/// A four-dimensional `f64` ray, under the same field obligation as `Ray2f`.
pub const Ray4d = RayOf(4, f64);

fn RayOf(comptime dim: usize, comptime T: type) type {
    const P = PointOf(dim, T);
    const V = VectorOf(dim, T);
    return struct {
        /// The end the ray starts from. Public and writable; `init` established
        /// that it is finite.
        origin: P,
        /// The unit direction the ray runs in. Public and writable; the caller
        /// keeps it unit.
        direction: V,

        const Self = @This();

        /// Builds the ray from `origin` along `direction`, which is
        /// normalized. `error.Degenerate` for a zero direction,
        /// `error.NonFiniteInput` for a non-finite origin or direction, and
        /// `error.NonFiniteResult` for a direction whose norm is past the top
        /// of the range.
        pub fn init(origin: P, direction: V) error{ NonFiniteInput, NonFiniteResult, Degenerate }!Self {
            const built = try lineBuild(P, origin, direction);
            return .{ .origin = built.origin, .direction = built.direction };
        }

        /// The point at `t`, clamped at `t = 0` so the ray never extends
        /// behind the end it starts from. `error.NonFiniteResult` when
        /// `t * direction` or the translation leaves the format.
        pub fn at(self: Self, t: T) error{ NonFiniteInput, NonFiniteResult }!P {
            try requireFiniteField(self);
            try requireFiniteScalar(T, t);
            if (t <= 0) return self.origin;
            return pointAtParameter(P, self.origin, self.direction, t);
        }

        /// The set of points common to two rays. Two rays share at most one
        /// point unless they leave the same end in the same direction, which
        /// is `contained` or `coincident`. An identical pair is
        /// `coincident`: a ray is a half-line, and the common set of a ray
        /// with its own copy is that half-line, never a whole `line`.
        /// `error.NonFiniteResult` when a derived point, parameter, or
        /// constructed result leaves the format.
        pub fn intersect(self: Self, other: Self, policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!IntersectionOf(dim, T) {
            try requireFiniteField(self);
            try requireFiniteField(other);
            return affineIntersectDispatch(dim, T, self.origin.toArray(), try pointAtArrayChecked(dim, T, self.origin.toArray(), self.direction.toArray(), 1), 0, infinityOf(T), other.origin.toArray(), try pointAtArrayChecked(dim, T, other.origin.toArray(), other.direction.toArray(), 1), 0, infinityOf(T), policy);
        }

        /// The set of points common to this ray and a line. A ray lying in the
        /// line reports `.contained = .first`, because the ray is the whole
        /// common set. `error.NonFiniteResult` as in `intersect`.
        pub fn intersectLine(self: Self, line: LineOf(dim, T), policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!IntersectionOf(dim, T) {
            try requireFiniteField(self);
            try requireFiniteField(line);
            return affineIntersectDispatch(dim, T, self.origin.toArray(), try pointAtArrayChecked(dim, T, self.origin.toArray(), self.direction.toArray(), 1), 0, infinityOf(T), line.origin.toArray(), try pointAtArrayChecked(dim, T, line.origin.toArray(), line.direction.toArray(), 1), -infinityOf(T), infinityOf(T), policy);
        }

        /// The set of points common to this ray and a segment.
        /// `error.NonFiniteResult` as in `intersect`.
        pub fn intersectSegment(self: Self, segment: SegmentOf(dim, T), policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!IntersectionOf(dim, T) {
            try requireFiniteField(self);
            try requireFiniteField(segment);
            return affineIntersectDispatch(dim, T, self.origin.toArray(), try pointAtArrayChecked(dim, T, self.origin.toArray(), self.direction.toArray(), 1), 0, infinityOf(T), segment.a.toArray(), segment.b.toArray(), 0, 1, policy);
        }

        /// True when the origin and the direction are finite. It says nothing
        /// about the direction being unit, which is the caller's invariant to
        /// keep rather than one this module re-checks.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(P, self.origin) and isFiniteQuiet(V, self.direction);
        }
    };
}

/// A two-dimensional `f32` segment.
///
/// A plain value of two ends. `init` established that both are finite; the
/// fields are public and writable and every operation below gates them again.
/// There is no unit invariant and no ordering invariant, so the ends may be
/// written in either direction and may coincide, which gives the degenerate
/// segment that is a single point.
pub const Segment2f = SegmentOf(2, f32);
/// A two-dimensional `f64` segment, under the same field obligation as
/// `Segment2f`.
pub const Segment2d = SegmentOf(2, f64);
/// A three-dimensional `f32` segment, under the same field obligation as
/// `Segment2f`.
pub const Segment3f = SegmentOf(3, f32);
/// A three-dimensional `f64` segment, under the same field obligation as
/// `Segment2f`.
pub const Segment3d = SegmentOf(3, f64);
/// A four-dimensional `f32` segment, under the same field obligation as
/// `Segment2f`.
pub const Segment4f = SegmentOf(4, f32);
/// A four-dimensional `f64` segment, under the same field obligation as
/// `Segment2f`.
pub const Segment4d = SegmentOf(4, f64);

fn SegmentOf(comptime dim: usize, comptime T: type) type {
    const P = PointOf(dim, T);
    const V = VectorOf(dim, T);
    return struct {
        /// The end at `t = 0`. Public and writable; `init` established that it
        /// is finite.
        a: P,
        /// The end at `t = 1`, under the same obligation as `a`.
        b: P,

        const Self = @This();

        /// Builds the segment from `a` to `b`. The two ends may coincide,
        /// which gives the degenerate segment that is a single point.
        /// `error.NonFiniteInput` when an end is an infinity or a NaN. Nothing
        /// is computed, so this reports no other failure.
        pub fn init(a: P, b: P) error{NonFiniteInput}!Self {
            try requireFiniteField(a);
            try requireFiniteField(b);
            return .{ .a = a, .b = b };
        }

        /// The point at `t`, clamped to `[0, 1]`, so the segment never
        /// extends past either end. The whole receiver is gated before the
        /// clamp, so an end a caller poisoned is refused even at `t = 0` or
        /// `t = 1` where the answer would be that end itself.
        /// `error.NonFiniteResult` when the difference or the interpolation
        /// leaves the format.
        pub fn at(self: Self, t: T) error{ NonFiniteInput, NonFiniteResult }!P {
            try requireFiniteField(self);
            try requireFiniteScalar(T, t);
            const clamped = clampUnitOf(T, t);
            if (clamped == 0) return self.a;
            if (clamped == 1) return self.b;
            const difference = try vectorBinaryChecked(V, toVectorOfPoint(P, self.b), toVectorOfPoint(P, self.a), .sub);
            return pointTranslateChecked(P, self.a, try vectorScaledChecked(V, difference, clamped));
        }

        /// The distance between the two ends, folded through the accepted
        /// `hypot` kernel so it does not overflow the way `lengthSquared`
        /// does. `error.NonFiniteResult` when the difference itself leaves the
        /// format, which two ends at opposite ends of the range reach.
        pub fn length(self: Self) error{ NonFiniteInput, NonFiniteResult }!T {
            try requireFiniteField(self);
            return vectorDistance(V, toVectorOfPoint(P, self.a), toVectorOfPoint(P, self.b));
        }

        /// The point of the segment closest to `point`, which is the
        /// orthogonal projection clamped to the segment. A degenerate segment
        /// answers its own single point. `error.NonFiniteResult` when the
        /// squared length of the segment, the offset, or the translation
        /// leaves the format.
        pub fn closestPointTo(self: Self, point: P) error{ NonFiniteInput, NonFiniteResult }!P {
            try requireFiniteField(self);
            try requireFiniteField(point);
            const from = toVectorOfPoint(P, self.a);
            const to = toVectorOfPoint(P, self.b);
            const difference = try vectorBinaryChecked(V, to, from, .sub);
            const denominator = try vectorDot(V, difference, difference);
            if (denominator == 0) return self.a;
            const offset = try vectorBinaryChecked(V, toVectorOfPoint(P, point), from, .sub);
            const t = clampUnitOf(T, (try vectorDot(V, offset, difference)) / denominator);
            if (t == 0) return self.a;
            if (t == 1) return self.b;
            return pointTranslateChecked(P, self.a, try vectorScaledChecked(V, difference, t));
        }

        /// The set of points common to two segments. The tags are `none` for
        /// disjoint segments, `point` for a crossing or a touching end,
        /// `segment` for a proper collinear overlap, `contained` when one
        /// segment lies inside the other, and `coincident` when the two are
        /// the same segment. `error.NonFiniteResult` when a derived point,
        /// parameter, or constructed result leaves the format.
        pub fn intersect(self: Self, other: Self, policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!IntersectionOf(dim, T) {
            try requireFiniteField(self);
            try requireFiniteField(other);
            return affineIntersectDispatch(dim, T, self.a.toArray(), self.b.toArray(), 0, 1, other.a.toArray(), other.b.toArray(), 0, 1, policy);
        }

        /// The set of points common to this segment and a ray. A segment lying
        /// wholly on the ray reports `.contained = .first`, because the
        /// segment is the whole common set. `error.NonFiniteResult` as in
        /// `intersect`.
        pub fn intersectRay(self: Self, ray: RayOf(dim, T), policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!IntersectionOf(dim, T) {
            try requireFiniteField(self);
            try requireFiniteField(ray);
            return affineIntersectDispatch(dim, T, self.a.toArray(), self.b.toArray(), 0, 1, ray.origin.toArray(), try pointAtArrayChecked(dim, T, ray.origin.toArray(), ray.direction.toArray(), 1), 0, infinityOf(T), policy);
        }

        /// The set of points common to this segment and a line.
        /// `error.NonFiniteResult` as in `intersect`.
        pub fn intersectLine(self: Self, line: LineOf(dim, T), policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!IntersectionOf(dim, T) {
            try requireFiniteField(self);
            try requireFiniteField(line);
            return affineIntersectDispatch(dim, T, self.a.toArray(), self.b.toArray(), 0, 1, line.origin.toArray(), try pointAtArrayChecked(dim, T, line.origin.toArray(), line.direction.toArray(), 1), -infinityOf(T), infinityOf(T), policy);
        }

        /// True when both ends are finite. This is the check a caller needs
        /// after writing the fields directly.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(P, self.a) and isFiniteQuiet(P, self.b);
        }
    };
}

// ---------------------------------------------------------------------------
// Array helpers shared by the geometric construction paths
//
// There are two kinds of helper here, and the difference is the contract
// rather than the arithmetic. The PLAIN ones compute and return, and the
// `fast` sign predicates are defined over exactly those evaluations, so they
// must not stop a policy from reporting the sign of its own arithmetic. The
// CHECKED ones are the same operations with the result gated, and every path
// that BUILDS a value uses those, because a constructed point, direction, or
// radius is a value the caller is about to hold.
// ---------------------------------------------------------------------------

fn arrayLength(comptime dim: usize, comptime T: type, a: [dim]T) T {
    var accumulator: T = 0;
    inline for (0..dim) |index| {
        accumulator = hypotOf(T, accumulator, a[index]);
    }
    return accumulator;
}

fn arrayLengthChecked(comptime dim: usize, comptime T: type, a: [dim]T) error{NonFiniteResult}!T {
    return requireFiniteResult(T, arrayLength(dim, T, a));
}

/// Scales `a` to unit length; a zero vector is unchanged. Callers must handle
/// degeneracy before using this as a direction. An overflowing norm returns
/// `error.NonFiniteResult`, not a silently zeroed direction.
fn arrayNormalizedChecked(comptime dim: usize, comptime T: type, a: [dim]T) error{NonFiniteResult}![dim]T {
    const length = try arrayLengthChecked(dim, T, a);
    if (length == 0) return a;
    var out: [dim]T = undefined;
    if (length < std.math.floatMin(T)) {
        var largest: T = 0;
        inline for (0..dim) |index| largest = maxOf(T, largest, absOf(T, a[index]));
        inline for (0..dim) |index| out[index] = divRaw(T, a[index], largest);
        const scaled_length = try arrayLengthChecked(dim, T, out);
        inline for (0..dim) |index| out[index] = divRaw(T, out[index], scaled_length);
    } else {
        // Keep one division, as in the original construction graph.
        inline for (0..dim) |index| out[index] = divRaw(T, a[index], length);
    }
    return requireFiniteRow(dim, T, out);
}

fn arrayScaleChecked(comptime dim: usize, comptime T: type, a: [dim]T, factor: T) error{NonFiniteResult}![dim]T {
    var out: [dim]T = undefined;
    inline for (0..dim) |index| {
        out[index] = scaleRaw(T, a[index], factor);
    }
    return requireFiniteRow(dim, T, out);
}

fn arrayAddChecked(comptime dim: usize, comptime T: type, a: [dim]T, b: [dim]T) error{NonFiniteResult}![dim]T {
    var out: [dim]T = undefined;
    inline for (0..dim) |index| {
        out[index] = addRaw(T, a[index], b[index]);
    }
    return requireFiniteRow(dim, T, out);
}

fn arraySub(comptime dim: usize, comptime T: type, a: [dim]T, b: [dim]T) [dim]T {
    var out: [dim]T = undefined;
    inline for (0..dim) |index| {
        out[index] = subRaw(T, a[index], b[index]);
    }
    return out;
}

fn arraySubChecked(comptime dim: usize, comptime T: type, a: [dim]T, b: [dim]T) error{NonFiniteResult}![dim]T {
    return requireFiniteRow(dim, T, arraySub(dim, T, a, b));
}

fn arrayDot(comptime dim: usize, comptime T: type, a: [dim]T, b: [dim]T) T {
    return dotOfArray(dim, T, a, b);
}

fn arrayDotChecked(comptime dim: usize, comptime T: type, a: [dim]T, b: [dim]T) error{NonFiniteResult}!T {
    return requireFiniteResult(T, arrayDot(dim, T, a, b));
}

/// A unit vector orthogonal to `v`, built by projecting out the coordinate
/// axis least aligned with `v` and normalizing. Choosing the least aligned
/// axis is what keeps the residual far from zero, so the normalization is
/// well conditioned. The caller uses this to complete a frame, never to
/// decide a sign, so it takes no policy; it does gate its own results,
/// because a frame with an infinite direction is not a frame.
fn anyPerpendicularChecked(comptime dim: usize, comptime T: type, v: [dim]T) error{NonFiniteResult}![dim]T {
    var best: usize = 0;
    var best_value: T = absOf(T, v[0]);
    inline for (1..dim) |index| {
        const value = absOf(T, v[index]);
        if (value < best_value) {
            best_value = value;
            best = index;
        }
    }
    var seed: [dim]T = undefined;
    inline for (0..dim) |index| {
        seed[index] = 0;
    }
    seed[best] = 1;
    const along = try arrayDotChecked(dim, T, seed, v);
    const squared = try arrayDotChecked(dim, T, v, v);
    return arrayNormalizedChecked(dim, T, try arraySubChecked(dim, T, seed, try arrayScaleChecked(dim, T, v, divRaw(T, along, squared))));
}

/// `v` with the component along each vector of `against` removed, in the
/// order given, then normalized. A slice rather than a fixed-length array
/// because the number of vectors depends on the dimension, and an `inline
/// for` over a slice is not what a comptime frame needs.
fn orthogonalizedChecked(comptime dim: usize, comptime T: type, v: [dim]T, against: []const [dim]T) error{NonFiniteResult}![dim]T {
    var current = v;
    for (against) |basis| {
        const projected = try arrayDotChecked(dim, T, current, basis);
        current = try arraySubChecked(dim, T, current, try arrayScaleChecked(dim, T, basis, projected));
    }
    return arrayNormalizedChecked(dim, T, current);
}

/// The point one parameter step along `direction` from `origin`, gated. Every
/// construction in this module that reports a point goes through here, so this
/// is the gate on the object a caller finally receives.
fn pointAtArrayChecked(comptime dim: usize, comptime T: type, origin: [dim]T, direction: [dim]T, t: T) error{NonFiniteResult}![dim]T {
    return arrayAddChecked(dim, T, origin, try arrayScaleChecked(dim, T, direction, t));
}

// ---------------------------------------------------------------------------
// Degree-two predicate forms shared by the balls and the flats
// ---------------------------------------------------------------------------

/// The sign of `dot(point - centre, point - centre) - radius * radius`,
/// degree two: negative inside the ball, zero on its boundary, positive
/// outside. Both policies read the ORIGINAL centre and point, so they differ
/// only in where the difference is formed. `error.NonFiniteResult` is the
/// `fast` policy saying its own evaluation came out a NaN.
fn ballContainsSign(comptime dim: usize, comptime T: type, center: [dim]T, radius: T, point: [dim]T, policy: PredicatePolicy) error{NonFiniteResult}!Sign {
    return switch (policy) {
        .fast => blk: {
            const v = arraySub(dim, T, point, center);
            break :blk try fastSign(T, subRaw(T, dotOfArray(dim, T, v, v), scaleRaw(T, radius, radius)));
        },
        .adaptive_exact => blk: {
            const centre = exactRow(dim, T, center);
            const probe = exactRow(dim, T, point);
            var deltas: [dim]Big = undefined;
            // SAFETY: Big.fromFloat writes metadata and active limbs before comparison.
            var scale: Big = undefined;
            Big.fromFloat(&scale, widen(T, radius));
            exactDeltasInto(dim, &deltas, &centre, &probe);
            break :blk squaredNormMinusSquareExact(dim, &deltas, &scale);
        },
    };
}

/// Whether each end of a segment lies inside a ball, over the ORIGINAL stored
/// coordinates. The centre and the radius are converted once and both ends
/// read them. `error.NonFiniteResult` is the `fast` policy saying its own
/// evaluation came out a NaN.
fn segmentInsideSigns(
    comptime dim: usize,
    comptime T: type,
    center: [dim]T,
    radius: T,
    a: [dim]T,
    b: [dim]T,
    policy: PredicatePolicy,
) error{NonFiniteResult}!EndSigns {
    return switch (policy) {
        .fast => blk: {
            const from_a = arraySub(dim, T, a, center);
            const from_b = arraySub(dim, T, b, center);
            break :blk .{
                .start = try fastSign(T, subRaw(T, dotOfArray(dim, T, from_a, from_a), scaleRaw(T, radius, radius))),
                .end = try fastSign(T, subRaw(T, dotOfArray(dim, T, from_b, from_b), scaleRaw(T, radius, radius))),
            };
        },
        .adaptive_exact => blk: {
            const anchor = exactRow(dim, T, center);
            // SAFETY: Big.fromFloat writes metadata and active limbs before comparison.
            var k: Big = undefined;
            Big.fromFloat(&k, widen(T, radius));
            var deltas: [dim]Big = undefined;
            const first = exactRow(dim, T, a);
            exactDeltasInto(dim, &deltas, &anchor, &first);
            const start = squaredNormMinusSquareExact(dim, &deltas, &k);
            const second = exactRow(dim, T, b);
            exactDeltasInto(dim, &deltas, &anchor, &second);
            break :blk .{ .start = start, .end = squaredNormMinusSquareExact(dim, &deltas, &k) };
        },
    };
}

/// The two tangency comparisons a pair of balls turns on, over the ORIGINAL
/// stored coordinates: the separation against the SUM of the radii, and
/// against their DIFFERENCE.
///
/// Both combinations and the separation come from the same four stored
/// numbers, so the exact path converts them once and reads both answers off
/// that one anchor. The radius combinations are exact values of their own,
/// because a rounded sum is the whole difference between two balls that cross
/// and two that touch. `error.NonFiniteResult` is the `fast` policy saying its
/// own evaluation came out a NaN.
fn ballTangencySigns(
    comptime dim: usize,
    comptime T: type,
    first_center: [dim]T,
    first_radius: T,
    second_center: [dim]T,
    second_radius: T,
    policy: PredicatePolicy,
) error{NonFiniteResult}!TangencySigns {
    return switch (policy) {
        .fast => blk: {
            const separation = arraySub(dim, T, second_center, first_center);
            const sum = addRaw(T, first_radius, second_radius);
            const difference = subRaw(T, first_radius, second_radius);
            break :blk .{
                .outer = try fastSign(T, subRaw(T, dotOfArray(dim, T, separation, separation), scaleRaw(T, sum, sum))),
                .inner = try fastSign(T, subRaw(T, dotOfArray(dim, T, separation, separation), scaleRaw(T, difference, difference))),
            };
        },
        .adaptive_exact => blk: {
            const first_row = exactRow(dim, T, first_center);
            const second_row = exactRow(dim, T, second_center);
            var separation: [dim]Big = undefined;
            // SAFETY: exactNormInto writes metadata and active limbs before comparison.
            var squared: Big = undefined;
            // SAFETY: exactRadiusPairInto writes both radii before squaring.
            var pair: ExactRadiusPair = undefined;
            // SAFETY: Big.mulInto writes metadata and active limbs before subtraction.
            var square: Big = undefined;
            exactDeltasInto(dim, &separation, &first_row, &second_row);
            exactNormInto(dim, &squared, &separation);
            exactRadiusPairInto(&pair, T, first_radius, second_radius);
            Big.mulInto(&square, &pair.sum, &pair.sum);
            Big.subInto(&square, &squared, &square);
            const outer = Big.sign(&square);
            Big.mulInto(&square, &pair.difference, &pair.difference);
            Big.subInto(&square, &squared, &square);
            break :blk .{ .outer = outer, .inner = Big.sign(&square) };
        },
    };
}

/// The sign of `p * p - (dot(w, w) - r * r)`, the ray-against-ball
/// discriminant, degree two, with `w = to - from` and `p = dot(w, direction)`.
///
/// The projection `p` is a product of two exact components, so it carries
/// `2^2148` and `p * p` carries `2^4296`, while the constant term carries
/// `2^2148`. `Big` holds integers, so the constant term is scaled UP by the
/// POSITIVE `2^2148` rather than the projection being scaled down, and the
/// comparison stays exact. `error.NonFiniteResult` is the `fast` policy
/// saying its own evaluation came out a NaN.
fn rayBallDiscriminantSign(
    comptime dim: usize,
    comptime T: type,
    from: [dim]T,
    to: [dim]T,
    direction: [dim]T,
    r: T,
    policy: PredicatePolicy,
) error{NonFiniteResult}!Sign {
    return switch (policy) {
        .fast => blk: {
            const w = arraySub(dim, T, to, from);
            const p = dotOfArray(dim, T, w, direction);
            break :blk try fastSign(T, subRaw(T, scaleRaw(T, p, p), subRaw(T, dotOfArray(dim, T, w, w), scaleRaw(T, r, r))));
        },
        .adaptive_exact => blk: {
            const from_row = exactRow(dim, T, from);
            const to_row = exactRow(dim, T, to);
            var deltas: [dim]Big = undefined;
            // SAFETY: Big.fromFloat writes metadata and active limbs before multiplication.
            var radius: Big = undefined;
            // SAFETY: exactDotInto writes metadata and active limbs before multiplication.
            var projection: Big = undefined;
            // SAFETY: exactNormInto writes metadata and active limbs before subtraction.
            var constant: Big = undefined;
            // SAFETY: Big.mulInto writes metadata and active limbs before subtraction.
            var square: Big = undefined;
            // SAFETY: Big.shiftedLeftInto writes metadata and active limbs before subtraction.
            var scaled: Big = undefined;
            // SAFETY: Big.mulInto writes metadata and active limbs before subtraction.
            var total: Big = undefined;
            exactDeltasInto(dim, &deltas, &from_row, &to_row);
            Big.fromFloat(&radius, widen(T, r));
            const direction_row = exactRow(dim, T, direction);
            exactDotInto(dim, &projection, &deltas, &direction_row);
            exactNormInto(dim, &constant, &deltas);
            Big.mulInto(&square, &radius, &radius);
            Big.subInto(&constant, &constant, &square);
            Big.shiftedLeftInto(&scaled, &constant, projection_scale);
            Big.mulInto(&total, &projection, &projection);
            Big.subInto(&total, &total, &scaled);
            break :blk Big.sign(&total);
        },
    };
}

/// The sign of `distance * distance - r * r`, the ball-against-flat
/// comparison, degree two, with `distance = dot(to - from, normal)`.
///
/// Same scale split as the ray discriminant: the distance carries `2^2148`
/// and its square carries `2^4296`, so the radius square is scaled up by the
/// POSITIVE `2^2148` before the subtraction. `error.NonFiniteResult` is the
/// `fast` policy saying its own evaluation came out a NaN.
fn flatBallSign(
    comptime dim: usize,
    comptime T: type,
    from: [dim]T,
    to: [dim]T,
    normal: [dim]T,
    r: T,
    policy: PredicatePolicy,
) error{NonFiniteResult}!Sign {
    return switch (policy) {
        .fast => blk: {
            const distance = dotOfArray(dim, T, arraySub(dim, T, to, from), normal);
            break :blk try fastSign(T, subRaw(T, scaleRaw(T, distance, distance), scaleRaw(T, r, r)));
        },
        .adaptive_exact => blk: {
            const from_row = exactRow(dim, T, from);
            const to_row = exactRow(dim, T, to);
            var deltas: [dim]Big = undefined;
            // SAFETY: exactDotInto writes metadata and active limbs before multiplication.
            var projection: Big = undefined;
            // SAFETY: Big.fromFloat writes metadata and active limbs before multiplication.
            var radius: Big = undefined;
            // SAFETY: Big.mulInto writes metadata and active limbs before shifting.
            var square: Big = undefined;
            // SAFETY: Big.shiftedLeftInto writes metadata and active limbs before subtraction.
            var scaled: Big = undefined;
            // SAFETY: Big.mulInto writes metadata and active limbs before subtraction.
            var total: Big = undefined;
            exactDeltasInto(dim, &deltas, &from_row, &to_row);
            const normal_row = exactRow(dim, T, normal);
            exactDotInto(dim, &projection, &deltas, &normal_row);
            Big.fromFloat(&radius, widen(T, r));
            Big.mulInto(&square, &radius, &radius);
            Big.shiftedLeftInto(&scaled, &square, projection_scale);
            Big.mulInto(&total, &projection, &projection);
            Big.subInto(&total, &total, &scaled);
            break :blk Big.sign(&total);
        },
    };
}

// ---------------------------------------------------------------------------
// A two-dimensional affine flat in four dimensions
// ---------------------------------------------------------------------------

/// A two-flat in four dimensions, which is what two non-parallel hyperplanes
/// meet in. The `API.md` family list for four dimensions has no plane type,
/// but the four-dimensional intersection union has to name the result of two
/// hyperplanes meeting, and the name of that shape is a two-flat: the
/// two-dimensional affine set they share.
///
/// A plain value of a point and two directions. It is built by
/// `Hyperplane4d.intersect` and read by the caller; there is no constructor,
/// and `init` does not exist, so the orthonormal invariant below is
/// established by the producer and is the caller's to keep if it rewrites
/// the fields. It is a SPAN, not a cut: `Hyperplane4d.intersect` reports one
/// only when the two hyperplanes are not parallel, and the shared set is the
/// whole two-flat.
pub const Flat2f = struct {
    /// A point on the flat, taken from the first hyperplane of the pair.
    origin: Point4f,
    /// The first unit direction spanning the flat. Both are unit and mutually
    /// orthogonal, so `u` and `v` span the common set exactly once.
    u: Vector4f,
    /// The second unit direction spanning the flat, orthogonal to `u`.
    v: Vector4f,
};

/// A two-flat in four dimensions, in `f64`, under the same field obligation as
/// `Flat2f`.
pub const Flat2d = struct {
    /// A point on the flat, taken from the first hyperplane of the pair.
    origin: Point4d,
    /// The first unit direction spanning the flat. Both are unit and mutually
    /// orthogonal, so `u` and `v` span the common set exactly once.
    u: Vector4d,
    /// The second unit direction spanning the flat, orthogonal to `u`.
    v: Vector4d,
};

fn Flat2Of(comptime T: type) type {
    return if (T == f32) Flat2f else Flat2d;
}

// ---------------------------------------------------------------------------
// Intersection unions
// ---------------------------------------------------------------------------

/// The result of any two-dimensional intersection, in `f32`.
///
/// The tags, and what each carries:
///
/// - `none`: the two operands have no common point. Every pair can produce
///   this one, which is why it is a tag and not the absence of a tag.
/// - `point`: the common set is a single point. A proper crossing, a touching
///   endpoint, or tangency.
/// - `two_points`: a boundary summary for two crossing balls. In two
///   dimensions these are the circumference crossings, not the closed-disc
///   lens. In three and four dimensions they are two selected boundary
///   points, not the full boundary or solid common set.
/// - `segment`: a closed common run, including a ray's entry-to-exit chord
///   through a ball and a proper overlap of collinear segments.
/// - `line`: the common set is a whole two-sided line, which happens exactly
///   when two lines are one line. A ray is not a line: an identical pair of
///   rays is `coincident`.
/// - `box`: the common set is an axis-aligned box that is neither operand.
///   Two boxes that overlap partially, including the touching case where the
///   box is degenerate. A box IS a set of points; the tag exists because the
///   intersection of two boxes is a box and no other name on this list does.
/// - `contained`: the common set is exactly one operand. A segment inside a
///   collinear segment, one circle inside another, one box inside another.
///   `ContainedSide` says which.
/// - `coincident`: the operands are the same set and that set is not an
///   infinite flat. Identical segments, identical rays, identical circles,
///   identical boxes.
///
/// There is no `plane` tag here and no `hyperplane` tag: two dimensions have
/// neither, and a union that carried a tag it could never produce would be
/// a lie about the geometry. There is no `circle` tag either, because a flat
/// that could cut a circle needs a third dimension, and `API.md` puts the
/// circle in three dimensions.
pub const Intersection2f = Intersection2(f32);
/// The result of any two-dimensional intersection, in `f64`, under the same
/// tag contract as `Intersection2f`.
pub const Intersection2d = Intersection2(f64);
/// The result of any three-dimensional intersection, in `f32`.
///
/// The two-dimensional tags, and three more:
///
/// - `plane`: the common set is the whole plane, which happens exactly when
///   two planes are the same plane. A coincident pair of INFINITE flats names
///   its own span rather than saying `coincident`.
/// - `circle`: the common set of a plane and a sphere that are not tangent.
///   The payload is a centre/radius record of the CUT, not the solid ball and
///   not the plane: `center` is the orthogonal foot of the sphere's centre on
///   the plane and `radius` is the radius of the section. A centre and a
///   radius do not encode a direction, so this payload says nothing about
///   WHICH plane it lies in; a caller that needs the curve as a set in space
///   keeps the plane it passed to `Sphere3d.intersectPlane` or
///   `Plane3d.intersectSphere`. The common set of the plane and the CLOSED
///   solid sphere is the disc this circle bounds, not the circle alone.
///
/// There is no `hyperplane` tag: three dimensions have no four-flat, so a
/// hyperplane tag could never be produced.
pub const Intersection3f = Intersection3(f32);
/// The result of any three-dimensional intersection, in `f64`, under the same
/// tag contract as `Intersection3f`.
pub const Intersection3d = Intersection3(f64);
/// The result of any four-dimensional intersection, in `f32`.
///
/// The three-dimensional tags, and three more:
///
/// - `hyperplane`: the common set is a whole three-flat, which happens
///   exactly when two hyperplanes are the same hyperplane.
/// - `hypersphere`: the common set of a hyperplane and a hypersphere that are
///   not tangent. The payload is the centre/radius record of the cut, read
///   exactly as `Intersection3f.circle` reads its own: `center` is the
///   orthogonal foot of the hypersphere's centre on the hyperplane and
///   `radius` is the radius of the section. It carries no direction, so a
///   caller that needs the section as a set in space keeps the hyperplane it
///   passed in.
///
/// The `plane` tag here names a TWO-flat, not a three-dimensional plane: in
/// four dimensions the tag names the dimension of the common affine span, so
/// two hyperplanes that are not parallel meet in a two-flat and say so. The
/// payload is `Flat2f`, a point plus two orthonormal directions spanning that
/// whole two-flat.
pub const Intersection4f = Intersection4(f32);
/// The result of any four-dimensional intersection, in `f64`, under the same
/// tag contract as `Intersection4f`.
pub const Intersection4d = Intersection4(f64);

fn IntersectionOf(comptime dim: usize, comptime T: type) type {
    return switch (dim) {
        2 => Intersection2(T),
        3 => Intersection3(T),
        4 => Intersection4(T),
        else => @compileError("zmath geometry dimensions are 2, 3, and 4"),
    };
}

fn Intersection2(comptime T: type) type {
    return union(enum) {
        /// No common point.
        none,
        /// One common point.
        point: Point2(T),
        /// Two common points.
        two_points: [2]Point2(T),
        /// A proper sub-segment of one of the operands.
        segment: SegmentOf(2, T),
        /// A whole two-sided line, which means the two lines are one line.
        line: LineOf(2, T),
        /// An axis-aligned box that is neither operand.
        box: Aabb2(T),
        /// One operand, whole.
        contained: ContainedSide,
        /// Both operands, the same set.
        coincident,
    };
}

fn Intersection3(comptime T: type) type {
    return union(enum) {
        /// No common point.
        none,
        /// One common point.
        point: Point3(T),
        /// Two common points.
        two_points: [2]Point3(T),
        /// A proper sub-segment of one of the operands.
        segment: SegmentOf(3, T),
        /// A whole two-sided line, which means the two lines are one line.
        line: LineOf(3, T),
        /// A whole plane, which means the two planes are one plane.
        plane: Plane3(T),
        /// The centre and radius of a plane cutting a sphere, read in the
        /// plane. The payload carries no plane, so the caller keeps its own.
        circle: Sphere3(T),
        /// An axis-aligned box that is neither operand.
        box: Aabb3(T),
        /// One operand, whole.
        contained: ContainedSide,
        /// Both operands, the same set.
        coincident,
    };
}

fn Intersection4(comptime T: type) type {
    return union(enum) {
        /// No common point.
        none,
        /// One common point.
        point: Point4(T),
        /// Two common points.
        two_points: [2]Point4(T),
        /// A proper sub-segment of one of the operands.
        segment: SegmentOf(4, T),
        /// A whole two-sided line, which means the two lines are one line.
        line: LineOf(4, T),
        /// A two-flat, which is what two non-parallel hyperplanes meet in.
        plane: Flat2Of(T),
        /// A whole three-flat, which means the two hyperplanes are one.
        hyperplane: Hyperplane4(T),
        /// The centre and radius of a hyperplane cutting a hypersphere, read
        /// in the hyperplane. The payload carries no hyperplane.
        hypersphere: Hypersphere4(T),
        /// An axis-aligned box that is neither operand.
        box: Aabb4(T),
        /// One operand, whole.
        contained: ContainedSide,
        /// Both operands, the same set.
        coincident,
    };
}

// ---------------------------------------------------------------------------
// Circle, sphere, hypersphere
// ---------------------------------------------------------------------------

fn ballBuild(comptime P: type, center: P, radius: componentType(P)) error{ NonFiniteInput, NegativeRadius }!struct { center: P, radius: componentType(P) } {
    try requireFiniteField(center);
    const T = componentType(P);
    try requireFiniteScalar(T, radius);
    if (radius < 0) return error.NegativeRadius;
    return .{ .center = center, .radius = radius };
}

/// A two-dimensional `f32` circle.
///
/// A plain value of a centre and a radius. It names the CLOSED DISC: `center`
/// and every point at or inside `radius`. `init` established that both are
/// finite and that `radius` is nonnegative; the fields are public and
/// writable, and every operation below re-checks the finiteness of both. The
/// nonnegative invariant is the caller's to keep, and it is what a negative
/// radius would break: no operation here reports `error.NegativeRadius`,
/// because `init` is where that is decided.
///
/// A zero radius is legal and names the single point `center`.
pub const Circle2f = Circle2(f32);
/// A two-dimensional `f64` circle, under the same field obligation as
/// `Circle2f`.
pub const Circle2d = Circle2(f64);
/// A three-dimensional `f32` sphere, under the same field obligation as
/// `Circle2f`.
pub const Sphere3f = Sphere3(f32);
/// A three-dimensional `f64` sphere, under the same field obligation as
/// `Circle2f`.
pub const Sphere3d = Sphere3(f64);
/// A four-dimensional `f32` hypersphere, under the same field obligation as
/// `Circle2f`.
pub const Hypersphere4f = Hypersphere4(f32);
/// A four-dimensional `f64` hypersphere, under the same field obligation as
/// `Circle2f`.
pub const Hypersphere4d = Hypersphere4(f64);

fn SphereOf(comptime dim: usize, comptime T: type) type {
    return switch (dim) {
        2 => Circle2(T),
        3 => Sphere3(T),
        4 => Hypersphere4(T),
        else => @compileError("zmath geometry dimensions are 2, 3, and 4"),
    };
}

fn Circle2(comptime T: type) type {
    return struct {
        /// The centre. Public and writable; `init` established that it is
        /// finite.
        center: Point2(T),
        /// The radius. Public and writable; `init` established that it is
        /// finite and nonnegative, and the caller keeps the nonnegative part.
        radius: T,

        const Self = @This();

        /// Builds the closed disc of that centre and radius. The radius must
        /// be finite and nonnegative: `error.NegativeRadius` rather than a sign
        /// flip, and `error.NonFiniteInput` for a non-finite centre or radius.
        /// Nothing is computed, so this reports no other failure.
        pub fn init(center: Point2(T), radius: T) error{ NonFiniteInput, NegativeRadius }!Self {
            const built = try ballBuild(Point2(T), center, radius);
            return .{ .center = built.center, .radius = built.radius };
        }

        /// True when `point` is inside or ON the disc. The boundary is
        /// inclusive, because a point at exactly the radius is on the set the
        /// circle names, and `adaptive_exact` is what makes that exact.
        /// `error.NonFiniteResult` is the `fast` policy saying its own
        /// evaluation came out a NaN; `adaptive_exact` answers for any finite
        /// centre, radius, and point.
        pub fn contains(self: Self, point: Point2(T), policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!bool {
            try requireFiniteField(self);
            try requireFiniteField(point);
            const sign = try ballContainsSign(2, T, self.center.toArray(), self.radius, point.toArray(), policy);
            return sign == .zero or sign == .negative;
        }

        /// The common set of two discs: `.none` when they are apart, `.point`
        /// when they touch, `.two_points` where their circumferences cross,
        /// `.contained` when one lies wholly inside the other, and `.coincident`
        /// when they are the same disc. `two_points` names where the two
        /// boundaries cross; the common set of the two closed discs is the lens
        /// between them, and this union has no tag for a lens.
        /// `error.NonFiniteResult` when a constructed point or radius leaves
        /// the format.
        pub fn intersect(self: Self, other: Self, policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!Intersection2(T) {
            try requireFiniteField(self);
            try requireFiniteField(other);
            return ballIntersect(2, T, self.center.toArray(), self.radius, other.center.toArray(), other.radius, policy);
        }

        /// The set of points common to this disc and a segment: nothing, one
        /// tangent point, the sub-segment inside the disc, or the whole segment
        /// when it lies inside. `error.NonFiniteResult` when a constructed
        /// point leaves the format.
        pub fn intersectSegment(self: Self, segment: SegmentOf(2, T), policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!Intersection2(T) {
            try requireFiniteField(self);
            try requireFiniteField(segment);
            return ballIntersectSegment(2, T, self.center.toArray(), self.radius, segment.a.toArray(), segment.b.toArray(), policy);
        }

        /// The closed common set of this disc and a ray: `none`, a tangent
        /// `point`, or a `segment` from entry to exit. For an inside origin,
        /// the segment starts at the origin. `error.NonFiniteResult` when
        /// a required construction leaves the format.
        pub fn intersectRay(self: Self, ray: RayOf(2, T), policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!Intersection2(T) {
            try requireFiniteField(self);
            try requireFiniteField(ray);
            return ballIntersectRay(2, T, self.center.toArray(), self.radius, ray.origin.toArray(), ray.direction.toArray(), policy);
        }

        /// True when the centre and the radius are finite. It says nothing
        /// about the radius being nonnegative, which is the caller's invariant
        /// to keep rather than one this module re-checks.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(Point2(T), self.center) and isFiniteOf(T, self.radius);
        }
    };
}

fn Sphere3(comptime T: type) type {
    return struct {
        /// The centre. Public and writable; `init` established that it is
        /// finite.
        center: Point3(T),
        /// The radius. Public and writable; `init` established that it is
        /// finite and nonnegative, and the caller keeps the nonnegative part.
        radius: T,

        const Self = @This();

        /// Builds the closed ball of that centre and radius. The radius must
        /// be finite and nonnegative: `error.NegativeRadius` rather than a sign
        /// flip, and `error.NonFiniteInput` for a non-finite centre or radius.
        /// Nothing is computed, so this reports no other failure.
        pub fn init(center: Point3(T), radius: T) error{ NonFiniteInput, NegativeRadius }!Self {
            const built = try ballBuild(Point3(T), center, radius);
            return .{ .center = built.center, .radius = built.radius };
        }

        /// True when `point` is inside or ON the ball. The boundary is
        /// inclusive. `error.NonFiniteResult` is the `fast` policy saying its
        /// own evaluation came out a NaN; `adaptive_exact` answers for any
        /// finite centre, radius, and point.
        pub fn contains(self: Self, point: Point3(T), policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!bool {
            try requireFiniteField(self);
            try requireFiniteField(point);
            const sign = try ballContainsSign(3, T, self.center.toArray(), self.radius, point.toArray(), policy);
            return sign == .zero or sign == .negative;
        }

        /// The common set of two balls: `.none` when they are apart, `.point`
        /// when they touch, `.two_points` where their surfaces cross,
        /// `.contained` when one lies wholly inside the other, and
        /// `.coincident` when they are the same ball. `two_points` names where
        /// the two surfaces cross; the common set of two closed balls is the
        /// lens between them, and this union has no tag for a lens.
        /// `error.NonFiniteResult` when a constructed point or radius leaves
        /// the format.
        pub fn intersect(self: Self, other: Self, policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!Intersection3(T) {
            try requireFiniteField(self);
            try requireFiniteField(other);
            return ballIntersect(3, T, self.center.toArray(), self.radius, other.center.toArray(), other.radius, policy);
        }

        /// The set of points common to this ball and a segment: nothing, one
        /// tangent point, the sub-segment inside the ball, or the whole segment
        /// when it lies inside. `error.NonFiniteResult` when a constructed
        /// point leaves the format.
        pub fn intersectSegment(self: Self, segment: SegmentOf(3, T), policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!Intersection3(T) {
            try requireFiniteField(self);
            try requireFiniteField(segment);
            return ballIntersectSegment(3, T, self.center.toArray(), self.radius, segment.a.toArray(), segment.b.toArray(), policy);
        }

        /// The closed common set of this ball and a ray: `none`, a tangent
        /// `point`, or a `segment` from entry to exit. For an inside origin,
        /// the segment starts at the origin. `error.NonFiniteResult` when
        /// a required construction leaves the format.
        pub fn intersectRay(self: Self, ray: RayOf(3, T), policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!Intersection3(T) {
            try requireFiniteField(self);
            try requireFiniteField(ray);
            return ballIntersectRay(3, T, self.center.toArray(), self.radius, ray.origin.toArray(), ray.direction.toArray(), policy);
        }

        /// The common set of this ball and a plane: nothing, the one tangent
        /// point, or the disc the plane cuts. The `.circle` payload is a centre
        /// and a radius read INSIDE that plane: `center` is the orthogonal foot
        /// of this ball's centre on the plane and `radius` is the radius of the
        /// section. A centre and a radius carry no direction, so the payload
        /// does not say which plane it lies in; a caller that needs the cut as
        /// a set in space keeps the plane it passed here.
        /// `error.NonFiniteResult` when a constructed point or radius leaves
        /// the format.
        pub fn intersectPlane(self: Self, plane: Plane3(T), policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!Intersection3(T) {
            try requireFiniteField(self);
            try requireFiniteField(plane);
            return ballIntersectFlat(3, T, self.center.toArray(), self.radius, plane.origin.toArray(), plane.normal.toArray(), policy);
        }

        /// True when the centre and the radius are finite. It says nothing
        /// about the radius being nonnegative, which is the caller's invariant
        /// to keep rather than one this module re-checks.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(Point3(T), self.center) and isFiniteOf(T, self.radius);
        }
    };
}

fn Hypersphere4(comptime T: type) type {
    return struct {
        /// The centre. Public and writable; `init` established that it is
        /// finite.
        center: Point4(T),
        /// The radius. Public and writable; `init` established that it is
        /// finite and nonnegative, and the caller keeps the nonnegative part.
        radius: T,

        const Self = @This();

        /// Builds the closed ball of that centre and radius. The radius must
        /// be finite and nonnegative: `error.NegativeRadius` rather than a sign
        /// flip, and `error.NonFiniteInput` for a non-finite centre or radius.
        /// Nothing is computed, so this reports no other failure.
        pub fn init(center: Point4(T), radius: T) error{ NonFiniteInput, NegativeRadius }!Self {
            const built = try ballBuild(Point4(T), center, radius);
            return .{ .center = built.center, .radius = built.radius };
        }

        /// True when `point` is inside or ON the ball. The boundary is
        /// inclusive. `error.NonFiniteResult` is the `fast` policy saying its
        /// own evaluation came out a NaN; `adaptive_exact` answers for any
        /// finite centre, radius, and point.
        pub fn contains(self: Self, point: Point4(T), policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!bool {
            try requireFiniteField(self);
            try requireFiniteField(point);
            const sign = try ballContainsSign(4, T, self.center.toArray(), self.radius, point.toArray(), policy);
            return sign == .zero or sign == .negative;
        }

        /// The common set of two balls: `.none` when they are apart, `.point`
        /// when they touch, `.two_points` where their surfaces cross,
        /// `.contained` when one lies wholly inside the other, and
        /// `.coincident` when they are the same ball. `two_points` names where
        /// the two surfaces cross; the common set of two closed balls is the
        /// lens between them, and this union has no tag for a lens.
        /// `error.NonFiniteResult` when a constructed point or radius leaves
        /// the format.
        pub fn intersect(self: Self, other: Self, policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!Intersection4(T) {
            try requireFiniteField(self);
            try requireFiniteField(other);
            return ballIntersect(4, T, self.center.toArray(), self.radius, other.center.toArray(), other.radius, policy);
        }

        /// The set of points common to this ball and a segment: nothing, one
        /// tangent point, the sub-segment inside the ball, or the whole segment
        /// when it lies inside. `error.NonFiniteResult` when a constructed
        /// point leaves the format.
        pub fn intersectSegment(self: Self, segment: SegmentOf(4, T), policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!Intersection4(T) {
            try requireFiniteField(self);
            try requireFiniteField(segment);
            return ballIntersectSegment(4, T, self.center.toArray(), self.radius, segment.a.toArray(), segment.b.toArray(), policy);
        }

        /// The closed common set of this ball and a ray: `none`, a tangent
        /// `point`, or a `segment` from entry to exit. For an inside origin,
        /// the segment starts at the origin. `error.NonFiniteResult` when
        /// a required construction leaves the format.
        pub fn intersectRay(self: Self, ray: RayOf(4, T), policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!Intersection4(T) {
            try requireFiniteField(self);
            try requireFiniteField(ray);
            return ballIntersectRay(4, T, self.center.toArray(), self.radius, ray.origin.toArray(), ray.direction.toArray(), policy);
        }

        /// The common set of this ball and a hyperplane: nothing, the one
        /// tangent point, or the three-ball the hyperplane cuts. The
        /// `.hypersphere` payload is a centre and a radius read INSIDE that
        /// hyperplane: `center` is the orthogonal foot of this ball's centre on
        /// the hyperplane and `radius` is the radius of the section. A centre
        /// and a radius carry no direction, so the payload does not say which
        /// hyperplane it lies in; a caller that needs the section as a set in
        /// space keeps the hyperplane it passed here.
        /// `error.NonFiniteResult` when a constructed point or radius leaves
        /// the format.
        pub fn intersectHyperplane(self: Self, plane: Hyperplane4(T), policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!Intersection4(T) {
            try requireFiniteField(self);
            try requireFiniteField(plane);
            return ballIntersectFlat(4, T, self.center.toArray(), self.radius, plane.origin.toArray(), plane.normal.toArray(), policy);
        }

        /// True when the centre and the radius are finite. It says nothing
        /// about the radius being nonnegative, which is the caller's invariant
        /// to keep rather than one this module re-checks.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(Point4(T), self.center) and isFiniteOf(T, self.radius);
        }
    };
}

// ---------------------------------------------------------------------------
// Plane, hyperplane
// ---------------------------------------------------------------------------

/// A three-dimensional `f32` plane.
///
/// A plain value of a point and a normal. `init` established that `origin` is
/// finite and that `normal` is a UNIT vector, and that unit invariant is the
/// caller's to keep: the fields are public, and every operation below gates
/// the finiteness of both without re-checking the length. The sign of the
/// normal decides which of the two sides is positive, so a caller that
/// rewrites it changes which side is which.
pub const Plane3f = Plane3(f32);
/// A three-dimensional `f64` plane, under the same field obligation as
/// `Plane3f`.
pub const Plane3d = Plane3(f64);
/// A four-dimensional `f32` hyperplane, under the same field obligation as
/// `Plane3f`.
pub const Hyperplane4f = Hyperplane4(f32);
/// A four-dimensional `f64` hyperplane, under the same field obligation as
/// `Plane3f`.
pub const Hyperplane4d = Hyperplane4(f64);

fn PlaneOf(comptime dim: usize, comptime T: type) type {
    return switch (dim) {
        3 => Plane3(T),
        4 => Hyperplane4(T),
        else => @compileError("zmath flats are three- and four-dimensional"),
    };
}

/// The shared constructor of a plane and a hyperplane: a point on the flat and
/// a normalized normal. It normalizes, so it inherits `vectorNormalize`'s
/// arithmetic and therefore its `error.NonFiniteResult`.
fn flatBuild(comptime P: type, point: P, normal: VectorOf(fieldCount(P), componentType(P))) error{ NonFiniteInput, NonFiniteResult, Degenerate }!struct { origin: P, normal: VectorOf(fieldCount(P), componentType(P)) } {
    const built = try lineBuild(P, point, normal);
    return .{ .origin = built.origin, .normal = built.direction };
}

fn Plane3(comptime T: type) type {
    return struct {
        /// A point on the plane. Public and writable; `init` established that
        /// it is finite.
        origin: Point3(T),
        /// The unit normal. Public and writable; the caller keeps it unit.
        normal: Vector3(T),

        const Self = @This();

        /// Builds the plane through `point` with `normal`, which is normalized.
        /// `error.Degenerate` for a zero normal, `error.NonFiniteInput` for a
        /// non-finite point or normal, and `error.NonFiniteResult` for a normal
        /// whose norm is past the top of the range.
        pub fn init(point: Point3(T), normal: Vector3(T)) error{ NonFiniteInput, NonFiniteResult, Degenerate }!Self {
            const built = try flatBuild(Point3(T), point, normal);
            return .{ .origin = built.origin, .normal = built.normal };
        }

        /// The signed distance from `point` to the plane, that is
        /// `dot(point - origin, normal)`. It is zero on the plane, positive on
        /// the side the normal points toward, and negative on the other side,
        /// and it equals the ordinary distance because `normal` is unit.
        /// `error.NonFiniteResult` when a term or the sum leaves the format,
        /// which a `point` at the opposite end of the range from `origin`
        /// reaches. This is the plain scalar evaluation, not a predicate: it
        /// takes no policy and returns no sign.
        pub fn signedDistance(self: Self, point: Point3(T)) error{ NonFiniteInput, NonFiniteResult }!T {
            try requireFiniteField(self);
            try requireFiniteField(point);
            return requireFiniteResult(T, flatDistanceFast(3, T, point.toArray(), self.origin.toArray(), self.normal.toArray()));
        }

        /// The set of points common to two planes: nothing, the whole plane
        /// when the two are one plane (`.plane`, because naming the span is
        /// more useful than saying `coincident`), or the line they meet in
        /// (`.line`). `error.NonFiniteResult` when a constructed point,
        /// direction, or parameter leaves the format.
        pub fn intersect(self: Self, other: Self, policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!Intersection3(T) {
            try requireFiniteField(self);
            try requireFiniteField(other);
            return flatPair3(T, self.origin.toArray(), self.normal.toArray(), other.origin.toArray(), other.normal.toArray(), policy);
        }

        /// The set of points common to this plane and a segment: nothing, the
        /// one point the plane cuts the segment at, or the whole segment when
        /// it lies in the plane (`.contained = .second`). `error.NonFiniteResult`
        /// when a constructed point or parameter leaves the format.
        pub fn intersectSegment(self: Self, segment: SegmentOf(3, T), policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!Intersection3(T) {
            try requireFiniteField(self);
            try requireFiniteField(segment);
            return flatIntersectSegment(3, T, self.origin.toArray(), self.normal.toArray(), segment.a.toArray(), segment.b.toArray(), policy);
        }

        /// The common set of this plane and a ball: nothing, the one tangent
        /// point, or the disc the plane cuts. The `.circle` payload is the
        /// centre and radius of that cut, read inside THIS plane, and it
        /// carries no direction, so the plane a caller needs to place it in
        /// space is the one it passed here. `error.NonFiniteResult` when a
        /// constructed point or radius leaves the format.
        pub fn intersectSphere(self: Self, sphere: Sphere3(T), policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!Intersection3(T) {
            try requireFiniteField(self);
            try requireFiniteField(sphere);
            return ballIntersectFlat(3, T, sphere.center.toArray(), sphere.radius, self.origin.toArray(), self.normal.toArray(), policy);
        }

        /// True when the origin and the normal are finite. It says nothing
        /// about the normal being unit, which is the caller's invariant to
        /// keep rather than one this module re-checks.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(Point3(T), self.origin) and isFiniteQuiet(Vector3(T), self.normal);
        }
    };
}

fn Hyperplane4(comptime T: type) type {
    return struct {
        /// A point on the hyperplane. Public and writable; `init` established
        /// that it is finite.
        origin: Point4(T),
        /// The unit normal. Public and writable; the caller keeps it unit.
        normal: Vector4(T),

        const Self = @This();

        /// Builds the hyperplane through `point` with `normal`, which is
        /// normalized. `error.Degenerate` for a zero normal,
        /// `error.NonFiniteInput` for a non-finite point or normal, and
        /// `error.NonFiniteResult` for a normal whose norm is past the top of
        /// the range.
        pub fn init(point: Point4(T), normal: Vector4(T)) error{ NonFiniteInput, NonFiniteResult, Degenerate }!Self {
            const built = try flatBuild(Point4(T), point, normal);
            return .{ .origin = built.origin, .normal = built.normal };
        }

        /// The signed distance from `point` to the hyperplane, read exactly as
        /// `Plane3d.signedDistance` reads its own. It takes no policy and
        /// returns no sign, and `error.NonFiniteResult` when a term or the sum
        /// leaves the format.
        pub fn signedDistance(self: Self, point: Point4(T)) error{ NonFiniteInput, NonFiniteResult }!T {
            try requireFiniteField(self);
            try requireFiniteField(point);
            return requireFiniteResult(T, flatDistanceFast(4, T, point.toArray(), self.origin.toArray(), self.normal.toArray()));
        }

        /// The set of points common to two hyperplanes: nothing, the whole
        /// three-flat when the two are one hyperplane (`.hyperplane`), or the
        /// two-flat they meet in (`.plane`, which in four dimensions names the
        /// dimension of the common span). `error.NonFiniteResult` when a
        /// constructed point, direction, or parameter leaves the format.
        pub fn intersect(self: Self, other: Self, policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!Intersection4(T) {
            try requireFiniteField(self);
            try requireFiniteField(other);
            return flatPair4(T, self.origin.toArray(), self.normal.toArray(), other.origin.toArray(), other.normal.toArray(), policy);
        }

        /// The set of points common to this hyperplane and a segment: nothing,
        /// the one point the hyperplane cuts the segment at, or the whole
        /// segment when it lies in the hyperplane (`.contained = .second`).
        /// `error.NonFiniteResult` when a constructed point or parameter
        /// leaves the format.
        pub fn intersectSegment(self: Self, segment: SegmentOf(4, T), policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!Intersection4(T) {
            try requireFiniteField(self);
            try requireFiniteField(segment);
            return flatIntersectSegment(4, T, self.origin.toArray(), self.normal.toArray(), segment.a.toArray(), segment.b.toArray(), policy);
        }

        /// The common set of this hyperplane and a ball: nothing, the one
        /// tangent point, or the three-ball the hyperplane cuts. The
        /// `.hypersphere` payload is the centre and radius of that cut, read
        /// inside THIS hyperplane, and it carries no direction, so the
        /// hyperplane a caller needs to place it in space is the one it passed
        /// here. `error.NonFiniteResult` when a constructed point or radius
        /// leaves the format.
        pub fn intersectSphere(self: Self, sphere: Hypersphere4(T), policy: PredicatePolicy) error{ NonFiniteInput, NonFiniteResult }!Intersection4(T) {
            try requireFiniteField(self);
            try requireFiniteField(sphere);
            return ballIntersectFlat(4, T, sphere.center.toArray(), sphere.radius, self.origin.toArray(), self.normal.toArray(), policy);
        }

        /// True when the origin and the normal are finite. It says nothing
        /// about the normal being unit, which is the caller's invariant to
        /// keep rather than one this module re-checks.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(Point4(T), self.origin) and isFiniteQuiet(Vector4(T), self.normal);
        }
    };
}

// ---------------------------------------------------------------------------
// Axis-aligned bounding box
// ---------------------------------------------------------------------------

fn aabbBuild(comptime P: type, minimum: P, maximum: P) error{ NonFiniteInput, InvertedBounds }!P {
    try requireFiniteField(minimum);
    try requireFiniteField(maximum);
    const lows = toArrayOf(P, minimum);
    const highs = toArrayOf(P, maximum);
    inline for (0..comptime fieldCount(P)) |index| {
        if (lows[index] > highs[index]) return error.InvertedBounds;
    }
    return minimum;
}

/// The finiteness gate for a box CONSUMED by an operation, as opposed to the
/// one a constructor ran. It has to cover BOTH corners: the fields are
/// public, so a caller can poison `maximum` after `init` accepted it, and
/// every operation below reads `maximum` in its arithmetic. Validating
/// `minimum` alone let an infinite `maximum` through and produced a finite
/// answer from it, which is exactly the propagation the module comment
/// forbids.
fn requireFiniteBox(comptime P: type, minimum: P, maximum: P) error{NonFiniteInput}!void {
    try requireFiniteField(minimum);
    try requireFiniteField(maximum);
}

/// A two-dimensional `f32` axis-aligned box.
///
/// A plain value of two corners. It names the CLOSED region between them, and
/// `init` established two invariants: both corners are finite, and no
/// `minimum` component is above the matching `maximum` component. The fields
/// are public and writable; every operation below re-checks both corners for
/// finiteness, but none of them re-checks the ordering, so that invariant is
/// the caller's to keep. A degenerate box, `minimum == maximum`, is legal and
/// names a single point.
pub const Aabb2f = Aabb2(f32);
/// A two-dimensional `f64` box, under the same field obligation as `Aabb2f`.
pub const Aabb2d = Aabb2(f64);
/// A three-dimensional `f32` box, under the same field obligation as `Aabb2f`.
pub const Aabb3f = Aabb3(f32);
/// A three-dimensional `f64` box, under the same field obligation as `Aabb2f`.
pub const Aabb3d = Aabb3(f64);
/// A four-dimensional `f32` box, under the same field obligation as `Aabb2f`.
pub const Aabb4f = Aabb4(f32);
/// A four-dimensional `f64` box, under the same field obligation as `Aabb2f`.
pub const Aabb4d = Aabb4(f64);

fn AabbOf(comptime dim: usize, comptime T: type) type {
    return switch (dim) {
        2 => Aabb2(T),
        3 => Aabb3(T),
        4 => Aabb4(T),
        else => @compileError("zmath geometry dimensions are 2, 3, and 4"),
    };
}

fn aabbContainsOf(comptime P: type, minimum: P, maximum: P, point: P) bool {
    const lows = toArrayOf(P, minimum);
    const highs = toArrayOf(P, maximum);
    const values = toArrayOf(P, point);
    var inside = true;
    inline for (0..comptime fieldCount(P)) |index| {
        inside = inside and values[index] >= lows[index] and values[index] <= highs[index];
    }
    return inside;
}

fn aabbOverlapsOf(comptime P: type, a_minimum: P, a_maximum: P, b_minimum: P, b_maximum: P) bool {
    const a_low = toArrayOf(P, a_minimum);
    const a_high = toArrayOf(P, a_maximum);
    const b_low = toArrayOf(P, b_minimum);
    const b_high = toArrayOf(P, b_maximum);
    var overlaps = true;
    inline for (0..comptime fieldCount(P)) |index| {
        overlaps = overlaps and a_low[index] <= b_high[index] and b_low[index] <= a_high[index];
    }
    return overlaps;
}

/// The product of the extents: the area in two dimensions, the volume in
/// three, and the hypervolume in four. One expression covers all three because
/// the measure a box "defines" is the product of its side lengths in every
/// dimension. Both a difference and a product can leave the format, so the
/// result is gated: a caller asking for the area of a box that spans the whole
/// range gets `error.NonFiniteResult`, not an infinity.
fn aabbMeasureOf(comptime P: type, minimum: P, maximum: P) error{NonFiniteResult}!componentType(P) {
    const T = componentType(P);
    const lows = toArrayOf(P, minimum);
    const highs = toArrayOf(P, maximum);
    var total = multRaw(T, highs[0] - lows[0], 1);
    inline for (1..comptime fieldCount(P)) |index| {
        total = multRaw(T, total, highs[index] - lows[index]);
    }
    return requireFiniteResult(T, total);
}

/// The full extent along each axis, `maximum - minimum`. A difference of two
/// finite values can leave the format, so it is gated.
fn aabbMeasuresOf(comptime P: type, minimum: P, maximum: P) error{NonFiniteResult}![fieldCount(P)]componentType(P) {
    const lows = toArrayOf(P, minimum);
    const highs = toArrayOf(P, maximum);
    var out: [fieldCount(P)]componentType(P) = undefined;
    inline for (0..comptime fieldCount(P)) |index| {
        out[index] = highs[index] - lows[index];
    }
    return requireFiniteRow(fieldCount(P), componentType(P), out);
}

/// The smallest box that holds `minimum`, `maximum`, and `point` on both
/// corners. `fmin` and `fmax` are the accepted C23 kernels, so a point that
/// is already inside leaves both corners exactly as they were.
fn aabbExpandedOf(comptime P: type, minimum: P, maximum: P, point: P) error{NonFiniteInput}!struct { minimum: P, maximum: P } {
    try requireFiniteField(point);
    const lows = toArrayOf(P, minimum);
    const highs = toArrayOf(P, maximum);
    const values = toArrayOf(P, point);
    var new_low: [fieldCount(P)]componentType(P) = undefined;
    var new_high: [fieldCount(P)]componentType(P) = undefined;
    inline for (0..comptime fieldCount(P)) |index| {
        new_low[index] = minOf(componentType(P), lows[index], values[index]);
        new_high[index] = maxOf(componentType(P), highs[index], values[index]);
    }
    return .{ .minimum = fromArrayOf(P, new_low), .maximum = fromArrayOf(P, new_high) };
}

/// The midpoint, `(minimum + maximum) / 2` per axis. Round the sum first.
/// Multiplication by 0.5 has the same IEEE result as division by 2 and avoids
/// compiler-rt division that flushes MIPS subnormal results.
/// Same-sign large corners can overflow the sum; the result gate rejects it.
fn aabbCenterOf(comptime P: type, minimum: P, maximum: P) error{NonFiniteResult}!P {
    const T = componentType(P);
    const lows = toArrayOf(P, minimum);
    const highs = toArrayOf(P, maximum);
    var out: [fieldCount(P)]componentType(P) = undefined;
    inline for (0..comptime fieldCount(P)) |index| {
        out[index] = scaleRaw(T, addRaw(T, lows[index], highs[index]), 0.5);
    }
    return fromArrayOf(P, try requireFiniteRow(fieldCount(P), T, out));
}

fn Aabb2(comptime T: type) type {
    const P = Point2(T);
    return struct {
        /// The lower corner. Public and writable; `init` established that it
        /// is finite and no component of it is above the matching component of
        /// `maximum`.
        minimum: P,
        /// The upper corner. Public and writable; the ordering against
        /// `minimum` is the caller's to keep.
        maximum: P,

        const Self = @This();

        /// Builds the box. Every `minimum` component must be at or below the
        /// matching `maximum` component: `error.InvertedBounds` rather than a
        /// silently swapped pair, and `error.NonFiniteInput` for a corner that
        /// is not finite. Nothing is computed, so this reports no other
        /// failure.
        pub fn init(minimum: P, maximum: P) error{ NonFiniteInput, InvertedBounds }!Self {
            _ = try aabbBuild(P, minimum, maximum);
            return .{ .minimum = minimum, .maximum = maximum };
        }

        /// True when `point` is inside or ON the box. The boundary is
        /// INCLUSIVE: a component exactly on a face is inside, because the
        /// face is part of the box and the range comparison rounds nothing.
        /// This is why `contains` takes no policy: there is no rounding here
        /// for a policy to change. Nothing is computed, so a non-finite
        /// receiver or point is the only failure.
        pub fn contains(self: Self, point: P) error{NonFiniteInput}!bool {
            try requireFiniteBox(P, self.minimum, self.maximum);
            try requireFiniteField(point);
            return aabbContainsOf(P, self.minimum, self.maximum, point);
        }

        /// The smallest box that holds both this one and `point`, on both
        /// corners. `fmin` and `fmax` are the accepted kernels, so a point
        /// already inside leaves both corners exactly as they were. Nothing is
        /// computed beyond those two selections, so this reports input failure
        /// only.
        pub fn expand(self: Self, point: P) error{NonFiniteInput}!Self {
            try requireFiniteBox(P, self.minimum, self.maximum);
            const grown = try aabbExpandedOf(P, self.minimum, self.maximum, point);
            return .{ .minimum = grown.minimum, .maximum = grown.maximum };
        }

        /// True when the two boxes have a common point. Touching counts:
        /// two boxes that share only a face do have that face in common. The
        /// comparison rounds nothing, so this takes no policy and computes
        /// nothing: a non-finite corner is the only failure.
        pub fn intersects(self: Self, other: Self) error{NonFiniteInput}!bool {
            try requireFiniteBox(P, self.minimum, self.maximum);
            try requireFiniteBox(P, other.minimum, other.maximum);
            return aabbOverlapsOf(P, self.minimum, self.maximum, other.minimum, other.maximum);
        }

        /// The midpoint of the box, read as `(minimum + maximum) / 2` on each
        /// axis. `error.NonFiniteResult` when that sum leaves the format, which
        /// a box spanning the whole range reaches.
        pub fn center(self: Self) error{ NonFiniteInput, NonFiniteResult }!P {
            try requireFiniteBox(P, self.minimum, self.maximum);
            return aabbCenterOf(P, self.minimum, self.maximum);
        }

        /// The size of the box along each axis, that is `maximum - minimum`.
        /// This is the full extent, not the half extent from the centre.
        /// `error.NonFiniteResult` when a difference leaves the format.
        pub fn extent(self: Self) error{ NonFiniteInput, NonFiniteResult }![2]T {
            try requireFiniteBox(P, self.minimum, self.maximum);
            return aabbMeasuresOf(P, self.minimum, self.maximum);
        }

        /// The area of the box, the product of its two extents.
        /// `error.NonFiniteResult` when an extent or the product leaves the
        /// format.
        pub fn area(self: Self) error{ NonFiniteInput, NonFiniteResult }!T {
            try requireFiniteBox(P, self.minimum, self.maximum);
            return aabbMeasureOf(P, self.minimum, self.maximum);
        }

        /// The set of points common to two boxes: `.none` when they are apart,
        /// the `.box` that is neither operand when they overlap partially, the
        /// `.contained` operand when one holds the other whole, and
        /// `.coincident` when the two are the same box. The answer is always
        /// the WHOLE intersection, so a box that meets its argument only on a
        /// face or a corner comes back as a degenerate `.box`, not as a point.
        /// `fmin` and `fmax` only select among corners the caller supplied, so
        /// nothing is computed here and a non-finite corner is the only
        /// failure.
        pub fn intersect(self: Self, other: Self) error{NonFiniteInput}!Intersection2(T) {
            try requireFiniteBox(P, self.minimum, self.maximum);
            try requireFiniteBox(P, other.minimum, other.maximum);
            return aabbIntersect(2, T, self.minimum, self.maximum, other.minimum, other.maximum);
        }

        /// True when both corners are finite. It says nothing about the
        /// ordering, which is the caller's invariant to keep rather than one
        /// this module re-checks. This is the check a caller needs after
        /// writing the fields directly.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(P, self.minimum) and isFiniteQuiet(P, self.maximum);
        }
    };
}

fn multRaw(comptime T: type, a: T, b: T) T {
    @setFloatMode(.strict);
    return a * b;
}

fn Aabb3(comptime T: type) type {
    const P = Point3(T);
    return struct {
        /// The lower corner. Public and writable; `init` established that it
        /// is finite and no component of it is above the matching component of
        /// `maximum`.
        minimum: P,
        /// The upper corner. Public and writable; the ordering against
        /// `minimum` is the caller's to keep.
        maximum: P,

        const Self = @This();

        /// Builds the box, rejecting an inverted one. `error.InvertedBounds`
        /// when a `minimum` component is above the matching `maximum`
        /// component, and `error.NonFiniteInput` for a corner that is not
        /// finite. Nothing is computed, so this reports no other failure.
        pub fn init(minimum: P, maximum: P) error{ NonFiniteInput, InvertedBounds }!Self {
            _ = try aabbBuild(P, minimum, maximum);
            return .{ .minimum = minimum, .maximum = maximum };
        }

        /// True when `point` is inside or ON the box. The boundary is
        /// INCLUSIVE: a component exactly on a face is inside, because the
        /// face is part of the box and the range comparison rounds nothing.
        /// This is why `contains` takes no policy. Nothing is computed, so a
        /// non-finite receiver or point is the only failure.
        pub fn contains(self: Self, point: P) error{NonFiniteInput}!bool {
            try requireFiniteBox(P, self.minimum, self.maximum);
            try requireFiniteField(point);
            return aabbContainsOf(P, self.minimum, self.maximum, point);
        }

        /// The smallest box that holds both this one and `point`, on both
        /// corners. `fmin` and `fmax` are the accepted kernels, so a point
        /// already inside leaves both corners exactly as they were. Nothing is
        /// computed beyond those two selections, so this reports input failure
        /// only.
        pub fn expand(self: Self, point: P) error{NonFiniteInput}!Self {
            try requireFiniteBox(P, self.minimum, self.maximum);
            const grown = try aabbExpandedOf(P, self.minimum, self.maximum, point);
            return .{ .minimum = grown.minimum, .maximum = grown.maximum };
        }

        /// True when the two boxes have a common point. Touching counts, and
        /// the comparison rounds nothing, so this takes no policy and computes
        /// nothing: a non-finite corner is the only failure.
        pub fn intersects(self: Self, other: Self) error{NonFiniteInput}!bool {
            try requireFiniteBox(P, self.minimum, self.maximum);
            try requireFiniteBox(P, other.minimum, other.maximum);
            return aabbOverlapsOf(P, self.minimum, self.maximum, other.minimum, other.maximum);
        }

        /// The midpoint of the box, read as `(minimum + maximum) / 2` on each
        /// axis. `error.NonFiniteResult` when that sum leaves the format.
        pub fn center(self: Self) error{ NonFiniteInput, NonFiniteResult }!P {
            try requireFiniteBox(P, self.minimum, self.maximum);
            return aabbCenterOf(P, self.minimum, self.maximum);
        }

        /// The size of the box along each axis, that is `maximum - minimum`,
        /// the full extent and not the half extent from the centre.
        /// `error.NonFiniteResult` when a difference leaves the format.
        pub fn extent(self: Self) error{ NonFiniteInput, NonFiniteResult }![3]T {
            try requireFiniteBox(P, self.minimum, self.maximum);
            return aabbMeasuresOf(P, self.minimum, self.maximum);
        }

        /// The volume of the box, the product of its three extents.
        /// `error.NonFiniteResult` when an extent or the product leaves the
        /// format.
        pub fn volume(self: Self) error{ NonFiniteInput, NonFiniteResult }!T {
            try requireFiniteBox(P, self.minimum, self.maximum);
            return aabbMeasureOf(P, self.minimum, self.maximum);
        }

        /// Twice the surface area of the box, `2 * (xy + xz + yz)`. The factor
        /// of two keeps the expression a product of the extents, so it is exact
        /// wherever the extents are; a caller who wants the area divides.
        /// `error.NonFiniteResult` when an extent or a product leaves the
        /// format.
        pub fn doubledSurfaceArea(self: Self) error{ NonFiniteInput, NonFiniteResult }!T {
            try requireFiniteBox(P, self.minimum, self.maximum);
            const sides = try aabbMeasuresOf(P, self.minimum, self.maximum);
            return requireFiniteResult(T, scaleRaw(T, addRaw(T, addRaw(T, multRaw(T, sides[0], sides[1]), multRaw(T, sides[0], sides[2])), multRaw(T, sides[1], sides[2])), 2));
        }

        /// The set of points common to two boxes: `.none` when they are apart,
        /// the `.box` that is neither operand when they overlap partially, the
        /// `.contained` operand when one holds the other whole, and
        /// `.coincident` when the two are the same box. The answer is always
        /// the WHOLE intersection, so a box that meets its argument only on a
        /// face, an edge, or a corner comes back as a degenerate `.box`, not as
        /// a lower-dimensional shape. Nothing is computed here beyond
        /// `fmin` and `fmax` over the caller's own corners, so a non-finite
        /// corner is the only failure.
        pub fn intersect(self: Self, other: Self) error{NonFiniteInput}!Intersection3(T) {
            try requireFiniteBox(P, self.minimum, self.maximum);
            try requireFiniteBox(P, other.minimum, other.maximum);
            return aabbIntersect(3, T, self.minimum, self.maximum, other.minimum, other.maximum);
        }

        /// True when both corners are finite. It says nothing about the
        /// ordering, which is the caller's invariant to keep rather than one
        /// this module re-checks. This is the check a caller needs after
        /// writing the fields directly.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(P, self.minimum) and isFiniteQuiet(P, self.maximum);
        }
    };
}

fn Aabb4(comptime T: type) type {
    const P = Point4(T);
    return struct {
        /// The lower corner. Public and writable; `init` established that it
        /// is finite and no component of it is above the matching component of
        /// `maximum`.
        minimum: P,
        /// The upper corner. Public and writable; the ordering against
        /// `minimum` is the caller's to keep.
        maximum: P,

        const Self = @This();

        /// Builds the box, rejecting an inverted one. `error.InvertedBounds`
        /// when a `minimum` component is above the matching `maximum`
        /// component, and `error.NonFiniteInput` for a corner that is not
        /// finite. Nothing is computed, so this reports no other failure.
        pub fn init(minimum: P, maximum: P) error{ NonFiniteInput, InvertedBounds }!Self {
            _ = try aabbBuild(P, minimum, maximum);
            return .{ .minimum = minimum, .maximum = maximum };
        }

        /// True when `point` is inside or ON the box. The boundary is
        /// INCLUSIVE: a component exactly on a face is inside, because the
        /// face is part of the box and the range comparison rounds nothing.
        /// This is why `contains` takes no policy. Nothing is computed, so a
        /// non-finite receiver or point is the only failure.
        pub fn contains(self: Self, point: P) error{NonFiniteInput}!bool {
            try requireFiniteBox(P, self.minimum, self.maximum);
            try requireFiniteField(point);
            return aabbContainsOf(P, self.minimum, self.maximum, point);
        }

        /// The smallest box that holds both this one and `point`, on both
        /// corners. `fmin` and `fmax` are the accepted kernels, so a point
        /// already inside leaves both corners exactly as they were. Nothing is
        /// computed beyond those two selections, so this reports input failure
        /// only.
        pub fn expand(self: Self, point: P) error{NonFiniteInput}!Self {
            try requireFiniteBox(P, self.minimum, self.maximum);
            const grown = try aabbExpandedOf(P, self.minimum, self.maximum, point);
            return .{ .minimum = grown.minimum, .maximum = grown.maximum };
        }

        /// True when the two boxes have a common point. Touching counts, and
        /// the comparison rounds nothing, so this takes no policy and computes
        /// nothing: a non-finite corner is the only failure.
        pub fn intersects(self: Self, other: Self) error{NonFiniteInput}!bool {
            try requireFiniteBox(P, self.minimum, self.maximum);
            try requireFiniteBox(P, other.minimum, other.maximum);
            return aabbOverlapsOf(P, self.minimum, self.maximum, other.minimum, other.maximum);
        }

        /// The midpoint of the box, read as `(minimum + maximum) / 2` on each
        /// axis. `error.NonFiniteResult` when that sum leaves the format.
        pub fn center(self: Self) error{ NonFiniteInput, NonFiniteResult }!P {
            try requireFiniteBox(P, self.minimum, self.maximum);
            return aabbCenterOf(P, self.minimum, self.maximum);
        }

        /// The size of the box along each axis, that is `maximum - minimum`,
        /// the full extent and not the half extent from the centre.
        /// `error.NonFiniteResult` when a difference leaves the format.
        pub fn extent(self: Self) error{ NonFiniteInput, NonFiniteResult }![4]T {
            try requireFiniteBox(P, self.minimum, self.maximum);
            return aabbMeasuresOf(P, self.minimum, self.maximum);
        }

        /// The four-dimensional hypervolume of the box, the product of its four
        /// extents. `error.NonFiniteResult` when an extent or the product
        /// leaves the format.
        pub fn hyperVolume(self: Self) error{ NonFiniteInput, NonFiniteResult }!T {
            try requireFiniteBox(P, self.minimum, self.maximum);
            return aabbMeasureOf(P, self.minimum, self.maximum);
        }

        /// The set of points common to two boxes: `.none` when they are apart,
        /// the `.box` that is neither operand when they overlap partially, the
        /// `.contained` operand when one holds the other whole, and
        /// `.coincident` when the two are the same box. The answer is always
        /// the WHOLE intersection, and the per-axis test decides it: a pair
        /// that agrees on three axes and not on the fourth has no common point
        /// at all, because a box that disagreed on an axis would have a lower
        /// corner above its upper one. Nothing is computed beyond `fmin` and
        /// `fmax` over the caller's own corners, so a non-finite corner is the
        /// only failure.
        pub fn intersect(self: Self, other: Self) error{NonFiniteInput}!Intersection4(T) {
            try requireFiniteBox(P, self.minimum, self.maximum);
            try requireFiniteBox(P, other.minimum, other.maximum);
            return aabbIntersect(4, T, self.minimum, self.maximum, other.minimum, other.maximum);
        }

        /// True when both corners are finite. It says nothing about the
        /// ordering, which is the caller's invariant to keep rather than one
        /// this module re-checks. This is the check a caller needs after
        /// writing the fields directly.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(P, self.minimum) and isFiniteQuiet(P, self.maximum);
        }
    };
}

/// The whole-set classification of two boxes.
///
/// Per axis the common interval is `[max(lows), min(highs)]`, and `fmin` and
/// `fmax` only select among corners the caller supplied, so the result is
/// always made of their own bits and this computes nothing that can overflow.
/// An axis whose common interval is empty makes the WHOLE intersection empty,
/// whatever the other axes say. When every axis is nonempty the common box is
/// compared with each operand: equal to one of them the answer is that
/// operand whole, equal to both the two are the same box, and equal to
/// neither the answer is a box that is neither operand.
fn aabbIntersect(comptime dim: usize, comptime T: type, a_minimum: PointOf(dim, T), a_maximum: PointOf(dim, T), b_minimum: PointOf(dim, T), b_maximum: PointOf(dim, T)) IntersectionOf(dim, T) {
    const P = PointOf(dim, T);
    const a_low = toArrayOf(P, a_minimum);
    const a_high = toArrayOf(P, a_maximum);
    const b_low = toArrayOf(P, b_minimum);
    const b_high = toArrayOf(P, b_maximum);
    var low: [dim]T = undefined;
    var high: [dim]T = undefined;
    var disjoint = false;
    inline for (0..dim) |index| {
        low[index] = maxOf(T, a_low[index], b_low[index]);
        high[index] = minOf(T, a_high[index], b_high[index]);
        disjoint = disjoint or low[index] > high[index];
    }
    if (disjoint) return .none;
    var equal_a = true;
    var equal_b = true;
    inline for (0..dim) |index| {
        equal_a = equal_a and low[index] == a_low[index] and high[index] == a_high[index];
        equal_b = equal_b and low[index] == b_low[index] and high[index] == b_high[index];
    }
    if (equal_a and equal_b) return .coincident;
    if (equal_a) return .{ .contained = .first };
    if (equal_b) return .{ .contained = .second };
    return .{ .box = .{ .minimum = fromArrayOf(P, low), .maximum = fromArrayOf(P, high) } };
}

// ---------------------------------------------------------------------------
// The affine intersection engine
//
// Every pair of lines, rays, and segments reduces to the same question, in
// any dimension: do two coplanar objects meet, and where. In two dimensions
// the question is asked directly. In three and four dimensions it is asked
// after an exact coplanarity gate and a projection into a plane, because two
// objects in a higher dimension meet only if they share a plane first.
// ---------------------------------------------------------------------------

const AffineKind = enum { none, point, segment, line, coincident, contained };

fn AffineResultOf(comptime T: type) type {
    return struct {
        kind: AffineKind,
        /// The near parameter along the first object, and for `.point` the
        /// only one.
        first: T,
        /// The far parameter along the first object, used by `.segment`.
        second: T,
        /// Which operand the whole intersection is, for `.contained`.
        side: ContainedSide,
    };
}

fn noneAffine(comptime T: type) AffineResultOf(T) {
    return .{ .kind = .none, .first = 0, .second = 0, .side = .first };
}

/// One affine object as the engine sees it: the two points that define its
/// line, the difference, and the closed range of parameters along that
/// difference. A segment is `[0, 1]`, a ray is `[0, +infinity)`, and a line
/// is the whole real line, so all three share one code path.
fn AffineObjectOf(comptime T: type) type {
    return struct {
        start: [2]T,
        finish: [2]T,
        direction: [2]T,
        lower: T,
        upper: T,
    };
}

/// The engine's view of one object. Its `direction` is a DERIVED value, and
/// the exact policy converts it to an integer later, so it is gated here: an
/// infinite difference would reach `fromFloat` as a value it cannot read and
/// would be answered as though it were zero.
fn affineObjectOf(comptime T: type, start: [2]T, finish: [2]T, lower: T, upper: T) error{NonFiniteResult}!AffineObjectOf(T) {
    return .{
        .start = start,
        .finish = finish,
        .direction = try arraySubChecked(2, T, finish, start),
        .lower = lower,
        .upper = upper,
    };
}

fn inRangeOf(comptime T: type, object: AffineObjectOf(T), t: T) bool {
    return t >= object.lower and t <= object.upper;
}

/// The sign of `cross(direction, point - start)`, degree two, over the
/// ORIGINAL stored coordinates. The offset is an exact value here, so a point
/// far from the line keeps the sign its own coordinates give it.
/// `error.NonFiniteResult` is the `fast` policy saying its own evaluation came
/// out a NaN.
fn lineSideSign2(comptime T: type, object: AffineObjectOf(T), point: [2]T, policy: PredicatePolicy) error{NonFiniteResult}!Sign {
    return switch (policy) {
        .fast => blk: {
            const offset = arraySub(2, T, point, object.start);
            break :blk try fastSign(T, subRaw(T, scaleRaw(T, object.direction[0], offset[1]), scaleRaw(T, object.direction[1], offset[0])));
        },
        .adaptive_exact => blk: {
            const direction = exactRow(2, T, object.direction);
            const start_row = exactRow(2, T, object.start);
            const probe = exactRow(2, T, point);
            var offset: [2]Big = undefined;
            exactDeltasInto(2, &offset, &start_row, &probe);
            break :blk exactDet2(&direction[0], &offset[1], &direction[1], &offset[0]);
        },
    };
}

/// True when `point` lies on the supporting line of `object` and its
/// parameter is inside the object's range. The collinearity test is the
/// exact predicate; the range test is a comparison of computed values, and a
/// parameter that left the format is reported rather than read as out of
/// range.
fn pointOnAffine2(comptime T: type, object: AffineObjectOf(T), point: [2]T, policy: PredicatePolicy) error{NonFiniteResult}!bool {
    if (try lineSideSign2(T, object, point, policy) != .zero) return false;
    const denominator = try arrayDotChecked(2, T, object.direction, object.direction);
    if (denominator == 0) return point[0] == object.start[0] and point[1] == object.start[1];
    const offset = try arraySubChecked(2, T, point, object.start);
    const parameter = try requireFiniteResult(T, divRaw(T, try arrayDotChecked(2, T, offset, object.direction), denominator));
    return inRangeOf(T, object, parameter);
}

/// The overlap of two collinear objects, projected onto the axis along which
/// the first object's direction is largest. That axis is the one the shared
/// direction is largest in, so the projection is monotone on both objects
/// and the one-dimensional overlap is the whole answer.
///
/// The projection goes through each object's PARAMETER RANGE, not through its
/// two defining points. A ray's two defining points are one unit apart, so
/// reading its extent off them would clip a ray to a segment, and a ray
/// overlapping a segment on the same line is the case a caller reaches first.
/// The range bounds are what make an infinite end an infinite coordinate, and
/// the four products that reach them are left ungated on purpose: an infinite
/// end of a range is the GEOMETRY, not an overflow, and the only answer that
/// reads one is `.line`, which carries no coordinate.
///
/// The last decision is the one that separates a line from a ray. When the
/// common set is the whole of BOTH objects the two are the same set, but what
/// that set IS depends on the first object's range: a range that runs both
/// ways is a whole line and answers `.line`, while `[0, +infinity)` is a
/// half-line and answers `.coincident`. A ray and its own copy are the same
/// half-line, and reporting that as a line would name points the two operands
/// do not share.
fn collinearOverlap(comptime T: type, a: AffineObjectOf(T), b: AffineObjectOf(T)) AffineResultOf(T) {
    var axis: usize = 0;
    var widest: T = absOf(T, a.direction[0]);
    inline for (1..2) |index| {
        const candidate = absOf(T, a.direction[index]);
        if (candidate > widest) {
            widest = candidate;
            axis = index;
        }
    }
    const a_at_lower = addRaw(T, a.start[axis], scaleRaw(T, a.direction[axis], a.lower));
    const a_at_upper = addRaw(T, a.start[axis], scaleRaw(T, a.direction[axis], a.upper));
    const b_at_lower = addRaw(T, b.start[axis], scaleRaw(T, b.direction[axis], b.lower));
    const b_at_upper = addRaw(T, b.start[axis], scaleRaw(T, b.direction[axis], b.upper));
    const a_low = minOf(T, a_at_lower, a_at_upper);
    const a_high = maxOf(T, a_at_lower, a_at_upper);
    const b_low = minOf(T, b_at_lower, b_at_upper);
    const b_high = maxOf(T, b_at_lower, b_at_upper);
    const low = maxOf(T, a_low, b_low);
    const high = minOf(T, a_high, b_high);
    if (low > high) return noneAffine(T);
    const span = a.direction[axis];
    const first = divRaw(T, subRaw(T, low, a.start[axis]), span);
    const second = divRaw(T, subRaw(T, high, a.start[axis]), span);
    if (low == high) return .{ .kind = .point, .first = first, .second = first, .side = .first };
    const whole_a = low == a_low and high == a_high;
    const whole_b = low == b_low and high == b_high;
    if (whole_a and whole_b) {
        // The two operands ARE the common set. Naming the span tells the caller
        // what that set is, and the span is a whole line only when the first
        // object's range runs both ways.
        const kind: AffineKind = if (a.lower < 0 and a.upper > 0) .line else .coincident;
        return .{ .kind = kind, .first = first, .second = second, .side = .first };
    }
    if (whole_a) return .{ .kind = .contained, .first = first, .second = second, .side = .first };
    if (whole_b) return .{ .kind = .contained, .first = first, .second = second, .side = .second };
    return .{ .kind = .segment, .first = first, .second = second, .side = .first };
}

/// The intersection of two coplanar two-dimensional affine objects.
fn affine2Intersect(comptime T: type, a: AffineObjectOf(T), b: AffineObjectOf(T), policy: PredicatePolicy) error{NonFiniteResult}!AffineResultOf(T) {
    // A degenerate segment is a single point, and then the whole question is
    // whether that point is on the other object.
    if (try arrayLengthChecked(2, T, a.direction) == 0) {
        if (try pointOnAffine2(T, b, a.start, policy)) {
            return .{ .kind = .point, .first = 0, .second = 0, .side = .first };
        }
        return noneAffine(T);
    }
    if (try arrayLengthChecked(2, T, b.direction) == 0) {
        if (try pointOnAffine2(T, a, b.start, policy)) {
            return .{ .kind = .point, .first = 1, .second = 1, .side = .second };
        }
        return noneAffine(T);
    }

    const a_first = try lineSideSign2(T, a, b.start, policy);
    const a_second = try lineSideSign2(T, a, b.finish, policy);
    const b_first = try lineSideSign2(T, b, a.start, policy);
    const b_second = try lineSideSign2(T, b, a.finish, policy);
    const a_collinear = a_first == .zero and a_second == .zero;
    const b_collinear = b_first == .zero and b_second == .zero;
    if (a_collinear or b_collinear) return collinearOverlap(T, a, b);

    // The supporting lines are not parallel, so they meet at exactly one
    // parameter pair, found by Cramer's rule on the two-by-two system. Every
    // intermediate is gated: a denominator or a numerator that left the
    // format would divide into a parameter that means nothing, and the range
    // tests below would read that as "off the end of the object".
    const denominator = try requireFiniteResult(T, subRaw(T, scaleRaw(T, a.direction[0], b.direction[1]), scaleRaw(T, a.direction[1], b.direction[0])));
    if (denominator == 0) return noneAffine(T);
    const gap = try arraySubChecked(2, T, b.start, a.start);
    const first = try requireFiniteResult(T, divRaw(T, try requireFiniteResult(T, subRaw(T, scaleRaw(T, gap[0], b.direction[1]), scaleRaw(T, gap[1], b.direction[0]))), denominator));
    const second = try requireFiniteResult(T, divRaw(T, try requireFiniteResult(T, subRaw(T, scaleRaw(T, gap[0], a.direction[1]), scaleRaw(T, gap[1], a.direction[0]))), denominator));
    if (inRangeOf(T, a, first) and inRangeOf(T, b, second)) {
        return .{ .kind = .point, .first = first, .second = second, .side = .first };
    }

    // The crossing is off the end of at least one object. The only place they
    // can still meet is a boundary of one of them, so each boundary is tested
    // against the other object.
    if (first < a.lower) {
        const candidate = try pointAtArrayChecked(2, T, a.start, a.direction, a.lower);
        if (try pointOnAffine2(T, b, candidate, policy)) {
            return .{ .kind = .point, .first = a.lower, .second = a.lower, .side = .first };
        }
    }
    if (first > a.upper) {
        const candidate = try pointAtArrayChecked(2, T, a.start, a.direction, a.upper);
        if (try pointOnAffine2(T, b, candidate, policy)) {
            return .{ .kind = .point, .first = a.upper, .second = a.upper, .side = .first };
        }
    }
    if (second < b.lower) {
        const candidate = try pointAtArrayChecked(2, T, b.start, b.direction, b.lower);
        if (try pointOnAffine2(T, a, candidate, policy)) {
            return .{ .kind = .point, .first = 0, .second = b.lower, .side = .second };
        }
    }
    if (second > b.upper) {
        const candidate = try pointAtArrayChecked(2, T, b.start, b.direction, b.upper);
        if (try pointOnAffine2(T, a, candidate, policy)) {
            return .{ .kind = .point, .first = 1, .second = b.upper, .side = .second };
        }
    }
    return noneAffine(T);
}

/// True when the two segments SHARE a flat in four dimensions, which is what
/// the four-dimensional affine intersections open on.
///
/// The three columns `a_finish - a_start`, `b_finish - b_start`, and
/// `b_start - a_start` are formed ONCE and the four three-by-four minors that
/// vanish exactly when those three vectors are linearly dependent all read
/// them, so the gate converts the four endpoints once rather than once per
/// minor. `fast` forms the same three differences in `T`, which is what `fast`
/// means, and reports `error.NonFiniteResult` when one of its own minors came
/// out a NaN.
fn coplanar4(
    comptime T: type,
    a_start: [4]T,
    a_finish: [4]T,
    b_start: [4]T,
    b_finish: [4]T,
    policy: PredicatePolicy,
) error{NonFiniteResult}!bool {
    return switch (policy) {
        .fast => blk: {
            const u = arraySub(4, T, a_finish, a_start);
            const v = arraySub(4, T, b_finish, b_start);
            const w = arraySub(4, T, b_start, a_start);
            break :blk !(try fastSign(T, fastMinor3x4(T, u, v, w, .{ 0, 1, 2 })) != .zero or
                try fastSign(T, fastMinor3x4(T, u, v, w, .{ 0, 1, 3 })) != .zero or
                try fastSign(T, fastMinor3x4(T, u, v, w, .{ 0, 2, 3 })) != .zero or
                try fastSign(T, fastMinor3x4(T, u, v, w, .{ 1, 2, 3 })) != .zero);
        },
        .adaptive_exact => blk: {
            const columns = exactMinor3Columns(T, a_start, a_finish, b_start, b_finish);
            // SAFETY: minor3ExactInto writes metadata and active limbs before sign extraction.
            var minor: Big = undefined;
            minor3ExactInto(&minor, 4, &columns, .{ 0, 1, 2 });
            if (Big.sign(&minor) != .zero) break :blk false;
            minor3ExactInto(&minor, 4, &columns, .{ 0, 1, 3 });
            if (Big.sign(&minor) != .zero) break :blk false;
            minor3ExactInto(&minor, 4, &columns, .{ 0, 2, 3 });
            if (Big.sign(&minor) != .zero) break :blk false;
            minor3ExactInto(&minor, 4, &columns, .{ 1, 2, 3 });
            break :blk Big.sign(&minor) == .zero;
        },
    };
}

/// The sign of the three-dimensional signed volume of the tetrahedron
/// `a`, `b`, `c`, `d`, which is zero exactly when the four points are
/// coplanar. That coplanarity is the gate the three-dimensional affine
/// intersections open on. `error.NonFiniteResult` is the `fast` policy saying
/// its own evaluation came out a NaN.
fn orient3Sign(comptime T: type, a: [3]T, b: [3]T, c: [3]T, d: [3]T, policy: PredicatePolicy) error{NonFiniteResult}!Sign {
    return switch (policy) {
        .fast => fastSign(T, crossDotFast(T, a, b, c, d)),
        .adaptive_exact => exactOrient3(T, a, b, c, d),
    };
}

/// A coordinate axis whose projection out of the whole basis is largest, so
/// the residual is far from zero and the normalization is well conditioned.
/// The axis is chosen against EVERY basis vector, not the first one alone:
/// an axis that is least aligned with the first vector can still lie wholly
/// inside the span of the others, and a seed inside the span projects to
/// exactly zero, which normalizes to the zero vector rather than a direction.
/// `basis` must be shorter than the dimension, so some axis is outside it.
fn seedOutsideBasis(comptime dim: usize, comptime T: type, basis: []const [dim]T) error{NonFiniteResult}![dim]T {
    var seed: [dim]T = undefined;
    inline for (0..dim) |index| {
        seed[index] = 0;
    }
    var best: usize = 0;
    var best_weight: T = std.math.inf(T);
    inline for (0..dim) |index| {
        // The squared projection of the axis onto the span, so the SMALLEST
        // weight names the axis that lies most outside it. A weight that left
        // the format would compare false against the running best and leave
        // the choice on the seed, which is a direction the span may contain.
        var weight: T = 0;
        for (basis) |vector| {
            weight = try requireFiniteResult(T, addRaw(T, weight, scaleRaw(T, vector[index], vector[index])));
        }
        if (weight < best_weight) {
            best_weight = weight;
            best = index;
        }
    }
    seed[best] = 1;
    return orthogonalizedChecked(dim, T, seed, basis);
}

/// The intersection of two coplanar affine objects in three or four
/// dimensions. The coplanarity gate is the exact predicate on the ORIGINAL
/// coordinates; the classification after the projection is the exact
/// two-dimensional predicate on the PROJECTED coordinates, which is what the
/// module comment promises and no more.
///
/// Every frame coordinate is gated. The projection feeds the exact
/// two-dimensional predicate, which converts those coordinates to integers,
/// and an infinite one would be read there as though it were zero.
fn affineIntersectHigh(
    comptime dim: usize,
    comptime T: type,
    a_start: [dim]T,
    a_finish: [dim]T,
    a_lower: T,
    a_upper: T,
    b_start: [dim]T,
    b_finish: [dim]T,
    b_lower: T,
    b_upper: T,
    policy: PredicatePolicy,
) error{NonFiniteResult}!AffineResultOf(T) {
    if (dim == 3) {
        if (try orient3Sign(T, a_start, a_finish, b_start, b_finish, policy) != .zero) return noneAffine(T);
    } else {
        if (!try coplanar4(T, a_start, a_finish, b_start, b_finish, policy)) return noneAffine(T);
    }

    // The gate is decided above, so the frame the projection needs is built
    // only once the two objects are known to share a flat.
    const u = try arraySubChecked(dim, T, a_finish, a_start);
    const v = try arraySubChecked(dim, T, b_finish, b_start);

    const first = try arrayNormalizedChecked(dim, T, u);
    if (try arrayLengthChecked(dim, T, first) == 0) return noneAffine(T);
    const residual = try arraySubChecked(dim, T, v, try arrayScaleChecked(dim, T, first, try arrayDotChecked(dim, T, v, first)));
    const second = if (try arrayLengthChecked(dim, T, residual) > 0)
        try arrayNormalizedChecked(dim, T, residual)
    else
        try anyPerpendicularChecked(dim, T, first);

    const first_object = try affineObjectOf(T, [2]T{
        try arrayDotChecked(dim, T, a_start, first),
        try arrayDotChecked(dim, T, a_start, second),
    }, [2]T{
        try arrayDotChecked(dim, T, a_finish, first),
        try arrayDotChecked(dim, T, a_finish, second),
    }, a_lower, a_upper);
    const second_object = try affineObjectOf(T, [2]T{
        try arrayDotChecked(dim, T, b_start, first),
        try arrayDotChecked(dim, T, b_start, second),
    }, [2]T{
        try arrayDotChecked(dim, T, b_finish, first),
        try arrayDotChecked(dim, T, b_finish, second),
    }, b_lower, b_upper);
    return affine2Intersect(T, first_object, second_object, policy);
}

/// Builds the concrete union from the engine's whole-set classification.
/// `collinearOverlap` distinguishes two-sided lines from coincident rays.
/// Every payload coordinate is checked before it reaches the caller.
fn affineToIntersection(
    comptime dim: usize,
    comptime T: type,
    a_start: [dim]T,
    a_direction: [dim]T,
    result: AffineResultOf(T),
) error{NonFiniteResult}!IntersectionOf(dim, T) {
    const P = PointOf(dim, T);
    switch (result.kind) {
        .none => return .none,
        .point => return .{ .point = fromArrayOf(P, try pointAtArrayChecked(dim, T, a_start, a_direction, result.first)) },
        .segment => return .{ .segment = .{
            .a = fromArrayOf(P, try pointAtArrayChecked(dim, T, a_start, a_direction, result.first)),
            .b = fromArrayOf(P, try pointAtArrayChecked(dim, T, a_start, a_direction, result.second)),
        } },
        .line => {
            if (try arrayLengthChecked(dim, T, a_direction) == 0) return .none;
            return .{ .line = .{
                .origin = fromArrayOf(P, a_start),
                .direction = fromArrayOf(VectorOf(dim, T), try arrayNormalizedChecked(dim, T, a_direction)),
            } };
        },
        // A segment or a ray has the same set as the other operand. The caller
        // already holds it, so this tag needs no computed payload.
        .coincident => return .coincident,
        .contained => return .{ .contained = result.side },
    }
}

/// The one route every affine intersection takes. The direction the payload
/// needs is derived here rather than taken from the operand, because the
/// engine's own objects are two points and a parameter range.
fn affineIntersectDispatch(
    comptime dim: usize,
    comptime T: type,
    a_start: [dim]T,
    a_finish: [dim]T,
    a_lower: T,
    a_upper: T,
    b_start: [dim]T,
    b_finish: [dim]T,
    b_lower: T,
    b_upper: T,
    policy: PredicatePolicy,
) error{NonFiniteResult}!IntersectionOf(dim, T) {
    const direction = try arraySubChecked(dim, T, a_finish, a_start);
    const result = if (dim == 2)
        try affine2Intersect(T, try affineObjectOf(T, a_start, a_finish, a_lower, a_upper), try affineObjectOf(T, b_start, b_finish, b_lower, b_upper), policy)
    else
        try affineIntersectHigh(dim, T, a_start, a_finish, a_lower, a_upper, b_start, b_finish, b_lower, b_upper, policy);
    return affineToIntersection(dim, T, a_start, direction, result);
}

// ---------------------------------------------------------------------------
// Ball intersections
// ---------------------------------------------------------------------------

/// The classification of two balls, and the point or the pair of points that
/// goes with it. Every coordinate it builds is gated, because each of them is
/// a value the caller receives.
fn ballIntersect(comptime dim: usize, comptime T: type, first_center: [dim]T, first_radius: T, second_center: [dim]T, second_radius: T, policy: PredicatePolicy) error{NonFiniteResult}!IntersectionOf(dim, T) {
    const P = PointOf(dim, T);
    const separation = try arraySubChecked(dim, T, second_center, first_center);
    // The two comparisons are the two tangency conditions, and their SIGNS run
    // the way the geometry does: the centres are farther apart than the radii
    // add when `outer` is positive, and closer than they differ when `inner`
    // is negative. Getting either sign backwards turns every classification
    // into its opposite, which is the mistake the tests below exist to catch.
    if (try arrayLengthChecked(dim, T, separation) == 0 and first_radius == second_radius) return .coincident;
    const tangency = try ballTangencySigns(dim, T, first_center, first_radius, second_center, second_radius, policy);
    if (tangency.outer == .positive) return .none;
    if (tangency.outer == .zero or tangency.inner == .zero) {
        // Tangency: one point on the line of centres, at distance `first`
        // from the first centre along the separation.
        const unit = try arrayNormalizedChecked(dim, T, separation);
        const touch = try arrayAddChecked(dim, T, first_center, try arrayScaleChecked(dim, T, unit, first_radius));
        return .{ .point = fromArrayOf(P, touch) };
    }
    if (tangency.inner == .negative) return .{ .contained = if (first_radius >= second_radius) .second else .first };
    return .{ .two_points = try ballCrossingPoints(dim, T, first_center, first_radius, second_center, second_radius) };
}

/// The two points where two balls cut each other. The geometry is fixed: the
/// chord lies in the plane through both centres, perpendicular to the line of
/// centres, at the distance the law of cosines puts it, and its half-length is
/// what `chordHeight` computes.
fn ballCrossingPoints(comptime dim: usize, comptime T: type, first_center: [dim]T, first_radius: T, second_center: [dim]T, second_radius: T) error{NonFiniteResult}![2]PointOf(dim, T) {
    const P = PointOf(dim, T);
    const separation = try arraySubChecked(dim, T, second_center, first_center);
    const distance = try arrayLengthChecked(dim, T, separation);
    const unit = try arrayScaleChecked(dim, T, separation, divRaw(T, 1, distance));
    const along = try requireFiniteResult(T, divRaw(T, try requireFiniteResult(T, addRaw(T, subRaw(T, scaleRaw(T, first_radius, first_radius), scaleRaw(T, second_radius, second_radius)), scaleRaw(T, distance, distance))), scaleRaw(T, 2, distance)));
    const height = try chordHeight(T, first_radius, second_radius, distance);
    const base = try arrayAddChecked(dim, T, first_center, try arrayScaleChecked(dim, T, unit, along));
    const side = try arrayScaleChecked(dim, T, try anyPerpendicularChecked(dim, T, unit), height);
    return .{
        fromArrayOf(P, try arrayAddChecked(dim, T, base, side)),
        fromArrayOf(P, try arraySubChecked(dim, T, base, side)),
    };
}

/// The half-length of the chord two balls cut, from the law of cosines.
///
/// The plain form of the third side is `radius * radius - along * along`, and
/// `along` is ROUNDED to the format before that subtraction, so a small ball
/// against the rim of a huge one leaves a half-chord of exactly zero and the
/// two reported points on top of each other. The factored form is the same
/// quantity with no square of a radius in it, because the product of the four
/// tangency margins IS `4 * distance * height * height`:
///
///     4 d h^2 = (r1 + r2 - d) (r1 + r2 + d) (d - r1 + r2) (d + r1 - r2)
///
/// Every margin is a difference of stored numbers, and each is taken against
/// the LARGER of the two it subtracts from, so no margin is ever read out of a
/// rounded `r1 + r2`. The two SMALL margins are each below `2 * distance` and
/// the two LARGE ones each below `3 * distance`, so dividing each PAIR by
/// `2 * distance` before multiplying leaves both quotients below
/// `5 * distance`, and the only value of the size of the answer is the final
/// product, which IS the half-chord. That is what keeps the calculation inside
/// the format at either end of the range. A margin that rounding puts below
/// zero leaves no chord, and its sign is dropped the way the plain form drops
/// the sign of its difference. Every product and quotient is gated, so a pair
/// of radii and a separation whose margins leave the format is
/// `error.NonFiniteResult` rather than a chord made of infinities.
fn chordHeight(comptime T: type, first_radius: T, second_radius: T, distance: T) error{NonFiniteResult}!T {
    const total = addRaw(T, first_radius, second_radius);
    const difference = subRaw(T, first_radius, second_radius);
    const outer_low = if (first_radius >= distance)
        addRaw(T, subRaw(T, first_radius, distance), second_radius)
    else if (second_radius >= distance)
        addRaw(T, subRaw(T, second_radius, distance), first_radius)
    else
        subRaw(T, total, distance);
    const outer_high = addRaw(T, total, distance);
    const inner_low = subRaw(T, distance, difference);
    const inner_high = addRaw(T, distance, difference);
    const divisor = scaleRaw(T, 2, distance);
    return requireFiniteResult(T, scaleRaw(
        T,
        sqrtOf(T, absOf(T, try requireFiniteResult(T, divRaw(T, try requireFiniteResult(T, scaleRaw(T, outer_low, inner_low)), divisor)))),
        sqrtOf(T, absOf(T, try requireFiniteResult(T, divRaw(T, try requireFiniteResult(T, scaleRaw(T, outer_high, inner_high)), divisor)))),
    ));
}

/// The common set of a segment and a ball. A segment whose BOTH ends are
/// inside lies in the ball, because the ball is convex, so it answers
/// `.contained = .second` and stops there. Otherwise the segment is framed as
/// the ray along it and the ray quadratic decides.
fn ballIntersectSegment(comptime dim: usize, comptime T: type, center: [dim]T, radius: T, a: [dim]T, b: [dim]T, policy: PredicatePolicy) error{NonFiniteResult}!IntersectionOf(dim, T) {
    const P = PointOf(dim, T);
    const inside = try segmentInsideSigns(dim, T, center, radius, a, b, policy);
    // The segment is the second operand of `intersectSegment`, and it is the
    // whole of the common set, so the side is `second`.
    if (inside.start != .positive and inside.end != .positive) return .{ .contained = .second };
    // Both ends are outside, so the rounded offset and the frame below are
    // built only for the case that reaches them.
    const from_a = try arraySubChecked(dim, T, a, center);
    const span = try arrayLengthChecked(dim, T, try arraySubChecked(dim, T, b, a));
    if (span == 0) return .none;
    // The unit frame the quadratic needs. It feeds `rayBallDiscriminantSign`,
    // whose exact branch converts it, so it is gated before it gets there.
    const unit = try arrayScaleChecked(dim, T, try arraySubChecked(dim, T, b, a), divRaw(T, 1, span));
    const along = try arrayDotChecked(dim, T, from_a, unit);
    const discriminant = try rayBallDiscriminantSign(dim, T, center, a, unit, radius, policy);
    if (discriminant == .negative) return .none;
    const root = try requireFiniteResult(T, sqrtOf(T, absOf(T, try requireFiniteResult(T, subRaw(T, scaleRaw(T, along, along), subRaw(T, try arrayDotChecked(dim, T, from_a, from_a), scaleRaw(T, radius, radius)))))));
    // The two roots of the ray quadratic are `-along -+ root`. The smaller one
    // enters the ball first, so it is the `low` end of the overlap.
    const near = try requireFiniteResult(T, divRaw(T, subRaw(T, negRaw(T, along), root), span));
    const far = try requireFiniteResult(T, divRaw(T, addRaw(T, negRaw(T, along), root), span));
    const low = maxOf(T, near, 0);
    const high = minOf(T, far, 1);
    if (low > high) return .none;
    const span_vector = try arraySubChecked(dim, T, b, a);
    if (low == high) return .{ .point = fromArrayOf(P, try pointAtArrayChecked(dim, T, a, span_vector, low)) };
    return .{ .segment = .{
        .a = fromArrayOf(P, try pointAtArrayChecked(dim, T, a, span_vector, low)),
        .b = fromArrayOf(P, try pointAtArrayChecked(dim, T, a, span_vector, high)),
    } };
}

/// The common set of a ray and a ball. The ray quadratic is
/// `t^2 + 2 * along * t + (dot(offset, offset) - r * r)` with the direction
/// already unit, so the roots are `-along -+ root`.
///
/// The three shapes follow from where the ORIGIN sits, because the ball is
/// the CLOSED solid ball and the ray is `[0, +infinity)`:
///
/// - `near >= 0`: the closed segment from entry to exit, or a single tangent
///   point when the endpoints coincide.
/// - `near < 0`: the origin is inside, so the common set is the whole run
///   from the origin to `far`, which is a SEGMENT and not a point. The
///   origin itself is a common point, so reporting only the exit would drop
///   it, and reporting only the origin would drop the run between them.
///   `far == 0` is the boundary case, where that run has zero length.
fn ballIntersectRay(comptime dim: usize, comptime T: type, center: [dim]T, radius: T, origin: [dim]T, direction: [dim]T, policy: PredicatePolicy) error{NonFiniteResult}!IntersectionOf(dim, T) {
    const P = PointOf(dim, T);
    const discriminant = try rayBallDiscriminantSign(dim, T, center, origin, direction, radius, policy);
    if (discriminant == .negative) return .none;
    // The discriminant is settled, so the rounded offset and projection that
    // place the two roots are built only for the case that has them.
    const offset = try arraySubChecked(dim, T, origin, center);
    const along = try arrayDotChecked(dim, T, offset, direction);
    const root = try requireFiniteResult(T, sqrtOf(T, absOf(T, try requireFiniteResult(T, subRaw(T, scaleRaw(T, along, along), subRaw(T, try arrayDotChecked(dim, T, offset, offset), scaleRaw(T, radius, radius)))))));
    // Both roots carry `-along`. The smaller is the near root, so it is
    // `-along - root`; the larger is `-along + root`.
    const near = try requireFiniteResult(T, negRaw(T, addRaw(T, along, root)));
    const far = try requireFiniteResult(T, addRaw(T, negRaw(T, along), root));
    // The whole chord lies behind the origin, which no point of the ray
    // reaches.
    if (far < 0) return .none;
    if (near >= 0) {
        if (near == far) return .{ .point = fromArrayOf(P, try pointAtArrayChecked(dim, T, origin, direction, near)) };
        return .{ .segment = .{
            .a = fromArrayOf(P, try pointAtArrayChecked(dim, T, origin, direction, near)),
            .b = fromArrayOf(P, try pointAtArrayChecked(dim, T, origin, direction, far)),
        } };
    }
    // The origin is inside the ball, so the run starts at the origin. A run
    // of zero length is a single point, which is the case of a ray that
    // starts on the surface and points away from the ball.
    if (far == 0) return .{ .point = fromArrayOf(P, origin) };
    return .{ .segment = .{
        .a = fromArrayOf(P, origin),
        .b = fromArrayOf(P, try pointAtArrayChecked(dim, T, origin, direction, far)),
    } };
}

/// The common set of a ball and a flat, and the cut when they are not
/// tangent. The payload carries the foot of the perpendicular from the centre
/// and the radius of the section, both read INSIDE the flat, which is why the
/// caller keeps the flat it passed: a centre and a radius say nothing about
/// which flat they sit in.
fn ballIntersectFlat(comptime dim: usize, comptime T: type, center: [dim]T, radius: T, flat_point: [dim]T, flat_normal: [dim]T, policy: PredicatePolicy) error{NonFiniteResult}!IntersectionOf(dim, T) {
    const sign = try flatBallSign(dim, T, flat_point, center, flat_normal, radius, policy);
    if (sign == .positive) return .none;
    // The signed distance is settled, so the rounded one that places the foot
    // and sizes the cut is built only for the case that uses them.
    const distance = try arrayDotChecked(dim, T, try arraySubChecked(dim, T, center, flat_point), flat_normal);
    const foot = try arraySubChecked(dim, T, center, try arrayScaleChecked(dim, T, flat_normal, distance));
    if (sign == .zero) return .{ .point = fromArrayOf(PointOf(dim, T), foot) };
    const cut_radius = try requireFiniteResult(T, sqrtOf(T, absOf(T, try requireFiniteResult(T, subRaw(T, scaleRaw(T, radius, radius), scaleRaw(T, distance, distance))))));
    const ball = SphereOf(dim, T){
        .center = fromArrayOf(PointOf(dim, T), foot),
        .radius = cut_radius,
    };
    return if (dim == 3) .{ .circle = ball } else .{ .hypersphere = ball };
}

// ---------------------------------------------------------------------------
// Flat intersections
// ---------------------------------------------------------------------------

/// The common set of a segment and a flat. Both ends on the flat means the
/// whole segment is (`.contained = .second`), one end on it means that end,
/// and the ends strictly opposite means the segment crosses at one parameter.
/// The parameter is the share of the two rounded signed distances, and it is
/// gated along with the point it builds.
fn flatIntersectSegment(comptime dim: usize, comptime T: type, flat_point: [dim]T, flat_normal: [dim]T, a: [dim]T, b: [dim]T, policy: PredicatePolicy) error{NonFiniteResult}!IntersectionOf(dim, T) {
    const P = PointOf(dim, T);
    const side = try segmentSideSigns(dim, T, flat_point, flat_normal, a, b, policy);
    const first = side.start;
    const second = side.end;
    if (first == .zero and second == .zero) return .{ .contained = .second };
    if (first == .zero) return .{ .point = fromArrayOf(P, a) };
    if (second == .zero) return .{ .point = fromArrayOf(P, b) };
    if ((first == .negative) == (second == .negative)) return .none;
    // The flat cuts the segment strictly inside, so the parameter is the
    // share of the two signed distances.
    const left = flatDistanceFast(dim, T, a, flat_point, flat_normal);
    const right = flatDistanceFast(dim, T, b, flat_point, flat_normal);
    const t = try requireFiniteResult(T, divRaw(T, try requireFiniteResult(T, left), try requireFiniteResult(T, subRaw(T, left, right))));
    return .{ .point = fromArrayOf(P, try pointAtArrayChecked(dim, T, a, try arraySubChecked(dim, T, b, a), t)) };
}

/// True when the two normals are parallel, decided on the two-by-two minors
/// of the matrix whose ROWS are the normals. A minor uses one component of
/// each normal against a DIFFERENT component of the other, which is what
/// `minor2Sign(first[i], second[j], first[j], second[i])` reads. Passing the
/// same index twice on each normal instead tests a relation between two
/// components of one normal and two of the other, which vanishes for every
/// coordinate-aligned pair and would call the two axis normals parallel.
fn normalsParallel3(comptime T: type, first: [3]T, second: [3]T, policy: PredicatePolicy) error{NonFiniteResult}!bool {
    return switch (policy) {
        .fast => try minor2Sign(T, first[0], second[1], first[1], second[0]) == .zero and
            try minor2Sign(T, first[1], second[2], first[2], second[1]) == .zero and
            try minor2Sign(T, first[2], second[0], first[0], second[2]) == .zero,
        // Each of the six components is read by two of the three minors, so
        // the exact path converts each one once and hands it to both.
        .adaptive_exact => blk: {
            const f = exactRow(3, T, first);
            const s = exactRow(3, T, second);
            break :blk exactDet2(&f[0], &s[1], &f[1], &s[0]) == .zero and
                exactDet2(&f[1], &s[2], &f[2], &s[1]) == .zero and
                exactDet2(&f[2], &s[0], &f[0], &s[2]) == .zero;
        },
    };
}

/// The four-dimensional gate, over all six index pairs. The same argument
/// order as the three-dimensional one, and all six minors vanish exactly when
/// the two rows are linearly dependent.
fn normalsParallel4(comptime T: type, first: [4]T, second: [4]T, policy: PredicatePolicy) error{NonFiniteResult}!bool {
    return switch (policy) {
        .fast => blk: {
            var parallel = true;
            inline for (0..4) |row| {
                inline for (row + 1..4) |column| {
                    if (try minor2Sign(T, first[row], second[column], first[column], second[row]) != .zero) {
                        parallel = false;
                    }
                }
            }
            break :blk parallel;
        },
        // Every component is read by three of the six minors, so the exact
        // path converts the eight coordinates once and hands them to all six.
        .adaptive_exact => blk: {
            const f = exactRow(4, T, first);
            const s = exactRow(4, T, second);
            var parallel = true;
            inline for (0..4) |row| {
                inline for (row + 1..4) |column| {
                    if (exactDet2(&f[row], &s[column], &f[column], &s[row]) != .zero) {
                        parallel = false;
                    }
                }
            }
            break :blk parallel;
        },
    };
}

/// The sign of `a * b - c * d`, degree two: the plain evaluation of one
/// two-by-two minor of the matrix whose ROWS are two normals. The exact
/// counterpart is `exactDet2` over cached anchors, and
/// `error.NonFiniteResult` is the plain evaluation coming out a NaN.
fn minor2Sign(comptime T: type, a: T, b: T, c: T, d: T) error{NonFiniteResult}!Sign {
    return fastSign(T, subRaw(T, scaleRaw(T, a, b), scaleRaw(T, c, d)));
}

/// The intersection of two planes in three dimensions.
///
/// When the normals are not parallel the meeting set is a line whose
/// direction is the cross product of the normals. The point on that line is
/// found by sliding along the component of the second normal that is
/// orthogonal to the first. That component is the right transversal to use
/// precisely because it is NOT orthogonal to the second normal: its dot
/// product with the second normal is `1 - (n1 . n2)^2`, which a non-parallel
/// pair of unit normals never makes zero. A transversal chosen to be
/// orthogonal to both normals instead would divide by zero on every pair
/// whose normals are coordinate-aligned, which is most of them.
fn flatPair3(comptime T: type, first_point: [3]T, first_normal: [3]T, second_point: [3]T, second_normal: [3]T, policy: PredicatePolicy) error{NonFiniteResult}!Intersection3(T) {
    if (try normalsParallel3(T, first_normal, second_normal, policy)) {
        if (try flatSideSign(3, T, second_point, first_point, first_normal, policy) != .zero) return .none;
        return .{ .plane = .{ .origin = fromArrayOf(Point3(T), first_point), .normal = fromArrayOf(Vector3(T), first_normal) } };
    }
    const direction = try arrayNormalizedChecked(3, T, cross3Raw(T, first_normal, second_normal));
    const transversal = try arrayNormalizedChecked(3, T, try arraySubChecked(3, T, second_normal, try arrayScaleChecked(3, T, first_normal, try arrayDotChecked(3, T, second_normal, first_normal))));
    const gap = try arrayDotChecked(3, T, try arraySubChecked(3, T, second_point, first_point), second_normal);
    const step = try requireFiniteResult(T, divRaw(T, gap, try arrayDotChecked(3, T, transversal, second_normal)));
    const meeting = try arrayAddChecked(3, T, first_point, try arrayScaleChecked(3, T, transversal, step));
    return .{ .line = .{
        .origin = fromArrayOf(Point3(T), meeting),
        .direction = fromArrayOf(Vector3(T), direction),
    } };
}

/// The intersection of two hyperplanes in four dimensions. Parallel normals
/// give the whole three-flat when the pair is one hyperplane and nothing when
/// it is not; otherwise the common set is a two-flat, whose payload is a point
/// on it plus two orthonormal directions spanning it.
fn flatPair4(comptime T: type, first_point: [4]T, first_normal: [4]T, second_point: [4]T, second_normal: [4]T, policy: PredicatePolicy) error{NonFiniteResult}!Intersection4(T) {
    if (try normalsParallel4(T, first_normal, second_normal, policy)) {
        if (try flatSideSign(4, T, second_point, first_point, first_normal, policy) != .zero) return .none;
        return .{ .hyperplane = .{ .origin = fromArrayOf(Point4(T), first_point), .normal = fromArrayOf(Vector4(T), first_normal) } };
    }
    // `in_plane_normal` is the part of the second normal orthogonal to the
    // first. It slides between the two hyperplanes, and it is orthogonal to
    // the first by construction, so the step lands on both.
    const in_plane_normal = try arrayNormalizedChecked(4, T, try arraySubChecked(4, T, second_normal, try arrayScaleChecked(4, T, first_normal, try arrayDotChecked(4, T, second_normal, first_normal))));
    const gap = try arrayDotChecked(4, T, try arraySubChecked(4, T, second_point, first_point), second_normal);
    const step = try requireFiniteResult(T, divRaw(T, gap, try arrayDotChecked(4, T, in_plane_normal, second_normal)));
    const meeting = try arrayAddChecked(4, T, first_point, try arrayScaleChecked(4, T, in_plane_normal, step));
    // The meeting set has direction span the two-dimensional orthogonal
    // complement of both normals, and in four dimensions that complement is
    // exactly two-dimensional, so two directions describe it.
    const first_direction = try seedOutsideBasis(4, T, &[_][4]T{ first_normal, in_plane_normal });
    const second_direction = try seedOutsideBasis(4, T, &[_][4]T{ first_normal, in_plane_normal, first_direction });
    return .{ .plane = .{
        .origin = fromArrayOf(Point4(T), meeting),
        .u = fromArrayOf(Vector4(T), first_direction),
        .v = fromArrayOf(Vector4(T), second_direction),
    } };
}
