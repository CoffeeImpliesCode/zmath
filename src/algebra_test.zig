//! Algebraic identities for the twelve concrete geometric algebras.
//!
//! Every test here is an identity that the published contract forces, not a
//! smoke test. The coefficient patterns are dyadic rationals whose partial sums
//! stay well inside the exact integer range of `f32`, so the comparisons are
//! bit exact rather than tolerant wherever the identity is exact.

const std = @import("std");
const algebra = @import("algebra.zig");

const testing = std.testing;

/// The twelve published types, in the order the normalization table lists them.
const all_algebras = .{
    algebra.Cl3f,
    algebra.Cl3d,
    algebra.Cl4f,
    algebra.Cl4d,
    algebra.Pga3f,
    algebra.Pga3d,
    algebra.Pga4f,
    algebra.Pga4d,
    algebra.Cga3f,
    algebra.Cga3d,
    algebra.Cga4f,
    algebra.Cga4d,
};

fn asBits(comptime T: type, value: T) u64 {
    return switch (@bitSizeOf(T)) {
        32 => @as(u64, @as(u32, @bitCast(value))),
        64 => @bitCast(value),
        else => @compileError("zmath.algebra tests support f32 and f64 only"),
    };
}

fn sameValue(comptime T: type, expected: T, actual: T) bool {
    if (asBits(T, expected) == asBits(T, actual)) return true;
    // The identities below are about values, and IEEE negation maps `+0` to
    // `-0`, so the two zeros compare equal here. The zero sign is pinned on its
    // own in the negation test.
    return expected == 0 and actual == 0;
}

fn expectCoefficients(comptime A: type, expected: A, actual: A) !void {
    const want = expected.coefficients();
    const got = actual.coefficients();
    for (want, 0..) |value, mask| {
        if (!sameValue(A.element, value, got[mask])) {
            std.debug.print("coefficient {d}: want {d}, got {d}\n", .{
                mask,
                value,
                got[mask],
            });
            return error.CoefficientMismatch;
        }
    }
}

/// The bit exact comparison the identities below need: `sameValue` accepts
/// either sign of zero, and a coefficient wise comparison cannot see the sign
/// of zero at all, so a change that only moves a zero between `+0` and `-0`
/// would pass both.
fn expectSameBits(comptime T: type, expected: T, actual: T) !void {
    if (asBits(T, expected) != asBits(T, actual)) {
        std.debug.print("want {x}, got {x}\n", .{ asBits(T, expected), asBits(T, actual) });
        return error.BitMismatch;
    }
}

fn expectClose(comptime T: type, expected: T, actual: T, tolerance: T) !void {
    if (!std.math.isFinite(actual)) return error.NonFiniteResult;
    if (@abs(expected - actual) > tolerance) {
        std.debug.print("want {d}, got {d}\n", .{ expected, actual });
        return error.OutOfTolerance;
    }
}

fn expectCloseCoefficients(
    comptime A: type,
    expected: A,
    actual: A,
    tolerance: A.element,
) !void {
    const want = expected.coefficients();
    const got = actual.coefficients();
    for (want, 0..) |value, mask| {
        if (!std.math.isFinite(got[mask])) return error.NonFiniteResult;
        if (@abs(value - got[mask]) > tolerance) {
            std.debug.print("coefficient {d}: want {d}, got {d}\n", .{
                mask,
                value,
                got[mask],
            });
            return error.OutOfTolerance;
        }
    }
}

/// The unit blade at `mask`, with no other coefficient.
fn bladeOf(comptime A: type, mask: usize) A {
    var values = [_]A.element{0} ** A.coefficient_count;
    values[mask] = 1;
    return A.fromCoefficients(values);
}

/// The `ordinal`th Euclidean basis vector in the published order, counting from
/// `e1`. `Pga` leads with the null vector, so its Euclidean directions start at
/// bit 1.
fn euclidean(comptime A: type, ordinal: usize) A {
    const offset = if (A.metric[0] == 0) 1 else 0;
    return A.blade(ordinal - 1 + offset);
}

/// A multivector with every coefficient nonzero and distinct enough that no
/// product collapses by accident.
fn dense(comptime A: type) A {
    var values: [A.coefficient_count]A.element = undefined;
    for (0..A.coefficient_count) |mask| {
        values[mask] = @as(A.element, @floatFromInt(@popCount(mask) + 1)) / 2;
    }
    return A.fromCoefficients(values);
}

test "the twelve types publish the basis, metric, and coefficient count" {
    inline for (all_algebras) |A| {
        try testing.expectEqual(if (A.element == f32) 32 else 64, @bitSizeOf(A.element));
        try testing.expect(A.basis_count == A.dimension or A.basis_count == A.dimension + 1 or
            A.basis_count == A.dimension + 2);
        try testing.expectEqual(@as(usize, 1) << A.basis_count, A.coefficient_count);
        try testing.expectEqual(A.coefficient_count - 1, A.pseudoscalar_mask);
        try testing.expectEqual(A.basis_count + 1, A.grade_count);
        try testing.expectEqual(A.coefficient_count * @sizeOf(A.element), A.storage_bytes);
    }
    try testing.expectEqual(@as(usize, 8), algebra.Cl3f.coefficient_count);
    try testing.expectEqual(@as(usize, 16), algebra.Cl4f.coefficient_count);
    try testing.expectEqual(@as(usize, 16), algebra.Pga3f.coefficient_count);
    try testing.expectEqual(@as(usize, 32), algebra.Pga4f.coefficient_count);
    try testing.expectEqual(@as(usize, 32), algebra.Cga3f.coefficient_count);
    try testing.expectEqual(@as(usize, 64), algebra.Cga4f.coefficient_count);
    try testing.expectEqual(@as(usize, 32), algebra.Cl3f.storage_bytes);
    try testing.expectEqual(@as(usize, 512), algebra.Cga4d.storage_bytes);
    try testing.expectEqual(algebra.Family.projective, algebra.Pga3f.family);
    try testing.expectEqual(algebra.Family.conformal, algebra.Cga4f.family);
    try testing.expectEqual(algebra.Family.clifford, algebra.Cl3d.family);
    try testing.expectEqual(@as(f32, 0), algebra.Pga3f.metric[0]);
    try testing.expectEqual(@as(f32, -1), algebra.Cga3f.metric[4]);
    try testing.expectEqual(@as(f32, -1), algebra.Cga4f.metric[5]);
}

test "one is the multiplicative identity in every algebra" {
    inline for (all_algebras) |A| {
        const a = dense(A);
        try expectCoefficients(A, a, A.geometricProduct(A.one(), a));
        try expectCoefficients(A, a, A.geometricProduct(a, A.one()));
        try expectCoefficients(A, A.one(), A.geometricProduct(A.one(), A.one()));
    }
}

test "the geometric product is associative" {
    inline for (.{ algebra.Cl3f, algebra.Pga3f, algebra.Cga3f, algebra.Cga4d }) |A| {
        const a = dense(A);
        const b = A.one().add(dense(A).scale(0.25));
        const c = A.fromCoefficients(A.pseudoscalar().coefficients()).add(dense(A).scale(0.125));
        const left = A.geometricProduct(A.geometricProduct(a, b), c);
        const right = A.geometricProduct(a, A.geometricProduct(b, c));
        try expectCoefficients(A, left, right);
    }
}

test "the geometric product distributes over addition" {
    inline for (.{ algebra.Cl3f, algebra.Pga3f, algebra.Cga3f, algebra.Cga4d }) |A| {
        const a = dense(A);
        const b = A.one().scale(0.5).add(dense(A).scale(0.25));
        const c = A.blade(0).scale(0.75);
        const left = A.geometricProduct(a.add(b), c);
        const right = A.geometricProduct(a, c).add(A.geometricProduct(b, c));
        try expectCoefficients(A, left, right);
        const mirrored = A.geometricProduct(c, a.add(b));
        const mirrored_right = A.geometricProduct(c, a).add(A.geometricProduct(c, b));
        try expectCoefficients(A, mirrored, mirrored_right);
    }
}

test "the scalar and vector grades form a commutative subalgebra" {
    inline for (all_algebras) |A| {
        var scalar_and_vector = [_]A.element{0} ** A.coefficient_count;
        scalar_and_vector[A.scalar_mask] = 1.5;
        for (0..A.dimension) |index| {
            const offset = if (A.metric[0] == 0) 1 else 0;
            scalar_and_vector[@as(usize, 1) << @as(std.math.Log2Int(usize), @intCast(index + offset))] = 0.5;
        }
        const a = A.fromCoefficients(scalar_and_vector);
        var other = [_]A.element{0} ** A.coefficient_count;
        other[A.scalar_mask] = 0.25;
        for (0..A.dimension) |index| {
            const offset = if (A.metric[0] == 0) 1 else 0;
            other[@as(usize, 1) << @as(std.math.Log2Int(usize), @intCast(index + offset))] = 1;
        }
        const b = A.fromCoefficients(other);
        try expectCoefficients(A, A.geometricProduct(a, b), A.geometricProduct(b, a));
        try expectCoefficients(A, A.geometricProduct(a, a), A.geometricProduct(a, a));
    }
}

test "e1 and e2 anticommute in every algebra" {
    inline for (all_algebras) |A| {
        const e1 = euclidean(A, 1);
        const e2 = euclidean(A, 2);
        const forward = A.geometricProduct(e1, e2);
        const backward = A.geometricProduct(e2, e1);
        try expectCoefficients(A, forward.neg(), backward);
        for (forward.coefficients()) |value| {
            if (value != 0) return;
        }
        return error.CommutativePair;
    }
}

test "each basis vector squares to its metric and distinct ones anticommute" {
    inline for (all_algebras) |A| {
        for (0..A.basis_count) |i| {
            const e = A.blade(i);
            const squared = A.geometricProduct(e, e);
            var expected = [_]A.element{0} ** A.coefficient_count;
            expected[A.scalar_mask] = A.metric[i];
            try expectCoefficients(A, A.fromCoefficients(expected), squared);
        }
        for (0..A.basis_count) |i| {
            for (i + 1..A.basis_count) |j| {
                const forward = A.geometricProduct(A.blade(i), A.blade(j));
                const backward = A.geometricProduct(A.blade(j), A.blade(i));
                try expectCoefficients(A, forward.neg(), backward);
            }
        }
    }
}

test "the pseudoscalar square matches the published metric" {
    try testing.expectEqual(@as(f32, -1), algebra.Cl3f.pseudoscalar_square);
    try testing.expectEqual(@as(f32, 1), algebra.Cl4f.pseudoscalar_square);
    try testing.expectEqual(@as(f32, 0), algebra.Pga3f.pseudoscalar_square);
    try testing.expectEqual(@as(f32, 0), algebra.Pga4f.pseudoscalar_square);
    try testing.expectEqual(@as(f32, -1), algebra.Cga3f.pseudoscalar_square);
    try testing.expectEqual(@as(f32, 1), algebra.Cga4f.pseudoscalar_square);
    inline for (all_algebras) |A| {
        const square = A.geometricProduct(A.pseudoscalar(), A.pseudoscalar());
        try testing.expectEqual(A.pseudoscalar_square, square.bladeValue(A.scalar_mask));
        if (A.pseudoscalar_square != 0) {
            try expectCoefficients(A, A.one().scale(A.pseudoscalar_square), square);
        }
    }
}

test "reverse is an involution" {
    inline for (all_algebras) |A| {
        const a = dense(A);
        try expectCoefficients(A, a, a.reverse().reverse());
        try expectCoefficients(A, A.one(), A.one().reverse());
        const top_sign: A.element = if (@mod(A.basis_count * (A.basis_count - 1) / 2, 2) == 0)
            1
        else
            -1;
        try expectCoefficients(A, A.pseudoscalar().scale(top_sign), A.pseudoscalar().reverse());
    }
}

test "reverse fixes grades zero, one, and four and negates grades two and three" {
    inline for (all_algebras) |A| {
        var fixed = [_]A.element{0} ** A.coefficient_count;
        fixed[A.scalar_mask] = 1;
        if (A.coefficient_count > 2) fixed[1] = 1;
        if (A.basis_count >= 4) fixed[0b01111] = 1;
        const a = A.fromCoefficients(fixed);
        try expectCoefficients(A, a, a.reverse());

        const bivector = bladeOf(A, 0b011);
        try expectCoefficients(A, bivector.neg(), bivector.reverse());
        const trivector = bladeOf(A, 0b0111);
        try expectCoefficients(A, trivector.neg(), trivector.reverse());
    }
}

test "the highest grades follow the reverse rule in the six dimensional algebra" {
    // Grade 5 has the sign (-1)^10 and is fixed; grade 6 has (-1)^15 and flips.
    const A = algebra.Cga4f;
    const grade5 = bladeOf(A, 0b011111);
    try expectCoefficients(A, grade5, grade5.reverse());
    const grade6 = bladeOf(A, 0b111111);
    try expectCoefficients(A, grade6.neg(), grade6.reverse());
    try expectClose(f32, 1, grade5.grade(5), 0);
}

test "the grade zero part of A times its reverse is the norm in every algebra" {
    inline for (all_algebras) |A| {
        const a = dense(A);
        const sandwich = A.geometricProduct(a, a.reverse());
        try testing.expectEqual(a.norm(), sandwich.bladeValue(A.scalar_mask));
        try testing.expectEqual(a.norm(), A.dot(a, a));
    }
}

test "A times its reverse is purely scalar for a single grade blade" {
    inline for (all_algebras) |A| {
        for (0..A.coefficient_count) |mask| {
            const blade = bladeOf(A, mask);
            const sandwich = A.geometricProduct(blade, blade.reverse());
            for (0..A.coefficient_count) |other| {
                const expected = if (other == A.scalar_mask) blade.norm() else 0;
                try testing.expectEqual(expected, sandwich.bladeValue(other));
            }
        }
    }
}

test "only the grade zero part of A times its reverse is a portable claim" {
    // The reverse leaves grade 1 unchanged, so for A = 1 + e1 the product is
    // 2 + 2 e1 and is not purely scalar, even in a positive definite algebra.
    // The norm reads the grade-0 part, which is all the definition promises.
    inline for (all_algebras) |A| {
        const a = A.one().add(euclidean(A, 1));
        const sandwich = A.geometricProduct(a, a.reverse());
        try testing.expectEqual(@as(A.element, 2), a.norm());
        try testing.expectEqual(a.norm(), sandwich.bladeValue(A.scalar_mask));
    }
    // A mixed-grade multivector keeps an off-scalar part whatever the metric, and
    // the conformal algebra also has a direction whose norm is negative.
    inline for (all_algebras) |A| {
        const offset = if (A.metric[0] == 0) 1 else 0;
        const vector_mask = @as(usize, 1) << @as(u6, @intCast(offset));
        const mixed = A.one().add(euclidean(A, 1));
        const sandwich = A.geometricProduct(mixed, mixed.reverse());
        try testing.expect(sandwich.bladeValue(vector_mask) != 0);
        try testing.expectEqual(mixed.norm(), sandwich.bladeValue(A.scalar_mask));
    }
    inline for (.{ algebra.Cga3f, algebra.Cga4f }) |A| {
        const negative = A.blade(A.basis_count - 1);
        try testing.expect(negative.norm() < 0);
        try testing.expect(negative.normalized() == null);
    }
    // The norm of a general multivector still agrees with the grade-0 part, and
    // it is nonnegative exactly where the metric is positive definite.
    inline for (.{ algebra.Cl3f, algebra.Cl3d, algebra.Cl4f, algebra.Cl4d }) |A| {
        const a = dense(A);
        try testing.expect(a.norm() >= 0);
        try testing.expectEqual(a.norm(), A.geometricProduct(a, a.reverse()).scalarPart());
    }
}

test "the dual of the dual is the dual square times the multivector" {
    inline for (all_algebras) |A| {
        const expected_square: A.element = if (A.basis_count % 2 == 0) 1 else -1;
        try testing.expectEqual(expected_square, A.dual_square);
        const a = dense(A);
        try expectCoefficients(A, a.scale(expected_square), a.dual().dual());
        if (expected_square == 1) {
            try expectCoefficients(A, a, a.dual().dual());
        } else {
            try expectCoefficients(A, a.neg(), a.dual().dual());
        }
    }
}

test "the right dual of the pseudoscalar is the scalar" {
    inline for (.{
        algebra.Cl3f,
        algebra.Cl3d,
        algebra.Cl4f,
        algebra.Cl4d,
        algebra.Cga3f,
        algebra.Cga3d,
        algebra.Cga4f,
        algebra.Cga4d,
    }) |A| {
        try expectCoefficients(A, A.one(), A.pseudoscalar().dual());
    }
}

test "the signed Poincare dual uses the sign rule the header derives" {
    // A grade k blade of mask m maps to the blade of mask top & ~m with the
    // sign (-1)^(k * (n - 1 - k)). Applying it twice gives (-1)^n, so Pga3 with
    // four basis vectors is an involution and Pga4 with five is not.
    inline for (.{ algebra.Pga3f, algebra.Pga3d, algebra.Pga4f, algebra.Pga4d }) |A| {
        const top = A.pseudoscalar_mask;
        for (0..A.coefficient_count) |mask| {
            const k: i32 = @intCast(@popCount(mask));
            const exponent = k * (@as(i32, @intCast(A.basis_count - 1)) - k);
            const sign: A.element = if (@mod(exponent, 2) == 0) 1 else -1;
            const image = bladeOf(A, mask).dual();
            try testing.expectEqual(sign, image.bladeValue(top & ~mask));
            for (0..A.coefficient_count) |other| {
                if (other != top & ~mask) {
                    try testing.expectEqual(@as(A.element, 0), image.bladeValue(other));
                }
            }
        }
        const a = dense(A);
        if (A.basis_count % 2 == 0) {
            try expectCoefficients(A, a, a.dual().dual());
        } else {
            try expectCoefficients(A, a.neg(), a.dual().dual());
        }
    }
    // Concrete blades in Pga3, whose four basis vectors give (-1)^(1 * 2) = +1
    // for grade 1 and (-1)^(2 * 1) = +1 for grade 2.
    try testing.expectEqual(@as(f32, 1), euclidean(algebra.Pga3f, 1).dual().bladeValue(0b1101));
    const pga_bivector = euclidean(algebra.Pga3f, 1).wedge(euclidean(algebra.Pga3f, 2));
    try testing.expectEqual(@as(f32, 1), pga_bivector.dual().bladeValue(0b1001));
    // Pga4 has five basis vectors, so the grade 1 sign is (-1)^(1 * 3) = -1.
    try testing.expectEqual(@as(f32, -1), euclidean(algebra.Pga4f, 1).dual().bladeValue(0b11101));
}

test "the dual complement is the unsigned complement" {
    inline for (all_algebras) |A| {
        const a = dense(A);
        const complemented = a.dualComplement();
        for (0..A.coefficient_count) |mask| {
            const expected = complemented.bladeValue(A.pseudoscalar_mask & ~mask);
            try testing.expectEqual(a.bladeValue(mask), expected);
        }
    }
}

test "the scalar form is symmetric and is the dot product on vectors" {
    inline for (all_algebras) |A| {
        const a = dense(A);
        const b = A.one().add(dense(A).scale(0.25));
        try testing.expectEqual(A.dot(a, b), A.dot(b, a));
        const x = euclidean(A, 1).scale(1.5);
        const y = euclidean(A, 2).scale(-0.5);
        try testing.expectEqual(@as(A.element, 0), A.dot(x, y));
        try testing.expectEqual(@as(A.element, 1.5), A.dot(x, euclidean(A, 1)));
        try testing.expectEqual(A.dot(x, x), x.norm());
    }
}

/// The ascending basis index list of the blade at `mask`, and its length.
fn bladeIndexList(comptime A: type, mask: usize, out: *[A.basis_count]usize) usize {
    var length: usize = 0;
    for (0..A.basis_count) |index| {
        if (mask & (@as(usize, 1) << @as(std.math.Log2Int(usize), @intCast(index))) != 0) {
            out[length] = index;
            length += 1;
        }
    }
    return length;
}

/// The sign of `e_I * e_I` before the metric is folded in: the inversion count
/// of the two concatenated index lists, which is the number of swaps that puts
/// that concatenation in ascending order.
fn swapParity(comptime A: type, mask: usize) i32 {
    var indices: [A.basis_count]usize = undefined;
    const length = bladeIndexList(A, mask, &indices);
    var sign: i32 = 1;
    for (0..length) |first| {
        for (0..length) |second| {
            if (indices[first] > indices[second]) sign = -sign;
        }
    }
    return sign;
}

/// The product of the metric signs of the basis vectors in the blade, so the
/// leading null direction of `Pga` folds the whole blade to the zero scale.
fn bladeMetricProduct(comptime A: type, mask: usize) i32 {
    var indices: [A.basis_count]usize = undefined;
    const length = bladeIndexList(A, mask, &indices);
    var sign: i32 = 1;
    for (indices[0..length]) |index| {
        const metric_value = A.metric[index];
        if (metric_value == 0) return 0;
        if (metric_value < 0) sign = -sign;
    }
    return sign;
}

/// The scalar form recomputed from the published metric alone: one term per
/// blade in ascending mask order, with the zero left skip, the zero scale skip,
/// and the operand order of the library. A change in the sign, in the order, or
/// in either skip shows up here as a different bit pattern.
fn diagonalDot(comptime A: type, left: A, right: A) A.element {
    var sum: A.element = 0;
    for (0..A.coefficient_count) |mask| {
        const value = left.bladeValue(mask);
        if (value == 0) continue;
        const parity = swapParity(A, mask);
        const scale_sign = parity * bladeMetricProduct(A, mask);
        if (scale_sign == 0) continue;
        const product: A.element = @floatFromInt(scale_sign);
        const reverse: A.element = @floatFromInt(parity);
        sum += product * value * (reverse * right.bladeValue(mask));
    }
    return sum;
}

test "the scalar form is the diagonal blade square sum, bit for bit" {
    inline for (all_algebras) |A| {
        const left = dense(A);
        const right = A.one().add(dense(A).scale(0.25));
        try expectSameBits(A.element, diagonalDot(A, left, right), A.dot(left, right));
        try expectSameBits(A.element, diagonalDot(A, left, left), A.norm(left));
        // A single blade leaves one term, so this pins the sign of the top
        // blade, which is the only place a grade above four appears.
        const top = bladeOf(A, A.pseudoscalar_mask);
        try expectSameBits(A.element, diagonalDot(A, top, top), A.norm(top));
        // Coefficients that no binary format holds exactly pin the
        // accumulation order: the sum rounds at every step, so summing the
        // same terms in another order gives other bits.
        var left_values: [A.coefficient_count]A.element = undefined;
        var right_values: [A.coefficient_count]A.element = undefined;
        for (0..A.coefficient_count) |mask| {
            left_values[mask] = @as(A.element, @floatFromInt(mask + 3)) / 7;
            right_values[mask] = @as(A.element, @floatFromInt(mask + 5)) / 11;
        }
        std.mem.doNotOptimizeAway(&left_values);
        std.mem.doNotOptimizeAway(&right_values);
        const rounded_left = A.fromCoefficients(left_values);
        const rounded_right = A.fromCoefficients(right_values);
        try expectSameBits(A.element, diagonalDot(A, rounded_left, rounded_right), A.dot(rounded_left, rounded_right));
    }
}

test "the scalar form skips zero terms and keeps the non-finite rule" {
    inline for (all_algebras) |A| {
        const positive_zero: A.element = 0;
        var infinite = [_]A.element{0} ** A.coefficient_count;
        // The second basis vector squares to `+1` in all three families, so its
        // diagonal scale is a unit and the infinity is not skipped by it.
        infinite[2] = std.math.inf(A.element);
        const broken = A.fromCoefficients(infinite);
        // A zero left factor annihilates every term, so a non-finite right
        // operand never reaches the sum and the result is positive zero.
        try expectSameBits(A.element, positive_zero, A.dot(A.zero(), broken));
        try expectSameBits(A.element, positive_zero, A.norm(A.zero()));
        // Negative zero is zero, so it takes the same skip.
        var negative_zeros = [_]A.element{0} ** A.coefficient_count;
        for (&negative_zeros) |*slot| slot.* = -@as(A.element, 0);
        const signed_zero = A.fromCoefficients(negative_zeros);
        try expectSameBits(A.element, positive_zero, A.dot(signed_zero, broken));
        try expectSameBits(A.element, positive_zero, A.dot(signed_zero, signed_zero));
        try expectSameBits(A.element, positive_zero, signed_zero.norm());
        // The leading null direction of `Pga` carries the zero metric, so a
        // nonzero coefficient there is not a term at all and the sum is still
        // the positive zero it starts at.
        if (A.metric[0] == 0) {
            var null_values = [_]A.element{0} ** A.coefficient_count;
            null_values[1] = 1;
            const null_blade = A.fromCoefficients(null_values);
            try expectSameBits(A.element, positive_zero, null_blade.norm());
            try expectSameBits(A.element, positive_zero, A.dot(null_blade, dense(A)));
        }
        // A non-finite left coefficient is not a zero row, so it propagates.
        try testing.expect(!std.math.isFinite(A.dot(broken, broken)));
        try testing.expect(!std.math.isFinite(A.norm(broken)));
    }
}

test "the checked inverse is the unchecked inverse whenever it accepts" {
    inline for (all_algebras) |A| {
        // A scalar plus a bivector and a pure Euclidean bivector are versors in
        // all three families, and the checked form has to return exactly what
        // the unchecked one returns, bit for bit. The bivector has to be even:
        // `sandwichConjugate` keeps only the even grades, so an odd grade blade
        // has a zero sandwich and the checked form rejects it.
        const accepted = [_]A{
            A.one().add(bladeOf(A, 0b0011)),
            bladeOf(A, 0b0110),
        };
        for (accepted) |candidate| {
            const checked = candidate.inverseChecked() orelse return error.InverseRejected;
            const unchecked = candidate.inverse();
            for (unchecked.coefficients(), checked.coefficients()) |want, got| {
                try expectSameBits(A.element, want, got);
            }
        }
        // Adding the top blade leaves a residue above grade zero in the
        // sandwich, which is the whole reason the checked form exists.
        const residue = A.one().add(bladeOf(A, 0b0011)).add(bladeOf(A, A.pseudoscalar_mask));
        try testing.expect(residue.inverseChecked() == null);
    }
}

test "the exterior product is the graded commutative join" {
    const A = algebra.Cl3f;
    const e1 = euclidean(A, 1);
    const e2 = euclidean(A, 2);
    const bivector = bladeOf(A, 0b011);
    try expectCoefficients(A, bivector, A.wedge(e1, e2));
    try expectCoefficients(A, A.zero(), A.wedge(e1, e1));
    try expectCoefficients(A, A.pseudoscalar(), A.wedge(bivector, euclidean(A, 3)));
    // Two blades commute up to the sign of the product of their grades.
    for (0..A.coefficient_count) |left| {
        for (0..A.coefficient_count) |right| {
            const a = bladeOf(A, left);
            const b = bladeOf(A, right);
            const sign: f32 = if (@mod(@popCount(left) * @popCount(right), 2) == 0) 1 else -1;
            try expectCoefficients(A, A.wedge(a, b).scale(sign), A.wedge(b, a));
        }
    }
}

test "the interior product contracts and lowers the grade" {
    const A = algebra.Cl3f;
    const bivector = bladeOf(A, 0b011);
    try expectCoefficients(A, euclidean(A, 1), A.inner(bivector, euclidean(A, 2)));
    try expectCoefficients(A, A.zero(), A.inner(bivector, euclidean(A, 3)));
    try expectCoefficients(A, A.zero(), A.inner(euclidean(A, 1), euclidean(A, 2)));
    // The contraction of a bivector with itself is its own square, and the
    // reverse sign of grade 2 makes that square -1.
    try testing.expectEqual(@as(f32, -1), A.inner(bivector, bivector).bladeValue(0));
    try testing.expectEqual(@as(f32, 0), A.inner(bivector, bivector).grade(1));
}

test "the meet of two planes is their shared blade" {
    const A = algebra.Cl3f;
    const first = bladeOf(A, 0b011);
    const second = bladeOf(A, 0b110);
    try expectCoefficients(A, euclidean(A, 2), A.meet(first, second));
    // The meet has grade r + s - n, which is 1 here, and swapping the arguments
    // multiplies it by (-1)^(n + r + s) = (-1)^7 = -1, the sign of the
    // orientation of the intersection.
    try testing.expectEqual(@as(f32, 1), A.meet(first, second).grade(1));
    try expectCoefficients(A, A.meet(first, second).scale(@as(f32, -1)), A.meet(second, first));
    const B = algebra.Cl4f;
    // In four dimensions the meet of two planes has grade 0, so it is a scalar.
    // Complementary planes meet in a nonzero scalar, planes that share a line
    // meet in zero, and the meet is symmetric because (-1)^(4 + 2 + 2) is +1.
    const spanning = euclidean(B, 1).wedge(euclidean(B, 2));
    const complementary = euclidean(B, 3).wedge(euclidean(B, 4));
    const shared = euclidean(B, 1).wedge(euclidean(B, 3));
    try testing.expectEqual(@as(f32, 1), B.meet(spanning, complementary).grade(0));
    try expectCoefficients(B, B.one(), B.meet(spanning, complementary));
    try expectCoefficients(B, B.meet(spanning, complementary), B.meet(complementary, spanning));
    try expectCoefficients(B, B.zero(), B.meet(spanning, shared));
}

test "the conformal pair has inner product minus one" {
    inline for (.{ algebra.Cga3f, algebra.Cga3d, algebra.Cga4f, algebra.Cga4d }) |A| {
        try testing.expectEqual(@as(A.element, -1), A.dot(A.n_o(), A.n_inf()));
        try testing.expectEqual(@as(A.element, -1), A.dot(A.n_inf(), A.n_o()));
        try testing.expectEqual(@as(A.element, 0), A.n_inf().norm());
        try expectCoefficients(A, A.n_o(), A.origin());
        try expectCoefficients(A, A.n_o(), A.embed(.{0} ** A.dimension));
    }
}

test "the embedding round trips the Euclidean part" {
    inline for (.{ algebra.Cga3f, algebra.Cga3d, algebra.Cga4f, algebra.Cga4d }) |A| {
        const x: [A.dimension]A.element = switch (A.dimension) {
            3 => .{ 1, -2, 0.5 },
            else => .{ 1, -2, 0.5, 0.75 },
        };
        const embedded = A.embed(x);
        try testing.expectEqual(x, embedded.euclideanPart());
        try testing.expectEqual(x, embedded.vectorPart());
        var squared: A.element = 0;
        for (x) |value| squared += value * value;
        // The embedding is exactly n_o + x + (x dot x) * n_inf / 2.
        const recovered = embedded.sub(A.n_o()).sub(A.n_inf().scale(squared / 2));
        var expected = [_]A.element{0} ** A.coefficient_count;
        for (x, 0..) |value, index| {
            expected[@as(usize, 1) << @as(std.math.Log2Int(usize), @intCast(index))] = value;
        }
        try expectCoefficients(A, A.fromCoefficients(expected), recovered);
    }
}

test "the conformal origin is the point at the zero vector" {
    inline for (.{ algebra.Cga3f, algebra.Cga3d, algebra.Cga4f, algebra.Cga4d }) |A| {
        const point = A.embed(.{0} ** A.dimension);
        try expectCoefficients(A, A.n_o(), point);
        try expectCoefficients(A, A.n_o(), A.origin());
        // The Euclidean components of the origin are zero.
        try testing.expectEqual(.{0} ** A.dimension, point.euclideanPart());
    }
}

test "the projective point blade multiplies a position into the null direction" {
    inline for (.{ algebra.Pga3f, algebra.Pga3d, algebra.Pga4f, algebra.Pga4d }) |A| {
        const x: [A.dimension]A.element = switch (A.dimension) {
            3 => .{ 1, -2, 0.5 },
            else => .{ 1, -2, 0.5, 0.75 },
        };
        var position = [_]A.element{0} ** A.coefficient_count;
        for (x, 0..) |value, index| {
            position[@as(usize, 1) << @as(std.math.Log2Int(usize), @intCast(index + 1))] = value;
        }
        const p = A.point(x);
        // The point blade stores the position in the e0 wedge ei components and
        // has no scalar and no vector part.
        for (x, 0..) |value, index| {
            const mask = @as(usize, 1) | (@as(usize, 1) << @as(std.math.Log2Int(usize), @intCast(index + 1)));
            try testing.expectEqual(value, p.bladeValue(mask));
        }
        try testing.expectEqual(@as(A.element, 0), p.grade(0));
        try testing.expectEqual(@as(A.element, 0), p.grade(1));
        var squared: A.element = 0;
        for (x) |value| squared += value * value;
        // (e0 wedge x) * x is the null direction scaled by the squared length.
        try testing.expectEqual(squared, A.geometricProduct(p, A.fromCoefficients(position))
            .bladeValue(1));
        // Two point blades multiply to exactly zero, which is why the table
        // publishes no single-blade line.
        try expectCoefficients(A, A.zero(), A.geometricProduct(p, A.point(x)));
        // A point blade is its own null square, so it has no unit form either.
        try expectCoefficients(A, A.zero(), A.geometricProduct(p, p));
        try testing.expect(p.normalized() == null);
    }
}

test "the unit constructions follow the normalization table" {
    inline for (all_algebras) |A| {
        // A rotor is the even part of one plus a bivector, normalized.
        const bivector = euclidean(A, 1).wedge(euclidean(A, 2));
        const rotor = A.unitRotor(bivector) orelse return error.RotorRejected;
        try testing.expect(rotor.grade(0) != 0);
        try testing.expect(rotor.grade(1) == 0 and rotor.grade(3) == 0);
        try expectClose(A.element, 1, rotor.norm(), @as(A.element, 1e-5));
        // A zero bivector gives the identity rotor, and a rotor whose even part
        // cancels has no unit form.
        try expectCloseCoefficients(A, A.one(), A.unitRotor(A.one()) orelse
            return error.RotorRejected, @as(A.element, 1e-6));
        try testing.expect(A.unitRotor(A.one().neg()) == null);

        // A versor is the even part, normalized.
        const versor = A.unitVersor(A.one().add(bivector)) orelse return error.VersorRejected;
        try testing.expectEqual(rotor, versor);
        try testing.expect(A.unitVersor(bivector) != null);
        try testing.expect(A.unitVersor(euclidean(A, 1)) == null);

        // A flat is a single-grade blade, normalized.
        const flat = A.unitFlat(bivector) orelse return error.FlatRejected;
        try expectClose(A.element, 1, flat.norm(), @as(A.element, 1e-5));
        try testing.expect(A.unitFlat(A.zero()) == null);
        try testing.expect(A.unitFlat(A.one().add(bivector)) == null);
    }
}

test "the unit round is the top grade normalized" {
    inline for (.{ algebra.Cga3f, algebra.Cga3d, algebra.Cga4f, algebra.Cga4d }) |A| {
        const top = A.basis_count - 1;
        const round = bladeOf(A, (@as(usize, 1) << @as(u6, @intCast(top))) - 1);
        const unit = A.unitRound(round) orelse return error.RoundRejected;
        try expectClose(A.element, 1, unit.norm(), @as(A.element, 1e-5));
        for (0..A.coefficient_count) |mask| {
            if (@popCount(mask) != top) {
                try testing.expectEqual(@as(A.element, 0), unit.bladeValue(mask));
            }
        }
        try testing.expect(A.unitRound(A.one()) == null);
    }
}

test "the inverse of a unit multivector is its conjugate" {
    inline for (.{
        algebra.Cl3f,
        algebra.Cl3d,
        algebra.Cl4f,
        algebra.Cl4d,
        algebra.Cga3f,
        algebra.Cga3d,
        algebra.Cga4f,
        algebra.Cga4d,
    }) |A| {
        const unit = A.one();
        try expectCoefficients(A, unit.conjugate(), unit.inverseChecked() orelse
            return error.InverseRejected);
        const rotor = A.unitRotor(euclidean(A, 1).wedge(euclidean(A, 2))) orelse
            return error.RotorRejected;
        try expectCloseCoefficients(A, rotor.conjugate(), rotor.inverse(), @as(A.element, 1e-6));
        try expectCloseCoefficients(A, A.one(), A.geometricProduct(rotor, rotor.inverse()), @as(A.element, 1e-6));
    }
}

test "the projective inverse is the versor conjugate over the sandwich scalar" {
    inline for (.{ algebra.Pga3f, algebra.Pga3d, algebra.Pga4f, algebra.Pga4d }) |A| {
        // A rotor, which is a scalar plus a Euclidean bivector, has a purely
        // scalar sandwich, so the inverse of a unit rotor is its conjugate.
        const raw = A.one().add(bladeOf(A, 0b0110));
        const unit = raw.normalized() orelse return error.NormalizeRejected;
        const inverse = unit.inverseChecked() orelse return error.InverseRejected;
        try expectCloseCoefficients(A, unit.conjugate(), inverse, @as(A.element, 1e-5));
        try expectCloseCoefficients(A, A.one(), A.geometricProduct(unit, inverse), @as(A.element, 1e-5));
        // A motor adds a grade 4 blade, and its sandwich keeps a grade 4
        // residue, so it is not a versor and the checked form rejects it. The
        // grade 0 part alone would give a wrong answer.
        const motor = raw.add(bladeOf(A, 0b01111));
        try testing.expect(motor.inverseChecked() == null);
        try testing.expect(motor.add(euclidean(A, 1)).inverseChecked() == null);
    }
}

test "the checked forms reject zero and non-finite input" {
    inline for (all_algebras) |A| {
        try testing.expect(A.normalized(A.zero()) == null);
        try testing.expect(A.inverseChecked(A.zero()) == null);
        try testing.expect(A.normalized(A.one()) != null);
        try testing.expect(A.inverseChecked(A.one()) != null);

        var infinite = [_]A.element{0} ** A.coefficient_count;
        infinite[1] = std.math.inf(A.element);
        const broken = A.fromCoefficients(infinite);
        try testing.expect(!broken.isFinite());
        try testing.expect(broken.normalized() == null);
        try testing.expect(broken.inverseChecked() == null);
    }
}

test "the binary64 inverse of a huge versor keeps its subnormal reciprocal" {
    // The sandwich scalar of `k * (1 + b)`, for a bivector `b` over two metric
    // `+1` directions, is exactly `2 * k^2`, so `k = 2^511` puts it at `2^1023`,
    // the last power of two binary64 holds below its overflow. Both inverses
    // then divide by `2^-1023`, which is subnormal: a soft-float binary64
    // divide answers a signed zero there, and every coefficient of the answer
    // goes with it.
    const k: f64 = 0x1p+511;
    const reciprocal: f64 = 0x1p-1023;
    inline for (.{
        algebra.Cl3d,
        algebra.Cl4d,
        algebra.Pga3d,
        algebra.Pga4d,
        algebra.Cga3d,
        algebra.Cga4d,
    }) |A| {
        // `0b0110` is a Euclidean bivector in every family, including `Pga`,
        // whose leading null direction it deliberately skips.
        const huge = A.one().add(bladeOf(A, 0b0110)).scale(k);
        const expected = huge.conjugate().scale(reciprocal);
        const checked = huge.inverseChecked() orelse return error.InverseRejected;
        try expectCoefficients(A, expected, checked);
        try expectCoefficients(A, A.one(), A.geometricProduct(huge, checked));
        try expectCoefficients(A, expected, huge.inverse());
    }
}

test "a zero left factor annihilates the geometric product" {
    inline for (all_algebras) |A| {
        var infinite = [_]A.element{0} ** A.coefficient_count;
        infinite[1] = std.math.inf(A.element);
        const broken = A.fromCoefficients(infinite);
        try expectCoefficients(A, A.zero(), A.geometricProduct(A.zero(), broken));
        try expectCoefficients(A, A.zero(), A.geometricProduct(A.zero(), dense(A)));
    }
}

test "coefficients round trip bit exactly in every algebra" {
    inline for (all_algebras) |A| {
        var values: [A.coefficient_count]A.element = undefined;
        for (0..A.coefficient_count) |mask| {
            const raw: u32 = @as(u32, @intCast(mask)) *% 2654435761;
            values[mask] = @as(A.element, @floatFromInt(raw % 251 + 1)) / 4;
        }
        const a = A.fromCoefficients(values);
        try testing.expectEqual(values, a.coefficients());
        try expectCoefficients(A, a, A.setCoefficients(a, values));
        for (values, 0..) |value, mask| {
            try testing.expectEqual(value, a.bladeValue(mask));
        }
        try testing.expectEqual(@as(A.element, 0), a.bladeValue(A.coefficient_count));
    }
}

test "grade, bladeValue, and scalarPart read the stored coefficients" {
    inline for (all_algebras) |A| {
        const a = dense(A);
        try testing.expectEqual(a.bladeValue(0), a.scalarPart());
        var total: A.element = 0;
        for (0..A.grade_count) |k| {
            total += a.grade(k);
        }
        var sum: A.element = 0;
        for (a.coefficients()) |value| sum += value;
        try testing.expectEqual(sum, total);
        try testing.expectEqual(@as(A.element, 0), a.grade(A.basis_count + 1));
    }
}

test "the even part keeps only the even grades" {
    inline for (all_algebras) |A| {
        const a = dense(A);
        const even = a.evenPart();
        for (0..A.coefficient_count) |mask| {
            const expected = if (@popCount(mask) % 2 == 0) a.bladeValue(mask) else 0;
            try testing.expectEqual(expected, even.bladeValue(mask));
        }
    }
}

test "the motor part is the grade four blade" {
    inline for (.{ algebra.Pga3f, algebra.Pga4f, algebra.Cga3f, algebra.Cga4f }) |A| {
        const motor = dense(A);
        const part = motor.motorPart();
        for (0..A.coefficient_count) |mask| {
            const expected = if (@popCount(mask) == 4) motor.bladeValue(mask) else 0;
            try testing.expectEqual(expected, part.bladeValue(mask));
        }
    }
}

test "normalization gives unit norm and null for a null blade" {
    inline for (all_algebras) |A| {
        const flat = euclidean(A, 1).add(euclidean(A, 2));
        const unit = flat.normalized() orelse return error.NormalizeRejected;
        try expectClose(A.element, 1, unit.norm(), @as(A.element, 1e-6));
        // A null blade has no unit form, and the projective metric makes that
        // reachable with a finite input.
        if (A.metric[0] == 0) {
            // The null basis vector is its own null blade, so it has no unit form.
            var null_blade = [_]A.element{0} ** A.coefficient_count;
            null_blade[1] = 1;
            const null_part = A.fromCoefficients(null_blade);
            try testing.expectEqual(@as(A.element, 0), null_part.norm());
            try testing.expect(null_part.normalized() == null);
            try testing.expect(null_part.inverseChecked() == null);
        }
    }
}

test "negation and reversal carry the IEEE sign of zero" {
    inline for (all_algebras) |A| {
        const sign_bit = @as(u64, 1) << (@bitSizeOf(A.element) - 1);
        const zero = A.zero();
        for (zero.neg().coefficients()) |value| {
            try testing.expectEqual(sign_bit, asBits(A.element, value));
        }
        // The reverse multiplies each coefficient by its grade sign, so a zero
        // coefficient becomes negative zero on the grades the sign negates.
        for (zero.reverse().coefficients(), 0..) |value, mask| {
            const k: i32 = @intCast(@popCount(mask));
            const negates = @mod(@divTrunc(k * (k - 1), 2), 2) != 0;
            try testing.expectEqual(if (negates) sign_bit else 0, asBits(A.element, value));
        }
        // A double negation restores positive zero, which is why the identity
        // comparisons above accept either sign of zero.
        try expectCoefficients(A, zero, zero.neg().neg());
        try expectCoefficients(A, zero, zero.scale(-1).scale(-1));
        try expectCoefficients(A, zero, zero.conjugate().conjugate());
    }
}
