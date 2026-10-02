// Correctly rounded inverse hyperbolic sine for binary16 values.
//
// Ported from CORE-MATH src/binary16/asinh/asinhf16.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa (Copyright (c) 2025 Paul Zimmermann).
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
// The caller must keep IEEE round-to-nearest, ties-to-even active.

const fma_impl = @import("scalar_fma.zig");

const sign_mask: u16 = 0x8000;
const magnitude_mask: u16 = 0x7fff;
const infinity_bits: u16 = 0x7c00;
const quiet_nan_bit: u16 = 0x0200;

const coefficients = [19][7]f32{
    .{
        0x1.000002p+0, -0x1.5554d6p-3, 0x1.32c4a4p-4, -0x1.5336f2p-5,
        0.0,           0.0,            0.0,
    },
    .{
        0x1.707d7ap-17, 0x1.ffe5ap-1,  0x1.86a284p-10, -0x1.608874p-3,
        0x1.1d36dp-7,   0x1.435f08p-4, -0x1.2861d6p-5,
    },
    .{
        0x1.2c296cp-10, 0x1.f9d57p-1,   0x1.b8aab8p-5, -0x1.30b6a4p-2,
        0x1.73af8ap-3,  -0x1.a1564ep-5, 0x1.6908p-8,
    },
    .{
        -0x1.171628p-6, 0x1.1382e6p+0,  -0x1.f8ac18p-4, -0x1.c6795cp-4,
        0x1.2d2dfep-4,  -0x1.38f626p-6, 0x1.f9d716p-10,
    },
    .{
        -0x1.29739cp-4, 0x1.43decap+0,  -0x1.91761cp-2,  0x1.82261ap-4,
        -0x1.eee226p-7, 0x1.746546p-10, -0x1.f02e6ap-15,
    },
    .{
        0x1.6722aep-3,  0x1.d6b2e4p-1,  -0x1.7dbe34p-3,  0x1.c10cbep-6,
        -0x1.502cc2p-9, 0x1.1e2eb2p-13, -0x1.a56434p-19,
    },
    .{
        0x1.6f7ecp-1,    0x1.05fe1p-1,  -0x1.c290eap-5,  0x1.12d0bep-8,
        -0x1.a5b282p-13, 0x1.6d57cp-18, -0x1.108528p-24,
    },
    .{
        0x1.5ed2bep+0,   0x1.0d3bp-2,    -0x1.d55ab8p-7,  0x1.20902cp-11,
        -0x1.bd0e6ap-17, 0x1.82e9d2p-23, -0x1.214f5cp-30,
    },
    .{
        0x1.06cee4p+1,   0x1.0f17c2p-3,  -0x1.da366ep-9,  0x1.24229ap-14,
        -0x1.c327d6p-21, 0x1.889372p-28, -0x1.25be4cp-36,
    },
    .{
        0x1.5f39aep+1,   0x1.0f8076p-4,  -0x1.db39e6p-11, 0x1.24d726p-17,
        -0x1.c44a96p-25, 0x1.899114p-33, -0x1.26784ap-42,
    },
    .{
        0x1.b7e21p+1,    0x1.0f9492p-5,  -0x1.db65a6p-13, 0x1.24f0f8p-20,
        -0x1.c46c0ep-29, 0x1.89a692p-38, -0x1.2681d4p-48,
    },
    .{
        0x1.084d12p+2,  0x1.0f923p-6,   -0x1.db5616p-15, 0x1.24df02p-23,
        -0x1.c4422p-33, 0x1.89751cp-43, -0x1.2653cp-54,
    },
    .{
        0x1.34afcep+2,   0x1.0f76f6p-7,  -0x1.daf3c8p-17, 0x1.24818cp-26,
        -0x1.c37d7ep-37, 0x1.889b7ep-48, -0x1.258ca6p-60,
    },
    .{
        0x1.610b6ep+2,  0x1.0f7ad8p-8, -0x1.db022p-19,  0x1.248fa8p-29,
        -0x1.c39c7p-41, 0x1.88bf6p-53, -0x1.25aeeep-66,
    },
    .{
        0x1.8d6a98p+2,   0x1.0f6f2ep-9,  -0x1.dad874p-21, 0x1.24689p-32,
        -0x1.c34ab6p-45, 0x1.886532p-58, -0x1.255d5ep-72,
    },
    .{
        0x1.b9c986p+2,   0x1.0f651cp-10, -0x1.dab52cp-23, 0x1.2447e4p-35,
        -0x1.c3079ep-49, 0x1.881ca6p-63, -0x1.251c46p-78,
    },
    .{
        0x1.e6269ap+2,   0x1.0f62aep-11, -0x1.daac6ep-25, 0x1.243f9ep-38,
        -0x1.c2f63cp-53, 0x1.88096cp-68, -0x1.250ac6p-84,
    },
    .{
        0x1.09425p+3,    0x1.0f5c6ap-12, -0x1.da96ep-27,  0x1.242c12p-41,
        -0x1.c2ceccp-57, 0x1.87df72p-73, -0x1.24e5eep-90,
    },
    .{
        0x1.1f707ep+3,   0x1.0f5d12p-13, -0x1.da990ap-29, 0x1.242df2p-44,
        -0x1.c2d272p-61, 0x1.87e32ep-78, -0x1.24e916p-96,
    },
};

pub fn asinhh(value: f16) f16 {
    @setFloatMode(.strict);
    const input_bits: u16 = @bitCast(value);
    const magnitude = input_bits & magnitude_mask;
    if (magnitude > infinity_bits) return @bitCast(input_bits | quiet_nan_bit);
    if (magnitude == infinity_bits) return value;

    const value_f: f32 = @floatCast(value);
    const value_f_bits: u32 = @bitCast(value_f);
    const absolute_bits = value_f_bits & 0x7fff_ffff;
    const exponent: i32 = @as(i32, @intCast(absolute_bits >> 23)) - 124;
    const sign = (input_bits & sign_mask) != 0;
    var t = value_f;

    if (exponent <= 0) {
        if (absolute_bits <= 0x3d38_a000) {
            // The exact result rounds to the input encoding in this region.
            return value;
        }
        const t_squared = t * t;
        const p = coefficients[0];
        const c5 = fma_impl.fma(p[3], t_squared, p[2]);
        var c1 = fma_impl.fma(p[1], t_squared, p[0]);
        c1 = fma_impl.fma(c5, t_squared * t_squared, c1);
        return @floatCast(t * c1);
    }

    if (sign) t = -t;
    const t_squared = t * t;
    const p = coefficients[@intCast(exponent)];
    var c4 = fma_impl.fma(p[5], t, p[4]);
    const c2 = fma_impl.fma(p[3], t, p[2]);
    var c0 = fma_impl.fma(p[1], t, p[0]);
    c4 = fma_impl.fma(p[6], t_squared, c4);
    c0 = fma_impl.fma(c2, t_squared, c0);
    var result = fma_impl.fma(c4, t_squared * t_squared, c0);

    switch (exponent) {
        7 => {
            if (absolute_bits == 0x4193_6000) result = 0x1.cdbffep+1;
        },
        10 => {
            if (absolute_bits == 0x436c_c000) result = 0x1.8a4002p+2;
        },
        11 => {
            if (absolute_bits == 0x43f1_0000) result = 0x1.b7bffap+2;
        },
        12 => {
            if (absolute_bits == 0x4458_8000) result = 0x1.dd4004p+2;
        },
        else => {},
    }
    if (sign) result = -result;
    return @floatCast(result);
}
