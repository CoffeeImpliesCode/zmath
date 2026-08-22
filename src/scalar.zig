const sign_mask_h: u16 = 0x8000;
const sign_mask_f: u32 = 0x8000_0000;
const sign_mask_d: u64 = 0x8000_0000_0000_0000;

const magnitude_mask_h: u16 = 0x7fff;
const magnitude_mask_f: u32 = 0x7fff_ffff;
const magnitude_mask_d: u64 = 0x7fff_ffff_ffff_ffff;

const infinity_bits_h: u16 = 0x7c00;
const infinity_bits_f: u32 = 0x7f80_0000;
const infinity_bits_d: u64 = 0x7ff0_0000_0000_0000;

// The binary16 exponential kernel and tables are ported from CORE-MATH
// src/binary16/exp/expf16.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa (MIT license).
const exph_data = @import("scalar_exph_data.zig");

pub fn exph(value: f16) f16 {
    const bits: u16 = @bitCast(value);
    const magnitude = bits & magnitude_mask_h;

    if (magnitude > infinity_bits_h) {
        const quiet_nan_bit_h: u16 = 0x0200;
        return @bitCast(bits | quiet_nan_bit_h);
    }
    if (magnitude == infinity_bits_h) {
        if ((bits & sign_mask_h) == 0) return value;
        return @bitCast(@as(u16, 0));
    }
    if (magnitude == 0) return @bitCast(@as(u16, 0x3c00));

    // Decompose x into x1 + x2, then use exp(x) = exp(x1) * exp(x2).
    const table1_index = bits >> 5;
    const table2_index = ((bits >> 10) << 5) | (bits & 0x001f);
    const exp_x1: f32 = @bitCast(exph_data.t1[table1_index]);
    const exp_x2: f32 = @bitCast(exph_data.t2[table2_index]);
    return @floatCast(exp_x1 * exp_x2);
}

pub fn absh(value: f16) f16 {
    const bits: u16 = @bitCast(value);
    return @bitCast(bits & magnitude_mask_h);
}

pub fn absf(value: f32) f32 {
    const bits: u32 = @bitCast(value);
    return @bitCast(bits & magnitude_mask_f);
}

pub fn absd(value: f64) f64 {
    const bits: u64 = @bitCast(value);
    return @bitCast(bits & magnitude_mask_d);
}

pub fn ulpDistanceh(a: f16, b: f16) u16 {
    const a_bits: u16 = @bitCast(a);
    const b_bits: u16 = @bitCast(b);

    if (isNanH(a_bits) or isNanH(b_bits)) return ~@as(u16, 0);

    const a_ordered = orderedBitsH(a_bits);
    const b_ordered = orderedBitsH(b_bits);
    if (a_ordered >= b_ordered) return a_ordered - b_ordered;
    return b_ordered - a_ordered;
}

pub fn ulpDistancef(a: f32, b: f32) u32 {
    const a_bits: u32 = @bitCast(a);
    const b_bits: u32 = @bitCast(b);

    if (isNanF(a_bits) or isNanF(b_bits)) return ~@as(u32, 0);

    const a_ordered = orderedBitsF(a_bits);
    const b_ordered = orderedBitsF(b_bits);
    if (a_ordered >= b_ordered) return a_ordered - b_ordered;
    return b_ordered - a_ordered;
}

pub fn ulpDistanced(a: f64, b: f64) u64 {
    const a_bits: u64 = @bitCast(a);
    const b_bits: u64 = @bitCast(b);

    if (isNanD(a_bits) or isNanD(b_bits)) return ~@as(u64, 0);

    const a_ordered = orderedBitsD(a_bits);
    const b_ordered = orderedBitsD(b_bits);
    if (a_ordered >= b_ordered) return a_ordered - b_ordered;
    return b_ordered - a_ordered;
}

fn orderedBitsH(bits: u16) u16 {
    if ((bits & sign_mask_h) != 0) return (~bits) + 1;
    return bits | sign_mask_h;
}

fn orderedBitsF(bits: u32) u32 {
    if ((bits & sign_mask_f) != 0) return (~bits) + 1;
    return bits | sign_mask_f;
}

fn orderedBitsD(bits: u64) u64 {
    if ((bits & sign_mask_d) != 0) return (~bits) + 1;
    return bits | sign_mask_d;
}

fn isNanH(bits: u16) bool {
    return (bits & magnitude_mask_h) > infinity_bits_h;
}

fn isNanF(bits: u32) bool {
    return (bits & magnitude_mask_f) > infinity_bits_f;
}

fn isNanD(bits: u64) bool {
    return (bits & magnitude_mask_d) > infinity_bits_d;
}

// Correctly rounded natural exponential for binary32 values.
//
// Copyright (c) 2023-2025 Alexei Sibidanov.
//
// This code is derived from CORE-MATH at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa and is used under the MIT License:
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
const expf_quiet_nan_bit: u32 = 0x0040_0000;

const expf_fast_coefficients = [4]f64{
    1.0,
    0x1.62e42fef4c4e7p-1,
    0x1.ebfd1b232f475p-3,
    0x1.c6b19384ecd93p-5,
};

const expf_accurate_coefficients = [6]f64{
    0x1.62e42fefa39efp-1,
    0x1.ebfbdff82c58fp-3,
    0x1.c6b08d702e0edp-5,
    0x1.3b2ab6fb92e5ep-7,
    0x1.5d886e6d54203p-10,
    0x1.430976b8ce6efp-13,
};

const expf_scale_bits = [64]u64{
    0x3ff0000000000000, 0x3ff02c9a3e778061, 0x3ff059b0d3158574,
    0x3ff0874518759bc8, 0x3ff0b5586cf9890f, 0x3ff0e3ec32d3d1a2,
    0x3ff11301d0125b51, 0x3ff1429aaea92de0, 0x3ff172b83c7d517b,
    0x3ff1a35beb6fcb75, 0x3ff1d4873168b9aa, 0x3ff2063b88628cd6,
    0x3ff2387a6e756238, 0x3ff26b4565e27cdd, 0x3ff29e9df51fdee1,
    0x3ff2d285a6e4030b, 0x3ff306fe0a31b715, 0x3ff33c08b26416ff,
    0x3ff371a7373aa9cb, 0x3ff3a7db34e59ff7, 0x3ff3dea64c123422,
    0x3ff4160a21f72e2a, 0x3ff44e086061892d, 0x3ff486a2b5c13cd0,
    0x3ff4bfdad5362a27, 0x3ff4f9b2769d2ca7, 0x3ff5342b569d4f82,
    0x3ff56f4736b527da, 0x3ff5ab07dd485429, 0x3ff5e76f15ad2148,
    0x3ff6247eb03a5585, 0x3ff6623882552225, 0x3ff6a09e667f3bcd,
    0x3ff6dfb23c651a2f, 0x3ff71f75e8ec5f74, 0x3ff75feb564267c9,
    0x3ff7a11473eb0187, 0x3ff7e2f336cf4e62, 0x3ff82589994cce13,
    0x3ff868d99b4492ed, 0x3ff8ace5422aa0db, 0x3ff8f1ae99157736,
    0x3ff93737b0cdc5e5, 0x3ff97d829fde4e50, 0x3ff9c49182a3f090,
    0x3ffa0c667b5de565, 0x3ffa5503b23e255d, 0x3ffa9e6b5579fdbf,
    0x3ffae89f995ad3ad, 0x3ffb33a2b84f15fb, 0x3ffb7f76f2fb5e47,
    0x3ffbcc1e904bc1d2, 0x3ffc199bdd85529c, 0x3ffc67f12e57d14b,
    0x3ffcb720dcef9069, 0x3ffd072d4a07897c, 0x3ffd5818dcfba487,
    0x3ffda9e603db3285, 0x3ffdfc97337b9b5f, 0x3ffe502ee78b3ff6,
    0x3ffea4afa2a490da, 0x3ffefa1bee615a27, 0x3fff50765b6e4540,
    0x3fffa7c1819e90d8,
};

pub fn expf(value: f32) f32 {
    const bits: u32 = @bitCast(value);
    const magnitude = bits & magnitude_mask_f;

    if (magnitude == 0) return 1.0;
    if (magnitude > infinity_bits_f) {
        return @bitCast(bits | expf_quiet_nan_bit);
    }
    if (magnitude == infinity_bits_f) {
        if ((bits & sign_mask_f) == 0) return @bitCast(infinity_bits_f);
        return 0.0;
    }

    return expfFinite(value, bits);
}

fn expfFinite(value: f32, bits: u32) f32 {
    @setFloatMode(.strict);

    const inverse_ln2: f64 = 0x1.71547652b82fep+0;
    const rounding_constant: f64 = 0x1.8p46;
    const value_wide: f64 = value;
    const scaled = inverse_ln2 * value_wide;
    const rounded = scaled + rounding_constant;
    const rounded_bits: u64 = @bitCast(rounded);
    const absolute_bits = bits << 1;

    if (absolute_bits > 0x8562_e42e or absolute_bits < 0x6f93_813e) {
        if (absolute_bits < 0x6f93_813e) {
            return @floatCast(1.0 + value_wide * (1.0 + value_wide * 0.5));
        }
        if (bits > 0xc2ce_8ec0) {
            var result = 0x1p-149 +
                (value_wide + 0x1.9d1d9fccf477p+6) * 0x1.71547652b82edp-150;
            if (result < 0x1p-151) result = 0x1p-151;
            return @floatCast(result);
        }
        if ((bits & sign_mask_f) == 0 and bits > 0x42b1_7217) {
            return @bitCast(infinity_bits_f);
        }
    }

    const correction = rounding_constant - rounded;
    var reduced = scaled + correction;
    const table_index: u6 = @truncate(rounded_bits);
    const exponent_bits = (rounded_bits >> 6) << 52;
    const scale_bits = expf_scale_bits[table_index] +% exponent_bits;
    const scale: f64 = @bitCast(scale_bits);
    var reduced_squared = reduced * reduced;
    var result =
        ((expf_fast_coefficients[0] + reduced * expf_fast_coefficients[1]) +
            reduced_squared *
                (expf_fast_coefficients[2] + reduced * expf_fast_coefficients[3])) *
        scale;
    var upper: f32 = @floatCast(result);
    const lower: f32 = @floatCast(result - result * 1.45e-10);

    if (upper != lower) {
        const inverse_ln2_high: f64 = 0x1.7154765p+0;
        const inverse_ln2_low: f64 = 0x1.5c17f0bbbe88p-31;
        reduced =
            (inverse_ln2_high * value_wide + correction) +
            inverse_ln2_low * value_wide;
        reduced_squared = reduced * reduced;
        const scaled_remainder = scale * reduced;
        result = scale +
            scaled_remainder *
                ((expf_accurate_coefficients[0] +
                    reduced * expf_accurate_coefficients[1]) +
                    reduced_squared *
                        ((expf_accurate_coefficients[2] +
                            reduced * expf_accurate_coefficients[3]) +
                            reduced_squared *
                                (expf_accurate_coefficients[4] +
                                    reduced * expf_accurate_coefficients[5])));
        upper = @floatCast(result);
    }

    return upper;
}
// Correctly rounded binary64 exponential, ported from CORE-MATH revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa.
// Copyright (c) 2022-2025 Alexei Sibidanov.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

const ExpdDoubleDouble = struct {
    high: f64,
    low: f64,
};

const expd_database: [51]f64 = .{
    0x1.fffffffffffffp-53,  0x1.ba07d73250de7p-14,  0x1.6a4d1af9cc989p-8,   0x1.5a75293a5dcdap-6,
    0x1.42ea46949b3c7p-5,   0x1.7c8bb0cf5d16p-5,    0x1.0948d39a41695p-3,   0x1.a065fefae814fp-3,
    0x1.f6e4c3ced7c72p-3,   0x1.1a0408712e00ap-2,   0x1.bcab27d05abdep-2,   0x1.005ae04256babp-1,
    0x1.273c188aa7b14p+2,   0x1.83d4bcdebb3f4p+2,   0x1.08f51434652c3p+4,   0x1.1d5c2daebe367p+4,
    0x1.c44ce0d716a1ap+4,   0x1.e07e71bfcf06fp+5,   0x1.f7216c4b435c9p+5,   0x1.54cd1fea7663ap+7,
    0x1.d6479eba7c971p+8,   -0x1.664716b68a409p-14, -0x1.a2fefefd580dfp-13, -0x1.ce3f638d0c742p-12,
    -0x1.ceff32831e2c2p-12, -0x1.33accae78b371p-11, -0x1.d792b60084f92p-11, -0x1.7fb235d76cce7p-8,
    -0x1.1ff9b8e8b38bep-7,  -0x1.54511e930898cp-7,  -0x1.5c5ed0ec83666p-6,  -0x1.8c56ff5326197p-6,
    -0x1.a4187f2ca71f9p-6,  -0x1.a8f783d749a8fp-4,  -0x1.bd44fdaed819fp-4,  -0x1.daf693d64fadap-4,
    -0x1.290ea09e36479p-3,  -0x1.8aeb636f3ce35p-3,  -0x1.d3f3799439415p-3,  -0x1.ea16274b0109bp-3,
    -0x1.22e24fa3d5cf9p-1,  -0x1.85068c07fbbf6p-1,  -0x1.bdc7955d1482cp-1,  -0x1.2a9cad9998262p+0,
    -0x1.cc37ef7de7501p+0,  -0x1.02393d5976769p+1,  -0x1.65061daf79a78p+1,  -0x1.e8bdbfcd9144ep+3,
    -0x1.8f80e06f3a04cp+4,  -0x1.59f038076039cp+6,  -0x1.981587ad4542fp+7,
};

// Double-double approximations of 2^(i / 2^6), for 0 <= i < 64.
const expd_table_0: [64]ExpdDoubleDouble = .{ .{ .high = 0x1p+0, .low = 0x0p+0 }, .{ .high = 0x1.02c9a3e778061p+0, .low = -0x1.19083535b085ep-56 }, .{ .high = 0x1.059b0d3158574p+0, .low = 0x1.d73e2a475b466p-55 }, .{ .high = 0x1.0874518759bc8p+0, .low = 0x1.186be4bb285p-57 }, .{ .high = 0x1.0b5586cf9890fp+0, .low = 0x1.8a62e4adc610ap-54 }, .{ .high = 0x1.0e3ec32d3d1a2p+0, .low = 0x1.03a1727c57b52p-59 }, .{ .high = 0x1.11301d0125b51p+0, .low = -0x1.6c51039449b3ap-54 }, .{ .high = 0x1.1429aaea92dep+0, .low = -0x1.32fbf9af1369ep-54 }, .{ .high = 0x1.172b83c7d517bp+0, .low = -0x1.19041b9d78a76p-55 }, .{ .high = 0x1.1a35beb6fcb75p+0, .low = 0x1.e5b4c7b4968e4p-55 }, .{ .high = 0x1.1d4873168b9aap+0, .low = 0x1.e016e00a2643cp-54 }, .{ .high = 0x1.2063b88628cd6p+0, .low = 0x1.dc775814a8494p-55 }, .{ .high = 0x1.2387a6e756238p+0, .low = 0x1.9b07eb6c70572p-54 }, .{ .high = 0x1.26b4565e27cddp+0, .low = 0x1.2bd339940e9dap-55 }, .{ .high = 0x1.29e9df51fdee1p+0, .low = 0x1.612e8afad1256p-55 }, .{ .high = 0x1.2d285a6e4030bp+0, .low = 0x1.0024754db41d4p-54 }, .{ .high = 0x1.306fe0a31b715p+0, .low = 0x1.6f46ad23182e4p-55 }, .{ .high = 0x1.33c08b26416ffp+0, .low = 0x1.32721843659a6p-54 }, .{ .high = 0x1.371a7373aa9cbp+0, .low = -0x1.63aeabf42eae2p-54 }, .{ .high = 0x1.3a7db34e59ff7p+0, .low = -0x1.5e436d661f5e2p-56 }, .{ .high = 0x1.3dea64c123422p+0, .low = 0x1.ada0911f09ebcp-55 }, .{ .high = 0x1.4160a21f72e2ap+0, .low = -0x1.ef3691c309278p-58 }, .{ .high = 0x1.44e086061892dp+0, .low = 0x1.89b7a04ef80dp-59 }, .{ .high = 0x1.486a2b5c13cdp+0, .low = 0x1.3c1a3b69062fp-56 }, .{ .high = 0x1.4bfdad5362a27p+0, .low = 0x1.d4397afec42e2p-56 }, .{ .high = 0x1.4f9b2769d2ca7p+0, .low = -0x1.4b309d25957e4p-54 }, .{ .high = 0x1.5342b569d4f82p+0, .low = -0x1.07abe1db13cacp-55 }, .{ .high = 0x1.56f4736b527dap+0, .low = 0x1.9bb2c011d93acp-54 }, .{ .high = 0x1.5ab07dd485429p+0, .low = 0x1.6324c054647acp-54 }, .{ .high = 0x1.5e76f15ad2148p+0, .low = 0x1.ba6f93080e65ep-54 }, .{ .high = 0x1.6247eb03a5585p+0, .low = -0x1.383c17e40b496p-54 }, .{ .high = 0x1.6623882552225p+0, .low = -0x1.bb60987591c34p-54 }, .{ .high = 0x1.6a09e667f3bcdp+0, .low = -0x1.bdd3413b26456p-54 }, .{ .high = 0x1.6dfb23c651a2fp+0, .low = -0x1.bbe3a683c88aap-57 }, .{ .high = 0x1.71f75e8ec5f74p+0, .low = -0x1.16e4786887a9ap-55 }, .{ .high = 0x1.75feb564267c9p+0, .low = -0x1.0245957316dd4p-54 }, .{ .high = 0x1.7a11473eb0187p+0, .low = -0x1.41577ee04993p-55 }, .{ .high = 0x1.7e2f336cf4e62p+0, .low = 0x1.05d02ba15797ep-56 }, .{ .high = 0x1.82589994cce13p+0, .low = -0x1.d4c1dd41532d8p-54 }, .{ .high = 0x1.868d99b4492edp+0, .low = -0x1.fc6f89bd4f6bap-54 }, .{ .high = 0x1.8ace5422aa0dbp+0, .low = 0x1.6e9f156864b26p-54 }, .{ .high = 0x1.8f1ae99157736p+0, .low = 0x1.5cc13a2e3976cp-55 }, .{ .high = 0x1.93737b0cdc5e5p+0, .low = -0x1.75fc781b57ebcp-57 }, .{ .high = 0x1.97d829fde4e5p+0, .low = -0x1.d185b7c1b85dp-54 }, .{ .high = 0x1.9c49182a3f09p+0, .low = 0x1.c7c46b071f2bep-56 }, .{ .high = 0x1.a0c667b5de565p+0, .low = -0x1.359495d1cd532p-54 }, .{ .high = 0x1.a5503b23e255dp+0, .low = -0x1.d2f6edb8d41e2p-54 }, .{ .high = 0x1.a9e6b5579fdbfp+0, .low = 0x1.0fac90ef7fd32p-54 }, .{ .high = 0x1.ae89f995ad3adp+0, .low = 0x1.7a1cd345dcc82p-54 }, .{ .high = 0x1.b33a2b84f15fbp+0, .low = -0x1.2805e3084d708p-57 }, .{ .high = 0x1.b7f76f2fb5e47p+0, .low = -0x1.5584f7e54ac3ap-56 }, .{ .high = 0x1.bcc1e904bc1d2p+0, .low = 0x1.23dd07a2d9e84p-55 }, .{ .high = 0x1.c199bdd85529cp+0, .low = 0x1.11065895048dep-55 }, .{ .high = 0x1.c67f12e57d14bp+0, .low = 0x1.2884dff483cacp-54 }, .{ .high = 0x1.cb720dcef9069p+0, .low = 0x1.503cbd1e949dcp-56 }, .{ .high = 0x1.d072d4a07897cp+0, .low = -0x1.cbc3743797a9cp-54 }, .{ .high = 0x1.d5818dcfba487p+0, .low = 0x1.2ed02d75b3706p-55 }, .{ .high = 0x1.da9e603db3285p+0, .low = 0x1.c2300696db532p-54 }, .{ .high = 0x1.dfc97337b9b5fp+0, .low = -0x1.1a5cd4f184b5cp-54 }, .{ .high = 0x1.e502ee78b3ff6p+0, .low = 0x1.39e8980a9cc9p-55 }, .{ .high = 0x1.ea4afa2a490dap+0, .low = -0x1.e9c23179c2894p-54 }, .{ .high = 0x1.efa1bee615a27p+0, .low = 0x1.dc7f486a4b6bp-54 }, .{ .high = 0x1.f50765b6e454p+0, .low = 0x1.9d3e12dd8a18ap-54 }, .{ .high = 0x1.fa7c1819e90d8p+0, .low = 0x1.74853f3a5931ep-55 } };

// Double-double approximations of 2^(i / 2^12), for 0 <= i < 64.
const expd_table_1: [64]ExpdDoubleDouble = .{
    .{ .high = 0x1p+0, .low = 0x0p+0 },                               .{ .high = 0x1.000b175effdc7p+0, .low = 0x1.ae8e38c59c72ap-54 },
    .{ .high = 0x1.00162f3904052p+0, .low = -0x1.7b5d0d58ea8f4p-58 }, .{ .high = 0x1.0021478e11ce6p+0, .low = 0x1.4115cb6b16a8ep-54 },
    .{ .high = 0x1.002c605e2e8cfp+0, .low = -0x1.d7c96f201bb2ep-55 }, .{ .high = 0x1.003779a95f959p+0, .low = 0x1.84711d4c35eap-54 },
    .{ .high = 0x1.0042936faa3d8p+0, .low = -0x1.0484245243778p-55 }, .{ .high = 0x1.004dadb113dap+0, .low = -0x1.4b237da2025fap-54 },
    .{ .high = 0x1.0058c86da1c0ap+0, .low = -0x1.5e00e62d6b30ep-56 }, .{ .high = 0x1.0063e3a559473p+0, .low = 0x1.a1d6cedbb948p-54 },
    .{ .high = 0x1.006eff583fc3dp+0, .low = -0x1.4acf197a00142p-54 }, .{ .high = 0x1.007a1b865a8cap+0, .low = -0x1.eaf2ea42391a6p-57 },
    .{ .high = 0x1.0085382faef83p+0, .low = 0x1.da93f90835f76p-56 },  .{ .high = 0x1.00905554425d4p+0, .low = -0x1.6a79084ab093cp-55 },
    .{ .high = 0x1.009b72f41a12bp+0, .low = 0x1.86364f8fbe8f8p-54 },  .{ .high = 0x1.00a6910f3b6fdp+0, .low = -0x1.82e8e14e3110ep-55 },
    .{ .high = 0x1.00b1afa5abcbfp+0, .low = -0x1.4f6b2a7609f72p-55 }, .{ .high = 0x1.00bcceb7707ecp+0, .low = -0x1.e1a258ea8f71ap-56 },
    .{ .high = 0x1.00c7ee448ee02p+0, .low = 0x1.4362ca5bc26f2p-56 },  .{ .high = 0x1.00d30e4d0c483p+0, .low = 0x1.095a56c919d02p-54 },
    .{ .high = 0x1.00de2ed0ee0f5p+0, .low = -0x1.406ac4e81a646p-57 }, .{ .high = 0x1.00e94fd0398ep+0, .low = 0x1.b5a6902767e08p-54 },
    .{ .high = 0x1.00f4714af41d3p+0, .low = -0x1.91b206085932p-54 },  .{ .high = 0x1.00ff93412315cp+0, .low = 0x1.427068ab22306p-55 },
    .{ .high = 0x1.010ab5b2cbd11p+0, .low = 0x1.c1d0660524e08p-54 },  .{ .high = 0x1.0115d89ff3a8bp+0, .low = -0x1.e7bdfb3204be8p-54 },
    .{ .high = 0x1.0120fc089ff63p+0, .low = 0x1.843aa8b9cbbc6p-55 },  .{ .high = 0x1.012c1fecd613bp+0, .low = -0x1.34104ee7edae8p-56 },
    .{ .high = 0x1.0137444c9b5b5p+0, .low = -0x1.2b6aeb6176892p-56 }, .{ .high = 0x1.01426927f5278p+0, .low = 0x1.a8cd33b8a1bb2p-56 },
    .{ .high = 0x1.014d8e7ee8d2fp+0, .low = 0x1.2edc08e5da99ap-56 },  .{ .high = 0x1.0158b4517bb88p+0, .low = 0x1.57ba2dc7e0c72p-55 },
    .{ .high = 0x1.0163da9fb3335p+0, .low = 0x1.b61299ab8cdb8p-54 },  .{ .high = 0x1.016f0169949edp+0, .low = -0x1.90565902c5f44p-54 },
    .{ .high = 0x1.017a28af25567p+0, .low = 0x1.70fc41c5c2d54p-55 },  .{ .high = 0x1.018550706ab62p+0, .low = 0x1.4b9a6e145d76cp-54 },
    .{ .high = 0x1.019078ad6a19fp+0, .low = -0x1.008eff5142bfap-56 }, .{ .high = 0x1.019ba16628de2p+0, .low = -0x1.77669f033c7dep-54 },
    .{ .high = 0x1.01a6ca9aac5f3p+0, .low = -0x1.09bb78eeead0ap-54 }, .{ .high = 0x1.01b1f44af9f9ep+0, .low = 0x1.371231477ece6p-54 },
    .{ .high = 0x1.01bd1e77170b4p+0, .low = 0x1.5e7626621eb5ap-56 },  .{ .high = 0x1.01c8491f08f08p+0, .low = -0x1.bc72b100828a4p-54 },
    .{ .high = 0x1.01d37442d507p+0, .low = -0x1.ce39cbbab8bbep-57 },  .{ .high = 0x1.01de9fe280ac8p+0, .low = 0x1.16996709da2e2p-55 },
    .{ .high = 0x1.01e9cbfe113efp+0, .low = -0x1.c11f5239bf536p-55 }, .{ .high = 0x1.01f4f8958c1c6p+0, .low = 0x1.e1d4eb5edc6b4p-55 },
    .{ .high = 0x1.020025a8f6a35p+0, .low = -0x1.afb99946ee3fp-54 },  .{ .high = 0x1.020b533856324p+0, .low = -0x1.8f06d8a148a32p-54 },
    .{ .high = 0x1.02168143b0281p+0, .low = -0x1.2bf310fc54eb6p-55 }, .{ .high = 0x1.0221afcb09e3ep+0, .low = -0x1.c95a035eb4176p-54 },
    .{ .high = 0x1.022cdece68c4fp+0, .low = -0x1.491793e46834cp-54 }, .{ .high = 0x1.02380e4dd22adp+0, .low = -0x1.3e8d0d9c4909p-56 },
    .{ .high = 0x1.02433e494b755p+0, .low = -0x1.314aa16278aa4p-54 }, .{ .high = 0x1.024e6ec0da046p+0, .low = 0x1.48daf888e965p-55 },
    .{ .high = 0x1.02599fb483385p+0, .low = 0x1.56dc8046821f4p-55 },  .{ .high = 0x1.0264d1244c719p+0, .low = 0x1.45b42356b9d46p-54 },
    .{ .high = 0x1.027003103b10ep+0, .low = -0x1.082ef51b61d7ep-56 }, .{ .high = 0x1.027b357854772p+0, .low = 0x1.2106ed0920a34p-56 },
    .{ .high = 0x1.0286685c9e059p+0, .low = -0x1.fd4cf26ea5d0ep-54 }, .{ .high = 0x1.02919bbd1d1d8p+0, .low = -0x1.09f8775e78084p-54 },
    .{ .high = 0x1.029ccf99d720ap+0, .low = 0x1.64cbba902ca28p-58 },  .{ .high = 0x1.02a803f2d170dp+0, .low = 0x1.4383ef231d206p-54 },
    .{ .high = 0x1.02b338c811703p+0, .low = 0x1.4a47a505b3a46p-54 },  .{ .high = 0x1.02be6e199c811p+0, .low = 0x1.e47120223468p-54 },
};

const expd_accurate_coefficients: [7]ExpdDoubleDouble = .{
    .{ .high = 0x1p+0, .low = 0x0p+0 },
    .{ .high = 0x1p-1, .low = 0x1.712f72ecec2cfp-99 },
    .{ .high = 0x1.5555555555555p-3, .low = 0x1.5555555554d07p-57 },
    .{ .high = 0x1.5555555555555p-5, .low = 0x1.55194d28275dap-59 },
    .{ .high = 0x1.1111111111111p-7, .low = 0x1.12faa0e1c0f7bp-63 },
    .{ .high = 0x1.6c16c16da6973p-10, .low = -0x1.4ba45ab25d2a3p-64 },
    .{ .high = 0x1.a01a019eb7f31p-13, .low = -0x1.9091d845ecd36p-67 },
};

fn expdRoundEvenFinite(value: f64) i64 {
    const truncated: i64 = @intFromFloat(value);
    const truncated_float: f64 = @floatFromInt(truncated);
    const fraction = value - truncated_float;
    const odd = (truncated & 1) != 0;

    if (fraction > 0x1p-1 or (fraction == 0x1p-1 and odd)) return truncated + 1;
    if (fraction < -0x1p-1 or (fraction == -0x1p-1 and odd)) return truncated - 1;
    return truncated;
}

fn expdFastTwoSum(x: f64, y: f64) ExpdDoubleDouble {
    const sum = x + y;
    const rounded_x = sum - x;
    return .{ .high = sum, .low = y - rounded_x };
}

fn expdFastSum(x_high: f64, x_low: f64, y_high: f64, y_low: f64) ExpdDoubleDouble {
    const sum = expdFastTwoSum(x_high, y_high);
    return .{ .high = sum.high, .low = (x_low + y_low) + sum.low };
}

fn expdMultiplyDoubleDouble(
    x_high: f64,
    x_low: f64,
    c_high: f64,
    c_low: f64,
) ExpdDoubleDouble {
    const high = c_high * x_high;
    const low = (c_high * x_low + c_low * x_high) + @mulAdd(f64, c_high, x_high, -high);
    return .{ .high = high, .low = low };
}

fn expdFastFmaDoubleWord(
    a_high: f64,
    a_low: f64,
    b_high: f64,
    b_low: f64,
    c_high: f64,
    c_low: f64,
) ExpdDoubleDouble {
    const high = @mulAdd(f64, a_high, b_high, c_high);
    const correction = c_high - high;
    const product_error = @mulAdd(f64, a_high, b_high, correction);
    const combined_error = product_error + c_low;
    const cross = a_high * b_low + combined_error;
    const low = a_low * b_high + cross;
    return .{ .high = high, .low = low };
}

fn expdEvaluateDoubleDoublePolynomial(
    x_high: f64,
    x_low: f64,
    coefficients: []const ExpdDoubleDouble,
) ExpdDoubleDouble {
    var index = coefficients.len - 1;
    var result = coefficients[index];
    while (index != 0) {
        index -= 1;
        const coefficient = coefficients[index];
        result = expdFastFmaDoubleWord(
            x_high,
            x_low,
            result.high,
            result.low,
            coefficient.high,
            coefficient.low,
        );
    }
    return result;
}

fn expdScaleNormal(value: f64, exponent: i64) f64 {
    const bits: u64 = @bitCast(value);
    const exponent_bits = @as(u64, @bitCast(exponent)) << 52;
    return @bitCast(bits +% exponent_bits);
}

fn expdToDenormal(value: f64) f64 {
    const bits: u64 = @bitCast(value);
    return @bitCast(bits & 0x000f_ffff_ffff_ffff);
}

fn expdDatabaseLookup(value: f64, fallback: f64) f64 {
    const value_bits: u64 = @bitCast(value);
    var lower: i32 = 0;
    var upper: i32 = @intCast(expd_database.len - 1);

    while (lower <= upper) {
        const middle = @divTrunc(lower + upper, 2);
        const middle_index: usize = @intCast(middle);
        const database_bits: u64 = @bitCast(expd_database[middle_index]);
        if (database_bits < value_bits) {
            lower = middle + 1;
        } else if (database_bits == value_bits) {
            const rounding_patterns = [2]u64{ 0x57f5fe2e5bde4075, 0x0000003c1f16b8ed };
            const sign_pattern: u64 = 333811522313371;
            const middle_shift: u6 = @intCast(middle);
            const delta_bits = (((sign_pattern >> middle_shift) & 1) << 63) | 0x3c90_0000_0000_0000;
            const delta: f64 = @bitCast(delta_bits);
            const pattern_index: usize = @intCast(middle >> 5);
            const pattern_shift: u6 = @intCast((middle << 1) & 63);
            const rounding_pattern = (rounding_patterns[pattern_index] >> pattern_shift) & 3;
            const fallback_bits: u64 = @bitCast(fallback);

            var offset: i64 = -1;
            while (offset <= 1) : (offset += 1) {
                const candidate_bits = if (offset < 0)
                    fallback_bits -% 1
                else if (offset > 0)
                    fallback_bits +% 1
                else
                    fallback_bits;
                if ((candidate_bits & 3) == rounding_pattern) {
                    const candidate: f64 = @bitCast(candidate_bits);
                    return candidate + delta;
                }
            }
            break;
        } else {
            upper = middle - 1;
        }
    }
    return fallback;
}

fn expdAccurate(value: f64) f64 {
    const value_bits: u64 = @bitCast(value);
    if (((value_bits >> 52) & 0x7ff) < 0x3c9) return 1.0 + value;

    const reduction_scale = 0x1.71547652b82fep+12;
    const reduction_integer = expdRoundEvenFinite(value * reduction_scale);
    const reduction: f64 = @floatFromInt(reduction_integer);
    const table_0_index: usize = @intCast((reduction_integer >> 6) & 0x3f);
    const table_1_index: usize = @intCast(reduction_integer & 0x3f);
    const exponent = reduction_integer >> 12;
    const table_0 = expd_table_0[table_0_index];
    const table_1 = expd_table_1[table_1_index];
    const table_product = expdMultiplyDoubleDouble(
        table_0.high,
        table_0.low,
        table_1.high,
        table_1.low,
    );

    const log2_high = 0x1.62e42ffp-13;
    const log2_low = 0x1.718432a1b0e26p-47;
    const log2_tail = 0x1.9ff0342542fc3p-102;
    const reduced = value - log2_high * reduction;
    var reduced_low = log2_low * reduction;
    const reduced_tail = log2_tail * reduction + @mulAdd(f64, log2_low, reduction, -reduced_low);
    const reduced_high = reduced + reduced_low;
    reduced_low = (reduced - reduced_high) + reduced_low + reduced_tail;

    var polynomial = expdEvaluateDoubleDoublePolynomial(
        reduced_high,
        reduced_low,
        &expd_accurate_coefficients,
    );
    polynomial = expdMultiplyDoubleDouble(
        reduced_high,
        reduced_low,
        polynomial.high,
        polynomial.low,
    );

    if (value_bits > 0xc086_232b_dd7a_bcd2) {
        const scale_bits = @as(u64, @intCast(1 - exponent)) << 52;
        const scale: f64 = @bitCast(scale_bits);
        const product = expdMultiplyDoubleDouble(
            polynomial.high,
            polynomial.low,
            table_product.high,
            table_product.low,
        );
        var result = expdFastSum(
            table_product.high,
            table_product.low,
            product.high,
            product.low,
        );
        const scaled_sum = expdFastTwoSum(scale, result.high);
        result.high = scaled_sum.high;
        result.low += scaled_sum.low;
        return expdToDenormal(result.high + result.low);
    }

    var result = ExpdDoubleDouble{ .high = 0.0, .low = 0.0 };
    if (table_product.high == 1.0) {
        const high_sum = expdFastTwoSum(table_product.high, polynomial.high);
        const low_sum = expdFastTwoSum(high_sum.low, polynomial.low);
        result = .{ .high = high_sum.high, .low = low_sum.high };

        var low_bits: u64 = @bitCast(result.low);
        if ((low_bits & 0x000f_ffff_ffff_ffff) == 0) {
            const error_bits: u64 = @bitCast(low_sum.low);
            if (((low_bits ^ error_bits) & sign_mask_d) == 0) {
                low_bits +%= 1;
            } else {
                low_bits -%= 1;
            }
            result.low = @bitCast(low_bits);
        }
    } else {
        const product = expdMultiplyDoubleDouble(
            polynomial.high,
            polynomial.low,
            table_product.high,
            table_product.low,
        );
        result = expdFastSum(
            table_product.high,
            table_product.low,
            product.high,
            product.low,
        );
    }

    result = expdFastTwoSum(result.high, result.low);
    const low_bits: u64 = @bitCast(result.low);
    const ambiguity = (low_bits +% 2) & 0x000f_ffff_ffff_ffff;
    if (ambiguity <= 2) result.high = expdDatabaseLookup(value, result.high);
    return expdScaleNormal(result.high, exponent);
}

pub fn expd(value: f64) f64 {
    const bits: u64 = @bitCast(value);
    const magnitude = bits & magnitude_mask_d;

    if (magnitude > infinity_bits_d) {
        return @bitCast(bits | 0x0008_0000_0000_0000);
    }
    if (magnitude == infinity_bits_d) {
        if ((bits & sign_mask_d) != 0) return 0.0;
        return value;
    }
    if (magnitude == 0) return 1.0;

    // exp(x) rounds to one in nearest-even mode for |x| <= 2^-54.
    if (magnitude <= 0x3c90_0000_0000_0000) return 1.0 + value;

    if (magnitude >= 0x4086_2e42_fefa_39f0) {
        if ((bits & sign_mask_d) == 0) return @bitCast(infinity_bits_d);
        if (magnitude >= 0x4087_4910_d52d_3052) return 0.0;
    }

    const reduction_scale = 0x1.71547652b82fep+12;
    const reduction_integer = expdRoundEvenFinite(value * reduction_scale);
    const reduction: f64 = @floatFromInt(reduction_integer);
    const table_0_index: usize = @intCast((reduction_integer >> 6) & 0x3f);
    const table_1_index: usize = @intCast(reduction_integer & 0x3f);
    const exponent = reduction_integer >> 12;
    const table_0 = expd_table_0[table_0_index];
    const table_1 = expd_table_1[table_1_index];
    const table_product = expdMultiplyDoubleDouble(
        table_0.high,
        table_0.low,
        table_1.high,
        table_1.low,
    );

    const log2_high = 0x1.62e42ffp-13;
    const log2_low = 0x1.718432a1b0e26p-47;
    const reduced = (value - log2_high * reduction) + log2_low * reduction;
    const reduced_squared = reduced * reduced;
    const coefficient_0 = 0x1p+0;
    const coefficient_1 = 0x1p-1;
    const coefficient_2 = 0x1.55555557e54ffp-3;
    const coefficient_3 = 0x1.55555553a12f4p-5;
    const polynomial = (coefficient_0 + reduced * coefficient_1) +
        reduced_squared * (coefficient_2 + reduced * coefficient_3);

    var high = table_product.high;
    const table_times_reduced = table_product.high * reduced;
    var low = table_product.low + table_times_reduced * polynomial;
    const error_bound = 1.64e-19;

    if (bits > 0xc086_232b_dd7a_bcd2) {
        const scale_bits = @as(u64, @intCast(1 - exponent)) << 52;
        const scale: f64 = @bitCast(scale_bits);
        const sum = expdFastTwoSum(scale, high);
        high = sum.high;
        low += sum.low;
        const upper = high + (low + error_bound);
        const lower = high + (low - error_bound);
        if (upper != lower) return expdAccurate(value);
        return expdToDenormal(lower);
    }

    const upper = high + (low + error_bound);
    const lower = high + (low - error_bound);
    if (upper != lower) return expdAccurate(value);
    return expdScaleNormal(lower, exponent);
}
