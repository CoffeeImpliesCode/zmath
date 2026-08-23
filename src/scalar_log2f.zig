// This file is part of the CORE-MATH project
// (https://core-math.gitlabpages.inria.fr/).
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
const sign_mask: u32 = 0x8000_0000;
const magnitude_mask: u32 = 0x7fff_ffff;
const infinity_bits: u32 = 0x7f80_0000;
const quiet_nan_bit: u32 = 0x0040_0000;
const canonical_nan_bits: u32 = 0x7fc0_0000;
const minus_infinity_bits: u32 = 0xff80_0000;

const ix = [_]f64{
    0x1.0000000000000p+0, 0x1.f81f820000000p-1, 0x1.f07c1f0000000p-1, 0x1.e9131ac000000p-1,
    0x1.e1e1e1e000000p-1, 0x1.dae6077000000p-1, 0x1.d41d41d000000p-1, 0x1.cd85689000000p-1,
    0x1.c71c71c000000p-1, 0x1.c0e0704000000p-1, 0x1.bacf915000000p-1, 0x1.b4e81b5000000p-1,
    0x1.af286bd000000p-1, 0x1.a98ef60000000p-1, 0x1.a41a41a000000p-1, 0x1.9ec8e95000000p-1,
    0x1.999999a000000p-1, 0x1.948b0fd000000p-1, 0x1.8f9c190000000p-1, 0x1.8acb90f000000p-1,
    0x1.8618618000000p-1, 0x1.8181818000000p-1, 0x1.7d05f41000000p-1, 0x1.78a4c81000000p-1,
    0x1.745d174000000p-1, 0x1.702e05c000000p-1, 0x1.6c16c17000000p-1, 0x1.6816817000000p-1,
    0x1.642c859000000p-1, 0x1.6058160000000p-1, 0x1.5c9882c000000p-1, 0x1.58ed231000000p-1,
    0x1.5555555000000p-1, 0x1.51d07eb000000p-1, 0x1.4e5e0a7000000p-1, 0x1.4afd6a0000000p-1,
    0x1.47ae148000000p-1, 0x1.446f865000000p-1, 0x1.4141414000000p-1, 0x1.3e22cbd000000p-1,
    0x1.3b13b14000000p-1, 0x1.3813814000000p-1, 0x1.3521cfb000000p-1, 0x1.323e34a000000p-1,
    0x1.2f684be000000p-1, 0x1.2c9fb4e000000p-1, 0x1.29e412a000000p-1, 0x1.27350b9000000p-1,
    0x1.2492492000000p-1, 0x1.21fb781000000p-1, 0x1.1f7047e000000p-1, 0x1.1cf06ae000000p-1,
    0x1.1a7b961000000p-1, 0x1.1811812000000p-1, 0x1.15b1e5f000000p-1, 0x1.135c811000000p-1,
    0x1.1111111000000p-1, 0x1.0ecf56c000000p-1, 0x1.0c97150000000p-1, 0x1.0a6810a000000p-1,
    0x1.0842108000000p-1, 0x1.0624dd3000000p-1, 0x1.0410410000000p-1, 0x1.0204081000000p-1,
    0x1.0000000000000p-1,
};

const lix = [_]f64{
    0x0.0p+0,              -0x1.6e7966ead8ac5p-6, -0x1.6bad38119a13ap-5, -0x1.0eb389ee9f560p-4,
    -0x1.663f6fc3a678dp-4, -0x1.bc8423d408321p-4, -0x1.08c588e79f8e9p-3, -0x1.32ae9e28fc362p-3,
    -0x1.5c01a3cde7f74p-3, -0x1.84c2bccf005a9p-3, -0x1.acf5e2c156d8fp-3, -0x1.d49ee4b90c47cp-3,
    -0x1.fbc16b67c143bp-3, -0x1.11307dc445fecp-2, -0x1.24407abf4dc03p-2, -0x1.37124cede831fp-2,
    -0x1.49a784a5bc715p-2, -0x1.5c01a3965cc38p-2, -0x1.6e221cc2bb868p-2, -0x1.800a564aa10b6p-2,
    -0x1.91bba8a906b7fp-2, -0x1.a33760adbb56ep-2, -0x1.b47ebf91d417cp-2, -0x1.c592faf028f8ep-2,
    -0x1.d6753e1a43e85p-2, -0x1.e726aa2157f61p-2, -0x1.f7a8567cd1cbdp-2, -0x1.03fda8a93ea1cp-1,
    -0x1.0c10500ed4feep-1, -0x1.140c9fb5a8f7fp-1, -0x1.1bf311daefb44p-1, -0x1.23c41d317edc1p-1,
    -0x1.2b80347f8250cp-1, -0x1.3327c6a752222p-1, -0x1.3abb3fb080128p-1, -0x1.423b07f5114e5p-1,
    -0x1.49a784b147150p-1, -0x1.5101187e9b871p-1, -0x1.5848226c6c7c3p-1, -0x1.5f7cff3de8f2ap-1,
    -0x1.66a008d8ede91p-1, -0x1.6db19694a04aap-1, -0x1.74b1fd6b5e715p-1, -0x1.7ba18f99ce2a5p-1,
    -0x1.82809d4d79ba9p-1, -0x1.894f749cf5047p-1, -0x1.900e615bac2f7p-1, -0x1.96bdad16f5160p-1,
    -0x1.9d5d9fe08baefp-1, -0x1.a3ee7f3e4a7ebp-1, -0x1.aa708f4de7fe0p-1, -0x1.b0e4125ca68fep-1,
    -0x1.b74948f9a72b0p-1, -0x1.bda071b77c9e2p-1, -0x1.c3e9ca4193f99p-1, -0x1.ca258dd397814p-1,
    -0x1.d053f6d543325p-1, -0x1.d6753dfedaa39p-1, -0x1.dc899aa874b31p-1, -0x1.e29142f2e8b3dp-1,
    -0x1.e88c6b41b14aep-1, -0x1.ee7b4718b4414p-1, -0x1.f45e08c87b090p-1, -0x1.fa34e117d8785p-1,
    -0x1.0000000000000p+0,
};

const b = [_]f64{
    0x1.7154765bab2b5p+0, -0x1.71574d6939f6cp-1, 0x1.ec60b584ca47fp-2,
};

const c = [_]f64{
    0x1.71547652b8314p+0, -0x1.71547652b7f67p-1, 0x1.ec709db872c63p-2, -0x1.715476b064cdfp-2,
    0x1.277c72c15441fp-2, -0x1.ec4ff35b64c31p-3,
};

fn log2fSpecial(value: f32, bits: u32) f32 {
    const doubled_magnitude = bits << 1;
    if (doubled_magnitude == 0) return @bitCast(minus_infinity_bits);
    if (bits == infinity_bits) return value;
    if (doubled_magnitude > 0xff00_0000) return @bitCast(bits | quiet_nan_bit);
    return @bitCast(canonical_nan_bits);
}

fn log2fFinite(input_bits: u32) f32 {
    @setFloatMode(.strict);

    var bits = input_bits;
    if (bits < (1 << 23)) {
        const leading_zeros: u32 = @clz(bits);
        const shift: u5 = @intCast(leading_zeros - 8);
        bits <<= shift;
        bits -%= (leading_zeros - 8) << 23;
    }

    const exponent: i32 = (@as(i32, @bitCast(bits)) >> 23) - 0x7f;
    const mantissa = bits & 0x007f_ffff;
    if (mantissa == 0) return @floatFromInt(exponent);

    const table_index: usize = @intCast((mantissa + (1 << 16)) >> 17);
    const encoded: u64 = (@as(u64, mantissa) << 29) | (@as(u64, 0x3ff) << 52);
    const normalized: f64 = @bitCast(encoded);
    const z = normalized * ix[table_index] - 1.0;
    const z_squared = z * z;
    const exponent_low: f64 = @as(f64, @floatFromInt(exponent)) - lix[table_index];
    const polynomial = (exponent_low + z * b[0]) + z_squared * (b[1] + z * b[2]);
    const lower: f32 = @floatCast(polynomial);
    const upper: f32 = @floatCast(polynomial + 0x1.661p-32);
    if (lower == upper) return lower;

    const coefficient_0 = c[0] + z * c[1];
    const coefficient_2 = c[2] + z * c[3];
    const coefficient_4 = c[4] + z * c[5];
    const accurate_polynomial = coefficient_0 + z_squared * (coefficient_2 + z_squared * coefficient_4);
    return @floatCast(exponent_low + z * accurate_polynomial);
}

// Correctly rounded binary32 log2, ported from CORE-MATH
// src/binary32/log2/log2f.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa.
// Copyright (c) 2022-2026 Alexei Sibidanov.
//
// This code is derived from CORE-MATH and used under the MIT License.

pub fn log2f(value: f32) f32 {
    const bits: u32 = @bitCast(value);
    const magnitude = bits & magnitude_mask;
    if (bits == 0 or (bits & sign_mask) != 0 or magnitude >= infinity_bits) {
        return log2fSpecial(value, bits);
    }
    return log2fFinite(bits);
}
