// Correctly rounded inverse hyperbolic tangent for binary16, binary32, and binary64.
//
// The finite kernels are ported from CORE-MATH at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa:
//   src/binary16/atanh/atanhf16.c (Copyright (c) 2025 Paul Zimmermann).
//   src/binary32/atanh/atanhf.c (Copyright (c) 2023-2025 Alexei Sibidanov).
//   src/binary64/atanh/atanh.c (Copyright (c) 2023-2025 Alexei Sibidanov).
//
// This code is derived from CORE-MATH and is used under the MIT License:
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

const sign_mask_h: u16 = 0x8000;
const magnitude_mask_h: u16 = 0x7fff;
const infinity_bits_h: u16 = 0x7c00;
const quiet_nan_bit_h: u16 = 0x0200;

const sign_mask_f: u32 = 0x8000_0000;
const magnitude_mask_f: u32 = 0x7fff_ffff;
const infinity_bits_f: u32 = 0x7f80_0000;
const quiet_nan_bit_f: u32 = 0x0040_0000;

const sign_mask_d: u64 = 0x8000_0000_0000_0000;
const magnitude_mask_d: u64 = 0x7fff_ffff_ffff_ffff;
const infinity_bits_d: u64 = 0x7ff0_0000_0000_0000;
const quiet_nan_bit_d: u64 = 0x0008_0000_0000_0000;

fn copySignH(magnitude: f16, bits: u16) f16 {
    const magnitude_bits: u16 = @bitCast(magnitude);
    return @bitCast((magnitude_bits & ~sign_mask_h) | (bits & sign_mask_h));
}

fn copySignF(magnitude: f32, bits: u32) f32 {
    const magnitude_bits: u32 = @bitCast(magnitude);
    return @bitCast((magnitude_bits & ~sign_mask_f) | (bits & sign_mask_f));
}

fn copySignD(magnitude: f64, bits: u64) f64 {
    const magnitude_bits: u64 = @bitCast(magnitude);
    return @bitCast((magnitude_bits & ~sign_mask_d) | (bits & sign_mask_d));
}

const tr = [_]f64{ 0x1.fc07f02p-1, 0x1.f44659ep-1, 0x1.ecc07b3p-1, 0x1.e573ac9p-1, 0x1.de5d6e4p-1, 0x1.d77b655p-1, 0x1.d0cb58fp-1, 0x1.ca4b305p-1, 0x1.c3f8f02p-1, 0x1.bdd2b8ap-1, 0x1.b7d6c3ep-1, 0x1.b20364p-1, 0x1.ac5701bp-1, 0x1.a6d01a7p-1, 0x1.a16d3f9p-1, 0x1.9c2d14fp-1, 0x1.970e4f8p-1, 0x1.920fb4ap-1, 0x1.8d3018dp-1, 0x1.886e5f1p-1, 0x1.83c977bp-1, 0x1.7f405fdp-1, 0x1.7ad2209p-1, 0x1.767dce4p-1, 0x1.724287fp-1, 0x1.6e1f76bp-1, 0x1.6a13cd1p-1, 0x1.661ec6ap-1, 0x1.623fa77p-1, 0x1.5e75bb9p-1, 0x1.5ac056bp-1, 0x1.571ed3cp-1, 0x1.5390949p-1, 0x1.5015015p-1, 0x1.4cab887p-1, 0x1.49539e4p-1, 0x1.460cbc8p-1, 0x1.42d6626p-1, 0x1.3fb014p-1, 0x1.3c995a4p-1, 0x1.3991c2cp-1, 0x1.3698df4p-1, 0x1.33ae45bp-1, 0x1.30d1901p-1, 0x1.2e025cp-1, 0x1.2b404adp-1, 0x1.288b013p-1, 0x1.25e2271p-1, 0x1.2345679p-1, 0x1.20b470cp-1, 0x1.1e2ef3bp-1, 0x1.1bb4a4p-1, 0x1.1945381p-1, 0x1.16e0689p-1, 0x1.1485f0ep-1, 0x1.12358e7p-1, 0x1.0fef011p-1, 0x1.0db20a9p-1, 0x1.0b7e6ecp-1, 0x1.0953f39p-1, 0x1.073260ap-1, 0x1.05197f8p-1, 0x1.03091b5p-1, 0x1.010101p-1 };
const tl = [_]f64{ 0x1.fe02a69106789p-9, 0x1.7b91b1155b11bp-7, 0x1.39e87ba1ebd6p-6, 0x1.b42dd713971bfp-6, 0x1.16536ee637ae1p-5, 0x1.51b073c96183fp-5, 0x1.8c345da019b21p-5, 0x1.c5e5492abc743p-5, 0x1.fec912fbbeabbp-5, 0x1.1b72ad33f67ap-4, 0x1.371fc1f6e8f74p-4, 0x1.526e5e5a1b438p-4, 0x1.6d60fe601d21dp-4, 0x1.87fa06438c911p-4, 0x1.a23bc223ab563p-4, 0x1.bc28673a58cd6p-4, 0x1.d5c216b8fbb91p-4, 0x1.ef0adcaec5936p-4, 0x1.040259530d041p-3, 0x1.1058bf8d24ad5p-3, 0x1.1c898c09d99fbp-3, 0x1.2895a13e286a3p-3, 0x1.347dd9a447d55p-3, 0x1.404308716a7e4p-3, 0x1.4be5f963b78a1p-3, 0x1.5767718015a6cp-3, 0x1.62c82f3a5c795p-3, 0x1.6e08eab13a1e4p-3, 0x1.792a55fe147a2p-3, 0x1.842d1d9928b17p-3, 0x1.8f11e873a62c7p-3, 0x1.99d958207e08bp-3, 0x1.a484090c1bb0ap-3, 0x1.af129324b786bp-3, 0x1.b9858970710fbp-3, 0x1.c3dd7a6ddad4dp-3, 0x1.ce1af0b65f3ebp-3, 0x1.d83e725022f3ep-3, 0x1.e2488197c6c26p-3, 0x1.ec399d3d68ccp-3, 0x1.f6123fac028acp-3, 0x1.ffd2e07e7f498p-3, 0x1.04bdf9e3b26d2p-2, 0x1.0986f4fa93521p-2, 0x1.0e4498651cc8cp-2, 0x1.12f719595efbcp-2, 0x1.179eabb0a99a1p-2, 0x1.1c3b81e933c25p-2, 0x1.20cdcd0e0ab6ep-2, 0x1.2555bcf50f7cbp-2, 0x1.29d37ff34b08bp-2, 0x1.2e47437640268p-2, 0x1.32b1338401d71p-2, 0x1.37117b5c147b6p-2, 0x1.3b6844a13fc23p-2, 0x1.3fb5b857f6f42p-2, 0x1.43f9fe2f7ce67p-2, 0x1.48353d11488dfp-2, 0x1.4c679b014ee3ap-2, 0x1.50913cc03686bp-2, 0x1.54b2468259498p-2, 0x1.58cadb57d7989p-2, 0x1.5cdb1dcaa1765p-2, 0x1.60e32f46788d9p-2 };
const ln2n = [_]f64{ 0x1.62e42fedb2a44p-2, 0x1.62e42feeab21ap-1, 0x1.0a2b23f33e789p+0, 0x1.62e42fef27604p+0, 0x1.bb9d3beb1048p+0, 0x1.0a2b23f37c97ep+1, 0x1.3687a9f1710bcp+1, 0x1.62e42fef657fap+1, 0x1.8f40b5ed59f38p+1, 0x1.bb9d3beb4e676p+1, 0x1.e7f9c1e942db4p+1, 0x1.0a2b23f39ba79p+2, 0x1.205966f295e18p+2, 0x1.3687a9f1901b7p+2, 0x1.4cb5ecf08a556p+2, 0x1.62e42fef848f5p+2, 0x1.791272ee7ec93p+2, 0x1.8f40b5ed79032p+2, 0x1.a56ef8ec733d1p+2, 0x1.bb9d3beb6d77p+2, 0x1.d1cb7eea67b0fp+2, 0x1.e7f9c1e961eaep+2, 0x1.fe2804e85c24dp+2, 0x1.0a2b23f3ab2f6p+3 };
const b = [_]f64{ 0x1.fffffffce5a6ap-2, -0x1.0001f81ec0ab8p-2, 0x1.555a0f53d79a5p-3 };
const atanhf_small_c = [_]f64{ 0x1.5555555555527p-2, 0x1.9999999ba4ee8p-3, 0x1.24922c280990ap-3, 0x1.c8236aae809c6p-4 };
const atanhf_refine_c = [_]f64{ 0x1p-1, -0x1.000000000001bp-2, 0x1.55555555555bap-3, -0x1.fffffff26d72ep-4, 0x1.99999989035p-4, -0x1.555c39cb9ee8p-4, 0x1.24992d8b014a1p-4 };
fn atanhfSpecial(bits: u32) f32 {
    const doubled = bits << 1;
    if (doubled == 0x7f00_0000) return copySignF(@bitCast(infinity_bits_f), bits);
    if (doubled > 0xff00_0000) return @bitCast(bits | quiet_nan_bit_f);
    return @bitCast(@as(u32, 0x7fc0_0000));
}

/// Return correctly rounded binary32 atanh with caller-selected nearest-even mode.
pub fn atanhf(value: f32) f32 {
    @setFloatMode(.strict);
    const bits: u32 = @bitCast(value);
    const doubled = bits << 1;
    if (doubled < 0x7a30_0000 or doubled >= 0x7f00_0000) {
        if (doubled >= 0x7f00_0000) return atanhfSpecial(bits);
        if (doubled < 0x7371_3744) {
            if (doubled == 0) return value;
            return @mulAdd(f32, value, 0x1p-25, value);
        }
        const z: f64 = @floatCast(value);
        const z2 = z * z;
        const z4 = z2 * z2;
        const r = atanhf_small_c[0] + z2 * atanhf_small_c[1] + z4 * (atanhf_small_c[2] + z2 * atanhf_small_c[3]);
        return @floatCast(z + (z * z2) * r);
    }

    const sign: f64 = if ((bits & sign_mask_f) == 0) 1.0 else -1.0;
    const e: u32 = doubled >> 24;
    const shift: u5 = @intCast(126 - e);
    const md = ((bits << 8) | (@as(u32, 1) << 31)) >> shift;
    var mn = 0 -% md;
    const leading_zeros: u32 = @clz(mn);
    const nz: u32 = leading_zeros + 1;
    mn <<= @intCast(nz);
    const jn: usize = @intCast(mn >> 26);
    const jd: usize = @intCast(md >> 26);
    const tn_bits: u64 = (@as(u64, mn) << 20) | (@as(u64, 1023) << 52);
    const td_bits: u64 = (@as(u64, md) << 20) | (@as(u64, 1023) << 52);
    const tn: f64 = @bitCast(tn_bits);
    const td: f64 = @bitCast(td_bits);
    const zn = tn * tr[jn] - 1.0;
    const zd = td * tr[jd] - 1.0;
    const zn2 = zn * zn;
    const zd2 = zd * zd;
    const rn = ((tl[jn] - ln2n[nz - 1]) + zn * b[0]) + zn2 * (b[1] + zn * b[2]);
    const rd = (tl[jd] + zd * b[0]) + zd2 * (b[1] + zd * b[2]);
    const r = sign * (rd - rn);
    var upper: f32 = @floatCast(r);
    const lower: f32 = @floatCast(r + sign * 0.226e-9);
    if (upper != lower) {
        const zn4 = zn2 * zn2;
        const zd4 = zd2 * zd2;
        const fn_value = zn * (((atanhf_refine_c[0] + zn * atanhf_refine_c[1]) + zn2 * (atanhf_refine_c[2] + zn * atanhf_refine_c[3])) +
            zn4 * ((atanhf_refine_c[4] + zn * atanhf_refine_c[5]) + zn2 * atanhf_refine_c[6]));
        var refined = fn_value + 0x1.0ca86c3898dp-50 * @as(f64, @floatFromInt(nz));
        refined += tl[jn];
        const en = @as(f64, @floatFromInt(nz)) * 0x1.62e42fefa3ap-2;
        const fd = zd * (((atanhf_refine_c[0] + zd * atanhf_refine_c[1]) + zd2 * (atanhf_refine_c[2] + zd * atanhf_refine_c[3])) +
            zd4 * ((atanhf_refine_c[4] + zd * atanhf_refine_c[5]) + zd2 * atanhf_refine_c[6]));
        refined = fd + tl[jd] - refined + en;
        upper = @floatCast(sign * refined);
    }
    return upper;
}

const atanhd_ch_zero = [_][2]f64{
    .{ 0x1.5555555555555p-2, 0x1.5555555555555p-56 }, // degree 3
    .{ 0x1.999999999999ap-3, -0x1.999999999611cp-57 }, // degree 5
    .{ 0x1.2492492492492p-3, 0x1.2492490f76b25p-57 }, // degree 7
    .{ 0x1.c71c71c71c71cp-4, 0x1.c71cd5c38a112p-58 }, // degree 9
    .{ 0x1.745d1745d1746p-4, -0x1.7556c4165f4cap-59 }, // degree 11
    .{ 0x1.3b13b13b13b14p-4, -0x1.b893c3b36052ep-59 }, // degree 13
    .{ 0x1.1111111111105p-4, 0x1.4e1afd723ed1fp-59 }, // degree 15
    .{ 0x1.e1e1e1e1e2678p-5, -0x1.f86ea96fb1435p-59 }, // degree 17
    .{ 0x1.af286bc9f90ccp-5, 0x1.1e51a6e54fde9p-60 }, // degree 19
    .{ 0x1.8618618c779b6p-5, -0x1.ab913de95c3bfp-61 }, // degree 21
    .{ 0x1.642c84aa383ebp-5, 0x1.632e747641b12p-59 }, // degree 23
    .{ 0x1.47ae2d205013cp-5, -0x1.0c9617e7bcff2p-60 }, // degree 25
    .{ 0x1.2f664d60473f9p-5, 0x1.3adb3e2b7f35ep-61 },
};
const atanhd_cl_zero = [_]f64{
    0x1.1a9a91fd692afp-5, // degree 29
    0x1.06dfbb35e7f44p-5, // degree 31
    0x1.037bed4d7588fp-5, // degree 33
    0x1.5aca6d6d720d6p-6, // degree 35
    0x1.99ea5700d53a5p-5,
};
const atanhd_c_small = [_]f64{ 0x1.999999999999ap-3, 0x1.2492492492244p-3, 0x1.c71c71c79715fp-4, 0x1.745d16f777723p-4, 0x1.3b13ca4174634p-4, 0x1.110c9724989bdp-4, 0x1.e2d17608a5b2ep-5, 0x1.a0b56308cba0bp-5, 0x1.fb6341208ad2ep-5 };
const atanhd_B = [_][2]i32{ .{ 301, 27565 }, .{ 7189, 24786 }, .{ 13383, 22167 }, .{ 18923, 19696 }, .{ 23845, 17361 }, .{ 28184, 15150 }, .{ 31969, 13054 }, .{ 35231, 11064 }, .{ 37996, 9173 }, .{ 40288, 7372 }, .{ 42129, 5657 }, .{ 43542, 4020 }, .{ 44546, 2457 }, .{ 45160, 962 }, .{ 45399, -468 }, .{ 45281, -1838 }, .{ 44821, -3151 }, .{ 44032, -4412 }, .{ 42929, -5622 }, .{ 41522, -6786 }, .{ 39825, -7905 }, .{ 37848, -8982 }, .{ 35602, -10020 }, .{ 33097, -11020 }, .{ 30341, -11985 }, .{ 27345, -12916 }, .{ 24115, -13816 }, .{ 20661, -14685 }, .{ 16989, -15526 }, .{ 13107, -16339 }, .{ 9022, -17126 }, .{ 4740, -17889 } };
const atanhd_r1 = [_]f64{ 0x1p+0, 0x1.f5076p-1, 0x1.ea4bp-1, 0x1.dfc98p-1, 0x1.d5818p-1, 0x1.cb72p-1, 0x1.c199cp-1, 0x1.b7f76p-1, 0x1.ae8ap-1, 0x1.a5504p-1, 0x1.9c492p-1, 0x1.93738p-1, 0x1.8ace6p-1, 0x1.8258ap-1, 0x1.7a114p-1, 0x1.71f76p-1, 0x1.6a09ep-1, 0x1.6247ep-1, 0x1.5ab08p-1, 0x1.5342cp-1, 0x1.4bfdap-1, 0x1.44e08p-1, 0x1.3dea6p-1, 0x1.371a8p-1, 0x1.306fep-1, 0x1.29e9ep-1, 0x1.2387ap-1, 0x1.1d488p-1, 0x1.172b8p-1, 0x1.11302p-1, 0x1.0b558p-1, 0x1.059bp-1, 0x1p-1 };
const atanhd_r2 = [_]f64{ 0x1p+0, 0x1.ffa74p-1, 0x1.ff4eap-1, 0x1.fef62p-1, 0x1.fe9dap-1, 0x1.fe452p-1, 0x1.fdeccp-1, 0x1.fd946p-1, 0x1.fd3c2p-1, 0x1.fce3ep-1, 0x1.fc8bcp-1, 0x1.fc33ap-1, 0x1.fbdbap-1, 0x1.fb83ap-1, 0x1.fb2bcp-1, 0x1.fad3ep-1, 0x1.fa7c2p-1, 0x1.fa246p-1, 0x1.f9ccap-1, 0x1.f975p-1, 0x1.f91d8p-1, 0x1.f8c6p-1, 0x1.f86e8p-1, 0x1.f8172p-1, 0x1.f7bfep-1, 0x1.f768ap-1, 0x1.f7116p-1, 0x1.f6ba4p-1, 0x1.f6632p-1, 0x1.f60c2p-1, 0x1.f5b52p-1, 0x1.f55e4p-1, 0x1.f5076p-1 };
const atanhd_l1 = [_][2]f64{ .{ 0x0p+0, 0x0p+0 }, .{ -0x1.532c1269e2038p-27, 0x1.62e5p-7 }, .{ 0x1.ce42d81b54e84p-27, 0x1.62e3cp-6 }, .{ -0x1.25826f815ec3dp-26, 0x1.0a2acp-5 }, .{ 0x1.0db1b1e7cee11p-26, 0x1.62e4ap-5 }, .{ -0x1.1f3a8c6c95003p-26, 0x1.bb9dcp-5 }, .{ -0x1.774cd4fb8c30dp-26, 0x1.0a2b2p-4 }, .{ 0x1.452e56c030a0ap-29, 0x1.3687fp-4 }, .{ 0x1.6b63c4966a79ap-28, 0x1.62e41p-4 }, .{ -0x1.b20a21ccb525ep-28, 0x1.8f40ap-4 }, .{ 0x1.4006cfb3d8f85p-26, 0x1.bb9d1p-4 }, .{ -0x1.cdb026b310c41p-26, 0x1.e7f9bp-4 }, .{ -0x1.69124fdc0f16dp-26, 0x1.0a2b08p-3 }, .{ -0x1.084656cdc2727p-26, 0x1.205958p-3 }, .{ -0x1.376fa8b0357fdp-26, 0x1.3687cp-3 }, .{ 0x1.e56ae55a47b4ap-28, 0x1.4cb5e8p-3 }, .{ 0x1.070ff8834eeb4p-26, 0x1.62e44p-3 }, .{ 0x1.623516109f4fep-26, 0x1.79129p-3 }, .{ -0x1.ec656b95fbdacp-29, 0x1.8f40bp-3 }, .{ 0x1.f0ca2e729f51p-28, 0x1.a56ed8p-3 }, .{ -0x1.7d260a858354ap-26, 0x1.bb9d68p-3 }, .{ 0x1.e7279075503d3p-27, 0x1.d1cb9p-3 }, .{ 0x1.39e1a0a503873p-27, 0x1.e7f9dp-3 }, .{ 0x1.cd86d7b87c3d6p-26, 0x1.fe27d8p-3 }, .{ 0x1.060ab88de341ep-26, 0x1.0a2b24p-2 }, .{ 0x1.20a860d3f939p-28, 0x1.154244p-2 }, .{ -0x1.dacee95fc2f1p-27, 0x1.205974p-2 }, .{ 0x1.45de3a86e0acap-26, 0x1.2b707p-2 }, .{ 0x1.c164cbfb991afp-27, 0x1.3687bp-2 }, .{ 0x1.d3f66b24225efp-26, 0x1.419ec4p-2 }, .{ 0x1.fc023efa144bap-26, 0x1.4cb5f8p-2 }, .{ 0x1.086a8af6f26cp-28, 0x1.57cd28p-2 }, .{ -0x1.05c610ca86c39p-30, 0x1.62e43p-2 } };
const atanhd_l2 = [_][2]f64{ .{ 0x0p+0, 0x0p+0 }, .{ -0x1.37e152a129e4ep-28, 0x1.632p-12 }, .{ -0x1.3f6c916b8be9cp-26, 0x1.63p-11 }, .{ 0x1.20505936739d5p-26, 0x1.0a24p-10 }, .{ -0x1.23e2e8cb541bap-26, 0x1.62dcp-10 }, .{ -0x1.acb7983ac4f5ep-32, 0x1.bbap-10 }, .{ 0x1.6f7c7689c63aep-28, 0x1.0a2ap-9 }, .{ 0x1.f5ca695b4c58bp-30, 0x1.368cp-9 }, .{ -0x1.c6c18bd953226p-27, 0x1.62e6p-9 }, .{ 0x1.7a516c34846bdp-26, 0x1.8f46p-9 }, .{ -0x1.f3b83dd8b853p-27, 0x1.bbap-9 }, .{ -0x1.c3459046e4e57p-31, 0x1.e8p-9 }, .{ 0x1.b5c7e34cb79f6p-38, 0x1.0a2cp-8 }, .{ -0x1.2487e9af9a692p-27, 0x1.205cp-8 }, .{ 0x1.f21bbc4ad79cep-26, 0x1.3687p-8 }, .{ -0x1.550ffc857b731p-29, 0x1.4cb7p-8 }, .{ 0x1.87458ec1b7b34p-27, 0x1.62e2p-8 }, .{ 0x1.103d4fe83ee81p-26, 0x1.7911p-8 }, .{ 0x1.810483d3b398cp-27, 0x1.8f44p-8 }, .{ -0x1.2085cb340608ep-27, 0x1.a573p-8 }, .{ 0x1.12698a119c42fp-26, 0x1.bb9dp-8 }, .{ -0x1.edb8c172b4c33p-26, 0x1.d1ccp-8 }, .{ -0x1.8b55b87a5e238p-26, 0x1.e7fep-8 }, .{ 0x1.be5e17763f78ap-26, 0x1.fe2bp-8 }, .{ -0x1.c2d496790073ep-30, 0x1.0a2a8p-7 }, .{ 0x1.6542f523abeecp-26, 0x1.1541p-7 }, .{ -0x1.b7fdbe5b193f8p-26, 0x1.205ap-7 }, .{ 0x1.fa4d42fe30c7cp-26, 0x1.2b7p-7 }, .{ 0x1.0d46ad04adc86p-26, 0x1.36888p-7 }, .{ -0x1.1c22d02d17c4cp-26, 0x1.419fp-7 }, .{ 0x1.a7d1e330dcccep-30, 0x1.4cb7p-7 }, .{ 0x1.187025e656ba3p-31, 0x1.57cdp-7 }, .{ -0x1.532c1269e2038p-27, 0x1.62e5p-7 } };
const atanhd_db = [_][3]f64{
    .{ 0x1.2dbb7b1c91363p-2, 0x1.36f33d51c264dp-2, 0x1p-56 },
    .{ 0x1.c493dc899e4a5p-2, 0x1.e611aa58ab608p-2, -0x1p-56 },
};
const atanhd_t1 = [_]f64{ 0x1p+0, 0x1.ea4afap-1, 0x1.d5818ep-1, 0x1.c199bep-1, 0x1.ae89f98p-1, 0x1.9c4918p-1, 0x1.8ace54p-1, 0x1.7a1147p-1, 0x1.6a09e68p-1, 0x1.5ab07ep-1, 0x1.4bfdad8p-1, 0x1.3dea65p-1, 0x1.306fe08p-1, 0x1.2387a7p-1, 0x1.172b84p-1, 0x1.0b5587p-1, 0x1p-1 };
const atanhd_t2 = [_]f64{ 0x1p+0, 0x1.fe9d968p-1, 0x1.fd3c228p-1, 0x1.fbdba38p-1, 0x1.fa7c18p-1, 0x1.f91d8p-1, 0x1.f7bfdbp-1, 0x1.f663278p-1, 0x1.f507658p-1, 0x1.f3ac948p-1, 0x1.f252b38p-1, 0x1.f0f9c2p-1, 0x1.efa1bfp-1, 0x1.ee4aaap-1, 0x1.ecf483p-1, 0x1.eb9f488p-1 };
const atanhd_t3 = [_]f64{ 0x1p+0, 0x1.ffe9d2p-1, 0x1.ffd3a58p-1, 0x1.ffbd798p-1, 0x1.ffa74e8p-1, 0x1.ff91248p-1, 0x1.ff7afb8p-1, 0x1.ff64d38p-1, 0x1.ff4eac8p-1, 0x1.ff38868p-1, 0x1.ff22618p-1, 0x1.ff0c3dp-1, 0x1.fef61ap-1, 0x1.fedff78p-1, 0x1.fec9d68p-1, 0x1.feb3b6p-1 };
const atanhd_t4 = [_]f64{ 0x1p+0, 0x1.fffe9dp-1, 0x1.fffd3ap-1, 0x1.fffbd78p-1, 0x1.fffa748p-1, 0x1.fff9118p-1, 0x1.fff7ae8p-1, 0x1.fff64cp-1, 0x1.fff4e9p-1, 0x1.fff386p-1, 0x1.fff2238p-1, 0x1.fff0c08p-1, 0x1.ffef5d8p-1, 0x1.ffedfa8p-1, 0x1.ffec98p-1, 0x1.ffeb35p-1 };
const atanhd_LL = [_][17][3]f64{
    .{
        .{ 0x0p+0, 0x0p+0, 0x0p+0 },
        .{ 0x1.62e432b24p-6, -0x1.745af34bb54b8p-42, -0x1.17e3ec05cde7p-97 },
        .{ 0x1.62e42e4a8p-5, 0x1.111a4eadf312p-44, 0x1.cff3027abb119p-93 },
        .{ 0x1.0a2b233f1p-4, -0x1.88ac4ec78af8p-42, 0x1.4fa087ca75dfdp-93 },
        .{ 0x1.62e43056cp-4, 0x1.6bd65e8b0b7p-46, -0x1.b18e160362c24p-95 },
        .{ 0x1.bb9d3cbd6p-4, 0x1.de14aa55ec2bp-42, -0x1.c6ac3f1862a6bp-94 },
        .{ 0x1.0a2b244dap-3, 0x1.94def487fea7p-42, -0x1.dead1a4581acfp-94 },
        .{ 0x1.3687aa9b78p-3, 0x1.9cec9a50db22p-43, 0x1.34a70684f8e0ep-93 },
        .{ 0x1.62e42fabap-3, -0x1.d69047a3aebp-44, -0x1.4e061f79144e2p-95 },
        .{ 0x1.8f40b56d28p-3, 0x1.de7d755fd2e2p-42, 0x1.bdc7ecf001489p-94 },
        .{ 0x1.bb9d3b61fp-3, 0x1.c14f1445b12p-46, 0x1.a1d78cbdc5b58p-93 },
        .{ 0x1.e7f9c11f08p-3, -0x1.6e3e0000dae7p-43, 0x1.6a4559fadde98p-94 },
        .{ 0x1.0a2b242ec4p-2, 0x1.bb7cf852a5fe8p-42, 0x1.a6aef11ee43bdp-93 },
        .{ 0x1.205966c764p-2, 0x1.ad3a5f214294p-45, 0x1.5cc344fa10652p-93 },
        .{ 0x1.3687a98aacp-2, 0x1.1623671842fp-45, -0x1.0b428fe1f9e43p-94 },
        .{ 0x1.4cb5ec93f4p-2, 0x1.3d50980ea513p-42, 0x1.67f0ea083b1c4p-93 },
        .{ 0x1.62e42fefa4p-2, -0x1.8432a1b0e264p-44, 0x1.803f2f6af40f3p-93 },
    },
    .{
        .{ 0x0p+0, 0x0p+0, 0x0p+0 },
        .{ 0x1.62e462b4p-10, 0x1.061d003b97318p-42, 0x1.d7faee66a2e1ep-93 },
        .{ 0x1.62e44c92p-9, 0x1.95a7bff5e239p-42, -0x1.f7e788a87135p-95 },
        .{ 0x1.0a2b1e33p-8, 0x1.2a3a1a65aa3ap-43, -0x1.54599c9605442p-93 },
        .{ 0x1.62e4367cp-8, -0x1.4a995b6d9ddcp-45, -0x1.56bb79b254f33p-100 },
        .{ 0x1.bb9d449ap-8, 0x1.8a119c42e9bcp-42, -0x1.8ecf7d8d661f1p-93 },
        .{ 0x1.0a2b1f19p-7, 0x1.8863771bd10a8p-42, 0x1.e9731de7f0155p-94 },
        .{ 0x1.3687ad11p-7, 0x1.e026a347ca1c8p-42, 0x1.fadc62522444dp-97 },
        .{ 0x1.62e436f28p-7, 0x1.25b84f71b70b8p-42, -0x1.fcb3f98612d27p-96 },
        .{ 0x1.8f40b7b38p-7, -0x1.62a0a4fd4758p-43, 0x1.3cb3c35d9f6a1p-93 },
        .{ 0x1.bb9d3abbp-7, -0x1.0ec48f94d786p-42, -0x1.6b47d410e4cc7p-93 },
        .{ 0x1.e7f9bb23p-7, 0x1.e4415cbc97ap-43, -0x1.3729fdb677231p-93 },
        .{ 0x1.0a2b22478p-6, -0x1.cb73f4505b03p-42, -0x1.1b3b3a3bc370ap-93 },
        .{ 0x1.2059691e8p-6, -0x1.abcc3412f264p-43, -0x1.fe6e998e48673p-95 },
        .{ 0x1.3687a768p-6, -0x1.43901e5c97a9p-42, 0x1.b54cdd52a5d88p-96 },
        .{ 0x1.4cb5eb5d8p-6, -0x1.8f106f00f13b8p-42, -0x1.8f793f5fce148p-93 },
        .{ 0x1.62e432b24p-6, -0x1.745af34bb54b8p-42, -0x1.17e3ec05cde7p-97 },
    },
    .{
        .{ 0x0p+0, 0x0p+0, 0x0p+0 },
        .{ 0x1.62e7bp-14, -0x1.868625640a68p-44, -0x1.34bf0db910f65p-93 },
        .{ 0x1.62e35f6p-13, -0x1.2ee3d96b696ap-43, 0x1.a2948cd558655p-94 },
        .{ 0x1.0a2b4b2p-12, 0x1.53edbcf1165p-47, -0x1.cfc26ccf6d0e4p-97 },
        .{ 0x1.62e4be1p-12, 0x1.783e334614p-52, -0x1.04b96da30e63ap-93 },
        .{ 0x1.bb9e085p-12, -0x1.60785f20acb2p-43, -0x1.f33369bf7dff1p-96 },
        .{ 0x1.0a2b94dp-11, 0x1.fd4b3a273353p-42, -0x1.685a35575eff1p-96 },
        .{ 0x1.368810f8p-11, 0x1.7ded26dc813p-47, -0x1.4c4d1abca79bfp-96 },
        .{ 0x1.62e47878p-11, 0x1.7d2bee9a1f63p-42, 0x1.860233b7ad13p-93 },
        .{ 0x1.8f40cb48p-11, -0x1.af034eaf471cp-42, 0x1.ae748822d57b7p-94 },
        .{ 0x1.bb9d094p-11, -0x1.7a223013a20fp-42, -0x1.1e499087075b6p-93 },
        .{ 0x1.e7fa32c8p-11, -0x1.b2e67b1b59bdp-43, -0x1.54a41eda30fa6p-93 },
        .{ 0x1.0a2b237p-10, -0x1.7ad97ff4ac7ap-44, 0x1.f932da91371ddp-93 },
        .{ 0x1.2059a338p-10, -0x1.96422d90df4p-44, -0x1.90800fbbf2ed3p-94 },
        .{ 0x1.36879824p-10, 0x1.0f9054001812p-44, 0x1.9567e01e48f9ap-93 },
        .{ 0x1.4cb602cp-10, -0x1.0d709a5ec0b5p-43, 0x1.253dfd44635d2p-94 },
        .{ 0x1.62e462b4p-10, 0x1.061d003b97318p-42, 0x1.d7faee66a2e1ep-93 },
    },
    .{
        .{ 0x0p+0, 0x0p+0, 0x0p+0 },
        .{ 0x1.63007cp-18, -0x1.db0e38e5aaaap-43, 0x1.259a7b94815b9p-93 },
        .{ 0x1.6300f6p-17, 0x1.2b1c75580438p-44, 0x1.78cabba01e3e4p-93 },
        .{ 0x1.0a2115p-16, -0x1.5ff223730759p-42, 0x1.8074feacfe49dp-95 },
        .{ 0x1.62e1ecp-16, -0x1.85d6f6487ce4p-45, 0x1.05485074b9276p-93 },
        .{ 0x1.bba301p-16, -0x1.af5d58a7c921p-43, -0x1.30a8c0fd2ff5fp-93 },
        .{ 0x1.0a32298p-15, 0x1.590faa0883bdp-43, 0x1.95e9bda999947p-93 },
        .{ 0x1.3682f1p-15, 0x1.f0224376efaf8p-42, -0x1.5843c0db50d1p-93 },
        .{ 0x1.62e3d8p-15, -0x1.142c13daed4ap-43, 0x1.c68a61183ce87p-93 },
        .{ 0x1.8f44dd8p-15, -0x1.aa489f399931p-43, 0x1.11c5c376854eap-94 },
        .{ 0x1.bb9601p-15, 0x1.9904d8b6a3638p-42, 0x1.8c89554493c8fp-93 },
        .{ 0x1.e7f744p-15, 0x1.5785ddbe7cba8p-42, 0x1.e7ff3cde7d70cp-94 },
        .{ 0x1.0a2c53p-14, -0x1.6d9e8780d0d5p-43, 0x1.ad9c178106693p-94 },
        .{ 0x1.205d134p-14, -0x1.214a2e893fccp-43, 0x1.548a9500c9822p-93 },
        .{ 0x1.3685e28p-14, 0x1.e23588646103p-43, 0x1.2a97b26da2d88p-94 },
        .{ 0x1.4cb6c18p-14, 0x1.2b7cfcea9e0d8p-42, -0x1.5095048a6b824p-93 },
        .{ 0x1.62e7bp-14, -0x1.868625640a68p-44, -0x1.34bf0db910f65p-93 },
    },
};
const atanhd_ch_refine = [_][2]f64{
    .{ 0x1p-1, 0x1.24b67ee516e3bp-111 }, .{ -0x1p-2, -0x1.932ce43199a8dp-110 }, .{ 0x1.5555555555555p-3, 0x1.55540c15cf91fp-57 },
};
const atanhd_cl_refine = [_]f64{ -0x1p-3, 0x1.9999999a0754fp-4, -0x1.55555555c3157p-4 };
const atanhd_c_refine = [_]f64{ -0x1p+0, 0x1.555555555553p+0, -0x1.fffffffffffap+0, 0x1.99999e33a6366p+1, -0x1.555559ef9525fp+2 };
fn fastTwoSum(x: f64, y: f64, low: *f64) f64 {
    const sum = x + y;
    const z = sum - x;
    low.* = y - z;
    return sum;
}

fn fastTwoSub(x: f64, y: f64, low: *f64) f64 {
    const sum = x - y;
    const z = x - sum;
    low.* = z - y;
    return sum;
}

fn addDoubleDouble(x_high: f64, x_low: f64, c_high: f64, c_low: f64, low: *f64) f64 {
    const sum = x_high + c_high;
    const delta = sum - x_high;
    low.* = ((c_high - delta) + (x_high + (delta - sum))) + (x_low + c_low);
    return sum;
}

fn mulDoubleDoubleAcc(x_high: f64, x_low: f64, c_high: f64, c_low: f64, low: *f64) f64 {
    const high_low = c_high * x_low;
    const low_high = c_low * x_high;
    const high_high = c_high * x_high;
    const high_error = @mulAdd(f64, c_high, x_high, -high_high);
    const middle = high_error + (low_high + high_low);
    return fastTwoSum(high_high, middle, low);
}

fn mulDoubleDouble(x_high: f64, x_low: f64, c_high: f64, c_low: f64, low: *f64) f64 {
    const high_low = c_high * x_low;
    const low_high = c_low * x_high;
    const high_high = c_high * x_high;
    const high_error = @mulAdd(f64, c_high, x_high, -high_high);
    low.* = high_error + (low_high + high_low);
    return high_high;
}

fn mulDoubleDoubleD(x_high: f64, x_low: f64, c: f64, low: *f64) f64 {
    const high = x_high * c;
    low.* = @mulAdd(f64, c, x_high, -high) + x_low * c;
    return high;
}

fn polyDd(x_high: f64, x_low: f64, comptime n: usize, coefficients: *const [n][2]f64, low: *f64) f64 {
    var index: usize = n - 1;
    var high: f64 = undefined;
    var coefficient_low: f64 = undefined;
    high = fastTwoSum(coefficients[index][0], low.*, &coefficient_low);
    coefficient_low += coefficients[index][1];
    while (index > 0) {
        index -= 1;
        high = mulDoubleDoubleAcc(x_high, x_low, high, coefficient_low, &coefficient_low);
        var temporary_low: f64 = undefined;
        const temporary_high = fastTwoSum(coefficients[index][0], high, &temporary_low);
        high = temporary_high;
        coefficient_low += temporary_low + coefficients[index][1];
    }
    low.* = coefficient_low;
    return high;
}

fn asAtanhZero(value: f64) f64 {
    const x_squared = value * value;
    const x_squared_low = @mulAdd(f64, value, value, -x_squared);
    const y2 = x_squared * (atanhd_cl_zero[0] + x_squared * (atanhd_cl_zero[1] + x_squared *
        (atanhd_cl_zero[2] + x_squared * (atanhd_cl_zero[3] + x_squared * atanhd_cl_zero[4]))));
    var y2_work = y2;
    var y1 = polyDd(x_squared, x_squared_low, 13, &atanhd_ch_zero, &y2_work);
    y1 = mulDoubleDoubleD(y1, y2_work, value, &y2_work);
    y1 = mulDoubleDoubleAcc(x_squared, x_squared_low, y1, y2_work, &y2_work);
    var y1_low: f64 = undefined;
    const y0 = fastTwoSum(value, y1, &y1_low);
    y1 = fastTwoSum(y1_low, y2_work, &y2_work);
    var y1_bits: u64 = @bitCast(y1);
    if ((y1_bits & (~@as(u64, 0) >> 12)) == 0) {
        const y2_bits: u64 = @bitCast(y2_work);
        if (((y2_bits ^ y1_bits) >> 63) != 0) y1_bits -%= 1 else y1_bits +%= 1;
        y1 = @bitCast(y1_bits);
    }
    return y0 + y1;
}

fn asAtanhDatabase(value: f64, fallback: f64) f64 {
    const absolute_value = @abs(value);
    const bits: u64 = @bitCast(value);
    if (absolute_value == atanhd_db[0][0]) return copySignD(atanhd_db[0][1] + atanhd_db[0][2], bits);
    if (absolute_value == atanhd_db[1][0]) return copySignD(atanhd_db[1][1] + atanhd_db[1][2], bits);
    return fallback;
}

fn asAtanhRefine(value: f64, zh: f64, zl: f64, a: f64) f64 {
    const bits: u64 = @bitCast(value);
    const value_sign: f64 = if ((bits & sign_mask_d) == 0) 1.0 else -1.0;
    const zh_bits: u64 = @bitCast(zh);
    const exponent_bits = zh_bits >> 52;
    const exponent: i64 = @as(i64, @intCast(exponent_bits)) - 0x3ff;
    var normalized_bits = zh_bits & (~@as(u64, 0) >> 12);
    normalized_bits |= @as(u64, 0x3ff) << 52;
    const normalized: f64 = @bitCast(normalized_bits);
    const v = a - @as(f64, @floatFromInt(exponent)) + 0x1.00008p+0;
    const v_bits: u64 = @bitCast(v);
    const table_index: u64 = (v_bits - (@as(u64, 0x3ff) << 52)) >> (52 - 16);
    const idx1: usize = @intCast(table_index >> 12);
    const idx2: usize = @intCast((table_index >> 8) & 0xf);
    const idx3: usize = @intCast((table_index >> 4) & 0xf);
    const idx4: usize = @intCast(table_index & 0xf);
    const l2h = 0x1.62e42fefa3ap-2;
    const l2l = -0x1.0ca86c3898dp-50;
    const l2ll = 0x1.f97b57a079ap-104;
    const exponent_f = @as(f64, @floatFromInt(exponent));
    const el2 = l2ll * exponent_f;
    const el1 = l2l * exponent_f;
    const el0 = l2h * exponent_f;
    var L: [3]f64 = undefined;
    L[0] = atanhd_LL[0][idx1][0] + atanhd_LL[1][idx2][0] + (atanhd_LL[2][idx3][0] + atanhd_LL[3][idx4][0]) + el0;
    L[1] = atanhd_LL[0][idx1][1] + atanhd_LL[1][idx2][1] + (atanhd_LL[2][idx3][1] + atanhd_LL[3][idx4][1]);
    L[2] = atanhd_LL[0][idx1][2] + atanhd_LL[1][idx2][2] + (atanhd_LL[2][idx3][2] + atanhd_LL[3][idx4][2]);
    const t12 = atanhd_t1[idx1] * atanhd_t2[idx2];
    const t34 = atanhd_t3[idx3] * atanhd_t4[idx4];
    var product_low: f64 = undefined;
    const th = t12 * t34;
    product_low = @mulAdd(f64, t12, t34, -th);
    const dh = th * normalized;
    const dl = @mulAdd(f64, th, normalized, -dh);
    var sh = product_low * normalized;
    var sl = @mulAdd(f64, product_low, normalized, -sh);
    var xl: f64 = undefined;
    var xh = fastTwoSum(dh - 1.0, dl, &xl);
    var zl_bits: u64 = @bitCast(zl);
    const exponent_shift: u64 = @bitCast(exponent << 52);
    zl_bits -%= exponent_shift;
    const zl_normalized: f64 = @bitCast(zl_bits);
    xl += th * zl_normalized;
    xh = addDoubleDouble(xh, xl, sh, sl, &xl);
    sl = xh * (atanhd_cl_refine[0] + xh * (atanhd_cl_refine[1] + xh * atanhd_cl_refine[2]));
    sh = polyDd(xh, xl, 3, &atanhd_ch_refine, &sl);
    sh = mulDoubleDoubleAcc(xh, xl, sh, sl, &sl);
    sh = addDoubleDouble(sh, sl, el1, el2, &sl);
    sh = addDoubleDouble(sh, sl, L[1], L[2], &sl);
    var v2: f64 = undefined;
    const v0 = fastTwoSum(L[0], sh, &v2);
    var v1 = fastTwoSum(v2, sl, &v2);
    var v1_bits: u64 = @bitCast(v1);
    if ((v1_bits & (~@as(u64, 0) >> 12)) == 0) {
        const v2_bits: u64 = @bitCast(v2);
        if (((v2_bits ^ v1_bits) >> 63) != 0) v1_bits -%= 1 else v1_bits +%= 1;
        v1 = @bitCast(v1_bits);
    }
    const v0_bits: u64 = @bitCast(v0);
    const er = (v1_bits +% 1) & (~@as(u64, 0) >> 12);
    const v0_exponent = (v0_bits >> 52) & 0x7ff;
    const v1_exponent = (v1_bits >> 52) & 0x7ff;
    const de = v0_exponent -% v1_exponent;
    const signed_v0 = v0 * value_sign;
    const signed_v1 = v1 * value_sign;
    const result = signed_v0 + signed_v1;
    if (de > 104 or er < 3) return asAtanhDatabase(value, result);
    return result;
}

const atanhd_main_c = [_]f64{ -0x1p+0, 0x1.555555555553p+0, -0x1.fffffffffffap+0, 0x1.99999e33a6366p+1, -0x1.555559ef9525fp+2 };

/// Return correctly rounded binary64 atanh with caller-selected nearest-even mode.
pub fn atanhd(value: f64) f64 {
    @setFloatMode(.strict);
    const bits: u64 = @bitCast(value);
    const absolute_bits = bits & magnitude_mask_d;
    if (absolute_bits == 0x3ff0_0000_0000_0000) return copySignD(@bitCast(infinity_bits_d), bits);
    if (absolute_bits > 0x3ff0_0000_0000_0000 and absolute_bits < infinity_bits_d)
        return @bitCast(@as(u64, 0x7ff8_0000_0000_0000));
    if (absolute_bits >= infinity_bits_d) {
        if (absolute_bits > infinity_bits_d) return @bitCast(bits | quiet_nan_bit_d);
        return @bitCast(@as(u64, 0x7ff8_0000_0000_0000));
    }

    if (absolute_bits < 0x3fd0_0000_0000_0000) {
        if (absolute_bits < 0x3e4d_12ed_0af1_a27f) {
            if (absolute_bits == 0) return value;
            return @mulAdd(f64, value, 0x1p-55, value);
        }
        const x_squared = value * value;
        const x4 = x_squared * x_squared;
        const x3 = x_squared * value;
        const x8 = x4 * x4;
        const dx2 = @mulAdd(f64, value, value, -x_squared);
        const dx3 = @mulAdd(f64, x_squared, value, -x3) + dx2 * value;
        const p = (atanhd_c_small[0] + x_squared * atanhd_c_small[1]) +
            x4 * (atanhd_c_small[2] + x_squared * atanhd_c_small[3]) +
            x8 * ((atanhd_c_small[4] + x_squared * atanhd_c_small[5]) +
                x4 * (atanhd_c_small[6] + x_squared * atanhd_c_small[7]) + x8 * atanhd_c_small[8]);
        const t = @mulAdd(f64, x_squared, p, 0x1.5555555555555p-56);
        var pl: f64 = undefined;
        var ph = fastTwoSum(0x1.5555555555555p-2, t, &pl);
        ph = mulDoubleDouble(ph, pl, x3, dx3, &pl);
        var tail_low: f64 = undefined;
        ph = fastTwoSum(value, ph, &tail_low);
        pl += tail_low;
        const eps = value * (x4 * 0x1.dp-53 + 0x1p-103);
        const lower = (ph + (pl - eps));
        const upper = (ph + (pl + eps));
        if (lower == upper) return lower;
        return asAtanhZero(value);
    }

    var pl: f64 = undefined;
    const absolute: f64 = @bitCast(absolute_bits);
    const ph = fastTwoSum(1.0, absolute, &pl);
    var ql: f64 = undefined;
    const qh = fastTwoSub(1.0, absolute, &ql);
    const iqh = 1.0 / qh;
    const th = ph * iqh;
    const reduction_low = @mulAdd(f64, ph, iqh, -th) + (pl + ph * (@mulAdd(f64, -qh, iqh, 1.0) - ql * iqh)) * iqh;

    const th_bits: u64 = @bitCast(th);
    const exponent_bits = th_bits >> 52;
    const exponent: i64 = @as(i64, @intCast(exponent_bits)) - 0x3ff;
    const normalized_th_bits = th_bits & (~@as(u64, 0) >> 12);
    const table_i: u64 = normalized_th_bits >> (52 - 5);
    const d = normalized_th_bits & (~@as(u64, 0) >> 17);
    const b_entry = atanhd_B[@intCast(table_i)];
    var j_numerator = normalized_th_bits + (@as(u64, @intCast(b_entry[0])) << 33);
    const d_product = d >> 16;
    if (b_entry[1] >= 0) j_numerator +%= @as(u64, @intCast(b_entry[1])) * d_product else j_numerator -%= @as(u64, @intCast(-b_entry[1])) * d_product;
    const j: u64 = j_numerator >> (52 - 10);
    const idx1_main: usize = @intCast(j >> 5);
    const idx2_main: usize = @intCast(j & 0x1f);
    const normalized_main: f64 = @bitCast(normalized_th_bits | (@as(u64, 0x3ff) << 52));
    const r = (0.5 * atanhd_r1[idx1_main]) * atanhd_r2[idx2_main];
    const dx = @mulAdd(f64, r, normalized_main, -0.5);
    const dx2 = dx * dx;
    const rx = r * normalized_main;
    const dx_low = @mulAdd(f64, r, normalized_main, -rx);
    const f = dx2 * ((atanhd_main_c[0] + dx * atanhd_main_c[1]) +
        dx2 * (atanhd_main_c[2] + dx * atanhd_main_c[3] + dx2 * atanhd_main_c[4]));
    const l2h = 0x1.62e42fefa3ap-2;
    const l2l = -0x1.0ca86c3898dp-50;
    var lh = (atanhd_l1[idx1_main][1] + atanhd_l2[idx2_main][1]) + l2h * @as(f64, @floatFromInt(exponent));
    var ll: f64 = undefined;
    lh = fastTwoSum(lh, rx - 0.5, &ll);
    ll += l2l * @as(f64, @floatFromInt(exponent)) + (atanhd_l1[idx1_main][0] + atanhd_l2[idx2_main][0]) + dx_low + 0.5 * reduction_low / th;
    ll += f;
    const value_sign: f64 = if ((bits & sign_mask_d) == 0) 1.0 else -1.0;
    lh *= value_sign;
    ll *= value_sign;
    const eps = 38e-24 + dx2 * 0x1p-49;
    const lower = lh + (ll - eps);
    const upper = lh + (ll + eps);
    if (lower == upper) return lower;
    var refined_low: f64 = undefined;
    const refined_high = fastTwoSum(th, reduction_low, &refined_low);
    return asAtanhRefine(value, refined_high, refined_low, 0x1.71547652b82fep+1 * @abs(lh + ll));
}

/// Return correctly rounded binary16 atanh with caller-selected nearest-even mode.
pub fn atanhh(value: f16) f16 {
    @setFloatMode(.strict);
    const bits: u16 = @bitCast(value);
    const absolute_bits = bits & magnitude_mask_h;
    if (absolute_bits > infinity_bits_h) return @bitCast(bits | quiet_nan_bit_h);
    if (absolute_bits == infinity_bits_h) return @bitCast(@as(u16, 0x7e00));
    if (absolute_bits == 0x3c00) return copySignH(@bitCast(infinity_bits_h), bits);
    const widened: f32 = @floatCast(value);
    return @floatCast(atanhf(widened));
}
