// Correctly-rounded exponential function biased by one for binary64 values.
//
// Copyright (c) 2024 Alexei Sibidanov.
//
// This file is part of the CORE-MATH project
// (https://core-math.gitlabpages.inria.fr/).
//
// Ported from src/binary64/expm1/expm1.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa.
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

const fma_impl = @import("scalar_fma.zig");

const sign_mask: u64 = 0x8000_0000_0000_0000;
const magnitude_mask: u64 = 0x7fff_ffff_ffff_ffff;
const infinity_bits: u64 = 0x7ff0_0000_0000_0000;
const quiet_nan_bit: u64 = 0x0008_0000_0000_0000;

const DoubleDouble = struct {
    low: f64,
    high: f64,
};

const database: [38]f64 = .{
    0x1.e923c188ea79bp-4,  0x1.1a0408712e00ap-2,  0x1.1c38132777b26p-2,  0x1.27f4980d511ffp-2,
    0x1.8172a0e02f90ep-2,  0x1.8bbe2fb45c151p-2,  0x1.bcab27d05abdep-2,  0x1.005ae04256babp-1,
    0x1.accfbe46b4efp-1,   0x1.d086543694c5ap-1,  0x1.273c188aa7b14p+2,  0x1.83d4bcdebb3f4p+2,
    0x1.08f51434652c3p+4,  0x1.1d5c2daebe367p+4,  0x1.c44ce0d716a1ap+4,  0x1.2ee70220fb1c5p+5,
    0x1.89d56a0c38e6fp+5,  0x1.7a60ee15e3e9dp+6,  0x1.1f0da93354198p+7,  0x1.54cd1fea7663ap+7,
    0x1.556c678d5e976p+7,  0x1.2da9e5e6af0bp+8,   0x1.9e7b643238a14p+8,  0x1.d6479eba7c971p+8,
    0x1.0bc04af1b09f5p+9,  -0x1.ab86cb1743b75p-4, -0x1.119aae6072d39p-2, -0x1.175693a03b59p-2,
    -0x1.474d4de7c14bbp-2, -0x1.789d025948efap-2, -0x1.82b5dfaf59b4cp-2, -0x1.9d871e078ebcep-2,
    -0x1.1397add4538acp-1, -0x1.22e24fa3d5cf9p-1, -0x1.dc2b5df1f7d3dp-1, -0x1.0a54d87783d6fp+0,
    -0x1.2a9cad9998262p+0, -0x1.e42a2abb1bf0fp+0,
};

const table_0: [64]DoubleDouble = .{
    .{ .low = 0x0p+0, .high = 0x1p+0 },
    .{ .low = -0x1.19083535b085ep-56, .high = 0x1.02c9a3e778061p+0 },
    .{ .low = 0x1.d73e2a475b466p-55, .high = 0x1.059b0d3158574p+0 },
    .{ .low = 0x1.186be4bb285p-57, .high = 0x1.0874518759bc8p+0 },
    .{ .low = 0x1.8a62e4adc610ap-54, .high = 0x1.0b5586cf9890fp+0 },
    .{ .low = 0x1.03a1727c57b52p-59, .high = 0x1.0e3ec32d3d1a2p+0 },
    .{ .low = -0x1.6c51039449b3ap-54, .high = 0x1.11301d0125b51p+0 },
    .{ .low = -0x1.32fbf9af1369ep-54, .high = 0x1.1429aaea92dep+0 },
    .{ .low = -0x1.19041b9d78a76p-55, .high = 0x1.172b83c7d517bp+0 },
    .{ .low = 0x1.e5b4c7b4968e4p-55, .high = 0x1.1a35beb6fcb75p+0 },
    .{ .low = 0x1.e016e00a2643cp-54, .high = 0x1.1d4873168b9aap+0 },
    .{ .low = 0x1.dc775814a8494p-55, .high = 0x1.2063b88628cd6p+0 },
    .{ .low = 0x1.9b07eb6c70572p-54, .high = 0x1.2387a6e756238p+0 },
    .{ .low = 0x1.2bd339940e9dap-55, .high = 0x1.26b4565e27cddp+0 },
    .{ .low = 0x1.612e8afad1256p-55, .high = 0x1.29e9df51fdee1p+0 },
    .{ .low = 0x1.0024754db41d4p-54, .high = 0x1.2d285a6e4030bp+0 },
    .{ .low = 0x1.6f46ad23182e4p-55, .high = 0x1.306fe0a31b715p+0 },
    .{ .low = 0x1.32721843659a6p-54, .high = 0x1.33c08b26416ffp+0 },
    .{ .low = -0x1.63aeabf42eae2p-54, .high = 0x1.371a7373aa9cbp+0 },
    .{ .low = -0x1.5e436d661f5e2p-56, .high = 0x1.3a7db34e59ff7p+0 },
    .{ .low = 0x1.ada0911f09ebcp-55, .high = 0x1.3dea64c123422p+0 },
    .{ .low = -0x1.ef3691c309278p-58, .high = 0x1.4160a21f72e2ap+0 },
    .{ .low = 0x1.89b7a04ef80dp-59, .high = 0x1.44e086061892dp+0 },
    .{ .low = 0x1.3c1a3b69062fp-56, .high = 0x1.486a2b5c13cdp+0 },
    .{ .low = 0x1.d4397afec42e2p-56, .high = 0x1.4bfdad5362a27p+0 },
    .{ .low = -0x1.4b309d25957e4p-54, .high = 0x1.4f9b2769d2ca7p+0 },
    .{ .low = -0x1.07abe1db13cacp-55, .high = 0x1.5342b569d4f82p+0 },
    .{ .low = 0x1.9bb2c011d93acp-54, .high = 0x1.56f4736b527dap+0 },
    .{ .low = 0x1.6324c054647acp-54, .high = 0x1.5ab07dd485429p+0 },
    .{ .low = 0x1.ba6f93080e65ep-54, .high = 0x1.5e76f15ad2148p+0 },
    .{ .low = -0x1.383c17e40b496p-54, .high = 0x1.6247eb03a5585p+0 },
    .{ .low = -0x1.bb60987591c34p-54, .high = 0x1.6623882552225p+0 },
    .{ .low = -0x1.bdd3413b26456p-54, .high = 0x1.6a09e667f3bcdp+0 },
    .{ .low = -0x1.bbe3a683c88aap-57, .high = 0x1.6dfb23c651a2fp+0 },
    .{ .low = -0x1.16e4786887a9ap-55, .high = 0x1.71f75e8ec5f74p+0 },
    .{ .low = -0x1.0245957316dd4p-54, .high = 0x1.75feb564267c9p+0 },
    .{ .low = -0x1.41577ee04993p-55, .high = 0x1.7a11473eb0187p+0 },
    .{ .low = 0x1.05d02ba15797ep-56, .high = 0x1.7e2f336cf4e62p+0 },
    .{ .low = -0x1.d4c1dd41532d8p-54, .high = 0x1.82589994cce13p+0 },
    .{ .low = -0x1.fc6f89bd4f6bap-54, .high = 0x1.868d99b4492edp+0 },
    .{ .low = 0x1.6e9f156864b26p-54, .high = 0x1.8ace5422aa0dbp+0 },
    .{ .low = 0x1.5cc13a2e3976cp-55, .high = 0x1.8f1ae99157736p+0 },
    .{ .low = -0x1.75fc781b57ebcp-57, .high = 0x1.93737b0cdc5e5p+0 },
    .{ .low = -0x1.d185b7c1b85dp-54, .high = 0x1.97d829fde4e5p+0 },
    .{ .low = 0x1.c7c46b071f2bep-56, .high = 0x1.9c49182a3f09p+0 },
    .{ .low = -0x1.359495d1cd532p-54, .high = 0x1.a0c667b5de565p+0 },
    .{ .low = -0x1.d2f6edb8d41e2p-54, .high = 0x1.a5503b23e255dp+0 },
    .{ .low = 0x1.0fac90ef7fd32p-54, .high = 0x1.a9e6b5579fdbfp+0 },
    .{ .low = 0x1.7a1cd345dcc82p-54, .high = 0x1.ae89f995ad3adp+0 },
    .{ .low = -0x1.2805e3084d708p-57, .high = 0x1.b33a2b84f15fbp+0 },
    .{ .low = -0x1.5584f7e54ac3ap-56, .high = 0x1.b7f76f2fb5e47p+0 },
    .{ .low = 0x1.23dd07a2d9e84p-55, .high = 0x1.bcc1e904bc1d2p+0 },
    .{ .low = 0x1.11065895048dep-55, .high = 0x1.c199bdd85529cp+0 },
    .{ .low = 0x1.2884dff483cacp-54, .high = 0x1.c67f12e57d14bp+0 },
    .{ .low = 0x1.503cbd1e949dcp-56, .high = 0x1.cb720dcef9069p+0 },
    .{ .low = -0x1.cbc3743797a9cp-54, .high = 0x1.d072d4a07897cp+0 },
    .{ .low = 0x1.2ed02d75b3706p-55, .high = 0x1.d5818dcfba487p+0 },
    .{ .low = 0x1.c2300696db532p-54, .high = 0x1.da9e603db3285p+0 },
    .{ .low = -0x1.1a5cd4f184b5cp-54, .high = 0x1.dfc97337b9b5fp+0 },
    .{ .low = 0x1.39e8980a9cc9p-55, .high = 0x1.e502ee78b3ff6p+0 },
    .{ .low = -0x1.e9c23179c2894p-54, .high = 0x1.ea4afa2a490dap+0 },
    .{ .low = 0x1.dc7f486a4b6bp-54, .high = 0x1.efa1bee615a27p+0 },
    .{ .low = 0x1.9d3e12dd8a18ap-54, .high = 0x1.f50765b6e454p+0 },
    .{ .low = 0x1.74853f3a5931ep-55, .high = 0x1.fa7c1819e90d8p+0 },
};

const table_1: [64]DoubleDouble = .{
    .{ .low = 0x0p+0, .high = 0x1p+0 },
    .{ .low = 0x1.ae8e38c59c72ap-54, .high = 0x1.000b175effdc7p+0 },
    .{ .low = -0x1.7b5d0d58ea8f4p-58, .high = 0x1.00162f3904052p+0 },
    .{ .low = 0x1.4115cb6b16a8ep-54, .high = 0x1.0021478e11ce6p+0 },
    .{ .low = -0x1.d7c96f201bb2ep-55, .high = 0x1.002c605e2e8cfp+0 },
    .{ .low = 0x1.84711d4c35eap-54, .high = 0x1.003779a95f959p+0 },
    .{ .low = -0x1.0484245243778p-55, .high = 0x1.0042936faa3d8p+0 },
    .{ .low = -0x1.4b237da2025fap-54, .high = 0x1.004dadb113dap+0 },
    .{ .low = -0x1.5e00e62d6b30ep-56, .high = 0x1.0058c86da1c0ap+0 },
    .{ .low = 0x1.a1d6cedbb948p-54, .high = 0x1.0063e3a559473p+0 },
    .{ .low = -0x1.4acf197a00142p-54, .high = 0x1.006eff583fc3dp+0 },
    .{ .low = -0x1.eaf2ea42391a6p-57, .high = 0x1.007a1b865a8cap+0 },
    .{ .low = 0x1.da93f90835f76p-56, .high = 0x1.0085382faef83p+0 },
    .{ .low = -0x1.6a79084ab093cp-55, .high = 0x1.00905554425d4p+0 },
    .{ .low = 0x1.86364f8fbe8f8p-54, .high = 0x1.009b72f41a12bp+0 },
    .{ .low = -0x1.82e8e14e3110ep-55, .high = 0x1.00a6910f3b6fdp+0 },
    .{ .low = -0x1.4f6b2a7609f72p-55, .high = 0x1.00b1afa5abcbfp+0 },
    .{ .low = -0x1.e1a258ea8f71ap-56, .high = 0x1.00bcceb7707ecp+0 },
    .{ .low = 0x1.4362ca5bc26f2p-56, .high = 0x1.00c7ee448ee02p+0 },
    .{ .low = 0x1.095a56c919d02p-54, .high = 0x1.00d30e4d0c483p+0 },
    .{ .low = -0x1.406ac4e81a646p-57, .high = 0x1.00de2ed0ee0f5p+0 },
    .{ .low = 0x1.b5a6902767e08p-54, .high = 0x1.00e94fd0398ep+0 },
    .{ .low = -0x1.91b206085932p-54, .high = 0x1.00f4714af41d3p+0 },
    .{ .low = 0x1.427068ab22306p-55, .high = 0x1.00ff93412315cp+0 },
    .{ .low = 0x1.c1d0660524e08p-54, .high = 0x1.010ab5b2cbd11p+0 },
    .{ .low = -0x1.e7bdfb3204be8p-54, .high = 0x1.0115d89ff3a8bp+0 },
    .{ .low = 0x1.843aa8b9cbbc6p-55, .high = 0x1.0120fc089ff63p+0 },
    .{ .low = -0x1.34104ee7edae8p-56, .high = 0x1.012c1fecd613bp+0 },
    .{ .low = -0x1.2b6aeb6176892p-56, .high = 0x1.0137444c9b5b5p+0 },
    .{ .low = 0x1.a8cd33b8a1bb2p-56, .high = 0x1.01426927f5278p+0 },
    .{ .low = 0x1.2edc08e5da99ap-56, .high = 0x1.014d8e7ee8d2fp+0 },
    .{ .low = 0x1.57ba2dc7e0c72p-55, .high = 0x1.0158b4517bb88p+0 },
    .{ .low = 0x1.b61299ab8cdb8p-54, .high = 0x1.0163da9fb3335p+0 },
    .{ .low = -0x1.90565902c5f44p-54, .high = 0x1.016f0169949edp+0 },
    .{ .low = 0x1.70fc41c5c2d54p-55, .high = 0x1.017a28af25567p+0 },
    .{ .low = 0x1.4b9a6e145d76cp-54, .high = 0x1.018550706ab62p+0 },
    .{ .low = -0x1.008eff5142bfap-56, .high = 0x1.019078ad6a19fp+0 },
    .{ .low = -0x1.77669f033c7dep-54, .high = 0x1.019ba16628de2p+0 },
    .{ .low = -0x1.09bb78eeead0ap-54, .high = 0x1.01a6ca9aac5f3p+0 },
    .{ .low = 0x1.371231477ece6p-54, .high = 0x1.01b1f44af9f9ep+0 },
    .{ .low = 0x1.5e7626621eb5ap-56, .high = 0x1.01bd1e77170b4p+0 },
    .{ .low = -0x1.bc72b100828a4p-54, .high = 0x1.01c8491f08f08p+0 },
    .{ .low = -0x1.ce39cbbab8bbep-57, .high = 0x1.01d37442d507p+0 },
    .{ .low = 0x1.16996709da2e2p-55, .high = 0x1.01de9fe280ac8p+0 },
    .{ .low = -0x1.c11f5239bf536p-55, .high = 0x1.01e9cbfe113efp+0 },
    .{ .low = 0x1.e1d4eb5edc6b4p-55, .high = 0x1.01f4f8958c1c6p+0 },
    .{ .low = -0x1.afb99946ee3fp-54, .high = 0x1.020025a8f6a35p+0 },
    .{ .low = -0x1.8f06d8a148a32p-54, .high = 0x1.020b533856324p+0 },
    .{ .low = -0x1.2bf310fc54eb6p-55, .high = 0x1.02168143b0281p+0 },
    .{ .low = -0x1.c95a035eb4176p-54, .high = 0x1.0221afcb09e3ep+0 },
    .{ .low = -0x1.491793e46834cp-54, .high = 0x1.022cdece68c4fp+0 },
    .{ .low = -0x1.3e8d0d9c4909p-56, .high = 0x1.02380e4dd22adp+0 },
    .{ .low = -0x1.314aa16278aa4p-54, .high = 0x1.02433e494b755p+0 },
    .{ .low = 0x1.48daf888e965p-55, .high = 0x1.024e6ec0da046p+0 },
    .{ .low = 0x1.56dc8046821f4p-55, .high = 0x1.02599fb483385p+0 },
    .{ .low = 0x1.45b42356b9d46p-54, .high = 0x1.0264d1244c719p+0 },
    .{ .low = -0x1.082ef51b61d7ep-56, .high = 0x1.027003103b10ep+0 },
    .{ .low = 0x1.2106ed0920a34p-56, .high = 0x1.027b357854772p+0 },
    .{ .low = -0x1.fd4cf26ea5d0ep-54, .high = 0x1.0286685c9e059p+0 },
    .{ .low = -0x1.09f8775e78084p-54, .high = 0x1.02919bbd1d1d8p+0 },
    .{ .low = 0x1.64cbba902ca28p-58, .high = 0x1.029ccf99d720ap+0 },
    .{ .low = 0x1.4383ef231d206p-54, .high = 0x1.02a803f2d170dp+0 },
    .{ .low = 0x1.4a47a505b3a46p-54, .high = 0x1.02b338c811703p+0 },
    .{ .low = 0x1.e47120223468p-54, .high = 0x1.02be6e199c811p+0 },
};

const small_table: [65]DoubleDouble = .{
    .{ .low = -0x1.797d4686c5393p-57, .high = -0x1.c5041854df7d4p-3 },
    .{ .low = -0x1.ea1cb9d163339p-55, .high = -0x1.b881a23aebb48p-3 },
    .{ .low = 0x1.f483a3e8cd60fp-55, .high = -0x1.abe60e1f21838p-3 },
    .{ .low = 0x1.dffd920f493dbp-56, .high = -0x1.9f3129931fabp-3 },
    .{ .low = -0x1.51bfdbb129094p-55, .high = -0x1.9262c1c3430ap-3 },
    .{ .low = 0x1.cd3e5225e2206p-55, .high = -0x1.857aa375db4e4p-3 },
    .{ .low = 0x1.e3a6bdaece8f9p-58, .high = -0x1.78789b0a5e0cp-3 },
    .{ .low = -0x1.daf2ae0c2d3d4p-55, .high = -0x1.6b5c7478983d8p-3 },
    .{ .low = -0x1.fd36226fadd44p-56, .high = -0x1.5e25fb4fde21p-3 },
    .{ .low = 0x1.d887cd0341abp-56, .high = -0x1.50d4fab639758p-3 },
    .{ .low = -0x1.676a52a1a618bp-55, .high = -0x1.43693d679612cp-3 },
    .{ .low = 0x1.9776b420ad283p-56, .high = -0x1.35e28db4ecd9cp-3 },
    .{ .low = 0x1.3d5fd7d70a5edp-56, .high = -0x1.2840b5836cf68p-3 },
    .{ .low = 0x1.a94ad2c8fa0bfp-58, .high = -0x1.1a837e4ba376p-3 },
    .{ .low = 0x1.6ad4c353465bp-61, .high = -0x1.0caab118a1278p-3 },
    .{ .low = -0x1.8bba170e59b65p-56, .high = -0x1.fd6c2d0e3d91p-4 },
    .{ .low = -0x1.e1e0a76cb0685p-55, .high = -0x1.e14aed893eefp-4 },
    .{ .low = 0x1.fe131f55e75f8p-55, .high = -0x1.c4f1331d22d4p-4 },
    .{ .low = -0x1.b5beee8bcee31p-55, .high = -0x1.a85e8c62d9c1p-4 },
    .{ .low = -0x1.7fe9b02c25e9bp-56, .high = -0x1.8b92870fa2b58p-4 },
    .{ .low = -0x1.32ae7bdaf1116p-55, .high = -0x1.6e8caff341fe8p-4 },
    .{ .low = 0x1.a6cfe58cbd73bp-56, .high = -0x1.514c92f634788p-4 },
    .{ .low = 0x1.8798de3138a56p-57, .high = -0x1.33d1bb17df2e8p-4 },
    .{ .low = -0x1.589321a7ef10bp-60, .high = -0x1.161bb26cbb59p-4 },
    .{ .low = -0x1.8d0e700fcfb65p-56, .high = -0x1.f0540438fd5cp-5 },
    .{ .low = 0x1.473ef07d5dd3bp-55, .high = -0x1.b3f864c08p-5 },
    .{ .low = -0x1.38e62149c16e2p-55, .high = -0x1.77239501304p-5 },
    .{ .low = -0x1.08bb6309bd394p-58, .high = -0x1.39d4a1a77e05p-5 },
    .{ .low = -0x1.bad3fd501a227p-55, .high = -0x1.f8152aee945p-6 },
    .{ .low = 0x1.3d27ac39ed253p-57, .high = -0x1.7b88f290230ep-6 },
    .{ .low = -0x1.b60bbd08aac55p-55, .high = -0x1.fc055004416cp-7 },
    .{ .low = -0x1.a00d03b3359dep-59, .high = -0x1.fe0154aaeed8p-8 },
    .{ .low = 0x0p+0, .high = 0x0p+0 },
    .{ .low = 0x1.861931c15e39bp-55, .high = 0x1.0100ab00222cp-7 },
    .{ .low = 0x1.7ab864b3e9045p-56, .high = 0x1.0202ad5778e4p-6 },
    .{ .low = 0x1.4e5659d75e95bp-56, .high = 0x1.84890d904374p-6 },
    .{ .low = 0x1.8e0bd083aba81p-56, .high = 0x1.040ac0224fd9p-5 },
    .{ .low = 0x1.45cc1cf959b1bp-60, .high = 0x1.465509d383ebp-5 },
    .{ .low = -0x1.eb6980ce14da7p-55, .high = 0x1.89246d053d18p-5 },
    .{ .low = 0x1.7324137d6c342p-56, .high = 0x1.cc79f4f5613ap-5 },
    .{ .low = -0x1.5272ff30eed1bp-59, .high = 0x1.082b577d34ed8p-4 },
    .{ .low = -0x1.1280f19dace1cp-55, .high = 0x1.2a5dd543ccc5p-4 },
    .{ .low = -0x1.d550af31c8ec3p-55, .high = 0x1.4cd4fc989cd68p-4 },
    .{ .low = 0x1.7923b72aa582dp-55, .high = 0x1.6f9157587069p-4 },
    .{ .low = -0x1.76c2e732457f1p-56, .high = 0x1.92937074e0cd8p-4 },
    .{ .low = 0x1.81f5c92a5200fp-55, .high = 0x1.b5dbd3f68122p-4 },
    .{ .low = 0x1.e8ac7a4d3206cp-55, .high = 0x1.d96b0eff0e79p-4 },
    .{ .low = -0x1.12db6f4bbe33bp-56, .high = 0x1.fd41afcba45e8p-4 },
    .{ .low = -0x1.8c4a5df1ec7e5p-58, .high = 0x1.10b022db7ae68p-3 },
    .{ .low = -0x1.bd4b1c37ea8a2p-57, .high = 0x1.22e3b09dc54d8p-3 },
    .{ .low = 0x1.5aeb9860044dp-55, .high = 0x1.353bc9fb00b2p-3 },
    .{ .low = -0x1.4c26602c63fdap-57, .high = 0x1.47b8b853aafecp-3 },
    .{ .low = -0x1.7f644c1f9d314p-55, .high = 0x1.5a5ac59b963ccp-3 },
    .{ .low = 0x1.f5aa8ec61fc2dp-55, .high = 0x1.6d223c5b10638p-3 },
    .{ .low = 0x1.7ab912c69ffebp-61, .high = 0x1.800f67b00d7b8p-3 },
    .{ .low = -0x1.b3564bc0ec9cdp-58, .high = 0x1.9322934f54148p-3 },
    .{ .low = 0x1.6a7062465be33p-55, .high = 0x1.a65c0b85ac1a8p-3 },
    .{ .low = -0x1.85718d2ff1bf4p-55, .high = 0x1.b9bc1d3910094p-3 },
    .{ .low = -0x1.045cb0c685e08p-55, .high = 0x1.cd4315e9e0834p-3 },
    .{ .low = -0x1.6e7fb859d5055p-62, .high = 0x1.e0f143b41a554p-3 },
    .{ .low = 0x1.51bbdee020603p-55, .high = 0x1.f4c6f5508ee5cp-3 },
    .{ .low = 0x1.e17611afc42c5p-57, .high = 0x1.04623d0b0f8c8p-2 },
    .{ .low = -0x1.1c5b2e8735a43p-56, .high = 0x1.0e7510fd7c564p-2 },
    .{ .low = -0x1.25fe139c4cffdp-55, .high = 0x1.189c1ecaeb084p-2 },
    .{ .low = -0x1.89843c4964554p-56, .high = 0x1.22d78f0fa061ap-2 },
};

fn roundEvenFinite(value: f64) i64 {
    @setFloatMode(.strict);

    // All callers bound |value| below 2^23. Under RNE, this fixed-exponent
    // bias rounds to nearest-even and its bit delta is the integer result.
    const rounding_constant: f64 = 0x1.8p52;
    const biased = value + rounding_constant;
    const biased_bits: u64 = @bitCast(biased);
    const rounding_bits: u64 = @bitCast(rounding_constant);
    return @bitCast(biased_bits -% rounding_bits);
}

fn fastTwoSum(x: f64, y: f64) DoubleDouble {
    const sum = x + y;
    const rounded_x = sum - x;
    return .{ .low = y - rounded_x, .high = sum };
}

fn fastSum(x_high: f64, x_low: f64, y_high: f64, y_low: f64) DoubleDouble {
    const sum = fastTwoSum(x_high, y_high);
    return .{ .low = (x_low + y_low) + sum.low, .high = sum.high };
}

fn multiplyDoubleDouble(x_high: f64, x_low: f64, c_high: f64, c_low: f64) DoubleDouble {
    const high = c_high * x_high;
    const low = (c_high * x_low + c_low * x_high) + fma_impl.fma64(c_high, x_high, -high);
    return .{ .low = low, .high = high };
}

fn multiplyDoubleWord(x_high: f64, x_low: f64, c: f64) DoubleDouble {
    const high = c * x_high;
    const low = c * x_low + fma_impl.fma64(c, x_high, -high);
    return .{ .low = low, .high = high };
}

fn evaluateDoubleDoublePolynomial(x_high: f64, x_low: f64, coefficients: []const DoubleDouble) DoubleDouble {
    var index = coefficients.len - 1;
    var result = coefficients[index];
    while (index != 0) {
        index -= 1;
        const coefficient = coefficients[index];
        result = multiplyDoubleDouble(x_high, x_low, result.high, result.low);
        const sum = fastSum(coefficient.high, coefficient.low, result.high, result.low);
        result = .{ .low = sum.low, .high = sum.high };
    }
    return result;
}

fn evaluateDoubleWordPolynomial(x: f64, coefficients: []const DoubleDouble, low: *f64) f64 {
    var index = coefficients.len - 1;
    const initial = fastTwoSum(coefficients[index].high, low.*);
    var high = initial.high;
    var current_low = coefficients[index].low + initial.low;
    while (index != 0) {
        index -= 1;
        const product = multiplyDoubleWord(high, current_low, x);
        const sum = fastSum(coefficients[index].high, coefficients[index].low, product.high, product.low);
        high = sum.high;
        current_low = sum.low;
    }
    low.* = current_low;
    return high;
}

fn scale(value: f64, exponent: i64) f64 {
    const bits: u64 = @bitCast(value);
    const exponent_bits = @as(u64, @bitCast(exponent)) << 52;
    return @bitCast(bits +% exponent_bits);
}

fn toDenormal(value: f64) f64 {
    const bits: u64 = @bitCast(value);
    return @bitCast(bits & 0x000f_ffff_ffff_ffff);
}

fn databaseLookup(value: f64, fallback: f64) f64 {
    const value_bits: u64 = @bitCast(value);
    var lower: isize = 0;
    var upper: isize = @intCast(database.len - 1);
    while (lower <= upper) {
        const middle = @divTrunc(lower + upper, 2);
        const middle_index: usize = @intCast(middle);
        const database_bits: u64 = @bitCast(database[middle_index]);
        if (database_bits < value_bits) {
            lower = middle + 1;
        } else if (database_bits == value_bits) {
            const rounding_patterns = [2]u64{ 0x76f58b0d65bd5553, 0xc06 };
            const sign_pattern: u64 = 0x300e81651c;
            const fallback_bits: u64 = @bitCast(fallback);
            const middle_shift: u6 = @intCast(middle);
            const delta_bits = ((sign_pattern >> middle_shift) << 63) |
                ((((fallback_bits >> 52) & 0x7ff) - 54) << 52);
            const delta: f64 = @bitCast(delta_bits);
            const pattern_index: usize = @intCast(middle >> 5);
            const pattern_shift: u6 = @intCast((middle << 1) & 63);
            const rounding_pattern = (rounding_patterns[pattern_index] >> pattern_shift) & 3;
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

const accurate_small_coefficients: [6]f64 = .{
    0x1.93974a8ca5354p-37,
    0x1.ae7f3e71e4908p-41,
    0x1.ae7f357341648p-45,
    0x1.952c7f96664cbp-49,
    0x1.686f8ce633aaep-53,
    0x1.2f49b2fbfb5b6p-57,
};
const accurate_small_polynomial: [11]DoubleDouble = .{
    .{ .low = 0x1.5555555555554p-57, .high = 0x1.5555555555555p-3 },
    .{ .low = 0x1.5555555555123p-59, .high = 0x1.5555555555555p-5 },
    .{ .low = 0x1.1111111118167p-63, .high = 0x1.1111111111111p-7 },
    .{ .low = -0x1.f49f49e220ceap-65, .high = 0x1.6c16c16c16c17p-10 },
    .{ .low = 0x1.a019eff6f919cp-73, .high = 0x1.a01a01a01a01ap-13 },
    .{ .low = 0x1.9fcff48a75b41p-76, .high = 0x1.a01a01a01a01ap-16 },
    .{ .low = -0x1.c14f73758cd7fp-73, .high = 0x1.71de3a556c734p-19 },
    .{ .low = 0x1.dfce97931018fp-76, .high = 0x1.27e4fb7789f5cp-22 },
    .{ .low = 0x1.c513da9e4c9c5p-80, .high = 0x1.ae64567f544e3p-26 },
    .{ .low = 0x1.ca00af84f2b6p-83, .high = 0x1.1eed8eff8d831p-29 },
    .{ .low = 0x1.f27ac6000898fp-87, .high = 0x1.6124613a86e8fp-33 },
};

noinline fn accurate(value: f64) f64 {
    const value_bits: u64 = @bitCast(value);
    if (((value_bits >> 52) & 0x7ff) < 0x3c9) return 1.0 + value;

    if (@abs(value) < 0.25) {
        var low = value * (accurate_small_coefficients[0] + value *
            (accurate_small_coefficients[1] + value *
                (accurate_small_coefficients[2] + value *
                    (accurate_small_coefficients[3] + value *
                        (accurate_small_coefficients[4] + value * accurate_small_coefficients[5])))));
        var high = evaluateDoubleWordPolynomial(value, &accurate_small_polynomial, &low);
        var product = multiplyDoubleWord(high, low, value);
        high = product.high;
        low = product.low;
        product = multiplyDoubleWord(high, low, value);
        high = product.high;
        low = product.low;
        product = multiplyDoubleWord(high, low, value);
        high = product.high;
        low = product.low;
        const half_x = 0.5 * value;
        const x2_high = value * half_x;
        const x2_low = fma_impl.fma64(value, half_x, -x2_high);
        const summed = fastSum(x2_high, x2_low, high, low);
        high = summed.high;
        low = summed.low;
        var first = fastTwoSum(value, high);
        var second = fastTwoSum(first.low, low);
        first = fastTwoSum(first.high, second.high);
        second = fastTwoSum(first.low, second.low);
        var rounded_bits: u64 = @bitCast(second.high);
        if ((rounded_bits & 0x000f_ffff_ffff_ffff) == 0) {
            if ((rounded_bits << 1) == 0) return databaseLookup(value, first.high);
            const low_bits: u64 = @bitCast(second.low);
            const same_sign = ((rounded_bits ^ low_bits) & sign_mask) == 0;
            const direction: u64 = if (same_sign) 1 else 0xffff_ffff_ffff_ffff;
            rounded_bits +%= direction;
            second.high = @bitCast(rounded_bits);
        }
        return first.high + second.high;
    }

    const coefficients: [7]DoubleDouble = .{
        .{ .low = 0, .high = 1.0 },
        .{ .low = 0x1.712f72ecec2cfp-99, .high = 0x1p-1 },
        .{ .low = 0x1.5555555554d07p-57, .high = 0x1.5555555555555p-3 },
        .{ .low = 0x1.55194d28275dap-59, .high = 0x1.5555555555555p-5 },
        .{ .low = 0x1.12faa0e1c0f7bp-63, .high = 0x1.1111111111111p-7 },
        .{ .low = -0x1.4ba45ab25d2a3p-64, .high = 0x1.6c16c16da6973p-10 },
        .{ .low = -0x1.9091d845ecd36p-67, .high = 0x1.a01a019eb7f31p-13 },
    };
    const reduction_scale = 0x1.71547652b82fep+12;
    const reduction_integer = roundEvenFinite(value * reduction_scale);
    const reduction: f64 = @floatFromInt(reduction_integer);
    const table_0_index: usize = @intCast((reduction_integer >> 6) & 0x3f);
    const table_1_index: usize = @intCast(reduction_integer & 0x3f);
    const exponent = reduction_integer >> 12;
    const first_table = table_0[table_0_index];
    const second_table = table_1[table_1_index];
    const table_product = multiplyDoubleDouble(
        first_table.high,
        first_table.low,
        second_table.high,
        second_table.low,
    );

    const log2_high = 0x1.62e42ffp-13;
    const log2_low = 0x1.718432a1b0e26p-47;
    const log2_tail = 0x1.9ff0342542fc3p-102;
    const reduced = value - log2_high * reduction;
    var reduced_low = log2_low * reduction;
    const reduced_tail = log2_tail * reduction + fma_impl.fma64(log2_low, reduction, -reduced_low);
    const reduced_high = reduced + reduced_low;
    reduced_low = (reduced - reduced_high) + reduced_low + reduced_tail;

    var polynomial = evaluateDoubleDoublePolynomial(reduced_high, reduced_low, &coefficients);
    polynomial = multiplyDoubleDouble(reduced_high, reduced_low, polynomial.high, polynomial.low);
    polynomial = multiplyDoubleDouble(
        polynomial.high,
        polynomial.low,
        table_product.high,
        table_product.low,
    );
    var result = fastSum(
        table_product.high,
        table_product.low,
        polynomial.high,
        polynomial.low,
    );

    const offset_bits = @as(u64, @intCast(2048 + 1023 - exponent)) << 52;
    const offset: f64 = @bitCast(offset_bits);
    var rounding_error: f64 = 0;
    if (exponent < 53) {
        const sum = fastTwoSum(offset, result.high);
        result.high = sum.high;
        rounding_error = sum.low;
    } else if (exponent < 104) {
        const sum = fastTwoSum(result.high, offset);
        result.high = sum.high;
        rounding_error = sum.low;
    }
    result.low += rounding_error;
    result = fastTwoSum(result.high, result.low);
    const low_bits: u64 = @bitCast(result.low);
    const database_probe = (low_bits +% 8) & 0x000f_ffff_ffff_ffff;
    if (database_probe <= 8) result.high = databaseLookup(value, result.high);
    return scale(result.high, exponent);
}

/// Correctly-rounded exponential function biased by one for binary64 values.
///
/// The caller must keep IEEE round-to-nearest, ties-to-even active.
pub fn expm1d(value: f64) f64 {
    @setFloatMode(.strict);
    const input_bits: u64 = @bitCast(value);
    const absolute_bits = input_bits & magnitude_mask;
    if (absolute_bits > infinity_bits) return @bitCast(input_bits | quiet_nan_bit);
    if (absolute_bits == infinity_bits) {
        if ((input_bits & sign_mask) == 0) return value;
        return @bitCast(@as(u64, 0xbff0_0000_0000_0000));
    }

    if (absolute_bits < 0x3fd0_0000_0000_0000) {
        if (absolute_bits < 0x3ca0_0000_0000_0000) {
            return value;
        }
        const scaled = 0x1p7 * value;
        const rounded = roundEvenFinite(scaled);
        const reduced = scaled - @as(f64, @floatFromInt(rounded));
        const reduced_squared = reduced * reduced;
        const table_value = small_table[@intCast(rounded + 32)];
        const coefficients = [6]f64{
            0x1p-7,
            0x1p-15,
            0x1.55555555551adp-24,
            0x1.555555555599cp-33,
            0x1.11111ad1ad69dp-42,
            0x1.6c16c168b1fb5p-52,
        };
        const high = reduced * coefficients[0];
        const low = reduced_squared * ((coefficients[1] + reduced * coefficients[2]) + reduced_squared *
            (coefficients[3] + reduced * (coefficients[4] + reduced * coefficients[5])));
        const epsilon = reduced_squared * 0x1.ap-65 + 0x1p-104;
        var sum = fastTwoSum(table_value.high, high);
        const low_sum = sum.low + table_value.low + low;
        const product = multiplyDoubleDouble(table_value.high, table_value.low, high, low);
        sum = fastSum(sum.high, low_sum, product.high, product.low);
        const upper = sum.high + (sum.low + epsilon);
        const lower = sum.high + (sum.low - epsilon);
        if (upper != lower) return accurate(value);
        return lower;
    }
    if (absolute_bits >= 0x4086_2e42_fefa_39f0) {
        if ((input_bits & sign_mask) == 0) return @bitCast(infinity_bits);
    }
    if (input_bits >= 0xc042_5e4f_7b27_37fa) {
        if (input_bits >= 0xc042_b708_8723_20e2) return @bitCast(@as(u64, 0xbff0_0000_0000_0000));
        return (0x1.25e4f7b2737fap+5 + value + 0x1.8486612173c69p-51) *
            0x1.71547652b82fep-54 - 0x1.fffffffffffffp-1;
    }

    const reduction_scale = 0x1.71547652b82fep+12;
    const reduction_integer = roundEvenFinite(value * reduction_scale);
    const reduction: f64 = @floatFromInt(reduction_integer);
    const table_0_index: usize = @intCast((reduction_integer >> 6) & 0x3f);
    const table_1_index: usize = @intCast(reduction_integer & 0x3f);
    const exponent = reduction_integer >> 12;
    const first = table_0[table_0_index];
    const second = table_1[table_1_index];
    const table_product = multiplyDoubleDouble(first.high, first.low, second.high, second.low);
    const reduced = (value - 0x1.62e42ffp-13 * reduction) + 0x1.718432a1b0e26p-47 * reduction;
    const reduced_squared = reduced * reduced;
    const polynomial = (1.0 + reduced * 0.5) + reduced_squared *
        (0x1.55555557e54ffp-3 + reduced * 0x1.55555553a12f4p-5);
    var high = table_product.high;
    const table_times_reduced = table_product.high * reduced;
    var low = table_product.low + table_times_reduced * polynomial;
    const epsilon = 1.64e-19 * table_product.high;
    const offset_bits = @as(u64, @intCast(2048 + 1023 - exponent)) << 52;
    const offset: f64 = @bitCast(offset_bits);
    var rounding_error: f64 = 0;
    if (exponent < 53) {
        const sum = fastTwoSum(offset, high);
        high = sum.high;
        rounding_error = sum.low;
    } else if (exponent < 75) {
        const sum = fastTwoSum(high, offset);
        high = sum.high;
        rounding_error = sum.low;
    }
    low += rounding_error;
    const final_sum = fastTwoSum(high, low);
    high = final_sum.high;
    low = final_sum.low;
    const upper = high + (low + epsilon);
    const lower = high + (low - epsilon);
    if (upper != lower) return accurate(value);
    return scale(lower, exponent);
}
