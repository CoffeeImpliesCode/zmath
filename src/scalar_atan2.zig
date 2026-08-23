// Correctly rounded atan2 for binary16, binary32, and binary64.
//
// The finite kernels are ported from CORE-MATH at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa:
//   src/binary16/atan2/atan2f16.c (Copyright (c) 2022-2025 Alexei Sibidanov and Paul Zimmermann).
//   src/binary32/atan2/atan2f.c (Copyright (c) 2022-2025 Alexei Sibidanov).
//   src/binary64/atan2/atan2.c (Copyright (c) 2022-2025 Alexei Sibidanov).
//
// This code is derived from CORE-MATH and is used under the MIT License:
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of the Software and associated documentation files (the "Software"), to deal
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

const atan2h_impl = @import("scalar_atan2h.zig");
const atan2f_impl = @import("scalar_atan2f.zig");
const atan2d_impl = @import("scalar_atan2d.zig");

/// Return correctly rounded binary16 atan2 with caller-selected nearest-even mode.
pub const atan2h = atan2h_impl.atan2h;

/// Return correctly rounded binary32 atan2 with caller-selected nearest-even mode.
pub const atan2f = atan2f_impl.atan2f;

/// Return correctly rounded binary64 atan2 with caller-selected nearest-even mode.
pub const atan2d = atan2d_impl.atan2d;
