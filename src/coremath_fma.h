/* Private fused multiply-add helper for the zmath C kernels.

   The kernels are compiled with -ffp-contract=off, so a fused multiply-add has
   to be requested explicitly. On a CPU that evaluates a fused multiply-add in
   one instruction the request is that instruction. On every other target the
   request goes to zmath_core_fmaf_soft for binary32 or zmath_core_fma_soft for
   binary64, which src/scalar_fma.zig exports and implements with the integer
   dyadic engine: exact product, exact sum with a sticky bit and a borrow, one
   round-to-nearest-even step. The library call that a C compiler emits for a
   target without the instruction is not correctly rounded for either width,
   which is what the software path replaces.

   The two widths do not fuse in the same set of instructions, so build.zig
   supplies ZMATH_FMA_NATIVE and ZMATH_FMA64_NATIVE from the target CPU policy
   that src/scalar_fma.zig uses. A core can therefore keep a fused binary32
   operation while its binary64 operations take the exact path.

   This header defines no libc symbol. A translation unit that wants a fused
   multiply-add calls zmath_fmaf or zmath_fma, so one project keeps one
   definition. */

#ifndef ZMATH_COREMATH_FMA_H
#define ZMATH_COREMATH_FMA_H

#ifndef ZMATH_FMA_NATIVE
#error "build.zig must supply the shared ZMATH_FMA_NATIVE target policy"
#endif

#ifndef ZMATH_FMA64_NATIVE
#error "build.zig must supply the shared ZMATH_FMA64_NATIVE target policy"
#endif

#if !ZMATH_FMA_NATIVE

/* Exported by src/scalar_fma.zig. The root module has to reference that file
   so that the symbol reaches every consumer of the C kernels. */
float zmath_core_fmaf_soft(float a, float b, float c);

#endif

#if !ZMATH_FMA64_NATIVE

double zmath_core_fma_soft(double a, double b, double c);

#endif

/* Correctly rounded fused multiply-add of three binary32 values, with a single
   rounding. The caller keeps IEEE round-to-nearest, ties-to-even active. */
static inline float zmath_fmaf(float a, float b, float c) {
#if ZMATH_FMA_NATIVE
  return __builtin_fmaf(a, b, c);
#else
  return zmath_core_fmaf_soft(a, b, c);
#endif
}

/* Correctly rounded fused multiply-add of three binary64 values, with a single
   rounding. The caller keeps IEEE round-to-nearest, ties-to-even active. */
static inline double zmath_fma(double a, double b, double c) {
#if ZMATH_FMA64_NATIVE
  return __builtin_fma(a, b, c);
#else
  return zmath_core_fma_soft(a, b, c);
#endif
}

#endif

