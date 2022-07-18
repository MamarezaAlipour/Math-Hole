// mathhole_avx2.hxx

#ifndef MATHHOLE_AVX2_HXX
#define MATHHOLE_AVX2_HXX

#include <cstdint> // std::int32_t
#include <immintrin.h>

namespace mathhole {
	namespace simd {

		//----------------------------------------------------------------------------//
		// Type aliases

		using v256f = __m256;
		using v256d = __m256d;

		//----------------------------------------------------------------------------//
		// Declarations

		//
		// Generic

		template <class Real>
		inline std::int32_t avx2_length();

		template <class Real>
		inline auto avx2_set_zero();

		template <class Real>
		inline auto avx2_set_scalar(Real a);

		//
		// Single precision

		// inline auto avx2_set_scalar(float a);
		inline v256f avx2_broadcast(float const* addr);
		inline v256f avx2_load(float const* addr);
		inline void avx2_store(float* addr, v256f a);
		inline v256f avx2_add(v256f a, v256f b);
		inline v256f avx2_sub(v256f a, v256f b);
		inline v256f avx2_mul(v256f a, v256f b);
		inline v256f avx2_div(v256f a, v256f b);
		inline v256f avx2_fma(v256f a, v256f b, v256f c);
		inline float avx2_reduce(v256f a);

		//
		// Double precision

		// inline auto avx2_set_scalar(double a);
		inline v256d avx2_broadcast(double const* addr);
		inline v256d avx2_load(double const* addr);
		inline void avx2_store(double* addr, v256d a);
		inline v256d avx2_add(v256d a, v256d b);
		inline v256d avx2_sub(v256d a, v256d b);
		inline v256d avx2_mul(v256d a, v256d b);
		inline v256d avx2_div(v256d a, v256d b);
		inline v256d avx2_fma(v256d a, v256d b, v256d c);
		inline double avx2_reduce(v256d a);

		//----------------------------------------------------------------------------//
		// Definitions

		//
		// Generic

		template <>
		constexpr std::int32_t avx2_length<float>() {
			return 8;
		}

		template <>
		constexpr std::int32_t avx2_length<double>() {
			return 4;
		}

		template <>
		inline auto avx2_set_zero<float>() {
			return _mm256_setzero_ps();
		}

		template <>
		inline auto avx2_set_zero<double>() {
			return _mm256_setzero_pd();
		}

		template <>
		inline auto avx2_set_scalar(float a) {
			return _mm256_set1_ps(a);
		}

		template <>
		inline auto avx2_set_scalar(double a) {
			return _mm256_set1_pd(a);
		}

		//
		// Single precision

		inline v256f avx2_set_scalar(float a) {
			return _mm256_set1_ps(a);
		}

		inline v256f avx2_broadcast(float const* addr) {
			return _mm256_broadcast_ss(addr);
		}

		inline v256f avx2_load(float const* addr) {
			return _mm256_loadu_ps(addr);
		}

		inline void avx2_store(float* addr, v256f a) {
			_mm256_storeu_ps(addr, a);
		}

		inline v256f avx2_add(v256f a, v256f b) {
			return _mm256_add_ps(a, b);
		}

		inline v256f avx2_sub(v256f a, v256f b) {
			return _mm256_sub_ps(a, b);
		}

		inline v256d avx2_set_scalar(double a) {
			return _mm256_set1_pd(a);
		}

		inline v256d avx2_set_zero() {
			return _mm256_setzero_pd();
		}

		inline v256d avx2_broadcast(double const* addr) {
			return _mm256_broadcast_sd(addr);
		}

		inline v256d avx2_load(double const* addr) {
			return _mm256_loadu_pd(addr);
		}

		inline void avx2_store(double* addr, v256d a) {
			_mm256_storeu_pd(addr, a);
		}

		inline v256d avx2_add(v256d a, v256d b) {
			return _mm256_add_pd(a, b);
		}

		inline v256d avx2_sub(v256d a, v256d b) {
			return _mm256_sub_pd(a, b);
		}

		inline v256d avx2_mul(v256d a, v256d b) {
			return _mm256_mul_pd(a, b);
		}

		inline v256d avx2_div(v256d a, v256d b) {
			return _mm256_div_pd(a, b);
		}

		inline v256d avx2_fma(v256d a, v256d b, v256d c) {
			return _mm256_fmadd_pd(a, b, c);
		}

		inline double avx2_reduce(v256d a) {
			__m128d low128 = _mm256_castpd256_pd128(a);
			__m128d high128 = _mm256_extractf128_pd(a, 1);
			low128 = _mm_add_pd(low128, high128);

			__m128d high64 = _mm_unpackhi_pd(low128, low128);
			return _mm_cvtsd_f64(_mm_add_sd(low128, high64));
		}

	} // namespace simd
} // namespace mathhole

#endif