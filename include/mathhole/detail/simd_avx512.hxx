// mathhole_avx512.hxx

#ifndef MATHHOLE_AVX512_HXX
#define MATHHOLE_AVX512_HXX

#include <cstdint> // std::int32_t
#include <immintrin.h>

namespace mathhole {
	namespace simd {

		//----------------------------------------------------------------------------//
		// Type aliases

		using v512f = __m512;
		using v512d = __m512d;

		//----------------------------------------------------------------------------//
		// Declarations

		//
		// Generic

		template <class Real>
		inline std::int32_t avx512_length();

		template <class Real>
		inline auto avx512_set_scalar(Real a);

		template <class Real>
		inline auto avx512_set_zero();

		//
		// Single precision

		inline v512f avx512_broadcast(float const* addr);
		inline v512f avx512_load(float const* addr);
		inline void avx512_store(float* addr, v512f a);
		inline v512f avx512_add(v512f a, v512f b);
		inline v512f avx512_sub(v512f a, v512f b);
		inline v512f avx512_mul(v512f a, v512f b);
		inline v512f avx512_div(v512f a, v512f b);
		inline v512f avx512_fma(v512f a, v512f b, v512f c);
		inline float avx512_reduce(v512f a);

		//
		// Double precision

		inline v512d avx512_broadcast(double const* addr);
		inline v512d avx512_load(double const* addr);
		inline void avx512_store(double* addr, v512d a);
		inline v512d avx512_add(v512d a, v512d b);
		inline v512d avx512_sub(v512d a, v512d b);
		inline v512d avx512_mul(v512d a, v512d b);
		inline v512d avx512_div(v512d a, v512d b);
		inline v512d avx512_fma(v512d a, v512d b, v512d c);
		inline double avx512_reduce(v512d a);

		//----------------------------------------------------------------------------//
		// Definitions

		//
		// Generic

		template <>
		constexpr std::int32_t avx512_length<float>() {
			return 16;
		}

		template <>
		constexpr std::int32_t avx512_length<double>() {
			return 8;
		}

		template <>
		inline auto avx512_set_scalar(float a) {
			return _mm512_set1_ps(a);
		}

		template <>
		inline auto avx512_set_scalar(double a) {
			return _mm512_set1_pd(a);
		}

		template <>
		inline auto avx512_set_zero<float>() {
			return _mm512_setzero_ps();
		}

		template <>
		inline auto avx512_set_zero<double>() {
			return _mm512_setzero_pd();
		}

		//
		// Single precision

		inline v512f avx512_broadcast(float const* addr) {
			__m128 a = _mm_broadcast_ss(addr);
			return _mm512_broadcastss_ps(a);
		}

		inline v512f avx512_load(float const* addr) {
			return _mm512_loadu_ps(addr);
		}

		inline void avx512_store(float* addr, v512f a) {
			_mm512_storeu_ps(addr, a);
		}

		inline v512f avx512_add(v512f a, v512f b) {
			return _mm512_add_ps(a, b);
		}

		inline v512f avx512_sub(v512f a, v512f b) {
			return _mm512_sub_ps(a, b);
		}

		inline v512f avx512_mul(v512f a, v512f b) {
			return _mm512_mul_ps(a, b);
		}

		inline v512f avx512_div(v512f a, v512f b) {
			return _mm512_div_ps(a, b);
		}

		inline v512f avx512_fma(v512f a, v512f b, v512f c) {
			return _mm512_fmadd_ps(a, b, c);
		}

		inline float avx512_reduce(v512f a) {
			return _mm512_reduce_add_ps(a);
		}

		//
		// Double precision

		inline v512d avx512_set_zero() {
			return _mm512_setzero_pd();
		}

		inline v512d avx512_broadcast(double const* addr) {
			__m128d a = _mm_load1_pd(addr);
			return _mm512_broadcastsd_pd(a);
		}

		inline v512d avx512_load(double const* addr) {
			return _mm512_loadu_pd(addr);
		}

		inline void avx512_store(double* addr, v512d a) {
			_mm512_storeu_pd(addr, a);
		}

		inline v512d avx512_add(v512d a, v512d b) {
			return _mm512_add_pd(a, b);
		}

		inline v512d avx512_sub(v512d a, v512d b) {
			return _mm512_sub_pd(a, b);
		}

		inline v512d avx512_mul(v512d a, v512d b) {
			return _mm512_mul_pd(a, b);
		}

		inline v512d avx512_div(v512d a, v512d b) {
			return _mm512_div_pd(a, b);
		}

		inline v512d avx512_fma(v512d a, v512d b, v512d c) {
			return _mm512_fmadd_pd(a, b, c);
		}

		inline double avx512_reduce(v512d a) {
			return _mm512_reduce_add_pd(a);
		}

	} // namespace simd
} // namespace mathhole

#endif