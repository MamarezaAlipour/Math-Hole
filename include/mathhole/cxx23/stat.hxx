// Copyright (c) 2022 Geek pww
// All rights reserved

#ifndef MATHHOLE_CXX23_STAT_HXX
#define MATHHOLE_CXX23_STAT_HXX

#include <cassert>
#include <cmath>   // std::sqrt
#include <cstdint> // std::int32_t
#include <experimental/mdspan>
#include <mathhole/detail/stat_impl.hxx>

namespace mathhole {
	namespace stat {

		//----------------------------------------------------------------------------//
		// Declarations

		template <class in_vector_t>
		inline auto arithmetic_mean(in_vector_t x);

		template <class in_vector_t>
		inline auto variance(in_vector_t x);

		template <class in_vector_t>
		inline auto standard_deviation(in_vector_t x);

		template <class in_vector_1_t, class in_vector_2_t>
		inline auto covariance(in_vector_1_t x, in_vector_2_t y);

		template <class in_vector_1_t, class in_vector_2_t>
		inline auto correlation(in_vector_1_t x, in_vector_2_t y);

		//----------------------------------------------------------------------------//
		// Definitions

		//
		// Arithmetic Mean

		template <class in_vector_t>
		inline auto arithmetic_mean(in_vector_t x) {
			static_assert(x.rank() == 1);
			return detail::arithmetic_mean_serial_impl(x.extent(0), x.data());
		}

		//
		// Variance

		template <class in_vector_t>
		inline auto variance(in_vector_t x) {
			static_assert(x.rank() == 1);
			return detail::variance_serial_impl(x.extent(0), x.data());
		}

		//
		// Standard Deviation

		template <class in_vector_t>
		inline auto standard_deviation(in_vector_t x) {
			static_assert(x.rank() == 1);
			return detail::standard_deviation_serial_impl(x.extent(0), x.data());
		}

		//
		// Covariance

		template <class in_vector_1_t, class in_vector_2_t>
		inline auto covariance(in_vector_1_t x, in_vector_2_t y) {
			static_assert(x.rank() == 1 or y.rank() == 1);
			assert(x.extent(0) == y.extent(0));
			return detail::covariance_serial_impl(x.extent(0), x.data(), y.data());
		}

		//
		// Correlation

		template <class in_vector_1_t, class in_vector_2_t>
		inline auto correlation(in_vector_1_t x, in_vector_2_t y) {
			static_assert(x.rank() == 1 or y.rank() == 1);
			assert(x.extent(0) == y.extent(0));
			return detail::correlation_serial_impl(x.extent(0), x.data(), y.data());
		}

	} // namespace stat
} // namespace mathhole

#endif