// Copyright (c) 2022 Parisa Khaleghi
// All rights reserved

#ifndef MATHHOLE_STAT_IMPL_HXX
#define MATHHOLE_STAT_IMPL_HXX

#include <cmath>   // std::sqrt
#include <cstdint> // std::int32_t
#include <mathhole/detail/stat_serial_impl.hxx>

namespace mathhole {
	namespace stat {
        //----------------------------------------------------------------------------//
		// Declarations

		template <class Real>
		inline Real arithmetic_mean(std::int32_t n, Real const x[]);

		template <class Real>
		inline Real variance(std::int32_t n, Real const x[]);

		template <class Real>
		inline Real standard_deviation(std::int32_t n, Real const x[]);

		template <class Real>
		inline Real covariance(std::int32_t n, Real const x[], Real const y[]);

		template <class Real>
		inline Real correlation(std::int32_t n, Real const x[], Real const y[]);

        //----------------------------------------------------------------------------//
		// Definitions

		//
		// Arithmetic Mean

		template <class Real>
		inline Real arithmetic_mean(std::int32_t n, Real const x[]) {
			return detail::arithmetic_mean_serial_impl(n, x);
		}

		//
		// Variance

		template <class Real>
		inline Real variance(std::int32_t n, Real const x[]) {
			return detail::variance_serial_impl(n, x);
		}

		//
		// Standard Deviation

		template <class Real>
		inline Real standard_deviation(std::int32_t n, Real const x[]) {
			return detail::standard_deviation_serial_impl(n, x);
		}

		//
		// Covariance

		template <class Real>
		inline Real covariance(std::int32_t n, Real const x[], Real const y[]) {
			return detail::covariance_serial_impl(n, x, y);
		}

		//
		// Correlation

		template <class Real>
		inline Real correlation(std::int32_t n, Real const x[], Real const y[]) {
			return detail::correlation_serial_impl(n, x, y);
		}
    }
}