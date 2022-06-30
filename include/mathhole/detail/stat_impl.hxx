// Copyright (c) 2022 Parisa Khaleghi
// All rights reserved

#ifndef MATHHOLE_STAT_IMPL_HXX
#define MATHHOLE_STAT_IMPL_HXX

#include <cmath>   // std::sqrt
#include <cstdint> // std::int32_t

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
		// Definitions (Soon)

        //
		// Arithmetic Mean

        //
		// Variance

        //
		// Standard Deviation

        //
		// Covariance

        //
		// Correlation
    }
}