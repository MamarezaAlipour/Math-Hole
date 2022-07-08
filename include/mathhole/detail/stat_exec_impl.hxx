// Copyright (c) 2022 Parisa Khaleghi
// All rights reserved

#ifndef MATHHOLE_STAT_IMPL_HXX
#define MATHHOLE_STAT_IMPL_HXX

#include <cmath>   // std::sqrt
#include <cstdint> // std::int32_t
#include <mathhole/detail/stat_serial_impl.hxx>
#include <mathhole/detail/stat_parallel_impl.hxx>
#include <mathhole/detail/stat_parallel_simd_impl.hxx>
#include <mathhole/cxx17/exec.hxx>

namespace mathhole {
	namespace stat {

		//----------------------------------------------------------------------------//
		// Declarations

		template <class Real>
		inline Real arithmetic_mean(std::int32_t n, Real const x[]);

		template <class ExecutionPolicy, class Real>
		inline Real arithmetic_mean(ExecutionPolicy&& exec, std::int32_t n, Real const x[]);

		template <class Real>
		inline Real variance(std::int32_t n, Real const x[]);

		template <class ExecutionPolicy, class Real>
		inline Real variance(ExecutionPolicy&& exec, std::int32_t n, Real const x[]);

		template <class Real>
		inline Real standard_deviation(std::int32_t n, Real const x[]);

		template <class ExecutionPolicy, class Real>
		inline Real standard_deviation(ExecutionPolicy&& exec, std::int32_t n, Real const x[]);

		template <class Real>
		inline Real covariance(std::int32_t n, Real const x[], Real const y[]);

		template <class ExecutionPolicy, class Real>
		inline Real covariance(ExecutionPolicy&& exec, std::int32_t n, Real const x[],
							   Real const y[]);

		template <class Real>
		inline Real correlation(std::int32_t n, Real const x[], Real const y[]);

		template <class ExecutionPolicy, class Real>
		inline Real correlation(ExecutionPolicy&& exec, std::int32_t n, Real const x[],
								Real const y[]);

		//----------------------------------------------------------------------------//
		// Definitions

		//
		// Arithmetic Mean
        template <class Real>
		inline Real arithmetic_mean(std::int32_t n, Real const x[]) {
			return detail::arithmetic_mean_serial_impl(n, x);
		}

        template <class ExecutionPolicy, class Real>
		inline Real arithmetic_mean(ExecutionPolicy&& /*exec*/, std::int32_t n, Real const x[]){
            // soon...
        }
    }
}