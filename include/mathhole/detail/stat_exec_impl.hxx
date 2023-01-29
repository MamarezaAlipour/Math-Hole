// Copyright (c) 2022 Geek pww
// All rights reserved

#ifndef MATHHOLE_STAT_IMPL_HXX
#define MATHHOLE_STAT_IMPL_HXX

#include <cmath>   // std::sqrt
#include <cstdint> // std::int32_t
#include <mathhole/detail/stat_serial_impl.hxx>
#include <mathhole/detail/stat_parallel_impl.hxx>
#include <mathhole/detail/stat_parallel_simd_impl.hxx>
#include <mathhole/detail/stat_simd_impl.hxx>
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
		inline Real arithmetic_mean(ExecutionPolicy&& /*exec*/, std::int32_t n, Real const x[]) {
			if constexpr (std::is_same_v<std::decay_t<ExecutionPolicy>,
										 mathhole::exec::sequential_policy>)
				return detail::arithmetic_mean_serial_impl(n, x);
			else if constexpr (std::is_same_v<std::decay_t<ExecutionPolicy>,
											  mathhole::exec::parallel_policy>)
				return detail::arithmetic_mean_parallel_impl(n, x);
			else if constexpr (std::is_same_v<std::decay_t<ExecutionPolicy>,
											  mathhole::exec::parallel_unsequential_policy>)
				return detail::arithmetic_mean_parallel_simd_impl(n, x);
			else if constexpr (std::is_same_v<std::decay_t<ExecutionPolicy>,
											  mathhole::exec::unsequential_policy>)
				return detail::arithmetic_mean_simd_impl(n, x);
			else
				static_assert(mathhole::exec::is_execution_policy_v<ExecutionPolicy>);
		}

		//
		// Variance

		template <class Real>
		inline Real variance(std::int32_t n, Real const x[]) {
			return detail::variance_serial_impl(n, x);
		}

		template <class ExecutionPolicy, class Real>
		inline Real variance(ExecutionPolicy&& /*exec*/, std::int32_t n, Real const x[]) {
			if constexpr (std::is_same_v<std::decay_t<ExecutionPolicy>,
										 mathhole::exec::sequential_policy>)
				return detail::variance_serial_impl(n, x);
			else if constexpr (std::is_same_v<std::decay_t<ExecutionPolicy>,
											  mathhole::exec::parallel_policy>)
				return detail::variance_parallel_impl(n, x);
			else if constexpr (std::is_same_v<std::decay_t<ExecutionPolicy>,
											  mathhole::exec::parallel_unsequential_policy>)
				return detail::variance_parallel_simd_impl(n, x);
			else if constexpr (std::is_same_v<std::decay_t<ExecutionPolicy>,
											  mathhole::exec::unsequential_policy>)
				return detail::variance_simd_impl(n, x);
			else
				static_assert(mathhole::exec::is_execution_policy_v<ExecutionPolicy>);
		}

		//
		// Standard Deviation

		template <class Real>
		inline Real standard_deviation(std::int32_t n, Real const x[]) {
			return detail::standard_deviation_serial_impl(n, x);
		}

		template <class ExecutionPolicy, class Real>
		inline Real standard_deviation(ExecutionPolicy&& /*exec*/, std::int32_t n, Real const x[]) {
			if constexpr (std::is_same_v<std::decay_t<ExecutionPolicy>,
										 mathhole::exec::sequential_policy>)
				return detail::standard_deviation_serial_impl(n, x);
			else if constexpr (std::is_same_v<std::decay_t<ExecutionPolicy>,
											  mathhole::exec::parallel_policy>)
				return detail::standard_deviation_parallel_impl(n, x);
			else if constexpr (std::is_same_v<std::decay_t<ExecutionPolicy>,
											  mathhole::exec::parallel_unsequential_policy>)
				return detail::standard_deviation_parallel_simd_impl(n, x);
			else if constexpr (std::is_same_v<std::decay_t<ExecutionPolicy>,
											  mathhole::exec::unsequential_policy>)
				return detail::standard_deviation_simd_impl(n, x);
			else
				static_assert(mathhole::exec::is_execution_policy_v<ExecutionPolicy>);
		}

		//
		// Covariance

		template <class Real>
		inline Real covariance(std::int32_t n, Real const x[], Real const y[]) {
			return detail::covariance_serial_impl(n, x, y);
		}

		template <class ExecutionPolicy, class Real>
		inline Real covariance(ExecutionPolicy&& /*exec*/, std::int32_t n, Real const x[],
							   Real const y[]) {
			if constexpr (std::is_same_v<std::decay_t<ExecutionPolicy>,
										 mathhole::exec::sequential_policy>)
				return detail::covariance_serial_impl(n, x, y);
			else if constexpr (std::is_same_v<std::decay_t<ExecutionPolicy>,
											  mathhole::exec::parallel_policy>)
				return detail::covariance_parallel_impl(n, x, y);
			else if constexpr (std::is_same_v<std::decay_t<ExecutionPolicy>,
											  mathhole::exec::parallel_unsequential_policy>)
				return detail::covariance_parallel_simd_impl(n, x, y);
			else if constexpr (std::is_same_v<std::decay_t<ExecutionPolicy>,
											  mathhole::exec::unsequential_policy>)
				return detail::covariance_simd_impl(n, x, y);
			else
				static_assert(mathhole::exec::is_execution_policy_v<ExecutionPolicy>);
		}

		//
		// Correlation

		template <class Real>
		inline Real correlation(std::int32_t n, Real const x[], Real const y[]) {
			return detail::correlation_serial_impl(n, x, y);
		}

		template <class ExecutionPolicy, class Real>
		inline Real correlation(ExecutionPolicy&& /*exec*/, std::int32_t n, Real const x[],
								Real const y[]) {
			if constexpr (std::is_same_v<std::decay_t<ExecutionPolicy>,
										 mathhole::exec::sequential_policy>)
				return detail::correlation_serial_impl(n, x, y);
			else if constexpr (std::is_same_v<std::decay_t<ExecutionPolicy>,
											  mathhole::exec::parallel_policy>)
				return detail::correlation_parallel_impl(n, x, y);
			else if constexpr (std::is_same_v<std::decay_t<ExecutionPolicy>,
											  mathhole::exec::parallel_unsequential_policy>)
				return detail::correlation_parallel_simd_impl(n, x, y);
			else if constexpr (std::is_same_v<std::decay_t<ExecutionPolicy>,
											  mathhole::exec::unsequential_policy>)
				return detail::correlation_simd_impl(n, x, y);
			else
				static_assert(mathhole::exec::is_execution_policy_v<ExecutionPolicy>);
		} 
    } // namespace stat
} // namespace mathhole

#endif