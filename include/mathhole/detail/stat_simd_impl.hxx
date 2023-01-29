// Copyright (c) 2022 Geek pww
// All rights reserved

#ifndef MATHHOLE_DETAIL_STAT_UNSEQ_HXX
#define MATHHOLE_DETAIL_STAT_UNSEQ_HXX

#include <cmath>
#include <cstdint>
#include <omp.h>

namespace mathhole {
	namespace stat {
		namespace detail {

			//----------------------------------------------------------------------------//
			// Statistics Unsequential Implementations

			//
			// Arithmetic Mean

			template <class Real>
			inline Real arithmetic_mean_simd_impl(std::int32_t n, const Real x[]) {
				// Compute and return the mean of the sequence
				Real sum {0.0};
				std::cout << "\n  arithmetic_mean_simd_impl()" << std::endl;
#pragma omp simd reduction(+ : sum)
				for (std::int32_t i = 0; i < n; ++i)
					sum += x[i];

				return sum / n;
			}

			//
			// Variance

			template <class Real>
			inline Real variance_simd_impl(std::int32_t n, const Real x[]) {
				// Compute the mean of the sequence
				Real sum {0.0};
				std::cout << "\n  variance_mean_simd_impl()" << std::endl;
#pragma omp simd reduction(+ : sum)
				for (std::int32_t i = 0; i < n; ++i)
					sum += x[i];

				Real mean = sum / n;

				// Compute and return the variance of the sequence
				sum = 0.0;

#pragma omp simd reduction(+ : sum)
				for (std::int32_t i = 0; i < n; ++i) {
					Real center = x[i] - mean;
					sum += center * center;
				}

				return sum / (n - 1);
			}

			//
			// Standard Deviation

			template <class Real>
			inline Real standard_deviation_simd_impl(std::int32_t n, const Real x[]) {
				// Compute the mean of the sequence
				Real sum {0.0};
				std::cout << "\n  standard_deviation_simd_impl()" << std::endl;
#pragma omp simd reduction(+ : sum)
				for (std::int32_t i = 0; i < n; ++i)
					sum += x[i];

				Real mean = sum / n;

				// Compute the variance of the sequence
				sum = 0.0;

#pragma omp simd reduction(+ : sum)
				for (std::int32_t i = 0; i < n; ++i) {
					Real center = x[i] - mean;
					sum += center * center;
				}

				Real var = sum / (n - 1);

				// Compute and return the standard deviation of the sequence
				return std::sqrt(var);
			}

			//
			// Covariance

			template <class Real>
			inline Real covariance_simd_impl(std::int32_t n, const Real x[], const Real y[]) {
				// Compute the mean of the two sequences
				Real sum1 {0.0};
				Real sum2 {0.0};
				std::cout << "\n  covariance_simd_impl()" << std::endl;
#pragma omp simd reduction(+ : sum1, sum2)
				for (std::int32_t i = 0; i < n; ++i) {
					sum1 += x[i];
					sum2 += y[i];
				}

				Real x_mean = sum1 / n;
				Real y_mean = sum2 / n;

				// Compute and return the covariance of two sequences
				sum1 = 0.0;

#pragma omp simd reduction(+ : sum1)
				for (std::int32_t i = 0; i < n; ++i) {
					Real x_center = x[i] - x_mean;
					Real y_center = y[i] - y_mean;
					sum1 += x_center * y_center;
				}

				return sum1 / (n - 1);
			}

			//
			// Correlation

			template <class Real>
			inline Real correlation_simd_impl(std::int32_t n, const Real x[], const Real y[]) {
				// Compute the mean of the two sequences
				Real sum1 {0.0};
				Real sum2 {0.0};
				std::cout << "\n  correlation_simd_impl()" << std::endl;
#pragma omp simd reduction(+ : sum1, sum2)
				for (std::int32_t i = 0; i < n; ++i) {
					sum1 += x[i];
					sum2 += y[i];
				}

				Real x_mean = sum1 / n;
				Real y_mean = sum2 / n;

				// Compute the variance of the two sequences
				sum1 = 0.0;
				sum2 = 0.0;

#pragma omp simd reduction(+ : sum1, sum2)
				for (std::int32_t i = 0; i < n; ++i) {
					Real x_center = x[i] - x_mean;
					Real y_center = y[i] - y_mean;
					sum1 += x_center * x_center;
					sum2 += y_center * y_center;
				}

				Real x_var = sum1 / (n - 1);
				Real y_var = sum2 / (n - 1);

				// Compute the standard deviation of two sequences
				Real x_std = std::sqrt(x_var);
				Real y_std = std::sqrt(y_var);

				// Compute and return the correlation between two sequences
				sum1 = 0.0;

#pragma omp simd reduction(+ : sum1)
				for (std::int32_t i = 0; i < n; ++i) {
					Real x_center = x[i] - x_mean;
					Real y_center = y[i] - y_mean;
					sum1 += x_center * y_center;
				}

				Real cov = sum1 / (n - 1);
				return cov / (x_std * y_std);
			}

		} // namespace detail
	}	  // namespace stat
} // namespace mathhole

#endif