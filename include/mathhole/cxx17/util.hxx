// Copyright (c) 2022 Geek pww
// All rights reserved

#ifndef MATHHOLE_CXX17_UTIL_HXX
#define MATHHOLE_CXX17_UTIL_HXX

#include <algorithm>   // std::max
#include <cmath>	   // std::abs
#include <cstdint>	   // std::int32_t
#include <iomanip>	   // std::setw
#include <iostream>	   // std::cout, std::endl
#include <string_view> // std::string_view

namespace mathhole {
	namespace util {

		//------------------------------------------------------------------------------
		// Declarations

		template <class Real>
		inline Real max_abs_difference(std::int32_t size, const Real a[], const Real b[]);

		template <class Real>
		inline void print_vector(std::string_view s, std::int32_t size, const Real a[]);

		template <class Real>
		inline void print_vector(std::string_view s, std::int32_t size, const Real a[],
								 std::int32_t precision, std::int32_t width);

		template <class Real>
		inline void print_matrix(std::string_view s, std::int32_t rows, std::int32_t cols,
								 const Real a[]);

		template <class Real>
		inline void print_matrix(std::string_view s, std::int32_t rows, std::int32_t cols,
								 const Real a[], std::int32_t precision, std::int32_t width);

		//------------------------------------------------------------------------------
		// Defintions

		//
		// Maximum Absolute Difference

		// Return the maximum absolute difference between two ranges
		template <class Real>
		Real max_abs_difference(std::int32_t size, const Real a[], const Real b[]) {
			Real diff {0.0};

			for (std::int32_t i = 0; i < size; ++i) {
				Real temp = std::abs(a[i] - b[i]);
				diff = std::max(temp, diff);
			}

			return diff;
		}

		//
		// Print Vector

		// Display all elements of a one-dimensional container
		template <class Real>
		void print_vector(std::string_view s, std::int32_t size, const Real a[],
						  std::int32_t precision, std::int32_t width) {
			std::cout << s << " = [" << std::endl;

			for (std::int32_t i = 0; i < size; ++i)
				std::cout << " " << std::setprecision(precision) << std::setw(width) << a[i];

			std::cout << "\n]" << std::endl;
		}

		template <class Real>
		void print_vector(std::string_view s, std::int32_t size, const Real a[]) {
			return print_vector(s, size, a, 6, 8);
		}

		//
		// Print Matrix

		// Display all elements of a two-dimensional container
		template <class Real>
		void print_matrix(std::string_view s, std::int32_t rows, std::int32_t cols, const Real a[],
						  std::int32_t precision, std::int32_t width) {
			std::cout << s << " = [" << std::endl;

			for (std::int32_t i = 0; i < rows; ++i) {
				for (std::int32_t j = 0; j < cols; ++j) {
					std::cout << " " << std::setprecision(precision) << std::setw(width)
							  << a[j * rows + i];
				}

				std::cout << std::endl;
			}

			std::cout << "]" << std::endl;
		}

		template <class Real>
		void print_matrix(std::string_view s, std::int32_t rows, std::int32_t cols,
						  const Real a[]) {
			return print_matrix(s, rows, cols, a, 6, 8);
		}

	} // namespace util
} // namespace mathhole

#endif