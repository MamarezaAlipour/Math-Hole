// Copyright (c) 2022 Parisa Khaleghi
// All rights reserved

#ifndef MATHHOLE_DETAIL_ROOT_HXX
#define MATHHOLE_DETAIL_ROOT_HXX

#include <cassert> // assert
#include <cmath>   // std::abs
#include <cstdint> // std::int32_t
#include <limits>  // std::numeric_limits
#include <tuple>   // std::tuple

namespace mathhole {
	namespace root {
		namespace detail {

			//----------------------------------------------------------------------------//
			// Root-Finding Implementation

			//
			// Bisection

			template <class UnaryFunction, class Real>
			inline std::tuple<Real, std::int32_t> bisection_impl(UnaryFunction&& f, Real a, Real b, Real tol, std::int32_t iter) {
				// Confirm a is the lower bound and b is the upper bound
				assert(a < b);

				Real fa = f(a);
				Real fb = f(b);

				// Confirm f(x0) and f(x1) have opposite signs
				assert(fa * fb <= 0.0);

				std::int32_t max_iter {iter};

				// Local variables used below
				Real c;
				Real fc;

				// Compute root
				iter = 0;
				while (iter < max_iter) {
					++iter;

					c = (a + b) / 2;
					fc = f(c);

					if (fb * fc > 0) {
						b = c;
						fb = fc; // Omit right half of current interval
					} else if (fa * fc > 0) {
						a = c;
						fa = fc; // Omit left half of current interval
					} else {
						a = c;
						b = c;
					}

					// Stopping criteria
					if (b - a < tol)
						break;

					c = (a + b) / 2;
				} // while loop

				return {c, iter};
			}

			//
			// Brent's Method

			template <class UnaryFunction, class Real>
			inline std::tuple<Real, std::int32_t> brent_impl(UnaryFunction&& f, Real x0, Real x1, Real tol, std::int32_t iter) {
				Real fx0 = f(x0);
				Real fx1 = f(x1);

				// Confirm f(x0) and f(x1) have opposite signs
				assert(fx0 * fx1 <= 0.0);

				std::int32_t max_iter {iter};
				Real eps = std::numeric_limits<Real>::epsilon();

				if (std::abs(fx0) < std::abs(fx1))
					std::swap(x0, x1);

				bool mflag {true};
				Real x2 {x0};

				// Local variables used below
				Real fx2;
				Real s;
				Real d;

				// Compute root
				iter = 0;
				while (iter < max_iter) {
					++iter;

					fx0 = f(x0);
					fx1 = f(x1);
					fx2 = f(x2);

					if ((fx0 != fx2) and (fx1 != fx2)) { // Inverse quadratic interpolation
						Real L0 = (x0 * fx1 * fx2) / ((fx0 - fx1) * (fx0 - fx2));
						Real L1 = (x1 * fx0 * fx2) / ((fx1 - fx0) * (fx1 - fx2));
						Real L2 = (x2 * fx1 * fx0) / ((fx2 - fx0) * (fx2 - fx1));
						s = L0 + L1 + L2;
					} else { // Secant method
						s = x1 - ((fx1 * (x1 - x0)) / (fx1 - fx0));
					}

					if (((s < ((fx1 * x0 + x1) / 4.0)) or (s > x1)) or
						((mflag == true) and (std::abs(s - x1) >= std::abs(x1 - x2) / 2.0)) or
						((mflag == false) and (std::abs(s - x1) >= std::abs(x2 - d) / 2.0)) or
						((mflag == true) and (std::abs(x1 - x2) < tol)) or
						((mflag == false) and (std::abs(x1 - d) < tol))) {
						s = (x0 + x1) / 2.0; // Bisection method
						mflag = true;
					} else {
						mflag = false;
					}

					Real fs = f(s);

					d = x2;
					x2 = x1;

					if (fx0 * fs < 0)
						x1 = s;
					else
						x0 = s;

					Real abserr = std::abs(x1 - x0);
					Real relerr = abserr / std::abs(x1 + eps);

					// Stopping criteria
					if (abserr < tol && relerr < tol)
						break;
					else if (std::abs(fx0) < std::abs(fx1))
						std::swap(x0, x1);
				} // while loop

				return {x0, iter};
			}

			//
			// Fixed-Point Method

			template <class UnaryFunction, class Real>
			inline std::tuple<Real, std::int32_t> fixed_point_impl(UnaryFunction&& g, Real x0, Real tol, std::int32_t iter) {
				Real eps = std::numeric_limits<Real>::epsilon();

				std::int32_t max_iter {iter};
				Real x_old {x0};

				// Local variable used below
				Real x_new;

				// Compute root
				iter = 0;
				while (iter < max_iter) {
					++iter;

					x_new = g(x_old);

					Real abserr = std::abs(x_new - x_old);
					Real relerr = abserr / std::abs(x_new + eps);

					// Stopping criteria
					if (abserr < tol && relerr < tol)
						break;

					x_old = x_new;
				} // while loop

				return {x_new, iter};
			}

			//
			// Newton's Method

			template <class UnaryFunction1, class UnaryFunction2, class Real>
			inline std::tuple<Real, std::int32_t> newton_impl(UnaryFunction1&& f, UnaryFunction2&& fp, Real x0, Real tol, std::int32_t iter) {
				Real eps = std::numeric_limits<Real>::epsilon();

				std::int32_t max_iter {iter};
				Real x_old {x0};

				// Local variable used below
				Real x_new;

				// Compute root
				iter = 0;
				while (iter < max_iter) {
					++iter;

					x_new = x_old - (f(x_old) / fp(x_old));

					Real abserr = std::abs(x_new - x_old);
					Real relerr = abserr / std::abs(x_new + eps);

					// Stopping criteria
					if (abserr < tol && relerr < tol)
						break;

					x_old = x_new;
				} // while loop

				return {x_new, iter};
			}

			//
			// Secant Method

			template <class UnaryFunction, class Real>
			inline std::tuple<Real, std::int32_t> secant_impl(UnaryFunction&& f, Real x0, Real x1, Real tol, std::int32_t iter) {
				Real eps = std::numeric_limits<Real>::epsilon();

				std::int32_t max_iter {iter};

				// Local variable used below
				Real x2;

				// Compute root
				iter = 0;
				while (iter < max_iter) {
					++iter;

					Real num = x0 * f(x1) - x1 * f(x0);
					Real den = f(x1) - f(x0);

					x2 = num / den;

					Real abserr = std::abs(x2 - x1);
					Real relerr = abserr / std::abs(x2 + eps);

					// Stopping criteria
					if (abserr < tol && relerr < tol)
						break;

					x0 = x1;
					x1 = x2;
				} // while loop

				return {x2, iter};
			}

		} // namespace detail
	}	  // namespace root
} // namespace mathhole

#endif