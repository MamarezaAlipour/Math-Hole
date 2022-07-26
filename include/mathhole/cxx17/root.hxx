// Copyright (c) 2022 Parisa Khaleghi
// All rights reserved

#ifndef MATHHOLE_CXX17_ROOT_HXX
#define MATHHOLE_CXX17_ROOT_HXX

#include <cstdint> // std::int32_t
#include <mathhole/detail/root_impl.hxx>
#include <tuple>   // std::tuple
#include <utility> // std::forward

namespace mathhole {
	namespace root {

		//----------------------------------------------------------------------------//
		// Declarations

		template <class UnaryFunction, class Real>
		inline std::tuple<Real, std::int32_t> bisection(UnaryFunction&& f, Real a, Real b);

		template <class UnaryFunction, class Real>
		inline std::tuple<Real, std::int32_t> bisection(UnaryFunction&& f, Real a, Real b, Real tol,
														std::int32_t maxiter);

		template <class UnaryFunction, class Real>
		inline std::tuple<Real, std::int32_t> brent(UnaryFunction&& f, Real x0, Real x1);

		template <class UnaryFunction, class Real>
		inline std::tuple<Real, std::int32_t> brent(UnaryFunction&& f, Real x0, Real x1, Real tol,
													std::int32_t maxiter);

		template <class UnaryFunction, class Real>
		inline std::tuple<Real, std::int32_t> fixed_point(UnaryFunction&& g, Real x0);

		template <class UnaryFunction, class Real>
		inline std::tuple<Real, std::int32_t> fixed_point(UnaryFunction&& g, Real x0, Real tol,
														  std::int32_t maxiter);

		template <class UnaryFunction1, class UnaryFunction2, class Real>
		inline std::tuple<Real, std::int32_t> newton(UnaryFunction1&& f, UnaryFunction2&& fp,
													 Real x0);

		template <class UnaryFunction1, class UnaryFunction2, class Real>
		inline std::tuple<Real, std::int32_t> newton(UnaryFunction1&& f, UnaryFunction2&& fp,
													 Real x0, Real tol, std::int32_t maxiter);

		template <class UnaryFunction, class Real>
		inline std::tuple<Real, std::int32_t> secant(UnaryFunction&& f, Real x0, Real x1);

		template <class UnaryFunction, class Real>
		inline std::tuple<Real, std::int32_t> secant(UnaryFunction&& f, Real x0, Real x1, Real tol,
													 std::int32_t maxiter);

		//----------------------------------------------------------------------------//
		// Definitions

		//
		// Bisection Method

		template <class UnaryFunction, class Real>
		inline std::tuple<Real, std::int32_t> bisection(UnaryFunction&& f, Real a, Real b, Real tol,
														std::int32_t maxiter) {
			return detail::bisection_impl(std::forward<UnaryFunction>(f), a, b, tol, maxiter);
		}

		template <class UnaryFunction>
		inline std::tuple<float, std::int32_t> bisection(UnaryFunction&& f, float a, float b) {
			float tol {1e-05};
			std::int32_t maxiter {100};
			return detail::bisection_impl(std::forward<UnaryFunction>(f), a, b, tol, maxiter);
		}

		template <class UnaryFunction>
		inline std::tuple<double, std::int32_t> bisection(UnaryFunction&& f, double a, double b) {
			double tol {1e-10};
			std::int32_t maxiter {100};
			return detail::bisection_impl(std::forward<UnaryFunction>(f), a, b, tol, maxiter);
		}

		//
		// Brent's Method

		template <class UnaryFunction, class Real>
		inline std::tuple<Real, std::int32_t> brent(UnaryFunction&& f, Real x0, Real x1, Real tol,
													std::int32_t maxiter) {
			return detail::brent_impl(std::forward<UnaryFunction>(f), x0, x1, tol, maxiter);
		}

		template <class UnaryFunction>
		inline std::tuple<float, std::int32_t> brent(UnaryFunction&& f, float x0, float x1) {
			float tol {1e-05};
			std::int32_t maxiter {100};
			return detail::brent_impl(std::forward<UnaryFunction>(f), x0, x1, tol, maxiter);
		}

		template <class UnaryFunction>
		inline std::tuple<double, std::int32_t> brent(UnaryFunction&& f, double x0, double x1) {
			double tol {1e-10};
			std::int32_t maxiter {100};
			return detail::brent_impl(std::forward<UnaryFunction>(f), x0, x1, tol, maxiter);
		}

		//
		// Fixed-Point Method

		template <class UnaryFunction, class Real>
		inline std::tuple<Real, std::int32_t> fixed_point(UnaryFunction&& g, Real x0, Real tol,
														  std::int32_t maxiter) {
			return detail::fixed_point_impl(std::forward<UnaryFunction>(g), x0, tol, maxiter);
		}

		template <class UnaryFunction>
		inline std::tuple<float, std::int32_t> fixed_point(UnaryFunction&& g, float x0) {
			float tol {1e-05};
			std::int32_t maxiter {100};
			return detail::fixed_point_impl(std::forward<UnaryFunction>(g), x0, tol, maxiter);
		}

		template <class UnaryFunction>
		inline std::tuple<double, std::int32_t> fixed_point(UnaryFunction&& g, double x0) {
			double tol {1e-10};
			std::int32_t maxiter {100};
			return detail::fixed_point_impl(std::forward<UnaryFunction>(g), x0, tol, maxiter);
		}

		//
		// Newton's Method

		template <class UnaryFunction1, class UnaryFunction2, class Real>
		inline std::tuple<Real, std::int32_t> newton(UnaryFunction1&& f, UnaryFunction2&& fp,
													 Real x0, Real tol, std::int32_t maxiter) {
			return detail::newton_impl(std::forward<UnaryFunction1>(f),
									   std::forward<UnaryFunction2>(fp), x0, tol, maxiter);
		}

		template <class UnaryFunction1, class UnaryFunction2>
		inline std::tuple<float, std::int32_t> newton(UnaryFunction1&& f, UnaryFunction2&& fp,
													  float x0) {
			float tol {1e-05};
			std::int32_t maxiter {100};
			return detail::newton_impl(std::forward<UnaryFunction1>(f),
									   std::forward<UnaryFunction2>(fp), x0, tol, maxiter);
		}

		template <class UnaryFunction1, class UnaryFunction2>
		inline std::tuple<double, std::int32_t> newton(UnaryFunction1&& f, UnaryFunction2&& fp,
													   double x0) {
			double tol {1e-10};
			std::int32_t maxiter {100};
			return detail::newton_impl(std::forward<UnaryFunction1>(f),
									   std::forward<UnaryFunction2>(fp), x0, tol, maxiter);
		}

		//
		// Secant Method

		template <class UnaryFunction, class Real>
		inline std::tuple<Real, std::int32_t> secant(UnaryFunction&& f, Real x0, Real x1, Real tol,
													 std::int32_t maxiter) {
			return detail::secant_impl(std::forward<UnaryFunction>(f), x0, x1, tol, maxiter);
		}

		template <class UnaryFunction>
		inline std::tuple<float, std::int32_t> secant(UnaryFunction&& f, float x0, float x1) {
			float tol {1e-05};
			std::int32_t maxiter {100};
			return detail::secant_impl(std::forward<UnaryFunction>(f), x0, x1, tol, maxiter);
		}

		template <class UnaryFunction>
		inline std::tuple<double, std::int32_t> secant(UnaryFunction&& f, double x0, double x1) {
			double tol {1e-10};
			std::int32_t maxiter {100};
			return detail::secant_impl(std::forward<UnaryFunction>(f), x0, x1, tol, maxiter);
		}

	} // namespace root
} // namespace mathhole

#endif