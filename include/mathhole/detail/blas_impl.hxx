// Copyright (c) 2022 Geek pww
// All rights reserved

#ifndef MATHHOLE_DETAIL_BLAS_HXX
#define MATHHOLE_DETAIL_BLAS_HXX

#include <cstdint> // std::int32_t

namespace mathhole {
	namespace blas {
		namespace detail {

			//----------------------------------------------------------------------------//
			// Declarations (for Reference BLAS)

			//
			// xAXPY (e.g. Linear Algebra Add)

			extern "C" void mathhole_saxpy_(std::int32_t const* n, float const* alpha,
											float const x[], std::int32_t const* incx, float y[],
											std::int32_t const* incy);

			extern "C" void mathhole_daxpy_(std::int32_t const* n, double const* alpha,
											double const x[], std::int32_t const* incx, double y[],
											std::int32_t const* incy);

			//
			// xNRM2 (e.g. Vector Euclidean Norm)

			extern "C" float mathhole_snrm2_(std::int32_t const* n, float const x[],
											 std::int32_t const* incx);

			extern "C" double mathhole_dnrm2_(std::int32_t const* n, double const x[],
											  std::int32_t const* incx);

			//
			// xGEMV (e.g. General Matrix Vector Produt)

			extern "C" void mathhole_sgemv_(char const* trans, std::int32_t const* m,
											std::int32_t const* n, float const* alpha,
											float const a[], std::int32_t const* lda,
											float const x[], std::int32_t const* incx,
											float const* beta, float y[], std::int32_t const* incy,
											std::int32_t length_trans);

			extern "C" void mathhole_dgemv_(char const* trans, std::int32_t const* m,
											std::int32_t const* n, double const* alpha,
											double const a[], std::int32_t const* lda,
											double const x[], std::int32_t const* incx,
											double const* beta, double y[],
											std::int32_t const* incy, std::int32_t length_trans);

			//
			// xTRSV (e.g. Triangular Matrix Vector Solve)

			extern "C" void mathhole_strsv_(char const* uplo, char const* trans, char const* diag,
											std::int32_t const* n, const float a[],
											std::int32_t const* lda, float x[],
											std::int32_t const* incx, std::int32_t length_uplo,
											std::int32_t length_trans, std::int32_t length_diag);

			extern "C" void mathhole_dtrsv_(char const* uplo, char const* trans, char const* diag,
											std::int32_t const* n, double const a[],
											std::int32_t const* lda, double x[],
											std::int32_t const* incx, std::int32_t length_uplo,
											std::int32_t length_trans, std::int32_t length_diag);

			//
			// xGEMM (e.g. General Matrix Product)

			extern "C" void mathhole_sgemm_(char const* transa, char const* transb,
											std::int32_t const* m, std::int32_t const* n,
											std::int32_t const* k, float const* alpha,
											float const a[], std::int32_t const* lda,
											float const b[], std::int32_t const* ldb,
											float const* beta, float c[], std::int32_t const* ldc,
											std::int32_t length_transa, std::int32_t length_transb);

			extern "C" void mathhole_dgemm_(char const* transa, char const* transb,
											std::int32_t const* m, std::int32_t const* n,
											std::int32_t const* k, double const* alpha,
											double const a[], std::int32_t const* lda,
											double const b[], std::int32_t const* ldb,
											double const* beta, double c[], std::int32_t const* ldc,
											std::int32_t length_transa, std::int32_t length_transb);

		} // namespace detail
	}	  // namespace blas
} // namespace mathhole

#endif