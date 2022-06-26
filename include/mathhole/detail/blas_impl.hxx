// Copyright (c) 2022 Parisa Khaleghi
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

			extern "C" void mathhole_saxpy_(const std::int32_t* n, const float* alpha,
											const float x[], const std::int32_t* incx, float y[],
											const std::int32_t* incy);

			extern "C" void mathhole_daxpy_(const std::int32_t* n, const double* alpha,
											const double x[], const std::int32_t* incx, double y[],
											const std::int32_t* incy);

			//
			// xNRM2 (e.g. Vector Euclidean Norm)

			extern "C" float mathhole_snrm2_(const std::int32_t* n, const float x[],
											 const std::int32_t* incx);

			extern "C" double mathhole_dnrm2_(const std::int32_t* n, const double x[],
											  const std::int32_t* incx);

			//
			// xGEMV (e.g. General Matrix Vector Produt)

			extern "C" void mathhole_sgemv_(const char* trans, const std::int32_t* m,
											const std::int32_t* n, const float* alpha,
											const float a[], const std::int32_t* lda,
											const float x[], const std::int32_t* incx,
											const float* beta, float y[], const std::int32_t* incy,
											std::int32_t length_trans);

			extern "C" void mathhole_dgemv_(const char* trans, const std::int32_t* m,
											const std::int32_t* n, const double* alpha,
											const double a[], const std::int32_t* lda,
											const double x[], const std::int32_t* incx,
											const double* beta, double y[],
											const std::int32_t* incy, std::int32_t length_trans);

			//
			// xTRSV (e.g. Triangular Matrix Vector Solve)

			extern "C" void mathhole_strsv_(const char* uplo, const char* trans, const char* diag,
											const std::int32_t* n, const float a[],
											const std::int32_t* lda, float x[],
											const std::int32_t* incx, std::int32_t length_uplo,
											std::int32_t length_trans, std::int32_t length_diag);

			extern "C" void mathhole_dtrsv_(const char* uplo, const char* trans, const char* diag,
											const std::int32_t* n, const double a[],
											const std::int32_t* lda, double x[],
											const std::int32_t* incx, std::int32_t length_uplo,
											std::int32_t length_trans, std::int32_t length_diag);

			//
			// xGEMM (e.g. General Matrix Product)

			extern "C" void mathhole_sgemm_(const char* transa, const char* transb,
											const std::int32_t* m, const std::int32_t* n,
											const std::int32_t* k, const float* alpha,
											const float a[], const std::int32_t* lda,
											const float b[], const std::int32_t* ldb,
											const float* beta, float c[], const std::int32_t* ldc,
											std::int32_t length_transa, std::int32_t length_transb);

			extern "C" void mathhole_dgemm_(const char* transa, const char* transb,
											const std::int32_t* m, const std::int32_t* n,
											const std::int32_t* k, const double* alpha,
											const double a[], const std::int32_t* lda,
											const double b[], const std::int32_t* ldb,
											const double* beta, double c[], const std::int32_t* ldc,
											std::int32_t length_transa, std::int32_t length_transb);

		} // namespace detail
	}	  // namespace blas
} // namespace mathhole

#endif