 !> \brief \b MATHHOLE_DTRSV
 !
 !  Definition:
 !  ===========
 !
 !       SUBROUTINE MATHHOLE_DTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
 !
 !       .. Scalar Arguments ..
 !       INTEGER INCX,LDA,N
 !       CHARACTER DIAG,TRANS,UPLO
 !       ..
 !       .. Array Arguments ..
 !       DOUBLE PRECISION A(LDA,*),X(*)
 !       ..
 !
 !
 !> \par Purpose:
 !  =============
 !>
 !> \verbatim
 !>
 !> MATHHOLE_DTRSV  solves one of the systems of equations
 !>
 !>    A*x = b,   or   A**T*x = b,
 !>
 !> where b and x are n element vectors and A is an n by n unit, or
 !> non-unit, upper or lower triangular matrix.
 !>
 !> No test for singularity or near-singularity is included in this
 !> routine. Such tests must be performed before calling this routine.
 !> \endverbatim
 !
 !  Arguments:
 !  ==========
 !
 !> \param[in] UPLO
 !> \verbatim
 !>          UPLO is CHARACTER*1
 !>           On entry, UPLO specifies whether the matrix is an upper or
 !>           lower triangular matrix as follows:
 !>
 !>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
 !>
 !>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
 !> \endverbatim
 !>
 !> \param[in] TRANS
 !> \verbatim
 !>          TRANS is CHARACTER*1
 !>           On entry, TRANS specifies the equations to be solved as
 !>           follows:
 !>
 !>              TRANS = 'N' or 'n'   A*x = b.
 !>
 !>              TRANS = 'T' or 't'   A**T*x = b.
 !>
 !>              TRANS = 'C' or 'c'   A**T*x = b.
 !> \endverbatim
 !>
 !> \param[in] DIAG
 !> \verbatim
 !>          DIAG is CHARACTER*1
 !>           On entry, DIAG specifies whether or not A is unit
 !>           triangular as follows:
 !>
 !>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
 !>
 !>              DIAG = 'N' or 'n'   A is not assumed to be unit
 !>                                  triangular.
 !> \endverbatim
 !>
 !> \param[in] N
 !> \verbatim
 !>          N is INTEGER
 !>           On entry, N specifies the order of the matrix A.
 !>           N must be at least zero.
 !> \endverbatim
 !>
 !> \param[in] A
 !> \verbatim
 !>          A is DOUBLE PRECISION array, dimension ( LDA, N )
 !>           Before entry with  UPLO = 'U' or 'u', the leading n by n
 !>           upper triangular part of the array A must contain the upper
 !>           triangular matrix and the strictly lower triangular part of
 !>           A is not referenced.
 !>           Before entry with UPLO = 'L' or 'l', the leading n by n
 !>           lower triangular part of the array A must contain the lower
 !>           triangular matrix and the strictly upper triangular part of
 !>           A is not referenced.
 !>           Note that when  DIAG = 'U' or 'u', the diagonal elements of
 !>           A are not referenced either, but are assumed to be unity.
 !> \endverbatim
 !>
 !> \param[in] LDA
 !> \verbatim
 !>          LDA is INTEGER
 !>           On entry, LDA specifies the first dimension of A as declared
 !>           in the calling (sub) program. LDA must be at least
 !>           max( 1, n ).
 !> \endverbatim
 !>
 !> \param[in,out] X
 !> \verbatim
 !>          X is DOUBLE PRECISION array, dimension at least
 !>           ( 1 + ( n - 1 )*abs( INCX ) ).
 !>           Before entry, the incremented array X must contain the n
 !>           element right-hand side vector b. On exit, X is overwritten
 !>           with the solution vector x.
 !> \endverbatim
 !>
 !> \param[in] INCX
 !> \verbatim
 !>          INCX is INTEGER
 !>           On entry, INCX specifies the increment for the elements of
 !>           X. INCX must not be zero.
 !>
 !>  Level 2 Blas routine.
 !>
 !> \endverbatim
 !
 !  Authors:
 !  ========
 !
 !> \author Geek pww
 !
 !> \date June 2022
 !
 !> \ingroup double_blas_level1
 !
 !  =====================================================================
       SUBROUTINE MATHHOLE_DTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
 !
 !     .. Scalar Arguments ..
       INTEGER INCX,LDA,N
       CHARACTER DIAG,TRANS,UPLO
 !     ..
 !     .. Array Arguments ..
       DOUBLE PRECISION A(LDA,*),X(*)
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       DOUBLE PRECISION ZERO
       parameter(zero=0.0d+0)
 !     ..
 !     .. Local Scalars ..
       DOUBLE PRECISION TEMP
       INTEGER I,INFO,IX,J,JX,KX
       LOGICAL NOUNIT
 !     ..
 !     .. External Functions ..
       LOGICAL MATHHOLE_LSAME
       EXTERNAL mathhole_lsame
 !     ..
 !     .. External Subroutines ..
       EXTERNAL mathhole_xerbla
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC max
 !     ..
 !
 !     Test the input parameters.
 !
       info = 0
       IF (.NOT.mathhole_lsame(uplo,'U') .AND. .NOT.mathhole_lsame(uplo,'L')) THEN
           info = 1
       ELSE IF (.NOT.mathhole_lsame(trans,'N') .AND. .NOT.mathhole_lsame(trans,'T') .AND. .NOT.mathhole_lsame(trans,'C')) THEN
           info = 2
       ELSE IF (.NOT.mathhole_lsame(diag,'U') .AND. .NOT.mathhole_lsame(diag,'N')) THEN
           info = 3
       ELSE IF (n.LT.0) THEN
           info = 4
       ELSE IF (lda.LT.max(1,n)) THEN
           info = 6
       ELSE IF (incx.EQ.0) THEN
           info = 8
       END IF
       IF (info.NE.0) THEN
           CALL mathhole_xerbla('MATHHOLE_DTRSV ',info)
           RETURN
       END IF
 !
 !     Quick return if possible.
 !
       IF (n.EQ.0) RETURN
 !
       nounit = mathhole_lsame(diag,'N')
 !
 !     Set up the start point in X if the increment is not unity. This
 !     will be  ( N - 1 )*INCX  too small for descending loops.
 !
       IF (incx.LE.0) THEN
           kx = 1 - (n-1)*incx
       ELSE IF (incx.NE.1) THEN
           kx = 1
       END IF
 !
 !     Start the operations. In this version the elements of A are
 !     accessed sequentially with one pass through A.
 !
       IF (mathhole_lsame(trans,'N')) THEN
 !
 !        Form  x := inv( A )*x.
 !
           IF (mathhole_lsame(uplo,'U')) THEN
               IF (incx.EQ.1) THEN
                   DO 20 j = n,1,-1
                       IF (x(j).NE.zero) THEN
                           IF (nounit) x(j) = x(j)/a(j,j)
                           temp = x(j)
                           DO 10 i = j - 1,1,-1
                               x(i) = x(i) - temp*a(i,j)
    10                     CONTINUE
                       END IF
    20             CONTINUE
               ELSE
                   jx = kx + (n-1)*incx
                   DO 40 j = n,1,-1
                       IF (x(jx).NE.zero) THEN
                           IF (nounit) x(jx) = x(jx)/a(j,j)
                           temp = x(jx)
                           ix = jx
                           DO 30 i = j - 1,1,-1
                               ix = ix - incx
                               x(ix) = x(ix) - temp*a(i,j)
    30                     CONTINUE
                       END IF
                       jx = jx - incx
    40             CONTINUE
               END IF
           ELSE
               IF (incx.EQ.1) THEN
                   DO 60 j = 1,n
                       IF (x(j).NE.zero) THEN
                           IF (nounit) x(j) = x(j)/a(j,j)
                           temp = x(j)
                           DO 50 i = j + 1,n
                               x(i) = x(i) - temp*a(i,j)
    50                     CONTINUE
                       END IF
    60             CONTINUE
               ELSE
                   jx = kx
                   DO 80 j = 1,n
                       IF (x(jx).NE.zero) THEN
                           IF (nounit) x(jx) = x(jx)/a(j,j)
                           temp = x(jx)
                           ix = jx
                           DO 70 i = j + 1,n
                               ix = ix + incx
                               x(ix) = x(ix) - temp*a(i,j)
    70                     CONTINUE
                       END IF
                       jx = jx + incx
    80             CONTINUE
               END IF
           END IF
       ELSE
 !
 !        Form  x := inv( A**T )*x.
 !
           IF (mathhole_lsame(uplo,'U')) THEN
               IF (incx.EQ.1) THEN
                   DO 100 j = 1,n
                       temp = x(j)
                       DO 90 i = 1,j - 1
                           temp = temp - a(i,j)*x(i)
    90                 CONTINUE
                       IF (nounit) temp = temp/a(j,j)
                       x(j) = temp
   100             CONTINUE
               ELSE
                   jx = kx
                   DO 120 j = 1,n
                       temp = x(jx)
                       ix = kx
                       DO 110 i = 1,j - 1
                           temp = temp - a(i,j)*x(ix)
                           ix = ix + incx
   110                 CONTINUE
                       IF (nounit) temp = temp/a(j,j)
                       x(jx) = temp
                       jx = jx + incx
   120             CONTINUE
               END IF
           ELSE
               IF (incx.EQ.1) THEN
                   DO 140 j = n,1,-1
                       temp = x(j)
                       DO 130 i = n,j + 1,-1
                           temp = temp - a(i,j)*x(i)
   130                 CONTINUE
                       IF (nounit) temp = temp/a(j,j)
                       x(j) = temp
   140             CONTINUE
               ELSE
                   kx = kx + (n-1)*incx
                   jx = kx
                   DO 160 j = n,1,-1
                       temp = x(jx)
                       ix = kx
                       DO 150 i = n,j + 1,-1
                           temp = temp - a(i,j)*x(ix)
                           ix = ix - incx
   150                 CONTINUE
                       IF (nounit) temp = temp/a(j,j)
                       x(jx) = temp
                       jx = jx - incx
   160             CONTINUE
               END IF
           END IF
       END IF
 !
       RETURN
 !
 !     End of MATHHOLE_DTRSV .
 !
       END