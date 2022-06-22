 !> \brief \b MATHHOLE_DGEMV
 !  Definition:
 !  ===========
 !
 !       SUBROUTINE MATHHOLE_DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
 !
 !       .. Scalar Arguments ..
 !       DOUBLE PRECISION ALPHA,BETA
 !       INTEGER INCX,INCY,LDA,M,N
 !       CHARACTER TRANS
 !       ..
 !       .. Array Arguments ..
 !       DOUBLE PRECISION A(LDA,*),X(*),Y(*)
 !       ..
 !
 !
 !> \par Purpose:
 !  =============
 !>
 !> \verbatim
 !>
 !> MATHHOLE_DGEMV  performs one of the matrix-vector operations
 !>
 !>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
 !>
 !> where alpha and beta are scalars, x and y are vectors and A is an
 !> m by n matrix.
 !> \endverbatim
 !
 !  Arguments:
 !  ==========
 !
 !> \param[in] TRANS
 !> \verbatim
 !>          TRANS is CHARACTER*1
 !>           On entry, TRANS specifies the operation to be performed as
 !>           follows:
 !>
 !>              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
 !>
 !>              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.
 !>
 !>              TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y.
 !> \endverbatim
 !>
 !> \param[in] M
 !> \verbatim
 !>          M is INTEGER
 !>           On entry, M specifies the number of rows of the matrix A.
 !>           M must be at least zero.
 !> \endverbatim
 !>
 !> \param[in] N
 !> \verbatim
 !>          N is INTEGER
 !>           On entry, N specifies the number of columns of the matrix A.
 !>           N must be at least zero.
 !> \endverbatim
 !>
 !> \param[in] ALPHA
 !> \verbatim
 !>          ALPHA is DOUBLE PRECISION.
 !>           On entry, ALPHA specifies the scalar alpha.
 !> \endverbatim
 !>
 !> \param[in] A
 !> \verbatim
 !>          A is DOUBLE PRECISION array, dimension ( LDA, N )
 !>           Before entry, the leading m by n part of the array A must
 !>           contain the matrix of coefficients.
 !> \endverbatim
 !>
 !> \param[in] LDA
 !> \verbatim
 !>          LDA is INTEGER
 !>           On entry, LDA specifies the first dimension of A as declared
 !>           in the calling (sub) program. LDA must be at least
 !>           max( 1, m ).
 !> \endverbatim
 !>
 !> \param[in] X
 !> \verbatim
 !>          X is DOUBLE PRECISION array, dimension at least
 !>           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
 !>           and at least
 !>           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
 !>           Before entry, the incremented array X must contain the
 !>           vector x.
 !> \endverbatim
 !>
 !> \param[in] INCX
 !> \verbatim
 !>          INCX is INTEGER
 !>           On entry, INCX specifies the increment for the elements of
 !>           X. INCX must not be zero.
 !> \endverbatim
 !>
 !> \param[in] BETA
 !> \verbatim
 !>          BETA is DOUBLE PRECISION.
 !>           On entry, BETA specifies the scalar beta. When BETA is
 !>           supplied as zero then Y need not be set on input.
 !> \endverbatim
 !>
 !> \param[in,out] Y
 !> \verbatim
 !>          Y is DOUBLE PRECISION array, dimension at least
 !>           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
 !>           and at least
 !>           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
 !>           Before entry with BETA non-zero, the incremented array Y
 !>           must contain the vector y. On exit, Y is overwritten by the
 !>           updated vector y.
 !> \endverbatim
 !>
 !> \param[in] INCY
 !> \verbatim
 !>          INCY is INTEGER
 !>           On entry, INCY specifies the increment for the elements of
 !>           Y. INCY must not be zero.
 !> \endverbatim
 !
 !  Authors:
 !  ========
 !
 !> \author Parisa Khaleghi
 !
 !> \date June 2022
 !
 !> \ingroup double_blas_level2
 !
 !> \par Further Details:
 !  =====================
 !>
 !> \verbatim
 !>
 !>  Level 2 Blas routine.
 !>  The vector and matrix arguments are not referenced when N = 0, or M = 0
 !>
 !  =====================================================================
       SUBROUTINE MATHHOLE_DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
 !
 !
 !     .. Scalar Arguments ..
       DOUBLE PRECISION ALPHA,BETA
       INTEGER INCX,INCY,LDA,M,N
       CHARACTER TRANS
 !     ..
 !     .. Array Arguments ..
       DOUBLE PRECISION A(LDA,*),X(*),Y(*)
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       DOUBLE PRECISION ONE,ZERO
       parameter(one=1.0d+0,zero=0.0d+0)
 !     ..
 !     .. Local Scalars ..
       DOUBLE PRECISION TEMP
       INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
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
       IF (.NOT.mathhole_lsame(trans,'N') .AND. .NOT.mathhole_lsame(trans,'T') .AND. .NOT.mathhole_lsame(trans,'C')) THEN
           info = 1
       ELSE IF (m.LT.0) THEN
           info = 2
       ELSE IF (n.LT.0) THEN
           info = 3
       ELSE IF (lda.LT.max(1,m)) THEN
           info = 6
       ELSE IF (incx.EQ.0) THEN
           info = 8
       ELSE IF (incy.EQ.0) THEN
           info = 11
       END IF
       IF (info.NE.0) THEN
           CALL mathhole_xerbla('MATHHOLE_DGEMV ',info)
           RETURN
       END IF
 !
 !     Quick return if possible.
 !
       IF ((m.EQ.0) .OR. (n.EQ.0) .OR. ((alpha.EQ.zero).AND. (beta.EQ.one))) RETURN
 !
 !     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
 !     up the start points in  X  and  Y.
 !
       IF (mathhole_lsame(trans,'N')) THEN
           lenx = n
           leny = m
       ELSE
           lenx = m
           leny = n
       END IF
       IF (incx.GT.0) THEN
           kx = 1
       ELSE
           kx = 1 - (lenx-1)*incx
       END IF
       IF (incy.GT.0) THEN
           ky = 1
       ELSE
           ky = 1 - (leny-1)*incy
       END IF
 !
 !     Start the operations. In this version the elements of A are
 !     accessed sequentially with one pass through A.
 !
 !     First form  y := beta*y.
 !
       IF (beta.NE.one) THEN
           IF (incy.EQ.1) THEN
               IF (beta.EQ.zero) THEN
                   DO 10 i = 1,leny
                       y(i) = zero
    10             CONTINUE
               ELSE
                   DO 20 i = 1,leny
                       y(i) = beta*y(i)
    20             CONTINUE
               END IF
           ELSE
               iy = ky
               IF (beta.EQ.zero) THEN
                   DO 30 i = 1,leny
                       y(iy) = zero
                       iy = iy + incy
    30             CONTINUE
               ELSE
                   DO 40 i = 1,leny
                       y(iy) = beta*y(iy)
                       iy = iy + incy
    40             CONTINUE
               END IF
           END IF
       END IF
       IF (alpha.EQ.zero) RETURN
       IF (mathhole_lsame(trans,'N')) THEN
 !
 !        Form  y := alpha*A*x + y.
 !
           jx = kx
           IF (incy.EQ.1) THEN
               DO 60 j = 1,n
                   temp = alpha*x(jx)
                   DO 50 i = 1,m
                       y(i) = y(i) + temp*a(i,j)
    50             CONTINUE
                   jx = jx + incx
    60         CONTINUE
           ELSE
               DO 80 j = 1,n
                   temp = alpha*x(jx)
                   iy = ky
                   DO 70 i = 1,m
                       y(iy) = y(iy) + temp*a(i,j)
                       iy = iy + incy
    70             CONTINUE
                   jx = jx + incx
    80         CONTINUE
           END IF
       ELSE
 !
 !        Form  y := alpha*A**T*x + y.
 !
           jy = ky
           IF (incx.EQ.1) THEN
               DO 100 j = 1,n
                   temp = zero
                   DO 90 i = 1,m
                       temp = temp + a(i,j)*x(i)
    90             CONTINUE
                   y(jy) = y(jy) + alpha*temp
                   jy = jy + incy
   100         CONTINUE
           ELSE
               DO 120 j = 1,n
                   temp = zero
                   ix = kx
                   DO 110 i = 1,m
                       temp = temp + a(i,j)*x(ix)
                       ix = ix + incx
   110             CONTINUE
                   y(jy) = y(jy) + alpha*temp
                   jy = jy + incy
   120         CONTINUE
           END IF
       END IF
 !
       RETURN
 !
 !     End of MATHHOLE_DGEMV .
 !
       END