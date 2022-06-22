 !> \brief \b MATHHOLE_DNRM2
 !
 !  Definition:
 !  ===========
 !
 !       DOUBLE PRECISION FUNCTION MATHHOLE_DNRM2(N,X,INCX)
 !
 !       .. Scalar Arguments ..
 !       INTEGER INCX,N
 !       ..
 !       .. Array Arguments ..
 !       DOUBLE PRECISION X(*)
 !       ..
 !
 !
 !> \par Purpose:
 !  =============
 !>
 !> \verbatim
 !>
 !> MATHHOLE_DNRM2 returns the euclidean norm of a vector via the function
 !> name, so that
 !>
 !>    MATHHOLE_DNRM2 := sqrt( x'*x )
 !> \endverbatim
 !
 !  Arguments:
 !  ==========
 !
 !> \param[in] N
 !> \verbatim
 !>          N is INTEGER
 !>         number of elements in input vector(s)
 !> \endverbatim
 !>
 !> \param[in] X
 !> \verbatim
 !>          X is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
 !> \endverbatim
 !>
 !> \param[in] INCX
 !> \verbatim
 !>          INCX is INTEGER
 !>         storage spacing between elements of DX
 !> \endverbatim
 !
 !  Authors:
 !  ========
 !
 !> \author Parisa Khaleghi
 !
 !> \date June 2022
 !
 !> \ingroup double_blas_level1
 !
 !>
 !  =====================================================================
       DOUBLE PRECISION FUNCTION MATHHOLE_DNRM2(N,X,INCX)
 !
 !     .. Scalar Arguments ..
       INTEGER incx,n
 !     .. Array Arguments ..
       DOUBLE PRECISION x(*)
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       DOUBLE PRECISION one,zero
       parameter(one=1.0d+0,zero=0.0d+0)
 !     ..
 !     .. Local Scalars ..
       DOUBLE PRECISION absxi,norm,scale,ssq
       INTEGER ix
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC abs,sqrt
 !     ..
       IF (n.LT.1 .OR. incx.LT.1) THEN
           norm = zero
       ELSE IF (n.EQ.1) THEN
           norm = abs(x(1))
       ELSE
           scale = zero
           ssq = one
 !        The following loop is equivalent to this call to the LAPACK
 !        auxiliary routine:
 !        CALL DLASSQ( N, X, INCX, SCALE, SSQ )
 !
           DO 10 ix = 1,1 + (n-1)*incx,incx
               IF (x(ix).NE.zero) THEN
                   absxi = abs(x(ix))
                   IF (scale.LT.absxi) THEN
                       ssq = one + ssq* (scale/absxi)**2
                       scale = absxi
                   ELSE
                       ssq = ssq + (absxi/scale)**2
                   END IF
               END IF
    10     CONTINUE
           norm = scale*sqrt(ssq)
       END IF
 !
       MATHHOLE_DNRM2 = norm
       RETURN
 !
 !     End of MATHHOLE_DNRM2.
 !
       END