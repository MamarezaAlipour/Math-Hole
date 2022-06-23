 !> \brief \b MATHHOLE_SNRM2
 !
 !  Definition:
 !  ===========
 !
 !       REAL FUNCTION MATHHOLE_SNRM2(N,X,INCX)
 !
 !       .. Scalar Arguments ..
 !       INTEGER INCX,N
 !       ..
 !       .. Array Arguments ..
 !       REAL X(*)
 !       ..
 !
 !
 !> \par Purpose:
 !  =============
 !>
 !> \verbatim
 !>
 !> MATHHOLE_SNRM2 returns the euclidean norm of a vector via the function
 !> name, so that
 !>
 !>    MATHHOLE_SNRM2 := sqrt( x'*x ).
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
 !>          X is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
 !> \endverbatim
 !>
 !> \param[in] INCX
 !> \verbatim
 !>          INCX is INTEGER
 !>         storage spacing between elements of SX
 !> \endverbatim
 !
 !  Authors:
 !  ========
 !
 !> \author Parisa Khaleghi
 !
 !> \date June 2022
 !
 !> \ingroup single_blas_level1
 !
 !>
 !  =====================================================================
       REAL FUNCTION MATHHOLE_SNRM2(N,X,INCX)
 !
 !     .. Scalar Arguments ..
       INTEGER incx,n
 !     ..
 !     .. Array Arguments ..
       REAL x(*)
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       REAL one,zero
       parameter(one=1.0e+0,zero=0.0e+0)
 !     ..
 !     .. Local Scalars ..
       REAL absxi,norm,scale,ssq
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
 !        CALL SLASSQ( N, X, INCX, SCALE, SSQ )
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
       MATHHOLE_SNRM2 = norm
       RETURN
 !
 !     End of MATHHOLE_SNRM2.
 !
       END