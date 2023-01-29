 !> \brief \b MATHHOLE_SAXPY
 !
 !  Definition:
 !  ===========
 !
 !       SUBROUTINE MATHHOLE_SAXPY(N,SA,SX,INCX,SY,INCY)
 !
 !       .. Scalar Arguments ..
 !       REAL SA
 !       INTEGER INCX,INCY,N
 !       ..
 !       .. Array Arguments ..
 !       REAL SX(*),SY(*)
 !       ..
 !
 !
 !> \par Purpose:
 !  =============
 !>
 !> \verbatim
 !>
 !>    MATHHOLE_SAXPY constant times a vector plus a vector.
 !>    uses unrolled loops for increments equal to one.
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
 !> \param[in] SA
 !> \verbatim
 !>          SA is REAL
 !>           On entry, SA specifies the scalar alpha.
 !> \endverbatim
 !>
 !> \param[in] SX
 !> \verbatim
 !>          SX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
 !> \endverbatim
 !>
 !> \param[in] INCX
 !> \verbatim
 !>          INCX is INTEGER
 !>         storage spacing between elements of SX
 !> \endverbatim
 !>
 !> \param[in,out] SY
 !> \verbatim
 !>          SY is REAL array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
 !> \endverbatim
 !>
 !> \param[in] INCY
 !> \verbatim
 !>          INCY is INTEGER
 !>         storage spacing between elements of SY
 !> \endverbatim
 !
 !  Authors:
 !  ========
 !
 !> \author Geek pww
 !
 !> \date June 2022
 !
 !> \ingroup single_blas_level1
 !
 !> \par Further Details:
 !  =====================
 !>
 !  =====================================================================
       SUBROUTINE MATHHOLE_SAXPY(N,SA,SX,INCX,SY,INCY)
 !
 !     .. Scalar Arguments ..
       REAL SA
       INTEGER INCX,INCY,N
 !     ..
 !     .. Array Arguments ..
       REAL SX(*),SY(*)
 !     ..
 !
 !  =====================================================================
 !
 !     .. Local Scalars ..
       INTEGER I,IX,IY,M,MP1
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC mod
 !     ..
       IF (n.LE.0) RETURN
       IF (sa.EQ.0.0) RETURN
       IF (incx.EQ.1 .AND. incy.EQ.1) THEN
 !
 !        code for both increments equal to 1
 !
 !
 !        clean-up loop
 !
          m = mod(n,4)
          IF (m.NE.0) THEN
             DO i = 1,m
                sy(i) = sy(i) + sa*sx(i)
             END DO
          END IF
          IF (n.LT.4) RETURN
          mp1 = m + 1
          DO i = mp1,n,4
             sy(i) = sy(i) + sa*sx(i)
             sy(i+1) = sy(i+1) + sa*sx(i+1)
             sy(i+2) = sy(i+2) + sa*sx(i+2)
             sy(i+3) = sy(i+3) + sa*sx(i+3)
          END DO
       ELSE
 !
 !        code for unequal increments or equal increments
 !          not equal to 1
 !
          ix = 1
          iy = 1
          IF (incx.LT.0) ix = (-n+1)*incx + 1
          IF (incy.LT.0) iy = (-n+1)*incy + 1
          DO i = 1,n
           sy(iy) = sy(iy) + sa*sx(ix)
           ix = ix + incx
           iy = iy + incy
          END DO
       END IF
       RETURN
       END