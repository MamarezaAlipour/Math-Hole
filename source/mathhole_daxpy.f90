 !> \brief \b MATHHOLE_DAXPY
 ! Definition:
 !  ===========
 !
 !       SUBROUTINE MATHHOLE_DAXPY(N,DA,DX,INCX,DY,INCY)
 !
 !       .. Scalar Arguments ..
 !       DOUBLE PRECISION DA
 !       INTEGER INCX,INCY,N
 !       ..
 !       .. Array Arguments ..
 !       DOUBLE PRECISION DX(*),DY(*)
 !       ..
 !
 !
 !> \par Purpose:
 !  =============
 !>
 !> \verbatim
 !>
 !>    MATHHOLE_DAXPY constant times a vector plus a vector.
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
 !> \param[in] DA
 !> \verbatim
 !>          DA is DOUBLE PRECISION
 !>           On entry, DA specifies the scalar alpha.
 !> \endverbatim
 !>
 !> \param[in] DX
 !> \verbatim
 !>          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
 !> \endverbatim
 !>
 !> \param[in] INCX
 !> \verbatim
 !>          INCX is INTEGER
 !>         storage spacing between elements of DX
 !> \endverbatim
 !>
 !> \param[in,out] DY
 !> \verbatim
 !>          DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
 !> \endverbatim
 !>
 !> \param[in] INCY
 !> \verbatim
 !>          INCY is INTEGER
 !>         storage spacing between elements of DY
 !> \endverbatim
 !
 !  Author:
 !  ========
 !
 !> \author Parisa Khaleghi
 !
 !> \date April 2022
 !
 !> \ingroup Partow Rayan Rastak
 !>
 !  =====================================================================


    SUBROUTINE MATHHOLE_DAXPY(N, DA, DX, INCX, DY, INCY)
    DOUBLE PRECISION DA
    INTEGER INCX, INCY, N
    DOUBLE PRECISION DX(*), DY(*)
    INTEGER I, IX, IY, M, MP1
    INTRINSIC mod

    IF (n.LE.0) RETURN
    IF (da.EQ.0d0) RETURN
    IF (incx.EQ.1 .AND. incy.EQ.1) THEN
        m = mod(n, 4)
        IF (m.NE.0) THEN
            DO i = 1, m
                dy(i) = dy(i) + da * dx(i)
            END DO
        END IF
        IF (n.LT.4) RETURN
        mp1 = m + 1
        DO i = mp1,n,4
            dy(i) = dy(i) + da * dx(i)
            dy(i + 1) = dy(i + 1) + da * dx(i + 1)
            dy(i + 2) = dy(i + 2) + da * dx(i + 2)
            dy(i + 3) = dy(i + 3) + da * dx(i + 3)
        END DO
    ELSE
        ix = 1
        iy = 1
        IF (incx.LT.0) ix = (-n+1) * incx + 1
        IF (incy.LT.0) iy = (-n+1) * incy + 1
        DO i = 1,n
            dy(iy) = dy(iy) + da * dx(ix)
            ix = ix + incx
            iy = iy + incy
        END DO
    END IF
    RETURN
    END