 !> \brief \b MATHHOLE_XERBLA
 !
 !  Definition:
 !  ===========
 !
 !       SUBROUTINE MATHHOLE_XERBLA( SRNAME, INFO )
 !
 !       .. Scalar Arguments ..
 !       CHARACTER*(*)      SRNAME
 !       INTEGER            INFO
 !       ..
 !
 !
 !> \par Purpose:
 !  =============
 !>
 !> \verbatim
 !>
 !> MATHHOLE_XERBLA  is an error handler for the LAPACK routines.
 !> It is called by an LAPACK routine if an input parameter has an
 !> invalid value.  A message is printed and execution stops.
 !>
 !> Installers may consider modifying the STOP statement in order to
 !> call system-specific exception-handling facilities.
 !> \endverbatim
 !
 !  Arguments:
 !  ==========
 !
 !> \param[in] SRNAME
 !> \verbatim
 !>          SRNAME is CHARACTER*(*)
 !>          The name of the routine which called MATHHOLE_XERBLA.
 !> \endverbatim
 !>
 !> \param[in] INFO
 !> \verbatim
 !>          INFO is INTEGER
 !>          The position of the invalid parameter in the parameter list
 !>          of the calling routine.
 !> \endverbatim
 !
 !  Authors:
 !  ========
 !
 !> \author Parisa Khaleghi
 !
 !> \date June 2022
 !
 !> \ingroup aux_blas
 !
 !  =====================================================================
       SUBROUTINE MATHHOLE_XERBLA( SRNAME, INFO )
 !
 !
 !     .. Scalar Arguments ..
       CHARACTER*(*)      SRNAME
       INTEGER            INFO
 !     ..
 !
 ! =====================================================================
 !
 !     .. Intrinsic Functions ..
       INTRINSIC          len_trim
 !     ..
 !     .. Executable Statements ..
 !
       WRITE( *, fmt = 9999 )srname( 1:len_trim( srname ) ), info
 !
       stop
 !
  9999 FORMAT( ' ** On entry to ', a, ' parameter number ', i2, ' had ', 'an illegal value' )
 !
 !     End of MATHHOLE_XERBLA
 !
       END