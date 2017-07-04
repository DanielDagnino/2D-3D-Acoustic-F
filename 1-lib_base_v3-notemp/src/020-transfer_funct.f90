  !*********************************************************************/
  !*****     Author: D. Dagnino     ************************************/
  !*********************************************************************/
  module m_dag_transfer_funct    
    
    interface butterworth
      module procedure butterworth, butterworth_r
    end interface
    
    interface butterworth0
      module procedure butterworth0, butterworth0_r
    end interface
    
    contains
    !*********************************************************************/
    !*********************************************************************/
    !*****     Double     ************************************************/
    !*********************************************************************/
    !*********************************************************************/
    elemental complex(r8) function butterworth( f, of, f0, f1, type_filt )
      use m_dag_data_kind
      use m_dag_mat
      implicit none
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r8
      real(prcsn), intent(in) :: f        ! Frequency.
      integer, intent(in)     :: of       ! Order filter.
      real(prcsn), intent(in) :: f0, f1   ! Filter frequency.
      integer, intent(in)     :: type_filt
      ! The variables which are generated inside the function.
      complex(prcsn) :: s
      integer :: k
      logical :: parity

      !---------------------------------------------------------------------/
      butterworth = cmplx( 1._prcsn, 0._prcsn )
      s = cmplx( 0._prcsn, f )

      if ( mod(of,2) == 0 ) then
        parity = .false.
      else
        parity = .true.
      end if

      !---------------------------------------------------------------------/
      ! low-pass filter
      if( type_filt == 1 ) then

        do k=1,of/2
          butterworth = butterworth * ( (s/f0)*(s/f0) - 2._prcsn*(s/f0)*cos(real(2*k+of-1,prcsn)*pi/real(2*of,prcsn)) + 1._prcsn )
        end do
        if ( parity ) butterworth = ((s/f0)+1._prcsn)*butterworth

      ! high-pass filter
      else if ( type_filt == 2 ) then

        do k=1,of/2
          butterworth = butterworth * ( (f1/s)*(f1/s) - 2._prcsn*(f1/s)*cos(real(2*k+of-1,prcsn)*pi/real(2*of,prcsn)) + 1._prcsn )
        end do
        if ( parity ) butterworth = ((f1/s)+1._prcsn)*butterworth

      ! band-pass filter
      else if ( type_filt == 3 ) then

        do k=1,of/2
          butterworth = butterworth * ( ((s+f0*f1/s)/(f1-f0))*((s+f0*f1/s)/(f1-f0)) - 2._prcsn*((s+f0*f1/s)/(f1-f0))*cos(real(2*k+of-1,prcsn)*pi/real(2*of,prcsn)) + 1._prcsn )
        end do
        if ( parity ) butterworth = (((s+f0*f1/s)/(f1-f0))+1._prcsn)*butterworth

      end if

      !---------------------------------------------------------------------/
      ! 
      butterworth = 1._prcsn/butterworth

    end function butterworth
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    elemental real(r8) function butterworth0( f, of, f0, f1, type_filt )
      use m_dag_data_kind
      implicit none
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r8
      real(prcsn), intent(in) :: f       ! Frequency.
      integer, intent(in)     :: of      ! Order filter.
      real(prcsn), intent(in) :: f0, f1  ! Filter frequency.
      integer, intent(in)     :: type_filt

      !---------------------------------------------------------------------/
      ! low-pass filter
      if ( type_filt == 1 ) then

        butterworth0 = 1._prcsn/sqrt( 1._prcsn + (f/f0)**(2*of) )

      ! high-pass filter
      else if ( type_filt == 2 ) then

        butterworth0 = 1._prcsn/sqrt( 1._prcsn + (f1/f)**(2*of) )

      ! band-pass filter
      else if ( type_filt == 3 ) then

        butterworth0 = 1._prcsn/sqrt( 1._prcsn + ((f+f0*f1/f)/(f1+f0))**(2*of) )

      end if

    end function butterworth0





    !*********************************************************************/
    !*********************************************************************/
    !*****     Real     **************************************************/
    !*********************************************************************/
    !*********************************************************************/
    elemental complex(r4) function butterworth_r( f, of, f0, f1, type_filt )
      use m_dag_data_kind
      use m_dag_mat
      implicit none
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r4
      real(prcsn), intent(in) :: f        ! Frequency.
      integer, intent(in)     :: of       ! Order filter.
      real(prcsn), intent(in) :: f0, f1   ! Filter frequency.
      integer, intent(in)     :: type_filt
      ! The variables which are generated inside the function.
      complex(prcsn) :: s
      integer :: k
      logical :: parity

      !---------------------------------------------------------------------/
      butterworth_r = cmplx( 1._prcsn, 0._prcsn )
      s = cmplx( 0._prcsn, f )

      if ( mod(of,2) == 0 ) then
        parity = .false.
      else
        parity = .true.
      end if

      !---------------------------------------------------------------------/
      ! low-pass filter
      if( type_filt == 1 ) then

        do k=1,of/2
          butterworth_r = butterworth_r * ( (s/f0)*(s/f0) - 2._prcsn*(s/f0)*cos(real(2*k+of-1,prcsn)*pi/real(2*of,prcsn)) + 1._prcsn )
        end do
        if ( parity ) butterworth_r = ((s/f0)+1._prcsn)*butterworth_r

      ! high-pass filter
      else if ( type_filt == 2 ) then

        do k=1,of/2
          butterworth_r = butterworth_r * ( (f1/s)*(f1/s) - 2._prcsn*(f1/s)*cos(real(2*k+of-1,prcsn)*pi/real(2*of,prcsn)) + 1._prcsn )
        end do
        if ( parity ) butterworth_r = ((f1/s)+1._prcsn)*butterworth_r

      ! band-pass filter
      else if ( type_filt == 3 ) then

        do k=1,of/2
          butterworth_r = butterworth_r * ( ((s+f0*f1/s)/(f1-f0))*((s+f0*f1/s)/(f1-f0)) - 2._prcsn*((s+f0*f1/s)/(f1-f0))*cos(real(2*k+of-1,prcsn)*pi/real(2*of,prcsn)) + 1._prcsn )
        end do
        if ( parity ) butterworth_r = (((s+f0*f1/s)/(f1-f0))+1._prcsn)*butterworth_r

      end if

      !---------------------------------------------------------------------/
      ! 
      butterworth_r = 1._prcsn/butterworth_r

    end function butterworth_r
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    elemental real(r4) function butterworth0_r( f, of, f0, f1, type_filt )
      use m_dag_data_kind
      implicit none
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r4
      real(prcsn), intent(in) :: f       ! Frequency.
      integer, intent(in)     :: of      ! Order filter.
      real(prcsn), intent(in) :: f0, f1  ! Filter frequency.
      integer, intent(in)     :: type_filt

      !---------------------------------------------------------------------/
      ! low-pass filter
      if ( type_filt == 1 ) then

        butterworth0_r = 1._prcsn/sqrt( 1._prcsn + (f/f0)**(2*of) )

      ! high-pass filter
      else if ( type_filt == 2 ) then

        butterworth0_r = 1._prcsn/sqrt( 1._prcsn + (f1/f)**(2*of) )

      ! band-pass filter
      else if ( type_filt == 3 ) then

        butterworth0_r = 1._prcsn/sqrt( 1._prcsn + ((f+f0*f1/f)/(f1+f0))**(2*of) )

      end if

    end function butterworth0_r
    
    
    
  !***********************************************************************/
  end module m_dag_transfer_funct