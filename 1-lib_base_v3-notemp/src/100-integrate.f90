  !*********************************************************************/
  !** This code has been done in the Barcelona Center for Subsurface 
  !** Imaging (BCSI).
  !** Authors: Daniel Dagnino.
  !*********************************************************************/
  module m_dag_integrate
    implicit none
    
    interface integrate_3o
      module procedure integrate_3o, integrate_3o_r
    end interface
    
    contains
    
    !*********************************************************************/
    !*********************************************************************/
    !*****     Double     ************************************************/
    !*********************************************************************/
    !*********************************************************************/
    pure subroutine integrate_3o( f, nx, dx, int_f, nx_int )
      use m_dag_data_kind
      use m_dag_interpol
      implicit none
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r8
      integer, intent(in) :: nx, nx_int
      real(prcsn), intent(in) :: dx
      real(prcsn), intent(in) :: f(:)
      real(prcsn), intent(inout) :: int_f(:)
      ! The variables which are generated inside the function.
      integer :: k
      real(prcsn) :: h_int, h1_o_3, h4_o_3, raux
      real(prcsn) :: integral(nx_int)
!dir$ assume_aligned f(1):64,int_f(1):64
      
      !---------------------------------------------------------------------/
      ! Interpolation.
      if ( nx /= nx_int ) then
        call interpol_quadratic_1D( f, nx, dx, int_f, nx_int )
      else
        int_f = f
      end if
      
      !---------------------------------------------------------------------/
      ! Init. parameters.
      h_int = dx*real(nx-1,prcsn)/real(nx_int-1,prcsn)
      h1_o_3 = h_int/3._prcsn
      h4_o_3 = 4._prcsn*h_int/3._prcsn
      
      !---------------------------------------------------------------------/
      ! Zero out.
      integral = 0._prcsn
      
      !---------------------------------------------------------------------/
      ! Ramping on: Assuming function=0 before the initial point.
      
      ! Interpolate integral at point 1.
      raux = 0._prcsn
!       raux = integral(-1)
!       raux = raux + h1_o_3*int_f(-1)
!       raux = raux + h4_o_3*int_f(0)
      raux = raux + h1_o_3*int_f(1)
      integral(1) = raux
      
      ! Interpolate integral at point 2.
      raux = 0._prcsn
!       raux = integral(0)
!       raux = raux + h1_o_3*int_f(0)
      raux = raux + h4_o_3*int_f(1)
      raux = raux + h1_o_3*int_f(2)
      integral(2) = raux
      
      !---------------------------------------------------------------------/
      ! Integrate 4,6,...,nx_int-0/1.
      do k=2,nx_int-2,2
        raux = integral(k)
        raux = raux + h1_o_3*int_f(k)
        raux = raux + h4_o_3*int_f(k+1)
        raux = raux + h1_o_3*int_f(k+2)
        integral(k+2) = raux
      end do
      
      ! Integrate 3,5,...,nx_int-1/0.
      do k=1,nx_int-2,2
        raux = integral(k)
        raux = raux + h1_o_3*int_f(k)
        raux = raux + h4_o_3*int_f(k+1)
        raux = raux + h1_o_3*int_f(k+2)
        integral(k+2) = raux
      end do
      
      !---------------------------------------------------------------------/
      ! Copy.
      int_f = integral
      
    end subroutine integrate_3o





    !*********************************************************************/
    !*********************************************************************/
    !*****     Real     **************************************************/
    !*********************************************************************/
    !*********************************************************************/
    pure subroutine integrate_3o_r( f, nx, dx, int_f, nx_int )
      use m_dag_data_kind
      use m_dag_interpol
      implicit none
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r4
      integer, intent(in) :: nx, nx_int
      real(prcsn), intent(in) :: dx
      real(prcsn), intent(in) :: f(:)
      real(prcsn), intent(inout) :: int_f(:)
      ! The variables which are generated inside the function.
      integer :: k
      real(prcsn) :: h_int, h1_o_3, h4_o_3, raux
      real(prcsn) :: integral(nx_int)
!dir$ assume_aligned f(1):32,int_f(1):32
      
      !---------------------------------------------------------------------/
      ! Interpolation.
      if ( nx /= nx_int ) then
        call interpol_quadratic_1D( f, nx, dx, int_f, nx_int )
      else
        int_f = f
      end if
      
      !---------------------------------------------------------------------/
      ! Init. parameters.
      h_int = dx*real(nx-1,prcsn)/real(nx_int-1,prcsn)
      h1_o_3 = h_int/3._prcsn
      h4_o_3 = 4._prcsn*h_int/3._prcsn
      
      !---------------------------------------------------------------------/
      ! Zero out.
      integral = 0._prcsn
      
      !---------------------------------------------------------------------/
      ! Ramping on: Assuming function=0 before the initial point.
      
      ! Interpolate integral at point 1.
      raux = 0._prcsn
!       raux = integral(-1)
!       raux = raux + h1_o_3*int_f(-1)
!       raux = raux + h4_o_3*int_f(0)
      raux = raux + h1_o_3*int_f(1)
      integral(1) = raux
      
      ! Interpolate integral at point 2.
      raux = 0._prcsn
!       raux = integral(0)
!       raux = raux + h1_o_3*int_f(0)
      raux = raux + h4_o_3*int_f(1)
      raux = raux + h1_o_3*int_f(2)
      integral(2) = raux
      
      !---------------------------------------------------------------------/
      ! Integrate 4,6,...,nx_int-0/1.
      do k=2,nx_int-2,2
        raux = integral(k)
        raux = raux + h1_o_3*int_f(k)
        raux = raux + h4_o_3*int_f(k+1)
        raux = raux + h1_o_3*int_f(k+2)
        integral(k+2) = raux
      end do
      
      ! Integrate 3,5,...,nx_int-1/0.
      do k=1,nx_int-2,2
        raux = integral(k)
        raux = raux + h1_o_3*int_f(k)
        raux = raux + h4_o_3*int_f(k+1)
        raux = raux + h1_o_3*int_f(k+2)
        integral(k+2) = raux
      end do
      
      !---------------------------------------------------------------------/
      ! Copy.
      int_f = integral
      
    end subroutine integrate_3o_r
    
  !***********************************************************************/
  end module m_dag_integrate