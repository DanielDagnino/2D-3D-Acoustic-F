  !*********************************************************************/
  !*****     Author: D. Dagnino     ************************************/
  !*********************************************************************/
  module m_lin_decay
    implicit none
    
    interface lin_decay
      module procedure lin_decay, lin_decay_r
    end interface
    
    contains
    !**********************************************************!
    !**********************************************************!
    !**********************************************************!
    elemental real(r8) function lin_decay( x, a, b, n )
      use m_dag_data_kind
      implicit none
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r8
      real(prcsn), intent(in) :: x, a, b, n
      ! The variables which are generated inside the function.
      real(prcsn) :: x_reg
      
      x_reg = (x-a)/(b-a)
      lin_decay = 0.5_prcsn * ( 1._prcsn + tanh(n*(-0.5_prcsn+x_reg))/tanh(0.5_prcsn*n) )
      
    end function lin_decay
    
    !***********************************************************************/
    elemental real(r4) function lin_decay_r( x, a, b, n )
      use m_dag_data_kind
      implicit none
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r4
      real(prcsn), intent(in) :: x, a, b, n
      ! The variables which are generated inside the function.
      real(prcsn) :: x_reg
      
      x_reg = (x-a)/(b-a)
      lin_decay_r = 0.5_prcsn * ( 1._prcsn + tanh(n*(-0.5_prcsn+x_reg))/tanh(0.5_prcsn*n) )
      
    end function lin_decay_r
    
    !***********************************************************************/
  end module m_lin_decay
  
  
  
  
  