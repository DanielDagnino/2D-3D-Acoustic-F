  !*********************************************************************/
  !** This code has been done in the Barcelona Center for Subsurface 
  !** Imaging (BCSI).
  !** Authors: Daniel Dagnino.
  !*********************************************************************/
  module m_dag_interpol_0
    implicit none
    
    contains
    !*********************************************************************/
    !*********************************************************************/
    !*****     Double     ************************************************/
    !*********************************************************************/
    !*********************************************************************/
    elemental real(r8) function interpol_1o( p, f0, f1 )
      use m_dag_data_kind
      implicit none
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r8
      real(prcsn), intent(in) :: p
      real(prcsn), intent(in) :: f0, f1      
      
      interpol_1o = (1._prcsn-p)*f0 + p*f1
      
    end function interpol_1o
    
    !*********************************************************************/
    elemental real(r8) function interpol_2o( p, fm1, f0, f1, f2 )
      use m_dag_data_kind
      implicit none
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r8
      real(prcsn), intent(in) :: p
      real(prcsn), intent(in) :: fm1, f0, f1, f2
      ! The variables which are generated inside the function.
      real(prcsn) ::  A, B, C, D
      real(prcsn), parameter :: o_o_6 = 1._prcsn/6._prcsn
      
      ! Abramowitz and Stegun. Pg 879
      A = -p*(p-1._prcsn)*(p-2._prcsn)
      B =  (p*p-1._prcsn)*(p-2._prcsn)
      C = -p*(p+1._prcsn)*(p-2._prcsn)
      D =  p*(p*p-1._prcsn)
      
      interpol_2o = o_o_6*( A*fm1 + D*f2 ) + 0.5_prcsn*( B*f0 + C*f1 )
      
    end function interpol_2o
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*****     Real     **************************************************/
    !*********************************************************************/
    !*********************************************************************/
    elemental real(r4) function interpol_1o_r( p, f0, f1 )
      use m_dag_data_kind
      implicit none
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r4
      real(prcsn), intent(in) :: p
      real(prcsn), intent(in) :: f0, f1      
      
      interpol_1o_r = (1._prcsn-p)*f0 + p*f1
      
    end function interpol_1o_r
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    elemental real(r4) function interpol_2o_r( p, fm1, f0, f1, f2 )
      use m_dag_data_kind
      implicit none
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r4
      real(prcsn), intent(in) :: p
      real(prcsn), intent(in) :: fm1, f0, f1, f2
      ! The variables which are generated inside the function.
      real(prcsn) ::  A, B, C, D
      real(prcsn), parameter :: o_o_6 = 1._prcsn/6._prcsn
      
      ! Abramowitz and Stegun. Pg 879
      A = -p*(p-1._prcsn)*(p-2._prcsn)
      B =  (p*p-1._prcsn)*(p-2._prcsn)
      C = -p*(p+1._prcsn)*(p-2._prcsn)
      D =  p*(p*p-1._prcsn)
      
      interpol_2o_r = o_o_6*( A*fm1 + D*f2 ) + 0.5_prcsn*( B*f0 + C*f1 )
      
    end function interpol_2o_r
    
  !***********************************************************************/
  end module m_dag_interpol_0
  
  
  
  
  