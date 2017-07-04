  
  
  !**********************************************************!
  module m_dag_det3
   implicit none
   
    interface det3
      module procedure det3, det3_r
    end interface
    
    !*********************************************************************/
    !*********************************************************************/
    contains
    !---------------------------------------------------------------------!
    elemental real(r4) function det3_r( a1, a2, a3, b1, b2, b3, c1, c2, c3 )
      use m_dag_data_kind
      implicit none
      ! The variables which are passed to the function.
      real(r4), intent(in)  :: a1, a2, a3, b1, b2, b3, c1, c2, c3
    
  !   | a1 a2 a3 |
  !   | b1 b2 b3 |
  !   | c1 c2 c3 |
    
      det3_r = a1*b2*c3 + b1*c2*a3 + c1*a2*b3 - ( c1*b2*a3 + c2*b3*a1 + c3*b1*a2 )
    
    end function det3_r
    
    !---------------------------------------------------------------------!
    elemental real(r8) function det3( a1, a2, a3, b1, b2, b3, c1, c2, c3 )
      use m_dag_data_kind
      implicit none
      ! The variables which are passed to the function.
      real(r8), intent(in)  :: a1, a2, a3, b1, b2, b3, c1, c2, c3
    
  !   | a1 a2 a3 |
  !   | b1 b2 b3 |
  !   | c1 c2 c3 |
    
      det3 = a1*b2*c3 + b1*c2*a3 + c1*a2*b3 - ( c1*b2*a3 + c2*b3*a1 + c3*b1*a2 )
    
    end function det3
    
    !---------------------------------------------------------------------!
  end module m_dag_det3
  
  
  
  
  