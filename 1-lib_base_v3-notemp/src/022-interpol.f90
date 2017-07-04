  !*********************************************************************/
  !** This code has been done in the Barcelona Center for Subsurface 
  !** Imaging (BCSI).
  !** Authors: Daniel Dagnino.
  !*********************************************************************/
  module m_dag_interpol_2
    use m_dag_interpol_0
    implicit none
    
    contains
    !*********************************************************************/
    !*********************************************************************/
    !*****     Double     ************************************************/
    !*********************************************************************/
    !*********************************************************************/
    pure subroutine interpol_linear_1D( func_ni, na_ni, da_ni, func, na, da, bound )
      use m_dag_data_kind
      implicit none
      
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r8
      integer, intent(in)   :: na, na_ni
      real(prcsn), intent(in)  :: da_ni
      real(prcsn), intent(in)  :: func_ni(:)
      real(prcsn), intent(out) :: func(:)
      real(prcsn), optional, intent(in) :: da
      integer, optional, intent(in)  :: bound      
      ! The variables which are generated inside the function.
      integer :: ia, bound_aux
      integer :: ina_ni
      real(prcsn) :: a, a_ni, size_a, p, da_aux
      logical :: fixed
!dir$ assume_aligned func_ni(1):64,func(1):64
      
      !---------------------------------------------------------------------/
      ! 
      fixed = .true.
      if ( na/=na_ni ) fixed = .false.
      
      if( present(da) )then
        if ( da*real(na-1,prcsn)/=da_ni*real(na_ni-1,prcsn) ) fixed = .false.
        da_aux = da
      else
        da_aux = (da_ni*real(na_ni-1,prcsn))/real(na-1,prcsn)
      end if
      
      !---------------------------------------------------------------------/
      ! 
      if( present(bound) )then
        bound_aux = bound
        size_a = da_ni*real(na_ni-1,prcsn)
      else
        bound_aux = 0
        size_a = huge(da_aux)
      end if
      
      !---------------------------------------------------------------------/
      if ( .not. fixed ) then
        
        do ia=1,na
          a = min( da_aux*real(ia-1,prcsn), size_a )
          ina_ni = floor( a/da_ni ) + 1
          if ( ina_ni == 0 ) ina_ni = 1
          if ( ina_ni >= na_ni ) ina_ni = na_ni-1
          a_ni  = real(ina_ni-1,prcsn)*da_ni
          
          if ( da_aux*real(ia-1,prcsn)<size_a .or. ( bound_aux == 1 ) ) then
            p = (a - a_ni)/da_ni
            func(ia) = interpol_1o( p, func_ni(ina_ni), func_ni(ina_ni+1) )
          else
            if ( bound == 2 ) func(ia) = 0._prcsn
          end if
        end do
        
      else
        func = func_ni
      end if
      
    end subroutine interpol_linear_1D
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    pure subroutine interpol_quadratic_1D( func_ni, na_ni, da_ni, func, na, da, bound )
      use m_dag_data_kind
      implicit none
      
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r8
      integer, intent(in)   :: na, na_ni
      real(prcsn), intent(in)  :: da_ni
      real(prcsn), intent(in)  :: func_ni(:)
      real(prcsn), intent(out) :: func(:)
      real(prcsn), optional, intent(in) :: da
      integer, optional, intent(in)  :: bound      
      ! The variables which are generated inside the function.
      integer :: ia, bound_aux
      integer :: ina_ni
      real(prcsn) :: a, a_ni, size_a, p, da_aux
      logical :: fixed
!dir$ assume_aligned func_ni(1):64,func(1):64
      
      !---------------------------------------------------------------------/
      ! 
      fixed = .true.
      if ( na/=na_ni ) fixed = .false.
      
      if( present(da) )then
        if ( da*real(na-1,prcsn)/=da_ni*real(na_ni-1,prcsn) ) fixed = .false.
        da_aux = da
      else
        da_aux = (da_ni*real(na_ni-1,prcsn))/real(na-1,prcsn)
      end if
      
      !---------------------------------------------------------------------/
      ! 
      if( present(bound) )then
        bound_aux = bound
        size_a = da_ni*real(na_ni-1,prcsn)
      else
        bound_aux = 0
        size_a = huge(da_aux)
      end if
      
      !---------------------------------------------------------------------/
      if ( .not. fixed ) then
        
        do ia=1,na
          a = min( da_aux*real(ia-1,prcsn), size_a )
          ina_ni = floor( a/da_ni ) + 1
          if ( ina_ni <= 1 ) then
            ina_ni = 2
          else if ( ina_ni >= na_ni-1 ) then
            ina_ni = ina_ni-2
          end if
          a_ni  = real(ina_ni-1,prcsn)*da_ni
          
          if ( da_aux*real(ia-1,prcsn)<size_a .or. ( bound_aux == 1 ) ) then
            ! Abramowitz and Stegun. Pg 879
            p = (a - a_ni)/da_ni
            func(ia) = interpol_2o( p, func_ni(ina_ni-1), func_ni(ina_ni), func_ni(ina_ni+1), func_ni(ina_ni+2) )
          else
            if ( bound == 2 ) func(ia) = 0._prcsn
          end if
        end do
        
      else
        func = func_ni
      end if
      
    end subroutine interpol_quadratic_1D
    
    
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*****     Real     **************************************************/
    !*********************************************************************/
    !*********************************************************************/
    pure subroutine interpol_linear_1D_r( func_ni, na_ni, da_ni, func, na, da, bound )
      use m_dag_data_kind
      implicit none
      
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r4
      integer, intent(in)   :: na, na_ni
      real(prcsn), intent(in)  :: da_ni
      real(prcsn), intent(in)  :: func_ni(:)
      real(prcsn), intent(out) :: func(:)
      real(prcsn), optional, intent(in) :: da
      integer, optional, intent(in)  :: bound      
      ! The variables which are generated inside the function.
      integer :: ia, bound_aux
      integer :: ina_ni
      real(prcsn) :: a, a_ni, size_a, p, da_aux
      logical :: fixed
!dir$ assume_aligned func_ni(1):32,func(1):32
      
      !---------------------------------------------------------------------/
      ! 
      fixed = .true.
      if ( na/=na_ni ) fixed = .false.
      
      if( present(da) )then
        if ( da*real(na-1,prcsn)/=da_ni*real(na_ni-1,prcsn) ) fixed = .false.
        da_aux = da
      else
        da_aux = (da_ni*real(na_ni-1,prcsn))/real(na-1,prcsn)
      end if
      
      !---------------------------------------------------------------------/
      ! 
      if( present(bound) )then
        bound_aux = bound
        size_a = da_ni*real(na_ni-1,prcsn)
      else
        bound_aux = 0
        size_a = huge(da_aux)
      end if
      
      !---------------------------------------------------------------------/
      if ( .not. fixed ) then
        
        do ia=1,na
          a = min( da_aux*real(ia-1,prcsn), size_a )
          ina_ni = floor( a/da_ni ) + 1
          if ( ina_ni == 0 ) ina_ni = 1
          if ( ina_ni >= na_ni ) ina_ni = na_ni-1
          a_ni  = real(ina_ni-1,prcsn)*da_ni
          
          if ( da_aux*real(ia-1,prcsn)<size_a .or. ( bound_aux == 1 ) ) then
            p = (a - a_ni)/da_ni
            func(ia) = interpol_1o_r( p, func_ni(ina_ni), func_ni(ina_ni+1) )
          else
            if ( bound == 2 ) func(ia) = 0._prcsn
          end if
        end do
        
      else
        func = func_ni
      end if
      
    end subroutine interpol_linear_1D_r
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    pure subroutine interpol_quadratic_1D_r( func_ni, na_ni, da_ni, func, na, da, bound )
      use m_dag_data_kind
      implicit none
      
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r4
      integer, intent(in)   :: na, na_ni
      real(prcsn), intent(in)  :: da_ni
      real(prcsn), intent(in)  :: func_ni(:)
      real(prcsn), intent(out) :: func(:)
      real(prcsn), optional, intent(in) :: da
      integer, optional, intent(in)  :: bound
      ! The variables which are generated inside the function.
      integer :: ia, bound_aux
      integer :: ina_ni
      real(prcsn) :: a, a_ni, size_a, p, da_aux
      logical :: fixed
!dir$ assume_aligned func_ni(1):32,func(1):32
      
      !---------------------------------------------------------------------/
      ! 
      fixed = .true.
      if ( na/=na_ni ) fixed = .false.
      
      if( present(da) )then
        if ( da*real(na-1,prcsn)/=da_ni*real(na_ni-1,prcsn) ) fixed = .false.
        da_aux = da
      else
        da_aux = (da_ni*real(na_ni-1,prcsn))/real(na-1,prcsn)
      end if
      
      !---------------------------------------------------------------------/
      ! 
      if( present(bound) )then
        bound_aux = bound
        size_a = da_ni*real(na_ni-1,prcsn)
      else
        bound_aux = 0
        size_a = huge(da_aux)
      end if
      
      !---------------------------------------------------------------------/
      if ( .not. fixed ) then
        
        do ia=1,na
          a = min( da_aux*real(ia-1,prcsn), size_a )
          ina_ni = floor( a/da_ni ) + 1
          if ( ina_ni <= 1 ) then
            ina_ni = 2
          else if ( ina_ni >= na_ni-1 ) then
            ina_ni = ina_ni-2
          end if
          a_ni  = real(ina_ni-1,prcsn)*da_ni
          
          if ( da_aux*real(ia-1,prcsn)<size_a .or. ( bound_aux == 1 ) ) then
            ! Abramowitz and Stegun. Pg 879
            p = (a - a_ni)/da_ni
            func(ia) = interpol_2o_r( p, func_ni(ina_ni-1), func_ni(ina_ni), func_ni(ina_ni+1), func_ni(ina_ni+2) )
          else
            if ( bound == 2 ) func(ia) = 0._prcsn
          end if
        end do
        
      else
        func = func_ni
      end if
      
    end subroutine interpol_quadratic_1D_r
    
    
    
    
    
  !***********************************************************************/
  end module m_dag_interpol_2
  
  
  
  
  