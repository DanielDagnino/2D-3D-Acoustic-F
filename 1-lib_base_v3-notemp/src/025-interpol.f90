  !*********************************************************************/
  !** This code has been done in the Barcelona Center for Subsurface 
  !** Imaging (BCSI).
  !** Authors: Daniel Dagnino.
  !*********************************************************************/
  module m_dag_interpol_dim
    use m_dag_interpol_0
    implicit none
    
    contains
    
    !*********************************************************************/
    !*********************************************************************/
    !*****     Double     ************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine interpol_linear_1D_2( i_dim, iib, ib, func_ni, na_ni, da_ni, func, na, da, bound )
      use m_dag_data_kind
      implicit none
      
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r8
      integer, intent(in) :: i_dim
      integer, target, intent(in) :: iib, ib
      integer, intent(in)   :: na, na_ni
      real(prcsn), intent(in)  :: da_ni
      real(prcsn), intent(in)  :: func_ni(:,:)
      real(prcsn), intent(out) :: func(:,:)
      real(prcsn), optional, intent(in) :: da
      integer, optional, intent(in)  :: bound
      
      ! The variables which are generated inside the function.
      integer, target :: ia
      integer, target :: iia
      integer, pointer :: ix, iy
      integer, pointer :: iix, iiy
      
      integer :: bound_aux
      integer :: ina_ni
      real(prcsn) :: a, a_ni, size_a, value1, value2, p, da_aux
      logical :: fixed
!dir$ assume_aligned func_ni(1,1):64,func(1,1):64
      
      !---------------------------------------------------------------------/
      ! 
      if ( i_dim == 1 ) then
        iix => iia
        iiy => iib
        
        ix => ia
        iy => ib
      end if
      
      if ( i_dim == 2 ) then
        iix => iib
        iiy => iia
        
        ix => ib
        iy => ia
      end if
      
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
            
            iia = ina_ni
            value1 = func_ni(iix,iiy)
            
            iia = ina_ni+1
            value2 = func_ni(iix,iiy)
            
            func(ix,iy) = interpol_1o( p, value1, value2 )
          else
            if ( bound == 2 ) func(ix,iy) = 0._prcsn
          end if
        end do
        
      else
        
        if ( i_dim == 1 ) then
          func(:,iy) = func_ni(:,iiy)
        else if ( i_dim == 2 ) then
          func(ix,:) = func_ni(iix,:)
        end if
        
      end if
      
    end subroutine interpol_linear_1D_2
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine interpol_quadratic_1D_2( i_dim, iib, ib, func_ni, na_ni, da_ni, func, na, da, bound )
      use m_dag_data_kind
      implicit none
      
      ! The variables whibh are passed to the function.
      integer, parameter :: prcsn = r8
      integer, intent(in) :: i_dim
      integer, target, intent(in) :: iib, ib
      integer, intent(in)   :: na, na_ni
      real(prcsn), intent(in)  :: da_ni
      real(prcsn), intent(in)  :: func_ni(:,:)
      real(prcsn), intent(out) :: func(:,:)
      real(prcsn), optional, intent(in) :: da
      integer, optional, intent(in)  :: bound
      
      ! The variables which are generated inside the function.
      integer, target :: ia
      integer, target :: iia
      integer, pointer :: ix, iy
      integer, pointer :: iix, iiy
      
      integer :: bound_aux
      integer :: ina_ni
      real(prcsn) :: a, a_ni, size_a, value1, value2, value3, value4, p, da_aux
      logical :: fixed
!dir$ assume_aligned func_ni(1,1):64,func(1,1):64
      
      !---------------------------------------------------------------------/
      ! 
      if ( i_dim == 1 ) then
        iix => iia
        iiy => iib
        
        ix => ia
        iy => ib
      end if
      
      if ( i_dim == 2 ) then
        iix => iib
        iiy => iia
        
        ix => ib
        iy => ia
      end if
      
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
          if ( ina_ni == 0 .or. ina_ni == 1 ) then
            ina_ni = 2
          else if ( ina_ni >= na_ni-1 ) then
            ina_ni = ina_ni-2
          end if
          a_ni  = real(ina_ni-1,prcsn)*da_ni
          
          if ( da_aux*real(ia-1,prcsn)<size_a .or. ( bound_aux == 1 ) ) then
            ! Abramowitz and Stegun. Pg 879
            p = (a - a_ni)/da_ni
            
            iia = ina_ni-1
            value1 = func_ni(iix,iiy)
            
            iia = ina_ni
            value2 = func_ni(iix,iiy)
            
            iia = ina_ni+1
            value3 = func_ni(iix,iiy)
            
            iia = ina_ni+2
            value4 = func_ni(iix,iiy)
            
            func(ix,iy) = interpol_2o( p, value1, value2, value3, value4 )
          else
            if ( bound == 2 ) func(ix,iy) = 0._prcsn
          end if
        end do
        
      else
      
        if ( i_dim == 1 ) then
          func(:,iy) = func_ni(:,iiy)
        else if ( i_dim == 2 ) then
          func(ix,:) = func_ni(iix,:)
        end if
      
      end if
      
    end subroutine interpol_quadratic_1D_2
    
    
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*****     Real     **************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine interpol_linear_1D_r_2( i_dim, iib, ib, func_ni, na_ni, da_ni, func, na, da, bound )
      use m_dag_data_kind
      implicit none
      
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r4
      integer, intent(in) :: i_dim
      integer, target, intent(in) :: iib, ib
      integer, intent(in)   :: na, na_ni
      real(prcsn), intent(in)  :: da_ni
      real(prcsn), intent(in)  :: func_ni(:,:)
      real(prcsn), intent(out) :: func(:,:)
      real(prcsn), optional, intent(in) :: da
      integer, optional, intent(in)  :: bound
      
      ! The variables which are generated inside the function.
      integer, target :: ia
      integer, target :: iia
      integer, pointer :: ix, iy
      integer, pointer :: iix, iiy
      
      integer :: bound_aux
      integer :: ina_ni
      real(prcsn) :: a, a_ni, size_a, value1, value2, p, da_aux
      logical :: fixed
!dir$ assume_aligned func_ni(1,1):32,func(1,1):32
      
      !---------------------------------------------------------------------/
      ! 
      if ( i_dim == 1 ) then
        iix => iia
        iiy => iib
        
        ix => ia
        iy => ib
      end if
      
      if ( i_dim == 2 ) then
        iix => iib
        iiy => iia
        
        ix => ib
        iy => ia
      end if
      
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
            
            iia = ina_ni
            value1 = func_ni(iix,iiy)
            
            iia = ina_ni+1
            value2 = func_ni(iix,iiy)
            
            func(ix,iy) = interpol_1o_r( p, value1, value2 )
          else
            if ( bound == 2 ) func(ix,iy) = 0._prcsn
          end if
        end do
        
      else
        
        if ( i_dim == 1 ) then
          func(:,iy) = func_ni(:,iiy)
        else if ( i_dim == 2 ) then
          func(ix,:) = func_ni(iix,:)
        end if
        
      end if
      
    end subroutine interpol_linear_1D_r_2
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine interpol_quadratic_1D_r_2( i_dim, iib, ib, func_ni, na_ni, da_ni, func, na, da, bound )
      use m_dag_data_kind
      implicit none
      
      ! The variables whibh are passed to the function.
      integer, parameter :: prcsn = r4
      integer, intent(in) :: i_dim
      integer, target, intent(in) :: iib, ib
      integer, intent(in)   :: na, na_ni
      real(prcsn), intent(in)  :: da_ni
      real(prcsn), intent(in)  :: func_ni(:,:)
      real(prcsn), intent(out) :: func(:,:)
      real(prcsn), optional, intent(in) :: da
      integer, optional, intent(in)  :: bound
      
      ! The variables which are generated inside the function.
      integer, target :: ia
      integer, target :: iia
      integer, pointer :: ix, iy
      integer, pointer :: iix, iiy
      
      integer :: bound_aux
      integer :: ina_ni
      real(prcsn) :: a, a_ni, size_a, value1, value2, value3, value4, p, da_aux
      logical :: fixed
!dir$ assume_aligned func_ni(1,1):32,func(1,1):32
      
      !---------------------------------------------------------------------/
      ! 
      if ( i_dim == 1 ) then
        iix => iia
        iiy => iib
        
        ix => ia
        iy => ib
      end if
      
      if ( i_dim == 2 ) then
        iix => iib
        iiy => iia
        
        ix => ib
        iy => ia
      end if
      
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
          if ( ina_ni == 0 .or. ina_ni == 1 ) then
            ina_ni = 2
          else if ( ina_ni >= na_ni-1 ) then
            ina_ni = ina_ni-2
          end if
          a_ni  = real(ina_ni-1,prcsn)*da_ni
          
          if ( da_aux*real(ia-1,prcsn)<size_a .or. ( bound_aux == 1 ) ) then
            ! Abramowitz and Stegun. Pg 879
            p = (a - a_ni)/da_ni
            
            iia = ina_ni-1
            value1 = func_ni(iix,iiy)
            
            iia = ina_ni
            value2 = func_ni(iix,iiy)
            
            iia = ina_ni+1
            value3 = func_ni(iix,iiy)
            
            iia = ina_ni+2
            value4 = func_ni(iix,iiy)
            
            func(ix,iy) = interpol_2o_r( p, value1, value2, value3, value4 )
          else
            if ( bound == 2 ) func(ix,iy) = 0._prcsn
          end if
        end do
        
      else
      
        if ( i_dim == 1 ) then
          func(:,iy) = func_ni(:,iiy)
        else if ( i_dim == 2 ) then
          func(ix,:) = func_ni(iix,:)
        end if
      
      end if
      
    end subroutine interpol_quadratic_1D_r_2
    
    
    
    
    
  !***********************************************************************/
  end module m_dag_interpol_dim
  
  
  
  
  