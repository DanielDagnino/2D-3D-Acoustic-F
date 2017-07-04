  !*********************************************************************/
  !** This code has been done in the Barcelona Center for Subsurface 
  !** Imaging (BCSI).
  !** Authors: Daniel Dagnino.
  !*********************************************************************/
  module m_dag_interpol_dim_2
    use m_dag_interpol_0
    implicit none
    
    contains
    !*********************************************************************/
    !*********************************************************************/
    !*****     Double     ************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine interpol_linear_2D_3( i_dim, iic, ic, func_ni, na_ni, nb_ni, da_ni, db_ni, func, na, nb, da, db, bound )
      use m_dag_data_kind
      implicit none
      
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r8
      integer, intent(in) :: i_dim
      integer, target, intent(in) :: iic, ic
      integer, intent(in)   :: na, nb, na_ni, nb_ni
      real(prcsn), intent(in)  :: da_ni, db_ni
      real(prcsn), intent(in)  :: func_ni(:,:,:)
      real(prcsn), intent(out) :: func(:,:,:)
      real(prcsn), optional, intent(in) :: da, db
      integer, optional, intent(in)  :: bound
      
      ! The variables which are generated inside the function.
      integer, target :: ia, ib
      integer, target :: iia, iib
      integer, pointer :: ix, iy, iz
      integer, pointer :: iix, iiy, iiz
      
      real(prcsn) :: a, b, a_ni, a_ni1, b_ni, b_ni1, raux, size_a, size_b, value, da_aux, db_aux
      integer :: bound_aux
      integer :: ina_ni, inb_ni
      logical :: fixed
!dir$ assume_aligned func_ni(1,1,1):64,func(1,1,1):64
      
      !---------------------------------------------------------------------/
      ! 
      if ( i_dim == 12 ) then
        iix => iia
        iiy => iib
        iiz => iic
        
        ix => ia
        iy => ib
        iz => ic
      end if
      
      if ( i_dim == 13 ) then
        iix => iia
        iiy => iic
        iiz => iib
        
        ix => ia
        iy => ic
        iz => ib
      end if
      
      if ( i_dim == 23 ) then
        iix => iic
        iiy => iia
        iiz => iib
        
        ix => ic
        iy => ia
        iz => ib
      end if
      
      !---------------------------------------------------------------------/
      ! 
      fixed = .true.
      if ( na/=na_ni .or. nb/=nb_ni ) fixed = .false.
      
      if( present(da) )then
        if ( da*real(na-1,prcsn)/=da_ni*real(na_ni-1,prcsn) .or. &
             db*real(nb-1,prcsn)/=db_ni*real(nb_ni-1,prcsn) ) fixed = .false.
        da_aux = da
        db_aux = db
      else
        da_aux = (da_ni*real(na_ni-1,prcsn))/real(na-1,prcsn)
        db_aux = (db_ni*real(nb_ni-1,prcsn))/real(nb-1,prcsn)
      end if
      
      !---------------------------------------------------------------------/
      ! 
      if( present(bound) )then
        bound_aux = bound
        size_a = da_ni*real(na_ni-1,prcsn)
        size_b = db_ni*real(nb_ni-1,prcsn)
      else
        bound_aux = 0
        size_a = huge(da_aux)
        size_b = huge(db_aux)
      end if
      
      !---------------------------------------------------------------------/
      if ( .not. fixed ) then
        
        raux = 1._prcsn/(da_ni*db_ni)
        do ia=1,na
          a = min( da_aux*real(ia-1,prcsn), size_a )
          ina_ni = floor( a/da_ni ) + 1
          if ( ina_ni == 0 ) ina_ni = 1
          if ( ina_ni >= na_ni ) ina_ni = na_ni-1
          a_ni  = real(ina_ni-1,prcsn)*da_ni
          a_ni1 = a_ni + da_ni
          
          do ib=1,nb
            b = min( db_aux*real(ib-1,prcsn), size_b )
            inb_ni = floor( b/db_ni ) + 1
            if ( inb_ni == 0 ) inb_ni = 1
            if ( inb_ni >= nb_ni ) inb_ni = nb_ni-1
            b_ni  = real(inb_ni-1,prcsn)*db_ni
            b_ni1 = b_ni + db_ni
            
            !-------------------------------/
            if ( da_aux*real(ia-1,prcsn)<size_a .and. db_aux*real(ib-1,prcsn)<size_b .or. ( bound_aux == 1 ) ) then
              
              iia = ina_ni
              iib = inb_ni
              value = (a_ni1-a)*(b_ni1-b)*func_ni(iix,iiy,iiz)
              
              iia = ina_ni+1
              iib = inb_ni
              value = value + (a-a_ni)*(b_ni1-b)*func_ni(iix,iiy,iiz)
              
              iia = ina_ni
              iib = inb_ni+1
              value = value + (a_ni1-a)*(b-b_ni)*func_ni(iix,iiy,iiz)
              
              iia = ina_ni+1
              iib = inb_ni+1
              value = value + (a-a_ni)*(b-b_ni)*func_ni(iix,iiy,iiz)
              
              func(ix,iy,iz) = raux*value
              
            else
              
              if ( bound == 2 ) func(ix,iy,iz) = 0._prcsn
              
            end if
            !-------------------------------/
              
          end do
          
        end do
        
      !---------------------------------------------------------------------/
      else
      
        if ( i_dim == 12 ) then
          func(:,:,iz) = func_ni(:,:,iiz)
        else if ( i_dim == 13 ) then
          func(:,iy,:) = func_ni(:,iiy,:)
        else if ( i_dim == 23 ) then
          func(ix,:,:) = func_ni(iix,:,:)
        end if
        
      end if
      
    end subroutine interpol_linear_2D_3





    !*********************************************************************/
    !*********************************************************************/
    !*****     Real     **************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine interpol_linear_2D_r_3( i_dim, iic, ic, func_ni, na_ni, nb_ni, da_ni, db_ni, func, na, nb, da, db, bound )
      use m_dag_data_kind
      implicit none
      
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r4
      integer, intent(in) :: i_dim
      integer, target, intent(in) :: iic, ic
      integer, intent(in)   :: na, nb, na_ni, nb_ni
      real(prcsn), intent(in)  :: da_ni, db_ni
      real(prcsn), intent(in)  :: func_ni(:,:,:)
      real(prcsn), intent(out) :: func(:,:,:)
      real(prcsn), optional, intent(in) :: da, db
      integer, optional, intent(in)  :: bound
      
      ! The variables which are generated inside the function.
      integer, target :: ia, ib
      integer, target :: iia, iib
      integer, pointer :: ix, iy, iz
      integer, pointer :: iix, iiy, iiz
      
      real(prcsn) :: a, b, a_ni, a_ni1, b_ni, b_ni1, raux, size_a, size_b, value, da_aux, db_aux
      integer :: bound_aux
      integer :: ina_ni, inb_ni
      logical :: fixed
!dir$ assume_aligned func_ni(1,1,1):32,func(1,1,1):32
      
      !---------------------------------------------------------------------/
      ! 
      if ( i_dim == 12 ) then
        iix => iia
        iiy => iib
        iiz => iic
        
        ix => ia
        iy => ib
        iz => ic
      end if
      
      if ( i_dim == 13 ) then
        iix => iia
        iiy => iic
        iiz => iib
        
        ix => ia
        iy => ic
        iz => ib
      end if
      
      if ( i_dim == 23 ) then
        iix => iic
        iiy => iia
        iiz => iib
        
        ix => ic
        iy => ia
        iz => ib
      end if
      
      !---------------------------------------------------------------------/
      ! 
      fixed = .true.
      if ( na/=na_ni .or. nb/=nb_ni ) fixed = .false.
      
      if( present(da) )then
        if ( da*real(na-1,prcsn)/=da_ni*real(na_ni-1,prcsn) .or. &
             db*real(nb-1,prcsn)/=db_ni*real(nb_ni-1,prcsn) ) fixed = .false.
        da_aux = da
        db_aux = db
      else
        da_aux = (da_ni*real(na_ni-1,prcsn))/real(na-1,prcsn)
        db_aux = (db_ni*real(nb_ni-1,prcsn))/real(nb-1,prcsn)
      end if
      
      !---------------------------------------------------------------------/
      ! 
      if( present(bound) )then
        bound_aux = bound
        size_a = da_ni*real(na_ni-1,prcsn)
        size_b = db_ni*real(nb_ni-1,prcsn)
      else
        bound_aux = 0
        size_a = huge(da_aux)
        size_b = huge(db_aux)
      end if
      
      !---------------------------------------------------------------------/
      if ( .not. fixed ) then
        
        raux = 1._prcsn/(da_ni*db_ni)
        do ia=1,na
          a = min( da_aux*real(ia-1,prcsn), size_a )
          ina_ni = floor( a/da_ni ) + 1
          if ( ina_ni == 0 ) ina_ni = 1
          if ( ina_ni >= na_ni ) ina_ni = na_ni-1
          a_ni  = real(ina_ni-1,prcsn)*da_ni
          a_ni1 = a_ni + da_ni
          
          do ib=1,nb
            b = min( db_aux*real(ib-1,prcsn), size_b )
            inb_ni = floor( b/db_ni ) + 1
            if ( inb_ni == 0 ) inb_ni = 1
            if ( inb_ni >= nb_ni ) inb_ni = nb_ni-1
            b_ni  = real(inb_ni-1,prcsn)*db_ni
            b_ni1 = b_ni + db_ni
            
            !-------------------------------/
            if ( da_aux*real(ia-1,prcsn)<size_a .and. db_aux*real(ib-1,prcsn)<size_b .or. ( bound_aux == 1 ) ) then
              
              iia = ina_ni
              iib = inb_ni
              value = (a_ni1-a)*(b_ni1-b)*func_ni(iix,iiy,iiz)
              
              iia = ina_ni+1
              iib = inb_ni
              value = value + (a-a_ni)*(b_ni1-b)*func_ni(iix,iiy,iiz)
              
              iia = ina_ni
              iib = inb_ni+1
              value = value + (a_ni1-a)*(b-b_ni)*func_ni(iix,iiy,iiz)
              
              iia = ina_ni+1
              iib = inb_ni+1
              value = value + (a-a_ni)*(b-b_ni)*func_ni(iix,iiy,iiz)
              
              func(ix,iy,iz) = raux*value
              
            else
              
              if ( bound == 2 ) func(ix,iy,iz) = 0._prcsn
              
            end if
            !-------------------------------/
              
          end do
          
        end do
        
      !---------------------------------------------------------------------/
      else
      
        if ( i_dim == 12 ) then
          func(:,:,iz) = func_ni(:,:,iiz)
        else if ( i_dim == 13 ) then
          func(:,iy,:) = func_ni(:,iiy,:)
        else if ( i_dim == 23 ) then
          func(ix,:,:) = func_ni(iix,:,:)
        end if
        
      end if
      
    end subroutine interpol_linear_2D_r_3
    
    
    
  !***********************************************************************/
  end module m_dag_interpol_dim_2
  
  
  
  
  