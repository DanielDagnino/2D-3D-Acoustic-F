  !*********************************************************************/
  !** This code has been done in the Barcelona Center for Subsurface 
  !** Imaging (BCSI).
  !** Authors: Daniel Dagnino.
  !*********************************************************************/
  module m_dag_interpol_3
    use m_dag_interpol_0
    implicit none
    
    contains
    
    !*********************************************************************/
    !*********************************************************************/
    !*****     Double     ************************************************/
    !*********************************************************************/
    !*********************************************************************/
    pure subroutine interpol_linear_3D( func_ni, nx_ni, ny_ni, nz_ni, dx_ni, dy_ni, dz_ni,  &
                                        func, nx, ny, nz, dx, dy, dz, bound )
      
      use m_dag_data_kind
      implicit none
      
      ! The variables which are passed to the function.
      integer, parameter  :: prcsn = r8
      integer, intent(in) :: nx_ni, ny_ni, nz_ni
      real(prcsn), intent(in)  :: dx_ni, dy_ni, dz_ni
      real(prcsn), intent(in)  :: func_ni(:,:,:)
      integer, intent(in)      :: nx, nz, ny
      real(prcsn), intent(out) :: func(:,:,:)
      real(prcsn), optional, intent(in) :: dx, dy, dz
      integer, optional, intent(in)     :: bound      
      ! The variables which are generated inside the function.
      integer :: i, j, k
      integer :: ix, iy, iz
      real(prcsn) :: x, y, z
      real(prcsn) :: x_ni, y_ni, z_ni
      real(prcsn) :: px, py, pz
      real(prcsn) :: size_x, size_y, size_z
      real(prcsn) :: dx_aux, dy_aux, dz_aux
      integer :: bound_aux
      logical :: fixed
!dir$ assume_aligned func_ni(1,1,1):64,func(1,1,1):64
      
      !*********************************************************************/
      ! 
      fixed = .true.
      if ( nx/=nx_ni .or. ny/=ny_ni .or. nz/=nz_ni ) fixed = .false.
      
      if( present(dx) )then
        if ( dx*real(nx-1)/=dx_ni*real(nx_ni-1) .or. &
             dy*real(ny-1)/=dy_ni*real(ny_ni-1) .or. &
             dz*real(nz-1)/=dz_ni*real(nz_ni-1) ) fixed = .false.
        dx_aux = dx
        dy_aux = dy
        dz_aux = dz
      else
        dx_aux = (dx_ni*real(nx_ni-1,prcsn))/real(nx-1,prcsn)
        dy_aux = (dy_ni*real(ny_ni-1,prcsn))/real(ny-1,prcsn)
        dz_aux = (dz_ni*real(nz_ni-1,prcsn))/real(nz-1,prcsn)
      end if
      
      !*********************************************************************/
      ! 
      if( present(bound) )then
        bound_aux = bound
        size_x = dx_ni*real(nx_ni-1,prcsn)
        size_y = dy_ni*real(ny_ni-1,prcsn)
        size_z = dz_ni*real(nz_ni-1,prcsn)
      else
        bound_aux = 0
        size_x = huge(dx_aux)
        size_y = huge(dy_aux)
        size_z = huge(dz_aux)
      end if
      
      !*********************************************************************/
      ! 
      if ( .not. fixed ) then
        
        !---------------------------------------------------------------------/
        ! 
        do i=1,nx
          x = min( dx_aux*real(i-1,prcsn), size_x )
          ix = floor( x/dx_ni ) + 1
          if ( ix == 0 ) then
            ix = 1
          else if ( ix >= nx_ni ) then
            ix = nx_ni-1
          end if
          x_ni = dx_ni*real(ix-1,prcsn)
          px = (x - x_ni)/dx_ni
          
          do j=1,ny
            y = min( dy_aux*real(j-1,prcsn), size_y )
            iy = floor( y/dy_ni ) + 1
            if ( iy == 0 ) then
              iy = 1
            else if ( iy >= ny_ni ) then
              iy = ny_ni-1
            end if
            y_ni = dy_ni*real(iy-1,prcsn)
            py = (y - y_ni)/dy_ni
            
            do k=1,nz
              z = min( dz_aux*real(k-1,prcsn), size_z )
              iz = floor( z/dz_ni ) + 1
              if ( iz == 0 ) then
                iz = 1
              else if ( iz >= nz_ni ) then
                iz = nz_ni-1
              end if
              z_ni = dz_ni*real(iz-1,prcsn)
              pz = (z - z_ni)/dz_ni
              
              !---------------------------------------------------------------------/      
              ! 
              func(i,j,k) = &
                (1._prcsn-px)*(1._prcsn-py)*(1._prcsn-pz) * func_ni(ix  ,iy  ,iz  ) + &
                          px *(1._prcsn-py)*(1._prcsn-pz) * func_ni(ix+1,iy  ,iz  ) + &
                (1._prcsn-px)*          py *(1._prcsn-pz) * func_ni(ix  ,iy+1,iz  ) + &
                          px *          py *(1._prcsn-pz) * func_ni(ix+1,iy+1,iz  ) + &
                (1._prcsn-px)*(1._prcsn-py)*          pz  * func_ni(ix  ,iy  ,iz+1) + &
                          px *(1._prcsn-py)*          pz  * func_ni(ix+1,iy  ,iz+1) + &
                (1._prcsn-px)*          py *          pz  * func_ni(ix  ,iy+1,iz+1) + &
                          px *          py *          pz  * func_ni(ix+1,iy+1,iz+1)
              !---------------------------------------------------------------------/ 
              
            end do
          end do
        end do
        
      !*********************************************************************/
      else
        func = func_ni
      end if
      
      !*********************************************************************/
      
    end subroutine interpol_linear_3D
    
    
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*****     Real     **************************************************/
    !*********************************************************************/
    !*********************************************************************/
    pure subroutine interpol_linear_3D_r( func_ni, nx_ni, ny_ni, nz_ni, dx_ni, dy_ni, dz_ni,  &
                                        func, nx, ny, nz, dx, dy, dz, bound )
      
      use m_dag_data_kind
      implicit none
      
      ! The variables which are passed to the function.
      integer, parameter  :: prcsn = r4
      integer, intent(in) :: nx_ni, ny_ni, nz_ni
      real(prcsn), intent(in)  :: dx_ni, dy_ni, dz_ni
      real(prcsn), intent(in)  :: func_ni(:,:,:)
      integer, intent(in)      :: nx, nz, ny
      real(prcsn), intent(out) :: func(:,:,:)
      real(prcsn), optional, intent(in) :: dx, dy, dz
      integer, optional, intent(in)     :: bound      
      ! The variables which are generated inside the function.
      integer :: i, j, k
      integer :: ix, iy, iz
      real(prcsn) :: x, y, z
      real(prcsn) :: x_ni, y_ni, z_ni
      real(prcsn) :: px, py, pz
      real(prcsn) :: size_x, size_y, size_z
      real(prcsn) :: dx_aux, dy_aux, dz_aux
      integer :: bound_aux
      logical :: fixed
!dir$ assume_aligned func_ni(1,1,1):32,func(1,1,1):32
      
      !*********************************************************************/
      ! 
      fixed = .true.
      if ( nx/=nx_ni .or. ny/=ny_ni .or. nz/=nz_ni ) fixed = .false.
      
      if( present(dx) )then
        if ( dx*real(nx-1)/=dx_ni*real(nx_ni-1) .or. &
             dy*real(ny-1)/=dy_ni*real(ny_ni-1) .or. &
             dz*real(nz-1)/=dz_ni*real(nz_ni-1) ) fixed = .false.
        dx_aux = dx
        dy_aux = dy
        dz_aux = dz
      else
        dx_aux = (dx_ni*real(nx_ni-1,prcsn))/real(nx-1,prcsn)
        dy_aux = (dy_ni*real(ny_ni-1,prcsn))/real(ny-1,prcsn)
        dz_aux = (dz_ni*real(nz_ni-1,prcsn))/real(nz-1,prcsn)
      end if
      
      !*********************************************************************/
      ! 
      if( present(bound) )then
        bound_aux = bound
        size_x = dx_ni*real(nx_ni-1,prcsn)
        size_y = dy_ni*real(ny_ni-1,prcsn)
        size_z = dz_ni*real(nz_ni-1,prcsn)
      else
        bound_aux = 0
        size_x = huge(dx_aux)
        size_y = huge(dy_aux)
        size_z = huge(dz_aux)
      end if
      
      !*********************************************************************/
      ! 
      if ( .not. fixed ) then
        
        !---------------------------------------------------------------------/
        ! 
        do i=1,nx
          x = min( dx_aux*real(i-1,prcsn), size_x )
          ix = floor( x/dx_ni ) + 1
          if ( ix == 0 ) then
            ix = 1
          else if ( ix >= nx_ni ) then
            ix = nx_ni-1
          end if
          x_ni = dx_ni*real(ix-1,prcsn)
          px = (x - x_ni)/dx_ni
          
          do j=1,ny
            y = min( dy_aux*real(j-1,prcsn), size_y )
            iy = floor( y/dy_ni ) + 1
            if ( iy == 0 ) then
              iy = 1
            else if ( iy >= ny_ni ) then
              iy = ny_ni-1
            end if
            y_ni = dy_ni*real(iy-1,prcsn)
            py = (y - y_ni)/dy_ni
            
            do k=1,nz
              z = min( dz_aux*real(k-1,prcsn), size_z )
              iz = floor( z/dz_ni ) + 1
              if ( iz == 0 ) then
                iz = 1
              else if ( iz >= nz_ni ) then
                iz = nz_ni-1
              end if
              z_ni = dz_ni*real(iz-1,prcsn)
              pz = (z - z_ni)/dz_ni
              
              !---------------------------------------------------------------------/      
              ! 
              func(i,j,k) = &
                (1._prcsn-px)*(1._prcsn-py)*(1._prcsn-pz) * func_ni(ix  ,iy  ,iz  ) + &
                          px *(1._prcsn-py)*(1._prcsn-pz) * func_ni(ix+1,iy  ,iz  ) + &
                (1._prcsn-px)*          py *(1._prcsn-pz) * func_ni(ix  ,iy+1,iz  ) + &
                          px *          py *(1._prcsn-pz) * func_ni(ix+1,iy+1,iz  ) + &
                (1._prcsn-px)*(1._prcsn-py)*          pz  * func_ni(ix  ,iy  ,iz+1) + &
                          px *(1._prcsn-py)*          pz  * func_ni(ix+1,iy  ,iz+1) + &
                (1._prcsn-px)*          py *          pz  * func_ni(ix  ,iy+1,iz+1) + &
                          px *          py *          pz  * func_ni(ix+1,iy+1,iz+1)
              !---------------------------------------------------------------------/ 
              
            end do
          end do
        end do
        
      !*********************************************************************/
      else
        func = func_ni
      end if
      
      !*********************************************************************/
      
    end subroutine interpol_linear_3D_r
    
    
    
    
    
  !***********************************************************************/
  end module m_dag_interpol_3
  
  
  
  
  