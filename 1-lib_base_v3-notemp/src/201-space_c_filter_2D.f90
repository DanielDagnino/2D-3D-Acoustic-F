  !*********************************************************************/
  !*****     Author: D. Dagnino     ************************************/
  !*********************************************************************/
  module m_dag_space_c_filter_2D
    
    interface space_c_filter_2D
      module procedure space_c_filter_2D, space_c_filter_2D_r
    end interface
    
    contains
    !*********************************************************************/
    !*********************************************************************/
    !*****     Double     ************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine space_c_filter_2D( image, lamd_x, lamd_y, nx, ny, dx, dy )
      use m_dag_data_kind
      implicit none
      
      ! The variables which are passed to the function.
      real(r8), intent(in)    :: lamd_x, lamd_y
      integer, intent(in)     :: nx, ny
      real(r8), intent(in)    :: dx, dy
      real(r8), intent(inout) :: image(:,:)
      
      ! The variables which are generated inside the function.
      integer  :: i, j, ix, iy, iix, iiy, rest, ix0, ix1, iy0, iy1
      real(r8) :: raux, norm
      real(r8) :: image_filt(nx,ny)
!dir$ assume_aligned image(1,1):64
      
      !---------------------------------------------------------------------/
      ! Filtering.
      
      iix = anint(0.5_r8*lamd_x/dx)
      iiy = anint(0.5_r8*lamd_y/dy)
      
!       iix = floor(0.5_r8*lamd_x/dx)
!       iiy = floor(0.5_r8*lamd_y/dy)
      
      if ( iix>floor(0.5_r8*dble(nx)/dx) ) iix = floor(0.5_r8*dble(nx)/dx)
      if ( iiy>floor(0.5_r8*dble(ny)/dy) ) iiy = floor(0.5_r8*dble(ny)/dy)
      
      norm = (2._r8*dble(iix)+1._r8)*(2._r8*dble(iiy)+1._r8)
      norm = 1._r8/norm
      
      !---------------------------------------------------------------------/
      if ( iix/=0 .or. iiy/=0 ) then
        
        do ix=1,nx
        do iy=1,ny
            
            ! 
            if ( ix-iix < 1 ) then
              ix0  = 1
              rest = 1 - (ix-iix)
              ix1 = ix+rest+iix
            else if ( ix+iix > nx ) then
              ix1  = nx
              rest = (ix+iix) - nx
              ix0 = ix-rest-iix
            else
              ix0 = ix-iix
              ix1 = ix+iix
            end if
            
              ! 
              if ( iy-iiy < 1 ) then
                iy0  = 1
                rest = 1 - (iy-iiy)
                iy1 = iy+rest+iiy
              else if ( iy+iiy > ny ) then
                iy1  = ny
                rest = (iy+iiy) - ny
                iy0 = iy-rest-iiy
              else
                iy0 = iy-iiy
                iy1 = iy+iiy
              end if
              
              ! 
              raux = 0._r8
              do i=ix0,ix1
              do j=iy0,iy1
                raux = raux + image(i,j)
              end do
              end do
                      
          image_filt(ix,iy) = raux
          
        end do
        end do
        
        image = norm*image_filt
        
      end if
      !---------------------------------------------------------------------/
      
    end subroutine space_c_filter_2D
    
    
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*****     Real     **************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine space_c_filter_2D_r( image, lamd_x, lamd_y, nx, ny, dx, dy )
      use m_dag_data_kind
      implicit none
      
      ! The variables which are passed to the function.
      real(r4), intent(in)    :: lamd_x, lamd_y
      integer, intent(in)     :: nx, ny
      real(r4), intent(in)    :: dx, dy
      real(r4), intent(inout) :: image(:,:)
      
      ! The variables which are generated inside the function.
      integer  :: i, j, ix, iy, iix, iiy, rest, ix0, ix1, iy0, iy1
      real(r4) :: raux, norm
      real(r4) :: image_filt(nx,ny)
!dir$ assume_aligned image(1,1):32
      
      !---------------------------------------------------------------------/
      ! Filtering.
      
      iix = floor(0.5_r4*lamd_x/dx)
      iiy = floor(0.5_r4*lamd_y/dy)
      if ( iix>floor(0.5_r4*real(nx)/dx) ) iix = floor(0.5_r4*real(nx)/dx)
      if ( iiy>floor(0.5_r4*real(ny)/dy) ) iiy = floor(0.5_r4*real(ny)/dy)
      
      norm = (2._r4*real(iix)+1._r4)*(2._r4*real(iiy)+1._r4)
      norm = 1._r4/norm
      
      if ( iix/=0 .or. iiy/=0 ) then
        
        do ix=1,nx
        do iy=1,ny
            
            ! 
            if ( ix-iix < 1 ) then
              ix0  = 1
              rest = 1 - (ix-iix)
              ix1 = ix+rest+iix
            else if ( ix+iix > nx ) then
              ix1  = nx
              rest = (ix+iix) - nx
              ix0 = ix-rest-iix
            else
              ix0 = ix-iix
              ix1 = ix+iix
            end if
            
              ! 
              if ( iy-iiy < 1 ) then
                iy0  = 1
                rest = 1 - (iy-iiy)
                iy1 = iy+rest+iiy
              else if ( iy+iiy > ny ) then
                iy1  = ny
                rest = (iy+iiy) - ny
                iy0 = iy-rest-iiy
              else
                iy0 = iy-iiy
                iy1 = iy+iiy
              end if
              
              ! 
              raux = 0._r4
              do i=ix0,ix1
              do j=iy0,iy1
                raux = raux + image(i,j)
              end do
              end do
!                       
          image_filt(ix,iy) = raux
          
        end do
        end do
        
        image = norm*image_filt
        
      end if
      !---------------------------------------------------------------------/
      
    end subroutine space_c_filter_2D_r
    
  !***********************************************************************/
  end module m_dag_space_c_filter_2D
  
  
  
  