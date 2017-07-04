  !*********************************************************************/
  !*****     Author: D. Dagnino     ************************************/
  !*********************************************************************/
  module m_dag_space_2_filter_2D_zp
    
    interface space_2_filter_2D_zp
      module procedure space_2_filter_2D_zp, space_2_filter_2D_zp_r
    end interface
    
    contains
    !*********************************************************************/
    !*********************************************************************/
    !*****     Double     ************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine space_2_filter_2D_zp( image, &
        freq1_filt_x, freq2_filt_x, p0_x, p1_x, type_filt_x, order_filt_x, zero_padding_x, &
        freq1_filt_z, freq2_filt_z, p0_z, p1_z, type_filt_z, order_filt_z, zero_padding_z, &
        nx, nz, dx, dz, xz, ext )
      use m_dag_data_kind
      use m_dag_transfer_funct
      use m_dag_space_2_filter_1D_win_zp
      use m_dag_space_2_filter_1d_win_zp
      implicit none
      
      ! The variables which are passed to the function.
      real(r8), intent(in)    :: freq1_filt_x, freq2_filt_x, freq1_filt_z, freq2_filt_z
      real(r8), intent(in)    :: p0_x, p1_x, p0_z, p1_z
      integer, intent(in)     :: type_filt_x, type_filt_z
      integer, intent(in)     :: order_filt_x, order_filt_z
      integer, intent(in)     :: zero_padding_x, zero_padding_z
      integer, intent(in)     :: nx, nz
      real(r8), intent(in)    :: dx, dz
      logical, intent(in)     :: xz, ext
      real(r8), intent(inout) :: image(:,:)
      
      ! The variables which are generated inside the function.
      integer :: i, j
!dir$ assume_aligned image(1,1):64
      
      !---------------------------------------------------------------------/
      ! Filtering.
      if ( xz .eqv. .false. ) then
        do i = 1,nx
          call space_2_filter_1D_win_zp( image(i,:), freq1_filt_z, freq2_filt_z, p0_z, p1_z, type_filt_z, order_filt_z, zero_padding_z, nz, dz, ext )
        end do
      end if  
      
      do j = 1,nz
        call space_2_filter_1D_win_zp( image(:,j), freq1_filt_x, freq2_filt_x, p0_x, p1_x, type_filt_x, order_filt_x, zero_padding_x, nx, dx, ext )
      end do
      
      if ( xz .eqv. .true. ) then
        do i = 1,nx
          call space_2_filter_1D_win_zp( image(i,:), freq1_filt_z, freq2_filt_z, p0_z, p1_z, type_filt_z, order_filt_z, zero_padding_z, nz, dz, ext )
        end do
      end if
    
    end subroutine space_2_filter_2D_zp
    
    
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*****     Real     **************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine space_2_filter_2D_zp_r( image, &
        freq1_filt_x, freq2_filt_x, p0_x, p1_x, type_filt_x, order_filt_x, zero_padding_x, &
        freq1_filt_z, freq2_filt_z, p0_z, p1_z, type_filt_z, order_filt_z, zero_padding_z, &
        nx, nz, dx, dz, xz, ext )
      use m_dag_data_kind
      use m_dag_transfer_funct
      use m_dag_space_2_filter_1d_win_zp
      implicit none
      
      ! The variables which are passed to the function.
      real(r4), intent(in)    :: freq1_filt_x, freq2_filt_x, freq1_filt_z, freq2_filt_z
      real(r4), intent(in)    :: p0_x, p1_x, p0_z, p1_z
      integer, intent(in)     :: type_filt_x, type_filt_z
      integer, intent(in)     :: order_filt_x, order_filt_z
      integer, intent(in)     :: zero_padding_x, zero_padding_z
      integer, intent(in)     :: nx, nz
      real(r4), intent(in)    :: dx, dz
      logical, intent(in)     :: xz, ext
      real(r4), intent(inout) :: image(:,:)
      
      ! The variables which are generated inside the function.
      integer :: i, j
!dir$ assume_aligned image(1,1):32
      
      !---------------------------------------------------------------------/
      ! Filtering.
      if ( xz .eqv. .false. ) then
        do i = 1,nx
          call space_2_filter_1D_win_zp_r( image(i,:), freq1_filt_z, freq2_filt_z, p0_z, p1_z, type_filt_z, order_filt_z, zero_padding_z, nz, dz, ext )
        end do
      end if  
      
      do j = 1,nz
        call space_2_filter_1D_win_zp_r( image(:,j), freq1_filt_x, freq2_filt_x, p0_x, p1_x, type_filt_x, order_filt_x, zero_padding_x, nx, dx, ext )
      end do
      
      if ( xz .eqv. .true. ) then
        do i = 1,nx
          call space_2_filter_1D_win_zp_r( image(i,:), freq1_filt_z, freq2_filt_z, p0_z, p1_z, type_filt_z, order_filt_z, zero_padding_z, nz, dz, ext )
        end do
      end if
    
    end subroutine space_2_filter_2D_zp_r

    
    
    
  !***********************************************************************/
  end module m_dag_space_2_filter_2D_zp