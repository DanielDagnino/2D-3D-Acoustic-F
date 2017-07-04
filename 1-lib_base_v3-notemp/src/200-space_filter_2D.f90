  !*********************************************************************/
  !*****     Author: D. Dagnino     ************************************/
  !*********************************************************************/
  module m_dag_space_filter_2D
    
    interface space_filter_2D
      module procedure space_filter_2D, space_filter_2D_r
    end interface
    
    contains
    !*********************************************************************/
    !*********************************************************************/
    !*****     Double     ************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine space_filter_2D( image, &
        freq1_filt_x, freq2_filt_x, type_filt_x, order_filt_x, zero_padding_x, &
        freq1_filt_z, freq2_filt_z, type_filt_z, order_filt_z, zero_padding_z, &
        nx, nz, dx, dz, xz, ext )
      use m_dag_data_kind
      use m_dag_transfer_funct
      use m_dag_space_filter_1D_win
      implicit none
      
      ! The variables which are passed to the function.
      real(r8), intent(in)    :: freq1_filt_x, freq2_filt_x, freq1_filt_z, freq2_filt_z
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
          call space_filter_1D_win( image(i,:), freq1_filt_z, freq2_filt_z, type_filt_z, order_filt_z, zero_padding_z, .true., nz, dz, ext )
        end do
      end if  
      
      do j = 1,nz
        call space_filter_1D_win( image(:,j), freq1_filt_x, freq2_filt_x, type_filt_x, order_filt_x, zero_padding_x, .true., nx, dx, ext )
      end do
      
      if ( xz .eqv. .true. ) then
        do i = 1,nx
          call space_filter_1D_win( image(i,:), freq1_filt_z, freq2_filt_z, type_filt_z, order_filt_z, zero_padding_z, .true., nz, dz, ext )
        end do
      end if
    
    end subroutine space_filter_2D
    
    
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*****     Real     **************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine space_filter_2D_r( image, &
        freq1_filt_x, freq2_filt_x, type_filt_x, order_filt_x, zero_padding_x, &
        freq1_filt_z, freq2_filt_z, type_filt_z, order_filt_z, zero_padding_z, &
        nx, nz, dx, dz, xz, ext )
      use m_dag_data_kind
      use m_dag_transfer_funct
      use m_dag_space_filter_1D_win
      implicit none
      
      ! The variables which are passed to the function.
      real(r4), intent(in)    :: freq1_filt_x, freq2_filt_x, freq1_filt_z, freq2_filt_z
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
          call space_filter_1D_win_r( image(i,:), freq1_filt_z, freq2_filt_z, type_filt_z, order_filt_z, zero_padding_z, .true., nz, dz, ext )
        end do
      end if  
      
      do j = 1,nz
        call space_filter_1D_win_r( image(:,j), freq1_filt_x, freq2_filt_x, type_filt_x, order_filt_x, zero_padding_x, .true., nx, dx, ext )
      end do
      
      if ( xz .eqv. .true. ) then
        do i = 1,nx
          call space_filter_1D_win_r( image(i,:), freq1_filt_z, freq2_filt_z, type_filt_z, order_filt_z, zero_padding_z, .true., nz, dz, ext )
        end do
      end if
    
    end subroutine space_filter_2D_r

    
    
    
  !***********************************************************************/
  end module m_dag_space_filter_2D