  !*********************************************************************/
  !*****     Author: Internet and addapted by D. Dagnino     ***********/
  !*********************************************************************/
  module m_dag_init_random_seed
    implicit none
  
    contains
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine init_random_seed( clock_shift )
      implicit none

      ! The variables which are passed to the function.
      integer, intent(in) :: clock_shift ! Useful to avoid same seed for different threats when parallel computation is used.

      ! The variables which are generated inside the function.
      integer :: i, n, clock
      integer, ALLOCATABLE :: seed(:)
      integer :: stal

      !---------------------------------------------------------------------/
      call RANDOM_SEED(size = n)
      ALLOCATE( seed(n), stat=stal ); if ( stal/=0 ) stop ' ***** AE-1: init_random_seed ***** '
      
      call SYSTEM_CLOCK(COUNT=clock)
      
      seed = clock*(10*clock_shift+1) + 37 * (/ (i - 1, i = 1, n) /)
      call RANDOM_SEED(PUT = seed)
      
      DEALLOCATE( seed, stat=stal ); if ( stal/=0 ) stop ' ***** AE-1d: init_random_seed ***** '
      
    end subroutine init_random_seed

  !***********************************************************************/
  end module m_dag_init_random_seed