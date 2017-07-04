  !**********************************************************!
  module m_dag_split_work
    implicit none
  
    !**********************************************************!
    !**********************************************************!
    !**********************************************************!
    contains
    
    !**********************************************************!
    pure subroutine split_work( n_shot, n_procs, n_sg )
      ! The variables which are passed to the function.
      integer, intent(in)  :: n_shot, n_procs
      integer, intent(out) :: n_sg(:)
      ! The variables which are generated inside the function.
      integer :: k
            
      n_sg = floor(dble(n_shot)/dble(n_procs))
      do k=1,mod(n_shot,n_procs)
        n_sg(k) = n_sg(k) + 1
      end do
      
    end subroutine split_work
    
    !**********************************************************!
  end module m_dag_split_work








