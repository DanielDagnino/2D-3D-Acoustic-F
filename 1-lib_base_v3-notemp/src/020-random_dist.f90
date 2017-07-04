  !*********************************************************************/
  !*****     Author: D. Dagnino     ************************************/
  !*********************************************************************/
  module m_dag_random_dist
    implicit none
    
    interface normal_distrib
      module procedure normal_distrib, normal_distrib_r
    end interface
    
    interface normal_var
      module procedure normal_var, normal_var_r
    end interface
    
    contains
    !*********************************************************************/
    !*********************************************************************/
    !*****     Double     ************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine normal_distrib( amp_norm_var, nt, norm_var )
      use m_dag_data_kind
      implicit none
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r8
      integer, intent(in)   :: nt
      real(prcsn), intent(in)  :: amp_norm_var   ! Amplitud of norm_var.
      real(prcsn), intent(out) :: norm_var(:)
      ! The variables which are generated inside the function.
      integer  :: k
      real(prcsn) :: rand_var1, rand_var2
!dir$ assume_aligned norm_var(1):64
      
      !---------------------------------------------------------------------/
      ! 
      do k=1,nt-1,2
        call normal_var( 0._prcsn, amp_norm_var, rand_var1, rand_var2 )
        norm_var(k)   = rand_var1
        norm_var(k+1) = rand_var2
      end do
      
      ! 
      if ( mod(nt,2) /= 0 ) then
        call normal_var( 0._prcsn, amp_norm_var, rand_var1, rand_var2 )
        norm_var(nt) = rand_var1
      end if
    
    end subroutine normal_distrib
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine normal_var( mean, sigma, rand_var1, rand_var2 )
      use m_dag_data_kind
      implicit none
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r8
      real(prcsn), intent(out) :: rand_var1, rand_var2
      real(prcsn), intent(in)  :: mean, sigma ! standard_deviation or sqrt(variance) or sqrt(<rand_var^2>-<rand_var>^2).
      ! The variables which are generated inside the function.
      real(prcsn) :: s
      
      ! Box-Muller Method.
      ! The Box-Muller method uses the technique of inverse transformation to turn two uniformly distributed randoms, U1 and U2, into two unit normal randoms, X and Y.
      s = 2._prcsn
      do while ( s >= 1._prcsn )
        call random_number(rand_var1)  ! U1=[0,1]
        call random_number(rand_var2)  ! U2=[0,1]
        rand_var1 = 2._prcsn*rand_var1-1._prcsn    ! V1=[-1,1]
        rand_var2 = 2._prcsn*rand_var2-1._prcsn    ! V2=[-1,1]
        s = rand_var1*rand_var1 + rand_var2*rand_var2
      end do
      
      rand_var1 = sqrt( -2._prcsn*log(s)/s ) * rand_var1
      rand_var2 = sqrt( -2._prcsn*log(s)/s ) * rand_var2
      
      ! The above is called the Polar Method and is fully described in the Ross book [Ros88].
      ! X and Y will be unit normal random variables (mean = 0 and sigma = 1), and can be easilly modified for different mean and variance.
      rand_var1 = mean + sigma * rand_var1
      rand_var2 = mean + sigma * rand_var2
    
    end subroutine normal_var
    




    !*********************************************************************/
    !*********************************************************************/
    !*****     Real     **************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine normal_distrib_r( amp_norm_var, nt, norm_var )
      use m_dag_data_kind
      implicit none
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r4
      integer, intent(in)   :: nt
      real(prcsn), intent(in)  :: amp_norm_var   ! Amplitud of norm_var.
      real(prcsn), intent(out) :: norm_var(:)
      ! The variables which are generated inside the function.
      integer  :: k
      real(prcsn) :: rand_var1, rand_var2
!dir$ assume_aligned norm_var(1):32
      
      !---------------------------------------------------------------------/
      ! 
      do k=1,nt-1,2
        call normal_var( 0._prcsn, amp_norm_var, rand_var1, rand_var2 )
        norm_var(k)   = rand_var1
        norm_var(k+1) = rand_var2
      end do
      
      ! 
      if ( mod(nt,2) /= 0 ) then
        call normal_var( 0._prcsn, amp_norm_var, rand_var1, rand_var2 )
        norm_var(nt) = rand_var1
      end if
    
    end subroutine normal_distrib_r
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine normal_var_r( mean, sigma, rand_var1, rand_var2 )
      use m_dag_data_kind
      implicit none
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r4
      real(prcsn), intent(out) :: rand_var1, rand_var2
      real(prcsn), intent(in)  :: mean, sigma ! standard_deviation or sqrt(variance) or sqrt(<rand_var^2>-<rand_var>^2).
      ! The variables which are generated inside the function.
      real(prcsn) :: s
      
      ! Box-Muller Method.
      ! The Box-Muller method uses the technique of inverse transformation to turn two uniformly distributed randoms, U1 and U2, into two unit normal randoms, X and Y.
      s = 2._prcsn
      do while ( s >= 1._prcsn )
        call random_number(rand_var1)  ! U1=[0,1]
        call random_number(rand_var2)  ! U2=[0,1]
        rand_var1 = 2._prcsn*rand_var1-1._prcsn    ! V1=[-1,1]
        rand_var2 = 2._prcsn*rand_var2-1._prcsn    ! V2=[-1,1]
        s = rand_var1*rand_var1 + rand_var2*rand_var2
      end do
      
      rand_var1 = sqrt( -2._prcsn*log(s)/s ) * rand_var1
      rand_var2 = sqrt( -2._prcsn*log(s)/s ) * rand_var2
      
      ! The above is called the Polar Method and is fully described in the Ross book [Ros88].
      ! X and Y will be unit normal random variables (mean = 0 and sigma = 1), and can be easilly modified for different mean and variance.
      rand_var1 = mean + sigma * rand_var1
      rand_var2 = mean + sigma * rand_var2
    
    end subroutine normal_var_r
    
  !***********************************************************************/
  end module m_dag_random_dist
