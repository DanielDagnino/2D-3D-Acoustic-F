  !*********************************************************************/
  !*****     Author: D. Dagnino     ************************************/
  !*********************************************************************/
  module m_dag_factor_number
    implicit none
    
    contains
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    pure integer function near_small_prime_factoriz( nt )
      implicit none
      ! The variables which are passed to the function.
      integer, intent(in) :: nt
      ! The variables which are generated inside the function.
      integer :: k
      integer, parameter :: n_best_number = 40
      integer :: best_number(n_best_number)
      
      !---------------------------------------------------------------------/
      k = 0
      
      k = k + 1
      best_number(k)  = 128   ! 2^7
      k = k + 1
      best_number(k)  = 256   ! 2^8
      k = k + 1
      best_number(k)  = 384   ! 2^7  * 3
      k = k + 1
      best_number(k)  = 512   ! 2^9
      k = k + 1
      best_number(k)  = 576  ! 2^6 * 3^2
      k = k + 1
      best_number(k)  = 640  ! 2^7 * 5
      k = k + 1
      best_number(k)  = 768   ! 2^8  * 3
      k = k + 1
      best_number(k)  = 896  ! 2^7 * 7
      k = k + 1
      best_number(k)  = 960  ! 2^6 * 3 * 5
      
      
      k = k + 1
      best_number(k)  = 1024  ! 2^10
      k = k + 1
      best_number(k)  = 1152  ! 2^7 * 3^2
      k = k + 1
      best_number(k)  = 1296  ! 2^4 * 3^4
      k = k + 1
      best_number(k)  = 1344  ! 2^6 * 3 * 7
      k = k + 1
      best_number(k)  = 1408  ! 2^7 * 11
      k = k + 1
      best_number(k)  = 1584  ! 2^4 * 3^2 * 11
      k = k + 1
      best_number(k)  = 1728  ! 2^6 * 3^3
      k = k + 1
      best_number(k)  = 1792  ! 2^8 * 7
      k = k + 1
      best_number(k)  = 1920  ! 2^7 * 3 * 5
      
      
      k = k + 1
      best_number(k)  = 2048  ! 2^11
      k = k + 1
      best_number(k)  = 2240  ! 2^6 * 5 * 7
      k = k + 1
      best_number(k)  = 2304  ! 2^8 * 3^2
      k = k + 1
      best_number(k)  = 2592  ! 2^5 * 3^4
      k = k + 1
      best_number(k)  = 2816  ! 2^8 * 11
      
      
      k = k + 1
      best_number(k)  = 3072  ! 2^10 * 3
      k = k + 1
      best_number(k)  = 4096  ! 2^12
      k = k + 1
      best_number(k) = 5120  ! 2^10 * 5
      k = k + 1
      best_number(k) = 6144  ! 2^11 * 3
      k = k + 1
      best_number(k) = 7168  ! 2^10 * 7
      k = k + 1
      best_number(k) = 8192  ! 2^13
      k = k + 1
      best_number(k) = 9216  ! 2^10 * 3^2
      k = k + 1
      best_number(k) = 10240 ! 2^11 * 5
      k = k + 1
      best_number(k) = 12288 ! 2^12 * 3
      k = k + 1
      best_number(k) = 14336 ! 2^11 * 7
      k = k + 1
      best_number(k) = 16384 ! 2^14
      k = k + 1
      best_number(k) = 18432 ! 2^11 * 3^2
      k = k + 1
      best_number(k) = 20480 ! 2^12 * 5
      k = k + 1
      best_number(k) = 24576 ! 2^13 * 3
      k = k + 1
      best_number(k) = 27648 ! 2^10 * 3^3
      k = k + 1
      best_number(k) = 30720 ! 2^11 * 3 * 5
      k = k + 1
      best_number(k) = 32768 ! 2^15
      
      ! 
      do k=1,n_best_number
        if ( best_number(k) >= nt ) exit
      end do
      
      near_small_prime_factoriz = best_number(k)
    
      !---------------------------------------------------------------------/
    end function near_small_prime_factoriz

  !***********************************************************************/
  end module m_dag_factor_number