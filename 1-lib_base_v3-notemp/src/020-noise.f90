  !*********************************************************************/
  !*****     Author: D. Dagnino     ************************************/
  !*********************************************************************/
  !*****     FFTPACK labrary is used in this code.     *****************/
  !*****     FFTPACK: P. Swarztrauber and D. Valent    *****************/
  !*********************************************************************/
  module m_dag_noise
    implicit none
    
    interface noise_gen
      module procedure noise_gen, noise_gen_r
    end interface
    
    contains
    !*********************************************************************/
    !*********************************************************************/
    !*****     Double     ************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine noise_gen( spectrum, nt, noise ) 
      use m_dag_data_kind
      use m_dag_random_dist
      implicit none
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r8
      integer, intent(in)   :: nt
      real(prcsn), intent(in)  :: spectrum(:)
      real(prcsn), intent(out) :: noise(:)
      ! The variables which are generated inside the function.
      integer :: k
      integer :: ind_middle
      integer :: lensav, lenwrk, iaux
      real(prcsn), allocatable :: wsave(:), work(:)
      integer :: stal
      complex(prcsn) :: cnoise(nt)
!dir$ assume_aligned spectrum(1):64,noise(1):64

      !---------------------------------------------------------------------/
      ! 
      if ( mod(nt,2) /= 0 ) then
        ind_middle = (nt-1)/2
      else
        ind_middle = nt/2
      end if
      
      !---------------------------------------------------------------------/
      ! White noise.
      call normal_distrib( 1._prcsn, nt, noise )
      
      cnoise = cmplx( noise )
      
      !---------------------------------------------------------------------/
      ! Spectrum weighting.
      lensav = 2*nt + ceiling(log(real(nt,prcsn))/log(2._prcsn)) + 4
      lenwrk = 2*nt
      
      allocate( wsave(lensav), work(lenwrk), stat=stal ); if ( stal/=0 ) stop ' ***** AE-1: noise_gen ***** '
      
      ! 
      call dcfft1i(nt,wsave,lensav,iaux)
      call dcfft1f(nt,1,cnoise,nt,wsave,lensav,work,lenwrk,iaux)
      
      cnoise(1) = cnoise(1)*spectrum(1)
      
      do k=2,ind_middle
        cnoise(k)      = cnoise(k)*spectrum(k)
        cnoise(nt-k+2) = conjg( cnoise(nt-k+2)*spectrum(nt-k+2) )
      end do
      
      if ( mod(nt,2) == 0 ) then
        k = ind_middle+1
        cnoise(k) = cnoise(k)*spectrum(k)
      end if
      
      ! 
      call dcfft1b(nt,1,cnoise,nt,wsave,lensav,work,lenwrk,iaux)
      
      ! 
      noise = real( cnoise )
      
      !---------------------------------------------------------------------/
      deallocate( wsave, work, stat=stal ); if ( stal/=0 ) stop ' ***** AE-1d: noise_gen ***** '
      
    end subroutine noise_gen





    !*********************************************************************/
    !*********************************************************************/
    !*****     Real     **************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine noise_gen_r( spectrum, nt, noise ) 
      use m_dag_data_kind
      use m_dag_random_dist
      implicit none
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r4
      integer, intent(in)   :: nt
      real(prcsn), intent(in)  :: spectrum(:)
      real(prcsn), intent(out) :: noise(:)
      ! The variables which are generated inside the function.
      integer :: k
      integer :: ind_middle
      integer :: lensav, lenwrk, iaux
      real(prcsn), allocatable :: wsave(:), work(:)
      integer :: stal
      complex(prcsn) :: cnoise(nt)
!dir$ assume_aligned spectrum(1):32,noise(1):32

      !---------------------------------------------------------------------/
      ! 
      if ( mod(nt,2) /= 0 ) then
        ind_middle = (nt-1)/2
      else
        ind_middle = nt/2
      end if
      
      !---------------------------------------------------------------------/
      ! White noise.
      call normal_distrib( 1._prcsn, nt, noise )
      
      cnoise = cmplx( noise )
      
      !---------------------------------------------------------------------/
      ! Spectrum weighting.
      lensav = 2*nt + ceiling(log(real(nt,prcsn))/log(2._prcsn)) + 4
      lenwrk = 2*nt
      
      allocate( wsave(lensav), work(lenwrk), stat=stal ); if ( stal/=0 ) stop ' ***** AE-1: noise_gen ***** '
      
      ! 
      call cfft1i(nt,wsave,lensav,iaux)
      call cfft1f(nt,1,cnoise,nt,wsave,lensav,work,lenwrk,iaux)
      
      cnoise(1) = cnoise(1)*spectrum(1)
      
      do k=2,ind_middle
        cnoise(k)      = cnoise(k)*spectrum(k)
        cnoise(nt-k+2) = conjg( cnoise(nt-k+2)*spectrum(nt-k+2) )
      end do
      
      if ( mod(nt,2) == 0 ) then
        k = ind_middle+1
        cnoise(k) = cnoise(k)*spectrum(k)
      end if
      
      ! 
      call cfft1b(nt,1,cnoise,nt,wsave,lensav,work,lenwrk,iaux)
      
      ! 
      noise = real( cnoise )
      
      !---------------------------------------------------------------------/
      deallocate( wsave, work, stat=stal ); if ( stal/=0 ) stop ' ***** AE-1d: noise_gen ***** '
      
    end subroutine noise_gen_r
    
  !***********************************************************************/
  end module m_dag_noise