
  !*********************************************************************/
  !*****     Author: D. Dagnino     ************************************/
  !*********************************************************************/
  module m_dag_filter_1D_zp
    
    interface filter_1D_zp
      module procedure filter_1D_zp, filter_1D_zp_r
    end interface
    
    contains
    !*********************************************************************/
    !*********************************************************************/
    !*****     Double     ************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine filter_1D_zp( trace, freq1_filt, freq2_filt, type_filt, order_filt, zero_padding, nt, dt )
      use m_dag_data_kind
      use m_dag_transfer_funct
      use m_dag_factor_number
      implicit none
     
      ! The variables which are passed to the function.
      real(r8), intent(in)    :: freq1_filt, freq2_filt
      integer, intent(in)     :: type_filt
      integer, intent(in)     :: order_filt
      integer, intent(in)     :: zero_padding
      integer, intent(in)     :: nt
      real(r8), intent(in)    :: dt
      real(r8), intent(inout) :: trace(:)
     
      ! The variables which are generated inside the function.
      integer     :: k
      integer     :: nt_pad
      complex(r8) :: cH, caux
      real(r8)    :: w, dw
      integer     :: lensav, lenwrk, iaux
      real(r8), allocatable :: wsave(:), work(:)
      complex(r8), allocatable :: ftmp(:)
      integer :: stal
      integer :: ind_middle
!dir$ assume_aligned trace(1):64
     
      !*********************************************************************/
!       nt_pad = (1+zero_padding)*nt
      if ( zero_padding==0 ) then
        nt_pad = nt
      else
        nt_pad = near_small_prime_factoriz( nt )*(1+zero_padding)
      end if
      
!       if ( mod(nt_pad,2) /= 0 ) stop '***** ERROR - filter_1D: mod(nt_pad,2) /= 0 *****'
     
      allocate( ftmp(nt_pad), stat=stal ); if ( stal/=0 ) stop ' ***** AE-1: filter_1D ***** '
     
      lensav = 2*nt_pad + ceiling(dlog(dble(nt_pad))/dlog(2._r8)) + 4
      lenwrk = 2*nt_pad
     
      allocate( wsave(lensav), work(lenwrk), stat=stal ); if ( stal/=0 ) stop ' ***** AE-2: filter_1D ***** '
      
      call dcfft1i(nt_pad,wsave,lensav,iaux)
      
      !*********************************************************************/
      ! 
      ftmp(1:nt) = dcmplx( trace )
      if ( nt_pad /= nt ) ftmp(nt+1:nt_pad) = 0._r8
      
      !*********************************************************************/
      !
      if ( mod(nt_pad,2) /= 0 ) then
        ind_middle = (nt_pad-1)/2
      else
        ind_middle = nt_pad/2
      end if
      
      dw   = 1._r8/( dt*dble(nt_pad) )
      
      !---------------------------------------------------------------------/
      call dcfft1f(nt_pad,1,ftmp,nt_pad,wsave,lensav,work,lenwrk,iaux)
      
      if ( type_filt == 1 ) then
        cH = butterworth( 0._r8, order_filt, freq1_filt, freq2_filt, type_filt )
        ftmp(1) = ftmp(1)*cH
      else
        ftmp(1) = 0._r8
      end if
      
      do k = 2,ind_middle
        w = dw*dble(k-1)
        cH = butterworth( w, order_filt, freq1_filt, freq2_filt, type_filt )
        ftmp(k)          = ftmp(k)*cH
        ftmp(nt_pad-k+2) = dconjg(ftmp(k))
      end do
      
      if ( mod(nt_pad,2) == 0 ) then
        k = ind_middle+1
        w = dw*dble(k-1)
        cH = butterworth( w, order_filt, freq1_filt, freq2_filt, type_filt )
        ftmp(k) = ftmp(k)*cH
      end if
      
      call dcfft1b(nt_pad,1,ftmp,nt_pad,wsave,lensav,work,lenwrk,iaux)
      
      !---------------------------------------------------------------------/
      ! 
      do k=1,ind_middle
        caux             = ftmp(nt_pad-k+1)
        ftmp(nt_pad-k+1) = ftmp(k)
        ftmp(k)          = caux
      end do
      
      !---------------------------------------------------------------------/
      call dcfft1f(nt_pad,1,ftmp,nt_pad,wsave,lensav,work,lenwrk,iaux)
      
      if ( type_filt == 1 ) then
        cH = butterworth( 0._r8, order_filt, freq1_filt, freq2_filt, type_filt )
        ftmp(1) = ftmp(1)*cH
      else
        ftmp(1) = 0._r8
      end if
      
      do k = 2,ind_middle
        w = dw*dble(k-1)
        cH = butterworth( w, order_filt, freq1_filt, freq2_filt, type_filt )
        ftmp(k)          = ftmp(k)*cH
        ftmp(nt_pad-k+2) = dconjg(ftmp(k))
      end do
      
      if ( mod(nt_pad,2) == 0 ) then
        k = ind_middle+1
        w = dw*dble(k-1)
        cH = butterworth( w, order_filt, freq1_filt, freq2_filt, type_filt )
        ftmp(k) = ftmp(k)*cH
      end if
      
      call dcfft1b(nt_pad,1,ftmp,nt_pad,wsave,lensav,work,lenwrk,iaux)
     
      !---------------------------------------------------------------------/
      ! 
      do k=1,nt
        trace(k) = dble( ftmp(nt_pad-k+1) )
      end do
      
      !---------------------------------------------------------------------/
      deallocate( ftmp, stat=stal ); if ( stal/=0 ) stop ' ***** AE-1d: filter_1D ***** '
      deallocate( wsave, work, stat=stal ); if ( stal/=0 ) stop ' ***** AE-2d: filter_1D ***** '
   
    end subroutine filter_1D_zp
   
   
   
   
   
    !*********************************************************************/
    !*********************************************************************/
    !*****     Real     **************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine filter_1D_zp_r( trace, freq1_filt, freq2_filt, type_filt, order_filt, zero_padding, nt, dt )
      use m_dag_data_kind
      use m_dag_transfer_funct
      use m_dag_factor_number
      implicit none
     
      ! The variables which are passed to the function.
      real(r4), intent(in)    :: freq1_filt, freq2_filt
      integer, intent(in)     :: type_filt
      integer, intent(in)     :: order_filt
      integer, intent(in)     :: zero_padding
      integer, intent(in)     :: nt
      real(r4), intent(in)    :: dt
      real(r4), intent(inout) :: trace(:)

      ! The variables which are generated inside the function.
      integer     :: k
      integer     :: nt_pad
      complex(r4) :: cH, caux
      real(r4)    :: w, dw
      integer     :: lensav, lenwrk, iaux
      real(r4), allocatable :: wsave(:), work(:)
      complex(r4), allocatable :: ftmp(:)
      integer :: stal
      integer :: ind_middle
!dir$ assume_aligned trace(1):32
     
      !*********************************************************************/
!       nt_pad = (1+zero_padding)*nt
      if ( zero_padding==0 ) then
        nt_pad = nt
      else
        nt_pad = near_small_prime_factoriz( nt )*(1+zero_padding)
      end if
      
!       if ( mod(nt_pad,2) /= 0 ) stop '***** ERROR - filter_1D: mod(nt_pad,2) /= 0 *****'
     
      allocate( ftmp(nt_pad), stat=stal ); if ( stal/=0 ) stop ' ***** AE-1: filter_1D ***** '
     
      lensav = 2*nt_pad + ceiling(log(real(nt_pad))/log(2._r4)) + 4
      lenwrk = 2*nt_pad
     
      allocate( wsave(lensav), work(lenwrk), stat=stal ); if ( stal/=0 ) stop ' ***** AE-2: filter_1D ***** '
      
      call cfft1i(nt_pad,wsave,lensav,iaux)
      
      !*********************************************************************/
      ! 
      ftmp(1:nt) = cmplx( trace )
      if ( nt_pad /= nt ) ftmp(nt+1:nt_pad) = 0._r4
      
      !*********************************************************************/
      !
      if ( mod(nt_pad,2) /= 0 ) then
        ind_middle = (nt_pad-1)/2
      else
        ind_middle = nt_pad/2
      end if
      
      dw   = 1._r4/( dt*real(nt_pad) )
      
      !---------------------------------------------------------------------/
      call cfft1f(nt_pad,1,ftmp,nt_pad,wsave,lensav,work,lenwrk,iaux)
      
      if ( type_filt == 1 ) then
        cH = butterworth( 0._r4, order_filt, freq1_filt, freq2_filt, type_filt )
        ftmp(1) = ftmp(1)*cH
      else
        ftmp(1) = 0._r4
      end if
      
      do k = 2,ind_middle
        w = dw*real(k-1)
        cH = butterworth( w, order_filt, freq1_filt, freq2_filt, type_filt )
        ftmp(k)          = ftmp(k)*cH
        ftmp(nt_pad-k+2) = conjg(ftmp(k))
      end do
      
      if ( mod(nt_pad,2) == 0 ) then
        k = ind_middle+1
        w = dw*real(k-1)
        cH = butterworth( w, order_filt, freq1_filt, freq2_filt, type_filt )
        ftmp(k) = ftmp(k)*cH
      end if
      
      call cfft1b(nt_pad,1,ftmp,nt_pad,wsave,lensav,work,lenwrk,iaux)
      
      !---------------------------------------------------------------------/
      ! 
      do k=1,ind_middle
        caux             = ftmp(nt_pad-k+1)
        ftmp(nt_pad-k+1) = ftmp(k)
        ftmp(k)          = caux
      end do
      
      !---------------------------------------------------------------------/
      call cfft1f(nt_pad,1,ftmp,nt_pad,wsave,lensav,work,lenwrk,iaux)
      
      if ( type_filt == 1 ) then
        cH = butterworth( 0._r4, order_filt, freq1_filt, freq2_filt, type_filt )
        ftmp(1) = ftmp(1)*cH
      else
        ftmp(1) = 0._r4
      end if
      
      do k = 2,ind_middle
        w = dw*real(k-1)
        cH = butterworth( w, order_filt, freq1_filt, freq2_filt, type_filt )
        ftmp(k)          = ftmp(k)*cH
        ftmp(nt_pad-k+2) = conjg(ftmp(k))
      end do
      
      if ( mod(nt_pad,2) == 0 ) then
        k = ind_middle+1
        w = dw*real(k-1)
        cH = butterworth( w, order_filt, freq1_filt, freq2_filt, type_filt )
        ftmp(k) = ftmp(k)*cH
      end if
      
      call cfft1b(nt_pad,1,ftmp,nt_pad,wsave,lensav,work,lenwrk,iaux)
     
      !---------------------------------------------------------------------/
      ! 
      do k=1,nt
        trace(k) = real( ftmp(nt_pad-k+1) )
      end do
      
      !---------------------------------------------------------------------/
      deallocate( ftmp, stat=stal ); if ( stal/=0 ) stop ' ***** AE-1d: filter_1D ***** '
      deallocate( wsave, work, stat=stal ); if ( stal/=0 ) stop ' ***** AE-2d: filter_1D ***** '
   
    end subroutine filter_1D_zp_r
   
  !***********************************************************************/
  end module m_dag_filter_1D_zp




