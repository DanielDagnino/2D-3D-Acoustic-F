  !*********************************************************************/
  !*****     Author: D. Dagnino     ************************************/
  !*********************************************************************/
  module m_dag_space_filter_1D_win
    
    interface space_filter_1D_win
      module procedure space_filter_1D_win, space_filter_1D_win_r
    end interface
    
    contains
    !*********************************************************************/
    !*********************************************************************/
    !*****     Double     ************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine space_filter_1D_win( trace, freq1_filt, freq2_filt, type_filt, order_filt, zero_padding, zero_delay, nt, dt, ext )
      use m_dag_data_kind
      use m_dag_transfer_funct
      use m_dag_mat
      implicit none
     
      ! The variables which are passed to the function.
      real(r8), intent(in)    :: freq1_filt, freq2_filt
      integer, intent(in)     :: type_filt
      integer, intent(in)     :: order_filt
      integer, intent(in)     :: zero_padding
      logical, intent(in)     :: zero_delay, ext
      integer, intent(in)     :: nt
      real(r8), intent(in)    :: dt
      real(r8), intent(inout) :: trace(:)
     
      ! The variables which are generated inside the function.
      integer     :: k
      integer     :: nt_pad
      real(r8)    :: H
      complex(r8) :: cH
      real(r8)    :: w, dw, alpha, inct
      integer     :: iinct
      integer     :: lensav, lenwrk, iaux
      real(r8), allocatable :: wsave(:), work(:)
      complex(r8), allocatable :: ftmp(:)
      real(r8), allocatable :: wind(:)
      integer :: stal
      integer :: ind_middle
      real(r8) :: tr1
!dir$ assume_aligned trace(1):64
      
      !---------------------------------------------------------------------/
      nt_pad = (1+zero_padding)*nt
      if ( mod(nt_pad,2) /= 0 ) stop '***** ERROR - space_filter_1D: mod(nt_pad,2) /= 0 *****'
     
      allocate( ftmp(nt_pad), stat=stal ); if ( stal/=0 ) stop ' ***** AE-1: space_filter_1D ***** '
      allocate( wind(nt_pad), stat=stal ); if ( stal/=0 ) stop ' ***** AE-3: space_filter_1D ***** '
     
      lensav = 2*nt_pad + ceiling(dlog(dble(nt_pad))/dlog(2._r8)) + 4
      lenwrk = 2*nt_pad
     
      allocate( wsave(lensav), work(lenwrk), stat=stal ); if ( stal/=0 ) stop ' ***** AE-2: space_filter_1D ***** '
     
      call dcfft1i(nt_pad,wsave,lensav,iaux)
     
      !---------------------------------------------------------------------/
      ! 
      inct = 0.5_r8*(nt_pad-nt)
      alpha = 2._r8*inct/dble(nt_pad)
      iinct = floor(inct)
      
      !---------------------------------------------------------------------/
      ! 
      tr1 = 0._r4
      if ( ext ) tr1 = trace(1)
      trace = trace-tr1
      if ( zero_padding /= 0 ) then
        ftmp(1:nt) = dcmplx( trace )
        ftmp(nt+1:nt+iinct) = dcmplx( trace(nt) )
        ftmp(1+nt+iinct:nt_pad) = dcmplx( 0._r8 )
      end if
      
      !---------------------------------------------------------------------/
      ! 
      if ( alpha > 1._r8 ) stop '***** ERROR - space_filter_1D: alpha > 1._r8 *****'
      if ( zero_padding == 0 ) stop '***** ERROR - space_filter_1D: zero_padding == 0 *****'
      
      ! 
      do k=1,nt_pad
      
        if ( k <= nt ) then
        
          wind(k) = 1._r8
!           write(*,*) wind(k)
          
        else
          
          wind(k) = 0.5_r8*( 1._r8 + dcos( pi*( 2._r8*dble(k-nt)/dble(nt_pad-nt+1)) ) )
!           write(*,*) wind(k)
          
        end if
        
      end do
      
!       write(*,*)
!       write(*,*)
! 
!       do k=1,nt_pad
!         if (k<=nt) then
!           write(*,*)  real( wind(k)*ftmp(k) ), real(ftmp(k)), trace(k)
!         else
!           write(*,*)  real( wind(k)*ftmp(k) ), real(ftmp(k))
!         end if
!       end do
!       
!       stop
      
      !---------------------------------------------------------------------/
      !
      ftmp = cmplx( wind*ftmp )
      
      !---------------------------------------------------------------------/
      !
      call dcfft1f(nt_pad,1,ftmp,nt_pad,wsave,lensav,work,lenwrk,iaux)
      
      !
      if ( mod(nt_pad,2) /= 0 ) then
        ind_middle = (nt_pad-1)/2
      else
        ind_middle = nt_pad/2
      end if
      
      dw   = 1._r8/( dt*dble(nt_pad) )
      
      if ( zero_delay ) then
      
        if ( type_filt == 1 ) then
          H = butterworth0( 0._r8, order_filt, freq1_filt, freq2_filt, type_filt )
          ftmp(1) = ftmp(1)*H
        else
          ftmp(1) = 0._r8
        end if
       
        do k = 2,ind_middle
          w = dw*dble(k-1)
          H = butterworth0( w, order_filt, freq1_filt, freq2_filt, type_filt )
          ftmp(k)          = ftmp(k)*H
          ftmp(nt_pad-k+2) = dconjg(ftmp(k))
        end do
        
        if ( mod(nt_pad,2) == 0 ) then
          k = ind_middle+1
          w = dw*dble(k-1)
          H = butterworth0( w, order_filt, freq1_filt, freq2_filt, type_filt )
          ftmp(k) = ftmp(k)*H
        end if
        
      else
      
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
        
      end if
     
      !
      call dcfft1b(nt_pad,1,ftmp,nt_pad,wsave,lensav,work,lenwrk,iaux)
     
      !
      trace = tr1 + dble( ftmp(1:nt) )
       
      !---------------------------------------------------------------------/
      deallocate( ftmp, stat=stal ); if ( stal/=0 ) stop ' ***** AE-1d: space_filter_1D ***** '
      deallocate( wind, stat=stal ); if ( stal/=0 ) stop ' ***** AE-3d: space_filter_1D ***** '
      deallocate( wsave, work, stat=stal ); if ( stal/=0 ) stop ' ***** AE-2d: space_filter_1D ***** '
   
    end subroutine space_filter_1D_win
   
   
   
   
   
    !*********************************************************************/
    !*********************************************************************/
    !*****     Real     **************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine space_filter_1D_win_r( trace, freq1_filt, freq2_filt, type_filt, order_filt, zero_padding, zero_delay, nt, dt, ext )
      use m_dag_data_kind
      use m_dag_transfer_funct
      use m_dag_mat
      implicit none
      
      ! The variables which are passed to the function.
      real(r4), intent(in)    :: freq1_filt, freq2_filt
      integer, intent(in)     :: type_filt
      integer, intent(in)     :: order_filt
      integer, intent(in)     :: zero_padding
      logical, intent(in)     :: zero_delay, ext
      integer, intent(in)     :: nt
      real(r4), intent(in)    :: dt
      real(r4), intent(inout) :: trace(:)
     
      ! The variables which are generated inside the function.
      integer     :: k
      integer     :: nt_pad
      real(r4)    :: H
      complex(r4) :: cH
      real(r4)    :: w, dw, alpha, inct
      integer     :: iinct
      integer     :: lensav, lenwrk, iaux
      real(r4), allocatable :: wsave(:), work(:)
      complex(r4), allocatable :: ftmp(:)
      real(r4), allocatable :: wind(:)
      integer :: stal
      integer :: ind_middle
      real(r4) :: tr1
!dir$ assume_aligned trace(1):32
     
      !---------------------------------------------------------------------/
      nt_pad = (1+zero_padding)*nt
      if ( mod(nt_pad,2) /= 0 ) stop '***** ERROR - space_filter_1D: mod(nt_pad,2) /= 0 *****'
     
      allocate( ftmp(nt_pad), stat=stal ); if ( stal/=0 ) stop ' ***** AE-1: space_filter_1D ***** '
      allocate( wind(nt_pad), stat=stal ); if ( stal/=0 ) stop ' ***** AE-3: space_filter_1D ***** '
     
      lensav = 2*nt_pad + ceiling(log(real(nt_pad))/log(2._r4)) + 4
      lenwrk = 2*nt_pad
     
      allocate( wsave(lensav), work(lenwrk), stat=stal ); if ( stal/=0 ) stop ' ***** AE-2: space_filter_1D ***** '
     
      call cfft1i(nt_pad,wsave,lensav,iaux)
     
      !---------------------------------------------------------------------/
      ! 
      inct = 0.5_r4*(nt_pad-nt)
      alpha = 2._r4*inct/real(nt_pad)
      iinct = floor(inct)
      
      !---------------------------------------------------------------------/
      ! 
      tr1 = 0._r4
      if ( ext ) tr1 = trace(1)
      trace = trace-tr1
      if ( zero_padding /= 0 ) then
        ftmp(1:nt) = cmplx( trace )
        ftmp(nt+1:nt+iinct) = cmplx( trace(nt) )
        ftmp(1+nt+iinct:nt_pad) = cmplx( 0._r4 )
      end if
      
      !---------------------------------------------------------------------/
      ! 
      if ( alpha > 1._r4 ) stop '***** ERROR - space_filter_1D: alpha > 1._r4 *****'
      if ( zero_padding == 0 ) stop '***** ERROR - space_filter_1D: zero_padding == 0 *****'
      
      ! 
      do k=1,nt_pad
      
        if ( k <= nt ) then
        
          wind(k) = 1._r4
          
        else
          
          wind(k) = 0.5_r4*( 1._r4 + cos( pi*( 2._r4*real(k-nt)/real(nt_pad-nt+1) ) ) )
          
        end if
        
      end do
      
      !---------------------------------------------------------------------/
      !
      ftmp = cmplx( wind*ftmp )
      
      !---------------------------------------------------------------------/
      !
      call cfft1f(nt_pad,1,ftmp,nt_pad,wsave,lensav,work,lenwrk,iaux)
      
      !
      if ( mod(nt_pad,2) /= 0 ) then
        ind_middle = (nt_pad-1)/2
      else
        ind_middle = nt_pad/2
      end if
      
      dw   = 1._r4/( dt*real(nt_pad) )
      
      if ( zero_delay ) then
      
        if ( type_filt == 1 ) then
          H = butterworth0( 0._r4, order_filt, freq1_filt, freq2_filt, type_filt )
          ftmp(1) = ftmp(1)*H
        else
          ftmp(1) = 0._r4
        end if
       
        do k = 2,ind_middle
          w = dw*real(k-1)
          H = butterworth0( w, order_filt, freq1_filt, freq2_filt, type_filt )
          ftmp(k)          = ftmp(k)*H
          ftmp(nt_pad-k+2) = conjg(ftmp(k))
        end do
        
        if ( mod(nt_pad,2) == 0 ) then
          k = ind_middle+1
          w = dw*real(k-1)
          H = butterworth0( w, order_filt, freq1_filt, freq2_filt, type_filt )
          ftmp(k) = ftmp(k)*H
        end if
        
      else
      
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
        
      end if
     
      !
      call cfft1b(nt_pad,1,ftmp,nt_pad,wsave,lensav,work,lenwrk,iaux)
     
      !
      trace = tr1 + real( ftmp(1:nt) )
       
      !---------------------------------------------------------------------/
      deallocate( ftmp, stat=stal ); if ( stal/=0 ) stop ' ***** AE-1d: space_filter_1D ***** '
      deallocate( wind, stat=stal ); if ( stal/=0 ) stop ' ***** AE-3d: space_filter_1D ***** '
      deallocate( wsave, work, stat=stal ); if ( stal/=0 ) stop ' ***** AE-2d: space_filter_1D ***** '
   
    end subroutine space_filter_1D_win_r
    
    
    
    

   
  !***********************************************************************/
  end module m_dag_space_filter_1D_win




