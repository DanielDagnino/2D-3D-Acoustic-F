  !*********************************************************************/
  !*****     Author: D. Dagnino     ************************************/
  !*********************************************************************/
  module m_dag_space_2_filter_1D_win_zp_extra
    
    interface space_2_filter_1D_win_zp_extra
      module procedure space_2_filter_1D_win_zp_extra, space_2_filter_1D_win_zp_extra_r
    end interface
    
    contains
    !*********************************************************************/
    !*********************************************************************/
    !*****     Double     ************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine space_2_filter_1D_win_zp_extra( trace, freq1_filt, freq2_filt, p0, p1, type_filt, order_filt, zero_padding, &
                                      nt, dt, ext, extra )
      use m_dag_data_kind
      use m_dag_transfer_funct
      use m_dag_mat
      use m_dag_factor_number
      implicit none
      
      ! The variables which are passed to the function.
      real(r8), intent(in)    :: freq1_filt, freq2_filt
      real(r8), intent(in)    :: p0, p1
      integer, intent(in)     :: type_filt
      integer, intent(in)     :: order_filt
      integer, intent(in)     :: zero_padding
      logical, intent(in)     :: ext
      integer, intent(in)     :: nt
      real(r8), intent(in)    :: dt
      real(r8), intent(inout) :: trace(:)
      real(r8), intent(in)    :: extra
    
      ! The variables which are generated inside the function.
      integer     :: k
      integer     :: nt_pad
      real(r8)    :: extra2
      complex(r8) :: cH, caux
      real(r8)    :: w, dw, alpha, inct
      integer     :: lensav, lenwrk, iaux
      real(r8), allocatable :: wsave(:), work(:)
      complex(r8), allocatable :: ftmp(:)
      real(r8), allocatable :: wind(:)
      integer :: stal
      integer :: ind_middle, it0, it1
      real(r8) :: tr1, tr2, T01, slp0, slp1
!dir$ assume_aligned trace(1):64
      
      !*********************************************************************/
!       nt_pad = (1+zero_padding)*nt
      if ( zero_padding==0 ) then
        nt_pad = nt
      else
        nt_pad = near_small_prime_factoriz( nt )*(1+zero_padding)
      end if
      
!       if ( mod(nt_pad,2) /= 0 ) stop '***** ERROR - space_2_filter_1D: mod(nt_pad,2) /= 0 *****'
    
      allocate( ftmp(nt_pad), stat=stal ); if ( stal/=0 ) stop ' ***** AE-1: space_2_filter_1D ***** '
      allocate( wind(nt_pad), stat=stal ); if ( stal/=0 ) stop ' ***** AE-3: space_2_filter_1D ***** '
    
      lensav = 2*nt_pad + ceiling(dlog(dble(nt_pad))/dlog(2._r8)) + 4
      lenwrk = 2*nt_pad
    
      allocate( wsave(lensav), work(lenwrk), stat=stal ); if ( stal/=0 ) stop ' ***** AE-2: space_2_filter_1D ***** '
      
      call dcfft1i(nt_pad,wsave,lensav,iaux)
      
      !*********************************************************************/
      !
      if ( mod(nt_pad,2) /= 0 ) then
        ind_middle = (nt_pad-1)/2
      else
        ind_middle = nt_pad/2
      end if
     
      !*********************************************************************/
      !
      inct = 0.5_r8*(nt_pad-nt)
      alpha = 2._r8*inct/dble(nt_pad)
     
      !*********************************************************************/
      !
      tr1 = trace(1)
      trace = trace-tr1
       
      ftmp(1:nt) = cmplx( trace )
      if ( nt_pad /= nt ) ftmp(nt+1:nt_pad) = 0._r8
      
      tr2 = trace(nt)
      do k=1,nt
        ftmp(k) = ftmp(k) - dble(k-1)*tr2/dble(nt-1)
      end do
     
      !*********************************************************************/
      !
      if ( alpha > 1._r8 ) stop '***** ERROR - space_2_filter_1D: alpha > 1._r8 *****'
      if ( zero_padding == 0 ) stop '***** ERROR - space_2_filter_1D: zero_padding == 0 *****'
     
      !
      do k=1,nt_pad
        if ( k <= nt ) then
          wind(k) = 1._r8
        else
          wind(k) = 0.5_r8*( 1._r8 + dcos( pi*( 2._r8*dble(k-nt)/dble(nt_pad-nt+1)) ) )
        end if
      end do
      
      !*********************************************************************/
      !
      T01 = dble(nt_pad-nt)
     
      it0 = floor(p0*dble(nt))     
      slp0 = 0._r8
      do k=1,it0
        slp0 = slp0 + ftmp(nt-k+1)-ftmp(nt-k+1-it0)
      end do
      slp0 = slp0/dble(it0)
      
      it1 = floor(p1*dble(nt))     
      slp1 = 0._r8
      do k=1,it1
        slp1 = slp1 + ftmp(k+it1)-ftmp(k)
      end do
      slp1 = slp1/dble(it1)
     
      !
      do k=1,nt_pad
        if ( k > nt .and. k <= ind_middle ) then
          ftmp(k) = slp1*dsin( 2._r8*pi*dble(k-nt)/T01 )
        else if ( k > ind_middle ) then
          ftmp(k) = slp0*dsin( 2._r8*pi*dble(k-nt_pad-1)/T01 )
        end if
      end do
      
      !*********************************************************************/
      !
      ftmp = cmplx( wind*ftmp )
      
      !*********************************************************************/
      !*********************************************************************/
      !
      dw   = 1._r8/( dt*dble(nt_pad) )
      
      !---------------------------------------------------------------------/
      call dcfft1f(nt_pad,1,ftmp,nt_pad,wsave,lensav,work,lenwrk,iaux)
      
      !----------------------/
      extra2 = 2._r8*pi*dble(int(anint( extra/dt )))*dw*dt
      do k = 1,nt_pad
        ftmp(k) = ftmp(k)*exp(-cmplx(0._r8,extra2*dble(k-1)))
      end do
      !----------------------/
      
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
        ftmp(nt_pad-k+2) = conjg(ftmp(k))
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
        ftmp(nt_pad-k+2) = conjg(ftmp(k))
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
      
      !*********************************************************************/
      !
      if ( ext ) then
       
        do k=1,nt
          trace(k) = trace(k) + dble(k-1)*tr2/dble(nt-1)
        end do
        
        trace = tr1 + trace
       
      end if
      
      !*********************************************************************/
      !*********************************************************************/
      deallocate( ftmp, stat=stal ); if ( stal/=0 ) stop ' ***** AE-1d: space_2_filter_1D ***** '
      deallocate( wind, stat=stal ); if ( stal/=0 ) stop ' ***** AE-3d: space_2_filter_1D ***** '
      deallocate( wsave, work, stat=stal ); if ( stal/=0 ) stop ' ***** AE-2d: space_2_filter_1D ***** '
  
    end subroutine space_2_filter_1D_win_zp_extra
  
  
  
  
  
    !*********************************************************************/
    !*********************************************************************/
    !*****     Real     **************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine space_2_filter_1D_win_zp_extra_r( trace, freq1_filt, freq2_filt, p0, p1, type_filt, order_filt, zero_padding, &
                                      nt, dt, ext, extra )
      use m_dag_data_kind
      use m_dag_mat
      use m_dag_transfer_funct
      use m_dag_mat
      use m_dag_factor_number
      implicit none
    
      ! The variables which are passed to the function.
      real(r4), intent(in)    :: freq1_filt, freq2_filt
      real(r4), intent(in)    :: p0, p1
      integer, intent(in)     :: type_filt
      integer, intent(in)     :: order_filt
      integer, intent(in)     :: zero_padding
      logical, intent(in)     :: ext
      integer, intent(in)     :: nt
      real(r4), intent(in)    :: dt
      real(r4), intent(inout) :: trace(:)
      real(r4), intent(in)    :: extra
      
      ! The variables which are generated inside the function.
      integer     :: k
      integer     :: nt_pad
      real(r4)    :: extra2
      complex(r4) :: cH, caux
      real(r4)    :: w, dw, alpha, inct
      integer     :: lensav, lenwrk, iaux
      real(r4), allocatable :: wsave(:), work(:)
      complex(r4), allocatable :: ftmp(:)
      real(r4), allocatable :: wind(:)
      integer :: stal
      integer :: ind_middle, it0, it1
      real(r4) :: tr1, tr2, T01, slp0, slp1
!dir$ assume_aligned trace(1):32
     
      !*********************************************************************/
!       nt_pad = (1+zero_padding)*nt
      if ( zero_padding==0 ) then
        nt_pad = nt
      else
        nt_pad = near_small_prime_factoriz( nt )*(1+zero_padding)
      end if
      
!       if ( mod(nt_pad,2) /= 0 ) stop '***** ERROR - space_2_filter_1D: mod(nt_pad,2) /= 0 *****'
    
      allocate( ftmp(nt_pad), stat=stal ); if ( stal/=0 ) stop ' ***** AE-1: space_2_filter_1D ***** '
      allocate( wind(nt_pad), stat=stal ); if ( stal/=0 ) stop ' ***** AE-3: space_2_filter_1D ***** '
    
      lensav = 2*nt_pad + ceiling(log(real(nt_pad))/log(2._r4)) + 4
      lenwrk = 2*nt_pad
    
      allocate( wsave(lensav), work(lenwrk), stat=stal ); if ( stal/=0 ) stop ' ***** AE-2: space_2_filter_1D ***** '
    
      call cfft1i(nt_pad,wsave,lensav,iaux)
     
      !*********************************************************************/
      !
      if ( mod(nt_pad,2) /= 0 ) then
        ind_middle = (nt_pad-1)/2
      else
        ind_middle = nt_pad/2
      end if
     
      !*********************************************************************/
      !
      inct = 0.5_r4*(nt_pad-nt)
      alpha = 2._r4*inct/real(nt_pad)
     
      !*********************************************************************/
      !
      tr1 = trace(1)
      trace = trace-tr1
       
      ftmp(1:nt) = cmplx( trace )
      if ( nt_pad /= nt ) ftmp(nt+1:nt_pad) = 0._r4
      
      tr2 = trace(nt)
      do k=1,nt
        ftmp(k) = ftmp(k) - real(k-1)*tr2/real(nt-1)
      end do
     
      !*********************************************************************/
      !
      if ( alpha > 1._r4 ) stop '***** ERROR - space_2_filter_1D: alpha > 1._r4 *****'
      if ( zero_padding == 0 ) stop '***** ERROR - space_2_filter_1D: zero_padding == 0 *****'
     
      !
      do k=1,nt_pad
        if ( k <= nt ) then
          wind(k) = 1._r4
        else
          wind(k) = 0.5_r4*( 1._r4 + dcos( pi*( 2._r4*real(k-nt)/real(nt_pad-nt+1)) ) )
        end if
      end do
      
      !*********************************************************************/
      !
      T01 = real(nt_pad-nt)
     
      it0 = floor(p0*real(nt))     
      slp0 = 0._r4
      do k=1,it0
        slp0 = slp0 + ftmp(nt-k+1)-ftmp(nt-k+1-it0)
      end do
      slp0 = slp0/real(it0)
      
      it1 = floor(p1*real(nt))     
      slp1 = 0._r4
      do k=1,it1
        slp1 = slp1 + ftmp(k+it1)-ftmp(k)
      end do
      slp1 = slp1/real(it1)
     
      !
      do k=1,nt_pad
        if ( k > nt .and. k <= ind_middle ) then
          ftmp(k) = slp1*sin( 2._r4*pi*real(k-nt)/T01 )
        else if ( k > ind_middle ) then
          ftmp(k) = slp0*sin( 2._r4*pi*real(k-nt_pad-1)/T01 )
        end if
      end do
      
      !*********************************************************************/
      !
      ftmp = cmplx( wind*ftmp )
     
      !*********************************************************************/
      !*********************************************************************/
      !
      dw   = 1._r4/( dt*real(nt_pad) )
      
      !---------------------------------------------------------------------/
      call cfft1f(nt_pad,1,ftmp,nt_pad,wsave,lensav,work,lenwrk,iaux)
      
      !----------------------/
      extra2 = 2._r4*pi*real(int(anint( extra/dt )))*dw*dt
      do k = 1,nt_pad
        ftmp(k) = ftmp(k)*exp(-cmplx(0._r4,extra2*real(k-1)))
      end do
      !----------------------/
      
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
    
      !*********************************************************************/
      !
      if ( ext ) then
       
        do k=1,nt
          trace(k) = trace(k) + real(k-1)*tr2/real(nt-1)
        end do
        
        trace = tr1 + trace
       
      end if
      
      !*********************************************************************/
      !*********************************************************************/
      deallocate( ftmp, stat=stal ); if ( stal/=0 ) stop ' ***** AE-1d: space_2_filter_1D ***** '
      deallocate( wind, stat=stal ); if ( stal/=0 ) stop ' ***** AE-3d: space_2_filter_1D ***** '
      deallocate( wsave, work, stat=stal ); if ( stal/=0 ) stop ' ***** AE-2d: space_2_filter_1D ***** '
  
    end subroutine space_2_filter_1D_win_zp_extra_r
  
  
  
  
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
   
   
   

  
  !***********************************************************************/
  end module m_dag_space_2_filter_1D_win_zp_extra





