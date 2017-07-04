  
  !*********************************************************************/
  !*****     Author: D. Dagnino     ************************************/
  !*********************************************************************/
  module m_dag_filter_1D_win
    
    interface filter_1D_win
      module procedure filter_1D_win, filter_1D_win_r
    end interface
    
    contains
    !*********************************************************************/
    !*********************************************************************/
    !*****     Double     ************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine filter_1D_win( trace, freq1_filt, freq2_filt, type_filt, order_filt, zero_padding, zero_delay, nt, dt, t_add )
      use m_dag_data_kind
      use m_dag_transfer_funct
      use m_dag_mat
      implicit none
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r8
      real(prcsn), intent(in) :: freq1_filt, freq2_filt, dt, t_add
      integer, intent(in)     :: type_filt, order_filt, zero_padding, nt
      logical, intent(in)     :: zero_delay
      real(prcsn), intent(inout) :: trace(:)
      ! The variables which are generated inside the function.
      integer     :: k, nt_pad, ind_middle
      integer     :: lensav, lenwrk, iaux
      complex(prcsn), allocatable :: ftmp(:)
      complex(prcsn) :: cH
      real(prcsn)    :: H, w, dw, alpha
      real(prcsn), allocatable :: wind(:)
      real(prcsn), allocatable :: wsave(:), work(:)
      integer :: stal
!dir$ assume_aligned trace(1):64
      
      !---------------------------------------------------------------------/
      ! 
      nt_pad = (1+zero_padding)*nt
      
      ! 
      if ( mod(nt_pad,2) /= 0 ) then
        ind_middle = (nt_pad-1)/2
      else
        ind_middle = nt_pad/2
      end if
      
      lensav = 2*nt_pad + ceiling(log(real(nt_pad,prcsn))/log(2._prcsn)) + 4
      lenwrk = 2*nt_pad
      
      allocate( wsave(lensav), work(lenwrk), stat=stal ); if ( stal/=0 ) stop ' ***** AE-2: filter_1D_win ***** '
      
      ! 
      call dcfft1i(nt_pad,wsave,lensav,iaux)

      !---------------------------------------------------------------------/
      ! 
      alpha = 2._prcsn*t_add/(dt*real(nt,prcsn))
      if ( alpha > 1._prcsn ) stop '***** ERROR - filter_1D_win: alpha > 1._prcsn *****'
      
      ! 
      do k=1,nt
        if ( k-1 <= floor(alpha*real(nt-1,prcsn)/2._prcsn) ) then
          wind(k) = 0.5_prcsn*( 1._prcsn + cos( pi*( 2._prcsn*real(k-1,prcsn)/(alpha*real(nt-1,prcsn)) -1._prcsn ) ) )
        else if ( k-1 >= floor(real(nt-1,prcsn)*(1._prcsn-alpha/2._prcsn)) ) then
          wind(k) = 0.5_prcsn*( 1._prcsn + cos( pi*( 2._prcsn*real(k-1,prcsn)/(alpha*real(nt-1,prcsn)) -2._prcsn/alpha +1._prcsn ) ) )
        else
          wind(k) = 1._prcsn
        end if
      end do
      
      !---------------------------------------------------------------------/
      ! 
      allocate( ftmp(nt_pad), stat=stal ); if ( stal/=0 ) stop ' ***** AE-1: filter_1D_win ***** '
      allocate( wind(nt), stat=stal ); if ( stal/=0 ) stop ' ***** AE-3: filter_1D_win ***** '
      
      ftmp(1:nt) = cmplx( wind*trace )
      if ( nt_pad /= nt ) ftmp(nt+1:nt_pad) = 0._prcsn
      
      !
      call dcfft1f(nt_pad,1,ftmp,nt_pad,wsave,lensav,work,lenwrk,iaux)
      
      !---------------------------------------------------------------------/
      ! 
      dw = 1._prcsn/( dt*real(nt_pad,prcsn) )
      
      ! 
      if ( zero_delay ) then
      
        if ( type_filt == 1 ) then
          H = butterworth0( 0._prcsn, order_filt, freq1_filt, freq2_filt, type_filt )
          ftmp(1) = ftmp(1)*H
        else
          ftmp(1) = 0._prcsn
        end if
        
        do k = 2,ind_middle
          w = dw*real(k-1,prcsn)
          H = butterworth0( w, order_filt, freq1_filt, freq2_filt, type_filt )
          ftmp(k)          = ftmp(k)*H
          ftmp(nt_pad-k+2) = conjg( ftmp(k) )
        end do
        
        if ( mod(nt_pad,2) == 0 ) then
          k = ind_middle+1
          w = dw*real(k-1,prcsn)
          H = butterworth0( w, order_filt, freq1_filt, freq2_filt, type_filt )
          ftmp(k) = ftmp(k)*H
        end if
        
      else
      
        if ( type_filt == 1 ) then
          cH = butterworth( 0._prcsn, order_filt, freq1_filt, freq2_filt, type_filt )
          ftmp(1) = ftmp(1)*cH
        else
          ftmp(1) = 0._prcsn
        end if
        
        do k = 2,ind_middle
          w = dw*real(k-1,prcsn)
          cH = butterworth( w, order_filt, freq1_filt, freq2_filt, type_filt )
          ftmp(k)          = ftmp(k)*cH
          ftmp(nt_pad-k+2) = conjg( ftmp(k) )
        end do
        
        if ( mod(nt_pad,2) == 0 ) then
          k = ind_middle+1
          w = dw*real(k-1,prcsn)
          cH = butterworth( w, order_filt, freq1_filt, freq2_filt, type_filt )
          ftmp(k) = ftmp(k)*cH
        end if
        
      end if
     
      !---------------------------------------------------------------------/
      !
      call dcfft1b(nt_pad,1,ftmp,nt_pad,wsave,lensav,work,lenwrk,iaux)
      
      !
      trace = real( ftmp(1:nt) )
      
      !---------------------------------------------------------------------/
      deallocate( ftmp, stat=stal ); if ( stal/=0 ) stop ' ***** AE-1d: filter_1D_win ***** '
      deallocate( wsave, work, stat=stal ); if ( stal/=0 ) stop ' ***** AE-2d: filter_1D_win ***** '
      deallocate( wind, stat=stal ); if ( stal/=0 ) stop ' ***** AE-3d: filter_1D_win ***** '
      
    end subroutine filter_1D_win
   
   
   
   
   
    !*********************************************************************/
    !*********************************************************************/
    !*****     Real     **************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine filter_1D_win_r( trace, freq1_filt, freq2_filt, type_filt, order_filt, zero_padding, zero_delay, nt, dt, t_add )
      use m_dag_data_kind
      use m_dag_transfer_funct
      use m_dag_mat
      implicit none
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r4
      real(prcsn), intent(in) :: freq1_filt, freq2_filt, dt, t_add
      integer, intent(in)     :: type_filt, order_filt, zero_padding, nt
      logical, intent(in)     :: zero_delay
      real(prcsn), intent(inout) :: trace(:)
      ! The variables which are generated inside the function.
      integer     :: k, nt_pad, ind_middle
      integer     :: lensav, lenwrk, iaux
      complex(prcsn), allocatable :: ftmp(:)
      complex(prcsn) :: cH
      real(prcsn)    :: H, w, dw, alpha
      real(prcsn), allocatable :: wind(:)
      real(prcsn), allocatable :: wsave(:), work(:)
      integer :: stal
!dir$ assume_aligned trace(1):32
      
      !---------------------------------------------------------------------/
      ! 
      nt_pad = (1+zero_padding)*nt
      
      ! 
      if ( mod(nt_pad,2) /= 0 ) then
        ind_middle = (nt_pad-1)/2
      else
        ind_middle = nt_pad/2
      end if
      
      lensav = 2*nt_pad + ceiling(log(real(nt_pad,prcsn))/log(2._prcsn)) + 4
      lenwrk = 2*nt_pad
      
      allocate( wsave(lensav), work(lenwrk), stat=stal ); if ( stal/=0 ) stop ' ***** AE-2: filter_1D_win ***** '
      
      ! 
      call cfft1i(nt_pad,wsave,lensav,iaux)

      !---------------------------------------------------------------------/
      ! 
      alpha = 2._prcsn*t_add/(dt*real(nt,prcsn))
      if ( alpha > 1._prcsn ) stop '***** ERROR - filter_1D_win: alpha > 1._prcsn *****'
      
      ! 
      do k=1,nt
        if ( k-1 <= floor(alpha*real(nt-1,prcsn)/2._prcsn) ) then
          wind(k) = 0.5_prcsn*( 1._prcsn + cos( pi*( 2._prcsn*real(k-1,prcsn)/(alpha*real(nt-1,prcsn)) -1._prcsn ) ) )
        else if ( k-1 >= floor(real(nt-1,prcsn)*(1._prcsn-alpha/2._prcsn)) ) then
          wind(k) = 0.5_prcsn*( 1._prcsn + cos( pi*( 2._prcsn*real(k-1,prcsn)/(alpha*real(nt-1,prcsn)) -2._prcsn/alpha +1._prcsn ) ) )
        else
          wind(k) = 1._prcsn
        end if
      end do
      
      !---------------------------------------------------------------------/
      ! 
      allocate( ftmp(nt_pad), stat=stal ); if ( stal/=0 ) stop ' ***** AE-1: filter_1D_win ***** '
      allocate( wind(nt), stat=stal ); if ( stal/=0 ) stop ' ***** AE-3: filter_1D_win ***** '
      
      ftmp(1:nt) = cmplx( wind*trace )
      if ( nt_pad /= nt ) ftmp(nt+1:nt_pad) = 0._prcsn
      
      !
      call cfft1f(nt_pad,1,ftmp,nt_pad,wsave,lensav,work,lenwrk,iaux)
      
      !---------------------------------------------------------------------/
      ! 
      dw = 1._prcsn/( dt*real(nt_pad,prcsn) )
      
      ! 
      if ( zero_delay ) then
      
        if ( type_filt == 1 ) then
          H = butterworth0( 0._prcsn, order_filt, freq1_filt, freq2_filt, type_filt )
          ftmp(1) = ftmp(1)*H
        else
          ftmp(1) = 0._prcsn
        end if
        
        do k = 2,ind_middle
          w = dw*real(k-1,prcsn)
          H = butterworth0( w, order_filt, freq1_filt, freq2_filt, type_filt )
          ftmp(k)          = ftmp(k)*H
          ftmp(nt_pad-k+2) = conjg( ftmp(k) )
        end do
        
        if ( mod(nt_pad,2) == 0 ) then
          k = ind_middle+1
          w = dw*real(k-1,prcsn)
          H = butterworth0( w, order_filt, freq1_filt, freq2_filt, type_filt )
          ftmp(k) = ftmp(k)*H
        end if
        
      else
      
        if ( type_filt == 1 ) then
          cH = butterworth( 0._prcsn, order_filt, freq1_filt, freq2_filt, type_filt )
          ftmp(1) = ftmp(1)*cH
        else
          ftmp(1) = 0._prcsn
        end if
        
        do k = 2,ind_middle
          w = dw*real(k-1,prcsn)
          cH = butterworth( w, order_filt, freq1_filt, freq2_filt, type_filt )
          ftmp(k)          = ftmp(k)*cH
          ftmp(nt_pad-k+2) = conjg( ftmp(k) )
        end do
        
        if ( mod(nt_pad,2) == 0 ) then
          k = ind_middle+1
          w = dw*real(k-1,prcsn)
          cH = butterworth( w, order_filt, freq1_filt, freq2_filt, type_filt )
          ftmp(k) = ftmp(k)*cH
        end if
        
      end if
     
      !---------------------------------------------------------------------/
      !
      call cfft1b(nt_pad,1,ftmp,nt_pad,wsave,lensav,work,lenwrk,iaux)
      
      !
      trace = real( ftmp(1:nt) )
      
      !---------------------------------------------------------------------/
      deallocate( ftmp, stat=stal ); if ( stal/=0 ) stop ' ***** AE-1d: filter_1D_win ***** '
      deallocate( wsave, work, stat=stal ); if ( stal/=0 ) stop ' ***** AE-2d: filter_1D_win ***** '
      deallocate( wind, stat=stal ); if ( stal/=0 ) stop ' ***** AE-3d: filter_1D_win ***** '
      
    end subroutine filter_1D_win_r
   
  !***********************************************************************/
  end module m_dag_filter_1D_win




