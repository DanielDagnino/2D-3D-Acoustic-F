  !*********************************************************************/
  !** This code has been done in the Barcelona Center for Subsurface 
  !** Imaging (BCSI).
  !** Authors: Daniel Dagnino.
  !*********************************************************************/
  !*****     FFTPACK labrary is used in this code.     *****************/
  !*****     Authors: P. Swarztrauber and D. Valent    *****************/
  !*********************************************************************/
  module m_dag_fourier
    
    interface fourier_trans_1D
      module procedure fourier_trans_1D, fourier_trans_1D_r
    end interface

    interface inv_fourier_trans_1D
      module procedure inv_fourier_trans_1D, inv_fourier_trans_1D_r, inv_fourier_trans_1D_c, inv_fourier_trans_1D_c_r
    end interface
    
    interface w_grid
      module procedure w_grid, w_grid_r
    end interface
    
    contains
    !*********************************************************************/
    !*********************************************************************/
    !*****     Double     ************************************************/
    !*********************************************************************/
    !*********************************************************************/
    pure subroutine w_grid( w, dt, nt_pad, dw )
      use m_dag_data_kind
      implicit none
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r8
      integer, intent(in)      :: nt_pad
      real(prcsn), intent(in)  :: dt
      real(prcsn), intent(out) :: w(:)      
      real(prcsn), intent(out) :: dw      
      ! The variables which are generated inside the function.
      integer :: k
      integer :: ind_middle
!dir$ assume_aligned w(1):64
      
      !---------------------------------------------------------------------/
      if ( mod(nt_pad,2) /= 0 ) then
        ind_middle = (nt_pad-1)/2
      else
        ind_middle = nt_pad/2
      end if
      
      !---------------------------------------------------------------------/
      dw = 1._prcsn/( dt*real(nt_pad,prcsn) )
      w(1) = 0._prcsn
      do k=2,ind_middle+1
        w(k)          =  real(k-1,prcsn)*dw
        w(nt_pad-k+2) = -real(k-1,prcsn)*dw
      end do
      
    end subroutine w_grid
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine fourier_trans_1D( trace, trace_Fourier, nt, nt_pad, w, dt )
      use m_dag_data_kind
      implicit none
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r8
      integer, intent(in)           :: nt
      integer, optional, intent(in) :: nt_pad
      complex(prcsn), intent(out)  :: trace_Fourier(:)
      real(prcsn), intent(inout) :: trace(:)
      real(prcsn), optional, intent(in)  :: dt
      real(prcsn), optional, intent(out) :: w(:)
      ! The variables which are generated inside the function.
      integer :: k
      integer :: ind_middle
      integer :: lensav, lenwrk, iaux, nt_pad_aux
      real(prcsn), allocatable :: wsave(:), work(:)
      real(prcsn) :: dw
      integer :: stal
!dir$ assume_aligned trace_Fourier(1):64,trace(1):64,w(1):64
      
      !---------------------------------------------------------------------/
      if ( present(nt_pad) ) then
        nt_pad_aux = nt_pad
      else
        nt_pad_aux = nt
      end if
      
      !---------------------------------------------------------------------/
      if ( mod(nt_pad_aux,2) /= 0 ) then
        ind_middle = (nt_pad_aux-1)/2
      else
        ind_middle = nt_pad_aux/2
      end if
      
      !---------------------------------------------------------------------/
      dw = 1._prcsn/( dt*real(nt_pad_aux,prcsn) )
      
      if ( present(w) ) then
        w(1) = 0._prcsn
        do k=2,ind_middle+1
          w(k)              =  real(k-1,prcsn)*dw
          w(nt_pad_aux-k+2) = -real(k-1,prcsn)*dw
        end do
      end if
      
      !---------------------------------------------------------------------/
      lensav = 2*nt_pad_aux + ceiling(log(real(nt_pad_aux,prcsn))/log(2._prcsn)) + 4
      lenwrk = 2*nt_pad_aux
      
      allocate( wsave(lensav), work(lenwrk), stat=stal ); if ( stal/=0 ) stop ' ***** AE-1: fourier_trans_1D ***** '
    
      call dcfft1i(nt_pad_aux,wsave,lensav,iaux)
      
      !---------------------------------------------------------------------/
      trace_Fourier(1:nt) = cmplx( trace )
      if ( nt/=nt_pad_aux ) trace_Fourier((nt+1):nt_pad_aux) = cmplx( 0._prcsn )
      
      call dcfft1f(nt_pad_aux,1,trace_Fourier,nt_pad_aux,wsave,lensav,work,lenwrk,iaux)
      
      deallocate( wsave, work, stat=stal ); if ( stal/=0 ) stop ' ***** AE-1: fourier_trans_1D ***** '
      
    end subroutine fourier_trans_1D
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine inv_fourier_trans_1D( trace_Fourier, trace, nt, nt_pad )
      use m_dag_data_kind
      implicit none
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r8
      complex(prcsn), intent(inout) :: trace_Fourier(:)
      real(prcsn), intent(out) :: trace(:)
      integer, intent(in)           :: nt
      integer, optional, intent(in) :: nt_pad
      ! The variables which are generated inside the function.
      integer :: lensav, lenwrk, iaux, nt_pad_aux
      real(prcsn), allocatable :: wsave(:), work(:)
      integer :: stal
      complex(prcsn), allocatable :: trace_Fourier_aux(:)
      
      !---------------------------------------------------------------------/
      if ( present(nt_pad) ) then
        nt_pad_aux = nt_pad
      else
        nt_pad_aux = nt
      end if
      
      !---------------------------------------------------------------------/
      lensav = 2*nt_pad_aux + ceiling(log(real(nt_pad_aux,prcsn))/log(2._prcsn)) + 4
      lenwrk = 2*nt_pad_aux
      
      allocate( wsave(lensav), work(lenwrk), stat=stal ); if ( stal/=0 ) stop ' ***** AE-1: inv_fourier_trans_1D ***** '
      
      call dcfft1i(nt_pad_aux,wsave,lensav,iaux)
      
      !---------------------------------------------------------------------/
      allocate( trace_Fourier_aux(nt_pad_aux), stat=stal ); if ( stal/=0 ) stop ' ***** AE-2: inv_fourier_trans_1D ***** '
      
      trace_Fourier_aux = trace_Fourier
      
      call dcfft1b(nt_pad_aux,1,trace_Fourier_aux,nt_pad_aux,wsave,lensav,work,lenwrk,iaux)
      
      trace = real( trace_Fourier_aux(1:nt), prcsn )
      
      deallocate( trace_Fourier_aux, stat=stal ); if ( stal/=0 ) stop ' ***** AE-2: inv_fourier_trans_1D ***** '
      deallocate( wsave, work, stat=stal ); if ( stal/=0 ) stop ' ***** AE-1: inv_fourier_trans_1D ***** '
      
    end subroutine inv_fourier_trans_1D
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine inv_fourier_trans_1D_c( trace_Fourier, trace, nt, nt_pad )
      use m_dag_data_kind
      implicit none
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r8
      complex(prcsn), intent(inout) :: trace_Fourier(:)
      complex(prcsn), optional, intent(out) :: trace(:)
      integer, intent(in)           :: nt
      integer, optional, intent(in) :: nt_pad
      ! The variables which are generated inside the function.
      integer :: lensav, lenwrk, iaux, nt_pad_aux
      real(prcsn), allocatable :: wsave(:), work(:)
      integer :: stal
      complex(prcsn), allocatable :: trace_Fourier_aux(:)
      
      !---------------------------------------------------------------------/
      if ( present(nt_pad) ) then
        nt_pad_aux = nt_pad
      else
        nt_pad_aux = nt
      end if
      
      !---------------------------------------------------------------------/
      lensav = 2*nt_pad_aux + ceiling(log(real(nt_pad_aux,prcsn))/log(2._prcsn)) + 4
      lenwrk = 2*nt_pad_aux
      
      allocate( wsave(lensav), work(lenwrk), stat=stal ); if ( stal/=0 ) stop ' ***** AE-1: inv_fourier_trans_1D_c ***** '
      
      call dcfft1i(nt_pad_aux,wsave,lensav,iaux)
      
      !---------------------------------------------------------------------/
      allocate( trace_Fourier_aux(nt_pad_aux), stat=stal ); if ( stal/=0 ) stop ' ***** AE-2: inv_fourier_trans_1D_c ***** '
      
      trace_Fourier_aux = trace_Fourier
      
      call dcfft1b(nt_pad_aux,1,trace_Fourier_aux,nt_pad_aux,wsave,lensav,work,lenwrk,iaux)
      
      if ( present(trace) ) then
        trace = trace_Fourier_aux(1:nt)
      else
        trace_Fourier = trace_Fourier_aux(1:nt)
      end if
      
      deallocate( trace_Fourier_aux, stat=stal ); if ( stal/=0 ) stop ' ***** AE-2: inv_fourier_trans_1D_c ***** '
      deallocate( wsave, work, stat=stal ); if ( stal/=0 ) stop ' ***** AE-1: inv_fourier_trans_1D_c ***** '
      
    end subroutine inv_fourier_trans_1D_c
    
    
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*****     Real     **************************************************/
    !*********************************************************************/
    !*********************************************************************/
    pure subroutine w_grid_r( w, dt, nt_pad, dw )
      use m_dag_data_kind
      implicit none
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r4
      integer, intent(in)      :: nt_pad
      real(prcsn), intent(in)  :: dt
      real(prcsn), intent(out) :: w(:)      
      real(prcsn), intent(out) :: dw      
      ! The variables which are generated inside the function.
      integer :: k
      integer :: ind_middle
!dir$ assume_aligned w(1):32
      
      !---------------------------------------------------------------------/
      if ( mod(nt_pad,2) /= 0 ) then
        ind_middle = (nt_pad-1)/2
      else
        ind_middle = nt_pad/2
      end if
      
      !---------------------------------------------------------------------/
      dw = 1._prcsn/( dt*real(nt_pad,prcsn) )
      w(1) = 0._prcsn
      do k=2,ind_middle+1
        w(k)          =  real(k-1,prcsn)*dw
        w(nt_pad-k+2) = -real(k-1,prcsn)*dw
      end do
    end subroutine w_grid_r
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine fourier_trans_1D_r( trace, trace_Fourier, nt, nt_pad, w, dt )
      use m_dag_data_kind
      implicit none
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r4
      integer, intent(in)           :: nt
      integer, optional, intent(in) :: nt_pad
      complex(prcsn), intent(out)  :: trace_Fourier(:)
      real(prcsn), intent(inout) :: trace(:)
      real(prcsn), optional, intent(in)  :: dt
      real(prcsn), optional, intent(out) :: w(:)
      ! The variables which are generated inside the function.
      integer :: k
      integer :: ind_middle
      integer :: lensav, lenwrk, iaux, nt_pad_aux
      real(prcsn), allocatable :: wsave(:), work(:)
      real(prcsn) :: dw
      integer :: stal
!dir$ assume_aligned trace_Fourier(1):32,trace(1):32,w(1):32
      
      !---------------------------------------------------------------------/
      if ( present(nt_pad) ) then
        nt_pad_aux = nt_pad
      else
        nt_pad_aux = nt
      end if
      
      !---------------------------------------------------------------------/
      if ( mod(nt_pad_aux,2) /= 0 ) then
        ind_middle = (nt_pad_aux-1)/2
      else
        ind_middle = nt_pad_aux/2
      end if
      
      !---------------------------------------------------------------------/
      dw = 1._prcsn/( dt*real(nt_pad_aux,prcsn) )
      
      if ( present(w) ) then
        w(1) = 0._prcsn
        do k=2,ind_middle+1
          w(k)              =  real(k-1,prcsn)*dw
          w(nt_pad_aux-k+2) = -real(k-1,prcsn)*dw
        end do
      end if
      
      !---------------------------------------------------------------------/
      lensav = 2*nt_pad_aux + ceiling(log(real(nt_pad_aux,prcsn))/log(2._prcsn)) + 4
      lenwrk = 2*nt_pad_aux
      
      allocate( wsave(lensav), work(lenwrk), stat=stal ); if ( stal/=0 ) stop ' ***** AE-1: fourier_trans_1D ***** '
    
      call cfft1i(nt_pad_aux,wsave,lensav,iaux)
      
      !---------------------------------------------------------------------/
      trace_Fourier(1:nt) = cmplx( trace )
      if ( nt/=nt_pad_aux ) trace_Fourier((nt+1):nt_pad_aux) = cmplx( 0._prcsn )
      
      call cfft1f(nt_pad_aux,1,trace_Fourier,nt_pad_aux,wsave,lensav,work,lenwrk,iaux)
      
      deallocate( wsave, work, stat=stal ); if ( stal/=0 ) stop ' ***** AE-1: fourier_trans_1D ***** '
      
    end subroutine fourier_trans_1D_r
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine inv_fourier_trans_1D_r( trace_Fourier, trace, nt, nt_pad )
      use m_dag_data_kind
      implicit none
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r4
      complex(prcsn), intent(inout) :: trace_Fourier(:)
      real(prcsn), intent(out) :: trace(:)
      integer, intent(in)           :: nt
      integer, optional, intent(in) :: nt_pad
      ! The variables which are generated inside the function.
      integer :: lensav, lenwrk, iaux, nt_pad_aux
      real(prcsn), allocatable :: wsave(:), work(:)
      integer :: stal
      complex(prcsn), allocatable :: trace_Fourier_aux(:)
      
      !---------------------------------------------------------------------/
      if ( present(nt_pad) ) then
        nt_pad_aux = nt_pad
      else
        nt_pad_aux = nt
      end if
      
      !---------------------------------------------------------------------/
      lensav = 2*nt_pad_aux + ceiling(log(real(nt_pad_aux,prcsn))/log(2._prcsn)) + 4
      lenwrk = 2*nt_pad_aux
      
      allocate( wsave(lensav), work(lenwrk), stat=stal ); if ( stal/=0 ) stop ' ***** AE-1: inv_fourier_trans_1D ***** '
      
      call cfft1i(nt_pad_aux,wsave,lensav,iaux)
      
      !---------------------------------------------------------------------/
      allocate( trace_Fourier_aux(nt_pad_aux), stat=stal ); if ( stal/=0 ) stop ' ***** AE-2: inv_fourier_trans_1D ***** '
      
      trace_Fourier_aux = trace_Fourier
      
      call cfft1b(nt_pad_aux,1,trace_Fourier_aux,nt_pad_aux,wsave,lensav,work,lenwrk,iaux)
      
      trace = real( trace_Fourier_aux(1:nt), prcsn )
      
      deallocate( trace_Fourier_aux, stat=stal ); if ( stal/=0 ) stop ' ***** AE-2: inv_fourier_trans_1D ***** '
      deallocate( wsave, work, stat=stal ); if ( stal/=0 ) stop ' ***** AE-1: inv_fourier_trans_1D ***** '
      
    end subroutine inv_fourier_trans_1D_r
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine inv_fourier_trans_1D_c_r( trace, trace_Fourier, nt, nt_pad )
      use m_dag_data_kind
      implicit none
      ! The variables which are passed to the function.
      integer, parameter :: prcsn = r4
      complex(prcsn), intent(inout) :: trace_Fourier(:)
      complex(prcsn), optional, intent(out) :: trace(:)
      integer, intent(in)           :: nt
      integer, optional, intent(in) :: nt_pad
      ! The variables which are generated inside the function.
      integer :: lensav, lenwrk, iaux, nt_pad_aux
      real(prcsn), allocatable :: wsave(:), work(:)
      integer :: stal
      complex(prcsn), allocatable :: trace_Fourier_aux(:)
      
      !---------------------------------------------------------------------/
      if ( present(nt_pad) ) then
        nt_pad_aux = nt_pad
      else
        nt_pad_aux = nt
      end if
      
      !---------------------------------------------------------------------/
      lensav = 2*nt_pad_aux + ceiling(log(real(nt_pad_aux,prcsn))/log(2._prcsn)) + 4
      lenwrk = 2*nt_pad_aux
      
      allocate( wsave(lensav), work(lenwrk), stat=stal ); if ( stal/=0 ) stop ' ***** AE-1: inv_fourier_trans_1D_c ***** '
      
      call cfft1i(nt_pad_aux,wsave,lensav,iaux)
      
      !---------------------------------------------------------------------/
      allocate( trace_Fourier_aux(nt_pad_aux), stat=stal ); if ( stal/=0 ) stop ' ***** AE-2: inv_fourier_trans_1D_c ***** '
      
      trace_Fourier_aux = trace_Fourier
      
      call cfft1b(nt_pad_aux,1,trace_Fourier_aux,nt_pad_aux,wsave,lensav,work,lenwrk,iaux)
      
      if ( present(trace) ) then
        trace = trace_Fourier_aux(1:nt)
      else
        trace_Fourier = trace_Fourier_aux(1:nt)
      end if
      
      deallocate( trace_Fourier_aux, stat=stal ); if ( stal/=0 ) stop ' ***** AE-2: inv_fourier_trans_1D_c ***** '
      deallocate( wsave, work, stat=stal ); if ( stal/=0 ) stop ' ***** AE-1: inv_fourier_trans_1D_c ***** '
      
    end subroutine inv_fourier_trans_1D_c_r
    
  !***********************************************************************/
  end module m_dag_fourier
  
  
  
  
  