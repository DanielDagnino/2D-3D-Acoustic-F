  !*********************************************************************/
  !*****     Author: HPC Stockolm summer school     ********************/
  !*****     Author: adapted by D. Dagnino     *************************/
  !*********************************************************************/
  module m_dag_time
    
    contains
    !*********************************************************************/
    !*********************************************************************/
    real(8) function rtc()
      implicit none
      integer :: icnt,irate
      real(8), save :: scaling
      logical, save :: scale = .true.
      
      call system_clock(icnt,irate)
      
      if ( scale ) then
        scaling = 1.d0/real(irate,8)
        scale = .false.
      end if
      
      rtc = icnt * scaling
      
    end function rtc
    
  !***********************************************************************/
  end module m_dag_time