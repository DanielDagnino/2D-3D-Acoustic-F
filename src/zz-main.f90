  !*********************************************************************/
  !** This code has been done in the Barcelona Center for Subsurface 
  !** Imaging (BCSI).
  !** Goal: Acoustic Full Wave Form Inversion with real data.
  !** Authors: Daniel Dagnino.
  !*********************************************************************/
  !** Short resume of this code/file:
  !** Main code. The code is basically (without writes) the following:
  !** 
  !** Read initial data.
  !** Split the work in the cores + QC.
  !** Read shot labels + QC.
  !** Read model + QC.
  !** Read frequency strategy.
  !** Call the main FWI algorithm.
  !** 
  !*********************************************************************/
  program main
    use m_mpi
    use m_dag_split_work
    use m_data_kind
    use m_geo
    use m_sg_type
    use m_sg, only: rmv_sg
    use m_sg_data, only: folder_strat, check_exit_folders_files
    use m_model
    use m_solv3Dahvc_allo_coll
    use m_solv3Dahvc_data_kind
    use m_solv3Dahvc_time
    use m_fwi
    use m_get_model_and_QC
    use m_get_initial_data
    use m_type_inversion
    use m_get_id_shots
    use m_get_freq_inv
    use m_filter, only:filter
    use m_Hessian, only: use_iH
    use m_inv_data, only: strat_inv
    use m_work, only: n_shot
    use m_get_strategy
#ifdef usempi
    use mpi
#endif
    implicit none
    integer, allocatable :: n_sg(:)
    type(shot_gather), allocatable :: shgat(:)
    integer :: stal
    
!     integer :: isg
    integer :: n_freq, isg
    real(sp), allocatable :: freq_0_k(:), freq_k(:)
    integer, allocatable :: iter_k(:)
    integer :: save_every_max
    real(sp) :: t1_inv_aux, t1_inv_max
    
    !**********************************************************!
    ! Compiler check.
#ifdef __INTEL_COMPILER
#elif __GFORTRAN__
#else
    if ( rank == 0 ) write(*,*) ' *** WARNING: We have never checked the code with this compiler !!! *** '
#endif
    
    !**********************************************************!
    ! If no MPI:
    n_procs = 1
    rank    = 0
    ierr    = 0
    mpi_fwi = 0
    
    ! Otherwise:
    ! Let the system do what it needs to start up MPI.
#ifdef usempi
    call MPI_INIT( ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, n_procs, ierr )
    call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
    ! Define the precision of the MPI variables.
    call prec_mpi_fwi( sp )
#endif
    
    !**********************************************************!
    ! MPI?
#ifdef usempi
    if ( rank == 0 ) write(*,*) '*** Multi core FWI ***********************************************************'
#else
    if ( rank == 0 ) write(*,*) '*** Single core FWI **********************************************************'
#endif
    
    !**********************************************************!
    ! Reading initial_data + check.
    if ( rank == 0 ) write(*,*) 'get_initial_data'
    call get_initial_data( n_shot, t1_inv_aux, dx )
    
    if ( rank == 0 ) call check_exit_folders_files()
    
    !**********************************************************!
    ! Read strategy info.
    if ( rank == 0 ) write(*,*) 'get_strategy'
    call get_strategy( folder_strat )
    
    ! Support of the source.
    p0_cut  = 0.05_sp   ! start sou at this accumulated weigth
    p1_cut  = 0.05_sp   ! end sou at this accumulated end weigth
    
    !**********************************************************!
#ifdef usempi
    if ( n_shot < n_procs ) stop " *** ERROR main: Dude, I am not going to waste time! n_shot < n_procs"
#endif
    
    ! Split the work in cores.
    allocate( n_sg(n_procs), stat=stal ); if ( stal/=0 ) stop 'AE shot'
    call split_work( n_shot, n_procs, n_sg )
    allocate( shgat(n_sg(rank+1)), stat=stal ); if ( stal/=0 ) stop 'AE shot'
    
    shgat(:)%t1_inv = t1_inv_aux
    
    !**********************************************************!
    ! Reading shots identifiers.
    if ( rank == 0 ) write(*,*) 'get_id_shots'
    call get_id_shots( folder_strat, n_shot, n_sg, shgat )
    
    ! Different t1_inv (and save_every=>only for OBS) for each shgat.
    different_time = .true.
    if ( different_time ) then
      
      ! Maximum save_every and t1_inv.
      save_every_max = save_every
      t1_inv_max = t1_inv_aux
      
      if ( geometry_type == 2 ) then
        do isg=1,n_sg(rank+1)
          ! Time to be inverted for each shgat.
          if ( shgat(isg)%id == 1 ) then
            shgat(isg)%t1_inv = 16.5
!             shgat(isg)%t1_inv = 15.0
          else if ( shgat(isg)%id == 2 ) then
            shgat(isg)%t1_inv = 18.0
!             shgat(isg)%t1_inv = 15.0
          else if ( shgat(isg)%id == 3 ) then
            shgat(isg)%t1_inv = 10.5
!             shgat(isg)%t1_inv = 10.5
          else if ( shgat(isg)%id == 4 ) then
            shgat(isg)%t1_inv = 13.5
!             shgat(isg)%t1_inv = 13.5
          else if ( shgat(isg)%id == 5 ) then
            shgat(isg)%t1_inv = 17.5
!             shgat(isg)%t1_inv = 15.0
          else if ( shgat(isg)%id == 6 ) then
            shgat(isg)%t1_inv = 15.5
!             shgat(isg)%t1_inv = 15.0
          else if ( shgat(isg)%id == 7 ) then
            shgat(isg)%t1_inv = 13.0
!             shgat(isg)%t1_inv = 13.0
          else if ( shgat(isg)%id == 8 ) then
            shgat(isg)%t1_inv = 16.5
!             shgat(isg)%t1_inv = 15.0
          else if ( shgat(isg)%id == 9 ) then
            shgat(isg)%t1_inv = 15.5
!             shgat(isg)%t1_inv = 15.0
          else if ( shgat(isg)%id == 11 ) then
            shgat(isg)%t1_inv = 21.0
!             shgat(isg)%t1_inv = 15.0
          else if ( shgat(isg)%id == 12 ) then
            shgat(isg)%t1_inv = 19.0
!             shgat(isg)%t1_inv = 15.0
          else
            stop ' ***** ERROR - main: It lacks a shgat(isg)%id ***** '
          end if
          
          ! Minimum field saving.
          save_every = max( save_every, ceiling(shgat(isg)%t1_inv*(dble(save_every_max)/t1_inv_max)) )
        end do
      end if
      
    end if
    
    !**********************************************************!
    ! Core who will write the outputs.
    ! It is not necessary the be the zero rank.
!     rank_wrt = n_procs-1
    rank_wrt = 0
!     if ( rank == rank_wrt ) isg_wrt = n_sg(1)
    if ( rank == rank_wrt ) isg_wrt = n_sg(rank_wrt+1)
    
    !**********************************************************!
    ! Reading frequency strategy.
    if ( rank == 0 ) write(*,*) 'get_freq_inv'
    call get_freq_inv( folder_strat, n_freq, freq_0_k, freq_k, iter_k )
    
    !**********************************************************!
    ! Reading model and performing some initial QC.
    if ( rank == 0 ) write(*,*) 'get_model_and_QC'
    call get_model_and_QC( cmpr, buoy, t1_inv_aux, n_freq, freq_k )
    
    ! From dnst to buoy.
    buoy = 1./buoy
    
    !**********************************************************!
    ! FWI method.
    if ( rank == 0 ) write(*,'(a,i3)') ' FWI strat_inv = ', strat_inv
    call fwi( n_sg, shgat, n_freq, freq_0_k, freq_k, iter_k )
    
    !**********************************************************!
    ! Deallocate.
    deallocate( freq_0_k, freq_k, iter_k, stat=stal ); if ( stal/=0 ) stop 'AE invA2b'
    deallocate( cmpr, buoy, stat=stal ); if ( stal/=0 ) stop 'dAE inv 3'
    call coll_deallo
    call rmv_sg( shgat, n_sg(rank+1) )
    deallocate( n_sg, stat=stal ); if ( stal/=0 ) stop 'dAE inv 2'
    
    !**********************************************************!
    ! Ending.
    if ( rank == 0 ) write(*,*) 'We have finiiiiiiiished :D !!!!!!!!!!'
    
    ! Shut down MPI.
#ifdef usempi
    call MPI_FINALIZE(ierr)
#endif
    
  end program main




