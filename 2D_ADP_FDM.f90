!
!
!
!     #######################################################
!     ##                                                   ##
!     ##               main program                        ##
!     ##                                                   ##
!     #######################################################
!
!
!
    program main

    use globals
!    use omp_lib
!    use mpi 
    implicit none
  
    real*8:: time0, timmax
    real*8:: start, finish,tstart,tfinish
    integer:: i, j, lpout, lkout, ms, ci, mstep
    integer:: nci, rc, uc, nthreads, nid,tempstart,tempend
    logical(4):: resultqq,systemqq
!   ===============================================
!              Initialization of MPI process
!   =====================================================
!print*,'start'
!      call MPI_INIT(ierror)
!      call MPI_COMM_SIZE(MPI_COMM_WORLD, nsize, ierror)
!      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
!      npart = nsize - 1
!      !$OMP PARALLEL private(nthreads)
!      nthreads = OMP_GET_NUM_THREADS()
!      nid = OMP_GET_THREAD_NUM()
!      if(rank.eq.0)then
!        if(nid.eq.1)write(*,*)'no. MPI: ', nsize,' no. omp: ',nthreads      
!      endif
!      !$OMP END PARALLEL

      
  
!   ===============================================
!              NUMBER OF SIMULATION CASES
!   =====================================================
	do ci=1, 1
!     start = MPI_Wtime()

!     if(rank.eq.0)print'("start cpu time=",f12.3,"<sec>")', start
   
!   =====================================================	
!                      READING DATA
!     =====================================================
!      if(rank.eq.0)print*,"start reading data"
      call rdat(time0, timmax, lpout, lkout,  ci)
!      if(rank.eq.0)print*,"finish reading data"
      

!     =====================================================
!					  Grid Generation
!     =====================================================
		ms = 2**div_max	  
		
		allocate(c_lev(im,jm), ind(im,jm,ms,ms))
		allocate(nei_info(maxcell,8), c_lev_id(0:maxcell))
		allocate(ccoords(maxcell,2))
        allocate(c_land(maxcell))

        ncell=0 ! number of cells

        write(*,*) 'call mesh_gen'
        call mesh_gen  ! mesh generation
        write(*,*)'end nei_info'
        
!       tempstart=istartarray(rank)
!       tempend=iendarray(rank)  
!     =====================================================        
!                    INITIALIZATION
!     =====================================================
      call initiald
      write(*,*)'end initiald'
!
!      call ascout_h(time, ci)
!      call ascout_u
!      call ascout_v
!
!      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
!     =====================================================
      time = time0
      mstep = 0
!     =====================================================
!       print*,'call radar'
!       call radar
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++
!                      LOOP START
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++
!     $$$$$$$$$$$$$$$$$$$$
!         FLUX
!     $$$$$$$$$$$$$$$$$$$$
1     continue
      call flux
!     ====================
!     ERROR TREATMENT
!     ====================
!!$omp barrier
!!$omp parallel do private(rc,uc)
        do nci=1, ncell
!      do nci=tempstart, tempend
          rc = nei_info(nci,2)
          uc = nei_info(nci,1)
          if(rc==0) cycle
          if(uc==0) cycle
          if(inf(nci) /= 1 .and. inf(nci) /= 2) cycle
          if(h(nci) >= th) cycle
          if(um(rc) > 0.0d0) um(rc) = 0.0d0
          if(um(nci) < 0.0d0)   um(nci) = 0.0d0
          if(vn(uc) > 0.0d0) vn(uc) = 0.0d0
          if(vn(nci) < 0.0d0)   vn(nci) = 0.0d0
        enddo
!     ====================
!         time step       
!     ====================
       time = time + dt
       mstep = mstep + 1
!
!     $$$$$$$$$$$$$$$$$$$$
!         WATER DEPTH    
!     $$$$$$$$$$$$$$$$$$$$
      call suisin
!      if(rank.eq.0) call suisinsw
!     =================================================================
!              ×µáÜªÎß©õó
!     =================================================================
      call mdvel
!     =================================================================
!                 N+2  -> N
!     =================================================================
      call forward
!     ====================
!         time step       
!     ====================
      time = time + dt
      mstep = mstep + 1
!     ====================
!          Rain data
!     ====================
!       if(mod(mstep,ltrain) == 0)then
!           deallocate(rain)
!           call radar
!       endif
!     =================================================================
!                           writing data
!     =================================================================
!      if(mod(mstep,lpout)==0)then
!        if(rank.eq.0)write(*,'(f14.4)') time
!          finish=MPI_Wtime()
!        if(rank.eq.0)print'("running cpu time=",f12.3," <sec>")', finish-start
!      endif

       if(mod(mstep,lkout) == 0)then
 
          call ascout_h(time,ci)
          ! call ascout_u(time,ci)
          ! call ascout_v(time,ci)

       endif

 
!     =================================================================
!              time judging
!     =================================================================
      
      if(time+dt <= timmax) goto 1

!     +++++++++++++++++++++++++++++++++++++++++++++++++++++
!                          loop end
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++
     
     ! if(rank.eq.0)close(111)

     !  if(rank.eq.0) then
     !   resultqq=systemqq('cp -f /tmp/ascout_h* ./results/ascout/')
     !   resultqq=systemqq('rm -f /tmp/shm/ascout_h*')
     !   write(*,*)'systemqq resultqq:',resultqq
     ! endif
!      call MPI_BARRIER(MPI_COMM_WORLD,ierror) 
      enddo

!	  call MPI_FINALIZE(ierror) 
	end program main
