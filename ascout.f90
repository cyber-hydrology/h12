!************************************
    subroutine ascout_h(t,ci)
!************************************
      use globals
      use mpi
      implicit none

      integer::i,j,k,ci,cind,nm
      real*8:: xnode, ynode, bs, t,total_h
      character(len=6)::ft
      integer:: tag(1:100)=(/(j,j=1,100)/), status(MPI_STATUS_SIZE)  ! MPI local variables
      
    if(rank.ne.0)then
        call MPI_SEND(h(istarttransfer(rank):iendtransfer(rank)),ichunktransfer(rank), &
        MPI_REAL8,0,tag(4), MPI_COMM_WORLD,ierror)
    else
        do i=1,npart
        call MPI_RECV(h(istarttransfer(i):iendtransfer(i)),ichunktransfer(i), &
        MPI_REAL8,i,tag(4), MPI_COMM_WORLD,status,ierror)
        enddo
    endif
    call MPI_BCAST(h,ncell,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
    
      if(rank.eq.0)then
      cind=0
      do i=1, im
        do j=1, jm
          if(bs_dem(i,j)<-9998) then
            h_grid(i,j)=0.0d0
          else
              cind=cind+1
              h_grid(i,j) = h(cind)
          endif
       enddo
      enddo
      endif
      
      if(rank.eq.0)print*,"ending h to h_grid"
     write(ft,'(i6)') int(t)+100000
      if(rank.eq.0)print*,"ending ft"
      if(rank.eq.0)then
      open(200, file='./results/ascout/ascout_h_'//ft//'.asc', action='write') 
      print*,"ending open asc file"
      write(200,'(a5,1x,i10)')'ncols',im    !i=x-direction
      write(200,'(a5,1x,i10)')'nrows',jm    !j=y-dirction
      write(200,'(a9,1x,f11.3)')'xllcorner',x_ori
      write(200,'(a9,1x,f11.3)')'yllcorner',y_ori
      write(200,'(a8,1x,f11.3)')'cellsize',dx
      write(200,'(a23)')'NODATA_value  -9999.000'
	  
	  do i=1, im
	    do j=1, jm
		 if(h_grid(i,j)<=0.0d0) h_grid(i,j)=-9999
		enddo
      enddo
	  
      do j=jm, 1, -1
        write(200,'(8614(f9.3,1x))') (h_grid(i,j), i=1, im)
      enddo

      close(200)
      print*,"ending writing asc file"
      endif
      
      end subroutine ascout_h

!************************************
    subroutine ascout_u(t,ci)
!************************************
      use globals
      use mpi
      implicit none

      integer::i,j,k,ci,cind,nm
      real*8:: xnode, ynode, bs, t,total_h
      character(len=6)::ft
      integer:: tag(1:100)=(/(j,j=1,100)/), status(MPI_STATUS_SIZE)  ! MPI local variables
      
    if(rank.ne.0)then
        call MPI_SEND(um(istarttransfer(rank):iendtransfer(rank)),ichunktransfer(rank), &
        MPI_REAL8,0,tag(4), MPI_COMM_WORLD,ierror)
    else
        do i=1,npart
        call MPI_RECV(um(istarttransfer(i):iendtransfer(i)),ichunktransfer(i), &
        MPI_REAL8,i,tag(4), MPI_COMM_WORLD,status,ierror)
        enddo
    endif
    call MPI_BCAST(um,ncell,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
    
      if(rank.eq.0)then
      cind=0
      do i=1, im
        do j=1, jm
          if(bs_dem(i,j)<-9998) then
            h_grid(i,j)=0.0d0
			um_grid(i,j)=0.0d0
          else
              cind=cind+1
              h_grid(i,j) = h(cind)
			  um_grid(i,j) = um(cind)
          endif
       enddo
      enddo
      endif
      
      if(rank.eq.0)print*,"ending u to u_grid"
      write(ft,'(i6)') int(t)+100000
      if(rank.eq.0)print*,"ending ft"
      if(rank.eq.0)then
      open(200, file='./results/ascout/ascout_u_'//ft//'.asc', action='write') 
      print*,"ending open asc file"
      write(200,'(a5,1x,i10)')'ncols',im    !i=x-direction
      write(200,'(a5,1x,i10)')'nrows',jm    !j=y-dirction
      write(200,'(a9,1x,f11.3)')'xllcorner',x_ori
      write(200,'(a9,1x,f11.3)')'yllcorner',y_ori
      write(200,'(a8,1x,f11.3)')'cellsize',dx
      write(200,'(a19)')'NODATA_value  -9999.0000'
	  do i=1,im
	    do j=1,jm
		  if(h_grid(i,j)<=0.0d0) then
		    um_grid(i,j)=-9999.0000d0
		  else
		    um_grid(i,j)=um_grid(i,j)/h_grid(i,j)
	          endif
		enddo
      enddo
      do j=jm, 1, -1
        write(200,'(8614(f10.4,1x))') (um_grid(i,j), i=1, im)
      enddo

      close(200)
      print*,"ending writing asc file"
      endif
      
      end subroutine ascout_u
	  
!************************************
    subroutine ascout_v(t,ci)
!************************************
      use globals
      use mpi
      implicit none

      integer::i,j,k,ci,cind,nm
      real*8:: xnode, ynode, bs, t,total_h
      character(len=6)::ft
      integer:: tag(1:100)=(/(j,j=1,100)/), status(MPI_STATUS_SIZE)  ! MPI local variables
      
    if(rank.ne.0)then
        call MPI_SEND(vn(istarttransfer(rank):iendtransfer(rank)),ichunktransfer(rank), &
        MPI_REAL8,0,tag(4), MPI_COMM_WORLD,ierror)
    else
        do i=1,npart
        call MPI_RECV(vn(istarttransfer(i):iendtransfer(i)),ichunktransfer(i), &
        MPI_REAL8,i,tag(4), MPI_COMM_WORLD,status,ierror)
        enddo
    endif
    call MPI_BCAST(vn,ncell,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
    
      if(rank.eq.0)then
      cind=0
      do i=1, im
        do j=1, jm
          if(bs_dem(i,j)<-9998) then
            h_grid(i,j)=0.0d0
			vn_grid(i,j)=0.0d0
          else
              cind=cind+1
              h_grid(i,j) = h(cind)
			  vn_grid(i,j) = vn(cind)
          endif
       enddo
      enddo
      endif
      
      if(rank.eq.0)print*,"ending v to v_grid"
      write(ft,'(i6)') int(t)+100000
      if(rank.eq.0)print*,"ending ft"
      if(rank.eq.0)then
      open(200, file='./results/ascout/ascout_v_'//ft//'.asc', action='write') 
      print*,"ending open asc file"
      write(200,'(a5,1x,i10)')'ncols',im    !i=x-direction
      write(200,'(a5,1x,i10)')'nrows',jm    !j=y-dirction
      write(200,'(a9,1x,f11.3)')'xllcorner',x_ori
      write(200,'(a9,1x,f11.3)')'yllcorner',y_ori
      write(200,'(a8,1x,f11.3)')'cellsize',dx
      write(200,'(a19)')'NODATA_value  -9999.0000'
	  do i=1,im
	    do j=1,jm
		  if(h_grid(i,j)<=0.0d0) then
		    vn_grid(i,j)=-9999.0000d0
		  else
		    vn_grid(i,j)=vn_grid(i,j)/h_grid(i,j)
	      endif
		enddo
      enddo
      do j=jm, 1, -1 
        write(200,'(8614(f10.4,1x))') (vn_grid(i,j), i=1, im)
      enddo

      close(200)
      print*,"ending writing asc file"
      endif
      
      end subroutine ascout_v
