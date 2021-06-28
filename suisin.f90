!     ====================================
!             CONTINUITY EQUATION
!     ====================================
	  subroutine suisin
      use globals
!      use omp_lib
!      use mpi
      implicit none

      integer:: nci, lc, rc, uc, dc, ruc, rdc, ldc, luc, it
      integer:: i, k, lev, levn, ndir, ir, jr, irh, j
      real*8:: ddx, ddy, xc, yc, um_rc, vn_uc,rfcd,rr, level
      real*8:: x_dist, y_dist, x_rest, y_rest
!      integer:: j,tag(1:100)=(/(j,j=1,100)/), status(MPI_STATUS_SIZE),tempstart,tempend
!      integer:: recv_request(100),send_request(100)
!      real*8:: htemp(mnhl*npart),qlmetemp(mnhl*npart)
!      integer:: send_request_h(cnode),recv_request_h(cnode)
!      integer:: send_request_qlme(cnode),recv_request_qlme(cnode)
!      integer:: sendingrank,receiverank,tag_h(1:500)=(/(j,j=1,500)/),tag_qlme(1:500)=(/(j,j=1,500)/),jstart,jend

!        if(rank.eq.0)then
!            do j=1,npart
!                receiverank = j
!                jstart=((j-1)*mnhl+1)
!                jend=((j-1)*mnhl+ic_no(j))
!                do i=1,ic_no(j)
!                  qlmetemp(jstart+i-1)=qlme(ic_location(j,i))
!                enddo
!                call MPI_ISEND(qlmetemp(jstart:jend),ic_no(j), &
!                MPI_REAL8,receiverank,tag_qlme(j), MPI_COMM_WORLD,send_request_qlme(j),ierror)   
!                call MPI_WAIT(send_request_qlme(j),status,ierror)
!            enddo
!        else
!            call MPI_IRECV(qlmetemp(1:ic_no(rank)),ic_no(rank), &
!            MPI_REAL8,0,tag_qlme(rank), MPI_COMM_WORLD,recv_request_qlme(rank),ierror)
!            call MPI_WAIT(recv_request_qlme(rank),status,ierror)
!                do i=1,ic_no(rank)
!                  qlme(ic_location(rank,i))=qlmetemp(i)
!                enddo
!        endif 

!      call MPI_BCAST(qlme,ncell,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
!      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      
!      tempstart=istartarray(rank)
!      tempend=iendarray(rank)

      
      err_ini=0.0d0
      
!!$omp barrier
!!$omp parallel do private(nci,lc,rc,uc,dc,ruc,rdc,ldc,luc,lev,levn,ndir, &
!ddx,ddy,xc,yc,um_rc,vn_uc,rfcd,x_dist,y_dist,x_rest,y_rest,ir,jrm,jr,rr,it,level,irh)&
!shared(h,ho,um,vn,nei_info,c_lev_id,um_fine,vn_fine,c_land,qlme,x_rad_ori,&
!y_rad_ori,rad_cell,baseo,dtlevel,dwlevel,time)
      do nci=1, ncell
          lc = nei_info(nci,4); rc = nei_info(nci,2)
          uc = nei_info(nci,1); dc = nei_info(nci,3)
          ruc = nei_info(nci,5); rdc = nei_info(nci,6)
          ldc = nei_info(nci,7); luc = nei_info(nci,8)
          
        if(lc == 0)then
          h(nci)=0.0d0
          goto 100
        elseif(rc == 0)then
          h(nci)=0.0d0
          goto 100
        elseif(uc==0)then
          h(nci)=0.0d0
          goto 100
        elseif(dc==0)then
          h(nci)=0.0d0
          goto 100
        elseif(ruc == 0)then
          h(nci)=0.0d0
          goto 100
        elseif(rdc == 0)then
          h(nci)=0.0d0
          goto 100
        elseif(ldc == 0)then
          h(nci)=0.0d0
          goto 100
        elseif(luc == 0)then
          h(nci)=0.0d0
          goto 100
        endif

        if(c_land(nci) ==10) rfcd = 1.0d0 ! water
        if(c_land(nci) ==20) rfcd = 1.0d0 ! impervious Surface
        if(c_land(nci) ==30) rfcd = 0.5d0 ! Soil & Barren
        if(c_land(nci) ==40) rfcd = 0.5d0 ! Trees & Forest
        if(c_land(nci) ==52) rfcd = 0.5d0 ! Shrubs
        if(c_land(nci) ==70) rfcd = 0.5d0 ! Grass & Herbaceous
        if(c_land(nci) ==80) rfcd = 0.7d0 ! Agriculture
        if(c_land(nci) ==91) rfcd = 1.0d0 ! Woody Wetlands
        if(c_land(nci) ==92) rfcd = 1.0d0 ! Emergent Wetlands

! =====================================================
!             x, y coordination calculation
! =====================================================
!        xc=ccoords(nci,1); yc=ccoords(nci,2)
        
!         x_dist=xc-x_rad_ori
!         y_dist=yc-y_rad_ori
 
!         x_rest=mod(x_dist,rad_cell)
!         y_rest=mod(y_dist,rad_cell)
 
!         ir=irm-int(y_dist-y_rest)/int(rad_cell)
!         jr=int(x_dist-x_rest)/int(rad_cell)+1
 
!         rr = rain(ir,jr)/1000.0d0/3600.d0*rfcd ! mm/hour -> m/sec ! it must be changed regarding radar data format.
!         if(rr<=0.0d0) rr=0.0d0

!          rr = 100.0d0 ! rainfall intencity
!          rr = rr/1000.0d0/3600.0d0*rfcd
           irh = int(time/dtrain) +1
           rr = rain_p(irh)/1000.0d0/3600.0d0*rfcd
!======================================================
!
         lev =c_lev_id(nci)
         levn=c_lev_id(nei_info(nci,1))
      if(lev-1==levn)then
         if(uc==0)then
           vn_uc=0.0d0
         else
           vn_uc=vn_fine(nci)
         endif
      else
         if(uc==0)then
           vn_uc=0.0d0
         else
           vn_uc=vn(nei_info(nci,1))
         endif
      endif
         ddy = dy/dble(2**lev)
         levn=c_lev_id(nei_info(nci,2))
      if(lev-1==levn)then
         if(rc==0)then
           um_rc=0.0d0
         else        
           um_rc=um_fine(nci)
         endif
      else
         if(rc==0)then
           um_rc=0.0d0
         else        
           um_rc=um(nei_info(nci,2))
         endif
      endif
           ddx = dx/dble(2**lev)

        h(nci) = ho(nci)-dt2*((um_rc-um(nci))/ddx+(vn_uc-vn(nci))/ddy-qlme(nci)-rr)
        h(nci) = max(h(nci),0.0d0)
!       if(time>100)print*,nci, h(nci)
!if(nci==647993) print*, time, h(nci-1),h(nci), h(nci+1)
!--------down stream boundary condition----------------------------
      if(time<=dtlevel)then
        level=dwlevel(1)
      else
        it=nint(time/dtlevel)
        level=dwlevel(it)+(dwlevel(it+1)-dwlevel(it))*(time-dtlevel*(dble(it)))/dtlevel
      endif
      if(rc==0.and.c_land(nci)==5)then
        h(nci)=max(level-baseo(nci),0.0d0)
      endif
100   enddo
!

!!$omp end parallel do
!!$omp barrier
      
!      if(rank.eq.0)tstart=omp_get_wtime()
!    if(rank.ne.npart)then
!        call MPI_ISEND(h(imoveupstart(rank):imoveupend(rank)),ighost, &
!        MPI_REAL8,rank+1,tag(2), MPI_COMM_WORLD,send_request(2),ierror)
        
!        call MPI_IRECV(h(imovedownstart(rank+1):imovedownend(rank+1)),ighost, &
!        MPI_REAL8,rank+1,tag(4), MPI_COMM_WORLD,recv_request(4),ierror)
        
!        call MPI_WAIT(send_request(2),status,ierror)
!        call MPI_WAIT(recv_request(4),status,ierror)
!        print*,'mpi_send from rank',rank,'ighost',ighost,imoveupend(rank)-imoveupstart(rank)+1, &
!        imoveupend(rank),imoveupstart(rank)
!    endif
!    if(rank.ne.0)then
!        call MPI_IRECV(h(imoveupstart(rank-1):imoveupend(rank-1)),ighost, &
!        MPI_REAL8,rank-1,tag(2), MPI_COMM_WORLD,recv_request(2),ierror)

!        call MPI_ISEND(h(imovedownstart(rank):imovedownend(rank)),ighost, &
!        MPI_REAL8,rank-1,tag(4), MPI_COMM_WORLD,send_request(4),ierror)
        
!        call MPI_WAIT(recv_request(2),status,ierror)
!        call MPI_WAIT(send_request(4),status,ierror)
!        print*,'mpi_recv at rank',rank,'ghost',ighost,imoveupend(rank-1)-imoveupstart(rank-1)+1, &
!        imoveupend(rank-1),imoveupstart(rank-1)
!    endif
!    call MPI_BCAST(h,ncell,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
!      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
!      if(rank.eq.0)tfinish=omp_get_wtime()
!      if(rank.eq.0)print'("mpi suisin =",f12.3," <sec>")', tfinish-tstart


!    if(rank.ne.0)then
!        call MPI_ISEND(h(istarttransfer(rank):iendtransfer(rank)),ichunktransfer(rank), &
!                       MPI_REAL8,0,tag_h(rank), MPI_COMM_WORLD,send_request_h(rank),ierror)        
!        call MPI_WAIT(send_request_h(rank),status,ierror)

!    else
!        do j=1,npart
!           sendingrank = j
!           call MPI_IRECV(h(istarttransfer(j):iendtransfer(j)),ichunktransfer(j), &
!                          MPI_REAL8,sendingrank,tag_h(j), MPI_COMM_WORLD,recv_request_h(j),ierror)
!           call MPI_WAIT(recv_request_h(j),status,ierror)
!        enddo
!    endif
!    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
!      if(rank.ne.0)then
!        do i=1,ic_no(rank)
!          htemp(i)=h(ic_location(rank,i))
!        enddo
!        call MPI_ISEND(htemp(1:ic_no(rank)),ic_no(rank), &
!                       MPI_REAL8,0,tag_h(rank), MPI_COMM_WORLD,send_request_h(rank),ierror)        
!        call MPI_WAIT(send_request_h(rank),status,ierror)
!      else
!        do j=1,npart
!           sendingrank = j
!           jstart=((j-1)*mnhl+1)
!           jend=((j-1)*mnhl+ic_no(j))
!           call MPI_IRECV(htemp(jstart:jend),ic_no(j), &
!                          MPI_REAL8,sendingrank,tag_h(j), MPI_COMM_WORLD,recv_request_h(j),ierror)
!           call MPI_WAIT(recv_request_h(j),status,ierror)
!           do k=1,ic_no(j)
!              h(ic_location(j,k))=htemp(jstart+k-1)
!           enddo
!        enddo
!      endif
!      call MPI_BARRIER(MPI_COMM_WORLD,ierror) 

 end subroutine suisin
      
      
