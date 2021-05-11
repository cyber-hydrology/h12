!     ******************************************************************
!                          2D INUNDATION CALCULATION
!     ******************************************************************
!     ====================================
!              MOMENTUM EQUATION
!     ====================================
      subroutine flux
      use globals
!      use omp_lib
!      use mpi
     
      implicit none
      
      real*8:: hhe, hhw, hhep, hhwp, hhan, sgnm, hh3
      real*8:: u13, u11uur, u11uul, umr, uml, u11, u12vvu, u12vvd, umu, umd, u12, sqx, ram
      real*8:: hhn, hhs, hhnp, hhsp
      real*8:: v13, v11uur, v11uul, vnr, vnl, v11, v12vvu, v12vvd, vnu, vnd, v12, sqy
      real*8:: ddx, ddy
      real*8:: h_lc, h_dc, baseo_lc, baseo_dc
      real*8:: uu1_rc, uu1_lc, uu1_rdc, uu1_dc
      real*8:: vv1_uc, vv1_luc, vv1_lc, vv1_dc
      real*8:: umo_rc, umo_lc, umo_uc, umo_dc
      real*8:: vno_rc, vno_lc, vno_uc, vno_dc
      real*8:: val_temp
!      integer:: nci, lc, rc, uc, dc, ruc, rdc, ldc, luc,lev, levn, ndir
      integer:: modex, modey, nm, int_temp
!     ------------------------------------
!             MPI local variables
!     ------------------------------------
      integer:: i,j,tag(1:100)=(/(j,j=1,100)/), status(MPI_STATUS_SIZE),tempstart,tempend
      integer:: send_request(100),recv_request(100)
      real*8:: tstart,tfinish

    tempstart=istartarray(rank)
    tempend=iendarray(rank)
!     ------------------------------------
!             UM, VN calculation
!     ------------------------------------
!$omp barrier
!$omp parallel do private(hhe,hhw,hhep,hhwp,hhan,sgnm,hh3,u13,u11uur,u11uul,umr,uml,u11,u12vvu, &
u12vvd,umu,umd,u12,sqx,ram,hhn,hhs,hhnp,hhsp,v13,v11uur,v11uul,vnr,vnl,v11,v12vvu,v12vvd,vnu, &
vnd,v12,sqy,ddx,ddy,h_lc,h_dc,baseo_lc,baseo_dc,uu1_rc,uu1_lc,uu1_rdc,uu1_dc,vv1_uc,vv1_luc, &
vv1_lc,vv1_dc,umo_rc,umo_lc,umo_uc,umo_dc,vno_rc,vno_lc,vno_uc,vno_dc,nci,lc,rc,uc,dc,ruc,rdc,ldc,luc, &
lev,levn,modex,modey,nm,val_temp, int_temp,rn) shared(h,baseo,um,vn,nei_info,umo,vno,vv1,uu1, &
inf,um_fine, vn_fine,c_lev_id, c_land)
!    do nci=1, ncell
      do nci=tempstart, tempend

        if(c_land(nci) ==10) rn = 1.0d0 ! water
        if(c_land(nci) ==20) rn = 1.0d0 ! impervious Surface
        if(c_land(nci) ==30) rn = 0.5d0 ! Soil & Barren
        if(c_land(nci) ==40) rn = 0.5d0 ! Trees & Forest
        if(c_land(nci) ==52) rn = 0.5d0 ! Shrubs
        if(c_land(nci) ==70) rn = 0.5d0 ! Grass & Herbaceous
        if(c_land(nci) ==80) rn = 0.7d0 ! Agriculture
        if(c_land(nci) ==91) rn = 1.0d0 ! Woody Wetlands
        if(c_land(nci) ==92) rn = 1.0d0 ! Emergent Wetlands

!     ------------------------------------
!                  X-DIRECTION
!     ------------------------------------
! x-, y- direction distance calculation

      lev =c_lev_id(nci)
      levn=c_lev_id(nei_info(nci,4))
      ddx = (dx/dble(2**lev)+dx/dble(2**levn))/2.0d0
      ddy = (dy/dble(2**lev)+dy/dble(2**c_lev_id(nei_info(nci,4))))/2.0d0
      nm=0
          lc = nei_info(nci,4); rc = nei_info(nci,2)
          uc = nei_info(nci,1); dc = nei_info(nci,3)
          ruc = nei_info(nci,5); rdc = nei_info(nci,6)
          ldc = nei_info(nci,7); luc = nei_info(nci,8)
!     ---------------------------------------------
!         find averaged neighbor cell information
!     ---------------------------------------------
!           
         if(uc==0)then
           vv1_uc=vv1_ugc(nci); umo_uc=0.0d0; vno_uc=0.0d0
         else
           vv1_uc=vv1(uc)
           umo_uc=umo(uc)
           vno_uc=vno(uc)
         endif
!
         if(rc==0)then
             uu1_rc=uu1_rgc(nci); umo_rc=0.0d0; vno_rc=0.0d0
         else
             uu1_rc=uu1(rc)
             umo_rc=umo(rc)
             vno_rc=vno(rc)
         endif
!
         if(dc==0)then
            h_dc=0.0d0; baseo_dc=0.0d0;	uu1_dc=0.0d0
            vv1_dc=0.0d0; umo_dc=0.0d0; vno_dc=0.0d0
         else        
            h_dc=h(dc)
            baseo_dc=baseo(dc)
            uu1_dc=uu1(dc)
            vv1_dc=vv1(dc)
            umo_dc=umo(dc)
            vno_dc=vno(dc)
         endif
!
          if(rdc==0)then
            uu1_rdc=0.0d0
          else        
            uu1_rdc=uu1(rdc)
          endif

          if(luc==0)then
            vv1_luc=0.0d0
          else        
            vv1_luc=vv1(luc)
          endif
!
1 continue
      nm=nm+1
!
         if(lc==0)then
            h_lc=0.0d0; baseo_lc=0.0d0
            uu1_lc=0.0d0; vv1_lc=0.0d0
            umo_lc=0.0d0; vno_lc=0.0d0
         else        
            h_lc=h(lc)
            baseo_lc=baseo(lc)
            uu1_lc=uu1(lc)
            vv1_lc=vv1(lc)
            umo_lc=umo(lc)
            vno_lc=vno(lc)
         endif      
!     ----------------------
!		 x-direction
!     ---------------------- 
           if(lc==0) goto 110
    	    if(inf(nci) == 0) goto 110
        	if(inf(lc) == 0) goto 110
        	if(h(lc) <= th .and. h(nci) <= th) goto 110
        	if(infx(nci) == 1) goto 110
             hhe = h(nci) + baseo(nci)
             hhw = h(lc) + baseo(lc)
             hhep = h(nci) - th
             hhwp = h(lc) - th
!     ----------------------
!       DRY BED TREATEMTN (1)
!     ----------------------
        	if(hhe < baseo(lc)) then
          	  if(h_lc > th) then
            	um(nci) = 0.544d0*h_lc*sqrt(gg*h_lc)
              else
            	um(nci) = 0.0d0
              endif
              modex =1
            	goto 120
        	elseif(hhw < baseo(nci)) then
          	  if(h(nci) > th) then
                um(nci) = -0.544d0*h(nci)*sqrt(gg*h(nci))
          	  else
            	um(nci) = 0.0d0
            endif
            	goto 120
!     ----------------------
!       DRY BED TREATEMTN (2)
!     ----------------------
        	elseif(hhep*hhwp < 0.0d0) then
          	  if(h(nci) > 0.0d0 .or. h_lc > 0.0d0) then
            	hhan = hhep - hhwp
            	sgnm = hhan/abs(hhan)
            	hh3 = max(h(nci)+baseo(nci), h_lc+baseo_lc) - max(baseo(nci), baseo_lc)
            	um(nci) = -sgnm*0.35d0*hh3*sqrt(2.0d0*gg*hh3)
            else
              um(nci) = 0.0d0
            endif
              goto 120
        	endif
!     ----------------------
!        GRAVITY TERM
!     ----------------------
        u13 = gg*hm(nci)*dt2/ddx*(h(nci) + baseo(nci) - h_lc - baseo_lc)
!     ----------------------
!       CONVECTION TERUM
!     ----------------------

          u11uur = 0.5d0*(uu1_rc + uu1(nci))
          u11uul = 0.5d0*(uu1(nci) + uu1_lc)
          umr = u11uur*(umo_rc + umo(nci))*0.5d0 + abs(u11uur)*(umo(nci) - umo_rc)*0.5d0
          uml = u11uul*(umo(nci) + umo_lc)*0.5d0 + abs(u11uul)*(umo_lc - umo(nci))*0.5d0
          u11 = dt2/ddx*(umr - uml)

        if(dc == 0) then
          u12 = 0.0d0
        elseif(uc == 0) then
          u12 = 0.0d0
        else
          u12vvu = 0.5d0*(vv1_uc + vv1_luc)
          if(luc==0) u12vvu=vv1_uc          ! modified
          u12vvd = 0.5d0*(vv1(nci) + vv1_lc)
          umu = u12vvu*(umo_uc + umo(nci))*0.5d0 + abs(u12vvu)*(umo(nci) - umo_uc)*0.5d0
          umd = u12vvd*(umo_dc + umo(nci))*0.5d0 + abs(u12vvd)*(umo_dc - umo(nci))*0.5d0
          u12 = dt2/ddy*(umu - umd)
        endif 
!     ----------------------
!        FRICTION TERM
!     ----------------------
        sqx = sqrt(uu1(nci)**2 + vva(nci)**2)
        ram = gg*rn**2*sqx/hm(nci)**1.333333d0
        if(hm(nci) <= th) ram = 0.0d0
!     ----------------------
!              UM
!     ----------------------
        um(nci) = ((1.0d0 - dt2*ram*0.5d0)*umo(nci) + (-u11 - u12 - u13))/(1.0d0 + dt2*ram*0.5d0)
        modex =4
		goto 120

110	  um(nci) = 0.0d0

120   continue
!if(nci==647993) print*, time, 'flux um',um(nci)
! re-calculation for fine mesh
!    if(levn==lev+1.and.nm==1)then
!      val_temp=um(nci)
!      um_fine(nei_info(nci,4))=val_temp
!      lc=nei_info(nei_info(nci,4),1)
!      goto 1
!    elseif(levn==lev+1.and.nm==2)then
!      val_temp=um(nci)
!      um_fine(nei_info(nei_info(nci,4),1))=val_temp
!      val_temp=(um_fine(nei_info(nci,4))+um_fine(nei_info(nei_info(nci,4),1)))*0.5d0
!      um(nci)=val_temp
!      lc = nei_info(nci,4)
!    endif
!=========================================================
      levn=c_lev_id(nei_info(nci,3))
!     ---------------------------------------------
!         find averaged neighbor cell information
!     ---------------------------------------------
         if(lc==0)then
            h_lc=0.0d0; baseo_lc=0.0d0
            uu1_lc=0.0d0; vv1_lc=0.0d0
            umo_lc=0.0d0; vno_lc=0.0d0
         else        
            h_lc=h(lc)
            baseo_lc=baseo(lc)
            uu1_lc=uu1(lc)
            vv1_lc=vv1(lc)
            umo_lc=umo(lc)
            vno_lc=vno(lc)
         endif
!
        nm=0
2 continue
        nm=nm+1
!
         if(dc==0)then
            h_dc=0.0d0; baseo_dc=0.0d0;	uu1_dc=0.0d0
            vv1_dc=0.0d0; umo_dc=0.0d0; vno_dc=0.0d0
         else        
            h_dc=h(dc)
            baseo_dc=baseo(dc)
            uu1_dc=uu1(dc)
            vv1_dc=vv1(dc)
            umo_dc=umo(dc)
            vno_dc=vno(dc)
         endif
!     ------------------------------------
!                  Y-DIRECTION
!     ------------------------------------
        if(dc == 0) goto 210
        if(inf(nci) == 0) goto 210
        if(inf(dc) == 0) goto 210
        if(h_dc <= th .and. h(nci) <= th) goto 210
        if(infy(nci) == 1) goto 210
          hhn = h(nci) + baseo(nci)
          hhs = h_dc + baseo_dc
          hhnp = h(nci) - th
          hhsp = h_dc - th
!     ----------------------
!      DRY BED TREATEMTN (1)
!     ----------------------
        if(hhn < baseo_dc) then
          if(h(dc) > th) then
            vn(nci) = 0.544d0*h_dc*sqrt(gg*h_dc)
          else
            vn(nci) = 0.0d0
          endif
          modey =5
            goto 220
        elseif(hhs < baseo(nci)) then
          if(h(nci) > th) then
            vn(nci) = -0.544d0*h(nci)*sqrt(gg*h(nci))
          else
            vn(nci) = 0.0d0
          endif
            goto 220
!     ----------------------
!      DRY BED TREATEMTN (2)
!     ----------------------
        elseif(hhnp*hhsp < 0.0d0) then
          if(h(nci) > 0.0d0 .or. h_dc > 0.0d0) then
            hhan = hhnp - hhsp
            sgnm = hhan/abs(hhan)
            hh3 = max(h(nci)+baseo(nci), h_dc+baseo_dc) - max(baseo(nci), baseo_dc)
            vn(nci) = -sgnm*0.35d0*hh3*sqrt(2.0d0*gg*hh3)
          else
            vn(nci) = 0.0d0
          endif

            goto 220               
        endif
!     ----------------------
!          GRAVITY TERM
!     ----------------------
        v13 = gg*hn(nci)*dt2/ddy*(h(nci) + baseo(nci) - h_dc - baseo_dc)
!     ----------------------
!        CONVECTION TERM
!     ----------------------
        if(lc == 0) then
          v11 = 0.0d0
        elseif(rc == 0) then
          v11 = 0.0d0
        else
          v11uur = 0.5d0*(uu1_rc + uu1_rdc)
          if(rdc==0) v11uur = uu1_rc  ! modified 
          v11uul = 0.5d0*(uu1(nci) + uu1_dc)
          vnr = v11uur*(vno_rc + vno(nci))*0.5d0 + abs(v11uur)*(vno(nci) - vno_rc)*0.5d0
          vnl = v11uul*(vno(nci) + vno_lc)*0.5d0 + abs(v11uul)*(vno_lc - vno(nci))*0.5d0
          v11 = dt2/ddx*(vnr - vnl)
        endif

          v12vvu = 0.5d0*(vv1_uc + vv1(nci))
          v12vvd = 0.5d0*(vv1(nci) + vv1_dc)
          vnu = v12vvu*(vno_uc + vno(nci))*0.5d0 + abs(v12vvu)*(vno(nci) - vno_uc)*0.5d0
          vnd = v12vvd*(vno_dc + vno(nci))*0.5d0 + abs(v12vvd)*(vno_dc - vno(nci))*0.5d0
          v12 = dt2/ddy*(vnu - vnd)
!$$$$$$         endif
!     ----------------------
!        FRICTION TERM
!     ----------------------
        sqy = sqrt(uua(nci)**2 + vv1(nci)**2)
        ram = gg*rn**2*sqy/hn(nci)**1.333333d0
        if(hn(nci) <= th) ram = 0.0d0
!     ----------------------
!               VN
!     ----------------------
        vn(nci) = ((1.0d0 - dt2*ram*0.5d0)*vno(nci) + (-v11 - v12 - v13))/(1.0d0 + dt2*ram*0.5d0)
        modey =8
		goto 220
  210   vn(nci) = 0.0d0

  220	continue
!if(nci==647993) print*, time, 'flux vn',um(nci)
! re-calculation for fine mesh
!    if(levn==lev+1.and.nm==1)then
!      val_temp=vn(nci)
!      vn_fine(nei_info(nci,3))=val_temp
!      dc = nei_info(nei_info(nci,3),2)
!      goto 2
!    elseif(levn==lev+1.and.nm==2)then
!      val_temp=vn(nci)
!      vn_fine(nei_info(nei_info(nci,3),2))=val_temp
!      val_temp=(vn_fine(nei_info(nci,3))+vn_fine(nei_info(nei_info(nci,3),2)))*0.5d0
!      vn(nci)=val_temp
!    endif
!		enddo
!$omp end parallel do
!$omp barrier
!      if(rank.eq.0)tstart=omp_get_wtime()
!      if(rank.ne.npart)then
!        call MPI_ISEND(um(imoveupstart(rank):imoveupend(rank)),ighost, &
!        MPI_REAL8,rank+1,tag(2), MPI_COMM_WORLD,send_request(2),ierror)
!        call MPI_ISEND(vn(imoveupstart(rank):imoveupend(rank)),ighost, &
!        MPI_REAL8,rank+1,tag(3), MPI_COMM_WORLD,send_request(3),ierror)
        
!        call MPI_IRECV(um(imovedownstart(rank+1):imovedownend(rank+1)),ighost, &
!        MPI_REAL8,rank+1,tag(4), MPI_COMM_WORLD,recv_request(4),ierror)
!        call MPI_IRECV(vn(imovedownstart(rank+1):imovedownend(rank+1)),ighost, &
!        MPI_REAL8,rank+1,tag(5), MPI_COMM_WORLD,recv_request(5),ierror)
        
!        call MPI_WAIT(send_request(2),status,ierror)
!        call MPI_WAIT(send_request(3),status,ierror)
!        call MPI_WAIT(recv_request(4),status,ierror)
!        call MPI_WAIT(recv_request(5),status,ierror)
!        print*,'recv at rank',rank,'ighost',ighost,imovedownend(rank+1)-imovedownstart(rank+1)+1, &
!        imovedownend(rank+1),imovedownstart(rank+1)
!      endif
!      if(rank.ne.0)then
!        call MPI_IRECV(um(imoveupstart(rank-1):imoveupend(rank-1)),ighost, &
!        MPI_REAL8,rank-1,tag(2), MPI_COMM_WORLD,recv_request(2),ierror)
!        call MPI_IRECV(vn(imoveupstart(rank-1):imoveupend(rank-1)),ighost, &
!        MPI_REAL8,rank-1,tag(3), MPI_COMM_WORLD,recv_request(3),ierror)
        
!        call MPI_ISEND(um(imovedownstart(rank):imovedownend(rank)),ighost, &
!        MPI_REAL8,rank-1,tag(4), MPI_COMM_WORLD,send_request(4),ierror)
!        call MPI_ISEND(vn(imovedownstart(rank):imovedownend(rank)),ighost, &
!        MPI_REAL8,rank-1,tag(5), MPI_COMM_WORLD,send_request(5),ierror)
        
!        call MPI_WAIT(recv_request(2),status,ierror)
!        call MPI_WAIT(recv_request(3),status,ierror)
!        call MPI_WAIT(send_request(4),status,ierror)
!        call MPI_WAIT(send_request(5),status,ierror)
!        print*,'send from rank',rank,'ighost',ighost,imovedownend(rank)-imovedownstart(rank)+1, &
!        imovedownend(rank),imovedownstart(rank)
!      endif
!      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
!      if(rank.eq.0)tfinish=omp_get_wtime()
!      if(rank.eq.0)print'("mpi flux =",f12.3," <sec>")', tfinish-tstart
!end subroutine flux

