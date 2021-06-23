        subroutine store_id_info

	use globals
!    use mpi    
	implicit none
	
	integer::i,j,k,is,js,ms,lev,in,jn,isn,jsn,cind,sldbound,ndir,ix,jy,nm,ic,irank,temp_no
	real*8:: x,y,dxl,dyl,xc,yc

	c_lev_id(0)=0
	
	cind=0
	do i=1, im
	  x=x_ori + (dble(i)-0.5d0)*dx
		do j=1, jm
        if(bs_dem(i,j)<-9998)cycle
		  y=y_ori + (dble(j)-0.5d0)*dy
		  lev=c_lev(i,j); ms=2**lev
		  dxl=dx/dble(2**lev); dyl=dy/dble(2**lev)
			do is=1, ms
			  xc=x+dxl/2.0d0+dble(is-1)*dxl
				do js=1, ms
				  yc=y+dyl/2.0d0+dble(js-1)*dyl
				  cind=cind+1
				  ind(i,j,is,js)=cind
				  c_lev_id(cind)=lev
				  ccoords(cind,1)=xc; ccoords(cind,2)=yc
                  c_land(cind)=mark_land(i,j)
				  ncell=ncell+1
				enddo
			enddo
		enddo
	enddo
      write(*,*)'ncell=', ncell
	
	!** store neighbour information (cell id)
	do i=1, im
		do j=1, jm
          if(bs_dem(i,j)<-9998)cycle
		  lev=c_lev(i,j); ms=2**lev
			do is=1, ms
				do js=1, ms
				  cind=ind(i,j,is,js)
					do ndir=1,8
					  call nei_find(ndir,i,j,is,js,in,jn,isn,jsn,sldbound)
						if(isn==0.or.jsn==0)then
						  nei_info(cind,ndir)=0
						else
						  nei_info(cind,ndir)=ind(in,jn,isn,jsn)
						endif
					enddo
				enddo
			enddo
		enddo
	enddo
      deallocate(ind)
      if(rank.eq.0) write(*,*)'end store neighborhood cell info'
!   -------------------------------------
!         Allocate variables
!   -------------------------------------
	  i=ncell
	  allocate(um(i),umo(i),vn(i),vno(i),uu(i),vv(i),uua(i),vva(i))
      allocate(um_fine(i), vn_fine(i))
      allocate(uu1(i),vv1(i),h(i),ho(i),baseo(i),mark(i))
      allocate(inf(i), infx(i), infy(i))
      allocate(hm(i), hn(i))
      allocate(cn_cell(i),cn_mh(i),qlme(i))
      allocate(cei(i))
 	 boundary data allocation
	  allocate(um_rgc(i), hm_rgc(i), vn_ugc(i), hn_ugc(i))
      allocate(uu1_rgc(i), vv1_ugc(i), uua_ugc(i), vva_rgc(i))
      if(rank.eq.0) write(*,*)'end allocation of variables'
!
!   -------------------------------------
!	     Calculation of baseo
!   -------------------------------------
	cind=0
    cnode=0
      do i=1, im
	  x=x_ori + (dble(i)-1.0d0)*dx
        do j=1, jm
!print*, i, j, bs_dem(i,j)
        if(bs_dem(i,j)<-9998)cycle
		  y=y_ori + (dble(j)-1.0d0)*dy
		  lev=c_lev(i,j); ms=2**lev
		  dxl=dx/dble(2**lev); dyl=dy/dble(2**lev)
            do is = 1, ms
              xc=x+dxl/2.0d0+dble(is-1)*dxl
              do js = 1, ms
                yc=y+dyl/2.0d0+dble(js-1)*dyl
                cind=cind+1
                ix=ceiling((xc-x_ori)/dxori)
                jy=ceiling((yc-y_ori)/dyori)
                nm=0
                baseo(cind)=0.0d0
                mark(cind)=0
                  if(mark_land(i,j)==4)then
                    do in=ix-int(dxl/dxori/2.0d0),ix+int(dxl/dxori/2.0d0)
                      do jn=jy-int(dyl/dyori/2.0d0),jy+int(dyl/dyori/2.0d0)
                        if(bs_ori_dem(in,jn)>-9998)then
                          baseo(cind) =baseo(cind)+ bs_ori_dem(in,jn)
                          nm=nm+1
                        endif
                      enddo
                    enddo
                    mark(cind) = 4
                    if(nm/=0) baseo(cind)=baseo(cind)/dble(nm)
                  endif
                  if(nm==0) baseo(cind) = bs_dem(i,j)
!print*, cind, baseo(cind)
                  cei(cind)=0
            enddo
          enddo
        enddo
      enddo	

!   -------------------------------------
!	     Calculation of baseo, find initial water depth position
!   -------------------------------------
	cind=0
    nm=0
      do i=1, im
	  x=x_ori + (dble(i)-1.0d0)*dx
        do j=1, jm
        if(bs_dem(i,j)<-9998)cycle
		  y=y_ori + (dble(j)-1.0d0)*dy
		  lev=c_lev(i,j); ms=2**lev
		  dxl=dx/dble(2**lev); dyl=dy/dble(2**lev)
            do is = 1, ms
              xc=x+dxl/2.0d0+dble(is-1)*dxl
              do js = 1, ms
                yc=y+dyl/2.0d0+dble(js-1)*dyl
                cind=cind+1
                cei(cind)=-999
                do k=1, nm_ini
                  if(i_ini(k)==i.and.j_ini(k)==j)then
                    nm=nm+1
                    cei(cind)=ele_ini(nm)
!                    print*, i, j, nm, cind, ele_ini(nm)
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
          
      if(rank.eq.0)print*,'end store_id_info'

!     =====================================================
        ichunk=int(ncell/nsize)
!        ighost=int(ichunk*0.05)
!        ighost=int(20000/(div_max*dx))  !default for 1m resolution at ver0.98
!        ighost=int(5000/1) !test value
!         ighost=int(jm/(dx/2**div_max))
!         ighost = int(jm*5)
!         ighost = int(jm*4)
         ighost = int(jm*2)

print*, 'jm=',jm,'dx=',dx, 'div_max=', div_max, 'ighost=',ighost

        allocate(ichunkarray(0:npart),istartarray(0:npart),iendarray(0:npart)) 
        allocate(istarttransfer(0:npart),iendtransfer(0:npart),ichunktransfer(0:npart))
        allocate(imoveupstart(0:npart),imoveupend(0:npart))
        allocate(imovedownstart(0:npart),imovedownend(0:npart))
        if(rank.eq.0)then
            ichunkarray(0)=ichunk+ighost
            ichunkarray(1:npart-1)=ichunk+2*ighost
            ichunkarray(npart)=ncell-ichunk*(npart)+ighost

            ichunktransfer(0:npart-1)=ichunk
            ichunktransfer(npart)=ncell-ichunk*(npart)
            
            istartarray(0)=1
            iendarray(0)=ichunk+ighost
            
            istarttransfer(0)=1
            iendtransfer(0)=ichunk
            
            imoveupstart(0)=ichunk-ighost+1
            imoveupend(0)=ichunk
            
            imovedownstart(0)=0
            imovedownend(0)=0
            do i=1,npart-1
                istartarray(i)=ichunk*(i)+1-ighost
                iendarray(i)=ichunk*(i+1)+ighost
                
                istarttransfer(i)=ichunk*(i)+1
                iendtransfer(i)=ichunk*(i+1)
                
                imoveupstart(i)=ichunk*(i+1)+1-ighost
                imoveupend(i)=ichunk*(i+1)
                
                imovedownstart(i)=ichunk*(i)+1
                imovedownend(i)=ichunk*(i)+ighost
            enddo
                istartarray(npart)=ncell-ichunkarray(npart)+1
                iendarray(npart)=ncell
                istarttransfer(npart)=ncell-ichunktransfer(npart)+1
                iendtransfer(npart)=ncell
                
                imoveupstart(npart)=0
                imoveupend(npart)=0
                
                imovedownstart(npart)=ichunk*(npart)+1
                imovedownend(npart)=ichunk*(npart)+ighost
            print*,'istartarray',istartarray
            print*,'iendarray',iendarray
            print*,'istarttransfer',istarttransfer
            print*,'iendtransfer',iendtransfer
            print*,'imoveupstart',imoveupstart
            print*,'imoveupend',imoveupend
            print*,'imovedownstart',imovedownstart
            print*,'imovedownend',imovedownend
        endif
!        call MPI_BCAST(istartarray,nsize,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
!        call MPI_BCAST(iendarray,nsize,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
        
!        call MPI_BCAST(ichunkarray,nsize,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
!        call MPI_BCAST(ichunktransfer,nsize,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
        
!        call MPI_BCAST(istarttransfer,nsize,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
!        call MPI_BCAST(iendtransfer,nsize,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
        
!        call MPI_BCAST(imoveupstart,nsize,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
!        call MPI_BCAST(imoveupend,nsize,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
!        call MPI_BCAST(imovedownstart,nsize,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
!        call MPI_BCAST(imovedownend,nsize,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
!==============================    
    allocate(ic_no(0:npart))
    allocate(ic_location(0:npart,int(mnhl/nsize*4)))
    do irank=0,npart    
       ic_no(irank)=0
    enddo
    
    do i=1,cnode
        ic=cn_cell(i)
        do irank=0,npart
            if(ic.ge.istarttransfer(irank) .and. ic.le.iendtransfer(irank))then
                ic_no(irank)=ic_no(irank)+1
                ic_location(irank,ic_no(irank))=ic
            endif
        enddo
    enddo 
!    write(*,*)'cnode',cnode
    temp_no=0
    do irank=0,npart
       temp_no=temp_no+ic_no(irank)
!       write(*,*)'ic_no',ic_no(irank)
!       if(rank.eq.0)write(*,*)'------rank, ic_no, ic_location',irank, ic_no(irank),ic_location(irank,1:ic_no(irank))
    enddo
!    write(*,*)'temp_no',temp_no
	endsubroutine store_id_info
