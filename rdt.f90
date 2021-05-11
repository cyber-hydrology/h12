!!     =======================================================
!                     READING INPUT DATA
!     =======================================================
      subroutine rdat(time0, timmax, lpout, lkout,ci)
      use globals
	  implicit none
	  
      integer:: i, j, k,lpout, lkout,ci,n,nm, nodata
      real*8:: tday, thour, tmin, tsec, dkout, dpout
      real*8:: time0, timmax, x_ini, y_ini
      character*50::fdem,froad,fland,fdwlevl,frain,fini,skip1,skip2
      character*50::fh,fho,fum,fumo,fvn,fvno,fuu,fvv
      character*50::fuua,fvva,fuu1,fvv1
      character(len=8)::obs_point

!     =====================================================
!                      READING DATA
!     =====================================================
! *****************READING SETUP DATA************************
      open(10, file='input/setup.dat', action='read')
      read(10,*) tday, thour, tmin, tsec
      read(10,*) dt
      read(10,*) dpout
      read(10,*) dkout
      read(10,*) dtrain
      read(10,*) dtlevel
      read(10,'(a50)') fdem
      read(10,'(a50)') fland
      read(10,'(a50)') fdwlevl
      read(10,'(a50)') frain
      close(10)
      print*,"finish reading setup file"
!		      
	! don't use open number 20, it's for rainfall data
      open(20,file=frain, action='read')
      open(30,file=fdem, action='read')
      open(31,file=fland, action='read')
      open(32,file=fdwlevl, action='read')
!----------------initial data--------------
! *****************READING RAINFALL DATA*************************
      print*, "START READING RAINFALL DATA"  
      read(20,*) nmrain ; allocate(rain_p(nmrain))
      read(20,*)
      read(20,*)
      do i=1, nmrain
        read(20,*) skip1, rain_p(i)
      enddo
      close(20)
      print*, "end reading rainfall data"
! *****************READING ELEVATION DATA************************
      print*,"start reading dem file"
      read(30,*)skip1,imori    !i=x-direction
      read(30,*)skip1,jmori    !j=y-dirction
      read(30,*)skip1,x_ori
      read(30,*)skip1,y_ori
      read(30,*)skip1,dxori
      read(30,*)skip1,nodata
      dyori=dxori
!      read(30,*)

      allocate(bs_ori_dem(imori,jmori))

      do j=jmori, 1, -1
        read(30,*) (bs_ori_dem(i,j), i=1, imori)
        do i=1, imori
          if(int(bs_ori_dem(i,j))==nodata)then
             bs_ori_dem(i,j)=-9999.0d0
          else
             bs_ori_dem(i,j)=bs_ori_dem(i,j)*0.001d0
!print*, i, j, bs_ori_dem(i,j)
          endif
        enddo
      enddo
write(*,*)'end reading bs_dem'
      close(30)
! *****************READING LANDCOVER DATA************************
      print*,"start reading landcover file"
      read(31,'(5/)')
      allocate(mark_ori_land(imori,jmori))

      do j=jmori, 1, -1
!        print*, j, jmori 
      read(31,*) (mark_ori_land(i,j), i=1, imori)
        do i=1, imori
           if(mark_ori_land(i,j)==nodata) mark_ori_land(i,j)=-9999.0d0
        enddo
      enddo
write(*,*)'end reading landcover data'
      close(31)
!
!       mark_ori_land(i,j) = 10 => water                  : cd=1.0
!       mark_ori_land(i,j) = 20 => impervious surface     : cd=1.0
!       mark_ori_land(i,j) = 30 => soil & barren          : cd=0.9
!       mark_ori_land(i,j) = 40 => tree & forest          : cd=0.9
!       mark_ori_land(i,j) = 52 => shrubs                 : cd=0.9
!       mark_ori_land(i,j) = 70 => Grass & Herbaceous     : cd=1.0
!       mark_ori_land(i,j) = 80 => Agriculture            : cd=1.0
!       mark_ori_land(i,j) = 91 => Woody Wetlands         : cd=1.0
!       mark_ori_land(i,j) = 92 => Emergent Wetlands      : cd=1.0

!
! **************reading downstream levl data*********************

      read(32,*)
      nm=0
      do 
        read(32,*,iostat=k)skip1
        if(k<0)exit
        nm=nm+1
      enddo
      allocate(dwlevel(nm))
! re-open
      rewind(32)
      read(32,*)
      do i=1, nm
        read(32,*)skip1, skip2, dwlevel(i)
      enddo
      close(32)
! *****************INCEREASE OF CELL SIZE************************
      call inc_cell

!     ===============================================
!               variable defination
!     ===============================================      
      c_lev_max = 0
      div_max=0
      maxcell = im*jm*2**div_max*2**div_max
!     ===============================================
      time0 =0.0d0
      timmax = 3600.0d0*24.0d0*tday + 3600.d0*thour + 60.d0*tmin + tsec
      dt2 = 2.0d0*dt
      th = 1.0d-4
      gg = 9.8d0
      pi = 4.0d0*atan(1.0d0)
!     =====================================================
!            MONITOR OUTPUT AND DISK SAVING TIME
!     =====================================================
      lpout = nint(dpout/dt)
      lkout = nint(dkout/dt)
!      ltrain = nint(dtrain/dt)
! *****************CAL SEEDING DATA************************
!      nps=1
!      npe=0
!      do i=1,imori
!        do j=1,jmori
!          if(mark_ori_land(i,j)==4) then  !road
!            npe=npe+1
!          endif
!        enddo
!      enddo
!      allocate(xs(npe), ys(npe))     
!      npe=0 
!      do i=1,imori
!        do j=1,jmori
!          if(mark_ori_land(i,j)==4) then
!            npe=npe+1
!            xs(npe) = x_ori + (dble(i)-0.5d0)*dxori   !cell center
!            ys(npe) = y_ori + (dble(j)-0.5d0)*dyori   !cell center
!         endif
!        enddo
!      enddo
!write(*,*)'end reading seeding point'
!      maxseed = npe
      end subroutine rdat
	
