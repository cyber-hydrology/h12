!     ====================================
!         PREPARING NEXT CALCULATION
!     ====================================
      subroutine mdvel
      use globals
      use omp_lib
      use mpi    
      implicit none

    real*8:: h_lc, h_n, h_dc, baseo_lc, baseo_dc
    real*8:: uu1_rc, uu1_lc, uu1_rdc, uu1_dc
    real*8:: vv1_uc, vv1_luc, vv1_lc, vv1_dc
    real*8:: umo_rc, umo_lc, umo_uc, umo_dc
    real*8:: vno_rc, vno_lc, vno_uc, vno_dc      
    integer:: nci, hh, lc, rc, uc, dc, ruc, rdc, ldc, luc
    integer:: ndir, i, j
    real*8:: phi,tstart,tfinish,dd1_sw_tmp
      integer:: recv_request(100),send_request(100),recv_request2(100),send_request2(100)
      integer:: tag(1:100)=(/(j,j=1,100)/), status(MPI_STATUS_SIZE),tempstart,tempend
!     =================================
!        hm, hn calculation
!     =================================
    tempstart=istartarray(rank)
    tempend=iendarray(rank)
!      if(rank.eq.0)tstart = omp_get_wtime()
!$omp barrier
!$omp parallel do private(h_lc,h_dc,lc,rc,uc,dc)
      do nci=tempstart,tempend
      lc = nei_info(nci,4); rc = nei_info(nci,2)
      uc = nei_info(nci,1); dc = nei_info(nci,3)
        if(dc==0)then
          h_dc=0.0d0
        else
          h_dc=h(nei_info(nci,3))              
        endif

        if(lc==0)then
          h_lc=0.0d0
        else        
          h_lc=h(nei_info(nci,4))
        endif
!======================================================
        if(lc/=0.and.dc/=0)then
          hm(nci) = 0.5d0*(h(nci) + h_lc)
          hn(nci) = 0.5d0*(h(nci) + h_dc)
        elseif(lc==0)then
          hm(nci) = h(nci)
            if(dc/=0)then
              hn(nci)=0.5d0*(h(nci) + h_dc)
            else
              hn(nci) = h(nci)
            endif
        elseif(dc==0)then
          if(lc/=0)hm(nci) = 0.5d0*(h(nci)+h_lc)
          hn(nci) = h(nci)
        endif   
!-------------------------------------------------------
!                  ghost cell
!-------------------------------------------------------
      if(rc==0) then
        hm_rgc(nci) = h(nci)
        um_rgc(nci) = 0.0d0
      endif
    
      if(uc==0) then
        hn_ugc(nci) = h(nci)
        vn_ugc(nci) = 0.0d0
      endif
    
    enddo
!$omp end parallel do
!$omp barrier
!
!     =================================
!        uu1, vv1 calculation
!     =================================
!$omp barrier
!$omp parallel do private(lc,rc,uc,dc)
	  do nci=tempstart,tempend
      lc = nei_info(nci,4); rc = nei_info(nci,2)
      uc = nei_info(nci,1); dc = nei_info(nci,3)
     
        if(hm(nci)>=th)then
          uu1(nci)=um(nci)/hm(nci)
        else
          uu1(nci)=0.0d0
        endif

        if(hn(nci)>=th)then
          vv1(nci)=vn(nci)/hn(nci)
        else
          vv1(nci)=0.0d0
        endif
!-------------------------------------------------------
!                  ghost cell
!-------------------------------------------------------
        if(rc==0)then
          if(hm_rgc(nci)>=th)then
            uu1_rgc(nci) = um_rgc(nci)/hm_rgc(nci)
          else
            uu1_rgc(nci) = 0.0d0
          endif
        endif

      if(uc==0)then
        if(hn_ugc(nci)>=th)then
          vv1_ugc(nci) = vn_ugc(nci)/hn_ugc(nci)
        else
          vv1_ugc(nci) = 0.0d0
        endif
      endif
    enddo
!$omp end parallel do
!$omp barrier
!     =================================
!        uu, vv calculation
!     =================================
!$omp barrier
!$omp parallel do private(vv1_uc,uu1_rc,lc,rc,uc,dc)
	  do nci=tempstart,tempend
      lc = nei_info(nci,4); rc = nei_info(nci,2)
      uc = nei_info(nci,1); dc = nei_info(nci,3)
        if(uc==0)then
          vv1_uc=0.0d0
        else
          vv1_uc=vv1(nei_info(nci,1))              
        endif
!
        if(rc==0)then
          uu1_rc=0.0d0
        else        
          uu1_rc=uu1(nei_info(nci,2))
        endif
!          
        if(h(nci)>=th)then
          if(rc/=0)then
            uu(nci)=(uu1_rc+uu1(nci))/2.0d0
          else
            uu(nci)=uu1(nci)
          endif
            
          if(uc/=0)then
            vv(nci)=(vv1_uc+vv1(nci))/2.0d0
          else
            vv(nci)=vv1(nci)
          endif
        else
          uu(nci)=0.0d0; vv(nci)=0.0d0
        endif
    enddo
!$omp end parallel do
!$omp barrier
!     =================================
!        uua, vva calculation
!     =================================
!$omp barrier
!$omp parallel do private(vv1_uc,uu1_rc,uu1_dc,vv1_lc,uu1_rdc,vv1_luc, &
lc,rc,uc,dc,ruc,rdc,ldc,luc)
!    do nci=1, ncell
	  do nci=tempstart,tempend
      lc = nei_info(nci,4); rc = nei_info(nci,2)
      uc = nei_info(nci,1); dc = nei_info(nci,3)
      ruc = nei_info(nci,5); rdc = nei_info(nci,6)
      ldc = nei_info(nci,7); luc = nei_info(nci,8)
!          
        if(uc==0)then
          vv1_uc=vv1_ugc(nci)
        else
          vv1_uc=vv1(nei_info(nci,1))              
        endif

        if(rc==0)then
          uu1_rc=uu1_rgc(nci)
        else        
          uu1_rc=uu1(nei_info(nci,2))
        endif

        if(dc==0)then
          uu1_dc=0.0d0
        else        
          uu1_dc=uu1(nei_info(nci,3))
        endif

        if(lc==0)then
          vv1_lc=0.0d0
        else        
          vv1_lc=vv1(nei_info(nci,4))
        endif

        if(rdc==0)then
          uu1_rdc=0.0d0
        else        
          uu1_rdc=uu1(nei_info(nci,6))
        endif

        if(luc==0)then
          vv1_luc=0.0d0
        else
          vv1_luc=vv1(nei_info(nci,8))
        endif       
!          
        if(dc/=0)then
          uua(nci)= 0.25d0*(uu1(nci) + uu1_rc + uu1_dc + uu1_rdc)
        else
          uua(nci)= 0.5d0*(uu1(nci)+uu1_rc)
        endif
        
        if(uc==0) uua_ugc(nci) =0.5d0*( uu1(nci)+ uu1_rc )
        if(lc/=0)then
          vva(nci) = 0.25d0*(vv1(nci) + vv1_uc + vv1_lc + vv1_luc)
        else
          vva(nci) = 0.5d0*(vv1(nci)+vv1_uc)
        endif
        if(rc==0) vva_rgc(nci) = 0.5d0*(vv1(nci)+vv1_uc)
          
    enddo
!$omp end parallel do
!$omp barrier
  
!    **************************************************
    end subroutine mdvel
