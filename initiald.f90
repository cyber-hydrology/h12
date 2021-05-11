!     =======================================================
!                 INITIALIZATION OF VARIABLES
!     =======================================================
      subroutine initiald
      use globals
      
      real*8:: hmn
      integer:: i, j, nci, lc, rc,uc, dc, ndir, lev, ms, cind

        rn = 0.03d0
      do i=1, ncell
        h(i)=0.0d0; ho(i)=0.0d0
        um(i)=0.0d0; umo(i)=0.0d0
        vn(i)=0.0d0; vno(i)=0.0d0
        uu(i)=0.0d0; uua(i)=0.0d0
        vv(i)=0.0d0; vva(i)=0.0d0
        uu1(i)=0.0d0; vv1(i)=0.0d0
        um_fine(i)=0.0d0
        vn_fine(i)=0.0d0
        inf(i)=1
        qlme(i)=0.0d0
      enddo
!
!   -------------------------------------
!       Calculation of infx and infy
!   -------------------------------------
	  do nci=1, ncell
        ! ndir:4=left, ndir:2=right 
        ! ndir:1=up, ndir:3=down
	! x-direction
    	lc = nei_info(nci,4); rc = nei_info(nci,2)
        uc = nei_info(nci,1); dc = nei_info(nci,3)
        
              if(lc==0)then
                hmn = h(nci)
                infx(nci)=1
              else
                hmn = 0.5d0*(h(nci)+h(lc))
                infx(nci) = abs(inf(nci)-inf(lc))
              endif
              hm(nci)=hmn
              
	! y-direction
!	        c1 = nei_info(nci,3); c2=nci
		
              if(dc==0)then
                hmn = h(nci)
                infy(nci) = 1
              else
                hmn = 0.5d0*(h(nci)+h(dc))
                infy(nci) = abs(inf(nci)-inf(dc))
              endif
              hn(nci)=hmn
      enddo

!     ===========================
!            CONDUIT
!     ===========================
!      do i = 1, iswr
!      do j = 1, jswr(i)
!        a_sw(i, j) = 0.0d0  ; ao_sw(i, j) = a_sw(i, j)
!        h_sw(i, j) = 0.0d0
!        q_sw(i, j) = 0.0d0  ; qo_sw(i, j) = 0.0d0 ; uu_sw(i, j) = 0.0d0
!        ai_sw(i, j) = 0.0d0 ; r_sw(i, j) = 0.0d0
!        ql_sw(i,j) = 0.0d0
!      enddo
!        q_sw(i, jswr(i)+1) = 0.0d0  ; qo_sw(i, jswr(i)+1) = 0.0d0
!      enddo
!     ===========================
!             MANHOLE
!     ===========================
!      do n = 1, mnhl
!        h_mh(n) = 0.0d0
!        ho_mh(n) = 0.0d0
!		a_mh(n) = 0.0d0
!         qsum_mh(n) = 0.0d0
!        hsum_mh(n) = 0.0d0
!      enddo
!
!     ===========================
!   accumulated exchange discharge
!     ===========================      
!      ac_qp = 0.0d0
!      ac_qm = 0.0d0
!      ac_cnq = 0.0d0
!      
!     ===========================
!               abox
!     ===========================
!    do i=1, iswr
!      do j=1, 2
!        if(mark_sw_check(i)==1)then
!          a_box(i,j) = 0.0d0
!          h_box(i,j) = 0.0d0
!        endif
!      enddo
!    enddo


!   ========================================================
!                 manhole elevation adjustment
!   ========================================================
!        do i=1, iswr
!       if(mark_sw_check(i)==1)then
!            mh = mhup_sw(i)
!            bs_box(i,1) = bs_sw(i,1) + slp_sw(i)*(dx_sw(i)/2.0d0 + sqrt(area_mh(mh)/pi))
!            if(bs_box(i,1) < bs_mh(mh)) bs_mh(mh) = bs_box(i,1)
            
!            mh = mhdw_sw(i)
!            bs_box(i,2) = bs_sw(i,jswr(i)) - slp_sw(i)*(dx_sw(i)/2.0d0 + sqrt(area_mh(mh)/pi))
!            if(bs_box(i,2) < bs_mh(mh)) bs_mh(mh) = bs_box(i,2)
!        endif
!        enddo

      end subroutine initiald
