	subroutine seed_find_div
	
	use globals
	implicit none
		
	integer:: i,j,is,k,countup, countdown, intup, intdown
	integer:: kk,mn, mb, nb, ib, jb, npos(4,2)
	real*8:: xc, yc, x1, x2, y1, y2, x, y, dxs, dys, ds
	
	do is=1, npe

    	x=xs(is)-x_ori; y=ys(is)-y_ori
        ib=ceiling(x/dx); jb=ceiling(y/dy)
        if(ib>im) ib=im
        if(jb>jm) jb=jm
        if(bs_dem(ib,jb)<-9998)cycle
	    c_lev(ib,jb)=div_max
	enddo
	
	
	endsubroutine seed_find_div
