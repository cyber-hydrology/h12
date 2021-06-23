	subroutine regularization
	
	use globals
	implicit none
	
	integer:: i,j,k
	integer:: lev, lev_nei, in, jn, c_reg, npos(8,2)
	
	npos=reshape((/0,1,0,-1,1,1,-1,-1,1,0,-1,0,1,-1,-1,1/),(/8,2/))
!	npos=reshape((/0,1,0,-1,1,1,-1,-1,-1,0,1,0,-1,0,0,-1/),(/8,2/))
	
	c_reg=1
	
	do while(c_reg/=0)
	  c_reg=0
		do i=1, im
			do j=1, jm
            if(bs_dem(i,j)<-9998)cycle
				do k=1, 8
				  lev=c_lev(i,j)
				  in=i+npos(k,1); jn=j+npos(k,2)

				  if(in>im.or.in<1.or.jn>jm.or.jn<1)cycle

				  lev_nei=c_lev(in,jn)
				  
		  		    if(lev_nei>lev+1)then
				      c_lev(i,j)=lev_nei-1
					  c_reg=1
					elseif(lev_nei<lev-1)then
					  c_lev(in,jn)=lev-1
					  c_reg=1
					endif

				enddo
			enddo
		enddo
	enddo
	
	endsubroutine regularization
					
