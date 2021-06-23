	subroutine sdp
	
	use globals
	implicit none
	
	integer:: k
	real*8::m_circle, dl_circle,theta
	
	m_circle=100.0d0
	dl_circle=2*pi/m_circle
	nps=1    !nps: number of starting point of seeding point
	npe=100  !npe: number of ending point of seeding point
	do k=nps, npe+1
	  theta=(k-1)*dl_circle
	  xs(k)=6.59d0+0.5d0*cos(theta)
	  ys(k)=1.96d0+0.5d0*sin(theta)
	enddo
	
	end subroutine sdp	  
