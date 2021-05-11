	subroutine mesh_gen
	
	use globals
	implicit none 

	integer :: i, j, ms, ic, jc
	
!	initialize c_lev
		
	do i=1, im
	  do j=1, jm
	    c_lev(i,j)=0
		ms=2**div_max
		do ic=1, ms 
		  do jc=1, ms 
		    ind(i,j,ic,jc)=0
		  enddo
		enddo
	  enddo
	enddo
write(*,*)'start seed find'
	call seed_find_div	! finding seed point in meshes
write(*,*)'end seed find'
	call regularization ! making surounded mesh for grade
write(*,*)'end regula'
	call store_id_info  ! finding neighbour mesh and storing id address
write(*,*)'end store_id'
	
	end subroutine mesh_gen
	