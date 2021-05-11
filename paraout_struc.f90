!************************************	
	subroutine paraout_struc(t,ci)
!************************************	    
      use globals
      implicit none

      integer::i,j, nm, ci
      real*8 xnode, ynode, bs, t
      character(len=6)::ft

      write(ft,'(i6)') int(time)+100000

     open(200, file='results/para_struc/paraout_struc_'//ft//'.vtk', action='write') 
      
	 write(200,'(a26)')'# vtk DataFile Version 2.0'
	 write(200,'(a9)')'Arlington'
	 write(200,'(a5)')'ASCII'
	 write(200,'(a25)')'DATASET STRUCTURED_POINTS'
	 write(200,'(a10,1x,i10,1x,i10,1x,a1)')'DIMENSIONS', im, jm, '1'
	 write(200,'(a18)')'ASPECT_RATIO 1 1 1'
	 write(200,'(a12)')'ORIGIN 0 0 0'
	 write(200,'(a10,1x,i10)')'POINT_DATA',im*jm
	 write(200,'(a21)')'SCALARS baseo float 1'
	 write(200,'(a20)')'LOOKUP_TABLE default'
       do j=1,jm
         do i=1, im  
           write(200,'(f11.2,1x)') bs_dem(i,j)
         enddo
     enddo
!$$$$$$      write(200,'(a30)')'SCALARS depth float 1'
!$$$$$$      write(200,'(a20)')'LOOKUP_TABLE default'     
!$$$$$$        do i=im, 1, -1     
!$$$$$$          do j=1, jm  
!$$$$$$            write(200,'(f8.3,1x)') h(i,j)
!$$$$$$          enddo
!$$$$$$      enddo
 
     close(200)
     endsubroutine paraout_struc