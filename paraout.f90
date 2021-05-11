!!************************************	
	subroutine paraout(t,ci)
!************************************	    
	use globals
      use mpi    
    implicit none
    
	integer:: i, j, is, js, maxnode, ci
    integer:: ms, lev
    real*8:: x, y, dxl, dyl, x1, y1, s0, x_end, t, z_end, zz
	character(len=5):: out_time
	character(len=3):: cases
      integer:: tag(1:100)=(/(j,j=1,100)/), status(MPI_STATUS_SIZE)    
!	========================================================    
    if(rank.ne.0)then
        call MPI_SEND(h(istarttransfer(rank):iendtransfer(rank)),ichunktransfer(rank), &
        MPI_REAL8,0,tag(4), MPI_COMM_WORLD,ierror)
        print*,'mpi_send h from rank',rank,'ichunktransfer',ichunktransfer(rank)
    else
        do i=1,npart
        call MPI_RECV(h(istarttransfer(i):iendtransfer(i)),ichunktransfer(i), &
        MPI_REAL8,i,tag(4), MPI_COMM_WORLD,status,ierror)
        print*,'mpi_recv h at rank',rank,'i,ichunktransfer',i,ichunktransfer(i)
        enddo
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
       if(rank.eq.0)then
!$$$$$$       write(out_time,'(i5)') nint(t*1000.d0)+10000
       write(out_time,'(i5)') nint(t)+10000
	  write(cases,'(i3)') ci+100
   
!      open(unit=30,file='results/para/surface/paraview_'//cases//out_time//'.vtk')
      open(unit=30,file='/tmp/paraview_'//cases//out_time//'.vtk')
		  write(30,'(a26)')'# vtk DataFile Version 2.0'
		  write(30,'(a25)')'Unstructured Grid Example'
		  write(30,'(a5)')'ASCII'
		  write(30,'(a25)')'DATASET UNSTRUCTURED_GRID'
          maxnode=0
          do i=1, im
            	do j=1, jm
                 if(bs_dem(i,j)<-9998)cycle
                lev=c_lev(i,j); Ms=2**lev
                    do is=1, ms
                      	do js=1, ms
           				 maxnode=maxnode+4
                        enddo
                    enddo
                enddo
           enddo            	
		  write(30,'(a6,1x,i10,1x,a5)')'POINTS',maxnode, 'float'
          
          do i=1, im
            x=x_ori+(dble(i)-0.5d0)*dx
                do j=1, jm
                 if(bs_dem(i,j)<-9998)cycle
                  y=y_ori+(dble(j)-0.5d0)*dy
                  lev=c_lev(i,j); ms=2**lev
                  dxl=dx/DBLE(2**lev); dyl=dy/DBLE(2**lev)
                    do is=1, ms
                      x1=X+DBLE(is-1)*dxl
                        do js=1, ms
                         y1=Y+DBLE(js-1)*dyl
                        write(30,'(2(F11.2,2x), a3)')x1,y1, '0.0'
                        write(30,'(2(F11.2,2x), a3)')x1,y1+dyl, '0.0'
                        write(30,'(2(F11.2,2x), a3)')x1+dxl,y1+dyl, '0.0'
                        write(30,'(2(F11.2,2x), a3)')x1+dxl,y1, '0.0'
                        enddo
                    enddo
                enddo
           enddo 
           write(30,'(a5,1x, I10, 1X, I10)')'CELLS', maxnode/4, 5*(maxnode/4)
           
            do i=1, maxnode, 4
                write(30, '(a1, 2x, 4(i10,2x))') '4', i-1, i, i+1, i+2
            enddo

	         WRITE(30,'(A10, 1X, I10)')'cell_types',  (maxnode/4)
	
			do i=1, maxnode, 4
              	write(30,'(a1)') '9'
            enddo

!$$$$$$ !-----------------------write cell data
             WRITE(30,'(A9, 1X, I10)')'cell_data',  (maxnode/4)
                          
!             WRITE(30,'(a21)') 'SCALARS baseo float 1'
!             WRITE(30,'(a18)') 'LOOKUP_TABLE baseo'
!             do i=1, ncell
!                WRITE(30,'(f15.3)') baseo(i)
!             enddo
             
             WRITE(30,'(a17)') 'SCALARS h float 1'
             WRITE(30,'(a17)') 'LOOKUP_TABLE h'
             do i=1, ncell
                WRITE(30,'(f15.3)') h(i)
             enddo              
!$$$$$$ 
!$$$$$$              WRITE(30,'(a23)') 'SCALARS H float 1'
!$$$$$$              WRITE(30,'(a20)') 'LOOKUP_TABLE infy'
!$$$$$$              do i=1, ncell
!$$$$$$                 WRITE(30,'(f15.3)') baseo(i)+h(i)
!$$$$$$              enddo
             
            CLOSE(30) 
        endif            
	end subroutine paraout

!	========================================
!		subroutine para_pipe    
!	========================================
    subroutine para_pipe(t,ci)
    use globals
      use mpi    
    implicit none
    
    integer:: i, j, k, nm,ci
    real*8::bs,t
    character(len=5):: out_time
    character(len=3):: cases

    if(rank.eq.0)then
    write(out_time,'(i5)') nint(t)+10000
    write(cases,'(i3)') ci+100

    open(unit=600,file='results/para/pipe/pipe_para_'//cases//out_time//'.vtk')
              
    write(600,'(a26)')'# vtk DataFile Version 2.0'
    write(600,'(a25)')'Unstructured Grid Example'
    write(600,'(a5)')'ASCII'
    write(600,'(a25)')'DATASET UNSTRUCTURED_GRID'
    write(600,'(a6,1x,i7,1x,a5)')'POINTS',num_pmax, 'float'
            
    nm=0
    do i=1, iswr
      do j=1, jswr(i)
        if(mark_sw_check(i)==1)then            
          if(j==1)then
            do k=1, 4
              nm=nm+1
              bs = bsdw_sw(i)+slp_sw(i)*jswr(i)*dx_sw(i)
               if(k>2)bs = bsdw_sw(i) + slp_sw(i)*(jswr(i)-1)*dx_sw(i)
	             write(600,'(3(F13.3,2x))') xx(nm), yy(nm), bs
            enddo
          else
            do k=1, 2
              nm=nm+1
              bs = bsdw_sw(i) + slp_sw(i)*(jswr(i)-j)*dx_sw(i)
              write(600,'(3(F13.3,2x))') xx(nm), yy(nm), bs
            enddo
          endif
        endif
      enddo
    enddo
            write(600,'(a5, I7, 1X, I7)')'CELLS', num_em, 5*num_em

           
     	nm=0
        do i=1, iswr
          do j=1, jswr(i)
			if(mark_sw_check(i)==1)then
        	    if(j==1)then
		         write(600, '(a1,2x,4(i7,2x))') '4', nm, nm+1, nm+2, nm+3
                 nm=nm+4
                else
		         write(600, '(a1,2x,4(i7,2x))') '4', nm-1, nm-2, nm, nm+1
                 nm=nm+2
	            endif
			endif
          enddo
        enddo
            


            write(600,'(A10, 1X, I7)')'cell_types',  num_em
    
        do i=1, num_em
            write(600,'(a1)') '9'
        enddo
!-----------------------write cell data         
             WRITE(600,'(A10, 1X, I7)')'cell_data',  num_em
             
             WRITE(600,'(a23)') 'SCALARS z float 1'
             WRITE(600,'(a20)') 'LOOKUP_TABLE z'
             do i=1, iswr
            if(mark_sw_check(i)==1)then                
               do j=1, jswr(i)
                WRITE(600,'(f15.3)') bs_sw(i,j)
               enddo
            endif
             enddo

             WRITE(600,'(a23)') 'SCALARS h float 1'
             WRITE(600,'(a20)') 'LOOKUP_TABLE h'
             do i=1, iswr
            if(mark_sw_check(i)==1)then                
               do j=1, jswr(i)
                WRITE(600,'(f15.3)') h_sw(i,j)
!                if(isnan(h_sw(i,j)))write(*,*)'para_pipe: time,i,j,h_sw(i,j)',time,i,j,h_sw(i,j)
               enddo
            endif
             enddo
             WRITE(600,'(a23)') 'SCALARS H float 1'
             WRITE(600,'(a20)') 'LOOKUP_TABLE H'
             do i=1, iswr
            if(mark_sw_check(i)==1)then                
               do j=1, jswr(i)
                WRITE(600,'(f15.3)') h_sw(i,j)+bs_sw(i,j)
               enddo
            endif
             enddo
             
             WRITE(600,'(a23)') 'SCALARS i integer 1'
             WRITE(600,'(a20)') 'LOOKUP_TABLE i'
        do i=1, iswr
        if(mark_sw_check(i)==1)then        
           do j=1, jswr(i)
                write(600,'(i5)') i
            enddo
        endif
        enddo

             WRITE(600,'(a23)') 'SCALARS x float 1'
             WRITE(600,'(a20)') 'LOOKUP_TABLE x'
        do i=1, iswr
        if(mark_sw_check(i)==1)then        
           do j=1, jswr(i)
                write(600,'(f15.3)') x_m(i,j)
            enddo
        endif
        enddo    

             WRITE(600,'(a23)') 'SCALARS y float 1'
             WRITE(600,'(a20)') 'LOOKUP_TABLE y'
        do i=1, iswr
        if(mark_sw_check(i)==1)then        
           do j=1, jswr(i)
                write(600,'(f15.3)') y_m(i,j)
            enddo
        endif
        
        enddo    

           CLOSE(600)             
       endif
       
    end subroutine para_pipe

!	========================================
!		subroutine para_mh  
!	========================================
    subroutine para_mh(t,ci)
    use globals
      use mpi     
    implicit none
    
    integer:: i,j,k,ci
    real*8::t
    character(len=5):: out_time
    character(len=3):: cases
    
    if(rank.eq.0)then
    write(out_time,'(i5)') nint(t)+10000
    write(cases,'(i3)')ci
    
    open(unit=610,file='results/para/mh/pipe_mh_'//cases//out_time//'.vtk')
    
    write(610,'(a26)')'# vtk DataFile Version 2.0'
    write(610,'(a25)')'Unstructured Grid Example'
    write(610,'(a5)')'ASCII'
    write(610,'(a25)')'DATASET UNSTRUCTURED_GRID'
    write(610,'(a6,1x,i7,1x,a5)')'POINTS',mh_no, 'float'
            
    do i=1, mnhl
      if(mark_mh_check(i)==1)then        
	    do j=1, 4              
          write(610,'(2(F13.3,2x), a3)') xx_mh(i,j), yy_mh(i,j), '0.0'
        enddo
      endif
    enddo

    write(610,'(a5, I7, 1X, I7)')'CELLS', mh_el, 5*mh_el

		k=0
		do i=1, mnhl
			if(mark_mh_check(i)==1) k=k+1        
		enddo
           
        do i=1, k*4, 4
            write(610,'(a1, 2x, 4(i7,2x))')'4', int(i)-1, int(i+1)-1, int(i+2)-1, int(i+3)-1
        enddo


            write(610,'(A10, 1X, I7)')'cell_types',  mh_el
    
        do i=1, mh_el
            write(610,'(a1)') '9'
        enddo
!-----------------------write cell data         
             WRITE(610,'(A10, 1X, I7)')'cell_data',  mh_el
             
             WRITE(610,'(a23)') 'SCALARS z float 1'
             WRITE(610,'(a20)') 'LOOKUP_TABLE z'
	         do i=1, mnhl
		        if(mark_mh_check(i)==1)then        
			         write(610,'(f15.3)') bs_mh(i)
		        endif     
        	enddo 

             WRITE(610,'(a23)') 'SCALARS h float 1'
             WRITE(610,'(a20)') 'LOOKUP_TABLE h'
	         do i=1, mnhl
		        if(mark_mh_check(i)==1)then        
			         write(610,'(f15.3)') h_mh(i)
		        endif     
        	enddo 
             WRITE(610,'(a23)') 'SCALARS H float 1'
             WRITE(610,'(a20)') 'LOOKUP_TABLE H'
	         do i=1, mnhl
		        if(mark_mh_check(i)==1)then        
			         write(610,'(f15.3)') h_mh(i)+bs_mh(i)
		        endif     
        	enddo               
              
             WRITE(610,'(a23)') 'SCALARS i integer 1'
             WRITE(610,'(a20)') 'LOOKUP_TABLE i'
	         do i=1, mnhl
		        if(mark_mh_check(i)==1)then        
			         write(610,'(i7)') i
		        endif     
        	enddo                  

           CLOSE(610)   
      endif           

    end subroutine para_mh  
