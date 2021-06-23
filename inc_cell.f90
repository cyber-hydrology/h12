!   ====================================================
!               Increase of Cell Size
 !   ====================================================
    subroutine inc_cell
    use globals
    implicit none

    integer::i,j,is,ie,js,je,ii,jj,nm
    real*8::maxcellsize,total_bs

    upscale = 0
    maxcellsize = dxori*(2**upscale)
    dx=maxcellsize; dy=dx

    im=int(dble(imori)*dxori/dx)
    jm=int(dble(jmori)*dyori/dy)

    allocate(bs_dem(im,jm),mark_land(im,jm),h_grid(im,jm))
    allocate(um_grid(im,jm),vn_grid(im,jm),ini_grid(im,jm))

    do i=1, im
      do j=1, jm
        is=int((dble(i-1)*dx)/dxori)+1; ie=int((dble(i)*dx)/dxori)
        js=int((dble(j-1)*dy)/dyori)+1; je=int((dble(j)*dy)/dyori)
        total_bs=0.0d0; nm=0
        do ii=is,ie
          do jj=js,je
            if(bs_ori_dem(ii,jj)>-9998)then
              total_bs=total_bs+bs_ori_dem(ii,jj)
                nm=nm+1
            endif
          enddo
        enddo

        if(nm>0)then
          bs_dem(i,j)=total_bs/dble(nm)
        else
          bs_dem(i,j) = -9999
        endif
            
        total_bs=0.0d0; nm=0
      
        do ii=is,ie
          do jj=js,je
            if(mark_ori_land(ii,jj)>-9998)then
              if(mark_ori_land(ii,jj)==10)then
                mark_land(i,j)=10
                nm=nm+1
                goto 100   ! to be modified for selecting mark_land
              elseif(mark_ori_land(ii,jj)==10)then
                mark_land(i,j)=10
                nm=nm+1
              elseif(mark_ori_land(ii,jj)==30)then
                mark_land(i,j)=30
                nm=nm+1
              elseif(mark_ori_land(ii,jj)==40)then
                mark_land(i,j)=40
                nm=nm+1
              elseif(mark_ori_land(ii,jj)==52)then
                mark_land(i,j)=52
                nm=nm+1
              elseif(mark_ori_land(ii,jj)==70)then
                mark_land(i,j)=70
                nm=nm+1
              elseif(mark_ori_land(ii,jj)==80)then
                mark_land(i,j)=80
                nm=nm+1
              elseif(mark_ori_land(ii,jj)==91)then
                mark_land(i,j)=91
                nm=nm+1
              elseif(mark_ori_land(ii,jj)==92)then
                mark_land(i,j)=92
                nm=nm+1
              endif
            endif
          enddo
100     enddo
       
        if(nm==0)mark_land(i,j) = -9999
      enddo
    enddo
!$$$$$$       mark_land(i,j) = 0 => bare ground    : cd=0.6
!$$$$$$       mark_land(i,j) = 1 => vegetation     : cd=0.5
!$$$$$$       mark_land(i,j) = 2 => building       : cd=0.9
!$$$$$$       mark_land(i,j) = 3 => parking lots   : cd=0.9
!$$$$$$       mark_land(i,j) = 4 => road           : cd=0.9
!$$$$$$       mark_land(i,j) = 5 => water          : cd=1.0
	
 	end subroutine inc_cell
