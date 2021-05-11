!!     ******************************************************************
      subroutine forward
      use globals! , only:umo,um,vno,vn,ho,h,um_fine,vn_fine,&
                 !       qlme,cn_cell,cnode,iswr,jswr,&
                 !       qo_sw,q_sw,ao_sw,a_sw,ql_sw,mnhl,&
                 !       ho_mh,h_mh,qsum_mh,hsum_mh,&
                 !       rank,istartarray,iendarray,ierror
      use omp_lib
      use mpi
      implicit none
      
      integer:: i,j,nci
    real*8:: tstart,tfinish
      integer:: tag(1:100)=(/(j,j=1,100)/), status(MPI_STATUS_SIZE),tempstart,tempend
!    call MPI_BARRIER(MPI_COMM_WORLD,ierror)   
    tempstart=istartarray(rank)
    tempend=iendarray(rank)
!    if(rank.eq.1)tstart = MPI_Wtime()
!$omp barrier
!$omp parallel do
    do nci=tempstart,tempend
      umo(nci)= um(nci); vno(nci)=vn(nci)
      ho(nci) = h(nci)
      um_fine(nci)=0.0d0
      vn_fine(nci)=0.0d0
      qlme(nci)=0.0d0
    enddo
!$omp end parallel do
!$omp barrier

    if(rank.eq.0)then
!$omp barrier
!$omp parallel do
    do i=1, cnode
      qlme(cn_cell(i))=0.0d0
    enddo
!$omp end parallel do
!$omp barrier

!!$omp barrier
!!$omp parallel do
!    do i = 1, iswr
!      do j = 1, jswr(i)
!        qo_sw(i, j) = q_sw(i, j)
!        ao_sw(i, j) = a_sw(i, j)
!        ql_sw(i, j) = 0.0d0
!      enddo
!        qo_sw(i, jswr(i) + 1) = q_sw(i, jswr(i) + 1)
!    enddo
!!$omp end parallel do
!!$omp barrier
!    do i = 1, mnhl
!      ho_mh(i) = h_mh(i)
!      qsum_mh(i) = 0.0d0
!      hsum_mh(i) = 0.0d0
!    enddo
    endif
    end subroutine forward
