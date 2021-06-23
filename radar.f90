!   **************************************
!   *              radar data            *
!   **************************************

    subroutine radar
    use globals
!    use omp_lib
!    use mpi
    implicit none

    
    integer::i,j,day,hour,minute
    real*8::rest, intime
    character(len=6)::cases,chour
    character(len=2)::cday
    character(len=50)::char
   
!    rest=mod(intime,3600.0d0)
    day=int((time)/(24.0d0*3600.0d0))
    hour=int((time-day*24.0d0*3600.0d0)/3600.0d0)  ! hour
    minute=int(time-3600.0d0*hour)/60
    
    write(cday,'(i2)')25+day
    write(chour,'(i6)')201700+hour
    print*,time, 'day=',day, 'hour=', hour, 'min=', minute
    
    open(20,file='input/Harvey-XMRG-UTM/xmrg08'//cday//chour//'z.asc',action='read')

    read(20,*)char, jrm
    read(20,*)char, irm
    read(20,*)char,x_rad_ori
    read(20,*)char,y_rad_ori
    read(20,*)char,rad_cell
    read(20,*)
    
    allocate(rain(irm,jrm))
    
    do i=1,irm
        read(20,*)(rain(i,j),j=1,jrm)
    enddo
    close(20)

    end subroutine radar

       
    
