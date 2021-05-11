    module globals
 
    implicit none
     
    ! general data
    real*8, parameter:: emsmall = 1.0d-12, g=9.81d0, tol_h=1.0d-5, hmin=1.0d-4
    integer,save:: im, jm, imori, jmori, ims, jms
    real*8,save:: gg, dt, dt2, time, th, pi, rn, dxori,dyori,dx, dy, dxr, dyr
    real*8,save:: xl, yl, x_man, y_man, slope, d_man, dtlevel
    real*8,save:: x_ori, y_ori, xr_ori, yr_ori

    real*8,allocatable:: rain(:,:), rain_p(:) !after get rainfall data
    real*8::dtrain,x_rad_ori, y_rad_ori, rad_cell
    integer::nmrain,irm,jrm

    ! flow data
    integer,allocatable:: inf(:), infx(:), infy(:), mark(:), c_land(:)
    real*8,allocatable:: baseo(:), h(:), ho(:), hm(:), hn(:)
    real*8,allocatable:: um(:), umo(:), uu(:), uua(:), uu1(:)
    real*8,allocatable:: um_fine(:), vn_fine(:)
    real*8,allocatable:: vn(:), vno(:), vv(:), vva(:), vv1(:),qlme(:)

    ! flow boundary data
    real*8,allocatable:: um_rgc(:), hm_rgc(:), uu1_rgc(:), vva_rgc(:) ! east boundary
    real*8,allocatable:: vn_ugc(:), hn_ugc(:), vv1_ugc(:), uua_ugc(:) ! north boundary
    real*8,allocatable:: dwlevel(:) !down stream sea level as boundary condition

    !** grid data
    integer:: nps, npe ! total number of seeding points
    integer:: maxcell, div_max, maxseed,upscale
    integer:: ncell, c_lev_max
    integer,allocatable:: c_lev(:,:), c_lev_id(:)
    integer,allocatable:: ind(:,:,:,:), nei_info(:,:)
    integer,allocatable::mark_ori_land(:,:),mark_land(:,:)
    real*8,allocatable:: xs(:), ys(:), ccoords(:,:)
    real*8,allocatable:: bs_rd(:,:), bs_dem(:,:), h_grid(:,:),bs_ori_dem(:,:)
	real*8,allocatable:: um_grid(:,:),vn_grid(:,:),ini_grid(:,:)

    ! manhole parameter
    real*8, allocatable:: x_mh(:), y_mh(:), bs_mh(:), area_mh(:), a_mh(:)
    real*8, allocatable:: h_mh(:), ho_mh(:), qsum_mh(:), hsum_mh(:)
    real*8, allocatable::  xx_mh(:,:), yy_mh(:,:)
    integer, allocatable:: isw_mh(:,:),inf_mh(:), ngrp(:)
    integer, allocatable:: check_mh(:)
    character(len=50),allocatable::id_mh(:)
    integer:: mnhl
    real*8:: ac_cnq, ac_qp, ac_qm

    ! pipe parameter
    real*8, allocatable:: dst_sw(:), rn_sw(:), bsdw_sw(:), slp_sw(:), dx_sw(:), bs_sw(:,:), dd1_sw(:), dd2_sw(:), bsl(:)
    real*8, allocatable:: a_sw(:,:), ao_sw(:,:), h_sw(:,:), ai_sw(:,:), r_sw(:,:), ql_sw(:,:)
    real*8, allocatable:: q_sw(:,:), qo_sw(:,:), uu_sw(:,:), x_sw(:,:), y_sw(:,:), xx(:), yy(:)
    real*8, allocatable:: hmm(:,:), a_box(:,:), bs_box(:,:), h_box(:,:), x_m(:,:), y_m(:,:)
    integer, allocatable:: inf_sw(:), igrp(:), idup_mh(:), iddw_mh(:), main_p(:), side_p(:),cnp_i(:), cnp_j(:)
    integer, allocatable:: idshp_sw(:), ipt_sw(:), mhup_sw(:), mhdw_sw(:), jswr(:), num(:),nm_maintoside(:,:)
    character(len=50),allocatable::id_sw(:)
    integer::nm_mp, nm_sp
    integer::mh_no, mh_el, pm_no, pm_el, num_pmax, num_em
    real*8:: h_spc(2000), a_spc(2000), r_spc(2000)
    integer:: iswr
    
 	!search
    integer, allocatable:: mark_mh(:,:), mark_sw(:,:), mark_sw_check(:), mark_mh_check(:)

    !Pipe mesh generation
    integer,allocatable::num_pg(:)

    !connection parameter
    integer::cnode
    integer,allocatable::cn_cell(:),cn_mh(:)

    !initial data parameter
    real*8,allocatable::ele_ini(:)
    integer,allocatable:: id_ini(:),i_ini(:),j_ini(:),cei(:)
    integer::nm_ini
    real*8::err_ini
    
    
    ! mpi variables
  integer(kind=4), public                              :: rank                ! unique number of each ensemble
  integer(kind=4), public                              :: nsize               ! size of mpi processes
  integer(kind=4), public                              :: npart               ! size of ensemble :: nsize - 1 (rank 0 is not used for ensemble simulation)
  integer(kind=4), public                              :: ierror              ! index for internal function of mpi
  integer(kind=4), public                              :: ichunk
  integer(kind=4), allocatable                         :: ichunkarray(:)
  integer(kind=4), public                              :: istart
  integer(kind=4), public                              :: iend
  integer(kind=4), public                              :: ighost
  integer(kind=4), allocatable                         :: istartarray(:)
  integer(kind=4), allocatable                         :: iendarray(:)
  integer(kind=4), allocatable                         :: istarttransfer(:)
  integer(kind=4), allocatable                         :: iendtransfer(:)
  integer(kind=4), allocatable                         :: ichunktransfer(:)
  integer(kind=4), allocatable                         :: imoveupstart(:)
  integer(kind=4), allocatable                         :: imoveupend(:)
  integer(kind=4), allocatable                         :: imovedownstart(:)
  integer(kind=4), allocatable                         :: imovedownend(:)
  
  integer(kind=4), allocatable                         :: ic_no(:)
  integer(kind=4), allocatable                         :: ic_location(:,:)
    end module globals
