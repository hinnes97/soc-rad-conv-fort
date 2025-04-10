module eos_dry_adj
  use interpolation, only: interpolate
  use netcdf
  use io, only : handle_err
  use params, only: dp, Finc, accelerate, inhibited, q0
  use tables, only: find_var_lin, satur, Tcrit, Pcrit
  use phys, only: Rstar, H2O
  use adjust_mod, only: gradient
  use atmosphere, only: cp_dry, mmw_dry, nqt, q_orig
  use q_sat_test, only : q_sat, dew_point_T
  
  implicit none

  integer :: nt, np
  real(dp), dimension(:), allocatable :: T_dat
  real(dp), dimension(:), allocatable :: P_dat

  real(dp), dimension(:,:), allocatable :: s_h2o
  real(dp), dimension(:,:), allocatable :: s_h2he

  real(dp), dimension(:,:), allocatable :: s_p_h2o
  real(dp), dimension(:,:), allocatable :: s_p_h2he

  real(dp), dimension(:,:), allocatable :: s_t_h2o
  real(dp), dimension(:,:), allocatable :: s_t_h2he

contains
  subroutine init_data
    ! Read in data for H2/He/H2O equations of state
    character(100) :: file_name
    integer :: i, n, id, err, status, ncid, id_t, id_p, ndims, dimids(2)
    integer :: id_sh2o, id_sh2he, id_s_p_h2o, id_s_t_h2o, id_s_t_h2he, id_s_p_h2he

    file_name = '/network/group/aopp/planetary/RTP026_INNES_MSN/soc-rad-conv-fort/fortran/s_h2he_h2o.nc'

    ! Open file
    status = nf90_open(file_name, 0, ncid)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_inq_dimid(ncid, 'logT', id_t)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_inq_dimid(ncid, 'logP', id_p)
    if (status /= nf90_noerr) call handle_err(status)
    
    status = nf90_inquire_dimension(ncid, id_t, len=nt)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_inquire_dimension(ncid, id_p, len=np)
    if (status /= nf90_noerr) call handle_err(status)

    allocate(T_dat(nt), P_dat(np),&
         s_h2o(np, nt), s_h2he(np, nt), s_p_h2o(np, nt), s_t_h2o(np, nt),&
         s_p_h2he(np, nt), s_t_h2he(np, nt))

    status = nf90_inq_varid(ncid, 'logT', id_t)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_inq_varid(ncid, 'logP', id_p)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_inq_varid(ncid, 'logs_h2o', id_sh2o)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_inq_varid(ncid, 'logs_h2he', id_sh2he)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_inq_varid(ncid, 'dlogS_dlogP_h2o', id_s_p_h2o)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_inq_varid(ncid, 'dlogS_dlogP_h2he', id_s_p_h2he)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_inq_varid(ncid, 'dlogS_dlogT_h2o', id_s_t_h2o)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_inq_varid(ncid, 'dlogS_dlogT_h2he', id_s_t_h2he)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_get_var(ncid, id_t, T_dat)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_get_var(ncid, id_p, P_dat)
    if (status /= nf90_noerr) call handle_err(status)

    ! TESTING
    status = nf90_inquire_variable(ncid, id_sh2o, ndims = ndims, dimids=dimids)


    status = nf90_get_var(ncid, id_sh2o, s_h2o)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_get_var(ncid, id_sh2he, s_h2he)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_get_var(ncid, id_s_p_h2o, s_p_h2o)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_get_var(ncid, id_s_p_h2he, s_p_h2he)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_get_var(ncid, id_s_t_h2o, s_t_h2o)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_get_var(ncid, id_s_t_h2he, s_t_h2he)
    if (status /= nf90_noerr) call handle_err(status)

  end subroutine init_data

  subroutine ad_grad(T, P, x_h2o, grad)
    real(dp), intent(in) :: T, P, x_h2o
    real(dp), intent(out) :: grad

    real(dp) :: log10s_h2o, log10s_h2he, gradp_h2o, gradt_h2o, gradp_h2he, gradt_h2he
    real(dp) :: num1, num2, den1, den2
   

    log10s_h2o = interpolate(np, P_dat, nt, T_dat, s_h2o, &
         log10(P), log10(T))

    log10s_h2he = interpolate(np, P_dat, nt, T_dat, s_h2he, &
         log10(P), log10(T))

    gradp_h2o = interpolate(np, P_dat, nt, T_dat, s_p_h2o, &
         log10(P), log10(T))

    gradt_h2o = interpolate(np, P_dat, nt, T_dat, s_t_h2o, &
         log10(P), log10(T))
    
    gradp_h2he = interpolate(np, P_dat, nt, T_dat, s_p_h2he, &
         log10(P), log10(T))

    gradt_h2he = interpolate(np, P_dat, nt, T_dat, s_t_h2he, &
         log10(P), log10(T))


    num1 = x_h2o * 10**(log10s_h2o) * gradp_h2o
    num2 = (1-x_h2o)*10**(log10s_h2he) * gradp_h2he
    den1 = x_h2o * 10**(log10s_h2o) * gradt_h2o
    den2 = (1 - x_h2o)*10**(log10s_h2he) * gradt_h2he

    grad = - (num1 + num2)/(den1 + den2)
    
  end subroutine ad_grad

  subroutine discern_ad_grad(T, P, x_h2o, mmw, cp, grad, grad_type, tstep)
    ! Discern between region where we use R/cp, and region where we use the full equation
    ! of state (depends on if we are in two phase region of water or not)

    real(dp), intent(in) :: T, P, x_h2o, mmw, cp
    real(dp), intent(out) :: grad, grad_type
    integer, intent(in) :: tstep
    real(dp) :: psat

    if (T .gt. Tcrit) then
       call ad_grad(T, P, x_h2o, grad)
       ! if (mod(tstep, 1000) .eq. 0) then
       !    call find_var_lin(T, satur, psat)
       !    write(*,*) 'Using interpolated grad, supercrit', T, P, x_h2o, grad
       ! endif
       grad_type = 1.0
       return
    endif
    
    call find_var_lin(T, satur, psat)

    if (p<psat .and. T .gt. 300) then
       call ad_grad(T, P, x_h2o, grad)
       ! if (mod(tstep, 1000) .eq. 0) then
       !    call find_var_lin(T, satur, psat)
       !    write(*,*) 'Using interpolated grad, p<psat', T, P, x_h2o, grad, psat
       ! endif
       grad_type = 1.0
    else
       grad = Rstar/mmw/cp
       grad_type = 2.
    endif

  end subroutine discern_ad_grad
    
  subroutine adjustment(p, delp, T, q, ktrop, grad, olr, mask, tstep, mmw, cps)
   
    !==========================================================================
    ! Input variables
    !==========================================================================    
     real(dp), intent(in)   , dimension(:) :: p,delp ! Pressure, pressure thickness
     real(dp), intent(in), dimension(:) :: mmw, cps
     
     real(dp), intent(inout), dimension(:) :: grad
     integer, intent(in) :: tstep

    real(dp), intent(in) :: olr
    !==========================================================================
    ! Output variables
    !==========================================================================
    integer, intent(inout) :: ktrop
    integer, intent(out) :: mask(:)
    
    !==========================================================================
    ! Mixed input/output
    !==========================================================================
    real(dp), intent(inout), dimension(:) :: T
    real(dp), intent(inout) :: q(:,:)
    !==========================================================================
    ! Local variables
    !==========================================================================
    ! Tune these parameters as necessasry
    real(dp), parameter          :: delta = 0.0001 ! Small number speeds up convergence

    integer, parameter       :: N_iter = 1000 ! Number of up-down iterations
    
    real(dp) :: pfact,  qcrit,temp
    real(dp) :: qsats(size(p),nqt), qcrits(size(p)), T_old(size(p))

    real(dp) :: grad_check(size(p)), grad_true(size(p)), grad_type(size(p))

    integer :: n,k
    integer :: info, counter
    integer :: npz, m
    logical :: condition, quit_adjust
    real(dp) :: f ! Helps global energy balance to be reached
    !==========================================================================
    ! Main body
    !==========================================================================
    
    npz = size(p)
    !ktrop = npz
    info = 0
    counter = 0
    n = 1
    !    write(*,*) lbound(mask), ubound(mask)
    mask = 0
    quit_adjust = .false.

    ! Uncomment this line and line at bottom for exit at tolerance
    !do while (.not. quit_adjust .and. n .lt. N_iter)

    grad_type = 0.0
    T_old = T
    do n=1,N_iter

       
       grad_check = 0.0_dp
       grad_true = 0.0_dp
       
       call q_sat(p, T, qsats)
       
       do k=npz-1,max(ktrop, 1),-1

          if (n.eq. 1 .and. tstep .gt. 1000 .and. accelerate) then
             f =  (Finc/olr)**(0.01_dp)
          else
             f = 1.
          endif

          if ( q(k+1,1) .gt. q0-1.e-10) then
             ! In dry part
             call discern_ad_grad(T(k+1), P(k+1), q(k+1, 1), mmw(k+1), cps(k+1), grad(k+1), grad_type(k+1), tstep)
             pfact =exp(grad(k+1)*log(p(k)/p(k+1)))
             condition = ( (T(k) - T(k+1)*pfact*(1+delta)) .lt. 0 )

             if (condition) then
                T(k+1) = (T(k)*delp(k) + T(k+1)*delp(k+1))/(delp(k+1) + pfact*delp(k))
                T(k+1) = f*T(k+1)
                T(k) = f*T(k+1)*pfact
                mask(k) = 1
                mask(k+1) = 1
             endif
             cycle
          endif

          ! Continue if the dry adjustment not called
          call gradient(p(k+1),T(k+1),cp_dry(k+1), mmw_dry(k+1), grad(k+1), temp)
          if (n .eq. 1 .and. mod(tstep, 1000) .eq. 0) write(*,*) 'gradients!'
          if (n .eq. 1 .and. mod(tstep, 1000) .eq. 0) then
             write(*,*) k, p(k+1), T(k+1), cp_dry(k), mmw_dry(k+1), grad(k+1)
          endif

          grad_type(k+1) = 3.0
          qcrit = 1./(1._dp - mmw_dry(k+1)/H2O%mmw) /temp
          
          qcrits(k) = qcrit
          ! if (qsats(k+1,1) .gt. 1) then
          !    call dew_point_T(p(k+1), T(k+1))

          !    if ((Finc/olr)**(0.01_dp) .lt. 1) then
          !       T(k+1) = T(k+1)*(Finc/olr)**(0.01_dp)
          !    else
          !       mask(k+1) = 1
          !       cycle
          !    endif
             
             
          ! endif

!       enddo
!    enddo
    
          pfact =exp(grad(k+1)*log(p(k)/p(k+1)))
              
          if (inhibited) then
             ! Do moisture inhibition
             condition = (  (T(k) - T(k+1)*pfact*(1 + delta))*(qcrit - q(k+1,1)) .lt. 0)
          else
             ! Regular criterion
             condition = ( (T(k) - T(k+1)*pfact*(1+delta)) .lt. 0 )
          endif
              
          if (condition .and. q(k+1,1) .gt. qsats(k+1,1) - 1.e-10) then
             T(k+1) = (T(k)*delp(k) + T(k+1)*delp(k+1))/(delp(k+1) + pfact*delp(k))
             T(k+1) = f*T(k+1)
             T(k) = f*T(k+1)*pfact
             q(k, 1) = qsats(k,1)
             q(k+1,1) = qsats(k,1)
             do m=2,nqt
                q(k,m) = (1 - q(k,1))/(1-q_orig(k,1))*q_orig(k,m)
             enddo
             mask(k) = 1
             mask(k+1) = 1
          endif

       enddo !k loop

       ! Code below if wanting to iterate until tolerance reached
!!$           quit_adjust = .true.
!!$           do k=npz-1,max(ktrop, 1),-1
!!$              if (mask(k) .and. mask(k+1)) then
!!$                 call gradient(p(k+1), T(k+1), grad_check(k+1))
!!$                 
!!$                 grad_true(k+1) = log(T(k+1)/T(k))/log(p(k+1)/p(k))
!!$              endif
!!$
!!$              if (abs((grad_true(k+1) - grad_check(k+1))/grad_check(k+1)) .gt. 1.e-5) then
!!$                 quit_adjust = .false.
!!$              endif
!!$           enddo
       if (maxval(abs(T - T_old)/T_old) .lt. 1.e-4) then
          exit
       endif
       T_old = T
       
    enddo ! do n=1,N_iter

    ! write(*,*) 'before after convection'
    ! write(*,*) 'k, T_before, T_after'
    if (mod(tstep,1000) .eq. 0) then
       do k=1,npz
          write(*,*) k, T(k), grad(k), grad_type(k), mask(k), q(k,1), qsats(k,1), sum(q(k,:))
       enddo
    endif
    ! stop
  end subroutine adjustment
end module eos_dry_adj
