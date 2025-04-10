module mat_step_mod
  use utils, only: linear_log_interp, bezier_interp
  use condense, only: set_q
  use atmosphere, only: nqt, nqr, q_orig, get_mmw, get_cp
  use params, only: dp, Nt, nf, ne, const, Finc, Fint, q0, surface, A_s, sb, surf_const, &
       moisture_scheme, sb, grav, cpair, del_time, accelerate, cp_s, rho_s, rdgas, &
       depth, U, C_d, conv_switch, sensible_heat, output_file, turb_diff

  use jacobian, only: calc_matrix, solve_matrix
  implicit none

  integer :: converge_count 
contains
  
  subroutine mat_step(Tf, pf, pe, file_name, q, Ts)
        !====================================================================================
    ! Description
    !====================================================================================
    !
    ! Controls timestepping. Single timestep controlled from single_step
    !
    !====================================================================================
    ! Input variables
    !====================================================================================
    
    real(dp), dimension(:), intent(in) :: pf ! Layer pressure
    real(dp), dimension(:), intent(in) :: pe ! Level pressure

    character(*), intent(in) :: file_name ! Output file name
    
    !====================================================================================
    ! Input/Output variables
    !====================================================================================
    
    real(dp), dimension(:), intent(inout) :: Tf  ! Layer temperatures
    real(dp), dimension(:,:), intent(inout) :: q   ! Specific humidity
    
    real(dp), intent(inout) :: Ts ! Surface temperature

    !====================================================================================
    ! Work variables
    !====================================================================================
    
    integer :: i,j,k
    integer :: ktrop 
    integer :: output_interval
    integer :: print_interval
    
    logical :: stop_switch = .false.   ! Controls when exit due to convergence occurs

    real(dp), dimension(ne) :: net_F  ! Net flux
    real(dp), dimension(ne) :: fup    ! Upwards IR flux
    real(dp), dimension(ne) :: fdn    ! Downwards IR flux
    real(dp), dimension(ne) :: s_up   ! Upwards SW flux
    real(dp), dimension(ne) :: s_dn   ! Downwards SW flux
    real(dp), dimension(ne) :: Te     ! Edge of layer temperature
    real(dp), dimension(ne) :: turb_flux

    real(dp), dimension(nf) :: flux_diff
    real(dp), dimension(nf) :: dflux ! Local convergence criterion
    real(dp), dimension(nf) :: dT    ! Change in temperature
    real(dp), dimension(nf) :: delp
    real(dp), dimension(nf) :: factor
    real(dp), dimension(nf) :: T_old
    real(dp) :: factor_surf, t_surf_old, dT_surf

    integer, dimension(nf) :: dry_mask
    
    real(dp) ::  sens ! Sensible heat flux
    real(dp) :: olr
    real(dp) :: dflux_surf
    character(len=300) :: plot_command
    real(dp) :: hrt(nf), t_consts(nf), delt(nf)
    integer :: tstep
    real(dp) :: jac(nf,nf), norm
    integer :: n
    logical :: stop_signal
    
    do i=1,nf
       delp(i) = pe(i+1) - pe(i)
    enddo
    call interp_to_edges(pf, pe, Tf, Te)

    call set_q(pf, Tf, q, ktrop, 1)
    call interp_to_edges(pf, pe, Tf, Te)

    stop_signal = .false.
    do n=1,1000
       ktrop = 10
       call calc_matrix(Tf, pf, pe, delp, 1.0_dp, Finc, Fint, olr, net_F, &
            Te, q, Ts, hrt, n, t_consts, jac, ktrop)

       call solve_matrix(jac, -hrt, delt)

       if (maxval(abs(delt)) .gt. 10.0) then
          norm = 0.0
          do k=1,nf
             norm = norm + delt(k)**2.
          enddo
          delt = delt/sqrt(norm) * 10.00
       endif
       Tf = Tf + delt
       !call set_q(pf, Tf, q, ktrop, 1)
       call interp_to_edges(pf, pe, Tf, Te)

       do i=1,nf
          write(*,*) i, Tf(i), delt(i), hrt(i)/sb/Tf(i)**4., q(k,1)
       enddo

       call check_convergence(Tf, delt, hrt, stop_signal)

       if (stop_signal) exit
    enddo

    write(*,*) 'CONVERGED'
    do k=1,nf
       write(*,*) k, pf(k), Tf(k), q(k,1)
    enddo
  end subroutine mat_step

  subroutine check_convergence(T, dT, dflux, stop_signal)
    real(dp), intent(in) :: T(:), dT(:), dflux(:)
    logical, intent(inout) :: stop_signal

    write(*,*) 'checking convergence', maxval(abs(dT)), maxval(abs(dflux/sb/T**4.))
    if (maxval(abs(dT)) .lt. 1.e-6 .and. maxval(abs(dflux/sb/T**4.)) .lt. 1.e-6)  then
       stop_signal = .True.
    endif

  end subroutine check_convergence
  subroutine interp_to_edges(pf, pe, Tf, Te)
    !====================================================================================
    ! Description
    !====================================================================================
    !
    ! Interpolate temperature to edge of layers
    !
    !====================================================================================
    ! Input Variables
    !====================================================================================
    
    real(dp), intent(in) :: pf(nf) ! Mid layer pressure
    real(dp), intent(in) :: pe(ne) ! Edge level pressure
    real(dp), intent(in) :: Tf(nf) ! Mid layer temperature
    
    !====================================================================================
    ! Output Variables
    !====================================================================================

    real(dp), intent(out) :: Te(ne) ! Edge of layer temperature

    !====================================================================================
    ! Work variables
    !====================================================================================

    integer :: i

    !====================================================================================
    ! Main body
    !====================================================================================

    do i = 2, nf-1
       call bezier_interp(pf(i-1:i+1), Tf(i-1:i+1), 3, pe(i), Te(i))
    end do

    call bezier_interp(pf(nf-2:nf), Tf(nf-2:nf), 3, pe(nf), Te(nf))
    call linear_log_interp(pe(1), pf(1), pf(2), Tf(1), Tf(2), Te(1))
    call linear_log_interp(pe(ne), pf(nf-1), pf(nf), Tf(nf-1), Tf(nf), Te(ne) )

  end subroutine interp_to_edges
  
end module mat_step_mod
