module condense

  use params, only : rdgas, q0, dp, nf
  use phys, only: T_TP => H2O_TriplePointT, P_TP => H2O_TriplePointP, &
       L_vap => H2O_L_vaporization_TriplePoint, Rstar, mu_v => H2O_MolecularWeight, &
       mu_d => H2He_solar_MolecularWeight
  use accurate_l, only : p_sat, L_calc, log_ps_pc, pc, dlogp_dlogt
  use tables, only : phase_grad, lheat, satur, find_var_lin, find_var_loglin
  implicit none

  
  real(dp) :: rvgas = 461.52
  real(dp) :: L_sub = 2.84e6 ! Defined at triple point
  real(dp) :: eps = mu_v/mu_d ! Ratio of water to hydrogen molecular weight

  real(dp) :: p_nought
  integer, parameter :: rk = kind(1.0D+00)
contains

  subroutine newton_iteration(f, dfdx, x1, sol, icode)

    interface
       real function f(x)
         use params, only: dp
         real(dp), intent(in) :: x
       end function f

       real function dfdx(x)
         use params, only: dp
         real(dp), intent(in) :: x
       end function dfdx
    end interface

    real(dp), intent(in) :: x1
    real(dp), intent(out) :: sol
    integer, intent(out) :: icode
    
    integer :: i, imax
    real(dp) :: x_old, x_new, deriv, x_test
    real(dp) :: eps = 1.e-15
    real(dp) :: tol = 1.e-7
    real(dp) :: lim = 1.e-5
    
    icode = 0

    i = 0
    imax = 10000

    x_old = x1
    do while (i .lt. imax)
       deriv = dfdx(x_old)
       if (deriv .lt. eps) then
          write(*,*) 'Bad derivative in Newton Method'
          icode = 1
          exit
       endif
       
       x_new = x_old - f(x_old)/dfdx(x_old)

       if ( abs(x_new/x_old - 1) .lt. tol .and. abs(f(x_new)) .lt. lim) then
          icode = 0
          sol = x_new
          exit
       endif

       i = i + 1
       x_test = x_old
       x_old = x_new
       
    end do

    if (i .eq. imax) then
       write(*,*) 'Reached max iterations without converging'
       icode = 2
       sol = x_new
    endif
    
    
  end subroutine newton_iteration
  
  real function find_root(T)

    real(dp), intent(in) :: T
    
    find_root = log_ps_pc(T) - log(p_nought/pc)
    
  end function find_root
  
  subroutine dew_point_T(pf, Tf)
    real(dp), intent(in) :: pf
    real(dp), intent(out) :: Tf

    !Tf = T_TP/(1 - Rstar*T_TP/mu_v/L_vap*log(pf/P_TP))

    ! New dew_point_T requires inversion of psat(T)
    integer :: info
    real(dp) :: x1

    if (pf .lt. 611.655) then
       Tf = T_TP/(1 - Rstar*T_TP/mu_v/L_vap*log(pf/P_TP))
    else
       
       x1 = T_TP/(1 - Rstar*T_TP/mu_v/L_vap*log(pf/P_TP))
       p_nought = pf
       
       call newton_iteration(find_root, dlogp_dlogt, x1, Tf, info)

       if (info .ne. 0) then
          write(*,*) 'Error in newton raphson iteration'
       endif
       
    endif
    
  end subroutine dew_point_T

  subroutine q_sat(p, T, q)
    real(dp), intent(in) :: p(:), T(:)
    real(dp), intent(out) :: q(:)

    integer :: k
    real(dp) :: psat

    do k=1,size(q)

       call sat_vp(p(k), T(k), psat)
       !psat = p_sat(T(k))
       q(k) = eps*psat/p(k)/(1 + (eps - 1)*psat/p(k))
    enddo
       
  end subroutine q_sat
  
  subroutine rain_out(p,T,q, qsat)
    real(dp), intent(in) :: p, T
    real(dp), intent(inout) :: q
    real(dp), intent(out) :: qsat
    real(dp) :: eps, sat_p

    eps = rdgas/rvgas

    ! Find saturation vapour pressure
    call sat_vp(p, T, sat_p)
    sat_p = p_sat(T)
    qsat = (eps*sat_p/p)/(1. + (eps - 1.)* sat_p/p )
    
    !where (q .gt. qsat)
    !   q = qsat
    !elsewhere
    !   q = q0
    !endwhere
    if (q .gt. qsat) then
       q = qsat
    else
       q = q0
    endif
    
  end subroutine rain_out

  subroutine sat_vp(p, T, sat_p)
    real(dp), intent(in) :: p
    real(dp), intent(in) :: T
    real(dp), intent(out) :: sat_p

    real(dp) :: L

!    if (T .gt. T_tp) then
       
!    else
!      L = L_sub
!    endif

       if (T .gt. 273.16) then
          call find_var_loglin(T, satur, sat_p)
          !sat_p = p_sat(T)
       else
          L = lheat(1)
          sat_p = p_tp*exp(-L/rvgas * (1./T - 1./T_tp) )
          !sat_p = satur(1)
          
       endif
       
       
    
  end subroutine sat_vp

  subroutine cold_trap(q, ktrop, p, T)
    real(dp), intent(inout) :: q(:), T(:)
    real(dp), intent(in) :: p(:)
    integer, intent(inout) :: ktrop

    integer :: i(1), j
    integer :: trop_index
    real(dp) q_min, T_smoothed(size(T)), q_smooth(size(T))

    !write(*,*) 'cold trap', adjust_mask(1)
    ! do j=1,nf
    !    if (adjust_mask(j)) then
    !       trop_index = j-1
    !       exit
    !    endif
    ! enddo
    q_min = 10000.

    do j=2,nf-1
       q_smooth(j) =  q(j-1)*0.25_dp + q(j)*0.5_dp + q(j+1)*0.25_dp
    enddo
    q_smooth(1) = q(1)*0.8 + q(2)*0.2
    q_smooth(nf) = q(nf)*0.8 + q(nf-1)*0.2

    !call q_sat(p, T_smoothed, q_smooth)
    ktrop=1
    do j=nf,1,-1
       if (q(j) .gt. 1) then
          q(j) = 1.
          call dew_point_T(p(j), T(j))
          cycle
       endif
         
       if (q(j) .lt. q_min) then
          q_min = q(j)
       else if (p(j) .lt. 9.*p(nf)/10.) then
          q(1:j) = q_min
          ktrop = j+1
          exit
       endif
       
    enddo
       
    ! i = minloc(q(1:ktrop))
    ! write(*,*) 'MINLOC', i(1), q(i(1))
    ! do j=1, i(1)
    !    q(j) = q(i(1))
    ! enddo
    
  end subroutine cold_trap
end module condense
