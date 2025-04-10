module q_sat_test
  use params, only: dp
  use tables, only: find_var_loglin, satur, lheat
  use atmosphere, only: mmw_dry, q_orig, nqt, nqr
  use accurate_l, only: dlogp_dlogt, pc, log_ps_pc
  use phys, only: H2O, Rstar
  
  implicit none

  real(dp) :: rvgas = 461.52
  real(dp) :: L_sub = 2.84e6 ! Defined at triple point
  !real(dp) :: eps = mu_v/mu_d ! Ratio of water to hydrogen molecular weight

  real(dp) :: p_nought
  integer :: a,b
  
contains
  real function find_root(T)

    real(dp), intent(in) :: T
    
    find_root = log_ps_pc(T) - log(p_nought/pc)
    
  end function find_root
  
  subroutine sat_vp(T, sat_p)
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
          sat_p = H2O%tp_p*exp(-L/rvgas * (1./T - 1./H2O%tp_t) )
          !sat_p = satur(1)
          
       endif
       
       
    
  end subroutine sat_vp
  
    subroutine q_sat(p, T, q)
    real(dp), intent(in) :: p(:), T(:)
    real(dp), intent(out) :: q(:,:)

    integer :: k,n
    real(dp) :: psat, q_old,eps

    do k=1,size(p)

       call sat_vp(T(k), psat)
       !psat = p_sat(T(k))
!q(k) = eps*psat/p(k)/(1 + (eps - 1)*psat/p(k))
       if (psat .gt. p(k) .or. T(k)>H2O%crit_T) then
          q(k,1) = 100
       else
! 1 Always water
          eps = H2O%mmw/mmw_dry(k)
          
          q(k,1) = eps*(psat/p(k))/(1 + (eps-1)*(psat/p(k)))

! Update dry species
          do n=2,nqt
             q(k,n) = (1. - q(k,1))/(1 - q_orig(k,1)) * q_orig(k,n)
          enddo
          
       endif
    enddo
       
  end subroutine q_sat

  
  subroutine dew_point_T(pf, Tf)
    real(dp), intent(in) :: pf
    real(dp), intent(out) :: Tf

    !Tf = T_TP/(1 - Rstar*T_TP/mu_v/L_vap*log(pf/P_TP))

    ! New dew_point_T requires inversion of psat(T)
    integer :: info
    real(dp) :: x1

    if (pf .lt. 611.655) then
       Tf = H2O%tp_t/(1 - Rstar*H2O%tp_t/H2O%mmw/H2O%L_vap_tp*log(pf/H2O%tp_p))
    else
       
       x1 = H2O%tp_t/(1 - Rstar*H2O%tp_t/H2O%mmw/H2O%L_vap_tp*log(pf/H2O%tp_p))
       p_nought = pf
       
       call newton_iteration(find_root, dlogp_dlogt, x1, Tf, info)

       if (info .ne. 0) then
          write(*,*) 'Error in newton raphson iteration'
       endif
       
    endif
    
  end subroutine dew_point_T

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

end module q_sat_test
