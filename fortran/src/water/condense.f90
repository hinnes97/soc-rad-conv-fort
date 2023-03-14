module condense

  use params, only : rdgas, q0, dp, nf, moisture_scheme
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
  integer :: a,b
  
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

       call sat_vp(T(k), psat)
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
    call sat_vp(T, sat_p)
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
          sat_p = p_tp*exp(-L/rvgas * (1./T - 1./T_tp) )
          !sat_p = satur(1)
          
       endif
       
       
    
  end subroutine sat_vp

  subroutine cold_trap(q, ktrop, p, T)
    real(dp), intent(inout) :: q(:), T(:)
    real(dp), intent(in) :: p(:)
    integer, intent(inout) :: ktrop

    integer ::  j
    real(dp) q_min


    q_min = 10000.

    ktrop=1
    do j=nf,1,-1

       if (moisture_scheme=='surface') then
          
          if (q(j) .gt. 1) then
             q(j) = 1.
             call dew_point_T(p(j), T(j))
             cycle
          endif


       
          if (q(j) .lt. q_min) then
             q_min = q(j)
       
          else if (p(j) .lt. 9.*p(nf)/10.) then
             ! Extra condition stops cold trap being too low down
             ! Liable to crashing if set too near the surface
             q(1:j) = q_min
             ktrop = j+1
             exit
          endif

       else if (moisture_scheme=='deep') then
          if (q(j) .lt. q_min) then
             q_min = q(j)

          else if (p(j) .lt. 9*p(nf)/10. .and. q(j) .lt. q0) then
             q(1:j) = q_min
             ktrop = j+1
          endif
          
       endif
    enddo
 
    
 
       
    ! i = minloc(q(1:ktrop))
    ! write(*,*) 'MINLOC', i(1), q(i(1))
    ! do j=1, i(1)
    !    q(j) = q(i(1))
    ! enddo
    
  end subroutine cold_trap

    subroutine set_q(p,T, q,ktrop)
    !==========================================================================
    ! Description
    !==========================================================================    
    ! Sets water vapour mixing ratio to saturation value, and then calculates
    ! the cold trap
    !
    ! moisture_scheme can currently have two values:
    ! - "surface" - assumes surface water ocean -- atmosphere must be saturated
    !               with water vapour and temperature limited by temp when q=1.
    !
    ! - "deep" - deep water content of q0, q is min of saturation or q0
    ! 
    !==========================================================================
    
    !==========================================================================
    ! Input variables
    !==========================================================================    
      real(dp), intent(in)   , dimension(:) :: p! Pressure, pressure thickness
      
    !==========================================================================
    ! Output variables
    !==========================================================================

    integer, intent(out) :: ktrop ! Tropopause defined by cold trap, don't do
                                  ! Moist convective adjustment above this index
    
    !==========================================================================
    ! Mixed input/output
    !==========================================================================
    real(dp), intent(inout), dimension(:) :: T,q ! Temperature and specific humid.
    

    integer :: k, npz    
    !==========================================================================
    ! Main body
    !==========================================================================

    ! Set to saturation
    call q_sat(p, T, q)
    call cold_trap(q, ktrop, p, T)

    
    ! Correct if above deep value
    if (moisture_scheme =='deep') then
       npz = size(q)
       do k=1,npz
          q(k) = min(q(k), q0)
       enddo
       ! Don't do moist convection if entire atmosphere is at q0
       if (q(ktrop) .gt. q0) ktrop = npz
    endif

   end subroutine set_q


end module condense
