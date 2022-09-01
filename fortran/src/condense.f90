module condense

  use params, only : rdgas, q0, dp, nf
  use phys, only: T_TP => H2O_TriplePointT, P_TP => H2O_TriplePointP, &
       L_vap => H2O_L_vaporization_TriplePoint, Rstar, mu_v => H2O_MolecularWeight

  implicit none

  
  real(dp) :: rvgas = 461.52
  real(dp) :: L_sub = 2.84e6 ! Defined at triple point
      real(dp) :: eps = 9. ! Ratio of water to hydrogen molecular weight
  
contains

  subroutine dew_point_T(pf, Tf)
    real(dp), intent(in) :: pf
    real(dp), intent(out) :: Tf

    Tf = T_TP/(1 - Rstar*T_TP/mu_v/L_vap*log(pf/P_TP))
    
  end subroutine dew_point_T

  subroutine q_sat(p, T, q)
    real(dp), intent(in) :: p(:), T(:)
    real(dp), intent(out) :: q(:)

    integer :: k
    real(dp) :: p_sat

    do k=1,size(q)
       call sat_vp(p(k), T(k), p_sat)
       q(k) = eps*p_sat/p(k)/(1 + (eps - 1)*p_sat/p(k))
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
       L = L_vap
!    else
!      L = L_sub
!    endif

    sat_p = p_tp*exp(-L/rvgas * (1./T - 1./T_tp) )
    
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
