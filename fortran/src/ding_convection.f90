module ding_convection
  use adjust_mod, only: gradient
  use phys, only : T_TP => H2O_TriplePointT, P_TP => H2O_TriplePointP, L_sub => H2O_L_sublimation, &
       L_vap => H2O_L_vaporization_TriplePoint, CP_v => H2O_cp, CP_d => H2He_solar_cp, &
       mu_d => H2He_solar_MolecularWeight, mu_v => H2O_MolecularWeight, Rstar
  use condense, only: q_sat, q_sat_single, r_sat_single
  use tables, only : find_var_lin, lheat, hv, cpv_new, cpl_new, Ttp, find_var_simplified
  use params, only: dp, grav
  implicit none

  real(dp) :: precision = 1.e-7
  integer :: n_iter_conv = 1000
  integer :: n_iter_newt = 100
contains

  subroutine adjustment(p, delp, T, q, qc, ktrop, grad, olr, mask, tstep)
    !==========================================================================
    ! Description
    !==========================================================================
    ! Performs moist adiabatic adjustment pairwise, conserving non-dilute moist
    ! enthalpy using Newton iteration, as in Ding & Pierrehumbert 2016

    !==========================================================================
    ! Input variables
    !==========================================================================    
    real(dp), intent(in)   , dimension(:) :: p! Pressure
    real(dp), intent(inout), dimension(:) :: grad,delp ! Pressure thickness can and will change
    integer, intent(in) :: tstep

    real(dp), intent(in) :: olr
    !==========================================================================
    ! Output variables
    !==========================================================================
    integer, intent(inout) :: ktrop
    logical, intent(out) :: mask(:)
    
    !==========================================================================
    ! Mixed input/output
    !==========================================================================
    real(dp), intent(inout), dimension(:) :: T,q ! Temperature and specific humid.
    real(dp), intent(inout), dimension(:) :: qc ! Condensate ratio

    !==========================================================================
    ! Local variables
    !==========================================================================
    real(dp) :: dlnTdlnp, temp, T_lift,rho_lift,q_lift, pfact, cp_v_loc,cp_c_loc
    real(dp) :: qsats(size(p)), rhos(size(p)), rc(size(p)), rv(size(p)), qt(size(p))
    real(dp) :: R_local(size(p)), cp_local(size(p))
    
    integer :: k, npz ,n
    
    !==========================================================================
    ! Main body
    !==========================================================================

    npz = size(p)
    qt = qc + q
    do n=1,N_iter_conv

       ! Find saturation vapour specific humidity
       call q_sat(p,T,qt,qsats)
       do k=1,npz
          if (T(k) .gt. Ttp) then
             call find_var_simplified(T(k), 'cp_v', cp_v_loc)
             call find_var_simplified(T(k), 'cp_l', cp_v_loc)
          else
             cp_v_loc = cpv_new(1)
             cp_c_loc = cpv_new(1)
          endif
          cp_local = cp_d*(1-q-qc) * cp_v_loc*q + cp_c_loc*qc
       enddo
       
       R_local = Rstar*(q/mu_v/(1-qc) + (1-q)/(1-qc)/mu_d)
       !cp_local = cp_d*(1-q) * cp_v * q
       rhos = p/R_local/T

       do k=npz-1,1,-1

          ! If saturated, gradient should be moist gradient
          if (q(k+1) .gt. 0.90*qsats(k+1)) then

             ! get moist gradient
             call gradient(p(k+1), T(k+1), dlnTdlnp,temp)

             ! Calculate lifted temperature from adiabatic expansion
             pfact = exp(dlnTdlnp*log(p(k)/p(k+1)))
             T_lift = T(k+1)*pfact
             
             ! Calculate q here as proportion of gas phase rv/(1+rv)
             call r_sat_single(p(k), T_lift, q_lift)
             q_lift = q_lift/(1 + q_lift)

             rho_lift = p(k)/T_lift/( Rstar*(q_lift/mu_v + (1-q_lift)/mu_d ))

             if (rho_lift .lt. rhos(k)) then
                ! Convection here! Conserve non-dilute moist enthalpy
                call conserve_moist_enthalpy(p(k), p(k+1), T(k), T(k+1), q(k), q(k+1), &
                                             delp(k), delp(k+1), qc(k), qc(k+1))
                mask(k) = .true.
                mask(k+1) = .true.
             endif
             
          else
             ! Do dry expansion with value of R_local
             dlnTdlnp = R_local(k+1)/cp_local(k+1)

             pfact = exp(dlnTdlnp*log(p(k)/p(k+1)))
             T_lift = T(k+1)*pfact

             rho_lift = p(k)/T_lift/R_local(k+1)

             if (rho_lift .lt. rhos(k)) then
                call conserve_dry_enthalpy(p(k), p(k+1), T(k), T(k+1), q(k), q(k+1), &
                                           delp(k), delp(k+1), qc(k), qc(k+1))
                mask(k) = .true.
                mask(k+1) = .true.
             endif
          endif ! if (q(k+1) .gt. qsats(k+1)) 
          
          
       enddo ! do k=npz-1,1,-1
    enddo ! do n=1,N_iter
    
  end subroutine adjustment

  subroutine conserve_dry_enthalpy(p1, p2, T1, T2, q1, q2, dp1, dp2, qc1, qc2)
    !==========================================================================
    ! Description
    !==========================================================================
    ! Calculates adjusted state with a Newton iteration by conserving dry enthalpy
    ! in the absence of condensation on adiabatic lifting

    !==========================================================================
    ! Input variables
    !==========================================================================    
    real(dp), intent(in) :: p1, p2 ! Pressure of upper/lower layer
    real(dp), intent(inout) :: T1, T2 ! Temperature of upper/lower layer
    real(dp), intent(inout) :: q1, q2,qc1,qc2 ! Specific humidity of lower/upper layer
    real(dp), intent(in)    :: dp1, dp2 ! Pressure thicknesses of the layers
    
    !==========================================================================
    ! Local variables
    !==========================================================================    
    real(dp) :: k_old, L1, L2, T2_guess, T1_new, dlnTdlnp, pfact,k_new, k1, k2

    real(dp) :: q1_new, q2_new, R_new, cpv_1, cpv_2, dT,k_diff, qc1_new, qc2_new

    integer :: n_iter_newt, n
    
    !==========================================================================
    ! Main body
    !==========================================================================    

    ! Enthalpy of the old state
    call calc_enthalpy(T1, q1, qc1,k1)
    call calc_enthalpy(T2, q2, qc2,k2)

    k_old = k1*dp1 + k2*dp2
    
    ! Calculate constant q of the after state by conserving mass
    q1_new = (q1*dp1 + q2*dp2)/(dp1 + dp2)
    q2_new = q1_new

    qc1_new = (qc1*dp1 + qc2*dp2)/(dp1+dp2)
    qc2_new = qc1_new
    
    T2_guess = T2
    do while ((abs(dT) .gt. precision) .and. (n .lt. n_iter_newt))
       ! Get R/cp of the lower layer
       
       call k_adj_dry(T2_guess, p1, p2, dp1, dp2, q1_new,  qc1_new, T1_new, k_new)
       call dk_adj_dry(T2_guess, p1,p2,dp1,dp2,q1_new, qc1_new, k_diff)

       dT = (k_new - k_old)/k_diff
       
       T2_guess = T2_guess  - dT
       n = n+1
    enddo

    call k_adj_dry(T2_guess, p1, p2, dp1, dp2, q1_new, q2_new, T1_new, k_new)
    
    q1 = q1_new
    q2 = q2_new
    T2 = T2_guess
    T1 = T1_new
  end subroutine conserve_dry_enthalpy

  subroutine k_adj_dry(T2_guess, p1, p2, dp1, dp2, q1, qc1, T1_new, k_new)
    real(dp), intent(in) :: T2_guess, p1, p2, dp1, dp2, q1, qc1
    real(dp), intent(out) :: T1_new, k_new

    real(dp) :: R_new, cpv_2,cpc_2, k1, k2, pfact,cp_new
    
    R_new = Rstar*(q1/(1-qc1)/mu_v + (1-q1)/(1-qc1)/mu_d)

    if (T2_guess .lt. Ttp) then
       cpv_2 = cpv_new(1)
       cpc_2 = cpv_new(1)
    else
       call find_var_simplified(T2_guess, 'cpv', cpv_2)
       call find_var_simplified(T2_guess, 'cpl', cpc_2)
    endif

    cp_new = q1*cpv_2 + qc1*cpc_2 + (1 - q1 - qc1)*cp_d
    pfact = exp(R_new/cpv_2 * log(p1/p2))
    T1_new = T2_guess*pfact

    call calc_enthalpy(T1_new, q1,qc1, k1)
    call calc_enthalpy(T2_guess, q1,qc1, k2)

    k_new = k1*dp1 + k2*dp2

  end subroutine k_adj_dry

  subroutine dk_adj_dry(T2_guess, p1, p2, dp1, dp2,q1, q2, dk_adj)
    real(dp), intent(in) :: T2_guess, p1, p2, dp1, dp2,q1,q2
    real(dp), intent(out) :: dk_adj

    real(dp) :: t, t2,  kp, km
    real(dp) :: eps = 1.e-8

    call k_adj_dry(T2_guess+eps/2, p1, p2, dp1, dp2, q1, q2, t, kp)
    call k_adj_dry(T2_guess-eps/2, p1, p2, dp1, dp2, q1, q2, t, km)

    dk_adj = (kp - km)/eps
  end subroutine dk_adj_dry
  
  subroutine conserve_moist_enthalpy(p1, p2, T1, T2, q1, q2, dp1,dp2, qc_1, qc_2)
    !==========================================================================
    ! Description
    !==========================================================================
    ! Calculates adjusted state with a Newton iteration by conserving moist enthalpy

    !==========================================================================
    ! Input variables
    !==========================================================================    
    real(dp), intent(in) :: p1, p2 ! Pressure of upper/lower layer
    real(dp), intent(inout) :: T1, T2 ! Temperature of upper/lower layer
    real(dp), intent(inout) :: q1, q2 ! Specific humidity of lower/upper layer
    real(dp), intent(in)    :: dp1, dp2 ! Pressure thicknesses of the layers
    real(dp), intent(out)   :: qc_1, qc_2 ! Vapour and condensate in adjusted state
    
    !==========================================================================
    ! Local variables
    !==========================================================================    
    real(dp) :: k_old, L1, L2, T2_guess, T1_guess, dlnTdlnp, pfact, q1_new, q2_new, k_new

    real(dp) :: eta, r1_new, r2_new, rc_2, rc_1, m_i, f1, k_diff, dT, T1_new, rv_1, rv_2
    real(dp) :: k1_old, k2_old
    integer :: n

    ! Calculate moist enthalpy of old state
    call calc_enthalpy(T1, q1, 0._dp, k1_old)
    call calc_enthalpy(T2, q2, 0._dp, k2_old)

    k_old = k1_old*dp1 + k2_old*dp2
    
    ! Do Newton iteration to find the appropriate T2 in the new state
    T2_guess = T2
    do while ((abs(dT) .gt. precision) .and. (n .lt. n_iter_newt)) 
       call k_adj(T2_guess, p1, p2, dp1, dp2,q1,q2, q1_new, q2_new, qc_1, qc_2, T1_new, k_new)
       call dk_adj(T2_guess, p1,p2,dp1,dp2, q1, q2, k_diff)

       dT = (k_new - k_old)/k_diff
       
       T2_guess = T2_guess  - dT
       n = n+1
    enddo

    call k_adj(T2_guess, p1, p2, dp1, dp2, q1,q2,q1_new, q2_new, qc_1, qc_2, T1_new, k_new)
    q1 = q1_new
    q2 = q2_new
    T2 = T2_guess
    T1 = T1_new

  end subroutine conserve_moist_enthalpy

  subroutine k_adj(T2_guess, p1, p2, dp1, dp2, q1_0,q2_0,qv_1,qv_2,qc_1, qc_2,T1_new,k_new)
    real(dp), intent(in):: T2_guess, p1, p2, dp1, dp2, q1_0,q2_0
    real(dp), intent(out) ::  qc_1, qc_2, k_new, T1_new,qv_1, qv_2


    real(dp) :: dlnTdlnp, T1_guess, rc_1, rc_2, rv_1, rv_2
    real(dp) :: k1, k2, eta, f1, m_i, pfact
    
    call gradient(p2, T2_guess, dlnTdlnp)
    
    pfact = exp(dlnTdlnp*log(p1/p2))
    T1_new = T2_guess*pfact

    call r_sat_single(p2, T2_guess,rv_2)
    call r_sat_single(p1, T1_new,rv_1)
    
    ! Ratio of upper condensate to lower condensate mixing ratio
    eta = (1 + rv_1)/(1+ rv_2)

    ! Initial mass of vapour
    m_i = q1_0*dp1 + q2_0*dp2

    ! Calculate rc_2 according to Ding and Pierrehumbert (note our lower + upper layer)
    ! indices are swapped wrt this paper
    rc_2 = (m_i*(1 + rv_2) - dp2*rv_2 - dp1*rv_1/eta)/(dp1 + dp2 - m_i)

    if (rc_2 .lt. 0) then
       ! Condensate is negative (not enough initial mass of water)
       qv_1 = rv_1/(1 + rv_1)
       qv_2 = rv_2/(1 + rv_2)

       f1 = m_i / (qv_2*dp2 + qv_1*dp1)

       qv_1 = qv_1 * f1
       qv_2 = qv_2 * f1

       qc_1 = 0._dp
       qc_2 = 0._dp
       
    else
       ! Condensate is positive
       rc_1 = rc_2 * eta

       qv_1 = rv_1/(1+rv_1 + rc_1)
       qc_1 = rc_1/(1+rv_1 + rc_1)

       qv_2 = rv_2/(1+rv_2 + rc_2)
       qc_2 = rc_2/(1+rv_2 + rc_2)
    endif


    ! Find enthalpy of the new structure
    ! Dry component = cp*T, moist component = tabulated

    call calc_enthalpy(T1_guess, qv_1, qc_1, k1)
    call calc_enthalpy(T2_guess, qv_2, qc_2, k2)

    k_new = k1*dp1 + k2*dp2
        
  end subroutine k_adj

  subroutine dk_adj(T2_guess, p1, p2, dp1, dp2, q1_0, q2_0, dk_a)
    real(dp), intent(in) :: T2_guess, p1, p2, dp1, dp2, q1_0, q2_0
    real(dp), intent(out) :: dk_a

    real(dp) :: qc_1, qc_2, qv_1, qv_2, k_out_up, k_out_dn, dqsat_dt, T_temp

    real(dp) :: eps = 1.e-8
    ! Estimate derivative of dk_adj wrt T2_guess

    call k_adj(T2_guess + eps/2, p1,p2,dp1,dp2,q1_0, q2_0, qv_1,qv_2,qc_1,qc_2,T_temp,k_out_up)
    call k_adj(T2_guess - eps/2, p1,p2,dp1,dp2,q1_0, q2_0, qv_1,qv_2,qc_1,qc_2,T_temp,k_out_dn)

    dk_a = (k_out_up - k_out_dn)/eps
  end subroutine dk_adj
    
  subroutine large_scale_cond(p, T, q, qc)
    real(dp), intent(in) :: p(:)
    real(dp), intent(inout) :: T(:), q(:), qc(:)

    ! Local variables
    integer:: k, n
    integer :: npz

    real(dp) :: qsat, psat, cp_local, dT , qt, q_guess, T_guess
    real(dp) :: kold, knew, rv, Tgp, qgp, Tgm,qgm, kp, km, diff, qc_guess


    npz = size(p)
    do k=1,npz
       qt = qc(k) + q(k)
       call q_sat_single(p(k), T(k), qt,qsat)
       
       if (q(k) .gt. qsat) then
          
          ! Do in terms of the enthalpy of beginning vs final state
          call calc_enthalpy(T(k), q(k), qc(k), kold)
          dT = 0.0_dp
          T_guess = T(k)
          q_guess = qsat
          do while ( (abs(dT) .gt. precision) .and. (n .lt. n_iter_newt))
             call q_sat_single(p(k), T_guess,qt, q_guess)
             qc_guess = qt - q_guess
             call calc_enthalpy(T_guess, q_guess, qc_guess, knew)
             
             Tgp = T_guess+1.e-8
             call q_sat_single(p(k), Tgp, qt,qgp)
             qc_guess = qt - qgp
             call calc_enthalpy(Tgp, qgp,qc_guess, kp)

             Tgm = T_guess - 1.e-8
             call q_sat_single(p(k), Tgm, qt, qgm)
             qc_guess = qt - qgm
             call calc_enthalpy(Tgm, qgm,qc_guess, km)

             diff = (kp-km)/(2.e-8)
             dT = (kold - knew)/diff
             
             T_guess = T_guess + dT
             n = n+1
          enddo
          
       
!!$       if (rv(k) .gt. rsat) then
!!$          ! Find adjustment that will bring atmosphere to saturation
!!$          cp_local = cp_l*(rv(k)+rc(k))/(1 + rv(k) + rc(k)) + cp_d/(1 + rv(k) + rc(k))
!!$
!!$          dT = 0.0_dp
!!$          n_newt = 100
!!$          T_guess = T(k)
!!$          q0 = rv(k) / (1 + rv(k) + rc(k))
!!$          do while ( (abs(dT) .gt. precision) .and. (n .lt. n_iter))
!!$             call r_sat_single(p(k), T_guess, rsat, psat)
!!$             
!!$             qsat_guess = rsat/(1 + rv(k) + rc(k)) ! Note rv+rc = const throughout to conserve mass
!!$
!!$             dqsat_dt = rsat*L*mu_v/Rstar/T_guess**2 * p(k)/(p(k)-psat)
!!$             dT = (L/cp_local * (q0 - qsat_guess) + (T_guess - T(k))) &
!!$                  / (1 + L/cp_local *dqsat_dt)
!!$
!!$             T_guess = T_guess + dT
!!$          enddo
!!$
          ! Update T(k) with the value from large scale condensation
          T(k) = T_guess
          call q_sat_single(p(k), T_guess,qt, q_guess)
          qc(k) = qt - q_guess
          q(k) = q_guess
       endif
       
    enddo
  end subroutine large_scale_cond

  subroutine calc_enthalpy_dry(T, q_v, q_c,kout)
    real(dp), intent(in)  :: T, q_v,q_c ! Temp, vapour and condensate ratio
    real(dp), intent(out) :: kout ! output enthalpy

    real(dp) :: h_d, h_v,h_c

    h_d = cp_d*T
    if (T .lt. Ttp) then
       h_v = hv(1) + cpv_new(1)*(T - Ttp)
       h_c = h_v - L_sub
    else
       call find_var_simplified(T, 'hv', h_v)
       call find_var_simplified(T, 'hl', h_c)
    endif

    kout = h_d*(1 - q_v - q_c) + h_v*q_v + h_c*q_c
  end subroutine calc_enthalpy_dry
  
  subroutine calc_enthalpy(T, q_v, q_c, kout)
    real(dp), intent(in)  :: T, q_v, q_c ! Temp, vapour and condensate ratio
    real(dp), intent(out) :: kout ! output enthalpy

    real(dp) :: h_d, h_v, h_c

    h_d = cp_d*T

    if (T .lt. Ttp) then
       h_v = hv(1) + cpv_new(1)*(T - Ttp)
       h_c = h_v - L_sub
    else
       call find_var_simplified(T, 'hv', h_v)
       call find_var_simplified(T, 'hl', h_c)
    endif

    ! Per unit mass, multiply by dp out of routine!
    kout = h_d*(1 - q_v - q_c) + h_c*q_c + h_v*q_v
  end subroutine calc_enthalpy
  
  subroutine rain_out(delp, T, q, qc, h_cond)
    !==========================================================================
    ! Description
    !==========================================================================
    ! Does rain out and conserves mass + enthalpy in doing so

    !==========================================================================
    ! Input variables
    !==========================================================================    
    real(dp), intent(in),dimension(:) :: T

    !==========================================================================
    ! Output variables
    !==========================================================================    
    real(dp), intent(out) :: h_cond

    !==========================================================================
    ! Input/Output variables
    !==========================================================================    
    real(dp), intent(inout) :: delp(:), q(:) ,qc(:)
    
    !==========================================================================
    ! Local variables
    !==========================================================================    
    integer :: k
    integer :: npz 

    real(dp) :: m_cond_tot, m_cond, rsat, T_guess, q0, rv
    real(dp) :: required, temp
    
    h_cond = 0.0_dp
    m_cond_tot = 0.0_dp
    npz = size(T)
    do k=1,npz
       
       ! Mass loss
       m_cond     = qc(k) * delp(k)/grav
       m_cond_tot = m_cond_tot + m_cond

       !=================================================================================
       ! Development notes
       !=================================================================================
       ! Re-evaporated condensate
       ! To do this in an enthalpy-conserving way, one would have to heat the condensate
       ! from source level to destination level, and then evaporate it into the source layer.
       ! This gives two equations in two variables -- mass of fallen condensate
       ! and the final temperature of layer below. Can be solved with a 2D newton iteration
       ! Easiest to do by going down layer by layer and testing for undersaturation, then
       ! selecting available condensate from layers above to transport to that layer and
       ! re-evaporate. What layer to select from though? Nearest layer, highest layer, layer
       ! with the highest mass of condensate? Probably best to do nearest layer since higher
       ! pressure thickness in log space will make it more likely to have enough available
       ! condensate. 
       !================================================================================

       if (T(k) .lt. Ttp) then
          h_cond = h_cond + m_cond * (cpv_new(1)*(T(k)-Ttp) - L_sub)
       else
          call find_var_simplified(T(k), 'hl', temp)
          h_cond = h_cond + m_cond * temp
       endif
       
       ! New q once condensate has rained out
       q(k) = q(k)/(1 - qc(k))

       delp(k) = delp(k) - m_cond*grav

       qc(k) = 0.0_dp
    end do

  end subroutine rain_out
  
end module ding_convection
