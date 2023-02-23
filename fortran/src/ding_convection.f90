module ding_convection
  use adjust, only: gradient
  use phys, only : T_TP => H2O_TriplePointT, P_TP => H2O_TriplePointP, L_sub => H2O_L_sublimation, &
       L_vap => H2O_L_vaporization_TriplePoint, CP_v => H2O_cp, CP_d => H2He_solar_cp, &
       mu_d => H2He_solar_MolecularWeight, mu_v => H2O_MolecularWeight, Rstar
  use condense, only: q_sat, q_sat_single, r_sat_single
  use tables, only : find_var_lin, lheat
  implicit none

  real(dp) :: precision = 1.e-7
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
    real(dp) :: dlnTdlnp, temp, T_lift,rho_lift,q_lift, p_fact, cp_v_loc
    real(dp) :: qsats(size(p)), rhos(size(p)), rc(size(p)), rv(size(p))
    real(dp) :: R_local(size(p)), cp_local(size(p))
    
    integer :: k, npz 
    
    !==========================================================================
    ! Main body
    !==========================================================================

    npz = size(p)

    do n=1,N_iter

       ! Find saturation vapour specific humidity
       call q_sat(p,T,qsats)
       do k=1,npz
          call find_var_simplified(T(k), 'cp_v', cp_v_loc)
          cp_local(k) = cp_d*(1-q) * cp_v_loc*q
       enddo
       
       R_local = Rstar*(q/mu_v + (1-q)/mu_d)
       !cp_local = cp_d*(1-q) * cp_v * q
       rhos = p/R_local/T

       do k=npz-1,1,-1

          ! If saturated, gradient should be moist gradient
          if (q(k+1) .gt. qsats(k+1)) then

             ! get moist gradient
             call gradient(p(k+1), T(k+1), dlnTdlnp,temp)

             ! Calculate lifted temperature from adiabatic expansion
             pfact = exp(dlnTdlnp*log(p(k)/p(k+1))
             T_lift = T(k+1)*pfact
             
             ! Calculate q here
             call q_sat_single(p(k), T_lift, q_lift)

             rho_lift = p(k)/T_lift/( Rstar*(q_lift/mu_v + (1-q_lift)/mu_d )

             if (rho_lift .lt. rhos(k)) then
                ! Convection here! Conserve non-dilute moist enthalpy
                call conserve_moist_enthalpy(p(k), p(k+1), T(k), T(k+1), q(k), q(k+1), &
                                             delp(k), delp(k+1), rv(k), rv(k+1), rc(k), rc(k+1))
                mask(k) = .true.
                mask(k+1) = .true.
             endif
             
          else
             ! Do dry expansion with value of R_local
             dlnTdlnp = R_local(k+1)/cp_local(k+1)

             pfact = exp(dlnTdlnp*log(p(k)/p(k+1))
             T_lift = T(k+1)*pfact

             rho_lift = p(k)/T_lift/R_local(k+1)

             if (rho_lift .lt. rhos(k)) then
                call conserve_dry_enthalpy(p(k), p(k+1), T(k), T(k+1), q(k), q(k+1), &
                                           delp(k), delp(k+1))
                mask(k) = .true.
                mask(k+1) = .true.

                rv(k) = q(k)/(1 - q(k))
                rv(k+1) = q(k+1)/(1 - q(k+1))
                rc(k) = 0._dp
                rc(k+1) = 0._dp
             endif
          endif ! if (q(k+1) .gt. qsats(k+1)) 
          
          
       enddo ! do k=npz-1,1,-1
    enddo ! do n=1,N_iter
    
  end subroutine adjustment

  subroutine conserve_dry_enthalpy(p1, p2, T1, T2, q1, q2, dp1, dp2)
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
    real(dp), intent(inout) :: q1, q2 ! Specific humidity of lower/upper layer
    real(dp), intent(in)    :: dp1, dp2 ! Pressure thicknesses of the layers
    
    !==========================================================================
    ! Local variables
    !==========================================================================    
    real(dp) :: k_old, L1, L2, T2_guess, T1_guess, dlnTdlnp, pfact,k_new

    real(dp) :: q1_new, q2_new, R_new

    integer :: n_newton, n
    
    !==========================================================================
    ! Main body
    !==========================================================================    

    ! Enthalpy of the old state
    cp_1 = q1*cp_v + (1 - q1)*cp_d
    cp_2 = q2*cp_v + (1 - q2)*cp_d

    k_old = cp_1*T1*dp1 + cp_2*T2*dp_2
    
    ! Calculate constant q of the after state by conserving mass
    q1_new = (q1*dp1 + q2*dp2)/(dp1 + dp2)
    q2_new = q1_new

    cp_1 = q1_new*cp_v + (1- q1_new)*cp_d
    cp_2 = cp_1
    R_new = R_star*(q1_new/mu_v + (1-q1_new)/mu_d)

    pfact = exp(R_new/cp_1 * log(p1/p2))

    T2 = k_old/cp_1 /(pfact*dp_1 + dp2)
    T1 = T2_guess*pfact
    q1 = q1_new
    q2 = q2_new
    
  end subroutine conserve_dry_enthalpy
  
  subroutine conserve_moist_enthalpy(p1, p2, T1, T2, q1, q2, dp1,dp2, rv_1, rv_2, rc_1, rc_2)
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
    real(dp), intent(out)   :: rv_1, rv_2, rc_1, rc_2 ! Vapour and condensate in adjusted state
    
    !==========================================================================
    ! Local variables
    !==========================================================================    
    real(dp) :: k_old, L1, L2, T2_guess, T1_guess, dlnTdlnp, pfact, q1_new, q2_new, k_new

    real(dp) :: eta, r1_new, r2_new, rc_2, rc_1, m_i, f1, k_diff, dT, T1_new

    integer :: n_newton, n
    ! Calculate moist enthalpy of old state
    cp_1 = c_l*q1 + cp_d*(1-q1)
    cp_2 = c_l*q2 + cp_d*(1-q2)

    
    if (T1 .gt. 273.16) then
       call find_var_lin(T1, lheat, L1)
    else
       L1 = lheat(1)
    endif
    
    if (T2 .gt. 273.16) then
       call find_var_lin(T2, lheat, L2)
    else
       L2 = lheat(1)
    endif
    
    k_old = cp_1*T1*dp1 + cp_2*T2*dp2 + L1*q1 + L2*q2

    ! Do Newton iteration to find the appropriate T2 in the new state
    T2_guess = T2
    do while ((abs(dT) .gt. precision) .and. (n .lt. n_newton)) 
       call k_adj(T2_guess, p1, p2, dp1, dp2, q1, q2, rv_1, rv_2,rc_1, rc_2, T1_new, k_new)
       call dk_adj(T2_guess, p1,p2,dp1,dp2, k_diff)

       dT = (k_new - k_old)/k_diff
       
       T2_guess = T2_guess  - dT
       n = n+1
    enddo
    
    T2 = T2_guess
    T1 = T1_new

  end subroutine conserve_moist_enthalpy

  subroutine k_adj(T2_guess, p1, p2, dp1, dp2, q1,q2, rv_1, rv_2, rc_1, rc_2,T1_new,k_new)
    real(dp), intent(in) T2_guess, p1, p2, q1, q2, dp1, dp2
    
    real(dp), intent(out) :: rc_1, rc_2, rv_1, rv_2, k_new, T1_new


    real(dp) :: dlnTdlnp, T1_guess, q1_new, q2_new
    
    call gradient(p2, T2, dlnTdlnp)
    
    pfact = exp(dlnTdlnp*log(p1/p2))
    T1_new = T2_guess*pfact

    call r_sat_single(p2,T2_guess,r2_new)
    call r_sat_single(p1,T1_new,r1_new)

    ! Ratio of upper condensate to lower condensate mixing ratio
    eta = (1 + rv_1)/(1+ rv_2)

    ! Initial mass of vapour
    m_i = q1*dp1 + q2*dp2

    ! Calculate rc_2 according to Ding and Pierrehumbert (note our lower + upper layer)
    ! indices are swapped wrt this paper
    rc_2 = (mi*(1 + rv_2) - dp2*rv_2 - dp1*rv_1/eta)/(dp1 + dp2 - m_i)

    if (rc_2 .lt. 0) then
       ! Condensate is negative (not enough initial mass of water)
       q1_new = rv_1/(1 + rv_1)
       q2_new = rv_2/(1 + rv_2)

       f1 = m_i / (q2_new*dp2 + q1_new*dp1)

       q1_new = q1_new * f1
       q2_new = q2_new * f1

       rv_1 = q1_new/(1 - q1_new)
       rv_2 = q2_new/(1 - q2_new)
       
       rc_1 = 0._dp
       rc_2 = 0._dp
       
    else
       ! Condensate is positive
       rc_1 = rc_2 * eta
    endif
    
    cp_1 = (rc_1 + rv_1)/(1 + rc_1 + rv_1)*cp_l + 1/(1 + rc_1 + rv_1)*cp_d
    cp_2 = (rc_2 + rv_2)/(1 + rc_2 + rv_2)*cp_l + 1/(1 + rc_2 + rv_2)*cp_d

    if (T1_new .gt. 273.16) then
       call find_var_lin(T1_guess, lheat, L1)
    else
       L1 = lheat(1)
    endif

    if (T2_guess .gt. 273.16) then
       call find_var_lin(T2_guess, lheat, L2)
    else
       L2 = lheat(1)
    endif

    k_new = cp_1*T1_guess + cp_2*T2_guess + L1*rv_1/(1 + rc_1 + rv_1) + L2*rv_2/(1 + rc_2 + rv_2)
    
  end subroutine k_adj

  subroutine dk_adj(T2_guess, p1, p2, dp1, dp2, dk_adj)
    real(dp), intent(in) :: T2_guess, p1, p2, dp1, dp2
    real(dp), intent(out) :: dk_adj

    real(dp) :: rc_1, rc_2, rv_1, rv_2, k_out_up, k_out_dn, dqsat_dt

    real(dp) :: eps = 1.e-8
    ! Estimate derivative of dk_adj wrt T2_guess

    call k_adj(T2_guess + eps/2, p1,p2,dp1,dp2,q1,q2,rv_1,rv_2,rc_1,rc_2,k_out_up)
    call k_adj(T2_guess - eps/2, p1,p2,dp1,dp2,q1,q2,rv_1,rv_2,rc_1,rc_2,k_out_dn)

    dk_adj = (k_out_up - k_out_dn)/eps
  end subroutine dk_adj
    
  subroutine large_scale_cond(p, T, q, rv, rc)
    real(dp), intent(in) :: p(:)
    real(dp), intent(inout) :: T(:), q(:)
    real(dp), intent(inout) :: rv, rc

    ! Local variables
    integer:: k, n, n_newt
    integer :: npz = size(p)

    real(dp) :: rsat, psat, cp_local, dT , q0, qsat_guess, T_guess
    do k=1,npz
       call r_sat_single(p(k), T(k), rsat)
       
       if (rv(k) .gt. rsat) then
          ! Find adjustment that will bring atmosphere to saturation
          cp_local = cp_l*(rv(k)+rc(k))/(1 + rv(k) + rc(k)) + cp_d/(1 + rv(k) + rc(k))

          dT = 0.0_dp
          n_newt = 100
          T_guess = T(k)
          q0 = rv(k) / (1 + rv(k) + rc(k))
          do while ( (abs(dT) .gt. precision) .and. (n .lt. n_iter))
             call r_sat_single(p(k), T_guess, rsat, psat)
             
             qsat_guess = rsat/(1 + rv(k) + rc(k)) ! Note rv+rc = const throughout to conserve mass

             dqsat_dt = rsat*L*mu_v/Rstar/T_guess**2 * p(k)/(p(k)-psat)
             dT = (L/cp_local * (q0 - qsat_guess) + (T_guess - T(k))) &
                  / (1 + L/cp_local *dqsat_dt)

             T_guess = T_guess + dT
          enddo

          ! Update T(k) with the value from large scale condensation
          T(k) = T_guess

          ! Step ensures (rc + rv)_after = (rc + rv)_before
          rc(k) = rc(k) + rv(k) - rsat
          rv(k) = rsat

          q(k) = rv(k)/(1 + rc(k) + rv(k))
    enddo
  end subroutine large_scale_cond

  subroutine rain_out(delp, T, q, rv, rc, h_cond)
    !==========================================================================
    ! Description
    !==========================================================================
    ! Does rain out and conserves mass + enthalpy in doing so

    !==========================================================================
    ! Input variables
    !==========================================================================    
    real(dp), intent(in),dimension(:) :: q, rv, rc

    !==========================================================================
    ! Output variables
    !==========================================================================    
    real(dp), intent(out) :: h_cond

    !==========================================================================
    ! Input/Output variables
    !==========================================================================    
    real(dp), intent(inout) :: delp(:), T(:)
    
    !==========================================================================
    ! Local variables
    !==========================================================================    
    integer :: k
    integer :: npz = size(p)

    real(dp) :: qc, m_cond_tot, m_cond, rsat, T_guess, q0
    real(dp) :: required
    
    h_cond = 0.0_dp
    m_cond_tot = 0.0_dp
    do k=1,npz
       ! q_cond before rain out
       qc = rc(k)/(1 + rc(k) + rv(k))

       ! Mass loss
       m_cond     = qc * delp(k)/grav
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

       h_cond = h_cond + m_cond * c_l * T(k)
       ! New q once condensate has rained out
       q(k) = rv(k)/(1 + rv(k))

       delp(k) = delp(k) - m_cond*grav

       
    end do

  end subroutine rain_out
  
end module ding_convection
