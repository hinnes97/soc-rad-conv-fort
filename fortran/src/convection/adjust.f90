! Hamish Innes 12/2021, University of Oxford

module adjust_mod
  
  use phys, only : T_TP => H2O_TriplePointT, P_TP => H2O_TriplePointP, L_sub => H2O_L_sublimation, &
       L_vap => H2O_L_vaporization_TriplePoint, CP_v => H2O_cp, CP_d => H2He_solar_cp, &
       mu_d => H2He_solar_MolecularWeight, mu_v => H2O_MolecularWeight, Rstar

  use params, only : dp, Finc, inhibited, accelerate
  use condense, only: q_sat, cold_trap, dew_point_T
  use tables, only : phase_grad, lheat, satur, find_var_lin, find_var_loglin
  
  implicit none

    
contains


  subroutine calc_q_and_grad(p, delp, T, q, mask, olr,ktrop)
    !==========================================================================
    ! Description
    !==========================================================================
    ! Performs moist adiabatic adjustment, conserving column-integrated moist
    ! enthalpy: \int(cp*T + Lq)*dp in the dilute limit, as in Manabe 1965. The
    ! adjustment is performed pairwise in layers, until convergence is reached.
    ! For two layers (labelled 1 and 2, with post-adjustment state having
    ! dashes), the conservation of moist enthalpy is:
    !
    ! cp*(T1*dp1+T2*dp2) + L*(q1(T1)*dp1+q2(T2)*dp2)
    !                    = cp*(T1'*dp1+T2'*dp2) + L*(q1'(T1')*dp1+q2'(T2')*dp2)
    !
    ! We can relate T1' to T2' by T1' = T2'*(p1/p2)**(dlnT/dlnp), and then since
    ! qi' is related to Ti' by the Clausius Clapyeron relation, we can solve for
    ! T2' using non-linear equation solver fsolve (More 1980).
    !
    ! Algorithm heavily based on Ray Pierrehumbert's dry convection code

    !==========================================================================
    ! Input variables
    !==========================================================================    
    real(dp), intent(in)   , dimension(:) :: p,delp ! Pressure, pressure thickness
    real(dp), intent(in) :: olr
    !==========================================================================
    ! Output variables
    !==========================================================================
    logical, intent(in),dimension(:) :: mask ! True where adiabatic
    integer, intent(out) :: ktrop
    !==========================================================================
    ! Mixed input/output
    !==========================================================================
    real(dp), intent(inout), dimension(:) :: T,q ! Temperature and specific humid.
    

    !==========================================================================
    ! Local variables
    !==========================================================================
    ! Tune these parameters as necessasry
    real(dp), parameter          :: delta = 0.000001 ! Small number speeds up convergence
    integer, parameter       :: N_iter = 1  ! Number of up-down iterations
    
    real(dp) :: qsat1, qsat2, pfact, grad2, qmin, grad_temp

    integer :: n,k
    integer :: info
    integer :: npz
    logical :: conv_switch = .true.

    real(dp) :: f ! Helps global energy balance to be reached
    real(dp) :: qcrit
    !==========================================================================
    ! Main body
    !==========================================================================
    
    npz = size(p)
    !ktrop = npz
    info = 0
    qmin = 10000.

    call q_sat(p, T, q)
    call cold_trap(q, ktrop, p, T)
    
   end subroutine calc_q_and_grad


   subroutine new_adjust(p, delp, T, q, ktrop, grad, olr, mask, tstep)
    !==========================================================================
    ! Description
    !==========================================================================
    ! Performs moist adiabatic adjustment, conserving column-integrated moist
    ! enthalpy: \int(cp*T + Lq)*dp in the dilute limit, as in Manabe 1965. The
    ! adjustment is performed pairwise in layers, until convergence is reached.
    ! For two layers (labelled 1 and 2, with post-adjustment state having
    ! dashes), the conservation of moist enthalpy is:
    !
    ! cp*(T1*dp1+T2*dp2) + L*(q1(T1)*dp1+q2(T2)*dp2)
    !                    = cp*(T1'*dp1+T2'*dp2) + L*(q1'(T1')*dp1+q2'(T2')*dp2)
    !
    ! We can relate T1' to T2' by T1' = T2'*(p1/p2)**(dlnT/dlnp), and then since
    ! qi' is related to Ti' by the Clausius Clapyeron relation, we can solve for
    ! T2' using non-linear equation solver fsolve (More 1980).
    !
    ! Algorithm heavily based on Ray Pierrehumbert's dry convection code

    !==========================================================================
    ! Input variables
    !==========================================================================    
       real(dp), intent(in)   , dimension(:) :: p,delp ! Pressure, pressure thickness
       real(dp), intent(inout), dimension(:) :: grad
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

    !==========================================================================
    ! Local variables
    !==========================================================================
    ! Tune these parameters as necessasry
    real(dp), parameter          :: delta = 0.0001 ! Small number speeds up convergence

    integer, parameter       :: N_iter = 1000 ! Number of up-down iterations
    
    real(dp) :: qsat1, qsat2, pfact, grad2, qmin, qcrit,temp
    real(dp) :: qsats(size(p)), qcrits(size(p))

    real(dp) :: grad_check(size(p)), grad_true(size(p))

    integer :: n,k
    integer :: info, counter
    integer :: npz
    logical :: conv_switch = .true.
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
    mask = .false.
    quit_adjust = .false.

    ! Uncomment this line and line at bottom for exit at tolerance
    !do while (.not. quit_adjust .and. n .lt. N_iter)
       
    do n=1,N_iter

       grad_check = 0.0_dp
       grad_true = 0.0_dp
       
       call q_sat(p, T, qsats)
       
       do k=npz-1,max(ktrop, 1),-1
          call gradient(p(k+1),T(k+1),grad(k), temp)
          qcrit = 1./(1._dp - mu_d/mu_v) /temp
          
          qcrits(k) = qcrit
          if (n.eq. 1 .and. tstep .gt. 1000 .and. accelerate) then
             f =  (Finc/olr)**(0.001_dp)
          else
             f = 1.
          endif
             
          if (qsats(k+1) .gt. 1) then
             call dew_point_T(p(k+1), T(k+1))

             if ((Finc/olr)**(0.01_dp) .lt. 1) then
                T(k+1) = T(k+1)*(Finc/olr)**(0.01_dp)
             else
                mask(k+1) = .true.
                cycle
             endif
             
             
          endif

!       enddo
!    enddo
    
              pfact =exp(grad(k)*log(p(k)/p(k+1)))
              
              if (inhibited) then
                 ! Do moisture inhibition
                 condition = (  (T(k) - T(k+1)*pfact*(1 + delta))*(qcrit - q(k+1)) .lt. 0)
              else
                 ! Regular criterion
                 condition = ( (T(k) - T(k+1)*pfact*(1+delta)) .lt. 0 )
              endif
              
              if (condition) then
                 T(k+1) = (T(k)*delp(k) + T(k+1)*delp(k+1))/(delp(k+1) + pfact*delp(k))
                 T(k+1) = f*T(k+1)
                 T(k) = f*T(k+1)*pfact
                 mask(k) = .true.
                 mask(k+1) = .true.
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
        enddo ! do n=1,N_iter
        
  end subroutine new_adjust

  subroutine gradient(p,T, dlnTdlnp, dlnpsat_dlnt)
    !==========================================================================
    ! Description
    !==========================================================================
    ! Calculates the moist adiabatic gradient d(log(T))/d(log(p)) as in
    ! Ding and Pierrehumbert 2016
    
    !==========================================================================
    ! Input variables
    !==========================================================================
    real(dp), intent(in) :: p, T !Pressure and temperature

    !==========================================================================
    ! Output variables
    !==========================================================================
    real(dp), intent(out) :: dlnTdlnp ! Moist adiabatic gradient

    real(dp),intent(inout),optional :: dlnpsat_dlnt
    !==========================================================================
    ! Local variables
    !==========================================================================
    real(dp) :: eps = mu_v/mu_d
    real(dp) :: L, psat, qsat, rsat, num, denom, temp, start, end,t2, ttt

    !==========================================================================
    ! Main body
    !==========================================================================
    
    if (T .lt. 273.16) then
      L = lheat(1)
   else
      call find_var_lin(T, lheat, L)
    endif
    !L = 2.5e6
    call sat(p, T, qsat, rsat, psat)
    
    num   = 1 + (L*mu_d/Rstar/T)*rsat
    !denom = 1 + ((cp_v/cp_d) + ((L*mu_v/Rstar/T) - 1)*(L/cp_d/T) )*rsat
    
    ! Changed for more accurate version
    if (T .gt. 273.16) then
       call find_var_lin(T, phase_grad,t2)
    else
       L = lheat(1)
       t2 = L*mu_v/Rstar/T
    endif
    

    if (present(dlnpsat_dlnt)) then
       dlnpsat_dlnt = t2
    endif
    
    denom = 1 + ((cp_v/cp_d) + t2/cp_d*(L/T - Rstar/mu_v)  )*rsat

    temp = Rstar/mu_d/cp_d * num/denom
    
    dlnTdlnp = 1. / ((psat/p)*L*mu_v/Rstar/T + (p - psat)/p/temp)
    
    
  end subroutine gradient
  
  subroutine sat(p,T, qsat, rsat, psat)
    !==========================================================================
    ! Description
    !==========================================================================

    ! Calculates the saturation vapour pressure, mass mixing ratio and mass
    ! concentration of water vapour according to simple Clausius Clapeyron
    ! relation. Change this function to some parametrisation if more accuracy
    ! required.

    !==========================================================================
    ! Input variables
    !==========================================================================
    
    real(dp), intent(in)  :: p, T ! Pressure, temperature

    !==========================================================================
    ! Output variables
    !==========================================================================
    
    real(dp), intent(out) :: qsat ! Mass concentration
    real(dp), intent(out), optional :: rsat, psat ! Mass mixing ratio, sat pressure

    !==========================================================================
    ! Local variables
    !==========================================================================
    real(dp)    :: L, eps, psat_temp, rsat_temp, start, end, p_sat_old,ttt

    !==========================================================================
    ! Main body
    !==========================================================================
    
    eps = mu_v/mu_d


    if (T .gt. 273.16) then
       call find_var_loglin(T, satur, psat_temp)
    else
       L = lheat(1)
       psat_temp = P_TP * exp(-L/Rstar*mu_v * (1./T - 1./T_TP))
    endif
    
    rsat_temp = psat_temp / (p-psat_temp) * eps
    qsat = rsat_temp / (1 + rsat_temp)

    if(present(psat)) psat = psat_temp
    if(present(rsat)) rsat = rsat_temp

  end subroutine sat


end module adjust_mod
