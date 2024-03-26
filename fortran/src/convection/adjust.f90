! Hamish Innes 12/2021, University of Oxford

module adjust_mod
  
  use phys, only : H2O, H2He_solar, Rstar
  use params, only : dp, Finc, inhibited, accelerate
  use condense, only: q_sat, cold_trap, dew_point_T
  use tables, only : phase_grad, lheat, satur, find_var_lin, find_var_loglin
  use atmosphere, only : mmw_dry, cp_dry, nqr, nqt
  implicit none

contains


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
    real(dp) :: qsats(size(p),nqt), qcrits(size(p))

    real(dp) :: grad_check(size(p)), grad_true(size(p))

    integer :: n,k
    integer :: info, counter
    integer :: npz
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
       
    do n=1,N_iter

       grad_check = 0.0_dp
       grad_true = 0.0_dp
       
       call q_sat(p, T, qsats)
       
       do k=npz-1,max(ktrop, 1),-1
          call gradient(p(k+1),T(k+1),cp_dry(k), mmw_dry(k), grad(k), temp)
          qcrit = 1./(1._dp - H2He_solar%mmw/H2O%mmw) /temp
          
          qcrits(k) = qcrit
          if (n.eq. 1 .and. tstep .gt. 1000 .and. accelerate) then
             f =  (Finc/olr)**(0.001_dp)
          else
             f = 1.
          endif
             
          if (qsats(k+1,1) .gt. 1) then
             call dew_point_T(p(k+1), T(k+1))

             if ((Finc/olr)**(0.01_dp) .lt. 1) then
                T(k+1) = T(k+1)*(Finc/olr)**(0.01_dp)
             else
                mask(k+1) = 1
                cycle
             endif
             
             
          endif

!       enddo
!    enddo
    
              pfact =exp(grad(k)*log(p(k)/p(k+1)))
              
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
        enddo ! do n=1,N_iter
        
  end subroutine new_adjust

  subroutine gradient(p,T, cpd, mmw_d,dlnTdlnp, dlnpsat_dlnt)
    !==========================================================================
    ! Description
    !==========================================================================
    ! Calculates the moist adiabatic gradient d(log(T))/d(log(p)) as in
    ! Ding and Pierrehumbert 2016
    
    !==========================================================================
    ! Input variables
    !==========================================================================
    real(dp), intent(in) :: p, T !Pressure and temperature
    real(dp), intent(in) :: cpd, mmw_d
    !==========================================================================
    ! Output variables
    !==========================================================================
    real(dp), intent(out) :: dlnTdlnp ! Moist adiabatic gradient

    real(dp),intent(inout),optional :: dlnpsat_dlnt
    
    !==========================================================================
    ! Local variables
    !==========================================================================

    real(dp) :: L, psat, qsat, rsat, num, denom, temp, t2
    
    !==========================================================================
    ! Main body
    !==========================================================================
    
    if (T .lt. 273.16) then
      L = lheat(1)
   else
      call find_var_lin(T, lheat, L)
    endif
    !L = 2.5e6
    call sat(p, T, cpd, mmw_d,qsat, rsat, psat)
    
    num   = 1 + (L*mmw_d/Rstar/T)*rsat
    !denom = 1 + ((cp_v/cp_d) + ((L*mu_v/Rstar/T) - 1)*(L/cp_d/T) )*rsat
    
    ! Changed for more accurate version
    if (T .gt. 273.16) then
       call find_var_lin(T, phase_grad,t2)
    else
       L = lheat(1)
       t2 = L*H2O%mmw/Rstar/T
    endif
    

    if (present(dlnpsat_dlnt)) then
       dlnpsat_dlnt = t2
    endif
    
    denom = 1 + ((H2O%cp/cpd) + t2/cpd*(L/T - Rstar/H2O%mmw)  )*rsat

    temp = Rstar/mmw_d/cpd * num/denom
    
    dlnTdlnp = 1. / ((psat/p)*L*H2O%mmw/Rstar/T + (p - psat)/p/temp)
 
  end subroutine gradient
  
  subroutine sat(p,T, cpd,mmw_d,qsat, rsat, psat)
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
    real(dp), intent(in) :: mmw_d, cpd
    !==========================================================================
    ! Output variables
    !==========================================================================
    
    real(dp), intent(out) :: qsat ! Mass concentration
    real(dp), intent(out), optional :: rsat, psat ! Mass mixing ratio, sat pressure

    !==========================================================================
    ! Local variables
    !==========================================================================
    real(dp)    :: L, rsat_temp, eps, psat_temp
    !==========================================================================
    ! Main body
    !==========================================================================
    
    eps = H2O%mmw/mmw_d


    if (T .gt. 273.16) then
       call find_var_loglin(T, satur, psat_temp)
    else
       L = lheat(1)
       psat_temp = H2O%TP_P * exp(-L/Rstar*H2O%mmw * (1./T - 1./H2O%tp_t))
    endif
    
    rsat_temp = psat_temp / (p-psat_temp) * eps
    qsat = rsat_temp / (1 + rsat_temp)

    if(present(psat)) psat = psat_temp
    if(present(rsat)) rsat = rsat_temp

  end subroutine sat


end module adjust_mod
