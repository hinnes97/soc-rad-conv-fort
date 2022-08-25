! Hamish Innes 12/2021, University of Oxford

module adjust_mod
  
  use phys, only : T_TP => H2O_TriplePointT, P_TP => H2O_TriplePointP, L_sub => H2O_L_sublimation, &
       L_vap => H2O_L_vaporization_TriplePoint, CP_v => H2O_cp, CP_d => H2He_solar_cp, &
       mu_d => H2He_solar_MolecularWeight, mu_v => H2O_MolecularWeight, Rstar

  use params, only : dp, Finc, inhibited
  use condense, only: q_sat, cold_trap, dew_point_T
  
  implicit none

  ! FOR FSOLVE PARAMS:
  integer, parameter :: rk = kind(1.0D+00)
  real(kind=rk) :: T1, T2, p1, p2, dp1, dp2, q1, q2, pfactor
  
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
    real(kind=rk), parameter :: tol = 0.00001 ! Tolerance for non-linear solver
    integer, parameter       :: N_iter = 1  ! Number of up-down iterations
    
    real(dp) :: qsat1, qsat2, pfact, grad2, qmin, grad_temp
    real(kind=rk) :: output(1), f_output(1)

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


     subroutine new_adjust(p, delp, T, q, ktrop, grad, olr, mask)
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
    real(kind=rk), parameter :: tol = 0.0001 ! Tolerance for non-linear solver
    integer, parameter       :: N_iter = 1000  ! Number of up-down iterations
    
    real(dp) :: qsat1, qsat2, pfact, grad2, qmin, qcrit
    real(dp) :: qsats(size(p))
    real(kind=rk) :: output(1), f_output(1)

    integer :: n,k
    integer :: info
    integer :: npz
    logical :: conv_switch = .true.
    logical :: condition
    real(dp) :: f ! Helps global energy balance to be reached
    !==========================================================================
    ! Main body
    !==========================================================================
    
    npz = size(p)
    !ktrop = npz
    info = 0
    !    write(*,*) lbound(mask), ubound(mask)
    mask = .false.
    
    
    do n=1,N_iter
       
       call q_sat(p, T, qsats)
       
       do k=npz-1,max(ktrop, 1),-1
          call gradient(p(k+1),T(k+1),grad(k))
          qcrit = (Rstar/mu_v)* T(k+1)/L_vap/(1._dp - mu_d/mu_v)
          if (n.eq. 1) then
             f = (Finc/olr)**(0.01_dp)
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
                 
                 !write(*,*) 'T(k+1), T(k) before', k,T(k+1), T(k)
                 T(k+1) = (T(k)*delp(k) + T(k+1)*delp(k+1))/(delp(k+1) + pfact*delp(k))
                 T(k+1) = f*T(k+1)
                 T(k) = f*T(k+1)*pfact
                 !write(*,*) 'T(k+1), T(k) after', k,T(k+1), T(k)
                 
                 mask(k) = .true.
                 mask(k+1) = .true.

              endif
             
 !          endif
          
          
        enddo

!!$        do k=max(ktrop, 1),npz-1
!!$           call gradient(p(k+1),T(k+1),grad(k))
!!$          qcrit = (Rstar/mu_v)* T(k+1)/L_vap/(1._dp - mu_d/mu_v)
!!$          f = (Finc/olr)**(0.01_dp)
!!$          f = 1.
!!$          if (qsats(k+1) .gt. 1) then
!!$             call dew_point_T(p(k+1), T(k+1))
!!$
!!$             if (f .lt. 1) then
!!$                T(k+1) = T(k+1)*f
!!$             else
!!$                mask(k+1) = .true.
!!$                cycle
!!$             endif
!!$             
!!$             
!!$          endif
!!$!       enddo
!!$!    enddo
!!$
!!$              pfact =exp(grad(k)*log(p(k)/p(k+1)))
!!$              
!!$!              !write(*,*) (T(k) .lt. T(k+1)*pfact*(1. + delta)), q(k), qsat1, q(k+1), qsat2
!!$              !if (((T(k) .lt. T(k+1)*pfact*(1. + delta)) .and. qsats(k+1) .lt. qcrit) .or. &
!!$              !     ((T(k) .gt. T(k+1)*pfact*(1-delta)) .and. (qsats(k+1) .gt. qcrit))) then
!!$
!!$              if ((T(k) .lt. T(k+1)*pfact*(1. + delta))) then
!!$                 
!!$                 !write(*,*) 'T(k+1), T(k) before', k,T(k+1), T(k)
!!$                 T(k+1) = (T(k)*delp(k) + T(k+1)*delp(k+1))/(delp(k+1) + pfact*delp(k))
!!$                 T(k+1) = f*T(k+1)
!!$                 T(k) = f*T(k+1)*pfact
!!$                 !write(*,*) 'T(k+1), T(k) after', k,T(k+1), T(k)
!!$                 
!!$                 mask(k) = .true.
!!$                 mask(k+1) = .true.
!!$
!!$              endif
!!$             
!!$ !          endif
!!$          
!!$          
!!$        enddo
!!$           
!!$
     enddo
     
!           write(*,*) f
           !do k=1,npz-1
              !write(*,*) log(T(k+1)/T(k))/log(p(k+1)/p(k)), grad(k)


  end subroutine new_adjust

  subroutine adjust(p, delp, T, q, mask, olr, ktrop)
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
    logical, intent(out),dimension(:) :: mask ! True where adiabatic
    integer, intent(inout) :: ktrop
    
    !==========================================================================
    ! Mixed input/output
    !==========================================================================
    real(dp), intent(inout), dimension(:) :: T,q ! Temperature and specific humid.

    !==========================================================================
    ! Local variables
    !==========================================================================
    ! Tune these parameters as necessasry
    real(dp), parameter          :: delta = 0.000001 ! Small number speeds up convergence
    real(kind=rk), parameter :: tol = 0.00001 ! Tolerance for non-linear solver
    integer, parameter       :: N_iter = 1  ! Number of up-down iterations
    
    real(dp) :: qsat1, qsat2, pfact, grad, grad2, qmin
    real(kind=rk) :: output(1), f_output(1)

    integer :: n,k
    integer :: info
    integer :: npz
    logical :: conv_switch = .true.

    real(dp) :: f ! Helps global energy balance to be reached
    !==========================================================================
    ! Main body
    !==========================================================================
    
    npz = size(p)
    !ktrop = npz
    info = 0
!    write(*,*) lbound(mask), ubound(mask)
    do n=1,N_iter
!!$
!!$       !Downwards pass
!!$        do k=1,npz-1
!!$           call sat(p(k), T(k), qsat1)
!!$           call sat(p(k+1), T(k+1), qsat2)
!!$
!!$           if ( (q(k) .gt. qsat1*(1.-delta)) .and. (q(k+1) .gt. qsat2*(1.-delta)) ) then
!!$
!!$              ! Equivalent to doing large-scale condensation without latent heating, remove
!!$!              ! if performed elsewhere
!!$
!!$              q(k)   = qsat1
!!$              q(k+1) = qsat2
!!$
!!$              call gradient(p(k+1), T(k+1), grad)
!!$
!!$              pfact = exp(grad*log(p(k)/p(k+1)))
!!$
!!$              ! Test for instability
!!$              if (T(k) .lt. T(k+1)*pfact*(1. + delta) ) then
!!$                 ! INSERT ROOT FINDER HERE FOR T(k+1)
!!$                 ! HACK: set module variables so that fsolve function can have parameters
!!$                 T1 = T(k)!real(T(k), rk)
!!$                 T2 = T(k+1)!real(T(k+1), rk)
!!$                 p1 = p(k)!real(p(k), rk)
!!$                 p2 = p(k+1)!real(p(k+1), rk)
!!$                 dp1 = delp(k)!real(dp(k), rk)
!!$                 dp2 = delp(k+1)!real(dp(k+1), rk)
!!$                 q1 = q(k)!real(q(k), rk)
!!$                 q2 = q(k+1)!real(q(k+1), rk)
!!$                 pfactor = pfact!real(pfact, rk)
!!$
!!$                 !Initial guess
!!$                 output(1) = 250.0D+00
!!$                 call find_my_root(1, output, f_output)
!!$
!!$                 ! Use fsolve to find root of non-linear equation
!!$                 call fsolve(find_my_root, 1, output, f_output, tol, info)
!!$
!!$                 if (info .ne. 1) then
!!$                    write(*,*) 'ERROR IN FSOLVE, CODE: ', info, 'level: ', k
!!$                    cycle
!!$                 endif
!!$                
!!$
!!$                 ! f factor from Malik 2019, tends to 1 when equilibrium reached
!!$                
!!$                 f = (Finc/olr)**(0.01_dp)
!!$
!!$                 if ( (output(1)*f - T(k+1) ).gt. 5.) then
!!$                    T(k+1) = T(k+1) + 5.
!!$                 else if ((output(1)*f - T(k+1) ) .lt. -5.) then
!!$                    T(k+1) = T(k+1) - 5.
!!$                 else
!!$                    T(k+1) = output(1)*f
!!$                 endif
!!$
!!$                 T(k)   = T(k+1)*pfact*f
!!$                 call sat(p(k), T(k), q(k))
!!$                 call sat(p(k+1),T(k+1), q(k+1))
!!$
!!$                 mask(k) = .true.
!!$                 mask(k+1) = .true.
!!$              else
!!$                 mask(k) = .false.
!!$                 mask(k+1) = .false.
!!$              endif
!!$
!!$          endif
!!$       enddo
       
       ! Upwards pass
       qmin = 10000.
        do k=npz-1,1,-1
           
           !if ( (q(k) .gt. qsat1*(1.-delta)) .and. (q(k+1) .gt. qsat2*(1.-delta)) ) then
           
              call sat(p(k), T(k), qsat1)
              call sat(p(k+1), T(k+1), qsat2)
              write(*,*) 'qsat1, qsat2', qsat1, qsat2, qmin
              if (qsat2 < qmin) then
                 qmin = qsat2
              else
                 q(1:k+1) = qmin
                 write(*,*) 'cold trap k', k+1
                 exit
              endif

              if (qsat1 < qmin) then
                 qmin = qsat1
              else
                 write(*,*) 'cold trap k', k
                 q(1:k) = qmin
                 exit
              endif
              
                 
             ! Equivalent to doing large-scale condensation without latent heating, remove
             ! if performed elsewhere

              if (qsat2 .gt. 1.) then
                 q(k+1) = 1.
                 call dew_point_T(p(k+1), T(k+1))
                 cycle
              endif
              
              q(k)   = qsat1
              q(k+1) = qsat2


              
             call gradient(p(k+1), T(k+1), grad)

             pfact =exp(grad*log(p(k)/p(k+1)))

             !write(*,*) (T(k) .lt. T(k+1)*pfact*(1. + delta)), q(k), qsat1, q(k+1), qsat2
             if ((T(k) .lt. T(k+1)*pfact*(1. + delta))) then
                ! INSERT ROOT FINDER HERE FOR T(k+1)
                ! HACK: set module variables so that fsolve function can have parameters
                !T1 = T(k)!real(T(k), rk)
                !T2 = T(k+1)!real(T(k+1), rk)
                !p1 = p(k)!, rk)
                !p2 = p(k+1)!, rk)
                !dp1 = delp(k)!, rk)
                !dp2 = delp(k+1)!, rk)
                !q1 = q(k)!, rk)
                !q2 = q(k+1)! rk)
                !pfactor = pfact!, rk)

                !Initial guess
                !output(1) = 250.0D+00
                !call find_my_root(1, output, f_output)

                !call fsolve(find_my_root, 1, output, f_output, tol, info)

                !if (info .ne. 1) then
                !   write(*,*) 'ERROR IN FSOLVE, CODE: ', info
                !   cycle
                !endif

                

                ! f factor from Malik 2019, tends to 1 when equilibrium reached
                f = (Finc/olr)**(0.01_dp)

                !if ( (output(1)*f - T(k+1) ).gt. 5.) then
                !   T(k+1) = T(k+1) + 5.
                !else if ((output(1)*f - T(k+1) ) .lt. -5.) then
                !   T(k+1) = T(k+1) - 5.
                !else
                !   T(k+1) = output(1)*f
                !endif

                T(k) = T(k+1)*pfact*f
!                write(*,*) 'T2, T1', T(k+1), T(k)
                call sat(p(k), T(k), q(k))
                !call sat(p(k+1),T(k+1), q(k+1))

                mask(k) = .true.
                mask(k+1) = .true.
!                conv_switch = .false.
             else
                mask(k) = .false.
                mask(k+1) = .false.
!                ktrop = min(k+1, ktrop)
!                if ( .not. conv_switch)  exit
                
             endif
!          endif
          
       enddo
       
          enddo
                   
          


  end subroutine adjust

  subroutine gradient(p,T, dlnTdlnp)
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

    !==========================================================================
    ! Local variables
    !==========================================================================
    real(dp) :: eps = mu_v/mu_d
    real(dp) :: L, psat, qsat, rsat, num, denom, temp

    !==========================================================================
    ! Main body
    !==========================================================================
    
    L = L_vap

    call sat(p, T, qsat, rsat, psat)

    num   = 1 + (L*mu_d/Rstar/T)*rsat
    denom = 1 + ((cp_v/cp_d) + ((L*mu_v/Rstar/T) - 1)*(L/cp_d/T) )*rsat

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
    real(dp)    :: L, eps, psat_temp, rsat_temp

    !==========================================================================
    ! Main body
    !==========================================================================
    
    eps = mu_v/mu_d

    L = L_vap


    psat_temp = P_TP * exp(-L/Rstar*mu_v * (1./T - 1./T_TP))
    
    rsat_temp = psat_temp / (p-psat_temp) * eps
    qsat = rsat_temp / (1 + rsat_temp)

    if(present(psat)) psat = psat_temp
    if(present(rsat)) rsat = rsat_temp

  end subroutine sat


  subroutine find_my_root(n, T, fvec)

    !==========================================================================
    ! Description
    !==========================================================================
    ! Function passed to fsolve, which will find its roots in T. Represents the
    ! conservation of moist enthalpy in the dilute limit

    !==========================================================================
    ! Input variables
    !==========================================================================
    integer, intent(in) :: n  ! Number of non-linear equations to be solved (1)
    real(kind=rk),    intent(in) :: T(n) ! Temperature to be solved for

    !==========================================================================
    ! Output variables
    !==========================================================================
    real(kind=rk),    intent(out) :: fvec(n) ! Value of function at T

    !==========================================================================
    ! Local variables
    !==========================================================================
    real(kind=rk) :: q1_new, q2_new, L
    real(dp) :: q1_temp, q2_temp

    !==========================================================================
    ! Main body
    !==========================================================================
    
    call sat(p1, T(1)*pfactor, q1_temp)
    call sat(p2, T(1), q2_temp)

    q1_new = q1_temp
    q2_new = q2_temp
    
    if (T2 .gt. real(T_TP, rk)) then
       L = real(L_vap, rk)
    else
       L = real(L_sub, rk)
    endif

    ! Moist enthalpy equation divided by cp*T, to be roughly order unity for solver
    fvec(1) = 1 - (T1*dp1 + T2*dp2 + L/cp_d*(q1 - q1_new)*dp1 + L/cp_d*(q2 - q2_new)*dp2)&
         / (dp2 + dp1*pfactor) / T(1)
    
  end subroutine find_my_root
end module adjust_mod
