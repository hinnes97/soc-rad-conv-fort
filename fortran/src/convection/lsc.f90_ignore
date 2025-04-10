module lsc_mod
  use condense, only: q_sat, q_sat_single, dqsat_dt_single, dew_point_T
  use atmosphere, only: get_cp, nqt, q_orig
  use tables, only: lheat, find_var_lin
  use params, only: dp
  

  implicit none
  
contains

  subroutine lsc(p, T, q, q_cond, dT_lsc, tstep, delp)
    real(dp), dimension(:), intent(in) :: p, delp
    real(dp), dimension(:), intent(in) :: T
    integer ,intent(in) :: tstep
    real(dp), dimension(:,:), intent(inout) :: q
    real(dp), dimension(:), intent(inout) :: q_cond
    real(dp), dimension(:), intent(out) :: dT_lsc

    integer, parameter :: N_iter = 100
    integer :: npz, k, n, l_rec, k_rec,m

    real(dp), dimension(size(q,1), size(q,2)) :: qsats
    real(dp), dimension(size(p)) :: T_dash
    real(dp) :: T0, T_new, qsat(nqt), T_guess, dqsatdt, cp, L, h_con, dq
    
    npz = size(p)
    T_dash = T
    call q_sat(p, T, qsats)

    !if (mod(tstep, 1) .eq. 0) then
    !   write(*,*) 'Before'
    !   do k=1,npz
    !      write(*,*) sum(q(k,:)), q(k,1), qsats(k,1)
    !   enddo
    !endif
    do k=10,npz
       if (q(k,1) - qsats(k,1) .gt. 0) then
          ! Condensation here
          call get_cp(q(k,:), cp)
          cp = 5500.
          call find_var_lin(T_dash(k), lheat, L)
          L = max(L, 0._dp)
          L = 2.6e6
          T0 = T(k)
          T_guess = T0
          do n=1,N_iter
             call q_sat_single(p(k), T_guess, k, qsat)
             call dqsat_dt_single(p(k), T_guess, k, dqsatdt, qsat)
             T_new = T_guess + ( cp*(T0 - T_guess) + L*(q(k,1) - qsat(1)))/(cp + L*dqsatdt)
             qsat(1) = qsat(1) + dqsatdt*(T_new - T_guess)
             if (abs(T_new - T_guess)/T_guess .lt. 1.e-8) exit
             T_guess = T_new
          enddo
          T_guess = T_new

          if(abs(T_guess - T(k)) .gt. 10._dp) then
             T_dash(k)  = min(max(T(k)-10.0_dp, T_guess), T(k) + 10.0_dp)
             q_cond(k) = q_cond(k) + cp*(T_dash(k) - T(k))/L
             dq = -cp*(T_dash(k) - T(k))/L
             write(*,*) 'lsc limited'
             if (q_cond(k) .lt. 0) write(*,*) 'q_cond .lt. 0 limited', k, T_guess - T(k)
          else
             T_dash(k) = T_guess
             q_cond(k) = q_cond(k) + q(k,1) - qsat(1)
             dq = qsat(1) - q(k,1)
             if (q_cond(k) .lt. 0) write(*,*) 'q_cond .lt. 0, normal', k, q(k,1), qsat(1), dqsatdt*(T_guess - T(k))
           endif

             
          !if (abs(dqsatdt*(T_guess-T(k))) .gt. 0.1_dp) then
          !   q_cond(k) = 0.1_dp
          !   T_guess = q_cond(k)*L/cp
          !endif
          !if (abs(T_guess - T(k)) .gt. 10.0_dp) write(*,*) 'gt 10 without limit!', k, T_guess - T(k)
          !T_dash(k) = max(min(T0+10.0_dp, T_guess), T0-10.0_dp)
          
          !T_dash(k) = min(max(T(k)-1.0_dp, T_guess), T(k) + 1.0_dp)
          
          !T_dash(k) = T_guess
          call q_sat_single(p(k), T_dash(k), k, qsat)
          ! if (mod(tstep,1000) .eq. 0 )  then
          !    write(*,*) 'lsc test', q(k,1), qsats(k,1), qsat(1), T0, T_dash(k)
          ! endif
          
          !q_cond(k)  = max(min(-dqsatdt*(T_guess - T(k)),0.01),-0.01)
          q(k,1) = q(k,1) + dq
          do m=2,nqt
             q(k,m) = (1. - q(k,1))/(1 - q_orig(k,1)) * q_orig(k,m)
          enddo

       endif
    enddo
    dT_lsc = T_dash - T

    if (mod(tstep, 1000) .eq. 0 )then
       h_con = 0.0
       L = 2.6e6
       do k=1,npz
          h_con = h_con + delp(k)*q_cond(k)*L
       enddo
       write(*,*) 'ENTHALPY OF CONDENSATE VS CHANGE IN TEMP', h_con, sum(delp*5500.*dT_lsc)
    endif

    if (mod(tstep, 1000) .eq. 0) then
       write(*,*) 'LSC AFTER -----------------------------------------'
       do k=1,npz
          if (dT_lsc(k) .gt. 0) then
             call q_sat_single(p(k), T_dash(k), k, qsat)
             write(*,*) k, q(k,1), qsat(1), q(k,1) - qsat(1), T(k), T_dash(k)
          endif
       enddo
    endif
    !stop
  end subroutine lsc

  subroutine condensation(p, T, q, delp, q_cond, dT_cond,tstep)
    real(dp), intent(in), dimension(:) :: p, delp, T
    real(dp), intent(inout), dimension(:) ::  q_cond
    integer, intent(in) :: tstep
    real(dp), intent(inout), dimension(:,:) :: q
    real(dp), intent(out), dimension(:) :: dT_cond

    integer :: k, npz, n, m
    integer, parameter :: n_iter = 100
    real(dp) ::  mcond, T0, T_guess, L, cp, T_new, dqsatdt, dew_t, h_con, h_tot
    real(dp), dimension(size(p)) :: T_dash, dmcond
    real(dp), dimension(size(q,2)) :: qsat
    real(dp), dimension(size(q,1), size(q,2)) :: qsats

    T_dash = T
    npz = size(p)
    mcond = 0.0_dp
    call q_sat(p, T, qsats)
    L = 2.6e6
    if (mod(tstep,1000) .eq. 0) then
       h_con = 0.0
       do k=1,npz
          h_con = h_con + L*q_cond(k)*delp(k)
       enddo
       write(*,*) 'INITIAL COND +TOT ENTHALPY', h_con, sum((L*(q(:,1)) + 5500.*T)*delp)
    endif

    dmcond = 0.0
    if (mod(tstep,1000) .eq. 0) then
    do k=1,npz
       write(*,*) 'q_cond', k, q_cond(k)
    enddo
    endif
       
    do k=10,npz
       if (q_cond(k) .gt. 0.0_dp) then
          mcond = mcond + q_cond(k)*delp(k)
          q_cond(k) = 0.0
       endif
       if (mod(tstep,1000) .eq. 0) write(*,*) 'BEFORE', k, sum(delp*q_cond) + mcond, sum(dmcond)
       if (qsats(k,1) .gt. q(k,1) .and. mcond > 0.0_dp) then
!          write(*,*) 'Re-evaporating!', k, mcond, qsats(k,1), q(k,1)
          call get_cp(q(k,:), cp)
          cp = 5500.
          call find_var_lin(T(k), lheat, L)
          L = max(L, 0._dp)
          L = 2.6e6
          ! Check if we can re-evaporate all of the condensate
           !write(*,*) 'mcondd', mcond, (qsats(k,1) - q(k,1))*delp(k)
          if (mcond .gt. 0.5_dp*(qsats(k,1) - q(k,1))*delp(k)) then
!              write(*,*) 'Trying to re-evaporate part'
             ! Do identical procedure to lsc - iterate to find equilibrium
             T0 = T_dash(k)
             T_guess = T0
             do n=1,N_iter
                call q_sat_single(p(k), T_guess, k, qsat)
                call dqsat_dt_single(p(k), T_guess, k, dqsatdt, qsat)
                T_new = T_guess + ( cp*(T0 - T_guess) + L*(q(k,1) - qsat(1)))/(cp + L*dqsatdt)
                qsat(1) = qsat(1) + dqsatdt*(T_new - T_guess)
                if (abs(T_new- T_guess)/T_guess .lt. 1.e-8) exit
                T_guess = T_new
             enddo
             T_guess = T_new

             if (abs(T_guess - T(k)) .gt. 10.0_dp) then
                T_dash(k) = min(max(T(k)-10._dp, T_guess), T(k)+10.0_dp)
                q(k,1) = q(k,1) - cp/L * (T_dash(k) - T(k))
                mcond = mcond + cp/L * (T_dash(k) - T(k))*delp(k)
                dmcond(k) = dmcond(k) + cp/L*(T_dash(k) - T(k))*delp(k)
                if (mod(tstep,1000) .eq. 0) then
                   write(*,*) '1', k, dmcond(k)*L, (T_dash(k) - T(k))*cp*delp(k)
                endif
             else
                T_dash(k) = T_guess
                mcond = mcond - (qsat(1) - q(k,1))*delp(k)
                dmcond(k) = -(qsat(1) - q(k,1))*delp(k)
                q(k,1) = qsat(1)
                if (mod(tstep,1000) .eq. 0) then
                   write(*,*) '2', k, dmcond(k)*L, (T_dash(k) - T(k))*cp*delp(k)
                endif
             endif

             do m=2,nqt
                q(k,m) = (1. - q(k,1))/(1 - q_orig(k,1)) * q_orig(k,m)
             enddo

             !mcond = mcond - (qsat(1) - q(k,1))*delp(k)
             !q(k,:) = qsat
!              write(*,*) 'mcond some', mcond, q(k,1)
          else
             ! If mcond < available "space", then attempt to evaporate all the stuff in this layer

             ! What happens of we are at psat>p
             if (qsats(k,1) .gt. 1) then
                ! Careful! Try and make this layer pure steam
                call dew_point_T(p(k), dew_t)
                if (T(k) - L/cp*mcond/delp(k) .lt. dew_t) then
                   T_dash(k) = dew_t
                   mcond = mcond - (T(k) - dew_t)*cp*delp(k)/L
                   dmcond(k) = - (T(k) - dew_t)*cp*delp(k)/L

                   if (mod(tstep,1000) .eq. 0) then
                      write(*,*) '3', k, dmcond(k)*L, (T_dash(k) - T(k))*cp*delp(k)
                   endif

                   q(k,1) = 1.0_dp
                   q(k,2:nqt) = 0.0_dp
                   !write(*,*) k, 'COND, MIGHT NOT CONSERVE ENTH',  (T(k) - dew_t)*cp*delp(k), 
                   cycle
                else
                   T_dash(k) = T(k) - L/cp*mcond/delp(k)
                   q(k,1) = q(k,1) + mcond/delp(k)
                   do m=2,nqt
                      q(k,m) = (1. - q(k,1))/(1 - q_orig(k,1)) * q_orig(k,m)
                   enddo
                   dmcond(k) = -mcond
                   mcond = 0.0_dp
                   if (mod(tstep,1000) .eq. 0) then
                      write(*,*) '4', k, dmcond(k)*L, (T_dash(k) - T(k))*cp*delp(k)
                   endif

                endif
                
             else
!             write(*,*) 'Re-evaporating all', k
                T_new = T(k)  - L/cp * mcond/delp(k)
                T_dash(k) = T_new
             q(k,1) = q(k,1) + mcond/delp(k)
             dmcond(k) = -mcond
             mcond = 0.0_dp
             if (mod(tstep,1000) .eq. 0) then
                write(*,*) '5', k, dmcond(k)*L, (T_dash(k) - T(k))*cp*delp(k)
             endif

             do m=2,nqt
                q(k,m) = (1. - q(k,1))/(1 - q_orig(k,1)) * q_orig(k,m)
             enddo
             
             ! Check if evaporation of all mcond has actually produced oversaturation again
             call q_sat_single(p(k), T_new, k, qsat)
!             write(*,*) 'qsat', T(k), T_new, qsat(1), L/cp* mcond/delp(k), q(k,1)
             if (q(k,1) > qsat(1)) then
                ! Recondense and move on to next layer
                T0 = T_new
                T_guess = T0
                if (mod(tstep,1000) .eq. 0) write(*,*) 'RECONDENSING', q(k,1), qsat(1)
                do n=1,N_iter
                   call q_sat_single(p(k), T_guess, k, qsat)
                   call dqsat_dt_single(p(k), T_guess, k, dqsatdt, qsat)
                   T_new = T_guess + ( cp*(T0 - T_guess) + L*(q(k,1) - qsat(1)))/(cp + L*dqsatdt)
                   qsat = qsat  + (T_new - T_guess)*dqsatdt
                   if (abs(T_new- T_guess)/T_guess .lt. 1.e-8) exit
!                   write(*,*) T_guess, t_new, qsat(1)
                   T_guess = T_new
                enddo
                
                T_guess = T_new
                if (mod(tstep,1000) .eq. 0) write(*,*) '7.5', q(k,1), qsat(1)
              if (abs(T_guess - T(k)) .gt. 10.0_dp) then
                T_dash(k) = min(max(T(k)-10._dp, T_guess), T(k)+10.0_dp)
                q(k,1) = q(k,1) - cp/L * (T_dash(k) - T(k))
                mcond = mcond + cp/L * (T_dash(k) - T(k))*delp(k)
                dmcond(k) = cp/L * (T_dash(k) - T(k))*delp(k)
                if (mod(tstep,1000) .eq. 0) then
                   write(*,*) '6', k, dmcond(k)*L, (T_dash(k) - T(k))*cp*delp(k)
                endif

             else
                T_dash(k) = T_guess
                if (qsat(1) .gt. 1.0_dp) then
                   write(*,*) 'COND, WEIRD q>1'
                   qsat(1) = 1.0_dp
                endif
                if (mod(tstep,1000) .eq. 0) write(*,*) -(qsat(1) - q(k,1)), mcond, mcond-(qsat(1) - q(k,1))*delp(k)
                mcond = mcond - (qsat(1) - q(k,1))*delp(k)
                dmcond(k)  = dmcond(k) - (qsat(1) - q(k,1))*delp(k)
                q(k,1) = qsat(1)
                if (mod(tstep,1000) .eq. 0) then
                   write(*,*) '7', k, dmcond(k)*L, (T_dash(k) - T(k))*cp*delp(k)
                endif

             endif

             do m=2,nqt
                q(k,m) = (1. - q(k,1))/(1 - q_orig(k,1)) * q_orig(k,m)
             enddo

                   
                ! T_guess = T_new
                ! T_dash(k) = T_guess
                ! call q_sat_single(p(k), T_guess, k, qsat)
                ! mcond = (q(k,1) - min(qsat(1), 1.0_dp))*delp(k)
                ! q(k,:) = qsat(:)
!                write(*,*) 'mcond oversat', mcond, k, T_guess
             
          endif
          endif
          endif
       endif
       if (mod(tstep,1000) .eq. 0) write(*,*) 'AFTER', k, sum(q_cond*delp) + mcond, sum(dmcond)
    enddo
    
    dT_cond = T_dash  - T
    if (mod(tstep,1000) .eq. 0) then
       write(*,*) 'ENTHALPY CONDENSATE VS T CHANGE', h_con, sum(5500*dT_cond*delp), sum(dmcond), sum(q_cond*delp), &
            sum((L*(q(:,1)) + 5500.*T_dash)*delp)
       do k=1,npz
          if (dT_cond(k) .lt. 0.) then
             write(*,*) 'CHECKING EACH LAYER', k, L*dmcond(k), 5500.*dT_cond(k)*delp(k)
          endif
       enddo
    endif
    q_cond = 0.0_dp

    if (mod(tstep,1000) .eq. 0) then
       call q_sat(p, T_dash, qsats)
       write(*,*) 'RAINOUT, q, qsat'
       do k=1,npz
          if (abs(dT_cond(k)) .gt. 0) then
             write(*,*) k, q(k,1), qsats(k,1)
          endif
       enddo
    endif
  end subroutine condensation
  
end module lsc_mod
