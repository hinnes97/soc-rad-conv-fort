module moist_adj_mod
  use condense, only : q_sat, dqsat_dt
  use adjust_mod, only: gradient
  use params, only : inhibited, dp
  use atmosphere, only: get_cp, q_orig, cp_dry, mmw_dry, nqt, get_cp
  use phys, only: H2O
  use tables, only: find_var_lin, phase_grad, lheat
  

  implicit none

contains

  subroutine moist_adj(p, delp, T, q, q_cond, grads, tstep, dT_moist)
    !=========================================================================
    ! Dummy variables
    !=========================================================================
    
    real(dp), intent(in),    dimension(:) :: p, delp, T    ! Pressure and pressure thickness
    real(dp), intent(inout), dimension(:) :: grads,  q_cond   ! Ad. gradient, temp
    real(dp), intent(inout), dimension(:,:) :: q
    real(dp), intent(out), dimension(:) :: dT_moist
    integer, intent(in) :: tstep

    !=========================================================================
    ! Local variables
    !=========================================================================

    integer :: npz,k,n,m, a,ltop, l_rec, k_rec
    integer, parameter :: N_iter=1000 ! Number of up-down iterations
    logical :: condition, super, moist

    real(dp) :: qsats(size(p), size(q,2)), kappa, qmin, f, mmw, dqsatdts(size(p)), dT_lifts(size(p)), &
         q_old(size(p)), T_new(size(p)), T_old(size(p))
    real(dp) :: qcrit, pfact, temp
    real(dp) :: qcrits(sizE(p))
    real(dp), parameter          :: delta = 0.0001 ! Small number speeds up convergence
    real(dp) :: h_bef, cp, L, T_lift, q_extra, num, denom, cps(size(p)), l_new, Ls(size(p))

    real(dp) ::dlnpsatdlnt, grad, h0, hend, m_cond, test
    ! Do moist adjustment all in one sweep, using analytic approximations for all the derivatives
    ! for speed. Should speed up moist convection a lot!

    ! Get vertical dimensions
    
    npz = size(p)
    nqt = size(q,2)


    l_rec = -1000
    k_rec = -1000
    ! Moist enthalpy before
    ! if (mod(tstep, 1) .eq.  0 ) then
    !    write(*,*) 'MOIST ADJ'
    !    h0  =0.0
    !    do k=1,npz
    !       call get_cp(q(k,:), cps(k))
    !       call find_var_lin(T(k), lheat, L)
    !       L = max(L, 0._dp)
    !       h0 = h0 + delp(k)*(cps(k)*T(k) + L*q(k,1))
    !       write(*,*) k, delp(k), cps(k), T(k), q(k,1), L
    !    enddo
    ! endif

    dT_moist = 0.0_dp
    q_old = q(:,1)

    T_new = T(:)
    ! Newton iteration loop
    m_cond = 0.0_dp
    outer: do n=1, N_iter
       T_old = T_new    
       ! Vertical loop (bottom up)
       h_bef = 0.0_dp
       num = 0.0_dp
       denom = 0.0_dp
       q_extra = 0.0_dp
       do k=npz,10, -1
          call q_sat(p, T_new, qsats)
          if (qsats(k,1) .gt. 1._dp) cycle
          call dqsat_dt(p, T, dqsatdts, qsats)
          do a=1,npz
             call get_cp(q(k,:), cps(k))
             cps(k) = 5500.
             call find_var_lin(T(a), lheat, Ls(a))
             Ls(a) = 2.6e6
             Ls(a) = max(Ls(a), 0._dp)
          enddo
          

          num = ((cps(k) + Ls(k)*dqsatdts(k))*T_new(k) + Ls(k)*(q(k,1) - qsats(k,1)))*delp(k)
          denom =(cps(k) + Ls(k)*dqsatdts(k)) * delp(k)
          q_extra = delp(k)*(q(k,1) - qsats(k,1))
          T_lift = T_new(k)
          dt_lifts(k) = 0.0_dp
          ! Lifting from layer k
          do a=k-1,10,-1
             
             q_extra = q_extra + (q(a,1) - qsats(a,1))*delp(a)
             ! Find lifted temperature on adiabat
             
             call gradient(p(a+1),T_new(a+1), cp_dry(a+1), mmw_dry(a+1), grad)
             grads(a+1) = grad
             call find_var_lin(T_new(a+1), phase_grad, dlnpsatdlnt)

             if (inhibited) then
                qcrit = 1./(1._dp - mmw_dry(a+1)/H2O%mmw)/dlnpsatdlnt
             else
                qcrit = 1.e10
             endif

             
             T_lift = T_lift*exp(grad*log(p(a)/p(a+1)))
             dt_lifts(a) = T_lift - T_new(k)

             
             if (T_lift .gt. T_new(a) .and. &
                  q_extra .gt. 0      .and. &
                  q(a+1,1) .lt. qcrit .and. &
                  qsats(a,1) .lt. 1.0_dp .and. &
                  qsats(a+1,1) .lt. 1.0_dp) then

!                write(*,*) 'MOIST CONV, k', k, qsats(k,1)
                !call get_cp(q(a,:), cp)
             
                !call find_var_lin(T(a), lheat, L)
             
                !L = max(L, 0._dp)
                ! Add to quantities needed for adjustment

                
                num = num + ((cps(a) + Ls(a)*dqsatdts(a))*(T_new(a) - dt_lifts(a))+ Ls(a)*(q(a,1) - qsats(a,1)))* delp(a)
                denom = denom + (cps(a) + Ls(a)*dqsatdts(a))*delp(a)
!                write(*,*) 'Hello', cps(a), Ls(a), dt_lifts(a), T_new(a), &
                !                     grad, qsats(a,1), q(a,1), dqsatdts(a)*(num/denom+dt_lifts(a) - T_new(a))
                ltop = a
             else

                ltop = a+1
                exit
             endif

          enddo

          ! Now we have layers from k to ltop to adjust
          if (ltop .lt. k) then
             k_rec = max(k_rec, k)
             l_rec = max(ltop, l_rec)

             h0 = 0.0
             hend = 0.0
             test = 0.0
             m_cond = 0.0
             do a=k,ltop,-1
                h0 = h0 + delp(a)*(q(a,1) + q_cond(a))
                !dT_moist(l) = T(k) + dt_lifts(l) - T(l)
                dT_moist(a) = num/denom + dt_lifts(a) - T_new(a)
                T_new(a) = num/denom + dt_lifts(a) 
                m_cond = m_cond + (q(a,1) - (qsats(a,1) + dqsatdts(a)*dT_moist(a)))*delp(a)
                !hend = hend + delp(a)*(q(a,1) - (qsats(a,1) + dqsatdts(a)*dT_moist(a)))*Ls(a)
                
                q(a,1) = qsats(a,1) + dqsatdts(a)*dT_moist(a)
                
                !h0  = h0 + delp(a)*cps(a)*(T_new(a) - T_old(a))
                !write(*,*) 'during', a, T_new(a) - T_old(a), q_cond(a)
                do m=2,nqt
                   q(a,m) = (1. - q(a,1))/(1. - q_orig(a,1)) * q_orig(a,m)
                enddo
!                write(*,*) a, q(a,1), T_new(a) - dT_moist(a), T_new(a)
             enddo
             !hend = hend + m_cond
             q_cond(ltop:k) = q_cond(ltop:k) + m_cond/sum(delp(ltop:k))
             do a=k,ltop,-1
                test = test + m_cond/sum(delp(ltop:k))*delp(a)
                hend = hend + (q(a,1)+q_cond(a))*delp(a)
             enddo
             if (mod(tstep,1000) .eq. 0) then
                write(*,*) tstep, 'MASS CONSV CHECK MOIST', h0,hend, m_cond, test
             endif
             !write(*,*) 'MOIST ENTH CHECK', h0, hend, ltop, k
          endif
          !if (k .eq. 1) exit outer 
       enddo
          
       if (maxval(abs(T_new - T_old)) .lt. 1.e-8 ) exit
    enddo outer

    dT_moist = T_new - T

    if (mod(tstep,1000) .eq. 0) then
       write(*,*) 'MOIST ADJ CHECKER'
       write(*,*) 'k, q, qsast, grad, real_grad'
       call q_sat(p, T_new, qsats)
       do k=2,npz
          call gradient(p(k),T_new(k), cp_dry(k), mmw_dry(k), grad)
          call find_var_lin(T_new(k), phase_grad, dlnpsatdlnt)
          
          write(*,*) k, q(k,1), qsats(k,1), grad, (log(T_new(k)) - log(T_new(k-1)))/(log(p(k)) - log(p(k-1))), &
               1./(1._dp - mmw_dry(k)/H2O%mmw)/dlnpsatdlnt, T(k), dT_moist(k)
          
          !write(*,*) p(k), T_new(k), cp_dry(k), mmw_dry(k)
       enddo
    endif

    if (any(abs(dT_moist) .gt. 0.) ) then
!       write(*,*) 'MOIST ENTH', l_rec, k_rec
       h0  =0.0
       hend=0.0
       do m=l_rec, k_rec
          !call get_cp(q(k,:), cps(k))
          !call find_var_lin(T(m), lheat, L)
          !call find_var_lin(T_new(m), lheat, L_new)
          !L = max(L, 0._dp)
          !L_new = max(L_new, 0.0_dp)
          h0  = h0 + delp(m)*cps(m)*(T_new(m) - T(m))
          hend = hend + Ls(m)*q_cond(m)*delp(m)

!          write(*,*) 'after', m, T_new(m) - T(m), q_cond(m)
!          write(*,*) k, T(m), T_new(m), q_old(m),q(m,1)
       enddo

       if (mod(tstep,1000) .eq. 0) then
          write(*,*) 'COMPARING MOIST ENTH TEMP VS COND', tstep, h0, hend, l_rec, k_rec
          write(*,*) sum((cps(l_rec:k_rec)*T_new(l_rec:k_rec) + Ls(l_rec:k_rec)*q(l_rec:k_rec,1))*delp(l_rec:k_rec)),  &
               sum((cps(l_rec:k_rec)*T(l_rec:k_rec) + Ls(l_rec:k_rec)*q_old(l_rec:k_rec))*delp(l_Rec:k_rec))
       
       write(*,*) 'q, qsat, grad, real_grad'
       call q_sat(p, T_new, qsats)
       do k=l_rec+1,k_rec
          write(*,*) k, q(k,1), qsats(k,1), grads(k), (log(T_new(k)) - log(T_new(k-1)))/(log(p(k)) - log(p(k-1)))
       enddo
    endif
     endif
  
   end subroutine moist_adj
end module moist_adj_mod
