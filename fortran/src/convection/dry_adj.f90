module dry_adj_mod
  
  use aqua_eos,   only: dim_output, interpolate_aqua_pt
  use phys, only: Rstar, H2O
  use params, only: dp, del_time, grav
  use atmosphere, only: get_mmw, get_cp, q_orig, mmw_dry, nqt
  
  implicit none
contains

  subroutine dry_adj(p, delp, T, q, dT_dry, tstep)
    !=========================================================================
    ! Dummy variables
    !=========================================================================

    integer, intent(in) :: tstep
    real(dp), intent(in),    dimension(:) :: p, delp    ! Pressure and pressure thickness
    real(dp), intent(in), dimension(:) :: T
    real(dp), intent(inout), dimension(:,:) :: q
    real(dp), intent(out),   dimension(:) :: dT_dry ! Change in dry temperature

    !=========================================================================
    ! Local variables
    !=========================================================================

    integer, parameter :: N_iter = 100
    integer :: npz,k,n,m, l,ltop
    real(dp) :: eos_pt(dim_output), grad, T_lift
    real(dp), dimension(size(p)) :: dt_lifts, T_new, T_old, cps, grads
    real(dp) :: num, denom, cp, mmw, mmw_up, grad_real
    real(dp) :: mass_tot, q_tot, hbef, haft, htest, h0, h00, htp

    npz = size(p)
    haft = 0.0_dp
    hbef = 0.0_dp
    h00 = 0.0
    do k=1,npz
       call get_cp(q(k,:), cp)
       cp = 5500.
       cps(k) = cp
       h00 = h00 + cps(k)*delp(k)*T(k)
    enddo
    
    if (mod(tstep, 1) .eq. 0) then
    !    write(*,*) 'INSIDE DRY ADJ, BEFORE'
    !    grads = 0.0
    !    do k=npz-1,1,-1
    !       grads(k) = (log(T(k+1)) - log(T(k)))/(log(p(k+1)) - log(p(k)))
    !       write(*,*) k, grads(k)
    !    enddo

       htp = 0.0_dp
       do k=1,npz
           htp = htp + delp(k)*q(k,1)
        enddo
    !    if (mod(tstep,1000) .eq. 0) write(*,*) 'enth calc 2324_PRE', T(23), t(24), T_new(23)
    !    if (mod(tstep,1000) .eq. 0) write(*,*) 'enth calc 2324_PRE', cps(23), cps(24), delp(23), delp(24) 

    !    write(*,*) '----------------------'
     endif
    T_new = T

    ! write(*,*) 'OUTSIDE DRY ADJ'
    ! do k=1,npz
    !    write(*,*) q(k,:)
    !     call get_mmw(q(k,:), mmw)
    !     write(*,*) k, mmw, T(k)*(1.e5/p(k))**(0.22)*mmw_dry(k)/mmw
    !  enddo

    outer: do n = 1,N_iter
       T_old = T_new
       do k=npz,10,-1

          if (q(k,1) .gt. 0.99_dp) then
             cp = H2O%cp
          else
             call get_cp(q(k,:), cp)
             call get_mmw(q(k,:), mmw)
             grad = Rstar/mmw/cp
             !grad = 0.22
          endif
          num = cps(k)*T_new(k)*delp(k)
          denom = cps(k)*delp(k)
          mass_tot = delp(k)
          q_tot = q(k,1)*delp(k)
          T_lift = T_new(k)
          dt_lifts(k) = 0.0_dp
          do l=k-1,10,-1
             if (q(l+1,1) .gt. 0.95_dp) then
                !eos_pt = interpolate_aqua_pt(p(l+1), T_new(l+1))
                grad = eos_pt(2)
                ! TODO: probably need Newton iteration to find enthalpy conserving adjustment
                cp = cps(k)!H2O%cp
                grad = 0.22
             else
                !cp = cps(k)
                call get_cp(q(l+1,:), cp)
                call get_mmw(q(l+1,:), mmw)
                grad = Rstar/mmw/cp
                !grad = 0.22
             endif

             
             T_lift = T_lift*exp(grad*log(p(l)/p(l+1)))
             dt_lifts(l) = T_lift - T_new(k)
             
             ! Use virtual temperature for comparison of temperature
             call get_mmw(q(l,:), mmw_up)
             !*mmw_dry(l+1)/mmw *mmw_dry(l)/mmw_up
             
             !if (T_lift*mmw_dry(l+1)/mmw .gt. T_new(l)*mmw_dry(l)/mmw_up) then
        !        write(*,*) 'DRY ADJ', l, T_lift*mmw_dry(l+1)/mmw, T_new(l)*mmw_dry(l)/mmw_up
             if (T_lift .gt. T_new(l)) then
                num = num + cps(l)*delp(l)*(T_new(l) - dt_lifts(l))
                denom = denom + cps(l)*delp(l)
                mass_tot = mass_tot + delp(l)
                q_tot = q_tot + q(l,1)*delp(l)
                ltop = l
             else
                ltop = l+1
                exit
             endif
          enddo
          !if (tstep .eq. 1000) write(*,*) k, ltop
          hbef = 0.0_dp
          if (ltop .lt. k) then
             do l=k,ltop,-1
                hbef = hbef + T(l)*delp(l)*cps(l)
                T_new(l) = min(max(num/denom + dt_lifts(l), T_new(l) - 50._dp), T_new(l) + 50._dp)
                
                ! Mix water perfectly through the dry layer
                q(l,1) = q_tot/mass_tot
                do m=2,nqt
                   q(l,m) = (1. - q(l,1))/(1 - q_orig(l,1)) * q_orig(l,m)
                enddo
             enddo

             haft = 0.0_dp
             htest = 0.0_dp
             h0 = 0.0_dp
             do l=1,npz
                h0 = h0 + cps(l)*T_new(l)*delp(l)
             enddo
             do l=k,ltop,-1
                haft = haft + cps(l)*T_new(l)*delp(l)
                htest = htest + cps(l)*(T_new(l) - T(l))*delp(l)
             enddo
!             if (mod(tstep,1000) .eq. 0) write(*,*) 'enth calc 2324', cps(23), cps(24), delp(23), delp(24) 
          endif
!          if (k .eq. 1) exit outer
       enddo
       if (maxval(abs(T_new - T_old)) < 1.e-8) exit
    enddo outer

    dT_dry = T_new - T
    h0 = 0.0
    haft = 0.0
    do k=1,npz
       
       haft = haft + cps(k)*dT_dry(k)*delp(k)/del_time/grav
       h0 = h0 + delp(k)*T_new(k)*cps(k)
       !if (mod(tstep,1000) .eq. 0) write(*,*) 'CHECK INSIDE', cps(k), delp(k), &
       !     dT_dry(k), del_time, grav
    enddo

if (mod(tstep, 1000) .eq. 0) then
       write(*,*) 'INSIDE DRY ADJ, AFTER'
       do k=1,npz
          write(*,'(I4, 6(1PE25.10))') k, T(k), T_new(k), &
               cps(k)*delp(k)*T_new(k) - cps(k)*delp(k)*T(k), cps(k)*delp(k)*dT_dry(k), dT_dry(k)
       enddo
       write(*,*) '----------------------'
       
       do k=1,npz-1
          call get_mmw(q(k,:), mmw_up)
          write(*,*) k, (log(T_new(k+1))  - log(T_new(k)))/(log(p(k+1)) - log(p(k))),&
               T(k)*(1.e5/p(k))**(0.22)*mmw_dry(k)/mmw_up
       enddo
       write(*,*) '------------------------'
    endif

    if (mod(tstep,1000) .eq. 0) then
       write(*,*) 'mass consv', htp
       htp = 0.0
       do k=1,npz
          htp = htp + cps(k)*T_new(k)*delp(k)
       enddo
       write(*,*) 'ENTHALPY CHANGE IN DRY_ADJ', haft, htp, h00, htp - h00
    endif
    
  end subroutine dry_adj
end module dry_adj_mod
