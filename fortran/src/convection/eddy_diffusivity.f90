module eddy_diff_mod
  use params, only: dp, grav, grav
  use atmosphere, only: get_mmw, get_cp, mmw_dry
  use phys, only: Rstar
  use adjust_mod, only: gradient
  use tables, only: phase_grad
  !use condense, only: q_sat
  use q_sat_test, only: q_sat
  !use lapack, only: DGSTV
  implicit none

contains
  subroutine get_kzz(conv_mask, Kzz)
    integer, intent(in) :: conv_mask(:)
    real(dp), intent(out) :: Kzz(:)

    integer :: npz, k
    npz = size(conv_mask)
    ! Test with constant value of Kzz
    Kzz = 0.0_dp
    do k=1,npz
       if (conv_mask(k) > 0) cycle
       Kzz(k) = 0.3
    enddo
    
    
  end subroutine get_kzz
 

  subroutine do_explicit_diff(Kzz, rho, delp, pf, Tf, q, flux_diff, turb_flux, tstep, dry_mask)
    
    ! Input variables
    real(dp), intent(in) :: Kzz(:),rho(:), pf(:)
    integer, intent(in) :: dry_mask(:)
    real(dp), intent(in)    :: delp(:) ! Layer thicknesses
    !real(dp), intent(in)  :: source(:) ! dT/dt of all explicit terms

    integer, intent(in) :: tstep
    ! Inout variables
    real(dp), intent(inout) :: Tf(:), q(:,:), flux_diff(:), turb_flux(:)

    ! Edge values of rho and delp (excluding first and last edges)
    real(dp), dimension(size(Tf)-1) :: rho_edge, delp_edge, Kzz_edge, th_edge,mmw_edge, cp_edge, pf_edge, Tf_edge, &
         grad_ad_edge, grad_edge, grad_mu_edge, qcrits_edge, q_edge, stabs

    ! Diffused quantity (potential temperature)
    real(dp), dimension(size(Tf)) :: theta, qcrits, grads
    real(dp), dimension(size(q,1),size(q,2)) :: qsats
    real(dp), dimension(size(Tf)) :: mmw, cp
    logical, dimension(size(Tf)):: moist
    !real(dp), dimension(size(Tf)) :: dT_dt
    real(dp), dimension(0:size(Tf)):: flux_arr
    real(dp) :: p0, dlnpsatdlnt, stab
    
    
    ! Local variables
    integer :: k,npz, nq, pi! Z index, N_z and N_tracer

    npz = size(pf)
    p0 = 1.e5

    call q_sat(pf, Tf, qsats)
    do k=1,npz
       call get_mmw(q(k,:), mmw(k))
       call get_cp(q(k,:), cp(k))
       if (q(k,1) .gt. qsats(k,1) - 1.e-10 ) then
          call gradient(pf(k), Tf(k), cp(k), mmw_dry(k), grads(k), dlnpsatdlnt)
          qcrits(k) = 1._dp/(1. - mmw_dry(k)/18.0_dp)/dlnpsatdlnt
       else
          ! TODO: change for supercritical
          grads(k) = Rstar/mmw(k)/cp(k)
       endif
    enddo

    
    delp_edge = 0.0_dp
    rho_edge = 0.0_dp
    Kzz_edge = 0.0_dp
    flux_arr = 0.0_dp
    
    do k=1,npz-1
       cp_edge(k) = 0.5_dp*(cp(k) + cp(k+1))
       pf_edge(k) = 0.5_dp*(pf(k+1) + pf(k))
       rho_edge(k) = 0.5_dp * (rho(k+1) + rho(k))
       Kzz_edge(k) = 0.5_dp*(Kzz(k+1) + Kzz(k))
       Tf_edge(k) = 0.5_dp*(Tf(k) + Tf(k+1))
       qcrits_edge(k) = 0.5_dp*(qcrits(k+1) + qcrits(k))
       q_edge(k) = 0.5_dp*(q(k+1,1) + q(k,1))
       grad_ad_edge(k) = 0.5_dp*(grads(k+1) +  grads(k))
       grad_edge(k) = log(Tf(k+1)/Tf(k))/log(pf(k+1)/pf(k))
       grad_mu_edge(k) = log(mmw(k+1)/mmw(k))/log(pf(k+1)/pf(k))

       if (q(k,1) .gt. qsats(k,1) - 1.e-10 .or. q(k+1,1) .gt. qsats(k+1,1) - 1.e-10) then
          stabs(k) = max(0._dp, (grad_ad_edge(k) - grad_edge(k)) * (qcrits_edge(k) - q_edge(k)))
          moist(k) = .true.
       else
          stabs(k) = max(0._dp, grad_ad_edge(k) - grad_edge(k) + grad_mu_edge(k))
          moist(k) = .false.
       endif
       flux_arr(k) = Kzz_edge(k)*rho_edge(k)**2._dp * grav * cp_edge(k) * Tf_edge(k)/pf_edge(k) * stabs(k)
       !flux_arr(k) = Kzz_edge(k)*rho_edge(k)**2._dp * grav *(theta(k+1) - theta(k))/delp_edge(k)
    enddo

    !dT_dt(1) = flux_arr(1) * exp(R(k)/cp(k)*log(pf(k)/p0))/delp(1)
    !dT_dt(npz) = -flux_arr(npz-1)*exp(R(npz)/cp(npz)*log(pf(npz)/p0))/delp(npz)

    do k=1,npz
       !dT_dt(k) = (flux_arr(k) - flux_arr(k-1))*exp(R(k)/cp(k)*log(p(k)/p0))&
       !     /cp(k)/delp(k)
       flux_diff(k) = (flux_arr(k) - flux_arr(k-1))
       turb_flux(k) = flux_arr(k-1)
    enddo
    turb_flux(npz+1) = flux_arr(npz)

    if (mod(tstep,500) .eq. 0) then
       write(*,*) 'k, rho(k), stabs(k), flux_arr(k)'
       do k=1,npz-1
          write(*,*) k,rho(k), stabs(k), grad_ad_edge(k), grad_edge(k), q_edge(k), qcrits(k), moist(k)
       enddo
       write(*,*) 'rho(npz)', rho(npz)
    endif
    
    ! Enthalpy fixer
  end subroutine do_explicit_diff
  
  subroutine do_turb_diff(Kzz, delt, rho, source, delp, pf, Tf, q, tstep)
    ! Use the prescription of Leconte 2024 to calculate the eddy diffusivity of vapour etc. within
    ! the layers. Can then decide whether I want to use eddy diffusivity for the heat and vapour?
    ! Leconte using it for the vapour and 'entropy' (read: theta_v?)

    
    ! Input variables
    real(dp), intent(in) :: Kzz(:), delt(:), rho(:), pf(:)
    real(dp), intent(in)    :: delp(:) ! Layer thicknesses
    real(dp), intent(in)  :: source(:) ! dT/dt of all explicit terms
    integer, intent(in) :: tstep
    
    ! Inout variables
    real(dp), intent(inout) :: Tf(:), q(:,:)
    
    ! Local variables
    integer :: k,npz, nq ! Z index, N_z and N_tracer

    ! Inputs for LAPACK tridiagonal solver
    integer :: nrhs, info
    ! Tridiagonal matrix components + source term (D)
    ! To solve equation A_k * q_{k-1} + B_k * q_k + C_k * q_{k+1} = D_k
    real(dp), dimension(size(Tf)) :: B
    real(dp), dimension(size(Tf) - 1) :: A, C
    real(dp), dimension(size(Tf), 1) :: D
    ! Edge values of rho and delp (excluding first and last edges)
    real(dp), dimension(size(Tf)-1) :: rho_edge, delp_edge, Kzz_edge

    ! Diffused quantity (potential temperature)
    real(dp), dimension(size(Tf)) :: theta
    real(dp), dimension(size(Tf)) :: mmw, cp

    real(dp) :: p0 ! Reference pressure for potential temperature

    p0 = 1.e5
    npz = size(Tf)
    nq = size(q,2)

    do k=1,npz
       call get_mmw(q(k,:), mmw(k))
       call get_cp(q(k,:), cp(k))
       !theta = T(k)*exp(R(k)/cp(k)*log(p0/pf(k)))
    enddo

    delp_edge = 0.0_dp
    rho_edge = 0.0_dp
    Kzz_edge = 0.0_dp
    do k=1,npz-1
       delp_edge(k) = pf(k+1) - pf(k)
       rho_edge(k) = 0.5_dp * (rho(k+1) + rho(k))
       Kzz_edge(k) = 0.5_dp*(Kzz(k+1) + Kzz(k))
    enddo

    A = 0.0_dp ; B= 0.0_dp ; C = 0.0_dp 
   do k=1,npz-1
      A(k) = - grav**2 * Kzz_edge(k) * rho_edge(k)**2/delp_edge(k)* exp(Rstar/mmw(k+1)/cp(k+1)*log(p0/pf(k+1)))
      C(k) = - grav**2 * Kzz_edge(k) * rho_edge(k)**2/delp_edge(k)* exp(Rstar/mmw(k)/cp(k)*log(p0/pf(k)))
   enddo

   B(1) = delp(1)/delt(1)*exp(Rstar/mmw(1)/cp(1)*log(p0/pf(1))) - C(1)
   B(npz) = delp(npz)/delt(npz)*exp(Rstar/mmw(npz)/cp(npz)*log(p0/pf(npz))) - A(npz-1)
    do k=2,npz-1
       B(k) = delp(k)/delt(k)*exp(Rstar/mmw(k)/cp(k)*log(p0/pf(k))) - A(k-1) - C(k)
    enddo
    
    D(1:npz,1) = (Tf(1:npz) +source(1:npz))* delp(1:npz)/delt(1:npz)*exp(Rstar/mmw(1:npz)/cp(1:npz)*log(p0/pf(1:npz)))
    ! do k=1,npz-1
    !    write(*,*) A(k), B(k), C(k), D(k,1)
    ! enddo
    ! write(*,*) '----------------------------'
    ! Solve tridiagonal matrix with LAPACK
    nrhs = 1
    call dgtsv(npz, nrhs, &
         A, B, C, D, &
         npz, info)

    if (info .ne. 0) then
       write(*,*) 'ERROR in diagonal solver'
       write(*,*) 'ERROR CODE: ', info
       stop
    endif

    ! Update T
    if (mod(tstep,500) .eq. 0) then
       do k=1,npz
          write(*,*) D(k,1) - Tf(k), source(k), Tf(k)
       enddo
    endif
    
          
       
    Tf(1:npz) = D(1:npz, 1)
!    do k=1,npz
!       write(*,*) Tf(k)
!    enddo
!    write(*,*) '-----------------------------------------------'
    
    !T(1:npz) = theta(1:npz)*exp(R(1:npz)/cp(1:npz)*log(p(1:npz)/p0)

    ! Try with the explicit solver instead
    Tf(1:npz) = theta(1:npz)
  end subroutine do_turb_diff
end module eddy_diff_mod
