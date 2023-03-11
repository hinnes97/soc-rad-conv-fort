!! EKH-Lee - Oct 2020
! Two-stream method following the Helios-r method without scattering (Kitzmann et al. 2020)
! Uses the method of short characteristics (Olson & Kunasz 1987) with linear interpolants.
! Properties: Fast, efficenct, flexible and stable
! Downsides: No scattering - unknown behaviour at very low optical depth (< 1e-6)
!!
module sc_split_mod
  use, intrinsic :: iso_fortran_env
  use params, only : sb, pi, twopi, Finc, Fint,  dp, kappa_sw, kappa_lw, surface, moist_rad, kappa_q
  use tau_mod, only: calc_tau
  implicit none

  !! Gauss quadrature variables - here you can comment in/out groups of mu values for testing

  !! single angle diffusion factor approximation - typically 1/1.66
  ! integer, parameter :: nmu = 1
  ! real(dp), dimension(nmu), parameter :: uarr = (/1.0_dp/1.66_dp/)
  ! real(dp), dimension(nmu), parameter :: w = (/1.0_dp/)
  ! real(dp), dimension(nmu), parameter :: wuarr = uarr * w


  !! Legendre quadrature for 2 nodes
  !integer, parameter :: nmu = 2
  !real(dp), dimension(nmu), parameter :: uarr = (/0.21132487_dp, 0.78867513_dp/)
  !real(dp), dimension(nmu), parameter :: w = (/0.5_dp, 0.5_dp/)
  !real(dp), dimension(nmu), parameter :: wuarr = uarr * w

  integer, parameter :: nmu = 1
  real(dp), dimension(nmu), parameter :: uarr = (/1.0_dp/)
  real(dp), dimension(nmu), parameter :: w = (/1.0_dp/)
  real(dp), dimension(nmu), parameter :: wuarr = uarr * w

  

  !! Lacis & Oinas (1991) numerical values - DOES NOT WORK CORRECTLY FOR SOME REASON
  ! integer, parameter :: nmu = 3
  ! real(dp), dimension(nmu), parameter :: uarr = (/0.1_dp, 0.5_dp, 1.0_dp/)
  ! real(dp), dimension(nmu), parameter :: wuarr = (/0.0433_dp, 0.5742_dp, 0.3825_dp/)

  !! Legendre quadrature for 4 nodes
  ! integer, parameter :: nmu = 4
  ! real(dp), dimension(nmu), parameter :: uarr = &
  !   & (/0.06943184_dp, 0.33000948_dp, 0.66999052_dp, 0.93056816_dp/)
  ! real(dp), dimension(nmu), parameter :: w = &
  !   & (/0.17392742_dp, 0.32607258_dp, 0.32607258_dp, 0.17392742_dp/)
  ! real(dp), dimension(nmu), parameter :: wuarr = uarr * w

  !! 5 point EGP quadrature values
  ! integer, parameter :: nmu = 5
  ! real(dp), dimension(nmu), parameter :: uarr = &
  !   &(/0.0985350858_dp, 0.3045357266_dp, 0.5620251898_dp, 0.8019865821_dp, 0.9601901429_dp/)
  ! real(dp), dimension(nmu), parameter :: wuarr = &
  !   & (/0.0157479145_dp, 0.0739088701_dp, 0.1463869871_dp, 0.1671746381_dp, 0.0967815902_dp/)

  public :: sc_split
  private :: lw_grey_updown_linear, sw_grey_down, linear_log_interp

contains

  subroutine sc_split(nlay, nlev, Te, pe, Tf, &
          &  net_F, mu_s,Ts, olr, q, lw_up, lw_down, sw_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev                         ! Number of layers, levels (lev = lay + 1)
    !real(dp), dimension(nlay), intent(in) :: Tl, pl           ! Temperature [K], pressure [pa] at layers
    real(dp), dimension(nlev), intent(in) :: Te, pe           ! pressure [pa] at levels
    real(dp), dimension(nlay), intent(in) :: Tf
    real(dp), intent(in) :: mu_s                        ! Incident flux [W m-2] and cosine zenith angle
    real(dp), intent(in) :: Ts ! Surface temp [K]
    real(dp), intent(in) :: q(:)

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: net_F, lw_down, lw_up, sw_down
    real(dp), intent(out) :: olr

    !! Work variables
    integer :: i
    real(dp), dimension(nlev) :: be
    real(dp), dimension(nlay) :: bf
    real(dp), dimension(nlev) :: sw_up
    real(dp), dimension(nlev) :: lw_net, sw_net, tau_IRe, tau_Ve, tau_dry, tau_wet
    real(dp), dimension(nlay) :: q_temp

    q_temp = 1.0
    tau_IRe = 0.0
    tau_Ve = 0.0
    tau_dry = 0.0
    tau_wet = 0.0

    if (moist_rad) then
       call calc_tau(kappa_q, pe, q, 0.5_dp, tau_wet)
       call calc_tau(kappa_lw, pe, 1-q, 0.5_dp, tau_dry)
       call calc_tau(kappa_sw, pe, q_temp, 0.5_dp, tau_Ve)
       tau_IRe = tau_wet + tau_dry
    else
       call calc_tau(kappa_lw, pe, q_temp, 0.5_dp, tau_IRe)
       call calc_tau(kappa_sw, pe, q_temp, 0.5_dp, tau_Ve)
    endif
    

    !do i=1,nlev
    !   tau_IRe(i) = 256*pe(i)/pe(nlev)
    !   tau_Ve(i) = 256*0.04*pe(i)/pe(nlev)
    !enddo

!    write(*,*) tau_IRe(nlev), tau_Ve(nlev)
    !! Shortwave fluxes
    sw_down(:) = 0.0_dp
    sw_up(:) = 0.0_dp
    ! Calculate sw flux
    if (mu_s > 0.0_dp) then
      call sw_grey_down(nlev, Finc, tau_Ve, sw_down, mu_s)
    end if

    !! Long wave two-stream fluxes
    ! Blackbody fluxes (note divide by pi for correct units)
    be(:) = sb * Te(:)**4
    bf(:) = sb* Tf(:)**4
    ! Calculate lw flux
    call lw_grey_updown_linear(nlay, nlev, be,bf, Ts, tau_IRe, lw_up, lw_down)

    ! Net fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)
    sw_net(:) = sw_up(:) - sw_down(:)
    net_F(:) = lw_net(:) + sw_net(:)

    ! Internal flux at lowest level
    !net_F(nlev) = Fint
    ! olr is upward longwave flux
    olr = lw_up(1)

    !write(*,*) 'TSRAD lwup(1), lwdown(1), swup(1), swdown(1)', lw_up(1), lw_down(1), sw_up(1), sw_down(1)
  end subroutine sc_split

  subroutine lw_grey_updown_linear(nlay, nlev, be, bf, Ts, tau_IRe, lw_up, lw_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: Ts
    real(dp), dimension(nlev), intent(in) :: be
    real(dp), dimension(nlay), intent(in) :: bf
    real(dp), dimension(nlev), intent(inout) :: tau_IRe

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: lw_up, lw_down

    !! Work variables and arrays
    integer :: k, m
    real(dp), dimension(nlay) :: dtau
    real(dp), dimension(2*nlay) :: dtau_half, edel
    real(dp), dimension(2*nlay + 1) :: b_comb, lw_up_g_half, lw_down_g_half
    real(dp) :: del, e0i, e1i, e1i_del
    real(dp), dimension(2*nlay) :: Am, Bm
    real(dp), dimension(nlev) :: lw_up_g, lw_down_g

    !! Calculate dtau in each layer
    do k = 1, nlay
       dtau(k) = tau_IRe(k+1) - tau_IRe(k)
       dtau_half(2*k-1) = 0.5*dtau(k)
       dtau_half(2*k) = 0.5*dtau(k)
       b_comb(2*k-1) = be(k)
       b_comb(2*k) = bf(k)
    end do
    b_comb(2*nlay+1) = be(nlev)

    ! Zero the total flux arrays
    lw_up(:) = 0.0_dp
    lw_down(:) = 0.0_dp

    !! Start loops to integrate in mu space
    do m = 1, nmu

      !! Prepare loop
      do k = 1, nlay*2
        ! Olson & Kunasz (1987) linear interpolant parameters
        del = dtau_half(k)/uarr(m)
        edel(k) = exp(-del)
        e0i = 1.0_dp - edel(k)
        e1i = del - e0i

        e1i_del = e1i/del

        Am(k) = e0i - e1i_del ! Am(k) = Gp(k), just indexed differently
        Bm(k) = e1i_del ! Bm(k) = Bp(k), just indexed differently
      end do

      !! Begin two-stream loops
      !! Peform downward loop first
      ! Top boundary condition - 0 intensity downward from top boundary
      lw_down_g(1) = 0.0_dp
      lw_down_g_half(1) = 0.0_dp
      do k = 1, nlay*2
        lw_down_g_half(k+1) = lw_down_g_half(k)*edel(k) + Am(k)*b_comb(k) + Bm(k)*b_comb(k+1) ! TS intensity
      end do

      !! Peform upward loop
      ! Lower boundary condition - planck function intensity upward from lowest level
      if (surface) then
         lw_up_g_half(2*nlay+1) = sb*Ts**4
      else
         !lw_up_g(nlev) = be(nlev)
         lw_up_g_half(2*nlay+1) = Fint/pi + lw_down_g(nlev)
      endif
      
      do k = 2*nlay, 1, -1
        lw_up_g_half(k) = lw_up_g_half(k+1)*edel(k) + Bm(k)*b_comb(k) + Am(k)*b_comb(k+1) ! TS intensity
      end do

      do k=1, nlev
         lw_down_g(k) = lw_down_g_half(2*k-1)
         lw_up_g(k) = lw_up_g_half(2*k-1)
      enddo
      
        
      !! Sum up flux arrays with Gauss weights and points for this mu stream
      lw_down(:) = lw_down(:) + lw_down_g(:) * wuarr(m)
      lw_up(:) = lw_up(:) + lw_up_g(:) * wuarr(m)
      
    end do

    ! Convert to flux by * 2pi
    lw_down(:) = lw_down(:)
    lw_up(:) =  lw_up(:)


  end subroutine lw_grey_updown_linear

  subroutine sw_grey_down(nlev, solar, solar_tau, sw_down, mu)
    implicit none

    integer, intent(in) :: nlev
    real(dp), intent(in) :: solar, mu
    real(dp), dimension(nlev), intent(in) :: solar_tau
    real(dp), dimension(nlev), intent(out) :: sw_down

    
    sw_down(:) = solar * mu * exp(-solar_tau(:)/mu)
    !write(*,*) solar, mu, solar_tau(1), sw_down(1)
  end subroutine sw_grey_down

  ! Perform linear interpolation in log10 space
  pure subroutine linear_log_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp) :: lxval, ly1, ly2, lx1, lx2
    real(dp), intent(out) :: yval
    real(dp) :: norm

    lxval = log10(xval)
    lx1 = log10(x1); lx2 = log10(x2)
    ly1 = log10(y1); ly2 = log10(y2)

    norm = 1.0_dp / (lx2 - lx1)

    yval = 10.0_dp**((ly1 * (lx2 - lxval) + ly2 * (lxval - lx1)) * norm)

  end subroutine linear_log_interp
  
end module sc_split_mod
