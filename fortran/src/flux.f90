module flux_mod

  use params, only: dp, rad_scheme
  use radiation_Kitzmann_noscatt, only: Kitzmann_TS_noscatt
  implicit none
contains
  subroutine get_fluxes(nf, ne, Tf, pf, Te, pe, tau_IR, tau_V, &
       net_F, mu_s, Finc, Fint, olr)

    !! Input variables
    integer, intent(in) :: nf, ne                         ! Number of layers, levels (lev = lay + 1)
    real(dp), dimension(:), intent(in) :: Tf, pf           ! Temperature [K], pressure [pa] at layers
    real(dp), dimension(:), intent(in) :: Te, pe           ! pressure [pa] at levels
    real(dp), dimension(:), intent(in) :: tau_V, tau_IR  ! V and IR band cumulative optical depth at levels
    real(dp), intent(in) :: Finc, mu_s                        ! Incident flux [W m-2] and cosine zenith angle
    real(dp), intent(in) :: Fint                              ! Internal flux [W m-2]
    
    !! Output variables
    real(dp), dimension(ne), intent(out) :: net_F
    real(dp), intent(out) :: olr

    select case(rad_scheme)
    case(1)
       call Kitzmann_TS_noscatt(nf, ne, Te, pe, tau_IR, tau_V, &
            net_F, mu_s, Finc, Fint, olr)
    end select    

  end subroutine get_fluxes
end module flux_mod
