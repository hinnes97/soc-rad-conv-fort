module flux_mod

  use params, only: dp, sb, invert_grid, moist_rad, surface, semi_grey_scheme, A_s, rad_scheme
  
  use socrates_interface_mod, only: run_socrates
  use radiation_Kitzmann_noscatt, only: Kitzmann_TS_noscatt
  use toon_mod, only : toon_driver
  use sc_split_mod, only: sc_split
  use ts_short_char_bezier_mod, only: ts_short_char_Bezier

  
  implicit none
contains
  subroutine get_fluxes(nf, ne, Tf, pf, Te, pe, &
       net_F, mu_s, Finc, Fint, olr, q, Ts, fup, fdn, s_dn, s_up)
    !! Input variables
    integer, intent(in) :: nf, ne                         ! Number of layers, levels (lev = lay + 1)
    real(dp), dimension(:), intent(in) :: Tf, pf           ! Temperature [K], pressure [pa] at layers
    real(dp), dimension(:), intent(in) :: Te, pe           ! pressure [pa] at levels
    real(dp), intent(in) :: Finc, mu_s                        ! Incident flux [W m-2] and cosine zenith angle
    real(dp), intent(in) :: Fint                              ! Internal flux [W m-2]
    real(dp), intent(in) :: q(:)
    real(dp), intent(in) :: Ts                            ! Surface temperature
    
    !! Output variables
    real(dp), dimension(ne), intent(out) :: net_F, fup, fdn
    real(dp), intent(out) :: olr, s_dn(:), s_up(:)

    integer :: i
    
    real(dp) :: rad_lat, rad_lon, t_surf_in, albedo_in, net_surf_sw_down, surf_lw_down
    real(dp), dimension(size(Tf)) :: q_in, temp_tend, h2_in, ch4_in, he_in
    

    select case(rad_scheme)

       case ('socrates')
    rad_lat = 0._dp
    rad_lon = 0._dp
    !q_in = 0.001_dp*9._dp
    !q_in = 0.1_dp
    !ch4_in =  0.1_dp
    !q_in = 0.01_dp
    !h2_in =  1._dp - q - ch4_in!1._dp - q_in
    !h2_in = 0.9_dp*h2_in
    !he_in = 0.1_dp*h2_in
    
    !h2_in =0.8_dp
    !he_in = h2_in*0.1_dp
    !h2_in = h2_in*0.9_dp

    ch4_in = 0.0_dp
    q_in = q
    h2_in = 1-q-ch4_in
    
    he_in =h2_in*0.25188_dp
    h2_in =h2_in*0.74812_dp
    
    albedo_in = A_s
    t_surf_in = Ts

!    q_in = 0.2
!    h2_in = 0.8
!    he_in = 0.0
    !write(*,*) '-----------------------------------------------------------------------------'
    call run_socrates(rad_lat, rad_lon, Tf, q_in, h2_in, ch4_in, he_in,t_surf_in, pf, pe, pf, pe, albedo_in, &
         temp_tend, net_surf_sw_down, surf_lw_down, net_F, fup, fdn, s_up, s_dn)
    olr = fup(1)

    case('short_char')
       !call Kitzmann_TS_noscatt(nf, ne, Te, pe, &
       !     net_F, mu_s, Ts, olr, q, fup, fdn, s_dn)
       !write(*,*) 'short char', net_F
       call sc_split(nf, ne, Te, pe, Tf, &
            net_F, mu_s, Ts, olr, q, fup, fdn, s_dn)
!write(*,*) 'sc_split', net_F

    case('short_char_elsie')
       call ts_short_char_Bezier(.true., surface, nf, ne, Ts, Tf, pf, pe, &
            net_F, mu_s, olr, fup, fdn, s_dn, s_up, 0.0_dp, q)
    case('toon')
       if ((moist_rad .and. surface)) then
          call toon_driver(Te, pe, net_F, olr, s_dn, fdn, fup, q, Ts)
       else if (moist_rad) then
          call toon_driver(Te, pe, net_F, olr, s_dn, fdn, fup, q)
       else if (surface) then
          call toon_driver(Te, pe, net_F, olr, s_dn, fdn, fup, Ts=Ts)
       else
          call toon_driver(Te, pe, net_F, olr, s_dn, fdn, fup)
       endif
    end select

       
  end subroutine get_fluxes
end module flux_mod
