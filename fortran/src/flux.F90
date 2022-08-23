module flux_mod

  use params, only: dp, sb, invert_grid, moist_rad, surface, semi_grey_scheme, A_s
!  use radiation_Kitzmann_noscatt, only: Kitzmann_TS_noscatt
  use condense, only : rain_out
  
#ifdef SOC
  use socrates_interface_mod, only: run_socrates
#elif defined SHORT_CHAR
  use radiation_Kitzmann_noscatt, only: Kitzmann_TS_noscatt
  use toon_mod, only : toon_driver
  use sc_split_mod, only: sc_split
#elif defined PICKET
  use k_Rosseland_mod, only: k_func_Freedman_local, gam_func_Parmentier, AB_func_Parmentier
  use short_char_ross, only: short_char_ross_driver
  use radiation_mod, only : radiation_interface
#elif defined TWOSTR
  use band_grey_mod, only : run_twostr
#endif
  
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
    
#ifdef SOC
    
    real(dp) :: rad_lat, rad_lon, t_surf_in, albedo_in, net_surf_sw_down, surf_lw_down
    real(dp), dimension(size(Tf)) :: q_in, temp_tend, h2_in, ch4_in, he_in
    
#elif defined PICKET
    real(dp), dimension(3) :: gam_V, Beta_V, A_Bond
    real(dp), dimension(2) :: beta
    real(dp) :: gam_1, gam_2, Tint, Tirr
    real(dp) :: met
    real(dp), dimension(2,nf) :: kIR_Ross
    real(dp), dimension(3,nf) :: kV_Ross        

#elif defined TWOSTR
    real(dp), allocatable :: delp(:)
#endif

#ifdef SOC
    
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

    ! Adjust so that bottom upwards flux = sigma T^4 in IR
    !net_F(ne) = net_F(ne) - fup(ne)
    !fup(ne) = sb*Te(ne)**4
    !net_F(ne) = net_F(ne) + fup(ne)
    
#elif defined SHORT_CHAR

    select case(semi_grey_scheme)
       
    case('short_char')
       !call Kitzmann_TS_noscatt(nf, ne, Te, pe, &
       !     net_F, mu_s, Ts, olr, q, fup, fdn, s_dn)
       !write(*,*) 'short char', net_F
       call sc_split(nf, ne, Te, pe, Tf, &
            net_F, mu_s, Ts, olr, q, fup, fdn, s_dn)
       !write(*,*) 'sc_split', net_F
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

       

#elif defined PICKET
    met =0.0_dp
    Tint = (Fint/sb)**0.25_dp
    Tirr = (Finc/sb)**0.25_dp
    call gam_func_Parmentier(Tint, Tirr, 2, 0._dp, 0._dp, gam_V, Beta_V, Beta, gam_1, gam_2, A_Bond)
    do i=1,nf
       call k_func_Freedman_local(Tf(i), pf(i)*10, met, kIR_Ross(1,i))
       kIR_Ross(1,i) = kIR_Ross(1,i)*0.1_dp
       kV_Ross(:,i) = kIR_Ross(1,i)*gam_V

       kIR_Ross(2,i) = kIR_Ross(1,i) * gam_2
       kIR_Ross(1,i) = kIR_Ross(1,i) * gam_1
    end do

    !A_Bond = 0.0_dp
    call short_char_ross_driver(nf,ne,Te,pe,net_F,1.0_dp,Finc, Fint,olr,kV_Ross,kIR_Ross,Beta_V, &
    Beta, A_Bond)!, &
    !kV_Ross, kIR_Ross, Beta_V, Beta, A_Bond)

    !call radiation_interface(pe,pf,Tf,Tf(nf),net_F, &
    !     kV_Ross, kIR_Ross, Beta_V, Beta, A_bond)

#elif defined TWOSTR

    if (invert_grid) then
       allocate(delp(nf+1))
       do i=2,nf
          delp(i) = pf(i) - pf(i-1)
       enddo
       delp(1) = pf(1) - pe(1)
       delp(nf+1) = pe(nf+1) - pf(nf)
    else
       allocate(delp(nf))
       do i=1,nf
          delp(i) = pe(i+1) - pe(i)
       enddo
    endif
    
    call run_twostr(nf, Tf, Te, pf, pe, delp, q, net_F, olr)

    deallocate(delp)
#endif
    
    
  end subroutine get_fluxes
end module flux_mod
