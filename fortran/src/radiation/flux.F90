module flux_mod

  use params, only: dp, sb, invert_grid, moist_rad, surface, semi_grey_scheme, A_s, grav
  use atmosphere, only : nqr, soc_indices, get_mmw, th_gases
  use phys, only : Rstar
#ifdef SOC
  !use socrates_interface_mod, only: run_socrates
  use soc_calc_mod, only: soc_calc
  use socrates_def_diag, only: StrDiag
  use gas_list_pcf
  use soc_init_mod, only : control_lw, control_sw, diag_sw, diag_lw, dimen_sw, dimen_lw
  use rad_pcf, only : ip_solar, ip_infra_red
  use omp_lib, only: omp_get_wtime
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
  subroutine get_fluxes(nf, ne, Tf, pf, Te, pe, delp, &
       net_F, mu_s, Finc, Fint, olr, q, Ts, fup, fdn, s_dn, s_up)
    !! Input variables
    integer, intent(in) :: nf, ne                         ! Number of layers, levels (lev = lay + 1)
    real(dp), dimension(:), intent(in) :: Tf, pf, delp   ! Temperature [K], pressure [pa] at layers
    real(dp), dimension(:), intent(in) :: Te, pe           ! pressure [pa] at levels
    real(dp), intent(in) :: Finc, mu_s                        ! Incident flux [W m-2] and cosine zenith angle
    real(dp), intent(in) :: Fint                              ! Internal flux [W m-2]
    real(dp), intent(inout) :: q(:,:)
    real(dp), intent(in) :: Ts                            ! Surface temperature
    
    !! Output variables
    real(dp), dimension(ne), intent(out) :: net_F, fup, fdn
    real(dp), intent(out) :: olr, s_dn(:), s_up(:)

    integer :: i, n, k, tim
    
#ifdef SOC
    
    real(dp) :: rad_lat, rad_lon, t_surf_in, albedo_in, net_surf_sw_down, surf_lw_down, test
    real(dp), dimension(size(Tf)) :: temp_tend, mass_1d, density_1d
    real(dp), dimension(size(Tf)) :: h2o_1d, h2_1d, he_1d, ch4_1d, co2_1d, co_1d, hcn_1d, &
         nh3_1d, c2h6_1d, n2_1d
    real(dp), dimension(1) :: mu_s_arr, ts_arr, insol_arr
    real(dp), dimension(2), target :: sup_out(1,0:nf), sdn_out(1,0:nf), fup_out(1,0:nf), fdn_out(1,0:nf)
    !type(StrDiag) :: diag_sw

    real(dp) :: start, finish
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

    ! ch4_in = 0.0_dp
    ! q_in = q(:,1)
    ! h2_in = 1-q(:,1)-ch4_in
    
    ! he_in =h2_in*0.25188_dp
    ! h2_in =h2_in*0.74812_dp
    
    albedo_in = A_s

!    q_in = 0.2
!    h2_in = 0.8
!    he_in = 0.0
    !write(*,*) '-----------------------------------------------------------------------------'
    !call run_socrates(rad_lat, rad_lon, Tf, q_in, h2_in, ch4_in, he_in,t_surf_in, pf, pe, pf, pe, albedo_in, &
    !     temp_tend, net_surf_sw_down, surf_lw_down, net_F, fup, fdn, s_up, s_dn)

    !call calc_soc(pf, pe, Tf, Te, q, fup, fdn, s_up, s_dn, soc_indices)
    !olr = fup(1) 
    
    do n = 1,nqr
       select case(soc_indices(n))
       !H2O
       case(IP_h2o)
          h2o_1d = q(:,n)
       case(IP_h2)
          h2_1d = q(:,n)
       case(IP_he)
          he_1d = q(:,n)
       case(IP_hcn)
          hcn_1d = q(:,n)
       case(IP_co)
          co_1d = q(:,n)
       case(IP_co2)
          co2_1d=q(:,n)
       case(IP_nh3)
          nh3_1d = q(:,n)
       case(IP_c2h6)
          c2h6_1d = q(:,n)
       case(IP_n2)
          n2_1d = q(:,n)
       case(IP_ch4)
          ch4_1d = q(:,n)
       end select
    enddo

    do k=1,nf
       call get_mmw(q(k,:), test)
       mass_1d(k) = delp(k)/grav
       density_1d(k) = pf(k) * test/Rstar/Tf(k)
    enddo
    mu_s_arr = mu_s
    ts_arr = Ts
    insol_arr = Finc/mu_s

    ! h2o_1d = 0.2; h2_1d = 0.2; he_1d = 0.2; co_1d = 0.2; co2_1d = 0.2
    ! hcn_1d = 0.0; c2h6_1d =0.0; nh3_1d = 0.0; ch4_1d = 0.0; n2_1d =0.0
    ! SW calculation
    diag_sw%flux_up => sup_out
    diag_sw%flux_down => sdn_out

    
    start = omp_get_wtime()

!    do tim=1,10
    q = 1.0
    call soc_calc(n_profile            = 1, &
                  n_layer              = nf, &
                  diag                 = diag_sw, &
                  control              = control_sw, &
                  dimen                = dimen_sw, &
                  p_layer_1d           = pf, &
                  t_layer_1d           = Tf, &
                  t_level_1d           = Te, &
                  t_ground             = ts_arr, &
                  mass_1d              = mass_1d, &
                  density_1d           = density_1d, &
                  h2o_1d               = h2o_1d, &
                  h2_1d                = h2_1d, &
                  he_1d                = he_1d, &
                  ch4_1d               = ch4_1d, &
                  co_1d                = co_1d, &
                  co2_1d               = co2_1d,&
                  hcn_1d               = hcn_1d,&
                  n2_1d                = n2_1d, &
                  c2h6_1d              = c2h6_1d, &
                  nh3_1d               = nh3_1d, &
                  cos_zenith_angle      = mu_s_arr, &
                  i_source             = ip_solar, &!, &
                  solar_irrad          = insol_arr, &
                  l_grey_albedo        = .true., &
                  grey_albedo          = albedo_in)
    s_up(1:ne) = sup_out(1,0:nf)
    s_dn(1:ne) = sdn_out(1,0:nf)

    ! open(unit=10, file='omp_tests/radout_sw.txt')
    ! do k=1,ne
    !    write(10,*) s_up(k), s_dn(k)
    ! enddo
    ! close(10)

    ! LW calculation

    diag_lw%flux_up => fup_out
    diag_lw%flux_down => fdn_out
    call soc_calc(n_profile            = 1, &
                   n_layer              = nf, &
                   diag                 = diag_lw, &
                   control              = control_lw, &
                   dimen             = dimen_lw, &
                   p_layer_1d           = pf, &
                   t_layer_1d           = Tf, &
                   t_level_1d           = Te, &
                   t_ground             = ts_arr, &
                   mass_1d              = mass_1d, &
                   density_1d           = density_1d, &
                   h2o_1d               = h2o_1d, &
                   h2_1d                = h2_1d, &
                   he_1d                = he_1d, &
                   ch4_1d               = ch4_1d, &
                   co_1d                = co_1d, &
                   co2_1d               = co2_1d,&
                   hcn_1d               = hcn_1d,&
                   n2_1d                = n2_1d, &
                   c2h6_1d              = c2h6_1d, &
                   nh3_1d               = nh3_1d, &
                   cos_zenith_angle     = mu_s_arr, &
                   i_source             = ip_infra_red)
    
    fup(1:ne) = fup_out(1,0:nf)
    fdn(1:ne) = fdn_out(1,0:nf)

    ! open(unit=10, file='omp_tests/radout_lw.txt')
    ! do k=1,ne
    !    write(10,*) fup(k), fdn(k)
    ! enddo
    ! close(10)
    ! enddo
    ! finish=omp_get_wtime()
    !  open(unit=10, file='omp_tests/timing.txt')
    !  write(10,*) finish-start
    !  close(10)
    !  stop
    net_F = fup + s_up - fdn - s_dn

    
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
