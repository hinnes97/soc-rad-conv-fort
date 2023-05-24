! Dummy module since we are not putting socrates under version control
module socrates_interface_mod
  ! 
  use netcdf
  use iso_fortran_env, only: real64
  implicit none

  integer, parameter ::  dp = real64
contains
  subroutine run_socrates(rad_lat, rad_lon, temp_in, q_in, h2_in, ch4_in, he_in, t_surf_in, p_full_in,&
    p_half_in, z_full_in, z_half_in, albedo_in, &
       temp_tend, net_surf_sw_down, surf_lw_down, net_flux, lw_up, lw_dn,sw_up, sw_dn)  
    real(dp), intent(in)                 :: t_surf_in, albedo_in
    real(dp), intent(in), dimension(:)   :: temp_in, p_full_in, q_in, z_full_in, h2_in, ch4_in,he_in
    real(dp), intent(in), dimension(:)  :: p_half_in, z_half_in
    real(dp), intent(inout), dimension(:) :: temp_tend
    real(dp), intent(out)   :: net_surf_sw_down, surf_lw_down
    real(dp), intent(in) :: rad_lat, rad_lon
    real(dp), intent(out) :: net_flux(:)
    real(dp), intent(out) :: lw_up(:), lw_dn(:), sw_up(:), sw_dn(:)
  end subroutine run_socrates
end module socrates_interface_mod
