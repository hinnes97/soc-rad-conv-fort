! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Set the diagnostics to be output from the Socrates runes interface
!
!------------------------------------------------------------------------------
module socrates_set_diag
implicit none
character(len=*), parameter, private :: ModuleName = 'SOCRATES_SET_DIAG'
contains

subroutine set_diag(diag, control, dimen, spectrum, &
  atm, cld, mcica_data, aer, bound, radout, &
  n_profile, n_layer, profile_list, n_layer_stride, n_tile, &
  layer_heat_capacity, layer_heat_capacity_1d, &
  l_invert, l_profile_last)

use socrates_def_diag, only: StrDiag

use def_control,  only: StrCtrl
use def_dimen,    only: StrDim
use def_spectrum, only: StrSpecData
use def_atm,      only: StrAtm
use def_cld,      only: StrCld
use def_mcica,    only: StrMcica
use def_aer,      only: StrAer
use def_bound,    only: StrBound
use def_out,      only: StrOut

use realtype_rd, only: RealK, RealExt
use rad_pcf, only: ip_cloud_off, ip_mcica, &
                   ip_clcmp_st_water, ip_clcmp_st_ice, &
                   ip_clcmp_cnv_water, ip_clcmp_cnv_ice

use socrates_cloud_level_diag, only: cloud_level_diag

implicit none

! Output diagnostic fields
type(StrDiag),   intent(inout) :: diag

! Control options:
type(StrCtrl),      intent(in) :: control

! Dimensions:
type(StrDim),       intent(in) :: dimen

! Spectral data:
type (StrSpecData), intent(in) :: spectrum

! Atmospheric properties:
type(StrAtm),       intent(in) :: atm

! Cloud properties:
type(StrCld),       intent(in) :: cld

! MCICA data:
type(StrMcica),     intent(in) :: mcica_data

! Aerosol properties:
type(StrAer),       intent(in) :: aer

! Boundary conditions:
type(StrBound),     intent(in) :: bound

! Output fields from core radiation code:
type(StrOut),       intent(in) :: radout

integer, intent(in) :: n_profile
!   Number of columns to operate on
integer, intent(in) :: n_layer
!   Number of layers for radiation
integer, intent(in), optional :: profile_list(:)
!   List of profiles to fill in output fields
integer, intent(in), optional :: n_layer_stride
!   Number of layers in input 1d arrays
integer, intent(in), optional :: n_tile
!   Number of surface tiles

real(RealExt), intent(in), optional :: layer_heat_capacity(:, :)
real(RealExt), intent(in), optional :: layer_heat_capacity_1d(:)
!   Heat capacity of layer

logical, intent(in), optional :: l_invert
!   Flag to invert fields in the vertical
logical, intent(in), optional :: l_profile_last
!   Loop over profiles is last in input fields and diagnostics


integer :: i, ii, l, ll, k, list(n_profile)
!   Loop variables
integer :: layer_offset, level_offset, stride_layer
!   Offset to loop counters to allow indexing in inverted order
logical :: l_last
!   Local flag to loop over profiles last in 1d fields and diagnostics
real(RealExt) :: flux_divergence(n_profile, n_layer)
!   Flux divergence across layer (Wm-2)


if (present(profile_list)) then
  list = profile_list(1:n_profile)
else
  do l=1, n_profile
    list(l) = l
  end do
end if

layer_offset = 0
level_offset = 0
if (present(l_invert)) then
  if (l_invert) then
    ! The layer is indexed using an inverted loop counter
    layer_offset = n_layer + 1
    level_offset = n_layer
  end if
end if

if (present(l_profile_last)) then
  l_last = l_profile_last
else
  l_last = .false.
end if

! Set the number of layers in the 1d arrays
if (present(n_layer_stride)) then
  stride_layer = n_layer_stride
else
  stride_layer = n_layer
end if


!------------------------------------------------------------------------------
! Heating rates
!------------------------------------------------------------------------------
if (associated(diag%heating_rate)) then
  if (control%l_flux_div) then
    do i=1, n_layer
      ii = abs(layer_offset-i)
      do l=1, n_profile
        flux_divergence(l, ii) = real( &
          sum(radout%flux_div(l, i, 1:control%n_channel)), RealExt)
      end do
    end do
  else
    do i=1, n_layer
      ii = abs(layer_offset-i)
      do l=1, n_profile
        flux_divergence(l, ii) = real( &
          sum(radout%flux_down(l, i-1, 1:control%n_channel)) - &
          sum(radout%flux_down(l, i,   1:control%n_channel)) + &
          sum(radout%flux_up(  l, i,   1:control%n_channel)) - &
          sum(radout%flux_up(  l, i-1, 1:control%n_channel)), RealExt)
      end do
    end do
  end if
  if (present(layer_heat_capacity)) then
    if (l_last) then
      do i=1, n_layer
        do l=1, n_profile
          diag%heating_rate(i, list(l)) = &
            flux_divergence(l, i) / layer_heat_capacity(i, list(l))
        end do
      end do
    else
      do i=1, n_layer
        do l=1, n_profile
          diag%heating_rate(list(l), i) = &
            flux_divergence(l, i) / layer_heat_capacity(list(l), i)
        end do
      end do
    end if
  else if (present(layer_heat_capacity_1d)) then
    if (l_last) then
      do l=1, n_profile
        do i=1, n_layer
          ll = stride_layer*(list(l)-1) + i
          diag%heating_rate(i, list(l)) = &
            flux_divergence(l, i) / layer_heat_capacity_1d(ll)
        end do
      end do
    else
      do i=1, n_layer
        do l=1, n_profile
          ll = n_profile*(i-1) + list(l)
          diag%heating_rate(list(l), i) = &
            flux_divergence(l, i) / layer_heat_capacity_1d(ll)
        end do
      end do
    end if
  else
    ! Just return the flux_divergence if no heat capacity supplied
    if (l_last) then
      do i=1, n_layer
        do l=1, n_profile
          diag%heating_rate(i, list(l)) = flux_divergence(l, i)
        end do
      end do
    else
      do i=1, n_layer
        do l=1, n_profile
          diag%heating_rate(list(l), i) = flux_divergence(l, i)
        end do
      end do
    end if
  end if
end if


!------------------------------------------------------------------------------
! Fluxes
!------------------------------------------------------------------------------
call sum_flux_channels(diag%flux_direct, radout%flux_direct)
call sum_flux_channels(diag%flux_down, radout%flux_down)
call sum_flux_channels(diag%flux_up, radout%flux_up)
call sum_flux_channels(diag%flux_direct_clear, radout%flux_direct_clear)
call sum_flux_channels(diag%flux_down_clear, radout%flux_down_clear)
call sum_flux_channels(diag%flux_up_clear, radout%flux_up_clear)

call sum_tile_channels(diag%flux_up_tile, radout%flux_up_tile)
call sum_tile_channels(diag%flux_up_blue_tile, radout%flux_up_blue_tile)
if (associated(diag%flux_direct_blue_surf)) then
  do l=1, n_profile
    diag%flux_direct_blue_surf(list(l)) &
      = real(radout%flux_direct_blue_surf(l), RealExt)
  end do
end if
if (associated(diag%flux_down_blue_surf)) then
  do l=1, n_profile
    diag%flux_down_blue_surf(list(l)) &
      = real(radout%flux_down_blue_surf(l), RealExt)
  end do
end if

!------------------------------------------------------------------------------
! Cloud diagnostics
!------------------------------------------------------------------------------
if (associated(diag%total_cloud_cover)) then
  if (control%i_cloud_representation == ip_cloud_off) then
    do l=1, n_profile
      diag%total_cloud_cover(list(l)) = 0.0_RealExt
    end do
  else
    if (control%i_inhom == ip_mcica) then
      do l=1, n_profile
        diag%total_cloud_cover(list(l)) &
          = real(cld%n_subcol_cld(l), RealExt) &
          / real(mcica_data%n_subcol_gen, RealExt)
      end do
    else
      do l=1, n_profile
        diag%total_cloud_cover(list(l)) &
          = real(radout%tot_cloud_cover(l), RealExt)
      end do
    end if
  end if
end if

if (associated(diag%total_cloud_fraction)) then
  if (l_last) then
    do l=1, n_profile
      do i=1, n_layer
        diag%total_cloud_fraction(i, list(l)) = 0.0_RealExt
      end do
    end do
  else
    do i=1, n_layer
      do l=1, n_profile
        diag%total_cloud_fraction(list(l), i) = 0.0_RealExt
      end do
    end do
  end if
  if (control%i_cloud_representation /= ip_cloud_off) then
    call set_cloud_prop(diag%total_cloud_fraction, cld%w_cloud)
  end if
end if

! Cloud component diagnostics
if (associated(diag%liq_dim))         diag%liq_dim         = 0.0_RealExt
if (associated(diag%liq_incloud_mmr)) diag%liq_incloud_mmr = 0.0_RealExt
if (associated(diag%liq_frac))        diag%liq_frac        = 0.0_RealExt
if (associated(diag%ice_dim))         diag%ice_dim         = 0.0_RealExt
if (associated(diag%ice_incloud_mmr)) diag%ice_incloud_mmr = 0.0_RealExt
if (associated(diag%ice_frac))        diag%ice_frac        = 0.0_RealExt
if (associated(diag%liq_conv_dim))    diag%liq_conv_dim    = 0.0_RealExt
if (associated(diag%liq_inconv_mmr))  diag%liq_inconv_mmr  = 0.0_RealExt
if (associated(diag%liq_conv_frac))   diag%liq_conv_frac   = 0.0_RealExt
if (associated(diag%ice_conv_dim))    diag%ice_conv_dim    = 0.0_RealExt
if (associated(diag%ice_inconv_mmr))  diag%ice_inconv_mmr  = 0.0_RealExt
if (associated(diag%ice_conv_frac))   diag%ice_conv_frac   = 0.0_RealExt
if (control%i_cloud_representation /= ip_cloud_off) then
  do k=1, cld%n_condensed
    select case (cld%type_condensed(k))
    case (ip_clcmp_st_water)
      call set_cloud_prop(diag%liq_dim, cld%condensed_dim_char(:,:,k))
      call set_cloud_prop(diag%liq_incloud_mmr, cld%condensed_mix_ratio(:,:,k))
      call set_cloud_prop(diag%liq_frac, cld%w_cloud, &
                          cld%frac_cloud(:, :, cld%i_cloud_type(k)))
    case (ip_clcmp_st_ice)
      call set_cloud_prop(diag%ice_dim, cld%condensed_dim_char(:,:,k))
      call set_cloud_prop(diag%ice_incloud_mmr, cld%condensed_mix_ratio(:,:,k))
      call set_cloud_prop(diag%ice_frac, cld%w_cloud, &
                          cld%frac_cloud(:, :, cld%i_cloud_type(k)))
    case (ip_clcmp_cnv_water)
      call set_cloud_prop(diag%liq_conv_dim, cld%condensed_dim_char(:,:,k))
      call set_cloud_prop(diag%liq_inconv_mmr, cld%condensed_mix_ratio(:,:,k))
      call set_cloud_prop(diag%liq_conv_frac, cld%w_cloud, &
                          cld%frac_cloud(:, :, cld%i_cloud_type(k)))
    case (ip_clcmp_cnv_ice)
      call set_cloud_prop(diag%ice_conv_dim, cld%condensed_dim_char(:,:,k))
      call set_cloud_prop(diag%ice_inconv_mmr, cld%condensed_mix_ratio(:,:,k))
      call set_cloud_prop(diag%ice_conv_frac, cld%w_cloud, &
                          cld%frac_cloud(:, :, cld%i_cloud_type(k)))
    end select
  end do
end if

! Observed effective radius diagnostics
if (associated(diag%cloud_top_liq_dim) .and. &
    associated(diag%cloud_top_liq_weight)) then
  call cloud_level_diag(control, dimen, atm, cld, .true., &
    list, diag%cloud_top_liq_dim, diag%cloud_top_liq_weight)
end if
if (associated(diag%cloud_top_warm_liq_dim) .and. &
    associated(diag%cloud_top_warm_liq_weight)) then
  call cloud_level_diag(control, dimen, atm, cld, .false., &
    list, diag%cloud_top_warm_liq_dim, diag%cloud_top_warm_liq_weight)
end if

!------------------------------------------------------------------------------
! Aerosol diagnostics
!------------------------------------------------------------------------------
call set_aerosol_prop(diag%aerosol_optical_depth, &
                      radout%aerosol_absorption_band, &
                      radout%aerosol_scattering_band)
call set_aerosol_prop(diag%aerosol_scat_optical_depth, &
                      radout%aerosol_scattering_band)
call set_aerosol_prop(diag%aerosol_asymmetry_scat, &
                      radout%aerosol_asymmetry_band)

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------
subroutine sum_flux_channels(field, field_channels)

  implicit none
  
  real(RealExt), intent(inout), pointer :: field(:, :)
  real(RealK), intent(in), allocatable :: field_channels(:, :, :)
  
  if (associated(field)) then
    if (l_last) then
      do i=0, n_layer
        ii = abs(level_offset-i)
        do l=1, n_profile
          field(ii, list(l)) &
            = real(sum(field_channels(l, i, 1:control%n_channel)), RealExt)
        end do
      end do
    else
      do i=0, n_layer
        ii = abs(level_offset-i)
        do l=1, n_profile
          field(list(l), ii) &
            = real(sum(field_channels(l, i, 1:control%n_channel)), RealExt)
        end do
      end do
    end if
  end if

end subroutine sum_flux_channels


!------------------------------------------------------------------------------
subroutine sum_tile_channels(field, field_channels)
  
  implicit none
  
  real(RealExt), intent(inout), pointer :: field(:, :)
  real(RealK), intent(in), allocatable :: field_channels(:, :, :)
  
  integer :: ll
  
  if (associated(field).and.present(n_tile)) then
    field = 0.0_RealExt
    if (control%l_tile) then
      if (l_last) then
        do i=1, n_tile
          do ll=1, bound%n_point_tile
            l = bound%list_tile(ll)
            field(i, list(l)) &
              = real(sum(field_channels(ll, i, 1:control%n_channel)), RealExt)
          end do
        end do
      else
        do i=1, n_tile
          do ll=1, bound%n_point_tile
            l = bound%list_tile(ll)
            field(list(l), i) &
              = real(sum(field_channels(ll, i, 1:control%n_channel)), RealExt)
          end do
        end do
      end if
    end if
  end if
  
end subroutine sum_tile_channels


!------------------------------------------------------------------------------
subroutine set_cloud_prop(field, cld_field1, cld_field2)

  implicit none

  real(RealExt), intent(inout), pointer :: field(:, :)
  real(RealK), intent(in) :: cld_field1(:, dimen%id_cloud_top:)
  real(RealK), intent(in), optional :: cld_field2(:, dimen%id_cloud_top:)
  integer :: i_lower, i_upper, i_dim

  if (associated(field)) then
    if (l_last) then
      i_dim = 1
    else
      i_dim = 2
    end if
    if (layer_offset == n_layer + 1) then
      ! Field output is inverted
      i_lower = max(n_layer + 1 - ubound(field, i_dim), dimen%id_cloud_top)
      i_upper = min(n_layer + 1 - lbound(field, i_dim), n_layer)
    else
      i_lower = max(lbound(field, i_dim), dimen%id_cloud_top)
      i_upper = min(ubound(field, i_dim), n_layer)
    end if
    ! Fill diagnostic between requested layers
    if (present(cld_field2)) then
      if (l_last) then
        do i=i_lower, i_upper
          ii = abs(layer_offset-i)
          do l=1, n_profile
            field(ii, list(l)) = real( &
              cld_field1(l, i) * cld_field2(l, i), RealExt)
          end do
        end do
      else
        do i=i_lower, i_upper
          ii = abs(layer_offset-i)
          do l=1, n_profile
            field(list(l), ii) = real( &
              cld_field1(l, i) * cld_field2(l, i), RealExt)
          end do
        end do
      end if
    else
      if (l_last) then
        do i=i_lower, i_upper
          ii = abs(layer_offset-i)
          do l=1, n_profile
            field(ii, list(l)) = real(cld_field1(l, i), RealExt)
          end do
        end do
      else
        do i=i_lower, i_upper
          ii = abs(layer_offset-i)
          do l=1, n_profile
            field(list(l), ii) = real(cld_field1(l, i), RealExt)
          end do
        end do
      end if
    end if
  end if

end subroutine set_cloud_prop


!------------------------------------------------------------------------------
subroutine set_aerosol_prop(field, aer_field1, aer_field2)

  implicit none

  real(RealExt), intent(inout), pointer :: field(:, :, :)
  real(RealK), intent(in), allocatable :: aer_field1(:, :, :)
  real(RealK), intent(in), allocatable, optional :: aer_field2(:, :, :)
  integer :: k_lower, k_upper, k_dim
  integer :: i_lower, i_upper, i_dim

  if (associated(field)) then
    if (l_last) then
      i_dim = 1
      k_dim = 2
    else
      i_dim = 2
      k_dim = 3
    end if
    k_lower = max(lbound(field, k_dim), 1)
    k_upper = min(ubound(field, k_dim), spectrum%basic%n_band)
    i_lower = max(lbound(field, i_dim), 1)
    i_upper = min(ubound(field, i_dim), n_layer)

    if (present(aer_field2)) then
      if (l_last) then
        do k=k_lower, k_upper
          do i=i_lower, i_upper
            ii = abs(layer_offset-i)
            do l=1, n_profile
              field(i, k, list(l)) = real(( aer_field1(l, ii, k) &
                + aer_field2(l, ii, k) ) * atm%mass(l, ii), RealExt)
            end do
          end do
        end do
      else
        do k=k_lower, k_upper
          do i=i_lower, i_upper
            ii = abs(layer_offset-i)
            do l=1, n_profile
              field(list(l), i, k) = real(( aer_field1(l, ii, k) &
                + aer_field2(l, ii, k) ) * atm%mass(l, ii), RealExt)
            end do
          end do
        end do
      end if
    else
      if (l_last) then
        do k=k_lower, k_upper
          do i=i_lower, i_upper
            ii = abs(layer_offset-i)
            do l=1, n_profile
              field(i, k, list(l)) = real( &
                aer_field1(l, ii, k) * atm%mass(l, ii), RealExt)
            end do
          end do
        end do
      else
        do k=k_lower, k_upper
          do i=i_lower, i_upper
            ii = abs(layer_offset-i)
            do l=1, n_profile
              field(list(l), i, k) = real( &
                aer_field1(l, ii, k) * atm%mass(l, ii), RealExt)
            end do
          end do
        end do
      end if
    end if
  end if

end subroutine set_aerosol_prop

end subroutine set_diag
end module socrates_set_diag
