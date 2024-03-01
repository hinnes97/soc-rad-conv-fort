! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! @brief Define the Socrates diagnostic output type

module socrates_def_diag

use realtype_rd, only: RealExt

implicit none

type :: StrDiag

real(RealExt), pointer :: heating_rate(:,:) => null()
! Heating rate, Ks-1 (n_profile, n_layer)

real(RealExt), pointer :: flux_direct(:,:) => null()
! Direct (unscattered) downwards flux, Wm-2 (n_profile, 0:n_layer)

real(RealExt), pointer :: flux_down(:,:) => null()
! Downwards flux, Wm-2 (n_profile, 0:n_layer)

real(RealExt), pointer :: flux_up(:,:) => null()
! Upwards flux, Wm-2 (n_profile, 0:n_layer)

real(RealExt), pointer :: flux_direct_clear(:,:) => null()
! Clear-sky direct (unscattered) downwards flux, Wm-2 (n_profile, 0:n_layer)

real(RealExt), pointer :: flux_down_clear(:,:) => null()
! Clear-sky downwards flux, Wm-2 (n_profile, 0:n_layer)

real(RealExt), pointer :: flux_up_clear(:,:) => null()
! Clear-sky upwards flux, Wm-2 (n_profile, 0:n_layer)

real(RealExt), pointer :: flux_up_tile(:,:) => null()
! Upwards flux on tiles, Wm-2 (n_profile, n_tile)

real(RealExt), pointer :: flux_up_blue_tile(:,:) => null()
! Upwards blue flux on tiles, Wm-2 (n_profile, n_tile)

real(RealExt), pointer :: flux_direct_blue_surf(:) => null()
! Direct blue flux at the surface, Wm-2 (n_profile)

real(RealExt), pointer :: flux_down_blue_surf(:) => null()
! Total downward blue flux at the surface, Wm-2 (n_profile)

real(RealExt), pointer :: total_cloud_cover(:) => null()
! Total cloud cover (n_profile)

real(RealExt), pointer :: total_cloud_fraction(:,:) => null()
! Total cloud fraction in layers (n_profile, n_layer)

real(RealExt), pointer :: liq_frac(:,:) => null()
! Liquid cloud fraction (n_profile, n_layer)

real(RealExt), pointer :: liq_conv_frac(:,:) => null()
! Liquid convective cloud fraction (n_profile, n_layer)

real(RealExt), pointer :: liq_incloud_mmr(:,:) => null()
! Liquid in-cloud mean mixing ratio (n_profile, n_layer)

real(RealExt), pointer :: liq_inconv_mmr(:,:) => null()
! Liquid convective in-cloud mean mixing ratio (n_profile, n_layer)

real(RealExt), pointer :: liq_dim(:,:) => null()
! Cloud droplet effective dimension (n_profile, n_layer)

real(RealExt), pointer :: liq_conv_dim(:,:) => null()
! Convective cloud droplet effective dimension (n_profile, n_layer)

real(RealExt), pointer :: ice_frac(:,:) => null()
! Ice cloud fraction (n_profile, n_layer)

real(RealExt), pointer :: ice_conv_frac(:,:) => null()
! Ice convective cloud fraction (n_profile, n_layer)

real(RealExt), pointer :: ice_incloud_mmr(:,:) => null()
! Ice in-cloud mean mixing ratio (n_profile, n_layer)

real(RealExt), pointer :: ice_inconv_mmr(:,:) => null()
! Ice convective in-cloud mean mixing ratio (n_profile, n_layer)

real(RealExt), pointer :: ice_dim(:,:) => null()
! Cloud ice-crystal effective dimension (n_profile, n_layer)

real(RealExt), pointer :: ice_conv_dim(:,:) => null()
! Convective cloud ice-crystal effective dimension (n_profile, n_layer)

real(RealExt), pointer :: cloud_top_liq_dim(:) => null()
! Cloud droplet effective radius at cloud top weighted by cloud fraction

real(RealExt), pointer :: cloud_top_liq_weight(:) => null()
! Weight for liquid cloud fraction at cloud top

real(RealExt), pointer :: cloud_top_warm_liq_dim(:) => null()
! Warm cloud droplet effective radius at cloud top weighted by cloud fraction

real(RealExt), pointer :: cloud_top_warm_liq_weight(:) => null()
! Weight for warm liquid cloud fraction at cloud top

real(RealExt), pointer :: aerosol_optical_depth(:,:,:) => null()
! Total aerosol optical depth (n_profile, n_layer, n_band)

real(RealExt), pointer :: aerosol_scat_optical_depth(:,:,:) => null()
! Total aerosol scattering optical depth (n_profile, n_layer, n_band)

real(RealExt), pointer :: aerosol_asymmetry_scat(:,:,:) => null()
! Total aerosol asymmetry weighted by scattering optical depth

end type StrDiag

end module socrates_def_diag
