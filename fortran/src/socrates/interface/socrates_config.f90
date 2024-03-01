module socrates_config_mod

use def_control,  only: StrCtrl
use def_spectrum, only: StrSpecData
USE filenamelength_mod, ONLY: filenamelength
USE missing_data_mod, ONLY: rmdi, imdi
USE realtype_rd, ONLY: RealK
use rad_pcf, only : ip_overlap_k_eqv_scl

implicit none


TYPE (StrCtrl), SAVE :: sw_control, lw_control
TYPE (StrSpecData), SAVE :: sw_spectrum, lw_spectrum


! Spectral region and bands
  INTEGER :: isolir                                               = imdi
!   Spectral region
  INTEGER :: first_band                                           = imdi
!   First band to use in the calculation
  INTEGER :: last_band                                            = imdi
!   Last band to use in the calculation


! Physical processes
  LOGICAL :: l_microphysics                                       = .FALSE.
!   Flag for microphysics
  LOGICAL :: l_gas                                                = .FALSE.
!   Flag for gaseous absorption
  LOGICAL :: l_rayleigh                                          = .FALSE.
!   Flag for Rayleigh scattering
  LOGICAL :: l_continuum                                          = .FALSE.
!   Flag for the continuum
  LOGICAL :: l_cont_gen                                           = .FALSE.
!   Flag for the generalised continua
  LOGICAL :: l_cloud                                              = .FALSE.
!   Flag for clouds
  LOGICAL :: l_drop                                               = .FALSE.
!   Flag for droplets
  LOGICAL :: l_ice                                                = .FALSE.
!   Flag for ice crystals
  LOGICAL :: l_aerosol                                            = .FALSE.
!   Flag for aerosols
  LOGICAL :: l_aerosol_mode                                       = .FALSE.
!   Flag for modal aerosols
  LOGICAL :: l_aerosol_ccn                                        = .FALSE.
!   Flag for aerosols as CCN
  LOGICAL :: l_solar_tail_flux                                    = .FALSE.
!   Flag for adding solar tail flux to LW ragion
  LOGICAL :: l_orog                                               = .FALSE.
!   Correct the direct solar flux at the surface for sloping terrain
  LOGICAL :: l_orog_fix                                           = .FALSE.
!   Fix a bug in the orographic correction for the column solver
  LOGICAL :: l_solvar                                             = .FALSE.
!   Time variation of the solar spectrum

! Gaseous absorption:
  INTEGER :: i_gas_overlap                                        = ip_overlap_k_eqv_scl
!   Treatment of gaseous overlaps
  INTEGER :: n_esft_red                                           = imdi
!   Number of reduced (resorted and rebinned) ESFT terms in
!   each band using random overlap with resorting and rebinning
  REAL (RealK) :: gpnt_split                                      = rmdi
!   g-coordinate point for splitting into two subintervals using
!   random overlap with resorting and rebinning
  INTEGER :: i_gas                                                = imdi
!   Gas to be considered (if only one gas)
  LOGICAL :: l_h2o                                                = .FALSE.
!   Flag for absorption by water vapour
  LOGICAL :: l_o2                                                 = .FALSE.
!   Flag for absorption by oxygen
  LOGICAL :: l_n2o                                                = .FALSE.
!   Flag for absorption by nitrous oxide
  LOGICAL :: l_ch4                                                = .FALSE.
!   Flag for absorption by methane
  LOGICAL :: l_cfc11                                              = .FALSE.
!   Flag for absorption by CFC11
  LOGICAL :: l_cfc12                                              = .FALSE.
!   Flag for absorption by CFC12
  LOGICAL :: l_cfc113                                             = .FALSE.
!   Flag for absorption by CFC113
  LOGICAL :: l_cfc114                                             = .FALSE.
!   Flag for absorption by CFC114
  LOGICAL :: l_hcfc22                                             = .FALSE.
!   Flag for absorption by HCFC22
  LOGICAL :: l_hfc125                                             = .FALSE.
!   Flag for absorption by HFC125
  LOGICAL :: l_hfc134a                                            = .FALSE.
!   Flag for absorption by HFC134A
  LOGICAL :: l_co                                                 = .FALSE.
!   Flag for absorption by carbon monoxide
  LOGICAL :: l_nh3                                                = .FALSE.
!   Flag for absorption by ammonia
  LOGICAL :: l_tio                                                = .FALSE.
!   Flag for absorption by titanium oxide
  LOGICAL :: l_vo                                                 = .FALSE.
!   Flag for absorption by vanadium oxide
  LOGICAL :: l_h2                                                 = .FALSE.
!   Flag for absorption by hydrogen
  LOGICAL :: l_he                                                 = .FALSE.
!   Flag for absorption by helium
  LOGICAL :: l_na                                                 = .FALSE.
!   Flag for absorption by sodium
  LOGICAL :: l_k                                                  = .FALSE.
!   Flag for absorption by potassium
  LOGICAL :: l_co2                                                = .FALSE.
!   Flag for absorption by carbon dioxide
  LOGICAL :: l_feh                                                = .FALSE.
!   Flag for absorption by iron hydride
  LOGICAL :: l_crh                                                = .FALSE.
!   Flag for absorption by chromium hydride
  LOGICAL :: l_li                                                 = .FALSE.
!   Flag for absorption by lithium
  LOGICAL :: l_rb                                                 = .FALSE.
!   Flag for absorption by rubidium
  LOGICAL :: l_cs                                                 = .FALSE.
!   Flag for absorption by cesium
  LOGICAL :: l_ph3                                                = .FALSE.
!   Flag for absorption by phosphine
  LOGICAL :: l_c2h2                                               = .FALSE.
!   Flag for absorption by acetylene
  LOGICAL :: l_hcn                                                = .FALSE.
!   Flag for absorption by hydrogen cyanide
  LOGICAL :: l_o3                                                 = .FALSE.
!   Flag for absorption by ozone
  LOGICAL :: l_so2                                                = .FALSE.
!   Flag for absorption by sulfur dioxide
  LOGICAL :: l_no2                                                = .FALSE.
!   Flag for absorption by nitrogen dioxide
  LOGICAL :: l_n2                                                 = .FALSE.
!   Flag for absorption by nitrogen
  LOGICAL :: l_ar                                                 = .FALSE.
!   Flag for absorption by argon
  LOGICAL :: l_o                                                  = .FALSE.
!   Flag for absorption by atomic oxygen
  LOGICAL :: l_n                                                  = .FALSE.
!   Flag for absorption by atomic nitrogen
  LOGICAL :: l_hno3                                               = .FALSE.
!   Flag for absorption by nitric acid
  LOGICAL :: l_no                                                 = .FALSE.
!   Flag for absorption by nitric oxide
  LOGICAL :: l_no3                                                = .FALSE.
!   Flag for absorption by nitrate radical
  LOGICAL :: l_n2o5                                               = .FALSE.
!   Flag for absorption by nitrogen pentoxide
  LOGICAL :: l_hono                                               = .FALSE.
!   Flag for absorption by nitrous acid
  LOGICAL :: l_ho2no2                                             = .FALSE.
!   Flag for absorption by peroxynitric acid
  LOGICAL :: l_h2o2                                               = .FALSE.
!   Flag for absorption by hydrogen peroxide
  LOGICAL :: l_c2h6                                               = .FALSE.
!   Flag for absorption by ethane
  LOGICAL :: l_ch3                                                = .FALSE.
!   Flag for absorption by methyl radical
  LOGICAL :: l_h2co                                               = .FALSE.
!   Flag for absorption by formaldehyde
  LOGICAL :: l_ho2                                                = .FALSE.
!   Flag for absorption by hydroperoxy radical
  LOGICAL :: l_all_gases                                          = .FALSE.
!   Flag for absorption by all gases in spectral file

! Properties of clouds:
  INTEGER :: i_cloud                                              = imdi
!   Cloud scheme
  INTEGER :: i_cloud_representation                               = imdi
!   Representation of clouds
  INTEGER :: i_st_water                                           = imdi
!   Type of water droplet in stratiform clouds
  INTEGER :: i_cnv_water                                          = imdi
  !   Type of water droplet in convective clouds
  INTEGER :: i_st_ice                                             = imdi
!   Type of ice crystal in stratiform clouds
  INTEGER :: i_cnv_ice                                            = imdi
!   Type of ice crystal in convective clouds
  INTEGER :: i_inhom                                              = imdi
!   Method of treating cloud water content variability
  INTEGER :: i_mcica_sampling                                     = imdi
!   Method of sampling sub-grid cloud using MCICA
  INTEGER :: i_overlap                                            = imdi
!   Method of treating cloud vertical overlap
  INTEGER :: i_drop_re                                            = imdi
!   Method of setting the effective radii of cloud droplets
  LOGICAL :: l_local_cnv_partition                                = .FALSE.
!   Flag to partition convective clouds between water and
!   ice using the local temperature
  LOGICAL :: l_global_cloud_top                                   = .FALSE.
!   Flag to use a global value for the topmost cloudy layer
!   (This is used to obtained bit-reproducible results
!   across different configurations of PEs on MPP systems)
  LOGICAL :: l_avg_phase_fnc                                      = .FALSE.
!   Use a grid-box average cloud phase function for sub-grid cloud

! Angular integration (including algorithmic options):
  INTEGER :: n_channel                                            = imdi
!   Number of channels in output
  INTEGER :: i_angular_integration                                = imdi
!   Method of angular integration
  INTEGER :: i_2stream                                            = imdi
!   Two-stream scheme
  INTEGER :: i_solver                                             = imdi
!   Two-stream solver
  INTEGER :: i_solver_clear                                       = imdi
!   Clear-sky solver
  INTEGER :: n_order_gauss                                        = imdi
!   Order of Gaussian quadrature
  INTEGER :: i_truncation                                         = imdi
!   Type of truncation for spherical harmonics
  INTEGER :: i_sph_algorithm                                      = imdi
!   Algorithm used for spherical harmonic calculations
  INTEGER :: n_order_phase_solar                                  = imdi
!   Order of truncation of the solar phase function
  INTEGER :: ls_global_trunc                                      = imdi
!   Global order of truncation
  INTEGER :: ms_min                                               = imdi
!   Minimum azimuthal order
  INTEGER :: ms_max                                               = imdi
!   Maximum azimuthal order
  INTEGER :: ls_brdf_trunc                                        = imdi
!   Order of truncation of BRDFs
  INTEGER :: n_order_forward                                      = imdi
!   Order of the term used to `define' the forward scattering fraction.
  INTEGER :: i_sph_mode                                           = imdi
!   Mode of operation of spherical harmonic code
  INTEGER :: i_scatter_method                                     = imdi
!   Method of treating scattering
  INTEGER :: i_solar_src                                          = imdi
!   Index of solar source function
!   i_solar_src = 1  Original Kurucz function
!   i_solar_src = 2  AER version of Kurucz
!   i_solar_src = 3  Kurucz function used in UK model
!   i_solar_src = 4  UK reduced version
!   i_solar_src = 5  Kurucz modtran
!   i_solar_src = 6  Labs Neckel
!                         Z. Sun
  INTEGER :: i_direct_tau                                         = imdi
!   Method of scaling optical depth for direct flux calculation
  LOGICAL :: l_ir_source_quad                                     = .FALSE.
!   Flag to use a quadratic source function in the IR
  LOGICAL :: l_rescale                                            = .FALSE.
!   Flag for rescaling
  LOGICAL :: l_spherical_solar                                    = .FALSE.
!   Flag to use spherical geometry for the direct beam path
  LOGICAL :: l_spherical_diffuse                                  = .FALSE.
!   Flag to use spherical geometry for the diffuse source terms
  LOGICAL :: l_henyey_greenstein_pf                               = .FALSE.
!   Flag to use Henyey-Greenstein phase functions
  LOGICAL :: l_lanczos                                            = .FALSE.
!   Flag to use Lanczos smoothing of solar phf
  LOGICAL :: l_euler_trnf                                         = .FALSE.
!   Flag to apply Euler's transformation to alternating series
  REAL (RealK) :: accuracy_adaptive                               = rmdi
!   Accuracy for adaptive truncation
  REAL (RealK) :: euler_factor                                    = rmdi
!   Factor applied to the last term of an alternating series
  REAL (RealK) :: half_angle                                      = rmdi
!   Half angle of FOV used for CSR fraction calculation 

! Miscallaneous options
  LOGICAL :: l_tile                                               = .FALSE.
!   Allow tiling of the surface
  LOGICAL :: l_tile_emissivity                                    = .FALSE.
!   Use tile emissivities to calculate the ground source function
  LOGICAL :: l_flux_ground                                        = .FALSE.
!   Use emission from the surface
  LOGICAL, ALLOCATABLE :: l_flux_tile(:)
!   Use emission from selected tiles
  LOGICAL :: l_extra_top                                          = .FALSE.
!   Flag to insert an extra layer into radiation above the
!   top of the model (this is sometimes desirable to ensure
!   proper radiative heating in the top resolved layer).
  LOGICAL :: l_rad_deg                                            = .FALSE.
!   Flag to apply spatial degradation in radiation: in this case
!   radiation quantities are interpolated to all grid-points,
!   whereas subsampling refers to selecting a portion of the area
!   on the PE and returning radiative quentatities only at those
!   points
  LOGICAL :: l_subsample                                          = .FALSE.
!   Flag to apply spatial subsampling (for satellite footprints)
!   in radiation.
  LOGICAL :: l_mixing_ratio                                       = .FALSE.
!   True if mixing ratios are with respect to dry mass
  LOGICAL :: l_map_sub_bands                                      = .FALSE.
!   Map sub-bands to channels rather than bands


! Band-by-band control options
  INTEGER, ALLOCATABLE :: i_scatter_method_band(:)
!   Method of treating scattering in each band
  INTEGER, ALLOCATABLE :: i_gas_overlap_band(:)
!   Gas overlap assumption in each band
  INTEGER, ALLOCATABLE :: map_channel(:)
!   Mapping of actual bands or sub-bands to the output channels
  REAL (RealK), ALLOCATABLE :: weight_band(:)
!   Weighting function for bands
  REAL (RealK), ALLOCATABLE :: weight_diag(:)
!   Weighting function for bands for diagnostic fluxes
  LOGICAL, ALLOCATABLE :: l_clear_band(:)
!   Calculate clear-sky fluxes for band


! Switches for diagnostic output
  LOGICAL :: l_clear                                              = .FALSE.
!   Calculate clear-sky fluxes
  LOGICAL :: l_flux_div                                           = .FALSE.
!   Calculate flux divergence
  LOGICAL :: l_blue_flux_surf                                     = .FALSE.
!   Calculate blue surface fluxes
  LOGICAL :: l_cloud_absorptivity                                 = .FALSE.
!   Calculate absorptivity of clouds (only infra-red)
  LOGICAL :: l_cloud_extinction                                   = .FALSE.
!   Calculate extinction of clouds (only solar)
  LOGICAL :: l_ls_cloud_absorptivity                              = .FALSE.
!   Calculate absorptivity of layer clouds (only infra-red)
  LOGICAL :: l_ls_cloud_extinction                                = .FALSE.
!   Calculate extinction of layer clouds (only solar)
  LOGICAL :: l_cnv_cloud_absorptivity                             = .FALSE.
!   Calculate absorptivity of conv.clouds (only infra-red)
  LOGICAL :: l_cnv_cloud_extinction                               = .FALSE.
!   Calculate extinction of conv.clouds (only solar)
  LOGICAL :: l_flux_direct_band                                   = .FALSE.
!   Calculate direct flux per band
  LOGICAL :: l_flux_direct_div_band                               = .FALSE.
!   Calculate direct flux divergence per band
  LOGICAL :: l_flux_direct_sph_band                               = .FALSE.
!   Calculate direct flux for spherical path per band
  LOGICAL :: l_flux_down_band                                     = .FALSE.
!   Calculate total downward flux per band
  LOGICAL :: l_flux_up_band                                       = .FALSE.
!   Calculate upward flux per band
  LOGICAL :: l_flux_div_band                                      = .FALSE.
!   Calculate flux divergence per band
  LOGICAL :: l_flux_direct_clear_band                             = .FALSE.
!   Calculate clear-sky direct flux per band
  LOGICAL :: l_flux_direct_clear_div_band                         = .FALSE.
!   Calculate clear-sky direct flux divergence per band
  LOGICAL :: l_flux_direct_clear_sph_band                         = .FALSE.
!   Calculate clear-sky direct flux for spherical path per band
  LOGICAL :: l_flux_down_clear_band                               = .FALSE.
!   Calculate clear-sky downward flux per band
  LOGICAL :: l_flux_up_clear_band                                 = .FALSE.
!   Calculate clear-sky upward flux per band
  LOGICAL :: l_flux_div_clear_band                                = .FALSE.
!   Calculate clear-sky flux divergence per band
  LOGICAL :: l_actinic_flux                                       = .FALSE.
!   Calculate actinic flux per channel
  LOGICAL :: l_actinic_flux_clear                                 = .FALSE.
!   Calculate clear-sky actinic flux per channel
  LOGICAL :: l_actinic_flux_band                                  = .FALSE.
!   Calculate actinic flux per band
  LOGICAL :: l_actinic_flux_clear_band                            = .FALSE.
!   Calculate clear-sky actinic flux per band
  LOGICAL :: l_photolysis_rate                                    = .FALSE.
!   Calculate photolysis rates per channel
  LOGICAL :: l_photolysis_rate_clear                              = .FALSE.
!   Calculate clear-sky photolysis rates per channel
  LOGICAL :: l_photolysis_div                                     = .FALSE.
!   Calculate flux divergence for photolysis
  LOGICAL :: l_photolysis_div_clear                               = .FALSE.
!   Calculate clear-sky flux divergence for photolysis
  LOGICAL :: l_aerosol_absorption_band                            = .FALSE.
!   Calculate total aerosol absorption per band
  LOGICAL :: l_aerosol_scattering_band                            = .FALSE.
!   Calculate total aerosol scattering per band
  LOGICAL :: l_aerosol_asymmetry_band                             = .FALSE.
!   Calculate total aerosol asymmetry (weighted by scattering) per band
  LOGICAL :: l_spherical_path_diag                                = .FALSE.
!   Output the direct beam path through spherical layers
  LOGICAL :: l_contrib_func                                       = .FALSE.
!   Output the contribution function
  LOGICAL :: l_contrib_func_band                                  = .FALSE.
!   Output the contribution function per band

! Satellite Data:

! Current provisions are largely intended for geostationary satellites.
! Note that in accordance with the UM's conventions these must be
! given in SI units; i.e. angles will be in radians etc.

  LOGICAL :: l_geostationary                                      = .FALSE.
!   Flag to signal that a geostationary satellite is assumed.
  CHARACTER (LEN=80) :: sat_desc                                  = ''
!   String for description of satellite
  REAL (RealK) :: sat_hgt                                         = rmdi
!   Height of the orbit above the Earth's surface
  REAL (RealK) :: sat_lon                                         = rmdi
!   Longitude of the (geostationary) satellite
  REAL (RealK) :: sat_lat                                         = rmdi
!   Latitude of the (geostationary) satellite (in practice, for
!   a geostationary satellite this must be 0.0)

! Viewing domain:
  REAL (RealK) :: max_view_lon                                    = rmdi
!   Maximum longitude of viewing domain
  REAL (RealK) :: min_view_lon                                    = rmdi
!   Minimum longitude of viewing domain
  REAL (RealK) :: max_view_lat                                    = rmdi
!   Maximum latitude of viewing domain
  REAL (RealK) :: min_view_lat                                    = rmdi
!   Minimum latitude of viewing domain

  ! Other stuff
      logical :: l_invert
    !   Flag to invert fields in the vertical
    logical :: l_profile_last
    !   Loop over profiles is last in input fields and diagnostics

    logical :: l_debug
    integer :: i_profile_debug
    !   Options for outputting debugging information

    ! Name of spectral file
    character(len=100) :: spectral_filename
    
  ! Define the namelist
namelist /socrates_rad_nml/ &
    isolir, first_band, last_band, &
    l_microphysics, l_gas, l_rayleigh, l_continuum, l_cont_gen, l_cloud, &
    l_drop, l_ice, l_aerosol, l_aerosol_mode, l_aerosol_ccn, &
    l_solar_tail_flux, l_orog, l_orog_fix, l_solvar, &
    i_gas_overlap, n_esft_red, gpnt_split, i_gas, &
    l_h2o, l_o2, l_n2o, l_ch4, l_cfc11, l_cfc12, l_cfc113, l_cfc114, l_hcfc22, &
    l_hfc125, l_hfc134a, l_co, l_nh3, l_tio, l_vo, l_h2, l_he, l_na, l_k, l_co2, &
    l_feh, l_crh, l_li, l_rb, l_cs, l_ph3, l_c2h2, l_hcn, l_o3, l_so2, l_no2, l_n2, &
    l_ar, l_o, l_n, l_hno3, l_no, l_no3, l_n2o5, l_hono, l_ho2no2, l_h2o2, l_c2h6, &
    l_ch3, l_h2co, l_ho2, &
    i_cloud, i_cloud_representation, i_st_water, i_cnv_water, i_st_ice, i_cnv_ice, &
    i_inhom, i_mcica_sampling, i_overlap, i_drop_re, l_local_cnv_partition, &
    l_global_cloud_top, l_avg_phase_fnc, &
    n_channel, i_angular_integration, i_2stream, i_solver, i_solver_clear, &
    n_order_gauss, i_truncation, i_sph_algorithm, n_order_phase_solar, ls_global_trunc, &
    ms_min, ms_max, ls_brdf_trunc, n_order_forward, i_sph_mode, i_scatter_method, &
    i_solar_src, i_direct_tau, l_ir_source_quad, l_rescale, l_spherical_solar, &
    l_spherical_diffuse, l_henyey_greenstein_pf, l_lanczos, l_euler_trnf, &
    accuracy_adaptive, euler_factor, half_angle, &
    l_tile, l_tile_emissivity, l_flux_ground, l_extra_top, l_rad_deg, l_subsample, &
    l_mixing_ratio, l_map_sub_bands, &
    i_scatter_method_band, i_gas_overlap_band, map_channel, weight_band, &
    weight_diag, l_clear_band, &
    l_clear, l_flux_div, l_blue_flux_surf, l_cloud_absorptivity, l_cloud_extinction, &
    l_ls_cloud_absorptivity, l_ls_cloud_extinction, l_cnv_cloud_absorptivity, &
    l_cnv_cloud_extinction, l_flux_direct_band, l_flux_direct_div_band, &
    l_flux_direct_sph_band, l_flux_down_band, l_flux_up_band, l_flux_div_band, &
    l_flux_direct_clear_band, l_flux_direct_clear_div_band, &
    l_flux_direct_clear_sph_band, l_flux_down_clear_band, l_flux_up_clear_band, &
    l_flux_div_clear_band, l_actinic_flux, l_actinic_flux_clear, l_actinic_flux_band, &
    l_actinic_flux_clear_band, l_photolysis_rate, l_photolysis_rate_clear, l_photolysis_div, &
    l_photolysis_div_clear, l_aerosol_absorption_band, l_aerosol_scattering_band, &
    l_aerosol_asymmetry_band, l_spherical_path_diag, l_contrib_func, l_contrib_func_band, &
    l_geostationary, sat_desc, sat_hgt, sat_lon, sat_lat, max_view_lon, min_view_lon, &
    max_view_lat, min_view_lat, &
    l_invert, l_profile_last, l_debug, i_profile_debug, &
    spectral_filename


end module socrates_config_mod
