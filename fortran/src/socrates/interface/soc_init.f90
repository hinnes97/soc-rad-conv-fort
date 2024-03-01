module soc_init_mod
  use socrates_config_mod
  use socrates_set_spectrum, only: set_spectrum
  use def_spectrum, only: StrSpecData
  use def_mcica,    only: StrMcica
  use def_control,  only: StrCtrl,  deallocate_control
  use def_dimen,    only: StrDim
  use def_atm,      only: StrAtm,   deallocate_atm
  use def_bound,    only: StrBound, deallocate_bound
  use def_cld,      only: StrCld,   deallocate_cld, deallocate_cld_prsc, &
                                  deallocate_cld_mcica
  use def_aer,      only: StrAer,   deallocate_aer, deallocate_aer_prsc
  use def_out,      only: StrOut,   deallocate_out
  use socrates_def_diag,     only : StrDiag
  
  use socrates_set_spectrum, only: spectrum_array_name, spectrum_array, &
       mcica_spectrum_name, mcica_data_array, &
       set_spectrum

  use socrates_set_control,   only: set_control
  use socrates_set_dimen,     only: set_dimen
  use socrates_set_atm,       only: set_atm
  use socrates_set_bound,     only: set_bound
  use socrates_set_cld,       only: set_cld
  use socrates_set_cld_dim,   only: set_cld_dim
  use socrates_set_cld_mcica, only: set_cld_mcica
  use socrates_set_aer,       only: set_aer
  use socrates_set_diag,      only: set_diag

  use realtype_rd, only: RealExt
  use ereport_mod, only: ereport
  use errormessagelength_mod, only: errormessagelength
  use rad_pcf, only: i_normal, i_err_fatal, ip_solar, ip_infra_red

  implicit none

  type(StrDiag) :: diag_sw, diag_lw
  type(StrDim) :: dimen_sw,dimen_lw
  type(StrCtrl)  :: control_sw, control_lw

  ! Spectral data:
  type(StrSpecData), pointer :: spec => null()
  type(StrSpecData), target  :: spectrum
  ! Mcica data:
  type(StrMcica), pointer :: mcica => null()
  type(StrMcica), target :: mcica_dummy


  integer :: n_aer_mode, n_tile, n_cloud_layer, n_profile
contains
  
  subroutine soc_init(nf, ne)
    integer, intent(in) :: nf, ne
    
    character(len=100) :: filename
    logical :: file_exist
    integer :: io, unit

    n_aer_mode = 1; n_tile = 1; n_cloud_layer = 0; n_profile=1
    
    filename='input.nml'
    
    inquire(file=filename, exist=file_exist)
    if ( file_exist ) then
       write(*,*) 'FILE EXISTS'
       open(newunit=unit, file=filename, iostat=io)
       if (io .ne. 0) write(*,*) 'OPEN ERROR: IOSTAT=',io
       rewind(unit)
       read (unit, socrates_rad_nml, iostat=io)
       if (io .ne. 0) write(*,*) 'READ ERROR: IOSTAT=',io
       !rewind(unit)
       !read (unit,radiation_nml, iostat=io)
       !if (io .ne. 0) write(*,*) 'READ ERROR: IOSTAT=',io
       close(unit)
       if (io .ne. 0) write(*,*) 'CLOSE ERROR: IOSTAT=',io
    endif

    ! LW
    call set_spectrum(spectrum=spectrum, &
                      spectral_file=spectral_filename,&
                      l_h2o=l_h2o, &
                      l_co2=l_co2, &
                      l_o3 = l_o3,&
                      l_o2 = l_o2, &
                      l_n2o=l_n2o, &
                      l_ch4=l_ch4, & 
                      l_so2=l_so2, &
                      l_cfc11=l_cfc11, &
                      l_cfc12=l_cfc12, &
                      l_cfc113=l_cfc113, &
                      l_cfc114=l_cfc114, &
                      l_hcfc22=l_hcfc22, &
                      l_hfc125=l_hfc125,&
                      l_hfc134a=l_hfc134a,&
                      l_co=l_co,&
                      l_nh3=l_nh3, &
                      l_tio=l_tio, &
                      l_vo=l_vo, &
                      l_h2= l_h2, &
                      l_he = l_he,&
                      l_na = l_na,  &
                      l_k = l_k,&
                      l_li = l_li,&
                      l_rb = l_rb, &
                      l_cs = l_cs, &
                      l_all_gases = l_all_gases)


        spec => spectrum
        mcica => mcica_dummy
        ! SW
        call set_control(control_sw, diag_sw, spec, &
             isolir                 = ip_solar, &
             l_rayleigh             = l_rayleigh, &
             l_orog                 = l_orog, &
             l_mixing_ratio         = l_mixing_ratio, &
             l_aerosol_mode         = l_aerosol_mode, &
             l_tile                 = l_tile, &
             l_flux_ground          = l_flux_ground, &
             l_flux_tile            = l_flux_tile, &
             n_tile                 = n_tile, &
             n_cloud_layer          = n_cloud_layer, &
             n_aer_mode             = n_aer_mode, &
             i_cloud_representation = i_cloud_representation, &
             i_overlap              = i_overlap, &
             i_inhom                = i_inhom, &
             i_mcica_sampling       = i_mcica_sampling, &
             i_st_water             = i_st_water, &
             i_cnv_water            = i_cnv_water, &
             i_st_ice               = i_st_ice, &
             i_cnv_ice              = i_cnv_ice, &
             i_drop_re              = i_drop_re, &
             l_set_defaults         = .true.)

        
        call set_control(control_lw, diag_lw, spec, &
             isolir                 = ip_infra_red, &
             l_rayleigh             = l_rayleigh, &
             l_orog                 = l_orog, &
             l_mixing_ratio         = l_mixing_ratio, &
             l_aerosol_mode         = l_aerosol_mode, &
             l_tile                 = l_tile, &
             l_flux_ground          = l_flux_ground, &
             l_flux_tile            = l_flux_tile, &
             n_tile                 = n_tile, &
             n_cloud_layer          = n_cloud_layer, &
             n_aer_mode             = n_aer_mode, &
             i_cloud_representation = i_cloud_representation, &
             i_overlap              = i_overlap, &
             i_inhom                = i_inhom, &
             i_mcica_sampling       = i_mcica_sampling, &
             i_st_water             = i_st_water, &
             i_cnv_water            = i_cnv_water, &
             i_st_ice               = i_st_ice, &
             i_cnv_ice              = i_cnv_ice, &
             i_drop_re              = i_drop_re, &
             l_set_defaults         = .true.)


        call set_dimen(dimen_lw, control_lw, n_profile, nf, &
             mcica_data    = mcica, &
             n_tile        = n_tile, &
             n_cloud_layer = n_cloud_layer, &
             n_aer_mode    = n_aer_mode )
        
        call set_dimen(dimen_sw, control_sw, n_profile, nf, &
             mcica_data    = mcica, &
             n_tile        = n_tile, &
             n_cloud_layer = n_cloud_layer, &
             n_aer_mode    = n_aer_mode )

   end subroutine soc_init
end module soc_init_mod
