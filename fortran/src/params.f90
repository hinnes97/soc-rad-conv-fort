module params

  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter:: dp = REAL64

  ! Physical constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: deg_to_rad = pi/180._dp

  ! Stefan boltzmann constant (W/m^2/K^4)
  real(dp), parameter :: sb = 5.670374419e-8_dp
  ! Gas constant (J/K/mol)
  real(dp), parameter :: rgas = 8.31446261815324_dp
  ! Speed of light (m/s)
  real(dp), parameter :: c = 299792458._dp
  ! Planck constant (m^2 kg/s)
  real(dp), parameter :: h = 6.62607004e-34
  ! Boltzmann constant (m^2 kg/s^2/K)
  real(dp), parameter :: kb = 1.38064852e-23

  ! Params read in via namelist

  ! -----------------------------------------------------------------------------
  !                           CONTROL
  !------------------------------------------------------------------------------
  ! Number of mid-levels
  integer :: nf = 100
  ! Number of edge levels (nf+1)
  integer :: ne 
  ! Whether to do matrix method or to timestep
  logical :: matrix_rt=.false.
  ! Whether to include a surface or not
  logical :: surface = .false.
  
  ! -----------------------------------------------------------------------------
  !                           INITIALISATION
  !------------------------------------------------------------------------------
  ! Log10 top pressure
  real(dp) :: log_top_p = 1._dp
  ! Log10 bottom pressure
  real(dp) :: log_bot_p = 6._dp
  ! TOA initial temperature
  real(dp) :: top_t = 200._dp
  ! BOA initial temperature
  real(dp) :: bot_t = 500._dp
  ! Type of pressure grid
  character(len=20) :: p_grid = 'log'
  ! Fraction of grid to be in hires trop
  integer :: frac = 3

  ! -----------------------------------------------------------------------------
  !                           RADIATION
  !------------------------------------------------------------------------------
  ! Incoming instellation
  real(dp) :: Finc = 1368.0_dp/4._dp
  ! Internal heat flux
  real(dp) :: Fint = sb*70._dp**4._dp

  ! -----------------------------------------------------------------------------
  !                           BAND GREY
  !------------------------------------------------------------------------------
  ! Directory where opacities are kept
  character(80) :: opacity_dir
  ! Fraction to alter SW bands by
  real(dp) :: sw_fac
  ! Fraction to alter LW bands by
  real(dp) :: lw_fac
  ! Whether to do grid inversion for stability
  logical :: invert_grid = .false.
  

  ! -----------------------------------------------------------------------------
  !                           SEMI-GREY (SHORT CHAR)
  !------------------------------------------------------------------------------
  ! Semi grey scheme selector (short_char or toon right now)
  character(30) :: semi_grey_scheme
  ! LW opacity in m^2/kg
  real(dp) :: kappa_LW = 1.e-3
  ! SW opacity in m^2/kg
  real(dp) :: kappa_SW = 1.e-5
  ! Moist radiation flag
  logical :: moist_rad = .false.
  ! Opacity for water vapour (m^2/kg)
  real(dp) :: kappa_q = 0.00789463795397635
  
  ! -----------------------------------------------------------------------------
  !                           TIMESTEPPING
  !------------------------------------------------------------------------------
  ! Number of timesteps
  integer :: Nt = 1000000
  ! Acceleration factor
  real(dp) :: const = 0.0001_dp
  ! Timestep
  real(dp) :: del_time
  ! Whether to use acceleration factor
  logical :: accelerate = .false.

  ! -----------------------------------------------------------------------------
  !                           MATRIX RT
  !------------------------------------------------------------------------------
  ! How much to weight prior temperature field
  real(dp) :: alpha = 1._dp
  ! Error in "measurements"
  real(dp) :: error_frac = 0.01
  ! How many iterations to do
  integer :: mat_iters = 20
  
  ! -----------------------------------------------------------------------------
  !                           CONVECTION
  !------------------------------------------------------------------------------
  ! Convection scheme (no_conv, dry_conv, moist_conv)
  character(80) :: conv_scheme
  ! Number of upwards/downwards passes during convection code
  integer :: passes
  ! Whether to include moisture inhibition
  logical :: inhibited = .false.
  
  ! -----------------------------------------------------------------------------
  !                           ATMOSPHERIC PARAMETERS
  !------------------------------------------------------------------------------
  ! Gas constant for dry air (J/kg/K)
  real(dp) :: rdgas = 3779._dp ! Solar metallicity
  ! Gravity (m/s^2)
  real(dp) :: grav = 12.43_dp
  ! Heat capcity (J/kg/K)
  real(dp) :: cpair = 3779._dp*7._dp/2._dp
  ! Dry lapse rate
  real(dp) :: Rcp
  

  ! -----------------------------------------------------------------------------
  !                           MOISTURE
  !------------------------------------------------------------------------------
  ! Moisture scheme (none, deep, surface)
  character(80) :: moisture_scheme
  ! Deep atmosphere moisture level
  real(dp) :: q0 = 0.01_dp
  
  ! -----------------------------------------------------------------------------
  !                           INPUT/OUTPUT
  !------------------------------------------------------------------------------

  ! Switch for whether input file is used on initialisation
  logical :: init_from_file
  ! Name of output file
  character(60) :: output_file
  ! Name of input file
  character(60) :: input_file

  ! -----------------------------------------------------------------------------
  !                           SURFACE
  !------------------------------------------------------------------------------

  ! Surface heat capacity (J/kg/K)
  real(dp) :: cp_s = 3989. ! Kaspi and Showman value
  ! Surface density [kg m^-3]
  real(dp) :: rho_s = 1035 ! Kaspi and Showman (2015) value
  ! Surface albedo
  real(dp) :: A_s = 0.3
  ! Const to multiply fluxes by
  real(dp) :: surf_const = 0.01
  ! Turbulent heat exchange coefficient
  real(dp) :: C_d = 0.001
  ! Depth of water ocean [m]
  real(dp) :: depth = 1.
  ! Characteristic surface wind speed
  real(dp) :: U = 10.! Pierrehumbert 2011 value
  
  namelist /control_nml/ nf, matrix_rt, surface
  namelist /initialisation_nml/ log_top_p, log_bot_p, bot_t, top_t, p_grid, frac
  namelist /io_nml/ init_from_file, input_file, output_file
  namelist /param_nml/ rdgas, grav, cpair, Rcp
  namelist /timestep_nml/ Nt, const, del_time, accelerate
  namelist /matrix_nml/ mat_iters, alpha, error_frac
  namelist /convection_nml/ conv_scheme, passes, inhibited
  namelist /radiation_nml/ Finc, Fint 
  namelist /band_grey_nml/ opacity_dir, invert_grid,sw_fac, lw_fac
  namelist /semi_grey_nml/ kappa_lw, kappa_sw, moist_rad, kappa_q, semi_grey_scheme
  namelist /moisture_nml/ moisture_scheme, q0
  namelist /surface_nml/ cp_s, A_s, surf_const, C_d, depth, U, rho_s
  
  
contains

  subroutine check_error(ios, filename, nml)
    integer, intent(in)  :: ios
    character(len=*), intent(in) :: filename, nml

    if (ios .gt. 0) then
       write(*,*) "Error reading ", trim(filename), "namelist ", trim(nml), "iostat=", ios
       stop
    endif

  end subroutine check_error
  
  subroutine read_constants()
    character(len=60) :: filename = "input.nml"
    logical :: exists
    integer :: f_unit = 1
    integer :: ios
  
    inquire(file=filename, exist=exists)
    if (.not. exists) then
       write(*,*) "file ", trim(filename), " doesn't exist"
       stop    
    else
       ! Read in all the different namelists

       open(f_unit, file=filename)

       ! Atmospheric parameters
       read(f_unit, param_nml, iostat=ios)
       rewind(f_unit)
       call check_error(ios, filename, "param_nml")

       ! Control namelist
       read(f_unit, control_nml, iostat=ios)
       rewind(f_unit)
       call check_error(ios, filename, "control_nml")

       ! Initialisation
       read(f_unit, initialisation_nml, iostat=ios)
       rewind(f_unit)
       call check_error(ios, filename, "initialisation_nml")

       ! Input/output
       read(f_unit, io_nml, iostat=ios)
       rewind(f_unit)
       call check_error(ios, filename, "io_nml")

       ! Timestep namelist
       read(f_unit, timestep_nml, iostat=ios)
       rewind(f_unit)
       call check_error(ios, filename, "timestep_nml")

       ! Matrix RT
       read(f_unit, matrix_nml, iostat=ios)
       rewind(f_unit)
       call check_error(ios, filename, "matrix_nml")

       ! Convection
       read(f_unit, convection_nml, iostat=ios)
       rewind(f_unit)
       call check_error(ios, filename, "convection_nml")

       ! Radiation
       read(f_unit, radiation_nml, iostat=ios)
       rewind(f_unit)
       call check_error(ios, filename, "radiation_nml")

       ! Band grey scheme
       read(f_unit, band_grey_nml, iostat=ios)
       rewind(f_unit)
       call check_error(ios, filename, "band_grey_nml")

       ! Semi-grey scheme
       read(f_unit, semi_grey_nml, iostat=ios)
       rewind(f_unit)
       call check_error(ios, filename, "semi_grey_nml")

       ! Moisture
       read(f_unit, moisture_nml, iostat=ios)
       rewind(f_unit)
       call check_error(ios, filename, "moisture_nml")

       ! Surface
       read(f_unit, surface_nml, iostat=ios)
       rewind(f_unit)
       call check_error(ios, filename, "surface_nml")

       close(f_unit)
    end if

    ne = nf + 1
    Rcp = rdgas/cpair

    
  end subroutine read_constants

  subroutine allocate_arrays(Tf, pf, pe, net_F, dT, Te, q, fup, fdn, s_dn, s_up)
    real(dp), dimension(:), allocatable, intent(inout) :: Tf, pf, pe,  net_F, dT,Te,q, &
         fup, fdn, s_dn, s_up
    
    allocate(Tf(nf), pf(nf), pe(ne),  net_F(ne), dT(nf), Te(ne), q(nf), &
         fup(ne), fdn(ne), s_dn(ne), s_up(ne))
    
  end subroutine allocate_arrays

  subroutine deallocate_arrays(Tf, pf, pe, net_F, dT,Te, s_dn, s_up)
    real(dp), dimension(:), allocatable, intent(inout) :: Tf, pf, pe,&
         & net_F, dT,Te, s_dn, s_up
    deallocate(Tf,pf,pe,net_F,dT,Te, s_dn, s_up)
    
  end subroutine deallocate_arrays

end module params
