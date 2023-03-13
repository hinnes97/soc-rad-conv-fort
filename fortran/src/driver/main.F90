program main
  use io
  use utils
  use params
  use timestep
  use condense, only: dew_point_T
  use init_pt_mod, only: init_pt

#ifdef SOC
  use socrates_interface_mod, only : socrates_init
#elif defined PICKET
  use radiation_mod, only : radiation_init
#elif defined TWOSTR
  use band_grey_mod, only : bg_init
#endif


  implicit none

  real(dp), dimension(:), allocatable :: Tf, pf ! Temp and pressure arrays
  real(dp), dimension(:), allocatable :: Te, pe     ! Edge pressure array

  real(dp), dimension(:), allocatable :: q
  real(dp) ::  start, end, Ts

  integer :: ncid


  ! Initialise parameters and allocate arrays
  call read_constants()
  call allocate_arrays(Tf, pf, pe, Te, q)

#ifdef SOC
  call socrates_init()
#elif defined PICKET
  call radiation_init(nf)
#elif defined TWOSTR
  call bg_init(nf)
#endif

  !Initialise output file
  call file_setup(output_file, nf, ne, ncid)

  if (init_from_file .eqv. .true.) then
     ! Read data from output file
     call read_initial_data(input_file, Tf, Te, q, pf, pe, Ts)
     write(*,*) 'INITIAL TS from file', Ts, Te(ne)
  else

     call init_pt(pf, pe, Tf, Te)
     Ts = Te(ne)
  endif
  
  ! Do timestepping
  call cpu_time(start)
  call step(Tf, pf, pe, output_file ,q, Ts)
  call cpu_time(end)
  
  call deallocate_arrays(Tf, pf, pe,Te)
  
end program main
