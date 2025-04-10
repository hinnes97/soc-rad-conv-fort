program main
  use io
  use utils
  use params
  use timestep
  use condense, only: dew_point_T
  use eos_dry_adj, only: init_data
  use init_pt_mod, only: init_pt
  use mat_step_mod, only: mat_step
  use atmosphere, only : read_abundances, init_atmos
  
#ifdef SOC
  !use socrates_interface_mod, only : soc_init
  use soc_init_mod, only : soc_init
#elif defined PICKET
  use radiation_mod, only : radiation_init
#elif defined TWOSTR
  use band_grey_mod, only : bg_init
#endif


  implicit none

  real(dp), dimension(:), allocatable :: Tf, pf ! Temp and pressure arrays
  real(dp), dimension(:), allocatable :: Te, pe     ! Edge pressure array

  real(dp), dimension(:,:), allocatable :: q
  real(dp) ::  start, end, Ts

  integer :: ncid, i


! Initialise parameters and allocate arrays
  call read_constants()
  call allocate_arrays(Tf, pf, pe, Te)

#ifdef SOC
  !call socrates_init()
  call soc_init(nf, ne)
#elif defined PICKET
  call radiation_init(nf)
#elif defined TWOSTR
  call bg_init(nf)
#endif

  !Initialise output file
  call file_setup(output_file, nf, ne, nqr, ncid)

  if (init_from_file .eqv. .true.) then
     ! Read data from output file
     call read_initial_data(input_file, Tf, Te, q, pf, pe, Ts)
     write(*,*) 'INITIAL TS from file', Ts, Te(ne)
  else

     call init_pt(pf, pe, Tf, Te)
     Ts = Te(ne)
  endif

! Read abundances from file
! Initialise atmosphere
  call init_atmos(pf, q)
  call init_data ! For dry adjust
  !call read_abundances(pf, q)
  
  ! Do timestepping
  call cpu_time(start)
  write(*,*) 'before step'
  do i=1,nf
     write(*,*) i, sum(q(i,:))
  enddo
  !call step(Tf, pf, pe, output_file ,q, Ts)
  call mat_step(Tf, pf, pe, output_file, q, Ts)
  call cpu_time(end)
  
  call deallocate_arrays(Tf, pf, pe,Te, q)
  
end program main
