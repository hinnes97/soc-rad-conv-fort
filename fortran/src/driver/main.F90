program main
  use io
  use utils
  use params
  use timestep
  use condense, only: dew_point_T

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
  real(dp), dimension(:), allocatable :: tau_IR, tau_V

  real(dp), dimension(:), allocatable :: net_F
  real(dp), dimension(:), allocatable :: dT
  real(dp), dimension(:), allocatable :: q
  real(dp) :: olr, start, end, Ts

  integer :: ncid

  integer :: i

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
     select case(p_grid)
     case('log')
        ! Logarithmically spaced pressure grid
        call logspace(log_top_p, log_bot_p, pe)
     case('hires_trop')
        ! Let section of the atmosphere between ps and ps/10 contain a larger proportion of points
        ! Linearly spaced in this region
        call logspace(log10(0.85)+log_bot_p, log_bot_p, pe((ne-ne/frac):ne))
        call logspace(log_top_p, log10(0.85) + log_bot_p, pe(1:(ne-ne/frac-1)), .false.)
     end select

     ! Initialise pf array from pe
     do i=1,nf
        pf(i) = (pe(i+1) - pe(i)) / (log(pe(i+1)) - log(pe(i)))
     end do

     ! Initialise temperature on a dry adiabat
     do i=1,nf
        Tf(i) = bot_t*(pf(i)/pe(ne))**(2./7.)
     enddo

     do i=1,ne
        Te(i) = bot_t*(pe(i)/pe(ne))**(2./7.)
     enddo

     ! Limit minimum temperature of atmosphere
     do i=1,nf
        Tf(i) = max(Tf(i), top_t)
     enddo

     do i=1,nf+1
        Te(i) = max(Te(i), top_t)
     enddo

     Ts = Te(ne)
  endif
  
  ! Do timestepping
  call cpu_time(start)
  call step(Tf, pf, pe, output_file ,q, Ts)
  call cpu_time(end)
  
  call deallocate_arrays(Tf, pf, pe,Te)
  
end program main
