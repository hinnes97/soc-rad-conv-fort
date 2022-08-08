program main
  use io
  use utils
  use params
  use timestep
  use matrix
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
  real(dp), dimension(:), allocatable :: q, fdn, fup, s_dn, s_up
  real(dp) :: olr, start, end, Ts

  integer :: ncid

  integer :: i

  ! Initialise parameters and allocate arrays
  write(*,*) 'BEFORE READ CONSTANTS'
  call read_constants()
  call allocate_arrays(Tf, pf, pe, net_F, dT, Te, q, fup, fdn, s_dn, s_up)

#ifdef SOC
  call socrates_init()
#elif defined PICKET
  call radiation_init(nf)
#elif defined TWOSTR
  call bg_init(nf)
#endif

  write(*,*) 'HERE'
  !Initialise output file
  call file_setup(output_file, nf, ne, ncid)

  if (init_from_file .eqv. .true.) then
     call read_initial_data(input_file, Tf, Te, q, Ts)
     call logspace(log_top_p, log_bot_p, pe)
     ! Initialise pf array from pe
     do i=1,nf
        pf(i) = (pe(i+1) - pe(i)) / (log(pe(i+1)) - log(pe(i)))
     end do
     !Ts = Te(ne)
     write(*,*) 'INITIAL TS from file', Ts, Te(ne)
     ! Smooth input
     !do i=2,nf-1
     !   Tf(i) = Tf(i-1)*0.25 + Tf(i)*0.5 + Tf(i+1)*0.25
     !bnddo

     !do i=2,nf-1
     !   call bezier_interp(pf(i-1:i+1), Tf(i-1:i+1), 3, pe(i), Te(i))
        !call linear_log_interp(pe(i), pf(i-1), pf(i), Tf(i-1), Tf(i), Te(i))
     !enddo
     !call bezier_interp(pf(nf-2:nf), Tf(nf-2:nf), 3, pe(nf), Te(nf))
     !call linear_log_interp(pe(1), pf(1), pf(2), Tf(1), Tf(2), Te(1))
     !call linear_log_interp(pe(ne), pf(nf-1), pf(nf), Tf(nf-1), Tf(nf), Te(ne) )

  else
     ! Initialise pressure and temperature arrays as log and lin spaced respectively
     select case(p_grid)
     case('log')
        call logspace(log_top_p, log_bot_p, pe)
     case('hires_trop')
        ! Let section of the atmosphere between ps and ps/10 contain a larger proportion of points
        ! Linearly spaced in this region
        !call linspace(9./10. * 10**log_bot_p, 10**(log_bot_p),pe((ne - ne/frac):ne))
        !call logspace(log_top_p, 9./10.*log_bot_p, pe(1:(ne-ne/frac-1)), .false.)
        
        call linspace(9./10.* 10.**log_bot_p, 10.**(log_bot_p), pe((ne - ne/frac):ne))
        call logspace(log_top_p, log10(0.9) + log_bot_p, pe(1:(ne-ne/frac-1)), .false.)
     end select

     ! Initialise pf array from pe
     do i=1,nf
        pf(i) = (pe(i+1) - pe(i)) / (log(pe(i+1)) - log(pe(i)))
     end do

     do i=1,nf
        Tf(i) = bot_t*(pf(i)/pe(ne))**(2./7.)

     enddo

     do i=1,ne
        Te(i) = bot_t*(pe(i)/pe(ne))**(2./7.)
     enddo

    !call linspace(top_t, bot_t, Tf)
    !call linspace(top_t, bot_t, Te)
     do i=1,nf
        if (Tf(i) .lt. top_t) then
           Tf(i) = top_t
        endif
     enddo

     do i=1,nf+1
        if (Te(i) .lt. top_t) then
           Te(i) = top_t
        endif
     enddo
     write(*,*) (Tf(i), i=1,nf)
     write(*,*) (Te(i), i=1,nf+1)
     !q = 1.

     Ts = Te(ne)
     write(*,*) 'Ts', Ts
  endif
  
  write(*,*) 'Before start'
  if (matrix_rt) then
     write(*,*) 'matrix method'
     ! Do matrix method
     call do_matrix(nf, ne, Tf, pf, Te, pe, 1.0_dp, Finc, Fint, olr,q, Ts)
  else
     ! Do timestepping
     write(*,*) 'Timestepping'
     call cpu_time(start)
     call step(Tf, pf, pe, net_F, dT, olr, output_file ,q, fup, fdn, s_dn, s_up,Ts)
     call cpu_time(end)
     write(*,*) 'TIME ELAPSED: ', end - start

     ! Interpolate to Te
     do i=2,nf-1
        call bezier_interp(pf(i-1:i+1), Tf(i-1:i+1), 3, pe(i), Te(i))
        !call linear_log_interp(pe(i), pf(i-1), pf(i), Tf(i-1), Tf(i), Te(i))
     enddo
     call bezier_interp(pf(nf-2:nf), Tf(nf-2:nf), 3, pe(nf), Te(nf))
     call linear_log_interp(pe(1), pf(1), pf(2), Tf(1), Tf(2), Te(1))
     call linear_log_interp(pe(ne), pf(nf-1), pf(nf), Tf(nf-1), Tf(nf), Te(ne) )
     Te(ne) = Ts
     write(*,*) Tf
  endif
  
  call dump_data(output_file, nf, ne, Tf, pf, pe, olr, Finc, Fint,Te, q, fup, fdn, s_dn, s_up, Ts)
  call deallocate_arrays(Tf, pf, pe, net_F, dT,Te, s_dn, s_up)
  
end program main
