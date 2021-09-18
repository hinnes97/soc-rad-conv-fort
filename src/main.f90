program main
  use radiation_Kitzmann_noscatt, only: Kitzmann_TS_noscatt
  !use, intrinsic  :: iso_fortran_env
  use io
  use utils
  use params

  implicit none

  ! Precision variable
  !integer, parameter :: dp = real64
  
  ! Set up arrays
  !integer, parameter :: nf = 100    ! Number of mid-levels
  !integer, parameter :: ne = nf + 1 ! Number of edge levels
  !real(dp), parameter :: log_top_p = 1._dp
  !real(dp), parameter :: log_bot_p = 6._dp
  !real(dp), parameter :: top_t = 200._dp
  !real(dp), parameter :: bot_t = 500._dp
  !real(dp), parameter :: tau_V_inf = 256._dp*0.04_dp
  !real(dp), parameter :: tau_IR_inf = 256._dp
  !real(dp), parameter :: Finc = 1368.0_dp/4._dp
  !real(dp), parameter :: Fint = sb*70._dp**4._dp
  !real(dp) :: const = 0.0001_dp

  !integer, parameter :: Nt = 1000000 ! Number of timesteps

  real(dp), dimension(:), allocatable :: Tf, pf ! Temp and pressure arrays
  real(dp), dimension(:), allocatable :: pe     ! Edge pressure array
  real(dp), dimension(:), allocatable :: tau_IR, tau_V

  real(dp), dimension(:), allocatable :: net_F
  real(dp), dimension(:), allocatable :: dT
  real(dp) :: olr

  integer :: ncid

  integer :: i,j

  ! Initialise parameters and allocate arrays
  call read_constants()
  call allocate_arrays(Tf, pf, pe, tau_IR, tau_V, net_F, dT)
  
  !Initialise output file
  call file_setup("output.nc", nf, ne, ncid)
  
  ! Initialise pressure and temperature arrays as log and lin spaced respectively
  call logspace(log_top_p, log_bot_p, pe)
  call linspace(top_t, bot_t, Tf)
  
  ! Initialise pf array from pe
  do i=1,nf
     pf(i) = (pe(i+1) - pe(i)) / (log(pe(i+1)) - log(pe(i)))
  end do

  write(*,*) maxval(pe), maxval(pf)
  
  ! Initialise tau arrays
  do i=1, ne
     tau_V(i) = tau_V_inf*(pe(i)/pe(ne))
     tau_IR(i)= tau_IR_inf*(pe(i)/pe(ne))
  end do

  do j =1,Nt
     if (mod(j, 100000) .eq. 0) then
        write(*,*) Tf
     end if
     
     call Kitzmann_TS_noscatt(nf, ne, Tf, pf, pe, tau_IR, tau_V, &
          net_F, 1._dp, Finc, Fint, olr)

     do i=1,nf
        dT(i) = const*(net_F(i+1) - net_F(i))
        if (dT(i)>5.0_dp) then
           dT(i) = 5.0_dp
        endif
        if (dT(i)<-5.0_dp) then
           dT(i) = -5.0_dp
        endif
     end do
     
     call Kitzmann_TS_noscatt(nf, ne, Tf+0.5_dp*dT, pf, pe, tau_IR, tau_V, &
          net_F, 1._dp, Finc, Fint, olr)

     do i=1,nf
        dT(i) = const*(net_F(i+1) - net_F(i))
        if (dT(i)>5.0_dp) then
           dT(i) = 5.0_dp
        endif
        if (dT(i)<-5.0_dp) then
           dT(i) = -5.0_dp
        endif
        Tf(i) = Tf(i) + dT(i)
        if (Tf(i) < 100._dp) then
           Tf(i) = 100._dp
        end if        
     end do

!     if (maxval(abs(dT)) < 0.1 .and. const<0.1) then
!        const = const*1.5_dp
!        write(*,*) const
!     end if
!     if (minval(abs(dT)) > 3) then
!        const = const/2._dp
!     end if
     
  end do 
  
  call dump_data(ncid, nf, ne, Tf, pf, pe, olr, tau_IR_inf, tau_V_inf, Finc, Fint)
  call close_file(ncid)
  call deallocate_arrays(Tf, pf, pe, tau_IR, tau_V, net_F, dT)
  
end program main
