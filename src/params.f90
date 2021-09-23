module params

  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter:: dp = REAL64

  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: sb = 5.670374419e-8_dp

  ! Params read in via namelist
  integer :: nf = 100    ! Number of mid-levels
  real(dp) :: log_top_p = 1._dp
  real(dp) :: log_bot_p = 6._dp
  real(dp) :: top_t = 200._dp
  real(dp) :: bot_t = 500._dp
  real(dp) :: tau_V_inf = 256._dp*0.04_dp
  real(dp) :: tau_IR_inf = 256._dp
  real(dp) :: Finc = 1368.0_dp/4._dp
  real(dp) :: Fint = sb*70._dp**4._dp
  real(dp) :: const = 0.0001_dp
  integer :: Nt = 1000000
  integer :: N_max = 20
  integer :: rad_scheme = 1
  integer :: ne
  
  namelist /param_nml/nf, log_top_p, log_bot_p, top_t, bot_t, tau_V_inf, tau_IR_inf, Finc, &
       Fint, const, Nt, N_max, rad_scheme
  
contains

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
       open(f_unit, file=filename)
       read(f_unit, param_nml, iostat=ios)

       if (ios .gt. 0) then
          write(*,*) "Error reading ", trim(filename), "iostat=", ios
          stop
       endif       
    end if

    ne = nf + 1
  end subroutine read_constants

  subroutine allocate_arrays(Tf, pf, pe, tau_IR, tau_V, net_F, dT)
    real(dp), dimension(:), allocatable, intent(inout) :: Tf, pf, pe, tau_IR, tau_V, net_F, dT
    
    allocate(Tf(nf))
    allocate(pf(nf))
    allocate(pe(ne))
    allocate(tau_IR(ne))
    allocate(tau_V(ne))
    allocate(net_F(ne))
    allocate(dT(nf))
    
  end subroutine allocate_arrays

  subroutine deallocate_arrays(Tf, pf, pe, tau_IR, tau_V, net_F, dT)
    real(dp), dimension(:), allocatable, intent(inout) :: Tf, pf, pe, tau_IR, tau_V, net_F, dT
    deallocate(Tf)
    deallocate(pf)
    deallocate(pe)
    deallocate(tau_IR)
    deallocate(tau_V)
    deallocate(net_F)
    deallocate(dT)
    
  end subroutine deallocate_arrays

end module params
