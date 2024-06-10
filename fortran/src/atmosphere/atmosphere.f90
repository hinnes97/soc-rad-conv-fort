module atmosphere

  use params, only: soc_index_file, dp,abundance_file, moisture_scheme,aqua_path
  use utils, only: linear_log_interp
  use phys
  use aqua_eos, only: load_table_pt
  implicit none

  integer :: nqr, nqt
  integer, allocatable :: soc_indices(:), th_indices(:)
  real(dp), allocatable :: mmw_dry(:), cp_dry(:), q_orig(:,:)
  type(gas),allocatable :: th_gases(:)
  type(gas), allocatable :: rad_gases(:)
  
  !real, allocatable :: mmw(:), cp(:), Rgas(:)
  
contains

  subroutine init_atmos(p, q)
    !real(dp), intent(in), allocatable, dimension(:) :: mmw, cp, Rgas
    real(dp), intent(in) :: p(:)
    real(dp), intent(inout), allocatable :: q(:,:)
    
    integer :: f_unit, ierr, i

    call init_gas_data
    
! Open file
    open(newunit=f_unit, iostat=ierr, action='read', file=soc_index_file)
    if (ierr .ne. 0) then
       write(*,*) 'Could not open ', soc_index_file
       stop
    endif

    read(f_unit, *) nqr, nqt
    allocate(soc_indices(nqr), rad_gases(nqr), th_indices(nqt), th_gases(nqt))
    
    read(f_unit,*) th_indices(1:nqt)
    
    close(f_unit)
    soc_indices(1:nqr) = th_indices(1:nqr)
    

    do i=1,nqt
       select case(th_indices(i))
       case(1)
          th_gases(i) = H2O
       case(2)
          th_gases(i) = CO2
       case(5)
          th_gases(i) = CO
       case(6)
          th_gases(i) = CH4
       case(11)
          th_gases(i) = NH3
       case(13)
          th_gases(i) = N2
       case(23)
          th_gases(i) = H2
       case(24)
          th_gases(i) = He
       case(35)
          th_gases(i) = HCN
       case(46)
          th_gases(i) = C2H6
       case default
          write(*,*) 'ERROR: incorrect species identifier'
          stop
       end select
    enddo

    call read_abundances(p,q)

    if (moisture_scheme == 'supercrit') call aqua_init
  end subroutine init_atmos

  subroutine aqua_init
    call load_table_pt(aqua_path)
  end subroutine aqua_init

   subroutine read_abundances(p, q)
     real(dp), intent(in) :: p(:)
     real(dp), intent(out),allocatable :: q(:,:)
     
     integer :: k, nqr,npz,ndat, io, unit, rand,n, i, test

     real(dp), allocatable :: qdat(:,:), pdat(:)
     character(len=100) :: formt
     
     !nqt = size(q, 2)
     npz = size(p)
     allocate(mmw_dry(npz), cp_dry(npz), q_orig(npz,nqt), q(npz,nqt))
     open(newunit=unit, file=abundance_file, iostat=io)

     read(unit,*) ndat, test
! Check that nqt == test
     write(*,*) ndat, test, nqt
     if (test .ne. nqt) then
        write(*,*) 'WARNING, abundance data length .ne. nqt'
        stop
     endif

     allocate(qdat(ndat,nqt), pdat(ndat))

     do k=1,ndat
        read(unit,*) pdat(k), qdat(k,:)
     enddo

     close(unit)

     
! Now interpolate all the qs in log space
     do k=1,npz
           i = minloc(abs(p(k) - pdat), dim=1)
           if (p(k) .lt. pdat(i) ) i=i-1
           i=max(min(i,ndat-1),1)
           do n=1,nqt
              ! Find nearest index in data
              !if (n .eq. 1) write(*,*), k, p(k), pdat(i), pdat(i+1)
              call linear_log_interp(p(k), pdat(i), pdat(i+1), qdat(i,n), qdat(i+1,n), q(k,n))
           enddo
     enddo


! Save the original abundance of species
     q_orig = q(:,:)
! Find the mmw and heat capacity of the dry component
     write(*,*) 'mmw_dry'
     do k=1,npz
        call get_dry_mmw(q(k,1:nqt), mmw_dry(k), th_gases)
        write(*,*) k, mmw_dry(k)
        call get_dry_cp(q(k,1:nqt), cp_dry(k), th_gases)
     enddo
     
   end subroutine read_abundances

  subroutine get_dry_mmw(q, mmw, gases)
    real(dp), intent(in) :: q(nqt)
    real(dp), intent(out) :: mmw
    type(gas),intent(in)  :: gases(nqt)

    integer :: i
    real(dp) :: mmw_r, weight

    mmw_r = 0.0_dp
    weight = 0.0_dp
    do i=1,nqt
       if (gases(i)%soc_index .ne. 1) then
          mmw_r = mmw_r + q(i)/gases(i)%mmw
          weight = weight + q(i)
       endif
    enddo
    mmw = 1._dp/mmw_r/weight
    

  end subroutine get_dry_mmw

  subroutine get_dry_cp(q, cp, gases)
    real(dp), intent(in) :: q(nqt)
    real(dp), intent(out) :: cp
    type(gas),intent(in)  :: gases(nqt)

    integer :: i
    real(dp) ::weight

    cp =  0.0_dp
    do i=1,nqt
       if (gases(i)%soc_index .ne. 1) then
          cp = cp + q(i)*gases(i)%cp
          weight = weight + q(i)
       endif
    enddo
    cp = cp/weight

  end subroutine get_dry_cp
  
  subroutine get_mmw(q, mmw)
    real(dp), intent(in) :: q(nqt)
    real(dp), intent(out) :: mmw
    
    integer :: i
    real(dp) ::  mmw_r

    mmw_r = 0.0_dp
    
    do i=1,nqt
       mmw_r = mmw_r + q(i)/th_gases(i)%mmw
    enddo
    mmw = 1._dp/mmw_r
    
  end subroutine get_mmw

  subroutine get_cp(q, cp)
    real(dp), intent(in) :: q(nqt)
    real(dp), intent(out) :: cp
    !type(gas), intent(in) :: gases(nqt)

    integer :: i

    cp = 0.0_dp
    do i=1,nqt
       cp = cp + q(i)*th_gases(i)%cp
    enddo
  end subroutine get_cp
  
end module atmosphere
