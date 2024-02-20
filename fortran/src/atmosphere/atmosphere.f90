module atmosphere

  use params, only: soc_index_file, therm_index_file, dp, p_sc, abundance_file
  use utils, only: linear_log_interp
  use phys
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
    real(dp), intent(inout):: q(:,:)
    
    integer :: f_unit, ierr, i

! Open file
    open(newunit=f_unit, iostat=ierr, action='read', file=soc_index_file)
    if (ierr .ne. 0) then
       write(*,*) 'Could not open ', soc_index_file
       stop
    endif

    read(f_unit, *) nqr
    allocate(soc_indices(nqr), rad_gases(nqr))
    
    do i=1,nqr
       read(f_unit,*) soc_indices(i)
    enddo
    close(f_unit)

    ! Open file
    open(newunit=f_unit, iostat=ierr, action='read', file=therm_index_file)
    if (ierr .ne. 0) then
       write(*,*) 'Could not open ', Therm_index_file
       stop
    endif

    read(f_unit, *) nqt
    allocate(th_indices(nqt), th_gases(nqt))

    do i=1,nqt
       read(f_unit,*) soc_indices(i)
    enddo
    close(f_unit)

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

       case default
          write(*,*) 'ERROR: incorrect species identifier'
          stop
       end select
    enddo

    call read_abundances(p,q)
  end subroutine init_atmos
  
   subroutine read_abundances(p, q)
     real(dp), intent(in) :: p(:)
     real(dp), intent(out) :: q(:,:)
     
     integer :: k, nqr,npz,ndat, io, unit, rand,n, i

     real(dp), allocatable :: qdat(:,:), pdat(:)
     
     nqr = size(q, 2)
     npz = size(p)

     allocate(mmw_dry(npz), cp_dry(npz), q_orig(npz,nqr))
     open(newunit=unit, file=abundance_file, iostat=io)

     ndat = 0
     do while (io .eq. 0)
        read(unit,*, iostat=io)
        ndat = ndat+1
     enddo

     allocate(qdat(ndat,nqr), pdat(ndat))
     rewind(unit)
     do k=1,ndat
        read(unit,*) rand,pdat(k), qdat(k,:)
     enddo

     close(unit)
     
! Now interpolate all the qs in log space
     do k=1,npz
        if (p(k) .gt. p_sc) then
           q(k,1) = 1.0_dp
           q(k,2:nqr) = 0.0_dp
        else
           i = minloc(abs(p(k) - pdat), dim=1)
           if (p(k) .lt. pdat(i) ) i=i-1
           i=max(min(i,ndat-1),1)
           do n=1,nqr
! Find nearest index in data
              call linear_log_interp(p(k), pdat(i), pdat(i+1), qdat(i,n), qdat(i+1,n), q(k,n))
           enddo
        endif
     enddo

! Save the original abundance of species
     q_orig = q(:,:)
! Find the mmw and heat capacity of the dry component
     do k=1,npz
        call get_dry_mmw(q(k,1:nqt), mmw_dry(k), th_gases)
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
  
  subroutine get_mmw(q, mmw, gases)
    real(dp), intent(in) :: q(nqt)
    real(dp), intent(out) :: mmw
    type(gas),intent(in)  :: gases(nqt)
    
    integer :: i
    real(dp) ::  mmw_r

    mmw_r = 0.0_dp
    do i=1,nqt
       mmw_r = mmw_r + q(i)/gases(i)%mmw
    enddo
    mmw = 1._dp/mmw_r
    
  end subroutine get_mmw

  subroutine get_cp(q, cp, gases)
    real(dp), intent(in) :: q(nqt)
    real(dp), intent(out) :: cp
    type(gas), intent(in) :: gases(nqt)

    integer :: i

    cp = 0.0_dp
    do i=1,nqt
       cp = cp + q(i)*gases(i)%cp
    enddo
  end subroutine get_cp
  
end module atmosphere
