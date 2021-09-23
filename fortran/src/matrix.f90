module matrix
  use params, only: dp, N_max
  use flux_mod, only: get_fluxes
  use utils, only: linear_log_interp
!  use lapack, only: dgetrf
  
  implicit none

  real(dp) :: pert_T = 0.1_dp
contains
  
  subroutine do_matrix(nf, ne, Tf, pf, Te, pe, tau_IR, tau_V, &
        mu_s, Finc, Fint, olr)

    integer, intent(in) :: nf, ne
    real(dp), intent(in) :: mu_s, Finc, Fint
    real(dp), intent(in), dimension(:) :: pf, pe, tau_IR, tau_V
    real(dp), intent(inout) :: Te(:), Tf(:)
    real(dp), intent(out) :: olr

    ! Work variables
    real(dp), dimension(ne) :: residual, del_T
    real(dp), dimension(ne,ne) :: mat
    integer :: n,i
    
    do n=1,N_max

       call calc_matrix(mat, nf, ne, Tf, pf, pe, tau_IR, tau_V, mu_s, Finc, Fint, olr, residual, Te)
       call solve_matrix(mat, del_T, residual, ne)

       Te = Te + del_T
       
       do i=1,ne
          if (Te(i) .lt.  100.0_dp) Te(i) = 100.0_dp
       end do
       
    end do
    
  end subroutine do_matrix

  subroutine calc_matrix(mat, nf, ne, Tf, pf, pe, tau_IR, tau_V, mu_s, Finc, Fint, olr, residual,Te)
    integer, intent(in) :: nf, ne
    real(dp), intent(in) :: Tf(:), pe(:), pf(:), tau_IR(:), tau_V(:)
    real(dp), intent(out) :: mat(:,:)
    real(dp), intent(in) :: mu_s, Finc, Fint
    real(dp), intent(out) :: olr
    real(dp), intent(inout) :: residual(:), Te(:)

    real(dp), dimension(ne) :: Tpert, flux, flux_pert
    integer :: i

    call get_fluxes(nf, ne, Tf, pf, Te, pe, tau_IR, tau_V, &
         flux, mu_s, Finc, Fint, olr)

    residual = flux - Fint
    
    do i=1,ne
       Tpert = Te
       Tpert(i) = Tpert(i) + pert_T
       call get_fluxes(nf, ne, Tpert(2:), pf, Tpert, pe, tau_IR, tau_V, &
            flux_pert, mu_s, Finc, Fint, olr)

       mat(:,i) = (flux_pert - flux)/pert_T
    end do
    
  end subroutine calc_matrix

  subroutine solve_matrix(mat, del_T, residual, ne)
    integer, intent(in) :: ne
    real(dp), intent(in) :: mat(ne,ne), residual(ne)
    real(dp), intent(out) :: del_T(ne)

    ! Work variables
    real(dp), dimension(ne,ne) :: factors,newmat
    integer :: ipiv(ne), info, lwork,nrhs
    real(dp), dimension(:), allocatable :: work
    real(dp),dimension(ne) ::  R, C
    character(len=1) :: equed
    real(dp) :: rcond
    real(dp), dimension(1) :: ferr, berr
    integer, dimension(ne) :: iwork
    real(dp), dimension(ne,1) :: input
    real(dp), dimension(ne,1) :: output
    nrhs = 1
    
    lwork = ne*4
    allocate(work(lwork))

    newmat = mat
    input(:,nrhs) = -residual
    Call dgesvx('E', 'N', ne, 1, newmat, ne, factors, ne, &
         ipiv, equed, r, c, input, &
         ne, output, ne, rcond, ferr, berr, work, iwork, &
         info)


    del_T = output(:,nrhs)

  end subroutine solve_matrix
end module matrix
