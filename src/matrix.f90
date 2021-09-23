module matrix
  use params, only: dp, N_max
  use flux_mod, only: get_fluxes
  use utils, only: linear_log_interp
!  use lapack, only: dgetrf
  
  implicit none

  real(dp) :: pert_T = 0.1_dp
contains
  
  subroutine do_matrix(nf, ne, Tf, pf, pe, tau_IR, tau_V, &
        mu_s, Finc, Fint, olr)

    integer, intent(in) :: nf, ne
    real(dp), intent(in) :: mu_s, Finc, Fint
    real(dp), intent(in) :: pf(:), pe(:), tau_IR(:), tau_V(:)
    real(dp), intent(inout) :: Tf(:)
    real(dp), intent(out) :: olr

    ! Work variables
    real(dp), dimension(ne) :: T, residual, del_T, Te
    real(dp), dimension(ne,ne) :: mat
    real(dp) :: Te2
    integer :: n,j,i
    
    ! Set up T
    call linear_log_interp(pe(2), pf(1), pf(2), Tf(1), Tf(2), Te2)
    T(1) = Tf(1) + (pe(1) - pe(2))/(pf(1) - pe(2)) * (Tf(1) - Te2)
    T(2:) = Tf
    
    do n=1,N_max

       call calc_matrix(T, mat, nf, ne, Tf, pf, pe, tau_IR, tau_V, mu_s, Finc, Fint, olr, residual, Te)
       call solve_matrix(mat, del_T, residual, ne)
       
       !do j=1,ne
       !   if(del_T(j)>20.0_dp) del_T(j) = 20.0_dp
       !   if(del_T(j)<-20.0_dp) del_T(j) = -20.0_dp
       !enddo
!       do j=1,ne
!          write(*,*) del_T(j)
!       end do
       
       !write(*,*) maxval(del_T), minval(del_T)
       T = T + del_T
       
       do j=1,ne
          if(T(j)<100.0_dp) T(j) = 100.0_dp
       end do
       
    end do
    Tf = T(2:)
    write(*,*) '----------------------------------'
    write(*,*) matmul(mat, -del_T)
    
    do j=1,ne
       write(*,*) mat(1,j)
    enddo
    write(*,*) '-.-.-.--.-.-.-.-.-.-.-.-.-.-.-.-.-.-.'
    do j=1,ne
       write(*,*) del_T(j)
    enddo

    open(10, file='mat.out')
    do i = 1, ne
        write(10,'((ES10.3))') (mat(i,j), j = 1, ne)
     end do
    close(unit=10)

    writE(*,*) 'TE', Te
    open(11, file='Te.out')
    do i = 1, ne
        write(11,'((ES10.3))') Te(i)
     end do
    close(unit=11)

  end subroutine do_matrix

  subroutine calc_matrix(T, mat, nf, ne, Tf, pf, pe, tau_IR, tau_V, mu_s, Finc, Fint, olr, residual,Te)
    integer, intent(in) :: nf, ne
    real(dp), intent(in) :: T(:), Tf(:), pe(:), pf(:), tau_IR(:), tau_V(:)
    real(dp), intent(out) :: mat(:,:)
    real(dp), intent(in) :: mu_s, Finc, Fint
    real(dp), intent(out) :: olr
    real(dp), intent(inout) :: residual(:), Te(:)

    real(dp), dimension(ne) :: Tpert, flux, flux_pert
    integer :: i, j, fu

    !fu = 1
    ! Temporary, read in Te to compare python and fortran
    !open(fu, file='Te.txt')
    !read(fu,*) (Te(i), i=1,ne)

    do j = 2, nf
       call linear_log_interp(pe(j), pf(j-1), pf(j), Tf(j-1), Tf(j), Te(j))
    end do
    
    ! Extrapolate to find Te at uppermost and lowest levels
    Te(1) = T(1)
    Te(ne) = Tf(nf) + (Tf(nf) - Tf(nf-1))/(pf(nf) - pf(nf-1)) *(pe(ne) - pf(nf))
    !Te(ne) = Tf(nf) + (pe(ne) - pe(nf))/(pf(nf) - pe(nf)) * (Tf(nf) - Te(nf))
    
    call get_fluxes(nf, ne, Tf, pf, Te, pe, tau_IR, tau_V, &
         flux, mu_s, Finc, Fint, olr)

    residual = flux - Fint
    
    !do i=1,ne
       !write(*,*) maxval(abs(residual)), minval(abs(residual))
!       write(*,*del_T(i)
       !write(*,*) residual(i)
    !end do
    
    do i=1,ne
       Tpert = T
       Tpert(i) = Tpert(i) + pert_T
       do j = 2, nf
          call linear_log_interp(pe(j), pf(j-1), pf(j), Tpert(j), Tpert(j+1), Te(j))
       end do
    
       ! Extrapolate to find Te at uppermost and lowest levels
       Te(1) = Tpert(1)
       !write(*,*) (T(j), Tpert(j), Te(j), j=1,ne)
       !Te(ne) = T(nf) + (pe(ne) - pe(nf))/(pf(nf) - pe(nf)) * (Tf(nf) - Te(nf))
       Te(ne) = Tpert(ne) + (Tpert(ne) - Tpert(nf))/(pf(nf) - pf(nf-1)) *(pe(ne) - pf(nf))

       call get_fluxes(nf, ne, Tpert(2:), pf, Te, pe, tau_IR, tau_V, &
            flux_pert, mu_s, Finc, Fint, olr)

       mat(:,i) = (flux_pert - flux)/pert_T
    end do
    !write(*,*) mat(:,1)
    
  end subroutine calc_matrix

  subroutine solve_matrix(mat, del_T, residual, ne)
    integer, intent(in) :: ne
    real(dp), intent(in) :: mat(ne,ne), residual(ne)
    real(dp), intent(out) :: del_T(ne)


!!$
!!$    Real (dp) :: rcond
!!$    Integer :: i, ifail, info, lda, ldaf, ldb, ldx, nrhs
!!$    Character (1) :: equed
!!$!     .. Local Arrays ..
!!$    Real (dp), Allocatable :: a(:, :), af(:, :), b(:, :), berr(:), &
!!$         c(:), ferr(:), r(:), work(:), x(:, :)
!!$    Integer, Allocatable :: ipiv(:), iwork(:)

    ! Work variables
    real(dp), dimension(ne,ne) :: factors,newmat, test
    integer :: ipiv(ne), info, lwork,nrhs
    real(dp), dimension(:), allocatable :: work
    real(dp),dimension(ne) ::  R, C
    integer :: i,j
    character(len=16) :: fact = 'Equilibration'
    character(len=16) :: trans ='No transpose'
    character(len=1) :: equed
    real(dp) :: rcond
    real(dp), dimension(1) :: ferr, berr
    integer, dimension(ne) :: iwork
    real(dp), dimension(ne,1) :: input, b
    real(dp), dimension(ne,1) :: output
    !real(dp), dimension(ne) :: output
    nrhs=1
    !nrhs = 1
    !lda = ne
    !ldaf = ne
    !ldb = ne
    !ldx = ne

    !Allocate (a(lda,ne), af(ldaf,ne), b(ldb,nrhs), berr(nrhs), c(ne), &
    !     ferr(nrhs), r(ne), work(4*ne), x(ldx,nrhs), ipiv(ne), iwork(ne))

    !a = mat
    !b(:,nrhs) = -residual
    !input(:,1) = -residual
    
    lwork = ne*4
    allocate(work(lwork))
    ! Store mat in mat_inv so lapack doesn't overwrite
    !factors = mat

    !del_T = -residual
    newmat = mat
    input(:,nrhs) = -residual
    Call dgesvx('E', 'N', ne, 1, newmat, ne, factors, ne, &
         ipiv, equed, r, c, input, &
         ne, output, ne, rcond, ferr, berr, work, iwork, &
         info)

    !b(:,1) = -residual
    !output(:,nrhs) = -residual
    !call dgesv(ne, nrhs,  newmat, ne, ipiv, b, ne, info)

    del_T = output(:,nrhs)
    write(*,*) matmul(mat, -output)
    write(*,*) '------------------------'
    write(*,*) matmul(mat, -del_T)
    !do i =1,ne
    !   write(*,*)  ouTput(I,1) - b(i,1)
    !end do
    !write(*,*) 1./rcond
    
    !output(:,1) = -residual
    !newmat = mat

    !write(*,*) output
    !call dgetrf(ne,ne,mat_inv,ne,ipiv,info)    
    !if (info .ne. 0) stop "Matrix is singular"

    !call dgetri(ne,mat_inv,ne,ipiv,work,lwork,info)
    !if (info .ne. 0) stop "Matrix inversion failed"
!!$
!!$    write(*,*) rcond, ferr, berr, info
!!$    open(10, file='mat.out')
!!$    do i = 1, ne
!!$        write(10,'((ES10.3))') (mat(i,j), j = 1, ne)
!!$     end do
!!$    close(unit=10)
!!$
!!$    open(11, file='del_T.out')
!!$    do i=1,ne
!!$       write(11,'((ES10.3))') del_T(i)
!!$    end do
!!$    close(unit=11)
!!$
!!$    open(12,file='res.out')
!!$    do i=1,ne
!!$       write(12, '((ES10.3))') residual(i)
!!$    end do
!!$
!!$    open(13,file='test.out')
!!$    test = matmul(-mat,del_T)
!!$    do i=1,ne
!!$       write(13, '((ES10.3))') test(i)
!!$    end do
!!$
!!$#endif    
    !write(*,*) maxval(del_T), minval(del_T)
    ! Multiply mat_inv by -residual for solution
    !del_T = matmul(mat_inv, -residual)
  end subroutine solve_matrix
end module matrix
