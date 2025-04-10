module jacobian
  use params, only: dp, grav,ne,nf, del_time, accelerate
  use flux_mod, only : get_fluxes
  use utils, only: linear_log_interp, interp_to_edges
  use mlt_mod, only : mlt_flux
  implicit none

  real(dp) :: pert_T = 0.01_dp
  !real(dp), allocatable :: T_old(:)
  !real(dp), allocatable :: mat(:,:)
  !real(dp), allocatable :: flux_old(:)
contains
  
  subroutine calc_matrix(Tf, pf, pe, delp,  mu_s, Finc, Fint, olr, flux,Te,q, Ts, hrt, tstep, tconsts, mat, ktrop)
    integer, intent(in) :: tstep, ktrop
    real(dp), intent(in) ::  pe(:), pf(:), delp(:)
  !  real(dp), intent(out) :: mat(:,:)
    real(dp), intent(in) :: mu_s, Finc, Fint
    real(dp), intent(out) :: olr
    real(dp), intent(inout) :: flux(:), Te(:), TF(:)
    real(dp), intent(inout) :: q(:,:)
    real(dp), intent(inout) :: Ts
    real(dp), intent(inout) :: hrt(:)
    real(dp), intent(out) :: tconsts(:)
    real(dp), intent(out) :: mat(:,:)
    
    real(dp), dimension(ne) ::  flux_pert, qpert, qsattemp, fup,fdn, s_up, Te_pert, test(nf)
    real(dp) :: s_dn(ne), hrt0, hrt1, delT(nf), hrt_test(nf), olr_tmp, Tpert(nf), mat_old(nf,nf), &
         T_older(nf)
    real(dp) :: F_mlt(nf+1), F_mlt_pert(nf+1), factor
    integer :: i,j
    
     call interp_to_edges(pf, pe, Tf, Te)
     call get_fluxes(nf, ne, Tf, pf, Te, pe,  &
         delp, flux, mu_s, Finc, Fint, olr,q, Ts, fup,fdn, s_dn, s_up)

     flux(ne) = 0.0
     call mlt_flux(pf, Tf, q, F_mlt, Te, pe, ktrop)
     !F_mlt = 0.0
     F_mlt(ne) = 0.0
     
     F_mlt = F_mlt * min(1.e-6 * 2**real(tstep, dp), 1._dp)
     write(*,*) 'MLT FLUX ------------------------'
     do i=1,nf
        write(*,*) i, F_mlt(i), F_mlt(i+1) - F_mlt(i), flux(i+1) - flux(i)
     enddo
  
     do i=1,nf
        hrt(i) = flux(i+1) - flux(i) + F_mlt(i+1) - F_mlt(i)
        Tpert = Tf
        Tpert(i) = Tpert(i) + pert_T*Tf(i)

        call interp_to_edges(pf, pe, Tpert, Te_pert)
        call get_fluxes(nf, ne, Tpert, pf, Te_pert, pe, &
             delp, flux_pert, mu_s, Finc, Fint, olr_tmp, q, Ts, fup,fdn, s_dn,s_up)

        call mlt_flux(pf, Tpert, q, F_mlt_pert, Te_pert, pe, ktrop)
        !F_mlt_pert = 0.0
        flux_pert(ne) = 0.0
        F_mlt_pert(ne) = 0.0
        F_mlt_pert = F_mlt_pert*min(1.e-6 * 2**real(tstep, dp), 1._dp)
        do j=1,nf
           hrt0 = (flux(j+1) - flux(j))  + (F_mlt(j+1) - F_mlt(j))!* grav/5500. /delp(j)
           hrt1 = (flux_pert(j+1) - flux_pert(j)) +(F_mlt_pert(j+1) - F_mlt_pert(j))!* grav/5500./delp(j)
           mat(j,i) = (hrt1 - hrt0)/(pert_T*Tf(i))
        enddo
        
!        hrt(i) = (flux(i+1)-flux(i))*grav/5500./delp(i)

     end do
 !    if (tstep .eq. 1) then
 !       allocate(T_old(nf), mat(nf,nf), flux_old(ne))
 !       mat = 0.0
 !       T_old = 0.0
 !       flux_old = 0.0
 !    endif
       
    
 !  if (any(abs(Tf - T_old) .gt. 10.0_dp) .and. accelerate) then
 !     call interp_to_edges(pf, pe, Tf, Te)
 !     call get_fluxes(nf, ne, Tf, pf, Te, pe,  &
 !         delp, flux, mu_s, Finc, Fint, olr,q, Ts, fup,fdn, s_dn, s_up)

 !     flux_old = flux
 !     !write(*,*) 'RECALCULATING MATRIX'
 !     T_older = T_old
 !     T_old = Tf
 !     !do i=1,nf
 !     !   call linear_log_interp(pf(i), pe(i), pe(i+1), Te(i), Te(i+1), Tf(i))
 !     !end do
 !     mat_old = mat
 !     do i=1,nf
 !        Tpert = Tf
 !        !qpert = q
       
 !        Tpert(i) = Tpert(i) + pert_T*Tf(i)
 ! !       call rain_out(pe(i), Tpert(i), qpert(i), qsattemp(i))

 !        call interp_to_edges(pf, pe, Tpert, Te_pert)
       
 ! !       call get_fluxes(nf, ne, Tf, pf, Tpert, pe, tau_IR, tau_V, &
 ! !            flux_pert, mu_s, Finc, Fint, olr, qpert,fup,fdn)
 !        call get_fluxes(nf, ne, Tpert, pf, Te_pert, pe, &
 !             delp, flux_pert, mu_s, Finc, Fint, olr_tmp, q, Ts, fup,fdn, s_dn,s_up)

 !        flux_pert(ne) = 0.0
 !        do j=1,nf
 !           hrt0 = (flux(j+1) - flux(j)) !* grav/5500. /delp(j)
 !           hrt1 = (flux_pert(j+1) - flux_pert(j)) !* grav/5500./delp(j)
 !           mat(j,i) = (hrt1 - hrt0)/(pert_T*Tf(i))
 !        enddo
        
 !        hrt(i) = (flux(i+1)-flux(i))*grav/5500./delp(i)

 !     end do

 !     write(*,*) '---------------- calculating fluxes ------------'
 !     do i=1,nf+1
 !        write(*,*) i,flux(i)
 !     enddo
 !     !    do j=1,nf
 !     !       test(i) = mat_old(i,j)*(Tf(j) - T_older(j))
 !     !    enddo
 !     !    write(*,*) test(i), hrt(i)
 !     ! enddo
 !   else
      
 !  ! call interp_to_edges(pf, pe, Tf, Te)
 !  ! call get_fluxes(nf, ne, Tf, pf, Te, pe,  &
 !  !      delp, flux, mu_s, Finc, Fint, olr,q, Ts, fup,fdn, s_dn, s_up)
     
      
 !    !flux(ne) =  Fint
 !      !flux = flux_old
 !      do i=1,nf
 !        do j=1,nf
 !           hrt(i) = hrt(i) + mat(i,j)*(Tf(j) - T_old(j))!(flux(i+1) - flux(i))/delp(i)*grav/8000.
 !        enddo
 !     enddo
 !      flux(ne) = Fint
 !      do i=nf,1,-1
 !         flux(i) = flux(i+1) + delp(i)*hrt(i)/grav*5500.
 !      enddo

      
 !      if (mod(tstep,500) .eq. 0) then
 !         write(*,*) 'new heating rate'
 !         do i=1,nf
 !            write(*,*) hrt(i)
 !         enddo
 !      endif
 !  endif
  
 !  do i=1,nf
 !     tconsts(i) = 1/abs(mat(i,i))
 !  enddo
 end subroutine calc_matrix

  subroutine solve_matrix(mat, resid, delT)
    real(dp), intent(in) :: mat(:,:)
    real(dp), intent(in) :: resid(:)
    real(dp), intent(inout) :: delT(:)

    integer :: M,N,LDA,IPIV(size(mat,1)),INFO
    real(dp), dimension(size(mat,1), size(mat,2)) :: fact

    fact = mat
    
    M = size(mat,1)
    N = size(mat,2)
    LDA = M

    call dgetrf(M,N,fact, LDA, IPIV, info)
    delT = resid
    
    call dgetrs('N',nf,1,fact,LDA,IPIV,delT,nf,info)
    
  end subroutine solve_matrix
end module jacobian
