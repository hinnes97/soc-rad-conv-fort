module supercrit_adjust
  use params,     only: p_sc, aqua_path,dp, accelerate, Finc, inhibited
  use aqua_eos,   only: dim_output, load_table_pt, interpolate_aqua_pt
  use adjust_mod, only: gradient
  use condense, only: q_sat
  use atmosphere, only : cp_dry, mmw_dry, q_orig, th_gases, get_mmw, get_cp
  use phys, only : H2O, Rstar
  implicit none
  
contains

  subroutine aqua_init
    call load_table_pt(aqua_path)
  end subroutine aqua_init
  
  subroutine adjustment(p, delp, T, q, ktrop, grad, olr, mask, tstep)
    
    !=========================================================================
    ! Dummy variables
    !=========================================================================
    
    real(dp), intent(in),    dimension(:) :: p, delp    ! Pressure and pressure thickness
    real(dp), intent(inout), dimension(:) :: grad, T   ! Ad. gradient, temp
    real(dp), intent(inout), dimension(:,:) :: q
    integer,  intent(inout)               :: ktrop      ! Index of tropopause (~cold trap)
    logical,  intent(out),   dimension(:) :: mask       ! Mask of whether convection occurring
    real(dp), intent(in) :: olr
    integer, intent(in) :: tstep
    !=========================================================================
    ! Local variables
    !=========================================================================

    integer :: npz,k,n,nqt,m
    integer, parameter :: N_iter=1000 ! Number of up-down iterations
    logical :: condition

    real(dp) :: qsats(size(p), size(q,2)), kappa, qmin, f, mmw,cp
    real(dp) :: eos_pt(dim_output), qcrit, pfact, temp
    real(dp), parameter          :: delta = 0.0001 ! Small number speeds up convergence
    !==========================================================================
    ! Main body
    !==========================================================================

    npz = size(p)
    nqt = size(q,2)
! Need to check for both dry and wet cases - how did I do this in the GCM?

    ! write(*,*) 'before  convection'
    ! do k=1,npz
    !    write(*,*) p(k), T(k), q(k,1)
    ! enddo
    ! write(*,*) '------------------------------------------------'
    mask = .false.
    do n=1,N_iter
       call q_sat(p, T, qsats)
       
       if (n.eq. 1 .and. tstep .gt. 1000 .and. accelerate) then
          f =  (Finc/olr)**(0.001_dp)
       else
          f = 1.
       endif

       do k=npz-1,max(ktrop,1), -1
          if (p(k) .gt. p_sc) then
             ! Do adjustment using the adiabat from AQUA
             eos_pt =  interpolate_aqua_pt(p(k), T(k))
             grad(k) = eos_pt(2)
             pfact =exp(grad(k)*log(p(k)/p(k+1)))
             condition = ( (T(k) - T(k+1)*pfact*(1+delta)) .lt. 0 )
          else
! Do dry/moist adjustment based on whether there is condensation (can also check composition too)
             if (q(k+1,1) .gt. qsats(k+1,1)-1.e-10) then
! Moist adjustment
                call gradient(p(k+1),T(k+1),mmw_dry(k+1), grad(k), temp)
                qcrit = 1./(1._dp - mmw_dry(k+1)/H2O%mmw) /temp
                pfact =exp(grad(k)*log(p(k)/p(k+1)))
                if (inhibited) then
                ! Do moisture inhibition
                   condition = (  (T(k) - T(k+1)*pfact*(1 + delta))*(qcrit - q(k+1,1)) .lt. 0)
                else
                  ! Regular criterion
                   condition = ( (T(k) - T(k+1)*pfact*(1+delta)) .lt. 0 )
                endif
                
             else
! Dry adjustment
                ! Calculate dry adiabatic index
                call get_mmw(q(k+1,:) , mmw)
                call get_cp(q(k+1,:), cp, th_gases)
                
                grad(k) = Rstar/mmw/cp
                pfact = exp(grad(k)*log(p(k)/p(k+1)))
                condition = ( (T(k) - T(k+1)*pfact*(1+delta)) .lt. 0 )                
             endif
          endif
          if (condition) then
             T(k+1) = (T(k)*delp(k) + T(k+1)*delp(k+1))/(delp(k+1) + pfact*delp(k))
             T(k+1) = f*T(k+1)
             T(k) = f*T(k+1)*pfact
             mask(k) = .true.
             mask(k+1) = .true.
          endif
       enddo

       ! Set to saturation
       call q_sat(p, T, qsats)
       do k=1,npz
          if (q(k,1) .gt. qsats(k,1)) then
             q(k,:)  = qsats(k,:)
          else
             q(k,:) = q_orig(k,:)
          endif
          q(k,1) = min(qsats(k,1), q_orig(k,1))
       enddo

       ! Cold trapping
       ktrop = 1
       qmin = 1000._dp
       do k = 1,npz
          if (q(k,1) .lt. qmin) then
             qmin = q(k,1)
          else if (p(k) .lt. 9*p(k)/10. .and. q(k,1) .lt. q_orig(k,1)) then
             q(1:k,1) = qmin
             ktrop = k+1
             do m=1,k
                q(m,2:) = (1. - q(m,1))/(1 - q_orig(m,1)) * q_orig(m,2:)
             enddo
          endif
       enddo
    enddo
    
    ! write(*,*) 'after convection'
    ! do k=1,npz
    !    write(*,*) p(k), T(k), q(k,1)
    ! enddo
    ! write(*,*) '------------------------------------------------' 
    
  end subroutine adjustment
end module supercrit_adjust
