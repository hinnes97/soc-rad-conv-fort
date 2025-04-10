module condense

  use params, only : rdgas, q0, dp, nf, moisture_scheme
  use phys, only: H2O, Rstar
  use accurate_l, only : p_sat, L_calc, log_ps_pc, pc, dlogp_dlogt
  use tables, only : phase_grad, lheat, satur, find_var_lin, find_var_loglin
  use atmosphere, only: mmw_dry, q_orig, nqt, nqr
  use q_sat_test, only: q_sat, dew_point_t, sat_vp, newton_iteration
  use utils, only: write_reala
  
  implicit none

  
  real(dp) :: rvgas = 461.52
  real(dp) :: L_sub = 2.84e6 ! Defined at triple point
  !real(dp) :: eps = mu_v/mu_d ! Ratio of water to hydrogen molecular weight

  real(dp) :: p_nought
  integer :: a,b
  
contains

  
  

  
  subroutine rain_out(p,T,q, qsat)
    real(dp), intent(in) :: p, T
    real(dp), intent(inout) :: q
    real(dp), intent(out) :: qsat
    real(dp) :: eps, sat_p

    eps = rdgas/rvgas

    ! Find saturation vapour pressure
    call sat_vp(T, sat_p)
    sat_p = p_sat(T)
    qsat = (eps*sat_p/p)/(1. + (eps - 1.)* sat_p/p )
    
    !where (q .gt. qsat)
    !   q = qsat
    !elsewhere
    !   q = q0
    !endwhere
    if (q .gt. qsat) then
       q = qsat
    else
       q = q0
    endif
    
  end subroutine rain_out


  subroutine cold_trap(q, ktrop, p, T)
    real(dp), intent(inout) :: q(:,:), T(:)
    real(dp), intent(in) :: p(:)
    integer, intent(inout) :: ktrop

    integer ::  j
    real(dp) q_min
    logical :: switch 

 !   write(*,*) 'BEGINNING OF COLD TRAP'
 !   switch = .false.
 !   call write_reala(q(:,1))
    
    q_min = 10000.

    if (moisture_scheme=='deep') q_min = q0
    ktrop=1
    do j=nf,1,-1

       if (moisture_scheme=='surface') then
          
          if (q(j,1) .gt. 1) then
             q(j,1) = 1.
!             write(*,*) 'cold trap dew point T'
             call dew_point_T(p(j), T(j))
             cycle
          endif


       
          if (q(j,1) .lt. q_min) then
             q_min = q(j,1)
       
          else if (p(j) .lt. 9.*p(nf)/10.) then
             ! Extra condition stops cold trap being too low down
             ! Liable to crashing if set too near the surface
             q(1:j,1) = q_min
             ktrop = j+1
             exit
          endif

       else if (moisture_scheme=='deep') then
!          write(*,*) 'cold trapping', j, q(j,1), q_min
          if (q(j,1) .lt. q_min) then
             q_min = q(j,1)
!             write(*,*) 'changing q_min', q_min, q(j,1)
          else if (p(j) .lt. 9*p(nf)/10 .and. q_min .lt. q0-1.e-20)  then 
             q(1:j,1) = q_min
             ktrop = j+1
             exit
!             write(*,*) 'trap set', ktrop, q(j,1), q(1,1)
          endif
          
       endif
    enddo
 
       
    ! i = minloc(q(1:ktrop))
    
    ! write(*,*) 'MINLOC', i(1), q(i(1))
    ! do j=1, i(1)
    !    q(j) = q(i(1))
    ! enddo
    
  end subroutine cold_trap

    subroutine set_q(p,T, q,ktrop, tstep)
    !==========================================================================
    ! Description
    !==========================================================================    
    ! Sets water vapour mixing ratio to saturation value, and then calculates
    ! the cold trap
    !
    ! moisture_scheme can currently have two values:
    ! - "surface" - assumes surface water ocean -- atmosphere must be saturated
    !               with water vapour and temperature limited by temp when q=1.
    !
    ! - "deep" - deep water content of q0, q is min of saturation or q0
    ! 
    !==========================================================================
    
    !==========================================================================
    ! Input variables
    !==========================================================================    
      real(dp), intent(in)   , dimension(:) :: p! Pressure, pressure thickness
      
    !==========================================================================
    ! Output variables
    !==========================================================================

    integer, intent(out) :: ktrop ! Tropopause defined by cold trap, don't do
                                  ! Moist convective adjustment above this index

    integer, intent(in) :: tstep
    !==========================================================================
    ! Mixed input/output
    !==========================================================================
    real(dp), intent(inout), dimension(:) :: T! Temperature 
    real(dp), intent(inout) :: q(:,:)

    integer :: k, npz, m
    real(dp) :: qsats(size(q,1), size(q,2))
    real(dp) :: qmin
    !==========================================================================
    ! Main body
    !==========================================================================

    ! Set to saturation

    call q_sat(p, T, q)
    call cold_trap(q, ktrop, p, T)
    
    
    ! Correct if above deep value
    if (moisture_scheme =='deep') then
       npz = size(q,1)
       do k=1,npz
          q(k,1) = min(q(k,1), q0)
       enddo
       ! Don't do moist convection if entire atmosphere is at q0
       if (q(ktrop,1) .gt. q0) ktrop = npz
       do m=1,npz
          q(m,2:) = (1. - q(m,1))/(1 - q_orig(m,1)) * q_orig(m,2:)
       enddo
       
    endif
    
    if (moisture_scheme == 'supercrit') then
       ! Set to saturation
       call q_sat(p, T, qsats)
       do k=1,npz
          if (qsats(k,1) .gt. 10.0_dp) then
             ! This happens above critical point or psat>p(k) - set to 1
             q(k,1) = 1.0_dp
          else
             q(k,:)  = qsats(k,:)
             !if (qsats(k,1) .lt. 1.e-10) write(*,*) 'q is 0', k
          endif
          !q(k,1) = min(qsats(k,1), q_orig(k,1))
       enddo

       ! Cold trapping
       ktrop = 1
       qmin = 1000._dp
       do k = npz,1,-1
          if (q(k,1) .lt. qmin*(1.+1.e-4)) then
             
             qmin = q(k,1)
          else if (p(k) .lt. 9*p(npz)/10) then
             q(1:k,1) = qmin
             ktrop = k+1
             do m=1,k
                q(m,2:) = (1. - q(m,1))/(1 - q_orig(m,1)) * q_orig(m,2:)
             enddo
             exit
          endif
       enddo

    endif

    if (mod(tstep,100) .eq. 0) then
       do k=1,npz
          write(*,*) k, T(k), q(k,1), qsats(k,1)
       enddo
    endif
    
   end subroutine set_q


end module condense
