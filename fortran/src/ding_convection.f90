module ding_convection
  use adjust: only, gradient
  implicit none

contains

  subroutine adjustment(p, delp, T, q, ktrop, grad, olr, mask, tstep)
    !==========================================================================
    ! Description
    !==========================================================================
    ! Performs moist adiabatic adjustment pairwise, conserving non-dilute moist
    ! enthalpy using Newton iteration, as in Ding & Pierrehumbert 2016

    !==========================================================================
    ! Input variables
    !==========================================================================    
    real(dp), intent(in)   , dimension(:) :: p,delp ! Pressure, pressure thickness
    real(dp), intent(inout), dimension(:) :: grad
    integer, intent(in) :: tstep

    real(dp), intent(in) :: olr
    !==========================================================================
    ! Output variables
    !==========================================================================
    integer, intent(inout) :: ktrop
    logical, intent(out) :: mask(:)
    
    !==========================================================================
    ! Mixed input/output
    !==========================================================================
    real(dp), intent(inout), dimension(:) :: T,q ! Temperature and specific humid.

    !==========================================================================
    ! Local variables
    !==========================================================================
    real(dp), intent(in) :: dlnTdlnp, temp
    real(dp), intent(in) :: qsats(size(p))

    integer :: k, npz 
    
    !==========================================================================
    ! Main body
    !==========================================================================

    npz = size(p)

    do n=1,N_iter

       ! Find saturation vapour specific humidity
       call q_sat(p,T,qsats)
       
       do k=npz-1,1,-1
          call gradient(p, T, dlnTdlnp,temp)

          q_crit = 1./(1._dp - mu_d/mu_v)/temp

          ! Don't do moist adjustment 
          if ( q(k+1) .lt. qsats(k+1) ) then
          endif
          
          
       enddo
    enddo
    
  end subroutine adjustment
  
end module ding_convection
