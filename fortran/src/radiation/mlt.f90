module mlt_mod

  use params, only : grav, q0, inhibited
  use adjust_mod, only: gradient
  use condense, only: q_sat
  use params, only: dp
  use phys, only : Rstar
  use atmosphere, only: cp_dry, mmw_dry
  implicit none
contains
  subroutine mlt_flux(p, T, q, flux, Te, pe, ktrop)
    integer, intent(in) :: ktrop
    real(dp), intent(in) :: p(:), T(:), q(:,:), Te(:), pe(:)
    real(dp), intent(out) :: flux(:)

    integer :: npz, k
    real(dp) :: qsats(size(q,1), size(q,2))
    real(dp) :: rho, cp(size(p))
    real(dp) :: grad, temp, pmid, Tmid, qmid, mmwdry_mid, R, grad_real, qcrit, cp_mid
    real(dp) :: gdiff, factor, scale_H, conv_V

    flux = 0.0
    npz = size(p)
    
    call q_sat(p, T, qsats)

    do k=npz-1,max(ktrop,1), -1
       grad_real = (log(T(k+1)) - log(T(k)))/(log(p(k+1)) - log(p(k)))
       pmid = (p(k+1) + p(k))*0.5
       Tmid = (T(k+1) + T(k))*0.5
       qmid = q(k+1,1)*0.5 + q(k,1)*0.5
       mmwdry_mid  = (mmw_dry(k) + mmw_dry(k+1))*0.5
       cp_mid = cp_dry(k+1)*0.5 + cp_dry(k)*0.5
       
       R = Rstar*((1-qmid)/mmwdry_mid + qmid/18.0)
       rho = pmid/Tmid/R
       
       scale_H = R*Tmid/grav
       
       
       if (q(k+1,1) .gt. qsats(k+1,1) - 1.e-10 .and. q(k+1,1) .lt. q0 ) then
          call gradient(pmid, Tmid, cp_mid, mmwdry_mid, grad, temp)
          qcrit = 1./(1._dp - mmw_dry(k)/18.0) /temp
          if (inhibited ) then
             factor = (qcrit - qmid)
          else
             factor= 1.0
          endif
          gdiff = max((grad_real - grad), 0.0)
       else if (q(k+1,1) .lt. qsats(k+1,1)- 1.e-10) then
          grad = 2./7.
          gdiff = max(grad_real - grad, 0.0)
       endif

       conv_V = sqrt(gdiff*grav*scale_H)
       flux(k) = 1./2. * rho * cp_mid * Tmid *gdiff * conv_V

       !write(*,*) 'mlt theory', k, rho, cp_mid, grav, Tmid, R, gdiff
       
       !flux(k) = 1./2.**2.5 * rho**1.5 * cp_dry(k)*R*Tmid**2/pmid**0.5*gdiff**1.5
    enddo

    flux(1)  = 0.0

    grad_real  = (log(Te(npz+1)) - log(T(npz)))/(log(pe(npz+1)) - log(p(npz)))
    if (q(npz,1) .gt. qsats(npz,1) - 1.e-10 .and. q(npz,1) .lt. q0 ) then
       call gradient(p(npz), T(npz), cp_dry(npz), mmw_dry(npz), grad, temp)
       qcrit = 1./(1._dp - mmw_dry(npz)/18.0) /temp
       if (inhibited) then
          factor = qcrit-qmid
       else
          factor = 1.0
       endif
       gdiff = max((grad_real - grad)*factor, 0.0)
    else if (q(npz,1) .lt. qsats(npz,1)- 1.e-10) then
       grad = 2./7.
       gdiff = max(grad_real - grad, 0.0)
    endif
    R = Rstar*((1-q(npz,1))/mmw_dry(npz) + q(npz,1)/18.0)
    scale_H = R*Te(npz+1)/grav
    rho = pe(npz+1)/Te(npz+1)/R
    conv_V = sqrt(gdiff*grav*scale_H)
    
    flux(npz+1) = 1./2. * rho * cp_mid * Tmid * gdiff * conv_V
    !flux(npz+1) = 1./2.**2.5 * rho**1.5*cp_dry(npz)*R*Te(npz+1)**2/pe(npz+1)**0.5*gdiff**1.5
  end subroutine mlt_flux
end module mlt_mod
