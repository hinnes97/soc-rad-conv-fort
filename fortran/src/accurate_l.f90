module accurate_l

  use params, only : dp
  implicit none


  ! Global data

  real(dp) :: pc = 22.064e6
  real(dp) :: Tc = 647.096

  real(dp) ::  a(6) = (/-7.85951783,1.84408259,-11.7866497,22.6807411,-15.9618719,1.80122502/)

contains
  
  real function log_ps_pc(T)
    
    real(dp), intent(in) :: T

    real(dp) :: th
    real(dp) :: temp
    real(dp) :: powers(6) = (/1., 1.5, 3., 3.5, 4., 7.5/)
    integer :: i

    th = 1. - T/Tc

    temp = 0.0
    
    do i=1,6
       temp = temp + a(i)*th**powers(i)
    enddo

    log_ps_pc  = temp * Tc/T
  end function log_ps_pc

  real function p_sat(T)
    real(dp), intent(in) :: T

    p_sat = pc * exp(log_ps_pc(T))
  end function p_sat

  real function dlogp_dlogt(T)
    real(dp), intent(in) :: T

    integer :: i

    real(dp) :: powers(6)  = (/0.,0.5,2.,2.5,3.,6.5/)
    real(dp) :: th
    real(dp) :: logps_pc
    real(dp) :: temp
    
    temp = 0.0

    logps_pc = log_ps_pc(T)
    th = 1 - T/Tc
    
    do i=1,6
       temp = temp + a(i)*(powers(i) + 1.)*th**powers(i)
    enddo

    dlogp_dlogt = -1 * (temp + log_ps_pc(T))
    
    
  end function dlogp_dlogt

  real function dpdt(T)
    real(dp), intent(in) :: T

    dpdt = p_sat(T)/T * dlogp_dlogt(T)
    
  end function dpdt
  
  real function alpha(T)

    real(dp), intent(in) :: T
    
    real(dp) :: d(6) = (/-1135.905627715,-5.65134998e-8,2690.66631,127.287297,&
         -135.003439,0.981825814/)
    real(dp) :: powers(6) = (/0.0, -19., 1., 4.5, 5., 54.5/)
    
    real(dp) :: th
    real(dp) :: temp
    
    integer :: i

    th = T/Tc
    temp = 0.0
    do i=1,6
       temp = temp + d(i)*th**powers(i)
    enddo

    alpha = temp
    
  end function alpha

  real function rho_liq(T)
    real(dp), intent(in) :: T

    real(dp) :: b(7)  = (/1., 1.99274064,1.09965342,-0.510839303,-1.75493479,&
         -45.5170352,-6.74694450e5/)
    real(dp) :: powers(7) = (/0., 1./3.,2./3.,5./3.,16./3., 43./3.,110./3. /)

    real(dp) :: rhoc = 322.
    real(dp) :: rho_rhoc = 0.0

    integer :: i
    real(dp) :: th

    th = 1-T/Tc

    do i=1,7
       rho_rhoc = rho_rhoc + b(i) * th**powers(i)
    enddo
    
    rho_liq = rhoc*rho_rhoc
    
  end function rho_liq
  
  real function rho_vap(T)
    real(dp), intent(in) :: T

    real(dp) :: th
    real(dp) :: rhoc = 322.
    real(dp) :: c(6) = (/-2.03150240,-2.68302940, -5.38626492,-17.2991605,&
         -44.7586581,-63.9201063/)
    real(dp) :: powers(6) = (/ 2./6., 4./6., 8./6., 18./6., 37./6., 71./6. /)

    real(dp) :: logrho_rhoc

    integer :: i
    
    th = 1.-T/Tc
    logrho_rhoc = 0.0
    
    do i=1,6
       logrho_rhoc = logrho_rhoc + c(i)*(th**powers(i))
    enddo

    rho_vap = rhoc*exp(logrho_rhoc)

    
  end function rho_vap

  
  real function h_liq(T)
    real(dp), intent(in) :: T

    h_liq = (1.e3*T/rho_liq(T) * dpdt(T)/1.e6)*1000.
    
  end function h_liq

  real function h_vap(T)
    real(dp), intent(in) :: T
    
    h_vap = (1.e3*T/rho_vap(T) * dpdt(T)/1.e6)*1000.
    
  end function h_vap
  
  real function L_calc(T)
    real(dp), intent(in) :: T

    L_calc = h_vap(T) - h_liq(T)
  end function L_calc
  
       
  
end module accurate_l
