MODULE phys
  IMPLICIT NONE
  real, parameter :: pi = 4.0 * atan(1.0)
!-------------Basic physical constants-------------------
!
  
!The following five are the most accurate 1986 values
REAL , PARAMETER :: h = 6.626075540E-34    !Planck's constant
REAL , PARAMETER :: c = 2.99792458E+08       !Speed of light in a vacuum
REAL , PARAMETER :: c2 = c**2
REAL , PARAMETER :: k = 1.38065812E-23     !Boltzman thermodynamic constant
REAL , PARAMETER :: kb = 1.36E-22          !Boltzmann constant in atm.molec^-1.cm^3.K^-1
REAL , PARAMETER :: sigma = 5.67051196E-08  !Stefan-Boltzman constant
REAL , PARAMETER :: G_cst = 6.67428E-11        !Gravitational constant (2006 measurements)

!-----------Thermodynamic constants------------
REAL , PARAMETER :: N_avogadro = 6.022136736E+23  !Avogadro's number
REAL , PARAMETER :: Rstar = 1000*k*N_avogadro   !Universal gas constant

  type gas
     ! Mean molecular weight
     real :: mmw
     ! Gas constant
     real :: Rgas
! Heat capacity @ const. pres
     ! TODO: Make this non-constant (how?)
     real :: cp
! Triple point temperature
     real :: tp_T
! Triple point pressure
     real :: tp_P
! Critical point temperature
     real :: crit_T
! Critical point pressure
     real :: crit_P
! Latent heat vaporization at triple point
     real :: L_vap_tp
! Latent heat vaporzation at boiling point
     real :: L_vap_bp
! Latent heat of fusion
     real :: L_fus
! Latent heat of sublimation
     real :: L_sub
! Gamma cp/cv
     real :: gamma
! Density of liquid at boiling point
     real :: rho_l_bp
! Density of liquid at triple point
     real :: rho_l_tp
! Density of solid phase
     real :: rho_s
! SOCRATES index
     integer :: soc_index
  end type gas

  
  type(gas) :: H2O
  type(gas) :: CH4
  type(gas) :: CO2
  type(gas) :: CO
  type(gas) :: N2
  type(gas) :: O2
  type(gas) :: H2
  type(gas) :: He
  type(gas) :: NH3
  type(gas) :: HCN
  type(gas) :: Earth_air
  type(gas) :: H2He_solar

contains
  

  subroutine init_gas_data
  ! H2O
  H2O%mmw      = 1.800000E+01
  H2O%Rgas     = Rstar/H2O%mmw
  H2O%cp       = 1.847000E+03
  H2O%tp_T     = 2.731500E+02
  H2O%tp_P     = 6.110000E+02
  H2O%crit_T   = 6.471000E+02
  H2O%crit_P   = 2.210000E+07
  H2O%L_vap_tp = 2.493000E+06
  H2O%L_vap_bp = 2.255000E+06 
  H2O%L_fus    = 3.340000E+05 
  H2O%L_sub    = 2.840000E+06
  H2O%gamma    = 1.331000E+00
  H2O%rho_l_bp = 9.584000E+02
  H2O%rho_l_tp = 9.998700E+02
  H2O%rho_s    = 9.170000E+02
  H2O%soc_index = 1
  
  ! CH4

  CH4%mmw      = 1.600000E+01
  CH4%Rgas     = Rstar/CH4%mmw
  CH4%cp       = 2.195000E+03
  CH4%tp_T     = 9.067000E+01
  CH4%tp_P     = 1.170000E+04
  CH4%crit_T   = 1.904400E+02
  CH4%crit_P   = 4.596000E+06
  CH4%L_vap_tp = 5.360000E+05
  CH4%L_vap_bp = 5.100000E+05
  CH4%L_fus    = 5.868000E+04
  CH4%L_sub    = 5.950000E+05
  CH4%gamma    = 1.305000E+00
  CH4%rho_l_bp = 4.502000E+02
  CH4%rho_l_tp = 0
  CH4%rho_s    = 5.093000E+02
  CH4%soc_index = 6
  
    ! CO2
  
  CO2%mmw      = 1.800000E+01
  CO2%Rgas     = Rstar/CO%mmw
  CO2%cp       = 1.847000E+03
  CO2%tp_T     = 2.165400E+02
  CO2%tp_P     = 5.185000E+05
  CO2%crit_T   = 3.042000E+02
  CO2%crit_P   = 7.382500E+06
  CO2%L_vap_tp = 3.970000E+05
  CO2%L_vap_bp = 0
  CO2%L_fus    = 1.960000E+05
  CO2%L_sub    = 5.930000E+05
  CO2%gamma    = 1.294000E+00
  CO2%rho_l_bp = 1.032000E+03
  CO2%rho_l_tp = 1.110000E+03
  CO2%rho_s    = 1.562000E+03
  CO2%soc_index = 2
  
    ! CO

  CO%mmw      = 2.800970E+01
  CO%Rgas     = Rstar/CO%mmw
  CO%cp       = 1.04000E+03
  CO%tp_T     = 6.795000E+01
  CO%tp_P     = 1.530000E+04
  CO%crit_T   = 1.134450E+02
  CO%crit_P   = 3.498750E+06
  CO%L_vap_tp = 2.142857E+05
  CO%L_vap_bp = 0
  CO%L_fus    = 0
  CO%L_sub    = 2.7142857E+05
  CO%gamma    = 50.241546E+00
  CO%rho_l_bp = 0
  CO%rho_l_tp = 0
  CO%rho_s    = 0
  CO%soc_index = 5

    ! N2

  N2%mmw      =  2.800000E+01
  N2%Rgas     = Rstar/N2%mmw
  N2%cp       = 1.037000E+03
  N2%tp_T     = 6.314000E+01
  N2%tp_P     = 1.253000E+04
  N2%crit_T   = 1.262000E+02
  N2%crit_P   = 3.400000E+06
  N2%L_vap_tp = 2.180000E+05
  N2%L_vap_bp = 1.980000E+05
  N2%L_fus    = 2.573000E+04
  N2%L_sub    = 2.437000E+05
  N2%gamma    = 1.403000E+00
  N2%rho_l_bp = 8.086000E+02
  N2%rho_l_tp = 0
  N2%rho_s    = 1.026000E+03
  N2%soc_index = 13
  
    ! O2

  O2%mmw      =  3.200000E+01
  O2%Rgas     = Rstar/O2%mmw
  O2%cp       = 9.160000E+02
  O2%tp_T     = 5.430000E+01
  O2%tp_P     =1.500000E+02
  O2%crit_T   =1.545400E+02
  O2%crit_P   = 5.043000E+06
  O2%L_vap_tp = 2.420000E+05
  O2%L_vap_bp = 2.130000E+05
  O2%L_fus    = 1.390000E+04
  O2%L_sub    = 2.560000E+05
  O2%gamma    =  1.393000E+00
  O2%rho_l_bp = 1.141000E+03
  O2%rho_l_tp = 1.307000E+03
  O2%rho_s    = 1.351000E+03
  O2%soc_index = 999999
  
    ! H2

  H2%mmw      = 2.000000E+00
  H2%Rgas     = Rstar/H2%mmw
  H2%cp       =  1.423000E+04
  H2%tp_T     = 1.395000E+01
  H2%tp_P     = 7.200000E+03
  H2%crit_T   = 3.320000E+01
  H2%crit_P   = 1.298000E+06
  H2%L_vap_tp = 0
  H2%L_vap_bp = 4.540000E+05
  H2%L_fus    = 5.820000E+04
  H2%L_sub    = 0
  H2%gamma    = 1.384000E+00
  H2%rho_l_bp = 7.097000E+01
  H2%rho_l_tp = 0
  H2%rho_s    = 8.800000E+01
  H2%soc_index = 23
  
    ! He

  He%mmw      = 4.000000E+00
  He%Rgas     = Rstar/He%mmw
  He%cp       = 5.196000E+03
  He%tp_T     = 2.170000E+00
  He%tp_P     = 5.070000E+03
  He%crit_T   = 5.100000E+00
  He%crit_P   = 2.280000E+05
  He%L_vap_tp = 0
  He%L_vap_bp = 2.030000E+04
  He%L_fus    = 0
  He%L_sub    = 0
  He%gamma    = 1.664000E+00
  He%rho_l_bp = 1.249600+02
  He%rho_l_tp = 0
  He%rho_s    = 2.000000E+02
  He%soc_index = 24
  
    ! NH3

  NH3%mmw      = 1.700000E+01
  NH3%Rgas     = Rstar/NH3%mmw
  NH3%cp       =  2.060000E+03
  NH3%tp_T     = 1.954000E+02
  NH3%tp_P     =6.100000E+03
  NH3%crit_T   =4.055000E+02
  NH3%crit_P   = 1.128000E+07
  NH3%L_vap_tp = 1.658000E+06
  NH3%L_vap_bp = 1.371000E+06
  NH3%L_fus    = 3.314000E+05
  NH3%L_sub    = 1.989000E+06
  NH3%gamma    = 1.309000E+00
  NH3%rho_l_bp = 6.820000E+02
  NH3%rho_l_tp = 7.342000E+02
  NH3%rho_s    = 8.226000E+02
  NH3%soc_index = 11
  
    ! HCN

  HCN%mmw      = 27.0253
  HCN%Rgas     = Rstar/HCN%mmw
  HCN%cp       = 1.332E+03
  HCN%tp_T     = 259.86
  HCN%tp_P     = 1.87E+04
  HCN%crit_T   = 456.85
  HCN%crit_P   = 5.39E+06
  HCN%L_vap_tp = 0
  HCN%L_vap_bp = 9.33E+05
  HCN%L_fus    = 3.11E+05
  HCN%L_sub    = 1.32E+06
  HCN%gamma    = 1.3
  HCN%rho_l_bp = 688.00
  HCN%rho_l_tp = 0
  HCN%rho_s    = 0
  HCN%soc_index = 35
  
  Earth_air%cp = 1004
  Earth_air%mmw = 28.97
  Earth_air%gamma = 1.4003

  H2He_solar%mmw =  2.288166
  H2He_solar%cp  =  1.195451608e4
  
end subroutine init_gas_data

!!$!Planck function (of frequency)
!!$    FUNCTION B(nu,T) RESULT(res)
!!$        IMPLICIT NONE
!!$        REAL , INTENT(IN) :: nu,T
!!$        REAL  :: u,res
!!$        u=min(h*nu/(k*T),200)
!!$        res=(2*h*nu**3/c**2)/(exp(u)-1)
!!$    END FUNCTION B
!!$    
!!$    FUNCTION Blam(lam,T) RESULT(res)
!!$        IMPLICIT NONE
!!$        REAL , INTENT(IN) :: lam,T
!!$        REAL  :: u,res
!!$        u=h*c/(lam*k*T)
!!$        res=(2*h*c**2/lam**5)/(exp(u)-1)
!!$    END FUNCTION Blam
!!$
!!$    FUNCTION dB(nu,T) RESULT(res)
!!$        IMPLICIT NONE
!!$        REAL , INTENT(IN) :: nu,T
!!$        REAL  :: u,res
!!$        u =min(h*nu/(k*T),200)
!!$        res=(2*h**2*nu**4/(k*c**2*T**2))*(exp(u)/(exp(u)-1)**2)
!!$        !print*, T
!!$    END FUNCTION dB
    
    !Saturation vapor pressure over ice (Smithsonian formula)
    !Input: Kelvin. Output: Pascal
    FUNCTION satvpi(T) RESULT(res)
        IMPLICIT NONE
        REAL , INTENT(IN) :: T
        REAL  :: esbasi,tbasi,aa,b,c,e,esice,res
        !Compute es over ice (valid between -153 c and 0 c)
        !see smithsonian meteorological tables page 350
   
        !Original source: GFDL climate model, circa 1995
        esbasi = 6107.1
        tbasi  = 273.16
   
        aa     = -9.09718 *(tbasi/T-1.0)
        b      = -3.56654 *log10(tbasi/T)
        c      = 0.876793*(1.0-T/tbasi)
        e      = log10(esbasi)
        esice  = 10**(aa+b+c+e)
        res =.1*esice  !Convert to Pascals
    END FUNCTION satvpi

    !Saturation vapor pressure over liquid water (Smithsonian formula)
    !Input: Kelvin. Output: Pascal
    FUNCTION satvpw(T) RESULT(res)
        IMPLICIT NONE
        REAL , INTENT(IN) :: T
        REAL  :: esbasw,tbasw,aa,b,c,d,e,esh2O,res
        !compute es over liquid water between -20c and freezing.
        !see smithsonian meteorological tables page 350.
   
        !Original source: GFDL climate model, circa 1995
        esbasw = 1013246.0
        tbasw  = 373.16

        aa     = -7.90298*(tbasw/T-1)
        b      =  5.02808*log10(tbasw/T)
        c      = -1.3816D-07*(  10**( ((1-T/tbasw)*11.344)-1 )  )
        d      =  8.1328D-03*(  10**( ((tbasw/T-1)*(-3.49149))-1)  )
        e      = log10(esbasw)
        esh2O  = 10**(aa+b+c+d+e)
        res = .1*esh2O  !Convert to Pascals
    END FUNCTION satvpw

    ! An alternate formula for saturation vapor pressure over liquid water
    FUNCTION satvpw_Heymsfield(T) RESULT(res)
        IMPLICIT NONE
        REAL , INTENT(IN) :: T
        REAL  :: ts,sr,ar,br,cr,dw,er,vp,res
        ts = 373.16
        sr = 3.0057166
        ! Vapor pressure over water. Heymsfield formula
        ar = ts/T
        br = 7.90298*(ar-1)
        cr  = 5.02808*log10(ar);
        dw = (1.3816D-07)*(10**(11.344*(1-1/ar))-1)
        er = 8.1328D-03*((10**(-(3.49149*(ar-1))) )-1)
        vp = 10**(cr-dw+er+sr-br)
        vp = vp*1.0D+02
        res = vp
    END FUNCTION satvpw_Heymsfield
      
    FUNCTION  satvpg(T) RESULT(res)
        IMPLICIT NONE
        REAL , INTENT(IN) :: T
        REAL  :: res
         
        !This is the saturation vapor pressure computation used in the
        !GFDL climate model.  It blends over from water saturation to
        !ice saturation as the temperature falls below 0C.
        IF ((T-273.16) .LT. -20) THEN
            res=satvpi(T)
        ENDIF
        IF ( ((T-273.16) .GE. -20) .AND. ((T-273.16) .LE. 0)) THEN
            res=0.05*(273.16-T)*satvpi(T) + 0.05*(T-253.16)*satvpw(T)
        ENDIF
        IF ((T-273.16) .GT. 0) THEN
            res=satvpw(T)
        ENDIF
    END FUNCTION satvpg

    !Saturation vapor pressure for any substance, computed using
    !the simplified form of Clausius-Clapeyron assuming the perfect
    !gas law and constant latent heat
    FUNCTION satvps(T,T0,e0,MolecularWeight,LatentHeat) RESULT(res)
        IMPLICIT NONE
        REAL , INTENT(IN) :: T,T0,e0,MolecularWeight,LatentHeat
        REAL  :: Rv,res
        Rv = Rstar/MolecularWeight 
        res = e0*exp(-(LatentHeat/Rv)*(1/T - 1/T0))
    END FUNCTION satvps
END MODULE phys
