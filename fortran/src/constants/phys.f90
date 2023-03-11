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

!Set properties of individual gases

!Thermodynamic properties of dry Earth air
REAL , PARAMETER :: cp_dry_Earth_air = 1004
REAL , PARAMETER :: MolecularWeight_dry_Earth_air = 28.97
REAL , PARAMETER :: gamma_dry_Earth_air = 1.4003
!--------------------------------------------------------    
!H2O
REAL , PARAMETER :: H2O_CriticalPointT = 6.471000E+02
REAL , PARAMETER :: H2O_CriticalPointP = 2.210000E+07
REAL , PARAMETER :: H2O_TriplePointT = 2.731500E+02
REAL , PARAMETER :: H2O_TriplePointP = 6.110000E+02
REAL , PARAMETER :: H2O_L_vaporization_BoilingPoint = 2.255000E+06
REAL , PARAMETER :: H2O_L_vaporization_TriplePoint = 2.493000E+06
REAL , PARAMETER :: H2O_L_fusion = 3.340000E+05
REAL , PARAMETER :: H2O_L_sublimation = 2.840000E+06
REAL , PARAMETER :: H2O_rho_liquid_BoilingPoint = 9.584000E+02
REAL , PARAMETER :: H2O_rho_liquid_TriplePoint = 9.998700E+02
REAL , PARAMETER :: H2O_rho_solid = 9.170000E+02
REAL , PARAMETER :: H2O_cp = 1.847000E+03
REAL , PARAMETER :: H2O_gamma = 1.331000E+00
REAL , PARAMETER :: H2O_MolecularWeight = 1.800000E+01
CHARACTER (LEN = 5), PARAMETER :: H2O_name = 'Water'
CHARACTER (LEN = 3), PARAMETER :: H2O_formula = 'H2O'
REAL , PARAMETER :: H2O_L_vaporization=2.493000E+06
REAL , PARAMETER :: H2O_rho_liquid=9.998700E+02
!--------------------------------------------------------
!CH4
REAL , PARAMETER :: CH4_CriticalPointT = 1.904400E+02
REAL , PARAMETER :: CH4_CriticalPointP = 4.596000E+06
REAL , PARAMETER :: CH4_TriplePointT = 9.067000E+01
REAL , PARAMETER :: CH4_TriplePointP = 1.170000E+04
REAL , PARAMETER :: CH4_L_vaporization_BoilingPoint = 5.100000E+05
REAL , PARAMETER :: CH4_L_vaporization_TriplePoint = 5.360000E+05
REAL , PARAMETER :: CH4_L_fusion = 5.868000E+04
REAL , PARAMETER :: CH4_L_sublimation = 5.950000E+05
REAL , PARAMETER :: CH4_rho_liquid_BoilingPoint = 4.502000E+02
REAL , PARAMETER :: CH4_rho_liquid_TriplePoint = 0
REAL , PARAMETER :: CH4_rho_solid = 5.093000E+02
REAL , PARAMETER :: CH4_cp = 2.195000E+03
REAL , PARAMETER :: CH4_gamma = 1.305000E+00
REAL , PARAMETER :: CH4_MolecularWeight = 1.600000E+01
CHARACTER (LEN = 7), PARAMETER :: CH4_name = 'Methane'
CHARACTER (LEN = 3), PARAMETER :: CH4_formula = 'CH4'
REAL , PARAMETER :: CH4_L_vaporization=5.360000E+05
REAL , PARAMETER :: CH4_rho_liquid=4.502000E+02
!--------------------------------------------------------
!CO2
REAL , PARAMETER :: CO2_CriticalPointT = 3.042000E+02
REAL , PARAMETER :: CO2_CriticalPointP = 7.382500E+06
REAL , PARAMETER :: CO2_TriplePointT = 2.165400E+02
REAL , PARAMETER :: CO2_TriplePointP = 5.185000E+05
REAL , PARAMETER :: CO2_L_vaporization_BoilingPoint = 0
REAL , PARAMETER :: CO2_L_vaporization_TriplePoint = 3.970000E+05
REAL , PARAMETER :: CO2_L_fusion = 1.960000E+05
REAL , PARAMETER :: CO2_L_sublimation = 5.930000E+05
REAL , PARAMETER :: CO2_rho_liquid_BoilingPoint = 1.032000E+03
REAL , PARAMETER :: CO2_rho_liquid_TriplePoint = 1.110000E+03
REAL , PARAMETER :: CO2_rho_solid = 1.562000E+03
REAL , PARAMETER :: CO2_cp = 8.200000E+02
REAL , PARAMETER :: CO2_gamma = 1.294000E+00
REAL , PARAMETER :: CO2_MolecularWeight = 4.400000E+01
CHARACTER (LEN = 14), PARAMETER :: CO2_name = 'Carbon Eioxide'
CHARACTER (LEN = 3), PARAMETER :: CO2_formula = 'CO2'
REAL , PARAMETER :: CO2_L_vaporization=3.970000E+05
REAL , PARAMETER :: CO2_rho_liquid=1.110000E+03
!--------------------------------------------------------
!CO
REAL , PARAMETER :: CO_CriticalPointT = 1.134450E+02
REAL , PARAMETER :: CO_CriticalPointP = 3.498750E+06
REAL , PARAMETER :: CO_TriplePointT = 6.795000E+01
REAL , PARAMETER :: CO_TriplePointP = 1.530000E+04
REAL , PARAMETER :: CO_L_vaporization_BoilingPoint = 0
REAL , PARAMETER :: CO_L_vaporization_TriplePoint = 2.142857E+05
REAL , PARAMETER :: CO_L_fusion = 0
REAL , PARAMETER :: CO_L_sublimation = 2.7142857E+05
REAL , PARAMETER :: CO_rho_liquid_BoilingPoint = 0
REAL , PARAMETER :: CO_rho_liquid_TriplePoint = 0
REAL , PARAMETER :: CO_rho_solid = 0
REAL , PARAMETER :: CO_cp = 1.04000E+03
REAL , PARAMETER :: CO_gamma = 50.241546E+00
REAL , PARAMETER :: CO_MolecularWeight = 2.800970E+01
CHARACTER (LEN = 14), PARAMETER :: CO_name = 'Carbon Monoxide'
CHARACTER (LEN = 3), PARAMETER :: CO_formula = 'CO'
REAL , PARAMETER :: CO_L_vaporization=2.142857E+05
REAL , PARAMETER :: CO_rho_liquid = 0
!--------------------------------------------------------
!N2
REAL , PARAMETER :: N2_CriticalPointT = 1.262000E+02
REAL , PARAMETER :: N2_CriticalPointP = 3.400000E+06
REAL , PARAMETER :: N2_TriplePointT = 6.314000E+01
REAL , PARAMETER :: N2_TriplePointP = 1.253000E+04
REAL , PARAMETER :: N2_L_vaporization_BoilingPoint = 1.980000E+05
REAL , PARAMETER :: N2_L_vaporization_TriplePoint = 2.180000E+05
REAL , PARAMETER :: N2_L_fusion = 2.573000E+04
REAL , PARAMETER :: N2_L_sublimation = 2.437000E+05
REAL , PARAMETER :: N2_rho_liquid_BoilingPoint = 8.086000E+02
REAL , PARAMETER :: N2_rho_liquid_TriplePoint = 0
REAL , PARAMETER :: N2_rho_solid = 1.026000E+03
REAL , PARAMETER :: N2_cp = 1.037000E+03
REAL , PARAMETER :: N2_gamma = 1.403000E+00
REAL , PARAMETER :: N2_MolecularWeight = 2.800000E+01
CHARACTER (LEN = 8), PARAMETER :: N2_name = 'Nitrogen'
CHARACTER (LEN = 2), PARAMETER :: N2_formula = 'N2'
REAL , PARAMETER :: N2_L_vaporization=2.180000E+05
REAL , PARAMETER :: N2_rho_liquid=8.086000E+02
!--------------------------------------------------------
!O2
REAL , PARAMETER :: O2_CriticalPointT = 1.545400E+02
REAL , PARAMETER :: O2_CriticalPointP = 5.043000E+06
REAL , PARAMETER :: O2_TriplePointT = 5.430000E+01
REAL , PARAMETER :: O2_TriplePointP = 1.500000E+02
REAL , PARAMETER :: O2_L_vaporization_BoilingPoint = 2.130000E+05
REAL , PARAMETER :: O2_L_vaporization_TriplePoint = 2.420000E+05
REAL , PARAMETER :: O2_L_fusion = 1.390000E+04
REAL , PARAMETER :: O2_L_sublimation = 2.560000E+05
REAL , PARAMETER :: O2_rho_liquid_BoilingPoint = 1.141000E+03
REAL , PARAMETER :: O2_rho_liquid_TriplePoint = 1.307000E+03
REAL , PARAMETER :: O2_rho_solid = 1.351000E+03
REAL , PARAMETER :: O2_cp = 9.160000E+02
REAL , PARAMETER :: O2_gamma = 1.393000E+00
REAL , PARAMETER :: O2_MolecularWeight = 3.200000E+01
CHARACTER (LEN = 6), PARAMETER :: O2_name = 'Oxygen'
CHARACTER (LEN = 2), PARAMETER :: O2_formula = 'O2'
REAL , PARAMETER :: O2_L_vaporization=2.420000E+05
REAL , PARAMETER :: O2_rho_liquid=1.307000E+03
!--------------------------------------------------------
!H2
REAL , PARAMETER :: H2_CriticalPointT = 3.320000E+01
REAL , PARAMETER :: H2_CriticalPointP = 1.298000E+06
REAL , PARAMETER :: H2_TriplePointT = 1.395000E+01
REAL , PARAMETER :: H2_TriplePointP = 7.200000E+03
REAL , PARAMETER :: H2_L_vaporization_BoilingPoint = 4.540000E+05
REAL , PARAMETER :: H2_L_vaporization_TriplePoint = 0
REAL , PARAMETER :: H2_L_fusion = 5.820000E+04
REAL , PARAMETER :: H2_L_sublimation = 0
REAL , PARAMETER :: H2_rho_liquid_BoilingPoint = 7.097000E+01
REAL , PARAMETER :: H2_rho_liquid_TriplePoint = 0
REAL , PARAMETER :: H2_rho_solid = 8.800000E+01
REAL , PARAMETER :: H2_cp = 1.423000E+04
REAL , PARAMETER :: H2_gamma = 1.384000E+00
REAL , PARAMETER :: H2_MolecularWeight = 2.000000E+00
CHARACTER (LEN = 8), PARAMETER :: H2_name = 'Hydrogen'
CHARACTER (LEN = 2), PARAMETER :: H2_formula = 'H2'
REAL , PARAMETER :: H2_L_vaporization=4.540000E+05
REAL , PARAMETER :: H2_rho_liquid=7.097000E+01
!--------------------------------------------------------
!He
REAL , PARAMETER :: He_CriticalPointT = 5.100000E+00
REAL , PARAMETER :: He_CriticalPointP = 2.280000E+05
REAL , PARAMETER :: He_TriplePointT = 2.170000E+00
REAL , PARAMETER :: He_TriplePointP = 5.070000E+03
REAL , PARAMETER :: He_L_vaporization_BoilingPoint = 2.030000E+04
REAL , PARAMETER :: He_L_vaporization_TriplePoint = 0
REAL , PARAMETER :: He_L_fusion = 0
REAL , PARAMETER :: He_L_sublimation = 0
REAL , PARAMETER :: He_rho_liquid_BoilingPoint = 1.249600+02
REAL , PARAMETER :: He_rho_liquid_TriplePoint = 0
REAL , PARAMETER :: He_rho_solid = 2.000000E+02
REAL , PARAMETER :: He_cp = 5.196000E+03
REAL , PARAMETER :: He_gamma = 1.664000E+00
REAL , PARAMETER :: He_MolecularWeight = 4.000000E+00
CHARACTER (LEN = 6), PARAMETER :: He_name = 'Helium'
CHARACTER (LEN = 2), PARAMETER :: He_formula = 'He'
REAL , PARAMETER :: He_L_vaporization=2.030000E+04
REAL , PARAMETER :: He_rho_liquid=1.249600E+02
!--------------------------------------------------------
!NH3
REAL , PARAMETER :: NH3_CriticalPointT = 4.055000E+02
REAL , PARAMETER :: NH3_CriticalPointP = 1.128000E+07
REAL , PARAMETER :: NH3_TriplePointT = 1.954000E+02
REAL , PARAMETER :: NH3_TriplePointP = 6.100000E+03
REAL , PARAMETER :: NH3_L_vaporization_BoilingPoint = 1.371000E+06
REAL , PARAMETER :: NH3_L_vaporization_TriplePoint = 1.658000E+06
REAL , PARAMETER :: NH3_L_fusion = 3.314000E+05
REAL , PARAMETER :: NH3_L_sublimation = 1.989000E+06
REAL , PARAMETER :: NH3_rho_liquid_BoilingPoint = 6.820000E+02
REAL , PARAMETER :: NH3_rho_liquid_TriplePoint = 7.342000E+02
REAL , PARAMETER :: NH3_rho_solid = 8.226000E+02
REAL , PARAMETER :: NH3_cp = 2.060000E+03
REAL , PARAMETER :: NH3_gamma = 1.309000E+00
REAL , PARAMETER :: NH3_MolecularWeight = 1.700000E+01
CHARACTER (LEN = 7), PARAMETER :: NH3_name = 'Ammonia'
CHARACTER (LEN = 3), PARAMETER :: NH3_formula = 'NH3'
REAL , PARAMETER :: NH3_L_vaporization=1.658000E+06
REAL , PARAMETER :: NH3_rho_liquid=7.342000E+02

! H2He solar mix
REAL, PARAMETER :: H2He_solar_MolecularWeight = 2.288166
REAL, PARAMETER :: H2He_solar_cp = 1.195451608e4

CONTAINS
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
