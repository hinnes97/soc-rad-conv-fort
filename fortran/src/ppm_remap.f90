module ppm_remap

  implicit none
  real, parameter :: r3 = 1./3., r23 = 2./3., r12 = 1./12.
contains

  subroutine remap(p_old, q_old, p_new, q_new)
    !--------------------------------------------------------------------------------------
    ! DESCRIPTION:
    ! Remaps from Lagrangian vertical coordinate (p_old) to new grid, p_new
    ! Uses the piecewise parabolic method from ExoFMS (Lin 2004). Code has largely been
    ! rewritten from the 3D FV3 version to work for a 1D grid
    !--------------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------
    ! INPUT
    !--------------------------------------------------------------------------------------
    real(dp), intent(in) :: p_old(:) ! Old pedge grid    
    real(dp), intent(in) :: p_new(:) ! New pedge grid to be remapped onto
    real(dp), intent(in) :: q_old(:) ! q defined on old grid
    
    !--------------------------------------------------------------------------------------
    ! OUTPUT
    !--------------------------------------------------------------------------------------
    real(dp), intent(out) :: q_new(:) ! Remapped variable

    !--------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    !--------------------------------------------------------------------------------------

    integer :: km = len(p_old) ! Length of old grid
    integer :: kn = len(p_new) ! Length of new grid
    
    integer :: k,l,m,k0 ! z indices
    
    
    real(dp) :: delp(km-1) ! old delp grid
    real(dp) :: delp_new(kn-1) ! new delp grid
    real(dp) :: q4(4,km) ! Interpolation parameters
    real(dp) :: pl, pr, qsum, diffp, esl
    
    do k=1,km
       delp(k) = p_old(k+1) - p_old(k)
       q4(1,k) = q(k)
    enddo

    do k=1,kn
       delp_new(k) = p_new(k+1) - p_new(k)
    enddo
    
    ! Compute vertical subgrid distribution
    call ppm_profile( q4, delp, km)

    ! Mapping
    k0=1
    do 555 k=1,kn
       do 100 l=k0,km
          ! locate the top edge: p_new(k)
          if(p_new(k) >= p_old(l) .and. p_new(k) <= p_old(l+1)) then
             pl = (p_new(i,k)-p_old(i,l)) / delp(i,l)
             if(p_new(i,k+1) <= p_old(i,l+1)) then
                ! entire new grid is within the original grid
                pr = (p_new(k+1)-p_old(l)) / delp(l)
                q_new(k) = q4(2,l) + 0.5*(q4(4,l)+q4(3,l)-q4(2,l))  &
                     *(pr+pl)-q4(4,l)*r3*(pr*(pr+pl)+pl**2)
                k0 = l
                goto 555
             else
                ! Fractional area...
                qsum = (p_old(l+1)-p_new(k))*(q4(2,l)+0.5*(q4(4,l)+   &
                          q4(3,l)-q4(2,l))*(1.+pl)-q4(4,l)*           &
                          (r3*(1.+pl*(1.+pl))))
                     do m=l+1,km
                        ! locate the bottom edge: pe2(i,k+1)
                        if(p_new(k+1) > p_old(m+1) ) then
                           ! Whole layer..
                           qsum = qsum + delp(m)*q4(1,m)
                        else
                           diffp = p_new(k+1)-p_old(m)
                           esl = diffp / delp(m)
                           qsum = qsum + dp*(q4(2,m)+0.5*esl*               &
                                (q4(3,m)-q4(2,m)+q4(4,m)*(1.-r23*esl)))
                           k0 = m
                           goto 123
                        endif
                     enddo
                  endif
                  goto 123
               endif
100         continue
123         q_new(k) = qsum / delp_new(k)
555      enddo


    end subroutine remap
    
    subroutine ppm_profile(a4, delp, km)
          
      ! !INPUT PARAMETERS:
      integer, intent(in):: km      ! vertical dimension
      integer, intent(in):: kord    ! Order (or more accurately method no.):
      ! 
      real , intent(in):: delp(km)     ! layer pressure thickness

      ! !INPUT/OUTPUT PARAMETERS:
      real , intent(inout):: a4(4,km)  ! Interpolated values

      ! DESCRIPTION:
      !
      !   Perform the piecewise parabolic reconstruction
      ! 
      ! !REVISION HISTORY: 
      ! S.-J. Lin   revised at GFDL 2007
      !-----------------------------------------------------------------------
      ! local arrays:
      real    dc(km)
      real    h2(km)
      real  delq(km)
      real   df2(km)
      real    d4(km)

      ! local scalars:
      integer i, k, km1, lmt, it
      real  fac
      real  a1, a2, c1, c2, c3, d1, d2
      real  qm, dq, lac, qmp, pmp

      km1 = km - 1
      it = i2 - i1 + 1

      do k=2,km
            delq(k-1) =   a4(1,k) - a4(1,k-1)
            d4(k) = delp(k-1) + delp(k)
      enddo

      do k=2,km1
            c1  = (delp(k-1)+0.5*delp(k))/d4(k+1)
            c2  = (delp(k+1)+0.5*delp(k))/d4(k)
            df2(k) = delp(k)*(c1*delq(k) + c2*delq(k-1)) /      &
                 (d4(k)+delp(k+1))
            dc(k) = sign( min(abs(df2(k)),              &
                 max(a4(1,k-1),a4(1,k),a4(1,k+1))-a4(1,k),  &
                 a4(1,k)-min(a4(1,k-1),a4(1,k),a4(1,k+1))), df2(k) )
      enddo

      !-----------------------------------------------------------
      ! 4th order interpolation of the provisional cell edge value
      !-----------------------------------------------------------

      do k=3,km1
            c1 = delq(k-1)*delp(k-1) / d4(k)
            a1 = d4(k-1) / (d4(k) + delp(k-1))
            a2 = d4(k+1) / (d4(k) + delp(k))
            a4(2,k) = a4(1,k-1) + c1 + 2./(d4(k-1)+d4(k+1)) *    &
                 ( delp(k)*(c1*(a1 - a2)+a2*dc(k-1)) -          &
                 delp(k-1)*a1*dc(k  ) )
      enddo


      ! Area preserving cubic with 2nd deriv. = 0 at the boundaries
      ! Top

      d1 = delp(1)
      d2 = delp(2)
      qm = (d2*a4(1,1)+d1*a4(1,2)) / (d1+d2)
      dq = 2.*(a4(1,2)-a4(1,1)) / (d1+d2)
      c1 = 4.*(a4(2,3)-qm-d2*dq) / ( d2*(2.*d2*d2+d1*(d2+3.*d1)) )
      c3 = dq - 0.5*c1*(d2*(5.*d1+d2)-3.*d1*d1)
      a4(2,2) = qm - 0.25*c1*d1*d2*(d2+3.*d1)
      ! Top edge:
      !-------------------------------------------------------
      a4(2,1) = d1*(2.*c1*d1**2-c3) + a4(2,2)
      !-------------------------------------------------------
      !        a4(2,1) = (12./7.)*a4(1,1)-(13./14.)*a4(1,2)+(3./14.)*a4(1,3)
      !-------------------------------------------------------
      ! No over- and undershoot condition
      a4(2,2) = max( a4(2,2), min(a4(1,1), a4(1,2)) )
      a4(2,2) = min( a4(2,2), max(a4(1,1), a4(1,2)) )
      dc(1) =  0.5*(a4(2,2) - a4(1,1))


      ! Enforce monotonicity of the "slope" within the top layer

      a4(2,1) = max(0., a4(2,1))
      a4(2,2) = max(0., a4(2,2))
      
      ! Bottom
      ! Area preserving cubic with 2nd deriv. = 0 at the surface
      d1 = delp(km)
      d2 = delp(km1)
      qm = (d2*a4(1,km)+d1*a4(1,km1)) / (d1+d2)
      dq = 2.*(a4(1,km1)-a4(1,km)) / (d1+d2)
      c1 = (a4(2,km1)-qm-d2*dq) / (d2*(2.*d2*d2+d1*(d2+3.*d1)))
      c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1*d1)
      a4(2,km) = qm - c1*d1*d2*(d2+3.*d1)
      ! Bottom edge:
      !-----------------------------------------------------
      a4(3,km) = d1*(8.*c1*d1**2-c3) + a4(2,km)
      !        dc(km) = 0.5*(a4(3,km) - a4(1,km))
      !-----------------------------------------------------
      !        a4(3,km) = (12./7.)*a4(1,km)-(13./14.)*a4(1,km-1)+(3./14.)*a4(1,km-2)
      ! No over- and under-shoot condition
      a4(2,km) = max( a4(2,km), min(a4(1,km), a4(1,km1)) )
      a4(2,km) = min( a4(2,km), max(a4(1,km), a4(1,km1)) )
      dc(km) = 0.5*(a4(1,km) - a4(2,km))

      ! Enforce constraint on the "slope" at the surface
      a4(2,km) = max(0.,a4(2,i,km))
      a4(3,km) = max(0.,a4(3,i,km))
      
      do k=1,km1
            a4(3,k) = a4(2,k+1)
      enddo

      !-----------------------------------------------------------
      ! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
      !-----------------------------------------------------------
      ! Top 2 and bottom 2 layers always use monotonic mapping
      do k=1,2
         a4(4,k) = 3.*(2.*a4(1,k) - (a4(2,k)+a4(3,k)))
         call ppm_limiters(dc(k), a4(1,k),  0)
      enddo

      if(kord >= 7) then
         !-----------------------
         ! Huynh's 2nd constraint
         !-----------------------
         do k=2,km1
               ! Method#1
               !           h2(k) = delq(k) - delq(k-1)
               ! Method#2 - better
               h2(k) = 2.*(dc(k+1)/delp(k+1) - dc(k-1)/delp(k-1))  &
                    / ( delp(k)+0.5*(delp(k-1)+delp(k+1)) )        &
                    * delp(k)**2 
               ! Method#3
!!!            h2(k) = dc(k+1) - dc(k-1)
         enddo

         fac = 1.5           ! original quasi-monotone

         do k=3,km-2
               ! Right edges
               !        qmp   = a4(1,k) + 2.0*delq(k-1)
               !        lac   = a4(1,k) + fac*h2(k-1) + 0.5*delq(k-1)
               !
               pmp   = 2.*dc(k)
               qmp   = a4(1,k) + pmp
               lac   = a4(1,k) + fac*h2(k-1) + dc(k)
               a4(3,k) = min(max(a4(3,k), min(a4(1,k), qmp, lac)),    &
                    max(a4(1,k), qmp, lac) )
               ! Left  edges
               !        qmp   = a4(1,k) - 2.0*delq(k)
               !        lac   = a4(1,k) + fac*h2(k+1) - 0.5*delq(k)
               !
               qmp   = a4(1,k) - pmp
               lac   = a4(1,k) + fac*h2(k+1) - dc(k)
               a4(2,k) = min(max(a4(2,k),  min(a4(1,k), qmp, lac)),   &
                    max(a4(1,k), qmp, lac))
               !-------------
               ! Recompute A6
               !-------------
               a4(4,k) = 3.*(2.*a4(1,k) - (a4(2,k)+a4(3,k)))

            ! Additional constraint to ensure positivity when kord=7
               call ppm_limiters(dc(k), a4(1,k), 2)
         enddo


      do k=km1,km
         a4(4,k) = 3.*(2.*a4(1,k) - (a4(2,k)+a4(3,k)))
         call ppm_limiters(dc(k), a4(1,k), 0)
      enddo

    end subroutine ppm_profile


    subroutine ppm_limiters(dm, a4, itot, lmt)

      ! !INPUT PARAMETERS:
      real , intent(in):: dm     ! the linear slope
      integer, intent(in) :: lmt       ! 0: Standard PPM constraint
      ! 1: Improved full monotonicity constraint (Lin)
      ! 2: Positive definite constraint
      ! 3: do nothing (return immediately)
      ! !INPUT/OUTPUT PARAMETERS:
      real , intent(inout) :: a4(4)   ! PPM array
      ! AA <-- a4(1,i)
      ! AL <-- a4(2,i)
      ! AR <-- a4(3,i)
      ! A6 <-- a4(4,i)
      ! !LOCAL VARIABLES:
      real  qmp
      real  da1, da2, a6da
      real  fmin
      integer i

      ! Developer: S.-J. Lin, NASA-GSFC
      ! Last modified: Apr 24, 2000

      if ( lmt == 3 ) return

      if(lmt == 0) then
         ! Standard PPM constraint
         if(dm(i) == 0.) then
            a4(2) = a4(1)
            a4(3) = a4(1)
            a4(4) = 0.
         else
            da1  = a4(3) - a4(2)
            da2  = da1**2
            a6da = a4(4)*da1
            if(a6da < -da2) then
               a4(4) = 3.*(a4(2)-a4(1))
               a4(3) = a4(2) - a4(4)
            elseif(a6da > da2) then
               a4(4) = 3.*(a4(3)-a4(1))
               a4(2) = a4(3) - a4(4)
            endif
         endif


      elseif (lmt == 1) then

         ! Improved full monotonicity constraint (Lin 2004)
         ! Note: no need to provide first guess of A6 <-- a4(4,i)
         qmp = 2.*dm
         a4(2) = a4(1)-sign(min(abs(qmp),abs(a4(2)-a4(1))), qmp)
         a4(3) = a4(1)+sign(min(abs(qmp),abs(a4(3)-a4(1))), qmp)
         a4(4) = 3.*( 2.*a4(1) - (a4(2)+a4(3)) )


      elseif (lmt == 2) then
         ! Positive definite constraint
         if( abs(a4(3)-a4(2)) < -a4(4) ) then
            fmin = a4(1)+0.25*(a4(3)-a4(2))**2/a4(4)+a4(4)*r12
            if( fmin < 0. ) then
               if(a4(1)<a4(3) .and. a4(1)<a4(2)) then
                  a4(3) = a4(1)
                  a4(2) = a4(1)
                  a4(4) = 0.
               elseif(a4(3) > a4(2)) then
                  a4(4) = 3.*(a4(2)-a4(1))
                  a4(3) = a4(2) - a4(4)
               else
                  a4(4) = 3.*(a4(3)-a4(1))
                  a4(2) = a4(3) - a4(4)
               endif
            endif
         endif
      endif

    end subroutine ppm_limiters
  
end module ppm_remap
