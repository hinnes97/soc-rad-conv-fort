module ppm_remap
  use params, only: dp
  implicit none
  real, parameter :: r3 = 1./3., r23 = 2./3., r12 = 1./12.
contains


  subroutine map1_q2(km,   pe1,   q1,            &
                    kn,   pe2,   q2,   dp2)

! !INPUT PARAMETERS:
      integer, intent(in) :: km                ! Original vertical dimension
      integer, intent(in) :: kn                ! Target vertical dimension

      real(dp), intent(in) ::  pe1(km+1)     ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the original vertical coordinate
      real(dp), intent(in) ::  pe2(kn+1)     ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the new vertical coordinate
      real(dp), intent(in) ::  q1(km) ! Field input
      real(dp), intent(in) ::  dp2(kn)
! !INPUT/OUTPUT PARAMETERS:
      real(dp), intent(inout):: q2(kn) ! Field output
! !LOCAL VARIABLES:
      real(dp)   dp1(km)
      real(dp)   q4(4,km)
      real(dp)   pl, pr, qsum, dp_k, esl

      integer i, k, l, m, k0

      do k=1,km
            dp1(k) = pe1(k+1) - pe1(k)
            q4(1,k) = q1(k)
      enddo

! Compute vertical subgrid distribution
!   if ( kord >7 ) then
!        call  cs_profile( q4, dp1, km, i1, i2, iv )
!   else
        call ppm_profile( q4, dp1, km )
!   endif

! Mapping
         k0 = 1
      do 555 k=1,kn
      do 100 l=k0,km
! locate the top edge: pe2(i,k)
      if(pe2(k) >= pe1(l) .and. pe2(k) <= pe1(l+1)) then
         pl = (pe2(k)-pe1(l)) / dp1(l)
         if(pe2(k+1) <= pe1(l+1)) then
! entire new grid is within the original grid
            pr = (pe2(k+1)-pe1(l)) / dp1(l)
            q2(k) = q4(2,l) + 0.5*(q4(4,l)+q4(3,l)-q4(2,l))  &
                       *(pr+pl)-q4(4,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
          else
! Fractional area...
            qsum = (pe1(l+1)-pe2(k))*(q4(2,l)+0.5*(q4(4,l)+   &
                    q4(3,l)-q4(2,l))*(1.+pl)-q4(4,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if(pe2(k+1) > pe1(m+1) ) then
                                                   ! Whole layer..
                    qsum = qsum + dp1(m)*q4(1,m)
                 else
                     dp_k = pe2(k+1)-pe1(m)
                    esl = dp_k / dp1(m)
                   qsum = qsum + dp_k*(q4(2,m)+0.5*esl*               &
                       (q4(3,m)-q4(2,m)+q4(4,m)*(1.-r23*esl)))
                   k0 = m
                   goto 123
                 endif
              enddo
              goto 123
          endif
      endif
100   continue
123   q2(k) = qsum / dp2(k)
555   continue

 end subroutine map1_q2

  subroutine ppm_profile(a4, delp, km)

! !INPUT PARAMETERS:
 integer, intent(in):: km      ! vertical dimension
 real(dp) , intent(in):: delp(km)     ! layer pressure thickness

! !INPUT/OUTPUT PARAMETERS:
 real(dp) , intent(inout):: a4(4,km)  ! Interpolated values

! DESCRIPTION:
!
!   Perform the piecewise parabolic reconstruction
! 
! !REVISION HISTORY: 
! S.-J. Lin   revised at GFDL 2007
!-----------------------------------------------------------------------
! local arrays:
      real(dp)    dc(km)
      real(dp)    h2(km)
      real(dp)  delq(km)
      real(dp)   df2(km)
      real(dp)    d4(km)

! local scalars:
      integer i, k, km1, lmt, it
      real(dp)  fac
      real(dp)  a1, a2, c1, c2, c3, d1, d2
      real(dp) qm, dq, lac, qmp, pmp

      km1 = km - 1
      do k=2,km
            delq(k-1) =   a4(1,k) - a4(1,k-1)
              d4(k  ) = delp(k-1) + delp(k)
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

!      if(km>8 .and. kord>3) call steepz(i1, i2, km, a4, df2, dc, delq, delp, d4)

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
!        a4(2,i,1) = (12./7.)*a4(1,i,1)-(13./14.)*a4(1,i,2)+(3./14.)*a4(1,i,3)
!-------------------------------------------------------
! No over- and undershoot condition
         a4(2,2) = max( a4(2,2), min(a4(1,1), a4(1,2)) )
         a4(2,2) = min( a4(2,2), max(a4(1,1), a4(1,2)) )
         dc(1) =  0.5*(a4(2,2) - a4(1,1))


! Enforce monotonicity of the "slope" within the top layer

 !     if( iv==0 ) then
            a4(2,1) = max(0., a4(2,1))
            a4(2,2) = max(0., a4(2,2))
!!$  !    elseif( iv==-1 ) then
!!$         do i=i1,i2
!!$            if ( a4(2,i,1)*a4(1,i,1) <= 0. ) then
!!$                 a4(2,i,1) = 0.
!!$            endif
!!$         enddo
!!$      endif


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
!        dc(i,km) = 0.5*(a4(3,i,km) - a4(1,i,km))
!-----------------------------------------------------
!        a4(3,i,km) = (12./7.)*a4(1,i,km)-(13./14.)*a4(1,i,km-1)+(3./14.)*a4(1,i,km-2)
! No over- and under-shoot condition
         a4(2,km) = max( a4(2,km), min(a4(1,km), a4(1,km1)) )
         a4(2,km) = min( a4(2,km), max(a4(1,km), a4(1,km1)) )
         dc(km) = 0.5*(a4(1,km) - a4(2,km))



! Enforce constraint on the "slope" at the surface

!!$#ifdef BOT_MONO
!!$      do i=i1,i2
!!$         a4(4,i,km) = 0
!!$         if( a4(3,i,km) * a4(1,i,km) <= 0. ) a4(3,i,km) = 0.
!!$         d1 = a4(1,i,km) - a4(2,i,km)
!!$         d2 = a4(3,i,km) - a4(1,i,km)
!!$         if ( d1*d2 < 0. ) then
!!$              a4(2,i,km) = a4(1,i,km)
!!$              a4(3,i,km) = a4(1,i,km)
!!$         else
!!$              dq = sign(min(abs(d1),abs(d2),0.5*abs(delq(i,km-1))), d1)
!!$              a4(2,i,km) = a4(1,i,km) - dq
!!$              a4(3,i,km) = a4(1,i,km) + dq
!!$         endif
!!$      enddo
!!$#else
 !     if( iv==0 ) then
             a4(2,km) = max(0.,a4(2,km))
             a4(3,km) = max(0.,a4(3,km))
!!$      elseif( iv==-1 ) then
!!$          do i=i1,i2
!!$             if ( a4(1,i,km)*a4(3,i,km) <= 0. ) then
!!$                  a4(3,i,km) = 0.
!!$             endif
!!$          enddo
!!$      endif
!!$!#endif

   do k=1,km1
         a4(3,k) = a4(2,k+1)
   enddo

!-----------------------------------------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
!-----------------------------------------------------------
! Top 2 and bottom 2 layers always use monotonic mapping
      do k=1,2
         a4(4,k) = 3.*(2.*a4(1,k) - (a4(2,k)+a4(3,k)))
         call ppm_limiters(dc(k), a4(1,k),0)
      enddo

!      if(kord >= 7) then
!-----------------------
! Huynh's 2nd constraint
!-----------------------
      do k=2,km1
! Method#1
!           h2(i,k) = delq(i,k) - delq(i,k-1)
! Method#2 - better
            h2(k) = 2.*(dc(k+1)/delp(k+1) - dc(k-1)/delp(k-1))  &
                     / ( delp(k)+0.5*(delp(k-1)+delp(k+1)) )        &
                     * delp(k)**2 
! Method#3
!!!            h2(i,k) = dc(i,k+1) - dc(i,k-1)
      enddo

      fac = 1.5           ! original quasi-monotone

      do k=3,km-2
! Right edges
!        qmp   = a4(1,i,k) + 2.0*delq(i,k-1)
!        lac   = a4(1,i,k) + fac*h2(i,k-1) + 0.5*delq(i,k-1)
!
         pmp   = 2.*dc(k)
         qmp   = a4(1,k) + pmp
         lac   = a4(1,k) + fac*h2(k-1) + dc(k)
         a4(3,k) = min(max(a4(3,k), min(a4(1,k), qmp, lac)),    &
                                        max(a4(1,k), qmp, lac) )
! Left  edges
!        qmp   = a4(1,i,k) - 2.0*delq(i,k)
!        lac   = a4(1,i,k) + fac*h2(i,k+1) - 0.5*delq(i,k)
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
         call ppm_limiters(dc(k), a4(1,k),2)
      enddo

!!$      else
!!$
!!$         lmt = kord - 3
!!$         lmt = max(0, lmt)
!!$         if (iv == 0) lmt = min(2, lmt)
!!$
!!$         do k=3,km-2
!!$            if( kord /= 4) then
!!$              do i=i1,i2
!!$                 a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
!!$              enddo
!!$            endif
!!$            if(kord/=6) call ppm_limiters(dc(i1,k), a4(1,i1,k), it, lmt)
!!$         enddo
!!$      endif

      do k=km1,km
            a4(4,k) = 3.*(2.*a4(1,k) - (a4(2,k)+a4(3,k)))
         call ppm_limiters(dc(k), a4(1,k), 0)
      enddo

 end subroutine ppm_profile


 subroutine ppm_limiters(dm, a4, lmt)

! !INPUT PARAMETERS:
      real(dp) , intent(in):: dm     ! the linear slope
      integer, intent(in) :: lmt       ! 0: Standard PPM constraint
                                       ! 1: Improved full monotonicity constraint (Lin)
                                       ! 2: Positive definite constraint
                                       ! 3: do nothing (return immediately)
! !INPUT/OUTPUT PARAMETERS:
      real(dp) , intent(inout) :: a4(4)   ! PPM array
                                           ! AA <-- a4(1,i)
                                           ! AL <-- a4(2,i)
                                           ! AR <-- a4(3,i)
                                           ! A6 <-- a4(4,i)
! !LOCAL VARIABLES:
      real(dp)  qmp
      real(dp)  da1, da2, a6da
      real(dp)  fmin
      integer i

! Developer: S.-J. Lin, NASA-GSFC
! Last modified: Apr 24, 2000

      if ( lmt == 3 ) return

      if(lmt == 0) then
! Standard PPM constraint
      if(dm == 0.) then
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
