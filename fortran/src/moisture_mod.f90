module moisture_mod

  use params, only: dp, moisture_scheme, q0, nf
  use condense, only: q_sat, cold_trap
  use phys, only: T_TP => H2O_TriplePointT, P_TP => H2O_TriplePointP, &
       L_vap => H2O_L_vaporization_TriplePoint, Rstar, mu_v => H2O_MolecularWeight
  implicit none

contains
  subroutine get_q(pf, Tf, pe, q, olr)
    real(dp), intent(in) :: pf(:), pe(:)
    real(dp), intent(in) :: olr
    real(dp), intent(inout) :: Tf(:)
    real(dp), intent(out) :: q(:)
    
    real(dp) :: qsat(nf), delp(nf), test(nf)
    integer  :: mask(nf) ! 0 if radiative, 1 adiabatic, 2 pure steam
    logical :: adjust_mask(nf)
    logical :: temp_mask(nf)
    integer :: k, klim, ktrop_temp(1), ktrop
    select case(moisture_scheme)

    case('none')
       q = 0.0
    case('deep')
!       do k=1,nf
!          call sat(pf(k), Tf(k), qsat(k))
!       enddo
       
       do k=1,nf
          q(k) = min(qsat(k), q0)
       enddo
!       call cold_trap(q, 1)

    case('surface')

       adjust_mask= .false.
       mask = 0
       
       do k=1,nf
          delp(k) = pe(k+1) - pe(k)
       enddo

!!$       call q_sat(pf, Tf, test)
!!$       do k=1,nf
!!$          call sat(pf(k), Tf(k), qsat(k))
!!$!          write(*,*) Tf(k), qsat(k), test(k)
!!$       enddo
!!$       
!!$       call cold_trap(qsat, ktrop)
!!$
!!$       !write(*,*) 'TESTING' ,qsat
!!$       q = qsat
!!$       klim=nf
!!$       do k=nf,1,-1
!!$          if ( qsat(k) .gt. 1. ) then
!!$             ! Put on a pure steam adiabat
!!$             call dew_point_T(pf(k), Tf(k))
!!$             q(k) = 1.
!!$             mask(k) = 2
!!$          else
!!$             exit
!!$          endif
!!$       enddo
!!$
       
       
       !call adjust(pf, delp, Tf, q, adjust_mask, olr, ktrop)
!!$       write(*,*) 'AFTER ADJUST', adjust_mask
!!$
!!$       do k=1,nf
!!$          if (adjust_mask(k)) then
!!$             mask(k) = 1
!!$          endif
!!$       enddo

!       write(*,*) 'AFTER ADJUST', Tf
       ! Now do cold trap
!       call cold_trap(q, ktrop)


    end select
    
  end subroutine get_q

  
end module moisture_mod
