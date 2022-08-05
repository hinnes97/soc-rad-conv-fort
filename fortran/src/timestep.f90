module timestep
  use params, only: dp, Nt, nf, ne, const, Finc, Fint, q0, surface, A_s, sb, surf_const, &
       moisture_scheme, sb, grav, cpair, del_time, accelerate, cp_s, rho_s, rdgas, &
       depth, U, C_d
  use flux_mod, only: get_fluxes
  use utils, only: linear_log_interp, bezier_interp
  use convection, only: dry_adjust
  use condense, only : rain_out, cold_trap,q_sat, dew_point_T
  use moisture_mod, only : get_q
  use adjust_mod, only : calc_q_and_grad, new_adjust
  use io, only: dump_data
  implicit none

contains
  subroutine step(Tf, pf, pe, net_F, dT, olr, file_name, q, fup, fdn, s_dn, s_up, Ts)
    character(*), intent(in) :: file_name
    real(dp), dimension(:), intent(inout) :: Tf
    real(dp), dimension(:), intent(out) :: dT
    real(dp), dimension(:), intent(in) :: pf
    real(dp), dimension(:), intent(in) :: pe
    real(dp), dimension(:), intent(out) :: net_F
    real(dp), intent(out) :: olr
    real(dp), intent(inout) :: q(:), fup(:), fdn(:), s_dn(:), s_up(:)
    real(dp), intent(inout) :: Ts

    integer :: i,j, l, ktrop, test
    real(dp), dimension(size(pe)) :: Te, q_half
    real(dp), dimension(size(Tf)) :: Tf_half, temp, dT_old, qsat,delp, grad
    logical, dimension(size(Tf)) :: dry_mask
    real(dp), dimension(size(Tf)) :: factor, dflux
    integer, dimension(size(Tf)) :: oscillate

    real(dp) :: Ts_half, dT_surf, time_const, dew_pt, Sens
    logical :: loop_cyc
    factor = const
    temp = 0
    oscillate = 0
    dT = 1._dp
    dry_mask = .false.

    olr = Finc ! Just for first get_q

    do i=1,nf
       delp(i) = pe(i+1) - pe(i)
    enddo
    
   
    !call get_q(pf, Tf, pe, q, olr)
    do i=1,nf
       write(*,*) 'Before', Tf(i), q(i) , pf(i)
    enddo
    call get_fluxes(nf, ne, Tf, pf, Te, pe, &
          net_F, 1.0_dp, Finc, Fint, olr, q, Ts, fup, fdn, s_dn, s_up)

    do l=1,1

        if (moisture_scheme == 'surface') then
           write(*,*) '-------------------------------------------------------'
           write(*,*) 'Recalculating q and gradient'
           write(*,*) '-------------------------------------------------------'
           call calc_q_and_grad(pf, delp, Tf, q, dry_mask, olr, grad, ktrop)
        endif
       

   Do j =1,Nt
       
     temp = Tf
     do i = 2, nf-1
       call bezier_interp(pf(i-1:i+1), Tf(i-1:i+1), 3, pe(i), Te(i))
       !call linear_log_interp(pe(i), pf(i-1), pf(i), Tf(i-1), Tf(i), Te(i))
     end do
     call bezier_interp(pf(nf-2:nf), Tf(nf-2:nf), 3, pe(nf), Te(nf))
    
    call linear_log_interp(pe(1), pf(1), pf(2), Tf(1), Tf(2), Te(1))
    call linear_log_interp(pe(ne), pf(nf-1), pf(nf), Tf(nf-1), Tf(nf), Te(ne) )

    ! do i=1,nf
    !    write(*,*) Tf(i), Te(i), q(i)
    ! enddo
    
     call get_fluxes(nf, ne, Tf, pf, Te, pe, &
          net_F, 1.0_dp, Finc, Fint, olr, q, Ts, fup, fdn, s_dn, s_up)
     !write(*,*) 'bottom boundary net flux', fup(ne)-fdn(ne), fup(1)
     do i=1,nf
        dT_old(i) = dT(i)

        if (accelerate) then
           time_const = const
        else
           time_const = grav/cpair/(pe(i+1) - pe(i))*del_time
        endif
        
        
        dT(i) = time_const*(net_F(i+1) - net_F(i))!*(pe(i+1) - pe(i))/pe(nf+1)
        !dT(i) = factor(i)*(net_F(i+1) - net_F(i))/(abs(net_F(i+1) - net_F(i))**0.9_dp)
        !*(pe(i+1) - pe(i))/pe(nf+1)
        if (dT(i)>5.0_dp) then
           dT(i) = 5.0_dp
        endif
        if (dT(i)<-5.0_dp) then
           dT(i) = -5.0_dp
        endif
     end do

     Tf_half = Tf + 0.5_dp*dT
     
     do i = 2, nf-1
       call bezier_interp(pf(i-1:i+1), Tf_half(i-1:i+1), 3, pe(i), Te(i))
       !call linear_log_interp(pe(i), pf(i-1), pf(i), Tf_half(i-1), Tf_half(i), Te(i))
    end do
    call bezier_interp(pf(nf-2:nf), Tf_half(nf-2:nf), 3, pe(nf), Te(nf))
    
    call linear_log_interp(pe(1), pf(1), pf(2), Tf(1), Tf(2), Te(1))
    !call linear_log_interp(pe(ne), pf(nf-1), pf(nf), Tf(nf-1), Tf(nf), Te(ne) )
    Te(ne) = Ts

    call get_fluxes(nf, ne, Tf_half, pf, Te, pe,  &
          net_F, 1.0_dp, Finc, Fint, olr, q, Ts, fup, fdn, s_dn, s_up)

    
    if (.not. surface) then 
       net_F(ne) = Fint
    endif
    

    do i=1,nf
       if (accelerate) then
          time_const = const
       else
          time_const = grav/cpair/(pe(i+1) - pe(i))*del_time
       endif
       
        dT(i) = time_const*(net_F(i+1) - net_F(i))

        if (dT(i)>5.0_dp) then
           dT(i) = 5.0_dp
        endif
        if (dT(i)<-5.0_dp) then
           dT(i) = -5.0_dp
        endif
        Tf(i) = Tf(i) + dT(i)

        if (Tf(i) < 50._dp) then
           Tf(i) = 50._dp
        else if (Tf(i) .gt. 10000. ) then
           Tf(i) = 10000._dp
        end if
     
     end do

     dT_surf = 0
     
     if (surface) then
        ! Do turbulent heat exchange between surface and lowest atmospheric layer
        Sens = cpair*pf(nf)/rdgas/Tf(nf) * C_d * U * (Tf(nf) - Ts)

        dT_surf = ( s_dn(ne)*(1-A_s) - sb*Ts**4 + fdn(ne)  + Sens)/cp_s/rho_s/depth * del_time
        Tf(nf) = Tf(nf)  - Sens*grav/cpair/(pe(ne) - pe(ne-1))*del_time
        Ts = Ts + dT_surf
        call dew_point_T(pe(ne), dew_pt)
        if (Ts .gt. dew_pt) then
           Ts = dew_pt
        endif

     endif


     if (moisture_scheme == 'surface') then
        call calc_q_and_grad(pf, delp, Tf, q, dry_mask, olr, grad, ktrop)
        call new_adjust(pf, delp, Tf, q, ktrop, grad, olr, dry_mask)
     endif
     
     do i = 2, nf-1
       call bezier_interp(pf(i-1:i+1), Tf_half(i-1:i+1), 3, pe(i), Te(i))
    end do

    call bezier_interp(pf(nf-2:nf), Tf_half(nf-2:nf), 3, pe(nf), Te(nf))

    call linear_log_interp(pe(1), pf(1), pf(2), Tf(1), Tf(2), Te(1))
    !call linear_log_interp(pe(ne), pf(nf-1), pf(nf), Tf(nf-1), Tf(nf), Te(ne) )
    Te(ne) = Ts
if (.true.) then
if (mod(j,1000) .eq. 0) then
    do i=1,nf
       if (dry_mask(i)) then
          dflux(i) = 0.
       else
          if (i.eq.nf) then
             dflux(i) = (net_F(i+1) - net_F(i) - Sens)/sb/Tf(i)**4
          else
             dflux(i) = (net_F(i+1) - net_F(i))/sb/Tf(i)**4
          endif
       endif
    end do

    do i=nf,1,-1
       if (dry_mask(i)) then
          test = i
       endif
    enddo
    !write(*,*) dry_mask
    write(*,*) j, 'max dflux', maxval(abs(dflux)), maxloc(abs(dflux)), maxval(abs(dT)), test, (net_F(1) - Fint)/Finc, Tf(nf), &
         net_F(1) - Fint, Ts
    if (maxval(abs(dflux))< 1.e-3 .and. abs((net_F(1) - Fint)/Finc) .lt. 1.e-2) then
       write(*,*) 'delta_Ts = ', dT_surf
       exit
    endif

    !write(*,*) dry_mask
    !write(*,*) q
    !    print*, 'near end', log(Tf(151)/Tf(150))/log(pf(151)/pf(150))
    call dump_data(file_name, nf,ne,Tf,pf,pe,olr,Finc,Fint,Te,q,fup,fdn,s_dn,s_up,Ts)
 endif
endif

 end do

 write(*,*) 'Surface temperature', Te(ne)

enddo

write(*,*) 'END T', Tf
write(*,*) 'end gradient'
do i=1,nf
   dflux(i) = (net_F(i+1) - net_F(i))/sb/Tf(i)**4
   print*, i, dflux(i), Tf(i), q(i)
enddo

call calc_q_and_grad(pf, delp, Tf, q, dry_mask, olr, grad, ktrop)
do i =1,nf-1
   
   write(*,*) 'grad', log(Tf(i+1)/Tf(i))/log(pf(i+1)/pf(i)), grad(i)
enddo

print*, 'end OLR', olr
!call calc_q_and_grad(pf, delp, Tf, q, dry_mask, olr, grad, ktrop)
!write(*,*) q
write(*,*) dry_mask
  end subroutine step
  
end module timestep
