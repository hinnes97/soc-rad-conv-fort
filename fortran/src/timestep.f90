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
  subroutine step(Tf, pf, pe, file_name, q, Ts)
    !====================================================================================
    ! Description
    !====================================================================================
    !
    ! Perform timestep
    !
    !====================================================================================
    ! Input variables
    !====================================================================================
    
    real(dp), dimension(:), intent(in) :: pf ! Layer pressure
    real(dp), dimension(:), intent(in) :: pe ! Level pressure

    character(*), intent(in) :: file_name ! Output file name
    
    !====================================================================================
    ! Input/Output variables
    !====================================================================================
    
    real(dp), dimension(:), intent(inout) :: Tf  ! Layer temperatures
    real(dp), dimension(:), intent(inout) :: q   ! Specific humidity
    
    real(dp), intent(inout) :: Ts ! Surface temperature

    

    integer :: i,j, l, ktrop, test, m
    real(dp), dimension(size(pe)) :: Te, q_half
    real(dp), dimension(size(Tf)) :: Tf_half, temp, dT_old, qsat,delp, grad, dT
    logical, dimension(size(Tf)) :: dry_mask
    real(dp), dimension(size(Tf)) :: factor,flux_diff, fsmooth_sum, fsmooth
    integer, dimension(size(Tf)) :: oscillate
    real(dp), dimension(size(Tf)) :: dflux
    
    real(dp), dimension(size(Te)) :: fup, fdn, s_up, s_dn, net_F
    
    real(dp) :: Ts_half, dT_surf, time_const, dew_pt, Sens, factor_surf, df_surf, t_surf_old, olr
    integer :: conv_start_step
    
    logical :: loop_cyc, conv_switch, sensible_heat
    conv_start_step = 2000
    factor = -50.
    factor_surf = -20.
    temp = 0
    oscillate = 0
    dT = 1._dp
    dry_mask = .false.
    conv_switch  = .false.
    sensible_heat = .false.

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


        if (moisture_scheme == 'surface') then
           write(*,*) '-------------------------------------------------------'
           write(*,*) 'Recalculating q and gradient'
           write(*,*) '-------------------------------------------------------'
           call calc_q_and_grad(pf, delp, Tf, q, dry_mask, olr, grad, ktrop)
        endif
       
   Do j =1,Nt
      !if (j .gt. 0 ) then
      !   accelerate = .false.
      !endif
      

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
        !dT_old(i) = dT(i)

        if (accelerate) then
           if (factor(i) .lt. 0) then
             factor(i) = 1.e-2
          endif

          if (mod(j,10000) .eq. 0) then
             factor(i) = 1.e-4
          endif
          
          time_const = factor(i) * pf(nf) / max((abs(net_F(i+1) - net_F(i))*1000. ), 1.e-30)**0.5 / (pe(ne) - pe(nf))

        else
           time_const = grav/cpair/(pe(i+1) - pe(i))*del_time
        endif
        
        
        dT(i) = time_const*(net_F(i+1) - net_F(i))!*(pe(i+1) - pe(i))/pe(nf+1)
!        write(*,*) dT(i)
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
    call linear_log_interp(pe(ne), pf(nf-1), pf(nf), Tf(nf-1), Tf(nf), Te(ne) )

    if (surface) then
       Te(ne) = Ts
    else
       Ts = Te(ne)
    endif
    
    
    call get_fluxes(nf, ne, Tf_half, pf, Te, pe,  &
          net_F, 1.0_dp, Finc, Fint, olr, q, Ts, fup, fdn, s_dn, s_up)

    
    if (.not. surface) then 
       net_F(ne) = Fint
    endif

    
    flux_diff(:) = 0.
    
    do i=1,nf
       if (accelerate) then
          if ((i .gt. 1) .and. (i .lt. nf)) then
          ! Calculate smoothing flux
             fsmooth(i) = 0.!1.e-3 * (0.5*(Tf(i-1) + Tf(i+1)) - Tf(i) )**7
             flux_diff(i) = flux_diff(i) + fsmooth(i)
          endif

          flux_diff(i) = flux_diff(i) + net_F(i+1) - net_F(i)


          if ((i .eq. nf) .and. surface) then
             
                     ! Do turbulent heat exchange between surface and lowest atmospheric layer
                     Sens = cpair*pf(nf)/rdgas/Tf(nf) * C_d * U * (Tf(nf) - Ts)
                     Sens = 0
                     df_surf = ( s_dn(ne)*(1-A_s) - fup(ne) + fdn(ne)  + Sens)
                     
                     if (factor_surf .lt. 0.) factor_surf = 1.e-3

                     dT_surf = factor_surf * pf(nf)/(abs(df_surf)*1000.)**0.5/(pe(ne) - pe(nf)) * df_surf
                     
                     if (abs(dT_surf) .gt. 5) then
                        dT_surf = 5* dT_surf/abs(dT_surf)
                     endif
                     
                     if (mod(j,6) .eq. 0) then
                        t_surf_old = ts
                     endif
                     if (mod(j,6) .eq. 5) then
                        if ( abs(Ts - t_surf_old) .lt. 3.*abs(dT_surf)) then
                           factor_surf = factor_surf/1.5
                        else
                           factor_surf = factor_surf*1.1
                        endif
                        if (factor_surf .lt. 1.e-10) factor_surf = 1.e-10
                     endif
                     
                     Ts = Ts + dT_surf
                     call dew_point_T(pe(ne), dew_pt)
                     if (Ts .gt. dew_pt) then
                        Ts = dew_pt
                     endif

                  endif
                  


                  !write(*,*) flux_diff(i), net_F(i+1) - net_f(i)
                  ! Require flux_diff not be zero
          time_const = factor(i) * pf(nf) / max((abs(flux_diff(i))*1000. ), 1.e-30)**0.5 / (pe(ne) - pe(nf))


          if (mod(j,6) .eq. 0) then
             dT_old(i) = Tf(i)
          endif

          if (mod(j,6) .eq. 5) then
             if ( abs(Tf(i) - dT_old(i)) .lt. 3.*abs(dT(i))) then
                factor(i) = factor(i)/1.5
             else
                factor(i) = factor(i)*1.1
             endif
          endif
             
       else
          time_const = grav/cpair/(pe(i+1) - pe(i))*del_time
          flux_diff(i) = flux_diff(i) + net_F(i+1) - net_F(i)
       endif


       dT(i) = time_const*flux_diff(i)

        if (dT(i)>5.0_dp) then
           dT(i) = 5.0_dp
        endif
        if (dT(i)<-5.0_dp) then
           dT(i) = -5.0_dp
        endif
        Tf(i) = Tf(i) + dT(i)

        if (Tf(i) < 1._dp) then
           Tf(i) = 1._dp
           factor(i) = 1.
        else if (Tf(i) .gt. 10000. ) then
           Tf(i) = 10000._dp
        end if

     end do

     fsmooth_sum = 0.0
     do i=1,nf
        do m=nf,i,-1
           fsmooth_sum(i) = fsmooth_sum(i) + fsmooth(m)
        enddo
     enddo
     
     
     if ((surface) .and. sensible_heat) then
        dT_surf = 0.
        ! Do turbulent heat exchange between surface and lowest atmospheric layer
        Sens = cpair*pf(nf)/rdgas/Tf(nf) * C_d * U * (Tf(nf) - Ts)
        !Sens = 0
        dT_surf = ( s_dn(ne)*(1-A_s) - fup(ne) + fdn(ne)  + Sens)/cp_s/rho_s/depth * del_time
        Tf(nf) = Tf(nf)  - Sens*grav/cpair/(pe(ne) - pe(ne-1))*del_time
        Ts = Ts + dT_surf
        call dew_point_T(pe(ne), dew_pt)
        if (Ts .gt. dew_pt) then
           Ts = dew_pt
        endif
     else
        Sens = 0.0
     endif


     if (conv_switch) then
        if (moisture_scheme == 'surface') then
           call calc_q_and_grad(pf, delp, Tf, q, dry_mask, olr, grad, ktrop)
           call new_adjust(pf, delp, Tf, q, ktrop, grad, olr+s_up(1), dry_mask)
        endif
     endif
     

     do i = 2, nf-1
       call bezier_interp(pf(i-1:i+1), Tf_half(i-1:i+1), 3, pe(i), Te(i))
    end do

    call bezier_interp(pf(nf-2:nf), Tf_half(nf-2:nf), 3, pe(nf), Te(nf))

    call linear_log_interp(pe(1), pf(1), pf(2), Tf(1), Tf(2), Te(1))
    call linear_log_interp(pe(ne), pf(nf-1), pf(nf), Tf(nf-1), Tf(nf), Te(ne) )

    if (surface)  then
       Te(ne) = Ts
    else
       Ts = Te(ne)
    endif
    
    
    if (.true.) then
       dflux = 0.

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
       enddo



    if (abs((net_F(1) - Fint)/Finc) .lt. 1.e-5 .and. maxval(abs(dflux))< 1.e-5) then
       if (.not. conv_switch) then
          write(*,*) 'CONVECTION STARTING'
          conv_switch = .true.
          conv_start_step = j
       else if (.not. sensible_heat .and. j - conv_start_step .gt. 1000) then
          write(*,*) 'SENSIBLE HEAT TRANSFER STARTING'
          sensible_heat = .true.
          accelerate = .false.
          conv_start_step = j
       else if (sensible_heat .and. j- conv_start_step .gt. 1000) then
          write(*,*) 'CONVERGED'
          exit
       endif
    endif
    

if (mod(j,100) .eq. 0) then
!   do i=1,nf
!      write(*,*) factor(i)
!   enddo
   write(*,*) '------------------------'
   write(*,*) factor_surf, abs(Ts - t_surf_old), 3.*abs(dT_surf)
   write(*,*) dry_mask
   write(*,*) '-------------------------'

    do i=nf,1,-1
       if (dry_mask(i)) then
          test = i
       endif
    enddo
    !write(*,*) dry_mask
    write(*,*) j, 'max dflux', maxval(abs(dflux)), maxloc(abs(dflux)), maxval(abs(dT)), test, (net_F(1) - Fint)/Finc, Tf(nf), &
         net_F(1) - Fint, Ts
    write(*,*) j, s_dn(1), s_up(1), fdn(1), fup(1), net_F(1)
    write(*,*) j, s_dn(ne), s_up(ne), fdn(ne), fup(ne), net_F(ne), df_surf
    write(*,*) sb*Ts**4 - fup(ne)

    !write(*,*) dry_mask
    !write(*,*) q
    !    print*, 'near end', log(Tf(151)/Tf(150))/log(pf(151)/pf(150))
    call dump_data(file_name, nf,ne,Tf,pf,pe,olr,Finc,Fint,Te,q,fup,fdn,s_dn,s_up,Ts)
 endif
endif

enddo ! do j=1,Nt

 write(*,*) 'Surface temperature', Te(ne)


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

call dump_data(file_name, nf,ne,Tf,pf,pe,olr,Finc,Fint,Te,q,fup,fdn,s_dn,s_up,Ts)
!call calc_q_and_grad(pf, delp, Tf, q, dry_mask, olr, grad, ktrop)
!write(*,*) q
write(*,*) dry_mask
  end subroutine step
  
end module timestep
