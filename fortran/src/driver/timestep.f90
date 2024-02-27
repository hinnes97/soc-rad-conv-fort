module timestep
  use params, only: dp, Nt, nf, ne, const, Finc, Fint, q0, surface, A_s, sb, surf_const, &
       moisture_scheme, sb, grav, cpair, del_time, accelerate, cp_s, rho_s, rdgas, &
       depth, U, C_d, conv_switch, sensible_heat
  use flux_mod, only: get_fluxes
  use utils, only: linear_log_interp, bezier_interp
  use convection, only: dry_adjust
  use condense, only : rain_out, cold_trap,q_sat, dew_point_T, set_q
  use adjust_mod, only : new_adjust
  use supercrit_adjust, only: adjustment
  use atmosphere, only : nqt, nqr
  use io, only: dump_data
  implicit none

  ! Variables needed across subroutines across timesteps
  real(dp) :: g_conv_old, dflux_surf_old
  real(dp) :: ts_conv_old, ts_old
  
contains
  subroutine step(Tf, pf, pe, file_name, q, Ts)
    !====================================================================================
    ! Description
    !====================================================================================
    !
    ! Controls timestepping. Single timestep controlled from single_step
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
    real(dp), dimension(:,:), intent(inout) :: q   ! Specific humidity
    
    real(dp), intent(inout) :: Ts ! Surface temperature

    !====================================================================================
    ! Work variables
    !====================================================================================
    
    integer :: i,j
    integer :: ktrop 
    integer :: output_interval
    integer :: print_interval
    
    logical :: stop_switch = .false.   ! Controls when exit due to convergence occurs

    real(dp), dimension(ne) :: net_F  ! Net flux
    real(dp), dimension(ne) :: fup    ! Upwards IR flux
    real(dp), dimension(ne) :: fdn    ! Downwards IR flux
    real(dp), dimension(ne) :: s_up   ! Upwards SW flux
    real(dp), dimension(ne) :: s_dn   ! Downwards SW flux
    real(dp), dimension(ne) :: Te     ! Edge of layer temperature
    
    real(dp), dimension(nf) :: dflux ! Local convergence criterion
    real(dp), dimension(nf) :: dT    ! Change in temperature
    real(dp), dimension(nf) :: delp
    real(dp), dimension(nf) :: factor
    real(dp), dimension(nf) :: T_old
    real(dp) :: factor_surf, t_surf_old, dT_surf

    logical, dimension(nf) :: dry_mask
    
    real(dp) ::  sens ! Sensible heat flux
    real(dp) :: olr
    real(dp) :: dflux_surf
    
    !====================================================================================
    ! Main Body
    !====================================================================================

    ! TODO: move to namelist
    print_interval = 1
    output_interval  = 100
    
    ! Initialise some of the variables
    factor = -50. ! Arbitary negative number
    factor_surf = -20. ! Arbitrary negative number
    dT = 0._dp
    dt_surf = 0._dp
    dry_mask = .false.
    g_conv_old = -33.

    ! Create delp array
    do i=1,nf
       delp(i) = pe(i+1) - pe(i)
    enddo


    if (moisture_scheme /= 'none' .and. moisture_scheme /= 'supercrit') call set_q(pf,Tf, q,ktrop)

    call get_fluxes(nf, ne, Tf, pf, Te, pe, &
          net_F, 1.0_dp, Finc, Fint, olr, q, Ts, fup, fdn, s_dn, s_up)

    call print_header
    
    do j =1,Nt

       ! Adapative timestepping
       if (accelerate) then
          call alter_timestep(factor, j, T_old, Tf, dT, factor_surf, t_surf_old, Ts, dt_surf)
       endif

       ! Advance forward one timestep
       call single_step(pf, pe, delp, Tf,Te, q,Ts, dT, net_F, fup, fdn, &
            s_up, s_dn, olr, dry_mask,  accelerate,&
            sens, j, factor, factor_surf, dT_surf, ktrop, dflux_surf)

       ! Check to see whether to enter next phase of model run
       call check_convergence(net_F, sens, Tf, pe, dry_mask,&
            stop_switch, j, dflux, Ts, dflux_surf)

       ! Output data to screen
       if (mod(j,print_interval) .eq. 0) then
          call print_update(j, dflux, net_F, dT, Ts, dry_mask, dflux_surf)
       endif

       ! Output data to file
       if (mod(j, output_interval) .eq. 0) then
          call dump_data(file_name, nf,ne,Tf,pf,pe,olr,Finc,Fint,&
                         Te,q,fup,fdn,s_dn,s_up,Ts, dry_mask)
       endif

       ! Exit if check_convergence says so
       if (stop_switch) exit
       
    enddo 

    call print_final_status(j, pf, Tf, Ts, q, dry_mask)
    
    call dump_data(file_name, nf,ne,Tf,pf,pe,olr,Finc,Fint,&
                   Te,q,fup,fdn,s_dn,s_up,Ts, dry_mask)

  end subroutine step

  subroutine single_step(pf, pe, delp, Tf, Te,q, Ts, dT, net_F, fup, fdn,&
                         s_up, s_dn, olr, dry_mask, accelerate, &
                         sens, tstep, factor, factor_surf, dt_surf, ktrop, dflux_surf)
    
    !====================================================================================
    ! Description
    !====================================================================================
    !
    ! Steps once forwards
    !
    !====================================================================================
    ! Input variables
    !====================================================================================
    
    real(dp), dimension(:), intent(in) :: pf ! Layer pressure
    real(dp), dimension(:), intent(in) :: pe ! Level pressure
    real(dp), dimension(:), intent(in) :: delp 

    logical, intent(in) :: accelerate    ! Controls adaptive timestepping

    integer, intent(in) :: tstep ! Current timestep
    integer, intent(inout) :: ktrop

    !====================================================================================
    ! Input/Output variables
    !====================================================================================
    
    real(dp), dimension(:), intent(inout) :: Tf  ! Layer temperatures
    real(dp), dimension(:), intent(inout) :: Te  ! Edge of layer temperature
    real(dp), dimension(:,:), intent(inout) :: q   ! Specific humidity
    
    real(dp), intent(inout) :: Ts ! Surface temperature
    real(dp) :: sens        ! Sensible heat flux at surface

    real(dp), dimension(:), intent(inout) :: net_F  ! Net flux
    real(dp), dimension(:), intent(inout) :: fup    ! Upwards IR flux
    real(dp), dimension(:), intent(inout) :: fdn    ! Downwards IR flux
    real(dp), dimension(:), intent(inout) :: s_up   ! Upwards SW flux
    real(dp), dimension(:), intent(inout) :: s_dn   ! Downwards SW flux

    real(dp), dimension(:), intent(inout) :: dT        ! Current timestep's change in temp
    real(dp), dimension(:), intent(inout) :: factor       
    real(dp), intent(inout) :: factor_surf
    logical, dimension(:), intent(inout) :: dry_mask ! Whether convection in layer

    real(dp), intent(inout) :: dt_surf

    real(dp), intent(inout) :: dflux_surf ! Surface flux convergence
    !====================================================================================
    ! Work variables
    !====================================================================================
    integer :: i
    
    real(dp), dimension(size(pf)) :: Tf_half   ! Temperature at half timestep
    real(dp), dimension(size(pf)) :: grad      ! Adiabatic lapse rate dlnT/dlnp
    real(dp), dimension(size(pf)) :: flux_diff ! Divergence of net flux
    
    real(dp) :: time_const  ! Multiplies net flux difference to give dT
    real(dp) :: df_surf     ! Net flux at surface
    real(dp) :: olr         ! Outgoing longwave radiation
    real(dp) :: dew_pt      ! Dew point temperature

    real(dp) :: minmax_dT, min_T, max_T

    minmax_dT = 5.0_dp
    min_T = 1.0_dp
    max_T = 10000._dp
    
    ! Interpolate Tf to edges for radiation scheme
    call interp_to_edges(pf, pe, Tf, Te)
    if (surface)  then
       Te(ne) = Ts
    else
       Ts = Te(ne)
    endif

    ! Calls radiation driver
    call get_fluxes(nf, ne, Tf, pf, Te, pe, &
         net_F, 1.0_dp, Finc, Fint, olr, q, Ts, fup, fdn, s_dn, s_up)
    do i=1,ne
       write(*,*) pe(i), Te(i), fup(i), fdn(i), s_dn(i), s_up(i)
    enddo
    stop


    ! Step half forwards
       do i=1,nf
          if (accelerate) then
             time_const = factor(i) * pf(nf) / max((abs(net_F(i+1) - net_F(i))*1000. ), 1.e-30)**0.9 / (pe(ne) - pe(nf))
          else
             time_const = grav/cpair/(pe(i+1) - pe(i))*del_time
          endif ! if accelerate

          dT(i) = time_const*(net_F(i+1) - net_F(i))

          ! Limit temperature increases to 5 degrees
          dT(i) = min(max(dT(i), -minmax_dT), minmax_dT)
          
       end do ! i=1,nf

       Tf_half = Tf + 0.5_dp*dT

       ! Interpolate this to edges
       call interp_to_edges(pf, pe, Tf_half, Te)

       if (surface) then
          Te(ne) = Ts
       else
          Ts = Te(ne)
       endif
    
       ! Call radiation scheme at half timestep
       call get_fluxes(nf, ne, Tf_half, pf, Te, pe,  &
            net_F, 1.0_dp, Finc, Fint, olr, q, Ts, fup, fdn, s_dn, s_up)

    
       if (.not. surface) then 
          net_F(ne) = Fint
       endif

    
       flux_diff(:) = 0.
       sens = 0.
       do i=1,nf
          if (accelerate) then
             ! Accelerated timestepping
             
             flux_diff(i) = flux_diff(i) + net_F(i+1) - net_F(i)

             if ((i .eq. nf) .and. surface) then
             
                ! Do turbulent heat exchange between surface and lowest atmospheric layer
                Sens = cpair*pf(nf)/rdgas/Tf(nf) * C_d * U * (Tf(nf) - Ts)
                
                df_surf = ( s_dn(ne)*(1-A_s) - fup(ne) + fdn(ne)  + Sens)
                dflux_surf = df_surf/sb/Ts**4
                
                flux_diff(i) = flux_diff(i)  - Sens
                dT_surf = factor_surf * pf(nf)/(abs(df_surf)*1000.)**0.9/(pe(ne) - pe(nf)) * (df_surf)

                ! Step forwards in time
                Ts = Ts + dT_surf

                ! If surface is liquid water, then temperature cannot exceed temperature
! where psat(T)>p_s.
                write(*,*) 'inside here', moisture_scheme, surface
                call dew_point_T(pe(ne), dew_pt)
                if (Ts .gt. dew_pt) then
                   Ts = dew_pt
                endif

             endif
                  


             ! Require flux_diff not be zero
             time_const = factor(i) * pf(nf) / max((abs(flux_diff(i))*1000. ), 1.e-30)**0.9 / (pe(ne) - pe(nf))


             
          else
             ! Real, physicsal timestep here
             time_const = grav/cpair/(pe(i+1) - pe(i))*del_time
             flux_diff(i) = flux_diff(i) + net_F(i+1) - net_F(i)

             if ((i.eq. nf) .and. surface) then
                Sens = cpair*pf(nf)/rdgas/Tf(nf) * C_d * U * (Tf(nf) - Ts)
                
                df_surf = ( s_dn(ne)*(1-A_s) - fup(ne) + fdn(ne)  + Sens)
                dflux_surf = df_surf/sb/Ts**4
                
                flux_diff(i) = flux_diff(i)  - Sens

                dT_surf = df_surf/depth/rho_s/cp_s * del_time

                ! Step forwards in time
                Ts = Ts + dT_surf

                ! If surface is liquid water, then temperature cannot exceed temperature
! where psat(T)>p_s.
                write(*,*) 'inside here 2', moisture_scheme, surface
                call dew_point_T(pe(ne), dew_pt)
                if (Ts .gt. dew_pt) then
                   Ts = dew_pt
                endif

             endif
             
          endif
          


          dT(i) = time_const*flux_diff(i)

          ! Again limit temperature changes
          dT(i) = min(max(dT(i), -minmax_dT), minmax_dT)

          ! Step forward full
          Tf(i) = Tf(i) + dT(i)

          ! Stop temperature going too low or high
          Tf(i) = max(min(max_T, Tf(i)), min_T)

       end do !i=1,nf

       ! Convective adjustment happens here
       if (conv_switch) then
          if (moisture_scheme /= 'none' .and. moisture_scheme /= 'supercrit') then
             write(*,*) 'before new_adjust'
             write(*,*) moisture_scheme, surface
             call set_q(pf,Tf, q, ktrop)
             call new_adjust(pf, delp, Tf, q, ktrop, grad, olr+s_up(1), dry_mask, tstep)
          else if (moisture_scheme == 'supercrit' ) then
             write(*,*) 'before supercrit adjust'
             !call adjustment(pf, delp, Tf, q, ktrop, grad, olr+s_up(1), dry_mask, tstep)
          endif
          
       endif

! Ensure final state is at saturation
       write(*,*) 'helloworld', moisture_scheme, surface
       if (moisture_scheme /='none' .and. moisture_scheme /='supercrit') call set_q(pf,Tf, q, ktrop)
       

       ! Interpolate convective adjustment to cell edges
       call interp_to_edges(pf, pe, Tf, Te)

       if (surface)  then
          Te(ne) = Ts
       else
          Ts = Te(ne)
       endif
       write(*,*) 'end of single step------------------------------'
     end subroutine single_step

     subroutine check_convergence(net_F, sens, Tf, pe,dry_mask, &
          stop_switch, tstep, dflux, Ts,dflux_surf)
       !====================================================================================
       ! Input variables
       !====================================================================================
       real(dp), intent(in) :: net_F(ne) ! Net fluxes
       
       real(dp), intent(in) :: sens, dflux_surf      ! Sensible heat

       real(dp), intent(in) :: pe(ne)
       
       logical, intent(in) :: dry_mask(nf)

       integer, intent(in) :: tstep

       !====================================================================================
       ! Input/Output Variables
       !====================================================================================
       
       real(dp), intent(inout) :: dflux(:)
       real(dp), intent(inout) :: Tf(nf)    ! Temperature

       real(dp), intent(inout) :: Ts
       !====================================================================================
       ! Input/Output Variables
       !====================================================================================

       logical, intent(out) :: stop_switch
       
       !====================================================================================
       ! Work variables
       !====================================================================================
       
       integer :: i
       real(dp) :: g_conv, dew_pt_surf
       
       dflux = 0.0
       
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

       if (abs((net_F(1) - Fint)/Finc) .lt. 1.e-5 .and. maxval(abs(dflux))< 1.e-5 .and. abs(dflux_surf) .lt. 1.e-5) then
!!$          if (.not. conv_switch) then
!!$             conv_switch = .true.
!!$             conv_start_step = tstep
!!$          else if ((.not. sensible_heat .and. tstep  - conv_start_step .gt. 100)) then
!!$             !write(*,*) 'Global and other flux satisfied'
!!$             kink = .false.
!!$             counter = 0
!!$             !do i=1,nf-3
!!$             !   if ( (Tf(i+1) - Tf(i))*(Tf(i+2)- Tf(i+1))*(Tf(i+3)-Tf(i+2)) .lt. 0) then
!!$             !      kink = .true.
!!$             !   endif
!!$             !enddo
!!$             
!!$             if (kink .and. counter .lt. 2) then
!!$                counter = counter + 1
!!$                call kink_smoother(pf,Tf)
!!$                call interp_to_edges(pf,pe,Tf,Te)
!!$             else
                !sensible_heat = .true.
                !accelerate = .false.
                !conv_start_step = tstep
                !stop_switch = .false.
                !write(*,*) 'MOVING ON', tstep

          if (mod(tstep,1000) .eq. 0) then
             ts_old = Ts
             write(*,*) 'RECORDING TS', ts_old
          else if (mod(tstep,1000) .eq. 999 ) then
             write(*,*) 'CHECKING TS', ts_old, Ts, abs(Ts - ts_old)/ts_old
             if ( abs(Ts - ts_old)/ts_old .lt. 1.e-5) then
                !sensible_heat = .true.
                !accelerate = .false.
                stop_switch = .true.
                !conv_start_step = tstep
             endif
          endif
                
           !endif

!!$          else if (sensible_heat .and. tstep - conv_start_step .gt. 1000 .and. abs(dflux_surf).lt. 1.e-5) then
!!$             stop_switch = .true.
!!$          endif
       endif

       ! Sometimes gets stuck in convection loop without properly converging due to kink.
       ! In this case, smooth temperature profile to aid convergence

       ! COMMENT FOR TESTING
!!$       if (tstep .eq. 2000 .and. .not. sensible_heat) then
!!$          ! think of what todo e
!!$          call kink_smoother(pf, Tf)
!!$          call interp_to_edges(pf,pe,Tf,Te)
!!$       endif
       
!!$       ! See what the max dT is and stop if at 5, restart with smaller timestep
!!$       if (sensible_heat .and. tstep-conv_start_step .gt. 100) then
!!$          if (mod(tstep,500) .eq. 0) then
!!$             ts_conv_old = maxval(abs(dT))
!!$          else if (mod(tstep,500) .eq. 499 ) then
!!$             if (ts_conv_old .gt. 4.999 .and. maxval(abs(dT)) .gt. 4.999) then
!!$                write(*,*) 'Program not converging, restart with smaller timestep'
!!$                stop 99
!!$             endif
!!$          endif
!!$       endif

       ! Check to see if the surface q is 1 and global convergence stuck
       if (sensible_heat .and. tstep-1000 .gt. 100) then
          ! Check if global convergence same as 500 ago
          if (mod(tstep,500) .eq. 0) then
             dflux_surf_old = dflux_surf
             g_conv_old = abs((net_F(1) - Fint)/Finc)
          else if (mod(tstep, 500) .eq. 499) then
             g_conv = abs((net_F(1) - Fint)/Finc)
             call dew_point_T(pe(ne), dew_pt_surf)
             if (abs((dflux_surf - dflux_surf_old)/dflux_surf_old) .lt. 1.e-4 .and. abs(Ts - dew_pt_surf) .lt. 1.e-5) then
                write(*,*) 'Surface now at q = 1, retry with lower Finc increment'
                stop 100
             endif
          endif
       endif

       if (tstep .gt. 10000 .and. abs(Ts - dew_pt_surf) .lt. 1.e-5) then
          stop 100
       endif

!!$       !Commenting for testing
!!$       if ( .not. sensible_heat .and. tstep - conv_start_step .gt. 100) then
!!$          ! Check to see if surface is stuck at dew point
!!$          if (mod(tstep,1000) .eq. 0) then
!!$                g_conv_old = Ts
!!$          else if (mod(tstep, 1000) .eq. 999) then
!!$             call dew_point_T(pe(ne), dew_pt_surf)
!!$             if (abs(g_conv_old - dew_pt_surf) .lt. 1.e-5 .and. abs(Ts - dew_pt_surf) .lt. 1.e-5) then
!!$                write(*,*) 'Surface now at q = 1, retry with pure timestepping'
!!$                stop 101
!!$             endif
!!$          endif
!!$       endif
       
       
     end subroutine check_convergence

     subroutine kink_smoother(pf, Tf)
       real(dp), intent(in) :: pf(nf)
       real(dp), intent(inout) :: Tf(nf)

       integer :: i,j

       do j=1,100
          do i=1,nf-3
             if ( (Tf(i+1) - Tf(i))*(Tf(i+2)- Tf(i+1))*(Tf(i+3)-Tf(i+2)) .lt. 0) then
                !Kink? Interpolate in between points to end points 
                call linear_log_interp(pf(i+1), pf(i), pf(i+3), Tf(i), Tf(i+3), Tf(i+1))
                call linear_log_interp(pf(i+2), pf(i), pf(i+3), Tf(i), Tf(i+3), Tf(i+2))
             endif
          enddo
       enddo
       
     end subroutine kink_smoother
     
  subroutine interp_to_edges(pf, pe, Tf, Te)
    !====================================================================================
    ! Description
    !====================================================================================
    !
    ! Interpolate temperature to edge of layers
    !
    !====================================================================================
    ! Input Variables
    !====================================================================================
    
    real(dp), intent(in) :: pf(nf) ! Mid layer pressure
    real(dp), intent(in) :: pe(ne) ! Edge level pressure
    real(dp), intent(in) :: Tf(nf) ! Mid layer temperature
    
    !====================================================================================
    ! Output Variables
    !====================================================================================

    real(dp), intent(out) :: Te(ne) ! Edge of layer temperature

    !====================================================================================
    ! Work variables
    !====================================================================================

    integer :: i

    !====================================================================================
    ! Main body
    !====================================================================================

    do i = 2, nf-1
       call bezier_interp(pf(i-1:i+1), Tf(i-1:i+1), 3, pe(i), Te(i))
    end do

    call bezier_interp(pf(nf-2:nf), Tf(nf-2:nf), 3, pe(nf), Te(nf))
    call linear_log_interp(pe(1), pf(1), pf(2), Tf(1), Tf(2), Te(1))
    call linear_log_interp(pe(ne), pf(nf-1), pf(nf), Tf(nf-1), Tf(nf), Te(ne) )

  end subroutine interp_to_edges

  subroutine print_header
    write(*,*) '==========================================================================================================='
    write(*,'(7A15)') ' ', 'Local', 'Global', ' ', ' ', 'Max dsurf', ''
    write(*,'(7A15)') 'Timestep', 'Convergence', 'Convergence', 'Max dT [K]', 'Ts [K]', 'maxloc(dflux)'
    write(*,*) '-----------------------------------------------------------------------------------------------------------'
  end subroutine print_header
  
  subroutine print_update(tstep, dflux, net_F, dT, Ts, dry_mask,  &
                         dflux_surf)
    !====================================================================================
    ! Description
    !====================================================================================
    !
    ! Print update of experiment to stdout
    !
    !====================================================================================
    ! Input variables
    !====================================================================================

    integer, intent(in) :: tstep ! Iteration number

    real(dp), intent(in) :: dflux(:) ! Convergence tester
    real(dp), intent(in) :: net_F(:) ! Net flux
    real(dp), intent(in) :: dT(:)    ! Temperature change

    logical, intent(in) :: dry_mask(:)
    
    real(dp), intent(in) :: Ts ! Surface temperature
    real(dp), intent(in) :: dflux_surf
    
    !====================================================================================
    ! Work variables
    !====================================================================================

    integer :: i,imax
    
    character(len=40) :: fmt
    
    !====================================================================================
    ! Main body
    !====================================================================================

    fmt = '(I15,3(1PE15.3),0PF15.3,1PE10.3,I5)'
    imax=nf
    do i=nf,1,-1
       if (dry_mask(i)) imax = i
    enddo
    write(*,fmt) tstep, maxval(abs(dflux)), (net_F(1) - Fint)/Finc, maxval(abs(dT)), Ts, &
         abs(dflux_surf), maxloc(abs(dflux))

  end subroutine print_update

  subroutine print_final_status(tstep, pf, Tf,Ts, q, dry_mask)
    !====================================================================================
    ! Description
    !====================================================================================
    !
    ! Print final status of experiment to stdout
    !
    !====================================================================================
    ! Input variables
    !====================================================================================

    integer, intent(in) :: tstep ! Iteration number

    real(dp), intent(in) :: Ts ! Surface temperature

    real(dp), intent(in) :: Tf(nf)
    real(dp), intent(in) :: q(nf,nqt)
    real(dp), intent(in) :: pf(nf)
    
    logical, intent(in) :: dry_mask(nf)
    
    !====================================================================================
    ! Work variables
    !====================================================================================

    integer :: i
    
    !====================================================================================
    ! Main body
    !====================================================================================

    write(*,*) '===================================================================================='
    write(*,'(2x,A)') 'CONVERGED'
    write(*,*) '------------------------------------------------------------------------------------'
    write(*,'(2x,A, I6)') 'Timestep: ', tstep
    write(*,*) '------------------------------------------------------------------------------------'
    write(*,'(4A17)') 'Pressure', 'Temperature', 'Specific', 'Convection'
    write(*,'(4A17)') '[Pa]', '[K]', 'Humidity', 'Mask'
    write(*,*) '------------------------------------------------------------------------------------'
    do i = 1,nf
       write(*,'(1PE17.3, 0PF17.3, 1PE17.3, L17)') pf(i), Tf(i), q(i,1), dry_mask(i)
    enddo
    write(*,*) '------------------------------------------------------------------------------------'
    write(*,'(2x,A,F10.3)') 'Surface Temperature [K]:', Ts 
    write(*,*) '===================================================================================='
  end subroutine print_final_status

  subroutine alter_timestep(factor, tstep, dT_old, Tf, dT, &
                            factor_surf, t_surf_old, Ts, dT_surf)
    !==========================================================================
    ! Description
    !==========================================================================    
    ! Uses adaptive timestepping (Malek 2017) to vary the timestep based
    ! on whether oscillations in the temperature are detected
    !==========================================================================

    !==========================================================================
    ! Input variables
    !==========================================================================    

    real(dp), intent(in) :: Tf(:), dT(:)
    real(dp), intent(in) :: Ts, dT_surf
    integer, intent(in) :: tstep

    !==========================================================================
    ! Mixed input/output
    !==========================================================================

    real(dp), intent(inout) :: factor(:)   ! Multiplies timestep for acceleration
    real(dp), intent(inout) :: factor_surf ! Separate factor for the surface
    real(dp), intent(inout) :: dT_old(:)   ! Memory of old dT for oscillation detection
    real(dp), intent(inout) :: t_surf_old  ! Same as above but for surface

    integer :: i

    ! Adaptive timestep for atmosphere
    do i=1,nf
       if (factor(i) .lt. 0) then
          factor(i) = 1.e-10
       endif
    
       if (factor_surf .lt. 0.) factor_surf = 1.e-10

       ! Set 'old' temperatures
       if (mod(tstep,6) .eq. 0) then
          t_surf_old = ts
          dT_old(i) = Tf(i)
       endif

       ! Compare with new temperatures
       if (mod(tstep,6) .eq. 5) then
          ! Condition for oscillation
          if (abs(Tf(i) - dT_old(i)) .lt. 3.*abs(dT(i))) then
             ! Reduce timestep
             factor(i) = factor(i)/1.5
          else
             ! Increase timestep
             factor(i) = factor(i) * 1.1
          endif
       endif
       
    
       if (factor(i) .lt. 1.e-10) factor(i) = 1.e-10
    enddo

    ! Same for surface (is this code redundant?)
    if (mod(tstep,6) .eq. 0) then
       t_surf_old = ts
    endif
    if (mod(tstep,6) .eq. 5) then
       if ( abs(Ts - t_surf_old) .lt. 3.*abs(dT_surf)) then
          factor_surf = factor_surf/1.5
       else
          factor_surf = factor_surf*1.1
       endif
       if (factor_surf .lt. 1.e-10) factor_surf = 1.e-10
    endif

  end subroutine alter_timestep
  
end module timestep
