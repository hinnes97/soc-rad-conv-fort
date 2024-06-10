module io
  use netcdf
  use iso_fortran_env
  use params, only: surface
  implicit none

  integer, parameter :: dp = real64

  private :: dp
contains
  
  subroutine handle_err(status)
    integer, intent(in) :: status

    if (status /= nf90_noerr) then
       write(*,*) trim(nf90_strerror(status))
       stop "Stopped"
    end if
    
  end subroutine handle_err

  subroutine file_setup(path, nf, ne, nqr, ncid)
    character (len=*), intent(in) :: path
    integer, intent(in) :: nf
    integer, intent(in) :: ne, nqr
    integer, intent(out) :: ncid

    integer :: status
    integer :: pf_dim_id, pe_dim_id, nqr_dim_id
    integer :: tf_id, pf_id, pe_id, olr_id, tau_ir_inf_id, tau_v_inf_id, te_id, q_id
    integer :: finc_id, fint_id, fdn_id, fup_id, sdn_id, Ts_id, sup_id, dm_id, fturb_id

    ! Create output file
    status = nf90_create(path, 0, ncid)
    if (status /= nf90_noerr) call handle_err(status)

    ! Create dimensions
    status = nf90_def_dim(ncid, "pfull", nf, pf_dim_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_def_dim(ncid, "pedge", ne, pe_dim_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_def_dim(ncid, "nqr", nqr, nqr_dim_id)
    ! Create variables
    status = nf90_def_var(ncid, "Tf", nf90_double, pf_dim_id, tf_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_def_var(ncid, "Te", nf90_double, pe_dim_id, te_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_def_var(ncid, "pfull", nf90_double, pf_dim_id, pf_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_def_var(ncid, "pedge", nf90_double, pe_dim_id, pe_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_def_var(ncid, "olr", nf90_double, olr_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_def_var(ncid, "tau_IR_inf", nf90_double, tau_ir_inf_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_def_var(ncid, "tau_V_inf", nf90_double, tau_v_inf_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_def_var(ncid, "Finc", nf90_double, finc_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_def_var(ncid, "Fint", nf90_double, fint_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_def_var(ncid, "q", nf90_double,(/pf_dim_id,nqr_dim_id/), q_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_def_var(ncid, "fdn", nf90_double,pe_dim_id, fdn_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_def_var(ncid, "fup", nf90_double,pe_dim_id, fup_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_def_var(ncid, "s_dn", nf90_double,pe_dim_id, sdn_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_def_var(ncid, "s_up", nf90_double,pe_dim_id, sup_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_def_var(ncid, "turb_flux", nf90_double, pe_dim_id, fturb_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_def_var(ncid, "Ts", nf90_double, Ts_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_def_var(ncid, "dry_mask", nf90_int,pf_dim_id, dm_id)
    if (status /= nf90_noerr) call handle_err(status)


    ! Create attributes
    status = nf90_put_att(ncid, pf_id, "Units", "Pa")
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, pe_id, "Units", "Pa")
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, tf_id, "Units", "K")
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, te_id, "Units", "K")
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, olr_id, "Units", "W/m^2")
    if (status /= nf90_noerr) call handle_err(status)    
    status = nf90_put_att(ncid, tau_ir_inf_id, "Units", "Dimensionless")
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, tau_v_inf_id, "Units", "Dimensionless")
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, finc_id, "Units", "W/m^2")
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, fint_id, "Units", "W/m^2")
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, fdn_id, "Units", "W/m^2")
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, fup_id, "Units", "W/m^2")
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, sdn_id, "Units", "W/m^2")
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, sup_id, "Units", "W/m^2")
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, Ts_id, "Units", "K")
    if (status /= nf90_noerr) call handle_err(status)


    
    status = nf90_put_att(ncid, pf_id, "Long name", "Mid level pressure")
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, pe_id, "Long name", "Interface pressure")
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, tf_id, "Long name", "Mid level temperature")
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, te_id, "Long name", "Interface temperature")
    if (status /= nf90_noerr) call handle_err(status)    
    status = nf90_put_att(ncid, olr_id, "Long name", "Outgoing longwave radiation")
    if (status /= nf90_noerr) call handle_err(status)    
    status = nf90_put_att(ncid, tau_ir_inf_id, "Long name", "IR optical depth")
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, tau_v_inf_id, "Long name", "SW optical depth")
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, finc_id, "Long name", "Incoming solar radiation")
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, fint_id, "Long name", "Internal heat flux")
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, fdn_id, "Long name", "Downwards LW flux")
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, fup_id, "Long name", "Upwards LW flux")
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, sdn_id, "Long name", "Downwards SW flux")
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, sup_id, "Long name", "Upwards SW flux")

    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, Ts_id, "Long name", "Surface temperature")
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_put_att(ncid, dm_id, "Long name", "Dry convection mask")
    if (status /= nf90_noerr) call handle_err(status)


    ! End file definition
    status = nf90_enddef(ncid)
    if (status /= nf90_noerr) call handle_err(status)
    
  end subroutine file_setup

  subroutine dump_data(file_name, nf, ne, Tf, pf, pe, olr, Finc, Fint, Te, q, fup, fdn, s_dn, &
       s_up, Ts, dry_mask, turb_flux)
    character(*), intent(in) :: file_name
    integer, intent(in) :: nf, ne
    real(dp), intent(in), dimension(nf) :: Tf,pf, q
    real(dp), intent(in), dimension(ne) :: te,pe, fdn, fup, s_dn, s_up, turb_flux
    real(dp), intent(in) :: olr, Finc, Fint, Ts

    integer, intent(in), dimension(nf) :: dry_mask
    integer :: dummy_id, status, k,ncid
    integer :: dry_mask_dummy(nf), i

    status = nf90_open(file_name, nf90_write, ncid)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_inq_varid(ncid, "Tf", dummy_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(ncid, dummy_id, Tf)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_inq_varid(ncid, "Te", dummy_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(ncid, dummy_id, Te)
    if (status /= nf90_noerr) call handle_err(status)


    status = nf90_inq_varid(ncid, "pfull", dummy_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(ncid, dummy_id, pf)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_inq_varid(ncid, "pedge", dummy_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(ncid, dummy_id, pe)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_inq_varid(ncid, "olr", dummy_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(ncid, dummy_id, olr)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_inq_varid(ncid, "Finc", dummy_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(ncid, dummy_id, Finc)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_inq_varid(ncid, "Fint", dummy_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(ncid, dummy_id, Fint)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_inq_varid(ncid, "q", dummy_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(ncid, dummy_id, q)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_inq_varid(ncid, "fdn", dummy_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(ncid, dummy_id, fdn)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_inq_varid(ncid, "fup", dummy_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(ncid, dummy_id, fup)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_inq_varid(ncid, "s_dn", dummy_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(ncid, dummy_id, s_dn)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_inq_varid(ncid, "s_up", dummy_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(ncid, dummy_id, s_up)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_inq_varid(ncid, "turb_flux", dummy_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(ncid, dummy_id, turb_flux)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_inq_varid(ncid, "Ts", dummy_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(ncid, dummy_id, Ts)
    if (status /= nf90_noerr) call handle_err(status)

    ! dry_mask_dummy = 0
    ! do i=1,nf
    !    if (dry_mask(i)) dry_mask_dummy(i) = 1
    ! enddo
    
    status = nf90_inq_varid(ncid, "dry_mask", dummy_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(ncid, dummy_id, dry_mask)
    if (status /= nf90_noerr) call handle_err(status)

    call close_file(ncid)

  end subroutine dump_data

  subroutine read_initial_data(file_name, Tf, Te, q, pf,pe, Ts)
    character(*), intent(in) :: file_name
    real(dp), dimension(:), intent(out) :: Tf, Te, pf, pe
    real(dp), dimension(:,:), intent(out) :: q
    real(dp), optional, intent(out) :: Ts
    integer :: ncid, status, id_tf, id_te, id_q,id_ts, k, id_pf, id_pe

    status = nf90_open(file_name, 0, ncid)
    if (status /= nf90_noerr) call handle_err(status)

    !status = nf90_inq_dimid(ncid, 'pfull', id_pf)
    !if (status /= nf90_noerr) call handle_err(status)
    !status = nf90_inq_dimid(ncid, 'pedge', id_pe)
    !if (status /= nf90_noerr) call handle_err(status)
        
    status = nf90_inq_varid(ncid, 'Tf', id_tf)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_inq_varid(ncid, 'Te', id_te)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_inq_varid(ncid, 'q', id_q)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_inq_varid(ncid, 'Ts', id_ts)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_inq_varid(ncid, 'pfull', id_pf)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_inq_varid(ncid, 'pedge', id_pe)
    if (status /= nf90_noerr) call handle_err(status)
    
    status = nf90_get_var(ncid, id_tf, Tf)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_get_var(ncid, id_te, Te)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_get_var(ncid, id_q, q)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_get_var(ncid, id_ts, Ts)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_get_var(ncid, id_pf, pf)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_get_var(ncid, id_pe, pe)
    if (status /= nf90_noerr) call handle_err(status)

    !do k=2, size(Tf)-1
    !   Tf(k) = 0.25*Tf(k-1) + 0.5*Tf(k) + 0.25*Tf(k+1)
    !enddo

    if (present(Ts)) then
       status = nf90_get_var(ncid, id_ts, Ts)
       if (status /= nf90_noerr) call handle_err(status)
    endif
!    status = nf90_get_var(ncid, id_pf, pf)
!    if (status /= nf90_noerr) call handle_err(status)
!    status = nf90_get_var(ncid, id_pe, pe)
!    if (status /= nf90_noerr) call handle_err(status)

    call close_file(ncid)
    
  end subroutine read_initial_data
  
  subroutine close_file(ncid)
    integer, intent(in) :: ncid

    integer :: status
    
    status = nf90_close(ncid)
    if (status /= nf90_noerr) call handle_err(status)
    
  end subroutine close_file

end module io
