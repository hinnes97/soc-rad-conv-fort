module io
  use netcdf
  use iso_fortran_env
  
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

  subroutine file_setup(path, nf, ne, ncid)
    character (len=*), intent(in) :: path
    integer, intent(in) :: nf
    integer, intent(in) :: ne
    integer, intent(out) :: ncid

    integer :: status
    integer :: pf_dim_id, pe_dim_id
    integer :: tf_id, pf_id, pe_id, olr_id, tau_ir_inf_id, tau_v_inf_id
    integer :: finc_id, fint_id

    ! Create output file
    status = nf90_create(path, 0, ncid)
    if (status /= nf90_noerr) call handle_err(status)

    ! Create dimensions
    status = nf90_def_dim(ncid, "pfull", nf, pf_dim_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_def_dim(ncid, "pedge", ne, pe_dim_id)
    if (status /= nf90_noerr) call handle_err(status)

    ! Create variables
    status = nf90_def_var(ncid, "Tf", nf90_double, pf_dim_id, tf_id)
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

    ! Create attributes
    status = nf90_put_att(ncid, pf_id, "Units", "Pa")
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, pe_id, "Units", "Pa")
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, tf_id, "Units", "K")
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
    
    status = nf90_put_att(ncid, pf_id, "Long name", "Mid level pressure")
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, pe_id, "Long name", "Interface pressure")
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, tf_id, "Long name", "Mid level temperature")
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

    ! End file definition
    status = nf90_enddef(ncid)
    if (status /= nf90_noerr) call handle_err(status)
    
  end subroutine file_setup

  subroutine dump_data(ncid, nf, ne, Tf, pf, pe, olr, tau_IR_inf, tau_V_inf, Finc, Fint)
    integer, intent(in) :: ncid
    integer, intent(in) :: nf, ne
    real(dp), intent(in), dimension(nf) :: Tf,pf
    real(dp), intent(in), dimension(ne) :: pe
    real(dp), intent(in) :: olr, tau_IR_inf, tau_V_inf, Finc, Fint
    
    integer :: dummy_id, status
    
    status = nf90_inq_varid(ncid, "Tf", dummy_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(ncid, dummy_id, Tf)
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

    status = nf90_inq_varid(ncid, "tau_IR_inf", dummy_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(ncid, dummy_id, tau_IR_inf)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_inq_varid(ncid, "tau_V_inf", dummy_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(ncid, dummy_id, tau_V_inf)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_inq_varid(ncid, "Finc", dummy_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(ncid, dummy_id, Finc)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_inq_varid(ncid, "Fint", dummy_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(ncid, dummy_id, Fint)
    if (status /= nf90_noerr) call handle_err(status)
    
  end subroutine dump_data

  subroutine close_file(ncid)
    integer, intent(in) :: ncid

    integer :: status
    
    status = nf90_close(ncid)
    if (status /= nf90_noerr) call handle_err(status)
    
  end subroutine close_file

end module io
