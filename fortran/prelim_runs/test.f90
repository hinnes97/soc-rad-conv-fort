program test
  use netcdf
  integer :: ierr, out
  ierr = nf90_open('init_condition.nc', nf90_write, out)

  if (ierr/=nf90_noerr) then
     write(*,*) trim (nf90_strerror(ierr))
  endif
  
  write(*,*) ierr
end program test


