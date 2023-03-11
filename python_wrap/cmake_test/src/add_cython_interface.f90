module add_cython_interface
  use add_cython, only: add
  use iso_c_binding, only: c_int
  
  implicit none

contains
  
  subroutine c_add(a,b,c) bind(c)

    integer(c_int), intent(in) :: a
    integer(c_int), intent(in) :: b
    integer(c_int), intent(out) :: c

    call add(a,b,c)
    
  end subroutine c_add
  
  
end module add_cython_interface
