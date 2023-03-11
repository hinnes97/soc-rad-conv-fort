module add_cython
  implicit none
contains
  subroutine add(a,b,c)
    integer, intent(in) :: a
    integer, intent(in) :: b
    integer, intent(out) :: c

    c = a + b
  end subroutine add
  
end module add_cython
