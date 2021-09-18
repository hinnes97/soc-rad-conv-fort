module utils
  use, intrinsic :: iso_fortran_env
  implicit none

  integer, parameter :: dp = real64

  private :: dp
  
contains
  
    subroutine linspace(start, end, arr)
    real(dp), intent(in) :: start
    real(dp), intent(in) :: end

    real(dp), dimension(:), intent(out) :: arr

    integer :: n 
    integer :: i

    n = size(arr)
    if (n==0) return
    if (n==1) then
       arr(1) = start
    end if

    do i=1,n
       arr(i) = start + (end - start) * (i-1)/(n-1)
    end do        
    
  end subroutine linspace

  subroutine logspace(start, end, arr)
    real(dp), intent(in) :: start, end
    real(dp), intent(inout)  :: arr(:)

    integer :: n
    integer :: i

    n = size(arr)
    call linspace(start, end, arr)

    do i=1,n
       arr(i) = 10**(arr(i))
    end do
    
  end subroutine logspace

end module utils
