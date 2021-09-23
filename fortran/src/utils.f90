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

  pure subroutine linear_log_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp) :: lxval, ly1, ly2, lx1, lx2
    real(dp), intent(out) :: yval
    real(dp) :: norm

    lxval = log10(xval)
    lx1 = log10(x1); lx2 = log10(x2)
    ly1 = log10(y1); ly2 = log10(y2)

    norm = 1.0_dp / (lx2 - lx1)

    yval = 10.0_dp**((ly1 * (lx2 - lxval) + ly2 * (lxval - lx1)) * norm)

  end subroutine linear_log_interp


end module utils
