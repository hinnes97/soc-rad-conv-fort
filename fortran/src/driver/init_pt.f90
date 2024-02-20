module init_pt_mod
  use params, only: log_bot_p, log_top_p, p_grid, frac, bot_t, top_t, dp, nqt
  use utils, only: logspace, linear_log_interp
  use atmosphere, only : th_gases, mmw_dry, cp_dry
  implicit none
  
contains

  subroutine init_pt(pf, pe, Tf, Te)

    real(dp), intent(out) :: pf(:), pe(:), Tf(:), Te(:)

    integer :: i, npz

    npz = size(pf)

    select case(p_grid)
     case('log')
        ! Logarithmically spaced pressure grid
        call logspace(log_top_p, log_bot_p, pe)
     case('hires_trop')
        ! Let section of the atmosphere between ps and ps/10 contain a larger proportion of points
        ! Linearly spaced in this region
        call logspace(log10(0.85)+log_bot_p, log_bot_p, pe((npz+1-(npz+1)/frac):(npz+1)))
        call logspace(log_top_p, log10(0.85) + log_bot_p, pe(1:(npz+1-(npz+1)/frac-1)), .false.)
     end select

     ! Initialise pf array from pe
     do i=1,npz
        pf(i) = (pe(i+1) - pe(i)) / (log(pe(i+1)) - log(pe(i)))
     end do

     ! Initialise temperature on a dry adiabat
     do i=1,npz
        Tf(i) = bot_t*(pf(i)/pe(npz+1))**(2./7.)
     enddo

     do i=1,npz+1
        Te(i) = bot_t*(pe(i)/pe(npz+1))**(2./7.)
     enddo

     ! Limit minimum temperature of atmosphere
     do i=1,npz
        Tf(i) = max(Tf(i), top_t)
     enddo

     do i=1,npz+1
        Te(i) = max(Te(i), top_t)
     enddo
   end subroutine init_pt

end module init_pt_mod
