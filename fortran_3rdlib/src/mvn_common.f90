module mvn_common
  implicit none
  private
  public :: dp, build_problem, checksum

  integer, parameter :: dp = kind(1.0d0)

contains

  subroutine build_problem(n, mu, sigma)
    integer, intent(in) :: n
    real(dp), allocatable, intent(out) :: mu(:), sigma(:,:)
    integer :: i, j

    allocate(mu(n), sigma(n, n))
    do i = 1, n
      mu(i) = 0.01_dp * real(i - 1, dp)
    end do
    do j = 1, n
      do i = 1, n
        sigma(i, j) = 0.25_dp ** abs(i - j)
      end do
    end do
  end subroutine build_problem

  pure function checksum(x) result(s)
    real(dp), intent(in) :: x(:)
    real(dp) :: s

    s = sum(x)
  end function checksum

end module mvn_common
