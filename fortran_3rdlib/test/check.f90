program check
  use mvn_common, only: dp, build_problem
  use stdlib_mvn, only: sample_stdlib
  implicit none

  integer, parameter :: n = 4
  integer, parameter :: nsamples = 100000
  real(dp), allocatable :: mu(:), sigma(:, :), x(:)
  real(dp) :: mean_err, cov_err

  call build_problem(n, mu, sigma)
  allocate(x(n))

  call random_seed()
  call accumulate_stdlib(mean_err, cov_err)

  write (*, '(a,es12.4)') 'stdlib max |sample mean - mu|   = ', mean_err
  write (*, '(a,es12.4)') 'stdlib max |sample cov - Sigma| = ', cov_err

  if (mean_err > 0.05_dp) error stop 'sample mean out of tolerance'
  if (cov_err > 0.05_dp) error stop 'sample covariance out of tolerance'

  write (*, '(a)') 'all checks passed'

contains

  subroutine accumulate_stdlib(mean_e, cov_e)
    real(dp), intent(out) :: mean_e, cov_e
    real(dp) :: s(n), c(n, n), d(n)
    integer :: k, a, b

    s = 0.0_dp
    c = 0.0_dp
    do k = 1, nsamples
      call sample_stdlib(mu, sigma, x)
      s = s + x
      do b = 1, n
        do a = 1, n
          c(a, b) = c(a, b) + x(a) * x(b)
        end do
      end do
    end do
    s = s / real(nsamples, dp)
    c = c / real(nsamples, dp)
    do b = 1, n
      do a = 1, n
        c(a, b) = c(a, b) - s(a) * s(b)
      end do
    end do
    d = abs(s - mu)
    mean_e = maxval(d)
    cov_e = maxval(abs(c - sigma))
  end subroutine accumulate_stdlib

end program check
