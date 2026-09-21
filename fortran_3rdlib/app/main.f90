program main
  use mvn_common, only: dp, build_problem, checksum
  use stdlib_mvn, only: sample_stdlib
  implicit none

  integer, parameter :: n = 4
  real(dp), allocatable :: mu(:), sigma(:, :), x(:)
  integer :: i

  call build_problem(n, mu, sigma)
  allocate(x(n))

  write (*, '(a,*(f10.6,1x))') 'mu     :', mu

  call sample_stdlib(mu, sigma, x)
  write (*, '(a,*(f10.6,1x))') 'stdlib :', x
  write (*, '(a,f12.6)') '  checksum = ', checksum(x)

  write (*, '(a)') ''
  write (*, '(a)') 'input Sigma:'
  do i = 1, n
    write (*, '(*(f8.4,1x))') sigma(i, :)
  end do
end program main
