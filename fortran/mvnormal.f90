module mvnormal_module
  use iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  real(real64), parameter :: symmetry_tolerance = 1.0e-12_real64
  real(real64), parameter :: pi = 4.0_real64 * atan(1.0_real64)

  public :: mvnormal_t
  public :: mvnormal_create

  type, public :: mvnormal_t
    integer, private :: n = 0
    real(real64), allocatable, private :: mean(:)
    real(real64), allocatable, private :: lower(:,:)
  contains
    procedure, public :: sample_into
    procedure, public :: sample
  end type mvnormal_t

contains

  function mvnormal_create(mean, covariance) result(distribution)
    real(real64), intent(in) :: mean(:)
    real(real64), intent(in) :: covariance(:,:)
    type(mvnormal_t) :: distribution
    integer :: n
    integer :: i
    integer :: j
    real(real64) :: scale

    n = size(mean)
    if (n == 0) then
      error stop 'MvNormal: mean must not be empty'
    end if

    if (size(covariance, 1) /= size(covariance, 2)) then
      error stop 'MvNormal: covariance must be square'
    end if
    if (size(covariance, 1) /= n) then
      error stop 'MvNormal: mean and covariance dimensions do not match'
    end if

    scale = max(1.0_real64, maxval(abs(covariance)))
    do i = 1, n
      if (.not. ieee_is_finite(mean(i))) then
        error stop 'MvNormal: mean must contain finite values'
      end if
      do j = 1, n
        if (.not. ieee_is_finite(covariance(i,j))) then
          error stop 'MvNormal: covariance must contain finite values'
        end if
      end do
    end do
    do i = 1, n
      do j = i + 1, n
        if (abs(covariance(i,j) - covariance(j,i)) &
            > symmetry_tolerance * scale) then
          error stop 'MvNormal: covariance must be symmetric'
        end if
      end do
    end do

    distribution%n = n
    allocate(distribution%mean(n), distribution%lower(n,n))
    distribution%mean = mean
    call cholesky_factor(covariance, distribution%lower)
  end function mvnormal_create

  subroutine cholesky_factor(covariance, lower)
    real(real64), intent(in) :: covariance(:,:)
    real(real64), intent(out) :: lower(:,:)
    integer :: n
    integer :: i
    integer :: j
    integer :: k
    real(real64) :: value

    n = size(covariance, 1)
    lower = 0.0_real64

    do i = 1, n
      do j = 1, i
        value = covariance(i,j)
        do k = 1, j - 1
          value = value - lower(i,k) * lower(j,k)
        end do

        if (i == j) then
          if (value <= 0.0_real64) then
            error stop 'MvNormal: covariance must be positive definite'
          end if
          lower(i,j) = sqrt(value)
        else
          if (lower(j,j) <= 0.0_real64) then
            error stop 'MvNormal: covariance must be positive definite'
          end if
          lower(i,j) = value / lower(j,j)
        end if
      end do
    end do
  end subroutine cholesky_factor

  subroutine sample_into(self, scratch, output)
    class(mvnormal_t), intent(in) :: self
    real(real64), intent(inout) :: scratch(:)
    real(real64), intent(out) :: output(:)
    integer :: i
    integer :: j
    real(real64) :: z1
    real(real64) :: z2

    if (self%n == 0) then
      error stop 'MvNormal: uninitialized distribution'
    end if
    if (size(scratch) /= self%n) then
      error stop 'MvNormal: scratch has the wrong dimension'
    end if
    if (size(output) /= self%n) then
      error stop 'MvNormal: output has the wrong dimension'
    end if

    do i = 1, self%n, 2
      call standard_normal_pair(z1, z2)
      scratch(i) = z1
      if (i + 1 <= self%n) scratch(i + 1) = z2
    end do

    output = self%mean
    do i = 1, self%n
      do j = 1, i
        output(i) = output(i) + self%lower(i,j) * scratch(j)
      end do
    end do
  end subroutine sample_into

  function sample(self) result(output)
    class(mvnormal_t), intent(in) :: self
    real(real64), allocatable :: output(:)
    real(real64), allocatable :: scratch(:)

    if (self%n == 0) then
      error stop 'MvNormal: uninitialized distribution'
    end if
    allocate(output(self%n), scratch(self%n))
    call self%sample_into(scratch, output)
  end function sample

  subroutine standard_normal_pair(first, second)
    real(real64), intent(out) :: first
    real(real64), intent(out) :: second
    real(real64) :: u1
    real(real64) :: u2
    real(real64) :: radius

    call random_number(u1)
    call random_number(u2)
    u1 = max(u1, tiny(1.0_real64))
    radius = sqrt(-2.0_real64 * log(u1))
    first = radius * cos(2.0_real64 * pi * u2)
    second = radius * sin(2.0_real64 * pi * u2)
  end subroutine standard_normal_pair

end module mvnormal_module
