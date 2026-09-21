module mvnormal_module
  use iso_fortran_env, only : real64, int64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  real(real64), parameter :: symmetry_tolerance = 1.0e-12_real64

  public :: mvnormal_t
  public :: mvnormal_create
  public :: normal_rng_t

  type, public :: normal_rng_t
    integer(int64), private :: state = 88172645463393265_int64
  contains
    procedure, public :: seed => normal_rng_seed
    procedure, private :: next_u64 => normal_rng_next_u64
    procedure, private :: uniform_open01 => normal_rng_uniform_open01
    procedure, public :: standard_normal_pair
  end type normal_rng_t

  type, public :: mvnormal_t
    integer, private :: n = 0
    real(real64), allocatable, private :: mean(:)
    real(real64), allocatable, private :: lower(:,:)
  contains
    procedure, public :: sample_into
    procedure, public :: sample_inplace
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

  subroutine sample_into(self, rng, scratch, output)
    class(mvnormal_t), intent(in) :: self
    type(normal_rng_t), intent(inout) :: rng
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
      call rng%standard_normal_pair(z1, z2)
      scratch(i) = z1
      if (i + 1 <= self%n) scratch(i + 1) = z2
    end do

    output = self%mean
    ! Fortran stores arrays column-major.  Accumulate one column of L at a
    ! time so lower(i,j) and output(i) are both traversed contiguously.
    do j = 1, self%n
      do i = j, self%n
        output(i) = output(i) + self%lower(i,j) * scratch(j)
      end do
    end do
  end subroutine sample_into

  subroutine sample_inplace(self, rng, output)
    class(mvnormal_t), intent(in) :: self
    type(normal_rng_t), intent(inout) :: rng
    real(real64), intent(inout) :: output(:)
    integer :: i
    integer :: j
    real(real64) :: z1
    real(real64) :: z2
    real(real64) :: z_j

    if (self%n == 0) then
      error stop 'MvNormal: uninitialized distribution'
    end if
    if (size(output) /= self%n) then
      error stop 'MvNormal: output has the wrong dimension'
    end if

    do i = 1, self%n, 2
      call rng%standard_normal_pair(z1, z2)
      output(i) = z1
      if (i + 1 <= self%n) output(i + 1) = z2
    end do

    ! Process columns from right to left.  output(j) still contains z(j)
    ! when column j is reached, and lower(i,j) is contiguous in this order.
    do j = self%n, 1, -1
      z_j = output(j)
      output(j) = self%mean(j)
      do i = j, self%n
        output(i) = output(i) + self%lower(i,j) * z_j
      end do
    end do
  end subroutine sample_inplace

  function sample(self, rng) result(output)
    class(mvnormal_t), intent(in) :: self
    type(normal_rng_t), intent(inout) :: rng
    real(real64), allocatable :: output(:)

    if (self%n == 0) then
      error stop 'MvNormal: uninitialized distribution'
    end if
    allocate(output(self%n))
    call self%sample_inplace(rng, output)
  end function sample

  subroutine normal_rng_seed(self, seed)
    class(normal_rng_t), intent(inout) :: self
    integer(int64), intent(in) :: seed

    if (seed == 0_int64) then
      self%state = 88172645463393265_int64
    else
      self%state = seed
    end if
  end subroutine normal_rng_seed

  function normal_rng_next_u64(self) result(value)
    class(normal_rng_t), intent(inout) :: self
    integer(int64) :: value

    value = self%state
    value = ieor(value, ishft(value, 13))
    value = ieor(value, ishft(value, -7))
    value = ieor(value, ishft(value, 17))
    self%state = value
  end function normal_rng_next_u64

  function normal_rng_uniform_open01(self) result(value)
    class(normal_rng_t), intent(inout) :: self
    integer(int64) :: raw
    real(real64) :: value

    raw = iand(self%next_u64(), int(z'7FFFFFFFFFFFFFFF', int64))
    value = (real(raw, real64) + 0.5_real64) * &
            (1.0_real64 / 9.223372036854775808e18_real64)
  end function normal_rng_uniform_open01

  subroutine standard_normal_pair(self, first, second)
    class(normal_rng_t), intent(inout) :: self
    real(real64), intent(out) :: first
    real(real64), intent(out) :: second
    real(real64) :: u
    real(real64) :: v
    real(real64) :: radius_squared
    real(real64) :: scale

    do
      u = 2.0_real64 * self%uniform_open01() - 1.0_real64
      v = 2.0_real64 * self%uniform_open01() - 1.0_real64
      radius_squared = u * u + v * v
      if (radius_squared > 0.0_real64 .and. radius_squared < 1.0_real64) exit
    end do
    scale = sqrt(-2.0_real64 * log(radius_squared) / radius_squared)
    first = u * scale
    second = v * scale
  end subroutine standard_normal_pair

end module mvnormal_module
