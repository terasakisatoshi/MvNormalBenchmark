module mvnormal_module
  use iso_fortran_env, only : real64, int64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  real(real64), parameter :: symmetry_tolerance = 1.0e-12_real64
  integer, parameter, public :: normal_algorithm_polar = 1
  integer, parameter, public :: normal_algorithm_ziggurat = 2
  integer, parameter :: ziggurat_layers = 256
  real(real64), parameter :: ziggurat_r = 3.6541528853610088_real64
  real(real64), parameter :: ziggurat_v = 0.004928673233974658_real64
  real(real64), parameter :: two_pow_63 = 9.223372036854775808e18_real64
  real(real64), save :: ziggurat_w(0:ziggurat_layers - 1)
  real(real64), save :: ziggurat_f(0:ziggurat_layers - 1)
  integer(int64), save :: ziggurat_k(0:ziggurat_layers - 1)
  logical, save :: ziggurat_initialized = .false.

  public :: mvnormal_t
  public :: mvnormal_create
  public :: normal_rng_t

  type, public :: normal_rng_t
    integer(int64), private :: state = 88172645463393265_int64
    integer, private :: algorithm = normal_algorithm_polar
    real(real64), private :: spare_normal = 0.0_real64
    logical, private :: has_spare = .false.
  contains
    procedure, public :: seed => normal_rng_seed
    procedure, public :: set_algorithm => normal_rng_set_algorithm
    procedure, public :: fill_standard_normals
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

    if (self%n == 0) then
      error stop 'MvNormal: uninitialized distribution'
    end if
    if (size(scratch) /= self%n) then
      error stop 'MvNormal: scratch has the wrong dimension'
    end if
    if (size(output) /= self%n) then
      error stop 'MvNormal: output has the wrong dimension'
    end if

    call normal_rng_fill_standard_normals(rng, scratch)

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

    call mvnormal_sample_inplace(self, rng, output)
  end subroutine sample_inplace

  subroutine mvnormal_sample_inplace(self, rng, output)
    type(mvnormal_t), intent(in) :: self
    type(normal_rng_t), intent(inout) :: rng
    real(real64), intent(inout) :: output(:)
    integer :: i
    integer :: j
    real(real64) :: z_j

    if (self%n == 0) then
      error stop 'MvNormal: uninitialized distribution'
    end if
    if (size(output) /= self%n) then
      error stop 'MvNormal: output has the wrong dimension'
    end if

    call normal_rng_fill_standard_normals(rng, output)

    ! Process columns from right to left.  output(j) still contains z(j)
    ! when column j is reached, and lower(i,j) is contiguous in this order.
    do j = self%n, 1, -1
      z_j = output(j)
      output(j) = self%mean(j)
      do i = j, self%n
        output(i) = output(i) + self%lower(i,j) * z_j
      end do
    end do
  end subroutine mvnormal_sample_inplace

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
    self%has_spare = .false.
  end subroutine normal_rng_seed

  subroutine normal_rng_set_algorithm(self, algorithm)
    class(normal_rng_t), intent(inout) :: self
    integer, intent(in) :: algorithm

    select case (algorithm)
    case (normal_algorithm_polar)
      continue
    case (normal_algorithm_ziggurat)
      call initialize_ziggurat_tables()
    case default
      error stop 'normal RNG: unknown algorithm'
    end select
    self%algorithm = algorithm
    self%has_spare = .false.
  end subroutine normal_rng_set_algorithm

  function normal_rng_next_u64(self) result(value)
    type(normal_rng_t), intent(inout) :: self
    integer(int64) :: value

    value = self%state
    value = ieor(value, ishft(value, 13))
    value = ieor(value, ishft(value, -7))
    value = ieor(value, ishft(value, 17))
    self%state = value
  end function normal_rng_next_u64

  function normal_rng_uniform_open01(self) result(value)
    type(normal_rng_t), intent(inout) :: self
    integer(int64) :: raw
    real(real64) :: value

    ! Match the other implementations: use the top 53 bits of the
    ! xorshift64 output to construct a double in (0, 1).
    raw = ishft(normal_rng_next_u64(self), -11)
    value = (real(raw, real64) + 0.5_real64) * &
            (1.0_real64 / 9.007199254740992e15_real64)
  end function normal_rng_uniform_open01

  subroutine fill_standard_normals(self, values)
    class(normal_rng_t), intent(inout) :: self
    real(real64), intent(out) :: values(:)

    call normal_rng_fill_standard_normals(self, values)
  end subroutine fill_standard_normals

  subroutine normal_rng_fill_standard_normals(self, values)
    type(normal_rng_t), intent(inout) :: self
    real(real64), intent(out) :: values(:)
    integer :: i
    real(real64) :: first
    real(real64) :: second

    if (self%algorithm == normal_algorithm_polar) then
      i = 1
      do while (i + 1 <= size(values))
        call normal_rng_polar_pair(self, first, second)
        values(i) = first
        values(i + 1) = second
        i = i + 2
      end do
      if (i <= size(values)) values(i) = normal_rng_standard_normal(self)
    else
      ! set_algorithm initializes the tables before this hot path is used.
      do i = 1, size(values)
        values(i) = normal_rng_ziggurat(self)
      end do
    end if
  end subroutine normal_rng_fill_standard_normals

  function normal_rng_standard_normal(self) result(value)
    type(normal_rng_t), intent(inout) :: self
    real(real64) :: value
    real(real64) :: second

    if (self%algorithm == normal_algorithm_ziggurat) then
      value = normal_rng_ziggurat(self)
      return
    end if
    if (self%has_spare) then
      value = self%spare_normal
      self%has_spare = .false.
      return
    end if
    call normal_rng_polar_pair(self, value, second)
    self%spare_normal = second
    self%has_spare = .true.
  end function normal_rng_standard_normal

  subroutine standard_normal_pair(self, first, second)
    class(normal_rng_t), intent(inout) :: self
    real(real64), intent(out) :: first
    real(real64), intent(out) :: second

    call normal_rng_standard_normal_pair(self, first, second)
  end subroutine standard_normal_pair

  subroutine normal_rng_standard_normal_pair(self, first, second)
    type(normal_rng_t), intent(inout) :: self
    real(real64), intent(out) :: first
    real(real64), intent(out) :: second

    first = normal_rng_standard_normal(self)
    second = normal_rng_standard_normal(self)
  end subroutine normal_rng_standard_normal_pair

  subroutine normal_rng_polar_pair(self, first, second)
    type(normal_rng_t), intent(inout) :: self
    real(real64), intent(out) :: first
    real(real64), intent(out) :: second
    real(real64) :: u
    real(real64) :: v
    real(real64) :: radius_squared
    real(real64) :: scale

    do
      u = 2.0_real64 * normal_rng_uniform_open01(self) - 1.0_real64
      v = 2.0_real64 * normal_rng_uniform_open01(self) - 1.0_real64
      radius_squared = u * u + v * v
      if (radius_squared > 0.0_real64 .and. radius_squared < 1.0_real64) exit
    end do
    scale = sqrt(-2.0_real64 * log(radius_squared) / radius_squared)
    first = u * scale
    second = v * scale
  end subroutine normal_rng_polar_pair

  subroutine initialize_ziggurat_tables()
    integer :: i
    real(real64) :: tail_density
    real(real64) :: q
    real(real64) :: x_next
    real(real64) :: x_current

    if (ziggurat_initialized) return
    tail_density = exp(-0.5_real64 * ziggurat_r * ziggurat_r)
    q = ziggurat_v / tail_density
    ziggurat_f = 0.0_real64
    ziggurat_w = 0.0_real64
    ziggurat_k = 0_int64
    ziggurat_f(0) = 1.0_real64
    ziggurat_f(ziggurat_layers - 1) = tail_density
    ziggurat_w(0) = q / two_pow_63
    ziggurat_w(ziggurat_layers - 1) = ziggurat_r / two_pow_63
    ziggurat_k(0) = int(ziggurat_r / q * two_pow_63, int64)
    ziggurat_k(1) = 0_int64
    x_next = ziggurat_r
    do i = ziggurat_layers - 2, 1, -1
      x_current = sqrt(-2.0_real64 * log(ziggurat_v / x_next + &
                                         exp(-0.5_real64 * x_next * x_next)))
      ziggurat_f(i) = exp(-0.5_real64 * x_current * x_current)
      ziggurat_w(i) = x_current / two_pow_63
      ziggurat_k(i + 1) = int(x_current / x_next * two_pow_63, int64)
      x_next = x_current
    end do
    ziggurat_initialized = .true.
  end subroutine initialize_ziggurat_tables

  function normal_rng_ziggurat(self) result(value)
    type(normal_rng_t), intent(inout) :: self
    real(real64) :: value
    integer(int64) :: bits
    integer(int64) :: magnitude
    integer :: index
    real(real64) :: sign
    real(real64) :: x
    real(real64) :: y
    real(real64) :: tail_x
    real(real64) :: tail_y

    do
      bits = normal_rng_next_u64(self)
      index = int(iand(bits, int(z'FF', int64)))
      if (btest(bits, 63)) then
        sign = -1.0_real64
      else
        sign = 1.0_real64
      end if
      magnitude = iand(bits, int(z'7FFFFFFFFFFFFFFF', int64))
      x = real(magnitude, real64) * ziggurat_w(index)
      if (magnitude < ziggurat_k(index)) then
        value = sign * x
        return
      end if

      if (index == 0) then
        do
          tail_x = -log(normal_rng_uniform_open01(self)) / ziggurat_r
          tail_y = -log(normal_rng_uniform_open01(self))
          if (2.0_real64 * tail_y >= tail_x * tail_x) then
            value = sign * (ziggurat_r + tail_x)
            return
          end if
        end do
      end if

      y = ziggurat_f(index) + normal_rng_uniform_open01(self) * &
          (ziggurat_f(index - 1) - ziggurat_f(index))
      if (y < exp(-0.5_real64 * x * x)) then
        value = sign * x
        return
      end if
    end do
  end function normal_rng_ziggurat

end module mvnormal_module
