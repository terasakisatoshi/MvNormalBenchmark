program mvnormal_benchmark
  use iso_fortran_env, only : real64, int32, int64, output_unit
  use mvnormal_module, only : mvnormal_t, mvnormal_create, normal_rng_t, &
                              normal_algorithm_polar, normal_algorithm_ziggurat
  implicit none

  integer(int64), parameter :: comparison_normal_seed = 1592598561_int64

  integer(int32) :: dimension
  integer(int32) :: samples
  integer(int32) :: repeats
  integer :: normal_algorithm
  integer :: argument_count
  integer :: argument_index
  integer :: i
  integer :: j
  integer :: repeat_index
  integer :: sample_index
  character(len=256) :: argument
  character(len=256) :: value
  character(len=16) :: algorithm_label
  real(real64), allocatable :: mean(:)
  real(real64), allocatable :: covariance(:,:)
  real(real64), allocatable :: output(:)
  type(mvnormal_t) :: distribution
  type(normal_rng_t) :: rng
  integer(int64) :: clock_start
  integer(int64) :: clock_finish
  integer(int64) :: clock_rate
  real(real64) :: setup_seconds
  real(real64) :: total_sample_seconds
  real(real64) :: repeat_sample_seconds
  real(real64) :: average_sample_seconds
  real(real64) :: minimum_sample_seconds
  real(real64) :: checksum

  dimension = 16_int32
  samples = 10000_int32
  repeats = 3_int32
  normal_algorithm = normal_algorithm_polar
  argument_count = command_argument_count()
  argument_index = 1

  do while (argument_index <= argument_count)
    call get_command_argument(argument_index, argument)
    select case (trim(argument))
    case ('--dim')
      call next_value(argument_index, value)
      call parse_integer(value, dimension)
    case ('--samples')
      call next_value(argument_index, value)
      call parse_integer(value, samples)
    case ('--repeats')
      call next_value(argument_index, value)
      call parse_integer(value, repeats)
    case ('--normal')
      call next_value(argument_index, value)
      call parse_normal(value, normal_algorithm)
    case default
      if (index(trim(argument), '--dim=') == 1) then
        value = argument(7:)
        call parse_integer(value, dimension)
      else if (index(trim(argument), '--samples=') == 1) then
        value = argument(11:)
        call parse_integer(value, samples)
      else if (index(trim(argument), '--repeats=') == 1) then
        value = argument(11:)
        call parse_integer(value, repeats)
      else if (index(trim(argument), '--normal=') == 1) then
        value = argument(10:)
        call parse_normal(value, normal_algorithm)
      else
        error stop 'usage: benchmark --dim N --samples N --repeats N'
      end if
    end select
    argument_index = argument_index + 1
  end do

  if (dimension <= 0 .or. samples <= 0 .or. repeats <= 0) then
    error stop '--dim, --samples, and --repeats must be positive'
  end if

  allocate(mean(dimension), covariance(dimension,dimension))
  allocate(output(dimension))
  do i = 1, dimension
    mean(i) = 0.01_real64 * real(i - 1, real64)
    do j = 1, dimension
      covariance(i,j) = scale(1.0_real64, -2 * abs(i - j))
    end do
  end do

  call rng%seed(comparison_normal_seed)
  call rng%set_algorithm(normal_algorithm)
  call system_clock(count_rate=clock_rate)

  call system_clock(count=clock_start)
  distribution = mvnormal_create(mean, covariance)
  call system_clock(count=clock_finish)
  setup_seconds = real(clock_finish - clock_start, real64) / &
                  real(clock_rate, real64)

  total_sample_seconds = 0.0_real64
  minimum_sample_seconds = huge(1.0_real64)
  checksum = 0.0_real64
  do repeat_index = 1, repeats
    call system_clock(count=clock_start)
    do sample_index = 1, samples
      call distribution%sample_inplace(rng, output)
      checksum = checksum + sum(output)
    end do
    call system_clock(count=clock_finish)
    repeat_sample_seconds = real(clock_finish - clock_start, real64) / &
                            real(clock_rate, real64)
    total_sample_seconds = total_sample_seconds + repeat_sample_seconds
    minimum_sample_seconds = min(minimum_sample_seconds, repeat_sample_seconds)
  end do
  average_sample_seconds = total_sample_seconds / &
                           real(repeats, real64)

  if (normal_algorithm == normal_algorithm_polar) then
    algorithm_label = 'fortran-polar'
  else
    algorithm_label = 'fortran-ziggurat'
  end if
  write(output_unit, '(A,",", I0, ",", I0, ",", I0, ",", ES0.8, ",", ES0.8, ",", ES0.8, ",", ES0.17)') &
       trim(algorithm_label), dimension, samples, repeats, setup_seconds, average_sample_seconds, &
       minimum_sample_seconds, checksum

contains

  subroutine next_value(current_index, result)
    integer, intent(inout) :: current_index
    character(len=*), intent(out) :: result

    if (current_index >= argument_count) then
      error stop 'missing value for command-line option'
    end if
    current_index = current_index + 1
    call get_command_argument(current_index, result)
  end subroutine next_value

  subroutine parse_integer(text, result)
    character(len=*), intent(in) :: text
    integer(int32), intent(out) :: result
    integer :: read_status

    read(text, *, iostat=read_status) result
    if (read_status /= 0) then
      error stop 'invalid integer command-line value'
    end if
  end subroutine parse_integer

  subroutine parse_normal(text, result)
    character(len=*), intent(in) :: text
    integer, intent(out) :: result

    select case (trim(text))
    case ('polar')
      result = normal_algorithm_polar
    case ('ziggurat')
      result = normal_algorithm_ziggurat
    case default
      error stop '--normal must be polar or ziggurat'
    end select
  end subroutine parse_normal

end program mvnormal_benchmark
