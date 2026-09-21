program mvnormal_benchmark
  use iso_fortran_env, only : real64, int32, output_unit
  use mvnormal_module, only : mvnormal_t, mvnormal_create
  implicit none

  integer(int32) :: dimension
  integer(int32) :: samples
  integer(int32) :: repeats
  integer :: argument_count
  integer :: argument_index
  integer :: i
  integer :: j
  integer :: repeat_index
  integer :: sample_index
  integer :: seed_size
  character(len=256) :: argument
  character(len=256) :: value
  real(real64), allocatable :: mean(:)
  real(real64), allocatable :: covariance(:,:)
  real(real64), allocatable :: scratch(:)
  real(real64), allocatable :: output(:)
  integer, allocatable :: seed(:)
  type(mvnormal_t) :: distribution
  real(real64) :: setup_start
  real(real64) :: setup_finish
  real(real64) :: sample_start
  real(real64) :: sample_finish
  real(real64) :: setup_seconds
  real(real64) :: total_sample_seconds
  real(real64) :: repeat_sample_seconds
  real(real64) :: average_sample_seconds
  real(real64) :: minimum_sample_seconds
  real(real64) :: checksum

  dimension = 16_int32
  samples = 10000_int32
  repeats = 3_int32
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
  allocate(scratch(dimension), output(dimension))
  do i = 1, dimension
    mean(i) = real(i - 1, real64) / real(max(1, dimension), real64)
    do j = 1, dimension
      covariance(i,j) = 0.25_real64 ** abs(i - j)
    end do
  end do

  call random_seed(size=seed_size)
  allocate(seed(seed_size))
  do i = 1, seed_size
    seed(i) = 12345 + 37 * i
  end do
  call random_seed(put=seed)

  call cpu_time(setup_start)
  distribution = mvnormal_create(mean, covariance)
  call cpu_time(setup_finish)
  setup_seconds = setup_finish - setup_start

  total_sample_seconds = 0.0_real64
  minimum_sample_seconds = huge(1.0_real64)
  checksum = 0.0_real64
  do repeat_index = 1, repeats
    call cpu_time(sample_start)
    do sample_index = 1, samples
      call distribution%sample_into(scratch, output)
      checksum = checksum + sum(output)
    end do
    call cpu_time(sample_finish)
    repeat_sample_seconds = sample_finish - sample_start
    total_sample_seconds = total_sample_seconds + repeat_sample_seconds
    minimum_sample_seconds = min(minimum_sample_seconds, repeat_sample_seconds)
  end do
  average_sample_seconds = total_sample_seconds / &
                           real(repeats, real64)

  write(output_unit, '("fortran,", I0, ",", I0, ",", I0, ",", ES0.8, ",", ES0.8, ",", ES0.8, ",", ES0.17)') &
       dimension, samples, repeats, setup_seconds, average_sample_seconds, &
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

end program mvnormal_benchmark
