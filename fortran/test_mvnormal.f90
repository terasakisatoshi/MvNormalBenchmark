program test_mvnormal
  use iso_fortran_env, only : real64, int64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mvnormal_module, only : mvnormal_t, mvnormal_create, normal_rng_t, &
                             normal_algorithm_polar, normal_algorithm_ziggurat
  implicit none
  integer :: algorithm, k
  integer, parameter :: dimensions(*) = [1, 2, 3, 4, 5, 6, 7, 8, 9, 17, 63, 64, 65, 128]

  do algorithm = normal_algorithm_polar, normal_algorithm_ziggurat
    call test_stream(algorithm)
    do k = 1, size(dimensions)
      call test_transform(dimensions(k), algorithm)
    end do
  end do
  print '(A)', 'Fortran tests passed (streams, odd/even dimensions, strided arrays, dense covariance).'

contains

  subroutine initialize_rng(rng, algorithm)
    type(normal_rng_t), intent(out) :: rng
    integer, intent(in) :: algorithm
    call rng%seed(1592598561_int64)
    call rng%set_algorithm(algorithm)
  end subroutine initialize_rng

  subroutine test_stream(algorithm)
    integer, intent(in) :: algorithm
    integer, parameter :: count = 20000
    type(normal_rng_t) :: bulk_rng, pair_rng, split_rng
    real(real64) :: bulk(count), pairs(count), split(count), empty(0)
    integer :: i, last, width

    call initialize_rng(bulk_rng, algorithm)
    call initialize_rng(pair_rng, algorithm)
    call initialize_rng(split_rng, algorithm)
    call bulk_rng%fill_standard_normals(bulk)
    do i = 1, count, 2
      call pair_rng%standard_normal_pair(pairs(i), pairs(i + 1))
    end do
    i = 1
    do while (i <= count)
      ! Interleave odd/even fills, empty fills and the scalar pair API.
      width = mod(i, 7) + 1
      last = min(count, i + width - 1)
      call split_rng%fill_standard_normals(split(i:last))
      call split_rng%fill_standard_normals(empty)
      i = last + 1
      if (i + 1 <= count) then
        call split_rng%standard_normal_pair(split(i), split(i + 1))
        i = i + 2
      end if
    end do
    if (.not. all(ieee_is_finite(bulk))) error stop 'nonfinite normal sample'
    if (any(bulk /= pairs)) error stop 'bulk and pair streams differ'
    if (any(bulk /= split)) error stop 'normal stream depends on fill boundaries'
  end subroutine test_stream

  subroutine test_transform(n, algorithm)
    integer, intent(in) :: n, algorithm
    type(mvnormal_t) :: distribution
    type(normal_rng_t) :: reference_rng, contiguous_rng, strided_rng, scratch_rng, allocating_rng
    real(real64) :: mean(n), lower(n,n), covariance(n,n), normals(n), expected(n)
    real(real64) :: output(n), strided(2*n), scratch(n), separate(n)
    real(real64), allocatable :: allocated(:)
    integer :: i, j, repeat

    ! A dense, signed triangular factor gives a known independent reference.
    lower = 0.0_real64
    do i = 1, n
      mean(i) = 0.01_real64 * real(i - 1, real64)
      lower(i,i) = 1.0_real64 + real(i, real64) / real(n, real64)
      do j = 1, i - 1
        lower(i,j) = real(mod(3*i + 7*j, 11) - 5, real64) / real(10*n, real64)
      end do
    end do
    covariance = matmul(lower, transpose(lower))
    distribution = mvnormal_create(mean, covariance)
    call initialize_rng(reference_rng, algorithm)
    call initialize_rng(contiguous_rng, algorithm)
    call initialize_rng(strided_rng, algorithm)
    call initialize_rng(scratch_rng, algorithm)
    call initialize_rng(allocating_rng, algorithm)
    strided = -999.0_real64
    do repeat = 1, 5
      call reference_rng%fill_standard_normals(normals)
      expected = mean + matmul(lower, normals)
      call distribution%sample_inplace(contiguous_rng, output)
      call distribution%sample_inplace(strided_rng, strided(1:2*n:2))
      call distribution%sample_into(scratch_rng, scratch, separate)
      allocated = distribution%sample(allocating_rng)
      call check_close(output, expected)
      call check_close(strided(1:2*n:2), expected)
      call check_close(separate, expected)
      call check_close(allocated, expected)
      if (any(strided(2:2*n:2) /= -999.0_real64)) error stop 'strided sample overwrote gaps'
    end do
  end subroutine test_transform

  subroutine check_close(actual, expected)
    real(real64), intent(in) :: actual(:), expected(:)
    if (.not. all(ieee_is_finite(actual))) error stop 'nonfinite multivariate sample'
    if (maxval(abs(actual - expected)) > 1.0e-12_real64 * max(1.0_real64, maxval(abs(expected)))) &
      error stop 'multivariate sample differs from dense reference'
  end subroutine check_close
end program test_mvnormal
