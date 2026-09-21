module stdlib_mvn
  use mvn_common, only: dp
  use stdlib_linalg, only: chol
  use stdlib_stats_distribution_normal, only: rvs_normal
  implicit none
  private
  public :: mvn_plan_t, sample_stdlib

  type :: mvn_plan_t
    real(dp), allocatable :: mu(:)
    real(dp), allocatable :: lower(:, :)
  contains
    procedure :: init => plan_init
    procedure :: sample => plan_sample
  end type mvn_plan_t

contains

  subroutine plan_init(self, mu, sigma)
    class(mvn_plan_t), intent(inout) :: self
    real(dp), intent(in) :: mu(:)
    real(dp), intent(in) :: sigma(:, :)

    self%mu = mu
    self%lower = chol(sigma)
  end subroutine plan_init

  subroutine plan_sample(self, out)
    class(mvn_plan_t), intent(in) :: self
    real(dp), intent(out) :: out(:)

    out = self%mu + matmul(self%lower, &
                           rvs_normal(0.0_dp, 1.0_dp, size(self%mu)))
  end subroutine plan_sample

  subroutine sample_stdlib(mu, sigma, out)
    real(dp), intent(in) :: mu(:)
    real(dp), intent(in) :: sigma(:, :)
    real(dp), intent(out) :: out(:)
    type(mvn_plan_t) :: plan

    call plan%init(mu, sigma)
    call plan%sample(out)
  end subroutine sample_stdlib

end module stdlib_mvn
