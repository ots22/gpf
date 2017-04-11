module m_gp
  use m_util
  use m_cov
  use m_noise

  private
  public BaseGP, nlog_lik, set_hyperparams

  type, abstract :: BaseGP
     ! noise hyperparameter (sigma^2)
     ! meaning depends on noise model 'noise_model'
     real(dp), dimension(:), allocatable :: nu 
     ! covariance hyperparameters
     ! meaning depends on covariance function `covariance'
     real(dp), dimension(:), allocatable :: theta
     ! inputs
     real(dp), dimension(:,:), allocatable :: x
     ! types of the observations
     integer, dimension(:), allocatable :: obs_type
     ! observations
     real(dp), dimension(:), allocatable :: t
     ! the covariance function
     class(cov_fn), allocatable :: covariance
     ! the noise model
     class(noise_model), allocatable :: noise_model
   contains
     procedure(log_lik), deferred :: log_lik
     procedure(update_matrices), deferred :: update_matrices
     procedure(predict), deferred :: predict
     procedure(write_out), deferred :: write_out
  end type BaseGP

  abstract interface
     function log_lik(this)
       use m_util, only: dp
       import BaseGP
       class(BaseGP), intent(in) :: this
       real(dp) log_lik
     end function log_lik

     subroutine update_matrices(this)
       import BaseGP
       class(BaseGP), intent(inout) :: this
     end subroutine update_matrices

     function predict(this, xnew, obs_type_new)
       use m_util, only: dp
       import BaseGP
       class(BaseGP), intent(in) :: this
       real(dp) predict
       real(dp), dimension(:), intent(in) :: xnew
       integer, optional, intent(in) :: obs_type_new
     end function predict

     subroutine write_out(this, filename)
       import BaseGP
       class(BaseGP), intent(in) :: this
       character(len=*), intent(in) :: filename
     end subroutine write_out
  end interface

contains

  subroutine set_hyperparams(gp, hypers)
    class(BaseGP), intent(inout) :: gp
    real(dp), dimension(:) :: hypers
    integer nnu
    integer ntheta
    nnu = size(gp%nu)
    ntheta = size(gp%theta)

    if ((any(hypers(1:nnu).ne.gp%nu)) &
         & .or.any(hypers(nnu+1:nnu+ntheta).ne.gp%theta)) then
       gp%nu = hypers(1:nnu)
       gp%theta = hypers(nnu+1:nnu+ntheta)
       call gp%update_matrices
    end if
  end subroutine set_hyperparams

  subroutine nlog_lik(val, n, hypers, grad, need_gradient, gp)
    class(BaseGP) :: gp
    integer :: n, need_gradient
    real(dp) :: val, hypers(n)
    real(dp), intent(inout) :: grad(n)

    call set_hyperparams(gp, hypers)
    val = gp%log_lik()

    call output_params(gp,val)
  end subroutine nlog_lik

  subroutine output_params(gp,val)
    class(BaseGP) :: gp
    real(dp) :: val
    integer, save :: u
    logical, save :: assigned = .false.
    if (.not.assigned) open(newunit=u, file="LOG_LIK_OPTIM"); assigned=.true.
    write (u,*) "log likelihood: noise = ", gp%nu, " theta = ", & 
         gp%theta, "log(likelihood) = ", val
  end subroutine output_params

end module m_gp
