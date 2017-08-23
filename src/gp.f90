! Base class for Gaussian process types

module m_gp
  use m_util
  use m_cov
  use m_noise

  private
  public BaseGP, nlog_lik, set_hyperparams

  type, abstract :: BaseGP
     ! The noise hyperparameter(s). This is passed to the noise model
     ! (`noise_model'), which determines its precise meaning.
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
     ! The log likelihood of the hyperparameters
     function log_lik(this)
       use m_util, only: dp
       import BaseGP
       class(BaseGP), intent(in) :: this
       real(dp) log_lik
     end function log_lik

     ! Helper routine to update the internal state.  Called when an
     ! observation or hyperparameter changes and the covariance matrix
     ! must be recomputed.
     subroutine update_matrices(this)
       import BaseGP
       class(BaseGP), intent(inout) :: this
     end subroutine update_matrices

     ! Make a prediction of the underlying function value at
     ! coordinate `xnew'.
     function predict(this, xnew, obs_type_new)
       use m_util, only: dp
       import BaseGP
       class(BaseGP), intent(in) :: this
       real(dp) predict
       real(dp), dimension(:), intent(in) :: xnew
       integer, optional, intent(in) :: obs_type_new
     end function predict

     ! Serialize to a file
     subroutine write_out(this, filename)
       import BaseGP
       class(BaseGP), intent(in) :: this
       character(len=*), intent(in) :: filename
     end subroutine write_out
  end interface

contains

  ! Set the hyperparameters of `gp' to `hypers', assuming they are
  ! ordered [nu(1:nnu), theta(1:ntheta)], where nu are the noise
  ! hyperparameters and theta are the covariance hyperparameters.
  ! Calls `update_matrices'.
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

  ! A routine with the required interface for the NLopt library, for
  ! maximizing the log-likelihood.
  subroutine nlog_lik(val, n, hypers, grad, need_gradient, gp)
    class(BaseGP) :: gp
    integer :: n, need_gradient
    real(dp) :: val, hypers(n)
    real(dp), intent(inout) :: grad(n)

    call set_hyperparams(gp, hypers)
    val = gp%log_lik()

    call output_params(gp,val)
  end subroutine nlog_lik

  ! Helper routine for nlog_lik: print out the current hyperparameters
  ! (associated with `gp') and their corresponding log-likelihood
  ! (`val').
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
