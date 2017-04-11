! train (and optionally optimize) a sparse GP from data files

program gp_in
  use m_util
  use m_gp
  use m_gp_sparse
  use m_gp_dense
  use m_gp_optim
  use m_cov_all
  use m_noise_all

  implicit none

!  type(SparseGP) :: gp
  class(BaseGP), allocatable :: gp
  
  ! number of training points, sparse points and dimension of the input
  integer :: n, nsparse, input_dimension
  ! unit number
  integer :: u
  ! loop counter
  integer :: i
  ! number of hyperparameters and noise parameters required
  integer :: ntheta, nnu

  real(dp), dimension(:), allocatable :: theta, nu, t, lbounds, ubounds
  real(dp), dimension(:,:), allocatable :: x
  integer, dimension(:), allocatable :: obs_type
  logical :: optimize

  integer :: optimize_max_iter
  real(dp) :: optimize_ftol

  character(len=max_name_len) :: covariance_function
  character(len=max_name_len) :: noise_model_name

  class(cov_fn), allocatable :: cf
  class(noise_model), allocatable :: nm

  namelist/DIMENSIONS/input_dimension, n, nsparse
  namelist/MODEL/covariance_function, noise_model_name
  namelist/HYPERPARAMETERS/nu, theta, lbounds, ubounds
  namelist/CONTROL/optimize, optimize_max_iter, optimize_ftol
  
  input_dimension = -1
  n = -1
  nsparse = -1
 
  read(*, nml=DIMENSIONS)
  
  if (input_dimension.le.0.or.n.le.0.or.nsparse.lt.-1) then
     write (*,*) "Invalid DIMENSIONS input block"
     stop 1
  end if
  if (nsparse.gt.n) then
     write (*,*) "More sparse points requested than total inputs (DIMENSIONS block)"
     stop 1
  end if

  read(*, nml=MODEL)

  call string_to_cov_fn(covariance_function, cf)
  call string_to_noise_model(noise_model_name, nm)

  nnu = nm%nparams_required(input_dimension)
  ntheta = cf%ntheta_required(input_dimension)

  allocate(real(dp) :: nu(nnu))
  allocate(real(dp) :: theta(ntheta))
  allocate(real(dp) :: lbounds(nnu + ntheta))
  allocate(real(dp) :: ubounds(nnu + ntheta))
  allocate(real(dp) :: t(n))
  allocate(real(dp) :: x(n,input_dimension))
  allocate(integer :: obs_type(n))
 
  nu = -1
  theta = -1
  lbounds(1) = 0.0
  lbounds(2:) = 0.01
  ubounds(:) = 100.0
  read(*, nml=HYPERPARAMETERS)

  optimize = .false.
  optimize_max_iter = 100
  optimize_ftol = 1.0d-2
  read(*, nml=CONTROL)

  open(newunit=u, file="DATA")
  
  read (u,*) (x(i,:), obs_type(i), t(i), i=1,n)

  if (nsparse.ge.n.or.nsparse.eq.-1) then
     allocate(gp, source=DenseGP(nu, theta, x, obs_type, t, cf, nm))
  else
     allocate(gp, source=SparseGP(nsparse, nu, theta, x, obs_type, t, cf, nm))
  end if

  if (optimize) then
     call log_lik_optim(nnu + ntheta, gp, lbounds, ubounds, optimize_max_iter, optimize_ftol)
  end if

  call gp%write_out("out.gp")

end program gp_in
