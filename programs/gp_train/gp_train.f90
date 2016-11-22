! train (and optionally optimize) a sparse GP from data files

program gp_in
  use m_util
  use m_gp
  use m_gp_optim
  
  implicit none

  type(SparseGP) :: gp
  
  real(dp) :: noise, a(1)
  ! number of training points, sparse points and dimension of the input
  integer :: n, nsparse, input_dimension
  ! unit number
  integer :: u
  ! loop counter
  integer :: i

  real(dp), dimension(:), allocatable :: theta, t, lbounds, ubounds
  real(dp), dimension(:,:), allocatable :: x
  integer, dimension(:), allocatable :: obs_type
  logical :: optimize

  integer :: optimize_max_iter
  real(dp) :: optimize_ftol

  namelist/DIMENSIONS/input_dimension,n,nsparse
  namelist/HYPERPARAMETERS/noise, theta, lbounds, ubounds
  namelist/CONTROL/optimize, optimize_max_iter, optimize_ftol
  
  input_dimension = -1
  n = -1
  nsparse = -1
 
  read(*, nml=DIMENSIONS)
  
  if (input_dimension.le.0.or.n.le.0.or.nsparse.le.0) then
     write (*,*) "Invalid DIMENSIONS input block"
     stop 1
  end if

  if (nsparse.gt.n) then
     write (*,*) "More sparse points requested than total inputs (DIMENSIONS block)"
     stop 1
  end if

  allocate(real(dp) :: theta(input_dimension))
  allocate(real(dp) :: lbounds(input_dimension+1))
  allocate(real(dp) :: ubounds(input_dimension+1))
  allocate(real(dp) :: t(n))
  allocate(real(dp) :: x(n,input_dimension))
  allocate(integer :: obs_type(n))
 
  noise = -1
  theta = -1
  lbounds(1) = 0.0
  lbounds(2:) = 0.01
  ubounds(:) = 100.0
  read(*, nml=HYPERPARAMETERS)

  if (noise.le.0.or.any(theta.le.0)) then
     write (*,*) "Invalid HYPERPARAMETERS input block"
     stop 1
  end if

  optimize = .false.
  optimize_max_iter = 100
  optimize_ftol = 1.0d-2
  read(*, nml=CONTROL)

  open(newunit=u, file="DATA")
  
  read (u,*) (x(i,:), obs_type(i), t(i), i=1,n)

  gp = SparseGP(nsparse, noise, theta, x, obs_type, t)

  if (optimize) then
     call log_lik_optim(input_dimension+1, gp, lbounds, ubounds, optimize_max_iter, optimize_ftol)
  end if

  call write_SparseGP("out.gp", gp) 
  
end program gp_in
