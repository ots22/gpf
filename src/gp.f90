module m_GP
  use m_util
  use m_cov_sqexp
  implicit none
  
  ! Projected process approximation to full GP
  type :: SparseGP
     ! number of sparse points
     integer :: m
     ! noise hyperparameter (sigma^2)
     real(dp) :: nu 
     ! hyperparameters
     real(dp), dimension(:), allocatable :: theta
     ! inputs
     real(dp), dimension(:,:), allocatable :: x
     ! types of the observations
     integer, dimension(:), allocatable :: obs_type
     ! observations
     real(dp), dimension(:), allocatable :: t
     
     ! reduced observations
     real(dp), dimension(:), allocatable :: Kmn_t, invQKmn_t
     ! sparsified covariance matrix and its inverse
     real(dp), dimension(:,:), allocatable :: invKmm, Q, invQ
  end type SparseGP

  interface SparseGP
     module procedure make_SparseGP
     module procedure read_SparseGP
  end interface SparseGP

contains

  subroutine alloc_SparseGP(gp, n, m, ntheta, d)
    type(SparseGP), intent(inout) :: gp
    integer, intent(in) :: n, m, ntheta, d

    if (m>n) then
       print *, "Error: Active data set larger than total data set."
       stop
    end if

    allocate(real(dp) :: gp%theta(ntheta))
    allocate(real(dp) :: gp%x(n, d))
    allocate(integer :: gp%obs_type(n))
    allocate(real(dp) :: gp%t(n))
    allocate(real(dp) :: gp%Kmn_t(m))
    allocate(real(dp) :: gp%invQKmn_t(m))
    allocate(real(dp) :: gp%invKmm(m, m))
    allocate(real(dp) :: gp%Q(m, m))
    allocate(real(dp) :: gp%invQ(m, m))
    
  end subroutine alloc_SparseGP

  subroutine update_matrices(gp)
    type(SparseGP), intent(inout) :: gp
    real(dp), dimension(:, :), allocatable :: Knm
    real(dp), dimension(:, :), allocatable :: noise
    integer :: i,j,m,n

    n = size(gp%t)
    m = gp%m
    
    allocate(real(dp) :: noise(m,m))
    allocate(real(dp) :: Knm(n,m))   

    noise = 0.0_dp
    do i=1,m
       noise(i,i) = merge(gp%nu, 0.0_dp, gp%obs_type(i).eq.0)
    end do

    do i=1,n
       do j=1,m
          Knm(i,j) = cov(gp%obs_type(i), gp%obs_type(j), &
               gp%x(i,:), gp%x(j,:), gp%theta)
          ! add a small diagonal term to stabilise the inversion
          if (i.eq.j) then 
             Knm(i,i) = Knm(i,i) + 1e-9
          end if
       end do
    end do

    gp%Kmn_t = matmul(gp%t, Knm)

    gp%Q = matmul(Knm(1:m,:), noise) + matmul(transpose(Knm), Knm)
    do i=1,m
       gp%Q(i,i) = gp%Q(i,i) + 1e-9
    end do
    
    gp%invQ = gp%Q
    call ninv(gp%invQ)

    gp%invQKmn_t = matmul(gp%invQ, gp%Kmn_t)

    gp%invKmm = Knm(1:m,:)
    call ninv(gp%invKmm)    
  end subroutine update_matrices
  
  function make_SparseGP(m, nu, theta, x, obs_type, t)
    type(SparseGP) :: make_SparseGP
    integer, intent(in) :: m
    real(dp), intent(in) :: nu
    real(dp), dimension(:), intent(in) :: theta
    real(dp), dimension(:,:), intent(in) :: x
    integer,  dimension(:), intent(in) :: obs_type
    real(dp), dimension(:), intent(in) :: t

    integer :: n
    n = size(t)

    call alloc_SparseGP(make_SparseGP, n, m, size(theta), size(x,2))

    make_SparseGP%m = m
    make_SparseGP%nu = nu
    make_SparseGP%theta = theta
    make_SparseGP%x = x
    make_SparseGP%obs_type = obs_type
    make_SparseGP%t = t
    
    call update_matrices(make_SparseGP)

  end function make_SparseGP

  subroutine write_SparseGP(filename, gp)
    type(SparseGP), intent(in) :: gp
    character(len=*), intent(in) :: filename
    integer :: u
    open(newunit=u, file=filename)
    write (u,'(I10)') size(gp%t), gp%m, size(gp%theta), size(gp%x,2)

    write (u,'(es24.15)') gp%nu, gp%theta, & 
         gp%x
    write (u,'(I4)') gp%obs_type
    write (u,'(es24.15)') gp%Kmn_t, gp%invQKmn_t, &
         gp%Q, gp%invQ
    close(u)

  end subroutine write_SparseGP
  
  function read_SparseGP(filename) result(gp)
    character(len=*), intent(in) :: filename
    type(SparseGP) :: gp
    integer n, ntheta, d, u
    open(newunit=u, file=filename)
    read (u,'(I10)') n, gp%m, ntheta, d

    ! check inputs read are sane

    call alloc_SparseGP(gp, n, gp%m, ntheta, d)
    
    read (u,'(es24.15)') gp%nu, gp%theta, & 
         gp%x
    read (u,'(I4)') gp%obs_type
    read (u,'(es24.15)') gp%Kmn_t, gp%invQKmn_t, &
         gp%Q, gp%invQ
    close(u)
  end function read_SparseGP

  pure function SparseGP_predict(gp, xnew, obs_type_new_arg)
    real(dp) :: SparseGP_predict
    type(SparseGP), intent(in) :: gp
    real(dp), dimension(:), intent(in) :: xnew
    integer, optional, intent(in) :: obs_type_new_arg
    
    integer :: obs_type_new, i

    real(dp), dimension(gp%m) :: km

    if (present(obs_type_new_arg)) then
       obs_type_new = obs_type_new_arg
    else
       obs_type_new = 0
    endif

    do i=1,gp%m
       km(i) = cov(obs_type_new, gp%obs_type(i), xnew, gp%x(i,:), gp%theta)
    end do
    
    SparseGP_predict = dot_product(km, gp%invQKmn_t)

  end function SparseGP_predict

  subroutine set_hyperparams(gp, hypers)
    type(SparseGP), intent(inout) :: gp
    real(dp), dimension(:) :: hypers

    if ((hypers(1).ne.gp%nu).or.all(hypers(2:).ne.gp%theta)) then
       gp%nu = hypers(1)
       gp%theta = hypers(2:)
       call update_matrices(gp)
    end if
  end subroutine set_hyperparams

  function log_lik(gp)
    type(SparseGP), intent(in) :: gp
    real(dp) :: log_lik
    integer m,n
    m = gp%m
    n = size(gp%t)
    log_lik = -0.5 * ((n-m)*log(gp%nu) + logdet(gp%Q) + (dot_product(gp%t, gp%t) &
         - dot_product(gp%Kmn_t, gp%invQKmn_t))/gp%nu)
  end function log_lik

  subroutine nlog_lik(val, n, hypers, grad, need_gradient, gp)
    type(SparseGP) :: gp
    integer :: n, need_gradient
    real(dp) :: val, hypers(n)
    real(dp), intent(inout) :: grad(n)

    call set_hyperparams(gp, hypers)
    val = log_lik(gp)

    call output_params(gp,val)
  end subroutine nlog_lik

  subroutine output_params(gp,val)
    type(SparseGP) :: gp
    real(dp) :: val
    integer, save :: u
    logical, save :: assigned = .false.
    if (.not.assigned) open(newunit=u, file="LOG_LIK_OPTIM"); assigned=.true.
    write (u,*) "log likelihood: noise = ", gp%nu, " theta = ", & 
         gp%theta, "log(likelihood) = ", val
  end subroutine output_params

end module m_GP
