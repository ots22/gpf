! # Projected-process approximation to full GP. 
!
! We loosely follow the notation of [1].  Given n observations, we
! construct an approximation to the covariance matrix of rank m.  To
! do this, we must choose m distinguished `support points' from amoung
! the training data.  Here, they are simply taken to be the first m
! inputs in the array of training points passed to the constructor.
!
! The full covariance matrix C is represented as
!
!   C = Knm * inv(Kmm) * Kmn
!
! where 'n' and 'm' are subscripts and represent the size along the
! indicated dimension (e.g. Kmn is an m by n matrix).  This is never
! computed explicitly.
!
! * Knm is the n by m covariance matrix of the inputs in the training
! set with each support point.
! 
! * Kmn is the transpose of Knm.
!
! * Kmm is the m by m covariance matrix of the support points.
!
! We also define
!
!   Q = nu**2 * Kmm + Kmn * Knm
! 
! where nu**2 is the variance of the input-independent noise, and I is
! the m by m identity matrix.  The matrix Q is useful when making
! predictions and computing the liklihood.  This is not to be confused
! with the \tilde{Q} term defined in reference [2].
!
! ## References
!
! [1] C. Rasmussen and C. Williams. Gaussian Processes for Machine
! Learning. Adaptative Computation and Machine Learning Series. MIT
! Press, 2006.
!
! [2] J. Quinonero and C. Rasmussen. Analysis of Some Methods for
! Reduced Rank Gaussian Process Regression, in Switching and Learning
! in Feedback Systems: European Summer School on Multi-Agent Control,
! Maynooth, Ireland, September 8-10, 2003, Revised Lectures and
! Selected Papers, Springer, 2005

module m_GP_sparse
  use m_gp
  use m_util
  use m_cov_all
  use m_noise_all
  implicit none

  private
  public SparseGP, read_SparseGP

  type, extends(BaseGP) :: SparseGP
     ! Number of sparse points (and the rank of the covariance matrix
     ! under this approximation)
     integer :: m
     
     ! See the above description for definitions of Kmn, Kmm, Q etc.
     
     ! The reduced observations. Kmn_t = Kmn * t, where t is the
     ! vector of training observations
     real(dp), dimension(:), allocatable :: Kmn_t
     ! inv(Q) * Kmn * t
     real(dp), dimension(:), allocatable :: invQKmn_t
     ! inv(Kmm)
     real(dp), dimension(:,:), allocatable :: invKmm
     ! Q as described above: the sparsified covariance matrix and its inverse
     real(dp), dimension(:,:), allocatable :: Q, invQ
   contains
     procedure log_lik
     procedure update_matrices
     procedure predict
     procedure write_out
  end type SparseGP

  interface SparseGP
     module procedure make_SparseGP
     module procedure read_SparseGP
  end interface SparseGP

contains

  subroutine alloc_SparseGP(gp, n, m, ntheta, d, CovFunction, NoiseModel)
    type(SparseGP), intent(inout) :: gp
    integer, intent(in) :: n, m, ntheta, d
    class(cov_fn) :: CovFunction
    class(noise_model) :: NoiseModel

    if (m>n) then
       print *, "Error: Active data set larger than total data set."
       stop
    end if

    allocate(real(dp) :: gp%nu(NoiseModel%nparams_required(d)))
    allocate(real(dp) :: gp%theta(ntheta))
    allocate(real(dp) :: gp%x(n, d))
    allocate(integer :: gp%obs_type(n))
    allocate(real(dp) :: gp%t(n))
    allocate(real(dp) :: gp%Kmn_t(m))
    allocate(real(dp) :: gp%invQKmn_t(m))
    allocate(real(dp) :: gp%invKmm(m, m))
    allocate(real(dp) :: gp%Q(m, m))
    allocate(real(dp) :: gp%invQ(m, m))
    allocate(gp%covariance, mold=CovFunction)
    allocate(gp%noise_model, mold=NoiseModel)
  end subroutine alloc_SparseGP

  subroutine update_matrices(this)
    class(SparseGP), intent(inout) :: this
    real(dp), dimension(:, :), allocatable :: Knm
    real(dp), dimension(:, :), allocatable :: noise
    integer :: i,j,m,n

    n = size(this%t)
    m = this%m
    
    allocate(real(dp) :: noise(m,m))
    allocate(real(dp) :: Knm(n,m))   

    noise = 0.0_dp
    do i=1,m
       noise(i,i) = this%noise_model%noise(this%obs_type(i), this%nu)
    end do

    do i=1,n
       do j=1,m
          Knm(i,j) = this%covariance%cov(this%obs_type(i), this%obs_type(j), &
               this%x(i,:), this%x(j,:), this%theta)
          ! add a small diagonal term to stabilise the inversion
          if (i.eq.j) then 
             Knm(i,i) = Knm(i,i) + 1e-9
          end if
       end do
    end do

    this%Kmn_t = matmul(this%t, Knm)

    this%Q = matmul(Knm(1:m,:), noise) + matmul(transpose(Knm), Knm)
    do i=1,m
       this%Q(i,i) = this%Q(i,i) + 1e-9
    end do

    this%invQ = this%Q
    call ninv(this%invQ)

    this%invQKmn_t = solve(this%Q, this%Kmn_t)
    this%invKmm = Knm(1:m,:)
    call ninv(this%invKmm)
  end subroutine update_matrices
  
  function make_SparseGP(m, nu, theta, x, obs_type, t, CovFunction, NoiseModel) result(gp)
    type(SparseGP) :: gp
    integer, intent(in) :: m
    real(dp), dimension(:), intent(in) :: nu
    real(dp), dimension(:), intent(in) :: theta
    real(dp), dimension(:,:), intent(in) :: x
    integer,  dimension(:), intent(in) :: obs_type
    real(dp), dimension(:), intent(in) :: t

    class(cov_fn) :: CovFunction
    class(noise_model) :: NoiseModel
    
    integer :: n, d
    n = size(t)
    d = size(x,2)

    if (size(theta) /= CovFunction%ntheta_required(d)) then
       print *, "size of theta does not match number of hyperparameters required by the covariance function"
       stop 1
    end if
    
    if (size(nu) /= NoiseModel%nparams_required(d)) then
       print *, "size of nu (noise params) does not match number required by the noise model"
       stop 1
    end if

    call alloc_SparseGP(gp, n, m, size(theta), d, CovFunction, NoiseModel)

    gp%m = m
    gp%nu = nu
    gp%theta = theta
    gp%x = x
    gp%obs_type = obs_type
    gp%t = t
    
    call update_matrices(gp)

  end function make_SparseGP

  subroutine write_out(this, filename)
    class(SparseGP), intent(in) :: this
    character(len=*), intent(in) :: filename
    character(len=max_name_len) :: cov_fn_name
    character(len=max_name_len) :: noise_model_name
    integer :: u ! unit number for output

    cov_fn_name = cov_fn_to_string(this%covariance)
    noise_model_name = noise_model_to_string(this%noise_model)

    open(newunit=u, file=filename)

    write (u,'(A)') "SparseGP"
    write (u,'(I10)') size(this%t), this%m, size(this%theta), size(this%nu), size(this%x,2)
    write (u,'(A)') trim(cov_fn_name)
    write (u,'(A)') trim(noise_model_name)
    write (u,'(es24.15)') this%nu, this%theta, this%x
    write (u,'(I4)') this%obs_type
    write (u,'(es24.15)') this%Kmn_t, this%invQKmn_t, this%Q, this%invQ

    close(u)
  end subroutine write_out
  
  function read_SparseGP(filename) result(gp)
    character(len=*), intent(in) :: filename
    type(SparseGP) :: gp
    integer n, ntheta, nnu, d, u
    character(len=max_name_len) :: label
    character(len=max_name_len) :: cov_fn_name
    character(len=max_name_len) :: noise_model_name
    class(cov_fn), allocatable :: CovFunction
    class(noise_model), allocatable :: NoiseModel
    open(newunit=u, file=filename)
    read (u,'(A)') label

    if (trim(label) /= "SparseGP") then
       print *, "read_SparseGP: Incompatible data file"
       stop 1
    end if

    read (u,'(I10)') n, gp%m, ntheta, nnu, d

    read (u,'(A)') cov_fn_name 
    call string_to_cov_fn(cov_fn_name, CovFunction)
    if (ntheta /= CovFunction%ntheta_required(d)) then
       print *, "read_SparseGP: ntheta does not match number required by the covariance function"
       stop 1
    end if

    read (u,'(A)') noise_model_name
    call string_to_noise_model(noise_model_name, NoiseModel)
    if (nnu /= NoiseModel%nparams_required(d)) then
       print *, "read_SparseGP: size of nu (noise params) does not match number required by the noise model"
       stop 1
    end if

    call alloc_SparseGP(gp, n, gp%m, ntheta, d, CovFunction, NoiseModel)
    
    read (u,'(es24.15)') gp%nu, gp%theta, & 
         gp%x
    read (u,'(I4)') gp%obs_type
    read (u,'(es24.15)') gp%Kmn_t, gp%invQKmn_t, &
         gp%Q, gp%invQ
    close(u)
  end function read_SparseGP

  pure function predict(this, xnew, obs_type_new)
    real(dp) :: predict
    class(SparseGP), intent(in) :: this 
    real(dp), dimension(:), intent(in) :: xnew
    integer, optional, intent(in) :: obs_type_new
    
    integer :: obs_type_new1, i

    ! covariance vector of the desired input with each support point
    real(dp), dimension(this%m) :: km

    ! Treatment of the optional obs_type_new. Default: a value
    ! observation.
    if (present(obs_type_new)) then
       obs_type_new1 = obs_type_new
    else
       obs_type_new1 = 0
    endif

    do i=1,this%m
       km(i) = this%covariance%cov(obs_type_new1, this%obs_type(i), xnew, this%x(i,:), this%theta)
    end do
    
    ! equation (8.26) in ref. [1]
    predict = dot_product(km, this%invQKmn_t)

  end function predict

  function log_lik(this)
    class(SparseGP), intent(in) :: this
    real(dp) :: log_lik

    integer m,n
    m = this%m
    n = size(this%t)

    ! See equations (24--25) in ref. [2]
    log_lik = -0.5 * ((n-m)*log(this%nu(1)) + logdet(this%Q) + (dot_product(this%t, this%t) &
         - dot_product(this%Kmn_t, this%invQKmn_t))/this%nu(1))
  end function log_lik
end module m_GP_sparse
