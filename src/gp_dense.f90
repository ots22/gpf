module m_gp_dense
  use m_gp
  use m_util
  use m_cov_all
  use m_noise_all
  implicit none

  private
  public DenseGP, read_DenseGP

  type, extends(BaseGP) :: DenseGP
     ! covariance matrix and inverse
     real(dp), dimension(:,:), allocatable :: C, invC
     ! precomuted product, used in the prediction
     real(dp), dimension(:), allocatable :: invCt
   contains
     procedure log_lik
     procedure update_matrices
     procedure predict
     procedure write_out
  end type DenseGP

  interface DenseGP
     module procedure make_DenseGP
     module procedure read_DenseGP
  end interface DenseGP

contains

  subroutine alloc_DenseGP(gp, n, ntheta, dims, CovFunction, NoiseModel)
    type(DenseGP), intent(inout) :: gp
    integer, intent(in) :: n, ntheta, dims
    class(cov_fn) :: CovFunction
    class(noise_model) :: NoiseModel
    
    allocate(real(dp) :: gp%nu(NoiseModel%nparams_required(dims)))
    allocate(real(dp) :: gp%theta(ntheta))
    allocate(real(dp) :: gp%x(n, dims))
    allocate(integer :: gp%obs_type(n))
    allocate(real(dp) :: gp%t(n))
    allocate(real(dp) :: gp%C(n,n))
    allocate(real(dp) :: gp%invC(n,n))
    allocate(real(dp) :: gp%invCt(n))
    allocate(gp%covariance, mold=CovFunction)
    allocate(gp%noise_model, mold=NoiseModel)
  end subroutine alloc_DenseGP

  function make_DenseGP(nu, theta, x, obs_type, t, CovFunction, NoiseModel) result(gp)
    type(DenseGP) :: gp
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

    call alloc_DenseGP(gp, n, size(theta), d, CovFunction, NoiseModel)

    gp%nu = nu
    gp%theta = theta
    gp%obs_type = obs_type
    gp%t = t
    gp%x = x
    call update_matrices(gp)    
  end function make_DenseGP

  subroutine write_out(this, filename)
    class(DenseGP), intent(in) :: this
    character(len=*), intent(in) :: filename
    character(len=max_name_len) :: cov_fn_name
    character(len=max_name_len) :: noise_model_name
    integer :: u

    cov_fn_name = cov_fn_to_string(this%covariance)
    noise_model_name = noise_model_to_string(this%noise_model)

    open(newunit=u, file=filename)

    write (u,'(A)') "DenseGP"
    write (u,'(I10)') size(this%t), size(this%theta), size(this%nu), size(this%x,2)
    write (u,'(A)') trim(cov_fn_name)
    write (u,'(A)') trim(noise_model_name)
    write (u,'(es24.15)') this%nu, this%theta, this%x
    write (u,'(I4)') this%obs_type
    write (u,'(es24.15)') this%C, this%invC, this%invCt

    close(u)
  end subroutine write_out

  function read_DenseGP(filename) result(gp)
    character(len=*), intent(in) :: filename
    type(DenseGP) :: gp
    integer n, ntheta, nnu, d, u
    character(len=max_name_len) :: label
    character(len=max_name_len) :: cov_fn_name
    character(len=max_name_len) :: noise_model_name
    class(cov_fn), allocatable :: cf
    class(noise_model), allocatable :: nm
    open(newunit=u, file=filename)
    read (u,'(A)') label

    if (trim(label) /= "DenseGP") then
       print *, "Incompatible data file"
       stop 1
    end if

    read (u,'(I10)') n, ntheta, nnu, d

    read (u,'(A)') cov_fn_name 
    call string_to_cov_fn(cov_fn_name, cf)
    if (ntheta /= cf%ntheta_required(d)) then
       print *, "ntheta does not match number required by the covariance function"
       stop 1
    end if

    read (u,'(A)') noise_model_name
    call string_to_noise_model(noise_model_name, nm)
    if (nnu /= nm%nparams_required(d)) then
       print *, "size of nu (noise params) does not match number required by the noise model"
       stop 1
    end if

    call alloc_DenseGP(gp, n, ntheta, d, cf, nm)
    
    read (u,'(es24.15)') gp%nu, gp%theta, & 
         gp%x
    read (u,'(I4)') gp%obs_type
    read (u,'(es24.15)') gp%C, gp%invC, gp%invCt
    close(u)
  end function read_DenseGP

  subroutine update_matrices(this)
    class(DenseGP), intent(inout) :: this
    integer :: i,j,n
    real(dp) noise
    ! always a small amount of noise to stabilize the inversion
    real(dp), parameter :: noise_stab = 1e-9_dp

    n = size(this%t)
    
    do i=1,n
       do j=1,n
          if (i.ne.j) then
             noise = 0.0_dp
          else
             noise = this%noise_model%noise(this%obs_type(i), this%nu) + noise_stab
          end if
          this%C(i,j) = this%covariance%cov(this%obs_type(i), this%obs_type(j), &
               this%x(i,:), this%x(j,:), this%theta) + noise
       end do
    end do
    
    ! replace the numerical inverse with a solve
    this%invC = this%C
    call ninv(this%invC)
    this%invCt = matmul(this%invC, this%t)
  end subroutine update_matrices

  function predict(this, xnew, obs_type_new)
    real(dp) predict
    class(DenseGP), intent(in) :: this
    real(dp), dimension(:), intent(in) :: xnew
    integer, optional, intent(in) :: obs_type_new
    
    integer :: obs_type_new1, i
    
    real(dp), dimension(size(this%t)) :: k
    if (present(obs_type_new)) then
       obs_type_new1 = obs_type_new
    else 
       obs_type_new1 = 0
    endif

    do i=1,size(this%t)
       k(i) = this%covariance%cov(obs_type_new1,this%obs_type(i),xnew,this%x(i,:),this%theta)
    end do

!    print *, obs_type_new1, this%obs_type(10), xnew, this%x(10,:), this%theta, k(10)

    predict = dot_product(k, this%invCt)
  end function predict

  function log_lik(this)
    class(DenseGP), intent(in) :: this
    real(dp) log_lik
    log_lik = -0.5_dp * (logdet(this%C) + dot_product(this%t, this%invCt))
  end function log_lik

end module m_gp_dense
