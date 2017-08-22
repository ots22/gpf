! contains a simple example of a gp to use in tests

module m_gp_example
  use m_gp
  use m_gp_sparse
  use m_gp_dense
  use m_util
  use m_cov_sqexp
  use m_noise_value_only
  implicit none
  
  type(SparseGP) :: gp
  type(DenseGP) :: gpdense
  integer, parameter :: N = 50
  real(dp) :: x(2*N,1), t(2*N)
  integer :: obs_type(2*N)

  public :: gp, gpDense, N, x, t, obs_type

contains

  subroutine gp_example_initialise   
    integer :: i, j, u
    real(dp) :: a
    type(cov_sqexp) :: cf
    type(noise_value_only) :: nm

    ! shuffled input vector
    do i=1,N
       call random_number(a)
       j = int(a*i)+1 ! random integer in [1,i]
       if (j.ne.i) x(i,:) = x(j,:)
       x(j,1) = 5.0*(i-1.0)/(N-1)
    end do

    do i=1,N
       ! value
       call random_number(a)
       t(i) = x(i,1)**2 ! + 0.0001 * a
       obs_type(i) = 0
       ! deriv
       x(i+N,1) = x(i,1)
       t(i+N) = 2*x(i,1)
       obs_type(i+N) = 1
    end do

    open(newunit=u,file="gp_example_training_data")
    write (u,'(2F20.10)') (x(i,1), t(i), i=1,N)
    close(u)

    gp = SparseGP(N-2, (/ 1.d-4 /), (/ 1.0_dp, 1.4_dp /), x(1:N,:), obs_type(1:N), t(1:N), cf, nm)

    gpDense = DenseGP( (/ 1.d-9 /), (/ 1.0_dp, 1.4_dp /), x(1:N,:), obs_type(1:N), t(1:N), cf, nm)

  end subroutine gp_example_initialise
end module m_gp_example
