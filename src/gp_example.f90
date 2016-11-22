! contains a simple example of a gp to use in tests

module m_gp_example
  use m_gp
  implicit none
  
  type(SparseGP) :: gp
  integer, parameter :: N = 50
  real(dp) :: x(2*N,1), t(2*N)
  integer :: obs_type(2*N)

  public :: gp, N, x, t, obs_type

contains

  subroutine gp_example_initialise   
    integer :: i, j, u
    real(dp) :: a

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

    gp = SparseGP(N-2, 1.d-4, (/ 1.4_dp /), x(1:N,:), obs_type(1:N), t(1:N))

  end subroutine gp_example_initialise
end module m_gp_example
