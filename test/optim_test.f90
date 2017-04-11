program optim_test
  use m_gp_optim
  use m_gp_example

  implicit none

  integer, parameter :: Nout = 1000
  integer :: i
  real(dp) :: a

  call gp_example_initialise

  call log_lik_optim(2, gpDense, [1e-10_DP, 0.1_dp], [10.0_dp, 5.0_dp], &
       100, 0.0001_dp)

  do i=1,N
     write (8,*) x(i,1), t(i)
  end do
  
  do i=0,Nout
     a = 5.0_dp*i/Nout
     write (7,*) a, gpDense%predict([a], 0)
  end do

end program optim_test
