program optim_test
  use m_gp_optim
  use m_gp_example

  implicit none

  integer, parameter :: Nout = 1000
  integer :: i
  real(dp) :: a

  call gp_example_initialise

  call log_lik_optim(2, gp, (/ 0.000001_dp, 0.1_dp /), (/ 10.0_dp, 5.0_dp /), &
       20, 0.001_dp)

  do i=0,N
     write (8,*) x(i,1), t(i)
  end do
  
  do i=0,Nout
     a = 5.0_dp*i/Nout
     write (7,*) a, SparseGP_predict(gp, (/ a /), 0)
  end do

end program optim_test
