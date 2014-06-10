program GP_test
  use m_GP
  use gp_example

  implicit none

  integer, parameter :: Nout = 1000
  integer :: i,u
  real(dp) :: a, startt, finisht

  call gp_example_initialise

  open(newunit=u, file="gp_test.out")

  call cpu_time(startt)
  do i=0,Nout
     a = 5.0_dp*i/Nout
     write (u,*) a, SparseGP_predict(gp, (/ a /), 0)
     !tmp = SparseGP_predict(gp, (/ a /), 1)
  end do
  call cpu_time(finisht)

  close(u)

  write (*,*) "prediction time: ", finisht-startt
end program GP_test
