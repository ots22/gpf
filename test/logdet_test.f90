program log_det_test
  use util
  implicit none
  
  real(dp), dimension(3,3) :: A
  real(dp) :: expected

  A = reshape((/ 1.0_dp, 5.0_dp, 20.0_dp, -2.0_dp, -0.5_dp, 8.0_dp, 7.0_dp, 1.0_dp, 0.0_dp /), (/ 3, 3 /))

  expected = 5.71042701737487_dp
  if (abs(logdet(A)-expected) > 1.0d-10) stop 1
  
end program log_det_test
