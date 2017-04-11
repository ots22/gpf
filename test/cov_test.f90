program cov_test
  use m_util
  use m_cov_sqexp

  implicit none

  integer, parameter :: Ndims = 2
  real(kind=dp), dimension(Ndims) :: x1, y1, r1

  x1 = (/ 1.0_dp, 1.0_dp /)
  y1 = (/ 2.0_dp, 1.4_dp /)
  r1 = (/ 1.0_dp, 0.3_dp /)

  call cov_test1
  call cov_test2
  call cov_test3
  call cov_test4

contains
  ! expected values in the routines below are from finite differences
  subroutine cov_test1
    type(cov_sqexp) cf
    real(dp) :: actual, expected
    actual = cf%cov_val(x1, y1, r1)       
    expected = 0.24935220877729633_dp
    if (abs(actual - expected) > 1d-15) stop 1
  end subroutine cov_test1

  subroutine cov_test2
    type(cov_sqexp) cf
    real(kind=dp), dimension(Ndims) :: actual, expected
    integer :: n
    actual = (/ (cf%dcov_x1(n, x1, y1, r1), n=1,2) /)
    expected = (/ .2493522088_dp, 1.108232039_dp /)
    if (maxval(abs(actual - expected)) > 1d-10) stop 2
  end subroutine cov_test2

  subroutine cov_test3
    type(cov_sqexp) cf
    real(kind=dp), dimension(Ndims) :: actual, expected
    integer :: n
    actual = (/ (cf%dcov_x2(n, x1, y1, r1), n=1,2) /)
    expected = (/ -0.2493522088_dp, -1.108232039_dp /)
    if (maxval(abs(actual - expected)) > 1d-10) stop 3
  end subroutine cov_test3

  subroutine cov_test4
    type(cov_sqexp) cf
    real(kind=dp), dimension(Ndims*Ndims) :: actual, expected
    integer :: n, m
    actual = (/ ((cf%d2cov_xx(n, m, x1, y1, r1), n=1,2), m=1,2) /)
    expected = (/ 0.0_dp, -1.10823209_dp, -1.10823209_dp, -2.1548956_dp /)
    if (maxval(abs(actual - expected)) > 1d-7) stop 4    
  end subroutine cov_test4
  
end program cov_test
