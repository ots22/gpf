program solve_test
  use m_util
  
  implicit none
  
  call solve_test1

contains
  
  subroutine solve_test1
    real(dp) :: actual(3), expected(3)
    real(dp) :: M(3,3), b(3)

    M(1,:) = [ 3.0,  2.0, -1.0]
    M(2,:) = [ 2.0, -2.0,  4.0]
    M(3,:) = [-1.0,  0.5, -1.0]

    b = [1.0, -2.0, 0.0]

    actual = solve(M, b)
    expected = [1.0, -2.0, -2.0]

    if (maxval(abs(actual - expected)) > 1D-10) stop 1
  end subroutine solve_test1
end program solve_test
