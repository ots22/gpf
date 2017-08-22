program inv_test
  use m_util
  integer, dimension(3) :: piv
  real(dp), dimension(3,3) :: A, invA_expected
  real(dp), allocatable :: work(:)
  real(dp) :: invA(3,3), invA_lapack(3,3), lwork_real(1)
  integer lwork

  A = reshape((/ 1.0, 2.0, 3.0, 5.0, 7.0, 9.0, 11.0, -1.0, -0.5 /), shape(A))

  invA_expected(1,:) = (/ -11.0, -203.0, 164.0 /)
  invA_expected(2,:) = (/   4.0,   67.0, -46.0 /)
  invA_expected(3,:) = (/   6.0,  -12.0,   6.0 /)
  invA_expected = invA_expected / 75

  invA_lapack = A
  call dgetrf(3, 3, invA_lapack, 3, piv, info)
  
  call dgetri(3, invA_lapack, 3, piv, lwork_real, -1, info)
  lwork = nint(lwork_real(1))
  allocate(real(dp) :: work(lwork))
  call dgetri(3, invA_lapack, 3, piv, work, lwork, info)
  
  invA = A
  call ninv(invA)
  
  if (maxval(abs(invA - invA_lapack)) > 1D-10) stop 1
  if (maxval(abs(invA - invA_expected)) > 1D-10) stop 2
end program inv_test
