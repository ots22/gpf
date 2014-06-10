program inv_test
  use util
  integer, dimension(3) :: piv
  real(dp), dimension(3,3) :: A, invA
  real(dp), dimension(1) :: lwork_real
  real(dp), dimension(:), allocatable :: work
  integer :: lwork

  A = reshape((/ 1.0, 2.0, 3.0, 5.0, 7.0, 9.0, 11.0, -1.0, -0.5 /), shape(A))

  invA = A
  call dgetrf(3, 3, invA, 3, piv, info)
  
  call dgetri(3, invA, 3, piv, lwork_real, -1, info)
  lwork = nint(lwork_real(1))
  print *, lwork
  allocate(real(dp) :: work(lwork))
  call dgetri(3, invA, 3, piv, work, lwork, info)

  print *, invA

  invA = A
  
  call ninv(invA)

  print *, invA


end program inv_test
