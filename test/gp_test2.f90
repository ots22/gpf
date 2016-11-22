program gp_test2
  use m_gp
  use m_gp_example

  implicit none

  type(SparseGP) :: gp2

  call gp_example_initialise

  call write_SparseGP('square.gp',gp)
  gp2 = read_SparseGP('square.gp')
  call write_SparseGP('square2.gp',gp2)

end program gp_test2
