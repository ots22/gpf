! Check that gp objects can be read in and written back out 
! without change

program gp_io_test
  use m_gp
  use m_gp_example

  implicit none

  type(SparseGP) :: gp2
  type(DenseGP) :: gpDense2

  call gp_example_initialise

  call gp%write_out('square.gp')

  gp2 = read_SparseGP('square.gp')
  call gp2%write_out('square2.gp')

  call gpDense%write_out('square.gpd')
  gpDense2 = read_DenseGP('square.gpd')
  call gpDense2%write_out('square2.gpd')

end program gp_io_test
