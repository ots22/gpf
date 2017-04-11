program gp_predict
  use m_util
  use m_gp
  use m_gp_sparse
  use m_gp_dense
  
  implicit none

  class(BaseGP), allocatable :: gp

  character(len=max_name_len) label
  integer u

  open(newunit=u, file='in.gp')
  read (u,'(A)') label
  close(u)
  
  select case (label)
  case('SparseGP')
     allocate(gp, source = SparseGP('in.gp'))
  case('DenseGP')
     allocate(gp, source = DenseGP('in.gp'))
  case default
     print *, "Incompatible data file (unrecognised GP type '", label, "')"
     stop 1
  end select

  call loop(size(gp%x,2))

contains

  subroutine loop(d)
    integer, intent(in) :: d
    real(dp) ::  x(d)
    integer obs_type
    do
       read (*,*) obs_type, x(:)
       write (*,*) obs_type, x(:), gp%predict(x, obs_type)
    end do
  end subroutine loop

end program gp_predict
