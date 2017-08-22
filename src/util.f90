module m_util
  implicit none
  integer, parameter ::  dp = selected_real_kind(14,100)

  integer, parameter :: max_name_len=100

contains
  ! Solve a system of linear equations A*x = b for a matrix A and
  ! vectors x and b.  Implemented as a wrapper around LAPACK dgesv.
  function solve(M,b) result(x)
    real(dp), intent(in) :: M(:,:)
    real(dp) :: b(:)
    real(dp), dimension(size(M,1),size(M,2)) :: A
    real(dp), dimension(size(b,1)) :: x
    integer ipiv(size(b,1)), N, info

    ! check M is square and b conforms
    if (size(M,1) /= size(M,2)) then
       write (0,*) "solve: Array M passed to solve should be square.  Got ", &
            size(M,1), "x", size(M,2)
       error stop
    end if

    N = size(M, 1)

    ! A and x are overwritten on output of dgesv, so copy
    A = M
    x = b
    call dgesv(N, 1, A, N, ipiv, x, N, info)

    ! check for success
    if (info /= 0) then
       write (0,*) "solve: dgesv returned an error code (", info, ")"
       error stop
    end if

  end function solve

  ! in place inverse
  subroutine ninv(A)
    real(dp), dimension(:,:), intent(inout) :: A
    integer, dimension(:), allocatable :: piv
    real(dp), dimension(1) :: lwork_real
    real(dp), dimension(:), allocatable :: work
    integer :: lwork, N, info

    N = size(A,1)

    allocate(integer :: piv(N))

    call dgetrf(N, N, A, N, piv, info)
    
    call dgetri(N, A, N, piv, lwork_real, -1, info)
    lwork = nint(lwork_real(1))
    allocate(real(dp) :: work(lwork))
    call dgetri(N, A, N, piv, work, lwork, info)
    
  end subroutine ninv

  function logdet(A)
    real(dp) :: logdet
    real(dp), dimension(:,:), intent(in) :: A
    real(dp), dimension(:,:), allocatable :: tmp
    integer, dimension(:), allocatable :: ipiv
    integer :: N, info, i

    N = size(A,1)

    allocate(real(dp) :: tmp(N, N))
    allocate(integer :: ipiv(N))

    tmp = A
    call dgetrf(N, N, tmp, N, ipiv, info)

    logdet = 0.0_dp
    do i=1,N
       logdet = logdet + log(abs(tmp(i,i)))
    end do

  end function logdet
  
end module m_util
