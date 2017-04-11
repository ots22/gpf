module m_cov
use m_util
implicit none

private
public cov_fn

  type, abstract :: cov_fn
   contains
     procedure(ntheta_required), deferred, nopass :: ntheta_required
     procedure(cov_val), deferred, nopass :: cov_val
     procedure(dcov_x1), deferred, nopass :: dcov_x1
     procedure(dcov_x2), deferred, nopass :: dcov_x2
     procedure(d2cov_xx), deferred, nopass :: d2cov_xx
     procedure :: cov
  end type cov_fn

  abstract interface
     pure function ntheta_required(dims)
       import cov_fn
       integer :: ntheta_required
       integer, intent(in) :: dims
     end function ntheta_required

     pure function cov_val(x,y,hypers)
       use m_util, only: dp
       import cov_fn
       real(dp) :: cov_val
       real(dp), intent(in), dimension(:) :: x, y, hypers
     end function cov_val
     
     pure function dcov_x1(n,x,y,hypers)
       use m_util, only: dp
       import cov_fn
       real(dp) :: dcov_x1
       real(dp), intent(in), dimension(:) :: x, y, hypers
       integer, intent(in) :: n
     end function dcov_x1
     
     pure function dcov_x2(n,x,y,hypers)
       use m_util, only: dp
       import cov_fn
       real(dp) :: dcov_x2
       real(dp), intent(in), dimension(:) :: x, y, hypers
       integer, intent(in) :: n
     end function dcov_x2
     
     pure function d2cov_xx(n,m,x,y,hypers)
       use m_util, only: dp
       import cov_fn
       real(dp) :: d2cov_xx
       real(dp), intent(in), dimension(:) :: x, y, hypers
       integer, intent(in) :: n, m
     end function d2cov_xx
  end interface

contains

  pure function cov(cf,n,m,x,y,hypers)
    class(cov_fn), intent(in) :: cf
    real(dp) :: cov
    real(dp), intent(in), dimension(:) :: x, y, hypers
    integer, intent(in) :: n, m
    if (n.eq.0 .and. m.eq.0) then
       cov = cf%cov_val(x,y,hypers)
    else if (n.eq.0) then
       cov = cf%dcov_x2(m,x,y,hypers)
    else if (m.eq.0) then
       cov = cf%dcov_x1(n,x,y,hypers)
    else
       cov = cf%d2cov_xx(n,m,x,y,hypers)
    end if
  end function cov

end module m_cov
