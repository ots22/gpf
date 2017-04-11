module m_cov_sqexp_param4
  use m_util
  use m_cov
  
  implicit none
  
  private
  public cov_sqexp_param4
  
  type, extends(cov_fn) :: cov_sqexp_param4
   contains
     procedure, nopass :: ntheta_required
     procedure, nopass :: cov_val
     procedure, nopass :: dcov_x1
     procedure, nopass :: dcov_x2
     procedure, nopass :: d2cov_xx
  end type cov_sqexp_param4

contains
  pure function ntheta_required(dims)
    integer ntheta_required
    integer, intent(in) :: dims
    ntheta_required = 4 ! scale parameter, and r parameters for diag, off-diag and energy
  end function ntheta_required

  pure function cov_val(x,y,hypers)
    real(dp) :: cov_val
    real(dp), intent(in), dimension(:) :: x, y, hypers
    real(dp) :: r(7), scale
    scale = hypers(1)
    r(1:3) = hypers(2)
    r(4:6) = hypers(3)
    r(7) = hypers(4)
    cov_val = scale * exp(-0.5*sum((x-y)**2/r**2))
  end function cov_val

  pure function dcov_x1(n,x,y,hypers)
    real(dp) :: dcov_x1
    real(dp), intent(in), dimension(:) :: x, y, hypers
    integer, intent(in) :: n
    real(dp) :: r(7)
    r(1:3) = hypers(2)
    r(4:6) = hypers(3)
    r(7) = hypers(4)
    dcov_x1 = -(x(n)-y(n))/r(n)**2 * cov_val(x,y,hypers)
  end function dcov_x1

  pure function dcov_x2(n,x,y,hypers)
    real(dp) :: dcov_x2
    real(dp), intent(in), dimension(:) :: x, y, hypers
    integer, intent(in) :: n
    real(dp) :: r(7)
    r(1:3) = hypers(2)
    r(4:6) = hypers(3)
    r(7) = hypers(4)
    dcov_x2 = (x(n)-y(n))/r(n)**2 * cov_val(x,y,hypers)
  end function dcov_x2

  pure function d2cov_xx(n,m,x,y,hypers)
    real(dp) :: d2cov_xx
    real(dp), intent(in), dimension(:) :: x, y, hypers
    integer, intent(in) :: n, m
    real(dp) :: r(7)
    r(1:3) = hypers(2)
    r(4:6) = hypers(3)
    r(7) = hypers(4)
    d2cov_xx = merge(cov_val(x,y,hypers)/r(n)**2, 0.0_dp, n.eq.m) &
         - (x(n)-y(n))/r(n)**2 * dcov_x2(m,x,y,hypers)
  end function d2cov_xx
end module m_cov_sqexp_param4
