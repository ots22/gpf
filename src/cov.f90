module m_cov_sqexp
use m_util
implicit none
contains
  pure function cov_val(x,y,r)
    real(dp) :: cov_val
    real(dp), intent(in), dimension(:) :: x, y, r
    cov_val = exp(-0.5*sum((x-y)**2/r**2))
  end function cov_val

  pure function dcov_x1(n,x,y,r)
    real(dp) :: dcov_x1
    real(dp), intent(in), dimension(:) :: x, y, r
    integer, intent(in) :: n
    dcov_x1 = -(x(n)-y(n))/r(n)**2 * cov_val(x,y,r)
  end function dcov_x1

  pure function dcov_x2(n,x,y,r)
    real(dp) :: dcov_x2
    real(dp), intent(in), dimension(:) :: x, y, r
    integer, intent(in) :: n
    dcov_x2 = (x(n)-y(n))/r(n)**2 * cov_val(x,y,r)
  end function dcov_x2

  pure function d2cov_xx(n,m,x,y,r)
    real(dp) :: d2cov_xx
    real(dp), intent(in), dimension(:) :: x, y, r
    integer, intent(in) :: n, m
    d2cov_xx = merge(cov_val(x,y,r)/r(n)**2, 0.0_dp, n.eq.m) &
         - (x(n)-y(n))/r(n)**2 * dcov_x2(m,x,y,r)
  end function d2cov_xx

  pure function cov(n,m,x,y,r)
    real(dp) :: cov
    real(dp), intent(in), dimension(:) :: x, y, r
    integer, intent(in) :: n, m
    if (n.eq.0 .and. m.eq.0) then
       cov = cov_val(x,y,r)
    else if (n.eq.0) then
       cov = dcov_x2(m,x,y,r)
    else if (m.eq.0) then
       cov = dcov_x1(n,x,y,r)
    else
       cov = d2cov_xx(n,m,x,y,r)
    end if
  end function cov

end module m_cov_sqexp
