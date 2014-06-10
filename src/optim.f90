module gp_optim
  use m_gp
  implicit none
  include 'nlopt.f'

contains
  subroutine log_lik_optim(ntheta, gp, lbounds, ubounds, iterations, ftol_rel)
    type(SparseGP) :: gp
    integer(kind=8) :: opt
    integer :: ntheta, iterations, ires
    real(dp) :: ftol_rel, maxf
    real(dp), dimension(ntheta) :: lbounds, ubounds
    real(dp), dimension(ntheta) :: hypers

    hypers(1) = gp%nu
    hypers(2:) = gp%theta

    call nlo_create(opt, NLOPT_LN_BOBYQA, size(lbounds))
    call nlo_set_lower_bounds(ires, opt, lbounds)
    call     check_error_code(ires)
    call nlo_set_upper_bounds(ires, opt, ubounds)
    call     check_error_code(ires)
    call nlo_set_max_objective(ires, opt, nlog_lik, gp)
    call     check_error_code(ires)
    call nlo_set_ftol_rel(ires, opt, ftol_rel)
    call     check_error_code(ires)
    call nlo_optimize(ires, opt, hypers, maxf)
    call     check_error_code(ires)
    call nlo_destroy(opt)
  end subroutine log_lik_optim

  subroutine check_error_code(ires)
    integer ires
    if (ires.le.0) then
       print *, "NLopt failed with error code ", ires
       stop 1
    end if   
  end subroutine check_error_code

end module gp_optim
