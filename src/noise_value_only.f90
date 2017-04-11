module m_noise_value_only
  use m_noise
  use m_util
  implicit none
  
  private
  public noise_value_only

  type, extends(noise_model) :: noise_value_only
   contains
     procedure, nopass :: nparams_required
     procedure, nopass :: noise
  end type noise_value_only

contains
  pure function nparams_required(dims)
    integer, intent(in) :: dims
    integer nparams_required
    ! single noise level to be applied to the target value (and zero
    ! to the derivatives)
    nparams_required = 1
  end function nparams_required

  pure function noise(obs_type, params)
    integer, intent(in) :: obs_type
    real(dp), intent(in) :: params(:)
    real(dp) noise
    noise = merge(params(1), 0.0_dp, obs_type.eq.0)
  end function noise
  
end module m_noise_value_only
