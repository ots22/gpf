module m_noise_param2
  use m_noise
  use m_util
  implicit none

  private
  public noise_param2

  type, extends(noise_model) :: noise_param2
   contains
     procedure, nopass :: nparams_required
     procedure, nopass :: noise
  end type noise_param2

contains

  pure function nparams_required(dims)
    integer, intent(in) :: dims
    integer nparams_required
    ! value and derivative noise parameters
    nparams_required = 2
  end function nparams_required

  pure function noise(obs_type, params)
    integer, intent(in) :: obs_type
    real(dp), intent(in) :: params(:)
    real(dp) noise

    noise = 0

    if (obs_type.eq.0) then
       noise = params(1)
    else if (obs_type.ge.1.and.obs_type.le.6) then
       noise = params(2)
    end if
  end function noise

end module m_noise_param2
