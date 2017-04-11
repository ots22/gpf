module m_noise
  use m_util
  implicit none
  
  private
  public noise_model
  
  type, abstract :: noise_model
   contains
     procedure(nparams_required), deferred, nopass :: nparams_required
     ! returns a vector of values for the value and each dimension of the gradient
     procedure(noise), deferred, nopass :: noise
  end type noise_model

  abstract interface
     pure function nparams_required(dims)
       import noise_model
       integer, intent(in) :: dims
       integer nparams_required       
     end function nparams_required

     pure function noise(obs_type, params)
       use m_util, only: dp
       import noise_model
       real(dp) noise
       integer, intent(in) :: obs_type
       real(dp), intent(in) :: params(:)
     end function noise
  end interface

end module m_noise
