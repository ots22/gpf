module m_noise_all
  use m_noise
  use m_noise_value_only
  use m_noise_param2

  use m_util, only: max_name_len

contains
  
  function noise_model_to_string(NoiseModel) result(noise_model_name)
    class(noise_model), intent(in) :: NoiseModel
    character(len=max_name_len) noise_model_name

    select type(nm => NoiseModel)
    type is (noise_value_only)
       noise_model_name = 'VAL'
    type is (noise_param2)
       noise_model_name = 'PARAM2'
    class default
       noise_model_name = 'UNKNOWN'
    end select
  end function noise_model_to_string

  subroutine string_to_noise_model(noise_model_name, NoiseModel)
    character(len=max_name_len), intent(in) :: noise_model_name
    class(noise_model), intent(out), allocatable :: NoiseModel

    select case (noise_model_name)
    case ('VAL')
       allocate(noise_value_only :: NoiseModel)
    case ('PARAM2')
       allocate(noise_param2 :: NoiseModel)
    case default
       print *, "unknown noise model, ", noise_model_name
       stop 1
    end select
  end subroutine string_to_noise_model

end module m_noise_all
