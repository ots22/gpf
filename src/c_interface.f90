module c_interface
  use m_gp
  use iso_c_binding
  implicit none
contains
  subroutine gp_projected_process_read(handle, fname_len, fname) bind(c)
    type(c_ptr), intent(out) :: handle
    integer(c_int), intent(in), value :: fname_len
    character(kind=c_char), intent(in) :: fname
    type(SparseGP), pointer :: gp
    allocate(gp)
    gp = read_SparseGP(fname(1:fname_len))
    handle = c_loc(gp)
  end subroutine gp_projected_process_read

  subroutine gp_projected_process_predict(output, handle, input_len, input, input_type) &
       bind(c)
    real(c_double), intent(out) :: output
    type(c_ptr), intent(in) :: handle
    integer(c_int), intent(in), value :: input_len
    real(c_double), intent(in) :: input(input_len)
    integer(c_int), intent(in), value :: input_type
    type(SparseGP), pointer :: gp
    call c_f_pointer(handle, gp)
    output = SparseGP_predict(gp, input, input_type)   
  end subroutine gp_projected_process_predict
  
end module c_interface
