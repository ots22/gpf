module c_interface
  use m_gp
  use m_gp_sparse 
  use iso_c_binding
  implicit none

  private
  public :: gp_projected_process_read, &
       & gp_projected_process_predict, &
       & gp_projected_process_destroy
  
  ! An approach based on passing around C pointers to the SparseGP
  ! object returned by C_LOC does not work, since SparseGP is not a
  ! c-interoperable type as defined by the standard (and does not
  ! compile with gfortran-6, despite it working with earlier versions
  ! of the compiler).  Instead, keep a simple pool of pointers, to
  ! which the `gpf_handle' type is an index.  This is rather
  ! simplistic, and only allows for MAX_GP_HANDLE GPs to be used over
  ! the entire lifetime of the program.
  integer, parameter :: MAX_GP_HANDLE=1024
  type SparseGP_ptr
     type(SparseGP), pointer :: ptr
  end type SparseGP_ptr
  type(SparseGP_ptr) :: gp(MAX_GP_HANDLE)

contains
  subroutine gp_projected_process_read(handle, fname_len, fname) bind(c)
    integer(c_int), intent(out) :: handle
    integer(c_int), intent(in), value :: fname_len
    character(kind=c_char), intent(in) :: fname
    integer(c_int), save :: current_gp_handle = 1
    
    if (current_gp_handle > MAX_GP_HANDLE) then
       write (*,*) "Error: out of space to store GP handles. &
            & Consider recompiling after increasing the value of MAX_GP_HANDLE."
    end if

    allocate(gp(current_gp_handle)%ptr)
    handle = current_gp_handle
    current_gp_handle = current_gp_handle + 1
    gp(current_gp_handle-1)%ptr = read_SparseGP(fname(1:fname_len))
  end subroutine gp_projected_process_read

  function gp_projected_process_predict(handle, input_len, input, input_type) result(output) &
       bind(c)
    real(c_double) :: output
    integer(c_int), intent(in), value :: handle
    integer(c_int), intent(in), value :: input_len
    real(c_double), intent(in) :: input(input_len)
    integer(c_int), intent(in), value :: input_type
    output = gp(handle)%ptr%predict(input, input_type)
  end function gp_projected_process_predict

  subroutine gp_projected_process_destroy(handle) bind(c)
    integer(c_int), intent(in), value :: handle
    deallocate(gp(handle)%ptr)
  end subroutine gp_projected_process_destroy
  
end module c_interface
