module amr
  use num_types
  use logger
  use p4est
  use mesh
  implicit none

  private
  public :: amr_refine_all

contains

  subroutine amr_refine_all(msh)
    ! argument list
    type(mesh_t), intent(inout) :: msh
    ! local variables
    integer(i4), allocatable, dimension(:) :: ref_mark
    integer(i4) :: level_max
    logical :: ifmod
    type(p4_mesh_restr_t) :: msh_rstr

    call neko_log%section("Mesh refinement")
    ! mark all the elements to be refined
    allocate(ref_mark(msh%nelv))
    ref_mark(:) = 1 ! THIS VALUE SHOULD BE TAKEN FORM p4est_wrap.h
    level_max = 3 ! THIS SHOULD BE TAKEN FROM RUNTIME PARAMETERS
    call p4_refine(ref_mark, level_max, ifmod, msh_rstr)

    ! Import mesh from p4est to neko
    if (ifmod) call p4_msh_get(msh)

    ! free memory
    deallocate(ref_mark)

    call neko_log%end_section()

    return
  end subroutine amr_refine_all

end module amr
