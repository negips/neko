module amr
  use mpi_f08
  use num_types
  use comm
  use logger
  use parameters
  use mxm_wrapper
  use speclib
  use p4est
  use mesh
  use fluid_method
  use fluid_plan1
  use fluid_plan4
  use fluid_pnpn
  implicit none

  private
  public :: amr_refine, amr_rcn_init, amr_rcn_free

  ! type for fields reconstruction (refinement/transfer/coarsening)
  type amr_rcn_t
     ! mesh
     type(mesh_t), pointer :: msh => null()
     ! function space
     type(space_t), pointer :: Xh => null()

     ! interpolation operators
     ! coarse to fine
     real(dp), allocatable, dimension(:,:,:) :: x_cr_to_fn, x_cr_to_fnT
     real(dp), allocatable, dimension(:,:,:) :: y_cr_to_fn, y_cr_to_fnT
     real(dp), allocatable, dimension(:,:,:) :: z_cr_to_fn, z_cr_to_fnT
     ! fine to coarse
     real(dp), allocatable, dimension(:,:,:) :: x_fn_to_cr, x_fn_to_crT
     real(dp), allocatable, dimension(:,:,:) :: y_fn_to_cr, y_fn_to_crT
     real(dp), allocatable, dimension(:,:,:) :: z_fn_to_cr, z_fn_to_crT
     ! point multiplicity; use for coarsening after children face summation
     ! for now I assume lx=ly=lz; in genera there shold be 3 face and edge arrays
     ! element
     real(dp), allocatable, dimension(:,:,:) :: el_mult
     ! face
     real(dp), allocatable, dimension(:,:) :: fc_mult
     ! edge
     real(dp), allocatable, dimension(:) :: ed_mult

     ! p4est <=> neko distribution mapping
     type(p4_msh_trs_t) :: msh_trs

     ! element reconstruction data
     type(p4_msh_rcn_t) :: msh_rcn

     ! work space
     real(dp), allocatable, dimension(:,:,:,:) :: tmp
  end type amr_rcn_t

  type(amr_rcn_t), save :: amr_rcn

  interface amr_refine
     module procedure amr_msh_refine, amr_fluid_refine
  end interface amr_refine

contains

  !> Initialise interpolation arrays for refinement/coarsening
  subroutine amr_rcn_init(msh, Xh)
    ! argument list
    type(mesh_t), target, intent(in) :: msh
    type(space_t), target, intent(in) :: Xh
    ! local variables
    integer(i4) :: il, jl, kl, nt2
    real(dp), allocatable, dimension(:) :: tmpl

    call amr_rcn_free(amr_rcn)

    ! for now GLL only
    if (Xh%t /= 1) call neko_error('Only GLL space supported for AMR for now.')

    ! which mesh space I point to
    amr_rcn%msh => msh
    amr_rcn%Xh => Xh

    ! work array
    allocate(tmpl(max(Xh%lx, Xh%ly, Xh%lz)))

    ! interpolation operators
    ! x-direction
    allocate(amr_rcn%x_cr_to_fn(Xh%lx, Xh%lx, 2), amr_rcn%x_cr_to_fnT(Xh%lx, Xh%lx, 2), &
         & amr_rcn%x_fn_to_cr(Xh%lx, Xh%lx, 2), amr_rcn%x_fn_to_crT(Xh%lx, Xh%lx, 2))
    amr_rcn%x_cr_to_fn(:,:,:) = 0.0_dp
    amr_rcn%x_cr_to_fnT(:,:,:) = 0.0_dp
    amr_rcn%x_fn_to_cr(:,:,:) = 0.0_dp
    amr_rcn%x_fn_to_crT(:,:,:) = 0.0_dp
    ! coarse -> fine
    ! negative
    do il = 1, Xh%lx
       tmpl(il) = 0.5_dp*(Xh%zg(il,1) - 1.0_dp)
    end do
    call igllm(amr_rcn%x_cr_to_fn, amr_rcn%x_cr_to_fnT, Xh%zg(:, 1), &
         &tmpl, Xh%lx, Xh%lx, Xh%lx, Xh%lx)
    ! positive; we use symmetry
    do jl = 1, Xh%lx
       do il = 1, Xh%lx
          amr_rcn%x_cr_to_fn(Xh%lx-il+1, Xh%lx-jl+1, 2) = amr_rcn%x_cr_to_fn(il, jl, 1)
          amr_rcn%x_cr_to_fnT(Xh%lx-il+1, Xh%lx-jl+1, 2) = amr_rcn%x_cr_to_fnT(il, jl, 1)
       end do
    end do
    ! fine -> coarse
    ! negative
    nt2 = Xh%lx/2 + mod(Xh%lx, 2)
    do il=1, nt2
       tmpl(il) = 2.0_dp*Xh%zg(il, 1) + 1.0_dp
    end do
    call igllm(amr_rcn%x_fn_to_cr, amr_rcn%x_fn_to_crT, Xh%zg(:, 1), &
         &tmpl, Xh%lx, nt2, Xh%lx, Xh%lx)
    ! positive; we use symmetry
    do jl = 1, Xh%lx
       do il = 1, nt2
          amr_rcn%x_fn_to_cr(Xh%lx-il+1, Xh%lx-jl+1, 2) = amr_rcn%x_fn_to_cr(il, jl, 1)
          amr_rcn%x_fn_to_crT(Xh%lx-il+1, Xh%lx-jl+1, 2) = amr_rcn%x_fn_to_crT(il, jl, 1)
       end do
    end do

    ! y-direction
    allocate(amr_rcn%y_cr_to_fn(Xh%ly, Xh%ly, 2), amr_rcn%y_cr_to_fnT(Xh%ly, Xh%ly, 2), &
         & amr_rcn%y_fn_to_cr(Xh%ly, Xh%ly, 2), amr_rcn%y_fn_to_crT(Xh%ly, Xh%ly, 2))
    amr_rcn%y_cr_to_fn(:,:,:) = 0.0_dp
    amr_rcn%y_cr_to_fnT(:,:,:) = 0.0_dp
    amr_rcn%y_fn_to_cr(:,:,:) = 0.0_dp
    amr_rcn%y_fn_to_crT(:,:,:) = 0.0_dp
    ! coarse -> fine
    ! negative
    do il = 1, Xh%ly
       tmpl(il) = 0.5_dp*(Xh%zg(il, 2) - 1.0_dp)
    end do
    call igllm(amr_rcn%y_cr_to_fn, amr_rcn%y_cr_to_fnT, Xh%zg(:, 2), &
         &tmpl, Xh%ly, Xh%ly, Xh%ly, Xh%ly)
    ! positive; we use symmetry
    do jl = 1, Xh%ly
       do il = 1, Xh%ly
          amr_rcn%y_cr_to_fn(Xh%ly-il+1, Xh%ly-jl+1, 2) = amr_rcn%y_cr_to_fn(il, jl, 1)
          amr_rcn%y_cr_to_fnT(Xh%ly-il+1, Xh%ly-jl+1, 2) = amr_rcn%y_cr_to_fnT(il, jl, 1)
       end do
    end do
    ! fine -> coarse
    ! negative
    nt2 = Xh%ly/2 + mod(Xh%ly, 2)
    do il = 1, nt2
       tmpl(il) = 2.0_dp*Xh%zg(il, 2) + 1.0_dp
    end do
    call igllm(amr_rcn%y_fn_to_cr, amr_rcn%y_fn_to_crT, Xh%zg(:, 2), &
         &tmpl, Xh%ly, nt2, Xh%ly, Xh%ly)
    ! positive; we use symmetry
    do jl = 1, Xh%ly
       do il = 1, nt2
          amr_rcn%y_fn_to_cr(Xh%ly-il+1, Xh%ly-jl+1, 2) = amr_rcn%y_fn_to_cr(il, jl, 1)
          amr_rcn%y_fn_to_crT(Xh%ly-il+1, Xh%ly-jl+1, 2) = amr_rcn%y_fn_to_crT(il, jl, 1)
       end do
    end do

    ! z-direction; IS THIS CORRECT? SHOULD I ALLOCATE THEM ANYHOW AND FILL WITH 1?
    if (msh%gdim == 3) then
       allocate(amr_rcn%z_cr_to_fn(Xh%lz, Xh%lz, 2), amr_rcn%z_cr_to_fnT(Xh%lz, Xh%lz, 2), &
            & amr_rcn%z_fn_to_cr(Xh%lz, Xh%lz, 2), amr_rcn%z_fn_to_crT(Xh%lz, Xh%lz, 2))
       amr_rcn%z_cr_to_fn(:,:,:) = 0.0_dp
       amr_rcn%z_cr_to_fnT(:,:,:) = 0.0_dp
       amr_rcn%z_fn_to_cr(:,:,:) = 0.0_dp
       amr_rcn%z_fn_to_crT(:,:,:) = 0.0_dp
       ! coarse -> fine
       ! negative
       do il = 1, Xh%lz
          tmpl(il) = 0.5_dp*(Xh%zg(il, 3) - 1.0_dp)
       end do
       call igllm(amr_rcn%z_cr_to_fn, amr_rcn%z_cr_to_fnT, Xh%zg(:, 3), &
            &tmpl, Xh%lz, Xh%lz, Xh%lz, Xh%lz)
       ! positive; we use symmetry
       do jl = 1, Xh%lz
          do il = 1, Xh%lz
             amr_rcn%z_cr_to_fn(Xh%lz-il+1, Xh%lz-jl+1, 2) = amr_rcn%z_cr_to_fn(il, jl, 1)
             amr_rcn%z_cr_to_fnT(Xh%lz-il+1, Xh%lz-jl+1, 2) = amr_rcn%z_cr_to_fnT(il, jl, 1)
          end do
       end do
       ! fine -> coarse
       ! negative
       nt2 = Xh%lz/2 + mod(Xh%lz, 2)
       do il = 1, nt2
          tmpl(il) = 2.0_dp*Xh%zg(il, 3) + 1.0_dp
       end do
       call igllm(amr_rcn%z_fn_to_cr, amr_rcn%z_fn_to_crT, Xh%zg(:, 3), &
            &tmpl, Xh%lz, nt2, Xh%lz, Xh%lz)
       ! positive; we use symmetry
       do jl = 1, Xh%lz
          do il = 1, nt2
             amr_rcn%z_fn_to_cr(Xh%lz-il+1, Xh%lz-jl+1, 2) = amr_rcn%z_fn_to_cr(il, jl, 1)
             amr_rcn%z_fn_to_crT(Xh%lz-il+1, Xh%lz-jl+1, 2) = amr_rcn%z_fn_to_crT(il, jl, 1)
          end do
       end do
    end if

    ! multiplicity arrays
    if ((Xh%lx /= Xh%ly).or.((Xh%lz /= 1).and.(Xh%lz /= Xh%lx))) &
         & call neko_error('Different polynomial ordes in various directions are not supported for AMR')
    allocate(amr_rcn%el_mult(Xh%lx, Xh%ly, Xh%lz), amr_rcn%fc_mult(Xh%lx, Xh%lx), &
         & amr_rcn%ed_mult(Xh%lx))
    amr_rcn%el_mult(:,:,:) = 1.0_dp
    amr_rcn%fc_mult(:,:) = 1.0_dp
    amr_rcn%ed_mult(:) = 1.0_dp
    ! X
    if (mod(Xh%lx, 2) == 1) then
       il = Xh%lx/2 + 1
       do kl = 1, Xh%lz
          do jl = 1, Xh%ly
             amr_rcn%el_mult(il, jl, kl) = amr_rcn%el_mult(il, jl, kl) + 1.0_dp
          end do
       end do
    end if
    ! Y
    if (mod(Xh%ly, 2) == 1) then
       jl = Xh%ly/2 + 1
       do kl = 1, Xh%lz
          do il = 1, Xh%lx
             amr_rcn%el_mult(il, jl, kl) = amr_rcn%el_mult(il, jl, kl) + 1.0_dp
          end do
       end do
       if (mod(Xh%lx, 2) == 1) then
          il = Xh%lx/2 + 1
          do kl = 1, Xh%lz
             amr_rcn%el_mult(il, jl, kl) = amr_rcn%el_mult(il, jl, kl) + 1.0_dp
          end do
       end if
    end if

    ! Z
    if (msh%gdim == 3) then
       if (mod(Xh%lz, 2) == 1) then
          kl = Xh%lz/2 + 1
          do jl = 1, Xh%ly
             do il = 1, Xh%lx
                amr_rcn%el_mult(il, jl, kl) = amr_rcn%el_mult(il, jl, kl) + 1.0_dp
             end do
          end do
          if (mod(Xh%lx, 2) == 1) then
             il = Xh%lx/2 + 1
             do jl = 1, XH%ly
                amr_rcn%el_mult(il, jl, kl) = amr_rcn%el_mult(il, jl, kl) + 1.0_dp
             end do
          end if
          if (mod(Xh%ly, 2) == 1) then
             jl = Xh%ly/2 + 1
             do il = 1, Xh%lx
                amr_rcn%el_mult(il, jl, kl) = amr_rcn%el_mult(il, jl, kl) + 1.0_dp
             end do
          end if
          if ((mod(Xh%lx, 2) == 1).and.(mod(Xh%ly,2 ) == 1)) then
             il = Xh%lx/2 + 1
             jl = Xh%ly/2 + 1
             amr_rcn%el_mult(il, jl, kl) = amr_rcn%el_mult(il, jl, kl) + 1.0_dp
          end if
       end if
    end if
    ! calculate inverse
    nt2 = Xh%lx*Xh%ly*Xh%lz
    call invcol1(amr_rcn%el_mult, nt2)

    ! to get proper J-1 on faces and edges for fast diagonalisation method
    ! I assume here LX1=LY1=LZ1, so only one array is needed
    !call ftovecl(amr_rcn%fc_mult,amr_rcn%el_mult,1,Xh%lx,Xh%ly,Xh%lz)
    !call etovec(amr_rcn%ed_mult,1,amr_rcn%el_mult,Xh%lx,Xh%ly,Xh%lz)
    ! THIS SHOLD BE CHECKED, BUT SHOULD BE FINE
    amr_rcn%fc_mult(:,:) = amr_rcn%el_mult(:,:, 1)
    amr_rcn%ed_mult(:) = amr_rcn%el_mult(:, 1, 1)

    ! work space
    allocate(amr_rcn%tmp(Xh%lx, Xh%ly, Xh%lz, 3))

    deallocate(tmpl)

    return
  end subroutine amr_rcn_init

  !> Free type arrays
  subroutine amr_rcn_free(ref)
    ! argument list
    type(amr_rcn_t), intent(inout) :: ref

    NULLIFY(ref%msh, ref%Xh)

    if (allocated(ref%x_cr_to_fn)) deallocate(ref%x_cr_to_fn)
    if (allocated(ref%x_cr_to_fnT)) deallocate(ref%x_cr_to_fnT)
    if (allocated(ref%x_fn_to_cr)) deallocate(ref%x_fn_to_cr)
    if (allocated(ref%x_fn_to_crT)) deallocate(ref%x_fn_to_crT)

    if (allocated(ref%y_cr_to_fn)) deallocate(ref%y_cr_to_fn)
    if (allocated(ref%y_cr_to_fnT)) deallocate(ref%y_cr_to_fnT)
    if (allocated(ref%y_fn_to_cr)) deallocate(ref%y_fn_to_cr)
    if (allocated(ref%y_fn_to_crT)) deallocate(ref%y_fn_to_crT)

    if (allocated(ref%z_cr_to_fn)) deallocate(ref%z_cr_to_fn)
    if (allocated(ref%z_cr_to_fnT)) deallocate(ref%z_cr_to_fnT)
    if (allocated(ref%z_fn_to_cr)) deallocate(ref%z_fn_to_cr)
    if (allocated(ref%z_fn_to_crT)) deallocate(ref%z_fn_to_crT)

    if (allocated(ref%el_mult)) deallocate(ref%el_mult)
    if (allocated(ref%fc_mult)) deallocate(ref%fc_mult)
    if (allocated(ref%ed_mult)) deallocate(ref%ed_mult)

    call ref%msh_trs%free()
    call ref%msh_rcn%free()

    if (allocated(ref%tmp)) deallocate(ref%tmp)
    
    return
  end subroutine amr_rcn_free

  !> @brief Map a single coarse element to a fine one
  !! @param[inout] ref     refinement data
  !! @param[in]    ch_pos  child position
  !! @param[in]    vc      coarse element vector
  !! @param[out]   vf      fine element vector
  subroutine amr_rcn_map_ctof(ref, ch_pos, vc, vf)
    ! argument list
    type(amr_rcn_t), intent(inout) :: ref
    integer(i4), dimension(3), intent(in) :: ch_pos
    real(dp), dimension(:,:,:), intent(in) :: vc
    real(dp), dimension(:,:,:), intent(out) :: vf
    ! local variables
    integer(i4) :: iz

    if (ref%msh%gdim == 3) then ! 3D
       call mxm(ref%x_cr_to_fn(:,:, ch_pos(1)), ref%Xh%lx, vc, &
            & ref%Xh%lx, ref%tmp(:,:,:, 1), ref%Xh%lyz)
       do iz = 1, ref%Xh%lz
          call mxm(ref%tmp(:,:, iz, 1), ref%Xh%lx, ref%y_cr_to_fnT(:,:, ch_pos(2)), &
               & ref%Xh%ly, ref%tmp(:,:, iz, 2), ref%Xh%ly)
       end do
       call mxm(ref%tmp(:,:,:, 2), ref%Xh%lxy, ref%z_cr_to_fnT(:,:, ch_pos(3)), &
            & ref%Xh%lz, vf, ref%Xh%lz)
    else ! 2D
       call mxm(ref%x_cr_to_fn(:,:, ch_pos(1)), ref%Xh%lx, vc, ref%Xh%lz, &
            & ref%tmp(:,:,1,1), ref%Xh%lyz)
       call mxm(ref%tmp(:,:,1,1), ref%Xh%lx, ref%y_cr_to_fnT(:,:, ch_pos(2)), &
            & ref%Xh%ly, vf, ref%Xh%ly)
    end if
    
    return
  end subroutine amr_rcn_map_ctof

  !> @brief Map a single fine element to a coarse one
  !! @param[inout] ref     refinement data
  !! @param[in]    ch_pos  child position
  !! @param[in]    vf      fine element vector
  !! @param[out]   vc      coarse element vector
  subroutine amr_rcn_map_ftoc(ref, ch_pos, vf, vc)
    ! argument list
    type(amr_rcn_t), intent(inout) :: ref
    integer(i4), dimension(3), intent(in) :: ch_pos
    real(dp), dimension(:,:,:), intent(in) :: vf
    real(dp), dimension(:,:,:), intent(out) :: vc
    ! local variables
    integer(i4) :: iz

    if (ref%msh%gdim == 3) then ! 3D
       call mxm(ref%x_fn_to_cr(:,:, ch_pos(1)), ref%Xh%lx, vf, &
            & ref%Xh%lx, ref%tmp(:,:,:, 1), ref%Xh%lyz)
       do iz = 1, ref%Xh%lz
          call mxm(ref%tmp(:,:, iz, 1), ref%Xh%lx, ref%y_fn_to_crT(:,:, ch_pos(2)), &
               & ref%Xh%ly, ref%tmp(:,:, iz, 2), ref%Xh%ly)
       end do
       call mxm(ref%tmp(:,:,:, 2), ref%Xh%lxy, ref%z_fn_to_crT(:,:, ch_pos(3)), &
            & ref%Xh%lz, vc, ref%Xh%lz)
    else ! 2D
       call mxm(ref%x_fn_to_cr(:,:, ch_pos(1)), ref%Xh%lx, vf, ref%Xh%lz, &
            & ref%tmp(:,:,1,1), ref%Xh%lyz)
       call mxm(ref%tmp(:,:,1,1), ref%Xh%lx, ref%y_fn_to_crT(:,:, ch_pos(2)), &
            & ref%Xh%ly, vc, ref%Xh%ly)
    end if
    
    return
  end subroutine amr_rcn_map_ftoc

  !> @brief Perform single field refinement operation
  !! @param[inout] ref     refinement data
  !! @param[inout] vcf     refined vector
  subroutine amr_rcn_refine_single(ref, vfc)
    ! argument list
    type(amr_rcn_t), intent(inout) :: ref
    real(dp), dimension(:,:,:,:), intent(inout) :: vfc
    ! local variables
    integer(i4) :: il, jl, itmp
    integer(i4), dimension(3) :: ch_pos

    ! I assume el_lst(1) gives position of the coarse block
    ! and final ch_pos() = 1,1,1
    ! loop over refined elements
    do il= 1, ref%msh_rcn%rfn_nr, ref%msh%npts
       ! copy coarse element to a temporary array
       ref%tmp(:,:,:, 3) = vfc(:,:,:, ref%msh_rcn%elgl_rfn(3, il))
!!$       test1 : block
!!$         integer :: kl,ll
!!$         if (pe_rank == 0) then
!!$         do kl = 1, ref%Xh%lx
!!$            do ll = 1, ref%Xh%lx
!!$               write(*,*) "OLD",ref%msh_rcn%elgl_rfn(3, il),ll,kl,ref%tmp(:,ll,kl, 3)
!!$            end do
!!$            write(*,*) ' '
!!$         end do
!!$         write(*,*) '==============================================',pe_rank, &
!!$              & ref%msh_rcn%elgl_rfn(3, il)
!!$         end if
!!$       end block test1
       ! loop over all the children
       do jl= 1, ref%msh%npts
          ! get child position
          itmp = ref%msh_rcn%elgl_rfn(3, il + jl - 1) ! new position in the array
          ch_pos(3) = (jl-1)/4 + 1 ! z position
          ch_pos(2) = mod((jl - 1)/2, 2) + 1 ! y position 
          ch_pos(1) = mod(jl - 1, 2) +1 ! x position
          ! refine
          call amr_rcn_map_ctof(ref, ch_pos, ref%tmp(:,:,:,3), vfc(:,:,:,itmp))
!!$          test2 : block
!!$            integer :: kl,ll
!!$            if (pe_rank == 0) then
!!$            do kl = 1, ref%Xh%lx
!!$               do ll = 1, ref%Xh%lx
!!$                  write(*,*) "NEW",itmp,ll,kl,vfc(:,ll,kl, itmp)
!!$               end do
!!$               write(*,*) ' '
!!$            end do
!!$            write(*,*) '==============================================',pe_rank, itmp,ch_pos(:)
!!$            endif
!!$          end block test2
       end do
    end do

    return
  end subroutine amr_rcn_refine_single

  !> @brief Perform single field coarsening operation
  !! @param[inout] ref     refinement data
  !! @param[inout] vcf     refined vector
  subroutine amr_rcn_coarsen_single(ref, vfc)
    ! argument list
    type(amr_rcn_t), intent(inout) :: ref
    real(dp), dimension(:,:,:,:), intent(inout) :: vfc
    ! local variables
    integer(i4) :: il, jl, itmp
    integer(i4), dimension(3) :: ch_pos

    ! I assume el_lst(1) gives position of the coarse block
    ! and final ch_pos() = 1,1,1
    ! loop over coarsened elements
    do il= 1, ref%msh_rcn%crs_nr
       ! loop over all the children
       do jl= 1, ref%msh%npts
          ! get child position
          itmp = ref%msh_rcn%elgl_crs(2, jl, il) ! new position in the array
          ch_pos(3) = (jl-1)/4 + 1 ! z position
          ch_pos(2) = mod((jl - 1)/2, 2) + 1 ! y position 
          ch_pos(1) = mod(jl - 1, 2) +1 ! x position
          ! coarsen
          call amr_rcn_map_ftoc(ref, ch_pos, vfc(:,:,:,itmp), ref%tmp(:,:,:,3))
          ! sum contribition
          if (jl == 1) then
             vfc(:,:,:,itmp) = ref%tmp(:,:,:,3)
          else
             vfc(:,:,:,itmp) =  vfc(:,:,:,itmp) + ref%tmp(:,:,:,3)
          end if
       end do
       vfc(:,:,:,itmp) =  vfc(:,:,:,itmp)*ref%el_mult(:,:,:)
    end do

    return
  end subroutine amr_rcn_coarsen_single

  !> @brief Perform refinement, transfer and coarsening of a single field
  !! @param[inout] ref     refinement data
  !! @param[inout] vcf     refined vector
  subroutine amr_rcn_refine_coarsen_single(ref, vfc)
    ! argument list
    type(amr_rcn_t), intent(inout) :: ref
    real(dp), dimension(:,:,:,:), intent(inout) :: vfc
    ! local variable
    integer(i4) :: il
    integer(i4), dimension(4) :: lshape
    real(dp), allocatable, dimension(:,:,:,:) :: svfc

    ! local refinement
    if (ref%msh_rcn%rfn_nr > 0) then
       if (mod(ref%msh_rcn%rfn_nr, ref%msh%npts) /= 0) &
            & call neko_error('Number of ref elem not multiply of nvert')
       call amr_rcn_refine_single(ref, vfc)
    end if

    ! place for transfer/sorting
    ! JUST TESTING VERSION; THIS HAS TO BE REWRITTEN
    lshape(:) = shape(vfc)
    allocate(svfc(lshape(1),lshape(2),lshape(3),lshape(4)))
    do il = 1, ref%msh%nelv
       if (ref%msh_rcn%elgl_map(3, il) /= pe_rank) &
            & call neko_error('No data transfer between ranks for now')
       svfc(:,:,:,ref%msh_rcn%elgl_map(1, il) - ref%msh%offset_el) = vfc(:,:,:,il)
    end do
    vfc(:,:,:,:) = svfc(:,:,:,:)
    deallocate(svfc)
    
!!$    test : block
!!$      integer(i4) :: il
!!$      do il = 1, ref%msh%nelv
!!$         write(*,*) 'TEST mapping', pe_rank, lshape(:), il,&
!!$              & ref%msh_rcn%elgl_map(1, il) - ref%msh%offset_el, ref%msh_rcn%elgl_map(3, il)
!!$      end do
!!$    end block test
    

    ! local coarsenig
    if (ref%msh_rcn%crs_nr > 0) then
       call amr_rcn_coarsen_single(ref, vfc)
    end if

    return
  end subroutine amr_rcn_refine_coarsen_single

  !> Refine all elements in the mesh
  subroutine amr_msh_refine(msh, param, ref_mark)
    ! argument list
    type(mesh_t), intent(inout) :: msh
    type(param_t), intent(in) :: param
    integer(i4), dimension(:), intent(inout) :: ref_mark
    ! local variables
    integer(i4) :: level_max, il, jl, ielo, itmp, ips, ich
    integer(i4) :: ierr, gnum, goff
    logical :: ifmod
    integer(i4), dimension(5,8) :: itmpv
    integer(i4), allocatable, dimension(:) :: el_gidx
    integer(i4), allocatable, dimension(:,:) :: lel_map

    call neko_log%section("Mesh refinement")
    ! sanity check
    if(msh%nelv /= size(ref_mark)) call neko_error('Inconsistent ref_mark size')
    ! check if any refinement is marked
    itmp = minval(ref_mark)
    call MPI_Allreduce(itmp, il, 1, MPI_INTEGER, MPI_MAX, NEKO_COMM, ierr)
    itmp = maxval(ref_mark)
    call MPI_Allreduce(itmp, jl, 1, MPI_INTEGER, MPI_MIN, NEKO_COMM, ierr)
    if (il == 0 .and. jl == 0) then
       call neko_log%message('No refinemnt mark set')
       call neko_log%end_section()
       return
    end if

    ! save old number of local elements
    amr_rcn%msh_rcn%nelvo = msh%nelv
    ! perform refinement/coarsening on p4est side
    level_max = param%amrlmax
    ! global elements numbers; THIS IS JUST A TEMPORARY HACK
    allocate(el_gidx(msh%nelv))
    do il = 1, msh%nelv
       el_gidx(il) = msh%offset_el + il
    end do
    call p4_refine(ref_mark, el_gidx, amr_rcn%msh_trs, level_max, ifmod, amr_rcn%msh_rcn)
    deallocate(el_gidx)

    if (ifmod)  then
       ! Import mesh from p4est to neko
       call p4_msh_get(msh)

!!$       ! PLACE FOR NEW MESH PARTITIONING
!!$
!!$       ! PLACE FOR TRANSFER/SORTING OF MAPPING DATA DATA
!!$       ! THIS SHOULD BE CHANGED AFTER ADDING COMMUNICATION AS IT WOULD INCLUDE COPY
!!$       ! collect information: old gidx, old nid, new gidx, new nid
!!$       ! the max size of tmp arrays could be max(nelv_new,nelv_old)*children_number
!!$       
!!$
!!$       ! Reshuffle data to get propper mapping on the destination rank
!!$       ! THIS SHOULD BE CHANGED AFTER ADDING COMMUNICATION
!!$       allocate(lel_map(3,max(msh%nelv,amr_rcn%msh_rcn%nelvo)*msh%npts)) ! at this point I do not know the proper size, so take a safe value
!!$       lel_map(:,:) = 0
!!$       lel_map(3,:) = -1
!!$       do il = 1, msh%nelv ! loop witn neko distribution
!!$          if (amr_rcn%msh_rcn%elgl_map(1,il) /= 0) then
!!$             lel_map(1, amr_rcn%msh_rcn%elgl_map(2,il)) = msh%offset_el + il ! this is just speciffic to a test case; in general wrong
!!$             lel_map(3, amr_rcn%msh_rcn%elgl_map(2,il)) = pe_rank ! this is just speciffic to a test case; in general wrong
!!$          end if
!!$       end do
!!$       call MOVE_ALLOC(lel_map,amr_rcn%msh_rcn%elgl_map)
!!$       ! PLACE TO UPDATE amr_rcn%msh_rcn%map_nr TO NEKO LOCAL VALUE
!!$
!!$       ! PLACE FOR TRANSFER/SORTING OF REFINEMENT DATA
!!$
!!$       ! Recalculate local position of refined element on the destination rank
!!$       ! THIS SHOULD BE CHANGED AFTER ADDING COMMUNICATION AS IT WOULD INCLUDE COPY
!!$       ! RIGHT NOW IT IS JUST LOCAL POSITION UPDATE
!!$       ! SORTING IS IMPORTANT
!!$       ! loop over refined elements
!!$       itmp = amr_rcn%msh_rcn%nelvo
!!$       do il = 1, amr_rcn%msh_rcn%rfn_nr, amr_rcn%msh%npts ! loop over parent owner
!!$          ips = itmp
!!$          ielo = amr_rcn%msh_rcn%elgl_rfn(3, il)
!!$          ! THIS COPY IS HERE AS THERE IS NO COMMUNICATION
!!$          itmpv(:, 1:amr_rcn%msh%npts) = &
!!$               & amr_rcn%msh_rcn%elgl_rfn(:, il:il + amr_rcn%msh%npts - 1)
!!$          ! loop over all the children
!!$          do jl = 1, amr_rcn%msh%npts
!!$             ! which child
!!$             ich = itmpv(5, jl)
!!$             ! local parent position sanity check
!!$             if (itmpv(3, jl) /= ielo) &
!!$                  & call neko_error('Children do not share parent')
!!$             ! new local position at the end of the array
!!$             if (ich == 0) then
!!$                amr_rcn%msh_rcn%elgl_rfn(3, il + ich) = ielo
!!$             else
!!$                itmp = itmp + 1
!!$                amr_rcn%msh_rcn%elgl_rfn(3, il + ich) = ips + ich
!!$             end if
!!$             ! new global element number
!!$             amr_rcn%msh_rcn%elgl_rfn(1, il + ich) = itmpv(1, jl)
!!$             ! parent local position in the array
!!$             amr_rcn%msh_rcn%elgl_rfn(2, il + ich) = ielo
!!$          end do
!!$       end do
!!$       ! save local number of elements after refinement
!!$       amr_rcn%msh_rcn%rfn_nr_a = itmp
!!$       ! PLACE TO UPDATE amr_rcn%msh_rcn%rfn_nr TO NEKO LOCAL VALUE
!!$
!!$       ! update global rank mapping
!!$       do il = 1, amr_rcn%msh_rcn%rfn_nr ! loop over parent owner
!!$          jl = amr_rcn%msh_rcn%elgl_rfn(3, il)
!!$          if (amr_rcn%msh_rcn%elgl_map(1, jl) == 0) then
!!$             amr_rcn%msh_rcn%elgl_map(1, jl) = amr_rcn%msh_rcn%elgl_rfn(1, il)
!!$             amr_rcn%msh_rcn%elgl_map(3, jl) = pe_rank ! this is just speciffic to a test case; in general wrong
!!$          else
!!$             call neko_error('Refinement transfer index already used; refinement.')
!!$          end if
!!$       end do
!!$
!!$       ! Provide global numbering of elements used for coarsening
!!$       call MPI_Allreduce(amr_rcn%msh_rcn%crs_nr, gnum, 1, MPI_INTEGER, MPI_SUM, & !????????
!!$            & NEKO_COMM, ierr)
!!$       ! get global offset
!!$       call MPI_Scan(amr_rcn%msh_rcn%crs_nr, goff, 1, MPI_INTEGER, MPI_SUM, &
!!$            & NEKO_COMM, ierr)
!!$       goff = goff - amr_rcn%msh_rcn%crs_nr
!!$       if (amr_rcn%msh_rcn%crs_nr > 0) then
!!$          itmp = msh%glb_nelv + (msh%npts - 1)*goff
!!$          do il = 1, amr_rcn%msh_rcn%crs_nr ! loop on p4est distribution
!!$             do jl = 2, msh%npts
!!$                itmp = itmp + 1
!!$                amr_rcn%msh_rcn%elgl_crs(1, jl, il) = itmp
!!$             end do
!!$          end do
!!$       end if
!!$
!!$       ! PLACE FOR COARSENING DATA TRANSFER AND SORTING; child owner
!!$
!!$       ! save local number of elements to be coarsened (child owner)
!!$       ! crs_nr_s = ???? ! THIS IS JUST TEST SPECIFFIC in general it should be result of communication
!!$       if (amr_rcn%msh_rcn%crs_nr == 0) then
!!$          amr_rcn%msh_rcn%crs_nr_s = 0
!!$       else
!!$          amr_rcn%msh_rcn%crs_nr_s = msh%nelv
!!$       end if
!!$       ! sanity check (sum of unchaged, refined parents and corsend children has to be equal nelvo)
!!$       if (amr_rcn%msh_rcn%nelvo /= amr_rcn%msh_rcn%crs_nr_s + &
!!$            & amr_rcn%msh_rcn%map_nr + &
!!$            & amr_rcn%msh_rcn%rfn_nr/msh%npts) call neko_error('Inconsiten')
!!$
!!$       ! update global rank mapping
!!$       ! THIS LOOP SHOULD LOOK DIFFERENT WITH COMMUNICATION
!!$       do il = 1, amr_rcn%msh_rcn%crs_nr ! loop over elements on child owner
!!$          do jl = 1, msh%npts
!!$             if (amr_rcn%msh_rcn%elgl_map(1, amr_rcn%msh_rcn%elgl_crs(3, jl, il)) == 0) then
!!$                amr_rcn%msh_rcn%elgl_map(1, jl) = amr_rcn%msh_rcn%elgl_crs(1, jl, il)
!!$                amr_rcn%msh_rcn%elgl_map(3, jl) = pe_rank ! this is just speciffic to a test case; in general wrong
!!$             else
!!$                call neko_error('Refinement transfer index already used; coarsening.')
!!$             end if
!!$          end do
!!$       end do
!!$
!!$       ! PLACE FOR COARSENING DATA TRANSFER AND SORTING; coarsened element owner
!!$
!!$       ! recalculate new childen position
!!$       ! THIS LOOP HAS TO BE CHANGED FOR COMMUNICATION
!!$       itmp = msh%nelv
!!$       do il = 1, amr_rcn%msh_rcn%crs_nr ! loop over children on parent owner
!!$          ips = itmp
!!$          !ielo = PARENT GLOBAL NUMBER
!!$          ! loop over all the children
!!$          do jl = 1, amr_rcn%msh%npts
!!$             ! sanity check; global parent position
!!$             ! NOT NEEDED FOR THIS EXAMPLE
!!$             ! which child
!!$             ! EQUAL TO JL IN THIS EXAMPLE
!!$             ich = jl
!!$             ! new global element number
!!$             ! DOESN'T CHANGE IN THIS EXAMPLE
!!$             ! local child position
!!$             if (ich == 1) then
!!$                amr_rcn%msh_rcn%elgl_crs(1, jl, il) = il ! this is not correct in general
!!$             else
!!$                itmp = itmp + 1
!!$                amr_rcn%msh_rcn%elgl_crs(1, jl, il) = ips + ich - 1
!!$             end if
!!$          end do
!!$       end do
!!$       ! save local number of elements before coarsening
!!$       amr_rcn%msh_rcn%crs_nr_b = itmp
!!$       ! PLACE TO UPDATE amr_rcn%msh_rcn%crs_nr TO NEKO LOCAL VALUE (parent owner)
          
    end if ! ifmod

!!$    testing : block
!!$      integer :: ierr
!!$      write(*,*) 'TESTING0', pe_rank, amr_rcn%msh_rcn%nelvo, amr_rcn%msh%nelv, &
!!$           & amr_rcn%msh_rcn%map_nr, amr_rcn%msh_rcn%rfn_nr, amr_rcn%msh_rcn%crs_nr,&
!!$           & amr_rcn%msh_rcn%rfn_nr_a, amr_rcn%msh_rcn%crs_nr_s, amr_rcn%msh_rcn%crs_nr_b, goff, gnum
!!$      call MPI_Barrier(NEKO_COMM, ierr)
!!$      !call neko_log%end_section()
!!$      return
!!$      !call neko_error('This is not an error.')
!!$    end block testing

    ! free memory

    call neko_log%end_section()

    return
  end subroutine amr_msh_refine

  ! Refine all elements in field
  subroutine amr_fluid_refine(msh, fld, param, ref_mark)
    ! argument list
    type(mesh_t), intent(inout) :: msh
    class(fluid_scheme_t), intent(inout) :: fld
    type(param_t), intent(in) :: param
    integer(i4), dimension(:), intent(inout) :: ref_mark
    ! local variables
    integer(i4) :: itmp
    real(dp), allocatable, dimension(:,:,:,:) :: tmpv

    ! ALL THIS SHOULD BE PART OF TYPE EXTENSION, BUT FOR TESTING IT SHOULD BE FINE
    select type(fld)
    type is(fluid_plan1_t)
       call neko_error('Nothing done for plan1')
    type is(fluid_plan4_t)
       call neko_error('Nothing done for plan4')
    type is(fluid_pnpn_t)
       ! perform refinement on p4est side and import the mesh
       call amr_msh_refine(msh, param, ref_mark)

       !test refinement on a single field
       allocate(tmpv(fld%Xh%lx, fld%Xh%ly, fld%Xh%lz,msh%nelv))
       tmpv(:,:,:,1:amr_rcn%msh_rcn%nelvo) = fld%dm_Xh%x(:,:,:,:)
       call amr_rcn_refine_coarsen_single(amr_rcn, tmpv)
       ! reinitialise variables
       ! degrees of freedom; THIS SHOULD BE POSSIBLY CHANGED TO SUBROUTINE TAKING VARIABLE
       !call fld%dm_Xh%resize()

  

       deallocate(tmpv)

       testing : block
         integer :: ierr
         write(*,*) 'TEST size', pe_rank, amr_rcn%msh_rcn%nelvo, msh%nelv, msh%glb_nelv
         call MPI_Barrier(NEKO_COMM, ierr)
         call neko_error('This is not an error.')
       end block testing
    end select

    return
  end subroutine amr_fluid_refine

end module amr
