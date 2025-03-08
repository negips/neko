! Copyright (c) 2020-2023, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> Implements type data_streamer_t.
module data_streamer
  use num_types, only: rp, c_rp
  use field, only: field_t
  use coefs, only: coef_t
  use utils, only: neko_warning
  use device
  use comm
  use neko_mpi_types
  use neko_config
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> Provides access to data streaming by interfacing with c++
  !! ADIOS2 subroutines.
  !! @details
  !! Adios2 is an API that allows for easy coupling of codes
  !! through data streaming and gives the posibility to perform
  !! other IO operations such as data compression, etc.
  !! This type wraps and interfaces the needed calls to allow
  !! the use of the c++ routines that ultimately expose the data
  !! from neko to any executable that counts with a proper reader.
  type, public :: data_streamer_t
     !> Define if the execution is asyncrhonous
     integer :: if_asynch
     !> global element numbers
     integer, allocatable :: lglel(:)
   contains
     !> Constructor
     procedure, pass(this) :: init => data_streamer_init
     !> Destructor
     procedure, pass(this) :: free => data_streamer_free
     !> Stream data
     procedure, pass(this) :: stream => data_streamer_stream

  end type data_streamer_t

contains

  !> Constructor
  !! Wraps the adios2 set-up.
  !! @param coef Type that contains geometrical information
  !! on the case.
  !! @param if_asynch Controls whether the asyncrhonous executions
  !! is to be enabled.
  subroutine data_streamer_init(this, coef, if_asynch)
    class(data_streamer_t), intent(inout) :: this
    type(coef_t), intent(inout) :: coef
    integer, intent(in) :: if_asynch
    integer :: nelb, nelb2, nelv, nelgv,npts,e

    !Allocate and initialize the global element number
    allocate(this%lglel(coef%msh%nelv))
    do e = 1, coef%msh%nelv
       this%lglel(e) = e + coef%msh%offset_el
    end do

    !Assign if the streaming is asynchronous
    this%if_asynch = if_asynch

    !Assign the set up parameters
    nelv  = coef%msh%nelv
    npts  = coef%Xh%lx*coef%Xh%ly*coef%Xh%lz
    nelgv = coef%msh%glb_nelv
    nelb  = coef%msh%offset_el
    ! Alternative way to get nelb:
    !nelb = elem_running_sum(nelv)
    !nelb = nelb - nelv

#ifdef HAVE_ADIOS2
    call fortran_adios2_setup(npts, nelv, nelb, nelgv, &
                nelgv, coef%dof%x, coef%dof%y, &
                coef%dof%z, if_asynch, NEKO_COMM)
#else
    call neko_warning('Is not being built with ADIOS2 support.')
    call neko_warning('Not able to use stream/compression functionality')
#endif


  end subroutine data_streamer_init

  !> Destructor
  !! wraps the adios2 finalize routine. Closes insitu writer
  subroutine data_streamer_free(this)
    class(data_streamer_t), intent(inout) :: this

    if (allocated(this%lglel))        deallocate(this%lglel)

#ifdef HAVE_ADIOS2
    call fortran_adios2_finalize()
#else
    call neko_warning('Is not being built with ADIOS2 support.')
    call neko_warning('Not able to use stream/compression functionality')
#endif

  end subroutine data_streamer_free

  !> streamer
  !! wraps the adios2 stream function.
  !! @param u velocity in x
  !! @param v velocity in y
  !! @param w velocity in z
  !! @param p pressure
  !! @param coef type
  subroutine data_streamer_stream(this, u, v, w, p, coef)
    class(data_streamer_t), intent(inout) :: this
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    integer :: nelv, npts

    nelv  = coef%msh%nelv
    npts  = coef%Xh%lx*coef%Xh%ly*coef%Xh%lz

    if (NEKO_BCKND_DEVICE .eq. 1) then
       ! Move the data to the CPU to be able to write it
       call device_memcpy(u%x, u%x_d, nelv*npts, DEVICE_TO_HOST, sync=.false.)
       call device_memcpy(v%x, v%x_d, nelv*npts, DEVICE_TO_HOST, sync=.false.)
       call device_memcpy(w%x, w%x_d, nelv*npts, DEVICE_TO_HOST, sync=.false.)
       call device_memcpy(p%x, p%x_d, nelv*npts, DEVICE_TO_HOST, sync=.true.)
    end if

#ifdef HAVE_ADIOS2
    call fortran_adios2_stream(this%lglel, p%x, u%x, v%x, w%x, coef%B, u%x)
#else
    call neko_warning('Is not being built with ADIOS2 support.')
    call neko_warning('Not able to use stream/compression functionality')
#endif

  end subroutine data_streamer_stream

  !> Supporting function to calculate the element number offset.
  !! returns the number of elements that the ranks previous to the
  !! present one have.
  !! @param nelv number of elements in current rank.
  function elem_running_sum(nelv) result(rbuff)
    integer, intent(in) :: nelv
    integer ::  ierr, xbuff, wbuff, rbuff

    xbuff = nelv  ! running sum
    wbuff = nelv  ! working buff
    rbuff = 0   ! recv buff

    call mpi_scan(xbuff, rbuff, 1, mpi_integer, mpi_sum, NEKO_COMM, ierr)
  end function elem_running_sum

#ifdef HAVE_ADIOS2

  !> Interface to adios2_setup in c++.
  !! @details This routine interfaces with c++ routine that set up adios2
  !! if streaming, the global array to pair writer and reader is opened.
  !! @param npts number of points per element
  !! @param nelv number of elements in this rank
  !! @param nelb number of elements in ranks before this one
  !! @param nelgv total number of elements in velocity mesh
  !! @param nelgt total number of elements in temperature mesh (not used)
  !! @param x coordinates in x direction
  !! @param y coordinates in y direction
  !! @param z coordinates in z direction
  !! @param asynch integer that indicates asynchronous execution
  !! @param comm simulation communicator
  subroutine fortran_adios2_setup(npts, nelv, nelb, nelgv, nelgt, x, y, &
                                  z, asynch, comm)
    use, intrinsic :: ISO_C_BINDING
    implicit none
    real(kind=rp), dimension(:,:,:,:), intent(inout) :: x
    real(kind=rp), dimension(:,:,:,:), intent(inout) :: y
    real(kind=rp), dimension(:,:,:,:), intent(inout) :: z
    integer, intent(in) :: npts, nelv, nelb,nelgv, nelgt, asynch
    type(MPI_COMM) :: comm

    interface
       !> C-definition is: void adios2_setup_(const int *nval,
       !! const int *nelvin,const int *nelb, const int *nelgv,
       !! const int *nelgt, const double *xml,const double *yml,
       !! const double *zml, const int *if_asynchronous,
       !! const int *comm_int)
       subroutine c_adios2_setup(npts, nelv, nelb, nelgv, nelgt, x, y, &
                             z, asynch, comm) bind(C,name="adios2_setup_")
         use, intrinsic :: ISO_C_BINDING
         import c_rp
         implicit none
         integer(kind=C_INT) :: npts
         integer(kind=C_INT) :: nelv
         integer(kind=C_INT) :: nelb
         integer(kind=C_INT) :: nelgv
         integer(kind=C_INT) :: nelgt
         real(kind=c_rp), intent(INOUT) :: x(*)
         real(kind=c_rp), intent(INOUT)  :: y(*)
         real(kind=c_rp), intent(INOUT)  :: z(*)
         integer(kind=C_INT) :: asynch
         type(*) :: comm
       end subroutine c_adios2_setup
    end interface

    call c_adios2_setup(npts, nelv, nelb, nelgv, nelgt, x, y, z, &
                        asynch, comm)
  end subroutine fortran_adios2_setup

  !> Interface to adios2_finalize in c++.
  !! closes any writer openned at initialization time
  subroutine fortran_adios2_finalize()
    use, intrinsic :: ISO_C_BINDING
    implicit none

    interface
       !> C-definition is: void adios2_finalize_()
       subroutine c_adios2_finalize() bind(C,name="adios2_finalize_")
         use, intrinsic :: ISO_C_BINDING
         implicit none
       end subroutine c_adios2_finalize
    end interface

    call c_adios2_finalize()
  end subroutine fortran_adios2_finalize

  !> Interface to adios2_stream in c++.
  !! @details This routine communicates the data to a global array that
  !! is accessed by a data processor. The operations do not write to disk.
  !! data is communicated with mpi.
  !! @param lglel global element number
  !! @param p pressure
  !! @param u velocity in x
  !! @param v velocity in y
  !! @param w velocity in z
  !! @param bm1 mass matrix
  !! @param t temperature / (Not really used in adios2 routine)
  subroutine fortran_adios2_stream(lglel, p, u, v, w, bm1, t)
    use, intrinsic :: ISO_C_BINDING
    implicit none
    integer, dimension(:), intent(inout) :: lglel
    real(kind=rp), dimension(:,:,:,:), intent(inout) :: p
    real(kind=rp), dimension(:,:,:,:), intent(inout) :: u
    real(kind=rp), dimension(:,:,:,:), intent(inout) :: v
    real(kind=rp), dimension(:,:,:,:), intent(inout) :: w
    real(kind=rp), dimension(:,:,:,:), intent(inout) :: bm1
    real(kind=rp), dimension(:,:,:,:), intent(inout) :: t

    interface
       !> C-definition is: void adios2_stream_(
       !! const int *lglel, const double *pr, const double *u,
       !! const double *v, const double *w, const double *mass1,
       !! const double *temp)
       subroutine c_adios2_stream(lglel, p, u, v, w, bm1, t) &
                                  bind(C,name="adios2_stream_")
         use, intrinsic :: ISO_C_BINDING
         import c_rp
         implicit none
         integer(kind=C_INT), intent(INOUT) :: lglel(*)
         real(kind=c_rp), intent(INOUT) :: p(*)
         real(kind=c_rp), intent(INOUT)  :: u(*)
         real(kind=c_rp), intent(INOUT)  :: v(*)
         real(kind=c_rp), intent(INOUT)  :: w(*)
         real(kind=c_rp), intent(INOUT)  :: bm1(*)
         real(kind=c_rp), intent(INOUT)  :: t(*)
       end subroutine c_adios2_stream
    end interface

    call c_adios2_stream(lglel, p, u, v, w, bm1, t)
  end subroutine fortran_adios2_stream
#endif

end module data_streamer
