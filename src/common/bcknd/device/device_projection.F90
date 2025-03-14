! Copyright (c) 2020-2021, The Neko Authors
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
!> Interface for device projection
!! @note Requires device MPI
module device_projection
  use num_types
  use utils
  use, intrinsic :: iso_c_binding
  implicit none

#ifdef HAVE_HIP
  interface
     subroutine hip_project_on(a_d, b_d, x_d_d, b_d_d, mult_d, x_d, j, n) &
          bind(c, name='hip_project_on')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d, x_d_d, b_d_d, mult_d, x_d
       integer(c_int) :: j, n
     end subroutine hip_project_on
  end interface

  interface
     subroutine hip_project_ortho(a_d, b_d, x_d_d, b_d_d, &
                                   w_d, xm_d, j, n, nrm) &
          bind(c, name='hip_project_ortho')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d, x_d_d, b_d_d, w_d
       type(c_ptr), value :: xm_d
       integer(c_int) :: j, n
       real(c_rp) :: nrm
     end subroutine hip_project_ortho
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_project_on(a_d, b_d, x_d_d, b_d_d, mult_d, x_d, j, n) &
          bind(c, name='cuda_project_on')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d, x_d_d, b_d_d, mult_d, x_d
       integer(c_int) :: j, n
     end subroutine cuda_project_on
  end interface

  interface
     subroutine cuda_project_ortho(a_d, b_d, x_d_d, b_d_d, &
                                   w_d, xm_d, j, n, nrm) &
          bind(c, name='cuda_project_ortho')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d, x_d_d, b_d_d, w_d
       type(c_ptr), value :: xm_d
       integer(c_int) :: j, n
       real(c_rp) :: nrm
     end subroutine cuda_project_ortho
  end interface
#endif

contains

  subroutine device_proj_on(alpha_d, b_d, x_d_d, b_d_d, mult_d, xbar_d, j, n)
    type(c_ptr), value :: alpha_d, b_d, x_d_d, b_d_d, mult_d, xbar_d
    integer(c_int) :: j, n
    integer :: ierr
#ifdef HAVE_HIP
    call hip_project_on(alpha_d, b_d, x_d_d, b_d_d, mult_d, xbar_d, j, n)
#elif HAVE_CUDA
    call cuda_project_on(alpha_d, b_d, x_d_d, b_d_d, mult_d, xbar_d, j, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_proj_on

  subroutine device_project_ortho(alpha_d, b_d, x_d_d, b_d_d, &
                                  w_d, xm_d, j, n, nrm)
    type(c_ptr), value :: alpha_d, b_d, x_d_d, b_d_d
    type(c_ptr), value :: w_d,  xm_d
    integer(c_int) :: j, n
    real(c_rp) :: nrm
    integer :: ierr
#ifdef HAVE_HIP
    call hip_project_ortho(alpha_d, b_d, x_d_d, b_d_d, w_d, xm_d, j, n, nrm)
#elif HAVE_CUDA
    call cuda_project_ortho(alpha_d, b_d, x_d_d, b_d_d, w_d, xm_d, j, n, nrm)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_project_ortho

end module device_projection
