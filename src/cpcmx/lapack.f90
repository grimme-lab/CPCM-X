! This file is part of CPCM-X.
! SPDX-Identifier: LGPL-3.0-or-later
!
! CPCM-X is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! CPCM-X is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with CPCM-X.  If not, see <https://www.gnu.org/licenses/>.

!> Convenience wrapper routines for common LAPACK functions
module crs_lapack
   use mctc_env, only : sp, dp
   implicit none
   private

   public :: wrap_getrf, wrap_getrs


   !> Computes an LU factorization of a general M-by-N matrix A
   !> using partial pivoting with row interchanges.
   !>
   !> The factorization has the form
   !>    A = P * L * U
   !> where P is a permutation matrix, L is lower triangular with unit
   !> diagonal elements (lower trapezoidal if m > n), and U is upper
   !> triangular (upper trapezoidal if m < n).
   interface wrap_getrf
      module procedure :: wrap_sgetrf
      module procedure :: wrap_dgetrf
   end interface wrap_getrf


   !> Solves a system of linear equations
   !>    A * X = B  or  A**T * X = B
   !> with a general N-by-N matrix A using the LU factorization computed
   !> by ?GETRF.
   interface wrap_getrs
      module procedure :: wrap_sgetrs
      module procedure :: wrap_dgetrs
   end interface wrap_getrs


   !> Computes an LU factorization of a general M-by-N matrix A
   !> using partial pivoting with row interchanges.
   !>
   !> The factorization has the form
   !>    A = P * L * U
   !> where P is a permutation matrix, L is lower triangular with unit
   !> diagonal elements (lower trapezoidal if m > n), and U is upper
   !> triangular (upper trapezoidal if m < n).
   interface lapack_getrf
      pure subroutine sgetrf(m, n, a, lda, ipiv, info)
         import :: sp
         real(sp), intent(inout) :: a(lda, *)
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine sgetrf
      pure subroutine dgetrf(m, n, a, lda, ipiv, info)
         import :: dp
         real(dp), intent(inout) :: a(lda, *)
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine dgetrf
   end interface lapack_getrf


   !> Solves a system of linear equations
   !>    A * X = B  or  A**T * X = B
   !> with a general N-by-N matrix A using the LU factorization computed
   !> by ?GETRF.
   interface lapack_getrs
      pure subroutine sgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: sp
         real(sp), intent(in) :: a(lda, *)
         integer, intent(in) :: ipiv(*)
         real(sp), intent(inout) :: b(ldb, *)
         character(len=1), intent(in) :: trans
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine sgetrs
      pure subroutine dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: dp
         real(dp), intent(in) :: a(lda, *)
         integer, intent(in) :: ipiv(*)
         real(dp), intent(inout) :: b(ldb, *)
         character(len=1), intent(in) :: trans
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine dgetrs
   end interface lapack_getrs


contains


subroutine wrap_sgetrf(amat, ipiv, info)
   real(sp), intent(inout) :: amat(:, :)
   integer, intent(out) :: ipiv(:)
   integer, intent(out) :: info
   integer :: m, n, lda
   lda = max(1, size(amat, 1))
   m = size(amat, 1)
   n = size(amat, 2)
   call lapack_getrf(m, n, amat, lda, ipiv, info)
end subroutine wrap_sgetrf


subroutine wrap_dgetrf(amat, ipiv, info)
   real(dp), intent(inout) :: amat(:, :)
   integer, intent(out) :: ipiv(:)
   integer, intent(out) :: info
   integer :: m, n, lda
   lda = max(1, size(amat, 1))
   m = size(amat, 1)
   n = size(amat, 2)
   call lapack_getrf(m, n, amat, lda, ipiv, info)
end subroutine wrap_dgetrf


subroutine wrap_sgetrs(amat, bmat, ipiv, info, trans)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(inout) :: bmat(:, :)
   integer, intent(in) :: ipiv(:)
   integer, intent(out) :: info
   character(len=1), intent(in), optional :: trans
   character(len=1) :: tra
   integer :: n, nrhs, lda, ldb
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   n = size(amat, 2)
   nrhs = size(bmat, 2)
   call lapack_getrs(tra, n, nrhs, amat, lda, ipiv, bmat, ldb, info)
end subroutine wrap_sgetrs


subroutine wrap_dgetrs(amat, bmat, ipiv, info, trans)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(inout) :: bmat(:, :)
   integer, intent(in) :: ipiv(:)
   integer, intent(out) :: info
   character(len=1), intent(in), optional :: trans
   character(len=1) :: tra
   integer :: n, nrhs, lda, ldb
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   n = size(amat, 2)
   nrhs = size(bmat, 2)
   call lapack_getrs(tra, n, nrhs, amat, lda, ipiv, bmat, ldb, info)
end subroutine wrap_dgetrs


end module crs_lapack
