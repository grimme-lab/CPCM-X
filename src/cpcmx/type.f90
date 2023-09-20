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

module type
    use mctc_env, only: wp
    implicit none

    private

    public :: parameter_type, molecule_data
    !> CRS parameter type
    type parameter_type
        real(wp) :: rav
        real(wp) :: aprime
        real(wp) :: fcorr
        real(wp) :: chb
        real(wp) :: shb
        real(wp) :: aeff
        real(wp) :: lambda
        real(wp) :: omega
        real(wp) :: eta
        real(wp) :: shift

    contains
        
        procedure :: copy_param_type
        generic :: assignment(=) => copy_param_type

    end type parameter_type


    !> Molecule data type
    type molecule_data
        !> COSMO sigma shieldings
        real(wp),dimension(:), allocatable :: su
        !> COSMO sigma shieldings (averaged)
        real(wp),dimension(:), allocatable :: sv, sv0, svt
        !> Area of the COSMO surface segments
        real(wp),dimension(:), allocatable :: area
        !> COSMO surface segment coordinates
        real(wp),dimension(:,:), allocatable :: xyz
        !> Potential of COSMO surface segments
        real(wp),dimension(:), allocatable :: pot
        !> Interacting potential of COSMO surface segments
        real(wp),dimension(:), allocatable :: poti
        !> Identification of COSMO surface segments
        integer, dimension(:), allocatable :: id
        !> COSMO sigma profile
        real(wp),dimension(:) :: sigma(0:50), sigma3(3,0:50)
        !> Element type of the molecule atoms
        character(len=2), dimension(:), allocatable :: element
        !> Coordinates of the molecule atoms
        real(wp), dimension(:,:), allocatable :: atom_xyz
        !> Bonding situation
        logical, dimension(:,:), allocatable :: is_bonded
        !> Volume of the molecule
        real(wp) :: volume
        !> COSMO energy of the molecule
        real(wp) :: energy
        !> Gas phase energy of the molecule
        real(wp) :: energy_gas
        !> Number of OH Bonds
        integer :: noh
        !> Number of NH Bonds
        integer :: nnh
        !> Hydrogen bonding groups
        character(len=2), allocatable, dimension(:) :: hb_group
        !> Belongs the atom to an ring?
        logical, dimension(:), allocatable :: ring
        !> Number of effective ring atoms
        integer :: near
    end type molecule_data

    contains

    subroutine copy_param_type(lhs,rhs)
        class(parameter_type),intent(out) :: lhs
        class(parameter_type), intent(in) :: rhs

        lhs%rav=rhs%rav
        lhs%aprime=rhs%aprime
        lhs%fcorr=rhs%fcorr
        lhs%chb=rhs%chb
        lhs%shb=rhs%shb
        lhs%aeff=rhs%aeff
        lhs%lambda=rhs%lambda
        lhs%omega=rhs%omega
        lhs%eta=rhs%eta
        lhs%shift=rhs%shift
    end subroutine copy_param_type


end module type