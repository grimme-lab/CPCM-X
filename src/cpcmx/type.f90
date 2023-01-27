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


end module type