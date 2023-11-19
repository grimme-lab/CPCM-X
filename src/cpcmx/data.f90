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

module data
    use mctc_env, only: wp
    use mctc_io, only: to_number
    use ieee_arithmetic, only: ieee_value, ieee_positive_inf
    use mctc_env, only: error_type, fatal_error
    private
    public :: AtomicMass,solvent_name,minnesota_eps,density

  !> Atomic masses in a.u.
    real(wp), parameter :: atomicMassNist(1:118) = [ &
    &   1.00794075_wp,   4.00260193_wp,   6.94003660_wp,   9.01218307_wp,&
    &  10.81102805_wp,  12.01073590_wp,  14.00670321_wp,  15.99940492_wp,&
    &  18.99840316_wp,  20.18004638_wp,  22.98976928_wp,  24.30505162_wp,&
    &  26.98153853_wp,  28.08549871_wp,  30.97376200_wp,  32.06478741_wp,&
    &  35.45293758_wp,  39.94779856_wp,  39.09830091_wp,  40.07802251_wp,&
    &  44.95590828_wp,  47.86674496_wp,  50.94146504_wp,  51.99613176_wp,&
    &  54.93804391_wp,  55.84514443_wp,  58.93319429_wp,  58.69334711_wp,&
    &  63.54603995_wp,  65.37778253_wp,  69.72306607_wp,  72.62755016_wp,&
    &  74.92159457_wp,  78.95938856_wp,  79.90352778_wp,  83.79800000_wp,&
    &  85.46766360_wp,  87.61664447_wp,  88.90584030_wp,  91.22364160_wp,&
    &  92.90637300_wp,  95.95978854_wp,  97.90721240_wp, 101.06494014_wp,&
    & 102.90549800_wp, 106.41532751_wp, 107.86814963_wp, 112.41155782_wp,&
    & 114.81808663_wp, 118.71011259_wp, 121.75978367_wp, 127.60312648_wp,&
    & 126.90447190_wp, 131.29276145_wp, 132.90545196_wp, 137.32689163_wp,&
    & 138.90546887_wp, 140.11573074_wp, 140.90765760_wp, 144.24159603_wp,&
    & 144.91275590_wp, 150.36635571_wp, 151.96437813_wp, 157.25213065_wp,&
    & 158.92535470_wp, 162.49947282_wp, 164.93032880_wp, 167.25908265_wp,&
    & 168.93421790_wp, 173.05415017_wp, 174.96681496_wp, 178.48497872_wp,&
    & 180.94787564_wp, 183.84177755_wp, 186.20670455_wp, 190.22485963_wp,&
    & 192.21605165_wp, 195.08445686_wp, 196.96656879_wp, 200.59916703_wp,&
    & 204.38341284_wp, 207.21690806_wp, 208.98039910_wp, 208.98243080_wp,&
    & 209.98714790_wp, 222.01757820_wp, 223.01973600_wp, 226.02541030_wp,&
    & 227.02775230_wp, 232.03805580_wp, 231.03588420_wp, 238.02891046_wp,&
    & 237.04817360_wp, 244.06420530_wp, 243.06138130_wp, 247.07035410_wp,&
    & 247.07030730_wp, 251.07958860_wp, 252.08298000_wp, 257.09510610_wp,&
    & 258.09843150_wp, 259.10103000_wp, 262.10961000_wp, 267.12179000_wp,&
    & 269.12791000_wp, 271.13393000_wp, 270.13336000_wp, 276.14846000_wp,&
    & 276.15159000_wp, 280.16131000_wp, 282.16912000_wp, 284.17416000_wp,&
    & 284.17873000_wp, 289.19042000_wp, 288.19274000_wp, 293.20449000_wp,&
    & 292.20746000_wp, 294.21392000_wp]

  !> Get atomic mass for a species
interface AtomicMass
    module procedure :: getAtomicMassSymbol
    module procedure :: getAtomicMassNumber
    module procedure :: getAtomicMassSymbolArray
    module procedure :: getAtomicMassNumberArray
 end interface AtomicMass    

contains

    !> Normalizes Input Solvent Name to Minnesota Solvent Name
    function solvent_name(solvent)
        use globals, only: to_lower
        !> Name of the solvent
        character(len=*), intent(in) :: solvent
        character(len=:), allocatable :: solvent_name

        solvent_name=to_lower(solvent)
        
        select case(solvent_name)
        case default !> Gives back empty solvent for error handling.
            solvent_name='' 
        case('2methylpyridine')
            solvent_name='2methylpyridine'
        case('4methyl2pentanone')
            solvent_name='4methyl2pentanone'
        case('aceticacid')
            solvent_name='aceticacid'
        case('acetonitrile')
            solvent_name='acetonitrile'
        case('acetophenone')
            solvent_name='acetophenone'
        case('aniline')
            solvent_name='aniline'
        case('anisole')
            solvent_name='anisole'
        case('benzene')
            solvent_name='benzene'
        case('benzonitrile')
            solvent_name='benzonitrile'
        case('benzylalcohol')
            solvent_name='benzylalcohol'
        case('bromobenzene')
            solvent_name='bromobenzene'
        case('bromoethane')
            solvent_name='bromoethane'
        case('bromoform')
            solvent_name='bromoform'
        case('bromooctane')
            solvent_name='bromooctane'
        case('butanol')
            solvent_name='butanol'
        case('butanone')
            solvent_name='butanone'
        case('butylacetate')
            solvent_name='butylacetate'
        case('butylbenzene')
            solvent_name='butylbenzene'
        case('carbondisulfide')
            solvent_name='carbondisulfide'
        case('carbontet','carbontetrachlorid','ccl4')
            solvent_name='carbontet' 
        case('chlorobenzene')
            solvent_name='chlorobenzene'
        case('chloroform','chcl3')
            solvent_name='chloroform'
        case('chlorohexane')
            solvent_name='chlorohexane'
        case('cyclohexane')
            solvent_name='cyclohexane'
        case('cyclohexanone')
            solvent_name='cyclohexanone'
        case('decalin')
            solvent_name='decalin'
        case('decane')
            solvent_name='decane'
        case('decanol')
            solvent_name='decanol'
        case('dibromoethane')
            solvent_name='dibromoethane'
        case('dibutylether')
            solvent_name='dibutylether'
        case('dichloroethane')
            solvent_name='dichloroethane'
        case('diethylether')
            solvent_name='diethylether'
        case('diisopropylether')
            solvent_name='diisopropylether'
        case('dimethylacetamide')
            solvent_name='dimethylacetamide'
        case('dimethylformamide')
            solvent_name='dimethylformamide'
        case('dimethylpyridine')
            solvent_name='dimethylpyridine'
        case('dimethylsulfoxide')
            solvent_name='dimethylsulfoxide'
        case('dodecane')
            solvent_name='dodecane'
        case('ethanol')
            solvent_name='ethanol'
        case('ethoxybenzene')
            solvent_name='ethoxybenzene'
        case('ethylacetate')
            solvent_name='ethylacetate'
        case('ethylbenzene')
            solvent_name='ethylbenzene'
        case('fluorobenzene')
            solvent_name='fluorobenzene'
        case('fluoroctane')
            solvent_name='fluoroctane'
        case('heptane')
            solvent_name='heptane'
        case('heptanol')
            solvent_name='heptanol'
        case('hexadecane')
            solvent_name='hexadecane'
        case('hexadecyliodide', 'iodohexadecane', '1-iodohexadecane')
            solvent_name='hexadecyliodide' 
        case('hexane')
            solvent_name='hexane'
        case('hexanol')
            solvent_name='hexanol'
        case('iodobenzene')
            solvent_name='iodobenzene'
        case('isobutanol')
            solvent_name='isobutanol'
        case('isooctane')
            solvent_name='isooctane'
        case('isopropanol')
            solvent_name='isopropanol'
        case('isopropylbenzene')
            solvent_name='isopropylbenzene'
        case('isopropyltoluene')
            solvent_name='isopropyltoluene'
        case('mcresol')
            solvent_name='mcresol'
        case('mesitylene')
            solvent_name='mesitylene'
        case('methoxyethanol')
            solvent_name='methoxyethanol'
        case('methylenechloride')
            solvent_name='methylenechloride'
        case('methylformamide')
            solvent_name='methylformamide'
        case('nitrobenzene')
            solvent_name='nitrobenzene'
        case('nitroethane')
            solvent_name='nitroethane'
        case('nitromethane')
            solvent_name='nitromethane'
        case('nonane')
            solvent_name='nonane'
        case('nonanol')
            solvent_name='nonanol'
        case('octane')
            solvent_name='octane'
        case('octanol')
            solvent_name='octanol'
        case('odichlorobenzene')
            solvent_name='odichlorobenzene'
        case('onitrotoluene')
            solvent_name='onitrotoluene'
        case('pentadecane')
            solvent_name='pentadecane'
        case('pentane')
            solvent_name='pentane'
        case('pentanol')
            solvent_name='pentanol'
        case('perfluorobenzene', 'hexaflurobenzene')
            solvent_name='perfluorobenzene' 
        case('phenylether')
            solvent_name='phenylether'
        case('propanol')
            solvent_name='propanol'
        case('pyridine')
            solvent_name='pyridine'
        case('secbutanol')
            solvent_name='secbutanol'
        case('secbutylbenzene')
            solvent_name='secbutylbenzene'
        case('tbutylbenzene')
            solvent_name='tbutylbenzene'
        case('tetrachloroethene','c2cl4')
            solvent_name='tetrachloroethene' 
        case('tetrahydrofuran')
            solvent_name='tetrahydrofuran'
        case('tetrahydrothiophenedioxide','sulfolan') 
            solvent_name='tetrahydrothiophenedioxide' 
        case('tetralin')
            solvent_name='tetralin'
        case('toluene')
            solvent_name='toluene'
        case('tributylphosphate')
            solvent_name='tributylphosphate'
        case('triethylamine')
            solvent_name='triethylamine'
        case('trimethylbenzene')
            solvent_name='trimethylbenzene'
        case('undecane')
            solvent_name='undecane'
        case('water','h2o')
            solvent_name='water'
        case('xylene')
            solvent_name='xylene'
        case('methanol')
            solvent_name='methanol'
        end select

    end function solvent_name
    !> Returns the density of a solvent in kg/m^3
    function density(solvent)
        !> Name of the solvent
        character(len=*), intent(in) :: solvent

        !> Density of the solvent
        real(wp) :: density
 
        select case(solvent_name(solvent))
        case('2methylpyridine')
            density=944.0_wp
        case('4methyl2pentanone')
            density=802.0_wp   
        case('aceticacid')
            density=1051.0_wp
        case('acetonitrile')
            density=787.0_wp
        case('acetophenone')
            density=1028.0_wp        
        case('aniline')
            density=1022.0_wp
        case('anisole')
            density=995.6_wp
        case('benzene')
            density=879.0_wp
        case('benzonitrile')
            density=1009.0_wp
        case('benzylalcohol')
            density=1041.0_wp
        case('bromobenzene')
            density=1490.0_wp
        case('bromoethane')
            density=1460.0_wp
        case('bromoform')
            density=2891.2_wp
        case('bromooctane')
            density=1118.0_wp
        case('butanol')
            density=810.0_wp
        case('butanone')
            density=806.0_wp
        case('butylacetate')
            density=921.2_wp
        case('butylbenzene')
            density=860.1_wp       
        case('carbondisulfide')
            density=1263.2_wp
        case('carbontet') !Carbontetrachloride
            density=1590.0_wp        
        case('chlorobenzene')
            density=1105.8_wp
        case('chloroform')
            density=1483.2_wp
        case('chlorohexane')
            density=879.0_wp
        case('cyclohexane')
            density=779.0_wp
        case('cyclohexanone')
            density=945.0_wp
        case('decalin')
            density=890.0_wp        
        case('decane')
            density=730.0_wp
        case('decanol')
            density=840.0_wp
        case('dibromoethane')
            density=2055.5_wp
        case('dibutylether')
            density=767.0_wp
        case('dichloroethane')
            density=1253.0_wp
        case('diethylether')
            density=714.0_wp       
        case('diisopropylether')
            density=724.0_wp       
        case('dimethylacetamide')
            density=943.0_wp
        case('dimethylformamide')
            density=950.0_wp
        case('dimethylpyridine')
            density=945.0_wp
        case('dimethylsulfoxide')
            density=1101.0_wp
        case('dodecane')
            density=748.7_wp        
        case('ethanol')
            density=790.0_wp
        case('ethoxybenzene')
            density=967.0_wp
        case('ethylacetate')
            density=902.0_wp
        case('ethylbenzene')
            density=867.0_wp
        case('fluorobenzene')
            density=1022.5_wp       
        case('fluoroctane')
            density=814.0_wp       
        case('heptane')
            density=683.8_wp
        case('heptanol')
            density=822.0_wp
        case('hexadecane')
            density=773.2_wp
        case('hexadecyliodide') !1-Iodohexadecane
            density=1121.0_wp
        case('hexane')
            density=659.0_wp
        case('hexanol')
            density=850.0_wp
        case('iodobenzene')
            density=1808.0_wp        
        case('isobutanol')
            density=802.0_wp
        case('isooctane')
            density=691.9_wp
        case('isopropanol')
            density=785.0_wp
        case('isopropylbenzene')
            density=866.0_wp
        case('isopropyltoluene')
            density=857.0_wp
        case('mcresol')
            density=1033.6_wp        
        case('mesitylene')
            density=860.0_wp        
        case('methoxyethanol')
            density=966.0_wp       
        case('methylenechloride')
            density=1322.0_wp       
        case('methylformamide')
            density=1011.0_wp
        case('nitrobenzene')
            density=1203.7_wp
        case('nitroethane')
            density=1050.0_wp
        case('nitromethane')
            density=1139.0_wp
        case('nonane')
            density=718.0_wp
        case('nonanol')
            density=827.0_wp
        case('octane')
            density=703.0_wp
        case('octanol')
            density=829.0_wp
        case('odichlorobenzene')
            density=1306.0_wp
        case('onitrotoluene')
            density=1162.2_wp
        case('pentadecane')
            density=768.5_wp
        case('pentane')
            density=626.0_wp
        case('pentanol')
            density=818.0_wp        
        case('perfluorobenzene') !Hexafluorobenzene
            density=1610.0_wp        
        case('phenylether')
            density=1070.0_wp       
        case('propanol')
            density=803.0_wp
        case('pyridine')
            density=983.0_wp
        case('secbutanol')
            density=810.0_wp
        case('secbutylbenzene')
            density=858.0_wp
        case('tbutylbenzene')
            density=866.9_wp
        case('tetrachloroethene') !C2Cl4
            density=1630.0_wp        
        case('tetrahydrofuran')
            density=888.0_wp
        case('tetrahydrothiophenedioxide') !Sulfolan
            density=1260.0_wp
        case('tetralin')
            density=974.0_wp        
        case('toluene')
            density=867.0_wp
        case('tributylphosphate')
            density=982.0_wp
        case('triethylamine')
            density=729.0_wp
        case('trimethylbenzene')
            density=894.4_wp
        case('undecane')
            density=740.1_wp        
        case('water','h2o')
            density=998.0_wp      
        case('xylene')
            density=880.0_wp
        case('methanol')
            density=792.0_wp
        case default
            density=998.0_wp
        end select

    end function density

    !> Get default dielectric constant from Minnesota Solvation Database
    function minnesota_eps(solvent,error) result(epsilon)
        character(len=*), intent(in) :: solvent
        type(error_type), allocatable :: error
        real(wp):: epsilon

        select case(solvent_name(solvent))
        case('infinity','inf')
            epsilon=ieee_value(epsilon,ieee_positive_inf)
        case('2methylpyridine')
            epsilon=9.9533_wp
        case('4methyl2pentanone')
            epsilon=12.8871
        case('aceticacid')
            epsilon=6.2528
        case('acetonitrile')
            epsilon=35.6881
        case('acetophenone')
            epsilon=17.44
        case('aniline')
            epsilon=6.8882
        case('anisole')
            epsilon=4.2247
        case('benzene')
            epsilon=2.2706
        case('benzonitrile')
            epsilon=25.592
        case('benzylalcohol')
            epsilon=12.4569
        case('bromobenzene')
            epsilon=5.3954
        case('bromoethane')
            epsilon=9.01
        case('bromoform')
            epsilon=4.2488
        case('bromooctane')
            epsilon=5.0244
        case('butanol')
            epsilon=17.3323
        case('butanone')
            epsilon=18.2457
        case('butylacetate')
            epsilon=4.9941
        case('butylbenzene')
            epsilon=2.36
        case('carbondisulfide')
            epsilon=2.6105
        case('carbontet')
            epsilon=2.228
        case('chlorobenzene')
            epsilon=5.6968
        case('chloroform')
            epsilon=4.7113
        case('chlorohexane')
            epsilon=5.9491
        case('cyclohexane')
            epsilon=2.0165
        case('cyclohexanone')
            epsilon=15.6186
        case('decalin')
            epsilon=2.196
        case('decane')
            epsilon=1.9846
        case('decanol')
            epsilon=7.5305
        case('dibromoethane')
            epsilon=4.9313
        case('dibutylether')
            epsilon=3.0473
        case('dichloroethane')
            epsilon=10.125
        case('diethylether')
            epsilon=4.24
        case('diisopropylether')
            epsilon=3.38
        case('dimethylacetamide')
            epsilon=37.7807
        case('dimethylformamide')
            epsilon=37.219
        case('dimethylpyridine')
            epsilon=7.1735
        case('dimethylsulfoxide')
            epsilon=46.826
        case('dodecane')
            epsilon=2.006
        case('ethanol')
            epsilon=24.852
        case('ethoxybenzene')
            epsilon=4.1797
        case('ethylacetate')
            epsilon=5.9867
        case('ethylbenzene')
            epsilon=2.4339
        case('fluorobenzene')
            epsilon=5.42
        case('fluoroctane')
            epsilon=3.89
        case('heptane')
            epsilon=1.9113
        case('heptanol')
            epsilon=11.321
        case('hexadecane')
            epsilon=2.0402
        case('hexadecyliodide')
            epsilon=3.5338
        case('hexane')
            epsilon=1.8819
        case('hexanol')
            epsilon=12.5102
        case('iodobenzene')
            epsilon=4.547
        case('isobutanol')
            epsilon=16.7766
        case('isooctane')
            epsilon=1.9358
        case('isopropanol')
            epsilon=19.2645
        case('isopropylbenzene')
            epsilon=2.3712
        case('isopropyltoluene')
            epsilon=2.2322
        case('mcresol')
            epsilon=12.44
        case('mesitylene')
            epsilon=2.265
        case('methoxyethanol')
            epsilon=17.2
        case('methylenechloride')
            epsilon=8.93
        case('methylformamide')
            epsilon=181.5619
        case('nitrobenzene')
            epsilon=34.8091
        case('nitroethane')
            epsilon=28.2896
        case('nitromethane')
            epsilon=36.5623
        case('nonane')
            epsilon=1.9605
        case('nonanol')
            epsilon=8.5991
        case('octane')
            epsilon=1.9406
        case('octanol')
            epsilon=9.8629
        case('odichlorobenzene')
            epsilon=9.9949
        case('onitrotoluene')
            epsilon=25.6692
        case('pentadecane')
            epsilon=2.0333
        case('pentane')
            epsilon=1.8371
        case('pentanol')
            epsilon=15.13
        case('perfluorobenzene')
            epsilon=2.029
        case('phenylether')
            epsilon=3.73
        case('propanol')
            epsilon=20.5237
        case('pyridine')
            epsilon=12.9776
        case('secbutanol')
            epsilon=15.9436
        case('secbutylbenzene')
            epsilon=2.3446
        case('tbutylbenzene')
            epsilon=2.3447
        case('tetrachloroethene')
            epsilon=2.268
        case('tetrahydrofuran')
            epsilon=7.4257
        case('tetrahydrothiophenedioxide')
            epsilon=43.9622
        case('tetralin')
            epsilon=2.771
        case('toluene')
            epsilon=2.3741
        case('tributylphosphate')
            epsilon=8.1781
        case('triethylamine')
            epsilon=2.3832
        case('trimethylbenzene')
            epsilon=2.3653
        case('undecane')
            epsilon=1.991
        case('water','h2o')
            epsilon=78.36_wp
        case('xylene')
            epsilon=2.3879
        case('benzene-water')
            epsilon=2.2706
        case('carbontet-water')
            epsilon=2.228
        case('chlorobenzene-water')
            epsilon=5.6968
        case('chloroform-water')
            epsilon=4.7113
        case('cyclohexane-water')
            epsilon=2.0165
        case('dibromoethane-water')
            epsilon=4.9313
        case('dibutylether-water')
            epsilon=3.0473
        case('dichloroethane-water')
            epsilon=10.125
        case('diethylether-water')
            epsilon=4.24
        case('ethylacetate-water')
            epsilon=5.9867
        case('heptane-water')
            epsilon=1.9113
        case('hexane-water')
            epsilon=1.8819
        case('nitrobenzene-water')
            epsilon=34.8091
        case('octanol-water')
            epsilon=9.8629
        case('methanol')
            epsilon=32.613
        case default
            Call fatal_error(error,'epsilon_water: unknown solvent '//solvent)
        end select
    end function minnesota_eps
    
    !> Get atomic mass for an array of Element Symbols
    function getAtomicMassSymbolArray(atoms) result (mass)
        !> Array of Element Symbols
        character(2), dimension(:), allocatable, intent(in) :: atoms
        !> Atomic Mass
        real(wp) :: mass
        !> Laufvariable
        integer :: i

        mass=0.0_wp
        do i=1,size(atoms)
            mass=mass+getAtomicMassSymbol(atoms(i))
        end do

    end function getAtomicMassSymbolArray


    !> Get atomic mass for an array of Element Symbols
    function getAtomicMassNumberArray(atoms) result (mass)
        !> Array of Element Numbers
        integer, dimension(:), allocatable, intent(in) :: atoms
        !> Atomic Mass
        real(wp) :: mass
        !> Laufvariable
        integer :: i

        mass=0.0_wp
        do i=1,size(atoms)
            mass=mass+getAtomicMassNumber(atoms(i))
        end do

    end function getAtomicMassNumberArray
        
  !> Get atomic mass for species with a given symbol
    elemental function getAtomicMassSymbol(symbol) result(mass)

    !> Element symbol
    character(len=*), intent(in) :: symbol

    !> atomic mass
    real(wp) :: mass

    mass = getAtomicMassNumber(to_number(symbol))

  end function getAtomicMassSymbol


  !> Get atomic mass for species with a given atomic number
  elemental function getAtomicMassNumber(number) result(mass)

    !> Atomic number
    integer, intent(in) :: number

    !> atomic mass
    real(wp) :: mass

    if (number > 0 .and. number <= size(atomicMassNist, dim=1)) then
      mass = atomicMassNist(number)
    else
      mass = -1.0_wp
    end if

  end function getAtomicMassNumber

end module data
