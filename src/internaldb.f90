module internaldb
use mctc_env, only: wp
use type, only: parameter_type
use globals, only: btoa

type(parameter_type), parameter :: xtb_water = parameter_type (&
&0.30880726_wp, &
&7540.06830745_wp, &
&-1.26119686_wp, &
&15813.08917624_wp, &
&0.00845974_wp, &
&16.58382447_wp, &
&-0.23040291_wp, &
&-0.13909842_wp, &
&-10.60951107_wp,&
&0.0_wp)

type(parameter_type), parameter :: xtb_other = parameter_type (&
&0.53725564_wp, &
&1977.88495090_wp, &
&1.61721497_wp, &
&2664.61735599_wp, &
&0.00806820_wp, &
&57.11402448_wp, &
&-0.12356319_wp, &
&-0.08426734_wp, &
&-8.41974122_wp,&
&0.0_wp)

character(len=200),dimension(41), parameter :: xtb_water_smd = [character(len=200) ::& 
&"#Zk - This are the default Parmeters for SMD/H2O", & 
&"H 	-77.70436564", & 
&"C 	-44.32921125", & 
&"N 	-25.97069034", & 
&"O 	-50.70232515", & 
&"F 	-6.43754644", & 
&"Cl	-36.05536752", & 
&"Br	-20.79864972", & 
&"S 	-18.13699434", & 
&"#Zkk Matrix", & 
&"H	00.00	76.79778940	00.00	00.00	00.00", & 
&"C	00.00	-1.94928804	00.00	00.00	00.00", & 
&"O	00.00	12.35744639	00.00	165.85557924	12.23268205", & 
&"N	00.00	-24.36906946	00.00	00.00	00.00", & 
&"P	00.00	00.00	00.00	00.00	00.00", & 
&"#rzkk Matrix", & 
&"H  00.00  1.84000000  00.00  1.55000000  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"C  1.55000000  1.84000000  1.84000000  1.84000000  1.84000000  2.20000000  2.20000000  2.10000000  2.30000000  2.60000000", & 
&"N  00.00  1.84000000  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"O  00.00  1.33000000  1.50000000  1.80000000  00.00  2.10000000  00.00  00.00  00.00  00.00", & 
&"F  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"P  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"S  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"Cl  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"Br  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"I  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"#drzkk Matrix", & 
&"H  00.00  0.30000000  00.00  0.30000000  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"C  0.30000000  0.30000000  0.30000000  0.30000000  0.30000000  0.30000000  0.30000000  0.30000000  0.30000000  0.30000000", & 
&"N  00.00  0.30000000  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"O  00.00  0.10000000  0.30000000  0.30000000  0.30000000  00.00  00.00  00.00  00.00  00.00", & 
&"F  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"P  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"S  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"Cl  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"Br  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"I  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"#NC3", & 
&"33.62391370", & 
&"1.22500000", & 
&"0.06500000"]

character(len=200),dimension(58), parameter :: xtb_other_smd = [character(len=200) ::& 
&"#Zk", & 
&"C	-3.70891582", & 
&"O	38.16986455", & 
&"N	11.67657766", & 
&"Cl	-26.23511892", & 
&"Br	-19.44296237", & 
&"S	-15.85282828", & 
&"H	-5.61002654", & 
&"P	-31.18737493", & 
&"#Zk2", & 
&"C	8.80543594", & 
&"O	-2.68778654", & 
&"#Zk3", & 
&"C	-5.80252796", & 
&"O	-87.53912863", & 
&"#Zkk Matrix", & 
&"H	00.00	-6.10646558	-72.27657574	00.00", & 
&"C	00.00	-15.05405921	00.00	-16.54455408", & 
&"O	00.00	-20.55673110	00.00	00.00", & 
&"N	00.00	00.00	00.00	00.00", & 
&"#Zkk2 Matrix", & 
&"O	00.00	65.41437683	00.00", & 
&"C	00.00	00.00	-2.81967860", & 
&"N	00.00	-32.31978386	00.00", & 
&"#Zkk3 Matrix", & 
&"O	-40.24403232	67.09645304", & 
&"N	00.00	00.00", & 
&"#rzkk Matrix", & 
&"H  00.00  1.55100494  00.00  1.53632796  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"C  1.54582282  1.83993169  1.83959886  2.07175304  1.84000000  2.24400000  2.20000000  2.10000000  2.30000000  2.60000000", & 
&"N  00.00  1.84000000  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"O  00.00  1.33181762  1.50000000  1.80000000  00.00  2.10000000  00.00  00.00  00.00  00.00", & 
&"F  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"P  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"S  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"Cl  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"Br  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"I  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"#drzkk Matrix", & 
&"H  00.00  0.30023238  00.00  0.30023238  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"C  0.30023238  0.30023238  0.30023238  0.30023238  0.30023238  0.30023238  0.30023238  0.30023238  0.30023238  0.30023238", & 
&"N  00.00  0.30023238  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"O  00.00  0.11286083  0.30023238  0.30023238  0.30023238  00.00  00.00  00.00  00.00  00.00", & 
&"F  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"P  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"S  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"Cl  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"Br  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"I  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00", & 
&"#NC3", & 
&"00.00", & 
&"1.22500000", & 
&"0.06500000", & 
&"#Solvent", & 
&"0.08563085", & 
&"-7.21289853", & 
&"-14.52469035", & 
&"-5.13736584"]

include "DB/xtb/2methylpyridine.fh"
include "DB/xtb/4methyl2pentanone.fh"
include "DB/xtb/aceticacid.fh"
include "DB/xtb/acetonitrile.fh"
include "DB/xtb/acetophenone.fh"
include "DB/xtb/aniline.fh"
include "DB/xtb/anisole.fh"
include "DB/xtb/benzene.fh"
include "DB/xtb/benzylalcohol.fh"
include "DB/xtb/bromobenzene.fh"
include "DB/xtb/bromoethane.fh"
include "DB/xtb/bromoform.fh"
include "DB/xtb/bromooctane.fh"
include "DB/xtb/butanol.fh"
include "DB/xtb/butanone.fh"
include "DB/xtb/butylacetate.fh"
include "DB/xtb/butylbenzene.fh"
include "DB/xtb/carbondisulfide.fh"
include "DB/xtb/carbontet.fh"
include "DB/xtb/chlorobenzene.fh"
include "DB/xtb/chloroform.fh"
include "DB/xtb/chlorohexane.fh"
include "DB/xtb/cyclohexane.fh"
include "DB/xtb/cyclohexanone.fh"
include "DB/xtb/decalin.fh"
include "DB/xtb/decane.fh"
include "DB/xtb/decanol.fh"
include "DB/xtb/dibromoethane.fh"
include "DB/xtb/dibutylether.fh"
include "DB/xtb/dichloroethane.fh"
include "DB/xtb/diethylether.fh"
include "DB/xtb/diisopropylether.fh"
include "DB/xtb/dimethylacetamide.fh"
include "DB/xtb/dimethylformamide.fh"
include "DB/xtb/dimethylpyridine.fh"
include "DB/xtb/dimethylsulfoxide.fh"
include "DB/xtb/dodecane.fh"
include "DB/xtb/ethanol.fh"
include "DB/xtb/ethoxybenzene.fh"
include "DB/xtb/ethylacetate.fh"
include "DB/xtb/ethylbenzene.fh"
include "DB/xtb/fluorobenzene.fh"
include "DB/xtb/fluoroctane.fh"
include "DB/xtb/heptane.fh"
include "DB/xtb/hexadecane.fh"
include "DB/xtb/hexadecyliodide.fh"
include "DB/xtb/hexane.fh"
include "DB/xtb/hexanol.fh"
include "DB/xtb/iodobenzene.fh"
include "DB/xtb/isobutanol.fh"
include "DB/xtb/isooctane.fh"
include "DB/xtb/isopropanol.fh"
include "DB/xtb/isopropylbenzene.fh"
include "DB/xtb/isopropyltoluene.fh"
include "DB/xtb/mcresol.fh"
include "DB/xtb/mesitylene.fh"
include "DB/xtb/methoxyethanol.fh"
include "DB/xtb/methylenechloride.fh"
include "DB/xtb/methylformamide.fh"
include "DB/xtb/nitrobenzene.fh"
include "DB/xtb/nitroethane.fh"
include "DB/xtb/nitromethane.fh"
include "DB/xtb/nonane.fh"
include "DB/xtb/nonanol.fh"
include "DB/xtb/octane.fh"
include "DB/xtb/octanol.fh"
include "DB/xtb/odichlorobenzene.fh"
include "DB/xtb/onitrotoluene.fh"
include "DB/xtb/pentadecane.fh"
include "DB/xtb/pentane.fh"
include "DB/xtb/pentanol.fh"
include "DB/xtb/perfluorobenzene.fh"
include "DB/xtb/phenylether.fh"
include "DB/xtb/propanol.fh"
include "DB/xtb/pyridine.fh"
include "DB/xtb/secbutanol.fh"
include "DB/xtb/secbutylbenzene.fh"
include "DB/xtb/tbutylbenzene.fh"
include "DB/xtb/tetrachloroethene.fh"
include "DB/xtb/tetrahydrofuran.fh"
include "DB/xtb/tetrahydrothiophenedioxide.fh"
include "DB/xtb/tetralin.fh"
include "DB/xtb/toluene.fh"
include "DB/xtb/tributylphosphate.fh"
include "DB/xtb/triethylamine.fh"
include "DB/xtb/trimethylbenzene.fh"
include "DB/xtb/undecane.fh"
include "DB/xtb/water.fh"


contains

subroutine internalcosmo(solvent,self,error)
   use type, only: molecule_data
   use mctc_env, only: error_type, fatal_error
   character(len=*), intent(in) :: solvent
   type(molecule_data) :: self
   type(error_type), intent(out), allocatable :: error

   select case(solvent)
   case default
      Call fatal_error(error,"Solvent "//solvent//" was not found in the internal database.")
      return
   case("2methylpyridine")
      self%id=sid_2methylpyridine_xtb
      self%area=sarea_2methylpyridine_xtb
      self%su=su_2methylpyridine_xtb
      self%atom_xyz=axyz_2methylpyridine_xtb
      self%element=aelement_2methylpyridine_xtb
      self%energy=energy_2methylpyridine_xtb
      self%xyz=sxyz_2methylpyridine_xtb
   case("4methyl2pentanone")
      self%id=sid_4methyl2pentanone_xtb
      self%area=sarea_4methyl2pentanone_xtb
      self%su=su_4methyl2pentanone_xtb
      self%atom_xyz=axyz_4methyl2pentanone_xtb
      self%element=aelement_4methyl2pentanone_xtb
      self%energy=energy_4methyl2pentanone_xtb
      self%xyz=sxyz_4methyl2pentanone_xtb
   case("aceticacid")
      self%id=sid_aceticacid_xtb
      self%area=sarea_aceticacid_xtb
      self%su=su_aceticacid_xtb
      self%atom_xyz=axyz_aceticacid_xtb
      self%element=aelement_aceticacid_xtb
      self%energy=energy_aceticacid_xtb
      self%xyz=sxyz_aceticacid_xtb
   case("acetonitrile")
      self%id=sid_acetonitrile_xtb
      self%area=sarea_acetonitrile_xtb
      self%su=su_acetonitrile_xtb
      self%atom_xyz=axyz_acetonitrile_xtb
      self%element=aelement_acetonitrile_xtb
      self%energy=energy_acetonitrile_xtb
      self%xyz=sxyz_acetonitrile_xtb
   case("acetophenone")
      self%id=sid_acetophenone_xtb
      self%area=sarea_acetophenone_xtb
      self%su=su_acetophenone_xtb
      self%atom_xyz=axyz_acetophenone_xtb
      self%element=aelement_acetophenone_xtb
      self%energy=energy_acetophenone_xtb
      self%xyz=sxyz_acetophenone_xtb
   case("aniline")
      self%id=sid_aniline_xtb
      self%area=sarea_aniline_xtb
      self%su=su_aniline_xtb
      self%atom_xyz=axyz_aniline_xtb
      self%element=aelement_aniline_xtb
      self%energy=energy_aniline_xtb
      self%xyz=sxyz_aniline_xtb
   case("anisole")
      self%id=sid_anisole_xtb
      self%area=sarea_anisole_xtb
      self%su=su_anisole_xtb
      self%atom_xyz=axyz_anisole_xtb
      self%element=aelement_anisole_xtb
      self%energy=energy_anisole_xtb
      self%xyz=sxyz_anisole_xtb
   case("benzene")
      self%id=sid_benzene_xtb
      self%area=sarea_benzene_xtb
      self%su=su_benzene_xtb
      self%atom_xyz=axyz_benzene_xtb
      self%element=aelement_benzene_xtb
      self%energy=energy_benzene_xtb
      self%xyz=sxyz_benzene_xtb
   case("benzylalcohol")
      self%id=sid_benzylalcohol_xtb
      self%area=sarea_benzylalcohol_xtb
      self%su=su_benzylalcohol_xtb
      self%atom_xyz=axyz_benzylalcohol_xtb
      self%element=aelement_benzylalcohol_xtb
      self%energy=energy_benzylalcohol_xtb
      self%xyz=sxyz_benzylalcohol_xtb
   case("bromobenzene")
      self%id=sid_bromobenzene_xtb
      self%area=sarea_bromobenzene_xtb
      self%su=su_bromobenzene_xtb
      self%atom_xyz=axyz_bromobenzene_xtb
      self%element=aelement_bromobenzene_xtb
      self%energy=energy_bromobenzene_xtb
      self%xyz=sxyz_bromobenzene_xtb
   case("bromoethane")
      self%id=sid_bromoethane_xtb
      self%area=sarea_bromoethane_xtb
      self%su=su_bromoethane_xtb
      self%atom_xyz=axyz_bromoethane_xtb
      self%element=aelement_bromoethane_xtb
      self%energy=energy_bromoethane_xtb
      self%xyz=sxyz_bromoethane_xtb
   case("bromoform")
      self%id=sid_bromoform_xtb
      self%area=sarea_bromoform_xtb
      self%su=su_bromoform_xtb
      self%atom_xyz=axyz_bromoform_xtb
      self%element=aelement_bromoform_xtb
      self%energy=energy_bromoform_xtb
      self%xyz=sxyz_bromoform_xtb
   case("bromooctane")
      self%id=sid_bromooctane_xtb
      self%area=sarea_bromooctane_xtb
      self%su=su_bromooctane_xtb
      self%atom_xyz=axyz_bromooctane_xtb
      self%element=aelement_bromooctane_xtb
      self%energy=energy_bromooctane_xtb
      self%xyz=sxyz_bromooctane_xtb
   case("butanol")
      self%id=sid_butanol_xtb
      self%area=sarea_butanol_xtb
      self%su=su_butanol_xtb
      self%atom_xyz=axyz_butanol_xtb
      self%element=aelement_butanol_xtb
      self%energy=energy_butanol_xtb
      self%xyz=sxyz_butanol_xtb
   case("butanone")
      self%id=sid_butanone_xtb
      self%area=sarea_butanone_xtb
      self%su=su_butanone_xtb
      self%atom_xyz=axyz_butanone_xtb
      self%element=aelement_butanone_xtb
      self%energy=energy_butanone_xtb
      self%xyz=sxyz_butanone_xtb
   case("butylacetate")
      self%id=sid_butylacetate_xtb
      self%area=sarea_butylacetate_xtb
      self%su=su_butylacetate_xtb
      self%atom_xyz=axyz_butylacetate_xtb
      self%element=aelement_butylacetate_xtb
      self%energy=energy_butylacetate_xtb
      self%xyz=sxyz_butylacetate_xtb
   case("butylbenzene")
      self%id=sid_butylbenzene_xtb
      self%area=sarea_butylbenzene_xtb
      self%su=su_butylbenzene_xtb
      self%atom_xyz=axyz_butylbenzene_xtb
      self%element=aelement_butylbenzene_xtb
      self%energy=energy_butylbenzene_xtb
      self%xyz=sxyz_butylbenzene_xtb
   case("carbondisulfide")
      self%id=sid_carbondisulfide_xtb
      self%area=sarea_carbondisulfide_xtb
      self%su=su_carbondisulfide_xtb
      self%atom_xyz=axyz_carbondisulfide_xtb
      self%element=aelement_carbondisulfide_xtb
      self%energy=energy_carbondisulfide_xtb
      self%xyz=sxyz_carbondisulfide_xtb
   case("carbontet")
      self%id=sid_carbontet_xtb
      self%area=sarea_carbontet_xtb
      self%su=su_carbontet_xtb
      self%atom_xyz=axyz_carbontet_xtb
      self%element=aelement_carbontet_xtb
      self%energy=energy_carbontet_xtb
      self%xyz=sxyz_carbontet_xtb
   case("chlorobenzene")
      self%id=sid_chlorobenzene_xtb
      self%area=sarea_chlorobenzene_xtb
      self%su=su_chlorobenzene_xtb
      self%atom_xyz=axyz_chlorobenzene_xtb
      self%element=aelement_chlorobenzene_xtb
      self%energy=energy_chlorobenzene_xtb
      self%xyz=sxyz_chlorobenzene_xtb
   case("chloroform")
      self%id=sid_chloroform_xtb
      self%area=sarea_chloroform_xtb
      self%su=su_chloroform_xtb
      self%atom_xyz=axyz_chloroform_xtb
      self%element=aelement_chloroform_xtb
      self%energy=energy_chloroform_xtb
      self%xyz=sxyz_chloroform_xtb
   case("chlorohexane")
      self%id=sid_chlorohexane_xtb
      self%area=sarea_chlorohexane_xtb
      self%su=su_chlorohexane_xtb
      self%atom_xyz=axyz_chlorohexane_xtb
      self%element=aelement_chlorohexane_xtb
      self%energy=energy_chlorohexane_xtb
      self%xyz=sxyz_chlorohexane_xtb
   case("cyclohexane")
      self%id=sid_cyclohexane_xtb
      self%area=sarea_cyclohexane_xtb
      self%su=su_cyclohexane_xtb
      self%atom_xyz=axyz_cyclohexane_xtb
      self%element=aelement_cyclohexane_xtb
      self%energy=energy_cyclohexane_xtb
      self%xyz=sxyz_cyclohexane_xtb
   case("cyclohexanone")
      self%id=sid_cyclohexanone_xtb
      self%area=sarea_cyclohexanone_xtb
      self%su=su_cyclohexanone_xtb
      self%atom_xyz=axyz_cyclohexanone_xtb
      self%element=aelement_cyclohexanone_xtb
      self%energy=energy_cyclohexanone_xtb
      self%xyz=sxyz_cyclohexanone_xtb
   case("decalin")
      self%id=sid_decalin_xtb
      self%area=sarea_decalin_xtb
      self%su=su_decalin_xtb
      self%atom_xyz=axyz_decalin_xtb
      self%element=aelement_decalin_xtb
      self%energy=energy_decalin_xtb
      self%xyz=sxyz_decalin_xtb
   case("decane")
      self%id=sid_decane_xtb
      self%area=sarea_decane_xtb
      self%su=su_decane_xtb
      self%atom_xyz=axyz_decane_xtb
      self%element=aelement_decane_xtb
      self%energy=energy_decane_xtb
      self%xyz=sxyz_decane_xtb
   case("decanol")
      self%id=sid_decanol_xtb
      self%area=sarea_decanol_xtb
      self%su=su_decanol_xtb
      self%atom_xyz=axyz_decanol_xtb
      self%element=aelement_decanol_xtb
      self%energy=energy_decanol_xtb
      self%xyz=sxyz_decanol_xtb
   case("dibromoethane")
      self%id=sid_dibromoethane_xtb
      self%area=sarea_dibromoethane_xtb
      self%su=su_dibromoethane_xtb
      self%atom_xyz=axyz_dibromoethane_xtb
      self%element=aelement_dibromoethane_xtb
      self%energy=energy_dibromoethane_xtb
      self%xyz=sxyz_dibromoethane_xtb
   case("dibutylether")
      self%id=sid_dibutylether_xtb
      self%area=sarea_dibutylether_xtb
      self%su=su_dibutylether_xtb
      self%atom_xyz=axyz_dibutylether_xtb
      self%element=aelement_dibutylether_xtb
      self%energy=energy_dibutylether_xtb
      self%xyz=sxyz_dibutylether_xtb
   case("dichloroethane")
      self%id=sid_dichloroethane_xtb
      self%area=sarea_dichloroethane_xtb
      self%su=su_dichloroethane_xtb
      self%atom_xyz=axyz_dichloroethane_xtb
      self%element=aelement_dichloroethane_xtb
      self%energy=energy_dichloroethane_xtb
      self%xyz=sxyz_dichloroethane_xtb
   case("diethylether")
      self%id=sid_diethylether_xtb
      self%area=sarea_diethylether_xtb
      self%su=su_diethylether_xtb
      self%atom_xyz=axyz_diethylether_xtb
      self%element=aelement_diethylether_xtb
      self%energy=energy_diethylether_xtb
      self%xyz=sxyz_diethylether_xtb
   case("diisopropylether")
      self%id=sid_diisopropylether_xtb
      self%area=sarea_diisopropylether_xtb
      self%su=su_diisopropylether_xtb
      self%atom_xyz=axyz_diisopropylether_xtb
      self%element=aelement_diisopropylether_xtb
      self%energy=energy_diisopropylether_xtb
      self%xyz=sxyz_diisopropylether_xtb
   case("dimethylacetamide")
      self%id=sid_dimethylacetamide_xtb
      self%area=sarea_dimethylacetamide_xtb
      self%su=su_dimethylacetamide_xtb
      self%atom_xyz=axyz_dimethylacetamide_xtb
      self%element=aelement_dimethylacetamide_xtb
      self%energy=energy_dimethylacetamide_xtb
      self%xyz=sxyz_dimethylacetamide_xtb
   case("dimethylformamide")
      self%id=sid_dimethylformamide_xtb
      self%area=sarea_dimethylformamide_xtb
      self%su=su_dimethylformamide_xtb
      self%atom_xyz=axyz_dimethylformamide_xtb
      self%element=aelement_dimethylformamide_xtb
      self%energy=energy_dimethylformamide_xtb
      self%xyz=sxyz_dimethylformamide_xtb
   case("dimethylpyridine")
      self%id=sid_dimethylpyridine_xtb
      self%area=sarea_dimethylpyridine_xtb
      self%su=su_dimethylpyridine_xtb
      self%atom_xyz=axyz_dimethylpyridine_xtb
      self%element=aelement_dimethylpyridine_xtb
      self%energy=energy_dimethylpyridine_xtb
      self%xyz=sxyz_dimethylpyridine_xtb
   case("dimethylsulfoxide")
      self%id=sid_dimethylsulfoxide_xtb
      self%area=sarea_dimethylsulfoxide_xtb
      self%su=su_dimethylsulfoxide_xtb
      self%atom_xyz=axyz_dimethylsulfoxide_xtb
      self%element=aelement_dimethylsulfoxide_xtb
      self%energy=energy_dimethylsulfoxide_xtb
      self%xyz=sxyz_dimethylsulfoxide_xtb
   case("dodecane")
      self%id=sid_dodecane_xtb
      self%area=sarea_dodecane_xtb
      self%su=su_dodecane_xtb
      self%atom_xyz=axyz_dodecane_xtb
      self%element=aelement_dodecane_xtb
      self%energy=energy_dodecane_xtb
      self%xyz=sxyz_dodecane_xtb
   case("ethanol")
      self%id=sid_ethanol_xtb
      self%area=sarea_ethanol_xtb
      self%su=su_ethanol_xtb
      self%atom_xyz=axyz_ethanol_xtb
      self%element=aelement_ethanol_xtb
      self%energy=energy_ethanol_xtb
      self%xyz=sxyz_ethanol_xtb
   case("ethoxybenzene")
      self%id=sid_ethoxybenzene_xtb
      self%area=sarea_ethoxybenzene_xtb
      self%su=su_ethoxybenzene_xtb
      self%atom_xyz=axyz_ethoxybenzene_xtb
      self%element=aelement_ethoxybenzene_xtb
      self%energy=energy_ethoxybenzene_xtb
      self%xyz=sxyz_ethoxybenzene_xtb
   case("ethylacetate")
      self%id=sid_ethylacetate_xtb
      self%area=sarea_ethylacetate_xtb
      self%su=su_ethylacetate_xtb
      self%atom_xyz=axyz_ethylacetate_xtb
      self%element=aelement_ethylacetate_xtb
      self%energy=energy_ethylacetate_xtb
      self%xyz=sxyz_ethylacetate_xtb
   case("ethylbenzene")
      self%id=sid_ethylbenzene_xtb
      self%area=sarea_ethylbenzene_xtb
      self%su=su_ethylbenzene_xtb
      self%atom_xyz=axyz_ethylbenzene_xtb
      self%element=aelement_ethylbenzene_xtb
      self%energy=energy_ethylbenzene_xtb
      self%xyz=sxyz_ethylbenzene_xtb
   case("fluorobenzene")
      self%id=sid_fluorobenzene_xtb
      self%area=sarea_fluorobenzene_xtb
      self%su=su_fluorobenzene_xtb
      self%atom_xyz=axyz_fluorobenzene_xtb
      self%element=aelement_fluorobenzene_xtb
      self%energy=energy_fluorobenzene_xtb
      self%xyz=sxyz_fluorobenzene_xtb
   case("fluoroctane")
      self%id=sid_fluoroctane_xtb
      self%area=sarea_fluoroctane_xtb
      self%su=su_fluoroctane_xtb
      self%atom_xyz=axyz_fluoroctane_xtb
      self%element=aelement_fluoroctane_xtb
      self%energy=energy_fluoroctane_xtb
      self%xyz=sxyz_fluoroctane_xtb
   case("heptane")
      self%id=sid_heptane_xtb
      self%area=sarea_heptane_xtb
      self%su=su_heptane_xtb
      self%atom_xyz=axyz_heptane_xtb
      self%element=aelement_heptane_xtb
      self%energy=energy_heptane_xtb
      self%xyz=sxyz_heptane_xtb
   case("hexadecane")
      self%id=sid_hexadecane_xtb
      self%area=sarea_hexadecane_xtb
      self%su=su_hexadecane_xtb
      self%atom_xyz=axyz_hexadecane_xtb
      self%element=aelement_hexadecane_xtb
      self%energy=energy_hexadecane_xtb
      self%xyz=sxyz_hexadecane_xtb
   case("hexadecyliodide")
      self%id=sid_hexadecyliodide_xtb
      self%area=sarea_hexadecyliodide_xtb
      self%su=su_hexadecyliodide_xtb
      self%atom_xyz=axyz_hexadecyliodide_xtb
      self%element=aelement_hexadecyliodide_xtb
      self%energy=energy_hexadecyliodide_xtb
      self%xyz=sxyz_hexadecyliodide_xtb
   case("hexane")
      self%id=sid_hexane_xtb
      self%area=sarea_hexane_xtb
      self%su=su_hexane_xtb
      self%atom_xyz=axyz_hexane_xtb
      self%element=aelement_hexane_xtb
      self%energy=energy_hexane_xtb
      self%xyz=sxyz_hexane_xtb
   case("hexanol")
      self%id=sid_hexanol_xtb
      self%area=sarea_hexanol_xtb
      self%su=su_hexanol_xtb
      self%atom_xyz=axyz_hexanol_xtb
      self%element=aelement_hexanol_xtb
      self%energy=energy_hexanol_xtb
      self%xyz=sxyz_hexanol_xtb
   case("iodobenzene")
      self%id=sid_iodobenzene_xtb
      self%area=sarea_iodobenzene_xtb
      self%su=su_iodobenzene_xtb
      self%atom_xyz=axyz_iodobenzene_xtb
      self%element=aelement_iodobenzene_xtb
      self%energy=energy_iodobenzene_xtb
      self%xyz=sxyz_iodobenzene_xtb
   case("isobutanol")
      self%id=sid_isobutanol_xtb
      self%area=sarea_isobutanol_xtb
      self%su=su_isobutanol_xtb
      self%atom_xyz=axyz_isobutanol_xtb
      self%element=aelement_isobutanol_xtb
      self%energy=energy_isobutanol_xtb
      self%xyz=sxyz_isobutanol_xtb
   case("isooctane")
      self%id=sid_isooctane_xtb
      self%area=sarea_isooctane_xtb
      self%su=su_isooctane_xtb
      self%atom_xyz=axyz_isooctane_xtb
      self%element=aelement_isooctane_xtb
      self%energy=energy_isooctane_xtb
      self%xyz=sxyz_isooctane_xtb
   case("isopropanol")
      self%id=sid_isopropanol_xtb
      self%area=sarea_isopropanol_xtb
      self%su=su_isopropanol_xtb
      self%atom_xyz=axyz_isopropanol_xtb
      self%element=aelement_isopropanol_xtb
      self%energy=energy_isopropanol_xtb
      self%xyz=sxyz_isopropanol_xtb
   case("isopropylbenzene")
      self%id=sid_isopropylbenzene_xtb
      self%area=sarea_isopropylbenzene_xtb
      self%su=su_isopropylbenzene_xtb
      self%atom_xyz=axyz_isopropylbenzene_xtb
      self%element=aelement_isopropylbenzene_xtb
      self%energy=energy_isopropylbenzene_xtb
      self%xyz=sxyz_isopropylbenzene_xtb
   case("isopropyltoluene")
      self%id=sid_isopropyltoluene_xtb
      self%area=sarea_isopropyltoluene_xtb
      self%su=su_isopropyltoluene_xtb
      self%atom_xyz=axyz_isopropyltoluene_xtb
      self%element=aelement_isopropyltoluene_xtb
      self%energy=energy_isopropyltoluene_xtb
      self%xyz=sxyz_isopropyltoluene_xtb
   case("mcresol")
      self%id=sid_mcresol_xtb
      self%area=sarea_mcresol_xtb
      self%su=su_mcresol_xtb
      self%atom_xyz=axyz_mcresol_xtb
      self%element=aelement_mcresol_xtb
      self%energy=energy_mcresol_xtb
      self%xyz=sxyz_mcresol_xtb
   case("mesitylene")
      self%id=sid_mesitylene_xtb
      self%area=sarea_mesitylene_xtb
      self%su=su_mesitylene_xtb
      self%atom_xyz=axyz_mesitylene_xtb
      self%element=aelement_mesitylene_xtb
      self%energy=energy_mesitylene_xtb
      self%xyz=sxyz_mesitylene_xtb
   case("methoxyethanol")
      self%id=sid_methoxyethanol_xtb
      self%area=sarea_methoxyethanol_xtb
      self%su=su_methoxyethanol_xtb
      self%atom_xyz=axyz_methoxyethanol_xtb
      self%element=aelement_methoxyethanol_xtb
      self%energy=energy_methoxyethanol_xtb
      self%xyz=sxyz_methoxyethanol_xtb
   case("methylenechloride")
      self%id=sid_methylenechloride_xtb
      self%area=sarea_methylenechloride_xtb
      self%su=su_methylenechloride_xtb
      self%atom_xyz=axyz_methylenechloride_xtb
      self%element=aelement_methylenechloride_xtb
      self%energy=energy_methylenechloride_xtb
      self%xyz=sxyz_methylenechloride_xtb
   case("methylformamide")
      self%id=sid_methylformamide_xtb
      self%area=sarea_methylformamide_xtb
      self%su=su_methylformamide_xtb
      self%atom_xyz=axyz_methylformamide_xtb
      self%element=aelement_methylformamide_xtb
      self%energy=energy_methylformamide_xtb
      self%xyz=sxyz_methylformamide_xtb
   case("nitrobenzene")
      self%id=sid_nitrobenzene_xtb
      self%area=sarea_nitrobenzene_xtb
      self%su=su_nitrobenzene_xtb
      self%atom_xyz=axyz_nitrobenzene_xtb
      self%element=aelement_nitrobenzene_xtb
      self%energy=energy_nitrobenzene_xtb
      self%xyz=sxyz_nitrobenzene_xtb
   case("nitroethane")
      self%id=sid_nitroethane_xtb
      self%area=sarea_nitroethane_xtb
      self%su=su_nitroethane_xtb
      self%atom_xyz=axyz_nitroethane_xtb
      self%element=aelement_nitroethane_xtb
      self%energy=energy_nitroethane_xtb
      self%xyz=sxyz_nitroethane_xtb
   case("nitromethane")
      self%id=sid_nitromethane_xtb
      self%area=sarea_nitromethane_xtb
      self%su=su_nitromethane_xtb
      self%atom_xyz=axyz_nitromethane_xtb
      self%element=aelement_nitromethane_xtb
      self%energy=energy_nitromethane_xtb
      self%xyz=sxyz_nitromethane_xtb
   case("nonane")
      self%id=sid_nonane_xtb
      self%area=sarea_nonane_xtb
      self%su=su_nonane_xtb
      self%atom_xyz=axyz_nonane_xtb
      self%element=aelement_nonane_xtb
      self%energy=energy_nonane_xtb
      self%xyz=sxyz_nonane_xtb
   case("nonanol")
      self%id=sid_nonanol_xtb
      self%area=sarea_nonanol_xtb
      self%su=su_nonanol_xtb
      self%atom_xyz=axyz_nonanol_xtb
      self%element=aelement_nonanol_xtb
      self%energy=energy_nonanol_xtb
      self%xyz=sxyz_nonanol_xtb
   case("octane")
      self%id=sid_octane_xtb
      self%area=sarea_octane_xtb
      self%su=su_octane_xtb
      self%atom_xyz=axyz_octane_xtb
      self%element=aelement_octane_xtb
      self%energy=energy_octane_xtb
      self%xyz=sxyz_octane_xtb
   case("octanol")
      self%id=sid_octanol_xtb
      self%area=sarea_octanol_xtb
      self%su=su_octanol_xtb
      self%atom_xyz=axyz_octanol_xtb
      self%element=aelement_octanol_xtb
      self%energy=energy_octanol_xtb
      self%xyz=sxyz_octanol_xtb
   case("odichlorobenzene")
      self%id=sid_odichlorobenzene_xtb
      self%area=sarea_odichlorobenzene_xtb
      self%su=su_odichlorobenzene_xtb
      self%atom_xyz=axyz_odichlorobenzene_xtb
      self%element=aelement_odichlorobenzene_xtb
      self%energy=energy_odichlorobenzene_xtb
      self%xyz=sxyz_odichlorobenzene_xtb
   case("onitrotoluene")
      self%id=sid_onitrotoluene_xtb
      self%area=sarea_onitrotoluene_xtb
      self%su=su_onitrotoluene_xtb
      self%atom_xyz=axyz_onitrotoluene_xtb
      self%element=aelement_onitrotoluene_xtb
      self%energy=energy_onitrotoluene_xtb
      self%xyz=sxyz_onitrotoluene_xtb
   case("pentadecane")
      self%id=sid_pentadecane_xtb
      self%area=sarea_pentadecane_xtb
      self%su=su_pentadecane_xtb
      self%atom_xyz=axyz_pentadecane_xtb
      self%element=aelement_pentadecane_xtb
      self%energy=energy_pentadecane_xtb
      self%xyz=sxyz_pentadecane_xtb
   case("pentane")
      self%id=sid_pentane_xtb
      self%area=sarea_pentane_xtb
      self%su=su_pentane_xtb
      self%atom_xyz=axyz_pentane_xtb
      self%element=aelement_pentane_xtb
      self%energy=energy_pentane_xtb
      self%xyz=sxyz_pentane_xtb
   case("pentanol")
      self%id=sid_pentanol_xtb
      self%area=sarea_pentanol_xtb
      self%su=su_pentanol_xtb
      self%atom_xyz=axyz_pentanol_xtb
      self%element=aelement_pentanol_xtb
      self%energy=energy_pentanol_xtb
      self%xyz=sxyz_pentanol_xtb
   case("perfluorobenzene")
      self%id=sid_perfluorobenzene_xtb
      self%area=sarea_perfluorobenzene_xtb
      self%su=su_perfluorobenzene_xtb
      self%atom_xyz=axyz_perfluorobenzene_xtb
      self%element=aelement_perfluorobenzene_xtb
      self%energy=energy_perfluorobenzene_xtb
      self%xyz=sxyz_perfluorobenzene_xtb
   case("phenylether")
      self%id=sid_phenylether_xtb
      self%area=sarea_phenylether_xtb
      self%su=su_phenylether_xtb
      self%atom_xyz=axyz_phenylether_xtb
      self%element=aelement_phenylether_xtb
      self%energy=energy_phenylether_xtb
      self%xyz=sxyz_phenylether_xtb
   case("propanol")
      self%id=sid_propanol_xtb
      self%area=sarea_propanol_xtb
      self%su=su_propanol_xtb
      self%atom_xyz=axyz_propanol_xtb
      self%element=aelement_propanol_xtb
      self%energy=energy_propanol_xtb
      self%xyz=sxyz_propanol_xtb
   case("pyridine")
      self%id=sid_pyridine_xtb
      self%area=sarea_pyridine_xtb
      self%su=su_pyridine_xtb
      self%atom_xyz=axyz_pyridine_xtb
      self%element=aelement_pyridine_xtb
      self%energy=energy_pyridine_xtb
      self%xyz=sxyz_pyridine_xtb
   case("secbutanol")
      self%id=sid_secbutanol_xtb
      self%area=sarea_secbutanol_xtb
      self%su=su_secbutanol_xtb
      self%atom_xyz=axyz_secbutanol_xtb
      self%element=aelement_secbutanol_xtb
      self%energy=energy_secbutanol_xtb
      self%xyz=sxyz_secbutanol_xtb
   case("secbutylbenzene")
      self%id=sid_secbutylbenzene_xtb
      self%area=sarea_secbutylbenzene_xtb
      self%su=su_secbutylbenzene_xtb
      self%atom_xyz=axyz_secbutylbenzene_xtb
      self%element=aelement_secbutylbenzene_xtb
      self%energy=energy_secbutylbenzene_xtb
      self%xyz=sxyz_secbutylbenzene_xtb
   case("tbutylbenzene")
      self%id=sid_tbutylbenzene_xtb
      self%area=sarea_tbutylbenzene_xtb
      self%su=su_tbutylbenzene_xtb
      self%atom_xyz=axyz_tbutylbenzene_xtb
      self%element=aelement_tbutylbenzene_xtb
      self%energy=energy_tbutylbenzene_xtb
      self%xyz=sxyz_tbutylbenzene_xtb
   case("tetrachloroethene")
      self%id=sid_tetrachloroethene_xtb
      self%area=sarea_tetrachloroethene_xtb
      self%su=su_tetrachloroethene_xtb
      self%atom_xyz=axyz_tetrachloroethene_xtb
      self%element=aelement_tetrachloroethene_xtb
      self%energy=energy_tetrachloroethene_xtb
      self%xyz=sxyz_tetrachloroethene_xtb
   case("tetrahydrofuran")
      self%id=sid_tetrahydrofuran_xtb
      self%area=sarea_tetrahydrofuran_xtb
      self%su=su_tetrahydrofuran_xtb
      self%atom_xyz=axyz_tetrahydrofuran_xtb
      self%element=aelement_tetrahydrofuran_xtb
      self%energy=energy_tetrahydrofuran_xtb
      self%xyz=sxyz_tetrahydrofuran_xtb
   case("tetrahydrothiophenedioxide")
      self%id=sid_tetrahydrothiophenedioxide_xtb
      self%area=sarea_tetrahydrothiophenedioxide_xtb
      self%su=su_tetrahydrothiophenedioxide_xtb
      self%atom_xyz=axyz_tetrahydrothiophenedioxide_xtb
      self%element=aelement_tetrahydrothiophenedioxide_xtb
      self%energy=energy_tetrahydrothiophenedioxide_xtb
      self%xyz=sxyz_tetrahydrothiophenedioxide_xtb
   case("tetralin")
      self%id=sid_tetralin_xtb
      self%area=sarea_tetralin_xtb
      self%su=su_tetralin_xtb
      self%atom_xyz=axyz_tetralin_xtb
      self%element=aelement_tetralin_xtb
      self%energy=energy_tetralin_xtb
      self%xyz=sxyz_tetralin_xtb
   case("toluene")
      self%id=sid_toluene_xtb
      self%area=sarea_toluene_xtb
      self%su=su_toluene_xtb
      self%atom_xyz=axyz_toluene_xtb
      self%element=aelement_toluene_xtb
      self%energy=energy_toluene_xtb
      self%xyz=sxyz_toluene_xtb
   case("tributylphosphate")
      self%id=sid_tributylphosphate_xtb
      self%area=sarea_tributylphosphate_xtb
      self%su=su_tributylphosphate_xtb
      self%atom_xyz=axyz_tributylphosphate_xtb
      self%element=aelement_tributylphosphate_xtb
      self%energy=energy_tributylphosphate_xtb
      self%xyz=sxyz_tributylphosphate_xtb
   case("triethylamine")
      self%id=sid_triethylamine_xtb
      self%area=sarea_triethylamine_xtb
      self%su=su_triethylamine_xtb
      self%atom_xyz=axyz_triethylamine_xtb
      self%element=aelement_triethylamine_xtb
      self%energy=energy_triethylamine_xtb
      self%xyz=sxyz_triethylamine_xtb
   case("trimethylbenzene")
      self%id=sid_trimethylbenzene_xtb
      self%area=sarea_trimethylbenzene_xtb
      self%su=su_trimethylbenzene_xtb
      self%atom_xyz=axyz_trimethylbenzene_xtb
      self%element=aelement_trimethylbenzene_xtb
      self%energy=energy_trimethylbenzene_xtb
      self%xyz=sxyz_trimethylbenzene_xtb
   case("undecane")
      self%id=sid_undecane_xtb
      self%area=sarea_undecane_xtb
      self%su=su_undecane_xtb
      self%atom_xyz=axyz_undecane_xtb
      self%element=aelement_undecane_xtb
      self%energy=energy_undecane_xtb
      self%xyz=sxyz_undecane_xtb
   case("water")
      self%id=sid_water_xtb
      self%area=sarea_water_xtb
      self%su=su_water_xtb
      self%atom_xyz=axyz_water_xtb
      self%element=aelement_water_xtb
      self%energy=energy_water_xtb
      self%xyz=sxyz_water_xtb 
   end select

end subroutine internalcosmo

end module internaldb
