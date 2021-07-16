# COSMO-SE

Building with Fortran Package Manager:
https://github.com/fortran-lang/fpm

fpm build &&
fpm install prefix [path without bin]

USES INPUT File, use with
./csm file.inp

Example Input File:

/path/to/sac/parameters/sac.param 
/path/to/SMD_Parameters/ 
sac #Keywords
#Comment Line
/path/to/solvent/cosmo/file.cosmo
/path/to/solute/cosmo/file.cosmo
smd_solvent smd_probe_radius
Temperature
0.995 0.005
