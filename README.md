# COSMO-SE

Building with Fortran Package Manager:
https://github.com/fortran-lang/fpm

fpm build
fpm install prefix [path without bin]

Setting Home Directory Global Variable:
export CSMHOME=path/to/source/directory

Searches for .param Files and .cosmo files in home directory and working directory.

USAGE: ./csm --c solute --s solvent [--sigma --model (sac,crs,sac2010,sac2013) --T temperatur(default=298.15)]
