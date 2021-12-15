# COSMO-X

## Introduction
This is an fully open source solvation model, based on the original 
conductor like screening model for realistic solvation (COSMO-RS) model by Klamt et al. in combination with the
universal solvation model based on solute electron density (SMD) by Marenich, Cramer and Truhlar.

While the parameters, that can be found in this repository are parametrized for r2scan-3c, the goal of this project is to allow reparametrization and
customization of this model for any functional or even semi-empirical methods.

## Building
You can use a release binary, but for compatibility, it is recommended to build this project by yourself from source. 
To do that, you need to clone the repository to your local environment.
```
git clone https://github.com/grimme-lab/COSMO-X.git
cd COSMO-X
```
You then have two options to build the project.

### Building with Fortran package Manager
You can use the Fortran Package Manager (https://github.com/fortran-lang/fpm) to build the project.
To install the project in your prefered path, just use 
```
fpm install -profile release -prefix [path]
```

More information about FPM can be found in the respective documentation.

### Building with Meson
You can also use Meson to build your project. To be able to do that, you need version 0.60 or newer.
To use the default backend, you also need ninja in version 1.7 or newer.
```
meson setup build --buildtype release
ninja -C build
```
You afterwards have to manually install the created binary in the path of your choice.

## Usage
You can either use a commandline version of COSMO-X or create a detailed input for your calculation.

### Using the commandline version (recommended)
To use the commandline version of COSMO-X, you have to set up an CSXHOME environment variable.
```
export CSXHOME=[path]
```
In the home directory, you need to place a configuration file, as well as the parameters and the COSMO database.
If you clone the repository, the commandline version should work out of the box, if you set the CSXHOME Variable to the cloned repository.
Otherwise, you can create a new sample configuration file and set up the home directory for yourself.
```
csx --newconfig
```
The commandline version needs a coord file for your solute, a working version of Turbomole, as well as an control file for the gas phase calculation.
You can then start a calculation (starting from the coord file) for example by
```
cefine -func r2scan-3c
csx --solvent water
```

### Using a detailed input
You can also set up an input file. COSMO-X needs cosmo files for the solute and the solvent. You can create a sample input by
```
csx --newinput
```
You then need to manually specify which parameters you want to use, as well as the respective cosmo files. You can start the calculation by
```
csx csx.input
```