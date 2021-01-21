# Comment
# Start of the makefile
test: main.o ElementParameters.o globals.o sort.o initialize.o sigma_av.o sac.o bonding.o
	gfortran -fbacktrace -o csm main.o ElementParameters.o globals.o sort.o initialize.o sigma_av.o sac.o bonding.o
main.o: main.f90 element_dict.mod globals.mod sort.mod sac.f90 initialize_cosmo.mod sigma_av.mod sac_mod.mod bonding.mod
	gfortran -c main.f90 
ElementParameters.o: dictionary.f90 linkedlist.f90 ElementParameters.f90
	gfortran -c ElementParameters.f90
element_dict.mod: dictionary.f90 linkedlist.f90 ElementParameters.f90
	gfortran -c ElementParameters.f90
globals.mod: globals.f90
	gfortran -c globals.f90
globals.o: globals.f90
	gfortran -c globals.f90
sort.o: sort.f90
	gfortran -c sort.f90
sort.mod: sort.f90
	gfortran -c sort.f90
initialize.o: initialize.f90 globals.mod
	gfortran -c initialize.f90
initialize_cosmo.mod: initialize.f90 globals.mod
	gfortran -c initialize.f90 
sigma_av.mod: sigma_av.f90
	gfortran -c sigma_av.f90
sigma_av.o : sigma_av.f90
	gfortran -c sigma_av.f90
sac.o : sac.f90
	gfortran -c sac.f90
sac_mod.mod : sac.f90
	gfortran -c sac.f90
bonding.mod : bonding.f90 globals.mod element_dict.mod
	gfortran -c bonding.f90
bonding.o : bonding.f90 globals.mod element_dict.mod
	gfortran -c bonding.f90

clean:
		rm main.o eledata_module.mod initialize_cosmo.mod initialize.o csm element_dict.mod globals.mod sort.mod ElementParameters.o globals.o sort.o sigma_av.o sigma_av.mod
# End of makefile
