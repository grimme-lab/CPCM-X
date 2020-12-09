# Comment
# Start of the makefile
test: main.o ElementParameters.o globals.o sort.o
	gfortran -o test main.o ElementParameters.o globals.o sort.o
main.o: main.f90 element_dict.mod globals.mod sort.mod
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

clean:
		rm main.o test element_dict.o element_dict.mod
# End of makefile
