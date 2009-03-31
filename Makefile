all:
	@echo 'Compiling and linking...'	
	@gfortran params.f90 Euler1Dlibs.f90 Euler1D.f90 tictoc.f90 -o Euler1D $(PGPLOT) -O3 -Wall
	@echo 'Done'
