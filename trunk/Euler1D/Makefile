PROGRAM  = Euler1D
SOURCES  = params.f90        \
	   Euler1Dlibs.f90   \
	   Euler1D.f90       \
           tictoc.f90
COMPILER =gfortran
FLAGS    =-O3 -Wall
PGPLOT   =-lX11 -lpgplot

$(PROGRAM): $(SOURCES)
	@echo 'Compiling and linking...'	
	$(COMPILER) $(SOURCES) -o $(PROGRAM) $(FLAGS) $(PGPLOT)
	@echo 'Done'
