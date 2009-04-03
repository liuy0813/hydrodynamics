PROGRAM  = Euler2D
COMPILER = gfortran
FLAGS    = -O3 -Wall
PGPLOT   = -lX11 -lpgplot

SOURCES = ./params.f90    \
	  ./Euler2D.f90   \
	  ./initconds.f90 \
	  ./output.f90    \
	  ./vars.f90      \
	  ./tictoc.f90

$(PROGRAM) : ${SOURCES}
	echo 'Compiling and linking...'	
	$(COMPILER) $(SOURCES) $(LIBS) -o $(PROGRAM) $(PGPLOT) $(FLAGS)
	echo 'Done'
