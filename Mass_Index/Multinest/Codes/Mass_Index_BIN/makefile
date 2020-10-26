NESTLIBDIR = ../
#LIBS = -L$(NESTLIBDIR) -lnest3 -L/usr/lib/lapack -llapack
LIBS = -L$(NESTLIBDIR) -lnest3  -llapack -lblas

FFLAGS = -O3  
FC=gfortran
#  FC=mpif90  -O3 -DMPI

OBJFILES = params.o  transit.o like.o nestwrap.o main.o

OBJFILES_MOD = params.o  transit_MOD.o like.o nestwrap.o main_MOD.o

F90FILES_MOD =  params.f90  transit_MOD.f90  like.f90  nestwrap.f90   main_MOD.f90

F90FILES =  params.f90  transit.f90  like.f90  nestwrap.f90   main.f90 

all: pred LRVplan

mod: pred_MOD MODALL

pred:
	$(FC) $(FFLAGS) -I$(NESTLIBDIR) -c $(F90FILES)

pred_MOD:  
	$(FC) $(FFLAGS) -I$(NESTLIBDIR) -c $(F90FILES_MOD)
 
LRVplan: $(OBJFILES)
	$(FC) -o ./mass_index $(OBJFILES) \
	$(FFLAGS) -L/usr/lib/ $(LIBS) 

MODALL: $(OBJFILES)
	$(FC) -o ./mass_index_MOD $(OBJFILES_MOD) \
	$(FFLAGS) -L/usr/lib/ $(LIBS) 

clean:
	rm -f *.o *.mod mass_index
