# COMPILER
CC = gcc                 
F77= gfortran            
FL = $(CC)#Linker        
#NOFOR_MAIN =
#CFLAGS_ADDON = 
#FFLAGS_ADDON = 
#LIB_ADDON =

#  Intel C compiler
#CC = icc
#F77= ifort
#FL = $(F77)#Linker 
#NOFOR_MAIN = -nofor_main
#CFLAGS_ADDON = -static-intel -Wunused-variable
#FFLAGS_ADDON = -static-intel #-warn unused
#LIB_ADDON = -Wl,-Bstatic 

#  mixed
#CC = gcc
#F77= ifort
#FL = $(F77)#Linker 
#NOFOR_MAIN = -nofor_main
#CFLAGS_ADDON = -static
#FFLAGS_ADDON = -static-intel
#LIB_ADDON = -Wl,-Bstatic 

default: all
all: SCF

r: c all #cp

#-DCUTOFF_MAX=1e6

CFLAGS = -O3 $(CFLAGS_ADDON)
FFLAGS = -O3 $(FFLAGS_ADDON) -fPIC 

LIB= -lm -lgfortran $(LIB_ADDON) \
     -lgsl -lgslcblas \
     ../../UMFPACK/Lib/libumfpack.a \
     ../../AMD/Lib/libamd.a \
     ../../UFconfig/xerbla/libcerbla.a

CPPFLAGS = -I./ \
	   -I/usr/include/x86_64-linux-gnu/ \
	   -I../../UMFPACK/Include \
	   -I../../AMD/Include \
	   -I../../UFconfig 


SRC = PR_EoS.o vis_n_therm.o Maxwell_Stefan_flux.o myUmfpack.o transport_LLE2.o gsl_optimizations2.o 

HEADER = 


SCF:	$(SRC) 
	'cp' *.[Hh] ../include/
	find . -name "*.o" > objs_file.tmp
	for FILE in `cat objs_file.tmp`; do ar -ru ../lib/libscfcore.a $$FILE; done
	ranlib ../lib/libscfcore.a
	rm -f objs_file.tmp 

%.o: %.f
	$(F77) -c $(FFLAGS) -o $@ $<
%.o: %.f90
	$(F77) -c $(FFLAGS) -o $@ $<
%.o: %.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) -o $@ $<

#-------------------------------------------------------------------------------
# Remove all but the files in the original distribution
#-------------------------------------------------------------------------------

purge: 
	- $(RM) *.o

c: purge
cm: c default

