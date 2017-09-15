# COMPILER
#CC = gcc                 
#F77= gfortran            
#FL = $(CC)#Linker        
#NOFOR_MAIN =
#CFLAGS_ADDON = 
#FFLAGS_ADDON = 
#LIB_ADDON =

#  Intel C compiler
CC = icc
F77= ifort
FL = $(F77)#Linker 
NOFOR_MAIN = -nofor_main
CFLAGS_ADDON = -static-intel #-Wunused-variable
FFLAGS_ADDON = -static-intel #-warn unused
LIB_ADDON = -Wl,-Bstatic 

#  mixed
#CC = gcc
#F77= ifort
#FL = $(F77)#Linker 
#NOFOR_MAIN = -nofor_main
#CFLAGS_ADDON = -static
#FFLAGS_ADDON = -static-intel
#LIB_ADDON = -Wl,-Bstatic 

default: all
all: testSCF

r: c all #cp

#-DCUTOFF_MAX=1e6

CFLAGS = -O3 $(CFLAGS_ADDON)
FFLAGS = -O3 $(FFLAGS_ADDON) -fPIC 

LIB= -lm -lgfortran $(LIB_ADDON) \
     -L../lib -lscfcore \
     -lgsl -lgslcblas \
     /home/phe/software/UMFPACK/Lib/libumfpack.a \
     /home/phe/software/AMD/Lib/libamd.a \
     /home/phe/software/UFconfig/xerbla/libcerbla.a  
			 
CPPFLAGS = -I./ \
	   -I/usr/include/x86_64-linux-gnu/ \
	   -I/home/phe/software/UMFPACK/Include \
	   -I/home/phe/software/AMD/Include \
	   -I/home/phe/software/UFconfig \
	   -I../include


SRC = main.o data.o

HEADER = CLEAN.H  INIT.H  PRINT_INIT.H  PRINT_INPUT.H  READ.H


testSCF: $(SRC) 
	$(FL) $(NOFOR_MAIN) $(CFLAGS) -o $@ $(SRC) $(LIB)

main.o: $(HEADER)

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

