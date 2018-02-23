#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include "MACROS.H"
#include "PR_EoS.h"
#include "myUmfpack.h"
#include "vis_n_therm.h"
#include "transport_LLE2.h"
#include "gsl_optimizations2.h"
#include "Maxwell_Stefan_flux.h"

#include "data.h"

int main(int argc, char *argv[])
{

  #include "READ.H"
  #include "PRINT_INPUT.H"
  #include "INIT.H"


  if (n==2)
  {
    int k = 0;
    double f_min;

    findequilibrium_new2_(
	&P, // ! pressure (Unit: Pa)
	&T0, // ! temperature (Unit: K)
	&k, // ! the major oil species' index
	&n, // ! number of species
	Pc,// ! vector of critical pressures
	Tc,// ! vector of critical temperatures
	w, // ! vector of acentric factors
	tk, // ! vector of binary interaction types
	x_1, // ! vector of mass fractions: oil-rich phase
	x_2, // ! vector of mass fractions: water-rich phase
	&f_min);// ! minimum evaluation value
  }
  else if (n==3)
  {
  }

  #include "CLEAN.H"
  return 0;
}
