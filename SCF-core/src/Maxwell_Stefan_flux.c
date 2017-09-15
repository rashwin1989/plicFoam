#include <math.h>
#include "MACROS.H"
#include "Maxwell_Stefan_flux.h"

void calc_MS_flux_interface(
    int n_flux_type,
    double dL, double dR, 
    double P,
    double T_s, 
    int n,
    double *MW,
    double *x_1, double *x_L, 
    double *x_2, double *x_R, 
    double *dlnphi_1,double *dlnphi_2,
    double *Dij_1,   double *Dij_2,
    double  V_1,     double  V_2,
    double *rhs_flux,
    double *flux_m_1,
    double *flux_m_2,
    LPT_UMFPACK flux_umf
    )
{
  int j, j1;
  double R_gas = 8.3144621;
  double Ta, Va, Z;

  static int n1st=1;
  static double *dx_tmp; 

  if (n1st) {
    n1st=0;
    _NEW_(dx_tmp,double,n);
  }

  Ta = T_s;

  // ---------------------------------------------
  // LEFT HAND SIDE OF THE INTERFACE
  Va = V_1;
  Z = P*Va/(R_gas*Ta);

  for (j=0; j<n; j++) 
  {
    dx_tmp[j] = (x_1[j]-x_L[j])/dL;
  }

  for (j=0; j<n; j++) 
  {
    rhs_flux[j] = - dx_tmp[j];
    if (n_flux_type>=0)
        for (j1=0;j1<n;j1++) 
	{
            int idx = j + j1*n;
            if (x_1[j]<1e-8) continue;

            rhs_flux[j] += - x_1[j]*dlnphi_1[idx]*dx_tmp[j1];
        }

    rhs_flux[j] *= 1./Va;         // 1/V[i] = c[i] mol/m^3
  }

  Maxwell_Stefan_mass_flux(Z,n,MW,x_1,Dij_1,rhs_flux,flux_m_1,flux_umf);

  // ---------------------------------------------
  // RIGHT HAND SIDE OF THE INTERFACE
  Va = V_2;
  Z = P*Va/(R_gas*Ta);

  for (j=0; j<n; j++) 
  {
    dx_tmp[j] = (x_R[j]-x_2[j])/dR;
  }

  for (j=0; j<n; j++) 
  {
    rhs_flux[j] = - dx_tmp[j];
    if (n_flux_type>=0)
        for (j1=0;j1<n;j1++) 
	{
            int idx = j + j1*n;
            if (x_2[j]<1e-8) continue;

            rhs_flux[j] += - x_2[j]*dlnphi_2[idx]*dx_tmp[j1];
        }

    rhs_flux[j] *= 1./Va;         // 1/V[i] = c[i] mol/m^3
  }

  // molar flux
  Maxwell_Stefan_mass_flux(Z,n,MW,x_2,Dij_2,rhs_flux,flux_m_2,flux_umf);
}

