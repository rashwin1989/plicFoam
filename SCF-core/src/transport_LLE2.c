#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "MACROS.H"
#include "PR_EoS.h"
#include "vis_n_therm.h"
#include "transport_LLE2.h"
#include "Maxwell_Stefan_flux.h"

void x2y(int n, double *MW, double *x, double *y){
    int i;
    double sum_MWx = 0.0, r1_sum;
    for (i=0; i<n; i++) sum_MWx += MW[i]*x[i];

    r1_sum = 1.0/sum_MWx;

    for (i=0; i<n; i++) y[i] = MW[i]*x[i] * r1_sum;    
}

void y2x(int n, double *MW, double *x, double *y){
    int i;
    double sum_y_MW = 0.0;
    for (i=0; i<n; i++) sum_y_MW += y[i]/MW[i];

    sum_y_MW = 1.0/sum_y_MW;

    for (i=0; i<n; i++) x[i] = y[i]/MW[i] * sum_y_MW;
}

void transport_LLE2_core(
    double P,   // system pressure (Pa)
    double T_s, // interface temperature (K)
    int n,      // number of species
    double *Pc, // critical pressure (Pa)
    double *Tc, // critical temperature (K)
    double *w,  // acentric factor
    double *kij,// matrix of BIPs
    double *X_1,// oil-rich phase
    double *X_2 // water-rich phase
    )
{  
  int bprint=0, binit=0;
  double s_min;  

  // fix x_1[i=3,...,n-2] (heavier oils in phase 1) and x_2[n-1] (water in phase 2)
  findequilibrium_fix_water_b_(&bprint,&binit,
      &P,&T_s,&n,Pc,Tc,w,kij,X_1,X_2,&s_min);
}


// evaluation function: sum_i(sum_{j!=i}(Ri-Rj)^2)
double transport_LLE2_eval_func(
    _PARAMETERS_
    )
{
  int j;
  double _tmp;

  static int n1st=1; 
  static double *rhs_flux, *y_1, *y_2, *dlnphi_1, *dlnphi_2;
  double V_1, V_2;

  if (n1st){
    n1st=0;
    _NEW_(rhs_flux,double,n);
    _NEW_(y_1,double,n);
    _NEW_(y_2,double,n);
    _NEW_(dlnphi_1,double,n*n);
    _NEW_(dlnphi_2,double,n*n);
  }

  // -1) latest, setup to search in positive domain
  for (j=0; j<n; j++) {
    if (x_1[j]<0) x_1[j]=-x_1[j];
    if (x_2[j]<0) x_2[j]=-x_2[j];
  }

  // LLE
  transport_LLE2_core(P,T_s,n,Pc,Tc,w,kij,x_1,x_2);

  // properties
  // phase-1  
  fugacities_n_its_derivatives3_(&P,&T_s,&n,Pc,Tc,w,x_1,kij,lnphi_1,dlnphi_1,&V_1);
  new_tlsm_diffusion_krishna_model_(&P,&T_s,&n,Pc,Tc,Vc,w,MW,kij,x_1,Dij_1);

  // phase-2  
  fugacities_n_its_derivatives3_(&P,&T_s,&n,Pc,Tc,w,x_2,kij,lnphi_2,dlnphi_2,&V_2);
  new_tlsm_diffusion_krishna_model_(&P,&T_s,&n,Pc,Tc,Vc,w,MW,kij,x_2,Dij_2);

  x2y(n, MW, x_1, y_1);
  x2y(n, MW, x_2, y_2);

  // MS Flux
  calc_MS_flux_interface(n_flux_type,dL,dR,
      P,T_s, n,MW, x_1,x_L,x_2,x_R,
      dlnphi_1,dlnphi_2,Dij_1,Dij_2,V_1,V_2,
      rhs_flux,flux_m_1,flux_m_2,flux_umf);

  // if the fugacity criteria is NOT fulfilled
  // return -3 as an error code
  {
    // criteria: x>=0
    for(j=0;j<n;j++){
      double xm = x_1[j];
      if (xm<0) {
        if (x_2[j]<xm) xm = x_2[j];
      }else{
        if (x_2[j]<0) xm = x_2[j];
      }

      if (xm<0) return x_2[n-1]*1e8;
    }

    // criteria: x<=1
    for(j=0;j<n;j++){
      double xm = x_1[j];
      if (xm>1) {
        if (x_2[j]>xm) xm = x_2[j];
      }else{
        if (x_2[j]>1) xm = x_2[j];
      }

      if (xm>1) return x_2[n-1]*1e8;
    }

  }

  _tmp = transport_evaluation_func(n,y_1,y_2,flux_m_1,flux_m_2);

  return _tmp;
}

double transport_evaluation_func(int n, 
                                  double *y_1, double *y_2, 
                                  double *J_1, double *J_2)
{
  int j;
  double _tmp, val, max;

  max = 0;
  for (j=0; j<n; j++) {
    _tmp = fabs(J_1[j]); if (max<_tmp) max = _tmp;
    _tmp = fabs(J_2[j]); if (max<_tmp) max = _tmp;
  }

  //-------------------------------------------------------------------
  // new eval func: when y_1 - y_2 ~ 0, this definition has singularity
  //-------------------------------------------------------------------
  val = 0;
  for (j=0; j<n; j++) 
  {
      int j1;
    
      for (j1=j+1; j1<n; j1++)
      {

        _tmp = (J_1[j1]-J_2[j1])*(y_1[j ]-y_2[j ]) 
              -(J_1[j ]-J_2[j ])*(y_1[j1]-y_2[j1]);
        val += _tmp * _tmp;

      }
    
  }

  //printf("sum:%le\n", val);
  //------------------------------------------------------------------

  return val;
}

void transport_LLE_core(
    double P,   // system pressure (Pa)
    double T_s, // interface temperature (K)
    int k,      // index of major oil species
    int n,      // number of species
    double *Pc, // critical pressure (Pa)
    double *Tc, // critical temperature (K)
    double *w,  // acentric factor
    double *kij,// matrix of BIPs
    double *x_1,// oil-rich phase
    double *x_2 // water-rich phase
    )
{
  int k_fortran;
  double s_min;

  k_fortran=k+1;

  printf("-------------------------------------------------------\n");
  printf("findequilibrium_new2_ FORTRAN function\n");
  printf("-------------------------------------------------------\n");
  printf("\n");

  findequilibrium_new2_(&P,&T_s,&k_fortran,&n,Pc,Tc,w,kij,x_1,x_2,&s_min);
}

// evaluation function: sum_i(sum_{j!=i}(Ri-Rj)^2)
double transport_LLE_eval_func(
    int k, // index of major oil species
    _PARAMETERS_
    )
{
  int j;
  double _tmp;

  static int n1st=1; 
  static double *rhs_flux, *y_1, *y_2, *dlnphi_1, *dlnphi_2;
  double V_1, V_2;

  if (n1st){
    n1st=0;
    _NEW_(rhs_flux,double,n);
    _NEW_(dlnphi_1,double,n*n);
    _NEW_(dlnphi_2,double,n*n);
    _NEW_(y_1,double,n);
    _NEW_(y_2,double,n);    
  }

  printf("-------------------------------------------------------\n");
  printf("transport_LLE_eval_func\n");
  printf("-------------------------------------------------------\n");
  printf("\n");

  printf("x_1: "); for(j=0; j<n; j++) printf("%.9f  ", x_1[j]);
  printf("\n");
  printf("x_2: "); for(j=0; j<n; j++) printf("%.9f  ", x_2[j]);
  printf("\n");

  // -1) latest, setup to search in positive domain
  for (j=0; j<n; j++) {
    if (x_1[j]<0) x_1[j]=-x_1[j];
    if (x_2[j]<0) x_2[j]=-x_2[j];
  }

  // LLE
  transport_LLE_core(P,T_s,k,n,Pc,Tc,w,kij,x_1,x_2);

  // properties
  // phase-1  
  fugacities_n_its_derivatives3_(&P,&T_s,&n,Pc,Tc,w,x_1,kij,lnphi_1,dlnphi_1,&V_1);
  new_tlsm_diffusion_krishna_model_(&P,&T_s,&n,Pc,Tc,Vc,w,MW,kij,x_1,Dij_1);

  // phase-2  
  fugacities_n_its_derivatives3_(&P,&T_s,&n,Pc,Tc,w,x_2,kij,lnphi_2,dlnphi_2,&V_2);
  new_tlsm_diffusion_krishna_model_(&P,&T_s,&n,Pc,Tc,Vc,w,MW,kij,x_2,Dij_2);

  x2y(n, MW, x_1, y_1);
  x2y(n, MW, x_2, y_2);

  // if the fugacity criteria is NOT fulfilled
  // return -3 as an error code
  {
    // criteria: x>=0
    for(j=0;j<n;j++){
      double xm = x_1[j];
      if (xm<0) {
        if (x_2[j]<xm) xm = x_2[j];
      }else{
        if (x_2[j]<0) xm = x_2[j];
      }

      if (xm<0) return -xm*1e5;
    }

    // criteria: x<=1
    for(j=0;j<n;j++){
      double xm = x_1[j];
      if (xm>1) {
        if (x_2[j]>xm) xm = x_2[j];
      }else{
        if (x_2[j]>1) xm = x_2[j];
      }

      if (xm>1) return (xm-1)*1e5;
    }
  }

  // MS Flux
  calc_MS_flux_interface(n_flux_type,dL,dR,
      P,T_s,
      n,MW, x_1,x_L,x_2,x_R,
      dlnphi_1,dlnphi_2,Dij_1,Dij_2,V_1,V_2,
      rhs_flux,flux_m_1,flux_m_2,flux_umf);

  _tmp = transport_evaluation_func(n,y_1,y_2,flux_m_1,flux_m_2);

  return _tmp;
}
