#ifndef       _GSL_OPTIMIZATIONS2_H_
#define       _GSL_OPTIMIZATIONS2_H_

#include <gsl/gsl_multimin.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct 
{
  int n;      /* number of species */ 
  double *Pc; /* vector of Pc (Pa) */
  double *Tc; /* vector of Tc (K)  */
  double *Vc; /* vector of Vc (cm3/mol) */
  double *w;  /* vector of acentric factor */
  double *MW; /* vector of MW (g/mol) */
  int *type_k;/* vector of PPR78 index */
  double *Tb; /* vector of normal boiling point (C) */
  double *SG; /* vector of SG */
  double *H8; /* vector of Cp_IG index, see vis_n_cond.f90 */

  int n_flux_type; /* flux type: 0 - ideal diffusion; 1 - non-ideal */
  double dL;       /* distance between interface and left side point */
  double dR;       /* distance between interface and right side point */
  double P;        /* pressure at interface (Pa) */
  double T_s;      /* temperature at interface (K) */

  double *x_L;     /* mole fraction at left  side point near the interface */
  double *x_R;     /* mole fraction at right side point near the interface */
  
  double *x_1;     /* mole fraction at left  side of the interface, denoted as 1 */
  double *x_2;     /* mole fraction at right side of the interface, denoted as 2 */

  double *lnphi_1; /* ln(phi) at 1 */
  double *lnphi_2; /* ln(phi) at 2 */
  double *Dij_1;   /* Dij at 1 (m^2/s) */
  double *Dij_2;   /* Dij at 2 (m^2/s) */
  double *H_1;     /* molar enthalpy at 1 (J/mol)*/
  double *H_2;     /* molar enthalpy at 2 (J/mol)*/
  double *flux_m_1;/* mass flux at 1 (kg/m^2-s)*/
  double *flux_m_2;/* mass flux at 2 (kg/m^2-s)*/
  LPT_UMFPACK flux_umf; /* pointer of an umfpack struct */
}
T_GSL_OPTIM,*LPT_GSL_OPTIM;

/*
  INPUT PARAMETERS defined in a macro called "_PARAMETERS_"
*/
    
#ifndef     _PARAMETERS_ 

#define     _PARAMETERS_  \
    int n,      \
    double *Pc, \
    double *Tc, \
    double *Vc, \
    double *w,  \
    double *MW, \
    int *type_k,\
    double *Tb, \
    double *SG, \
    double *H8, \
    \
    int n_flux_type, \
    double dL, \
    double dR, \
    double P,   \
    double T_s, \
    \
    double *x_L, \
    double *x_R, \
    double *x_1, \
    double *x_2, \
    double *lnphi_1, \
    double *lnphi_2, \
    double *Dij_1,   \
    double *Dij_2,   \
    double *H_1,     \
    double *H_2,     \
    double *flux_m_1,\
    double *flux_m_2,\
    LPT_UMFPACK flux_umf

#endif

int    my2_gsl_find_transport_LLE(_PARAMETERS_);

int    my2_gsl_find_transport_LLE2(_PARAMETERS_);


double transport_LLE_f_(const gsl_vector *v, void *params);

double transport_LLE2_f_(const gsl_vector *v, void *params);

#ifdef __cplusplus
}
#endif

#endif
