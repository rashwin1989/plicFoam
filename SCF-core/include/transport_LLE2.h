#ifndef         _TRANSPORT_LLE2_H_
#define         _TRANSPORT_LLE2_H_

#include "myUmfpack.h"

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
    double dL,  \
    double dR,  \
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

///////////////////////////////////////////////////////////////////
//  light oil-heavy oil phase equilibrium
//
// fixing x_1[j] (j!=k, jW) to
// obtain x_1 & x_2 using LLE flash calulations
void transport_LLE2_core(
    double P,   // system pressure (Pa)
    double T_s, // interface temperature (K)
    int n,      // number of species
    double *Pc, // critical pressure (Pa)
    double *Tc, // critical temperature (K)
    double *w,  // acentric factor
    int *type_k,// index for BIP
    double *x_1,// oil-rich phase
    double *x_2 // water-rich phase
    );

// evaluation function: sum_i(sum_{j!=i}(Ri-Rj)^2)
double transport_LLE2_eval_func(
    _PARAMETERS_);
//
///////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////
//  general oil-water phase behavior
//
void transport_LLE_core(
    double P,   // system pressure (Pa)
    double T_s, // interface temperature (K)
    int k,      // index heaviest HC
    int n,      // number of species
    double *Pc, // critical pressure (Pa)
    double *Tc, // critical temperature (K)
    double *w,  // acentric factor
    int *type_k,// index for BIP
    double *x_1,// oil-rich phase
    double *x_2 // water-rich phase
    );

// evaluation function: sum_i(sum_{j!=i}(Ri-Rj)^2)
double transport_LLE_eval_func(
    int k, // index of major oil species
    _PARAMETERS_);
//
///////////////////////////////////////////////////////////////////

// to evaluate the condition of the transport constraints
double transport_evaluation_func(int n, 
                                  double *y_1, double *y_2, 
                                  double *J_1, double *J_2);

#endif
