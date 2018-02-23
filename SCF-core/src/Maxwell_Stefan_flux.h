#ifndef _MAXWELL_STEFAN_FLUX_H_
#define _MAXWELL_STEFAN_FLUX_H_

#include <stdio.h>
#include "myUmfpack.h"

#ifdef __cplusplus
extern "C" {
#endif

/* calc_MS_flux_interface_combination() 
 * computes the mass fluxes on both sides of the interface
 */
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
    );

#ifdef __cplusplus
}
#endif

#endif
