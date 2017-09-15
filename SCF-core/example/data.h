#include <stdio.h>
#include "myUmfpack.h"

//-------------------------------------------------------------------------------
// Data for species
//-------------------------------------------------------------------------------
extern
int     n,             // num of species
       *idx_species,   // idx of species according to database petro.dat+
       *idx_groups,    // idx of PPR78 groups according to database groups.dat+
       *tk;            // type of BIP
extern
double *Pc, *Tc, *Vc,  // critical properties (Pa, K, cm^3/mol)
       *w,             // acentric factor
       *MW,            // molecular weight (g/mol)
       *dm,            // dipole moment (Debye)
       *k,             // association parameter
       *H8,            // molar enthalpy @ T=.8Tc (Btu/lb), which can be 0
       *Tb,            // boiling point (C)
       *SG;            // specific gravity
//-------------------------------------------------------------------------------

//-------------------------------------------------------------------------------
// Data for grid points: 
//          grid points are marked using 'x'
//          'o' is not used explicitly
//-------------------------------------------------------------------------------
//         r=0                         r=R
//          |                           |
//       x  o  x  o  x ....... x  o  x  o  x 
//       |     |     |         |     |     |
// lo_ghost=0  1     2       NR-3   NR-2  hi_ghost=NR-1
//            i_lo                  i_hi
//-------------------------------------------------------------------------------
extern
int NR,                         // number of total grid points on r-axis
    i_lo, i_lo_ghost,           // index of r_min (r_min=0)
    i_hi, i_hi_ghost;           // index of r_max
extern
double dr,                      // delta r
       R,                       // maximum radial distance
       *r;                      // radial locations r[NR] for 'x' points

// x[NR][n]
extern int *phase_type;         // 0: gas; 1: liquid; -1: only phase
extern
double **x, **y, **x_n, **y_n,  // molar & mass fractions (current and last time)
       **rhoY,   **rhoY_n,      // rho*y
       // fluxes are calculated at the center of a cell, i.e. in the mid of 2 grid points
       // flux_m[i][n] is between i and i+1, i has range 0 - NR-1
       **flux_m,                // mass flux current from dx/dr
       **flux_m_,               // mass flux from Maxwell
        *flux_h,                // heat flux: -k*dT/dr 
       **rhs_flux,              // rhs of Maxwell-Stefan equation
       **lnphi,                 // ln(phi)[NR grid][n species]
       **phi,                   // phi    [NR grid][n species]
       **xphi,                  // x*phi  [NR grid][n species]
       **advect,                // rho*y*u*r^2
       **dlnphidT,              // d(lnphi)/dT [NR grid][n species]

       **D,                     // mass diffusivity (m^2/s)
       **Cp_IG,                 // ideal gas Cp
       **Hpar,                  // partial molar enthalpy
       **Hdep,                  // molar enthalpy departure
       // Dij[NR][n*n]
       **Dij,                   // mass diffusivity (m^2/s) index = i+j*n FORTRAN FORMAT
       **Dij_mc,                // multicomponent Dij using MTGL
       **Dij_Fick,              // multicomponent Dij in Fick's law: J=-c*D*grad(x)
       **Dij_new,               // Wesselingh & Krishna Model
       **dxj_dr,                // d(x_j)/dr [NR][n]
       **dxj_dr_,               // d(x_j)/dr [NR][n]
       **dlnphi_dxj,            // d(lnphi_i)/dx_j [NR][n*n]
       *dlnphi_dxj_tmp,         // d(lnphi_i)/dx_j [n*n]
       *dlnphidT_tmp,           // d(lnphi_i)/dT [n]
       // rho[NR]
       *rho, *rho_n, *rho_s,    // density (kg/m^3): rho_s is directly from mass transfer eq.
       *u,                      // velocity (m/s)
       *u_f,                    // (u+u_n)/2
       *T,                      // temperature (K)
       *p,                      // pressure (Pa)
       *cond,                   // thermal conductivity (W/m/K)
       *V,                      // molar volume (m^3/mol)
       *MWave,                  // average MW (g/mol)
       *Cp,                     // heat capacity (J/mol/K)
       *Cv,                     // heat capacity (J/mol/K)
       *Cv_IG,                  // ideal gas Cv
       *vis,                    // viscosity Pa*sec
       *u_n, *T_n,              // *_n for storage of last time step
       *am, *bm,                // PR-EOS paramters

       *drhodT,                 // drho/dT
      **drhodY,                 // drho/dy[j]

      **S_IG,                   // ideal gas entropy
      **H_IG,                   // ideal gas enthalpy
      **G_IG,                   // ideal gas Gibbs energy
      **dG_IG,                  // ideal gas Gibbs energy
      **G_par,                  // partial Gibbs energy, i.e. chemical potential

// NOT directly used DATA; BUT requied by FORTRAN CODE
// Vpar[n] 
       *Vpar, *lnphi_tmp, dVdT, G, *coef_ab, *x_tmp, *dx_tmp,
       *Dij_tmp, *kij_tmp; // [n*n] 
extern int G_only;
extern int bKijSet;    // to set kij's state
//-------------------------------------------------------------------------------

//-------------------------------------------------------------------------------
// Data for crude - oil interface
//-------------------------------------------------------------------------------
extern int    i_s, i_s_n;              // index of interface location
extern int phase_type_1, phase_type_2;
extern
double r_s,                     // interface location
       r_s_n,                   //   ..        ..      of LAST TIME STEP
       u_s,                     // interface velocity
       u_s_n, u_s_f,            // u_s_f=(u_s+u_s_n)/2
      *x_1, *xi1,               // species fraction on left (xi1: initial guess)
      *x_2, *xi2,               // species fraction on right
      *x_a, *y_a,                // average x and y for face point
      *x_1_n,*x_2_n,
      *y_1, *rhoY_1,
      *y_2, *rhoY_2,
      *flux_m_1, *flux_m_2, // one-side difference for transport LLE
       flux_h_1,  flux_h_2, // heat flux at interface
      *flux_m_A, *flux_m_B, // central difference for mass transfer
      *dydr_1,   *dydr_2,
       MWave_1, MWave_2,
       T_1, T_2, u_1, u_2,
       cond_1, cond_2,
       Cp_1, Cp_2,
      *D_1, *D_2,
      *Dij_1, *Dij_2,           // mass diffusivity (m^2/s) index = i+j*n FORTRAN FORMAT
      *Dij_mc_1, *Dij_mc_2,     // multicomponent Dij
      *h_1, *h_2,
      *Hdep_1,*Hdep_2,
       rho_1,   rho_2,
       rho_1_n, rho_2_n,
       rho1,    rho2,
      *lnphi_1,      *lnphi_2,       // length n
      *phi_1,        *phi_2,         // ln(x*phi)
      *xphi_1,       *xphi_2,        // ln(x*phi)
      *K, // fugacity ratio
      *RM_ni, *M_ai_yi, // coefs for transport relations (see TRANSPROT_LLE.H)
      *x_1_max,         // maximum fractions for each binary set (oil_i + water)
      *x_1_min,         // minimum fractions for each binary set (oil_i + water)
      *x_1_dx,          // delta molar fraction for search
      *dlnphi_dxj_1, *dlnphi_dxj_2,  // length nxn, for d(lnf_i)/dx_j : idx = i + j*n
      *dfdx,            // [n] for derivatives of evaluation func of transport constraints
       am_1, am_2,   
       bm_1, bm_2,
       V_1, V_2,
       vis_1, vis_2,
       V_tmp,Cp_tmp,Cv_tmp,*Cp_IG_tmp,Cv_IG_tmp,*Vpar_tmp,
       T_s,                     // inteface temperature
       dL,                      // dist from interface to right grid point
       dR;                      // dist from interface to left grid point
extern int    n_miscible;              // 1: YES; 0: NO
extern double du1_s;
//-------------------------------------------------------------------------------

// chemical potentials: 
extern double **mu_IG, *mu_IG_1, *mu_IG_2, 
              **dMu,   *dMu_1,   *dMu_2;

// data for upwind HJ-WENO5 using ghost points
/*
// 1. indices
int i_lo_gb, // = 0        // fixed
    i_hi_gb, // = NR  + 6  // fixed
    i_lo_fb_L, // = 4        // fixed
    i_hi_fb_L, // = i_s + 3           // moving with interface
    i_lo_fb_R, // = i_s + 4           // moving with interface
    i_hi_fb_R; // = NR  + 2  // fixed
// 2. data
double **x_L,     **x_R,     // [n][3+NR+4=NR+7]
       **lnphi_L, **lnphi_R, //   the same as above
       **flux_L,  **flux_R,  //   ...
       **rhoyu_L, **rhoyu_R, //   ...
        *T_L,      *T_R,     // [NR+7]
        *D1;                 // [NR+7] temporary vector
double **dxdr_L,     **dxdr_R,     
       **dlnphidr_L, **dlnphidr_R, 
       **drhoyudr_L, **drhoyudr_R, 
        *dTdr_L,      *dTdr_R;     
        */

// Runge-Kutta data
extern double dt_org,
       **RHS_Y_RK1, *RHS_T_RK1,
       **RHS_Y_RK2, *RHS_T_RK2,
       **RHS_Y_RK3, *RHS_T_RK3,
       **RHS_Y_RK4, *RHS_T_RK4,
       **RHS_Y[4],  *RHS_T[4],  RHS_rs[4];
extern int iRK;

extern double mass_rate;  // mass flux rate at the interface

// data for implicit interface
#ifndef _IMPLICIT_INTERFACE_
#define _IMPLICIT_INTERFACE_
typedef struct implicit_interface{
  int found;
  int n_occur; // occurrence times
  int i_s, i_s_n;
  double jump;  // jump evaluation, see utility_func.c
  double r_s, r_s_n;
  double u_s;
  struct implicit_interface * next;
}T_IMP_INTERFACE,*LPT_IMP_INTERFACE;
#endif

extern int n_imp_interface;
extern LPT_IMP_INTERFACE imp_interface;

//-------------------------------------------------------------------------------
// Data for initial conditions & boundary conditions
//-------------------------------------------------------------------------------
extern
double T0, // @ r=0, t=0
       T1, // @ r=R, t=any time
       P,  // system pressure (Pa)
       *x1, *x2;   // molar fractions of pure oil-phase and pure water phase
//-------------------------------------------------------------------------------

//-------------------------------------------------------------------------------
// Data for time marching
//-------------------------------------------------------------------------------
extern
double t, t_max, dt,
       max_u;
//-------------------------------------------------------------------------------

//-------------------------------------------------------------------------------
// UMFPACK Data for Maxwell-Stefan mass flux solver
//-------------------------------------------------------------------------------
extern T_UMFPACK flux_umf;
//-------------------------------------------------------------------------------

//-------------------------------------------------------------------------------
// Control variables
//-------------------------------------------------------------------------------
extern int restart;  // 1: YES;   0: NO, initialize system with 0
// if restart=0: u=0, [xW=0, T=T0](r<=r_s); [xW=1, T=T1](r>r_s);
// if restart=1: read files to initialize u & T, x
extern int nfreq;    // output frequency
extern double Courant; // Courant Number
extern int bConvection; // 1: YES; 0: NO convection in Mass Transfer and Energy Equations
extern int n_flux_type;   // 1: Maxwell-Stefan flux; 0: Fick's law
//-------------------------------------------------------------------------------

// Multiple interfaces
// 1. struct
// moved to interface.h
//#ifndef _T_INTERFACE_
//#define _T_INTERFACE_
//#endif

// 2. Number of interfaces and their data structure linked list
//extern int nInterfaces;
//extern LPT_INTERFACE interfaces, ifc_curr;

// a buffer distance only for initial temperature profile
// within this distance from the droplet, the lower temperature will
// be set to prevent the crash of the simulation
extern double r0_buffer;

// number of expected interfaces used in the output data
extern int nExpectInterfaces;

// Fully miscible conditions
extern int nFullyMiscible;

// using constant delta time
extern int bConstantDt;

// maximum GSL loop number
extern int n_GSL_max_loop;

extern int nDirichlet_U1;

// output file
extern FILE * f_output;

// to check conservation of species and enthalpy
// 1 - 1st checkup point after species & energy conserv. eqs.
// 2 - 2nd checkup point using control volume analysis
extern int    checkpoint; // = 0 - none; 1 - first; 2 - second
extern double **n1,       // number of moles in each cell
              *n1s, *n2s, // number of moles in the entire domain
	      *n3s, *n4s, // check points # 3, 4
	       *H1,       // enthalpy in each cell
	       H1s,  H2s, // enthalpy in the entire domain
	       M1s,  M2s; // sum of total mass
extern double total_enthalpy_mixing;

extern double D_interface; // interfacial mass diffusivity [m^2/s]
