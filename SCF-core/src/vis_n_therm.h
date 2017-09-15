#ifndef                 _PHE_VIS_N_THERM_H_
#define                 _PHE_VIS_N_THERM_H_

void vis_n_cond_(double *P,       // P (Pa)  scalar
                 double *T,       // T (K)   scalar
                 int    *n,       // n       scalar
                 double *Pc,      // Pc (Pa) vector
                 double *Tc,      // Tc (K)  vector
                 double *Vc,      // Vc (cm^3/mol)   vector
                 double *w,       // acentric factor vector
                 double *MW,      // MW (g/mol)      vector
                 double *k,       // association f.  vector
                 double *dm,      // dipole m.  (D)  vector
                 double *x,       // molar frac.     vector
                 double *Cv_m,    // Cv of mixture   scalar    (J/mol/K)
                 double *V,       // molar volume    scalar    (m^3/mol)
                 double *cond_m,  // thermal cond.   scalar    OUTPUT  (W/m/K)
                 double *vis_m);  // viscosity       scalar    OUTPUT  (Pa*s)

void thermo_properties_(double *P, double *T,                   // INPUT
    int *n,double *Pc,double *Tc,double *w,double *MW,double *x,// INPUT
    double *Tb,double *SG,double *H8,int *type_k,               // INPUT
    double *V,double *CP,double *CV,       // OUTPUT scalar 
    double *CP_IG,                         // OUTPUT vector
    double *Hpar,double *Hpur,double *Vpar,// OUTPUT vector
    double *dVdT,                          // OUTPUT scalar: dV/dT
    double *G,                             // OUTPUT scalar: excess Gibbs energy
    double *lnphi,                         // OUTPUT vector: ln(fuga_coeff)
    double *am, double *bm,                // OUTPUT vector: PR_EOS parameter a & b
    int *G_only);                          // INPUT: 1-onlyGibbs; 0-everything

void tlsm_diffusion_ij_( 
    double *rho,             // density      scalar
    double *T,               // temperature  scalar
    int    *n,               // # of species scalar
    double *Pc,              // vector of Pc (Pa)
    double *Tc,              // vector of Tc (K)
    double *Vc,              // vector of Vc (cm^3/mol)
    double *MW,              // vector of molecular weights
    double *x,               // vector of molar fractions
    double *D,               // vector of mass diffusivity
    double *Dij);            // matrix of binary mass diffusivity 
                             //   in vector form: idx = i + j*n

void tlsm_diffusion_and_wesselingh_krishna_model_( 
    double *rho,             // density      scalar
    double *T,               // temperature  scalar
    int    *n,               // # of species scalar
    double *Pc,              // vector of Pc (Pa)
    double *Tc,              // vector of Tc (K)
    double *Vc,              // vector of Vc (cm^3/mol)
    double *MW,              // vector of molecular weights
    double *x,               // vector of molar fractions
    double *Dij);            // matrix of binary mass diffusivity 
                             //   in vector form: idx = i + j*n
                             //
void tlsm_diffusion_trace_ij_( 
    double *rho,             // density      scalar
    double *T,               // temperature  scalar
    int    *n,               // # of species scalar
    double *Pc,              // vector of Pc (Pa)
    double *Tc,              // vector of Tc (K)
    double *Vc,              // vector of Vc (cm^3/mol)
    double *MW,              // vector of molecular weights
    double *x,               // vector of molar fractions
    double *Dij);            // matrix of binary mass diffusivity 
                             //   in vector form: idx = i + j*n
void new_tlsm_diffusion_krishna_model_( 
    double *P,               // pressure     scalar
    double *T,               // temperature  scalar
    int    *n,               // # of species scalar
    double *Pc,              // vector of Pc (Pa)
    double *Tc,              // vector of Tc (K)
    double *Vc,              // vector of Vc (cm^3/mol)
    double *w,               // vector of acentric factors
    int    *tk,              // vector of BIP types
    double *coef_ab,         // vector of coef_ab
    double *MW,              // vector of molecular weights
    double *x,               // vector of molar fractions
    double *Dij);            // matrix of binary mass diffusivity 
                             //   in vector form: idx = i + j*n
                             
void tlsm_diffusion_trace_new_( 
    double *P,               // pressure     scalar
    double *T,               // temperature  scalar
    int    *n,               // # of species scalar
    double *Pc,              // vector of Pc (Pa)
    double *Tc,              // vector of Tc (K)
    double *Vc,              // vector of Vc (cm^3/mol)
    double *w,               // vector of acentric factors
    int    *tk,              // vector of BIP types
    double *coef_ab,         // vector of coef_ab
    double *MW,              // vector of molecular weights
    double *Dij);            // matrix of binary mass diffusivity 
                             //   in vector form: idx = i + j*n
                             
void species_lle_( //
    double *P, //! pressure (Pa)         scalar
    double *T, //! temperature (K)       scalar
    int *n, //! # of species             scalar
    double *Pc,// ! vector of Pc (Pa)
    double *Tc,// ! vector of Tc (K)
    double *Vc,// ! vector of Vc (cm^3/mol)
    double *w, // ! vector of acentric factor
    double *MW,// ! molecular weight (g/mol)
    int *type_k, //! type of kij
    double *xL,//! molar fractions of left  side of interface
    double *xR,//! molar fractions of right side of interface
    double *c0, //! original mixing fractions: c0*xL + (1-c0)*xR
    double *x1,//! equilibrium molar fractions for left side  ! OUTPUT 
    double *x2,//! equilibrium molar fractions for right side ! OUTPUT 
    int *n_miscible); //! >0: miscible; <=0: immiscible ! OUTPUT scalar

void species_lle2_( //
    double *P, //! pressure (Pa)         scalar
    double *T, //! temperature (K)       scalar
    int *n, //! # of species             scalar
    double *Pc,// ! vector of Pc (Pa)
    double *Tc,// ! vector of Tc (K)
    double *Vc,// ! vector of Vc (cm^3/mol)
    double *w, // ! vector of acentric factor
    double *MW,// ! molecular weight (g/mol)
    int *type_k, //! type of kij
    double *xL,//! molar fractions of left  side of interface
    double *xR,//! molar fractions of right side of interface
    double *c0, //! original mixing fractions: c0*xL + (1-c0)*xR
    double *x1,//! equilibrium molar fractions for left side  ! OUTPUT 
    double *x2,//! equilibrium molar fractions for right side ! OUTPUT 
    int *n_miscible); //! >0: miscible; <=0: immiscible ! OUTPUT scalar

void findequilibrium_fix_hc_ratios_( 
    double *P, //! pressure (Pa)         scalar
    double *T, //! temperature (K)       scalar
    int *n, //! # of species             scalar
    double *Pc,// ! vector of Pc (Pa)
    double *Tc,// ! vector of Tc (K)
    double *w, // ! vector of acentric factor
    int *type_k, //! type of kij
    double *xL,//! molar fractions of left  side of interface
    double *xR,//! molar fractions of right side of interface
    double *x1,//! equilibrium molar fractions for left side  ! OUTPUT 
    double *x2,//! equilibrium molar fractions for right side ! OUTPUT 
    int *n_miscible); //! >0: miscible; <=0: immiscible ! OUTPUT scalar

void ideal_gas_chem_potential_(
    double *P, double *T, int *n, double *MW, double *x, double *mu);

void ig_cp_h_( // note that this is for single species; NOT for a vector of species
        double *Tb,double *Tc,double *SG,double *H8,double *MW,
        double *T,double *H_IG,double *CP_IG);

#endif

