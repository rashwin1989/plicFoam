#ifndef                    _PR_EOS_H_
#define                    _PR_EOS_H_

#ifdef __cplusplus
extern "C" {
#endif

/* species' information:

  I. Type of BIP: type_k
  1. PPR78:  idx_groups[j] > 0
    type_k[j] = - (idx_groups[j]*n_pseudospecies_db + idx_species[j]) < 0
  (-type_k[j]) is the line number in groups.dat+ starting from 0

  2. Guang's BIPs: idx_groups[j] < 0
    type_k[j] = -idx_groups[j] > 0
    there are only three species
            water toulene n-decane 
    type_k   1       2        3

  II. Type of ab: coef_ab
  This is an old coefficient aiming to solve asphaltene aggregations
  Its value equals to the molecular weight ratios of an asphaltene over its monomer
  When coef_ab[j]<0, it will be treated as a normal species (USE THIS)
 */

// most important function: calculate_kij to set kij before any change T state
    void calculate_kij_( 
        int    *bSet,          // >0: set kij, <=0: retrieve data
        double *T,             // temperature (Unit: K)
        int    *n,             // number of species
        double *Pc,            // vector of critical pressures
        double *Tc,            // vector of critical temperatures
        double *w,             // vector of acentric factors
        int    *type_k,        // vector of binary interaction types
        double *kij);          // vector of BIP [n*n]

    void calculate_kij_from_table_
    (
        int *bSet,
        double *T,
        int *n,
        double *kij        
    );

    void fugacities2_
    ( 
        double *P,         // pressure (Unit: Pa)
        double *T,         // temperature (Unit: K)
        int    *n,         // number of species
        double *Pc,        // vector of critical pressures
        double *Tc,        // vector of critical temperatures
        double *w,         // vector of acentric factors
        double *x,         // vector of mole fractions
        double *kij,       // matrix of BIPs
        double *fuga,      // vector of fugacity coefficients (output)
        double *Gex        // excess Gibbs free energy (output)
    );

    void fugacities3_
    ( 
        double *P,         // pressure (Unit: Pa)
        double *T,         // temperature (Unit: K)
        int    *n,         // number of species
        double *Pc,        // vector of critical pressures
        double *Tc,        // vector of critical temperatures
        double *w,         // vector of acentric factors
        double *x,         // vector of mole fractions
        double *kij,       // matrix of BIPs
        double *lnphi      // vector of fugacity coefficients (output)
    );

    void fugacities_n_its_derivatives_( 
        double *P,             // pressure (Unit: Pa)
        double *T,             // temperature (Unit: K)
        int    *n,             // number of species
        double *Pc,            // vector of critical pressures
        double *Tc,            // vector of critical temperatures
        double *w,             // vector of acentric factors
        double *x,             // vector of mass fractions
        int    *type_k,        // vector of binary interaction types
        double *coef_ab,       // vector of a, b coefficients
        double *lnphi,         // vector of fugacity coefficients (output)
        double *dlnphi_dxj,    // matrix of dlnphi_i/dx_j(i,j): [idx=i+j*n]
        double *V);            // molar volume (m^3/mol), scalar (output)
    void fugacities_n_its_derivatives_general_( 
        double *P,             // pressure (Unit: Pa)
        double *T,             // temperature (Unit: K)
        int    *n,             // number of species
        double *Pc,            // vector of critical pressures
        double *Tc,            // vector of critical temperatures
        double *w,             // vector of acentric factors
        double *x,             // vector of mass fractions
        int    *type_k,        // vector of binary interaction types
        double *lnphi,         // vector of fugacity coefficients (output)
        double *dlnphi_dxj,    // matrix of dlnphi_i/dx_j(i,j): [idx=i+j*n]
        double *V,             // molar volume (m^3/mol), scalar (output)
        double *am,
        double *dV_dxi, 
        double *dam_dxi, 
        double *d2am_dxidxj);

    void fugacities_n_its_derivatives2_
    ( 
        double *P,         // pressure (Unit: Pa)
        double *T,         // temperature (Unit: K)
        int    *n,         // number of species
        double *Pc,        // vector of critical pressures
        double *Tc,        // vector of critical temperatures
        double *w,         // vector of acentric factors
        double *x,         // vector of mole fractions
        double *kij,       // matrix of BIPs
        double *lnphi,     // vector of fugacity coefficients (output)
        double *dlnphi_dxj // matrix of dlnphi_i/dx_j(i,j): (output)
    );

    void fugacities_n_its_derivatives3_
    ( 
        double *P,          // pressure (Unit: Pa)
        double *T,          // temperature (Unit: K)
        int    *n,          // number of species
        double *Pc,         // vector of critical pressures
        double *Tc,         // vector of critical temperatures
        double *w,          // vector of acentric factors
        double *x,          // vector of mole fractions
        double *kij,        // matrix of BIPs
        double *lnphi,      // vector of fugacity coefficients (output)
        double *dlnphi_dxj, // matrix of dlnphi_i/dx_j(i,j): (output)
        double *V           // molar volume (Unit: m3/mol) (output)
    );

    void fugacities_( 
        double *P, // pressure (Unit: Pa)
        double *T, // temperature (Unit: K)
        int    *n, // number of species
        double *Pc,// vector of critical pressures
        double *Tc,// vector of critical temperatures
        double *w, // vector of acentric factors
        double *x, // vector of mass fractions
        int    *type_k, // vector of binary interaction types
        double *coef_ab,// vector of a, b coefficients
        double *fuga,   // vector of fugacity coefficients (output results)
        double *Gex,    // value of excess Gibbs energy (output results)
        double *delta); // matrix of binary interaction parameters

    void pressure_pr_eos_( //
        double *rho, // ! density (Unit: kg/m^3)
        double *T,   // ! temperature (Unit: K)
        int    *n,   // ! number of species
        double *Pc,  // ! vector of critical pressures
        double *Tc,  // ! vector of critical temperatures
        double *w,   // ! vector of acentric factors
        double *M,   // ! vector of molecular weights
        double *x,   // ! vector of mass fractions
        int    *type_k, // ! vector of binary interaction types
        double *coef_ab,// ! vector of a, b coefficients
        double *P);  // ! pressure (Unit: Pa) (OUTPUT)

    void density_pr_eos_( //
        double *P,   // ! pressure (Unit: Pa) 
        double *T,   // ! temperature (Unit: K)
        int    *n,   // ! number of species
        double *Pc,  // ! vector of critical pressures
        double *Tc,  // ! vector of critical temperatures
        double *w,   // ! vector of acentric factors
        double *M,   // ! vector of molecular weights
        double *x,   // ! vector of mass fractions
        int    *type_k, // ! vector of binary interaction types
        double *coef_ab,// ! vector of a, b coefficients
        double *rho);// ! density (Unit: kg/m^3) (OUTPUT)

    void density_pr_eos2_( //
        double *P,   // ! pressure (Unit: Pa) 
        double *T,   // ! temperature (Unit: K)
        double *x,   // ! vector of mass fractions
        int    *n,   // ! number of species
        double *Pc,  // ! vector of critical pressures
        double *Tc,  // ! vector of critical temperatures
        double *w,   // ! vector of acentric factors
        double *M,   // ! vector of molecular weights
        double *kij, // ! matrix of BIPs
        double *rho);// ! density (Unit: kg/m^3) (OUTPUT)

    // ONLY FOR Binary LLE
    void findequilibrium_fugacity_( 
        double *P, // ! pressure (Unit: Pa)
        double *T, // ! temperature (Unit: K)
        int    *n, // ! number of species
        double *Pc,// ! vector of critical pressures
        double *Tc,// ! vector of critical temperatures
        double *w, // ! vector of acentric factors
        int    *type_k, // ! vector of binary interaction types
        double *coef_ab,// ! vector of a, b coefficients
        double *x_a, // ! vector of mass fractions: alpha phase
        double *x_b);// ! vector of mass fractions: beta phase

    void findequilibrium_search_( 
        double *P, // ! pressure (Unit: Pa)
        double *T, // ! temperature (Unit: K)
        int    *k, // ! the major oil species' index
        int    *n, // ! number of species
        double *Pc,// ! vector of critical pressures
        double *Tc,// ! vector of critical temperatures
        double *w, // ! vector of acentric factors
        int    *type_k, // ! vector of binary interaction types
        double *coef_ab,// ! vector of a, b coefficients
        double *x_a, // ! vector of mass fractions: alpha phase
        double *x_b, // ! vector of mass fractions: beta phase
        double *s_min);// ! minimum evaluation value

    void findequilibrium_new_( 
        int *bPrint, // ! 1: print debugging; 0: no prints
        double *P, // ! pressure (Unit: Pa)
        double *T, // ! temperature (Unit: K)
        int    *k, // ! the major oil species' index
        int    *n, // ! number of species
        double *Pc,// ! vector of critical pressures
        double *Tc,// ! vector of critical temperatures
        double *w, // ! vector of acentric factors
        int    *type_k, // ! vector of binary interaction types
        double *coef_ab,// ! vector of a, b coefficients
        double *x_a, // ! vector of mass fractions: alpha phase
        double *x_b, // ! vector of mass fractions: beta phase
        double *s_min);// ! minimum evaluation value

    void findequilibrium_new2_( 
        double *P, // ! pressure (Unit: Pa)
        double *T, // ! temperature (Unit: K)
        int    *k, // ! the major oil species' index
        int    *n, // ! number of species
        double *Pc,// ! vector of critical pressures
        double *Tc,// ! vector of critical temperatures
        double *w, // ! vector of acentric factors
        double *kij, // ! matrix of BIPs
        double *x_a, // ! vector of mass fractions: alpha phase
        double *x_b, // ! vector of mass fractions: beta phase
        double *s_min);// ! minimum evaluation value

    void findequilibrium_fix_water_b_( 
        int *bPrint, // ! 1: print debugging; 0: no prints
        int *bInitial, // ! 1: do initialization; 0: no init.
        double *P, // ! pressure (Unit: Pa)
        double *T, // ! temperature (Unit: K)
        int    *n, // ! number of species
        double *Pc,// ! vector of critical pressures
        double *Tc,// ! vector of critical temperatures
        double *w, // ! vector of acentric factors
        double *kij, // ! matrix of BIPs
        double *x_a, // ! vector of mass fractions: alpha phase
        double *x_b, // ! vector of mass fractions: beta phase
        double *s_min);// ! minimum evaluation value
    
    void phase_stability2_
    (
        int *stable,//! OUTPUT: 1 - stable; 0 - unstable
        double *lnK,//! initial lnK values for LLE (corresponding to test phase)
        double *P,//! pressure (Unit: Pa)
        double *T,//! temperature (Unit: K)
        int *n,//! number of species
        double *Pc,//! vector of critical pressures
        double *Tc,//! vector of critical temperatures
        double *w,//! vector of acentric factors
        double *kij,//! matrix of BIPs
        double *z//! vector of inital mole fractions
    );

    void species_lle4_
    ( 
        double *P,// ! pressure (Pa)
        double *T,// ! temperature (K)
        int *n,// ! # of species
        double *Pc,// ! vector of Pc (Pa)
        double *Tc,// ! vector of Tc (K)
        double *w,// ! vector of acentric factor
        double *kij,// ! matrix of BIPs
        double *xL,// ! specified initial mole fractions oil phase
        double *xR,// ! specified initial mole fractions oil phase
        double *c1,
        double *c0,//! original mixing fractions: c1*xL+(1-c1)*xR resulting ==> c0*x1 + (1-c0)*x2
        double *x1,//! equilibrium mole fractions oil phase (OUTPUT)
        double *x2,//! equilibrium mole fractions water phase (OUTPUT)
        int *n_miscible//! >0: miscible; <=0: immiscible (OUTPUT)
    );

    void binarylle_( 
        double *P,//! pressure (Pa)
        double *T,//! temperature (K)
        int *n,//! # of species
        double *Pc,// ! vector of Pc (Pa)
        double *Tc,// ! vector of Tc (K)
        double *w,// ! vector of acentric factor
        double *kij,// ! matrix of BIPs
        double *x1,//! equilibrium mole fractions oil phase (OUTPUT)
        double *x2,//! equilibrium mole fractions water phase (OUTPUT)
        int *n_miscible//! >0: miscible; <=0: immiscible (OUTPUT)
    );

    void ln_fugacities_( //
        double *P, // ! pressure (Unit: Pa)
        double *T, // ! temperature (Unit: K)
        int    *n, // ! number of species
        double *Pc,// ! vector of critical pressures
        double *Tc,// ! vector of critical temperatures
        double *w, // ! vector of acentric factors
        double *x, // ! vector of mass fractions
        int    *type_k, // ! vector of binary interaction types
        double *coef_ab,// vector of a, b coefficients
        double *lnphi);   // ! vector of log(fugacity coefficients) (output results)

    void molar_volume_( //
        double *P, // ! pressure (Unit: Pa)
        double *T, // ! temperature (Unit: K)
        int    *n, // ! number of species
        double *Pc,// ! vector of critical pressures
        double *Tc,// ! vector of critical temperatures
        double *w, // ! vector of acentric factors
        double *x, // ! vector of mass fractions
        int    *type_k, // ! vector of binary interaction types
        double *coef_ab,// vector of a, b coefficients
        double *V);     // ! molar volume (output results)

    void molar_volume2_( //
        double *P, // ! pressure (Unit: Pa)
        double *T, // ! temperature (Unit: K)
        int    *n, // ! number of species
        double *Pc,// ! vector of critical pressures
        double *Tc,// ! vector of critical temperatures
        double *w, // ! vector of acentric factors
        double *kij, // ! matrix of BIPs
        double *x, // ! vector of mass fractions
        double *V);     // ! molar volume (output results)

    void pr_phase_( //
        double *P, // ! pressure (Unit: Pa)
        double *T, // ! temperature (Unit: K)
        int    *n, // ! number of species
        double *Pc,// ! vector of critical pressures
        double *Tc,// ! vector of critical temperatures
        double *w, // ! vector of acentric factors
        double *x, // ! vector of mass fractions
        int    *type_k,   // ! vector of binary interaction types
        double *coef_ab,  // vector of a, b coefficients
        int    *phase);   // ! value of phase type: 0: Vapor, 1: Liquid, 2: only phase

    void get_dkdt_( //
        double *P, // ! pressure (Unit: Pa)
        double *T, // ! temperature (Unit: K)
        int    *n, // ! number of species
        double *Pc,// ! vector of critical pressures
        double *Tc,// ! vector of critical temperatures
        double *w, // ! vector of acentric factors
        double *x, // ! vector of mass fractions
        int    *type_k,   // ! vector of binary interaction types
        double *coef_ab,  // vector of a, b coefficients
        double *kij,    // ! matrix of BIP       (output results)
        double *dkdT);  // ! matrix of d(kij)/dT (output results)
    void get_dadt_( //
        double *P, // ! pressure (Unit: Pa)
        double *T, // ! temperature (Unit: K)
        int    *n, // ! number of species
        double *Pc,// ! vector of critical pressures
        double *Tc,// ! vector of critical temperatures
        double *w, // ! vector of acentric factors
        double *x, // ! vector of mass fractions
        int    *type_k,   // ! vector of binary interaction types
        double *coef_ab,  // vector of a, b coefficients
        double *a,  double *dadT, 
        double *ami,double *damidT, 
        double *am, double *damdT);  
    void get_dvdt_( //
        double *P, // ! pressure (Unit: Pa)
        double *T, // ! temperature (Unit: K)
        int    *n, // ! number of species
        double *Pc,// ! vector of critical pressures
        double *Tc,// ! vector of critical temperatures
        double *w, // ! vector of acentric factors
        double *x, // ! vector of mass fractions
        int    *type_k,   // ! vector of binary interaction types
        double *coef_ab,  // vector of a, b coefficients
        double *V,  double *dVdT);

        #ifdef __cplusplus
}
#endif

#endif
