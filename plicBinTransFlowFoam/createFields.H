//------------------------------------------------------------------//
//Define all label, scalar, vector (OpenFOAM format) variables
//Input parameters, tolerances
label nSpecies, nOCorr, nPCorr, nNonOrthCorr, ALPHA_BOUND_ITERS_MAX, Y_BOUND_ITERS_MAX;
scalar ALPHA_MAX_BOUND, ALPHA_MIN_BOUND, ALPHA_BOUND_ERR, SMALLEST_ALPHA_FOR_YADV, SMALLEST_ALPHA_FOR_YDIFF, T1_TMP_VAL, T0_TMP_VAL,MIN_ALPHA_DIFF, /*ALPHA_2PH_MIN, */ALPHA_2PH_MIN_IST, ALPHA_2PH_MAX_IST;
//Temp variables
label iOCorr, iPCorr, iNonOrthCorr, faceI, celli;
label nCells, nFaces, nIntFaces, n_flatFld_diff, n_flatFld_ist;
scalar deltaT, dt, deltaT_diff, maxAlpha1, minAlpha1, div_alpha1, div_phi;
vector minAlpha1CellC, maxAlpha1CellC;
word fieldName;
fileName outputFile;
//Done define all label, scalar, vector (OpenFOAM format) variables
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Define all int, bool, double (standard C++) variables
//Input parameters, tolerances
int MAX_ITER_T, n_flux_type, MAX_ITER_Ts;
bool alpha_debug, redist_alpha, adv_debug, adv_debug2, redist_C, grad_debug, diff_debug, diff_debug2, diff_debug_T, ist_debug, ist_debug2, momentumPredictor, advFlux_debug, advFlux_debug2, pEqn_debug, U_debug, isothermal, advFlux_debugH, interp_debug, ist_ON, transProp_debug, useCombinedHeatMassDiffusion;
//bool calcKij;
double T_TOL, H_TOL, R_gas, Ts_TOL, MASS_FRAC_TOL;
//Temp variables
int n, i, j, ijIndex, maxAlpha1Cell, minAlpha1Cell, nT_kij, iT, idx, iread, n_miscible;
double T_0, P_thermo_tmp, pSCoeff, Y1i_cellI, Y1t_cellI, Y0i_cellI, Y0t_cellI ,Y1n_cellI, Y0n_cellI, C1Tot_cellI, C0Tot_cellI, C1i_cellI, C0i_cellI, alpha1_cellI, alpha0_cellI, rho_tmp, mu_tmp, Ta_kij, Tb_kij, T_cellI_tmp, h_cellI_tmp, rhoh_cellI_tmp;
int nErf, iErf;
double erf_a, erf_b;
double *erfInv_table;
bool useErf;
bool calcDiff;
FILE *f;
//Done define all int, bool, double (standard C++) variables
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Define all int, bool, double (standard C++) arrays/pointers
//Input and parameter arrays
double *Pc, *Tc, *Vc, *w, *MW, *Tb, *SG, *H8, *k, *dm;
double *kij_T, *kij;
//Temp variable arrays
double *y_cellI_tmp, *x_cellI_tmp, *Dij1_tmp, *Dij0_tmp;
double *xs1_tmp, *xs0_tmp, *ys1_tmp, *ys0_tmp;
//Done define all int, bool, double (standard C++) arrays/pointers
//------------------------------------------------------------------//

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//------------------------------------------------------------------//
//Define and load the OpenFOAM input dinctionaries
//to read parameters from
IOdictionary transPropDict
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        p_rgh.db(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const dictionary& solDict = mesh.solutionDict();
//Done define and load the OpenFOAM input dinctionaries
//to read parameters from
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Reading tolerances, max iters and switches
isothermal = transPropDict.lookupOrDefault("isothermal", true);
useCombinedHeatMassDiffusion = transPropDict.lookupOrDefault("useCombinedHeatMassDiffusion", false);
// iters, correctors
nOCorr = solDict.lookupOrDefault("nOuterCorrector", 1);
nPCorr = solDict.lookupOrDefault("nPressureCorrector", 2);
nNonOrthCorr = solDict.lookupOrDefault("nNonOrthogonalCorrector", 1);
MAX_ITER_T = transPropDict.lookupOrDefault("MAX_ITER_T", 108);
if(isothermal)
{
    MAX_ITER_Ts = 1;
    Ts_TOL = 1000;
}
else
{
    MAX_ITER_Ts = transPropDict.lookupOrDefault("MAX_ITER_Ts", 108);
    Ts_TOL = transPropDict.lookupOrDefault("Ts_TOL", 1E-3);
}
ALPHA_BOUND_ITERS_MAX = 12;
Y_BOUND_ITERS_MAX = 12;
// tolerances
ALPHA_BOUND_ERR = transPropDict.lookupOrDefault("ALPHA_BOUND_ERR", 1E-15);
SMALLEST_ALPHA_FOR_YADV = transPropDict.lookupOrDefault("SMALLEST_ALPHA_FOR_YADV", 1E-15);
T_TOL = transPropDict.lookupOrDefault("T_TOL", 1E-3);
H_TOL = transPropDict.lookupOrDefault("H_TOL", 0.1);
//MASS_FRAC_TOL = transPropDict.lookupOrDefault("MASS_FRAC_TOL", 1E-3);
SMALLEST_ALPHA_FOR_YDIFF = transPropDict.lookupOrDefault("SMALLEST_ALPHA_FOR_YDIFF", 1E-15);
MIN_ALPHA_DIFF = transPropDict.lookupOrDefault("MIN_ALPHA_DIFF", 1E-09);
//ALPHA_2PH_MIN = transPropDict.lookupOrDefault("ALPHA_2PH_MIN", 1E-09);
//ALPHA_2PH_MAX = 1 - ALPHA_2PH_MIN;
ALPHA_2PH_MIN_IST = transPropDict.lookupOrDefault("ALPHA_2PH_MIN_IST", 1E-06);
ALPHA_2PH_MAX_IST = 1 - ALPHA_2PH_MIN_IST;
//A_INTFC_2PH_MIN_COEFF = transPropDict.lookupOrDefault("A_INTFC_2PH_MIN_COEFF", 1E-03);
//scalar avgSf = (mesh.magSf()).average().value();
//A_INTFC_2PH_MIN = A_INTFC_2PH_MIN_COEFF*avgSf;
// debug and other bool switches
alpha_debug = transPropDict.lookupOrDefault("alpha_debug", false);
redist_alpha = transPropDict.lookupOrDefault("redist_alpha", false);
adv_debug = transPropDict.lookupOrDefault("adv_debug", false);
adv_debug2 = transPropDict.lookupOrDefault("adv_debug2", false);
redist_C = transPropDict.lookupOrDefault("redist_C", false);
grad_debug = transPropDict.lookupOrDefault("grad_debug", false);
//grad_debug2 = transPropDict.lookupOrDefault("grad_debug2", false);
diff_debug = transPropDict.lookupOrDefault("diff_debug", false);
diff_debug2 = transPropDict.lookupOrDefault("diff_debug2", false);
diff_debug_T = transPropDict.lookupOrDefault("diff_debug_T", false);
ist_debug = transPropDict.lookupOrDefault("ist_debug", false);
ist_debug2 = transPropDict.lookupOrDefault("ist_debug2", false);
momentumPredictor = solDict.lookupOrDefault("momentumPredictor", false);
advFlux_debug = transPropDict.lookupOrDefault("advFlux_debug", false);
advFlux_debug2 = transPropDict.lookupOrDefault("advFlux_debug2", false);
advFlux_debugH = transPropDict.lookupOrDefault("advFlux_debugH", false);
pEqn_debug = transPropDict.lookupOrDefault("pEqn_debug", false);
U_debug = transPropDict.lookupOrDefault("U_debug", false);
interp_debug = transPropDict.lookupOrDefault("interp_debug", false);
transProp_debug = transPropDict.lookupOrDefault("transProp_debug", false);
useErf = transPropDict.lookupOrDefault("useErf", false);
calcDiff = transPropDict.lookupOrDefault("calcDiff", true);
ist_ON = transPropDict.lookupOrDefault("ist_ON", true);
// int switches
//n_flux_type = transPropDict.lookupOrDefault("n_flux_type", 1);
//Done reading tolerances, max iters and switches
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Read and set no. of species
//arrays and lists need to be initialized with this size
nSpecies = transPropDict.lookupOrDefault("nSpecies", 2);
Info<< "Number of species " << nSpecies << endl;
n = nSpecies;
R_gas = transPropDict.lookupOrDefault("R_GAS", 8.31412);
//Done read and set no. of species
//arrays and lists need to be initialized with this size
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Allocate memory for arrays
_NNEW2_(Pc, double, n);
_NNEW2_(Tc, double, n);
_NNEW2_(Vc, double, n);
_NNEW2_(w, double, n);
_NNEW2_(MW, double, n);
_NNEW2_(Tb, double, n);
_NNEW2_(SG, double, n);
_NNEW2_(H8, double, n);
_NNEW2_(k, double, n);
_NNEW2_(dm, double, n);
_NNEW2_(y_cellI_tmp, double, n);
_NNEW2_(x_cellI_tmp, double, n);
_NNEW2_(Dij1_tmp, double, n*n);
_NNEW2_(Dij0_tmp, double, n*n);
_NNEW2_(xs1_tmp, double, n);
_NNEW2_(xs0_tmp, double, n);
_NNEW2_(ys1_tmp, double, n);
_NNEW2_(ys0_tmp, double, n);
//Done allocate memory for arrays
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Read input parameters for BIP table from input file
f = fopen("constant/kij_input.dat+", "r");
iread = fscanf(f, "%lf %lf %d", &Ta_kij, &Tb_kij, &nT_kij);
fclose(f);

_NNEW2_(kij_T, double, nT_kij*n*n);
_NNEW2_(kij, double, n*n);

f = fopen("constant/kij_T.dat+", "r");
for(iT=0; iT<nT_kij; iT++)
{
    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
        {
            idx = iT*n*n + i*n + j;
            iread = fscanf(f, "%lf", &kij_T[idx]);
        }
    }
}
fclose(f);
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Read thermo input parameters from files
int n_species_db;
char str_tmp[50];
int n_pseudospecies_db, n_purespecies_db;
int *idx_species, *idx_groups;

_NNEW2_(idx_species, int, n);
_NNEW2_(idx_groups, int, n);

f = fopen("constant/species.dat+", "r");

iread = fscanf(f, "%d %d", &n_pseudospecies_db, &n_purespecies_db);
if(iread <= 0) Info<< "Input file reading error-----------------" << endl;
n_species_db = n_pseudospecies_db + n_purespecies_db;

Info<< "n_pseudoSpecies_db = " << n_pseudospecies_db << "  n_purespecies_db = " << n_purespecies_db << endl;

fclose(f);

f = fopen("constant/input.dat+", "r");

iread=fscanf(f, "%s %d", str_tmp, &n);  // num of species
if(iread <= 0) Info<< "Input file reading error-----------------" << endl;

iread=fscanf(f, "%s", str_tmp);         // index of species
for (i=0; i<n; i++) iread=fscanf(f, "%d", &idx_species[i]);
if(iread <= 0) Info<< "Input file reading error-----------------" << endl;

iread=fscanf(f, "%s", str_tmp);         // index of PPR78 groups
for (i=0; i<n; i++) iread=fscanf(f, "%d", &idx_groups[i]);
if(iread <= 0) Info<< "Input file reading error-----------------" << endl;

Info<< "Species indices in database:  ";
for(i=0; i<n; i++) Info<< idx_species[i] << "  ";
Info<< endl;

Info<< "Species group indices in database:  ";
for(i=0; i<n; i++) Info<< idx_groups[i] << "  ";
Info<< endl;

fclose(f);

f = fopen("constant/petro.dat+","r");

for (j=0;j<n;j++) 
{
    for (i=0;i<n_species_db;i++) 
    {
        double tmp;
        double xPNA[3];

        iread=fscanf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
        &tmp, &tmp, &Tb[j], &Tc[j], &Pc[j], &w[j], &MW[j], &Vc[j], 
        &SG[j], &H8[j], &xPNA[0], &xPNA[1], &xPNA[2]);

        if(iread <= 0) Info<< "Input file reading error-----------------" << endl;

        Pc[j] *= 1e5;  // bar --> Pa        

        //x1[j] *= xPNA[ idx_groups[j] ];

        dm[j] = 0.0; 
        k [j] = 0.0;

        if (i==idx_species[j]) 
        {
            // SPECIAL DATA POINTS BUILT IN THE CODE
            // 1. water
            if (fabs(Tc[j]-647)<0.5 && fabs(Pc[j]-220.6e5)<0.5e5) 
            {
                Vc[j] = 57.1; // cm^3/mol
                dm[j] = 1.85; // Debye
                k [j] = 0.076;

            }                                
            // 2. toluene
            else if (fabs(Tc[j]-591.8)<0.5 && fabs(Pc[j]-41.1e5)<0.5e5)
            {
                    Vc[j] = 316.0;
                    dm[j] = 0.36; // Debye

            }// 3. n-decane
            else if (fabs(Tc[j]-617.7)<0.5 && fabs(Pc[j]-21.1e5)<0.5e5)
            {
                Vc[j] = 600.0;
            }// 4. n-C30
            else if (fabs(Tc[j]-844.0)<0.5 && fabs(Pc[j]- 8.0e5)<0.5e5) {
                Vc[j] = 1805.0;

                // 5. n-C50
            }else if (fabs(Tc[j]-1073.6)<0.5&& fabs(Pc[j]-3.51e5)<0.5e5) {
                Vc[j] = 2999;

                // 6. Benzene-C10
            }else if (fabs(Tc[j]-753.0)<0.5 && fabs(Pc[j]-17.7e5)<0.5e5) {
                Vc[j] = 813;

                // 7. Benzene-C30
            }else if (fabs(Tc[j]-965.6)<0.5 && fabs(Pc[j]- 5.3e5)<0.5e5) {
                Vc[j] = 1943.5;

                // 8. Naphthalene-C10
            }else if (fabs(Tc[j]-859.0)<0.5 && fabs(Pc[j]-15.8e5)<0.5e5) {
                Vc[j] = 1070;

                // 9. Naphthalene-C12
            }else if (fabs(Tc[j]-854.6)<0.5 && fabs(Pc[j]-13.0e5)<0.5e5) {
                Vc[j] = 1081.5;

                // 10. Benzene
            }else if (fabs(Tc[j]-562.0)<0.5 && fabs(Pc[j]-49.0e5)<0.5e5) {
                Vc[j] = 256.0;

                // 11. o-xylene
            }else if (fabs(Tc[j]-630.3)<0.5 && fabs(Pc[j]-37.3e5)<0.5e5) {
                Vc[j] = 370.0;

                // 12. p-xylene
            }else if (fabs(Tc[j]-616.2)<0.5 && fabs(Pc[j]-35.1e5)<0.5e5) {
                Vc[j] = 378.0;

                // 13. 1,3,5-trimethylbenzene
            }else if (fabs(Tc[j]-637.3)<0.5 && fabs(Pc[j]-31.3e5)<0.5e5) {
                Vc[j] = 433.0;

                // 14. naphthalene
            }else if (fabs(Tc[j]-748.4)<0.5 && fabs(Pc[j]-40.5e5)<0.5e5) {
                Vc[j] = 407.0;

                // 15. 1-methylnaphthalene
            }else if (fabs(Tc[j]-772.0)<0.5 && fabs(Pc[j]-36.0e5)<0.5e5) {
                Vc[j] = 465;

                // 16. anthracene
            }else if (fabs(Tc[j]-873.0)<0.5 && fabs(Pc[j]-29.0e5)<0.5e5) {
                Vc[j] = 554;

                // 17. 1,2-diphenylethane
            }else if (fabs(Tc[j]-780.0)<0.5 && fabs(Pc[j]-26.5e5)<0.5e5) {
                Vc[j] = 616;

                // 18. pyrene
            }else if (fabs(Tc[j]-936.0)<0.5 && fabs(Pc[j]-26.1e5)<0.5e5) {
                Vc[j] = 660;

                // 19. n-C16
            }else if (fabs(Tc[j]-723.0)<0.5 && fabs(Pc[j]-14.0e5)<0.5e5) {
                Vc[j] = 969.2;

                // 20. trans-decalin
            }else if (fabs(Tc[j]-687.0)<0.5 && fabs(Pc[j]-32.0e5)<0.5e5) {
                Vc[j] = 480.0;

                // 21. butylbenzene
            }else if (fabs(Tc[j]-660.5)<0.5 && fabs(Pc[j]-28.9e5)<0.5e5) {
                Vc[j] = 497.0;

                // 22. hexylbenzene
            }else if (fabs(Tc[j]-698.0)<0.5 && fabs(Pc[j]-23.8e5)<0.5e5) {
                Vc[j] = 593.0;
            }
            // WHEN EQUALS, WE FOUND THE DATA 
            // SO QUIT THE READING FOR SPECIES I
            break;
        }
    }
    rewind(f);
}

fclose(f);
//Done read thermo input parameters from files
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Initialize UMFPACK
T_UMFPACK flux_umf;
initialUmfpack(&flux_umf);
flux_memoryMatVec(&flux_umf, n);
flux_compRow_setRowCol(&flux_umf, n);
//Done initialize UMFPACK
//------------------------------------------------------------------//
Info<< "Done initialize UMFPACK" << endl;
//------------------------------------------------------------------//
//Set up output streams linked to output files
outputFile = "plic_dbg";
OFstream osPlic(runTime.path()/outputFile);
osPlic.precision(5);
osPlic.setf(std::ios::scientific);
//osPlic.setf(std::ios::showpos);
osPlic.setf(std::ios::left);

outputFile = "alpha_dbg";
OFstream osAlpha(runTime.path()/outputFile);
osAlpha.precision(5);
osAlpha.setf(std::ios::scientific);
//osAlpha.setf(std::ios::showpos);
osAlpha.setf(std::ios::left);

outputFile = "adv_dbg";
OFstream osAdv(runTime.path()/outputFile);
osAdv.precision(5);
osAdv.setf(std::ios::scientific);
//osAdv.setf(std::ios::showpos);
osAdv.setf(std::ios::left);

outputFile = "advH_dbg";
OFstream osAdvH(runTime.path()/outputFile);
osAdv.precision(5);
osAdv.setf(std::ios::scientific);
//osAdv.setf(std::ios::showpos);
osAdv.setf(std::ios::left);

outputFile = "grad_dbg";
OFstream osGrad(runTime.path()/outputFile);
osGrad.precision(5);
osGrad.setf(std::ios::scientific);
//osGrad.setf(std::ios::showpos);
osGrad.setf(std::ios::left);

outputFile = "diff_dbg";
OFstream osDiff(runTime.path()/outputFile);
osDiff.precision(5);
osDiff.setf(std::ios::scientific);
//osDiff.setf(std::ios::showpos);
osDiff.setf(std::ios::left);

outputFile = "ist_dbg";
OFstream osIst(runTime.path()/outputFile);
osIst.precision(5);
osIst.setf(std::ios::scientific);
//osIst.setf(std::ios::showpos);
osIst.setf(std::ios::left);

outputFile = "pEqn_dbg";
OFstream osPEqn(runTime.path()/outputFile);
osPEqn.precision(5);
osPEqn.setf(std::ios::scientific);
//osPEqn.setf(std::ios::showpos);
osPEqn.setf(std::ios::left);

outputFile = "U_dbg";
OFstream osU(runTime.path()/outputFile);
osU.precision(5);
osU.setf(std::ios::scientific);
//osU.setf(std::ios::showpos);
osU.setf(std::ios::left);

outputFile = "transProp_dbg";
OFstream osTransProp(runTime.path()/outputFile);
osTransProp.precision(5);
osTransProp.setf(std::ios::scientific);
//osTransProp.setf(std::ios::showpos);
osTransProp.setf(std::ios::left);

//Done set up output streams linked to output files
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Set mesh related variables and fields
nCells = mesh.nCells();
nFaces = mesh.nFaces();
nIntFaces = mesh.nInternalFaces();

const labelList& own = mesh.owner();
const labelList& nei = mesh.neighbour();
const vectorField& meshCi = mesh.C().internalField();
const scalarField& meshV = mesh.V();
//Done set mesh related variables and fields
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Read and initialize fundamental variable fields
P_thermo_tmp = transPropDict.lookupOrDefault("P_thermo", 30.0E+06);

Info<< "Initializing field P_thermo\n" << endl;
volScalarField P_thermo
(
    IOobject
    (
        "P_thermo",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    P_thermo_tmp
);

Info<< "Initializing field Pf_thermo\n" << endl;
surfaceScalarField Pf_thermo
(
    IOobject
    (
        "Pf_thermo",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    P_thermo_tmp
);

Info<< "Reading field alpha1\n" << endl;
volScalarField alpha1
(
    IOobject
    (
        "alpha1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Initializing field alpha1_old\n" << endl;
volScalarField alpha1_old
(
    IOobject
    (
        "alpha1_old",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1,
    alpha1.boundaryField().types()
);

Info<< "Calculating field alpha0\n" << endl;
volScalarField alpha0("alpha0", scalar(1) - alpha1);
volScalarField alpha0_old = alpha0;

volScalarField limitedAlpha1
(
    "limitedAlpha1", 
    min(max(alpha1, scalar(0)), scalar(1))
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"
U.correctBoundaryConditions();
phi = linearInterpolate(U) & mesh.Sf();

surfaceScalarField alpha1f = fvc::interpolate(limitedAlpha1, phi, "alpha1f");
surfaceScalarField alpha0f = scalar(1) - alpha1f;
//Done read and initialize fundamental variable fields
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//T initialization
Info<< "Reading field T1" << endl;
volScalarField T1
(
    IOobject
    (
        "T1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
T1.correctBoundaryConditions();

Info<< "Reading field T0" << endl;
volScalarField T0
(
    IOobject
    (
        "T0",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
T0.correctBoundaryConditions();

T_0 = T1.internalField()[0];
//Done T initialization
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Thermo and transport properties fields initialization
Info<< "Initializing transportProperties\n" << endl;

List<scalar> Y1_0(n);
List<scalar> Y0_0(n);

scalar Y1i_tmp, Y0i_tmp; 
scalar Y1Tot_tmp = 0;
scalar Y0Tot_tmp = 0;
for(i=0; i<(n-1); i++)
{    
    fieldName = "Y1_0[" + Foam::name(i) + "]";
    Y1i_tmp = transPropDict.lookupOrDefault(fieldName, 1.0);
    Y1_0[i] = Y1i_tmp;
    Y1Tot_tmp += Y1i_tmp;

    fieldName = "Y0_0[" + Foam::name(i) + "]";
    Y0i_tmp = transPropDict.lookupOrDefault(fieldName, 0.0);
    Y0_0[i] = Y0i_tmp;
    Y0Tot_tmp += Y0i_tmp;
}
Y1_0[n-1] = scalar(1.0) - Y1Tot_tmp;
Y1_0[n-1] = max(Y1_0[n-1], 0.0);
Y0_0[n-1] = scalar(1.0) - Y0Tot_tmp;
Y0_0[n-1] = max(Y0_0[n-1], 0.0);

//Since isothermal case, calculate BIP matrix from input table once for T_0
plicFuncs::calc_kij_from_table(T_0, n, Ta_kij, Tb_kij, nT_kij, kij_T, kij);
//calcKij = false;
Info<< "BIP matrix:  ";
for(i=0; i<n*n; i++) Info<< kij[i] << "  ";
Info<< endl;

//initialization of values for thermal conductivities
const dimensionedScalar& lambda1_0 = transPropDict.lookup("lambda1");
const dimensionedScalar& lambda0_0 = transPropDict.lookup("lambda0");
//initialization of values for specific heats
const dimensionedScalar& Cp1_0 = transPropDict.lookup("Cp1");
const dimensionedScalar& Cp0_0 = transPropDict.lookup("Cp0");
// rho and mu are calculated (below)

//phase-1
//set pure oil phase initial mole fractions
x_cellI_tmp[0] = 1; x_cellI_tmp[1] = 0;
y_cellI_tmp[0] = 1; y_cellI_tmp[1] = 0;

plicFuncs::calc_rho(P_thermo_tmp,T_0,x_cellI_tmp,n,Pc,Tc,w,MW,kij,rho_tmp);
dimensionedScalar rho1_0("rho1_0", dimDensity, rho_tmp);

plicFuncs::calc_mu(P_thermo_tmp,T_0,y_cellI_tmp,n,Pc,Tc,Vc,w,MW,Tb,SG,H8,k,dm,kij,mu_tmp);
dimensionedScalar mu1_0("mu1_0", dimPressure*dimTime, mu_tmp);

//new_tlsm_diffusion_krishna_model_(&P_thermo_tmp,&T_0,&n,Pc,Tc,Vc,w,MW,kij,x_cellI_tmp,Dij1_tmp);

//phase-0
//set pure water phase initial mole fractions
x_cellI_tmp[0] = 0; x_cellI_tmp[1] = 1;
y_cellI_tmp[0] = 0; y_cellI_tmp[1] = 1;

plicFuncs::calc_rho(P_thermo_tmp,T_0,x_cellI_tmp,n,Pc,Tc,w,MW,kij,rho_tmp);
dimensionedScalar rho0_0("rho0_0", dimDensity, rho_tmp);

plicFuncs::calc_mu(P_thermo_tmp,T_0,y_cellI_tmp,n,Pc,Tc,Vc,w,MW,Tb,SG,H8,k,dm,kij,mu_tmp);
dimensionedScalar mu0_0("mu0_0", dimPressure*dimTime, mu_tmp);

//new_tlsm_diffusion_krishna_model_(&P_thermo_tmp,&T_0,&n,Pc,Tc,Vc,w,MW,kij,x_cellI_tmp,Dij0_tmp);

volScalarField rho1
(
    IOobject
    (
        "rho1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    rho1_0
);

volScalarField v1
(
    IOobject
    (
        "v1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    1e-3*MW[0]/rho1_0
);

volScalarField rho0
(
    IOobject
    (
        "rho0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    rho0_0
);

volScalarField v0
(
    IOobject
    (
        "v0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    1e-3*MW[1]/rho0_0
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1*rho1 + alpha0*rho0,
    alpha1.boundaryField().types()
);

volScalarField mu1
(
    IOobject
    (
        "mu1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    mu1_0
);

volScalarField mu0
(
    IOobject
    (
        "mu0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    mu0_0
);

volScalarField mu1_thermo = mu1;
volScalarField mu0_thermo = mu0;

surfaceScalarField mu1f(linearInterpolate(mu1));
surfaceScalarField mu0f(linearInterpolate(mu0));

volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1*mu1 + alpha0*mu0,
    alpha1.boundaryField().types()
);

// Mixture viscosity at cell faces
surfaceScalarField muf(alpha1f*mu1f + alpha0f*mu0f);

// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

volScalarField lambda1
(
    IOobject
    (
        "lambda1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    lambda1_0
);

volScalarField lambda0
(
    IOobject
    (
        "lambda0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    lambda0_0
);
// Interpolation for lambdas done in diff_grad_interp.H

volScalarField Cp1
(
    IOobject
    (
        "Cp1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    Cp1_0
);

volScalarField Cp0
(
    IOobject
    (
        "Cp0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    Cp0_0
);

volScalarField Dh1 = lambda1/rho1/Cp1;
volScalarField Dh0 = lambda0/rho0/Cp0;

volScalarField dgdt
(
    pos(alpha0)*fvc::div(phi)/max(alpha0, scalar(0.0001))
);

int Nij = n*n;
PtrList<dimensionedScalar> D1_0(n);
PtrList<dimensionedScalar> D0_0(n);
PtrList<volScalarField> Dij1(Nij);
PtrList<volScalarField> Dij0(Nij);
PtrList<volScalarField> D1(n);
PtrList<volScalarField> D0(n);
PtrList<surfaceScalarField> D1f(n);
PtrList<surfaceScalarField> D0f(n);

//Done thermo and transport properties fields initialization
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Species concentrations Y, x, C initialization
PtrList<volScalarField> Y1(n);
PtrList<volScalarField> x1(n);
PtrList<volScalarField> c1(n);
PtrList<volScalarField> C1(n);
PtrList<volScalarField> C1_old(n);
PtrList<volScalarField> c1_old(n);
PtrList<volScalarField> Y0(n);
PtrList<volScalarField> x0(n);
PtrList<volScalarField> c0(n);
PtrList<volScalarField> C0(n);
PtrList<volScalarField> C0_old(n);
PtrList<volScalarField> c0_old(n);

PtrList<volScalarField> C_phAvg(n);

for(i=0; i<(n-1); i++)
{
    fieldName = "Y1" + Foam::name(i);
    Info<< "Reading field " << fieldName << endl;
    Y1.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
    Y1[i].correctBoundaryConditions();

    fieldName = "x1" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    x1.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            Y1[i]
        )
    );

    fieldName = "c1" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    c1.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rho1*Y1[i],
            Y1[i].boundaryField().types()
        )
    );

    fieldName = "C1" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    C1.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            alpha1*rho1*Y1[i],
            Y1[i].boundaryField().types()
        )
    );

    fieldName = "C1_old" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    C1_old.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            C1[i],
            Y1[i].boundaryField().types()
        )
    );

    fieldName = "c1_old" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    c1_old.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rho1*Y1[i],
            Y1[i].boundaryField().types()
        )
    );

    fieldName = "Y0" + Foam::name(i);
    Info<< "Reading field " << fieldName << endl;
    Y0.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
    Y0[i].correctBoundaryConditions();

    fieldName = "x0" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    x0.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            Y0[i]
        )
    );

    fieldName = "c0" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    c0.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rho0*Y0[i],
            Y0[i].boundaryField().types()
        )
    );
    
    fieldName = "C0" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    C0.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            alpha0*rho0*Y0[i],
            Y0[i].boundaryField().types()
        )
    );

    fieldName = "C0_old" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    C0_old.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            C0[i],
            Y0[i].boundaryField().types()
        )
    );

    fieldName = "c0_old" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    c0_old.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rho0*Y0[i],
            Y0[i].boundaryField().types()
        )
    );

    fieldName = "C_phAvg" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    C_phAvg.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            C0[i] + C1[i],
            Y1[i].boundaryField().types()
        )
    );
}

volScalarField Yt
(
    IOobject
    (
        "Yt",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField Y1t = Yt;
volScalarField Y0t = Yt;

for(i=0; i<(n-1); i++)
{
    Y1t -= Y1[i];
    Y0t -= Y0[i];
}

int endIndex = n-1;

fieldName = "Y1" + Foam::name(endIndex);
Info<< "Reading field " << fieldName << endl;
Y1.set
(
    endIndex,
    new volScalarField
    (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Y1t
    )
);
Y1[endIndex].correctBoundaryConditions();

fieldName = "x1" + Foam::name(endIndex);
Info<< "Initializing field " << fieldName << endl;
x1.set
(
    endIndex,
    new volScalarField
    (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Y1[endIndex]
    )
);

fieldName = "c1" + Foam::name(endIndex);
Info<< "Initializing field " << fieldName << endl;
c1.set
(
    endIndex,
    new volScalarField
    (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho1*Y1[endIndex],
        Y1[endIndex].boundaryField().types()
    )
);

fieldName = "C1" + Foam::name(endIndex);
Info<< "Initializing field " << fieldName << endl;
C1.set
(
    endIndex,
    new volScalarField
    (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1*rho1*Y1[endIndex],
        Y1[endIndex].boundaryField().types()
    )
);

fieldName = "C1_old" + Foam::name(endIndex);
Info<< "Initializing field " << fieldName << endl;
C1_old.set
(
    endIndex,
    new volScalarField
    (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        C1[endIndex],
        Y1[endIndex].boundaryField().types()
    )
);

fieldName = "c1_old" + Foam::name(endIndex);
Info<< "Initializing field " << fieldName << endl;
c1_old.set
(
    endIndex,
    new volScalarField
    (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho1*Y1[endIndex],
        Y1[endIndex].boundaryField().types()
    )
);

fieldName = "Y0" + Foam::name(endIndex);
Info<< "Reading field " << fieldName << endl;
Y0.set
(
    endIndex,
    new volScalarField
    (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Y0t
    )
);
Y0[endIndex].correctBoundaryConditions();

fieldName = "x0" + Foam::name(endIndex);
Info<< "Initializing field " << fieldName << endl;
x0.set
(
    endIndex,
    new volScalarField
    (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Y0[endIndex]
    )
);

fieldName = "c0" + Foam::name(endIndex);
Info<< "Initializing field " << fieldName << endl;
c0.set
(
    endIndex,
    new volScalarField
    (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho0*Y0[endIndex],
        Y0[endIndex].boundaryField().types()
    )
);

fieldName = "C0" + Foam::name(endIndex);
Info<< "Initializing field " << fieldName << endl;
C0.set
(
    endIndex,
    new volScalarField
    (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha0*rho0*Y0[endIndex],
        Y0[endIndex].boundaryField().types()
    )
);

fieldName = "C0_old" + Foam::name(endIndex);
Info<< "Initializing field " << fieldName << endl;
C0_old.set
(
    endIndex,
    new volScalarField
    (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        C0[endIndex],
        Y0[endIndex].boundaryField().types()
    )
);

fieldName = "c0_old" + Foam::name(endIndex);
Info<< "Initializing field " << fieldName << endl;
c0_old.set
(
    endIndex,
    new volScalarField
    (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho0*Y0[endIndex],
        Y0[endIndex].boundaryField().types()
    )
);

fieldName = "C_phAvg" + Foam::name(endIndex);
Info<< "Initializing field " << fieldName << endl;
C_phAvg.set
(
    endIndex,
    new volScalarField
    (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        C0[endIndex] + C1[endIndex],
        Y1[endIndex].boundaryField().types()
    )
);
//Done species concentrations Y, x, C initialization
//------------------------------------------------------------------//


//------------------------------------------------------------------//
//Enthalpy initialization

Info<< "Initializing field h1" << endl;
volScalarField h1
(
    IOobject
    (
        "h1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionSet(0, 2, -2, 0, 0, 0, 0)
);

Info<< "Initializing field rhoh1" << endl;
volScalarField rhoh1
(
    IOobject
    (
        "rhoh1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho1*h1,
    T1.boundaryField().types()
);

Info<< "Initializing field H1" << endl;
volScalarField H1
(
    IOobject
    (
        "H1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha1*rhoh1,
    T1.boundaryField().types()
);
volScalarField H1_old
(
    IOobject
    (
        "H1_old",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    H1,
    T1.boundaryField().types()
);

Info<< "Initializing field h0" << endl;
volScalarField h0
(
    IOobject
    (
        "h0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionSet(0, 2, -2, 0, 0, 0, 0)
);

Info<< "Initializing field rhoh0" << endl;
volScalarField rhoh0
(
    IOobject
    (
        "rhoh0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho0*h0,
    T0.boundaryField().types()
);

Info<< "Initializing field H0" << endl;
volScalarField H0
(
    IOobject
    (
        "H0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha0*rhoh0,
    T0.boundaryField().types()
);
volScalarField H0_old
(
    IOobject
    (
        "H0_old",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    H0,
    T0.boundaryField().types()
);

PtrList<volScalarField> hpar1(n);
PtrList<volScalarField> hpar0(n);

for(i=0; i<n; i++)
{
    fieldName = "hpar1" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    hpar1.set
        (
            i,
        new volScalarField
            (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 2, -2, 0, 0, 0, 0)
            )
        );

    fieldName = "hpar0" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    hpar0.set
        (
            i,
        new volScalarField
            (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 2, -2, 0, 0, 0, 0)
            )
        );
}
//Done Enthalpy initialization
//------------------------------------------------------------------//


//------------------------------------------------------------------//
//Max and min values for scalar fields
List<scalar> C1max(n);
List<scalar> C0max(n);
List<scalar> C1min(n);
List<scalar> C0min(n);
List<scalar> Y1max(n);
List<scalar> Y0max(n);
List<scalar> Y1min(n);
List<scalar> Y0min(n);

for(i=0; i<n; i++)
{
    C1max[i] = max(C1[i]).value();
    C0max[i] = max(C0[i]).value();
    C1min[i] = min(C1[i]).value();
    C0min[i] = min(C0[i]).value();

    Y1max[i] = 1;
    Y0max[i] = 1;
    Y1min[i] = 0;
    Y0min[i] = 0;
}

double T1min, T1max, H1min, H1max;
T1min = transPropDict.lookupOrDefault("T1MIN", 500);
T1max = transPropDict.lookupOrDefault("T1MAX", 800);
T1_TMP_VAL = 0.5*(T1min + T1max);

double T0min, T0max, H0min, H0max;
T0min = transPropDict.lookupOrDefault("T0MIN", 500);
T0max = transPropDict.lookupOrDefault("T0MAX", 800);
T0_TMP_VAL = 0.5*(T0min + T0max);

//Done max and min values for scalar fields
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//PLIC interface class set up and interface properties like normal, curvature 
plic interface(mesh, alpha1, alpha0, U, phi, rho1, rho0, rho, osPlic);

volScalarField alpha1_smooth
(
    IOobject
    (
        "alpha1_smooth",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::average(alpha1),
    alpha1.boundaryField().types()
);

dimensionedScalar deltaN
(
    "deltaN",
    1e-8/pow(average(mesh.V()), 1.0/3.0)
);

surfaceScalarField nHatf
(
    IOobject
    (
        "nHatf",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("nHatf", dimArea, 0.0)
);

dimensionedScalar sigma(transPropDict.lookup("sigma"));

volScalarField K
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("K", dimless/dimLength, 0.0)
);
//Done PLIC interface class set up and interface properties like normal, curvature
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Obtain references to phase centroids, phase face centres, phase face areas, face phase states
const scalarField& Af_ph1_own = interface.Af_ph1_own();
const scalarField& Af_ph0_own = interface.Af_ph0_own();
const vectorField& Cf_ph1_own = interface.Cf_ph1_own();
const vectorField& Cf_ph0_own = interface.Cf_ph0_own();

const scalarField& Af_ph1_nei = interface.Af_ph1_nei();
const scalarField& Af_ph0_nei = interface.Af_ph0_nei();
const vectorField& Cf_ph1_nei = interface.Cf_ph1_nei();
const vectorField& Cf_ph0_nei = interface.Cf_ph0_nei();

vectorField Cf_ph1(nFaces);
vectorField Cf_ph0(nFaces);

const vectorField& C_ph1 = interface.C_ph1();
const vectorField& C_ph0 = interface.C_ph0();

const volVectorField& C_intfc = interface.C_intfc();
const volScalarField& A_intfc = interface.A_intfc();
const volVectorField& nHat = interface.nHat();
//Done obtain references to phase centroids, phase face centres, phase face areas, face phase states
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Gravity terms and total pressure
#include "readGravitationalAcceleration.H"

Info<< "Calculating field g.h\n" << endl;
volScalarField gh("gh", g & mesh.C());
surfaceScalarField ghf("ghf", g & mesh.Cf());

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);
 
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    mesh.solutionDict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}
//Done gravity terms and total pressure
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//alpha, species, density advection
surfaceScalarField& phiAlpha1(interface.phiAlpha1());
volScalarField surfInt_phiAlpha1(fvc::surfaceIntegrate(phiAlpha1));
scalarField& surfInt_phiAlpha1Cells = surfInt_phiAlpha1.internalField();
//volScalarField surfInt_phiAlpha0(fvc::surfaceIntegrate((phi- phiAlpha1)));
//scalarField& surfInt_phiAlpha0Cells = surfInt_phiAlpha0.internalField();
centredCPCCellToCellStencilObject cell_stencil(mesh);

labelList unbounded_Y_in_cell(nCells, 0);

PtrList<surfaceScalarField> advFlux_Y1(n-1);
PtrList<surfaceScalarField> advFlux_Y0(n-1);
List<scalarField> dC1_redistCells(n-1);
List<scalarField> dC0_redistCells(n-1);

for(i=0; i<(n-1); i++)
{
    fieldName = "advFlux_Y1" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    advFlux_Y1.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("advFlux_Y1", dimVol*dimDensity/dimTime, 0)
        )
    );

    fieldName = "advFlux_Y0" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    advFlux_Y0.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("advFlux_Y0", dimVol*dimDensity/dimTime, 0)
        )
    );

    dC1_redistCells[i].setSize(nCells);
    dC0_redistCells[i].setSize(nCells);
}
//Done alpha, species, density advection
//------------------------------------------------------------------//


//------------------------------------------------------------------//
//Enthalpy advection
fieldName = "advFlux_rhoh1";
Info<< "Initializing field " << fieldName << endl;
surfaceScalarField advFlux_rhoh1
(
    IOobject
    (
        fieldName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("advFlux_rhoh1", dimVol*dimPressure/dimTime, 0)
);

fieldName = "advFlux_rhoh0";
Info<< "Initializing field " << fieldName << endl;
surfaceScalarField advFlux_rhoh0
(
    IOobject
    (
        fieldName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("advFlux_rhoh0", dimVol*dimPressure/dimTime, 0)
);
//Done enthalpy advection
//------------------------------------------------------------------//


//------------------------------------------------------------------//
//Species,enthalpy diffusion
centredCFCCellToCellStencilObject diffCellStencilObject(mesh);
const labelListList& diffCellStencil = diffCellStencilObject.stencil();
const mapDistribute& diffMap = diffCellStencilObject.map();
n_flatFld_diff = diffMap.constructSize();

// flat field for phase centroid position vectors
List<vector> C_ph1_flatFld_diff(n_flatFld_diff, vector::one);
List<vector> C_ph0_flatFld_diff(n_flatFld_diff, vector::one);

// flat field for scalar variables
List<scalar> alpha1_flatFld_diff(n_flatFld_diff, 0.0);
List<List<scalar> > Y1_flatFld_diff(n); 
List<List<scalar> > Y0_flatFld_diff(n);
List<scalar> Y_flatFld_diff(n_flatFld_diff, 0.0);
for(i=0; i<n; i++)
{
    Y1_flatFld_diff[i] = Y_flatFld_diff;    
    Y0_flatFld_diff[i] = Y_flatFld_diff;        
}
List<scalar> T1_flatFld_diff(n_flatFld_diff, T_0);
List<scalar> T0_flatFld_diff(n_flatFld_diff, T_0);

// labelList to indicate phase state of each mesh face
labelList face_phaseState_diff(nFaces);

// weights for interpolation of fields to faces
surfaceScalarField w_ph1
(
    IOobject
    (
        "w_ph1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("w_ph1", dimless, 1)
);

surfaceScalarField w_ph0
(
    IOobject
    (
        "w_ph0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("w_ph0", dimless, 1)
);

// fields interpolated to faces
surfaceScalarField rho1f(linearInterpolate(rho1));
surfaceScalarField rho0f(linearInterpolate(rho0));

surfaceScalarField T1f(linearInterpolate(T1));
surfaceScalarField T0f(linearInterpolate(T0));
surfaceScalarField v1f(linearInterpolate(v1));
surfaceScalarField v0f(linearInterpolate(v0));

PtrList<surfaceScalarField> x1f(n);
PtrList<surfaceScalarField> x0f(n);
PtrList<surfaceScalarField> hpar1f(n);
PtrList<surfaceScalarField> hpar0f(n);
for(i=0; i<n; i++)
{
    fieldName = "x1f" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    x1f.set
        (
            i,
        new surfaceScalarField
            (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            linearInterpolate(x1[i])            
            )
        );

    fieldName = "x0f" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    x0f.set
        (
            i,
        new surfaceScalarField
            (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            linearInterpolate(x0[i])            
            )
        );

    fieldName = "hpar1f" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    hpar1f.set
        (
            i,
        new surfaceScalarField
            (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 2, -2, 0, 0, 0, 0)
            )
        );    

    fieldName = "hpar0f" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    hpar0f.set
        (
            i,
        new surfaceScalarField
            (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 2, -2, 0, 0, 0, 0)
            )
        );
}

surfaceScalarField lambda1f
(
    IOobject
    (
        "lambda1f",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    linearInterpolate(lambda1)
);

surfaceScalarField lambda0f
(
    IOobject
    (
        "lambda0f",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    linearInterpolate(lambda0)
);

// gradient fields and diffusive flux fields
PtrList<surfaceScalarField> gradf_Y1(n);
PtrList<surfaceScalarField> gradf_Y0(n);

PtrList<surfaceScalarField> diffFlux_Y1(n);
PtrList<surfaceScalarField> diffFlux_Y0(n);

PtrList<volScalarField> diffTerm_Y1(n);
PtrList<volScalarField> diffTerm_Y0(n);

dimensionedScalar diffTerm_zero("diffTerm_zero", C1[0].dimensions(), 0.0);

for(i=0; i<n; i++)
{
    fieldName = "gradf_Y1" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    gradf_Y1.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("gradf_Y1", dimless/dimLength, 0)
        )
    );

    fieldName = "gradf_Y0" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    gradf_Y0.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("gradf_Y0", dimless/dimLength, 0)
        )
    );

    fieldName = "diffFlux_Y1" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    diffFlux_Y1.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("diffFlux_Y1", dimVol*dimDensity/dimTime, 0)
        )
    );

    fieldName = "diffTerm_Y1" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    diffTerm_Y1.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            diffTerm_zero
        )
    );

    fieldName = "diffFlux_Y0" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    diffFlux_Y0.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("diffFlux_Y0", dimVol*dimDensity/dimTime, 0)
        )
    );

    fieldName = "diffTerm_Y0" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    diffTerm_Y0.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            diffTerm_zero
        )
    );
}

//Done species diffusion
//------------------------------------------------------------------//


//------------------------------------------------------------------//
//T, enthalpy diffusion

fieldName = "gradf_T1";
Info<< "Initializing field " << fieldName << endl;
surfaceScalarField gradf_T1
(
    IOobject
    (
        fieldName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("gradf_T1", dimTemperature/dimLength, 0)
);

fieldName = "gradf_T0";
Info<< "Initializing field " << fieldName << endl;
surfaceScalarField gradf_T0
(
    IOobject
    (
        fieldName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("gradf_T0", dimTemperature/dimLength, 0)
);

fieldName = "diffFlux_h1";
Info<< "Initializing field " << fieldName << endl;
surfaceScalarField diffFlux_h1
(
    IOobject
    (
        fieldName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("diffFlux_h1", dimVol*dimPressure/dimTime, 0)
);

fieldName = "diffFlux_h0";
Info<< "Initializing field " << fieldName << endl;
surfaceScalarField diffFlux_h0
(
    IOobject
    (
        fieldName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("diffFlux_h0", dimVol*dimPressure/dimTime, 0)
);

fieldName = "diffFlux_Jh1";
Info<< "Initializing field " << fieldName << endl;
surfaceScalarField diffFlux_Jh1
(
    IOobject
    (
        fieldName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("diffFlux_Jh1", dimVol*dimPressure/dimTime, 0)
);

fieldName = "diffFlux_Jh0";
Info<< "Initializing field " << fieldName << endl;
surfaceScalarField diffFlux_Jh0
(
    IOobject
    (
        fieldName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("diffFlux_Jh0", dimVol*dimPressure/dimTime, 0)
);

volScalarField surfInt_diffFlux_h1(fvc::surfaceIntegrate(diffFlux_h1));
volScalarField surfInt_diffFlux_h0(fvc::surfaceIntegrate(diffFlux_h0));
volScalarField surfInt_diffFlux_Jh1(fvc::surfaceIntegrate(diffFlux_Jh1));
volScalarField surfInt_diffFlux_Jh0(fvc::surfaceIntegrate(diffFlux_Jh0));

scalarField& surfInt_diffFlux_h1Cells = surfInt_diffFlux_h1.internalField();
scalarField& surfInt_diffFlux_h0Cells = surfInt_diffFlux_h0.internalField();
scalarField& surfInt_diffFlux_Jh1Cells = surfInt_diffFlux_Jh1.internalField();
scalarField& surfInt_diffFlux_Jh0Cells = surfInt_diffFlux_Jh0.internalField();

//Done T, enthalpy diffusion
//------------------------------------------------------------------//


//------------------------------------------------------------------//
//Interfacial conditions and source terms
const labelListList& istCellStencil = cell_stencil.stencil();
const mapDistribute& istMap = cell_stencil.map();
n_flatFld_ist = istMap.constructSize();

labelList calc_ist_for_cell(nCells, 0);

// flat field for phase centroid position vectors
List<vector> C_ph1_flatFld_ist(n_flatFld_ist, vector::one);
List<vector> C_ph0_flatFld_ist(n_flatFld_ist, vector::one);

// flat field for scalar variables
List<scalar> alpha1_flatFld_ist(n_flatFld_ist, 0.0);
List<List<scalar> > Y1_flatFld_ist(n); 
List<List<scalar> > Y0_flatFld_ist(n);
List<scalar> Y_flatFld_ist(n_flatFld_ist, 0.0);
for(i=0; i<n; i++)
{
    Y1_flatFld_ist[i] = Y_flatFld_ist;    
    Y0_flatFld_ist[i] = Y_flatFld_ist;        
}
List<scalar> T1_flatFld_ist(n_flatFld_ist, T1min);
List<scalar> T0_flatFld_ist(n_flatFld_ist, T0min);

PtrList<volScalarField> Ys1(n);
PtrList<volScalarField> Ys0(n);
PtrList<volScalarField> Js1(n);
PtrList<volScalarField> Js0(n);
PtrList<volScalarField> grads1(n);
PtrList<volScalarField> grads0(n);
PtrList<volScalarField> mS1(n);
PtrList<volScalarField> mS0(n);

// Binary LLE calculation to set equilibrium mass fractions
n_miscible = 0;
binarylle_(&P_thermo_tmp,&T_0,&n,Pc,Tc,w,kij,xs1_tmp,xs0_tmp,&n_miscible);

plicFuncs::x2y(n,MW,xs1_tmp,ys1_tmp);
plicFuncs::x2y(n,MW,xs0_tmp,ys0_tmp);

Info<< "n_miscible = " << n_miscible << endl;
Info<< "Equilibrium mass fractions: " << endl;
Info<< "ys1:  ";
for(i=0; i<n; i++) Info<< ys1_tmp[i] << "  ";
Info<< endl;
Info<< "ys0:  ";
for(i=0; i<n; i++) Info<< ys0_tmp[i] << "  ";
Info<< endl;

if(n_miscible==1)
{    
    Info<< "FATAL ERROR: Fully miscible species pair at given T0..." << endl;
    exit(-1);
}

if(calcDiff)
{
    new_tlsm_diffusion_krishna_model_(&P_thermo_tmp,&T_0,&n,Pc,Tc,Vc,w,MW,kij,xs1_tmp,Dij1_tmp);
    new_tlsm_diffusion_krishna_model_(&P_thermo_tmp,&T_0,&n,Pc,Tc,Vc,w,MW,kij,xs0_tmp,Dij0_tmp);
}

for(i=0; i<n; i++)
{    
    fieldName = "Ys1" + Foam::name(i);    
    Info<< "Initializing field " << fieldName << endl;
    Ys1.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(fieldName+"_0", dimless, ys1_tmp[i])
        )
    );

    fieldName = "Ys0" + Foam::name(i);    
    Info<< "Initializing field " << fieldName << endl;
    Ys0.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(fieldName+"_0", dimless, ys0_tmp[i])
        )
    );

    fieldName = "Js1" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    Js1.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(fieldName, dimMass/dimTime, 0)
        )
    );

    fieldName = "Js0" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    Js0.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(fieldName, dimMass/dimTime, 0)
        )
    );

    fieldName = "grads1" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    grads1.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(fieldName, dimless/dimLength, 0)
        )
    );

    fieldName = "grads0" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    grads0.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(fieldName, dimless/dimLength, 0)
        )
    );

    fieldName = "mS1" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    mS1.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(fieldName, dimMass/dimVol/dimTime, 0)
        )
    );

    fieldName = "mS0" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    mS0.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(fieldName, dimMass/dimVol/dimTime, 0)
        )
    );
}

fieldName = "JsTot";
Info<< "Initializing field " << fieldName << endl;
volScalarField JsTot
(
    IOobject
    (
        fieldName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(fieldName, dimMass/dimTime, 0)
);

fieldName = "mS1Tot";
Info<< "Initializing field " << fieldName << endl;
volScalarField mS1Tot
(
    IOobject
    (
        fieldName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(fieldName, dimMass/dimVol/dimTime, 0)
);

fieldName = "mS0Tot";
Info<< "Initializing field " << fieldName << endl;
volScalarField mS0Tot
(
    IOobject
    (
        fieldName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(fieldName, dimMass/dimVol/dimTime, 0)
);

fieldName = "alphaS1";
Info<< "Initializing field " << fieldName << endl;
volScalarField alphaS1
(
    IOobject
    (
        fieldName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(fieldName, dimless/dimTime, 0)
);

fieldName = "alphaS0";
Info<< "Initializing field " << fieldName << endl;
volScalarField alphaS0
(
    IOobject
    (
        fieldName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(fieldName, dimless/dimTime, 0)
);

fieldName = "pS";
Info<< "Initializing field " << fieldName << endl;
volScalarField pS
(
    IOobject
    (
        fieldName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(fieldName, dimless/dimTime, 0)
);

fieldName = "delta1";
Info<< "Initializing field " << fieldName << endl;
volScalarField delta1
(
    IOobject
    (
        fieldName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(fieldName, dimLength, 0)
);

fieldName = "delta0";
Info<< "Initializing field " << fieldName << endl;
volScalarField delta0
(
    IOobject
    (
        fieldName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(fieldName, dimLength, 0)
);

fieldName = "Ts";
Info<< "Initializing field " << fieldName << endl;
volScalarField Ts
(
    IOobject
    (
        fieldName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1*T1 + alpha0*T0
);

fieldName = "Qs";
Info<< "Initializing field " << fieldName << endl;
volScalarField Qs
(
    IOobject
    (
        fieldName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(fieldName, dimPressure/dimTime, 0)
);

labelList n_iters_Ts(nCells);
labelList status_transLLE(nCells);
labelList cell_had_intfc(nCells);

for(i=0; i<n; i++)
{
    if(calcDiff)
    {
        fieldName = "D1" + Foam::name(i);
        Info<< "Initializing diffusivity " << fieldName << endl;
        D1_0.set
            (
                i,
                new dimensionedScalar
                (
                    "D1_0", dimLength*dimLength/dimTime, Dij1_tmp[1]
                )
            );

        fieldName = "D0" + Foam::name(i);
        Info<< "Initializing diffusivity " << fieldName << endl;
        D0_0.set
            (
                i,
                new dimensionedScalar
                (
                    "D0_0", dimLength*dimLength/dimTime, Dij0_tmp[1]
                )
            );
    }
    else
    {
        fieldName = "D1" + Foam::name(i);
        Info<< "Reading diffusivity " << fieldName << endl;
        D1_0.set
            (
                i,
                new dimensionedScalar
                (
                    transPropDict.lookup("D1_0")
                )
            );

        fieldName = "D0" + Foam::name(i);
        Info<< "Reading diffusivity " << fieldName << endl;
        D0_0.set
            (
                i,
                new dimensionedScalar
                (
                    transPropDict.lookup("D0_0")
                )
            );
    }
}

for(i=0; i<n; i++)
{    
    for(j=0; j<n; j++)
    {
        ijIndex = n*i + j;

        fieldName = "Dij1" + Foam::name(i) + Foam::name(j);
        Info<< "Initializing field " << fieldName << endl;
        Dij1.set
        (
            ijIndex,
            new volScalarField
            (
                IOobject
                (
                    fieldName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                D1_0[i]
            )
        );
        
        fieldName = "Dij0" + Foam::name(i) + Foam::name(j);
        Info<< "Initializing field " << fieldName << endl;
        Dij0.set
        (
            ijIndex,
            new volScalarField
            (
                IOobject
                (
                    fieldName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                D0_0[i]
            )
        );        
    }

    fieldName = "D1" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    D1.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            D1_0[i]
        )
    );

    fieldName = "D0" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    D0.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            D0_0[i]
        )
    );

    fieldName = "D1f" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    D1f.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            D1_0[i]
        )
    );

    fieldName = "D0f" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    D0f.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            D0_0[i]
        )
    );        
}
//------------------------------------------------------------------//
/*
//Read input parameters for error function inverse table from input file
f = fopen("constant/erfInv_input.dat+", "r");
iread = fscanf(f, "%lf %lf %d", &erf_a, &erf_b, &nErf);
fclose(f);

_NNEW2_(erfInv_table, double, nErf);

f = fopen("constant/erfInv_table.dat+", "r");
for(iErf=0; iErf<nErf; iErf++)
{
    iread = fscanf(f, "%lf", &erfInv_table[iErf]);    
}
fclose(f);
*/
//------------------------------------------------------------------//

//Done interfacial conditions variables and source terms
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Reference to internal fields for easy access
scalarField& alpha1Cells = alpha1.internalField();
scalarField& alpha0Cells = alpha0.internalField();
scalarField& alpha1OldCells = alpha1_old.internalField();
//scalarField& alpha0OldCells = alpha0_old.internalField();
scalarField& rho1Cells = rho1.internalField();
scalarField& rho0Cells = rho0.internalField();
scalarField& mS1TotCells = mS1Tot.internalField();
scalarField& alphaS1Cells = alphaS1.internalField();
scalarField& alphaS0Cells = alphaS0.internalField();
scalarField& pSCells = pS.internalField();
scalarField& H1Cells = H1.internalField();
scalarField& H0Cells = H0.internalField();
scalarField& H1OldCells = H1_old.internalField();
scalarField& H0OldCells = H0_old.internalField();
scalarField& rhoh1Cells = rhoh1.internalField();
scalarField& rhoh0Cells = rhoh0.internalField();
scalarField& h1Cells = h1.internalField();
scalarField& h0Cells = h0.internalField();
scalarField& T1Cells = T1.internalField();
scalarField& T0Cells = T0.internalField();
scalarField& TsCells = Ts.internalField();
scalarField& QsCells = Qs.internalField();
scalarField& lambda1Cells = lambda1.internalField();
scalarField& lambda0Cells = lambda0.internalField();
scalarField& Cp1Cells = Cp1.internalField();
scalarField& Cp0Cells = Cp0.internalField();
scalarField& Dh1Cells = Dh1.internalField();
scalarField& Dh0Cells = Dh0.internalField();
//Done reference to internal fields for easy access
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Initialize dependent variable fields correctly
plicFuncs::correct_x_from_Y(n, MW, Y1, x1);
plicFuncs::correct_x_from_Y(n, MW, Y0, x0);

//plicFuncs::correct_rho(P_thermo_tmp,T1,x1,n,Pc,Tc,w,MW,kij_T,Ta_kij,Tb_kij,nT_kij,calcKij,kij,rho1);
//plicFuncs::correct_rho(P_thermo_tmp,T0,x0,n,Pc,Tc,w,MW,kij_T,Ta_kij,Tb_kij,nT_kij,calcKij,kij,rho0);

for(i=0; i<n; i++)
{
    volScalarField& Y1i = Y1[i];
    volScalarField& Y0i = Y0[i];
    volScalarField& c1i = c1[i];
    volScalarField& c0i = c0[i];
    volScalarField& C1i = C1[i];
    volScalarField& C0i = C0[i];
    volScalarField& C1Oldi = C1_old[i];
    volScalarField& C0Oldi = C0_old[i];
    volScalarField& c1Oldi = c1_old[i];
    volScalarField& c0Oldi = c0_old[i];

    c1i = rho1*Y1i;
    C1i = alpha1*c1i;
    C1Oldi = C1i;
    c1Oldi = c1i;

    c0i = rho0*Y0i;
    C0i = alpha0*c0i;
    C0Oldi = C0i;
    c0Oldi = c0i;

    Y1i.correctBoundaryConditions();
    plicFuncs::correct_boundaryField_C(Y1i, rho1, alpha1, c1i, C1i);

    Y0i.correctBoundaryConditions();
    plicFuncs::correct_boundaryField_C(Y0i, rho0, alpha0, c0i, C0i);

    forAll(C1i.boundaryField(), patchI)
    {
        forAll(C1i.boundaryField()[patchI], fcI)
        {
            C1Oldi.boundaryField()[patchI][fcI] = C1i.boundaryField()[patchI][fcI];
            C0Oldi.boundaryField()[patchI][fcI] = C0i.boundaryField()[patchI][fcI];
            c1Oldi.boundaryField()[patchI][fcI] = c1i.boundaryField()[patchI][fcI];
            c0Oldi.boundaryField()[patchI][fcI] = c0i.boundaryField()[patchI][fcI];
        }
    }
}


    plicFuncs::correct_h(P_thermo_tmp, T1, x1, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, h1);
    rhoh1 = rho1*h1;
    H1 = alpha1*rhoh1;
    H1_old = H1;

    plicFuncs::correct_h(P_thermo_tmp, T0, x0, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, h0);
    rhoh0 = rho0*h0;
    H0 = alpha0*rhoh0;
    H0_old = H0;

plicFuncs::correct_boundaryField_h_rhoh_H(P_thermo_tmp, T1, x1, rho1, alpha1, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, h1, rhoh1, H1);

plicFuncs::correct_boundaryField_h_rhoh_H(P_thermo_tmp, T0, x0, rho0, alpha0, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, h0, rhoh0, H0);

#include "correct_thermo_trans_prop.H"

#include "diff_grad_interp.H"

rho = alpha1*rho1 + alpha0*rho0;
mu = alpha1*mu1 + alpha0*mu0;
//Done initialize dependent variable fields correctly
//------------------------------------------------------------------//

for(i=0; i<n; i++)
{
    C1max[i] = max(alpha1*rho1*Y1max[i]).value();
    C0max[i] = max(alpha0*rho0*Y0max[i]).value();
    C1min[i] = min(alpha1*rho1*Y1max[i]).value();
    C0min[i] = min(alpha0*rho0*Y0max[i]).value();
}

//------------------------------------------------------------------//
//Clipping mass conservation errors
scalar alpha1_avg;
scalar alpha1_clip_avg;
scalar alpha1_clip_err;
scalar alpha1_clip_cum_err = 0;

List<scalar> C1_avg(n);
List<scalar> C1_clip_avg(n);
List<scalar> C1_clip_err(n);
List<scalar> C1_clip_cum_err(n);

List<scalar> C0_avg(n);
List<scalar> C0_clip_avg(n);
List<scalar> C0_clip_err(n);
List<scalar> C0_clip_cum_err(n);

for(i=0; i<n; i++)
{
    C1_clip_cum_err[i] = 0;
    C0_clip_cum_err[i] = 0;
}
//------------------------------------------------------------------//

/*
volScalarField p_rghEqnComp1 =
    fvc::ddt(rho1)
    + fvc::div(phi, rho1, "div(phi,rho)") - fvc::Sp(fvc::div(phi), rho1);

volScalarField p_rghEqnComp0 =
    fvc::ddt(rho0)
    + fvc::div(phi, rho0, "div(phi,rho)") - fvc::Sp(fvc::div(phi), rho0);

volScalarField p_rghEqnComp =
    (max(alpha1, scalar(0))/rho1)*p_rghEqnComp1
    + (max(alpha0, scalar(0))/rho0)*p_rghEqnComp0;
    */
Info << "Done createFields.H" << endl;
