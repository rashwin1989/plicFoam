//------------------------------------------------------------------//
//Define all label, scalar, vector (OpenFOAM format) variables
//Input parameters, tolerances
label nSpecies, ALPHA_BOUND_ITERS_MAX, Y_BOUND_ITERS_MAX;
scalar ALPHA_MAX_BOUND, ALPHA_MIN_BOUND, ALPHA_BOUND_ERR, SMALLEST_ALPHA_FOR_YADV, T1_TMP_VAL, MIN_ALPHA_DIFF;
//Temp variables
label faceI;
label nCells, nFaces, n_flatFld_diff, n_flatFld_ist;
scalar deltaT, dt, maxAlpha1, minAlpha1, div_alpha1, div_phi;
scalar C1Tot_cellI, C0Tot_cellI;
vector minAlpha1CellC, maxAlpha1CellC;
fileName outputFile;
//Done define all label, scalar, vector (OpenFOAM format) variables
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Define all int, bool, double (standard C++) variables
//Input parameters, tolerances
int MAX_ITER_T;
bool alpha_debug, redist_alpha, adv_debug, adv_debug2, grad_debug, diff_debug, diff_debug2;
double T_TOL, H_TOL, T1_TMP_VAL, T0_TMP_VAL;
//Temp variables
int n, i, j, ijIndex, nCells, iOCorr, nOCorr, maxAlpha1Cell, minAlpha1Cell;
double T_cellI_tmp, h_cellI_tmp, P_thermo_tmp, P_thermo;
//Done define all int, bool, double (standard C++) variables
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Define all int, bool, double (standard C++) arrays/pointers
//Input and parameter arrays
int *tk;
double *Pc, *Tc, *Vc, *w, *MW, *coef_ab, *Tb, *SG, *H8, *k, *dm;
//Temp variable arrays
double *y_cellI_tmp, *x_cellI_tmp;
//Done define all int, bool, double (standard C++) arrays/pointers
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Define and load the OpenFOAM input dinctionaries
//to read parameters from
IOdictionary transPropDict
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh.db(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

IOdictionary solDict
(
    IOobject
    (
        "fvSolution",
        mesh.solutionDict(),
        mesh.db(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
//Done define and load the OpenFOAM input dinctionaries
//to read parameters from
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Read and set no. of species
//arrays and lists need to be initialized with this size
nSpecies = transPropDict.lookupOrDefault("nSpecies", 2);
Info<< "Number of species " << nSpecies << endl;
n = nSpecies;
//Done read and set no. of species
//arrays and lists need to be initialized with this size
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Alloacte memory for arrays
_NEW_(Pc, double, n);
_NEW_(Tc, double, n);
_NEW_(Vc, double, n);
_NEW_(w, double, n);
_NEW_(MW, double, n);
_NEW_(tk, int, n);
_NEW_(coef_ab, double, n);
_NEW_(Tb, double, n);
_NEW_(SG, double, n);
_NEW_(H8, double, n);
_NEW_(k, double, n);
_NEW_(dm, double, n);
_NEW_(y_cellI_tmp, double, n);
_NEW_(x_cellI_tmp, double, n);
//Done alloacte memory for arrays
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Set up output streams linked to output files
outputFile = "alpha_dbg";
OFstream osAlpha(runTime.path()/outputFile);
outputFile = "adv_dbg";
OFstream osAdv(runTime.path()/outputFile);
outputFile = "grad_dbg";
OFstream osGrad(runTime.path()/outputFile);
outputFile = "diff_dbg";
OFstream osDiff(runTime.path()/outputFile);
outputFile = "ist_dbg";
OFstream osIst(runTime.path()/outputFile);
//Done set up output streams linked to output files
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Set mesh related variables and fields
nCells = mesh.nCells();
nFaces = mesh.nFaces();

const labelList& own = mesh.owner();
const labelList& nei = mesh.neighbour();
const vectorField& meshCi = mesh.C().internalField();
const cellList& meshCells = mesh.cells();
const scalarField& meshV = mesh.V();
const vectorField& meshSf = mesh.Sf();
const scalarField& meshMagSf = mesh.magSf();
//Done set mesh related variables and fields
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Read and initialize fundamental variable fields
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
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    alpha1,
    alpha1.boundaryField().types()
);

Info<< "Calculating field alpha0\n" << endl;
volScalarField alpha0("alpha0", scalar(1) - alpha1);

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
//Done read and initialize fundamental variable fields
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Thermo and transport properties fields initialization
Info<< "Reading transportProperties\n" << endl;

//initialization of values for phase densities
const dimensionedScalar& rho1_0 = transPropDict.lookup("rho1");
const dimensionedScalar& rho0_0 = transPropDict.lookup("rho0");
//initialization of values for phase viscosities
const dimensionedScalar& mu1_0 = transPropDict.lookup("mu1");
const dimensionedScalar& mu0_0 = transPropDict.lookup("mu0");
//initialization of values for thermal conductivities
const dimensionedScalar& lambda1_0 = transPropDict.lookup("lambda1");
const dimensionedScalar& lambda0_0 = transPropDict.lookup("lambda0");
//initialization of values for mass diffusivities
const dimensionedScalar& D1_0 = transPropDict.lookup("D1");
const dimensionedScalar& D0_0 = transPropDict.lookup("D0");

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
    1e-3*MW1_0/rho1_0
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
    1e-3*MW0_0/rho0_0
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
surfaceScalarField muf
(
    IOobject
    (
        "muf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(mu)
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

int Nij = n*n;
PtrList<dimensionedScalar> D1_0(n);
PtrList<dimensionedScalar> D0_0(n);
PtrList<volScalarField> D1(Nij);
PtrList<volScalarField> D0(Nij);
PtrList<surfaceScalarField> D1f(Nij);
PtrList<surfaceScalarField> D0f(Nij);

for(i=0; i<n; i++)
{
    fieldName = "D1" + Foam::name(i);
    Info<< "Reading diffusivity " << fieldName << endl;
    D1_0.set
    (
        i,
        new dimensionedScalar
        (
            transPropDict.lookup(fieldName)
        )
    );

    fieldName = "D0" + Foam::name(i);
    Info<< "Reading diffusivity " << fieldName << endl;
    D0_0.set
    (
        i,
        new dimensionedScalar
        (
            transPropDict.lookup(fieldName)
        )
    );
}

for(i=0; i<n; i++)
{    
    for(j=0; j<n; j++)
    {
        ijIndex = n*i + j;

        fieldName = "D1" + Foam::name(i) + Foam::name(j);
        Info<< "Initializing field " << fieldName << endl;
        D1.set
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

        fieldName = "D1f" + Foam::name(i) + Foam::name(j);
        Info<< "Initializing field " << fieldName << endl;
        D1f.set
        (
            ijIndex,
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
        
        fieldName = "D0" + Foam::name(i) + Foam::name(j);
        Info<< "Initializing field " << fieldName << endl;
        D0.set
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

        fieldName = "D0f" + Foam::name(i) + Foam::name(j);
        Info<< "Initializing field " << fieldName << endl;
        D0f.set
        (
            ijIndex,
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
}
//Done thermo and transport properties fields initialization
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Species concentrations Y, x, C initialization
PtrList<volScalarField> Y1(n);
PtrList<volScalarField> x1(n);
PtrList<volScalarField> c1(n);
PtrList<volScalarField> C1(n);
PtrList<volScalarField> C1_old(n);
PtrList<volScalarField> Y0(n);
PtrList<volScalarField> x0(n);
PtrList<volScalarField> c0(n);
PtrList<volScalarField> C0(n);
PtrList<volScalarField> C0_old(n);

PtrList<volScalarField> C_phAvg(n);

word fieldName;
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
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
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
                IOobject::AUTO_WRITE
            ),
            C1[i],
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
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
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
                IOobject::AUTO_WRITE
            ),
            C0[i],
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

volScalarField Y1t = 0.0*Y1[0];
volScalarField Y0t = 0.0*Y0[0];

for(i=0; i<(n-1); i++)
{
    Y1t += Y1[i];
    Y0t += Y0[i];
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
        Y1[0]
    )
);
Y1[endIndex] = scalar(1.0) - Y1t;
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
        Y0[0]
    )
);
Y0[endIndex] = scalar(1.0) - Y0t;
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
//T, Enthalpy initialization
Info<< "Reading field T1\n" << endl;
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

Info<< "Initializing field h1\n" << endl;
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
    dimensionedScalar("h_0", [0 2 -2 0 0 0 0], 0),
    T1.boundaryField().types()
);

Info<< "Initializing field rhoh1\n" << endl;
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

Info<< "Initializing field H1\n" << endl;
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

Info<< "Reading field T0\n" << endl;
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

Info<< "Initializing field h0\n" << endl;
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
    dimensionedScalar("h_0", [0 2 -2 0 0 0 0], 0),
    T0.boundaryField().types()
);

Info<< "Initializing field rhoh0\n" << endl;
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

Info<< "Initializing field H0\n" << endl;
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
            dimensionedScalar("hpar_0", [0 2 -2 0 0 0 0], 0)
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
            dimensionedScalar("hpar_0", [0 2 -2 0 0 0 0], 0)
        )
    );
}
//Done T, Enthalpy initialization
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

    Y1max[i] = max(Y1[i]).value();
    Y0max[i] = max(Y0[i]).value();
    Y1min[i] = min(Y1[i]).value();
    Y0min[i] = min(Y0[i]).value();
}

double T1min, T1max, H1min, H1max;
T1min = transPropDict.lookupOrDefault("T1MIN", 300);
T1max = transPropDict.lookupOrDefault("T1MAX", 800);

double T0min, T0max, H0min, H0max;
T0min = transPropDict.lookupOrDefault("T0MIN", 300);
T0max = transPropDict.lookupOrDefault("T0MAX", 800);
//Done max and min values for scalar fields
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//PLIC interface class set up and interface properties like normal, curvature 
plic interface(mesh, alpha1, U, phi, rho1, rho0, rho);

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
    mesh.solutionDict().subDict("PIMPLE"),
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
//alpha, species advection
surfaceScalarField& phiAlpha1(interface.phiAlpha1());
centredCPCCellToCellStencilObject cell_stencil(mesh);

scalarField& alpha1Cells = alpha1.internalField();
scalarField& alpha0Cells = alpha0.internalField();
scalarField& alpha1OldCells = alpha1_old.internalField();

scalarField& rho1Cells = rho1.internalField();
scalarField& rho0Cells = rho0.internalField();

PtrList<surfaceScalarField> advFlux_Y1(n-1);
PtrList<surfaceScalarField> advFlux_Y0(n-1);

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
}
//Done alpha, species advection
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Enthalpy advection
fieldName = "advFlux_h1";
Info<< "Initializing field " << fieldName << endl;
surfaceScalarField advFlux_h1
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
    dimensionedScalar("advFlux_h1", dimVol*dimPressure/dimTime, 0)
);

fieldName = "advFlux_h0";
Info<< "Initializing field " << fieldName << endl;
surfaceScalarField advFlux_h0
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
    dimensionedScalar("advFlux_h0", dimVol*dimPressure/dimTime, 0)
);
//Done enthalpy advection
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Species, enthalpy diffusion
centredCFCCellToCellStencilObject diffCellStencilObject(mesh);
const labelListList& diffCellStencil = diffCellStencilObject.stencil();
const mapDistribute& diffMap = diffCellStencilObject.map();
n_flatFld_diff = diffMap.constructSize();

// flat field for phase centroid position vectors
List<vector> C_ph1_flatFld_diff(n_flatFld_diff, vector::one);
List<vector> C_ph0_flatFld_diff(n_flatFld_diff, vector::one);

// flat field for scalar variables
List<scalar> alpha1_flatFld_diff(n_flatFld_diff, 0.0);
List<List<scalar> > x1_flatFld_diff(n); 
List<List<scalar> > x0_flatFld_diff(n);
List<scalar> x_flatFld_diff(n_flatFld_diff, 0.0);
for(i=0; i<n; i++)
{
    x1_flatFld_diff[i] = x_flatFld_diff;    
    x0_flatFld_diff[i] = x_flatFld_diff;        
}
List<scalar> T1_flatFld_diff(n_flatFld_diff, T1min);
List<scalar> T0_flatFld_diff(n_flatFld_diff, T0min);

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
surfaceScalarField T1f(fvc::linearInterpolate(T1));
surfaceScalarField T0f(fvc::linearInterpolate(T0));
surfaceScalarField v1f(fvc::linearInterpolate(v1));
surfaceScalarField v0f(fvc::linearInterpolate(v0));

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
            dimensionedScalar("hparf_0", [0 2 -2 0 0 0 0], 0)
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
            dimensionedScalar("hparf_0", [0 2 -2 0 0 0 0], 0)
        )
    );
}

surfaceScalarField rho1f
(
    IOobject
    (
        "rho1f",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    linearInterpolate(rho1)
);

surfaceScalarField rho0f
(
    IOobject
    (
        "rho0f",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    linearInterpolate(rho0)
);

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
PtrList<surfaceScalarField> gradf_x1(n);
PtrList<surfaceScalarField> gradf_x0(n);

PtrList<surfaceScalarField> diffFlux_Y1(n);
PtrList<surfaceScalarField> diffFlux_Y0(n);

for(i=0; i<n; i++)
{
    fieldName = "gradf_x1" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    gradf_x1.set
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
            dimensionedScalar("gradf_x1", dimless/dimLength, 0)
        )
    );

    fieldName = "gradf_x0" + Foam::name(i);
    Info<< "Initializing field " << fieldName << endl;
    gradf_x0.set
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
            dimensionedScalar("gradf_x0", dimless/dimLength, 0)
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
}

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
//Done species diffusion
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Interfacial conditions and source terms
const labelListList& istCellStencil = cell_stencil.stencil();
const mapDistribute& istMap = cell_stencil.map();
n_flatFld_ist = istMap.constructSize();

// flat field for phase centroid position vectors
List<vector> C_ph1_flatFld_ist(n_flatFld_ist, vector::one);
List<vector> C_ph0_flatFld_ist(n_flatFld_ist, vector::one);

// flat field for scalar variables
List<scalar> alpha1_flatFld_ist(n_flatFld_ist, 0.0);
List<List<scalar> > x1_flatFld_ist(n); 
List<List<scalar> > x0_flatFld_ist(n);
List<scalar> x_flatFld_ist(n_flatFld_ist, 0.0);
for(i=0; i<n; i++)
{
    x1_flatFld_ist[i] = x_flatFld_ist;    
    x0_flatFld_ist[i] = x_flatFld_ist;        
}
List<scalar> T1_flatFld_ist(n_flatFld_ist, T1min);
List<scalar> T0_flatFld_ist(n_flatFld_ist, T0min);

PtrList<dimensionedScalar> xs1_0(n);
PtrList<dimensionedScalar> xs0_0(n);
PtrList<dimensionedScalar> Ys1_0(n);
PtrList<dimensionedScalar> Ys0_0(n);
PtrList<volScalarField> xs1(n);
PtrList<volScalarField> xs0(n);
PtrList<volScalarField> Ys1(n);
PtrList<volScalarField> Ys0(n);
PtrList<volScalarField> Js1(n);
PtrList<volScalarField> Js0(n);
PtrList<volScalarField> mS1(n);

for(i=0; i<n; i++)
{
    fieldName = "xs1" + Foam::name(i);
    Info<< "Reading interface equilibrium mole fraction " << fieldName << endl;
    xs1_0.set
    (
        i,
        new dimensionedScalar
        (
            transPropDict.lookup(fieldName)
        )
    );

    Info<< "Initializing field " << fieldName << endl;
    xs1.set
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
            xs1_0[i]            
        )
    );

    fieldName = "Ys1" + Foam::name(i);
    Info<< "Reading interface equilibrium mass fraction " << fieldName << endl;
    Ys1_0.set
    (
        i,
        new dimensionedScalar
        (
            transPropDict.lookup(fieldName)
        )
    );

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
            Ys1_0[i]            
        )
    );

    fieldName = "xs0" + Foam::name(i);
    Info<< "Reading interface equilibrium mole fraction " << fieldName << endl;
    xs0_0.set
    (
        i,
        new dimensionedScalar
        (
            transPropDict.lookup(fieldName)
        )
    );

    Info<< "Initializing field " << fieldName << endl;
    xs0.set
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
            xs0_0[i]            
        )
    );

    fieldName = "Ys0" + Foam::name(i);
    Info<< "Reading interface equilibrium mass fraction " << fieldName << endl;
    Ys0_0.set
    (
        i,
        new dimensionedScalar
        (
            transPropDict.lookup(fieldName)
        )
    );

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
            Ys0_0[i]            
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
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(fieldName, dimMass/dimVol/dimTime, 0)
        )
    );
}

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
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(fieldName, dimMass/dimVol/dimTime, 0)
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
        IOobject::NO_WRITE
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
        IOobject::NO_WRITE
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
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(fieldName, dimless/dimTime, 0)
);

//Done interfacial conditions variables and source terms
//------------------------------------------------------------------//

// tolerances and max iters
ALPHA_BOUND_ITERS_MAX = 12;
Y_BOUND_ITERS_MAX = 12;
scalar ALPHA_BOUND_ERR = transPropDict.lookupOrDefault("ALPHA_BOUND_ERR", 1E-15);
scalar SMALLEST_ALPHA_FOR_YADV = transPropDict.lookupOrDefault("SMALLEST_ALPHA_FOR_YADV", 1E-15);
scalar SMALLEST_ALPHA_FOR_YDIFF = transPropDict.lookupOrDefault("SMALLEST_ALPHA_FOR_YDIFF", 1E-15);
scalar MIN_ALPHA_DIFF = transPropDict.lookupOrDefault("MIN_ALPHA_DIFF", 1E-09);
scalar MIN_PHI_DEBUG = transPropDict.lookupOrDefault("MIN_PHI_DEBUG", 1E-12);
scalar ALPHA_2PH_MIN = transPropDict.lookupOrDefault("ALPHA_2PH_MIN", 1E-09);

// debug bool switches
bool redist_alpha;
redist_alpha = transPropDict.lookupOrDefault("redist_alpha", false);
bool phi_interp_debug;
phi_interp_debug = transPropDict.lookupOrDefault("phi_interp_debug", false);
bool alpha_debug;
alpha_debug = transPropDict.lookupOrDefault("alpha_debug", false);
/*
bool adv_debug;
adv_debug = transPropDict.lookupOrDefault("adv_debug", false);
*/
bool adv_debug2;
adv_debug2 = transPropDict.lookupOrDefault("adv_debug2", false);

bool ist_debug;
ist_debug = transPropDict.lookupOrDefault("ist_debug", false);

bool ist_debug2;
ist_debug2 = transPropDict.lookupOrDefault("ist_debug2", false);

bool diff_debug;
diff_debug = transPropDict.lookupOrDefault("diff_debug", false);

bool diff_debug2;
diff_debug2 = transPropDict.lookupOrDefault("diff_debug2", false);
