int iOCorr, nOCorr, i, nCells, maxAlpha1Cell, minAlpha1Cell, n;
label faceI, ALPHA_BOUND_ITERS_MAX, Y_BOUND_ITERS_MAX, nSpecies;
bool alpha_debug, redist_alpha, adv_debug, adv_debug2;
scalar deltaT, dt, maxAlpha1, minAlpha1, ALPHA_MAX_BOUND, ALPHA_MIN_BOUND, ALPHA_BOUND_ERR, div_alpha1, div_phi, SMALLEST_ALPHA_FOR_YADV, T1_TMP_VAL;
scalar C1Tot_cellI, C0Tot_cellI;
vector minAlpha1CellC, maxAlpha1CellC;
fileName outputFile;

double T_cellI_tmp, h_cellI_tmp;
double *Pc, *Tc, *w, *MW, *Tb, *SG, *H8;
double *y_cellI_tmp, *x_cellI_tmp;

_NEW_(Pc, double, n);
_NEW_(Tc, double, n);
_NEW_(w, double, n);
_NEW_(MW, double, n);
_NEW_(Tb, double, n);
_NEW_(SG, double, n);
_NEW_(H8, double, n);
_NEW_(y_cellI_tmp, double, n);
_NEW_(x_cellI_tmp, double, n);

outputFile = "alpha_dbg";
OFstream osAlpha(mesh.time().path()/outputFile);
outputFile = "adv_dbg";
OFstream osAdv(mesh.time().path()/outputFile);
outputFile = "ist_dbg";
OFstream osIst(mesh.time().path()/outputFile);
outputFile = "diff_dbg";
OFstream osDiff(mesh.time().path()/outputFile);

nCells = mesh.nCells();

const vectorField& meshCi = mesh.C().internalField();
const cellList& meshCells = mesh.cells();
const scalarField& meshV = mesh.V();

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

Info<< "Reading transportProperties\n" << endl;

IOdictionary transPropDict
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        U.db(),
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
        U.db(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// initialization values for phase densities
const dimensionedScalar& rho1_0 = transPropDict.lookup("rho1");
const dimensionedScalar& rho0_0 = transPropDict.lookup("rho0");

const dimensionedScalar& mu1_0 = transPropDict.lookup("mu1");
const dimensionedScalar& mu0_0 = transPropDict.lookup("mu0");

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

// Construct PLIC interface object from alpha1 distribution, U field etc. 
plic interface(mesh, alpha1, U, phi, rho1, rho0, rho);

surfaceScalarField& phiAlpha1(interface.phiAlpha1());

scalarField& alpha1Cells = alpha1.internalField();
scalarField& alpha1OldCells = alpha1_old.internalField();

centredCPCCellToCellStencilObject cell_stencil(mesh);

//interfaceProperties intfcProp(alpha1, U, transPropDict);

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

// gravity and pressure
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

// Y, C, species initialization
nSpecies = transPropDict.lookupOrDefault("nSpecies", 2);
Info<< "Number of species " << nSpecies << endl;
n = nSpecies;

PtrList<volScalarField> Y1(nSpecies);
PtrList<volScalarField> x1(nSpecies);
PtrList<volScalarField> c1(nSpecies);
PtrList<volScalarField> C1(nSpecies);
PtrList<volScalarField> C1_old(nSpecies);
PtrList<volScalarField> Y0(nSpecies);
PtrList<volScalarField> x0(nSpecies);
PtrList<volScalarField> c0(nSpecies);
PtrList<volScalarField> C0(nSpecies);
PtrList<volScalarField> C0_old(nSpecies);

PtrList<volScalarField> C_phAvg(nSpecies);

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
                mesh.time().timeName(),
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
                mesh.time().timeName(),
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
                mesh.time().timeName(),
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
                mesh.time().timeName(),
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
                mesh.time().timeName(),
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
                mesh.time().timeName(),
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
                mesh.time().timeName(),
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
                mesh.time().timeName(),
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
                mesh.time().timeName(),
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
                mesh.time().timeName(),
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
                mesh.time().timeName(),
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
            mesh.time().timeName(),
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
            mesh.time().timeName(),
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
            mesh.time().timeName(),
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
            mesh.time().timeName(),
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
            mesh.time().timeName(),
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
            mesh.time().timeName(),
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
            mesh.time().timeName(),
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
            mesh.time().timeName(),
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
            mesh.time().timeName(),
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
            mesh.time().timeName(),
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
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        C0[endIndex] + C1[endIndex],
        Y1[endIndex].boundaryField().types()
    )
);

// Enthalpy initialization
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

// Y advection
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

// Enthalpy advection

scalarField& rho1Cells = rho1.internalField();
scalarField& rho0Cells = rho0.internalField();
scalarField& alpha1Cells = alpha1.internalField();
scalarField& alpha0Cells = alpha0.internalField();

// Interfacial species transfer
PtrList<dimensionedScalar> Ys1_0(n);
PtrList<dimensionedScalar> Ys0_0(n);
PtrList<volScalarField> Ys1(n);
PtrList<volScalarField> Ys0(n);
PtrList<volScalarField> Js1(n);
PtrList<volScalarField> Js0(n);
PtrList<volScalarField> mS1(n);
PtrList<volScalarField> mS0(n);

for(i=0; i<n; i++)
{
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
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            Ys1_0[i]            
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
                mesh.time().timeName(),
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
                mesh.time().timeName(),
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
                mesh.time().timeName(),
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
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(fieldName, dimMass/dimVol, 0)
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
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(fieldName, dimMass/dimVol, 0)
        )
    );
}

fieldName = "mS1Tot";
Info<< "Initializing field " << fieldName << endl;
volScalarField mS1Tot
(
    IOobject
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(fieldName, dimMass/dimVol, 0)
);

fieldName = "mS0Tot";
Info<< "Initializing field " << fieldName << endl;
volScalarField mS0Tot
(
    IOobject
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(fieldName, dimMass/dimVol, 0)
);

fieldName = "alphaS1";
Info<< "Initializing field " << fieldName << endl;
volScalarField alphaS1
(
    IOobject
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(fieldName, dimless, 0)
);

fieldName = "alphaS0";
Info<< "Initializing field " << fieldName << endl;
volScalarField alphaS0
(
    IOobject
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(fieldName, dimless, 0)
);

fieldName = "pS";
Info<< "Initializing field " << fieldName << endl;
volScalarField pS
(
    IOobject
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(fieldName, dimless/dimTime, 0)
);

const labelListList& istCellStencil = cell_stencil.stencil();
const mapDistribute& istMap = cell_stencil.map();
label n_flatFld_ist = istMap.constructSize();

List<scalar> alpha1_flatFld_ist(n_flatFld_ist, 0.0);

List<List<scalar> > Y1_flatFld_ist(nSpecies); 
List<List<scalar> > Y0_flatFld_ist(nSpecies);
List<scalar> Y_flatFld_ist(n_flatFld_ist, 0.0);     

for(label i=0; i<nSpecies; i++)
{
    Y1_flatFld_ist[i] = Y_flatFld_ist;    
    Y0_flatFld_ist[i] = Y_flatFld_ist;        
}

// flat field for phase centroid position vectors
List<vector> C_ph1_flatFld_ist(n_flatFld_ist, vector::one);
List<vector> C_ph0_flatFld_ist(n_flatFld_ist, vector::one);


// Y diffusion

surfaceScalarField rho1f
(
    IOobject
    (
        "rho1f",
        mesh.time().timeName(),
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
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    linearInterpolate(rho0)
);

PtrList<dimensionedScalar> D1_0(nSpecies);
PtrList<dimensionedScalar> D0_0(nSpecies);

PtrList<volScalarField> D1(nSpecies);
PtrList<volScalarField> D0(nSpecies);

PtrList<surfaceScalarField> D1f(nSpecies);
PtrList<surfaceScalarField> D0f(nSpecies);

for(label i=0; i<nSpecies; i++)
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

    Info<< "Initializing field " << fieldName << endl;
    D1.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                fieldName,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            D1_0[i]            
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
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            linearInterpolate(D1[i])
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

    Info<< "Initializing field " << fieldName << endl;
    D0.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                fieldName,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            D0_0[i]            
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
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            linearInterpolate(D0[i])
        )
    );
}

PtrList<surfaceScalarField> gradf_Y1(nSpecies);
PtrList<surfaceScalarField> gradf_Y0(nSpecies);

PtrList<surfaceScalarField> diffFlux_Y1(nSpecies);
PtrList<surfaceScalarField> diffFlux_Y0(nSpecies);

for(label i=0; i<nSpecies; i++)
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

centredCFCCellToCellStencilObject diffCellStencilObject(mesh);
const labelListList& diffCellStencil = diffCellStencilObject.stencil();
const mapDistribute& diffMap = diffCellStencilObject.map();
label n_flatFld_diff = diffMap.constructSize();

label nFaces = mesh.nFaces();
labelList face_phaseState_diff(nFaces);

List<scalar> alpha1_flatFld_diff(n_flatFld_diff, 0.0);

List<List<scalar> > Y1_flatFld_diff(nSpecies); 
List<List<scalar> > Y0_flatFld_diff(nSpecies);
List<scalar> Y_flatFld_diff(n_flatFld_diff, 0.0);     

for(label i=0; i<nSpecies; i++)
{
    Y1_flatFld_diff[i] = Y_flatFld_diff;    
    Y0_flatFld_diff[i] = Y_flatFld_diff;        
}

// flat field for phase centroid position vectors
List<vector> C_ph1_flatFld_diff(n_flatFld_diff, vector::one);
List<vector> C_ph0_flatFld_diff(n_flatFld_diff, vector::one);


//const labelList& own = mesh.owner();
//const labelList& nei = mesh.neighbour();
//const vectorField& meshC = mesh.C();
//const vectorField& meshSf = mesh.Sf();
//const scalarField& meshMagSf = mesh.magSf();

//------------------------------------------------------------------//
//Obtain references to phase centroids, phase face centres, phase face areas, face phase states
const scalarField& Af_ph1_own = interface.Af_ph1_own();
const scalarField& Af_ph0_own = interface.Af_ph0_own();
const vectorField& Cf_ph1_own = interface.Cf_ph1_own();
const vectorField& Cf_ph0_own = interface.Cf_ph0_own();
//const labelList& face_phaseState_own = interface.face_phaseState_own();
const scalarField& Af_ph1_nei = interface.Af_ph1_nei();
const scalarField& Af_ph0_nei = interface.Af_ph0_nei();
const vectorField& Cf_ph1_nei = interface.Cf_ph1_nei();
const vectorField& Cf_ph0_nei = interface.Cf_ph0_nei();
//const labelList& face_phaseState_nei = interface.face_phaseState_nei();

vectorField Cf_ph1(nFaces);
vectorField Cf_ph0(nFaces);

const vectorField& C_ph1 = interface.C_ph1();
const vectorField& C_ph0 = interface.C_ph0();

//Done obtain references to phase centroids, phase face centres, phase face areas, face phase states
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
