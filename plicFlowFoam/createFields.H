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
alpha1.correctBoundaryConditions();

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
//incompTwoPhaseMixture twoPhaseProperties(alpha1);

IOdictionary transPropDict
(
    IOobject
    (
        "transportProperties",
        U.time().constant(),
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

centredCPCCellToCellStencilObject cell_stencil(mesh);

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

// tolerances and max iters
label ALPHA_BOUND_ITERS_MAX = 108;
scalar ALPHA_BOUND_ERR = transPropDict.lookupOrDefault("ALPHA_BOUND_ERR", 1E-15);
scalar MIN_PHI_DEBUG = transPropDict.lookupOrDefault("MIN_PHI_DEBUG", 1E-12);

// debug bool switches
bool debug;
debug = transPropDict.lookupOrDefault("debug", false);
bool unbound_debug;
unbound_debug = transPropDict.lookupOrDefault("unbound_debug", false);
bool phi_interp_debug;
phi_interp_debug = transPropDict.lookupOrDefault("phi_interp_debug", false);

