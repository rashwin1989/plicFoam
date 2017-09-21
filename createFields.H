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

const dimensionedScalar& rho1_0 = transPropDict.lookup("rho1");
const dimensionedScalar& rho0_0 = transPropDict.lookup("rho0");

volScalarField rho1
(
    IOobject
    (
        "rho1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
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
        IOobject::AUTO_WRITE
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
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha0*rho0,
    alpha1.boundaryField().types()
);

// Construct PLIC interface object from alpha1 distribution, U field etc. 
plic interface(mesh, alpha1, U, phi, rho1, rho0, rho);

surfaceScalarField& phiAlpha1(interface.phiAlpha1());

centredCPCCellToCellStencilObject cell_stencil(mesh);

// tolerances and max iters
label ALPHA_BOUND_ITERS_MAX = 108;
scalar ALPHA_BOUND_ERR = transPropDict.lookupOrDefault("ALPHA_BOUND_ERR", 1E-15);

// debug bool switches
bool debug;
debug = transPropDict.lookupOrDefault("debug", false);
bool unbound_debug;
unbound_debug = transPropDict.lookupOrDefault("unbound_debug", false);
bool phi_interp_debug;
phi_interp_debug = transPropDict.lookupOrDefault("phi_interp_debug", false);

// twoD bool switch
bool twoD;
twoD = runTime.controlDict().lookupOrDefault("twoD", true);

/*const mapDistribute& map = cell_stencil.map();

List<scalar> alpha1Flat(map.constructSize());

forAll(alpha1, cellI)
{
    alpha1Flat[cellI] = alpha1[cellI];
}

forAll(alpha1.boundaryField(), patchI)
{
    const fvPatchScalarField& palpha1 = alpha1.boundaryField()[patchI];

    label nCompact =
        palpha1.patch().start()
        -alpha1.mesh().nInternalFaces()
        +alpha1.mesh().nCells();

    forAll(palpha1, i)
    {
        alpha1Flat[nCompact++] = palpha1[i];  
    }
}*/

    
