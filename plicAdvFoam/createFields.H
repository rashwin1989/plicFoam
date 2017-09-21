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

// initialization values for phase densities
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
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1*rho1 + alpha0*rho0,
    alpha1.boundaryField().types()
);

// Construct PLIC interface object from alpha1 distribution, U field etc. 
plic interface(mesh, alpha1, U, phi, rho1, rho0, rho);

surfaceScalarField& phiAlpha1(interface.phiAlpha1());

centredCPCCellToCellStencilObject cell_stencil(mesh);

// Y, C, species initialization
label nSpecies;
nSpecies = transPropDict.lookupOrDefault("nSpecies", 2);
Info<< "Number of species " << nSpecies << endl;

PtrList<volScalarField> Y1(nSpecies);
PtrList<volScalarField> C1(nSpecies);
PtrList<volScalarField> Y0(nSpecies);
PtrList<volScalarField> C0(nSpecies);

word fieldName;
for(label i=0; i<nSpecies; i++)
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
}

// Y advection
PtrList<surfaceScalarField> advFlux_Y1(nSpecies);
PtrList<surfaceScalarField> advFlux_Y0(nSpecies);

for(label i=0; i<nSpecies; i++)
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
                IOobject::AUTO_WRITE
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
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("advFlux_Y0", dimVol*dimDensity/dimTime, 0)
        )
    );
}

// Y diffusion
/*
PtrList<dimensionedScalar> D1(nSpecies);
PtrList<dimensionedScalar> D0(nSpecies);

word DName;
for(label i=0; i<nSpecies; i++)
{
    DName = "D1" + Foam::name(i);
    Info<< "Reading diffusivity " << DName << endl;
    D1.set
    (
        i,
        new dimensionedScalar
        (
            transPropDict.lookup(DName)
        )
    );

    DName = "D0" + Foam::name(i);
    Info<< "Reading diffusivity " << DName << endl;
    D0.set
    (
        i,
        new dimensionedScalar
        (
            transPropDict.lookup(DName)
        )
    );
}

// scalarField diffWeightsOwn_ph1(mesh.nFaces());
// scalarField diffWeightsNei_ph1(mesh.nFaces());
// scalarField diffWeightsOwn_ph0(mesh.nFaces());
// scalarField diffWeightsNei_ph0(mesh.nFaces());

surfaceScalarField diffFlux_ph1
(
    IOobject
    (
        "diffFlux_ph1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("diffFlux_ph1", dimVol/dimTime, 0)
);

surfaceScalarField diffFlux_ph0
(
    IOobject
    (
        "diffFlux_ph0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("diffFlux_ph0", dimVol/dimTime, 0)
);

centredCFCCellToCellStencilObject diffCellStencilObject(mesh);
const labelListList& diffCellStencil = diffCellStencilObject.stencil();


const labelList& own = mesh.owner();
const labelList& nei = mesh.neighbour();
//const vectorField& meshC = mesh.C();
const vectorField& meshSf = mesh.Sf();
const scalarField& meshMagSf = mesh.magSf();

const mapDistribute& diffMap = diffCellStencilObject.map();
label n_flatFld_diff = diffMap.constructSize();

List<scalar> alpha1_flatFld_diff(n_flatFld_diff, 0.0);

List<List<scalar> > Y1_flatFld_diff(nSpecies); 
List<List<scalar> > Y0_flatFld_diff(nSpecies);
List<scalar> Y_flatFld_diff(n_flatFld_diff, 0.0);     

for(label i=0; i<nSpecies; i++)
{
    Y1_flatFld_diff[i] = Y_flatFld_diff;    
    Y0_flatFld_diff[i] = Y_flatFld_diff;        
}

List<vector> C_ph1_flatFld_diff(n_flatFld_diff, vector::one);
List<vector> C_ph0_flatFld_diff(n_flatFld_diff, vector::one);
*/

// tolerances and max iters
label ALPHA_BOUND_ITERS_MAX = 108;
scalar ALPHA_BOUND_ERR = transPropDict.lookupOrDefault("ALPHA_BOUND_ERR", 1E-15);
scalar SMALLEST_ALPHA_FOR_YADV = transPropDict.lookupOrDefault("SMALLEST_ALPHA_FOR_YADV", 1E-15);

// debug bool switches
bool debug;
debug = transPropDict.lookupOrDefault("debug", false);
bool unbound_debug;
unbound_debug = transPropDict.lookupOrDefault("unbound_debug", false);
bool phi_interp_debug;
phi_interp_debug = transPropDict.lookupOrDefault("phi_interp_debug", false);
bool adv_debug;
adv_debug = transPropDict.lookupOrDefault("adv_debug", false);
bool adv_debug2;
adv_debug2 = transPropDict.lookupOrDefault("adv_debug2", false);

// 2-D bool switch
bool twoD;
twoD = runTime.controlDict().lookupOrDefault("twoD", true);
