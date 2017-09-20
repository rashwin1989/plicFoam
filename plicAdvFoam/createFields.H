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

const dimensionedScalar& rho1 = transPropDict.lookup("rho1");
const dimensionedScalar& rho0 = transPropDict.lookup("rho0");

//bool debug = false;
//debug = transPropDict.lookupOrDefault("debug", false);

//bool twoD = true;
//twoD = runTime.controlDict().lookupOrDefault("twoD", true);

//bool unbound_debug = false;
//unbound_debug = transPropDict.lookupOrDefault("unbound_debug", false);
//bool phi_interp_debug = false;
//phi_interp_debug = transPropDict.lookupOrDefault("phi_interp_debug", false);

//label ALPHA_BOUND_ITERS_MAX = 108;

//scalar ALPHA_BOUND_ERR = transPropDict.lookupOrDefault("ALPHA_BOUND_ERR", 1E-15);

bool diff_debug = false;
diff_debug = runTime.controlDict().lookupOrDefault("diff_debug", false);

bool diff_debug2 = false;
diff_debug2 = runTime.controlDict().lookupOrDefault("diff_debug2", false);

label nSpecies = 2;
nSpecies = transPropDict.lookupOrDefault("nSpecies", 2);
Info<< "Number of species " << nSpecies << endl;

PtrList<volScalarField> Y1(nSpecies);
PtrList<volScalarField> C1(nSpecies);
PtrList<volScalarField> Y0(nSpecies);
PtrList<volScalarField> C0(nSpecies);
PtrList<dimensionedScalar> D1(nSpecies);
PtrList<dimensionedScalar> D0(nSpecies);

for(label i=0; i<nSpecies; i++)
{
    word fieldName = "Y1" + Foam::name(i);
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
}
for(label i=0; i<nSpecies; i++)
{
    word fieldName = "C1" + Foam::name(i);
    Info<< "Reading field " << fieldName << endl;

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
            alpha1*Y1[i],
            Y1[i].boundaryField().types()
        )
    );
}

for(label i=0; i<nSpecies; i++)
{
    word fieldName = "Y0" + Foam::name(i);
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
}
for(label i=0; i<nSpecies; i++)
{
    word fieldName = "C0" + Foam::name(i);
    Info<< "Reading field " << fieldName << endl;

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
            alpha0*Y0[i],
            Y0[i].boundaryField().types()
        )
    );
}

for(label i=0; i<nSpecies; i++)
{
    word DName = "D1" + Foam::name(i);
    Info<< "Reading diffusivity " << DName << endl;

    D1.set
    (
        i,
        new dimensionedScalar
        (
            transPropDict.lookup(DName)
        )
    );
}

for(label i=0; i<nSpecies; i++)
{
    word DName = "D0" + Foam::name(i);
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
rho.oldTime();

// Construct PLIC interface object from alpha1 distribution, U field etc. 
plic interface(mesh, alpha1, U, phi, rho);

//surfaceScalarField& phiAlpha1(interface.phiAlpha1());

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
