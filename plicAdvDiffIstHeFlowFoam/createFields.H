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

interfaceProperties intfcProp(alpha1, U, transPropDict);

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
label nSpecies;
nSpecies = transPropDict.lookupOrDefault("nSpecies", 2);
Info<< "Number of species " << nSpecies << endl;

PtrList<volScalarField> Y1(nSpecies);
PtrList<volScalarField> C1(nSpecies);
PtrList<volScalarField> Y0(nSpecies);
PtrList<volScalarField> C0(nSpecies);

PtrList<volScalarField> C_phAvg(nSpecies);

word fieldName;
for(label i=0; i<(nSpecies-1); i++)
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
    Y0[i].correctBoundaryConditions();

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

for(label i=0; i<(nSpecies - 1); i++)
{
    Y1t += Y1[i];
    Y0t += Y0[i];
}

label endSp = nSpecies - 1;

fieldName = "Y1" + Foam::name(endSp);
Info<< "Reading field " << fieldName << endl;
Y1.set
(
    endSp,
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
Y1[endSp] = scalar(1.0) - Y1t;
Y1[endSp].correctBoundaryConditions();

fieldName = "C1" + Foam::name(endSp);
Info<< "Initializing field " << fieldName << endl;
C1.set
(
    endSp,
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
        alpha1*rho1*Y1[endSp],
        Y1[endSp].boundaryField().types()
    )
);

fieldName = "Y0" + Foam::name(endSp);
Info<< "Reading field " << fieldName << endl;
Y0.set
(
    endSp,
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
Y0[endSp] = scalar(1.0) - Y0t;
Y0[endSp].correctBoundaryConditions();

fieldName = "C0" + Foam::name(endSp);
Info<< "Initializing field " << fieldName << endl;
C0.set
(
    endSp,
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
        alpha0*rho0*Y0[endSp],
        Y0[endSp].boundaryField().types()
    )
);

fieldName = "C_phAvg" + Foam::name(endSp);
Info<< "Initializing field " << fieldName << endl;
C_phAvg.set
(
    endSp,
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
        C0[endSp] + C1[endSp],
        Y1[endSp].boundaryField().types()
    )
);


List<scalar> C1max(nSpecies);
List<scalar> C0max(nSpecies);
List<scalar> C1min(nSpecies);
List<scalar> C0min(nSpecies);
List<scalar> Y1max(nSpecies);
List<scalar> Y0max(nSpecies);
List<scalar> Y1min(nSpecies);
List<scalar> Y0min(nSpecies);

for(label i=0; i<nSpecies; i++)
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

// Y advection
PtrList<surfaceScalarField> advFlux_Y1(nSpecies - 1);
PtrList<surfaceScalarField> advFlux_Y0(nSpecies - 1);

for(label i=0; i<(nSpecies - 1); i++)
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

scalarField& rho1Cells = rho1.internalField();
scalarField& rho0Cells = rho0.internalField();
scalarField& alpha1Cells = alpha1.internalField();
scalarField& alpha0Cells = alpha0.internalField();

// Interfacial species transfer
PtrList<dimensionedScalar> Ys1_0(nSpecies);
PtrList<dimensionedScalar> Ys0_0(nSpecies);
PtrList<volScalarField> Ys1(nSpecies);
PtrList<volScalarField> Ys0(nSpecies);
PtrList<volScalarField> Js1(nSpecies);
PtrList<volScalarField> Js0(nSpecies);
PtrList<volScalarField> mS1(nSpecies);
PtrList<volScalarField> mS0(nSpecies);
List<scalar> He(nSpecies);

for(label i=0; i<nSpecies; i++)
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

    fieldName = "He" + Foam::name(i);
    Info<< "Reading Henry's Law Co-efficient " << fieldName << endl;
    He[i] = transPropDict.lookupOrDefault(fieldName, 0.1);    
}

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


const labelList& own = mesh.owner();
const labelList& nei = mesh.neighbour();
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
label ALPHA_BOUND_ITERS_MAX = 12;
label Y_BOUND_ITERS_MAX = 12;
scalar ALPHA_BOUND_ERR = transPropDict.lookupOrDefault("ALPHA_BOUND_ERR", 1E-15);
scalar SMALLEST_ALPHA_FOR_YADV = transPropDict.lookupOrDefault("SMALLEST_ALPHA_FOR_YADV", 1E-15);
scalar SMALLEST_ALPHA_FOR_YDIFF = transPropDict.lookupOrDefault("SMALLEST_ALPHA_FOR_YDIFF", 1E-15);
scalar MIN_ALPHA_DIFF = transPropDict.lookupOrDefault("MIN_ALPHA_DIFF", 1E-09);
scalar MIN_PHI_DEBUG = transPropDict.lookupOrDefault("MIN_PHI_DEBUG", 1E-12);
scalar ALPHA_2PH_MIN = transPropDict.lookupOrDefault("ALPHA_2PH_MIN", 1E-09);
scalar MIN_ALPHA_IST = transPropDict.lookupOrDefault("MIN_ALPHA_IST", 1E-09);
scalar A_INTFC_2PH_MIN_COEFF = transPropDict.lookupOrDefault("A_INTFC_2PH_MIN_COEFF", 1E-12);
scalar avgSf = (mesh.magSf()).average().value();
scalar A_INTFC_2PH_MIN = A_INTFC_2PH_MIN_COEFF*avgSf;

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

fileName outputFile("alpha_dbg");
OFstream osAlpha(mesh.time().path()/outputFile);

outputFile = "adv_dbg";
OFstream osAdv(mesh.time().path()/outputFile);

outputFile = "ist_dbg";
OFstream osIst(mesh.time().path()/outputFile);

outputFile = "diff_dbg";
OFstream osDiff(mesh.time().path()/outputFile);
