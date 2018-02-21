scalar IFTCoNum = 0.0;
scalar meanIFTCoNum = 0.0;

scalar pi = Foam::constant::mathematical::pi;

volScalarField rhoAvg = 0.5*(rho1+rho0);

scalarField avgMagSf
(
    fvc::average(mesh.magSf())().internalField()
);

IFTCoNum = Foam::sqrt(2*pi*sigma.value()*gMax(avgMagSf*avgMagSf*avgMagSf/mesh.V().field()/mesh.V().field()/mesh.V().field()/rhoAvg.field()))*runTime.deltaTValue();

meanIFTCoNum = Foam::sqrt(2*pi*sigma.value()*pow((gSum(avgMagSf)/gSum(mesh.V().field())), 3)/gSum(rhoAvg.field()))*runTime.deltaTValue();

Info<< "IFT Courant Number mean: " << meanIFTCoNum
    << " max: " << IFTCoNum << endl;
