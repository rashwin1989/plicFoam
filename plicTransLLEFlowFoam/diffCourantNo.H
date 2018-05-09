
scalar diffCoNum = 0.0;
scalar meanDiffCoNum = 0.0;

dimensionedScalar DMax = gMax(D1[0]);
for(i=0; i<D1.size(); i++)
{
    dimensionedScalar DMaxI = gMax(D1[i]);
    if(DMaxI.value() > DMax.value())
    {
        DMax = DMaxI;
    }
}
for(i=0; i<D0.size(); i++)
{
    dimensionedScalar DMaxI = gMax(D0[i]);
    if(DMaxI.value() > DMax.value())
    {
        DMax = DMaxI;
    }    
}

forAll(Dh1Cells, cellI)
{
    if(alpha1Cells[cellI] > 1E-6)
    {
        Dh1Cells[cellI] = lambda1Cells[cellI]/rho1Cells[cellI]/Cp1Cells[cellI];
    }
    else
    {
        Dh1Cells[cellI] = (lambda1_0/rho1_0/Cp1_0).value();
    }
}

forAll(Dh0Cells, cellI)
{
    if(alpha0Cells[cellI] > 1E-6)
    {
        Dh0Cells[cellI] = lambda0Cells[cellI]/rho0Cells[cellI]/Cp0Cells[cellI];
    }
    else
    {
        Dh0Cells[cellI] = (lambda0_0/rho0_0/Cp0_0).value();
    }
}

dimensionedScalar DhMax = max(gMax(Dh1), gMax(Dh0));

DMax = max(DMax, DhMax);

surfaceScalarField DSfDelta(DMax*mesh.magSf()*mesh.magSf());

scalarField sumDSfDelta
(
    fvc::surfaceSum(DSfDelta)().internalField()
);

diffCoNum = gMax(sumDSfDelta/mesh.V().field()/mesh.V().field())*runTime.deltaTValue();

meanDiffCoNum = (gSum(sumDSfDelta)/gSum(mesh.V().field())/gSum(mesh.V().field()))*runTime.deltaTValue();

Info<< "Diffusion Courant Number mean: " << meanDiffCoNum
    << " max: " << diffCoNum << endl;
