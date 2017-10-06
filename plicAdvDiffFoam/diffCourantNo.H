
scalar diffCoNum = 0.0;
scalar meanDiffCoNum = 0.0;

dimensionedScalar DMax = gMax(D1[0]);
forAll(D1,i)
{
    dimensionedScalar DMaxI = gMax(D1[i]);
    if(DMaxI.value() > DMax.value())
    {
        DMax = DMaxI;
    }
}
forAll(D0,i)
{
    dimensionedScalar DMaxI = gMax(D0[i]);
    if(DMaxI.value() > DMax.value())
    {
        DMax = DMaxI;
    }    
}

surfaceScalarField DSfDelta(DMax*mesh.magSf()*mesh.magSf());

scalarField sumDSfDelta
(
    fvc::surfaceSum(DSfDelta)().internalField()
);

diffCoNum = 0.5*gMax(sumDSfDelta/mesh.V().field()/mesh.V().field())*runTime.deltaTValue();

meanDiffCoNum = 0.5*(gSum(sumDSfDelta)/gSum(mesh.V().field()))*runTime.deltaTValue();

Info<< "Diffusion Courant Number mean: " << meanDiffCoNum
    << " max: " << diffCoNum << endl;
