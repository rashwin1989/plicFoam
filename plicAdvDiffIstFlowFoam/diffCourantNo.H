
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

surfaceScalarField DSfDelta(DMax*mesh.magSf()*mesh.magSf());

scalarField sumDSfDelta
(
    fvc::surfaceSum(DSfDelta)().internalField()
);

diffCoNum = gMax(sumDSfDelta/mesh.V().field()/mesh.V().field())*runTime.deltaTValue();

meanDiffCoNum = (gSum(sumDSfDelta)/gSum(mesh.V().field())/gSum(mesh.V().field()))*runTime.deltaTValue();

Info<< "Diffusion Courant Number mean: " << meanDiffCoNum
    << " max: " << diffCoNum << endl;
