
scalar diffCoNum = 0.0;
scalar meanDiffCoNum = 0.0;

dimensionedScalar DMax(D1[0]);
forAll(D1,i)
{
    if(D1[i].value() > DMax.value())
    {
        DMax = D1[i];
    }
}
forAll(D0,i)
{
    if(D0[i].value() > DMax.value())
    {
        DMax = D0[i];
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
