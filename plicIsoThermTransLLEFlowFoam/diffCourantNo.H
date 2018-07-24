
scalar diffCoNum = 0.0;
scalar meanDiffCoNum = 0.0;

dimensionedScalar DMax = gMax(Dij1[1]);
for(i=0; i<n; i++)
{
    for(j=i+1; j<n; j++)
    {
        ijIndex = i*n + j;
        dimensionedScalar DMaxI = gMax(Dij1[ijIndex]);
        if(DMaxI.value() > DMax.value())
        {
            DMax = DMaxI;
        }
    }
}

for(i=0; i<n; i++)
{
    for(j=i+1; j<n; j++)
    {
        ijIndex = i*n + j;
        dimensionedScalar DMaxI = gMax(Dij0[ijIndex]);
        if(DMaxI.value() > DMax.value())
        {
            DMax = DMaxI;
        }
    }
}


surfaceScalarField DSfSf(DMax*mesh.magSf()*mesh.magSf());

scalarField sumDSfSf
(
    fvc::surfaceSum(DSfSf)().internalField()
);

diffCoNum = gMax(sumDSfSf/mesh.V().field()/mesh.V().field())*runTime.deltaTValue();

meanDiffCoNum = gSum(sumDSfSf/mesh.V().field())/gSum(mesh.V().field())*runTime.deltaTValue();

Info<< "Diffusion Courant Number mean: " << meanDiffCoNum
    << " max: " << diffCoNum << endl;
