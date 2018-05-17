{    
    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("pcorr", p_rgh.dimensions(), 0.0),
        pcorrTypes
    );

    dimensionedScalar rAUf("rAUf", dimTime/rho.dimensions(), 1.0);

    adjustPhi(phi, U, pcorr);

    for(iNonOrthCorr=0; iNonOrthCorr<nNonOrthCorr; iNonOrthCorr++)
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAUf, pcorr) == fvc::div(phi)
        );

        pcorrEqn.setReference(pRefCell, pRefValue);
        Info<< "Non-orthogonal pressure corrector " << iNonOrthCorr+1 << endl;
        pcorrEqn.solve();

        if (iNonOrthCorr == nNonOrthCorr-1)
        {
            phi -= pcorrEqn.flux();
            U = fvc::reconstruct(phi);
            U.correctBoundaryConditions();
        }
    }

    #include "continuityErrs.H"
}
