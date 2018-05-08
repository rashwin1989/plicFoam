{    
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn.H();   

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(HbyA) & mesh.Sf())
      + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
    );
    
    adjustPhi(phiHbyA, U, p_rgh);

    surfaceScalarField phig
    (
        (
            fvc::interpolate(sigma*K)*fvc::snGrad(alpha1)
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh.boundaryField(),
        (
            phiHbyA.boundaryField()
          - (mesh.Sf().boundaryField() & U.boundaryField())
        )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
    );

    tmp<fvScalarMatrix> p_rghEqnComp1;
    tmp<fvScalarMatrix> p_rghEqnComp0;

    p_rghEqnComp1 =
        fvc::ddt(rho1)
        + fvc::div(phi, rho1) - fvc::Sp(fvc::div(phi), rho1);

    p_rghEqnComp0 =
        fvc::ddt(rho0)
        + fvc::div(phi, rho0) - fvc::Sp(fvc::div(phi), rho0);

    for(iNonOrthCorr=0; iNonOrthCorr<nNonOrthCorr; iNonOrthCorr++)
    {
        fvScalarMatrix p_rghEqnIncomp
        (
             fvc::div(phiHbyA) 
           - fvm::laplacian(rAUf, p_rgh)
           - pS
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        solve
        (
            (
                (max(alpha1, scalar(0))/rho1)*p_rghEqnComp1()
              + (max(alpha0, scalar(0))/rho0)*p_rghEqnComp0()
            )
          + p_rghEqnIncomp,
            mesh.solver("p_rgh"))
        );

        if (iNonOrthCorr == (nNonOrthCorr-1))
        {
            phi = phiHbyA + p_rghEqnIncomp.flux();

            U = HbyA + rAU*fvc::reconstruct((phig + p_rghEqnIncomp.flux())/rAUf);            
            U.correctBoundaryConditions();
        }
    }

    #include "continuityErrs.H"

    Info<< "Min(U) = " << min(U).value()
        << "  Max(U) = " << max(U).value()
        << endl;
    
    Info<< "Min(p_rgh) = " << min(p_rgh).value()
        << "  Max(p_rgh) = " << max(p_rgh).value()
        << endl;

    p == p_rgh + rho*gh;

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
}
