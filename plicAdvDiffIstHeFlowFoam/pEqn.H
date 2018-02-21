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

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
             fvc::div(phiHbyA) - fvm::laplacian(rAUf, p_rgh)             
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + p_rghEqn.flux();            

            U = HbyA + rAU*fvc::reconstruct((phig + p_rghEqn.flux())/rAUf);            
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

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();
    surfaceVectorField U_interp(linearInterpolate(U));
 
    if(phi_interp_debug)
    {
        Info<< endl;
        forAll(phi.internalField(), faceI)
        {
            scalar curPhi = phi.internalField()[faceI];
            vector curSf = mesh.Sf()[faceI];
            scalar curMagSf = mesh.magSf()[faceI];
            scalar curPhiU = curPhi/curMagSf;         
            scalar curU_interp = (U_interp.internalField()[faceI] & curSf)/curMagSf;

            if(mag(curPhi) > MIN_PHI_DEBUG && mag((curPhiU - curU_interp)/(curU_interp + VSMALL)) > 0.5)
            {
                Info<< "Face " << faceI << nl
                    << "Flux = " << curPhi << "  Face vel. = " << curPhiU << "  Interp. U = " << curU_interp << nl
                    << "Own U = " << U.internalField()[own[faceI]] << "  Nei U = " << U.internalField()[nei[faceI]] << endl;
            }            
        }
        Info<< endl;
    }
}
