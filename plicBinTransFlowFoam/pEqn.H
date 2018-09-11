{
    if(pEqn_debug)
    {
        plicFuncs::print_line(osPEqn, 100);
        osPEqn<< "Solution time = " << runTime.timeName() << nl
            << "OCorr: " << iOCorr+1 << "  PCorr: " << iPCorr+1 
            << endl;
        plicFuncs::print_line(osPEqn, 100);
        osPEqn<< endl;
    }

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
    /*
    p_rghEqnComp1 =
        fvc::ddt(rho1)
        + fvc::div(phi, rho1, "div(phi,rho)") - fvc::Sp(fvc::div(phi), rho1);

    p_rghEqnComp0 =
        fvc::ddt(rho0)
        + fvc::div(phi, rho0, "div(phi,rho)") - fvc::Sp(fvc::div(phi), rho0);

    p_rghEqnComp =
        (max(alpha1, scalar(0))/rho1)*p_rghEqnComp1
        + (max(alpha0, scalar(0))/rho0)*p_rghEqnComp0;

    p_rghEqnComp1_2 =
        fvc::ddt(rho1)
        + fvc::div(phiHbyA, rho1, "div(phi,rho)") - fvc::Sp(fvc::div(phiHbyA), rho1);

    p_rghEqnComp0_2 =
        fvc::ddt(rho0)
        + fvc::div(phiHbyA, rho0, "div(phi,rho)") - fvc::Sp(fvc::div(phiHbyA), rho0);

    p_rghEqnComp_2 =
        (max(alpha1, scalar(0))/rho1)*p_rghEqnComp1_2
        + (max(alpha0, scalar(0))/rho0)*p_rghEqnComp0_2;
        */

    volVectorField U_0 = U;
    surfaceScalarField phi_0 = phi;

    for(iNonOrthCorr=0; iNonOrthCorr<nNonOrthCorr; iNonOrthCorr++)
    {
        fvScalarMatrix p_rghEqn
        (
             fvc::div(phiHbyA) 
           - fvm::laplacian(rAUf, p_rgh)
           - pS
             //+ p_rghEqnComp             
        );

        //p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver("p_rgh"));        

        if (iNonOrthCorr == nNonOrthCorr-1)
        {
            /*
            dgdt =
            (
                pos(alpha0)*p_rghEqnComp0/rho0
              - pos(alpha1)*p_rghEqnComp1/rho1
            );
                */
            phi = phiHbyA + p_rghEqn.flux();

            U = HbyA + rAU*fvc::reconstruct((phig + p_rghEqn.flux())/rAUf);            
            U.correctBoundaryConditions();
        }

        if(pEqn_debug)
        {
            volScalarField div_phiHbyA = fvc::div(phiHbyA);
            volScalarField div_rAUf_gradp = fvc::laplacian(rAUf, p_rgh);
            volScalarField pEqn_val = div_phiHbyA - div_rAUf_gradp - pS;
            volScalarField pEqn_val_2 = fvc::div(phi) - pS;

            plicFuncs::print_line(osPEqn, 100);
            osPEqn<< setw(6) << "celli" << " " << setw(12) << "pS" << " " << setw(12) << "div(phiHbyA)" << " " << setw(12) << "div(rAUf,grad(p_rgh))" << " " << setw(12) << "pEqn" << " " << setw(12) << "pEqn2" << endl;
            plicFuncs::print_line(osPEqn, 100);

            for(celli=0; celli<nCells; celli++)
            {
                osPEqn<< setw(6) << celli << " " << setw(12) << pS[celli] << " " << setw(12) << div_phiHbyA[celli] << " " << setw(12) << div_rAUf_gradp[celli] << " " << setw(12) << pEqn_val[celli] << " " << setw(12) << pEqn_val_2[celli] << endl;
            }
            plicFuncs::print_line(osPEqn, 100);
            osPEqn<< endl;

            surfaceScalarField pFlux = p_rghEqn.flux();
            surfaceScalarField phiCorr = phig + p_rghEqn.flux();
            plicFuncs::print_line(osPEqn, 100);
            osPEqn<< setw(6) << "faceI" << " " << setw(8) << "phi_0" << " " << setw(8) << "phiHbyA" << " " << setw(8) << "pFlux" << " " << setw(8) << "phi" << " " << setw(8) << "phig" << " " << setw(8) << "rAUf" << " " << setw(8) << "phiCorr" << endl;
            plicFuncs::print_line(osPEqn, 100);

            for(faceI=0; faceI<nIntFaces; faceI++)
            {
                osPEqn<< setw(6) << faceI << " " << setw(8) << phi_0[faceI] << " " << setw(8) << phiHbyA[faceI] << " " << setw(8) << pFlux[faceI] << " " << setw(8) << phi[faceI] << " " << setw(8) << phig[faceI] << " " << setw(8) << rAUf[faceI] << " " << setw(8) << phiCorr[faceI] << endl;
            }
            plicFuncs::print_line(osPEqn, 100);
            osPEqn<< endl;
            
            volVectorField UCorr = fvc::reconstruct((phig + p_rghEqn.flux())/rAUf);
            plicFuncs::print_line(osPEqn, 100);
            osPEqn<< setw(6) << "celli" << " " << setw(12) << "U_0" << " " << setw(12) << "HbyA" << " " << setw(12) << "rAU" << " " << setw(12) << "UCorr" << setw(12) << " " << "U" << endl;
            plicFuncs::print_line(osPEqn, 100);

            for(celli=0; celli<nCells; celli++)
            {
                osPEqn<< setw(6) << celli << " " << setw(12) << U_0[celli] << " " << setw(12) << HbyA[celli] << " " << setw(12) << rAU[celli] << " " << setw(12) << UCorr[celli] << " " << setw(12) << U[celli] << endl;
            }
            plicFuncs::print_line(osPEqn, 100);
            osPEqn<< endl;
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

    /*
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
        */   
}
