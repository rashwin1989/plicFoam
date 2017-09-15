
nSpecies = Y1.size();

volScalarField Y1t = 0.0*Y1[0];
volScalarField Y0t = 0.0*Y0[0];

for(label i=0; i<(nSpecies - 1); i++)
{

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();
    const vectorField& meshC = mesh.C();
    volScalarField& Y1i = Y1[i];
    volScalarField& Y0i = Y0[i];
    scalar& D1i = D1[i];
    scalar& D0i = D0[i];
    const scalarField& magSf_ph1 = intfc.Af_ph1();
    const scalarField& magSf_ph0 = intfc.Af_ph0();

    scalar MIN_ALPHA_DIFF = 1E-06;
    scalar MAX_ALPHA_DIFF = 1 - MIN_ALPHA_DIFF;

    // Calculate own and nei cell weights for diffusion fluxes at faces

    // Internal faces
    const scalarField& Y1iCells = Y1i.internalField();
    const scalarField& Y0iCells = Y0i.internalField();
    const scalarField& alpha1Cells = alpha1.internalField();

    for(label faceI=0; faceI<mesh.nInternalFaces(); faceI++)
    {
        scalar curMagSf = mesh.magSf()[faceI];
        vector nf = mesh.Sf()[faceI]/curMagSf;
        scalar curMagSf_ph1 = magSf_ph1[faceI];
        scalar curAlpha1f = curMagSf_ph1/curMagSf;
        const label& faceOwn = own[faceI];
        const label& faceNei = nei[faceI];
        label donor = own[faceI];

        scalar wOwn = 1;
        scalar wNei = 1;

        if(curAlpha1f < MIN_ALPHA_DIFF)
        {
            Foam::plicFuncs::calcDiffWeights(mesh,intfc,faceI,Y2i,alpha1,diffCellStencil,0.1*MIN_ALPHA_DIFF,0,wOwn,wNei);
            diffWeightsOwn_ph0[faceI] = wOwn*D0i*magSf_ph0[faceI];
            diffWeightsNei_ph0[faceI] = wNei*D0i*magSf_ph0[faceI];    
            diffFlux_ph0[faceI] = diffWeightsOwn_ph0[faceI]*Y0iCells[faceOwn] + diffWeightsNei_ph0[faceI]*Y0iCells[faceNei];

            diffWeightsOwn_ph1[faceI] = 0;
            diffWeightsNei_ph1[faceI] = 0;        
            diffFlux_ph1[faceI] = 0;
        }
        else if(curAlpha1f > MAX_ALPHA_DIFF)
        {
            Foam::plicFuncs::calcDiffWeights(mesh,intfc,faceI,Y1i,alpha1,diffCellStencil,0.1*MIN_ALPHA_DIFF,1,wOwn,wNei);
            diffWeightsOwn_ph1[faceI] = wOwn*D1i*magSf_ph1[faceI];
            diffWeightsNei_ph1[faceI] = wNei*D1i*magSf_ph1[faceI];
            diffFlux_ph1[faceI] = diffWeightsOwn_ph1[faceI]*Y1iCells[faceOwn] + diffWeightsNei_ph1[faceI]*Y1iCells[faceNei];
        
            diffWeightsOwn_ph0[faceI] = 0;
            diffWeightsNei_ph0[faceI] = 0;        
            diffFlux_ph0[faceI] = 0;
        }
        else
        {
            Foam::plicFuncs::calcDiffWeights(mesh,intfc,faceI,Y1i,alpha1,diffCellStencil,0.1*MIN_ALPHA_DIFF,1,wOwn,wNei);
            diffWeightsOwn_ph1[faceI] = wOwn*D1i*magSf_ph1[faceI];
            diffWeightsNei_ph1[faceI] = wNei*D1i*magSf_ph1[faceI];
            diffFlux_ph1[faceI] = diffWeightsOwn_ph1[faceI]*Y1iCells[faceOwn] + diffWeightsNei_ph1[faceI]*Y1iCells[faceNei];

            Foam::plicFuncs::calcDiffWeights(mesh,intfc,faceI,Y2i,alpha1,diffCellStencil,0.1*MIN_ALPHA_DIFF,0,wOwn,wNei);
            diffWeightsOwn_ph0[faceI] = wOwn*D0i*magSf_ph0[faceI];
            diffWeightsNei_ph0[faceI] = wNei*D0i*magSf_ph0[faceI];    
            diffFlux_ph0[faceI] = diffWeightsOwn_ph0[faceI]*Y0iCells[faceOwn] + diffWeightsNei_ph0[faceI]*Y0iCells[faceNei];
        }

        if(diffFlux_ph1[faceI] > 0)
        {
            donor = faceOwn;
        }
        else
        {
            donor = faceNei;
        }
        scalar maxDiffFlux = Y1iCells[donor]*alpha1Cells[donor]*mesh.V()[donor]/runTime.deltaT().value();
        if(mag(diffFlux_ph1[faceI]) > maxDiffFlux)
        {
            if(diff_debug)
            {
                Info<< "Face " << faceI << "  diff flux ph 1: " << diffFlux_ph1[faceI] << nl
                    << "Own " << faceOwn << "  Nei " << faceNei << "  Donor " << donor << nl
                    << "MaxDiffFlux: " << maxDiffFlux << nl 
                    << "alpha1: " << alpha1Cells[donor] << "  Y1: " <<  Y1iCells[donor] << "  deltaT: " << runTime.deltaT().value() << "  V: " <<  mesh.V()[donor] << nl 
                    << endl;
            }
        }

        if(diffFlux_ph0[faceI] > 0)
        {
            donor = faceOwn;
        }
        else
        {
            donor = faceNei;
        }
        maxDiffFlux = Y0iCells[donor]*(1 - alpha1Cells[donor])*mesh.V()[donor]/runTime.deltaT().value();
        if(mag(diffFlux_ph0[faceI]) > maxDiffFlux)
        {
            if(diff_debug)
            {
                Info<< "Face " << faceI << "  diff flux ph 0: " << diffFlux_ph0[faceI] << nl
                    << "Own " << faceOwn << "  Nei " << faceNei << "  Donor " << donor << nl
                    << "MaxDiffFlux: " << maxDiffFlux << nl 
                    << "alpha0: " << (1 - alpha1Cells[donor]) << "  Y0: " <<  Y0iCells[donor] << "  deltaT: " << runTime.deltaT().value() << "  V: " <<  mesh.V()[donor] << nl 
                    << endl;
            }
        }
    } 

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches,patchI)
    {
        const polyPatch& pp = patches[patchI];
        const vectorField& pCf = mesh.Cf().boundaryField()[patchI];
        const fvPatchScalarField& pY1i = Y1i.boundaryField()[patchI];
        const fvPatchScalarField& pY0i = Y0i.boundaryField()[patchI];
        label faceI = pp.start();

        if(isA<zeroGradientFvPatchScalarField>(pY1i))
        {
            forAll(pY1i, fcI)
            {
                diffWeightsOwn_ph1[faceI] = 0;
                diffWeightsNei_ph1[faceI] = 0;        
                diffFlux_ph1[faceI] = 0;

                diffWeightsOwn_ph0[faceI] = 0;
                diffWeightsNei_ph0[faceI] = 0;        
                diffFlux_ph0[faceI] = 0;

                faceI++;
            }
        }
        else if(isA<fixedValueFvPatchScalarField>(pY1i))
        {
            forAll(pY1i, fcI)
            {
                scalar curMagSf = mesh.magSf()[faceI];
                vector nf = mesh.Sf()[faceI]/curMagSf;
                scalar curMagSf_ph1 = magSf_ph1[faceI];
                scalar curAlpha1f = curMagSf_ph1/curMagSf;
                vector curCf = pCf[fcI];
                const label& faceOwn = own[faceI];            

                if(curAlpha1f < MIN_ALPHA_DIFF)
                {                
                    vector COwn = intfc.C_ph0()[faceOwn];
                    scalar curDelta = nf & (curCf - COwn);
                    diffFlux_ph0[faceI] = D0i*magSf_ph0[faceI]*(Y0iCells[faceOwn] - pY0i[fcI])/curDelta;

                    diffFlux_ph1[faceI] = 0;
                }
                else if(curAlpha1f > MAX_ALPHA_DIFF)
                {
                    vector COwn = intfc.C_ph0()[faceOwn];
                    scalar curDelta = nf & (curCf - COwn);
                    diffFlux_ph1[faceI] = D1i*magSf_ph1[faceI]*(Y1iCells[faceOwn] - pY1i[fcI])/curDelta;

                    diffFlux_ph0[faceI] = 0;
                }
                else
                {
                    vector COwn = intfc.C_ph0()[faceOwn];
                    scalar curDelta = nf & (curCf - COwn);
                    diffFlux_ph1[faceI] = D1i*magSf_ph1[faceI]*(Y1iCells[faceOwn] - pY1i[fcI])/curDelta;

                    vector COwn = intfc.C_ph0()[faceOwn];
                    scalar curDelta = nf & (curCf - COwn);
                    diffFlux_ph0[faceI] = D0i*magSf_ph0[faceI]*(Y0iCells[faceOwn] - pY0i[fcI])/curDelta;
                }

                faceI++;
            }
        }
        else
        {
            forAll(pY1i, fcI)
            {
                diffWeightsOwn_ph1[faceI] = 0;
                diffWeightsNei_ph1[faceI] = 0;        
                diffFlux_ph1[faceI] = 0;

                diffWeightsOwn_ph0[faceI] = 0;
                diffWeightsNei_ph0[faceI] = 0;        
                diffFlux_ph0[faceI] = 0;

                faceI++;
            }
        }
    }

    fvScalarMatrix Y1iEqn
        (
            fvm::ddt(alpha1, Y1i) + fvc::div(diffFlux_ph1)
        );

    Y1iEqn.solve(mesh.solver("Yi"));

    Y1i.correctBoundaryConditions();

    Info<< "Phase-1 Yi= "
        << (alpha1*Y1i).weightedAverage(mesh.Vsc()).value()
        << "  Min(Y1i) = " << min(Y1i).value()
        << "  Max(Y1i) = " << max(Y1i).value()
        << nl
        << "Y1i undershoot = " << min(Y1i).value()
        << "  Y1i overshoot = " << (max(Y1i).value() - 1)
        << endl;

    Y1i.max(0.0);
    Y1i.min(1.0);

    Info<< "Phase-1 Yi after correction= "
        << (alpha1*Y1i).weightedAverage(mesh.Vsc()).value()
        << "  Min(Y1i) = " << min(Y1i).value()
        << "  Max(Y1i) = " << max(Y1i).value()
        << endl;


    fvScalarMatrix Y0iEqn
        (
            fvm::ddt(alpha1, Y0i) + fvc::div(diffFlux_ph1)
        );

    Y0iEqn.solve(mesh.solver("Yi"));

    Y0i.correctBoundaryConditions();

    Info<< "Phase-1 Yi= "
        << (alpha1*Y0i).weightedAverage(mesh.Vsc()).value()
        << "  Min(Y0i) = " << min(Y0i).value()
        << "  Max(Y0i) = " << max(Y0i).value()
        << nl
        << "Y0i undershoot = " << min(Y0i).value()
        << "  Y0i overshoot = " << (max(Y0i).value() - 1)
        << endl;

    Y0i.max(0.0);
    Y0i.min(1.0);

    Info<< "Phase-1 Yi after correction= "
        << (alpha1*Y0i).weightedAverage(mesh.Vsc()).value()
        << "  Min(Y0i) = " << min(Y0i).value()
        << "  Max(Y0i) = " << max(Y0i).value()
        << endl;

    Y1t += Y1i;
    Y0t += Y0i;
}

Y1[nSpecies - 1] = scalar(1.0) - Y1t;
Y1[nSpecies - 1].max(0.0);
Y1[nSpecies - 1].min(1.0);

Y0[nSpecies - 1] = scalar(1.0) - Y0t;
Y0[nSpecies - 1].max(0.0);
Y0[nSpecies - 1].min(1.0);
