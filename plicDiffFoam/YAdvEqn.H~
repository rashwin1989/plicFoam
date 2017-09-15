
nSpecies = Y1.size();

volScalarField Y1t = 0.0*Y1[0];
volScalarField Y0t = 0.0*Y0[0];

scalarField& alpha1Cells = alpha1.internalField();
scalarField& alpha0Cells = alpha0.internalField();

plicFuncs::makeFlatFld(alpha1, diffMap, alpha1_flatFld_diff);
plicFuncs::makeFlatFld(interface.C_ph1(), diffMap, C_ph1_flatFld_diff);
plicFuncs::makeFlatFld(interface.C_ph0(), diffMap, C_ph0_flatFld_diff);

if(diff_debug2)
{
Info<< "alpha1_flatFld" << nl
    << alpha1_flatFld_diff << nl
    << endl;
}

const labelList& face_phaseState = interface.face_phaseState();

for(label i=0; i<(nSpecies - 1); i++)
{
    if(diff_debug2)
    {
        Info<< "-----------------------------------" << nl
            << "              Field Y[" << i << "]" << nl
            << "-----------------------------------" << nl
            << endl;
    }

    volScalarField& Y1i = Y1[i];
    volScalarField& Y0i = Y0[i];
    volScalarField& C1i = C1[i];
    volScalarField& C0i = C0[i];

    scalarField& Y1iCells = Y1i.internalField();
    scalarField& Y0iCells = Y0i.internalField();
    scalarField& C1iCells = C1i.internalField();
    scalarField& C0iCells = C0i.internalField();

    List<scalar>& Y1i_flatFld_diff = Y1_flatFld_diff[i];
    List<scalar>& Y0i_flatFld_diff = Y0_flatFld_diff[i];
    plicFuncs::makeFlatFld(Y1i, diffMap, Y1i_flatFld_diff);
    plicFuncs::makeFlatFld(Y0i, diffMap, Y0i_flatFld_diff);

    scalar& D1i = D1[i].value();
    scalar& D0i = D0[i].value();   
    const scalarField& magSf_ph1 = interface.Af_ph1();
    const scalarField& magSf_ph0 = interface.Af_ph0();

    scalar MIN_ALPHA_DIFF = 1E-06;
    scalar MAX_ALPHA_DIFF = 1 - MIN_ALPHA_DIFF;    

    if(diff_debug2)
    {
        Info<< "Diffusion flux calculation" << nl
            << nl
            << "Internal faces" << nl
            << endl;
    }

    //Internal faces    
    for(label faceI=0; faceI<mesh.nInternalFaces(); faceI++)
    {
        scalar curMagSf = meshMagSf[faceI];
        scalar curMagSf_ph1 = magSf_ph1[faceI];
        scalar curMagSf_ph0 = magSf_ph0[faceI];
        scalar curAlpha1f = curMagSf_ph1/curMagSf;
        scalar curAlpha0f = curMagSf_ph0/curMagSf;
        vector nf = meshSf[faceI]/curMagSf;
        const label& faceOwn = own[faceI];
        const label& faceNei = nei[faceI];
        label donor = own[faceI];
        label rec = nei[faceI];

        if(diff_debug2)
        {
            Info<< "Face: " << faceI << nl
                << "Own: " << faceOwn << "  Nei: " << faceNei << nl
                << "Own Y_ph1: " << Y1i_flatFld_diff[faceOwn] << "  Nei Y_ph1: " << Y1i_flatFld_diff[faceNei] << nl
                << "Own Y_ph0: " << Y0i_flatFld_diff[faceOwn] << "  Nei Y_ph0: " << Y0i_flatFld_diff[faceNei] << nl
                << "mag(Sf) = " << curMagSf << "  mag(Sf)_ph1 =  " << curMagSf_ph1 << "  alpha1f = " << curAlpha1f << nl
                << endl;
        }

        scalar wOwn = 1;
        scalar wNei = 1;

        if(curAlpha1f < MIN_ALPHA_DIFF && curAlpha0f < MIN_ALPHA_DIFF)
        {
            diffFlux_ph0[faceI] = 0;
            diffFlux_ph1[faceI] = 0;
        }
        else if(curAlpha1f < MIN_ALPHA_DIFF)
        {
            plicFuncs::calcTwoSidedFaceGradWeights
            (
                faceI,
                faceOwn,
                faceNei,
                nf,
                Y0i_flatFld_diff,
                alpha1_flatFld_diff,
                C_ph0_flatFld_diff,
                diffCellStencil,
                0.1*MIN_ALPHA_DIFF,
                0,
                wOwn,
                wNei,
                diff_debug2
            );
            
            diffFlux_ph0[faceI] = -D0i*magSf_ph0[faceI]*(wOwn*Y0i_flatFld_diff[faceOwn] + wNei*Y0i_flatFld_diff[faceNei]);
            
            diffFlux_ph1[faceI] = 0;
        }//end if(curAlpha1f < MIN_ALPHA_DIFF)
        else if(curAlpha1f > MAX_ALPHA_DIFF)
        {
            plicFuncs::calcTwoSidedFaceGradWeights
            (
                faceI,
                faceOwn,
                faceNei,
                nf,
                Y1i_flatFld_diff,
                alpha1_flatFld_diff,
                C_ph1_flatFld_diff,
                diffCellStencil,
                0.1*MIN_ALPHA_DIFF,
                1,
                wOwn,
                wNei,
                diff_debug2
            );
            
            diffFlux_ph1[faceI] = -D1i*magSf_ph1[faceI]*(wOwn*Y1i_flatFld_diff[faceOwn] + wNei*Y1i_flatFld_diff[faceNei]);
            
            diffFlux_ph0[faceI] = 0;
        }//end else if(curAlpha1f > MAX_ALPHA_DIFF)
        else
        {
            plicFuncs::calcTwoSidedFaceGradWeights
            (
                faceI,
                faceOwn,
                faceNei,
                nf,
                Y0i_flatFld_diff,
                alpha1_flatFld_diff,
                C_ph0_flatFld_diff,
                diffCellStencil,
                0.1*MIN_ALPHA_DIFF,
                0,
                wOwn,
                wNei,
                diff_debug2
            );
            
            diffFlux_ph0[faceI] = -D0i*magSf_ph0[faceI]*(wOwn*Y0i_flatFld_diff[faceOwn] + wNei*Y0i_flatFld_diff[faceNei]);

            plicFuncs::calcTwoSidedFaceGradWeights
            (
                faceI,
                faceOwn,
                faceNei,
                nf,
                Y1i_flatFld_diff,
                alpha1_flatFld_diff,
                C_ph1_flatFld_diff,
                diffCellStencil,
                0.1*MIN_ALPHA_DIFF,
                1,
                wOwn,
                wNei,
                diff_debug2
            );
            
            diffFlux_ph1[faceI] = -D1i*magSf_ph1[faceI]*(wOwn*Y1i_flatFld_diff[faceOwn] + wNei*Y1i_flatFld_diff[faceNei]);
        }//end if(curAlpha1f < MIN_ALPHA_DIFF)

        if(diff_debug2)
        {
            Info<< nl
                << "diffusion flux ph1 for Y" << i << " = " << diffFlux_ph1[faceI] << nl
                << "diffusion flux ph0 for Y" << i << " = " << diffFlux_ph0[faceI] << nl
                << endl;
        }

        //check diffFlux > maxDiffFlux causing unboundedness
        //ph1
        if(diffFlux_ph1[faceI] > 0)
        {
            donor = faceOwn;
            rec = faceNei;
        }
        else
        {
            donor = faceNei;
            rec = faceOwn;
        }
        scalar maxDiffFlux = Y1iCells[donor]*alpha1Cells[donor]*mesh.V()[donor]/runTime.deltaTValue();
        if(mag(diffFlux_ph1[faceI]) > maxDiffFlux)
        {
            if(diff_debug)
            {
                Info<< "Face " << faceI << "  diff flux ph 1: " << diffFlux_ph1[faceI] << nl
                    << "Af_ph1: " << magSf_ph1[faceI]/meshMagSf[faceI] << "  Af_ph0: " << magSf_ph0[faceI]/meshMagSf[faceI] << "  phaseState: " << face_phaseState[faceI] << nl
                    << "Own " << faceOwn << "  Nei " << faceNei << "  Donor " << donor << nl
                    << "MaxDiffFlux: " << maxDiffFlux << nl 
                    << "Don alpha1: " << alpha1Cells[donor] << "  Y1: " <<  Y1iCells[donor] << "  deltaT: " << runTime.deltaTValue() << "  V: " <<  mesh.V()[donor] << nl 
                    << "Rec alpha1: " << alpha1Cells[rec] << "  Y1: " <<  Y1iCells[rec] << "  deltaT: " << runTime.deltaTValue() << "  V: " <<  mesh.V()[rec] << nl 
                    << endl;
            }
        }

        //ph0
        if(diffFlux_ph0[faceI] > 0)
        {
            donor = faceOwn;
            rec = faceNei;
        }
        else
        {
            donor = faceNei;
            rec = faceOwn;
        }
        maxDiffFlux = Y0iCells[donor]*(1 - alpha1Cells[donor])*mesh.V()[donor]/runTime.deltaTValue();
        if(mag(diffFlux_ph0[faceI]) > maxDiffFlux)
        {
            if(diff_debug)
            {
                Info<< "Face " << faceI << "  diff flux ph 0: " << diffFlux_ph0[faceI] << nl
                    << "Af_ph1: " << magSf_ph1[faceI]/meshMagSf[faceI] << "  Af_ph0: " << magSf_ph0[faceI]/meshMagSf[faceI] << "  phaseState: " << face_phaseState[faceI] << nl
                    << "Own " << faceOwn << "  Nei " << faceNei << "  Donor " << donor << nl
                    << "MaxDiffFlux: " << maxDiffFlux << nl 
                    << "Don alpha0: " << (1 - alpha1Cells[donor]) << "  Y0: " <<  Y0iCells[donor] << "  deltaT: " << runTime.deltaTValue() << "  V: " <<  mesh.V()[donor] << nl 
                    << "Rec alpha0: " << (1 - alpha1Cells[rec]) << "  Y0: " <<  Y0iCells[rec] << "  deltaT: " << runTime.deltaTValue() << "  V: " <<  mesh.V()[rec] << nl 
                    << endl;
            }
        }
        //end check diffFlux > maxDiffFlux
    }//end for(label faceI=0; faceI<mesh.nInternalFaces(); faceI++)

    //end internal faces

    if(diff_debug2)
    {
        Info<< "Boundary faces" << nl
            << endl;
    }

    //Boundary faces
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const label nBnd = mesh.nFaces() - mesh.nInternalFaces();

    List<scalar> ownYi_ph1(nBnd);
    List<scalar> ownAlpha_ph1(nBnd);
    List<scalar> ownBeta_ph1(nBnd);
    List<scalar> ownMagt1_ph1(nBnd);
    List<scalar> ownMagt2_ph1(nBnd);
    List<scalar> ownD_ph1(nBnd);
    List<scalar> neiYi_ph1(nBnd);
    List<scalar> neiAlpha_ph1(nBnd);
    List<scalar> neiBeta_ph1(nBnd);
    List<scalar> neiMagt1_ph1(nBnd);
    List<scalar> neiMagt2_ph1(nBnd);
    List<scalar> neiD_ph1(nBnd);

    List<scalar> ownYi_ph0(nBnd);
    List<scalar> ownAlpha_ph0(nBnd);
    List<scalar> ownBeta_ph0(nBnd);
    List<scalar> ownMagt1_ph0(nBnd);
    List<scalar> ownMagt2_ph0(nBnd);
    List<scalar> ownD_ph0(nBnd);
    List<scalar> neiYi_ph0(nBnd);
    List<scalar> neiAlpha_ph0(nBnd);
    List<scalar> neiBeta_ph0(nBnd);
    List<scalar> neiMagt1_ph0(nBnd);
    List<scalar> neiMagt2_ph0(nBnd);
    List<scalar> neiD_ph0(nBnd);

    forAll(patches,patchI)
    {
        const polyPatch& pp = patches[patchI];

        if(pp.coupled())
        {
            label faceI = pp.start();

            forAll(pp, fcI)
            {
                label bndFaceI = faceI - mesh.nInternalFaces();
                scalar curMagSf = meshMagSf[faceI];
                scalar curMagSf_ph1 = magSf_ph1[faceI];
                scalar curAlpha1f = curMagSf_ph1/curMagSf;
                vector nf = meshSf[faceI]/curMagSf;
                const label& faceOwn = own[faceI];
                const labelList& ownCells = diffCellStencil[faceOwn];

                ownYi_ph1 = Y1iCells[faceOwn];
                neiYi_ph1 = Y1iCells[faceOwn];
                ownYi_ph0 = Y0iCells[faceOwn];
                neiYi_ph0 = Y0iCells[faceOwn];

                scalar alpha;
                scalar beta;
                scalar magt1;
                scalar magt2;
                scalar d;

                if(curAlpha1f < MIN_ALPHA_DIFF)
                {                    
                    plicFuncs::calcCellGradWeights
                    (
                        faceOwn,
                        nf,
                        Y0i_flatFld_diff,
                        alpha1_flatFld_diff,
                        C_ph0_flatFld_diff,
                        ownCells,
                        0.1*MIN_ALPHA_DIFF,
                        0,
                        alpha,
                        beta,
                        magt1,
                        magt2,
                        d,
                        diff_debug2
                    );
                    ownAlpha_ph0[bndFaceI] = alpha;
                    ownBeta_ph0[bndFaceI] = beta;
                    ownMagt1_ph0[bndFaceI] = magt1;
                    ownMagt2_ph0[bndFaceI] = magt2;
                    ownD_ph0[bndFaceI] = d;
                    neiAlpha_ph0[bndFaceI] = alpha;
                    neiBeta_ph0[bndFaceI] = beta;
                    neiMagt1_ph0[bndFaceI] = magt1;
                    neiMagt2_ph0[bndFaceI] = magt2;
                    neiD_ph0[bndFaceI] = d;                    
                }//end if(curAlpha1f < MIN_ALPHA_DIFF)
                else if(curAlpha1f > MAX_ALPHA_DIFF)
                {
                    plicFuncs::calcCellGradWeights
                    (
                        faceOwn,
                        nf,
                        Y1i_flatFld_diff,
                        alpha1_flatFld_diff,
                        C_ph1_flatFld_diff,
                        ownCells,
                        0.1*MIN_ALPHA_DIFF,
                        1,
                        alpha,
                        beta,
                        magt1,
                        magt2,
                        d,
                        diff_debug2
                    );
                    ownAlpha_ph1[bndFaceI] = alpha;
                    ownBeta_ph1[bndFaceI] = beta;
                    ownMagt1_ph1[bndFaceI] = magt1;
                    ownMagt2_ph1[bndFaceI] = magt2;
                    ownD_ph1[bndFaceI] = d;
                    neiAlpha_ph1[bndFaceI] = alpha;
                    neiBeta_ph1[bndFaceI] = beta;
                    neiMagt1_ph1[bndFaceI] = magt1;
                    neiMagt2_ph1[bndFaceI] = magt2;
                    neiD_ph1[bndFaceI] = d;
                }//end else if(curAlpha1f > MAX_ALPHA_DIFF)
                else
                {                    
                    plicFuncs::calcCellGradWeights
                    (
                        faceOwn,
                        nf,
                        Y0i_flatFld_diff,
                        alpha1_flatFld_diff,
                        C_ph0_flatFld_diff,
                        ownCells,
                        0.1*MIN_ALPHA_DIFF,
                        0,
                        alpha,
                        beta,
                        magt1,
                        magt2,
                        d,
                        diff_debug2
                    );
                    ownAlpha_ph0[bndFaceI] = alpha;
                    ownBeta_ph0[bndFaceI] = beta;
                    ownMagt1_ph0[bndFaceI] = magt1;
                    ownMagt2_ph0[bndFaceI] = magt2;
                    ownD_ph0[bndFaceI] = d;
                    neiAlpha_ph0[bndFaceI] = alpha;
                    neiBeta_ph0[bndFaceI] = beta;
                    neiMagt1_ph0[bndFaceI] = magt1;
                    neiMagt2_ph0[bndFaceI] = magt2;
                    neiD_ph0[bndFaceI] = d;

                    plicFuncs::calcCellGradWeights
                    (
                        faceOwn,
                        nf,
                        Y1i_flatFld_diff,
                        alpha1_flatFld_diff,
                        C_ph1_flatFld_diff,
                        ownCells,
                        0.1*MIN_ALPHA_DIFF,
                        1,
                        alpha,
                        beta,
                        magt1,
                        magt2,
                        d,
                        diff_debug2
                    );
                    ownAlpha_ph1[bndFaceI] = alpha;
                    ownBeta_ph1[bndFaceI] = beta;
                    ownMagt1_ph1[bndFaceI] = magt1;
                    ownMagt2_ph1[bndFaceI] = magt2;
                    ownD_ph1[bndFaceI] = d;
                    neiAlpha_ph1[bndFaceI] = alpha;
                    neiBeta_ph1[bndFaceI] = beta;
                    neiMagt1_ph1[bndFaceI] = magt1;
                    neiMagt2_ph1[bndFaceI] = magt2;
                    neiD_ph1[bndFaceI] = d;
                }//end if(curAlpha1f < MIN_ALPHA_DIFF)

                faceI++;
            }//end forAll(pp, fcI)
        }//end if(pp.coupled())
    }//end forAll(patches,patchI)

    syncTools::swapBoundaryFaceList(mesh, neiYi_ph0);
    syncTools::swapBoundaryFaceList(mesh, neiAlpha_ph0);
    syncTools::swapBoundaryFaceList(mesh, neiBeta_ph0);
    syncTools::swapBoundaryFaceList(mesh, neiMagt1_ph0);
    syncTools::swapBoundaryFaceList(mesh, neiMagt2_ph0);
    syncTools::swapBoundaryFaceList(mesh, neiD_ph0);

    syncTools::swapBoundaryFaceList(mesh, neiYi_ph1);
    syncTools::swapBoundaryFaceList(mesh, neiAlpha_ph1);
    syncTools::swapBoundaryFaceList(mesh, neiBeta_ph1);
    syncTools::swapBoundaryFaceList(mesh, neiMagt1_ph1);
    syncTools::swapBoundaryFaceList(mesh, neiMagt2_ph1);
    syncTools::swapBoundaryFaceList(mesh, neiD_ph1);

    //Info<< "Done syncing coupled info" << endl;

    forAll(patches,patchI)
    {
        const polyPatch& pp = patches[patchI];        
        const fvPatchScalarField& pY1i = Y1i.boundaryField()[patchI];
        fvsPatchScalarField pdiffFlux_ph0 = diffFlux_ph0.boundaryField()[patchI];
        fvsPatchScalarField pdiffFlux_ph1 = diffFlux_ph1.boundaryField()[patchI];
        label faceI = pp.start();
        scalar faceGrad = 0;

        if(pp.coupled())
        {
            forAll(pY1i, fcI)
            {
                label bndFaceI = faceI - mesh.nInternalFaces();
                scalar curMagSf = meshMagSf[faceI];
                scalar curMagSf_ph1 = magSf_ph1[faceI];
                scalar curAlpha1f = curMagSf_ph1/curMagSf;

                if(curAlpha1f < MIN_ALPHA_DIFF)
                {                
                    plicFuncs::calcFaceGradFromWeights
                    (
                        ownYi_ph0[bndFaceI],
                        neiYi_ph0[bndFaceI],
                        ownAlpha_ph0[bndFaceI],
                        ownBeta_ph0[bndFaceI],
                        ownMagt1_ph0[bndFaceI],
                        ownMagt2_ph0[bndFaceI],
                        ownD_ph0[bndFaceI],
                        neiAlpha_ph0[bndFaceI],
                        neiBeta_ph0[bndFaceI],
                        neiMagt1_ph0[bndFaceI],
                        neiMagt2_ph0[bndFaceI],
                        neiD_ph0[bndFaceI],
                        faceGrad,
                        diff_debug2
                    );
                    pdiffFlux_ph0[fcI] = -D0i*magSf_ph0[faceI]*faceGrad;

                    pdiffFlux_ph1[fcI] = 0;
                }//end if(curAlpha1f < MIN_ALPHA_DIFF)
                else if(curAlpha1f > MAX_ALPHA_DIFF)
                {
                    plicFuncs::calcFaceGradFromWeights
                    (
                        ownYi_ph1[bndFaceI],
                        neiYi_ph1[bndFaceI],
                        ownAlpha_ph1[bndFaceI],
                        ownBeta_ph1[bndFaceI],
                        ownMagt1_ph1[bndFaceI],
                        ownMagt2_ph1[bndFaceI],
                        ownD_ph1[bndFaceI],
                        neiAlpha_ph1[bndFaceI],
                        neiBeta_ph1[bndFaceI],
                        neiMagt1_ph1[bndFaceI],
                        neiMagt2_ph1[bndFaceI],
                        neiD_ph1[bndFaceI],
                        faceGrad,
                        diff_debug2
                    );
                    pdiffFlux_ph1[fcI] = -D1i*magSf_ph1[faceI]*faceGrad;

                    pdiffFlux_ph0[fcI] = 0;
                }//end else if(curAlpha1f > MAX_ALPHA_DIFF)
                else
                {
                    plicFuncs::calcFaceGradFromWeights
                    (
                        ownYi_ph0[bndFaceI],
                        neiYi_ph0[bndFaceI],
                        ownAlpha_ph0[bndFaceI],
                        ownBeta_ph0[bndFaceI],
                        ownMagt1_ph0[bndFaceI],
                        ownMagt2_ph0[bndFaceI],
                        ownD_ph0[bndFaceI],
                        neiAlpha_ph0[bndFaceI],
                        neiBeta_ph0[bndFaceI],
                        neiMagt1_ph0[bndFaceI],
                        neiMagt2_ph0[bndFaceI],
                        neiD_ph0[bndFaceI],
                        faceGrad,
                        diff_debug2
                    );
                    pdiffFlux_ph0[fcI] = -D0i*magSf_ph0[faceI]*faceGrad;

                    plicFuncs::calcFaceGradFromWeights
                    (
                        ownYi_ph1[bndFaceI],
                        neiYi_ph1[bndFaceI],
                        ownAlpha_ph1[bndFaceI],
                        ownBeta_ph1[bndFaceI],
                        ownMagt1_ph1[bndFaceI],
                        ownMagt2_ph1[bndFaceI],
                        ownD_ph1[bndFaceI],
                        neiAlpha_ph1[bndFaceI],
                        neiBeta_ph1[bndFaceI],
                        neiMagt1_ph1[bndFaceI],
                        neiMagt2_ph1[bndFaceI],
                        neiD_ph1[bndFaceI],
                        faceGrad,
                        diff_debug2
                    );
                    pdiffFlux_ph1[fcI] = -D1i*magSf_ph1[faceI]*faceGrad;
                }//end if(curAlpha1f < MIN_ALPHA_DIFF)
                
                faceI++;
            }//end forAll(pY1i, fcI)
        }//end if(pp.coupled())
        else if(isA<zeroGradientFvPatchScalarField>(pY1i))
        {            
            forAll(pY1i, fcI)
            {                
                pdiffFlux_ph1[fcI] = 0;                
                pdiffFlux_ph0[fcI] = 0;
                
                faceI++;
            }
        }//end if(isA<zeroGradientFvPatchScalarField>(pY1i))
        else if(isA<fixedValueFvPatchScalarField>(pY1i))
        {
            forAll(pY1i, fcI)
            {                
                scalar curMagSf = meshMagSf[faceI];
                scalar curMagSf_ph1 = magSf_ph1[faceI];
                scalar curAlpha1f = curMagSf_ph1/curMagSf;
                vector nf = meshSf[faceI]/curMagSf;
                const label& faceOwn = own[faceI];
                const labelList& ownCells = diffCellStencil[faceOwn];
                scalar ownCellGrad = 0;
                
                if(curAlpha1f < MIN_ALPHA_DIFF)
                {                
                    plicFuncs::calcCellGrad
                    (
                        faceOwn,
                        nf,
                        Y0i_flatFld_diff,
                        alpha1_flatFld_diff,
                        C_ph0_flatFld_diff,
                        ownCells,
                        0.1*MIN_ALPHA_DIFF,
                        0,
                        ownCellGrad,
                        diff_debug2
                    );
                    pdiffFlux_ph0[fcI] = -D0i*magSf_ph0[faceI]*ownCellGrad;

                    pdiffFlux_ph1[fcI] = 0;
                }//end if(curAlpha1f < MIN_ALPHA_DIFF)
                else if(curAlpha1f > MAX_ALPHA_DIFF)
                {
                    plicFuncs::calcCellGrad
                    (
                        faceOwn,
                        nf,
                        Y1i_flatFld_diff,
                        alpha1_flatFld_diff,
                        C_ph1_flatFld_diff,
                        ownCells,
                        0.1*MIN_ALPHA_DIFF,
                        1,
                        ownCellGrad,
                        diff_debug2
                    );
                    pdiffFlux_ph1[fcI] = -D1i*magSf_ph1[faceI]*ownCellGrad;

                    pdiffFlux_ph0[fcI] = 0;
                }//end else if(curAlpha1f > MAX_ALPHA_DIFF)
                else
                {
                    plicFuncs::calcCellGrad
                    (
                        faceOwn,
                        nf,
                        Y0i_flatFld_diff,
                        alpha1_flatFld_diff,
                        C_ph0_flatFld_diff,
                        ownCells,
                        0.1*MIN_ALPHA_DIFF,
                        0,
                        ownCellGrad,
                        diff_debug2
                    );
                    pdiffFlux_ph0[fcI] = -D0i*magSf_ph0[faceI]*ownCellGrad;

                    plicFuncs::calcCellGrad
                    (
                        faceOwn,
                        nf,
                        Y1i_flatFld_diff,
                        alpha1_flatFld_diff,
                        C_ph1_flatFld_diff,
                        ownCells,
                        0.1*MIN_ALPHA_DIFF,
                        1,
                        ownCellGrad,
                        diff_debug2
                    );
                    pdiffFlux_ph1[fcI] = -D1i*magSf_ph1[faceI]*ownCellGrad;
                }//end if(curAlpha1f < MIN_ALPHA_DIFF)

                faceI++;
            }//end forAll(pY1i, fcI)
        }//end else if(isA<fixedValueFvPatchScalarField>(pY1i))
        else
        {
            forAll(pY1i, fcI)
            {                
                pdiffFlux_ph1[fcI] = 0;                
                pdiffFlux_ph0[fcI] = 0;

                faceI++;
            }
        }//end if(pp.coupled())        
    }//end forAll(patches,patchI)

    //end boundary faces

    if(diff_debug)
    {
        Info<<"-------------------------------------" << nl
            << "Diffusive fluxes at faces" << nl
            << "------------------------------------" << nl
            << endl;

        Info<< "-----------------------------------------------------------------------------------------------------" << nl
            << " Face    face_phaseState        alpha1f                alpha0f        diffFlux_ph1       diffFlux_ph0" << nl
            << "-----------------------------------------------------------------------------------------------------" << nl
            << endl;
        
        for(label faceI=0; faceI<mesh.nInternalFaces(); faceI++)
        {
            Info<< "   " << faceI << "              " << face_phaseState[faceI] << "                " << magSf_ph1[faceI]/meshMagSf[faceI] << "                " << magSf_ph0[faceI]/meshMagSf[faceI] << "                " << diffFlux_ph1[faceI] << "                " << diffFlux_ph0[faceI] << endl;
        }

        Info<< "------------------------------------------------------------------------------------------------------" << nl
            << endl;
    }


    /*fvScalarMatrix Y1iEqn
    (
        fvm::ddt(alpha1, Y1i) - fvc::div(diffFlux_ph1)
    );

        Y1iEqn.solve(mesh.solver("Yi"));*/

    volScalarField surfInt_diffFlux_ph1(fvc::surfaceIntegrate(diffFlux_ph1));    
    scalarField& surfInt_diffFlux_ph1Cells = surfInt_diffFlux_ph1.internalField();

    forAll(C1iCells,cellI)
    {
        if(diff_debug2)
        {
            Info<< "Solving Y1[" << i << "] eqn in cell " << cellI << nl
                << "Y1i = " << Y1iCells[cellI] << "  alpha1 = " << alpha1Cells[cellI] << nl
                << "div(diffFlux) = " << surfInt_diffFlux_ph1Cells[cellI] << "  div(diffFlux)*dt = " << surfInt_diffFlux_ph1Cells[cellI]*runTime.deltaTValue()
                << endl;
        }
        
        //C1iCells[cellI] = alpha1Cells[cellI]*Y1iCells[cellI];
        C1iCells[cellI] += -surfInt_diffFlux_ph1Cells[cellI]*runTime.deltaTValue();

        if(alpha1Cells[cellI] < 1E-15)
        {
            Y1iCells[cellI] = 0;
        }
        else
        {
            Y1iCells[cellI] = C1iCells[cellI]/alpha1Cells[cellI];
        }

        if(diff_debug2)
        {
            Info<< "C1i = " << C1iCells[cellI] << nl
                << "Y1i = " << Y1iCells[cellI] << nl
                << endl;
        }
    }

    C1i.correctBoundaryConditions();
    Y1i.correctBoundaryConditions();

    Info<< "Phase-1 Yi= "
        << (C1i).weightedAverage(mesh.Vsc()).value()
        << "  Min(Y1i) = " << min(Y1i).value()
        << "  Max(Y1i) = " << max(Y1i).value()
        << nl
        << "Y1i undershoot = " << min(min(Y1i).value(), 0)
        << "  Y1i overshoot = " << max((max(Y1i).value() - 1), 0)
        << endl;

    //Y1i.max(0.0);
    //Y1i.min(1.0);

    //C1i = alpha1*Y1i;

    Info<< "Phase-1 Yi after correction= "
        << (C1i).weightedAverage(mesh.Vsc()).value()
        << "  Min(Y1i) = " << min(Y1i).value()
        << "  Max(Y1i) = " << max(Y1i).value()
        << endl;


    /*fvScalarMatrix Y0iEqn
    (
        fvm::ddt(alpha0, Y0i) - fvc::div(diffFlux_ph0)
    );

        Y0iEqn.solve(mesh.solver("Yi"));*/

volScalarField surfInt_diffFlux_ph0(fvc::surfaceIntegrate(diffFlux_ph0));
    scalarField& surfInt_diffFlux_ph0Cells = surfInt_diffFlux_ph0.internalField();

    forAll(C0iCells,cellI)
    {
        if(diff_debug2)
        {
            Info<< "Solving Y0[" << i << "] eqn in cell " << cellI << nl
                << "Y0i = " << Y0iCells[cellI] << "  alpha0 = " << alpha0Cells[cellI] << nl
                << "div(diffFlux) = " << surfInt_diffFlux_ph0Cells[cellI] << "  div(diffFlux)*dt = " << surfInt_diffFlux_ph0Cells[cellI]*runTime.deltaTValue()
                << endl;
        }
        
        //C0iCells[cellI] = alpha0Cells[cellI]*Y0iCells[cellI];
        C0iCells[cellI] += -surfInt_diffFlux_ph0Cells[cellI]*runTime.deltaTValue();

        if(alpha0Cells[cellI] < 1E-15)
        {
            Y0iCells[cellI] = 0;
        }
        else
        {
            Y0iCells[cellI] = C0iCells[cellI]/alpha0Cells[cellI];
        }

        if(diff_debug2)
        {
            Info<< "C0i = " << C0iCells[cellI] << nl
                << "Y0i = " << Y0iCells[cellI] << nl
                << endl;
        }
    }

    C0i.correctBoundaryConditions();
    Y0i.correctBoundaryConditions();

    Info<< "Phase-0 Yi= "
        << (C0i).weightedAverage(mesh.Vsc()).value()
        << "  Min(Y0i) = " << min(Y0i).value()
        << "  Max(Y0i) = " << max(Y0i).value()
        << nl
        << "Y0i undershoot = " << min(min(Y0i).value(), 0)
        << "  Y0i overshoot = " << max((max(Y0i).value() - 1), 0)
        << endl;

    //Y0i.max(0.0);
    //Y0i.min(1.0);

    //C0i = alpha0*Y0i;

    Info<< "Phase-0 Yi after correction= "
        << (C0i).weightedAverage(mesh.Vsc()).value()
        << "  Min(Y0i) = " << min(Y0i).value()
        << "  Max(Y0i) = " << max(Y0i).value()
        << endl;

    Y1t += Y1i;
    Y0t += Y0i;
}

Y1[nSpecies - 1] = scalar(1.0) - Y1t;
Y1[nSpecies - 1].max(0.0);
Y1[nSpecies - 1].min(1.0);
C1[nSpecies - 1] = alpha1*Y1[nSpecies - 1];

Y0[nSpecies - 1] = scalar(1.0) - Y0t;
Y0[nSpecies - 1].max(0.0);
Y0[nSpecies - 1].min(1.0);
C0[nSpecies - 1] = alpha0*Y0[nSpecies - 1];

Info<< "Done YEqns" << endl;
