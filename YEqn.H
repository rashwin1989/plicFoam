
//for(label i=0; i < C.size(); i++)
//{

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

volScalarField rho1Alpha1(rho1*alpha1);
surfaceScalarField& phiRho1Y1i = phiRho1Y11;

fvScalarMatrix Y1iEqn
(
    fvm::ddt(rho1Alpha1, Y1i) + fvc::div(phiRho1Y1i)
);

Y1iEqn.solve();

volScalarField C1i(rho1Alpha1*Y1i);

Info<< "Total C = "
    << C1i.weightedAverage(mesh.Vsc()).value()
    << "  Min(Y) = " << min(Y1i).value()
    << "  Max(Y) = " << max(Y1i).value()
    << nl
    << "Y undershoot = " << min(Y1i).value()
    << "  Y overshoot = " << (max(Y1i).value() - 1)
    << endl;

scalarField& YCells = Y1i.internalField();
label nCells = mesh.nCells();

scalar maxY = 1;
scalar minY = 0;
label maxY1iCell = -1;
label minY1iCell = -1;

scalar Y_MAX_BOUND = 1 + Y_BOUND_ERR;
scalar Y_MIN_BOUND = -Y_BOUND_ERR;

const labelUList& owner = mesh.owner();
const labelUList& neighbour = mesh.neighbour();
const scalarField& meshV = mesh.V();

if(unbound_debug)
{
    Info<< endl;
    forAll(YCells, cellI)
    {
        if(YCells[cellI] < Y_MIN_BOUND || YCells[cellI] > Y_MAX_BOUND)
        {
            Info<< "Unbounded Y in cell " << cellI << " Y = " << YCells[cellI] << endl;
            if(YCells[cellI] < Y_MIN_BOUND)
            {
                Info<< "Y undershoot = " << YCells[cellI] << endl;
            }
            else
            {
                Info<< "Y overshoot = " << (YCells[cellI] - 1) << endl;
            }

            scalar div_phiRhoY = 0;
            scalar div_phi = 0;

            forAll(owner, faceI)
            {
                if(owner[faceI] == cellI)
                {
                    Info<< "Face " << faceI << " owner  phiRhoY = " << phiRho1Y1i[faceI] << "  phi = " << phi[faceI] << endl;
                    div_phiRhoY += phiRho1Y1i[faceI];
                    div_phi += phi[faceI];
                }
                else if(neighbour[faceI] == cellI)
                {
                    Info<< "Face " << faceI << " neighbour  phiRhoY = " << -phiRho1Y1i[faceI] << "  phi = " << -phi[faceI] << endl;
                    div_phiRhoY -= phiRho1Y11[faceI];
                    div_phi -= phi[faceI];
                }
            }

            forAll(mesh.boundary(), patchI)
            {
                const labelUList& pFaceCells = mesh.boundary()[patchI].faceCells();
                const fvsPatchScalarField& pphiAlpha1 = phiAlpha1.boundaryField()[patchI];
                const fvsPatchScalarField& pphi = phi.boundaryField()[patchI];

                forAll(mesh.boundary()[patchI], faceI)
                {
                    if(pFaceCells[faceI] == cellI)
                    {
                        Info<< "Patch " << patchI << " face " << faceI << "  phiAlpha1 = " << pphiAlpha1[faceI] << "  phi = " << pphi[faceI] << endl;
                        div_alpha1 += pphiAlpha1[faceI];
                        div_phi += pphi[faceI];
                    }
                }
            }

            Info<< "Cell alpha1 old time = " << alpha1.oldTime().internalField()[cellI] << nl
                << "alpha1 face flux integrate = " << div_alpha1 << nl
                << "face flux integrate = " << div_phi << endl;

            div_alpha1 /= meshV[cellI];
            div_phi /= meshV[cellI];

            Info<< "divergence term in eqn = " << div_alpha1 << nl
                << "div of velocity field = " << div_phi << nl
                << "alpha1 solution update = " << div_alpha1*runTime.deltaT().value() << nl            
                << "Cell alpha1 new time = " << alphaCells[cellI] << nl
                << endl;
        }
    }
}

forAll(alphaCells, cellI)
{
    if(alphaCells[cellI] < ALPHA_MIN_BOUND)
    {
        if(debug)
        {
            Info<< "Correcting alpha1 in cell " << cellI << nl
                << "Cell alpha1:  " << alphaCells[cellI] << nl << endl;
        }

        const labelList& curCellCells = cell_stencil.stencil()[cellI];
        scalar minAlpha = 1;
        label minAlphaCell = cellI;
        label nIters = 0;
        bool allNeiDone = true;
        for(label i=0; i<curCellCells.size(); i++)
        {
            label curCell = curCellCells[i];

            if(!(curCell==cellI) && (curCell < nCells) && (alphaCells[curCell] > 0))
            {
                allNeiDone = false;
                break;
            }
        }

        if(!allNeiDone)
        {
            do
            {
                if(debug)
                {
                    Info<< "alpha1 correction iteration no: " << nIters+1 << endl;
                }
            
                allNeiDone = true;
                minAlpha = 1;
                for(label i=0; i<curCellCells.size(); i++)
                {
                    label curCell = curCellCells[i];

                    if(!(curCell==cellI) && (curCell < nCells) && (alphaCells[curCell] > 0))
                    {
                        allNeiDone = false;
                        if(alphaCells[curCell] < minAlpha)
                        {
                            minAlphaCell = curCell;
                            minAlpha = alphaCells[curCell];
                        }
                    }
                }            

                if(debug)
                {
                    Info<< "Nei cell with minimum non-zero alpha1: " << minAlphaCell << nl
                        << "Nei cell alpha1: " << alphaCells[minAlphaCell] << endl;
                }

                scalar tAlpha = alphaCells[minAlphaCell] + alphaCells[cellI];
                alphaCells[minAlphaCell] = max(tAlpha, 0);
                alphaCells[cellI] = min(tAlpha, 0);

                if(debug)
                {
                    Info<< "New nei cell alpha1: " << alphaCells[minAlphaCell] << nl
                        << "New cell alpha1: " << alphaCells[cellI] << nl << endl;
                }
            
                nIters++;
            }while(alphaCells[cellI] < ALPHA_MIN_BOUND && !allNeiDone && nIters < ALPHA_BOUND_ITERS_MAX);        
        }        
    }

    if(alphaCells[cellI] > ALPHA_MAX_BOUND)
    {
        if(debug)
        {
            Info<< "Correcting alpha1 in cell " << cellI << nl
                << "Cell alpha1:  " << alphaCells[cellI] << nl << endl;
        }

        const labelList& curCellCells = cell_stencil.stencil()[cellI];
        scalar maxAlpha = 0;
        label maxAlphaCell = cellI;
        label nIters = 0;
        bool allNeiDone = true;
        for(label i=0; i<curCellCells.size(); i++)
        {
            label curCell = curCellCells[i];

            if(!(curCell==cellI) && (curCell < nCells) && (alphaCells[curCell] < 1))
            {
                allNeiDone = false;
                break;
            }
        }
        
        if(!allNeiDone)
        {
            do
            {
                if(debug)
                {
                    Info<< "alpha1 correction iteration no: " << nIters+1 << endl;
                }

                maxAlpha = 0;
                allNeiDone = true;
                for(label i=0; i<curCellCells.size(); i++)
                {
                    label curCell = curCellCells[i];

                    if(!(curCell==cellI) && (curCell < nCells) && (alphaCells[curCell] < 1))
                    {
                        allNeiDone = false;
                        if(alphaCells[curCell] > maxAlpha)
                        {
                            maxAlphaCell = curCell;
                            maxAlpha = alphaCells[curCell];
                        }
                    }
                }

                if(debug)
                {
                    Info<< "Nei cell with maximum alpha1 below 1: " << maxAlphaCell << nl
                        << "Nei cell alpha1: " << alphaCells[maxAlphaCell] << endl;
                }
            
                scalar tAlpha = alphaCells[maxAlphaCell] + alphaCells[cellI] - 1;
                alphaCells[maxAlphaCell] = min(tAlpha, 1);
                alphaCells[cellI] = max(tAlpha, 1);

                if(debug)
                {
                    Info<< "New nei cell alpha1: " << alphaCells[maxAlphaCell] << nl
                        << "New cell alpha1: " << alphaCells[cellI] << nl << endl;
                }

                nIters++;
            }while(alphaCells[cellI] > ALPHA_MAX_BOUND && !allNeiDone && nIters < ALPHA_BOUND_ITERS_MAX);        
        }
    }

    if(alphaCells[cellI] < minAlpha1)
    {
        minAlpha1 = alphaCells[cellI];
        minAlpha1Cell = cellI;
    }

    if(alphaCells[cellI] > maxAlpha1)
    {
        maxAlpha1 = alphaCells[cellI];
        maxAlpha1Cell = cellI;
    }
}

const vectorField& meshCi = mesh.C().internalField();
vector minAlpha1CellC = vector::zero;
vector maxAlpha1CellC = vector::zero;
if(minAlpha1Cell >= 0)
{
    minAlpha1CellC = meshCi[minAlpha1Cell];
}
if(maxAlpha1Cell >= 0)
{
    maxAlpha1CellC = meshCi[maxAlpha1Cell];
}

alpha1.correctBoundaryConditions();

Info<< "Phase-1 volume fraction after correction= "
    << alpha1.weightedAverage(mesh.Vsc()).value()
    << "  Min(alpha1) = " << min(alpha1).value()
    << "  Max(alpha1) = " << max(alpha1).value()
    << nl
    << "alpha1 undershoot = " << min(alpha1).value()
    << "  alpha1 overshoot = " << (max(alpha1).value() - 1)
    << endl;

/*Info<< "alpha1 extrema locations: " << nl
<< "Min alpha1 cell: " << minAlpha1Cell << " Min cell loc: " << minAlpha1CellC << " Min alpha1 value: " << minAlpha1 << nl 
    << "Max alpha1 cell: " << maxAlpha1Cell  << " Max cell loc: " << maxAlpha1CellC << " Max alpha1 value: " << maxAlpha1 << endl;*/

alpha1.max(0.0);
alpha1.min(1.0);

Info<< "Phase-1 volume fraction after correction= "
    << alpha1.weightedAverage(mesh.Vsc()).value()
    << "  Min(alpha1) = " << min(alpha1).value()
    << "  Max(alpha1) = " << max(alpha1).value()
    << endl;

//}
