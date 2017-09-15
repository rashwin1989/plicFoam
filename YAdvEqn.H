
//for(label i=0; i < C.size(); i++)
//{


volScalarField& Y1i = Y1[i];
volScalarField& Y0i = Y0[i];
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
