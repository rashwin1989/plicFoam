
if(alpha_debug)
{
    plicFuncs::print_line(osAlpha, 80);
    osAlpha<< "Solution Time = " << runTime.timeName() << "  alpha1 advection" << endl;
    plicFuncs::print_line(osAlpha, 80);
    osAlpha<< endl;    
}

scalarField& alphaCells = alpha1.internalField();
scalarField& alphaOldCells = alpha1_old.internalField();

volScalarField surfInt_phiAlpha1(fvc::surfaceIntegrate(phiAlpha1));
scalarField& surfInt_phiAlpha1Cells = surfInt_phiAlpha1.internalField();

forAll(alphaCells, cellI)
{
    alphaCells[cellI] = alphaOldCells[cellI] - surfInt_phiAlpha1Cells[cellI]*dt;
}

Info<< "Phase-1 volume fraction = "
    << alpha1.weightedAverage(mesh.Vsc()).value()
    << "  Min(alpha1) = " << min(alpha1).value()
    << "  Max(alpha1) = " << max(alpha1).value()
    << nl
    << "alpha1 undershoot = " << min(alpha1).value()
    << "  alpha1 overshoot = " << (max(alpha1).value() - 1)
    << endl;

nCells = mesh.nCells();

maxAlpha1 = 1;
minAlpha1 = 0;
maxAlpha1Cell = -1;
minAlpha1Cell = -1;

ALPHA_MAX_BOUND = 1 + ALPHA_BOUND_ERR;
ALPHA_MIN_BOUND = -ALPHA_BOUND_ERR;

if(alpha_debug)
{
    osAlpha<< endl;
    forAll(alphaCells, cellI)
    {
        if(alphaCells[cellI] < ALPHA_MIN_BOUND || alphaCells[cellI] > ALPHA_MAX_BOUND)
        {
            osAlpha<< "Unbounded alpha in cell " << cellI << " alpha1 = " << alphaCells[cellI] << endl;
            if(alphaCells[cellI] < ALPHA_MIN_BOUND)
            {
                osAlpha<< "alpha1 undershoot = " << alphaCells[cellI] << endl;
            }
            else
            {
                osAlpha<< "alpha1 overshoot = " << (alphaCells[cellI] - 1) << endl;
            }

            scalar div_alpha1 = 0;
            scalar div_phi = 0;

            forAll(owner, faceI)
            {
                if(owner[faceI] == cellI)
                {
                    osAlpha<< "Face " << faceI << " owner  phiAlpha1 = " << phiAlpha1[faceI] << "  phi = " << phi[faceI] << endl;
                    div_alpha1 += phiAlpha1[faceI];
                    div_phi += phi[faceI];
                }
                else if(neighbour[faceI] == cellI)
                {
                    osAlpha<< "Face " << faceI << " neighbour  phiAlpha1 = " << -phiAlpha1[faceI] << "  phi = " << -phi[faceI] << endl;
                    div_alpha1 -= phiAlpha1[faceI];
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
                        osAlpha<< "Patch " << patchI << " face " << faceI << "  phiAlpha1 = " << pphiAlpha1[faceI] << "  phi = " << pphi[faceI] << endl;
                        div_alpha1 += pphiAlpha1[faceI];
                        div_phi += pphi[faceI];
                    }
                }
            }

            osAlpha<< "Cell alpha1 old time = " << alpha1.oldTime().internalField()[cellI] << nl
                << "alpha1 face flux integrate = " << div_alpha1 << nl
                << "face flux integrate = " << div_phi << endl;

            div_alpha1 /= meshV[cellI];
            div_phi /= meshV[cellI];

            osAlpha<< "divergence term in eqn = " << div_alpha1 << nl
                << "div of velocity field = " << div_phi << nl
                << "alpha1 solution update = " << div_alpha1*runTime.deltaT().value() << nl            
                << "Cell alpha1 new time = " << alphaCells[cellI] << nl
                << endl;
        }
    }
}

if(redist_alpha)
{
    forAll(alphaCells, cellI)
    {
        if(alphaCells[cellI] < ALPHA_MIN_BOUND)
        {
            if(alpha_debug)
            {
                osAlpha<< "Correcting alpha1 in cell " << cellI << nl
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
                    if(alpha_debug)
                    {
                        osAlpha<< "alpha1 correction iteration no: " << nIters+1 << endl;
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

                    if(alpha_debug)
                    {
                        osAlpha<< "Nei cell with minimum non-zero alpha1: " << minAlphaCell << nl
                            << "Nei cell alpha1: " << alphaCells[minAlphaCell] << endl;
                    }

                    scalar tAlpha = alphaCells[minAlphaCell] + alphaCells[cellI];
                    alphaCells[minAlphaCell] = max(tAlpha, 0);
                    alphaCells[cellI] = min(tAlpha, 0);

                    if(alpha_debug)
                    {
                        osAlpha<< "New nei cell alpha1: " << alphaCells[minAlphaCell] << nl
                            << "New cell alpha1: " << alphaCells[cellI] << nl << endl;
                    }
            
                    nIters++;
                }while(alphaCells[cellI] < ALPHA_MIN_BOUND && !allNeiDone && nIters < ALPHA_BOUND_ITERS_MAX);        
            }        
        }

        if(alphaCells[cellI] > ALPHA_MAX_BOUND)
        {
            if(alpha_debug)
            {
                osAlpha<< "Correcting alpha1 in cell " << cellI << nl
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
                    if(alpha_debug)
                    {
                        osAlpha<< "alpha1 correction iteration no: " << nIters+1 << endl;
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

                    if(alpha_debug)
                    {
                        osAlpha<< "Nei cell with maximum alpha1 below 1: " << maxAlphaCell << nl
                            << "Nei cell alpha1: " << alphaCells[maxAlphaCell] << endl;
                    }
            
                    scalar tAlpha = alphaCells[maxAlphaCell] + alphaCells[cellI] - 1;
                    alphaCells[maxAlphaCell] = min(tAlpha, 1);
                    alphaCells[cellI] = max(tAlpha, 1);

                    if(alpha_debug)
                    {
                        osAlpha<< "New nei cell alpha1: " << alphaCells[maxAlphaCell] << nl
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

    alpha1.correctBoundaryConditions();

    Info<< "Phase-1 volume fraction after redistribution = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << nl
        << "alpha1 undershoot = " << min(alpha1).value()
        << "  alpha1 overshoot = " << (max(alpha1).value() - 1)
        << endl;
}

if(alpha_debug)
{
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

    osAlpha<< "alpha1 extrema locations: " << nl
        << "Min alpha1 cell: " << minAlpha1Cell << " Min cell loc: " << minAlpha1CellC << " Min alpha1 value: " << minAlpha1 << nl 
        << "Max alpha1 cell: " << maxAlpha1Cell  << " Max cell loc: " << maxAlpha1CellC << " Max alpha1 value: " << maxAlpha1 << nl
        << endl;
}

alpha1.max(0.0);
alpha1.min(1.0);

Info<< "Phase-1 volume fraction after clipping = "
    << alpha1.weightedAverage(mesh.Vsc()).value()
    << "  Min(alpha1) = " << min(alpha1).value()
    << "  Max(alpha1) = " << max(alpha1).value()
    << endl;

alpha1.correctBoundaryConditions();

limitedAlpha1 = min(max(alpha1, scalar(0)), scalar(1));

alpha0 = scalar(1) - limitedAlpha1;
