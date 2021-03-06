
plicFuncs::print_line(osAlpha, 100);
osAlpha<< "Solution Time = " << runTime.timeName() << endl;
plicFuncs::print_line(osAlpha, 100);
osAlpha<< endl;

volScalarField surfInt_phiAlpha1(fvc::surfaceIntegrate(phiAlpha1));
scalarField& surfInt_phiAlpha1Cells = surfInt_phiAlpha1.internalField();

forAll(alpha1Cells, cellI)
{
    alpha1Cells[cellI] = alpha1OldCells[cellI] - surfInt_phiAlpha1Cells[cellI]*dt;    
}

Info<< "alpha1:  avg = "
    << alpha1.weightedAverage(mesh.V()).value()
    << "  min = " << min(alpha1).value()
    << "  undershoot = " << min(min(alpha1).value(), 0)
    << "  max = " << max(alpha1).value()
    << "  overshoot = " << max((max(alpha1).value() - 1), 0)
    << endl;

maxAlpha1 = 1;
minAlpha1 = 0;
maxAlpha1Cell = -1;
minAlpha1Cell = -1;

ALPHA_MIN_BOUND = -ALPHA_BOUND_ERR;
ALPHA_MAX_BOUND = 1 + ALPHA_BOUND_ERR;

if(alpha_debug)
{
    osAlpha<< endl;
    forAll(alpha1Cells, cellI)
    {
        if(alpha1Cells[cellI] < ALPHA_MIN_BOUND || alpha1Cells[cellI] > ALPHA_MAX_BOUND)
        {
            osAlpha<< "Unbounded alpha in cell " << cellI << " alpha1 = " << alpha1Cells[cellI] << endl;
            if(alpha1Cells[cellI] < ALPHA_MIN_BOUND)
            {
                osAlpha<< "alpha1 undershoot = " << alpha1Cells[cellI] << endl;
            }
            else
            {
                osAlpha<< "alpha1 overshoot = " << (alpha1Cells[cellI] - 1) << endl;
            }

            div_alpha1 = 0;
            div_phi = 0;

            forAll(own, faceI)
            {
                if(own[faceI] == cellI)
                {
                    osAlpha<< "Face " << faceI << " own  phiAlpha1 = " << phiAlpha1[faceI] << "  phi = " << phi[faceI] << endl;
                    div_alpha1 += phiAlpha1[faceI];
                    div_phi += phi[faceI];
                }
                else if(nei[faceI] == cellI)
                {
                    osAlpha<< "Face " << faceI << " nei  phiAlpha1 = " << -phiAlpha1[faceI] << "  phi = " << -phi[faceI] << endl;
                    div_alpha1 -= phiAlpha1[faceI];
                    div_phi -= phi[faceI];
                }
            }

            forAll(alpha1.boundaryField(), patchI)
            {
                const labelUList& pFaceCells = mesh.boundary()[patchI].faceCells();
                const fvsPatchScalarField& pphiAlpha1 = phiAlpha1.boundaryField()[patchI];
                const fvsPatchScalarField& pphi = phi.boundaryField()[patchI];

                forAll(pphi, faceI)
                {
                    if(pFaceCells[faceI] == cellI)
                    {
                        osAlpha<< "Patch " << patchI << " face " << faceI << "  phiAlpha1 = " << pphiAlpha1[faceI] << "  phi = " << pphi[faceI] << endl;
                        div_alpha1 += pphiAlpha1[faceI];
                        div_phi += pphi[faceI];
                    }
                }
            }

            osAlpha<< "alpha1 old time = " << alpha1.oldTime().internalField()[cellI] << "  alpha1Old = " << alpha1OldCells[cellI] << nl
                << "div(phi)*V = " << div_phi
                << "  div(phi, alpha1)*V = " << div_alpha1 << endl;

            div_alpha1 /= meshV[cellI];
            div_phi /= meshV[cellI];
            osAlpha<< "div(phi) = " << div_phi
                << "  div(phi, alpha1) = " << div_alpha1 << nl
                << "pS = " << pSCells[cellI] << nl
                << "div(phi, alpha1)*dt = " << div_alpha1*dt
                << "  alpha1 new time = " << alpha1Cells[cellI] << nl
                << endl;
        }
    }
}

if(redist_alpha)
{
    plicFuncs::redistribute_alpha_field(mesh, alpha1Cells, cell_stencil.stencil(), ALPHA_MIN_BOUND, ALPHA_BOUND_ITERS_MAX, minAlpha1Cell, maxAlpha1Cell, minAlpha1, maxAlpha1, alpha_debug, osAlpha);

    alpha1.correctBoundaryConditions();

    Info<< "alpha1 after redistribution:  avg = "
    << alpha1.weightedAverage(mesh.V()).value()
    << "  min = " << min(alpha1).value()
    << "  undershoot = " << min(min(alpha1).value(), 0)
    << "  max = " << max(alpha1).value()
    << "  overshoot = " << max((max(alpha1).value() - 1), 0)
    << endl;
}

if(alpha_debug)
{    
    minAlpha1CellC = vector::zero;
    maxAlpha1CellC = vector::zero;
    if(minAlpha1Cell >= 0)
    {
        minAlpha1CellC = meshCi[minAlpha1Cell];
    }
    if(maxAlpha1Cell >= 0)
    {
        maxAlpha1CellC = meshCi[maxAlpha1Cell];
    }

    osAlpha<< "alpha1 extrema locations: " << nl
        << "Min alpha1 cell: " << minAlpha1Cell << " Min cell loc: " << minAlpha1CellC << " Min alpha1 value = " << minAlpha1 << nl 
        << "Max alpha1 cell: " << maxAlpha1Cell  << " Max cell loc: " << maxAlpha1CellC << " Max alpha1 value = " << maxAlpha1 << nl
        << endl;
}

alpha1.max(0.0);
alpha1.min(1.0);

Info<< "alpha1 after clipping:  avg = "
    << alpha1.weightedAverage(mesh.V()).value()    
    << endl;

alpha1.correctBoundaryConditions();

limitedAlpha1 = min(max(alpha1, scalar(0)), scalar(1));

alpha0 = scalar(1) - limitedAlpha1;
