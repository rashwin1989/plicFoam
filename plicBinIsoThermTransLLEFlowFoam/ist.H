
if(ist_debug)
{
    plicFuncs::print_line(osIst, 100);
    osIst<< "Solution Time = " << runTime.timeName() << "  Interfacial transfer" << endl;
    plicFuncs::print_line(osIst, 100);
    osIst<< endl;
}

dt = deltaT;

//-----------------------------------------------------------------------------------------------------//
//Make flat fields for alpha1 and phase centroids required for gradient calculation
plicFuncs::makeFlatFld(alpha1, istMap, alpha1_flatFld_ist);
plicFuncs::makeFlatFld(interface.C_ph1(), istMap, C_ph1_flatFld_ist);
plicFuncs::makeFlatFld(interface.C_ph0(), istMap, C_ph0_flatFld_ist);

for(i=0; i<n; i++)
{
    const volScalarField& Y1i = Y1[i];
    const volScalarField& Y0i = Y0[i];
    List<scalar>& Y1i_flatFld_ist = Y1_flatFld_ist[i];
    List<scalar>& Y0i_flatFld_ist = Y0_flatFld_ist[i];
    plicFuncs::makeFlatFld(Y1i, istMap, Y1i_flatFld_ist);
    plicFuncs::makeFlatFld(Y0i, istMap, Y0i_flatFld_ist);
}
//Done make flat fields for alpha1 and phase centroids required for gradient calculation
//-----------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------//
//Calculate interfacial species mass fluxes
plicFuncs::calc_Js
(
    mesh,
    istCellStencil,
    Y1_flatFld_ist,
    Y0_flatFld_ist,
    alpha1_flatFld_ist,
    alpha1, 
    rho1,
    rho0,
    D1,
    D0,
    C_ph1_flatFld_ist,
    C_ph0_flatFld_ist,
    C_intfc,
    A_intfc,
    nHat,
    Ys1,
    Ys0,
    Js1,
    Js0,
    grads1,
    grads0,
    delta1,
    delta0,
    nSpecies,
    ALPHA_2PH_MIN_IST,
    erf_a,
    erf_b,
    nErf,
    erfInv_table,
    useErf,
    Y1_0,
    Y0_0,
    ist_debug2,
    osIst
);
//Done calculate interfacial species mass fluxes
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Calculate interfacial species transfer source terms and alpha source terms
plicFuncs::calc_mS_alphaS
(
    mesh,
    C1,
    C0,
    Y1,
    Y0,
    alpha1,
    rho1,
    rho0,
    Ys1,
    Ys0,
    Js1,
    Js0,
    nSpecies,
    ALPHA_2PH_MIN_IST,
    mS1,
    mS0,
    mS1Tot,
    mS0Tot,
    alphaS1,
    alphaS0,
    dt,
    ist_debug,
    osIst
);
//Done calculate interfacial species transfer source terms and alpha source terms
//------------------------------------------------------------------//


//-----------------------------------------------------------------------------------------------------//
//Apply interfacial species transfer and alpha source terms, calculate source term for pressure equation
forAll(pSCells, cellI)
{
    pSCoeff = 1.0/rho1Cells[cellI] - 1.0/rho0Cells[cellI];
    pSCells[cellI] = pSCoeff*mS1TotCells[cellI];
}

forAll(alpha1Cells, cellI)
{
    alpha1_cellI = alpha1Cells[cellI];
    
    if(alpha1_cellI > ALPHA_2PH_MIN_IST && alpha1_cellI < ALPHA_2PH_MAX_IST)
    {        
        calc_ist_for_cell[cellI] = 1;
    }
    else
    {
        calc_ist_for_cell[cellI] = 0;
    }
}
/*
for(faceI=0; faceI<nIntFaces; faceI++)
{
    label faceOwn = own[faceI];
    label faceNei = nei[faceI];

    if((alpha1Cells[faceOwn] >= ALPHA_2PH_MAX_IST && alpha1Cells[faceNei] <= ALPHA_2PH_MIN_IST) || (alpha1Cells[faceNei] >= ALPHA_2PH_MAX_IST && alpha1Cells[faceOwn] <= ALPHA_2PH_MIN_IST))
    {
        calc_ist_for_cell[faceOwn] = 1;
        calc_ist_for_cell[faceNei] = 1;
    } 
}

const polyBoundaryMesh& patches = mesh.boundaryMesh();

forAll(alpha1.boundaryField(), patchI)
{
    const polyPatch& pp = patches[patchI];
    const fvPatchScalarField& pAlpha1 = alpha1.boundaryField()[patchI];

    if(pp.coupled())
    {
        faceI = pp.start();
        const scalarField& alpha1NeiFld = alpha1.boundaryField()[patchI].patchNeighbourField();

        forAll(pAlpha1, fcI)
        {
            label faceOwn = own[faceI];
            scalar alpha1Own = alpha1Cells[faceOwn];
            scalar alpha1Nei = alpha1NeiFld[fcI];
            
            if((alpha1Own >= ALPHA_2PH_MAX_IST && alpha1Nei <= ALPHA_2PH_MIN_IST) || (alpha1Nei >= ALPHA_2PH_MAX_IST && alpha1Own <= ALPHA_2PH_MIN_IST))
            {
                calc_ist_for_cell[faceOwn] = 1;
            }
        }
    }
}
*/
forAll(alpha1Cells, cellI)
{
    if(calc_ist_for_cell[cellI] == 1)
    {
        if(ist_debug)
        {    
            plicFuncs::print_line(osIst, 100);
            osIst<< "Cell " << cellI << endl;
            plicFuncs::print_line(osIst, 100);
            osIst<< "Before interfacial species transfer" << endl;
            plicFuncs::print_line(osIst, 100);
            osIst<< setw(12) << "alpha1" << "  " << setw(12) << "alpha0" << "  " << setw(16) << "alpha1+alpha0" << "  " << setw(12) << "rho1" << "  " << setw(12) << "rho0" << endl;
            plicFuncs::print_line(osIst, 100);
            osIst<< setw(12) << alpha1Cells[cellI] << "  " << setw(12) << alpha0Cells[cellI] << "  " << setw(16) << alpha1Cells[cellI] + alpha0Cells[cellI] << "  " << setw(12) << rho1Cells[cellI] << "  " << setw(12) << rho0Cells[cellI] << endl;
            plicFuncs::print_line(osIst, 100);
            osIst<< setw(7) << "Species" << "  " << setw(12) << "C1" << "  " << setw(12) << "C0" << "  " << setw(12) << "Y1" << "  " << setw(12) << "Y0" << endl;
            plicFuncs::print_line(osIst, 100);
            for(i=0; i<n; i++)
            {
                osIst<< setw(7) << i << "  " << setw(12) << C1[i].internalField()[cellI] << "  " << setw(12) << C0[i].internalField()[cellI] << "  " << setw(12) << Y1[i].internalField()[cellI] << "  " << setw(12) << Y0[i].internalField()[cellI] << endl;
            }
            plicFuncs::print_line(osIst, 100);
            osIst<< "Interface conditions" << endl;
            plicFuncs::print_line(osIst, 100);            
            osIst<< setw(7) << "Species" << "  " << setw(12) << "ys1" << "  " << setw(12) << "ys0" << endl;
            plicFuncs::print_line(osIst, 100);
            for(i=0; i<n; i++)
            {
                osIst<< setw(7) << i << "  " << setw(12) << Ys1[i].internalField()[cellI] << "  " << setw(12) << Ys0[i].internalField()[cellI] << endl;
            }
            plicFuncs::print_line(osIst, 100);
            osIst<< "Interfacial transfer source terms" << endl;
            plicFuncs::print_line(osIst, 100);
            osIst<< setw(12) << "mS1Tot" << "  " << setw(12) << "alphaS1" << "  " << setw(12) << "alphaS0" << "  " << setw(12) << "pS" << "  " << setw(12) << "dt" << endl;
            plicFuncs::print_line(osIst, 100);
            osIst<< setw(12) << mS1TotCells[cellI] << "  " << setw(12) << alphaS1Cells[cellI] << "  " << setw(12) << alphaS0Cells[cellI] << "  " << setw(12) << pSCells[cellI] << "  " << setw(12) << dt << endl;
            plicFuncs::print_line(osIst, 100);
            osIst<< setw(7) << "Species" << "  " << setw(12) << "mS1" << endl;
            plicFuncs::print_line(osIst, 100);
            for(i=0; i<n; i++)
            {
                osIst<< setw(7) << i << "  " << setw(12) << mS1[i].internalField()[cellI] << endl;
            }
            plicFuncs::print_line(osIst, 100);
        }
    }

    
    //ph-1
    //volume fraction
    alpha1_cellI = alpha1Cells[cellI] + alphaS1Cells[cellI]*dt;
    if(alpha1_cellI < 0)
    {
        osIst<< "Time = " << runTime.timeName() << nl
            << "alpha1 undershoot in cell " << cellI << "  alpha1 = " << alpha1_cellI << endl;
    }
    if(alpha1_cellI > 1)
    {
        osIst<< "Time = " << runTime.timeName() << nl
            << "alpha1 overshoot in cell " << cellI << "  alpha1 = " << alpha1_cellI << endl;
    }    
    alpha1_cellI = max(alpha1_cellI, 0.0);
    alpha1Cells[cellI] = alpha1_cellI;

    C1Tot_cellI = 0.0;
    //species
    for(i=0; i<n; i++)
    {
        C1i_cellI = C1[i].internalField()[cellI] + mS1[i].internalField()[cellI]*dt;
        // if(C1i_cellI < 0)
        // {
        //     osIst<< "Time = " << runTime.timeName() << nl
        //         << "C1[" << i << "] -ve in cell " << cellI << "  C1[" << i << "] = " << C1i_cellI << endl;
        // }
        
        C1i_cellI = max(C1i_cellI, 0.0);
        C1[i].internalField()[cellI] = C1i_cellI;
        C1Tot_cellI += C1i_cellI;
    }

    if(C1Tot_cellI < SMALL)
    {
        C1Tot_cellI += SMALL;
        C1[0].internalField()[cellI] += SMALL;
    }

    Y1t_cellI = 0.0;

    for(i=0; i<(n-1); i++)
    {   
        Y1i_cellI = C1[i].internalField()[cellI]/C1Tot_cellI;
        Y1[i].internalField()[cellI] = Y1i_cellI;
        c1[i].internalField()[cellI] = rho1Cells[cellI]*Y1i_cellI;
        Y1t_cellI += Y1i_cellI;
        y_cellI_tmp[i] = Y1i_cellI;
    }

    Y1n_cellI = max(0.0, 1.0 - Y1t_cellI);
    y_cellI_tmp[n-1] = Y1n_cellI;
    Y1[n-1].internalField()[cellI] = Y1n_cellI;
    c1[n-1].internalField()[cellI] = rho1Cells[cellI]*Y1n_cellI;    
    
    plicFuncs::y2x(n, MW, y_cellI_tmp, x_cellI_tmp);

    for(i=0; i<n; i++)
    {
        x1[i].internalField()[cellI] = x_cellI_tmp[i];
    }

    //ph-0
    //volume fraction
    alpha0_cellI = alpha0Cells[cellI] + alphaS0Cells[cellI]*dt;
    if(alpha0_cellI < 0)
    {
        osIst<< "Time = " << runTime.timeName() << nl
            << "alpha0 undershoot in cell " << cellI << "  alpha0 = " << alpha0_cellI << endl;
    }
    if(alpha0_cellI > 1)
    {
        osIst<< "Time = " << runTime.timeName() << nl
            << "alpha0 overshoot in cell " << cellI << "  alpha0 = " << alpha0_cellI << endl;
    }
    if(alpha0_cellI > 1.5)
    {
        plicFuncs::print_line(osIst, 100);
        osIst<< "Time = " << runTime.timeName() << endl;
        plicFuncs::print_line(osIst, 100);
        plicFuncs::print_line(osIst, 100);
        osIst<< "Cell " << cellI << endl;
        plicFuncs::print_line(osIst, 100);
        osIst<< "Before interfacial species transfer" << endl;
        plicFuncs::print_line(osIst, 100);
        osIst<< setw(14) << "alpha1" << setw(14) << "alpha0" << setw(14) << "alpha1+alpha0" << setw(14) << "rho1" << setw(14) << "rho0" << endl;
        plicFuncs::print_line(osIst, 100);
        osIst<< setw(14) << alpha1Cells[cellI] << setw(14) << alpha0Cells[cellI] << setw(16) << alpha1Cells[cellI] + alpha0Cells[cellI] << setw(14) << rho1Cells[cellI] << setw(14) << rho0Cells[cellI] << endl;
        plicFuncs::print_line(osIst, 100);
        osIst<< setw(9) << "Species" << setw(14) << "C1" << setw(14) << "C0" << setw(14) << "Y1" << setw(14) << "Y0" << endl;
        plicFuncs::print_line(osIst, 100);
        for(i=0; i<n; i++)
        {
            osIst<< setw(9) << i << setw(14) << C1[i].internalField()[cellI] << setw(14) << C0[i].internalField()[cellI] << setw(14) << Y1[i].internalField()[cellI] << setw(14) << Y0[i].internalField()[cellI] << endl;
        }
        plicFuncs::print_line(osIst, 100);
        osIst<< "Interface conditions" << endl;
        plicFuncs::print_line(osIst, 100);            
        osIst<< setw(9) << "Species" << setw(14) << "ys1" << setw(14) << "ys0" << endl;
        plicFuncs::print_line(osIst, 100);
        for(i=0; i<n; i++)
        {
            osIst<< setw(9) << i << setw(14) << Ys1[i].internalField()[cellI] << setw(14) << Ys0[i].internalField()[cellI] << endl;
        }
        plicFuncs::print_line(osIst, 100);
        osIst<< "Interfacial transfer source terms" << endl;
        plicFuncs::print_line(osIst, 100);
        osIst<< setw(14) << "mS1Tot" << setw(14) << "alphaS1" << setw(14) << "alphaS0" << setw(14) << "pS" << setw(14) << "dt" << endl;
        plicFuncs::print_line(osIst, 100);
        osIst<< setw(14) << mS1TotCells[cellI] << setw(14) << alphaS1Cells[cellI] << setw(14) << alphaS0Cells[cellI] << setw(14) << pSCells[cellI] << setw(14) << dt << endl;
        plicFuncs::print_line(osIst, 100);
        osIst<< "D1 = " << D1[0].internalField()[cellI] << "  D0 = " << D1[0].internalField()[cellI] << "  A_intfc = " << A_intfc.internalField()[cellI] << endl;
        plicFuncs::print_line(osIst, 100);
        osIst<< setw(9) << "Species" << setw(14) << "mS1" << setw(14) << "Js1" << setw(14) << "Js0" << endl;
        plicFuncs::print_line(osIst, 100);
        for(i=0; i<n; i++)
        {
            osIst<< setw(9) << i << setw(14) << mS1[i].internalField()[cellI] << setw(14) << Js1[i].internalField()[cellI] << setw(14) << Js0[i].internalField()[cellI] << endl;
        }
        plicFuncs::print_line(osIst, 100);        
    }
    alpha0_cellI = max(alpha0_cellI, 0.0);
    alpha0Cells[cellI] = alpha0_cellI;
    //alpha0Cells[cellI] = scalar(1.0) - alpha1Cells[cellI];

    C0Tot_cellI = 0.0;
    //species
    for(i=0; i<n; i++)
    {
        C0i_cellI = C0[i].internalField()[cellI] + mS0[i].internalField()[cellI]*dt;
        // if(C0i_cellI < 0)
        // {
        //     osIst<< "Time = " << runTime.timeName() << nl
        //         << "C0[" << i << "] -ve in cell " << cellI << "  C0[" << i << "] = " << C0i_cellI << endl;
        // }
        C0i_cellI = max(C0i_cellI, 0.0);
        C0[i].internalField()[cellI] = C0i_cellI;
        C0Tot_cellI += C0i_cellI;
    }

    if(C0Tot_cellI < SMALL)
    {
        C0Tot_cellI += SMALL;
        C0[n-1].internalField()[cellI] += SMALL;
    }

    Y0t_cellI = 0.0;
    for(i=0; i<(n-1); i++)
    {   
        Y0i_cellI = C0[i].internalField()[cellI]/C0Tot_cellI;
        Y0[i].internalField()[cellI] = Y0i_cellI;
        c0[i].internalField()[cellI] = rho0Cells[cellI]*Y0i_cellI;
        Y0t_cellI += Y0i_cellI;
        y_cellI_tmp[i] = Y0i_cellI;
    }

    Y0n_cellI = max(0.0, 1.0 - Y0t_cellI);
    y_cellI_tmp[n-1] = Y0n_cellI;
    Y0[n-1].internalField()[cellI] = Y0n_cellI;
    c0[n-1].internalField()[cellI] = rho0Cells[cellI]*Y0n_cellI;

    plicFuncs::y2x(n, MW, y_cellI_tmp, x_cellI_tmp);

    for(i=0; i<n; i++)
    {
        x0[i].internalField()[cellI] = x_cellI_tmp[i];
    }

    
    if(calc_ist_for_cell[cellI] == 1)
    {
        if(ist_debug)
        {            
            osIst<< "After interfacial species transfer" << endl;
            plicFuncs::print_line(osIst, 100);
            osIst<< setw(12) << "alpha1" << "  " << setw(12) << "alpha0" << "  " << setw(16) << "alpha1+alpha0" << "  " << setw(12) << "rho1" << "  " << setw(12) << "rho0" << endl;
            plicFuncs::print_line(osIst, 100);
            osIst<< setw(12) << alpha1Cells[cellI] << "  " << setw(12) << alpha0Cells[cellI] << "  " << setw(16) << alpha1Cells[cellI] + alpha0Cells[cellI] << "  " << setw(12) << rho1Cells[cellI] << "  " << setw(12) << rho0Cells[cellI] << endl;
            plicFuncs::print_line(osIst, 100);
            osIst<< setw(7) << "Species" << "  " << setw(12) << "C1" << "  " << setw(12) << "C0" << "  " << setw(12) << "Y1" << "  " << setw(12) << "Y0" << endl;
            plicFuncs::print_line(osIst, 100);
            for(i=0; i<n; i++)
            {
                osIst<< setw(7) << i << "  " << setw(12) << C1[i].internalField()[cellI] << "  " << setw(12) << C0[i].internalField()[cellI] << "  " << setw(12) << Y1[i].internalField()[cellI] << "  " << setw(12) << Y0[i].internalField()[cellI] << endl;
            }        
            plicFuncs::print_line(osIst, 100);        
        }
    }
}

Info<< "alpha1 avg = " << alpha1.weightedAverage(mesh.V()).value()
    << "  alpha1 min = " << min(alpha1).value()
    << "  alpha1 undershoot = " << min(min(alpha1).value(), 0)
    << "  alpha1 max = " << max(alpha1).value()
    << "  alpha1 overshoot = " << max((max(alpha1).value() - 1), 0)
    << endl;

Info<< "alpha0 avg = " << alpha0.weightedAverage(mesh.V()).value()
    << "  alpha0 min = " << min(alpha0).value()
    << "  alpha0 undershoot = " << min(min(alpha0).value(), 0)
    << "  alpha0 max = " << max(alpha0).value()
    << "  alpha0 overshoot = " << max((max(alpha0).value() - 1), 0)
    << endl;

alpha1.correctBoundaryConditions();
alpha0.correctBoundaryConditions();

for(i=0; i<n; i++)
{    
    volScalarField& Y1i = Y1[i];
    volScalarField& Y0i = Y0[i];
    volScalarField& c1i = c1[i];
    volScalarField& c0i = c0[i];
    volScalarField& C1i = C1[i];
    volScalarField& C0i = C0[i];

    Y1i.correctBoundaryConditions();
    plicFuncs::correct_boundaryField_C(Y1i, rho1, alpha1, c1i, C1i);

    Y0i.correctBoundaryConditions();
    plicFuncs::correct_boundaryField_C(Y0i, rho0, alpha0, c0i, C0i);
}
for(i=0; i<(n-1); i++)
{
    volScalarField& Y1i = Y1[i];
    volScalarField& Y0i = Y0[i];
    volScalarField& C1i = C1[i];
    volScalarField& C0i = C0[i];

    Info<< "C1[" << i << "] avg = " << C1i.weightedAverage(mesh.V()).value()
        << "  C1[" << i << "] min = " << min(C1i).value()
        << "  C1[" << i << "] max = " << max(C1i).value()
        << nl
        << "Y1[" << i << "] min = " << min(Y1i).value()
        << "  Y1[" << i << "] max = " << max(Y1i).value()
        << endl;

    Info<< "C0[" << i << "] avg = " << C0i.weightedAverage(mesh.V()).value()
        << "  C0[" << i << "] min = " << min(C0i).value()
        << "  C0[" << i << "] max = " << max(C0i).value()
        << nl
        << "Y0[" << i << "] min = " << min(Y0i).value()
        << "  Y0[" << i << "] max = " << max(Y0i).value()
        << endl;
}
//Done apply interfacial species transfer and alpha source terms, calculate source term for pressure equation
//-----------------------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------------------//
//Reset the allowable limits for concentrations and mass fractions in the phases


for(i=0; i<n; i++)
{
    C1max[i] = max(C1[i]).value();
    C0max[i] = max(C0[i]).value();
    C1min[i] = min(C1[i]).value();
    C0min[i] = min(C0[i]).value();
}


Info<< "Done Interfacial Species Transfer Source Terms Calculation" << endl;
