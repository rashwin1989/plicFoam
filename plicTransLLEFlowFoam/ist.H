
if(ist_debug)
{
    plicFuncs::print_line(osIst, 80);
    osIst<< "Solution Time = " << runTime.timeName() << "  Interfacial transfer" << endl;
    plicFuncs::print_line(osIst, 80);
    osIst<< endl;
}

dt = deltaT;

//-----------------------------------------------------------------------------------------------------//
//Make flat fields for alpha1 and phase centroids required for gradient calculation
plicFuncs::makeFlatFld(alpha1, istMap, alpha1_flatFld_ist);
plicFuncs::makeFlatFld(interface.C_ph1(), istMap, C_ph1_flatFld_ist);
plicFuncs::makeFlatFld(interface.C_ph0(), istMap, C_ph0_flatFld_ist);
plicFuncs::makeFlatFld(T1, istMap, T1_flatFld_ist);
plicFuncs::makeFlatFld(T0, istMap, T0_flatFld_ist);

for(label i=0; i<nSpecies; i++)
{
    const volScalarField& x1i = x1[i];
    const volScalarField& x0i = x0[i];
    List<scalar>& x1i_flatFld_ist = x1_flatFld_ist[i];
    List<scalar>& x0i_flatFld_ist = x0_flatFld_ist[i];
    plicFuncs::makeFlatFld(x1i, istMap, x1i_flatFld_ist);
    plicFuncs::makeFlatFld(x0i, istMap, x0i_flatFld_ist);
}
//Done make flat fields for alpha1 and phase centroids required for gradient calculation
//-----------------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------------//
//Calculate interfacial species concentrations and mass fluxes
//Calculate interfacial species transfer source terms and alpha source terms
plicFuncs::calc_Xs_Ys_Js_mS_alphaS(mesh, Pc, Tc, Vc, w, MW, tk, coef_ab, Tb, SG, H8, k, dm, istCellStencil, x1_flatFld_ist, x0_flatFld_ist, T1_flatFld_ist, T0_flatFld_ist, alpha1_flatFld_ist, alpha1, rho1, rho0, T1, T0, Ts, P_thermo, C1, C0, x1, x0, Y1, Y0, C_ph1_flatFld_ist, C_ph0_flatFld_ist, C_intfc, A_intfc, nHat, dt, Ys1, Ys0, xs1, xs0, Js1, Js0, mS1, JsTot, mS1Tot, alphaS1, alphaS0, Qs, n_iters_Ts, status_transLLE, cell_had_intfc, n_flux_type, flux_umf, nSpecies, ALPHA_2PH_MIN, A_INTFC_2PH_MIN, Ts_TOL, MAX_ITERS_Ts, MASS_FRAC_TOL, ist_debug2, osIst);
//Done calculate interfacial species concentrations and mass fluxes
//Done calculate interfacial species transfer source terms and alpha source terms
//-----------------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------------//
//Apply interfacial species transfer and alpha source terms, calculate source term for pressure equation
forAll(pSCells, cellI)
{
    pSCoeff = 1.0/rho1Cells[cellI] - 1.0/rho0Cells[cellI];
    pSCells[cellI] = pSCoeff*mS1TotCells[cellI];
}

forAll(alpha1Cells, cellI)
{
    if(ist_debug)
    {    
        plicFuncs::print_line(osIst, 80);
        osIst<< "Cell " << cellI << endl;
        plicFuncs::print_line(osIst, 80);
        osIst<< "Before interfacial species transfer" << endl;
        plicFuncs::print_line(osIst, 80);
        osIst<< setw(12) << "alpha1" << "  " << setw(12) << "alpha0" << "  " << setw(16) << "alpha1+alpha0" << "  " << setw(12) << "rho1" << "  " << setw(12) << "rho0" << "  " << setw(12) << "T1" << "  " << setw(12) << "T0" << endl;
        plicFuncs::print_line(osIst, 80);
        osIst<< setw(12) << alpha1Cells[cellI] << "  " << setw(12) << alpha0Cells[cellI] << "  " << setw(16) << alpha1Cells[cellI] + alpha0Cells[cellI] << "  " << setw(12) << rho1Cells[cellI] << "  " << setw(12) << rho0Cells[cellI] << "  " << setw(12) << T1Cells[cellI] << "  " << setw(12) << T0Cells[cellI] << endl;
        plicFuncs::print_line(osIst, 80);
        osIst<< setw(7) << "Species" << "  " << setw(12) << "C1" << "  " << setw(12) << "C0" << "  " << setw(12) << "Y1" << "  " << setw(12) << "Y0" << endl;
        plicFuncs::print_line(osIst, 80);
        for(i=0; i<n; i++)
        {
            osIst<< setw(7) << i << "  " << setw(12) << C1[i].internalField()[cellI] << "  " << setw(12) << C0[i].internalField()[cellI] << "  " << setw(12) << Y1[i].internalField()[cellI] << "  " << setw(12) << Y0[i].internalField()[cellI] << endl;
        }
        plicFuncs::print_line(osIst, 80);
        osIst<< "Interface conditions" << endl;
        plicFuncs::print_line(osIst, 80);
        osIst<< "Ts = " << TsCells[cellI] << endl;
        plicFuncs::print_line(osIst, 80);
        osIst<< setw(7) << "Species" << "  " << setw(12) << "xs1" << "  " << setw(12) << "xs0" << "  " << setw(12) << "ys1" << "  " << setw(12) << "ys0" << endl;
        plicFuncs::print_line(osIst, 80);
        for(i=0; i<n; i++)
        {
            osIst<< setw(7) << i << "  " << setw(12) << xs1[i].internalField()[cellI] << "  " << setw(12) << xs0[i].internalField()[cellI] << "  " << setw(12) << Ys1[i].internalField()[cellI] << "  " << setw(12) << Ys0[i].internalField()[cellI] << endl;
        }
        plicFuncs::print_line(osIst, 80);
        osIst<< "Interfacial transfer source terms" << endl;
        plicFuncs::print_line(osIst, 80);
        osIst<< setw(12) << "mS1Tot" << "  " << setw(12) << "alphaS1" << "  " << setw(12) << "alphaS0" << "  " << setw(12) << "Qs" << "  " << setw(12) << "pS" << "  " << setw(12) << "dt" << endl;
        plicFuncs::print_line(osIst, 80);
        osIst<< setw(12) << mS1TotCells[cellI] << "  " << setw(12) << alphaS1Cells[cellI] << "  " << setw(12) << alphaS0Cells[cellI] << "  " << setw(12) << QsCells[cellI] << "  " << setw(12) << pSCells[cellI] << "  " << setw(12) << dt << endl;
        plicFuncs::print_line(osIst, 80);
        osIst<< setw(7) << "Species" << "  " << setw(12) << "mS1" << endl;
        plicFuncs::print_line(osIst, 80);
        for(i=0; i<n; i++)
        {
            osIst<< setw(7) << i << "  " << setw(12) << mS1[i].internalField()[cellI] << endl;
        }
        plicFuncs::print_line(osIst, 80);
    }
}

forAll(alpha1Cells, cellI)
{
    //ph-1
    //volume fraction
    alpha1Cells[cellI] += alphaS1Cells[cellI]*dt;    

    //species
    for(i=0; i<(n-1); i++)
    {
        C1[i].internalField()[cellI] += mS1[i].internalField()[cellI]*dt;              
    }

    Y1t_cellI = 0.0;

    for(i=0; i<(n-1); i++)
    {   
        if(alpha1Cells[cellI] < SMALLEST_ALPHA_FOR_YADV)
        {
            Y1i_cellI = 0;
        }
        else
        {
            Y1i_cellI = C1[i].internalField()[cellI]/alpha1Cells[cellI]/rho1Cells[cellI];
        }

        if(Y1i_cellI < 0)
        {
            Y1i_cellI = 0;            
        }
        
        if(Y1i_cellI > 1)
        {
            Y1i_cellI = 1;
        }
        
        Y1[i].internalField()[cellI] = Y1i_cellI;
        c1[i].internalField()[cellI] = rho1Cells[cellI]*Y1i_cellI;
        C1[i].internalField()[cellI] = alpha1Cells[cellI]*rho1Cells[cellI]*Y1i_cellI;

        Y1t_cellI += Y1i_cellI;
        y_cellI_tmp[i] = Y1i_cellI;
    }

    Y1n_cellI = max(0.0, 1.0 - Y1t_cellI);
    y_cellI_tmp[n-1] = Y1n_cellI;
    Y1[n-1].internalField()[cellI] = Y1n_cellI;
    c1[n-1].internalField()[cellI] = rho1Cells[cellI]*Y1n_cellI;
    C1[n-1].internalField()[cellI] = alpha1Cells[cellI]*rho1Cells[cellI]*Y1n_cellI;

    plicFuncs::y2x(n, MW, y_cellI_tmp, x_cellI_tmp);

    for(i=0; i<n; i++)
    {
        x1[i].internalField()[cellI] = x_cellI_tmp[i];
    }

    //enthalpy and temperature
    H1Cells[cellI] += QsCells[cellI]*dt;
    
    if(alpha1Cells[cellI] < SMALLEST_ALPHA_FOR_YADV)
    {
        T1Cells[cellI] = T1_TMP_VAL;
    }
    else
    {
        rhoh_cellI_tmp = H1Cells[cellI]/alpha1Cells[cellI];
        h_cellI_tmp = rhoh_cellI_tmp/rho1Cells[cellI];
        T_cellI_tmp = T1Cells[cellI];        
        plicFuncs::calc_T_from_h(x_cellI_tmp, P_thermo, h_cellI_tmp, n, Pc, Tc, w, MW, tk, Tb, SG, H8, T_cellI_tmp, MAX_ITER_T, T_TOL, H_TOL);
        if(T_cellI_tmp < T1min)
        {
            T_cellI_tmp = T1min;                
        }
        
        if(T_cellI_tmp > T1max)
        {
            T_cellI_tmp = T1max;        
        }        
    }
    plicFuncs::calc_h_from_T(x_cellI_tmp, P_thermo, T_cellI_tmp, n, Pc, Tc, w, MW, tk, Tb, SG, H8, h_cellI_tmp);
    
    T1Cells[cellI] = T_cellI_tmp;
    h1Cells[cellI] = h_cellI_tmp;
    rhoh1Cells[cellI] = rho1Cells[cellI]*h1_cellI_tmp;
    H1Cells[cellI] = alpha1Cells[cellI]*rhoh1Cells[cellI];

    //ph-0
    //volume fraction
    alpha0Cells[cellI] += alphaS0Cells[cellI]*dt;    

    //species
    for(i=0; i<(n-1); i++)
    {
        C0[i].internalField()[cellI] += -mS1[i].internalField()[cellI]*dt;              
    }

    Y0t_cellI = 0.0;

    for(i=0; i<(n-1); i++)
    {   
        if(alpha0Cells[cellI] < SMALLEST_ALPHA_FOR_YADV)
        {
            Y0i_cellI = 0;
        }
        else
        {
            Y0i_cellI = C0[i].internalField()[cellI]/alpha0Cells[cellI]/rho0Cells[cellI];
        }

        if(Y0i_cellI < 0)
        {
            Y0i_cellI = 0;            
        }
        
        if(Y0i_cellI > 1)
        {
            Y0i_cellI = 1;
        }
        
        Y0[i].internalField()[cellI] = Y0i_cellI;
        c0[i].internalField()[cellI] = rho0Cells[cellI]*Y0i_cellI;
        C0[i].internalField()[cellI] = alpha0Cells[cellI]*rho0Cells[cellI]*Y0i_cellI;

        Y0t_cellI += Y0i_cellI;
        y_cellI_tmp[i] = Y0i_cellI;
    }

    Y0n_cellI = max(0.0, 1.0 - Y0t_cellI);
    y_cellI_tmp[n-1] = Y0n_cellI;
    Y0[n-1].internalField()[cellI] = Y0n_cellI;
    c0[n-1].internalField()[cellI] = rho0Cells[cellI]*Y0n_cellI;
    C0[n-1].internalField()[cellI] = alpha0Cells[cellI]*rho0Cells[cellI]*Y0n_cellI;

    plicFuncs::y2x(n, MW, y_cellI_tmp, x_cellI_tmp);

    for(i=0; i<n; i++)
    {
        x0[i].internalField()[cellI] = x_cellI_tmp[i];
    }

    //enthalpy and temperature
    H0Cells[cellI] += -QsCells[cellI]*dt;
    
    if(alpha0Cells[cellI] < SMALLEST_ALPHA_FOR_YADV)
    {
        T0Cells[cellI] = T0_TMP_VAL;
    }
    else
    {
        rhoh_cellI_tmp = H0Cells[cellI]/alpha0Cells[cellI];
        h_cellI_tmp = rhoh_cellI_tmp/rho0Cells[cellI];
        T_cellI_tmp = T0Cells[cellI];        
        plicFuncs::calc_T_from_h(x_cellI_tmp, P_thermo, h_cellI_tmp, n, Pc, Tc, w, MW, tk, Tb, SG, H8, T_cellI_tmp, MAX_ITER_T, T_TOL, H_TOL);
        if(T_cellI_tmp < T0min)
        {
            T_cellI_tmp = T0min;                
        }
        
        if(T_cellI_tmp > T0max)
        {
            T_cellI_tmp = T0max;        
        }        
    }
    plicFuncs::calc_h_from_T(x_cellI_tmp, P_thermo, T_cellI_tmp, n, Pc, Tc, w, MW, tk, Tb, SG, H8, h_cellI_tmp);
    
    T0Cells[cellI] = T_cellI_tmp;
    h0Cells[cellI] = h_cellI_tmp;
    rhoh0Cells[cellI] = rho0Cells[cellI]*h_cellI_tmp;
    H0Cells[cellI] = alpha0Cells[cellI]*rhoh0Cells[cellI];
}

alpha1.correctBoundaryConditions();

for(i=0; i<n; i++)
{    
    Y1[i].correctBoundaryConditions();
    plicFuncs::correct_boundaryField_C(Y1[i], rho1, alpha1, c1[i], C1[i]);

    Y0[i].correctBoundaryConditions();
    plicFuncs::correct_boundaryField_C(Y0[i], rho0, alpha0, c0[i], C0[i]);

    Info<< "Phase-1 C[" << i << "] = "
        << (C1[i]).weightedAverage(mesh.Vsc()).value()
        << nl
        << "Min(C1i) = " << min(C1[i]).value()
        << "  Max(C1i) = " << max(C1[i]).value()
        << nl
        << "Min(Y1i) = " << min(Y1[i]).value()
        << "  Max(Y1i) = " << max(Y1[i]).value()
        << endl;

    Info<< "Phase-0 C[" << i << "] = "
        << (C0[i]).weightedAverage(mesh.Vsc()).value()
        << nl
        << "Min(C0i) = " << min(C0[i]).value()
        << "  Max(C0i) = " << max(C0[i]).value()
        << nl
        << "Min(Y0i) = " << min(Y0[i]).value()
        << "  Max(Y0i) = " << max(Y0[i]).value()
        << endl;        
}

T1.correctBoundaryConditions();
plicFuncs::correct_boundaryField_h_rhoh_H(n,Pc,Tc,w,MW,tk,Tb,SG,H8,P_thermo,x1,T1,rho1,alpha1,h1,rhoh1,H1);

T0.correctBoundaryConditions();
plicFuncs::correct_boundaryField_h_rhoh_H(n,Pc,Tc,w,MW,tk,Tb,SG,H8,P_thermo,x0,T0,rho0,alpha0,h0,rhoh0,H0);

forAll(alpha1Cells, cellI)
{
    if(ist_debug)
    {    
        plicFuncs::print_line(osIst, 80);
        osIst<< "Cell " << cellI << endl;
        plicFuncs::print_line(osIst, 80);
        osIst<< "After interfacial species transfer" << endl;
        plicFuncs::print_line(osIst, 80);
        osIst<< setw(12) << "alpha1" << "  " << setw(12) << "alpha0" << "  " << setw(16) << "alpha1+alpha0" << "  " << setw(12) << "rho1" << "  " << setw(12) << "rho0" << "  " << setw(12) << "T1" << "  " << setw(12) << "T0" << endl;
        plicFuncs::print_line(osIst, 80);
        osIst<< setw(12) << alpha1Cells[cellI] << "  " << setw(12) << alpha0Cells[cellI] << "  " << setw(16) << alpha1Cells[cellI] + alpha0Cells[cellI] << "  " << setw(12) << rho1Cells[cellI] << "  " << setw(12) << rho0Cells[cellI] << "  " << setw(12) << T1Cells[cellI] << "  " << setw(12) << T0Cells[cellI] << endl;
        plicFuncs::print_line(osIst, 80);
        osIst<< setw(7) << "Species" << "  " << setw(12) << "C1" << "  " << setw(12) << "C0" << "  " << setw(12) << "Y1" << "  " << setw(12) << "Y0" << endl;
        plicFuncs::print_line(osIst, 80);
        for(i=0; i<n; i++)
        {
            osIst<< setw(7) << i << "  " << setw(12) << C1[i].internalField()[cellI] << "  " << setw(12) << C0[i].internalField()[cellI] << "  " << setw(12) << Y1[i].internalField()[cellI] << "  " << setw(12) << Y0[i].internalField()[cellI] << endl;
        }        
        plicFuncs::print_line(osIst, 80);        
    }
}

//Done apply interfacial species transfer and alpha source terms, calculate source term for pressure equation
//-----------------------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------------------//
//Reset the allowable limits for concentrations and mass fractions in the phases
/*
for(label i=0; i<nSpecies; i++)
{
    C1max[i] = max(C1[i]).value();
    C0max[i] = max(C0[i]).value();
    C1min[i] = min(C1[i]).value();
    C0min[i] = min(C0[i]).value();

    Y1max[i] = max(Y1[i]).value();
    Y0max[i] = max(Y0[i]).value();
    Y1min[i] = min(Y1[i]).value();
    Y0min[i] = min(Y0[i]).value();
}
    */

Info<< "Done Interfacial Species Transfer" << endl;
