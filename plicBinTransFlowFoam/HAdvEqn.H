
if(adv_debug2)
{
    plicFuncs::print_line(osAdv, 80);
    osAdv<< "Solution Time = " << runTime.timeName() << "  H Adv" << endl;
    plicFuncs::print_line(osAdv, 80);
    osAdv<< endl;
}

dt = deltaT;

//Phase-1
volScalarField surfInt_advFlux_h1(fvc::surfaceIntegrate(advFlux_rhoh1));
scalarField& surfInt_advFlux_h1Cells = surfInt_advFlux_h1.internalField();

if(adv_debug2)
{
    plicFuncs::print_advFluxFld(mesh, advFlux_rhoh1, alpha1, H1, h1, T1, "h1", "T1", osAdv);
    plicFuncs::print_advFluxIntData(mesh, advFlux_rhoh1, surfInt_advFlux_h1Cells, dt, alpha1, H1, h1, T1, Af_ph1_own, Af_ph1_nei, "ph1", "alpha1", "H1", "h1", "T1", osAdv);

    plicFuncs::print_line(osAdv, 80);
    osAdv<< setw(6) << "cellI" << "  " << setw(12) << "H1" << "  " << setw(12) << "rhoh1" << "  " << setw(12) << "h1" << "  " << setw(12) << "T1" << endl;
    plicFuncs::print_line(osAdv, 80);
}

forAll(H1Cells,cellI)
{                    
    H1Cells[cellI] = H1OldCells[cellI] - surfInt_advFlux_h1Cells[cellI]*dt;

    if(alpha1Cells[cellI] < SMALLEST_ALPHA_FOR_YADV)
    {
        T1Cells[cellI] = T1_TMP_VAL;
    }
    else
    {
        rhoh1Cells[cellI] = H1Cells[cellI]/alpha1Cells[cellI];
        h1Cells[cellI] = rhoh1Cells[cellI]/rho1Cells[cellI];

        h_cellI_tmp = h1Cells[cellI];
        T_cellI_tmp = T1Cells[cellI];
        for(i=0; i<n; i++)
        {
            x_cellI_tmp[i] = x1[i].internalField()[cellI];
        }

        plicFuncs::calc_T_from_h(x_cellI_tmp, P_thermo_tmp, h_cellI_tmp, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, T_cellI_tmp, MAX_ITER_T, T_TOL, H_TOL, T1min, T1max);
        T1Cells[cellI] = T_cellI_tmp;

    }

    if(adv_debug2)
    {
        osAdv<< setw(6) << cellI << "  " << setw(12) << H1Cells[cellI] << "  " << setw(12) << rhoh1Cells[cellI] << "  " << setw(12) << h1Cells[cellI] << "  " << setw(12) << T1Cells[cellI] << endl;
    }
}

if(adv_debug2)
{
    plicFuncs::print_line(osAdv, 80);
}

Info<< "Phase-1 H = ";
Info << (H1).weightedAverage(mesh.Vsc()).value();
Info << nl;
Info << "Min(H1) = " << min(H1).value();
Info << "  Max(H1) = " << max(H1).value();
Info << nl;
Info << "Min(T1) = " << min(T1).value();
Info << "  Max(T1) = " << max(T1).value();
Info << nl;
Info << "H1 undershoot = " << min((min(H1).value() - H1min), 0);
Info << "  H1 overshoot = " << max((max(H1).value() - H1max), 0);
Info << nl;
Info << "T1 undershoot = " << min((min(T1).value() - T1min), 0);
Info << "  T1 overshoot = " << max((max(T1).value() - T1max), 0);
Info << endl;
        
forAll(T1Cells,cellI)
{        
    if(T1Cells[cellI] < T1min)
    {
        T1Cells[cellI] = T1min;                
    }
        
    if(T1Cells[cellI] > T1max)
    {
        T1Cells[cellI] = T1max;        
    }   

    h_cellI_tmp = h1Cells[cellI];
    T_cellI_tmp = T1Cells[cellI];
    for(i=0; i<n; i++)
    {
        x_cellI_tmp[i] = x1[i].internalField()[cellI];
    }
    plicFuncs::calc_h_from_T(x_cellI_tmp, P_thermo_tmp, T_cellI_tmp, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, h_cellI_tmp);

    h1Cells[cellI] = h_cellI_tmp;
    rhoh1Cells[cellI] = rho1Cells[cellI]*h_cellI_tmp;
    H1Cells[cellI] = alpha1Cells[cellI]*rhoh1Cells[cellI];
}

Info << "Done" << endl;

Info<< "Phase-1 H after clipping = "
<< (H1).weightedAverage(mesh.Vsc()).value()
<< nl
<< "Min(H1) = " << min(H1).value()
<< "  Max(H1) = " << max(H1).value()
<< nl
<< "Min(T1) = " << min(T1).value()
<< "  Max(T1) = " << max(T1).value()        
<< endl;

T1.correctBoundaryConditions();
plicFuncs::correct_boundaryField_h_rhoh_H(P_thermo_tmp, T1, x1, rho1, alpha1, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, h1, rhoh1, H1);

//Phase-0
volScalarField surfInt_advFlux_h0(fvc::surfaceIntegrate(advFlux_rhoh0));
scalarField& surfInt_advFlux_h0Cells = surfInt_advFlux_h0.internalField();

if(adv_debug2)
{
    plicFuncs::print_advFluxFld(mesh, advFlux_rhoh0, alpha0, H0, h0, T0, "h0", "T0", osAdv);
    plicFuncs::print_advFluxIntData(mesh, advFlux_rhoh0, surfInt_advFlux_h0Cells, dt, alpha0, H0, h0, T0, Af_ph0_own, Af_ph0_nei, "ph0", "alpha0", "H0", "h0", "T0", osAdv);

    plicFuncs::print_line(osAdv, 80);
    osAdv<< setw(6) << "cellI" << "  " << setw(12) << "H0" << "  " << setw(12) << "rhoh0" << "  " << setw(12) << "h0" << "  " << setw(12) << "T0" << endl;
    plicFuncs::print_line(osAdv, 80);
}

forAll(H0Cells,cellI)
{                    
    H0Cells[cellI] = H0OldCells[cellI] - surfInt_advFlux_h0Cells[cellI]*dt;

    if(alpha0Cells[cellI] < SMALLEST_ALPHA_FOR_YADV)
    {
        T0Cells[cellI] = T0_TMP_VAL;
    }
    else
    {
        rhoh0Cells[cellI] = H0Cells[cellI]/alpha0Cells[cellI];
        h0Cells[cellI] = rhoh0Cells[cellI]/rho0Cells[cellI];

        h_cellI_tmp = h0Cells[cellI];
        T_cellI_tmp = T0Cells[cellI];
        for(i=0; i<n; i++)
        {
            x_cellI_tmp[i] = x0[i].internalField()[cellI];
        }

        plicFuncs::calc_T_from_h(x_cellI_tmp, P_thermo_tmp, h_cellI_tmp, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, T_cellI_tmp, MAX_ITER_T, T_TOL, H_TOL, T0min, T0max);

        T0Cells[cellI] = T_cellI_tmp;

    }

    if(adv_debug2)
    {
        osAdv<< setw(6) << cellI << "  " << setw(12) << H0Cells[cellI] << "  " << setw(12) << rhoh0Cells[cellI] << "  " << setw(12) << h0Cells[cellI] << "  " << setw(12) << T0Cells[cellI] << endl;
    }
}

if(adv_debug2)
{
    plicFuncs::print_line(osAdv, 80);
}

Info<< "Phase-0 H = "
<< (H0).weightedAverage(mesh.Vsc()).value()
<< nl
<< "Min(H0) = " << min(H0).value()
<< "  Max(H0) = " << max(H0).value()
<< nl
<< "Min(T0) = " << min(T0).value()
<< "  Max(T0) = " << max(T0).value()
<< nl
<< "H0 undershoot = " << min((min(H0).value() - H0min), 0)
<< "  H0 overshoot = " << max((max(H0).value() - H0max), 0)
<< nl
<< "T0 undershoot = " << min((min(T0).value() - T0min), 0)
<< "  T0 overshoot = " << max((max(T0).value() - T0max), 0)
<< endl;
        
forAll(T0Cells,cellI)
{        
    if(T0Cells[cellI] < T0min)
    {
        T0Cells[cellI] = T0min;                
    }
        
    if(T0Cells[cellI] > T0max)
    {
        T0Cells[cellI] = T0max;        
    }   

    h_cellI_tmp = h0Cells[cellI];
    T_cellI_tmp = T0Cells[cellI];
    for(i=0; i<n; i++)
    {
        x_cellI_tmp[i] = x0[i].internalField()[cellI];
    }
    plicFuncs::calc_h_from_T(x_cellI_tmp, P_thermo_tmp, T_cellI_tmp, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, h_cellI_tmp);

    h0Cells[cellI] = h_cellI_tmp;
    rhoh0Cells[cellI] = rho0Cells[cellI]*h_cellI_tmp;
    H0Cells[cellI] = alpha0Cells[cellI]*rhoh0Cells[cellI];
}

Info<< "Phase-0 H after clipping = "
<< (H0).weightedAverage(mesh.Vsc()).value()
<< nl
<< "Min(H0) = " << min(H0).value()
<< "  Max(H0) = " << max(H0).value()
<< nl
<< "Min(T0) = " << min(T0).value()
<< "  Max(T0) = " << max(T0).value()        
<< endl;

T0.correctBoundaryConditions();
plicFuncs::correct_boundaryField_h_rhoh_H(P_thermo_tmp, T0, x0, rho0, alpha0, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, h0, rhoh0, H0);

Info<< "Done HAdvEqns" << endl;
