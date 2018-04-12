
if(diff_debug)
{
    plicFuncs::print_line(osDiff, 80);
    osDiff<< "Solution Time = " << runTime.timeName() << "  H Diff" << endl;
    plicFuncs::print_line(osDiff, 80);
    osDiff<< endl;
}

//set fractional time-step for diffusion
dt = 0.5*deltaT;

//------------------------------------------------------------------//
//Calculate enthalpy diffusive fluxes for both phases
plicFuncs::calc_2ph_diffFluxes_T
(
    mesh,    
    R_gas,
    n,
    Pc,
    Tc,
    w,
    MW,
    tk,
    Tb,
    SG,
    H8,
    P_thermo_tmp,
    lambda1f,
    lambda0f,
    gradf_T1,
    gradf_T0,
    Af_ph1_own,
    Af_ph0_own,
    Af_ph1_nei,
    Af_ph0_nei,
    face_phaseState_diff,
    rho1,
    alpha1,
    x1,
    H1,
    T1,
    rho0,
    alpha0,
    x0,
    H0,
    T0,
    dt,
    diffFlux_h1,
    diffFlux_h0,
    MAX_ITER_T,
    T_TOL,
    H_TOL,
    T1min,
    T1max,
    T0min,
    T0max,
    diff_debug2,
    osDiff    
);
//Done calculate enthalpy diffusive fluxes for both phases
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Calculate enthalpy carried by species diffusive fluxes for both phases
diffFlux_Jh1 = 0.0*diffFlux_h1;
diffFlux_Jh0 = 0.0*diffFlux_h0;
for(i=0; i<n; i++)
{
    diffFlux_Jh1 += diffFlux_Y1[i]*hpar1f[i];
    diffFlux_Jh0 += diffFlux_Y0[i]*hpar0f[i];
}
//Done calculate enthalpy carried by species diffusive fluxes for both phases
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Solve enthalpy diffusion equations in both phases

scalarField& T1Cells = T1.internalField();
scalarField& T0Cells = T0.internalField();
scalarField& H1Cells = H1.internalField();
scalarField& H0Cells = H0.internalField();
scalarField& rhoh1Cells = rhoh1.internalField();
scalarField& rhoh0Cells = rhoh0.internalField();
scalarField& h1Cells = h1.internalField();
scalarField& h0Cells = h0.internalField();

//Phase-1
volScalarField surfInt_diffFlux_h1(fvc::surfaceIntegrate(diffFlux_h1));
scalarField& surfInt_diffFlux_h1Cells = surfInt_diffFlux_h1.internalField();

volScalarField surfInt_diffFlux_Jh1(fvc::surfaceIntegrate(diffFlux_Jh1));
scalarField& surfInt_diffFlux_Jh1Cells = surfInt_diffFlux_Jh1.internalField();

if(diff_debug2)
{
    plicFuncs::print_diffFluxFld(mesh, diffFlux_h1, alpha1, H1, h1, T1, "h1", "T1", osDiff);
    plicFuncs::print_diffFluxIntData(mesh, diffFlux_h1, surfInt_diffFlux_h1Cells, dt, alpha1, H1, h1, T1, Af_ph1_own, Af_ph1_nei, "ph1", "alpha1", "H1", "h1", "T1", osDiff);

    plicFuncs::print_line(osDiff, 80);
    osDiff<< setw(6) << "cellI" << "  " << setw(12) << "H1" << "  " << setw(12) << "rhoh1" << "  " << setw(12) << "h1" << "  " << setw(12) << "T1" << endl;
    plicFuncs::print_line(osDiff, 80);
}

forAll(H1Cells,cellI)
{                    
    H1Cells[cellI] += -(surfInt_diffFlux_h1Cells[cellI] + surfInt_diffFlux_Jh1Cells[cellI])*dt;

    if(alpha1Cells[cellI] < SMALLEST_ALPHA_FOR_YDIFF)
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
        plicFuncs::calc_T_from_h(x_cellI_tmp, P_thermo, h_cellI_tmp, n, Pc, Tc, w, MW, tk, Tb, SG, H8, T_cellI_tmp, MAX_ITER_T, T_TOL, H_TOL);
        T1Cells[cellI] = T_cellI_tmp;
    }

    if(diff_debug2)
    {
        osDiff<< setw(6) << cellI << "  " << setw(12) << H1Cells[cellI] << "  " << setw(12) << rhoh1Cells[cellI] << "  " << setw(12) << h1Cells[cellI] << "  " << setw(12) << T1Cells[cellI] << endl;
    }
}

if(diff_debug2)
{
    plicFuncs::print_line(osDiff, 80);
}

Info<< "Phase-1 H = "
<< (H1).weightedAverage(mesh.Vsc()).value()
<< nl
<< "Min(H1) = " << min(H1).value()
<< "  Max(H1) = " << max(H1).value()
<< nl
<< "Min(T1) = " << min(T1).value()
<< "  Max(T1) = " << max(T1).value()
<< nl
<< "H1 undershoot = " << min((min(H1).value() - H1min), 0)
<< "  H1 overshoot = " << max((max(H1).value() - H1max), 0)
<< nl
<< "T1 undershoot = " << min((min(T1).value() - T1min), 0)
<< "  T1 overshoot = " << max((max(T1).value() - T1max), 0)
<< endl;
        
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
    plicFuncs::calc_h_from_T(x_cellI_tmp, P_thermo, T_cellI_tmp, n, Pc, Tc, w, MW, tk, Tb, SG, H8, h_cellI_tmp);

    h1Cells[cellI] = h_cellI_tmp;
    rhoh1Cells[cellI] = rho1Cells[cellI]*h1_cellI_tmp;
    H1Cells[cellI] = alpha1Cells[cellI]*rhoh1Cells[cellI];
}

Info<< "Phase-1 H after clipping = "
<< (H1).weightedAverage(mesh.Vsc()).value()
<< nl
<< "Min(H1) = " << min(H1).value()
<< "  Max(H1) = " << max(H1).value()
<< nl
<< "Min(T1) = " << min(T1).value()
<< "  Max(T1) = " << max(T1).value()        
<< endl;


//Phase-0
volScalarField surfInt_diffFlux_h0(fvc::surfaceIntegrate(diffFlux_h0));
scalarField& surfInt_diffFlux_h0Cells = surfInt_diffFlux_h0.internalField();

volScalarField surfInt_diffFlux_Jh0(fvc::surfaceIntegrate(diffFlux_Jh0));
scalarField& surfInt_diffFlux_Jh0Cells = surfInt_diffFlux_Jh0.internalField();

if(diff_debug2)
{
    plicFuncs::print_diffFluxFld(mesh, diffFlux_h0, alpha0, H0, h0, T0, "h0", "T0", osDiff);
    plicFuncs::print_diffFluxIntData(mesh, diffFlux_h0, surfInt_diffFlux_h0Cells, dt, alpha0, H0, h0, T0, Af_ph0_own, Af_ph0_nei, "ph0", "alpha0", "H0", "h0", "T0", osDiff);

    plicFuncs::print_line(osDiff, 80);
    osDiff<< setw(6) << "cellI" << "  " << setw(12) << "H0" << "  " << setw(12) << "rhoh0" << "  " << setw(12) << "h0" << "  " << setw(12) << "T0" << endl;
    plicFuncs::print_line(osDiff, 80);
}

forAll(H0Cells,cellI)
{                    
    H0Cells[cellI] += -(surfInt_diffFlux_h0Cells[cellI] + surfInt_diffFlux_Jh0Cells[cellI])*dt;

    if(alpha0Cells[cellI] < SMALLEST_ALPHA_FOR_YDIFF)
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
        plicFuncs::calc_T_from_h(x_cellI_tmp, P_thermo, h_cellI_tmp, n, Pc, Tc, w, MW, tk, Tb, SG, H8, T_cellI_tmp, MAX_ITER_T, T_TOL, H_TOL);
        T0Cells[cellI] = T_cellI_tmp;
    }

    if(diff_debug2)
    {
        osDiff<< setw(6) << cellI << "  " << setw(12) << H0Cells[cellI] << "  " << setw(12) << rhoh0Cells[cellI] << "  " << setw(12) << h0Cells[cellI] << "  " << setw(12) << T0Cells[cellI] << endl;
    }
}

if(diff_debug2)
{
    plicFuncs::print_line(osDiff, 80);
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
    plicFuncs::calc_h_from_T(x_cellI_tmp, P_thermo, T_cellI_tmp, n, Pc, Tc, w, MW, tk, Tb, SG, H8, h_cellI_tmp);

    h0Cells[cellI] = h_cellI_tmp;
    rhoh0Cells[cellI] = rho0Cells[cellI]*h0_cellI_tmp;
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

H0.correctBoundaryConditions();
rhoh0.correctBoundaryConditions();
h0.correctBoundaryConditions();
T0.correctBoundaryConditions();

//Done solve enthalpy diffusion equations in both phases
//------------------------------------------------------------------//
Info<< "Done HDiffEqns" << endl;
