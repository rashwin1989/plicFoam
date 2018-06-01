
if(diff_debug)
{
    plicFuncs::print_line(osDiff, 80);
    osDiff<< "Solution Time = " << runTime.timeName() << "  Y Diff" << endl;
    plicFuncs::print_line(osDiff, 80);
    osDiff<< endl;
}

//set fractional time-step for diffusion
dt = 0.5*deltaT;

Y1t = 0.0*Y1[0];
Y0t = 0.0*Y0[0];

//------------------------------------------------------------------//
//Calculate species diffusive fluxes for both phases
for(i=0; i<n; i++)
{
    const volScalarField& Y1i = Y1[i];
    const volScalarField& Y0i = Y0[i];
    const surfaceScalarField& D1fi = D1f[i];
    const surfaceScalarField& D0fi = D0f[i];
    const surfaceScalarField& gradf_Y1i = gradf_Y1[i];
    const surfaceScalarField& gradf_Y0i = gradf_Y0[i];

    surfaceScalarField& diffFlux_Y1i = diffFlux_Y1[i];
    surfaceScalarField& diffFlux_Y0i = diffFlux_Y0[i];

    plicFuncs::calc_2ph_diffFluxes_Yi_Fick
    (
        mesh,
        rho1f,
        rho0f,
        D1fi,
        D0fi,
        gradf_Y1i,
        gradf_Y0i,
        Af_ph1_own,
        Af_ph0_own,
        Af_ph1_nei,
        Af_ph0_nei,
        face_phaseState_diff,
        rho1,
        alpha1,
        Y1i,
        rho0,
        alpha0,
        Y0i,
        diffFlux_Y1i,
        diffFlux_Y0i,
        Y1min[i],
        Y1max[i],
        Y0min[i],
        Y0max[i],
        i,
        dt,
        diff_debug2,
        osDiff
    );
}
//Done calculate species diffusive fluxes for both phases
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Solve species diffusion equations in both phases
for(i=0; i<(n-1); i++)
{
    word YiName = "Y["+name(i)+"]";
    if(diff_debug)
    {
        plicFuncs::print_line(osDiff, 80);
        osDiff<< setw(40) << YiName << endl;
        plicFuncs::print_line(osDiff, 80);
        osDiff<< endl;
    }

    volScalarField& Y1i = Y1[i];
    volScalarField& Y0i = Y0[i];
    volScalarField& c1i = c1[i];
    volScalarField& c0i = c0[i];
    volScalarField& C1i = C1[i];
    volScalarField& C0i = C0[i];    

    scalarField& Y1iCells = Y1i.internalField();
    scalarField& Y0iCells = Y0i.internalField();
    scalarField& c1iCells = c1i.internalField();
    scalarField& c0iCells = c0i.internalField();
    scalarField& C1iCells = C1i.internalField();
    scalarField& C0iCells = C0i.internalField();    

    word C1iName = "C1["+name(i)+"]";
    word c1iName = "c1["+name(i)+"]";
    word Y1iName = "Y1["+name(i)+"]";
    word C0iName = "C0["+name(i)+"]";
    word c0iName = "c0["+name(i)+"]";
    word Y0iName = "Y0["+name(i)+"]";

    const surfaceScalarField& diffFlux_Y1i = diffFlux_Y1[i];
    const surfaceScalarField& diffFlux_Y0i = diffFlux_Y0[i];

    //--------------------------------------------------------------//
    //ph1
    volScalarField surfInt_diffFlux_Y1i(fvc::surfaceIntegrate(diffFlux_Y1i));    
    scalarField& surfInt_diffFlux_Y1iCells = surfInt_diffFlux_Y1i.internalField();

    if(diff_debug)
    {
        plicFuncs::print_diffFluxFld(mesh, diffFlux_Y1i, alpha1, C1i, c1i, Y1i, c1iName, Y1iName, osDiff);
        plicFuncs::print_diffFluxIntData(mesh, diffFlux_Y1i, surfInt_diffFlux_Y1iCells, dt, alpha1, C1i, c1i, Y1i, Af_ph1_own, Af_ph1_nei, "ph1", "alpha1", C1iName, c1iName, Y1iName, osDiff);

        plicFuncs::print_line(osDiff, 80);
        osDiff<< setw(6) << "cellI" << "  " << setw(12) << C1iName << "  " << setw(12) << Y1iName << endl;
        plicFuncs::print_line(osDiff, 80);
    }

    forAll(C1iCells,cellI)
    {                
        C1iCells[cellI] += -surfInt_diffFlux_Y1iCells[cellI]*dt;

        if(alpha1Cells[cellI] < SMALLEST_ALPHA_FOR_YDIFF)
        {
            Y1iCells[cellI] = 0;
        }
        else
        {
            Y1iCells[cellI] = C1iCells[cellI]/alpha1Cells[cellI]/rho1Cells[cellI];
        }

        if(diff_debug)
        {
            osDiff<< setw(6) << cellI << "  " << setw(12) << C1iCells[cellI] << "  " << setw(12) << Y1iCells[cellI] << endl;
        }
    }

    if(diff_debug)
    {
        plicFuncs::print_line(osDiff, 80);
    }

    Info<< "Phase-1 Ci = "
        << (C1i).weightedAverage(mesh.Vsc()).value()
        << nl
        << "Min(C1i) = " << min(C1i).value()
        << "  Max(C1i) = " << max(C1i).value()
        << nl
        << "Min(Y1i) = " << min(Y1i).value()
        << "  Max(Y1i) = " << max(Y1i).value()
        << nl
        << "C1i undershoot = " << min((min(C1i).value() - C1min[i]), 0)
        << "  C1i overshoot = " << max((max(C1i).value() - C1max[i]), 0)
        << nl
        << "Y1i undershoot = " << min((min(Y1i).value() - Y1min[i]), 0)
        << "  Y1i overshoot = " << max((max(Y1i).value() - Y1max[i]), 0)
        << endl;    

    forAll(C1iCells,cellI)
    {        
        if(Y1iCells[cellI] < Y1min[i])
        {
            Y1iCells[cellI] = Y1min[i];            
        }
        
        if(Y1iCells[cellI] > Y1max[i])
        {
            Y1iCells[cellI] = Y1max[i];            
        }

        C1iCells[cellI] = rho1Cells[cellI]*alpha1Cells[cellI]*Y1iCells[cellI];
        c1iCells[cellI] = rho1Cells[cellI]*Y1iCells[cellI];
    }

    Info<< "Phase-1 Ci after clipping = "
        << (C1i).weightedAverage(mesh.Vsc()).value()
        << nl
        << "Min(C1i) = " << min(C1i).value()
        << "  Max(C1i) = " << max(C1i).value()
        << nl
        << "Min(Y1i) = " << min(Y1i).value()
        << "  Max(Y1i) = " << max(Y1i).value()        
        << endl;

    Y1i.correctBoundaryConditions();
    plicFuncs::correct_boundaryField_C(Y1i, rho1, alpha1, c1i, C1i);
    //Done ph1
    //--------------------------------------------------------------//

    //--------------------------------------------------------------//
    //ph0
    volScalarField surfInt_diffFlux_Y0i(fvc::surfaceIntegrate(diffFlux_Y0i));    
    scalarField& surfInt_diffFlux_Y0iCells = surfInt_diffFlux_Y0i.internalField();

    if(diff_debug)
    {
        plicFuncs::print_diffFluxFld(mesh, diffFlux_Y0i, alpha0, C0i, c0i, Y0i, c0iName, Y0iName, osDiff);
        plicFuncs::print_diffFluxIntData(mesh, diffFlux_Y0i, surfInt_diffFlux_Y0iCells, dt, alpha0, C0i, c0i, Y0i, Af_ph0_own, Af_ph0_nei, "ph0", "alpha0", C0iName, c0iName, Y0iName, osDiff);

        plicFuncs::print_line(osDiff, 80);
        osDiff<< setw(6) << "cellI" << "  " << setw(12) << C0iName << "  " << setw(12) << Y0iName << endl;
        plicFuncs::print_line(osDiff, 80);
    }

    forAll(C0iCells,cellI)
    {                
        C0iCells[cellI] += -surfInt_diffFlux_Y0iCells[cellI]*dt;

        if(alpha0Cells[cellI] < SMALLEST_ALPHA_FOR_YDIFF)
        {
            Y0iCells[cellI] = 0;
        }
        else
        {
            Y0iCells[cellI] = C0iCells[cellI]/alpha0Cells[cellI]/rho0Cells[cellI];
        }

        if(diff_debug)
        {
            osDiff<< setw(6) << cellI << "  " << setw(12) << C0iCells[cellI] << "  " << setw(12) << Y0iCells[cellI] << endl;
        }
    }

    if(diff_debug)
    {
        plicFuncs::print_line(osDiff, 80);
    }

    Info<< "Phase-0 Ci = "
        << (C0i).weightedAverage(mesh.Vsc()).value()
        << nl
        << "Min(C0i) = " << min(C0i).value()
        << "  Max(C0i) = " << max(C0i).value()
        << nl
        << "Min(Y0i) = " << min(Y0i).value()
        << "  Max(Y0i) = " << max(Y0i).value()
        << nl
        << "C0i undershoot = " << min((min(C0i).value() - C0min[i]), 0)
        << "  C0i overshoot = " << max((max(C0i).value() - C0max[i]), 0)
        << nl
        << "Y0i undershoot = " << min((min(Y0i).value() - Y0min[i]), 0)
        << "  Y0i overshoot = " << max((max(Y0i).value() - Y0max[i]), 0)
        << endl;    

    forAll(C0iCells,cellI)
    {        
        if(Y0iCells[cellI] < Y0min[i])
        {
            Y0iCells[cellI] = Y0min[i];            
        }
        
        if(Y0iCells[cellI] > Y0max[i])
        {
            Y0iCells[cellI] = Y0max[i];            
        }

        C0iCells[cellI] = rho0Cells[cellI]*alpha0Cells[cellI]*Y0iCells[cellI];
        c0iCells[cellI] = rho0Cells[cellI]*Y0iCells[cellI];
    }

    Info<< "Phase-0 Ci after clipping = "
        << (C0i).weightedAverage(mesh.Vsc()).value()
        << nl
        << "Min(C0i) = " << min(C0i).value()
        << "  Max(C0i) = " << max(C0i).value()
        << nl
        << "Min(Y0i) = " << min(Y0i).value()
        << "  Max(Y0i) = " << max(Y0i).value()        
        << endl;

    Y0i.correctBoundaryConditions();
    plicFuncs::correct_boundaryField_C(Y0i, rho0, alpha0, c0i, C0i);
    //Done ph0
    //--------------------------------------------------------------//

    Y1t += Y1i;
    Y0t += Y0i;
}
//Done solve species diffusion equations in both phases
//------------------------------------------------------------------//

Y1[n-1] = scalar(1.0) - Y1t;
Y1[n-1].max(0.0);
Y1[n-1].min(1.0);
C1[n-1] = rho1*alpha1*Y1[n-1];
c1[n-1] = rho1*Y1[n-1];
Y1[n-1].correctBoundaryConditions();
plicFuncs::correct_boundaryField_C(Y1[n-1], rho1, alpha1, c1[n-1], C1[n-1]);

Y0[n-1] = scalar(1.0) - Y0t;
Y0[n-1].max(0.0);
Y0[n-1].min(1.0);
C0[n-1] = rho0*alpha0*Y0[n-1];
c0[n-1] = rho0*Y0[n-1];
Y0[n-1].correctBoundaryConditions();
plicFuncs::correct_boundaryField_C(Y0[n-1], rho0, alpha0, c0[n-1], C0[n-1]);

plicFuncs::correct_x_from_Y(n, MW, Y1, x1);
plicFuncs::correct_x_from_Y(n, MW, Y0, x0);

Info<< "Done YDiffEqns" << endl;
