
dt = deltaT;

Y1t = 0.0*Y1[0];
Y0t = 0.0*Y0[0];

for(i=0; i<(n-1); i++)
{
    if(adv_debug2)
    {
        plicFuncs::print_line(osAdv, 80);
        osAdv<< "Solution Time = " << runTime.timeName() << "  Y[" << i << "] Adv" << endl;
        plicFuncs::print_line(osAdv, 80);
        osAdv<< endl;    
    }

    volScalarField& Y1i = Y1[i];
    volScalarField& Y0i = Y0[i];
    volScalarField& c1i = c1[i];
    volScalarField& c0i = c0[i];
    volScalarField& C1i = C1[i];
    volScalarField& C0i = C0[i];
    volScalarField& C1Oldi = C1_old[i];
    volScalarField& C0Oldi = C0_old[i];

    scalarField& Y1iCells = Y1i.internalField();
    scalarField& Y0iCells = Y0i.internalField();
    scalarField& c1iCells = c1i.internalField();
    scalarField& c0iCells = c0i.internalField();
    scalarField& C1iCells = C1i.internalField();
    scalarField& C0iCells = C0i.internalField();
    scalarField& C1OldiCells = C1Oldi.internalField();
    scalarField& C0OldiCells = C0Oldi.internalField();

    word C1iName = "C1["+name(i)+"]";
    word c1iName = "c1["+name(i)+"]";
    word Y1iName = "Y1["+name(i)+"]";
    word C0iName = "C0["+name(i)+"]";
    word c0iName = "c0["+name(i)+"]";
    word Y0iName = "Y0["+name(i)+"]";

    surfaceScalarField& advFlux_Y1i = advFlux_Y1[i];
    surfaceScalarField& advFlux_Y0i = advFlux_Y0[i];

    //Phase-1
    volScalarField surfInt_advFlux_Y1i(fvc::surfaceIntegrate(advFlux_Y1i));
    scalarField& surfInt_advFlux_Y1iCells = surfInt_advFlux_Y1i.internalField();

    if(adv_debug2)
    {
        plicFuncs::print_advFluxFld(mesh, advFlux_Y1i, alpha1, C1i, c1i, Y1i, c1iName, Y1iName, osAdv);
        plicFuncs::print_advFluxIntData(mesh, advFlux_Y1i, surfInt_advFlux_Y1iCells, dt, alpha1, C1i, c1i, Y1i, Af_ph1_own, Af_ph1_nei, "ph1", "alpha1", C1iName, c1iName, Y1iName, osAdv);

        plicFuncs::print_line(osAdv, 80);
        osAdv<< setw(6) << "cellI" << "  " << setw(12) << C1iName << "  " << setw(12) << Y1iName << endl;
        plicFuncs::print_line(osAdv, 80);
    }

    forAll(C1iCells,cellI)
    {                        
        C1iCells[cellI] = C1OldiCells[cellI] - surfInt_advFlux_Y1iCells[cellI]*dt;

        if(alpha1Cells[cellI] < SMALLEST_ALPHA_FOR_YADV)
        {
            Y1iCells[cellI] = 0;
        }
        else
        {
            Y1iCells[cellI] = C1iCells[cellI]/alpha1Cells[cellI]/rho1Cells[cellI];
        }

        if(adv_debug2)
        {
            osAdv<< setw(6) << cellI << "  " << setw(12) << C1iCells[cellI] << "  " << setw(12) << Y1iCells[cellI] << endl;
        }
    }

    if(adv_debug2)
    {
        plicFuncs::print_line(osAdv, 80);
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

    if(adv_debug2)
    {
        plicFuncs::print_line(osAdv, 80);
        osAdv<< "Redistributing " << C1iName << " field" << endl;
        plicFuncs::print_line(osAdv, 80);
        osAdv<< endl;
    }
    
    if(redist_C)
    {
        plicFuncs::redistribute_Ci_field(mesh, C1iCells, Y1iCells, alpha1Cells, rho1Cells, cell_stencil.stencil(), Y1min[i], Y1max[i], SMALLEST_ALPHA_FOR_YADV, Y_BOUND_ITERS_MAX, adv_debug2, osAdv);

        Info<< "Phase-1 Ci after redistribution = "
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
    }
   
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

    C1i.correctBoundaryConditions();
    c1i.correctBoundaryConditions();
    Y1i.correctBoundaryConditions();

    //Phase-0
    volScalarField surfInt_advFlux_Y0i(fvc::surfaceIntegrate(advFlux_Y0i));
    scalarField& surfInt_advFlux_Y0iCells = surfInt_advFlux_Y0i.internalField();

    if(adv_debug2)
    {
        plicFuncs::print_advFluxFld(mesh, advFlux_Y0i, alpha0, C0i, c0i, Y0i, c0iName, Y0iName, osAdv);
        plicFuncs::print_advFluxIntData(mesh, advFlux_Y0i, surfInt_advFlux_Y0iCells, dt, alpha0, C0i, c0i, Y0i, Af_ph0_own, Af_ph0_nei, "ph0", "alpha0", C0iName, c0iName, Y0iName, osAdv);

        plicFuncs::print_line(osAdv, 80);
        osAdv<< setw(6) << "cellI" << "  " << setw(12) << C0iName << "  " << setw(12) << Y0iName << endl;
        plicFuncs::print_line(osAdv, 80);
    }

    forAll(C0iCells,cellI)
    {                
        C0iCells[cellI] = C0OldiCells[cellI] - surfInt_advFlux_Y0iCells[cellI]*dt;

        if(alpha0Cells[cellI] < SMALLEST_ALPHA_FOR_YADV)
        {
            Y0iCells[cellI] = 0;
        }
        else
        {
            Y0iCells[cellI] = C0iCells[cellI]/alpha0Cells[cellI]/rho0Cells[cellI];
        }
        
        if(adv_debug2)
        {
            osAdv<< setw(6) << cellI << "  " << setw(12) << C0iCells[cellI] << "  " << setw(12) << Y0iCells[cellI] << endl;
        }
    }

    if(adv_debug2)
    {
        plicFuncs::print_line(osAdv, 80);
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

    if(adv_debug2)
    {
        plicFuncs::print_line(osAdv, 80);
        osAdv<< "Redistributing " << C0iName << " field" << endl;
        plicFuncs::print_line(osAdv, 80);
        osAdv<< endl;
    }

    if(redist_C)
    {
        plicFuncs::redistribute_Ci_field(mesh, C0iCells, Y0iCells, alpha0Cells, rho0Cells, cell_stencil.stencil(), Y0min[i], Y0max[i], SMALLEST_ALPHA_FOR_YADV, Y_BOUND_ITERS_MAX, adv_debug2, osAdv);

        Info<< "Phase-0 Ci after redistribution = "
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
    }     
     
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
        
    C0i.correctBoundaryConditions();
    c0i.correctBoundaryConditions();
    Y0i.correctBoundaryConditions();

    Y1t += Y1i;
    Y0t += Y0i;
}

Y1[nSpecies - 1] = scalar(1.0) - Y1t;
Y1[nSpecies - 1].max(0.0);
Y1[nSpecies - 1].min(1.0);
C1[nSpecies - 1] = rho1*alpha1*Y1[nSpecies - 1];
c1[nSpecies - 1] = rho1*Y1[nSpecies - 1];

Y0[nSpecies - 1] = scalar(1.0) - Y0t;
Y0[nSpecies - 1].max(0.0);
Y0[nSpecies - 1].min(1.0);
C0[nSpecies - 1] = rho0*alpha0*Y0[nSpecies - 1];
c0[nSpecies - 1] = rho0*Y0[nSpecies - 1];

forAll(alpha1Cells, cellI)
{
    for(i=0; i<n; i++)
    {
        y_cellI_tmp[i] = Y1[i].internalField()[cellI];
    }
    plicFuncs::y2x(n, MW, y_cellI_tmp, x_cellI_tmp);
    for(i=0; i<n; i++)
    {
        x1[i].internalField()[cellI] = x_cellI_tmp;
    }

    for(i=0; i<n; i++)
    {
        y_cellI_tmp[i] = Y0[i].internalField()[cellI];
    }
    plicFuncs::y2x(n, MW, y_cellI_tmp, x_cellI_tmp);
    for(i=0; i<n; i++)
    {
        x0[i].internalField()[cellI] = x_cellI_tmp;
    }
}

Info<< "Done YAdvEqns" << endl;
