if(adv_debug)
{
    plicFuncs::print_line(osAdv, 100);
    osAdv<< "Solution Time = " << runTime.timeName() << "  Y Adv" << endl;
    plicFuncs::print_line(osAdv, 100);
    osAdv<< endl;
}

dt = deltaT;

Y1t = 0.0*Y1[0];
Y0t = 0.0*Y0[0];

forAll(alpha1Cells, cellI)
{
    unbounded_Y_in_cell[cellI] = 0;
}

for(i=0; i<(n-1); i++)
{
    volScalarField& Y1i = Y1[i];
    volScalarField& Y0i = Y0[i];
    volScalarField& c1i = c1[i];
    volScalarField& c0i = c0[i];
    volScalarField& C1i = C1[i];
    volScalarField& C0i = C0[i];
    //volScalarField& C1Oldi = C1_old[i];
    //volScalarField& C0Oldi = C0_old[i];

    scalarField& Y1iCells = Y1i.internalField();
    scalarField& Y0iCells = Y0i.internalField();
    scalarField& c1iCells = c1i.internalField();
    scalarField& c0iCells = c0i.internalField();
    scalarField& C1iCells = C1i.internalField();
    scalarField& C0iCells = C0i.internalField();
    //scalarField& C1OldiCells = C1Oldi.internalField();
    //scalarField& C0OldiCells = C0Oldi.internalField();

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

    plicFuncs::print_line(osAdv, 100);
    osAdv<< "Advecting " << C1iName << " field" << endl;
    plicFuncs::print_line(osAdv, 100);
    osAdv<< endl;

    if(adv_debug2)
    {        
        plicFuncs::print_advFluxFld(mesh, advFlux_Y1i, alpha1, C1i, c1i, Y1i, c1iName, Y1iName, osAdv);
        plicFuncs::print_advFluxIntData(mesh, advFlux_Y1i, surfInt_advFlux_Y1iCells, dt, alpha1, C1i, c1i, Y1i, Af_ph1_own, Af_ph1_nei, "ph1", "alpha1", C1iName, c1iName, Y1iName, osAdv);
    }
    if(adv_debug)
    {
        plicFuncs::print_line(osAdv, 100);
        osAdv<< setw(6) << "cellI" << "  " << setw(12) << "alpha1" << "  " << setw(12) << "rho1" << "  " << setw(12) << Y1iName << "  " << setw(12) << c1iName << "  " << setw(12) << C1iName << endl;
        plicFuncs::print_line(osAdv, 100);
    }

    forAll(C1iCells,cellI)
    {                        
        C1iCells[cellI] = C1iCells[cellI] - surfInt_advFlux_Y1iCells[cellI]*dt;

        if(alpha1Cells[cellI] < SMALLEST_ALPHA_FOR_YADV)
        {
            Y1iCells[cellI] = Y1_TMP[i];
        }
        else
        {
            Y1iCells[cellI] = C1iCells[cellI]/alpha1Cells[cellI]/rho1Cells[cellI];
        }

        if(adv_debug)
        {
            osAdv<< setw(6) << cellI << "  " << setw(12) << alpha1Cells[cellI] << "  " << setw(12) << rho1Cells[cellI] << "  " << setw(12) << Y1iCells[cellI] << "  " << setw(12) << c1iCells[cellI] << "  " << setw(12) << C1iCells[cellI] << endl;
        }
    }

    if(adv_debug)
    {
        plicFuncs::print_line(osAdv, 100);
        osAdv<< endl;
    }

    C1_avg[i] = (C1i).weightedAverage(mesh.V()).value();

    Info<< "C1[" << i << "]:"
        << "  avg = " << C1_avg[i]
        << "  min = " << min(C1i).value()
        << "  max = " << max(C1i).value()
        << nl
        << "Y1[" << i << "]:"
        << "  min = " << min(Y1i).value()
        << "  max = " << max(Y1i).value()        
        << endl;
    
    if(redist_C)
    {
        if(adv_debug)
        {
            plicFuncs::print_line(osAdv, 100);
            osAdv<< "Redistributing " << C1iName << " field" << endl;
            plicFuncs::print_line(osAdv, 100);
            osAdv<< endl;
        }

        plicFuncs::redistribute_Ci_field(mesh, C1iCells, Y1iCells, alpha1Cells, rho1Cells, cell_stencil.stencil(), Y1min[i], Y1max[i], SMALLEST_ALPHA_FOR_YADV, Y_BOUND_ITERS_MAX, adv_debug, osAdv);

        C1_avg[i] = (C1i).weightedAverage(mesh.V()).value();

        Info<< "Redistributed C1[" << i << "]:"
            << "  avg = " << C1_avg[i]
            << "  min = " << min(C1i).value()
            << "  max = " << max(C1i).value()
            << nl
            << "Y1[" << i << "]:"
            << "  min = " << min(Y1i).value()
            << "  max = " << max(Y1i).value()        
            << endl;
    }
   
    forAll(C1iCells,cellI)
    {        
        if(Y1iCells[cellI] < Y1min[i])
        {
            Y1iCells[cellI] = Y1min[i];
            unbounded_Y_in_cell[cellI] = 1;
        }
        
        if(Y1iCells[cellI] > Y1max[i])
        {
            Y1iCells[cellI] = Y1max[i];
            unbounded_Y_in_cell[cellI] = 1;
        }   
        
        c1iCells[cellI] = rho1Cells[cellI]*Y1iCells[cellI];
        C1iCells[cellI] = alpha1Cells[cellI]*c1iCells[cellI];
    }

    C1_clip_avg[i] = (C1i).weightedAverage(mesh.V()).value();
    C1_clip_err[i] = mag(C1_clip_avg[i] - C1_avg[i]);
    C1_clip_cum_err[i] += C1_clip_err[i];

    Info<< "Clipped C1[" << i << "]:"        
        << "  avg = " << C1_clip_avg[i]
        << "  min = " << min(C1i).value()
        << "  max = " << max(C1i).value()
        << "  err = " << C1_clip_err[i]
        << "  cum_err = " << C1_clip_cum_err[i]
        << endl;

    Y1i.correctBoundaryConditions();
    plicFuncs::correct_boundaryField_C(Y1i, rho1, alpha1, c1i, C1i);

    //Phase-0
    volScalarField surfInt_advFlux_Y0i(fvc::surfaceIntegrate(advFlux_Y0i));
    scalarField& surfInt_advFlux_Y0iCells = surfInt_advFlux_Y0i.internalField();

    plicFuncs::print_line(osAdv, 100);
    osAdv<< "Advecting " << C0iName << " field" << endl;
    plicFuncs::print_line(osAdv, 100);
    osAdv<< endl;

    if(adv_debug2)
    {
        plicFuncs::print_advFluxFld(mesh, advFlux_Y0i, alpha0, C0i, c0i, Y0i, c0iName, Y0iName, osAdv);
        plicFuncs::print_advFluxIntData(mesh, advFlux_Y0i, surfInt_advFlux_Y0iCells, dt, alpha0, C0i, c0i, Y0i, Af_ph0_own, Af_ph0_nei, "ph0", "alpha0", C0iName, c0iName, Y0iName, osAdv);
    }
    if(adv_debug)
    {
        plicFuncs::print_line(osAdv, 100);
        osAdv<< setw(6) << "cellI" << "  " << setw(12) << "alpha0" << "  " << setw(12) << "rho0" << "  " << setw(12) << Y0iName << "  " << setw(12) << c0iName << "  " << setw(12) << C0iName << endl;
        plicFuncs::print_line(osAdv, 100);
    }

    forAll(C0iCells,cellI)
    {                
        C0iCells[cellI] = C0iCells[cellI] - surfInt_advFlux_Y0iCells[cellI]*dt;

        if(alpha0Cells[cellI] < SMALLEST_ALPHA_FOR_YADV)
        {
            Y0iCells[cellI] = Y0_TMP[i];
        }
        else
        {
            Y0iCells[cellI] = C0iCells[cellI]/alpha0Cells[cellI]/rho0Cells[cellI];
        }
        
        if(adv_debug)
        {
            osAdv<< setw(6) << cellI << "  " << setw(12) << alpha0Cells[cellI] << "  " << setw(12) << rho0Cells[cellI] << "  " << setw(12) << Y0iCells[cellI] << "  " << setw(12) << c0iCells[cellI] << "  " << setw(12) << C0iCells[cellI] << endl;
        }
    }

    if(adv_debug)
    {
        plicFuncs::print_line(osAdv, 100);
        osAdv<< endl;
    }

    C0_avg[i] = (C0i).weightedAverage(mesh.V()).value();

    Info<< "C0[" << i << "]:"
        << "  avg = " << C0_avg[i]
        << "  min = " << min(C0i).value()
        << "  max = " << max(C0i).value()
        << nl
        << "Y0[" << i << "]:"
        << "  min = " << min(Y0i).value()
        << "  max = " << max(Y0i).value()        
        << endl;

    if(redist_C)
    {
        if(adv_debug)
        {
            plicFuncs::print_line(osAdv, 100);
            osAdv<< "Redistributing " << C0iName << " field" << endl;
            plicFuncs::print_line(osAdv, 100);
            osAdv<< endl;
        }

        plicFuncs::redistribute_Ci_field(mesh, C0iCells, Y0iCells, alpha0Cells, rho0Cells, cell_stencil.stencil(), Y0min[i], Y0max[i], SMALLEST_ALPHA_FOR_YADV, Y_BOUND_ITERS_MAX, adv_debug, osAdv);

        C0_avg[i] = (C0i).weightedAverage(mesh.V()).value();

        Info<< "Redistributed C0[" << i << "]:"
            << "  avg = " << C0_avg[i]
            << "  min = " << min(C0i).value()
            << "  max = " << max(C0i).value()
            << nl
            << "Y0[" << i << "]:"
            << "  min = " << min(Y0i).value()
            << "  max = " << max(Y0i).value()        
            << endl;
    }
     
    forAll(C0iCells,cellI)
    {        
        if(Y0iCells[cellI] < Y0min[i])
        {
            Y0iCells[cellI] = Y0min[i];
            unbounded_Y_in_cell[cellI] = 1;
        }
        
        if(Y0iCells[cellI] > Y0max[i])
        {
            Y0iCells[cellI] = Y0max[i];
            unbounded_Y_in_cell[cellI] = 1;
        }   

        c0iCells[cellI] = rho0Cells[cellI]*Y0iCells[cellI];
        C0iCells[cellI] = alpha0Cells[cellI]*c0iCells[cellI];
    }
    
    C0_clip_avg[i] = (C0i).weightedAverage(mesh.V()).value();
    C0_clip_err[i] = mag(C0_clip_avg[i] - C0_avg[i]);
    C0_clip_cum_err[i] += C0_clip_err[i];

    Info<< "Clipped C0[" << i << "]:"        
        << "  avg = " << C0_clip_avg[i]
        << "  min = " << min(C0i).value()
        << "  max = " << max(C0i).value()
        << "  err = " << C0_clip_err[i]
        << "  cum_err = " << C0_clip_cum_err[i]
        << endl;

    Y0i.correctBoundaryConditions();
    plicFuncs::correct_boundaryField_C(Y0i, rho0, alpha0, c0i, C0i);

    Y1t += Y1i;
    Y0t += Y0i;
}

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

Info<< "Done YAdvEqns" << endl;
