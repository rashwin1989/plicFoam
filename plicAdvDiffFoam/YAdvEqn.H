
Y1t = 0.0*Y1[0];
Y0t = 0.0*Y0[0];

for(label i=0; i<(nSpecies - 1); i++)
{
    volScalarField& Y1i = Y1[i];
    volScalarField& Y0i = Y0[i];
    volScalarField& C1i = C1[i];
    volScalarField& C0i = C0[i];

    scalarField& Y1iCells = Y1i.internalField();
    scalarField& Y0iCells = Y0i.internalField();
    scalarField& C1iCells = C1i.internalField();
    scalarField& C0iCells = C0i.internalField();

    surfaceScalarField& advFlux_Y1i = advFlux_Y1[i];
    surfaceScalarField& advFlux_Y0i = advFlux_Y0[i];
    

    //Phase-1
    volScalarField surfInt_advFlux_Y1i(fvc::surfaceIntegrate(advFlux_Y1i));
    scalarField& surfInt_advFlux_Y1iCells = surfInt_advFlux_Y1i.internalField();

    forAll(C1iCells,cellI)
    {
        if(adv_debug2)
        {
            Info<< "Solving Y1[" << i << "] eqn in cell " << cellI << nl
                << "C1i = " << C1iCells[cellI] << "  Y1i = " << Y1iCells[cellI] << "  alpha1 = " << alpha1Cells[cellI] << nl
                << "div(advFlux) = " << surfInt_advFlux_Y1iCells[cellI] << "  div(advFlux)*dt = " << surfInt_advFlux_Y1iCells[cellI]*runTime.deltaTValue()
                << endl;
        }
                
        C1iCells[cellI] += -surfInt_advFlux_Y1iCells[cellI]*runTime.deltaTValue();

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
            Info<< "C1i = " << C1iCells[cellI] << nl
                << "Y1i = " << Y1iCells[cellI] << nl
                << endl;
        }
    }

    Info<< "Phase-1 Yi= "
        << (C1i).weightedAverage(mesh.Vsc()).value()
        << nl
        << "Min(C1i) = " << min(C1i).value()
        << "  Max(C1i) = " << max(C1i).value()
        << nl
        << "Min(Y1i) = " << min(Y1i).value()
        << "  Max(Y1i) = " << max(Y1i).value()
        << nl
        << "C1i undershoot = " << min(min(C1i).value(), 0)
        << "  C1i overshoot = " << max((max(C1i).value() - rho1_0.value()), 0)
        << nl
        << "Y1i undershoot = " << min(min(Y1i).value(), 0)
        << "  Y1i overshoot = " << max((max(Y1i).value() - 1), 0)
        << endl;
        
    forAll(C1iCells,cellI)
    {        
        if(Y1iCells[cellI] < 0)
        {
            Y1iCells[cellI] = 0;
            C1iCells[cellI] = rho1Cells[cellI]*alpha1Cells[cellI]*Y1iCells[cellI];
        }
        
        if(Y1iCells[cellI] > 1)
        {
            Y1iCells[cellI] = 1;
            C1iCells[cellI] = rho1Cells[cellI]*alpha1Cells[cellI]*Y1iCells[cellI];
        }        
    }

    Info<< "Phase-1 Yi after correction= "
        << (C1i).weightedAverage(mesh.Vsc()).value()
        << nl
        << "Min(C1i) = " << min(C1i).value()
        << "  Max(C1i) = " << max(C1i).value()
        << nl
        << "Min(Y1i) = " << min(Y1i).value()
        << "  Max(Y1i) = " << max(Y1i).value()        
        << endl;

    C1i.correctBoundaryConditions();
    Y1i.correctBoundaryConditions();

    //Phase-0
    volScalarField surfInt_advFlux_Y0i(fvc::surfaceIntegrate(advFlux_Y0i));
    scalarField& surfInt_advFlux_Y0iCells = surfInt_advFlux_Y0i.internalField();

    forAll(C0iCells,cellI)
    {
        if(adv_debug2)
        {
            Info<< "Solving Y0[" << i << "] eqn in cell " << cellI << nl
                << "Y0i = " << Y0iCells[cellI] << "  alpha0 = " << alpha0Cells[cellI] << nl
                << "div(advFlux) = " << surfInt_advFlux_Y0iCells[cellI] << "  div(advFlux)*dt = " << surfInt_advFlux_Y0iCells[cellI]*runTime.deltaTValue()
                << endl;
        }
                
        C0iCells[cellI] += -surfInt_advFlux_Y0iCells[cellI]*runTime.deltaTValue();

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
            Info<< "C0i = " << C0iCells[cellI] << nl
                << "Y0i = " << Y0iCells[cellI] << nl
                << endl;
        }
    }    
    
    Info<< "Phase-0 Yi= "
        << (C0i).weightedAverage(mesh.Vsc()).value()
        << "  Min(Y0i) = " << min(Y0i).value()
        << "  Max(Y0i) = " << max(Y0i).value()
        << nl
        << "Y0i undershoot = " << min(min(Y0i).value(), 0)
        << "  Y0i overshoot = " << max((max(Y0i).value() - 1), 0)
        << endl;
          
    forAll(C0iCells,cellI)
    {        
        if(Y0iCells[cellI] < 0)
        {
            Y0iCells[cellI] = 0;
            C0iCells[cellI] = rho0Cells[cellI]*alpha0Cells[cellI]*Y0iCells[cellI];

            if(adv_debug2)
            {
                Info<< "Y0[" << i << "] < 0 in cell " << cellI << endl;
            }
        }
        
        if(Y0iCells[cellI] > 1)
        {
            Y0iCells[cellI] = 1;
            C0iCells[cellI] = rho0Cells[cellI]*alpha0Cells[cellI]*Y0iCells[cellI];

            if(adv_debug2)
            {
                Info<< "Y0[" << i << "] > 1 in cell " << cellI << endl;
            }
        }        
    }
    
    Info<< "Phase-0 Yi after correction= "
        << (C0i).weightedAverage(mesh.Vsc()).value()
        << "  Min(Y0i) = " << min(Y0i).value()
        << "  Max(Y0i) = " << max(Y0i).value()        
        << endl;
        
    C0i.correctBoundaryConditions();
    Y0i.correctBoundaryConditions();

    Y1t += Y1i;
    Y0t += Y0i;
}

Y1[nSpecies - 1] = scalar(1.0) - Y1t;
//Y1[nSpecies - 1].max(0.0);
//Y1[nSpecies - 1].min(1.0);
C1[nSpecies - 1] = rho1*alpha1*Y1[nSpecies - 1];

Y0[nSpecies - 1] = scalar(1.0) - Y0t;
//Y0[nSpecies - 1].max(0.0);
//Y0[nSpecies - 1].min(1.0);
C0[nSpecies - 1] = rho0*alpha0*Y0[nSpecies - 1];

Info<< "Done YAdvEqns" << endl;
