
nSpecies = Y1.size();

volScalarField Y1t = 0.0*Y1[0];
volScalarField Y0t = 0.0*Y0[0];

scalarField& alpha1Cells = alpha1.internalField();
scalarField& alpha0Cells = alpha0.internalField();

const labelList& face_phaseState = interface.face_phaseState();

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
    Foam::plic::calc_2ph_advFlux(Y1i, Y0i, advFlux_Y1i, advFlux_Y0i);

    //Phase-1
    volScalarField surfInt_advFlux_Y1i(fvc::surfaceIntegrate(advFlux_Y1i));
    scalarField& surfInt_advFlux_Y1iCells = surfInt_advFlux_Y1i.internalField();

    forAll(C1iCells,cellI)
    {
        if(adv_debug2)
        {
            Info<< "Solving Y1[" << i << "] eqn in cell " << cellI << nl
                << "Y1i = " << Y1iCells[cellI] << "  alpha1 = " << alpha1Cells[cellI] << nl
                << "div(advFlux) = " << surfInt_advFlux_Y1iCells[cellI] << "  div(advFlux)*dt = " << surfInt_advFlux_Y1iCells[cellI]*runTime.deltaTValue()
                << endl;
        }
                
        C1iCells[cellI] += -surfInt_advFlux_Y1iCells[cellI]*runTime.deltaTValue();

        if(alpha1Cells[cellI] < 1E-15)
        {
            Y1iCells[cellI] = 0;
        }
        else
        {
            Y1iCells[cellI] = C1iCells[cellI]/alpha1Cells[cellI]/rho1;
        }

        if(adv_debug2)
        {
            Info<< "C1i = " << C1iCells[cellI] << nl
                << "Y1i = " << Y1iCells[cellI] << nl
                << endl;
        }
    }

    C1i.correctBoundaryConditions();
    Y1i.correctBoundaryConditions();

    Info<< "Phase-1 Yi= "
        << (C1i).weightedAverage(mesh.Vsc()).value()
        << "  Min(Y1i) = " << min(Y1i).value()
        << "  Max(Y1i) = " << max(Y1i).value()
        << nl
        << "Y1i undershoot = " << min(min(Y1i).value(), 0)
        << "  Y1i overshoot = " << max((max(Y1i).value() - 1), 0)
        << endl;
        
    forAll(C1iCells,cellI)
    {        
        if(Y1iCells[cellI] < 0)
        {
            Y1iCells[cellI] = 0;
            C1iCells[cellI] = rho1*alpha1Cells[cellI]*Y1iCells[cellI];
        }
        
        if(Y1iCells[cellI] > 1)
        {
            Y1iCells[cellI] = 1;
            C1iCells[cellI] = rho1*alpha1Cells[cellI]*Y1iCells[cellI];
        }        
    }

    Info<< "Phase-1 Yi after correction= "
        << (C1i).weightedAverage(mesh.Vsc()).value()
        << "  Min(Y1i) = " << min(Y1i).value()
        << "  Max(Y1i) = " << max(Y1i).value()        
        << endl;

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

        if(alpha0Cells[cellI] < 1E-15)
        {
            Y0iCells[cellI] = 0;
        }
        else
        {
            Y0iCells[cellI] = C0iCells[cellI]/alpha0Cells[cellI]/rho0;
        }

        if(adv_debug2)
        {
            Info<< "C0i = " << C0iCells[cellI] << nl
                << "Y0i = " << Y0iCells[cellI] << nl
                << endl;
        }
    }

    C0i.correctBoundaryConditions();
    Y0i.correctBoundaryConditions();

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
            C0iCells[cellI] = rho0*alpha0Cells[cellI]*Y0iCells[cellI];
        }
        
        if(Y0iCells[cellI] > 1)
        {
            Y0iCells[cellI] = 1;
            C0iCells[cellI] = rho0*alpha0Cells[cellI]*Y0iCells[cellI];
        }        
    }

    Info<< "Phase-0 Yi after correction= "
        << (C0i).weightedAverage(mesh.Vsc()).value()
        << "  Min(Y0i) = " << min(Y0i).value()
        << "  Max(Y0i) = " << max(Y0i).value()        
        << endl;

    Y1t += Y1i;
    Y0t += Y0i;
}

Y1[nSpecies - 1] = scalar(1.0) - Y1t;
Y1[nSpecies - 1].max(0.0);
Y1[nSpecies - 1].min(1.0);
C1[nSpecies - 1] = rho1*alpha1*Y1[nSpecies - 1];

Y0[nSpecies - 1] = scalar(1.0) - Y0t;
Y0[nSpecies - 1].max(0.0);
Y0[nSpecies - 1].min(1.0);
C0[nSpecies - 1] = rho0*alpha0*Y0[nSpecies - 1];

Info<< "Done YAdvEqns" << endl;
