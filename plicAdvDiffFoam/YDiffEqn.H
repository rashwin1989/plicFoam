
Y1t = 0.0*Y1[0];
Y0t = 0.0*Y0[0];

//------------------------------------------------------------------//
//Make flat fields for alpha1 and phase centroids required for gradient calculation
plicFuncs::makeFlatFld(alpha1, diffMap, alpha1_flatFld_diff);
plicFuncs::makeFlatFld(interface.C_ph1(), diffMap, C_ph1_flatFld_diff);
plicFuncs::makeFlatFld(interface.C_ph0(), diffMap, C_ph0_flatFld_diff);
//Done make flat fields for alpha1 and phase centroids required for gradient calculation
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Obtain references to phase centroids, phase face centres, phase face areas, face phase states
const scalarField& Af_ph1 = interface.Af_ph1();
const scalarField& Af_ph0 = interface.Af_ph0();
const vectorField& Cf_ph1 = interface.Cf_ph1();
const vectorField& Cf_ph0 = interface.Cf_ph0();
const vectorField& C_ph1 = interface.C_ph1();
const vectorField& C_ph0 = interface.C_ph0();
const labelList& face_phaseState = interface.face_phaseState();
//Done obtain references to phase centroids, phase face centres, phase face areas, face phase states
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Linear interpolate density and diffusivity to faces in both phases
plicFuncs::calc_2ph_linearInterpolation_weights(mesh, C_ph1, Cf_ph1, w_ph1);
plicFuncs::calc_2ph_linearInterpolation_weights(mesh, C_ph0, Cf_ph0, w_ph0);

plicFuncs::linearInterpolate_2ph(rho1, mesh, w_ph1, rho1f);
plicFuncs::linearInterpolate_2ph(rho0, mesh, w_ph0, rho0f);

for(label i=0; i<nSpecies; i++)
{
    const volScalarField& D1i = D1[i];
    surfaceScalarField& D1fi = D1f[i];
    plicFuncs::linearInterpolate_2ph(D1i, mesh, w_ph1, D1fi);

    const volScalarField& D0i = D0[i];
    surfaceScalarField& D0fi = D0f[i];
    plicFuncs::linearInterpolate_2ph(D0i, mesh, w_ph0, D0fi);
}
//Done linear interpolate density and diffusivity to faces in both phases
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Calculate species mass fraction gradients for both phases
for(label i=0; i<nSpecies; i++)
{
    const volScalarField& Y1i = Y1[i];
    const volScalarField& Y0i = Y0[i];

    List<scalar>& Y1i_flatFld_diff = Y1_flatFld_diff[i];
    List<scalar>& Y0i_flatFld_diff = Y0_flatFld_diff[i];
    plicFuncs::makeFlatFld(Y1i, diffMap, Y1i_flatFld_diff);
    plicFuncs::makeFlatFld(Y0i, diffMap, Y0i_flatFld_diff);

    surfaceScalarField& gradf_Y1i = gradf_Y1[i];
    surfaceScalarField& gradf_Y0i = gradf_Y0[i];

    plicFuncs::calc_2ph_gradf
    (
        mesh,
        diffCellStencil,
        Y1i,
        Y0i,
        Y1i_flatFld_diff,
        Y0i_flatFld_diff,
        alpha1_flatFld_diff,
        C_ph1_flatFld_diff,
        C_ph0_flatFld_diff,
        Af_ph1,
        Af_ph0,
        gradf_Y1i,
        gradf_Y0i,
        i,
        MIN_ALPHA_DIFF,
        diff_debug2
    );
}
//Done calculate species mass fraction gradients for both phases
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Calculate species diffusive fluxes for both phases
for(label i=0; i<nSpecies; i++)
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
        Af_ph1,
        Af_ph0,
        face_phaseState,
        rho1,
        alpha1,
        Y1i,
        rho0,
        alpha0,
        Y0i,
        diffFlux_Y1i,
        diffFlux_Y0i,
        i,
        diff_debug2
    );
}
//Done calculate species diffusive fluxes for both phases
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Solve species diffusion equations in both phases
for(label i=0; i<(nSpecies - 1); i++)
{
    if(diff_debug2)
    {
        Info<< "-----------------------------------" << nl
            << "              Field Y[" << i << "]" << nl
            << "-----------------------------------" << nl
            << endl;
    }

    volScalarField& Y1i = Y1[i];
    volScalarField& Y0i = Y0[i];
    volScalarField& C1i = C1[i];
    volScalarField& C0i = C0[i];

    scalarField& Y1iCells = Y1i.internalField();
    scalarField& Y0iCells = Y0i.internalField();
    scalarField& C1iCells = C1i.internalField();
    scalarField& C0iCells = C0i.internalField();

    const surfaceScalarField& diffFlux_Y1i = diffFlux_Y1[i];
    const surfaceScalarField& diffFlux_Y0i = diffFlux_Y0[i];

    //--------------------------------------------------------------//
    //ph1
    volScalarField surfInt_diffFlux_Y1i(fvc::surfaceIntegrate(diffFlux_Y1i));    
    scalarField& surfInt_diffFlux_Y1iCells = surfInt_diffFlux_Y1i.internalField();

    forAll(C1iCells,cellI)
    {
        if(diff_debug2)
        {
            Info<< "Solving Y1[" << i << "] eqn in cell " << cellI << nl
                << "Y1i = " << Y1iCells[cellI] << "  alpha1 = " << alpha1Cells[cellI] << nl
                << "div(diffFlux) = " << surfInt_diffFlux_Y1iCells[cellI] << "  div(diffFlux)*dt = " << surfInt_diffFlux_Y1iCells[cellI]*runTime.deltaTValue()
                << endl;
        }
                
        C1iCells[cellI] += -surfInt_diffFlux_Y1iCells[cellI]*runTime.deltaTValue();

        if(alpha1Cells[cellI] < SMALLEST_ALPHA_FOR_YDIFF)
        {
            Y1iCells[cellI] = 0;
        }
        else
        {
            Y1iCells[cellI] = C1iCells[cellI]/alpha1Cells[cellI]/rho1Cells[cellI];
        }

        if(diff_debug2)
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
    //done ph1
    //--------------------------------------------------------------//

    //--------------------------------------------------------------//
    //ph0
    volScalarField surfInt_diffFlux_Y0i(fvc::surfaceIntegrate(diffFlux_Y0i));    
    scalarField& surfInt_diffFlux_Y0iCells = surfInt_diffFlux_Y0i.internalField();

    forAll(C0iCells,cellI)
    {
        if(diff_debug2)
        {
            Info<< "Solving Y0[" << i << "] eqn in cell " << cellI << nl
                << "Y0i = " << Y0iCells[cellI] << "  alpha0 = " << alpha0Cells[cellI] << nl
                << "div(diffFlux) = " << surfInt_diffFlux_Y0iCells[cellI] << "  div(diffFlux)*dt = " << surfInt_diffFlux_Y0iCells[cellI]*runTime.deltaTValue()
                << endl;
        }
                
        C0iCells[cellI] += -surfInt_diffFlux_Y0iCells[cellI]*runTime.deltaTValue();

        if(alpha0Cells[cellI] < SMALLEST_ALPHA_FOR_YDIFF)
        {
            Y0iCells[cellI] = 0;
        }
        else
        {
            Y0iCells[cellI] = C0iCells[cellI]/alpha0Cells[cellI]/rho0Cells[cellI];
        }

        if(diff_debug2)
        {
            Info<< "C0i = " << C0iCells[cellI] << nl
                << "Y0i = " << Y0iCells[cellI] << nl
                << endl;
        }
    }

    Info<< "Phase-0 Yi= "
        << (C0i).weightedAverage(mesh.Vsc()).value()
        << nl
        << "Min(C0i) = " << min(C0i).value()
        << "  Max(C0i) = " << max(C0i).value()
        << nl
        << "Min(Y0i) = " << min(Y0i).value()
        << "  Max(Y0i) = " << max(Y0i).value()
        << nl
        << "C0i undershoot = " << min(min(C0i).value(), 0)
        << "  C0i overshoot = " << max((max(C0i).value() - rho0_0.value()), 0)
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
        }
        
        if(Y0iCells[cellI] > 1)
        {
            Y0iCells[cellI] = 1;
            C0iCells[cellI] = rho0Cells[cellI]*alpha0Cells[cellI]*Y0iCells[cellI];
        }        
    }

    Info<< "Phase-0 Yi after correction= "
        << (C0i).weightedAverage(mesh.Vsc()).value()
        << nl
        << "Min(C0i) = " << min(C0i).value()
        << "  Max(C0i) = " << max(C0i).value()
        << nl
        << "Min(Y0i) = " << min(Y0i).value()
        << "  Max(Y0i) = " << max(Y0i).value()        
        << endl;

    C0i.correctBoundaryConditions();
    Y0i.correctBoundaryConditions();
    //done ph0
    //--------------------------------------------------------------//

    Y1t += Y1i;
    Y0t += Y0i;
}
//Done solve species diffusion equations in both phases
//------------------------------------------------------------------//

Y1[nSpecies - 1] = scalar(1.0) - Y1t;
Y1[nSpecies - 1].max(0.0);
Y1[nSpecies - 1].min(1.0);
C1[nSpecies - 1] = rho1*alpha1*Y1[nSpecies - 1];

Y0[nSpecies - 1] = scalar(1.0) - Y0t;
Y0[nSpecies - 1].max(0.0);
Y0[nSpecies - 1].min(1.0);
C0[nSpecies - 1] = rho0*alpha0*Y0[nSpecies - 1];

Info<< "Done YEqns" << endl;
