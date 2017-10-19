
if(diff_debug)
{
    osDiff<< nl
          << "Solution Time = " << runTime.timeName()
          << nl << endl;
}

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
//Determine phase state at faces for diffusion flux calculation
plicFuncs::calc_face_phaseState_diff(mesh, Y1[0], Af_ph1_own, Af_ph0_own, Af_ph1_nei, Af_ph0_nei, MIN_ALPHA_DIFF, face_phaseState_diff, diff_debug, osDiff);
//Done determine phase state at faces for diffusion flux calculation
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Linear interpolate density and diffusivity to faces in both phases
plicFuncs::calc_2ph_Cf(Cf_ph1_own, Cf_ph0_own, Cf_ph1_nei, Cf_ph0_nei, Af_ph1_own, Af_ph0_own, Af_ph1_nei, Af_ph0_nei, Cf_ph1, Cf_ph0);

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
        face_phaseState_diff,
        gradf_Y1i,
        gradf_Y0i,
        i,
        MIN_ALPHA_DIFF,
        diff_debug2,
        osDiff
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
        diff_debug2,
        osDiff
    );
}
//Done calculate species diffusive fluxes for both phases
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Solve species diffusion equations in both phases
for(label i=0; i<(nSpecies - 1); i++)
{
    if(diff_debug)
    {
        osDiff<< "-----------------------------------" << nl
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
        if(diff_debug)
        {
            const cell& curCell = meshCells[cellI];

            osDiff<< "Solving Y1[" << i << "] diff eqn in cell " << cellI << nl
                << "C1i = " << C1iCells[cellI] << "  Y1i = " << Y1iCells[cellI] << "  alpha1 = " << alpha1Cells[cellI] << "  alpha1 oldTime = " << alpha1.oldTime().internalField()[cellI] << endl;

            forAll(curCell, faceI)
            {
                label curFaceLbl = curCell[faceI];
                label curOwn = own[curFaceLbl];

                osDiff<< "cell face label: " << curFaceLbl 
                    << nl
                    <<"face own: " << curOwn << "  alpha1: " << alpha1Cells[curOwn] << "  C1i: " << C1iCells[curOwn] << "  Y1i: " << Y1iCells[curOwn] 
                    << endl;
                if(curFaceLbl < mesh.nInternalFaces())
                {
                    label curNei = nei[curFaceLbl];
                    osDiff<< "face nei: " << curNei << "  alpha1: " << alpha1Cells[curNei] << "  C1i: " << C1iCells[curNei] << "  Y1i: " << Y1iCells[curNei]
                        << nl
                        << "Af_ph1_own: " << Af_ph1_own[curFaceLbl] << "  Af_ph1_nei: " << Af_ph1_nei[curFaceLbl] << "  gradf_Y1i: " << gradf_Y1[i][curFaceLbl] << "  diffFlux_Y1i: " << diffFlux_Y1i[curFaceLbl] 
                        << endl;
                }        
                else
                {
                    osDiff<< "Af_ph1_own: " << Af_ph1_own[curFaceLbl] << "  Af_ph1_nei: " << Af_ph1_nei[curFaceLbl]
                        << endl;
                }    
                
            }

            osDiff<< "cell V: " << meshV[cellI] << endl;

            osDiff<< "div(diffFlux) = " << surfInt_diffFlux_Y1iCells[cellI] << "  div(diffFlux)*dt = " << surfInt_diffFlux_Y1iCells[cellI]*runTime.deltaTValue()
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

        if(diff_debug)
        {
            osDiff<< "C1i = " << C1iCells[cellI] << nl
                << "Y1i = " << Y1iCells[cellI] << nl
                << endl;
        }
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
            //C1iCells[cellI] = rho1Cells[cellI]*alpha1Cells[cellI]*Y1iCells[cellI];
        }
        
        if(Y1iCells[cellI] > Y1max[i])
        {
            Y1iCells[cellI] = Y1max[i];
            //C1iCells[cellI] = rho1Cells[cellI]*alpha1Cells[cellI]*Y1iCells[cellI];
        }

        C1iCells[cellI] = rho1Cells[cellI]*alpha1Cells[cellI]*Y1iCells[cellI];
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
    Y1i.correctBoundaryConditions();
    //done ph1
    //--------------------------------------------------------------//

    //--------------------------------------------------------------//
    //ph0
    volScalarField surfInt_diffFlux_Y0i(fvc::surfaceIntegrate(diffFlux_Y0i));    
    scalarField& surfInt_diffFlux_Y0iCells = surfInt_diffFlux_Y0i.internalField();

    forAll(C0iCells,cellI)
    {
        if(diff_debug)
        {
            const cell& curCell = meshCells[cellI];

            osDiff<< "Solving Y0[" << i << "] diff eqn in cell " << cellI << nl
                << "C0i = " << C0iCells[cellI] << "  Y0i = " << Y0iCells[cellI] << "  alpha0 = " << alpha0Cells[cellI] << "  alpha0 oldTime = " << alpha0.oldTime().internalField()[cellI] << endl;

            forAll(curCell, faceI)
            {
                label curFaceLbl = curCell[faceI];
                label curOwn = own[curFaceLbl];

                osDiff<< "cell face label: " << curFaceLbl 
                    << nl
                    << "face own: " << curOwn << "  alpha0: " << alpha0Cells[curOwn] << "  C0i: " << C0iCells[curOwn] << "  Y0i: " << Y0iCells[curOwn] 
                    << endl;
                if(curFaceLbl < mesh.nInternalFaces())
                {
                    label curNei = nei[curFaceLbl];
                    osDiff<< "face nei: " << curNei << "  alpha0: " << alpha0Cells[curNei] << "  C0i: " << C0iCells[curNei] << "  Y0i: " << Y0iCells[curNei]
                        << nl
                        << "Af_ph0_own: " << Af_ph0_own[curFaceLbl] << "  Af_ph0_nei: " << Af_ph0_nei[curFaceLbl] << "  gradf_Y0i: " << gradf_Y0[i][curFaceLbl] << "  diffFlux_Y0i: " << diffFlux_Y0i[curFaceLbl] 
                        << endl;                        
                }       
                else
                {
                    osDiff<< "Af_ph0_own: " << Af_ph0_own[curFaceLbl] << "  Af_ph0_nei: " << Af_ph0_nei[curFaceLbl]
                        << endl;
                }
            }

            osDiff<< "cell V: " << meshV[cellI] << endl;

            osDiff<< "div(diffFlux) = " << surfInt_diffFlux_Y0iCells[cellI] << "  div(diffFlux)*dt = " << surfInt_diffFlux_Y0iCells[cellI]*runTime.deltaTValue()
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

        if(diff_debug)
        {
            osDiff<< "C0i = " << C0iCells[cellI] << nl
                << "Y0i = " << Y0iCells[cellI] << nl
                << endl;
        }
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
            //C0iCells[cellI] = rho0Cells[cellI]*alpha0Cells[cellI]*Y0iCells[cellI];
        }
        
        if(Y0iCells[cellI] > Y0max[i])
        {
            Y0iCells[cellI] = Y0max[i];
            //C0iCells[cellI] = rho0Cells[cellI]*alpha0Cells[cellI]*Y0iCells[cellI];
        }

        C0iCells[cellI] = rho0Cells[cellI]*alpha0Cells[cellI]*Y0iCells[cellI];
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

Info<< "Done YDiffEqns" << endl;
