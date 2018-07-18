
if(grad_debug)
{
    plicFuncs::print_line(osGrad, 100);
    osGrad<< "Solution Time = " << runTime.timeName() << nl 
        <<"Gradients and Interpolation for diffusion" << endl;
    plicFuncs::print_line(osGrad, 100);
    osGrad<< endl;
}

//------------------------------------------------------------------//
//Make flat fields for alpha1, phase centroids, x, T required for gradient calculation
plicFuncs::makeFlatFld(alpha1, diffMap, alpha1_flatFld_diff);
plicFuncs::makeFlatFld(interface.C_ph1(), diffMap, C_ph1_flatFld_diff);
plicFuncs::makeFlatFld(interface.C_ph0(), diffMap, C_ph0_flatFld_diff);
plicFuncs::makeFlatFld(T1, diffMap, T1_flatFld_diff);
plicFuncs::makeFlatFld(T0, diffMap, T0_flatFld_diff);
for(i=0; i<n; i++)
{
    const volScalarField& Y1i = Y1[i];
    const volScalarField& Y0i = Y0[i];
    List<scalar>& Y1i_flatFld_diff = Y1_flatFld_diff[i];
    List<scalar>& Y0i_flatFld_diff = Y0_flatFld_diff[i];
    plicFuncs::makeFlatFld(Y1i, diffMap, Y1i_flatFld_diff);
    plicFuncs::makeFlatFld(Y0i, diffMap, Y0i_flatFld_diff);
}
//Done make flat fields for alpha1, phase centroids, x, T required for gradient calculation
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Determine phase state at faces for diffusion flux calculation
plicFuncs::calc_face_phaseState_diff(mesh, Y1[0], Af_ph1_own, Af_ph0_own, Af_ph1_nei, Af_ph0_nei, MIN_ALPHA_DIFF, face_phaseState_diff, grad_debug, osGrad);
//Done determine phase state at faces for diffusion flux calculation
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Linear interpolate T, v, lambda, mu, diffusivity, x to faces in both phases
plicFuncs::calc_2ph_Cf(Cf_ph1_own, Cf_ph0_own, Cf_ph1_nei, Cf_ph0_nei, Af_ph1_own, Af_ph0_own, Af_ph1_nei, Af_ph0_nei, Cf_ph1, Cf_ph0);

plicFuncs::calc_2ph_linearInterpolation_weights(mesh, C_ph1, Cf_ph1, w_ph1);
plicFuncs::calc_2ph_linearInterpolation_weights(mesh, C_ph0, Cf_ph0, w_ph0);

plicFuncs::linearInterpolate_2ph(rho1, mesh, w_ph1, rho1f);
plicFuncs::linearInterpolate_2ph(rho0, mesh, w_ph0, rho0f);
plicFuncs::linearInterpolate_2ph(mu1, mesh, w_ph1, mu1f, interp_debug, osGrad);
plicFuncs::linearInterpolate_2ph(mu0, mesh, w_ph0, mu0f, interp_debug, osGrad);

for(i=0; i<n; i++)
{        
    const volScalarField& D1i = D1[i];
    const volScalarField& D0i = D0[i];
    surfaceScalarField& D1fi = D1f[i];
    surfaceScalarField& D0fi = D0f[i];
    plicFuncs::linearInterpolate_2ph(D1i, mesh, w_ph1, D1fi);
    plicFuncs::linearInterpolate_2ph(D0i, mesh, w_ph1, D0fi);
}
//Done linear interpolate T, v, lambda, mu, diffusivity, x to faces in both phases
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Calculate species mole fraction and temperature gradients for both phases
plicFuncs::calc_2ph_gradf
(
    mesh,
    diffCellStencil,
    Y1,
    Y0,
    T1,
    T0,
    Y1_flatFld_diff,
    Y0_flatFld_diff,
    T1_flatFld_diff,
    T0_flatFld_diff,
    alpha1_flatFld_diff,
    C_ph1_flatFld_diff,
    C_ph0_flatFld_diff,
    Cf_ph1_own,
    Cf_ph1_nei,
    Cf_ph0_own,
    Cf_ph0_nei,
    face_phaseState_diff,
    gradf_Y1,
    gradf_Y0,
    gradf_T1,
    gradf_T0,
    n,
    MIN_ALPHA_DIFF,
    grad_debug,
    osGrad
);
//Done calculate species mole fraction and temperature gradients for both phases
//------------------------------------------------------------------//

Info<< "Done gradients and face interpolations for diffusion" << endl;
