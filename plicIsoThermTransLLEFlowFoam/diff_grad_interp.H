
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
    const volScalarField& x1i = x1[i];
    const volScalarField& x0i = x0[i];
    List<scalar>& x1i_flatFld_diff = x1_flatFld_diff[i];
    List<scalar>& x0i_flatFld_diff = x0_flatFld_diff[i];
    plicFuncs::makeFlatFld(x1i, diffMap, x1i_flatFld_diff);
    plicFuncs::makeFlatFld(x0i, diffMap, x0i_flatFld_diff);
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

plicFuncs::linearInterpolate_2ph(T1, mesh, w_ph1, T1f);
plicFuncs::linearInterpolate_2ph(T0, mesh, w_ph0, T0f);
plicFuncs::linearInterpolate_2ph(v1, mesh, w_ph1, v1f);
plicFuncs::linearInterpolate_2ph(v0, mesh, w_ph0, v0f);
plicFuncs::linearInterpolate_2ph(mu1, mesh, w_ph1, mu1f);
plicFuncs::linearInterpolate_2ph(mu0, mesh, w_ph0, mu0f);

for(i=0; i<n; i++)
{
    const volScalarField& x1i = x1[i];
    const volScalarField& x0i = x0[i];
    surfaceScalarField& x1fi = x1f[i];
    surfaceScalarField& x0fi = x0f[i];
    plicFuncs::linearInterpolate_2ph(x1i, mesh, w_ph1, x1fi);
    plicFuncs::linearInterpolate_2ph(x0i, mesh, w_ph1, x0fi);    

    for(j=0; j<n; j++)
    {
        idx = i*n + j;
        const volScalarField& D1ij = Dij1[idx];
        const volScalarField& D0ij = Dij0[idx];
        surfaceScalarField& D1fij = Dij1f[idx];
        surfaceScalarField& D0fij = Dij0f[idx];
        plicFuncs::linearInterpolate_2ph(D1ij, mesh, w_ph1, D1fij);
        plicFuncs::linearInterpolate_2ph(D0ij, mesh, w_ph1, D0fij);
    }
}

if(grad_debug2)
{        
    plicFuncs::print_line(osGrad, 120);
    plicFuncs::print_line(osGrad, 120);
    osGrad<< "Cell values ph-1: alpha1, y, x, v, rho" << endl;
    plicFuncs::print_line(osGrad, 120);
    osGrad<< setw(6) << "cellI" << "  " << setw(8) << "alpha1" << "  ";
    for(i=0; i<n; i++)
    {
        word y1iName = "y1["+name(i)+"]";
        osGrad<< setw(8) << y1iName << "  ";
    }
    for(i=0; i<n; i++)
    {
        word x1iName = "x1["+name(i)+"]";
        osGrad<< setw(8) << x1iName << "  ";
    }
    osGrad<< setw(8) << "v1" << "  " << setw(8) << "rho1" << endl;
    plicFuncs::print_line(osGrad, 120);
    forAll(alpha1Cells, cellI)
    {
        osGrad<< setw(6) << cellI << "  " << setw(8) << alpha1Cells[cellI] << "  ";
        for(i=0; i<n; i++)
        {
            osGrad<< setw(8) << Y1[i].internalField()[cellI] << "  ";
        }
        for(i=0; i<n; i++)
        {
            osGrad<< setw(8) << x1[i].internalField()[cellI] << "  ";
        }
        osGrad<< setw(8) << v1.internalField()[cellI] << "  " << setw(8) << rho1.internalField()[cellI] << endl;
    }
    plicFuncs::print_line(osGrad, 120);
    plicFuncs::print_line(osGrad, 120);
    osGrad<< endl;

    plicFuncs::print_line(osGrad, 120);
    plicFuncs::print_line(osGrad, 120);
    osGrad<< "Cell values ph-1: alpha1, D" << endl;
    plicFuncs::print_line(osGrad, 120);
    osGrad<< setw(6) << "cellI" << "  " << setw(8) << "alpha1" << "  ";
    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
        {
            word D1iName = "D1("+name(i)+","+name(j)+")";
            osGrad<< setw(8) << D1iName << "  ";
        }
    }    
    osGrad<< endl;
    plicFuncs::print_line(osGrad, 120);
    forAll(alpha1Cells, cellI)
    {
        osGrad<< setw(6) << cellI << "  " << setw(8) << alpha1Cells[cellI] << "  ";
        for(i=0; i<n; i++)
        {
            for(j=0; j<n; j++)
            {
                idx = i*n + j;
                osGrad<< setw(8) << Dij1[idx].internalField()[cellI] << "  ";
            }
        }
        osGrad<< endl;
    }
    plicFuncs::print_line(osGrad, 120);
    plicFuncs::print_line(osGrad, 120);
    osGrad<< endl;


    plicFuncs::print_line(osGrad, 120);
    plicFuncs::print_line(osGrad, 120);
    osGrad<< "Face values ph-1: alpha1, y, x, v, rho" << endl;
    plicFuncs::print_line(osGrad, 120);
    osGrad<< setw(6) << "faceI" << "  " << setw(8) << "alpha1" << "  ";    
    for(i=0; i<n; i++)
    {
        word x1iName = "x1["+name(i)+"]";
        osGrad<< setw(8) << x1iName << "  ";
    }
    osGrad<< setw(8) << "v1" << endl;
    plicFuncs::print_line(osGrad, 120);
    for(faceI=0; faceI<mesh.nInternalFaces(); faceI++)
    {
        osGrad<< setw(6) << faceI << "  " << setw(8) << alpha1f[faceI] << "  ";
        for(i=0; i<n; i++)
        {
            osGrad<< setw(8) << x1f[i][faceI] << "  ";
        }
        osGrad<< setw(8) << v1f[faceI] << endl;
    }
    plicFuncs::print_line(osGrad, 120);
    plicFuncs::print_line(osGrad, 120);
    osGrad<< endl;

    plicFuncs::print_line(osGrad, 120);
    plicFuncs::print_line(osGrad, 120);
    osGrad<< "Face values ph-1: alpha1, D" << endl;
    plicFuncs::print_line(osGrad, 120);
    osGrad<< setw(6) << "faceI" << "  " << setw(8) << "alpha1" << "  ";
    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
        {
            word D1iName = "D1("+name(i)+","+name(j)+")";
            osGrad<< setw(8) << D1iName << "  ";
        }
    }    
    osGrad<< endl;
    plicFuncs::print_line(osGrad, 120);
    for(faceI=0; faceI<mesh.nInternalFaces(); faceI++)
    {
        osGrad<< setw(6) << faceI << "  " << setw(8) << alpha1f[faceI] << "  ";
        for(i=0; i<n; i++)
        {
            for(j=0; j<n; j++)
            {
                idx = i*n + j;
                osGrad<< setw(8) << Dij1f[idx][faceI] << "  ";
            }
        }
        osGrad<< endl;
    }
    plicFuncs::print_line(osGrad, 120);
    plicFuncs::print_line(osGrad, 120);
    osGrad<< endl;
}
//Done linear interpolate T, v, lambda, mu, diffusivity, x to faces in both phases
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Calculate species mole fraction and temperature gradients for both phases
plicFuncs::calc_2ph_gradf
(
    mesh,
    diffCellStencil,
    x1,
    x0,
    T1,
    T0,
    x1_flatFld_diff,
    x0_flatFld_diff,
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
    gradf_x1,
    gradf_x0,
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
