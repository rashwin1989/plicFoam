
if(ist_debug)
{
    osIst<< nl
         << "Solution Time = " << runTime.timeName()
         << nl << endl;
}

//------------------------------------------------------------------//
//Make flat fields for alpha1 and phase centroids required for gradient calculation
plicFuncs::makeFlatFld(alpha1, istMap, alpha1_flatFld_ist);
plicFuncs::makeFlatFld(interface.C_ph1(), istMap, C_ph1_flatFld_ist);
plicFuncs::makeFlatFld(interface.C_ph0(), istMap, C_ph0_flatFld_ist);

for(label i=0; i<nSpecies; i++)
{
    List<scalar>& Y1i_flatFld_ist = Y1_flatFld_ist[i];
    List<scalar>& Y0i_flatFld_ist = Y0_flatFld_ist[i];
    plicFuncs::makeFlatFld(Y1i, istMap, Y1i_flatFld_ist);    
    plicFuncs::makeFlatFld(Y1i, istMap, Y0i_flatFld_ist);    
}
//Done make flat fields for alpha1 and phase centroids required for gradient calculation
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Calculate interfacial species mass fluxes
plicFuncs::calc_Js(mesh, cell_stencil, Y1_flatFld_ist, Y0_flatFld_ist, alpha1_flatFld_ist, rho1, rho0, D1, D0, C_ph1_flatFld_ist, C_ph0_flatFld_ist, C_intfc, A_intfc, nHat, Ys1, Ys0, Js1, Js0, nSpecies, ALPHA_2PH_MIN, ist_debug, osIst);
//Done calculate interfacial species mass fluxes
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Calculate interfacial species transfer source terms and alpha source terms
plicFuncs::calc_mS_alphaS(mesh, C1, C0, alpha1, rho1, rho0, Js1, Js0, nSpecies, ALPHA_2PH_MIN, mS1, mS0, mS1Tot, mS0Tot, alphaS1, alphaS0, ist_debug, osIst);
//Done calculate interfacial species transfer source terms and alpha source terms
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Apply interfacial species transfer and alpha source terms, calculate source term for pressure equation

//------------------------------------------------------------------//
//ph1
forAll(alpha1.internalField(), cellI)

//Done apply interfacial species transfer and alpha source terms, calculate source term for pressure equation
//------------------------------------------------------------------//


Info<< "Done Interfacial Species Transfer" << endl;
