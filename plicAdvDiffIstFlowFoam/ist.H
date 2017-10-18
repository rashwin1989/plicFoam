
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
plicFuncs::calc_Js(mesh, cell_stencil, Y1_flatFld_ist, Y0_flatFld_ist, alpha1_flatFld_ist, rho1, rho0, D1, D0, C_ph1_flatFld_ist, C_ph0_flatFld_ist, C_intfc, A_intfc, nHat, Ys1, Ys0, Js1, Js0, nSpecies, ALPHA_2PH_MIN, ist_debug2, osIst);
//Done calculate interfacial species mass fluxes
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Calculate interfacial species transfer source terms and alpha source terms
plicFuncs::calc_mS_alphaS(mesh, C1, C0, alpha1, rho1, rho0, Js1, Js0, nSpecies, ALPHA_2PH_MIN, mS1, mS0, mS1Tot, mS0Tot, alphaS1, alphaS0, ist_debug2, osIst);
//Done calculate interfacial species transfer source terms and alpha source terms
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Apply interfacial species transfer and alpha source terms, calculate source term for pressure equation
scalar C1Tot_cellI;
scalar C0Tot_cellI;
forAll(alpha1.internalField(), cellI)
{
    if(ist_debug)
    {        
        osIst<< "--------------------------------------------------------------------------" << nl
            << "Cell: " << cellI << nl
            << "--------------------------------------------------------------------------" << nl
            << "Before interfacial species transfer" << nl
            << "alpha1 = " << alpha1.internalField()[cellI] << "  alpha0 = " << alpha0.internalField()[cellI]
            << "  (alpha1+alpha0) = " << (alpha1.internalField()[cellI] + alpha1.internalField()[cellI]) << nl
            << "rho1 = " << rho1.internalField()[cellI] << "rho0 = " << rho0.internalField()[cellI] << endl;
        for(label i=0; i<nSpecies; i++)
        {
            osIst<< "C1[" << i << "] = " << C1[i].internalField()[cellI] << "  Y1[" << i << "] = " << Y1[i].internalField()[cellI]
                << "  C0[" << i << "] = " << C0[i].internalField()[cellI] << "  Y0[" << i << "] = " << Y0[i].internalField()[cellI]
                << endl;
        }
        osIst<< "--------------------------------------------------------------------------" << nl
            << "Interfacial source terms" << nl
            << "--------------------------------------------------------------------------" << nl
            << "mS1Tot = " << mS1Tot.internalField()[cellI] << "mS0Tot = " << mS0Tot.internalField()[cellI] << nl
            << "alphaS1 = " << alphaS1.internalField()[cellI] << "alphaS0 = " << alphaS0.internalField()[cellI] << endl;
        for(label i=0; i<nSpecies; i++)
        {
            osIst<< "mS1[" << i << "] = " << mS1[i].internalField()[cellI] << "  mS0[" << i << "] = " << mS0[i].internalField()[cellI]
                << endl;
        }
    }

    C1Tot_cellI = 0;
    C0Tot_cellI = 0;

    for(label i=0; i<nSpecies; i++)
    {
        C1[i].internalField()[cellI] += mS1[i].internalField()[cellI];
        C0[i].internalField()[cellI] += mS0[i].internalField()[cellI];
        C1Tot_cellI += C1[i].internalField()[cellI];
        C0Tot_cellI += C0[i].internalField()[cellI];
    }

    alpha1.internalField()[cellI] += alphaS1.internalField()[cellI];
    alpha0.internalField()[cellI] += alphaS0.internalField()[cellI];

    scalar pSCoeff = 1/rho1.internalField()[cellI] - 1/rho0.internalField()[cellI];
    pS.internalField()[cellI] = pSCoeff*mS1Tot.internalField()[cellI]/runTime.deltaTValue();

    if(C1Tot_cellI < SMALL)
    {
        C1Tot_cellI += SMALL;
    }
    if(C0Tot_cellI < SMALL)
    {
        C0Tot_cellI += SMALL;
    }

    for(label i=0; i<nSpecies; i++)
    {        
        Y1[i].internalField()[cellI] = C1[i].internalField()[cellI]/C1Tot_cellI;
        Y0[i].internalField()[cellI] = C0[i].internalField()[cellI]/C0Tot_cellI;               
    }

    if(ist_debug)
    {        
        osIst<< "pS = " << pS.internalField()[cellI] << nl
            << "--------------------------------------------------------------------------" << nl
            << "After interfacial species transfer" << nl
            << "alpha1 = " << alpha1.internalField()[cellI] << "  alpha0 = " << alpha0.internalField()[cellI]
            << "  (alpha1+alpha0) = " << (alpha1.internalField()[cellI] + alpha1.internalField()[cellI]) << nl
            << "rho1 = " << rho1.internalField()[cellI] << "rho0 = " << rho0.internalField()[cellI] << endl;
        for(label i=0; i<nSpecies; i++)
        {
            osIst<< "C1[" << i << "] = " << C1[i].internalField()[cellI] << "  Y1[" << i << "] = " << Y1[i].internalField()[cellI]
                << "  C0[" << i << "] = " << C0[i].internalField()[cellI] << "  Y0[" << i << "] = " << Y0[i].internalField()[cellI]
                << endl;
        }        
        osIst<< "--------------------------------------------------------------------------" << nl
            << "--------------------------------------------------------------------------" << nl
            << endl;
    }
}

//Done apply interfacial species transfer and alpha source terms, calculate source term for pressure equation
//------------------------------------------------------------------//


Info<< "Done Interfacial Species Transfer" << endl;
