
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
    plicFuncs::makeFlatFld(Y1[i], istMap, Y1i_flatFld_ist);    
    plicFuncs::makeFlatFld(Y0[i], istMap, Y0i_flatFld_ist);    
}
//Done make flat fields for alpha1 and phase centroids required for gradient calculation
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Calculate interfacial species mass fluxes
//plicFuncs::calc_Js(mesh, istCellStencil, Y1_flatFld_ist, Y0_flatFld_ist, alpha1_flatFld_ist, alpha1, rho1, rho0, D1, D0, C_ph1_flatFld_ist, C_ph0_flatFld_ist, interface.C_intfc(), interface.A_intfc(), interface.nHat(), Ys1, Ys0, Js1, Js0, nSpecies, ALPHA_2PH_MIN, ist_debug2, osIst);
//Done calculate interfacial species mass fluxes
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Calculate interfacial species transfer source terms and alpha source terms
plicFuncs::calc_mS_He(mesh, istCellStencil, Y1_flatFld_ist, Y0_flatFld_ist, alpha1_flatFld_ist, C1, C0, Y1, Y0, alpha1, rho1, rho0, D1, D0, C_ph1_flatFld_ist, C_ph0_flatFld_ist, interface.C_intfc(), interface.A_intfc(), interface.nHat(), He, Ys1, Ys0, Js1, Js0, nSpecies, MIN_ALPHA_IST, A_INTFC_2PH_MIN, mS1, mS0, ist_debug2, osIst);
//Done calculate interfacial species transfer source terms and alpha source terms
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Apply interfacial species transfer and alpha source terms, calculate source term for pressure equation
scalar dt = runTime.deltaTValue();
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
            << "  (alpha1+alpha0) = " << (alpha1.internalField()[cellI] + alpha0.internalField()[cellI]) << nl
            << "rho1 = " << rho1.internalField()[cellI] << "  rho0 = " << rho0.internalField()[cellI] << endl;
        for(label i=0; i<nSpecies; i++)
        {
            osIst<< "C1[" << i << "] = " << C1[i].internalField()[cellI] << "  Y1[" << i << "] = " << Y1[i].internalField()[cellI]
                << "  C0[" << i << "] = " << C0[i].internalField()[cellI] << "  Y0[" << i << "] = " << Y0[i].internalField()[cellI]
                << endl;
        }
        osIst<< "--------------------------------------------------------------------------" << nl
            << "Interfacial source terms" << nl
            << "--------------------------------------------------------------------------" << endl;            
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
        C1[i].internalField()[cellI] += mS1[i].internalField()[cellI]*dt;
        C1[i].internalField()[cellI] = max(C1[i].internalField()[cellI], 0);
        C0[i].internalField()[cellI] += mS0[i].internalField()[cellI]*dt;
        C0[i].internalField()[cellI] = max(C0[i].internalField()[cellI], 0);
        C1Tot_cellI += C1[i].internalField()[cellI];
        C0Tot_cellI += C0[i].internalField()[cellI];
    }

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
        osIst<< "--------------------------------------------------------------------------" << nl
            << "After interfacial species transfer" << nl            
            << "rho1 = " << rho1.internalField()[cellI] << "  rho0 = " << rho0.internalField()[cellI] << endl;
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

for(label i=0; i<nSpecies; i++)
{
    C1[i].correctBoundaryConditions();
    Y1[i].correctBoundaryConditions();

    C0[i].correctBoundaryConditions();
    Y0[i].correctBoundaryConditions();

    Info<< "Phase-1 C[" << i << "] = "
        << (C1[i]).weightedAverage(mesh.Vsc()).value()
        << nl
        << "Min(C1i) = " << min(C1[i]).value()
        << "  Max(C1i) = " << max(C1[i]).value()
        << nl
        << "Min(Y1i) = " << min(Y1[i]).value()
        << "  Max(Y1i) = " << max(Y1[i]).value()
        << endl;

    Info<< "Phase-0 C[" << i << "] = "
        << (C0[i]).weightedAverage(mesh.Vsc()).value()
        << nl
        << "Min(C0i) = " << min(C0[i]).value()
        << "  Max(C0i) = " << max(C0[i]).value()
        << nl
        << "Min(Y0i) = " << min(Y0[i]).value()
        << "  Max(Y0i) = " << max(Y0[i]).value()
        << endl;        
}

//Done apply interfacial species transfer and alpha source terms, calculate source term for pressure equation
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//Reset the allowable limits for concentrations and mass fractions in the phases
for(label i=0; i<nSpecies; i++)
{
    C1max[i] = max(C1[i]).value();
    C0max[i] = max(C0[i]).value();
    C1min[i] = min(C1[i]).value();
    C0min[i] = min(C0[i]).value();

    Y1max[i] = max(Y1[i]).value();
    Y0max[i] = max(Y0[i]).value();
    Y1min[i] = min(Y1[i]).value();
    Y0min[i] = min(Y0[i]).value();
}

Info<< "Done Interfacial Species Transfer" << endl;
