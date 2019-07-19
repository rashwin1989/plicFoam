// Ist file for 2-component, non-isothermal interfacial mass transfer

if(ist_debug)
{
    plicFuncs::print_line(osIst, 100);
    osIst<< "Solution Time = " << runTime.timeName() << "  Interfacial transfer" << endl;
    plicFuncs::print_line(osIst, 100);
    osIst<< endl;
}

dt = deltaT;

//-----------------------------------------------------------------------------------------------------//
//Make flat fields for alpha1 and phase centroids required for gradient calculation
plicFuncs::makeFlatFld(alpha1, istMap, alpha1_flatFld_ist);
plicFuncs::makeFlatFld(interface.C_ph1(), istMap, C_ph1_flatFld_ist);
plicFuncs::makeFlatFld(interface.C_ph0(), istMap, C_ph0_flatFld_ist);
plicFuncs::makeFlatFld(T1, istMap, T1_flatFld_ist);
plicFuncs::makeFlatFld(T0, istMap, T0_flatFld_ist);

for(i=0; i<n; i++)
{
    const volScalarField& Y1i = Y1[i];
    const volScalarField& Y0i = Y0[i];
    List<scalar>& Y1i_flatFld_ist = Y1_flatFld_ist[i];
    List<scalar>& Y0i_flatFld_ist = Y0_flatFld_ist[i];
    plicFuncs::makeFlatFld(Y1i, istMap, Y1i_flatFld_ist);
    plicFuncs::makeFlatFld(Y0i, istMap, Y0i_flatFld_ist);
}
//Done make flat fields for alpha1 and phase centroids required for gradient calculation
//-----------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------//
//Calculate interfacial species mass fluxes

plicFuncs::calc_Ts
(
    mesh,
    istCellStencil,
    Y1_flatFld_ist,
    Y0_flatFld_ist,
    T1_flatFld_ist,
    T0_flatFld_ist,
    alpha1_flatFld_ist,
    alpha1,
    lambda1,
    lambda0,
    C_ph1_flatFld_ist,
    C_ph0_flatFld_ist,
    C_intfc,
    A_intfc,
    nHat,
    nSpecies,
    T1,
    T0,
    ALPHA_2PH_MIN_IST,
    erf_a,
    erf_b,
    nErf,
    erfInv_table,
    useErf,
    ist_debug2,
    Hs1,
    Hs0,
    Ts,
    delta1,
    delta0,
    osIst
);

Info << "Ts max = " << gMax(Ts) << " Ts min = " << gMin(Ts) << endl;
Info << "Hs1 max = " << gMax(Hs1) << " Hs1 min = " << gMin(Hs1) << endl;

// Apply Hs with limiter - T0MAX and T1MIN in constant/transportProperties
forAll(alpha1Cells, cellI)
{
    alpha1_cellI = alpha1Cells[cellI];
    scalar Hs1_limiter, Hs1_max, Hs1_min, h_tmp, H1_tmp, H0_tmp, T_tmp;
        if(alpha1_cellI > ALPHA_2PH_MIN_IST && alpha1_cellI < ALPHA_2PH_MAX_IST)
        {
            Hs1_limiter = 1;
            if(Hs1Cells[cellI] > 0)
            {
                // Receiver limitation
                T_tmp = min(T1max, T0Cells[cellI]);
                for(i=0; i<n; i++)
                {
                    x_cellI_tmp[i] = x1[i].internalField()[cellI];
                }
                plicFuncs::calc_h_from_T(x_cellI_tmp, P_thermo_tmp, T_tmp, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, h_tmp);
                H1_tmp = h_tmp*rho1Cells[cellI]*alpha1Cells[cellI];
                Hs1_max = H1_tmp - H1Cells[cellI];
                Hs1_limiter = min(Hs1_limiter, Hs1_max/Hs1Cells[cellI]/dt);
                // Donorer limitation
                T_tmp = max(T1Cells[cellI], T0min);
                for(i=0; i<n; i++)
                {
                    x_cellI_tmp[i] = x0[i].internalField()[cellI];
                }
                plicFuncs::calc_h_from_T(x_cellI_tmp, P_thermo_tmp, T_tmp, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, h_tmp);
                H0_tmp = h_tmp*rho0Cells[cellI]*(1 - alpha1Cells[cellI]);
                Hs1_max = H0Cells[cellI] - H0_tmp;
                Hs1_limiter = min(Hs1_limiter, Hs1_max/Hs1Cells[cellI]/dt);
                // Corrected value
                Hs1Cells[cellI] = Hs1Cells[cellI]*Hs1_limiter;
            } 
            else if(Hs1Cells[cellI] < 0)
            {
                // Receiver limitation
                T_tmp = min(T0max, T1Cells[cellI]);
                for(i=0; i<n; i++)
                {
                    x_cellI_tmp[i] = x0[i].internalField()[cellI];
                }
                plicFuncs::calc_h_from_T(x_cellI_tmp, P_thermo_tmp, T_tmp, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, h_tmp);
                H0_tmp = h_tmp*rho0Cells[cellI]*(1 - alpha1Cells[cellI]);
                Hs1_min = -(H0_tmp - H0Cells[cellI]);
                Hs1_limiter = min(Hs1_limiter, Hs1_min/Hs1Cells[cellI]/dt);
                // Donorer limitation
                T_tmp = max(T0Cells[cellI], T1min);
                for(i=0; i<n; i++)
                {
                    x_cellI_tmp[i] = x1[i].internalField()[cellI];
                }
                plicFuncs::calc_h_from_T(x_cellI_tmp, P_thermo_tmp, T_tmp, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, h_tmp);
                H1_tmp = h_tmp*rho1Cells[cellI]*alpha1Cells[cellI];
                Hs1_min = -(H1Cells[cellI] - H1_tmp);
                Hs1_limiter = min(Hs1_limiter, Hs1_min/Hs1Cells[cellI]/dt);
                // Corrected value
                Hs1Cells[cellI] = Hs1Cells[cellI]*Hs1_limiter;
            } // if(Hs1Cells[cellI] > 0)

            Info << "Hs-limiter = " << Hs1_limiter << endl;

            // Apply H-source terms
            H1Cells[cellI] += Hs1Cells[cellI]*dt;
            H0Cells[cellI] -= Hs1Cells[cellI]*dt;
            // Update temperature, h, rhoh fields
            h_tmp = H1Cells[cellI]/rho1Cells[cellI]/alpha1Cells[cellI];
            h1Cells[cellI] = h_tmp;
            rhoh1Cells[cellI] = h_tmp*rho1Cells[cellI];
                for(i=0; i<n; i++)
                {
                    x_cellI_tmp[i] = x1[i].internalField()[cellI];
                }
        plicFuncs::calc_T_from_h(x_cellI_tmp, P_thermo_tmp, h_tmp, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, T_tmp, MAX_ITER_T, T_TOL, H_TOL, T1min, T1max);
        T1Cells[cellI] = T_tmp;

        h_tmp = H0Cells[cellI]/rho0Cells[cellI]/(1 - alpha1Cells[cellI]);
            h0Cells[cellI] = h_tmp;
            rhoh0Cells[cellI] = h_tmp*rho0Cells[cellI];
                for(i=0; i<n; i++)
                {
                    x_cellI_tmp[i] = x0[i].internalField()[cellI];
                }
        plicFuncs::calc_T_from_h(x_cellI_tmp, P_thermo_tmp, h_tmp, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, T_tmp, MAX_ITER_T, T_TOL, H_TOL, T0min, T0max);
        T0Cells[cellI] = T_tmp;
        } // if(alpha1_cellI > ALPHA_2PH_MIN && alpha1_cellI < ALPHA_2PH_MAX)
} // forAll(alpha1Cells, cellI)
// Done Applying Hs with limiter - T0MAX and T1MIN in constant/transportProperties


//plicFuncs::calc_Js_T


//plicFuncs::calc_Ts


