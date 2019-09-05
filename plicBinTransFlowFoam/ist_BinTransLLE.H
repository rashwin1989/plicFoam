// Ist file for 2-component, non-isothermal interfacial mass transfer

if(ist_debug)
{
    plicFuncs::print_line(osIst, 100);
    osIst<< "Solution Time = " << runTime.timeName() << "  Interfacial transfer" << endl;
    plicFuncs::print_line(osIst, 100);
    osIst<< endl;
}

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


dt = deltaT;

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
    Hs1,
    Hs0,
    Ts,
    delta1,
    delta0,
    ist_debug2,
    osIst
);

Info << "Ts max = " << gMax(Ts) << " Ts min = " << gMin(Ts) << endl;
Info << "Hs1 max = " << gMax(Hs1) << " Hs1 min = " << gMin(Hs1) << endl;

// Hs1 limiter
forAll(alpha1Cells,cellI)
{
    alpha1_cellI = alpha1Cells[cellI];
    if(alpha1_cellI > ALPHA_2PH_MIN_IST && alpha1_cellI < ALPHA_2PH_MAX_IST)
    {
    V_cellI = mesh.V()[cellI];
    alpha1_cellI = alpha1Cells[cellI];
    rho1_cellI = rho1Cells[cellI];
    rho0_cellI = rho0Cells[cellI];
    H1_cellI = H1Cells[cellI];
    H0_cellI = H0Cells[cellI];
    T1_cellI = T1Cells[cellI];
    T0_cellI  = T0Cells[cellI];
    //    n = nSpecies;
    for(i=0;i<n;i++)
    {
        x1_cellI[i] = x1[i].internalField()[cellI];
        x0_cellI[i] = x0[i].internalField()[cellI];
    }
    scalar limiter = 1;
    Hs1_cellI = Hs1Cells[cellI];
    plicFuncs::calc_intfcDiffFlux_limiter_T
        (
            alpha1_cellI,
            V_cellI,
            rho1_cellI,
            x1_cellI,
            T1_cellI,
            H1_cellI,
            rho0_cellI,
            x0_cellI,
            T0_cellI,
            H0_cellI,
            P_thermo_tmp, 
            n, 
            Pc,
            Tc,
            w,
            MW,
            Tb,
            SG,
            H8,
            kij_T,
            Ta_kij,
            Tb_kij,
            nT_kij,
            kij, 
            dt,
            Hs1_cellI,
            limiter,
            MAX_ITER_T,
            T_TOL,
            H_TOL,
            T1min,
            T1max,
            T0min,
            T0max
        );
    //    Info << "My limiter 2 = " << limiter << endl;
    //    Hs1Cells[cellI] = Hs1Cells[cellI]*limiter;
    //}
    //}

    //Info << "T-intfc conduction my limiter done" << endl;

// Apply Hs with limiter - T0MAX and T1MIN in constant/transportProperties
//forAll(alpha1Cells, cellI)
//{
    alpha1_cellI = alpha1Cells[cellI];
    scalar Hs1_limiter, Hs1_max, Hs1_min, h_tmp, H1_tmp, H0_tmp, T_tmp;
    //        if(alpha1_cellI > ALPHA_2PH_MIN_IST && alpha1_cellI < ALPHA_2PH_MAX_IST)
    //        {
                       
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
                Hs1_limiter = max(Hs1_limiter,0);
            } // if(Hs1Cells[cellI] > 0)

                // Corrected value
            if(Hs1_limiter < limiter)
            {
                Info << "Hs1_limiter < limiter. limiter. = " << limiter << " Hs1_limiter = " << Hs1_limiter << " alpha = " << alpha1_cellI << endl;
                limiter = Hs1_limiter;
            }
            Hs1Cells[cellI] = Hs1Cells[cellI]*limiter;

            //            Info << "Hs-limiter = " << Hs1_limiter << endl;
    

            // Apply H-source terms
            H1Cells[cellI] += Hs1Cells[cellI]*dt;
            H0Cells[cellI] -= Hs1Cells[cellI]*dt;

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


dt = deltaT;

#include "correct_thermo_trans_prop.H"
#include "diff_grad_interp.H" // for hpar

// Set Ts_display field for post-processing
forAll(TsCells, cellI)
{
    alpha1_cellI = alpha1Cells[cellI];
    if(alpha1_cellI > 0.01 && alpha1_cellI < 0.99)
    {
        Ts_display.internalField()[cellI] = TsCells[cellI];
    }
    else
    {
        Ts_display.internalField()[cellI] = 0;
    }
}

//Done interfacial conduction
//-----------------------------------------------------------------------------------------------------//


//------------------------------------------------------------------//
//Mass and Accompanying Heat transfer


if(ist_BinTransLLE_mass_Transfer)
{

dt = deltaT;

// Set Ys1 and Ys0 fields
forAll(TsCells, cellI)
{
    alpha1_cellI = alpha1Cells[cellI];
    scalar Ts_cellI = TsCells[cellI];
    if(alpha1_cellI > ALPHA_2PH_MIN_IST && alpha1_cellI < ALPHA_2PH_MAX_IST)
    {
        plicFuncs::calc_kij_from_table(Ts_cellI, n, Ta_kij, Tb_kij, nT_kij, kij_T, kij);
        binarylle_(&P_thermo_tmp,&Ts_cellI,&n,Pc,Tc,w,kij,xs1_tmp,xs0_tmp,&n_miscible);
        plicFuncs::x2y(n,MW,xs1_tmp,ys1_tmp);
        plicFuncs::x2y(n,MW,xs0_tmp,ys0_tmp);
        for(i=0;i<n;i++)
        {
            Ys1[i].internalField()[cellI] = ys1_tmp[i];
            Ys0[i].internalField()[cellI] = ys0_tmp[i];
        }
    }
}

Info << "Done Set Ys1 and Ys0 fields using Ts" << endl;

plicFuncs::calc_Js
(
    mesh,
    istCellStencil,
    Y1_flatFld_ist,
    Y0_flatFld_ist,
    alpha1_flatFld_ist,
    alpha1, 
    rho1,
    rho0,
    D1,
    D0,
    C_ph1_flatFld_ist,
    C_ph0_flatFld_ist,
    C_intfc,
    A_intfc,
    nHat,
    Ys1,
    Ys0,
    Js1,
    Js0,
    grads1,
    grads0,
    delta1,
    delta0,
    nSpecies,
    ALPHA_2PH_MIN_IST,
    erf_a,
    erf_b,
    nErf,
    erfInv_table,
    useErf,
    Y1_0,
    Y0_0,
    ist_debug2,
    osIst
);

Info << "Done calc_Js" << endl;

plicFuncs::calc_mS_Hs_alphaS
(
    mesh,
    C1,
    C0,
    Y1,
    Y0,
    T1,
    T0,
    alpha1,
    rho1,
    rho0,
    Ys1,
    Ys0,
    Js1,
    Js0,
    hpar1,
    hpar0,
    nSpecies,
    ALPHA_2PH_MIN_IST,
    mS1,
    mS0,
    mS1Tot,
    mS0Tot,
    Hs1,
    Hs0,
    alphaS1,
    alphaS0,
    dt,
    ist_debug2,
    osIst
);

Info << "Done calc_mS_Hs_alphaS" << endl;

forAll(pSCells, cellI)
{
    pSCoeff = 1.0/rho1Cells[cellI] - 1.0/rho0Cells[cellI];
    pSCells[cellI] = pSCoeff*mS1TotCells[cellI];
}

forAll(alpha1Cells, cellI)
{
    alpha1_cellI = alpha1Cells[cellI];
    
    if(alpha1_cellI > ALPHA_2PH_MIN_IST && alpha1_cellI < ALPHA_2PH_MAX_IST)
    {        
        calc_ist_for_cell[cellI] = 1;
    }
    else
    {
        calc_ist_for_cell[cellI] = 0;
    }
}

//Info << "Done calculating pSCells" << endl; 

forAll(alpha1Cells, cellI)
{

 if(calc_ist_for_cell[cellI] == 1)
    {
        if(ist_debug)
        {    
            plicFuncs::print_line(osIst, 100);
            osIst<< "Cell " << cellI << endl;
            plicFuncs::print_line(osIst, 100);
            osIst<< "Before interfacial species transfer" << endl;
            plicFuncs::print_line(osIst, 100);
            osIst<< setw(12) << "alpha1" << "  " << setw(12) << "alpha0" << "  " << setw(16) << "alpha1+alpha0" << "  " << setw(12) << "rho1" << "  " << setw(12) << "rho0" << endl;
            plicFuncs::print_line(osIst, 100);
            osIst<< setw(12) << alpha1Cells[cellI] << "  " << setw(12) << alpha0Cells[cellI] << "  " << setw(16) << alpha1Cells[cellI] + alpha0Cells[cellI] << "  " << setw(12) << rho1Cells[cellI] << "  " << setw(12) << rho0Cells[cellI] << endl;
            plicFuncs::print_line(osIst, 100);
            osIst<< setw(7) << "Species" << "  " << setw(12) << "C1" << "  " << setw(12) << "C0" << "  " << setw(12) << "Y1" << "  " << setw(12) << "Y0" << endl;
            plicFuncs::print_line(osIst, 100);
            for(i=0; i<n; i++)
            {
                osIst<< setw(7) << i << "  " << setw(12) << C1[i].internalField()[cellI] << "  " << setw(12) << C0[i].internalField()[cellI] << "  " << setw(12) << Y1[i].internalField()[cellI] << "  " << setw(12) << Y0[i].internalField()[cellI] << endl;
            }
            plicFuncs::print_line(osIst, 100);
            osIst<< "Interface conditions" << endl;
            plicFuncs::print_line(osIst, 100);            
            osIst<< setw(7) << "Species" << "  " << setw(12) << "ys1" << "  " << setw(12) << "ys0" << endl;
            plicFuncs::print_line(osIst, 100);
            for(i=0; i<n; i++)
            {
                osIst<< setw(7) << i << "  " << setw(12) << Ys1[i].internalField()[cellI] << "  " << setw(12) << Ys0[i].internalField()[cellI] << endl;
            }
            plicFuncs::print_line(osIst, 100);
            osIst<< "Interfacial transfer source terms" << endl;
            plicFuncs::print_line(osIst, 100);
            osIst<< setw(12) << "mS1Tot" << "  " << setw(12) << "alphaS1" << "  " << setw(12) << "alphaS0" << "  " << setw(12) << "pS" << "  " << setw(12) << "dt" << endl;
            plicFuncs::print_line(osIst, 100);
            osIst<< setw(12) << mS1TotCells[cellI] << "  " << setw(12) << alphaS1Cells[cellI] << "  " << setw(12) << alphaS0Cells[cellI] << "  " << setw(12) << pSCells[cellI] << "  " << setw(12) << dt << endl;
            plicFuncs::print_line(osIst, 100);
            osIst<< setw(7) << "Species" << "  " << setw(12) << "mS1" << endl;
            plicFuncs::print_line(osIst, 100);
            for(i=0; i<n; i++)
            {
                osIst<< setw(7) << i << "  " << setw(12) << mS1[i].internalField()[cellI] << endl;
            }
            plicFuncs::print_line(osIst, 100);
        }
    }

    //ph-1
    //volume fraction
    alpha1_cellI = alpha1Cells[cellI] + alphaS1Cells[cellI]*dt;
    if(alpha1_cellI < 0)
    {
        osIst<< "Time = " << runTime.timeName() << nl
            << "alpha1 undershoot in cell " << cellI << "  alpha1 = " << alpha1_cellI << endl;
    }
    if(alpha1_cellI > 1)
    {
        osIst<< "Time = " << runTime.timeName() << nl
            << "alpha1 overshoot in cell " << cellI << "  alpha1 = " << alpha1_cellI << endl;
    }    
    alpha1_cellI = max(alpha1_cellI, 0.0);
    alpha1Cells[cellI] = alpha1_cellI;

    //species
    C1Tot_cellI = 0.0;
    for(i=0; i<n; i++)
    {
        C1i_cellI = C1[i].internalField()[cellI] + mS1[i].internalField()[cellI]*dt;
        // if(C1i_cellI < 0)
        // {
        //     osIst<< "Time = " << runTime.timeName() << nl
        //         << "C1[" << i << "] -ve in cell " << cellI << "  C1[" << i << "] = " << C1i_cellI << endl;
        // }
        
        C1i_cellI = max(C1i_cellI, 0.0);
        C1[i].internalField()[cellI] = C1i_cellI;
        C1Tot_cellI += C1i_cellI;
    }

    if(C1Tot_cellI < SMALL)
    {
        C1Tot_cellI += SMALL;
        C1[0].internalField()[cellI] += SMALL;
    }

    Y1t_cellI = 0.0;

    for(i=0; i<(n-1); i++)
    {   
        Y1i_cellI = C1[i].internalField()[cellI]/C1Tot_cellI;
        Y1[i].internalField()[cellI] = Y1i_cellI;
        c1[i].internalField()[cellI] = rho1Cells[cellI]*Y1i_cellI;
        Y1t_cellI += Y1i_cellI;
        y_cellI_tmp[i] = Y1i_cellI;
    }

    Y1n_cellI = max(0.0, 1.0 - Y1t_cellI);
    y_cellI_tmp[n-1] = Y1n_cellI;
    Y1[n-1].internalField()[cellI] = Y1n_cellI;
    c1[n-1].internalField()[cellI] = rho1Cells[cellI]*Y1n_cellI;    
    
    plicFuncs::y2x(n, MW, y_cellI_tmp, x_cellI_tmp);

    for(i=0; i<n; i++)
    {
        x1[i].internalField()[cellI] = x_cellI_tmp[i];
    }

    //    Info << "Phase 1: Done alpha and mS" << endl; 

    // Enthalpy
    if(alpha1_cellI > ALPHA_2PH_MIN_IST)
    {
        H1_cellI = H1Cells[cellI] +  Hs1.internalField()[cellI]*dt;
        h_cellI_tmp = H1_cellI/rho1Cells[cellI]/alpha1_cellI;
        plicFuncs::calc_T_from_h(x_cellI_tmp, P_thermo_tmp, h_cellI_tmp, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, T_cellI_tmp, MAX_ITER_T, T_TOL, H_TOL, T1min, T1max);

    if(T_cellI_tmp < T1min)
    {
        Info << "T_cellI = " << T_cellI_tmp << endl;
        T_cellI_tmp = T1min;
    plicFuncs::calc_h_from_T(x_cellI_tmp, P_thermo_tmp, T_cellI_tmp, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, h_cellI_tmp);
    }
        
    if(T_cellI_tmp > T1max)
    {
        Info << "T_cellI = " << T_cellI_tmp << endl;
        T_cellI_tmp = T1max;        
    plicFuncs::calc_h_from_T(x_cellI_tmp, P_thermo_tmp, T_cellI_tmp, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, h_cellI_tmp);
    }   
    h1Cells[cellI] = h_cellI_tmp;
    rhoh1Cells[cellI] = rho1Cells[cellI]*h_cellI_tmp;
    H1Cells[cellI] = rhoh1Cells[cellI]*alpha1_cellI;
    T1Cells[cellI] = T_cellI_tmp;
    } // if(alpha1_cellI > ALPHA_2PH_MIN_IST)

    //    Info << "Done Phase 1 in ist_BinTransLLE.H" << endl; 

    //ph-0
    //volume fraction
    alpha0_cellI = alpha0Cells[cellI] + alphaS0Cells[cellI]*dt;
    if(alpha0_cellI < 0)
    {
        osIst<< "Time = " << runTime.timeName() << nl
            << "alpha0 undershoot in cell " << cellI << "  alpha0 = " << alpha0_cellI << endl;
    }
    if(alpha0_cellI > 1)
    {
        osIst<< "Time = " << runTime.timeName() << nl
            << "alpha0 overshoot in cell " << cellI << "  alpha0 = " << alpha0_cellI << endl;
    }
    if(alpha0_cellI > 1.5)
    {
        plicFuncs::print_line(osIst, 100);
        osIst<< "Time = " << runTime.timeName() << endl;
        plicFuncs::print_line(osIst, 100);
        plicFuncs::print_line(osIst, 100);
        osIst<< "Cell " << cellI << endl;
        plicFuncs::print_line(osIst, 100);
        osIst<< "Before interfacial species transfer" << endl;
        plicFuncs::print_line(osIst, 100);
        osIst<< setw(14) << "alpha1" << setw(14) << "alpha0" << setw(14) << "alpha1+alpha0" << setw(14) << "rho1" << setw(14) << "rho0" << endl;
        plicFuncs::print_line(osIst, 100);
        osIst<< setw(14) << alpha1Cells[cellI] << setw(14) << alpha0Cells[cellI] << setw(16) << alpha1Cells[cellI] + alpha0Cells[cellI] << setw(14) << rho1Cells[cellI] << setw(14) << rho0Cells[cellI] << endl;
        plicFuncs::print_line(osIst, 100);
        osIst<< setw(9) << "Species" << setw(14) << "C1" << setw(14) << "C0" << setw(14) << "Y1" << setw(14) << "Y0" << endl;
        plicFuncs::print_line(osIst, 100);
        for(i=0; i<n; i++)
        {
            osIst<< setw(9) << i << setw(14) << C1[i].internalField()[cellI] << setw(14) << C0[i].internalField()[cellI] << setw(14) << Y1[i].internalField()[cellI] << setw(14) << Y0[i].internalField()[cellI] << endl;
        }
        plicFuncs::print_line(osIst, 100);
        osIst<< "Interface conditions" << endl;
        plicFuncs::print_line(osIst, 100);            
        osIst<< setw(9) << "Species" << setw(14) << "ys1" << setw(14) << "ys0" << endl;
        plicFuncs::print_line(osIst, 100);
        for(i=0; i<n; i++)
        {
            osIst<< setw(9) << i << setw(14) << Ys1[i].internalField()[cellI] << setw(14) << Ys0[i].internalField()[cellI] << endl;
        }
        plicFuncs::print_line(osIst, 100);
        osIst<< "Interfacial transfer source terms" << endl;
        plicFuncs::print_line(osIst, 100);
        osIst<< setw(14) << "mS1Tot" << setw(14) << "alphaS1" << setw(14) << "alphaS0" << setw(14) << "pS" << setw(14) << "dt" << endl;
        plicFuncs::print_line(osIst, 100);
        osIst<< setw(14) << mS1TotCells[cellI] << setw(14) << alphaS1Cells[cellI] << setw(14) << alphaS0Cells[cellI] << setw(14) << pSCells[cellI] << setw(14) << dt << endl;
        plicFuncs::print_line(osIst, 100);
        osIst<< "D1 = " << D1[0].internalField()[cellI] << "  D0 = " << D1[0].internalField()[cellI] << "  A_intfc = " << A_intfc.internalField()[cellI] << endl;
        plicFuncs::print_line(osIst, 100);
        osIst<< setw(9) << "Species" << setw(14) << "mS1" << setw(14) << "Js1" << setw(14) << "Js0" << endl;
        plicFuncs::print_line(osIst, 100);
        for(i=0; i<n; i++)
        {
            osIst<< setw(9) << i << setw(14) << mS1[i].internalField()[cellI] << setw(14) << Js1[i].internalField()[cellI] << setw(14) << Js0[i].internalField()[cellI] << endl;
        }
        plicFuncs::print_line(osIst, 100);        
    }
    alpha0_cellI = max(alpha0_cellI, 0.0);
    alpha0Cells[cellI] = alpha0_cellI;
    //alpha0Cells[cellI] = scalar(1.0) - alpha1Cells[cellI];

    //species
    C0Tot_cellI = 0.0;
    for(i=0; i<n; i++)
    {
        C0i_cellI = C0[i].internalField()[cellI] + mS0[i].internalField()[cellI]*dt;
        // if(C0i_cellI < 0)
        // {
        //     osIst<< "Time = " << runTime.timeName() << nl
        //         << "C0[" << i << "] -ve in cell " << cellI << "  C0[" << i << "] = " << C0i_cellI << endl;
        // }
        C0i_cellI = max(C0i_cellI, 0.0);
        C0[i].internalField()[cellI] = C0i_cellI;
        C0Tot_cellI += C0i_cellI;
    }

    if(C0Tot_cellI < SMALL)
    {
        C0Tot_cellI += SMALL;
        C0[n-1].internalField()[cellI] += SMALL;
    }

    Y0t_cellI = 0.0;
    for(i=0; i<(n-1); i++)
    {   
        Y0i_cellI = C0[i].internalField()[cellI]/C0Tot_cellI;
        Y0[i].internalField()[cellI] = Y0i_cellI;
        c0[i].internalField()[cellI] = rho0Cells[cellI]*Y0i_cellI;
        Y0t_cellI += Y0i_cellI;
        y_cellI_tmp[i] = Y0i_cellI;
    }

    Y0n_cellI = max(0.0, 1.0 - Y0t_cellI);
    y_cellI_tmp[n-1] = Y0n_cellI;
    Y0[n-1].internalField()[cellI] = Y0n_cellI;
    c0[n-1].internalField()[cellI] = rho0Cells[cellI]*Y0n_cellI;

    plicFuncs::y2x(n, MW, y_cellI_tmp, x_cellI_tmp);

    for(i=0; i<n; i++)
    {
        x0[i].internalField()[cellI] = x_cellI_tmp[i];
    }

    //    Info << "Phase 0: Done alpha and mS" << endl; 

    // Enthalpy
    if(alpha0_cellI > ALPHA_2PH_MIN_IST)
    {
        H0_cellI = H0Cells[cellI] +  Hs0Cells[cellI]*dt;
        h_cellI_tmp = H0_cellI/rho0Cells[cellI]/alpha0_cellI;
        plicFuncs::calc_T_from_h(x_cellI_tmp, P_thermo_tmp, h_cellI_tmp, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, T_cellI_tmp, MAX_ITER_T, T_TOL, H_TOL, T0min, T0max);

    if(T_cellI_tmp < T0min)
    {
        Info << "T_cellI = " << T_cellI_tmp << endl;
        T_cellI_tmp = T0min;
    plicFuncs::calc_h_from_T(x_cellI_tmp, P_thermo_tmp, T_cellI_tmp, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, h_cellI_tmp);
    }
        
    if(T_cellI_tmp > T0max)
    {
        Info << "T_cellI = " << T_cellI_tmp << endl;
        T_cellI_tmp = T0max;        
    plicFuncs::calc_h_from_T(x_cellI_tmp, P_thermo_tmp, T_cellI_tmp, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, h_cellI_tmp);
    }   
    h0Cells[cellI] = h_cellI_tmp;
    rhoh0Cells[cellI] = rho0Cells[cellI]*h_cellI_tmp;
    H0Cells[cellI] = rhoh0Cells[cellI]*alpha0_cellI;    
    T0Cells[cellI] = T_cellI_tmp;
    } // if(alpha1_cellI < ALPHA_2PH_MAX_IST)

    //    Info << "Done Phase 0 in ist_BinTransLLE.H" << endl; 

    if(calc_ist_for_cell[cellI] == 1)
    {
        if(ist_debug)
        {            
            osIst<< "After interfacial species transfer" << endl;
            plicFuncs::print_line(osIst, 100);
            osIst<< setw(12) << "alpha1" << "  " << setw(12) << "alpha0" << "  " << setw(16) << "alpha1+alpha0" << "  " << setw(12) << "rho1" << "  " << setw(12) << "rho0" << endl;
            plicFuncs::print_line(osIst, 100);
            osIst<< setw(12) << alpha1Cells[cellI] << "  " << setw(12) << alpha0Cells[cellI] << "  " << setw(16) << alpha1Cells[cellI] + alpha0Cells[cellI] << "  " << setw(12) << rho1Cells[cellI] << "  " << setw(12) << rho0Cells[cellI] << endl;
            plicFuncs::print_line(osIst, 100);
            osIst<< setw(7) << "Species" << "  " << setw(12) << "C1" << "  " << setw(12) << "C0" << "  " << setw(12) << "Y1" << "  " << setw(12) << "Y0" << endl;
            plicFuncs::print_line(osIst, 100);
            for(i=0; i<n; i++)
            {
                osIst<< setw(7) << i << "  " << setw(12) << C1[i].internalField()[cellI] << "  " << setw(12) << C0[i].internalField()[cellI] << "  " << setw(12) << Y1[i].internalField()[cellI] << "  " << setw(12) << Y0[i].internalField()[cellI] << endl;
            }        
            plicFuncs::print_line(osIst, 100);        
        }
    }
} // forAll(alpha1Cells, cellI)

T1.correctBoundaryConditions();
plicFuncs::correct_boundaryField_h_rhoh_H(P_thermo_tmp, T1, x1, rho1, alpha1, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, h1, rhoh1, H1);
T0.correctBoundaryConditions();
plicFuncs::correct_boundaryField_h_rhoh_H(P_thermo_tmp, T0, x0, rho0, alpha0, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, h0, rhoh0, H0);



Info<< "alpha1 avg = " << alpha1.weightedAverage(mesh.V()).value()
    << "  alpha1 min = " << min(alpha1).value()
    << "  alpha1 undershoot = " << min(min(alpha1).value(), 0)
    << "  alpha1 max = " << max(alpha1).value()
    << "  alpha1 overshoot = " << max((max(alpha1).value() - 1), 0)
    << endl;

Info<< "alpha0 avg = " << alpha0.weightedAverage(mesh.V()).value()
    << "  alpha0 min = " << min(alpha0).value()
    << "  alpha0 undershoot = " << min(min(alpha0).value(), 0)
    << "  alpha0 max = " << max(alpha0).value()
    << "  alpha0 overshoot = " << max((max(alpha0).value() - 1), 0)
    << endl;

alpha1.correctBoundaryConditions();
alpha0.correctBoundaryConditions();

for(i=0; i<n; i++)
{    
    volScalarField& Y1i = Y1[i];
    volScalarField& Y0i = Y0[i];
    volScalarField& c1i = c1[i];
    volScalarField& c0i = c0[i];
    volScalarField& C1i = C1[i];
    volScalarField& C0i = C0[i];

    Y1i.correctBoundaryConditions();
    plicFuncs::correct_boundaryField_C(Y1i, rho1, alpha1, c1i, C1i);

    Y0i.correctBoundaryConditions();
    plicFuncs::correct_boundaryField_C(Y0i, rho0, alpha0, c0i, C0i);
}
for(i=0; i<(n-1); i++)
{
    volScalarField& Y1i = Y1[i];
    volScalarField& Y0i = Y0[i];
    volScalarField& C1i = C1[i];
    volScalarField& C0i = C0[i];

    Info<< "C1[" << i << "] avg = " << C1i.weightedAverage(mesh.V()).value()
        << "  C1[" << i << "] min = " << min(C1i).value()
        << "  C1[" << i << "] max = " << max(C1i).value()
        << nl
        << "Y1[" << i << "] min = " << min(Y1i).value()
        << "  Y1[" << i << "] max = " << max(Y1i).value()
        << endl;

    Info<< "C0[" << i << "] avg = " << C0i.weightedAverage(mesh.V()).value()
        << "  C0[" << i << "] min = " << min(C0i).value()
        << "  C0[" << i << "] max = " << max(C0i).value()
        << nl
        << "Y0[" << i << "] min = " << min(Y0i).value()
        << "  Y0[" << i << "] max = " << max(Y0i).value()
        << endl;
}

//Reset the allowable limits for concentrations and mass fractions in the phases
for(i=0; i<n; i++)
{
    C1max[i] = max(C1[i]).value();
    C0max[i] = max(C0[i]).value();
    C1min[i] = min(C1[i]).value();
    C0min[i] = min(C0[i]).value();
}




Info<< "Done Interfacial Species and Heat Transfer" << endl;

//Done apply interfacial species transfer and alpha source terms, HS source terms, calculate source term for pressure equation
//-----------------------------------------------------------------------------------------------------------//

} // if(ist_BinTransLLE_mass_Transfer)

//dt = deltaT;

// plicFuncs::calc_Ts












//-----------------------------------------------------------------------------------------------------------//
Info<< "Done Interfacial Species Transfer Source Terms Calculation" << endl;


