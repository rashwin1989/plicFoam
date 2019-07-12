if(!isothermal || variableTransportProp)
{
plicFuncs::correct_thermo_trans_prop
(
    mesh,
    P_thermo_tmp,    
    T1,
    x1,
    Y1,
    n,
    Pc,
    Tc,
    Vc,
    w,
    MW,    
    Tb,
    SG,
    H8,
    k,
    dm,
    kij_T,
    Ta_kij,
    Tb_kij,
    nT_kij,
    kij,
    v1,   
    mu1_thermo,
    lambda1_thermo,
    Cp1_thermo,
    rho1_thermo,
    Dij1,
    transProp_debug,
    osTransProp
);

plicFuncs::correct_thermo_trans_prop
(
    mesh,
    P_thermo_tmp,    
    T0,
    x0,
    Y0,
    n,
    Pc,
    Tc,
    Vc,
    w,
    MW,    
    Tb,
    SG,
    H8,
    k,
    dm,
    kij_T,
    Ta_kij,
    Tb_kij,
    nT_kij,
    kij,
    v0,   
    mu0_thermo,
    lambda0_thermo,
    Cp0_thermo,
    rho0_thermo,
    Dij0,
    transProp_debug,
    osTransProp
);
}    
if(variableTransportProp)
{
    plicFuncs::correct_D_from_Dij(n, x1, Dij1, D1);
    plicFuncs::correct_D_from_Dij(n, x0, Dij0, D0);
    // diffusivityFactor
    for(i=0; i<n; i++)
    {
        D1[i] *= diffusivityFactor;
        D0[i] *= diffusivityFactor;
    }
}

Info<< "D1[0]:  min = " << gMin(D1[0]) << "  max = " << gMax(D1[0]) << endl;
Info<< "D0[0]:  min = " << gMin(D0[0]) << "  max = " << gMax(D0[0]) << endl;

if(variableTransportProp)
{
Info<< "mu1_thermo:  min = " << gMin(mu1_thermo) << "  max = " << gMax(mu1_thermo) << endl;
Info<< "mu0_thermo:  min = " << gMin(mu0_thermo) << "  max = " << gMax(mu0_thermo) << endl;

mu1 = mu1_thermo;
mu0 = mu0_thermo;

Info<< "lambda1_thermo:  min = " << gMin(lambda1_thermo) << "  max = " << gMax(lambda1_thermo) << endl;
Info<< "lambda0_thermo:  min = " << gMin(lambda0_thermo) << "  max = " << gMax(lambda0_thermo) << endl;

lambda1 = lambda1_thermo;
lambda0 = lambda0_thermo;
// diffusivityFactor
lambda1 *= diffusivityFactor;
lambda0 *= diffusivityFactor;
}

if(!isothermal || variableTransportProp)
{
Info<< "Cp1_thermo:  min = " << gMin(Cp1_thermo) << "  max = " << gMax(Cp1_thermo) << endl;
Info<< "Cp0_thermo:  min = " << gMin(Cp0_thermo) << "  max = " << gMax(Cp0_thermo) << endl;

Cp1 = Cp1_thermo; // Required to be updates in non-isothermal cases to set diffCourant number for deltat
Cp0 = Cp0_thermo; // Required to be updates in non-isothermal cases to set diffCourant number for deltat 
}

Info << "Diffusivity factor applied = " << diffusivityFactor << endl;

Dh1 = lambda1/rho1/Cp1;
Dh0 = lambda0/rho0/Cp0;

Le1 = Dh1/D1[0];
Le0 = Dh0/D0[0];

Info<< "Done correcting thermo and transport props" << endl;
