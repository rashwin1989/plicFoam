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
    kij,
    v1,   
    mu1_thermo,
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
    kij,
    v0,   
    mu0_thermo,
    Dij0,
    transProp_debug,
    osTransProp
);

//plicFuncs::correct_D_from_Dij(n, x1, Dij1, D1);
//plicFuncs::correct_D_from_Dij(n, x0, Dij0, D0);

for(i=0; i<n; i++)
{
    for(j=i+1; j<n; j++)
    {
        ijIndex = i*n + j;
        Info<< "D1[" << i << "," << j << "]:  min = " << gMin(Dij1[ijIndex]) << "  max = " << gMax(Dij1[ijIndex]) << endl;
        Info<< "D0[" << i << "," << j << "]:  min = " << gMin(Dij0[ijIndex]) << "  max = " << gMax(Dij0[ijIndex]) << endl;
    }
}

Info<< "mu1_thermo:  min = " << gMin(mu1_thermo) << "  max = " << gMax(mu1_thermo) << endl;
Info<< "mu0_thermo:  min = " << gMin(mu0_thermo) << "  max = " << gMax(mu0_thermo) << endl;

Info<< "Done correcting thermo and transport props" << endl;
