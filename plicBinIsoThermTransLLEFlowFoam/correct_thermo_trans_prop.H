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

plicFuncs::correct_D_from_Dij(n, x1, Dij1, D1);
plicFuncs::correct_D_from_Dij(n, x0, Dij0, D0);

Info<< "Done correcting thermo and transport props" << endl;
