plicFuncs::correct_thermo_trans_prop
(
    mesh,    
    n,
    Pc,
    Tc,
    Vc,
    w,
    MW,
    tk,
    coef_ab,
    Tb,
    SG,
    H8,
    k,
    dm,
    P_thermo_tmp,
    x1,
    T1,
    v1,    
    mu1,
    Dij1
);

plicFuncs::correct_thermo_trans_prop
(
    mesh,    
    n,
    Pc,
    Tc,
    Vc,
    w,
    MW,
    tk,
    coef_ab,
    Tb,
    SG,
    H8,
    k,
    dm,
    P_thermo_tmp,
    x0,
    T0,
    v0,    
    mu0,
    Dij0
);

plicFuncs::correct_D_from_Dij(x1, Dij1, D1, n);
plicFuncs::correct_D_from_Dij(x0, Dij0, D0, n);

Info<< "Done correcting thermo and transport props" << endl;
