plicFuncs::correct_rho(n, Pc, Tc, w, MW, tk, coef_ab, P_thermo_tmp, x1, T1, rho1);
plicFuncs::correct_rho(n, Pc, Tc, w, MW, tk, coef_ab, P_thermo_tmp, x0, T0, rho0);

//Linear interpolate density to faces in both phases
plicFuncs::calc_2ph_Cf(Cf_ph1_own, Cf_ph0_own, Cf_ph1_nei, Cf_ph0_nei, Af_ph1_own, Af_ph0_own, Af_ph1_nei, Af_ph0_nei, Cf_ph1, Cf_ph0);

plicFuncs::calc_2ph_linearInterpolation_weights(mesh, C_ph1, Cf_ph1, w_ph1);
plicFuncs::calc_2ph_linearInterpolation_weights(mesh, C_ph0, Cf_ph0, w_ph0);

plicFuncs::linearInterpolate_2ph(rho1, mesh, w_ph1, rho1f);
plicFuncs::linearInterpolate_2ph(rho0, mesh, w_ph0, rho0f);
