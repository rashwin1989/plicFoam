{
  int i, j;

  // initialize data -- clean initialization
  // Grid data
  i_lo_ghost = 0;    i_lo = 1; 
  i_hi_ghost = NR-1; i_hi = NR-2; 
  dr = R / (double) (NR-2);

  // Grid locations
  r[0] = -0.5*dr;
  for (i=1; i<NR; i++) r[i] = r[i-1] + dr;

  // physical properties on grid points
  // initial data for two phases only
  // multiple phases are available in RESTART.H
  for (i=0; i<NR; i++) {
    if (r[i]<=r_s) {for (j=0; j<n; j++) x[i][j] = x1[j]; }//T[i] = T0;
    else           {for (j=0; j<n; j++) x[i][j] = x2[j]; }//T[i] = T1;

    //x2y(n, MW, x[i], y[i]);

    u[i] = 0.0;
  }

  // interfacial properties
  u_s = 0;

  for (i=0; i<NR; i++) {
    //myHeaviside(EPSILON*dr, r_s+r0_buffer, T0, T1, 1, &r[i], &T[i]);
    if (r[i]<=r_s) T[i] = T0;
    else T[i]=T1;

    // HACKING for testing
    //T[i] = T0+(T1-T0)/(double)(NR-1)*(double)i;
  }

  t = 0;


}
