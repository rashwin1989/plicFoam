{
  int i;

  // input data
  printf("\n--------------------------------------------------\n");
  printf(  "-----Input data: input.dat+-----------------------\n");
  printf("n_species: %d\n", n);
  printf("species:   ");
  for (i=0;i<n;i++) printf("%d ", idx_species[i]); printf("\n");
  printf("groups:   ");
  for (i=0;i<n;i++) printf("%d ", idx_groups[i]);  printf("\n");

  printf("X1:       ");
  for (i=0;i<n;i++) printf("%6.4lf ", x1[i]);  printf("\n");
  printf("X2:       ");
  for (i=0;i<n;i++) printf("%6.4lf ", x2[i]);  printf("\n");

  printf("n_grid: %d\n", NR);
  printf("R:      %lf\n", R);
  printf("rs0:    %lf\n", r_s);
  printf("T0:     %lf\n", T0);
  printf("T1:     %lf\n", T1);
  printf("P:      %le\n", P );
  printf("restart: %d\n", restart);
  printf("n_freq:  %d\n", nfreq);
  printf("Courant:  %lf\n", Courant);
  printf("bConvect:  %d\n", bConvection);
  printf("nFluxType:  %d\n", n_flux_type);
  printf("nExpectInterfaces:  %d\n", nExpectInterfaces);
  printf("nFullyMiscible:  %d\n", nFullyMiscible);
  printf("bConstantDt:  %d\n", bConstantDt);
  printf("GSLmaxLoop:   %d\n", n_GSL_max_loop);
  printf("Dirichlet:    %d\n", nDirichlet_U1);
  printf( "---------------------------------------------------\n");

  // species data
  printf("\n--------------------------------------------------\n");
  printf(  "-----Species' properties--------------------------\n");
  for (i=0; i<n; i++){
    printf("%5.1f,%6.1f,%4.1f,%6.1f,%6.1f,%4.1f,%5.1f,%5.3f,%4.2f,%6.1f\n",
        Pc[i]*1e-5,Tc[i],w[i],MW[i],Tb[i],SG[i],H8[i],k[i],dm[i],Vc[i]);
  }
  printf( "---------------------------------------------------\n");
}
