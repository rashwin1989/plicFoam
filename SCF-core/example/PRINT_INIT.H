{
  int i;
  FILE *f = fopen("debug_init.dat", "w");
  for (i=0; i<NR; i++) {
    fprintf(f, "%le %le %le %le %le %le %le %le %le %le %le %le %le\n", 
        r[i],T[i],V[i],rho[i],Cp[i],cond[i],x[i][0],x[i][n-1],Hpar[i][0],Hpar[i][n-1],D[i][0],D[i][10],D[i][n-1]);
    //  1    2    3    4      5     6       7       8         9          10           11      12       13
  }
  fclose(f);
}
