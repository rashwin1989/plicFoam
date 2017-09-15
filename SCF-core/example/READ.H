{
  // I. READ INPUT FILE
  int i, j, n_species_db, iread;
  char c[50];
  FILE *f;

  int n_pseudospecies_db, n_purespecies_db;
  f = fopen("species.dat+", "r");
  iread = fscanf(f, "%d %d", &n_pseudospecies_db, &n_purespecies_db);
  n_species_db = n_pseudospecies_db + n_purespecies_db;
  fclose(f);

  f = fopen("input.dat+", "r");

  iread=fscanf(f, "%s %d", c, &n);  // num of species

  _NEW_(idx_species, int, n);
  _NEW_(idx_groups, int, n);
  _NEW_(x1, double, n);
  _NEW_(x2, double, n);
  _NEW_(xi1, double, n);
  _NEW_(xi2, double, n);

  iread=fscanf(f, "%s", c);         // index of species
  for (i=0; i<n; i++) iread=fscanf(f, "%d", &idx_species[i]);

  iread=fscanf(f, "%s", c);         // index of PPR78 groups
  for (i=0; i<n; i++) iread=fscanf(f, "%d", &idx_groups[i]);

  iread=fscanf(f, "%s", c);         // left hand side phase molar fractions
  for (i=0; i<n; i++) iread=fscanf(f, "%lf", &x1[i]);
  iread=fscanf(f, "%s", c);         // right hand side phase molar fractions
  for (i=0; i<n; i++) iread=fscanf(f, "%lf", &x2[i]);

  iread=fscanf(f, "%s %d", c, &NR);  // num of grid points
  iread=fscanf(f, "%s %lf",c, &R);   // radial distance

  iread=fscanf(f, "%s %lf",c, &r_s); // initial droplet radius
  iread=fscanf(f, "%s %lf",c, &r0_buffer); // initial droplet buffer distance for T

  iread=fscanf(f, "%s %lf",c, &T0);  // initial oil temperature
  iread=fscanf(f, "%s %lf",c, &T1);  // initial water temperature
  iread=fscanf(f, "%s %lf",c, &P );  // uniform pressure (Pa)

  iread=fscanf(f, "%s %lf",c, &dt);     // delta time (sec)
  iread=fscanf(f, "%s %lf",c, &t_max);  // maximum time (sec)

  iread=fscanf(f, "%s %d", c, &restart); // restart tag
  iread=fscanf(f, "%s %d", c, &nfreq);   // output frequency
  iread=fscanf(f, "%s %lf", c, &Courant);   // output frequency
  iread=fscanf(f, "%s %d", c, &bConvection);   // Convection YES or NO
  iread=fscanf(f, "%s %d", c, &n_flux_type);   // 0: Fick's law; >0: Maxwell-Stefan

  iread=fscanf(f, "%s %d", c, &nExpectInterfaces);   // expected interface number

  iread=fscanf(f, "%s %d", c, &nFullyMiscible);   
  iread=fscanf(f, "%s %d", c, &bConstantDt);   
  
  iread=fscanf(f, "%s %d", c, &n_GSL_max_loop);   
  iread=fscanf(f, "%s %d", c, &nDirichlet_U1);   

  if (nfreq<=0) nfreq=1000; // default value
  if (Courant>1) Courant=1e-2; // default value
  //if (n_flux_type>0) n_flux_type=1; else n_flux_type=0;
  if (n_GSL_max_loop<=0) n_GSL_max_loop = 100;

  iread=fscanf(f, "%s", c);         // left hand side phase molar fractions
  for (i=0; i<n; i++) iread=fscanf(f, "%lf", &xi1[i]);
  iread=fscanf(f, "%s", c);         // right hand side phase molar fractions
  for (i=0; i<n; i++) iread=fscanf(f, "%lf", &xi2[i]);

  fclose(f);

  // II. ALLOCATE MEMORY
  // 1. species
  _NEW_(tk, int, n);
  _NEW_(Pc, double, n);
  _NEW_(Tc, double, n);
  _NEW_(Vc, double, n);
  _NEW_(w,  double, n);
  _NEW_(MW, double, n);
  _NEW_(dm, double, n);
  _NEW_(k,  double, n);
  _NEW_(H8, double, n);
  _NEW_(Tb, double, n);
  _NEW_(SG, double, n);
  // 2. grid
  _NEW_(r, double, NR);
  _NEW_(phase_type, int, NR);
  _NEW2_(x,double, NR, n); _NEW2_(x_n, double, NR, n);
  _NEW2_(y,double, NR, n); _NEW2_(y_n, double, NR, n);
  _NEW2_(rhoY,double, NR, n); _NEW2_(rhoY_n, double, NR, n);
  _NEW2_(D, double, NR, n);   _NEW2_(Dij, double, NR, n*n);  
  _NEW2_(dlnphi_dxj, double, NR, n*n); 
  _NEW2_(dxj_dr, double, NR, n); 
  _NEW2_(dxj_dr_,double, NR, n); 
  _NEW_ (dlnphi_dxj_tmp, double, n*n); 
  _NEW_ (dlnphidT_tmp, double, n); 
  _NEW2_(flux_m, double, NR, n);   
  _NEW_ (flux_h, double, NR);   
  _NEW2_(flux_m_,double, NR, n);   
  _NEW2_(lnphi, double, NR, n);   
  _NEW2_(phi, double, NR, n);   
  _NEW2_(dlnphidT, double, NR, n);   
  _NEW2_(xphi, double, NR, n);   
  _NEW2_(advect,double, NR, n);   
  _NEW2_(rhs_flux,double, NR, n); 
  _NEW2_(Cp_IG,double,NR, n);
  _NEW2_(Hpar,double, NR, n);
  _NEW2_(Hdep,double, NR, n);
  _NEW_(rho, double, NR);  _NEW_(rho_n, double, NR);_NEW_(rho_s, double, NR);
  _NEW_(u,   double, NR);  _NEW_(u_n,   double, NR);_NEW_(u_f,   double, NR);
  _NEW_(T,   double, NR);  _NEW_(T_n,   double, NR);
  _NEW_(cond,double, NR);
  _NEW_(p,   double, NR);
  _NEW_(am,  double, NR);
  _NEW_(bm,  double, NR);
  _NEW_ (drhodT,double, NR);
  _NEW2_(drhodY,double, NR, n);
  _NEW_(V,   double, NR);  _NEW_(MWave, double, NR);
  _NEW_(Cp,  double, NR);
  _NEW_(Cv,  double, NR);
  _NEW_(Cv_IG,double, NR);
  _NEW_(vis,double, NR);
  _NEW2_(S_IG,double, NR, n);
  _NEW2_(H_IG,double, NR, n);
  _NEW2_(G_IG,double, NR, n);
  _NEW2_(dG_IG,double,NR, n);
  _NEW2_(G_par,double, NR, n);
  _NEW2_(n1, double, NR, n);
  _NEW_(n1s, double, n);
  _NEW_(n2s, double, n);
  _NEW_(n3s, double, n);
  _NEW_(n4s, double, n);
  _NEW_(H1,  double, NR);
  // 3. interface
  _NEW_(x_1, double, n); _NEW_(y_1, double, n); _NEW_(rhoY_1, double, n);
  _NEW_(x_2, double, n); _NEW_(y_2, double, n); _NEW_(rhoY_2, double, n);
  _NEW_(x_a, double, n); _NEW_(y_a, double, n);
  _NEW_(x_1_n, double, n);
  _NEW_(x_2_n, double, n);
  _NEW_(flux_m_1, double, n);  _NEW_(dydr_1, double, n);
  _NEW_(flux_m_2, double, n);  _NEW_(dydr_2, double, n);
  _NEW_(flux_m_A, double, n);
  _NEW_(flux_m_B, double, n);
  //_NEW_(xi1, double, n);
  //_NEW_(xi2, double, n);
  _NEW_(h_1, double, n); _NEW_(D_1, double, n); _NEW_(Dij_1, double, n*n);
  _NEW_(h_2, double, n); _NEW_(D_2, double, n); _NEW_(Dij_2, double, n*n);
  _NEW_(Hdep_1,double, n);
  _NEW_(Hdep_2,double, n);
  _NEW_(Cp_IG_tmp, double, n);
  _NEW_(lnphi_1, double, n); 
  _NEW_(lnphi_2, double, n);
  _NEW_(phi_1, double, n); 
  _NEW_(phi_2, double, n);
  _NEW_(xphi_1, double, n); 
  _NEW_(xphi_2, double, n);
  _NEW_(K, double, n); 
  _NEW_(RM_ni, double, n); 
  _NEW_(M_ai_yi, double, n); 
  _NEW_(x_1_max, double, n); 
  _NEW_(x_1_min, double, n); 
  _NEW_(x_1_dx,  double, n); 
  _NEW_(dlnphi_dxj_1, double, n*n); 
  _NEW_(dlnphi_dxj_2, double, n*n);
  _NEW_(dfdx, double, n);
  // 4. initial
  //_NEW_(x1, double, n);
  //_NEW_(x2, double, n);
  // 5. unused; but required 
  _NEW_(Vpar,  double, n);
  _NEW_(lnphi_tmp, double, n);
  _NEW_(coef_ab, double, n);
  for (j=0;j<n;j++) coef_ab[j]=-1;
  _NEW_(Dij_tmp, double, n*n);
  _NEW_(kij_tmp, double, n*n);
  _NEW_(x_tmp, double, n);
  _NEW_(dx_tmp, double, n);
  _NEW_(Vpar_tmp, double, n); 
  // 6. chemical potentials
  _NEW2_(mu_IG, double, NR, n); _NEW2_(dMu, double, NR, n); 
  _NEW_(mu_IG_1, double, n);    _NEW_(dMu_1, double, n);
  _NEW_(mu_IG_2, double, n);    _NEW_(dMu_2, double, n);

  _NEW2_(Dij_mc,  double, NR, n*n);  
  _NEW_ (Dij_mc_1,double, n*n);
  _NEW_ (Dij_mc_2,double, n*n);

  _NEW2_(Dij_Fick,  double, NR, n*n);  
  //_NEW2_(Dij_new,  double, NR, n*n);  

  // 7. upwind HJ-WENO5 data
  /*
  i_lo_gb = 0      ; // fixed
  i_hi_gb = NR + 6; // fixed
  i_lo_fb_L = 4      ; // fixed
  i_hi_fb_R = NR + 2; // fixed
  _NEW2_(x_L,    double, n, NR+7); _NEW2_(x_R,    double, n, NR+7);
  _NEW2_(lnphi_L,double, n, NR+7); _NEW2_(lnphi_R,double, n, NR+7);
  _NEW2_(rhoyu_L,double, n, NR+7); _NEW2_(rhoyu_R,double, n, NR+7);
  _NEW_ (T_L,    double, NR+7);    _NEW_ (T_R,    double, NR+7);  
  _NEW2_(flux_L, double, n, NR+7); _NEW2_(flux_R, double, n, NR+7);
  _NEW2_(dxdr_L,    double, n, NR+7); _NEW2_(dxdr_R,    double, n, NR+7);
  _NEW2_(dlnphidr_L,double, n, NR+7); _NEW2_(dlnphidr_R,double, n, NR+7);
  _NEW2_(drhoyudr_L,double, n, NR+7); _NEW2_(drhoyudr_R,double, n, NR+7);
  _NEW_ (dTdr_L,    double, NR+7);    _NEW_ (dTdr_R,    double, NR+7);  
  _NEW_ (D1,    double, NR+7);

  // zero all *_L and *_R
  for (i=0;i<NR+7;i++) {
    for (j=0;j<n;j++) {
      x_L[j][i]=0; lnphi_L[j][i]=0; flux_L[j][i]=0; rhoyu_L[j][i]=0;
      x_R[j][i]=0; lnphi_R[j][i]=0; flux_R[j][i]=0; rhoyu_R[j][i]=0;
    }
    T_L[i] = 0;
    T_R[i] = 0;
  }
  */

  // 8. Runga-Kutta data
  _NEW2_(RHS_Y_RK1,double,NR,n);_NEW_(RHS_T_RK1,double,NR);
  _NEW2_(RHS_Y_RK2,double,NR,n);_NEW_(RHS_T_RK2,double,NR);
  _NEW2_(RHS_Y_RK3,double,NR,n);_NEW_(RHS_T_RK3,double,NR);
  _NEW2_(RHS_Y_RK4,double,NR,n);_NEW_(RHS_T_RK4,double,NR);

  RHS_Y[0] = RHS_Y_RK1; RHS_T[0] = RHS_T_RK1;
  RHS_Y[1] = RHS_Y_RK2; RHS_T[1] = RHS_T_RK2;
  RHS_Y[2] = RHS_Y_RK3; RHS_T[2] = RHS_T_RK3;
  RHS_Y[3] = RHS_Y_RK4; RHS_T[3] = RHS_T_RK4;

  // III. READ SPECIES

  f = fopen("petro.dat+","r");

  for (j=0;j<n;j++) {
    for (i=0;i<n_species_db;i++) {
      double tmp;
      double xPNA[3];

      //[x' y' Tb' Tc' Pc' w' MW' Vc' SG' H_8' xP' xN' xA'];
      iread=fscanf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
          &tmp, &tmp, &Tb[j], &Tc[j], &Pc[j], &w[j], &MW[j], &Vc[j], 
          &SG[j], &H8[j], &xPNA[0], &xPNA[1], &xPNA[2]);

      //tmp = 5.811+4.919*w[j];
      //Vc[j] = 83.14*Tc[j]/Pc[j]/(3.72+.26*(tmp-7)); // cm^3/mol

      Pc[j] *= 1e5;  // bar --> Pa

      if (idx_groups[j]>=0) {
      // the first condition
      // we use PPR78 to calculate BIPs
        tk[j] = -(n_pseudospecies_db*idx_groups[j]+i+1);
      }
      else {
      // in the following
      // Guang's BIPs are used
      // NOTE that there are only 3 species
      // water, toluene, n-decane
        tk[j] = -idx_groups[j];
      }

      //x1[j] *= xPNA[ idx_groups[j] ];

      dm[j] = 0.0; 
      k [j] = 0.0;

      if (i==idx_species[j]) {
        // SPECIAL DATA POINTS BUILT IN THE CODE
        // 1. water
        if (fabs(Tc[j]-647)<0.5 && fabs(Pc[j]-220.6e5)<0.5e5) {
          Vc[j] = 57.1; // cm^3/mol
          dm[j] = 1.85; // Debye
          k [j] = 0.076;

        }
        
        if (H8[j]<0)
        // 2. toluene
        if (fabs(Tc[j]-591.8)<0.5 && fabs(Pc[j]-41.1e5)<0.5e5) {
          Vc[j] = 316.0;
          dm[j] = 0.36; // Debye

        // 3. n-decane
        }else if (fabs(Tc[j]-617.7)<0.5 && fabs(Pc[j]-21.1e5)<0.5e5) {
          Vc[j] = 600.0;

        // 4. n-C30
        }else if (fabs(Tc[j]-844.0)<0.5 && fabs(Pc[j]- 8.0e5)<0.5e5) {
          Vc[j] = 1805.0;

        // 5. n-C50
        }else if (fabs(Tc[j]-1073.6)<0.5&& fabs(Pc[j]-3.51e5)<0.5e5) {
          Vc[j] = 2999;

        // 6. Benzene-C10
        }else if (fabs(Tc[j]-753.0)<0.5 && fabs(Pc[j]-17.7e5)<0.5e5) {
          Vc[j] = 813;

        // 7. Benzene-C30
        }else if (fabs(Tc[j]-965.6)<0.5 && fabs(Pc[j]- 5.3e5)<0.5e5) {
          Vc[j] = 1943.5;

        // 8. Naphthalene-C10
        }else if (fabs(Tc[j]-859.0)<0.5 && fabs(Pc[j]-15.8e5)<0.5e5) {
          Vc[j] = 1070;

        // 9. Naphthalene-C12
        }else if (fabs(Tc[j]-854.6)<0.5 && fabs(Pc[j]-13.0e5)<0.5e5) {
          Vc[j] = 1081.5;

        // 10. Benzene
        }else if (fabs(Tc[j]-562.0)<0.5 && fabs(Pc[j]-49.0e5)<0.5e5) {
          Vc[j] = 256.0;

        // 11. o-xylene
        }else if (fabs(Tc[j]-630.3)<0.5 && fabs(Pc[j]-37.3e5)<0.5e5) {
          Vc[j] = 370.0;

        // 12. p-xylene
        }else if (fabs(Tc[j]-616.2)<0.5 && fabs(Pc[j]-35.1e5)<0.5e5) {
          Vc[j] = 378.0;

        // 13. 1,3,5-trimethylbenzene
        }else if (fabs(Tc[j]-637.3)<0.5 && fabs(Pc[j]-31.3e5)<0.5e5) {
          Vc[j] = 433.0;

        // 14. naphthalene
        }else if (fabs(Tc[j]-748.4)<0.5 && fabs(Pc[j]-40.5e5)<0.5e5) {
          Vc[j] = 407.0;

        // 15. 1-methylnaphthalene
        }else if (fabs(Tc[j]-772.0)<0.5 && fabs(Pc[j]-36.0e5)<0.5e5) {
          Vc[j] = 465;

        // 16. anthracene
        }else if (fabs(Tc[j]-873.0)<0.5 && fabs(Pc[j]-29.0e5)<0.5e5) {
          Vc[j] = 554;

        // 17. 1,2-diphenylethane
        }else if (fabs(Tc[j]-780.0)<0.5 && fabs(Pc[j]-26.5e5)<0.5e5) {
          Vc[j] = 616;

        // 18. pyrene
        }else if (fabs(Tc[j]-936.0)<0.5 && fabs(Pc[j]-26.1e5)<0.5e5) {
          Vc[j] = 660;

        // 19. n-C16
        }else if (fabs(Tc[j]-723.0)<0.5 && fabs(Pc[j]-14.0e5)<0.5e5) {
          Vc[j] = 969.2;

        // 20. trans-decalin
        }else if (fabs(Tc[j]-687.0)<0.5 && fabs(Pc[j]-32.0e5)<0.5e5) {
          Vc[j] = 480.0;

        // 21. butylbenzene
        }else if (fabs(Tc[j]-660.5)<0.5 && fabs(Pc[j]-28.9e5)<0.5e5) {
          Vc[j] = 497.0;

        // 22. hexylbenzene
        }else if (fabs(Tc[j]-698.0)<0.5 && fabs(Pc[j]-23.8e5)<0.5e5) {
          Vc[j] = 593.0;
        }

        // WHEN EQUALS, WE FOUND THE DATA 
        // SO QUIT THE READING FOR SPECIES I
        break;
      }
    }
    rewind(f);
  }
  fclose(f);

  // init. pressure
  for (i=0; i<NR; i++) p[i] = P;

  // initialize umfpack
  initialUmfpack(&flux_umf);
  flux_memoryMatVec(&flux_umf, n);
  flux_compRow_setRowCol(&flux_umf, n);
}
