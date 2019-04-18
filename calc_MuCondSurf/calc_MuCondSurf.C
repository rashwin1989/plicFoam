/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    calcRho

Description
    Generate lookup table for density

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "PR_EoS.h"
#include "MACROS2.H"
#include "plic.H"
#include "plicFuncs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    double *Pc, *Tc, *Vc, *w, *MW, *Tb, *SG, *H8, *k, *dm;
    double Ta_kij, Tn_kij, P, T, T_tmp, V, CvIG, mu_tmp, cond_tmp,rho_tmp,Cp_tmp,v_tmp,h_tmp;
    int n, nT_kij, nT, nX, i, j, iT, iX, idx, n_miscible;
    double *kij, *kij_T;
    double Tmin, Tmax, dT, Xmin, Xmax, dX;
    double *x, *x_eqm0, *x_eqm1;
    bool is_2phase;

    int n_species_db, iread;
    char str_tmp[50];
    FILE *f;
    //    FILE *ft;
    //    FILE *fx;
    int n_pseudospecies_db, n_purespecies_db;
    int *idx_species, *idx_groups;

    f = fopen("constant/BIP_nT.dat+", "r");
    iread = fscanf(f, "%lf %lf %d", &Ta_kij, &Tn_kij, &nT_kij);
    fclose(f);

    f = fopen("constant/species.dat+", "r");

    iread = fscanf(f, "%d %d", &n_pseudospecies_db, &n_purespecies_db);
    if(iread <= 0) Info<< "Input file reading error-----------------" << endl;
    n_species_db = n_pseudospecies_db + n_purespecies_db;

    Info<< "n_pseudoSpecies_db = " << n_pseudospecies_db << "  n_purespecies_db = " << n_purespecies_db << endl;

    fclose(f);

    //Read thermo input parameters from files       
    f = fopen("constant/input.dat+", "r");

    iread=fscanf(f, "%s %d", str_tmp, &n);  // num of species
    if(iread <= 0) Info<< "Input file reading error-----------------" << endl;


    //Allocate memory for arrays
    _NNEW2_(kij, double, n*n);
    _NNEW2_(kij_T, double, nT_kij*n*n);
    _NNEW2_(Pc, double, n);
    _NNEW2_(Tc, double, n);
    _NNEW2_(Vc, double, n);
    _NNEW2_(w, double, n);
    _NNEW2_(MW, double, n);
    _NNEW2_(Tb, double, n);
    _NNEW2_(SG, double, n);
    _NNEW2_(H8, double, n);
    _NNEW2_(k, double, n);
    _NNEW2_(dm, double, n);
    _NNEW2_(idx_species, int, n);
    _NNEW2_(idx_groups, int, n);
    _NNEW2_(x, double, n);
    _NNEW2_(x_eqm0, double, n);
    _NNEW2_(x_eqm1, double, n);

    iread=fscanf(f, "%s", str_tmp);         // index of species
    for (i=0; i<n; i++) iread=fscanf(f, "%d", &idx_species[i]);
    if(iread <= 0) Info<< "Input file reading error-----------------" << endl;

    iread=fscanf(f, "%s", str_tmp);         // index of PPR78 groups
    for (i=0; i<n; i++) iread=fscanf(f, "%d", &idx_groups[i]);
    if(iread <= 0) Info<< "Input file reading error-----------------" << endl;

    Info<< "Species indices in database:  ";
    for(i=0; i<n; i++) Info<< idx_species[i] << "  ";
    Info<< endl;

    Info<< "Species group indices in database:  ";
    for(i=0; i<n; i++) Info<< idx_groups[i] << "  ";
    Info<< endl;

    fclose(f);

    f = fopen("constant/petro.dat+","r");

    for (j=0;j<n;j++) 
    {
        for (i=0;i<n_species_db;i++) 
        {
            double tmp;
            double xPNA[3];

            //[x' y' Tb' Tc' Pc' w' MW' Vc' SG' H_8' xP' xN' xA'];
            iread=fscanf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &tmp, &tmp, &Tb[j], &Tc[j], &Pc[j], &w[j], &MW[j], &Vc[j], 
            &SG[j], &H8[j], &xPNA[0], &xPNA[1], &xPNA[2]);

            if(iread <= 0) Info<< "Input file reading error-----------------" << endl;

            Pc[j] *= 1e5;  // bar --> Pa
            dm[j] = 0.0; // Dipole moment
            k [j] = 0.0;

            if (i==idx_species[j]) 
            {
                // SPECIAL DATA POINTS BUILT IN THE CODE
                // 1. water
                if (fabs(Tc[j]-647)<0.5 && fabs(Pc[j]-220.6e5)<0.5e5) 
                {
                    Vc[j] = 57.1; // cm^3/mol
                    dm[j] = 1.85; // Debye
                    k [j] = 0.076;

                }                                
                // 2. toluene
                else if (fabs(Tc[j]-591.8)<0.5 && fabs(Pc[j]-41.1e5)<0.5e5)
                {
                    Vc[j] = 316.0;
                    dm[j] = 0.36; // Debye

                }// 3. n-decane
                else if (fabs(Tc[j]-617.7)<0.5 && fabs(Pc[j]-21.1e5)<0.5e5)
                {
                    Vc[j] = 600.0;
                }// 4. n-C30
                else if (fabs(Tc[j]-844.0)<0.5 && fabs(Pc[j]- 8.0e5)<0.5e5) {
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

    f = fopen("constant/kij_T.dat+", "r");

    for(iT=0; iT<nT_kij; iT++)
    {
        for(i=0; i<n; i++)
        {
            for(j=0; j<n; j++)
            {
                idx = iT*n*n + i*n + j;
                iread = fscanf(f, "%lf", &kij_T[idx]);
            }
        }
    }
    
    fclose(f);

    f = fopen("constant/solver_inputs.dat+","r");
    iread = fscanf(f, "%d", &nT);
    iread = fscanf(f, "%lf", &Tmin);
    iread = fscanf(f, "%lf", &Tmax);
    iread = fscanf(f, "%d", &nX);
    iread = fscanf(f, "%lf", &Xmin);
    iread = fscanf(f, "%lf", &Xmax);
    iread = fscanf(f, "%lf", &P);
    Info << "Tmin = " << Tmin << tab <<
            "Tmax = " << Tmax << tab <<
            "nT = "   << nT   << tab <<
            "Xmin = " << Xmin << tab <<
            "Xmax = " << Xmax << tab <<
            "nx = "   << nX   << tab <<
            "P = "    << P    << endl;

    dT = (Tmax - Tmin)/(nT -1);
    dX = (Xmax - Xmin)/(nX -1);
    fclose(f);

    f = fopen("prop.dat+", "w");
    T = Tmin;
    for(iT=0;iT<nT;iT++)
    {
        //        Info << "T = " << T << " K, " << endl;
        plicFuncs::calc_kij_from_table(T, n, Ta_kij, Tn_kij, nT_kij, kij_T, kij);
        x[0] = Xmin;
        x[1] = 1- x[0];
        for(iX=0;iX<nX;iX++)
        {
            //  Info << "xWater = " << x[1] << endl;

            calc_v_cvig_(&P,&T,x,&n,Pc,Tc,w,MW,Tb,SG,H8, kij, &V,&CvIG);
            vis_n_cond_(&P,&T,&n,Pc,Tc,Vc,w,MW,k,dm,x,&CvIG,&V,&cond_tmp,&mu_tmp);

            Info << "Mu = " << mu_tmp << " and K = " << cond_tmp << "at T = " << T << " and X_oil = " << x[0] << endl;
            

            density_pr_eos2_(&P, &T, x, &n, Pc, Tc, w, MW, kij, &rho_tmp);
            Info<< " rho = " << rho_tmp << endl;

            calc_v_cp_h_(&P, &T, x, &n, Pc, Tc, w, MW, Tb, SG, H8, kij, &v_tmp, &Cp_tmp, &h_tmp);
            Info<< " Cp = " << Cp_tmp << endl;

            fprintf(f,"%lf    %lf   %lf   %lf   %lf   %lf" "\n", Cp_tmp, mu_tmp, cond_tmp, rho_tmp, T, x[1]);
            x[0] = min(1, x[0] + dX);
            x[1] = 1- x[0];
        }
        T = T + dT;
    }
    fclose(f);

    f = fopen("is2ph.dat+", "w");
    T = Tmin;
    for(iT=0;iT<nT;iT++)
    {
        //        Info << "T = " << T << " K, " << endl;
        plicFuncs::calc_kij_from_table(T, n, Ta_kij, Tn_kij, nT_kij, kij_T, kij);
        x[0] = Xmin;
        x[1] = 1- x[0];
        for(iX=0;iX<nX;iX++)
        {
            //  Info << "xWater = " << x[1] << endl;
            //            plicFuncs::calc_h_from_T(x, P, T, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tn_kij, nT_kij, kij, h_tmp);
            binarylle_(&P, &T, &n, Pc, Tc, w, kij,x_eqm0, x_eqm1, &n_miscible );
            // h_par vector in J/mol [oil water]
            
            if(n_miscible > 0)
            {
                is_2phase = false;
            }
            else if(x_eqm0[1] < x[1] && x_eqm1[1] > x[1])
            {
                is_2phase = true;
            }
            else
            {
                is_2phase = false;
            }

            fprintf(f,"%d    %lf    %lf    %lf    %lf " "\n", is_2phase, T, x_eqm0[1], x[1], x_eqm1[1]);
            x[0] = min(1, x[0] + dX);
            x[1] = 1- x[0];
        }
        T = T + dT;
    }
    fclose(f);

    _DDELETE2_(idx_species);
    _DDELETE2_(idx_groups);

    _DDELETE2_(kij);
    _DDELETE2_(kij_T);
    _DDELETE2_(Pc);
    _DDELETE2_(Tc);
    _DDELETE2_(Vc);
    _DDELETE2_(w);
    _DDELETE2_(MW);
    _DDELETE2_(dm);
    _DDELETE2_(k);
    _DDELETE2_(H8);
    _DDELETE2_(Tb);
    _DDELETE2_(SG);
    _DDELETE2_(x);
    _DDELETE2_(x_eqm0);
    _DDELETE2_(x_eqm1);
   
    Info<< "ExecutionTime = "
        << runTime.elapsedCpuTime()
        << " s\n\n" << endl;     

    Info<< "End\n" << endl;    

    return 0;
}


// ************************************************************************* //
