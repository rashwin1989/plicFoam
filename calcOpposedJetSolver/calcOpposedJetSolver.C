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
    calcIstAxial

Description
    Calculates total normalized C10 inside droplet at each time

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

#include "fvc.H"
#include "OFstream.H"
#include "fileName.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"            

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    string outputFileName = "T_Y00_profiles.dat+";
    fileName outputFile1(outputFileName);
    OFstream os(runTime.path()/outputFile1);

    int i;
    const vector zeroVec(0,0,0);
    const double sqrt2            = 1.41421356;
    const scalar rho0             = 180;
    const scalar rho1             = 800;
    const scalarField& V          = mesh.V();    
    const vectorField& coordinate = mesh.C();    
    // Parameters of the Uniform mesh
    scalar dx;
    scalar dy;
    scalar dz;
    scalar yTop;
    scalar yBottom;
    scalar xLeft;
    scalar xRight;
    scalar tSteady1; // Time at which plots are made - lb
    scalar tSteady2; // Time at which plots are made - ub
    dx       = runTime.controlDict().lookupOrDefault("dx",8e-6);    
    dy       = runTime.controlDict().lookupOrDefault("dy",8e-6);    
    dz       = runTime.controlDict().lookupOrDefault("dz",0.001);    
    yTop     = runTime.controlDict().lookupOrDefault("yTop",0.002);    
    yBottom  = runTime.controlDict().lookupOrDefault("yBottom",0.0);    
    xLeft    = runTime.controlDict().lookupOrDefault("xLeft",0.0);    
    xRight   = runTime.controlDict().lookupOrDefault("xRight",0.002);    
    tSteady1  = runTime.controlDict().lookupOrDefault("tSteady1",0.79);    
    tSteady2  = runTime.controlDict().lookupOrDefault("tSteady2",0.81);
    
    scalar nCellsInX = (xRight - xLeft)/dx;
    scalar nCellsInY = (yTop - yBottom)/dy;
    if(nCellsInX - floor(nCellsInX) != 0)  // Want integer # of cells.Can perhaps be done in a better way!
    { 
           Info << "nCellsInX = " << nCellsInX << ". Resolution not entered consistently in controlDict." << endl;
    }
    label nCellsX     = floor(nCellsInX);
    label nCellsY     = floor(nCellsInY);
    //    label n = nCellsX;
    Info << nCellsX << endl;
    double xCoordinate[nCellsX];
    double Temperature[nCellsX];
    double Y[nCellsX];
    double ThDiff[nCellsX];
    double Diff[nCellsX];
    double Le[nCellsX];
    double velX[nCellsX];
    double velY[nCellsX];
    int phaseStability[nCellsX];
    // Initializing dynamic Arrays
    for (i=0; i < nCellsX; i++)  
    {
        xCoordinate[i]              = 0;
        Temperature[i]                        = 0;
        Y[i]                        = 0;
        velX[i]                     = 0;
        velY[i]                     = 0;
        phaseStability[i]           = 0;
        ThDiff[i]                   = 0;
        Diff[i]                     = 0;
        Le[i]                       = 0;
    }

    //----------------------------------------
    //Begin thermodynamic set-up

    double *Pc, *Tc, *Vc, *w, *MW, *Tb, *SG, *H8, *k, *dm;
    double Ta_kij, Tn_kij, P_thermo, T, h_tmp, T_tmp;
    int n, nT_kij, nT, nX, j, iT, iX, idx;
    double *kij, *kij_T;
    double *x;

    int n_species_db, iread;
    char str_tmp[50];
    FILE *f;
    FILE *ft;
    FILE *fx;
    int n_pseudospecies_db, n_purespecies_db;
    int *idx_species, *idx_groups;

    P_thermo = 30000000; // System pressure: 30MPa
    int n_miscible; // For binaryLLE output; >0: miscible; <=0: immiscible
    double *x1, *x2; // For binaryLLE output

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
    _NNEW2_(x1, double, n);
    _NNEW2_(x2, double, n);

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

    // Done thermodynamic set-up
    //----------------------------------------

    // Starting time loop
    while (runTime.run())
    {        
        runTime++;        
        scalar t = runTime.value(); 
        Info << "t = " << t << endl;
        if ((t>tSteady1)&&(t<tSteady2))   // Read fields at this time only
        {

            Info << "Entered (t == tSteady) region" << endl;
        // Reading fields

        volScalarField T0
        (
            IOobject
            (
                "T0",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        ); 

        volScalarField Y00
        (
            IOobject
            (
                "Y00",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        ); 

        volScalarField Dh0
        (
            IOobject
            (
                "Dh0",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        ); 

        volScalarField D00
        (
            IOobject
            (
                "D00",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        ); 

        volScalarField Le0
        (
            IOobject
            (
                "Le0",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        ); 

        volVectorField U
            (
               IOobject
            (
             "U",
             runTime.timeName(),
             mesh,
             IOobject::MUST_READ
            ), 
            mesh
            );

        // Done Reading Fields

        /*        
        // Calculating Derived Fields

         volVectorField vorticity
            (
            IOobject
            (
                "vorticity",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            fvc::curl(U)
            );

        volScalarField magVorticity
            (
            IOobject
            (
                "magVorticity",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mag(vorticity)
            );

         volScalarField QCriterion
            (
            IOobject
            (
                "QCriterion",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            0.5*mag(skew(fvc::grad(U)))*mag(skew(fvc::grad(U))) - mag(symm(fvc::grad(U)))*mag(symm(fvc::grad(U)))
            );

         volScalarField magStrainRate
            (
            IOobject
            (
                "magStrainRate",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            sqrt2*mag(symm(fvc::grad(U)))
            );

            // Done Calculating Derived Fields
*/

           // Populating arrays from read fields - y averaging
            forAll(T0.internalField(), cellI)
            { 
                scalar xCoordinate = coordinate[cellI].x();
                scalar yCoordinate = coordinate[cellI].y(); 
                scalar xIndex = (xCoordinate - xLeft)/dx;
                label xi = floor(xIndex);
                /* Add quantities for all y-coordinate steps, normalized before printing */
                Temperature[xi]  = Temperature[xi] + T0[cellI];
                Y[xi]            = Y[xi] + Y00[cellI];
                ThDiff[xi]       = ThDiff[xi] + Dh0[cellI];
                Diff[xi]         = Diff[xi] + D00[cellI];
                Le[xi]           = Le[xi] + Le0[cellI];
                velX[xi]         = velX[xi] + U[cellI].x();
                velY[xi]         = velY[xi] + U[cellI].y();
            }   // Done Calculating y-Averaged Quantities

        Info << "Time = " << runTime.timeName() << nl << endl;
        }   
    }   // Time loop over

    // Normalize quantities by # of cells in Y
    for (i=0; i < nCellsX; i++)
    {
        Temperature[i]  = Temperature[i]/nCellsY;
        Y[i]  = Y[i]/nCellsY;
        ThDiff[i]  = ThDiff[i]/nCellsY;
        Diff[i]  = Diff[i]/nCellsY;
        Le[i]  = Le[i]/nCellsY;
        velX[i]  = velX[i]/nCellsY;
        velY[i]  = velY[i]/nCellsY;
        T = Temperature[i];
        Info << T << endl;
        Info << "Entering Binary LLE" << endl;
        plicFuncs::calc_kij_from_table(T, n, Ta_kij, Tn_kij, nT_kij, kij_T, kij);
        binarylle_(&P_thermo,&T,&n,Pc,Tc,w,kij,x1,x2,&n_miscible);
        Info << "Binary LLE Done" << endl;
        if(n_miscible > 0) phaseStability[i] = 1;        
        if(n_miscible <= 0) 
        {
            double y_oil = Y[i];
            double mw_oil = MW[0];
            double mw_water = MW[1];
            double x_oil = y_oil/mw_oil/(y_oil/mw_oil + (1-y_oil)/mw_water);
            if(x_oil > x1[0]) phaseStability[i] = 1;
            else if(x_oil < x2[0]) phaseStability[i] = 1;
            else phaseStability[i] = -1;
        }
        scalar xCoordinate = i*dx + dx/2;
        os<< xCoordinate << tab <<  Temperature[i]<< tab << Y[i] << tab << ThDiff[i] << tab << Diff[i] << tab << Le[i] << tab << velX[i] <<tab << velY[i]<< tab << phaseStability[i] << endl;
    }

   
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
    _DDELETE2_(idx_species);
    _DDELETE2_(idx_groups);
    _DDELETE2_(x);
    _DDELETE2_(x1);
    _DDELETE2_(x2);

    Info<< "ExecutionTime = "
        << runTime.elapsedCpuTime()
        << " s\n\n" << endl;     

    Info<< "End\n" << endl;    

    return 0;
}


// ************************************************************************* //
