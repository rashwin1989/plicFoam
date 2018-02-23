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
    plicAdvDiffIstFlowFoam

Description
    Solver for 2 incompressible, immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing with 
    piecewise linear interface reconstruction. 
    Advection and diffusion of species within each bulk phase.
    Interfacial species mass transfer calculation.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "plic.H"
#include "plicFuncs.H"
//#include "interfaceProperties.H"
#include "IOdictionary.H"
#include "centredCPCCellToCellStencilObject.H"
#include "centredCFCCellToCellStencilObject.H"
#include "pimpleControl.H"
#include "fixedFluxPressureFvPatchScalarField.H"

/*
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include "MACROS.H"
#include "PR_EoS.h"
#include "myUmfpack.h"
#include "vis_n_therm.h"
#include "transport_LLE2.h"
#include "gsl_optimizations2.h"
#include "Maxwell_Stefan_flux.h"
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{    
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"        
   
    pimpleControl pimple(mesh);
 
    #include "readTimeControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createPrghCorrTypes.H"
    #include "correctPhi.H"
    #include "CourantNo.H"
    #include "diffCourantNo.H"
    #include "setInitialDeltaT.H"

    /*
    double Z;
    double *MW;
    double *x;
    double *Dij;
    double *rhs_flux;
    double *flux_m;
    T_UMFPACK flux_umf;

    _NEW_(MW, double, nSpecies);    
    _NEW_(x, double, nSpecies);    
    _NEW_(Dij, double, nSpecies*nSpecies);
    _NEW_(rhs_flux, double, nSpecies);
    _NEW_(flux_m, double, nSpecies);    

    // initialize umfpack
    initialUmfpack(&flux_umf);
    flux_memoryMatVec(&flux_umf, nSpecies);
    flux_compRow_setRowCol(&flux_umf, nSpecies);
    
    int iread;
    FILE *f;
    f = fopen("MW.dat+", "r");
    for(int i=0; i<nSpecies; i++) iread=fscanf(f, "%lf", &MW[i]);
    fclose(f);

    for(int i=0; i<nSpecies; i++) x[i] = Y1[i].internalField()[0]/MW[i];

    for(int i=0; i<nSpecies; i++)
    {
        for(int j=0; j<nSpecies; j++)
        {
            Dij[i+nSpecies*j] = D1_0[0].value();
        }
    }

    for(int i=0; i<nSpecies; i++) rhs_flux[i] = -gradf_Y1[i].internalField()[0];

    Z = 0.5;

    Maxwell_Stefan_mass_flux(Z, nSpecies, MW, x, Dij, rhs_flux, flux_m, &flux_umf);
    */


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        //--------------------------------------------------------------------//
        // At start of each time-step, read time controls like max Courant No.
        // Make sure to use a "CourantNo.H" file which calculates both Adv
        // and Diff Courant numbers. Use a "setDeltaT.H" file which scales dt
        // according to both Adv and Diff Courant numbers
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "diffCourantNo.H"
        #include "setDeltaT.H"
        //--------------------------------------------------------------------//

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        //--------------------------------------------------------------------//
        // code to check and debug big differences between face flux field
        // and velocity field linearly interpolated to mesh faces
        const labelList& own = mesh.owner();
        const labelList& nei = mesh.neighbour();
        surfaceVectorField U_interp(linearInterpolate(U));
 
        if(phi_interp_debug)
        {            
            forAll(phi.internalField(),faceI)
            {
                scalar curPhi = phi.internalField()[faceI];
                vector curSf = mesh.Sf()[faceI];
                scalar curMagSf = mesh.magSf()[faceI];
                scalar curPhiU = curPhi/curMagSf;         
                scalar curU_interp = (U_interp.internalField()[faceI] & curSf)/curMagSf;

                if(mag((curPhiU - curU_interp)/(curU_interp + VSMALL)) > 0.5)
                {
                    Info<< "Face " << faceI << nl
                        << "Flux = " << curPhi << "  Face vel. = " << curPhiU 
                        << "  Interp. U = " << curU_interp << nl
                        << "Own U = " << U.internalField()[own[faceI]] 
                        << "  Nei U = " << U.internalField()[nei[faceI]] << endl;
                }            
            }            
        }
        //--------------------------------------------------------------------//

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if(pimple.firstIter())
            {
                Info<< "Calculating two-phase advective fluxes" << endl;
            
                interface.calc_2ph_advFluxes(Y1, Y0, advFlux_Y1, advFlux_Y0);
     
                Info<< "ExecutionTime = "
                    << runTime.elapsedCpuTime()
                    << " s" << endl; 

                #include "alpha1Eqn.H"        

                Info<< "ExecutionTime = "
                    << runTime.elapsedCpuTime()
                    << " s" << nl << endl; 

                #include "YAdvEqn.H"

                Info<< "ExecutionTime = "
                    << runTime.elapsedCpuTime()
                    << " s" << nl << endl;

                interface.intfc_correct();                

                Info<< "ExecutionTime = "
                    << runTime.elapsedCpuTime()
                    << " s" << nl << endl;

                #include "ist.H"

                Info<< "ExecutionTime = "
                    << runTime.elapsedCpuTime()
                    << " s" << nl << endl;

                interface.intfc_correct();                

                Info<< "ExecutionTime = "
                    << runTime.elapsedCpuTime()
                    << " s" << nl << endl;

                #include "YDiffEqn.H"

                Info<< "ExecutionTime = "
                    << runTime.elapsedCpuTime()
                    << " s" << nl << endl;

                //intfcProp.correct();
                #include "curvature.H"

                for(label i=0; i<nSpecies; i++)
                {
                    C_phAvg[i] = C0[i] + C1[i];
                }
            }

            #include "UEqn.H"

            while(pimple.correct())
            {
                #include "pEqn.H"
                
                Info<< "ExecutionTime = "
                    << runTime.elapsedCpuTime()
                    << " s" << endl; 
            }

            Info<< nl << endl; 
        }

        runTime.write();        
    }    

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
