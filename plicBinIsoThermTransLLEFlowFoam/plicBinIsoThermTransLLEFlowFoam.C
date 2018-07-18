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
    plicBinIsoThermTransLLEFlowFoam

Description
    Solver for binary isothermal 2 compressible (density varying with x), 
    partially miscible fluids using a VOF (volume of fluid) 
    phase-fraction based interface capturing with 
    piecewise linear interface reconstruction. 
    Advection and diffusion of species within each bulk phase.
    Interfacial species mass transfer calculation consistent 
    with phase equilibrium and transport constraints at interface.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "plic.H"
#include "plicFuncs.H"
//#include "interfaceProperties.H"
#include "IOdictionary.H"
#include "centredCPCCellToCellStencilObject.H"
#include "centredCFCCellToCellStencilObject.H"
//#include "pimpleControl.H"
#include "fixedFluxPressureFvPatchScalarField.H"

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

#include "MACROS2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{    
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"       
 
    #include "readTimeControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createPrghCorrTypes.H"
    #include "correctPhi.H"
    #include "CourantNo.H"
    #include "diffCourantNo.H"
    #include "setInitialDeltaT.H"

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
        
        deltaT = runTime.deltaTValue();

        for(iOCorr=0; iOCorr<nOCorr; iOCorr++)
        {
            Info<< "Outer corrector " << iOCorr+1 << endl;
            
            for(i=0; i<n; i++)
            {
                diffTerm_Y1[i] = diffTerm_zero;
                diffTerm_Y0[i] = diffTerm_zero;
            }
                       
            #include "ist.H"

            Info<< "ExecutionTime = "
                << runTime.elapsedCpuTime()
                << " s" << nl << endl;

            interface.intfc_correct();
            Info<< "Done interface reconstruction" << endl;

            Info<< "ExecutionTime = "
                << runTime.elapsedCpuTime()
                << " s" << nl << endl;

            //alpha1 == alpha1_old;
            //alpha1.correctBoundaryConditions();

            #include "curvature.H"

            dt = deltaT;
            /*
            #include "correctRho.H"

            Info<< "ExecutionTime = "
                << runTime.elapsedCpuTime()
                << " s" << nl << endl;
                */
            rho = alpha1*rho1 + alpha0*rho0;
            rhoPhi = phiAlpha1*(rho1f - rho0f) + phi*rho0f;

            #include "UEqn.H"

            for(iPCorr=0; iPCorr<nPCorr; iPCorr++)
            {
                Info<< "Pressure corrector " << iPCorr+1 << endl;

                #include "pEqn.H"
                
                Info<< "ExecutionTime = "
                    << runTime.elapsedCpuTime()
                    << " s" << endl; 
            }
            Info<< endl;

            dt = deltaT;
            Info<< "Calculating two-phase advective fluxes" << endl;
            interface.calc_2ph_advFluxes(c1, c0, dt, advFlux_Y1, advFlux_Y0, advFlux_debug, advFlux_debug2, osAdv);
     
            Info<< "ExecutionTime = "
                << runTime.elapsedCpuTime()
                << " s" << nl << endl;            

            #include "alpha1Eqn.H"

            Info<< "ExecutionTime = "
                << runTime.elapsedCpuTime()
                << " s" << nl << endl;
            
            #include "YAdvEqn.H"

            Info<< "ExecutionTime = "
                << runTime.elapsedCpuTime()
                << " s" << nl << endl;
            
            interface.intfc_correct();
            Info<< "Done interface reconstruction" << endl;

            Info<< "ExecutionTime = "
                << runTime.elapsedCpuTime()
                << " s" << nl << endl;

            #include "correct_thermo_trans_prop.H"

            Info<< "ExecutionTime = "
                << runTime.elapsedCpuTime()
                << " s" << nl << endl;

            #include "diff_grad_interp.H"

            Info<< "ExecutionTime = "
                << runTime.elapsedCpuTime()
                << " s" << nl << endl;

            #include "YDiffEqn.H"            

            Info<< "ExecutionTime = "
                << runTime.elapsedCpuTime()
                << " s" << nl << endl;
        }

        for(i=0; i<n; i++)
        {
            C_phAvg[i] = C0[i] + C1[i];
        }

        #include "copyOldFields.H"

        Info<< nl << endl; 

        runTime.write();        
    }

    #include "CLEAN.H"

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
