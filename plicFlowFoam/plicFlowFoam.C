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
    plicFlowFoam

Description
    Solver for 2 incompressible, immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing with 
    piecewise linear interface reconstruction.        

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "plic.H"
#include "plicFuncs.H"
#include "interfaceProperties.H"
#include "IOdictionary.H"
#include "centredCPCCellToCellStencilObject.H"
#include "pimpleControl.H"
#include "fixedFluxPressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"        

    pimpleControl pimple(mesh);

    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "readTimeControls.H"
    #include "createPrghCorrTypes.H"
    #include "correctPhi.H"
    #include "CourantNo.H"
    #include "IFTCourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "IFTCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {   
            if(pimple.firstIter())
            {
                Info<< "Calculating two-phase fluxes" << endl;

                interface.calc_face_phaseFluxes();

                Info<< "ExecutionTime = "
                    << runTime.elapsedCpuTime()
                    << " s" << endl;                 
                
                #include "alpha1Eqn.H"

                Info<< "ExecutionTime = "
                    << runTime.elapsedCpuTime()
                    << " s" << endl; 

                interface.intfc_correct();

                Info<< "ExecutionTime = "
                    << runTime.elapsedCpuTime()
                    << " s" << endl;

                intfcProp.correct();
            }

            #include "UEqn.H"

            while(pimple.correct())
            {
                #include "pEqn.H"
                
                Info<< "ExecutionTime = "
                    << runTime.elapsedCpuTime()
                    << " s" << endl; 
            }

        }

        runTime.write();

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl; 
    }    

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
