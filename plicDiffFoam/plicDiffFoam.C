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
    plicDiffFoam

Description
    Solver for intra-phase diffusion within two immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing algorithm with 
    piecewise linear interface reconstruction and diffusion fluxes between
    arbitrary-shaped half cells

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "plic.H"
#include "plicFuncs.H"
#include "IOdictionary.H"
//#include "centredCPCCellToCellStencilObject.H"
#include "centredCFCCellToCellStencilObject.H"
#include "syncTools.H"

#include "myUmfpack.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"        

    #include "readTimeControls.H"
    #include "createFields.H"
    #include "diffCourantNo.H"
    #include "setInitialDeltaT.H"

    T_UMFPACK umf;
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "diffCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        interface.intfc_correct();

        //Foam::plicFuncs::display_field_vector(interface.C_ph0());
        //Foam::plicFuncs::display_field_vector(interface.C_ph1());

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s" << endl; 

        #include "YDiffEqn.H"
       
        runTime.write();

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl; 
    }    

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
