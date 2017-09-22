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
    plicTFoam

Description
    Solver for 2 incompressible, immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing with 
    piecewise linear interface reconstruction.
    Time-varying velocity field set

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "plic.H"
#include "plicFuncs.H"
#include "IOdictionary.H"
#include "centredCPCCellToCellStencilObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"        

    #include "readTimeControls.H"
    //#include "initContinuityErrs.H"
    #include "createFields.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "U_shear.H"

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

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
                        << "Flux = " << curPhi << "  Face vel. = " << curPhiU << "  Interp. U = " << curU_interp << nl
                        << "Own U = " << U.internalField()[own[faceI]] << "  Nei U = " << U.internalField()[nei[faceI]] << endl;
                }            
            }            
        }

        interface.calc_2ph_advFluxes(Y1, Y0, advFlux_Y1, advFlux_Y0);
     
        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s" << endl; 

        #include "alpha1Eqn.H"        

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s" << endl; 

        #include "YAdvEqn.H"

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s" << endl; 

        #include "U_shear.H"
       
        interface.intfc_correct();

        runTime.write();

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl; 
    }    

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
