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

Class
    Foam::centredFPCCellToFaceStencil

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef centredFPCCellToFaceStencil_H
#define centredFPCCellToFaceStencil_H

#include "extendedCentredCellToFaceStencil.H"
#include "FPCCellToFaceStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class centredFPCCellToFaceStencil Declaration
\*---------------------------------------------------------------------------*/

class centredFPCCellToFaceStencil
:   
    public extendedCentredCellToFaceStencil
{

public:

    ClassName("centredFPCCellToFaceStencil");

    // Constructors

        //- Construct from uncompacted face stencil
        explicit centredFPCCellToFaceStencil
        (
            const fvMesh& mesh
        )
        :            
            extendedCentredCellToFaceStencil(FPCCellToFaceStencil(mesh))
        {
            if (extendedCellToFaceStencil::debug)
            {
                Info<< "Generated FPC centred stencil "
                    << nl << endl;
                writeStencilStats(Info, stencil(), map());
            }
        }


    //- Destructor
    virtual ~centredFPCCellToFaceStencil()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
