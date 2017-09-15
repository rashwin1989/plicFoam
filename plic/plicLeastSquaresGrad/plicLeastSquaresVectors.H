/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::fv::plicLeastSquaresVectors

Description
    Least-squares gradient scheme vectors

See Also
    Foam::fv::plicLeastSquaresGrad

SourceFiles
    plicLeastSquaresVectors.C

\*---------------------------------------------------------------------------*/

#ifndef plicLeastSquaresVectors_H
#define plicLeastSquaresVectors_H

#include "centredCPCCellToCellStencil.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class plicLeastSquaresVectors Declaration
\*---------------------------------------------------------------------------*/

class plicLeastSquaresVectors
{
    // Private data

        //- cell stencil
        centredCPCCellToCellStencil cellStencil_;

        //- mesh cell centres field with extrapolated boundary cells
        volVectorField C_;

        //- Least-squares gradient vectors
        List<List<vector> > vectors_;


    // Private Member Functions

        //- Calculate Least-squares gradient vectors
        void calcplicLeastSquaresVectors();


public:

    // Declare name of the class and its debug switch
    TypeName("plicLeastSquaresVectors");


    // Constructors

        //- Construct given an fvMesh and the minimum determinant criterion
        plicLeastSquaresVectors
        (
            const fvMesh&
        );


    //- Destructor
    virtual ~plicLeastSquaresVectors();


    // Member functions

        void correct(const volVectorField&);

        //- Return reference to the stencil
        const centredCPCCellToCellStencil& stencil() const
        {
            return cellStencil_;
        }

        //- Return reference to the least square vectors
        const List<List<vector> >& vectors() const
        {
            return vectors_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
