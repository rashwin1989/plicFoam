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
    Foam::plicLeastSquaresGrad

Description
    Gradient calculated using weighted least-squares on an cell-point-cell stencil.

SourceFiles
    plicLeastSquaresGrad.C
    plicLeastSquaresVectors.H
    plicLeastSquaresVectors.C

\*---------------------------------------------------------------------------*/

#ifndef plicLeastSquaresGrad_H
#define plicLeastSquaresGrad_H

#include "volFields.H"
#include "plicLeastSquaresVectors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


/*---------------------------------------------------------------------------*\
                      Class plicLeastSquaresGrad Declaration
\*---------------------------------------------------------------------------*/

class plicLeastSquaresGrad
{
    // Private Data

    plicLeastSquaresVectors lsv_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        plicLeastSquaresGrad(const plicLeastSquaresGrad&);

        //- Disallow default bitwise assignment
        void operator=(const plicLeastSquaresGrad&);


public:

    //- Runtime type information
    TypeName("plicLeastSquaresGrad");


    // Constructors

        //- Construct from Istream
        plicLeastSquaresGrad(const fvMesh& mesh)
        :
        lsv_(mesh)
        {}

    // Destructor
        virtual ~plicLeastSquaresGrad();


    // Member Functions

        void correctLeastSquaresVectors(const volVectorField& meshC)
        {
            lsv_.correct(meshC);
        }

        //- Return the gradient of the given field       
        tmp<volVectorField> calcGrad
        (
            const volScalarField&,            
            const word&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
