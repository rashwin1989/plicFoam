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

\*---------------------------------------------------------------------------*/

#include "plicLeastSquaresVectors.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(plicLeastSquaresVectors, 0);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::plicLeastSquaresVectors::plicLeastSquaresVectors
(
    const fvMesh& mesh
)
:
    cellStencil_(mesh),
    C_
    (
        IOobject
        (
            "lsvMeshC",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    vectors_(mesh.nCells())
{
    const volVectorField& C = mesh.C();
    forAll(C_.internalField(), cellI)
    {
        C_.internalField()[cellI] = C.internalField()[cellI];
    }

    forAll(C_.boundaryField(), patchI)
    {
        fvPatchVectorField& pC_ = C_.boundaryField()[patchI];
        const fvPatchVectorField& pC = C.boundaryField()[patchI];
        
        forAll(pC_, faceI)
        {
            pC_[faceI] = pC[faceI];
        }
    }

    calcplicLeastSquaresVectors();
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::plicLeastSquaresVectors::~plicLeastSquaresVectors()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::plicLeastSquaresVectors::correct(const volVectorField& meshC)
{
    forAll(C_.internalField(), cellI)
    {
        C_.internalField()[cellI] = meshC.internalField()[cellI];
    }

    forAll(C_.boundaryField(), patchI)
    {
        fvPatchVectorField& pC_ = C_.boundaryField()[patchI];
        const fvPatchVectorField& pC = meshC.boundaryField()[patchI];
        
        forAll(pC_, faceI)
        {
            pC_[faceI] = pC[faceI];
        }
    }

    calcplicLeastSquaresVectors();
}


void Foam::plicLeastSquaresVectors::calcplicLeastSquaresVectors()
{
    if (debug)
    {
        Info<< "plicLeastSquaresVectors::calcplicLeastSquaresVectors() :"
            << "Calculating least square gradient vectors"
            << endl;
    }

    const fvMesh& mesh = C_.mesh();
    const centredCPCCellToCellStencil& stencil = this->stencil();

    stencil.collectData(C_, vectors_);

    // Create the base form of the dd-tensor
    // including components for the "empty" directions
    symmTensor dd0(sqr((Vector<label>::one - mesh.geometricD())/2));

    forAll (vectors_, i)
    {
        List<vector>& lsvi = vectors_[i];
        symmTensor dd(dd0);

        // The current cell is 0 in the stencil
        // Calculate the deltas and sum the weighted dd
        for (label j=1; j<lsvi.size(); j++)
        {
            lsvi[j] = lsvi[j] - lsvi[0];
            scalar magSqrLsvi = magSqr(lsvi[j]);
            dd += sqr(lsvi[j])/magSqrLsvi;
            lsvi[j] /= magSqrLsvi;
        }

        // Invert dd
        dd = inv(dd);

        // Remove the components corresponding to the empty directions
        dd -= dd0;

        // Finalize the gradient weighting vectors
        lsvi[0] = vector::zero;
        for (label j=1; j<lsvi.size(); j++)
        {
            lsvi[j] = dd & lsvi[j];
            lsvi[0] -= lsvi[j];
        }
    }

    if (debug)
    {
        Info<< "plicLeastSquaresVectors::calcplicLeastSquaresVectors() :"
            << "Finished calculating least square gradient vectors"
            << endl;
    }
}


// ************************************************************************* //
