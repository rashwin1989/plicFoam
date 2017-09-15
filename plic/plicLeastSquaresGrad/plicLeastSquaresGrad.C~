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

#include "plicLeastSquaresGrad.H"
#include "plicLeastSquaresVectors.H"
#include "gaussGrad.H"
#include "fvMesh.H"
#include "zeroGradientFvPatchField.H"
#include "centredCPCCellToCellStencil.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(plicLeastSquaresGrad, 0);
}


// * * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * //

Foam::plicLeastSquaresGrad::~plicLeastSquaresGrad()
{}


Foam::tmp<Foam::volVectorField>
Foam::plicLeastSquaresGrad::calcGrad
(
    const volScalarField& vtf,
    const word& name
) const
{
    const fvMesh& mesh = vtf.mesh();

    tmp<volVectorField> tlsGrad
    (
        new volVectorField
        (
            IOobject
            (
                name,
                vtf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector
            (
                "zero",
                vtf.dimensions()/dimLength,
                vector::zero
            ),
            zeroGradientFvPatchField<vector>::typeName
        )
    );
    volVectorField& lsGrad = tlsGrad();
    Field<vector>& lsGradIf = lsGrad;

    const centredCPCCellToCellStencil& stencil = lsv_.stencil();
    const List<List<label> >& stencilAddr = stencil.stencil();
    const List<List<vector> >& lsvs = lsv_.vectors();

    // Construct flat version of vtf
    // including all values referred to by the stencil
    List<scalar> flatVtf(stencil.map().constructSize());

    // Insert internal values
    forAll(vtf, celli)
    {
        flatVtf[celli] = vtf[celli];
    }

    // Insert boundary values
    forAll(vtf.boundaryField(), patchi)
    {
        const fvPatchScalarField& ptf = vtf.boundaryField()[patchi];

        label nCompact =
            ptf.patch().start()
          - mesh.nInternalFaces()
          + mesh.nCells();

        forAll(ptf, i)
        {
            flatVtf[nCompact++] = ptf[i];
        }
    }

    // Do all swapping to complete flatVtf
    stencil.map().distribute(flatVtf);

    // Accumulate the cell-centred gradient from the
    // weighted least-squares vectors and the flattened field values
    forAll(stencilAddr, celli)
    {
        const labelList& compactCells = stencilAddr[celli];
        const List<vector>& lsvc = lsvs[celli];

        forAll(compactCells, i)
        {
            lsGradIf[celli] += lsvc[i]*flatVtf[compactCells[i]];
        }
    }

    // Correct the boundary conditions
    lsGrad.correctBoundaryConditions();
    Foam::fv::gaussGrad<scalar>::correctBoundaryConditions(vtf, lsGrad);

    return tlsGrad;
}


// ************************************************************************* //
