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

\*---------------------------------------------------------------------------*/

#include "FPCCellToFaceStencil.H"
#include "syncTools.H"
#include "emptyPolyPatch.H"
#include "dummyTransform.H"
#include "faceList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Calculates per point the neighbour data (= pointCells)
void Foam::FPCCellToFaceStencil::calcPointBoundaryData
(
    const boolList& isValidBFace,
    const labelList& boundaryPoints,
    Map<labelList>& neiGlobal
) const
{
    neiGlobal.resize(2*boundaryPoints.size());

    labelHashSet pointGlobals;

    forAll(boundaryPoints, i)
    {
        label pointI = boundaryPoints[i];

        neiGlobal.insert
        (
            pointI,
            calcFaceCells
            (
                isValidBFace,
                mesh().pointFaces()[pointI],
                pointGlobals
            )
        );
    }

    syncTools::syncPointMap(mesh(), neiGlobal, unionEqOp(), dummyTransform());
}


// Calculates per face the point connected data (= cell or boundary in global
// numbering).
void Foam::FPCCellToFaceStencil::calcFaceStencil
(
    labelListList& faceStencil
) const
{
    const polyBoundaryMesh& patches = mesh().boundaryMesh();
    const label nBnd = mesh().nFaces()-mesh().nInternalFaces();
    const labelList& own = mesh().faceOwner();
    const labelList& nei = mesh().faceNeighbour();



    // Determine neighbouring global cell
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList neiGlobalCell(nBnd);
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            label faceI = pp.start();

            forAll(pp, i)
            {
                neiGlobalCell[faceI-mesh().nInternalFaces()] =
                    globalNumbering().toGlobal(own[faceI]);
                faceI++;
            }
        }
    }
    syncTools::swapBoundaryFaceList(mesh(), neiGlobalCell);



    // Determine on coupled points the point cells on the other side
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Calculate points on coupled patches
    labelList boundaryPoints
    (
        allCoupledFacesPatch()().meshPoints()
    );

    // Mark boundary faces to be included in stencil (i.e. not coupled or empty)
    boolList isValidBFace;
    validBoundaryFaces(isValidBFace);

    // Swap pointCells for coupled points. Note: use Map (from pt labels to list 
    // of cell labels) since we've got
    // syncTools::syncPointMap for those. 
    Map<labelList> neiGlobal;
    calcPointBoundaryData
    (
        isValidBFace,
        boundaryPoints,
        neiGlobal
    );



    // Construct stencil in global numbering
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    faceStencil.setSize(mesh().nFaces());

    // Do coupled points first

    forAll(boundaryPoints, i)
    {
        label pointI = boundaryPoints[i];

        const labelList& pGlobals = neiGlobal[pointI];

        // Distribute to all pointFaces
        const labelList& pFaces = mesh().pointFaces()[pointI];

        forAll(pFaces, j)
        {
            label faceI = pFaces[j];

            // Insert pGlobals into faceStencil.
            merge(-1, -1, pGlobals, faceStencil[faceI]);
        }
    }
    neiGlobal.clear();


    // Do remaining points by looping over all faces

    // Work arrays
    DynamicList<label> fPointsSet;
    DynamicList<label> pFacesSet;
    labelHashSet faceStencilSet;

    const faceList& meshFaces = mesh().faces();

    for (label faceI = 0; faceI < mesh().nInternalFaces(); faceI++)
    {
        label globalOwn = globalNumbering().toGlobal(own[faceI]);
        label globalNei = globalNumbering().toGlobal(nei[faceI]);

        // Convert any existing faceStencil (from coupled points) into
        // set and operate on this.

        faceStencilSet.clear();

        // Insert all but global owner and neighbour
        forAll(faceStencil[faceI], i)
        {
            label globalI = faceStencil[faceI][i];
            if (globalI != globalOwn && globalI != globalNei)
            {
                faceStencilSet.insert(globalI);
            }
        }
        faceStencil[faceI].clear();

        // Collect all point connected (internal) cells
        const face& fPoints = meshFaces[faceI];

        forAll(fPoints, i)
        {
            label pointI = fPoints[i];

            insertFaceCells
            (
                globalOwn,
                globalNei,
                isValidBFace,
                mesh().pointFaces()[pointI],
                faceStencilSet
            );
        }

        // Extract, guarantee owner first, neighbour second.
        faceStencil[faceI].setSize(faceStencilSet.size()+2);
        label n = 0;
        faceStencil[faceI][n++] = globalOwn;
        faceStencil[faceI][n++] = globalNei;
        forAllConstIter(labelHashSet, faceStencilSet, iter)
        {
            if (iter.key() == globalOwn || iter.key() == globalNei)
            {
                FatalErrorIn("FPCCellToFaceStencil::calcFaceStencil(..)")
                    << "problem:" << faceStencilSet
                    << abort(FatalError);
            }
            faceStencil[faceI][n++] = iter.key();
        }
    }
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        label faceI = pp.start();

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label globalOwn = globalNumbering().toGlobal(own[faceI]);
                label globalNei = neiGlobalCell[faceI-mesh().nInternalFaces()];

                // Convert any existing faceStencil (from coupled points) into
                // set and operate on this.

                faceStencilSet.clear();

                // Insert all but global owner and neighbour
                forAll(faceStencil[faceI], i)
                {
                    label globalI = faceStencil[faceI][i];
                    if (globalI != globalOwn && globalI != globalNei)
                    {
                        faceStencilSet.insert(globalI);
                    }
                }
                faceStencil[faceI].clear();

                // Collect all point connected (internal) cells
                const face& fPoints = meshFaces[faceI];

                forAll(fPoints, i)
                {
                    label pointI = fPoints[i];

                    insertFaceCells
                    (
                        globalOwn,
                        globalNei,
                        isValidBFace,
                        mesh().pointFaces()[pointI], 
                        faceStencilSet
                    );
                }

                // Extract, guarantee owner first, neighbour second.
                faceStencil[faceI].setSize(faceStencilSet.size()+2);
                label n = 0;
                faceStencil[faceI][n++] = globalOwn;
                faceStencil[faceI][n++] = globalNei;
                forAllConstIter(labelHashSet, faceStencilSet, iter)
                {
                    if (iter.key() == globalOwn || iter.key() == globalNei)
                    {
                        FatalErrorIn
                        (
                            "FPCCellToFaceStencil::calcFaceStencil(..)"
                        )   << "problem:" << faceStencilSet
                            << abort(FatalError);
                    }
                    faceStencil[faceI][n++] = iter.key();
                }

                if (n != faceStencil[faceI].size())
                {
                    FatalErrorIn("problem") << "n:" << n
                        << " size:" << faceStencil[faceI].size()
                        << abort(FatalError);
                }

                faceI++;
            }
        }
        else if (!isA<emptyPolyPatch>(pp))
        {
            forAll(pp, i)
            {
                label globalOwn = globalNumbering().toGlobal(own[faceI]);

                // Convert any existing faceStencil (from coupled points) into
                // set and operate on this.

                faceStencilSet.clear();

                // Insert all but global owner and neighbour
                forAll(faceStencil[faceI], i)
                {
                    label globalI = faceStencil[faceI][i];
                    if (globalI != globalOwn)
                    {
                        faceStencilSet.insert(globalI);
                    }
                }
                faceStencil[faceI].clear();

                // Collect all point connected (internal) cells
                const face& fPoints = meshFaces[faceI];

                forAll(fPoints, i)
                {
                    label pointI = fPoints[i];

                    insertFaceCells
                    (
                        globalOwn,
                        -1,
                        isValidBFace,
                        mesh().pointFaces()[pointI],
                        faceStencilSet
                    );
                }

                // Extract, guarantee owner first, neighbour second.
                faceStencil[faceI].setSize(faceStencilSet.size()+1);
                label n = 0;
                faceStencil[faceI][n++] = globalOwn;
                forAllConstIter(labelHashSet, faceStencilSet, iter)
                {
                    if (iter.key() == globalOwn)
                    {
                        FatalErrorIn
                        (
                            "FPCCellToFaceStencil::calcFaceStencil(..)"
                        )   << "problem:" << faceStencilSet
                            << abort(FatalError);
                    }
                    faceStencil[faceI][n++] = iter.key();
                }

                faceI++;
            }
        }
    }


    for (label faceI = 0; faceI < mesh().nInternalFaces(); faceI++)
    {
        label globalOwn = globalNumbering().toGlobal(own[faceI]);
        if (faceStencil[faceI][0] != globalOwn)
        {
            FatalErrorIn("FPCCellToFaceStencil::calcFaceStencil(..)")
                << "problem:" << faceStencil[faceI]
                << " globalOwn:" << globalOwn
                << abort(FatalError);
        }
        label globalNei = globalNumbering().toGlobal(nei[faceI]);
        if (faceStencil[faceI][1] != globalNei)
        {
            FatalErrorIn("FPCCellToFaceStencil::calcFaceStencil(..)")
                << "problem:" << faceStencil[faceI]
                << " globalNei:" << globalNei
                << abort(FatalError);
        }
    }


    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label faceI = pp.start()+i;

                label globalOwn = globalNumbering().toGlobal(own[faceI]);
                if (faceStencil[faceI][0] != globalOwn)
                {
                    FatalErrorIn("FPCCellToFaceStencil::calcFaceStencil(..)")
                        << "problem:" << faceStencil[faceI]
                        << " globalOwn:" << globalOwn
                        << abort(FatalError);
                }
                label globalNei = neiGlobalCell[faceI-mesh().nInternalFaces()];
                if (faceStencil[faceI][1] != globalNei)
                {
                    FatalErrorIn("FPCCellToFaceStencil::calcFaceStencil(..)")
                        << "problem:" << faceStencil[faceI]
                        << " globalNei:" << globalNei
                        << abort(FatalError);
                }
            }
        }
        else if (!isA<emptyPolyPatch>(pp))
        {
            forAll(pp, i)
            {
                label faceI = pp.start()+i;

                label globalOwn = globalNumbering().toGlobal(own[faceI]);
                if (faceStencil[faceI][0] != globalOwn)
                {
                    FatalErrorIn("FPCCellToFaceStencil::calcFaceStencil(..)")
                        << "problem:" << faceStencil[faceI]
                        << " globalOwn:" << globalOwn
                        << abort(FatalError);
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FPCCellToFaceStencil::FPCCellToFaceStencil(const polyMesh& mesh)
:
    cellToFaceStencil(mesh)
{
    // Calculate per face the (point) connected cells (in global numbering)
    labelListList faceStencil;
    calcFaceStencil(faceStencil);

    // Transfer to *this
    transfer(faceStencil);
}


// ************************************************************************* //
