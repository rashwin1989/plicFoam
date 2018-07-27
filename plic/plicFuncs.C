/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Author
    Ashwin Raghavan

\*---------------------------------------------------------------------------*/

#include "plicFuncs.H"
#include "OFstream.H"
#include "fileName.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace plicFuncs
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template <class Type>
void display_field
(
    const GeometricField<Type, fvPatchField, volMesh>& fld
)
{
    Info<< "//========================================================================\\" << endl << endl
        << "                             " << fld.name() << " field" << endl << endl
        << "\\========================================================================//" << endl
        << endl       
        << "============================================================================" << endl
        << "                              Internal field" << endl
        << "============================================================================" << endl
        << endl
        << "----------------------------------------------------------------------------" << endl
        << "    Cell index                        Field value" << endl
        << "----------------------------------------------------------------------------" << endl;

    forAll(fld.internalField(), cellI)
    {
        Info<< "        " << cellI << "                           " << fld.internalField()[cellI] << endl;
    }

    Info<< endl << endl
        << "============================================================================" << endl
        << "                              Boundary field" << endl
        << "============================================================================" << endl
        << endl;
    
    const fvMesh& mesh = fld.mesh();
    wordList patchNames(mesh.boundaryMesh().names());

    forAll(fld.boundaryField(), patchI)
    {
        const fvPatchField<Type>& pfld = fld.boundaryField()[patchI];
        const fvPatch& pp = mesh.boundary()[patchI];

        label nf = pp.start();            
        
        Info<< "----------------------------------------------------------------------------" << endl
            << "                       " << patchNames[patchI] << endl
            << "----------------------------------------------------------------------------" << endl << endl
            << "----------------------------------------------------------------------------" << endl
            << "Patch face index    Global face index                Field value" << endl
            << "----------------------------------------------------------------------------" << endl;

        forAll(pfld, i)
        {
            Info<< "        " << i << "                  " << nf << "                         " << pfld[i] << endl;
            nf++;
        }
        Info<< endl;
    }
    
    Info<< "//========================================================================\\" << endl << endl
        << "                         End of " << fld.name() << " field" << endl << endl
        << "\\========================================================================//" << endl
        << endl;    
}


template <class Type>
void write_field
(
    const GeometricField<Type, fvPatchField, volMesh>& fld
)
{
    const fvMesh& mesh = fld.mesh();

    fileName outputFile(fld.name()+"_db");
    OFstream os(mesh.time().path()/mesh.time().timeName()/outputFile);

    os<< "//==========================================================================\\" << endl << endl
      << "                            " << fld.name() << " field" << endl << endl
      << "\\==========================================================================//" << endl
      << endl       
      << "==============================================================================" << endl
      << "                                Internal field" << endl
      << "==============================================================================" << endl
      << endl
      << "------------------------------------------------------------------------------" << endl
      << "    Cell index                        Field value" << endl
      << "------------------------------------------------------------------------------" << endl;

    forAll(fld.internalField(), cellI)
    {
        os<< "        " << cellI << "                           " << fld.internalField()[cellI] << endl;
    }

    os<< endl << endl
      << "==============================================================================" << endl
      << "                                Boundary field" << endl
      << "==============================================================================" << endl
      << endl;
        
    wordList patchNames(mesh.boundaryMesh().names());

    forAll(fld.boundaryField(), patchI)
    {
        const fvPatchField<Type>& pfld = fld.boundaryField()[patchI];
        const fvPatch& pp = mesh.boundary()[patchI];

        label nf = pp.start();            
        
        os<< "------------------------------------------------------------------------------" << endl
          << "                        " << patchNames[patchI] << endl
          << "------------------------------------------------------------------------------" << endl << endl
          << "------------------------------------------------------------------------------" << endl
          << "Patch face index    Global face index                Field value" << endl
          << "------------------------------------------------------------------------------" << endl;

        forAll(pfld, i)
        {
            os<< "        " << i << "                  " << nf << "                         " << pfld[i] << endl;
            nf++;
        }
        os<< endl;
    }
    
    os<< "//==========================================================================\\" << endl << endl
      << "                          End of " << fld.name() << " field" << endl << endl
      << "\\==========================================================================//" << endl
      << endl;         
}


template <class Type>
void print_field
(
    const GeometricField<Type, fvPatchField, volMesh>& fld,
    OFstream& os
)
{
    const fvMesh& mesh = fld.mesh();    

    os<< "//==========================================================================\\" << endl << endl
      << "                            " << fld.name() << " field" << endl << endl
      << "\\==========================================================================//" << endl
      << endl       
      << "==============================================================================" << endl
      << "                                Internal field" << endl
      << "==============================================================================" << endl
      << endl
      << "------------------------------------------------------------------------------" << endl
      << "    Cell index                        Field value" << endl
      << "------------------------------------------------------------------------------" << endl;

    forAll(fld.internalField(), cellI)
    {
        os<< "        " << cellI << "                           " << fld.internalField()[cellI] << endl;
    }

    os<< endl << endl
      << "==============================================================================" << endl
      << "                                Boundary field" << endl
      << "==============================================================================" << endl
      << endl;
        
    wordList patchNames(mesh.boundaryMesh().names());

    forAll(fld.boundaryField(), patchI)
    {
        const fvPatchField<Type>& pfld = fld.boundaryField()[patchI];
        const fvPatch& pp = mesh.boundary()[patchI];

        label nf = pp.start();            
        
        os<< "------------------------------------------------------------------------------" << endl
          << "                        " << patchNames[patchI] << endl
          << "------------------------------------------------------------------------------" << endl << endl
          << "------------------------------------------------------------------------------" << endl
          << "Patch face index    Global face index                Field value" << endl
          << "------------------------------------------------------------------------------" << endl;

        forAll(pfld, i)
        {
            os<< "        " << i << "                  " << nf << "                         " << pfld[i] << endl;
            nf++;
        }
        os<< endl;
    }
    
    os<< "//==========================================================================\\" << endl << endl
      << "                          End of " << fld.name() << " field" << endl << endl
      << "\\==========================================================================//" << endl
      << endl;         
}


template <class Type>
void display_point_field
(
    const GeometricField<Type, pointPatchField, pointMesh>& ptfld
)
{
    Info<< "//========================================================================\\" << endl << endl
        << "                      " << ptfld.name() << " field" << endl << endl
        << "\\========================================================================//" << endl
        << endl       
        << "============================================================================" << endl
        << "                             Internal field" << endl
        << "============================================================================" << endl
        << endl
        << "----------------------------------------------------------------------------" << endl
        << "    Point index                        Field value" << endl
        << "----------------------------------------------------------------------------" << endl;

    forAll(ptfld.internalField(), pointI)
    {
        Info<< "         " << pointI << "                           " << ptfld.internalField()[pointI] << endl;
    }

    Info<< endl << endl
        << "============================================================================" << endl
        << "                             Boundary field" << endl
        << "============================================================================" << endl
        << endl;
    
    const pointMesh& mesh = ptfld.mesh();
    wordList patchNames(mesh().boundaryMesh().names());

    forAll(ptfld.boundaryField(), patchI)
    {
        Field<Type> pptfld(ptfld.boundaryField()[patchI].patchInternalField());    
        
        Info<< "----------------------------------------------------------------------------" << endl
            << "                        " << patchNames[patchI] << endl
            << "----------------------------------------------------------------------------" << endl << endl
            << "----------------------------------------------------------------------------" << endl
            << "Patch point index                Field value" << endl
            << "----------------------------------------------------------------------------" << endl;

        forAll(pptfld, i)
        {
            Info<< "         " << i << "                         " << pptfld[i] << endl;
        }
        Info<< endl;
    }
    
    Info<< "//========================================================================\\" << endl << endl
        << "                      End of " << ptfld.name() << " field" << endl << endl
        << "\\========================================================================//" << endl
        << endl;    
}


template <class Type>
void write_point_field
(
    const GeometricField<Type, pointPatchField, pointMesh>& ptfld,
    const fvMesh& mesh
)
{
    fileName outputFile(ptfld.name()+"_db");
    OFstream os(mesh.time().path()/mesh.time().timeName()/outputFile);

    os<< "//==========================================================================\\" << endl << endl
      << "                           " << ptfld.name() << " field" << endl << endl
      << "\\==========================================================================//" << endl
      << endl       
      << "==============================================================================" << endl  
      << "                                Internal field" << endl
      << "==============================================================================" << endl
      << endl
      << "------------------------------------------------------------------------------" << endl
      << "    Point index                        Field value" << endl
      << "------------------------------------------------------------------------------" << endl;

    forAll(ptfld.internalField(), pointI)
    {
        os<< "         " << pointI << "                           " << ptfld.internalField()[pointI] << endl;
    }

    os<< endl << endl
      << "==============================================================================" << endl
      << "                                Boundary field" << endl
      << "==============================================================================" << endl
      << endl;
        
    wordList patchNames(mesh.boundaryMesh().names());

    forAll(ptfld.boundaryField(), patchI)
    {
        Field<Type> pptfld(ptfld.boundaryField()[patchI].patchInternalField());    
        
        os<< "------------------------------------------------------------------------------" << endl
          << "                       " << patchNames[patchI] << endl
          << "------------------------------------------------------------------------------" << endl << endl
          << "------------------------------------------------------------------------------" << endl
          << "Patch point index                Field value" << endl
          << "------------------------------------------------------------------------------" << endl;

        forAll(pptfld, i)
        {
            os<< "         " << i << "                         " << pptfld[i] << endl;
        }
        os<< endl;
    }
    
    os<< "//==========================================================================\\" << endl << endl
      << "                      End of " << ptfld.name() << " field" << endl << endl
      << "\\==========================================================================//" << endl
      << endl;    
}


void display_cell
(
    const cell& curCell,
    const faceList& faces,
    const pointField& points
)
{
    Info<< "==================================================================================================================================" << endl;
    Info<< "  Cell  n.o.f. = " << curCell.size() << endl;
    Info<< "----------------------------------------------------------------------------------------------------------------------------------" << endl;
    Info<< " Face index                                  Face points          " << endl;
    Info<< "----------------------------------------------------------------------------------------------------------------------------------" << endl;
    forAll(curCell, faceI)
    {
        Info<< "    " << curCell[faceI] << "        ";
        face curFace = faces[curCell[faceI]];
        forAll(curFace, pointI)
        {
            Info<< points[curFace[pointI]] << "  "; 
        }
        Info<< endl;
    }
    Info<< "==================================================================================================================================" << endl << endl;
}


void display_cellInfo
(
    const cellInfo& curCellInfo    
)
{
    const pointField& points = curCellInfo.points();
    const faceList& faces = curCellInfo.faces();

    Info<< "==================================================================================================================================" << endl;
    Info<< "  Cell  n.o.f. = " << faces.size() << endl;
    Info<< "----------------------------------------------------------------------------------------------------------------------------------" << endl;
    Info<< " Face index                                  Face points          " << endl;
    Info<< "----------------------------------------------------------------------------------------------------------------------------------" << endl;

    forAll(faces, faceI)
    {
        Info<< "    " << faceI  << "        ";
        face curFace = faces[faceI];
        forAll(curFace, pointI)
        {
            Info<< points[curFace[pointI]] << "  "; 
        }
        Info<< endl;
    }
    Info<< "==================================================================================================================================" << endl << endl;
}


void display_face
(
    const face& curFace,
    const pointField& points
)
{
    Info<< "==================================================================" << endl;
    Info<< "  Face  n.o.p. = " << curFace.size() << endl;
    Info<< "------------------------------------------------------------------" << endl;
    Info<< " Point index                          Point          " << endl;
    Info<< "------------------------------------------------------------------" << endl;
    forAll(curFace, pointI)
    {
        Info<< "      " << curFace[pointI] << "                          " << points[curFace[pointI]] << endl;        
    }
    Info<< "==================================================================" << endl << endl;
}
 

void display_pointField
(
    const pointField& points
)
{
    Info<< "==================================================================" << endl;
    Info<< "  pointField  n.o.p. = " << points.size() << endl;
    Info<< "------------------------------------------------------------------" << endl;
    Info<< " Point index                          Point          " << endl;
    Info<< "------------------------------------------------------------------" << endl;
    forAll(points, i)
    {
        Info<< "      " << i << "                          " << points[i] << endl;   
    }
    Info<< "==================================================================" << endl << endl;
}


void display_Plane
(
    const Plane& curPlane
)
{
    Info<< "==================================================================" << endl;
    Info<< "  Plane " << endl;
    Info<< "------------------------------------------------------------------" << endl;
    Info<< "  Plane normal:  " << curPlane.normal() << "    Ref pt:  " << curPlane.refPoint() << endl;
    Info<< "==================================================================" << endl << endl;
}


void display_labelList
(
    const labelList& lbls
)
{
    Info<< "==================================================================" << endl;
    Info<< "  labelList  n.o.l = " << lbls.size() << endl;
    Info<< "------------------------------------------------------------------" << endl;
    Info<< "     Index                          Label          " << endl;
    Info<< "------------------------------------------------------------------" << endl;
    forAll(lbls, i)
    {
        Info<< "      " << i << "                           " << lbls[i] << endl;   
    }
    Info<< "==================================================================" << endl << endl;
}


void write_labelList
(
    const labelList& lbls,
    const fvMesh& mesh,
    const word& name
)
{
    fileName outputFile(name+"_db");
    OFstream os(mesh.time().path()/mesh.time().timeName()/outputFile);

    os<< "==================================================================" << endl;
    os<< "  labelList  n.o.l = " << lbls.size() << endl;
    os<< "------------------------------------------------------------------" << endl;
    os<< "     Index                          Label          " << endl;
    os<< "------------------------------------------------------------------" << endl;
    forAll(lbls, i)
    {
        os<< "      " << i << "                           " << lbls[i] << endl;   
    }
    os<< "==================================================================" << endl << endl;
}


void write_boolList
(
    const boolList& lbls,
    const fvMesh& mesh,
    const word& name
)
{
    fileName outputFile(name+"_db");
    OFstream os(mesh.time().path()/mesh.time().timeName()/outputFile);

    os<< "==================================================================" << endl;
    os<< "  boolList  n.o.l = " << lbls.size() << endl;
    os<< "------------------------------------------------------------------" << endl;
    os<< "     Index                          Label          " << endl;
    os<< "------------------------------------------------------------------" << endl;
    forAll(lbls, i)
    {
        os<< "      " << i << "                           " << lbls[i] << endl;   
    }
    os<< "==================================================================" << endl << endl;
}


void write_stencil
(
    const labelListList& stencil,
    const fvMesh& mesh,
    const word& stencilName
)
{
    fileName outputFile(stencilName);
    OFstream os(mesh.time().path()/mesh.time().constant()/outputFile);

    os<< "//========================================================================\\" << nl
        << "                         " << stencilName << nl
        << "\\========================================================================//" << nl
        << nl
        << "No. of mesh cells:            " << mesh.nCells() << nl
        << "No. of mesh bfaces:           " << (mesh.nFaces() - mesh.nInternalFaces()) << nl
        << "Total no. of cells + bfaces:  " << (mesh.nCells() + mesh.nFaces() - mesh.nInternalFaces()) << nl
        << nl
        << "----------------------------------------------------------------------------" << nl
        << "    Face                                     Stencil cells" << nl
        << "----------------------------------------------------------------------------" << endl;

    forAll(stencil, faceI)
    {
        const labelList& curStencil = stencil[faceI];
        os<< "     " << faceI << "                               ";
        forAll(curStencil, cellI)
        {
            os<< curStencil[cellI] << "  ";
        }
        os<< endl;
    }
    os<< "----------------------------------------------------------------------------" << nl << nl
        << "//========================================================================\\" << nl
        << "                        Done " << stencilName << nl
        << "\\========================================================================//" << nl
        << endl;
        
}


template<class Type>
void write_flatFld
(
    const List<Type>& flatFld,
    const GeometricField<Type, fvPatchField, volMesh>& fld
)
{
    const fvMesh& mesh = fld.mesh();

    fileName outputFile(fld.name()+"_flatFld");
    OFstream os(mesh.time().path()/mesh.time().timeName()/outputFile);

    os<< "//==========================================================================\\" << nl 
        << nl
        << "                        " << fld.name() << " flat field" << nl << nl
        << "\\==========================================================================//" << nl
        << nl       
        << "==============================================================================" << nl
        << "                                Internal field" << nl
        << "==============================================================================" << nl
        << nl
        << "------------------------------------------------------------------------------" << nl
        << "cell index    flatFld index          field value          flatFld value" << nl
        << "------------------------------------------------------------------------------" << endl;

    label nFlat = 0;

    forAll(fld, cellI)
    {
        os<< "    " << cellI << "          " << cellI << "                " << fld.internalField()[cellI]  << "                " << flatFld[nFlat] << endl;
        nFlat++;
    }

    os<< nl << nl
        << "==============================================================================" << nl
        << "                                Boundary field" << nl
        << "==============================================================================" << nl
        << endl;
        
    wordList patchNames(mesh.boundaryMesh().names());

    forAll(fld.boundaryField(), patchI)
    {
        const fvPatchField<Type>& pfld = fld.boundaryField()[patchI];
        const fvPatch& pp = mesh.boundary()[patchI];

        label nf = pp.start();            
        
        os<< "------------------------------------------------------------------------------" << nl
            << "                        " << patchNames[patchI] << nl
            << "------------------------------------------------------------------------------" << nl << nl
            << "------------------------------------------------------------------------------" << nl
            << "face index    flatFld index          field value          flatFld value" << nl
            << "------------------------------------------------------------------------------" << endl;

        forAll(pfld, i)
        {
            os<< "    " << nf << "            " << nFlat << "                    " << pfld[i] << "                " << flatFld[nFlat] << endl;
            nf++;
            nFlat++;
        }
        os<< endl;
    }
    os<< nl
        << "==============================================================================" << nl
        << "                                Non-local field data" << nl
        << "==============================================================================" << nl
        << nl
        << "------------------------------------------------------------------------------" << nl
        << "flatFld index            flatFld value" << nl
        << "------------------------------------------------------------------------------" << endl;
        
    
    for(label i=nFlat; i<flatFld.size(); i++)
    {
        os<< "    " << i << "                        " << flatFld[i] << endl;
    }
    
    os<< nl  
        << "//==========================================================================\\" << nl << nl
        << "                        End of " << fld.name() << " flat field" << nl << nl
        << "\\==========================================================================//" << nl
        << endl;
}

template<class Type>
void write_surfaceField
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sfld,
    const fvMesh& mesh
)
{
    fileName outputFile(sfld.name()+"_db");
    OFstream os(mesh.time().path()/mesh.time().timeName()/outputFile);

    os<< "//==========================================================================\\" << endl << endl
      << "                           " << sfld.name() << " field" << endl << endl
      << "\\==========================================================================//" << endl
      << endl       
      << "==============================================================================" << endl  
      << "                                Internal field" << endl
      << "==============================================================================" << endl
      << endl
      << "------------------------------------------------------------------------------" << endl
      << "    Point index                        Field value" << endl
      << "------------------------------------------------------------------------------" << endl;

    forAll(sfld.internalField(), pointI)
    {
        os<< "         " << pointI << "                           " << sfld.internalField()[pointI] << endl;
    }

    os<< endl << endl
      << "==============================================================================" << endl
      << "                                Boundary field" << endl
      << "==============================================================================" << endl
      << endl;
        
    wordList patchNames(mesh.boundaryMesh().names());

    forAll(sfld.boundaryField(), patchI)
    {
        fvsPatchField<Type> psfld(sfld.boundaryField()[patchI]);    
        
        os<< "------------------------------------------------------------------------------" << endl
          << "                       " << patchNames[patchI] << endl
          << "------------------------------------------------------------------------------" << endl << endl
          << "------------------------------------------------------------------------------" << endl
          << "Patch face index                Field value" << endl
          << "------------------------------------------------------------------------------" << endl;

        forAll(psfld, i)
        {
            os<< "         " << i << "                         " << psfld[i] << endl;
        }
        os<< endl;
    }
    
    os<< "//==========================================================================\\" << endl << endl
      << "                      End of " << sfld.name() << " field" << endl << endl
      << "\\==========================================================================//" << endl
      << endl;    
}


template<class Type>
void print_surfaceField
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sfld,
    const fvMesh& mesh,
    OFstream& os
)
{    
    os<< "//==========================================================================\\" << endl << endl
      << "                           " << sfld.name() << " field" << endl << endl
      << "\\==========================================================================//" << endl
      << endl       
      << "==============================================================================" << endl  
      << "                                Internal field" << endl
      << "==============================================================================" << endl
      << endl
      << "------------------------------------------------------------------------------" << endl
      << "    Point index                        Field value" << endl
      << "------------------------------------------------------------------------------" << endl;

    forAll(sfld.internalField(), pointI)
    {
        os<< "         " << pointI << "                           " << sfld.internalField()[pointI] << endl;
    }

    os<< endl << endl
      << "==============================================================================" << endl
      << "                                Boundary field" << endl
      << "==============================================================================" << endl
      << endl;
        
    wordList patchNames(mesh.boundaryMesh().names());

    forAll(sfld.boundaryField(), patchI)
    {
        fvsPatchField<Type> psfld(sfld.boundaryField()[patchI]);    
        
        os<< "------------------------------------------------------------------------------" << endl
          << "                       " << patchNames[patchI] << endl
          << "------------------------------------------------------------------------------" << endl << endl
          << "------------------------------------------------------------------------------" << endl
          << "Patch face index                Field value" << endl
          << "------------------------------------------------------------------------------" << endl;

        forAll(psfld, i)
        {
            os<< "         " << i << "                         " << psfld[i] << endl;
        }
        os<< endl;
    }
    
    os<< "//==========================================================================\\" << endl << endl
      << "                      End of " << sfld.name() << " field" << endl << endl
      << "\\==========================================================================//" << endl
      << endl;    
}


template<class Type>
void display_surfaceField
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sfld,
    const fvMesh& mesh
)
{   
    Info<< "//==========================================================================\\" << endl << endl
        << "                           " << sfld.name() << " field" << endl << endl
        << "\\==========================================================================//" << endl
        << endl       
        << "==============================================================================" << endl  
        << "                                Internal field" << endl
        << "==============================================================================" << endl
        << endl
        << "------------------------------------------------------------------------------" << endl
        << "    Point index                        Field value" << endl
        << "------------------------------------------------------------------------------" << endl;

    forAll(sfld.internalField(), pointI)
    {
        Info<< "         " << pointI << "                           " << sfld.internalField()[pointI] << endl;
    }

    Info<< endl << endl
        << "==============================================================================" << endl
        << "                                Boundary field" << endl
        << "==============================================================================" << endl
        << endl;
        
    wordList patchNames(mesh.boundaryMesh().names());

    forAll(sfld.boundaryField(), patchI)
    {
        fvsPatchField<Type> psfld(sfld.boundaryField()[patchI]);    
        
        Info<< "------------------------------------------------------------------------------" << endl
            << "                       " << patchNames[patchI] << endl
            << "------------------------------------------------------------------------------" << endl << endl
            << "------------------------------------------------------------------------------" << endl
            << "Patch face index                Field value" << endl
            << "------------------------------------------------------------------------------" << endl;

        forAll(psfld, i)
        {
            Info<< "         " << i << "                         " << psfld[i] << endl;
        }
        Info<< endl;
    }
    
    Info<< "//==========================================================================\\" << endl << endl
        << "                      End of " << sfld.name() << " field" << endl << endl
        << "\\==========================================================================//" << endl
        << endl;    
}


point centre
(
    const face& curFace,
    const pointField& points
)
{
    bool debug_ = false;

    if(debug_)
    {
        Info<< "Calculating face centre..." << endl;
    }

    // Calculate the centre by breaking the face into triangles and                                                            
    // area-weighted averaging their centres                                                                                   
    const label nPoints = curFace.size();
    // If the face is a triangle, do a direct calculation                                                                      
    if (nPoints == 3)
    {
        return
            (1.0/3.0)
            *(
                points[curFace[0]]
                + points[curFace[1]]
                + points[curFace[2]]
            );
    }

    point centrePoint = point::zero;
    for(label pI=0; pI<nPoints; pI++)
    {
        centrePoint += points[curFace[pI]];
    }
    centrePoint /= nPoints;

    if(debug_)
    {
        Info<< "Rough centre point:  " << centrePoint << endl;
    }
    scalar sumA = 0;
    vector sumAc = vector::zero;

    for(label pI=0; pI<nPoints; pI++)
    {
        const point& nextPoint = points[curFace[(pI + 1) % nPoints]];
        // Calculate 3*triangle centre                                                                                         
        const vector ttc
            (
                points[curFace[pI]]
                + nextPoint
                + centrePoint
            );
        // Calculate 2*triangle area                                                                                           
        const scalar ta = mag
            (
                (points[curFace[pI]] - centrePoint)
                ^ (nextPoint - centrePoint)
            );

        sumA += ta;
        sumAc += ta*ttc;
    }

    if (sumA > VSMALL)
    {
        return sumAc/(3.0*sumA);
    }
    else
    {
        return centrePoint;
    }
}


label findCellInFaceDir
(
    const labelList& cells,
    const List<vector>& C,
    const vector& Cp,
    const vector& nf,
    const label& C1_lbl,
    bool& foundCell,
    bool debug,
    OFstream& os
)
{
    foundCell = false;
    scalar cosThetaMax = 0;
    label cell_lbl = 0;

    if(debug)
    {
        os<< "Finding cell in face normal direction" << nl 
          << "nf: " << nf
          << endl;
    }

    for(label cellI=0; cellI<cells.size(); cellI++)
    {
        label curCell = cells[cellI];
        vector Ci = C[curCell];
        vector CCi = Ci - Cp;

        if(debug)
        {
            os<< "Cell " << cellI << " in stencil: " << curCell << nl
              << "C = " << Cp << "  Ci = " << Ci << nl
              << "CCi" << CCi
              << endl;
        }

        if(curCell != C1_lbl)
        {            
            CCi /= mag(CCi);
            scalar cosTheta = nf & CCi;

            if(debug)
            {
                os<< "costheta = " << cosTheta << endl;
            }

            if(cosTheta > cosThetaMax)
            {
                cell_lbl = curCell;
                cosThetaMax = cosTheta;
                foundCell = true;
            }
        }
    }

    if(debug)
    {
        os<< "Cell in nf direction: " << cell_lbl << nl
          << endl;
    }

    return cell_lbl;
}


label findCellInFaceOrthDir
(
    const labelList& cells,
    const List<vector>& C,
    const vector& Cp,
    const vector& C1,
    const vector& nf,
    const label& C1_lbl,
    bool& foundCell,
    bool debug,
    OFstream& os
)
{
    foundCell = false;
    scalar magCosThetaMin = 1;
    label cell_lbl = 0;

    if(debug)
    {
        os<< "Finding cell in face normal orthogonal direction" << nl 
            << "nf: " << nf
            << endl;
    }

    scalar CC1_CCi_dir;
    //vector Ci;
    vector CC1 = C1 - Cp;
    vector nfXCC1 = nf ^ CC1;
    nfXCC1 /= mag(CC1);
    vector CCi;

    // first try to find cell C2 so that vector nf lies between 
    // vector CC1 and CC2 ==> CC1_CCi_dir < 0
    for(label cellI=0; cellI<cells.size(); cellI++)
    {
        label curCell = cells[cellI];
        CCi = C[curCell] - Cp;

        CC1_CCi_dir = nfXCC1 & (nf ^ CCi);
        CC1_CCi_dir /= mag(CCi);

        if(debug)
        {
            os<< "Cell " << cellI << " in stencil: " << curCell << nl
                << "C = " << Cp << "  C1 = " << C1 << "  Ci = " << C[curCell] << nl
                << "CC1 = " << CC1 << "  CCi" << CCi << "  CC1_CCi_dir = " << CC1_CCi_dir
                << endl;
        }

        if(curCell != C1_lbl && CC1_CCi_dir < 0)
        {
            //CCi /= mag(CCi);
            //scalar cosTheta = nf & CCi;
            scalar cosTheta = (CC1 & CCi)/mag(CCi)/mag(CC1);

            if(debug)
            {
                os<< "costheta = " << cosTheta << endl;
            }

            if(mag(cosTheta) < magCosThetaMin)
            {
                cell_lbl = curCell;
                magCosThetaMin = mag(cosTheta);
                foundCell = true;
            }
        }
    }

    // if can't find cell C2 so that vector nf lies between 
    // vector CC1 and CC2 ==> CC1_CCi_dir < 0, find cell
    // closest ortogonal to CC1 direction and use backward difference
    // in calcCellGradWeights()

    if(!foundCell)
    {
        for(label cellI=0; cellI<cells.size(); cellI++)
        {
            label curCell = cells[cellI];
            CCi = C[curCell] - Cp;            

            if(debug)
            {
                os<< "Cell " << cellI << " in stencil: " << curCell << nl
                    << "C = " << Cp << "  C1 = " << C1 << "  Ci = " << C[curCell] << nl                    
                    << endl;
            }

            if(curCell != C1_lbl)
            {
                //CCi /= mag(CCi);
                //scalar cosTheta = nf & CCi;
                scalar cosTheta = (CC1 & CCi)/mag(CCi)/mag(CC1);

                if(debug)
                {
                    os<< "costheta = " << cosTheta << endl;
                }

                if(mag(cosTheta) < magCosThetaMin)
                {
                    cell_lbl = curCell;
                    magCosThetaMin = mag(cosTheta);
                    foundCell = true;
                }
            }
        }
    }

    if(debug)
    {
        os<< "Cell in nf orthogonal direction: " << cell_lbl << nl
            << endl;
    }

    return cell_lbl;
}


void calcCellGradWeights
(
    const label& curCell_lbl,
    const vector& nf,
    const List<scalar>& Y,
    const List<scalar>& alpha1,
    const List<vector>& C,
    const labelList& curCellsAll,
    const scalar& MIN_ALPHA_DIFF,
    const label& phaseLbl,
    scalar& alpha,
    scalar& beta,
    scalar& magt1,
    scalar& magt2,
    scalar& d,
    bool debug, 
    OFstream& os
)
{        
    if(debug)
    {
        os<< "Calculating cell grad weights in cell " << curCell_lbl << nl
            << endl;
    }

    scalar MAX_ALPHA_DIFF = 1 - MIN_ALPHA_DIFF;
    
    labelList curCells(curCellsAll.size());
    label n_ph = 0;

    if(debug)
    {
        os<< "Reducing cell stencil for phase " << phaseLbl << nl
            << "Full cell stencil" << nl
            << curCellsAll << nl
            << endl;
    }

    if(phaseLbl == 1)
    {
        for(label cellI=0; cellI<curCellsAll.size(); cellI++)
        {            
            label cellI_lbl = curCellsAll[cellI];

            if(debug)
            {
                os<< "cellI: " << cellI << "  cellI_lbl: " << cellI_lbl << nl
                    << "alpha1 = " << alpha1[cellI_lbl] << "  MIN_ALPHA_DIFF = " << MIN_ALPHA_DIFF << nl
                    << endl;
            }

            if(alpha1[cellI_lbl] > MIN_ALPHA_DIFF && cellI_lbl != curCell_lbl)
            {
                curCells[n_ph++] = cellI_lbl;
            }            
        }        
    }
    else
    {
        for(label cellI=0; cellI<curCellsAll.size(); cellI++)
        {
            label cellI_lbl = curCellsAll[cellI];

            if(debug)
            {
                os<< "cellI: " << cellI << "  cellI_lbl: " << cellI_lbl << nl
                    << "alpha1 = " << alpha1[cellI_lbl] << "  MAX_ALPHA_DIFF = " << MAX_ALPHA_DIFF << nl
                    << endl;
            }

            if(alpha1[cellI_lbl] < MAX_ALPHA_DIFF && cellI_lbl != curCell_lbl)
            {
                curCells[n_ph++] = curCellsAll[cellI];
            }            
        }
    }
    
    curCells.setSize(n_ph);

    if(debug)
    {
        os<< "Cell reduced stencil" << nl
            << curCells << nl
            << endl;
    }

    // suffix 1: direction closest to nf
    // suffix 2: direction closest orthogonal to nf in 2-D
    // further improvements needed for 3-D calculation
    vector Cp = C[curCell_lbl];
    bool foundCell1 = false;
    label C1_lbl = findCellInFaceDir(curCells,C,Cp,nf,-1,foundCell1,debug,os);
    if(foundCell1)
    {
        vector C1 = C[C1_lbl];
        scalar Yp = Y[curCell_lbl];
        scalar Y1 = Y[C1_lbl];

        vector t1 = C1 - Cp;
        magt1 = mag(t1);
        if(debug)
        {
            os<< "Cell: " << curCell_lbl << "  Cp = " << Cp << "  nf = " << nf << nl
                << "C1_lbl: " << C1_lbl << "  C1 = " << C1 << " t1 = " << t1 << "  mag(t1) = " << magt1
                << endl;
        }
        if(magt1 < SMALL){magt1 += SMALL;}
        scalar costheta1 = (nf & t1)/magt1;
        if(mag(costheta1) > 1){costheta1 = 1;}
        scalar theta1 = acos(costheta1);

        if(debug)
        {
            os<< "Cell: " << curCell_lbl << "  Cp = " << Cp << "  nf = " << nf << nl
                << "C1_lbl: " << C1_lbl << "  C1 = " << C1 << " t1 = " << t1 << "  mag(t1) = " << magt1
                << "costheta1 = " << costheta1 << "  theta1 = " << theta1
                << endl;
        }

        if(theta1 > 1E-3)
        {
            bool foundCell2 = false;
            label C2_lbl = findCellInFaceOrthDir(curCells,C,Cp,C1,nf,C1_lbl,foundCell2,debug,os);

            if(foundCell2)
            {
                vector C2 = C[C2_lbl];
                scalar Y2 = Y[C2_lbl];
                scalar Y2_1 = Y2;
                vector t2 = C2 - Cp;
                magt2 = mag(t2);
                if(magt2 < SMALL){magt2 += SMALL;}    
                scalar magt1t2 = magt1*magt2;
                if(magt1t2 < SMALL){magt1t2 += SMALL;}

                if(debug)
                {
                    os<< "C2_lbl: " << C2_lbl << "  C2 = " << C2 << " t1 = " << t2 << "  mag(t2) = " << magt2 << "  mag(t1*t2) = " << magt1t2
                        << endl;
                }

                scalar theta2_sign = -((nf ^ t1) & (nf ^ t2))/magt1t2;                

                if(debug)
                {
                    os<< "mag(t1*t2) = " << magt1t2 << "  theta2_sign = " << theta2_sign;
                }

                scalar theta2;
                scalar theta2_1;
                scalar costheta2 = (nf & t2)/magt2;
                if(debug)
                {
                    os<< "  costheta2 = " << costheta2; 
                }
                if(mag(costheta2) > 1){costheta2 /= mag(costheta2);}
                if(theta2_sign >= 0)
                {
                    theta2 = acos(costheta2);
                    theta2_1 = theta2;
                }
                else
                {
                    theta2 = constant::mathematical::pi - acos(costheta2);
                    theta2_1 = acos(costheta2);
                    Y2 = Yp + (Yp - Y2_1);
                }

                if(debug)
                {
                    os
                        << "theta2 = " << theta2 << "  theta2_1 = " << theta2_1 << nl                    
                            << "Yp = " << Yp << "  Y1 = " << Y1 << "  Y2 = " << Y2 << "  Y2_1 = " << Y2_1 
                            << endl;
                }

                scalar sintheta12 = sin(theta1 + theta2);
                if(sintheta12 < SMALL){sintheta12 += SMALL;}

                alpha = sin(theta2)/sintheta12;
                beta = sin(theta1)/sintheta12;    
                d = beta*Y2/magt2;
            }
            else
            {
                magt2 = mag(C[curCellsAll[0]] - C[curCellsAll[1]]);
                alpha = 1;
                beta = 0;
                d = 1.0E6/magt2;
            }
        }
        else
        {
            magt2 = mag(C[curCellsAll[0]] - C[curCellsAll[1]]);
            alpha = 1;
            beta = 0;
            d = 0;
        }
    }
    else
    {
        os<< "Cell 1 not found! Fatal error!" << nl
            << endl;

        magt1 = mag(C[curCellsAll[0]] - C[curCellsAll[1]]);
        magt2 = mag(C[curCellsAll[0]] - C[curCellsAll[1]]);
        alpha = 1;
        beta = 0;
        d = 1.0E6/magt2;
    }

    if(debug)
    {
        os<< "alpha = " << alpha << "  beta = " << beta << "  d = " << d << nl
            << endl;
    }
}


void calcCellGradWeights
(
    const label& curCell_lbl,
    const vector& nf,
    const List<List<scalar> >& Y,
    const List<scalar>& T,
    const List<scalar>& alpha1,
    const List<vector>& C,
    const labelList& curCellsAll,
    const scalar& MIN_ALPHA_DIFF,
    const label& phaseLbl,
    scalar& alpha,
    scalar& beta,
    scalar& magt1,
    scalar& magt2,
    List<scalar>& Yp,
    List<scalar>& Y1,
    List<scalar>& Y2,
    scalar& Tp,
    scalar& T1,
    scalar& T2,
    bool& foundCell1,
    bool& foundCell2,
    const label& n,
    bool debug, 
    OFstream& os
)
{
    label i, n_ph, cellI, cellI_lbl, C1_lbl, C2_lbl;    
    scalar MAX_ALPHA_DIFF, costheta1, theta1, magt1t2, theta2_sign, theta2, theta2_1, costheta2;
    scalar sintheta12, T2_1;
    vector Cp, C1, C2, t1, t2;
    labelList curCells(curCellsAll.size());    
    List<scalar> Y2_1(n);    

    if(debug)
    {
        os<< "Calculating cell grad weights in cell " << curCell_lbl << nl
            << endl;
    }

    MAX_ALPHA_DIFF = 1 - MIN_ALPHA_DIFF;        
    n_ph = 0;

    if(debug)
    {
        os<< "Reducing cell stencil for phase " << phaseLbl << nl
            << "Full cell stencil" << nl
            << curCellsAll << nl
            << endl;
    }

    if(phaseLbl == 1)
    {
        for(cellI=0; cellI<curCellsAll.size(); cellI++)
        {            
            cellI_lbl = curCellsAll[cellI];

            if(debug)
            {
                os<< "cellI: " << cellI << "  cellI_lbl: " << cellI_lbl << nl
                    << "alpha1 = " << alpha1[cellI_lbl] << "  MIN_ALPHA_DIFF = " << MIN_ALPHA_DIFF << nl
                    << endl;
            }

            if(alpha1[cellI_lbl] > MIN_ALPHA_DIFF && cellI_lbl != curCell_lbl)
            {
                curCells[n_ph++] = cellI_lbl;
            }            
        }        
    }
    else
    {
        for(cellI=0; cellI<curCellsAll.size(); cellI++)
        {
            cellI_lbl = curCellsAll[cellI];

            if(debug)
            {
                os<< "cellI: " << cellI << "  cellI_lbl: " << cellI_lbl << nl
                    << "alpha1 = " << alpha1[cellI_lbl] << "  MAX_ALPHA_DIFF = " << MAX_ALPHA_DIFF << nl
                    << endl;
            }

            if(alpha1[cellI_lbl] < MAX_ALPHA_DIFF && cellI_lbl != curCell_lbl)
            {
                curCells[n_ph++] = curCellsAll[cellI];
            }            
        }
    }
    
    curCells.setSize(n_ph);

    if(debug)
    {
        os<< "Cell reduced stencil" << nl
            << curCells << nl
            << endl;
    }

    // suffix 1: direction closest to nf
    // suffix 2: direction closest orthogonal to nf in 2-D
    // further improvements needed for 3-D calculation
    Cp = C[curCell_lbl];
    Tp = T[curCell_lbl];
    for(i=0; i<n; i++)
    {
        Yp[i] = Y[i][curCell_lbl];
    }
    foundCell1 = false;
    C1_lbl = findCellInFaceDir(curCells,C,Cp,nf,-1,foundCell1,debug,os);
    if(foundCell1)
    {
        C1 = C[C1_lbl];
        for(i=0; i<n; i++)
        {
            Y1[i] = Y[i][C1_lbl];
        }        
        T1 = T[C1_lbl];

        t1 = C1 - Cp;
        magt1 = mag(t1);
        if(debug)
        {
            os<< "Cell: " << curCell_lbl << "  Cp = " << Cp << "  nf = " << nf << nl
                << "C1_lbl: " << C1_lbl << "  C1 = " << C1 << " t1 = " << t1 << "  mag(t1) = " << magt1
                << endl;
        }
        if(magt1 < SMALL){magt1 += SMALL;}
        costheta1 = (nf & t1)/magt1;
        if(mag(costheta1) > 1){costheta1 = 1;}
        theta1 = acos(costheta1);

        if(debug)
        {
            os<< "Cell: " << curCell_lbl << "  Cp = " << Cp << "  nf = " << nf << nl
                << "C1_lbl: " << C1_lbl << "  C1 = " << C1 << " t1 = " << t1 << "  mag(t1) = " << magt1
                << "costheta1 = " << costheta1 << "  theta1 = " << theta1
                << endl;
        }

        if(theta1 > 1E-3)
        {
            foundCell2 = false;
            C2_lbl = findCellInFaceOrthDir(curCells,C,Cp,C1,nf,C1_lbl,foundCell2,debug,os);

            if(foundCell2)
            {
                C2 = C[C2_lbl];
                for(i=0; i<n; i++)
                {
                    Y2[i] = Y[i][C2_lbl];
                    Y2_1[i] = Y2[i];
                }
                T2 = T[C2_lbl];
                T2_1 = T2;
                t2 = C2 - Cp;
                magt2 = mag(t2);
                if(magt2 < SMALL){magt2 += SMALL;}    
                magt1t2 = magt1*magt2;
                if(magt1t2 < SMALL){magt1t2 += SMALL;}

                if(debug)
                {
                    os<< "C2_lbl: " << C2_lbl << "  C2 = " << C2 << " t1 = " << t2 << "  mag(t2) = " << magt2 << "  mag(t1*t2) = " << magt1t2
                        << endl;
                }

                theta2_sign = -((nf ^ t1) & (nf ^ t2))/magt1t2;                

                if(debug)
                {
                    os<< "mag(t1*t2) = " << magt1t2 << "  theta2_sign = " << theta2_sign;
                }
                
                costheta2 = (nf & t2)/magt2;
                if(debug)
                {
                    os<< "  costheta2 = " << costheta2; 
                }
                if(mag(costheta2) > 1){costheta2 /= mag(costheta2);}
                if(theta2_sign >= 0)
                {
                    theta2 = acos(costheta2);
                    theta2_1 = theta2;
                }
                else
                {
                    theta2 = constant::mathematical::pi - acos(costheta2);
                    theta2_1 = acos(costheta2);
                    for(i=0; i<n; i++)
                    {
                        Y2[i] = Yp[i] + (Yp[i] - Y2_1[i]);
                    }
                    T2 = Tp + (Tp - T2_1);
                }

                if(debug)
                {
                    os<< "theta2 = " << theta2 << "  theta2_1 = " << theta2_1 << nl                    
                        << setw(7) << "Species" << "  " << setw(10) << "Yp" << "  " << setw(10) << "Y1" << "  " << setw(10) << "Y2" << "  " <<setw(10) << "Y2_1" << endl;
                    for(i=0; i<n; i++)
                    {
                        os<< setw(7) << i << "  " << setw(10) << Yp[i] << "  " << setw(10) << Y1[i] << "  " << setw(10) << Y2[i] << "  " <<setw(10) << Y2_1[i] << endl;
                    }
                    os<< "Tp = " << Tp << "  T1 = " << T1 << "  T2 = " << T2 << "  T2_1 = " << T2_1 << endl;
                }

                sintheta12 = sin(theta1 + theta2);
                if(sintheta12 < SMALL){sintheta12 += SMALL;}

                alpha = sin(theta2)/sintheta12;
                beta = sin(theta1)/sintheta12;                
            }
            else
            {
                magt2 = mag(C[curCellsAll[0]] - C[curCellsAll[1]]);
                alpha = 1;
                beta = 0;
                for(i=0; i<n; i++)
                {
                    Y2[i] = Y1[i];
                }
                T2 = T1;
            }            
        }
        else
        {
            foundCell2 = false;
            magt2 = mag(C[curCellsAll[0]] - C[curCellsAll[1]]);
            alpha = 1;
            beta = 0;
            for(i=0; i<n; i++)
            {
                Y2[i] = Y1[i];
            }
            T2 = T1;
        }
    }
    else
    {
        os<< "Cell 1 not found! Fatal error!" << nl
            << endl;
        foundCell1 = false;
        foundCell2 = false;
        magt1 = mag(C[curCellsAll[0]] - C[curCellsAll[1]]);
        magt2 = mag(C[curCellsAll[0]] - C[curCellsAll[1]]);
        alpha = 0;
        beta = 0;
        for(i=0; i<n; i++)
        {
            Y1[i] = Yp[i];
            Y2[i] = Yp[i];            
        }
        T1 = Tp;
        T2 = Tp;
    }

    if(debug)
    {
        os<< "alpha = " << alpha << "  beta = " << beta << nl
            << endl;
    }
}


void calcCellGrad
(
    const label& curCell_lbl,
    const vector& nf,
    const List<scalar>& Y,
    const List<scalar>& alpha1,
    const List<vector>& C,
    const labelList& curCellsAll,
    const scalar& MIN_ALPHA_DIFF,
    const label& phaseLbl,
    scalar& cellGrad,
    bool debug, 
    OFstream& os
)
{    
    scalar alpha; scalar beta; scalar magt1; scalar magt2;

    if(debug)
    {
        os<< "Calculating cell grad in cell " << curCell_lbl << nl
            << endl;
    }
    
    scalar MAX_ALPHA_DIFF = 1 - MIN_ALPHA_DIFF;
    
    labelList curCells(curCellsAll.size());
    label n_ph = 0;

    if(phaseLbl == 1)
    {
        for(label cellI=0; cellI<curCellsAll.size(); cellI++)
        {
            label cellI_lbl = curCellsAll[cellI];
            if(alpha1[cellI_lbl] > MIN_ALPHA_DIFF && cellI_lbl != curCell_lbl)
            {
                curCells[n_ph++] = cellI_lbl;
            }            
        }        
    }
    else
    {
        for(label cellI=0; cellI<curCellsAll.size(); cellI++)
        {
            label cellI_lbl = curCellsAll[cellI];
            if(alpha1[cellI_lbl] < MAX_ALPHA_DIFF && cellI_lbl != curCell_lbl)
            {
                curCells[n_ph++] = curCellsAll[cellI];
            }            
        }
    }
    
    curCells.setSize(n_ph);

    if(debug)
    {
        os<< "Cell reduced stencil" << nl
            << curCells << endl;    
    }

    // suffix 1: direction closest to nf
    // suffix 2: direction closest orthogonal to nf in 2-D
    // further improvements needed for 3-D calculation
     vector Cp = C[curCell_lbl];
    bool foundCell1 = false;
    label C1_lbl = findCellInFaceDir(curCells,C,Cp,nf,-1,foundCell1,debug,os);
    if(foundCell1)
    {
        vector C1 = C[C1_lbl];
        scalar Yp = Y[curCell_lbl];
        scalar Y1 = Y[C1_lbl];

        vector t1 = C1 - Cp;
        magt1 = mag(t1);
        if(debug)
        {
            os<< "Cell: " << curCell_lbl << "  Cp = " << Cp << "  nf = " << nf << nl
                << "C1_lbl: " << C1_lbl << "  C1 = " << C1 << " t1 = " << t1 << "  mag(t1) = " << magt1
                << endl;
        }
        if(magt1 < SMALL){magt1 += SMALL;}
        scalar costheta1 = (nf & t1)/magt1;
        if(mag(costheta1) > 1){costheta1 = 1;}
        scalar theta1 = acos(costheta1);

        if(debug)
        {
            os<< "Cell: " << curCell_lbl << "  Cp = " << Cp << "  nf = " << nf << nl
                << "C1_lbl: " << C1_lbl << "  C1 = " << C1 << " t1 = " << t1 << "  mag(t1) = " << magt1
                << "costheta1 = " << costheta1 << "  theta1 = " << theta1
                << endl;
        }

        if(theta1 > 1E-3)
        {
            bool foundCell2 = false;
            label C2_lbl = findCellInFaceOrthDir(curCells,C,Cp,C1,nf,C1_lbl,foundCell2,debug,os);

            if(foundCell2)
            {
                vector C2 = C[C2_lbl];
                scalar Y2 = Y[C2_lbl];
                scalar Y2_1 = Y2;
                vector t2 = C2 - Cp;
                magt2 = mag(t2);
                if(magt2 < SMALL){magt2 += SMALL;}    
                scalar magt1t2 = magt1*magt2;
                if(magt1t2 < SMALL){magt1t2 += SMALL;}

                if(debug)
                {
                    os<< "C2_lbl: " << C2_lbl << "  C2 = " << C2 << " t1 = " << t2 << "  mag(t2) = " << magt2 << "  mag(t1*t2) = " << magt1t2
                        << endl;
                }

                scalar theta2_sign = -((nf ^ t1) & (nf ^ t2))/magt1t2;                

                if(debug)
                {
                    os<< "mag(t1*t2) = " << magt1t2 << "  theta2_sign = " << theta2_sign;
                }

                scalar theta2;
                scalar theta2_1;
                scalar costheta2 = (nf & t2)/magt2;
                if(debug)
                {
                    os<< "  costheta2 = " << costheta2; 
                }
                if(mag(costheta2) > 1){costheta2 /= mag(costheta2);}
                if(theta2_sign >= 0)
                {
                    theta2 = acos(costheta2);
                    theta2_1 = theta2;
                }
                else
                {
                    theta2 = constant::mathematical::pi - acos(costheta2);
                    theta2_1 = acos(costheta2);
                    Y2 = Yp + (Yp - Y2_1);
                }

                if(debug)
                {
                    os
                        << "theta2 = " << theta2 << "  theta2_1 = " << theta2_1 << nl                    
                            << "Yp = " << Yp << "  Y1 = " << Y1 << "  Y2 = " << Y2 << "  Y2_1 = " << Y2_1 
                            << endl;
                }

                scalar sintheta12 = sin(theta1 + theta2);
                if(sintheta12 < SMALL){sintheta12 += SMALL;}

                alpha = sin(theta2)/sintheta12;
                beta = sin(theta1)/sintheta12;

                cellGrad = alpha*(Y1 - Yp)/magt1 + beta*(Y2 - Yp)/magt2;
            }
            else
            {
                alpha = 1;
                beta = 0;

                cellGrad = (Y1 - Yp)/magt1;
            }
        }
        else
        {
            alpha = 1;
            beta = 0;

            cellGrad = (Y1 - Yp)/magt1;
        }
    }
    else
    {
        os<< "Cell 1 not found! Fatal error!" << nl
            << endl;

        alpha = 1;
        beta = 0;

        cellGrad = 0;
    }    

    if(debug)
    {
        os<< "alpha = " << alpha << "  beta = " << beta << "  cell grad = " << cellGrad << nl
            << endl;
    }
}


void calcCellGrad
(
    const label& curCell_lbl,
    const vector& nf,
    const List<List<scalar> >& Y,
    const List<scalar>& T,
    const List<scalar>& alpha1,
    const List<vector>& C,
    const labelList& curCellsAll,
    const scalar& MIN_ALPHA_DIFF,
    const label& phaseLbl,
    List<scalar>& cellGrad_Y,
    scalar& cellGrad_T,
    const label& n,
    bool debug, 
    OFstream& os
)
{    
    label i, n_ph, cellI, cellI_lbl, C1_lbl, C2_lbl;
    bool foundCell1, foundCell2;
    scalar MAX_ALPHA_DIFF, costheta1, theta1, magt1t2, theta2_sign, theta2, theta2_1, costheta2;
    scalar sintheta12, alpha, beta, magt1, magt2, Tp, T1, T2, T2_1;
    vector Cp, C1, C2, t1, t2;
    labelList curCells(curCellsAll.size());
    List<scalar> Yp(n);
    List<scalar> Y1(n);
    List<scalar> Y2(n);
    List<scalar> Y2_1(n);

    if(debug)
    {
        os<< "Calculating cell grad weights in cell " << curCell_lbl << nl
            << endl;
    }

    MAX_ALPHA_DIFF = 1 - MIN_ALPHA_DIFF;        
    n_ph = 0;

    if(debug)
    {
        os<< "Reducing cell stencil for phase " << phaseLbl << nl
            << "Full cell stencil" << nl
            << curCellsAll << nl
            << endl;
    }

    if(phaseLbl == 1)
    {
        for(cellI=0; cellI<curCellsAll.size(); cellI++)
        {            
            cellI_lbl = curCellsAll[cellI];

            if(debug)
            {
                os<< "cellI: " << cellI << "  cellI_lbl: " << cellI_lbl << nl
                    << "alpha1 = " << alpha1[cellI_lbl] << "  MIN_ALPHA_DIFF = " << MIN_ALPHA_DIFF << nl
                    << endl;
            }

            if(alpha1[cellI_lbl] > MIN_ALPHA_DIFF && cellI_lbl != curCell_lbl)
            {
                curCells[n_ph++] = cellI_lbl;
            }            
        }        
    }
    else
    {
        for(cellI=0; cellI<curCellsAll.size(); cellI++)
        {
            cellI_lbl = curCellsAll[cellI];

            if(debug)
            {
                os<< "cellI: " << cellI << "  cellI_lbl: " << cellI_lbl << nl
                    << "alpha1 = " << alpha1[cellI_lbl] << "  MAX_ALPHA_DIFF = " << MAX_ALPHA_DIFF << nl
                    << endl;
            }

            if(alpha1[cellI_lbl] < MAX_ALPHA_DIFF && cellI_lbl != curCell_lbl)
            {
                curCells[n_ph++] = curCellsAll[cellI];
            }            
        }
    }
    
    curCells.setSize(n_ph);

    if(debug)
    {
        os<< "Cell reduced stencil" << nl
            << curCells << nl
            << endl;
    }

    // suffix 1: direction closest to nf
    // suffix 2: direction closest orthogonal to nf in 2-D
    // further improvements needed for 3-D calculation
    Cp = C[curCell_lbl];
    for(i=0; i<n; i++)
    {
        Yp[i] = Y[i][curCell_lbl];
    }
    Tp = T[curCell_lbl];
    foundCell1 = false;
    C1_lbl = findCellInFaceDir(curCells,C,Cp,nf,-1,foundCell1,debug,os);
    if(foundCell1)
    {
        C1 = C[C1_lbl];
        for(i=0; i<n; i++)
        {
            Y1[i] = Y[i][C1_lbl];
        }
        T1 = T[C1_lbl];

        t1 = C1 - Cp;
        magt1 = mag(t1);
        if(debug)
        {
            os<< "Cell: " << curCell_lbl << "  Cp = " << Cp << "  nf = " << nf << nl
                << "C1_lbl: " << C1_lbl << "  C1 = " << C1 << " t1 = " << t1 << "  mag(t1) = " << magt1
                << endl;
        }
        if(magt1 < SMALL){magt1 += SMALL;}
        costheta1 = (nf & t1)/magt1;
        if(mag(costheta1) > 1){costheta1 = 1;}
        theta1 = acos(costheta1);

        if(debug)
        {
            os<< "Cell: " << curCell_lbl << "  Cp = " << Cp << "  nf = " << nf << nl
                << "C1_lbl: " << C1_lbl << "  C1 = " << C1 << " t1 = " << t1 << "  mag(t1) = " << magt1
                << "costheta1 = " << costheta1 << "  theta1 = " << theta1
                << endl;
        }

        if(theta1 > 1E-3)
        {
            foundCell2 = false;
            C2_lbl = findCellInFaceOrthDir(curCells,C,Cp,C1,nf,C1_lbl,foundCell2,debug,os);

            if(foundCell2)
            {
                C2 = C[C2_lbl];
                for(i=0; i<n; i++)
                {
                    Y2[i] = Y[i][C2_lbl];
                    Y2_1[i] = Y2[i];
                }
                T2 = T[C2_lbl];
                T2_1 = T2;
                t2 = C2 - Cp;
                magt2 = mag(t2);
                if(magt2 < SMALL){magt2 += SMALL;}    
                magt1t2 = magt1*magt2;
                if(magt1t2 < SMALL){magt1t2 += SMALL;}

                if(debug)
                {
                    os<< "C2_lbl: " << C2_lbl << "  C2 = " << C2 << " t1 = " << t2 << "  mag(t2) = " << magt2 << "  mag(t1*t2) = " << magt1t2
                        << endl;
                }

                theta2_sign = -((nf ^ t1) & (nf ^ t2))/magt1t2;                

                if(debug)
                {
                    os<< "mag(t1*t2) = " << magt1t2 << "  theta2_sign = " << theta2_sign;
                }
                
                costheta2 = (nf & t2)/magt2;
                if(debug)
                {
                    os<< "  costheta2 = " << costheta2; 
                }
                if(mag(costheta2) > 1){costheta2 /= mag(costheta2);}
                if(theta2_sign >= 0)
                {
                    theta2 = acos(costheta2);
                    theta2_1 = theta2;
                }
                else
                {
                    theta2 = constant::mathematical::pi - acos(costheta2);
                    theta2_1 = acos(costheta2);
                    for(i=0; i<n; i++)
                    {
                        Y2[i] = Yp[i] + (Yp[i] - Y2_1[i]);
                    }
                    T2 = Tp + (Tp - T2_1);
                }

                if(debug)
                {
                    os<< "theta2 = " << theta2 << "  theta2_1 = " << theta2_1 << nl                    
                        << setw(7) << "Species" << "  " << setw(10) << "Yp" << "  " << setw(10) << "Y1" << "  " << setw(10) << "Y2" << "  " <<setw(10) << "Y2_1" << endl;
                    for(i=0; i<n; i++)
                    {
                        os<< setw(7) << i << "  " << setw(10) << Yp[i] << "  " << setw(10) << Y1[i] << "  " << setw(10) << Y2[i] << "  " <<setw(10) << Y2_1[i] << endl;
                    }
                    os<< "Tp = " << Tp << "  T1 = " << T1 << "  T2 = " << T2 << "  T2_1 = " << T2_1 << endl;
                }

                sintheta12 = sin(theta1 + theta2);
                if(sintheta12 < SMALL){sintheta12 += SMALL;}

                alpha = sin(theta2)/sintheta12;
                beta = sin(theta1)/sintheta12;
                for(i=0; i<n; i++)
                {
                    cellGrad_Y[i] = alpha*(Y1[i] - Yp[i])/magt1 + beta*(Y2[i] - Yp[i])/magt2;
                }
                cellGrad_T = alpha*(T1 - Tp)/magt1 + beta*(T2 - Tp)/magt2;
            }
            else
            {                
                alpha = 1;
                beta = 0;
                for(i=0; i<n; i++)
                {
                    cellGrad_Y[i] = (Y1[i] - Yp[i])/magt1;
                }
                cellGrad_T = (T1 - Tp)/magt1;
            }
        }
        else
        {
            alpha = 1;
            beta = 0;
            for(i=0; i<n; i++)
            {
                cellGrad_Y[i] = (Y1[i] - Yp[i])/magt1;
            }
            cellGrad_T = (T1 - Tp)/magt1;
        }
    }
    else
    {
        os<< "Cell 1 not found! Fatal error!" << nl
            << endl;
        
        alpha = 0;
        beta = 0;
        for(i=0; i<n; i++)
        {
            cellGrad_Y[i] = 0;
        }
        cellGrad_T = 0;
    }

    if(debug)
    {
        os<< "alpha = " << alpha << "  beta = " << beta << nl 
            << setw(7) << "Species" << "  " << setw(10) << "cellGrad_Y" << endl;
        for(i=0; i<n; i++)
        {
            os<< setw(7) << i << "  " << setw(10) << cellGrad_Y[i] << endl;
        }
        os<< "cellGrad_T = " << cellGrad_T << endl;
        os<< endl;
    }
}


void calcTwoSidedFaceGradWeights
(
    const label& faceI,
    const label& own,
    const label& nei,
    const vector& nf,
    const List<scalar>& Y,
    const List<scalar>& alpha1,
    const List<vector>& C,
    const labelListList& diffCellStencil,
    const scalar& MIN_ALPHA_DIFF,
    const label& phaseLbl,
    scalar& wOwn,
    scalar& wNei,    
    bool debug, 
    OFstream& os
)
{
    if(debug)
    {
        os<< "Calculating two-sided grad weights for face " << faceI << "  Phase: " << phaseLbl << nl
            << "Own: " << own << "  Nei: " << nei << nl
            << endl;
    }

    if(debug)
    {
        os<< "Own side calculation" << nl 
            << endl;
    }
    scalar alphap = 1;
    scalar betap = 1;
    scalar magt1p = 1;
    scalar magt2p = 1;
    scalar dp = 1;
    const labelList& ownCells = diffCellStencil[own];
    calcCellGradWeights(own,nf,Y,alpha1,C,ownCells,MIN_ALPHA_DIFF,phaseLbl,alphap,betap,magt1p,magt2p,dp,debug,os);

    if(debug)
    {
        os<< "Nei side calculation" << nl
            << endl;
    }
    scalar alpham = 1;
    scalar betam = 1;
    scalar magt1m = 1;
    scalar magt2m = 1;
    scalar dm = 1;
    const labelList& neiCells = diffCellStencil[nei];
    calcCellGradWeights(nei,-nf,Y,alpha1,C,neiCells,MIN_ALPHA_DIFF,phaseLbl,alpham,betam,magt1m,magt2m,dm,debug,os);

    scalar mup;
    scalar mum;
    scalar dmdp;
    if(mag(dm) < SMALL && mag(dp) < SMALL)
    {
        mup = 0.5;
    }
    else
    {
        dmdp = dm + dp;
        if(mag(dmdp) < SMALL)
        {
            dmdp = mag(dmdp) + SMALL;
        }
        mup = dm/dmdp;
        mup = mag(mup);
        mup = min(mup, 1);
    }
    mum = 1 - mup;

    wOwn = -(mup*(alphap/magt1p + betap/magt2p) + mum*alpham/magt1m);
    wNei = mum*(alpham/magt1m + betam/magt2m) + mup*alphap/magt1p;

    if(debug)
    {
        os<< nl
            << "mup = " << mup << "  mum = " << mum << nl
            << "wOwn = " << wOwn << "  wNei = " << wNei << nl
            << endl;
    }
}


void calcFaceGradFromWeights
(
    const scalar& Yp,
    const scalar& Ym,
    const scalar& alphap,
    const scalar& betap,
    const scalar& magt1p,
    const scalar& magt2p,
    const scalar& dp,
    const scalar& alpham,
    const scalar& betam,
    const scalar& magt1m,
    const scalar& magt2m,
    const scalar& dm,    
    scalar& faceGrad, 
    bool debug, 
    OFstream& os
)
{
    if(debug)
    {
        os<< "Calculating face grad from weights" << nl            
            << endl;
    }

    scalar mup;
    scalar mum;
    scalar dmdp;
    if(mag(dm) < SMALL && mag(dp) < SMALL)
    {
        mup = 0.5;
    }
    else
    {
        dmdp = dm + dp;
        if(mag(dmdp) < SMALL)
        {
            dmdp = mag(dmdp) + SMALL;
        }
        mup = dm/dmdp;
        mup = mag(mup);
        mup = min(mup, 1);
    }
    
    mum = 1 - mup;

    scalar wOwn = -(mup*(alphap/magt1p + betap/magt2p) + mum*alpham/magt1m);
    scalar wNei = mum*(alpham/magt1m + betam/magt2m) + mup*alphap/magt1p;

    faceGrad = wOwn*Yp + wNei*Ym;

    if(debug)
    {
        os<< nl
            << "mup = " << mup << "  mum = " << mum << nl
            << "wOwn = " << wOwn << "  wNei = " << wNei << nl
            << "YOwn = " << Yp << "  YNei = " << Ym << nl
            << "face grad = " << faceGrad << nl
            << endl;
    }
}


void calcFaceGradFromWeights
(
    const vector& nf,
    const vector& Cf,
    const vector& Cp,    
    const scalar& alphap,
    const scalar& betap,
    const scalar& magt1p,
    const scalar& magt2p,
    const List<scalar>& Yp,
    const List<scalar>& Y1p,
    const List<scalar>& Y2p,
    const scalar& Tp,
    const scalar& T1p,
    const scalar& T2p,
    const bool& foundCell1p,
    const bool& foundCell2p,
    const vector& Cm,
    const scalar& alpham,
    const scalar& betam,
    const scalar& magt1m,
    const scalar& magt2m,
    const List<scalar>& Ym,
    const List<scalar>& Y1m,
    const List<scalar>& Y2m,
    const scalar& Tm,
    const scalar& T1m,
    const scalar& T2m,
    const bool& foundCell1m,
    const bool& foundCell2m,
    List<scalar>& faceGrad_Y,
    scalar& faceGrad_T,
    const label& n,
    bool debug, 
    OFstream& os
)
{
    label i;
    scalar mup, mum, dp, dm, dpdm, gradTOwn, gradTNei;
    List<scalar> gradYOwn(n);
    List<scalar> gradYNei(n);

    if(debug)
    {
        os<< "Calculating face grad from weights" << nl            
            << endl;
    }

    if(!foundCell1p && !foundCell1m)
    {
        mup = 0;
        mum = 0;
    }
    else
    {
        if(!foundCell1p)
        {
            mup = 0;
            mum = 1 - mup;
        }
        else
        {
            if(!foundCell1m)
            {
                mup = 1;
                mum = 1 - mup;
            }
            else
            {
                dp = mag(nf & (Cf - Cp));
                dm = mag(nf & (Cf - Cm));
                if(dp < SMALL && dm < SMALL)
                {
                    mup = 0.5;
                    mum = 1 - mup;
                }
                else
                {
                    dpdm = dp + dm;
                    if(dpdm < SMALL) dpdm += SMALL;
                    mup = dm/dpdm;
                    mup = min(mup, 1);
                    mum = 1 - mup;
                }
            }
        }
    }

    for(i=0; i<n; i++)
    {
        gradYOwn[i] = alphap*(Y1p[i] - Yp[i])/magt1p + betap*(Y2p[i] - Yp[i])/magt2p;
        gradYNei[i] = alpham*(Y1m[i] - Ym[i])/magt1m + betam*(Y2m[i] - Ym[i])/magt2m;
        faceGrad_Y[i] = mup*gradYOwn[i] - mum*gradYNei[i];
    }
    gradTOwn = alphap*(T1p - Tp)/magt1p + betap*(T2p - Tp)/magt2p;
    gradTNei = alpham*(T1m - Tm)/magt1m + betam*(T2m - Tm)/magt2m;
    faceGrad_T = mup*gradTOwn - mum*gradTNei;

    if(debug)
    {
        os<< "mup = " << mup << "  mum = " << mum << nl
            << setw(7) << "Species" << "  " << setw(16) << "gradYOwn" << "  " << setw(16) << "gradYNei" << "  " << setw(16) << "faceGradY"
            << endl;
        for(i=0; i<n; i++)
        {
            os<< setw(7) << i << "  " << setw(16) << gradYOwn[i] << "  " << setw(16) << gradYNei[i] << "  " << setw(16) << faceGrad_Y[i]
                << endl;
        }
        os<< "gradTOwn = " << gradTOwn << "  gradTNei = " << gradTNei << "  faceGradT = " << faceGrad_T << nl
            << endl;
    }
}


template<class Type>
void makeFlatFld
(
    const GeometricField<Type, fvPatchField, volMesh>& fld,
    const mapDistribute& map,
    List<Type>& flatFld
)
{
    label nFlatFld = map.constructSize();
    flatFld.resize(nFlatFld);

    forAll(fld, cellI)
    {
        flatFld[cellI] = fld[cellI];
    }

    forAll(fld.boundaryField(), patchI)
    {
        const fvPatchField<Type>& pfld = fld.boundaryField()[patchI];

        label nCompact = 
            pfld.patch().start() 
            -fld.mesh().nInternalFaces()
            +fld.mesh().nCells();

        forAll(pfld, i)
        {
            flatFld[nCompact++] = pfld[i];
        }
    }

    // Do all swapping
    map.distribute(flatFld);
}


void calc_2ph_gradf
(    
    const fvMesh& mesh,
    const labelListList& diffCellStencil,
    const volScalarField& Y1i,
    const volScalarField& Y0i,
    const List<scalar>& Y1i_flatFld_diff,
    const List<scalar>& Y0i_flatFld_diff,
    const List<scalar>& alpha1_flatFld_diff,
    const List<vector>& C_ph1_flatFld_diff,
    const List<vector>& C_ph0_flatFld_diff,    
    const labelList& face_phaseState_diff,
    surfaceScalarField& gradf_Y1i,
    surfaceScalarField& gradf_Y0i,
    const label& i,
    const scalar& MIN_ALPHA_DIFF,
    const bool debug,
    OFstream& os
)
{   
    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();
    const surfaceVectorField& meshSf = mesh.Sf();
    const surfaceScalarField& meshMagSf = mesh.magSf();    

    if(debug)
    {
        os<< "Gradient calculation" << nl
            << nl
            << "Internal faces" << nl
            << endl;
    }

    scalar curMagSf;    
    label curPhaseState;
    vector nf;
    label faceOwn; 
    label faceNei;     
    scalar wOwn; 
    scalar wNei;

    //Internal faces    
    const scalarField& Y1iCells = Y1i.internalField();
    const scalarField& Y0iCells = Y0i.internalField();

    for(label faceI=0; faceI<mesh.nInternalFaces(); faceI++)
    {
        curMagSf = meshMagSf[faceI];        
        nf = meshSf[faceI]/curMagSf;
        faceOwn = own[faceI];
        faceNei = nei[faceI];
        curPhaseState = face_phaseState_diff[faceI];

        if(debug)
        {
            os<< "Face: " << faceI << "  mag(Sf) = " << curMagSf << nl                
                << "phase state for gradient calculation: " << curPhaseState << nl
                << "Own: " << faceOwn << "  Nei: " << faceNei << nl
                << "Own Y_ph1 = " << Y1i_flatFld_diff[faceOwn] << "  Nei Y_ph1 = " << Y1i_flatFld_diff[faceNei] << nl
                << "Own Y_ph0 = " << Y0i_flatFld_diff[faceOwn] << "  Nei Y_ph0 = " << Y0i_flatFld_diff[faceNei] << nl
                
                << endl;
        }

        wOwn = 1;
        wNei = 1;

        if(curPhaseState == 3)
        {
            gradf_Y0i[faceI] = 0;
            gradf_Y1i[faceI] = 0;
        }//end if(curPhaseState == 3)
        else if(curPhaseState == 0)
        {
            plicFuncs::calcTwoSidedFaceGradWeights
            (
                faceI,
                faceOwn,
                faceNei,
                nf,
                Y0i_flatFld_diff,
                alpha1_flatFld_diff,
                C_ph0_flatFld_diff,
                diffCellStencil,
                0.1*MIN_ALPHA_DIFF,
                0,
                wOwn,
                wNei,
                debug, 
                os
            );
            
            gradf_Y0i[faceI] = wOwn*Y0i_flatFld_diff[faceOwn] + wNei*Y0i_flatFld_diff[faceNei];            
            gradf_Y1i[faceI] = 0;
        }//end if(curPhaseState == 0)        
        else if(curPhaseState == 1)
        {
            plicFuncs::calcTwoSidedFaceGradWeights
            (
                faceI,
                faceOwn,
                faceNei,
                nf,
                Y1i_flatFld_diff,
                alpha1_flatFld_diff,
                C_ph1_flatFld_diff,
                diffCellStencil,
                0.1*MIN_ALPHA_DIFF,
                1,
                wOwn,
                wNei,
                debug, 
                os
            );
            
            gradf_Y1i[faceI] = wOwn*Y1i_flatFld_diff[faceOwn] + wNei*Y1i_flatFld_diff[faceNei];            
            gradf_Y0i[faceI] = 0;
        }//end if(curPhaseState == 1)
        else
        {
            plicFuncs::calcTwoSidedFaceGradWeights
            (
                faceI,
                faceOwn,
                faceNei,
                nf,
                Y0i_flatFld_diff,
                alpha1_flatFld_diff,
                C_ph0_flatFld_diff,
                diffCellStencil,
                0.1*MIN_ALPHA_DIFF,
                0,
                wOwn,
                wNei,
                debug, 
                os
            );
            
            gradf_Y0i[faceI] = wOwn*Y0i_flatFld_diff[faceOwn] + wNei*Y0i_flatFld_diff[faceNei];

            plicFuncs::calcTwoSidedFaceGradWeights
            (
                faceI,
                faceOwn,
                faceNei,
                nf,
                Y1i_flatFld_diff,
                alpha1_flatFld_diff,
                C_ph1_flatFld_diff,
                diffCellStencil,
                0.1*MIN_ALPHA_DIFF,
                1,
                wOwn,
                wNei,
                debug, 
                os
            );
            
            gradf_Y1i[faceI] = wOwn*Y1i_flatFld_diff[faceOwn] + wNei*Y1i_flatFld_diff[faceNei];
        }//end if(curPhaseState == 2)

        if(debug)
        {
            os<< nl
                << "gradient ph1 for Y" << i << " = " << gradf_Y1i[faceI] << nl
                << "gradient ph0 for Y" << i << " = " << gradf_Y0i[faceI] << nl
                << endl;
        }        
    }//end for(label faceI=0; faceI<mesh.nInternalFaces(); faceI++)

    //end internal faces

    if(debug)
    {
        os<< "Boundary faces" << nl
            << endl;
    }

    //Boundary faces
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const label nBnd = mesh.nFaces() - mesh.nInternalFaces();

    List<scalar> ownYi_ph1(nBnd);
    List<scalar> ownAlpha_ph1(nBnd);
    List<scalar> ownBeta_ph1(nBnd);
    List<scalar> ownMagt1_ph1(nBnd);
    List<scalar> ownMagt2_ph1(nBnd);
    List<scalar> ownD_ph1(nBnd);
    List<scalar> neiYi_ph1(nBnd);
    List<scalar> neiAlpha_ph1(nBnd);
    List<scalar> neiBeta_ph1(nBnd);
    List<scalar> neiMagt1_ph1(nBnd);
    List<scalar> neiMagt2_ph1(nBnd);
    List<scalar> neiD_ph1(nBnd);

    List<scalar> ownYi_ph0(nBnd);
    List<scalar> ownAlpha_ph0(nBnd);
    List<scalar> ownBeta_ph0(nBnd);
    List<scalar> ownMagt1_ph0(nBnd);
    List<scalar> ownMagt2_ph0(nBnd);
    List<scalar> ownD_ph0(nBnd);
    List<scalar> neiYi_ph0(nBnd);
    List<scalar> neiAlpha_ph0(nBnd);
    List<scalar> neiBeta_ph0(nBnd);
    List<scalar> neiMagt1_ph0(nBnd);
    List<scalar> neiMagt2_ph0(nBnd);
    List<scalar> neiD_ph0(nBnd);

    scalar alpha;
    scalar beta;
    scalar magt1;
    scalar magt2;
    scalar d;
    label bndFaceI;

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        const fvsPatchVectorField& pSf = meshSf.boundaryField()[patchI];
        const fvsPatchScalarField& pMagSf = meshMagSf.boundaryField()[patchI];

        if(pp.coupled())
        {
            if(debug)
            {
                os<< "---------------------------------------------------------------------------------" << nl
                    << "Calculation of weights for own and nei cell gradient for coupled patch " << patchI << nl
                    << "---------------------------------------------------------------------------------" << nl
                    << endl;
            }

            label faceI = pp.start();

            forAll(pp, fcI)
            {
                bndFaceI = faceI - mesh.nInternalFaces();                

                nf = pSf[fcI]/pMagSf[fcI];
                faceOwn = own[faceI];
                const labelList& ownCells = diffCellStencil[faceOwn];
                curPhaseState = face_phaseState_diff[faceI];

                ownYi_ph1 = Y1iCells[faceOwn];
                neiYi_ph1 = Y1iCells[faceOwn];
                ownYi_ph0 = Y0iCells[faceOwn];
                neiYi_ph0 = Y0iCells[faceOwn];                

                if(curPhaseState == 0)
                {                    
                    plicFuncs::calcCellGradWeights
                    (
                        faceOwn,
                        nf,
                        Y0i_flatFld_diff,
                        alpha1_flatFld_diff,
                        C_ph0_flatFld_diff,
                        ownCells,
                        0.1*MIN_ALPHA_DIFF,
                        0,
                        alpha,
                        beta,
                        magt1,
                        magt2,
                        d,
                        debug, 
                        os
                    );
                    ownAlpha_ph0[bndFaceI] = alpha;
                    ownBeta_ph0[bndFaceI] = beta;
                    ownMagt1_ph0[bndFaceI] = magt1;
                    ownMagt2_ph0[bndFaceI] = magt2;
                    ownD_ph0[bndFaceI] = d;
                    neiAlpha_ph0[bndFaceI] = alpha;
                    neiBeta_ph0[bndFaceI] = beta;
                    neiMagt1_ph0[bndFaceI] = magt1;
                    neiMagt2_ph0[bndFaceI] = magt2;
                    neiD_ph0[bndFaceI] = d;

                    ownAlpha_ph1[bndFaceI] = 0;
                    ownBeta_ph1[bndFaceI] = 0;
                    ownMagt1_ph1[bndFaceI] = 1;
                    ownMagt2_ph1[bndFaceI] = 1;
                    ownD_ph1[bndFaceI] = 1;
                    neiAlpha_ph1[bndFaceI] = 0;
                    neiBeta_ph1[bndFaceI] = 0;
                    neiMagt1_ph1[bndFaceI] = 1;
                    neiMagt2_ph1[bndFaceI] = 1;
                    neiD_ph1[bndFaceI] = 1;
                }//end if(curPhaseState == 0)
                else if(curPhaseState == 1)
                {
                    ownAlpha_ph0[bndFaceI] = 0;
                    ownBeta_ph0[bndFaceI] = 0;
                    ownMagt1_ph0[bndFaceI] = 1;
                    ownMagt2_ph0[bndFaceI] = 1;
                    ownD_ph0[bndFaceI] = 1;
                    neiAlpha_ph0[bndFaceI] = 0;
                    neiBeta_ph0[bndFaceI] = 0;
                    neiMagt1_ph0[bndFaceI] = 1;
                    neiMagt2_ph0[bndFaceI] = 1;
                    neiD_ph0[bndFaceI] = 1;

                    plicFuncs::calcCellGradWeights
                    (
                        faceOwn,
                        nf,
                        Y1i_flatFld_diff,
                        alpha1_flatFld_diff,
                        C_ph1_flatFld_diff,
                        ownCells,
                        0.1*MIN_ALPHA_DIFF,
                        1,
                        alpha,
                        beta,
                        magt1,
                        magt2,
                        d,
                        debug, 
                        os
                    );
                    ownAlpha_ph1[bndFaceI] = alpha;
                    ownBeta_ph1[bndFaceI] = beta;
                    ownMagt1_ph1[bndFaceI] = magt1;
                    ownMagt2_ph1[bndFaceI] = magt2;
                    ownD_ph1[bndFaceI] = d;
                    neiAlpha_ph1[bndFaceI] = alpha;
                    neiBeta_ph1[bndFaceI] = beta;
                    neiMagt1_ph1[bndFaceI] = magt1;
                    neiMagt2_ph1[bndFaceI] = magt2;
                    neiD_ph1[bndFaceI] = d;
                }//end else if(curPhaseState == 0)
                else if(curPhaseState == 2)
                {                    
                    plicFuncs::calcCellGradWeights
                    (
                        faceOwn,
                        nf,
                        Y0i_flatFld_diff,
                        alpha1_flatFld_diff,
                        C_ph0_flatFld_diff,
                        ownCells,
                        0.1*MIN_ALPHA_DIFF,
                        0,
                        alpha,
                        beta,
                        magt1,
                        magt2,
                        d,
                        debug, 
                        os
                    );
                    ownAlpha_ph0[bndFaceI] = alpha;
                    ownBeta_ph0[bndFaceI] = beta;
                    ownMagt1_ph0[bndFaceI] = magt1;
                    ownMagt2_ph0[bndFaceI] = magt2;
                    ownD_ph0[bndFaceI] = d;
                    neiAlpha_ph0[bndFaceI] = alpha;
                    neiBeta_ph0[bndFaceI] = beta;
                    neiMagt1_ph0[bndFaceI] = magt1;
                    neiMagt2_ph0[bndFaceI] = magt2;
                    neiD_ph0[bndFaceI] = d;

                    plicFuncs::calcCellGradWeights
                    (
                        faceOwn,
                        nf,
                        Y1i_flatFld_diff,
                        alpha1_flatFld_diff,
                        C_ph1_flatFld_diff,
                        ownCells,
                        0.1*MIN_ALPHA_DIFF,
                        1,
                        alpha,
                        beta,
                        magt1,
                        magt2,
                        d,
                        debug, 
                        os
                    );
                    ownAlpha_ph1[bndFaceI] = alpha;
                    ownBeta_ph1[bndFaceI] = beta;
                    ownMagt1_ph1[bndFaceI] = magt1;
                    ownMagt2_ph1[bndFaceI] = magt2;
                    ownD_ph1[bndFaceI] = d;
                    neiAlpha_ph1[bndFaceI] = alpha;
                    neiBeta_ph1[bndFaceI] = beta;
                    neiMagt1_ph1[bndFaceI] = magt1;
                    neiMagt2_ph1[bndFaceI] = magt2;
                    neiD_ph1[bndFaceI] = d;
                }//end if(curPhaseState == 2)
                else
                {
                    ownAlpha_ph0[bndFaceI] = 0;
                    ownBeta_ph0[bndFaceI] = 0;
                    ownMagt1_ph0[bndFaceI] = 1;
                    ownMagt2_ph0[bndFaceI] = 1;
                    ownD_ph0[bndFaceI] = 1;
                    neiAlpha_ph0[bndFaceI] = 0;
                    neiBeta_ph0[bndFaceI] = 0;
                    neiMagt1_ph0[bndFaceI] = 1;
                    neiMagt2_ph0[bndFaceI] = 1;
                    neiD_ph0[bndFaceI] = 1;

                    ownAlpha_ph1[bndFaceI] = 0;
                    ownBeta_ph1[bndFaceI] = 0;
                    ownMagt1_ph1[bndFaceI] = 1;
                    ownMagt2_ph1[bndFaceI] = 1;
                    ownD_ph1[bndFaceI] = 1;
                    neiAlpha_ph1[bndFaceI] = 0;
                    neiBeta_ph1[bndFaceI] = 0;
                    neiMagt1_ph1[bndFaceI] = 1;
                    neiMagt2_ph1[bndFaceI] = 1;
                    neiD_ph1[bndFaceI] = 1;
                }//end if(curPhaseState == 3)

                faceI++;
            }//end forAll(pp, fcI)
        }//end if(pp.coupled())
    }//end forAll(patches,patchI)

    syncTools::swapBoundaryFaceList(mesh, neiYi_ph0);
    syncTools::swapBoundaryFaceList(mesh, neiAlpha_ph0);
    syncTools::swapBoundaryFaceList(mesh, neiBeta_ph0);
    syncTools::swapBoundaryFaceList(mesh, neiMagt1_ph0);
    syncTools::swapBoundaryFaceList(mesh, neiMagt2_ph0);
    syncTools::swapBoundaryFaceList(mesh, neiD_ph0);

    syncTools::swapBoundaryFaceList(mesh, neiYi_ph1);
    syncTools::swapBoundaryFaceList(mesh, neiAlpha_ph1);
    syncTools::swapBoundaryFaceList(mesh, neiBeta_ph1);
    syncTools::swapBoundaryFaceList(mesh, neiMagt1_ph1);
    syncTools::swapBoundaryFaceList(mesh, neiMagt2_ph1);
    syncTools::swapBoundaryFaceList(mesh, neiD_ph1);

    if(debug)
    {
        os<< "Done calculation of weights for own and nei cell gradient for coupled patches" << nl
            << "---------------------------------------------------------------------------------" << nl 
            << endl;
    }

    scalar faceGrad = 0;

    forAll(Y1i.boundaryField(), patchI)
    {
        const polyPatch& pp = patches[patchI];        
        const fvPatchScalarField& pY1i = Y1i.boundaryField()[patchI];
        fvsPatchScalarField& pgradf_Y0i = gradf_Y0i.boundaryField()[patchI];
        fvsPatchScalarField& pgradf_Y1i = gradf_Y1i.boundaryField()[patchI];
        const fvsPatchVectorField& pSf = meshSf.boundaryField()[patchI];
        const fvsPatchScalarField& pMagSf = meshMagSf.boundaryField()[patchI];
        label faceI = pp.start();        

        if(pp.coupled())
        {
            forAll(pY1i, fcI)
            {
                bndFaceI = faceI - mesh.nInternalFaces();                                

                curPhaseState = face_phaseState_diff[faceI];

                if(curPhaseState == 0)
                {                
                    plicFuncs::calcFaceGradFromWeights
                    (
                        ownYi_ph0[bndFaceI],
                        neiYi_ph0[bndFaceI],
                        ownAlpha_ph0[bndFaceI],
                        ownBeta_ph0[bndFaceI],
                        ownMagt1_ph0[bndFaceI],
                        ownMagt2_ph0[bndFaceI],
                        ownD_ph0[bndFaceI],
                        neiAlpha_ph0[bndFaceI],
                        neiBeta_ph0[bndFaceI],
                        neiMagt1_ph0[bndFaceI],
                        neiMagt2_ph0[bndFaceI],
                        neiD_ph0[bndFaceI],
                        faceGrad,
                        debug, 
                        os
                    );
                    pgradf_Y0i[fcI] = faceGrad;
                    pgradf_Y1i[fcI] = 0;
                }//end if(curPhaseState == 0)
                else if(curPhaseState == 1)
                {
                    plicFuncs::calcFaceGradFromWeights
                    (
                        ownYi_ph1[bndFaceI],
                        neiYi_ph1[bndFaceI],
                        ownAlpha_ph1[bndFaceI],
                        ownBeta_ph1[bndFaceI],
                        ownMagt1_ph1[bndFaceI],
                        ownMagt2_ph1[bndFaceI],
                        ownD_ph1[bndFaceI],
                        neiAlpha_ph1[bndFaceI],
                        neiBeta_ph1[bndFaceI],
                        neiMagt1_ph1[bndFaceI],
                        neiMagt2_ph1[bndFaceI],
                        neiD_ph1[bndFaceI],
                        faceGrad,
                        debug, 
                        os
                    );
                    pgradf_Y1i[fcI] = faceGrad;
                    pgradf_Y0i[fcI] = 0;
                }//end if(curPhaseState == 1)
                else if(curPhaseState == 2)
                {
                    plicFuncs::calcFaceGradFromWeights
                    (
                        ownYi_ph0[bndFaceI],
                        neiYi_ph0[bndFaceI],
                        ownAlpha_ph0[bndFaceI],
                        ownBeta_ph0[bndFaceI],
                        ownMagt1_ph0[bndFaceI],
                        ownMagt2_ph0[bndFaceI],
                        ownD_ph0[bndFaceI],
                        neiAlpha_ph0[bndFaceI],
                        neiBeta_ph0[bndFaceI],
                        neiMagt1_ph0[bndFaceI],
                        neiMagt2_ph0[bndFaceI],
                        neiD_ph0[bndFaceI],
                        faceGrad,
                        debug, 
                        os
                    );
                    pgradf_Y0i[fcI] = faceGrad;

                    plicFuncs::calcFaceGradFromWeights
                    (
                        ownYi_ph1[bndFaceI],
                        neiYi_ph1[bndFaceI],
                        ownAlpha_ph1[bndFaceI],
                        ownBeta_ph1[bndFaceI],
                        ownMagt1_ph1[bndFaceI],
                        ownMagt2_ph1[bndFaceI],
                        ownD_ph1[bndFaceI],
                        neiAlpha_ph1[bndFaceI],
                        neiBeta_ph1[bndFaceI],
                        neiMagt1_ph1[bndFaceI],
                        neiMagt2_ph1[bndFaceI],
                        neiD_ph1[bndFaceI],
                        faceGrad,
                        debug, 
                        os
                    );
                    pgradf_Y1i[fcI] = faceGrad;
                }//end if(curPhaseState == 2)
                else
                {
                    pgradf_Y0i[fcI] = 0;
                    pgradf_Y1i[fcI] = 0;
                }//end if(curPhaseState == 3)
                
                faceI++;
            }//end forAll(pY1i, fcI)
        }//end if(pp.coupled())
        else if(isA<zeroGradientFvPatchScalarField>(pY1i))
        {            
            forAll(pY1i, fcI)
            {                
                pgradf_Y1i[fcI] = 0;                
                pgradf_Y0i[fcI] = 0;
                
                faceI++;
            }
        }//end if(isA<zeroGradientFvPatchScalarField>(pY1i))
        else if(isA<fixedValueFvPatchScalarField>(pY1i))
        {
            forAll(pY1i, fcI)
            {
                nf = pSf[fcI]/pMagSf[fcI];
                faceOwn = own[faceI];
                const labelList& ownCells = diffCellStencil[faceOwn];
                curPhaseState = face_phaseState_diff[faceI];                                
                scalar ownCellGrad = 0;
                
                if(curPhaseState == 0)
                {                
                    plicFuncs::calcCellGrad
                    (
                        faceOwn,
                        nf,
                        Y0i_flatFld_diff,
                        alpha1_flatFld_diff,
                        C_ph0_flatFld_diff,
                        ownCells,
                        0.1*MIN_ALPHA_DIFF,
                        0,
                        ownCellGrad,
                        debug, 
                        os
                    );
                    pgradf_Y0i[fcI] = ownCellGrad;
                    pgradf_Y1i[fcI] = 0;
                }//end if(curPhaseState == 0)
                else if(curPhaseState == 1)
                {
                    plicFuncs::calcCellGrad
                    (
                        faceOwn,
                        nf,
                        Y1i_flatFld_diff,
                        alpha1_flatFld_diff,
                        C_ph1_flatFld_diff,
                        ownCells,
                        0.1*MIN_ALPHA_DIFF,
                        1,
                        ownCellGrad,
                        debug, 
                        os
                    );
                    pgradf_Y1i[fcI] = ownCellGrad;
                    pgradf_Y0i[fcI] = 0;
                }//end if(curPhaseState == 1)
                else if(curPhaseState == 2)
                {
                    plicFuncs::calcCellGrad
                    (
                        faceOwn,
                        nf,
                        Y0i_flatFld_diff,
                        alpha1_flatFld_diff,
                        C_ph0_flatFld_diff,
                        ownCells,
                        0.1*MIN_ALPHA_DIFF,
                        0,
                        ownCellGrad,
                        debug, 
                        os
                    );
                    pgradf_Y0i[fcI] = ownCellGrad;

                    plicFuncs::calcCellGrad
                    (
                        faceOwn,
                        nf,
                        Y1i_flatFld_diff,
                        alpha1_flatFld_diff,
                        C_ph1_flatFld_diff,
                        ownCells,
                        0.1*MIN_ALPHA_DIFF,
                        1,
                        ownCellGrad,
                        debug, 
                        os
                    );
                    pgradf_Y1i[fcI] = ownCellGrad;
                }//end if(curPhaseState == 2)
                else
                {
                    pgradf_Y0i[fcI] = 0;
                    pgradf_Y1i[fcI] = 0;
                }//end if(curPhaseState == 3)

                faceI++;
            }//end forAll(pY1i, fcI)
        }//end if(isA<fixedValueFvPatchScalarField>(pY1i))
        else
        {
            forAll(pY1i, fcI)
            {                
                pgradf_Y1i[fcI] = 0;                
                pgradf_Y0i[fcI] = 0;

                faceI++;
            }
        }//end if(pp.coupled())        
    }//end forAll(patches,patchI)

    //end boundary faces
}


void calc_2ph_gradf
(    
    const fvMesh& mesh,
    const labelListList& diffCellStencil,
    const PtrList<volScalarField>& Y1,
    const PtrList<volScalarField>& Y0,
    const volScalarField& T1,
    const volScalarField& T0,
    const List<List<scalar> >& Y1_flatFld_diff,
    const List<List<scalar> >& Y0_flatFld_diff,
    const List<scalar>& T1_flatFld_diff,
    const List<scalar>& T0_flatFld_diff,
    const List<scalar>& alpha1_flatFld_diff,
    const List<vector>& C_ph1_flatFld_diff,
    const List<vector>& C_ph0_flatFld_diff,
    const vectorField& Cf_ph1_own,
    const vectorField& Cf_ph1_nei,
    const vectorField& Cf_ph0_own,
    const vectorField& Cf_ph0_nei,
    const labelList& face_phaseState_diff,
    PtrList<surfaceScalarField>& gradf_Y1,
    PtrList<surfaceScalarField>& gradf_Y0,
    surfaceScalarField& gradf_T1,
    surfaceScalarField& gradf_T0,
    const label& n,
    const scalar& MIN_ALPHA_DIFF,
    const bool debug,
    OFstream& os
)
{   
    label faceI, curPhaseState, faceOwn, faceNei, i, nBnd, bndFaceI;
    bool foundCell1p, foundCell2p, foundCell1m, foundCell2m;
    scalar curMagSf, alpha1Own, gradf_T1_faceI, gradf_T0_faceI;
    scalar alphap, betap, magt1p, magt2p, Tp, T1p, T2p;
    scalar alpham, betam, magt1m, magt2m, Tm, T1m, T2m;
    vector nf, Cf_ph1_faceI, Cf_ph0_faceI, Cp_ph1, Cp_ph0, Cm_ph1, Cm_ph0;

    List<scalar> gradf_Y1_faceI(n);
    List<scalar> gradf_Y0_faceI(n);
    List<scalar> Yp(n);
    List<scalar> Y1p(n);
    List<scalar> Y2p(n);
    List<scalar> Ym(n);
    List<scalar> Y1m(n);
    List<scalar> Y2m(n);

    nBnd = mesh.nFaces() - mesh.nInternalFaces();
    List<vector> Cp_ph1_bnd(nBnd);
    List<scalar> alphap_ph1_bnd(nBnd);
    List<scalar> betap_ph1_bnd(nBnd);
    List<scalar> magt1p_ph1_bnd(nBnd);
    List<scalar> magt2p_ph1_bnd(nBnd);
    List<List<scalar> > Yp_ph1_bnd(n);
    List<List<scalar> > Y1p_ph1_bnd(n);
    List<List<scalar> > Y2p_ph1_bnd(n);
    List<scalar> Tp_ph1_bnd(nBnd);
    List<scalar> T1p_ph1_bnd(nBnd);
    List<scalar> T2p_ph1_bnd(nBnd);
    List<bool> foundCell1p_ph1_bnd(nBnd);
    List<bool> foundCell2p_ph1_bnd(nBnd);

    List<vector> Cp_ph0_bnd(nBnd);
    List<scalar> alphap_ph0_bnd(nBnd);
    List<scalar> betap_ph0_bnd(nBnd);
    List<scalar> magt1p_ph0_bnd(nBnd);
    List<scalar> magt2p_ph0_bnd(nBnd);
    List<List<scalar> > Yp_ph0_bnd(n);
    List<List<scalar> > Y1p_ph0_bnd(n);
    List<List<scalar> > Y2p_ph0_bnd(n);
    List<scalar> Tp_ph0_bnd(nBnd);
    List<scalar> T1p_ph0_bnd(nBnd);
    List<scalar> T2p_ph0_bnd(nBnd);
    List<bool> foundCell1p_ph0_bnd(nBnd);
    List<bool> foundCell2p_ph0_bnd(nBnd);

    List<vector> Cm_ph1_bnd(nBnd);
    List<scalar> alpham_ph1_bnd(nBnd);
    List<scalar> betam_ph1_bnd(nBnd);
    List<scalar> magt1m_ph1_bnd(nBnd);
    List<scalar> magt2m_ph1_bnd(nBnd);
    List<List<scalar> > Ym_ph1_bnd(n);
    List<List<scalar> > Y1m_ph1_bnd(n);
    List<List<scalar> > Y2m_ph1_bnd(n);
    List<scalar> Tm_ph1_bnd(nBnd);
    List<scalar> T1m_ph1_bnd(nBnd);
    List<scalar> T2m_ph1_bnd(nBnd);
    List<bool> foundCell1m_ph1_bnd(nBnd);
    List<bool> foundCell2m_ph1_bnd(nBnd);

    List<vector> Cm_ph0_bnd(nBnd);
    List<scalar> alpham_ph0_bnd(nBnd);
    List<scalar> betam_ph0_bnd(nBnd);
    List<scalar> magt1m_ph0_bnd(nBnd);
    List<scalar> magt2m_ph0_bnd(nBnd);
    List<List<scalar> > Ym_ph0_bnd(n);
    List<List<scalar> > Y1m_ph0_bnd(n);
    List<List<scalar> > Y2m_ph0_bnd(n);
    List<scalar> Tm_ph0_bnd(nBnd);
    List<scalar> T1m_ph0_bnd(nBnd);
    List<scalar> T2m_ph0_bnd(nBnd);
    List<bool> foundCell1m_ph0_bnd(nBnd);
    List<bool> foundCell2m_ph0_bnd(nBnd);

    for(i=0; i<n; i++)
    {
        Yp_ph1_bnd[i].setSize(nBnd);
        Yp_ph0_bnd[i].setSize(nBnd);
        Ym_ph1_bnd[i].setSize(nBnd);
        Ym_ph0_bnd[i].setSize(nBnd);
        Y1p_ph1_bnd[i].setSize(nBnd);
        Y1p_ph0_bnd[i].setSize(nBnd);
        Y1m_ph1_bnd[i].setSize(nBnd);
        Y1m_ph0_bnd[i].setSize(nBnd);
        Y2p_ph1_bnd[i].setSize(nBnd);
        Y2p_ph0_bnd[i].setSize(nBnd);
        Y2m_ph1_bnd[i].setSize(nBnd);
        Y2m_ph0_bnd[i].setSize(nBnd);
    }    

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();
    const surfaceVectorField& meshSf = mesh.Sf();
    const surfaceScalarField& meshMagSf = mesh.magSf();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if(debug)
    {
        os<< "Gradient calculation" << nl
            << nl
            << "Internal faces" << nl
            << endl;
    }

    //Internal faces
    for(faceI=0; faceI<mesh.nInternalFaces(); faceI++)
    {
        curMagSf = meshMagSf[faceI];        
        nf = meshSf[faceI]/curMagSf;
        faceOwn = own[faceI];
        faceNei = nei[faceI];
        curPhaseState = face_phaseState_diff[faceI];
        Cf_ph1_faceI = 0.5*(Cf_ph1_own[faceI] + Cf_ph1_nei[faceI]);
        Cf_ph0_faceI = 0.5*(Cf_ph0_own[faceI] + Cf_ph0_nei[faceI]);
        Cp_ph1 = C_ph1_flatFld_diff[faceOwn];
        Cp_ph0 = C_ph0_flatFld_diff[faceOwn];
        Cm_ph1 = C_ph1_flatFld_diff[faceNei];
        Cm_ph0 = C_ph0_flatFld_diff[faceNei];

        if(debug)
        {
            os<< "Face: " << faceI << "  mag(Sf) = " << curMagSf << nl                
                << "phase state for gradient calculation: " << curPhaseState << nl
                << "Own: " << faceOwn << "  Nei: " << faceNei << endl;
            print_line(os, 80);
            os<< setw(7) << "Species" << "  " << setw(16) << "Y1Own" << "  " << setw(16) << "Y0Own" << "  " << setw(16) << "Y1Nei" << "  " << setw(16) << "Y0Nei" << endl;
            print_line(os, 80);
            for(i=0; i<n; i++)
            {
                os<< setw(7) << i << "  " << setw(16) << Y1_flatFld_diff[i][faceOwn] << "  " << setw(16) << Y0_flatFld_diff[i][faceOwn] << "  " << setw(16) << Y1_flatFld_diff[i][faceNei] << "  " << setw(16) << Y0_flatFld_diff[i][faceNei] << endl;
            }
            print_line(os, 80);
            os<< "T1Own = " << T1_flatFld_diff[faceOwn] << "  T0Own = " << T0_flatFld_diff[faceOwn] << "T1Nei = " << T1_flatFld_diff[faceNei] << "  T0Nei = " << T0_flatFld_diff[faceNei] << endl;
            print_line(os, 80);
            os<< endl;
        }

        if(curPhaseState == 3)
        {
            for(i=0; i<n; i++)
            {
                gradf_Y1_faceI[i] = 0;
                gradf_Y0_faceI[i] = 0;
            }
            gradf_T1_faceI = 0;
            gradf_T0_faceI = 0;
        }//end if(curPhaseState == 3)
        else if(curPhaseState == 0)
        {
            //ph-0
            //own
            calcCellGradWeights(faceOwn, nf, Y0_flatFld_diff, T0_flatFld_diff, alpha1_flatFld_diff, C_ph0_flatFld_diff, diffCellStencil[faceOwn], 0.1*MIN_ALPHA_DIFF, 0, alphap, betap, magt1p, magt2p, Yp, Y1p, Y2p, Tp, T1p, T2p, foundCell1p, foundCell2p, n, debug, os);
            //nei
            calcCellGradWeights(faceNei, -nf, Y0_flatFld_diff, T0_flatFld_diff, alpha1_flatFld_diff, C_ph0_flatFld_diff, diffCellStencil[faceNei], 0.1*MIN_ALPHA_DIFF, 0, alpham, betam, magt1m, magt2m, Ym, Y1m, Y2m, Tm, T1m, T2m, foundCell1m, foundCell2m, n, debug, os);

            calcFaceGradFromWeights(nf, Cf_ph0_faceI, Cp_ph0, alphap, betap, magt1p, magt2p, Yp, Y1p, Y2p, Tp, T1p, T2p, foundCell1p, foundCell2p, Cm_ph0, alpham, betam, magt1m, magt2m, Ym, Y1m, Y2m, Tm, T1m, T2m, foundCell1m, foundCell2m, gradf_Y0_faceI, gradf_T0_faceI, n, debug, os);
            //ph-1
            for(i=0; i<n; i++)
            {
                gradf_Y1_faceI[i] = 0;
            }
            gradf_T1_faceI = 0;
        }//end if(curPhaseState == 0)        
        else if(curPhaseState == 1)
        {
            //ph-1
            //own
            calcCellGradWeights(faceOwn, nf, Y1_flatFld_diff, T1_flatFld_diff, alpha1_flatFld_diff, C_ph1_flatFld_diff, diffCellStencil[faceOwn], 0.1*MIN_ALPHA_DIFF, 1, alphap, betap, magt1p, magt2p, Yp, Y1p, Y2p, Tp, T1p, T2p, foundCell1p, foundCell2p, n, debug, os);
            //nei
            calcCellGradWeights(faceNei, -nf, Y1_flatFld_diff, T1_flatFld_diff, alpha1_flatFld_diff, C_ph1_flatFld_diff, diffCellStencil[faceNei], 0.1*MIN_ALPHA_DIFF, 1, alpham, betam, magt1m, magt2m, Ym, Y1m, Y2m, Tm, T1m, T2m, foundCell1m, foundCell2m, n, debug, os);

            calcFaceGradFromWeights(nf, Cf_ph1_faceI, Cp_ph1, alphap, betap, magt1p, magt2p, Yp, Y1p, Y2p, Tp, T1p, T2p, foundCell1p, foundCell2p, Cm_ph1, alpham, betam, magt1m, magt2m, Ym, Y1m, Y2m, Tm, T1m, T2m, foundCell1m, foundCell2m, gradf_Y1_faceI, gradf_T1_faceI, n, debug, os);
            //ph-0
            for(i=0; i<n; i++)
            {
                gradf_Y0_faceI[i] = 0;
            }
            gradf_T0_faceI = 0;
        }//end if(curPhaseState == 1)
        else
        {
            //ph-1
            //own
            calcCellGradWeights(faceOwn, nf, Y1_flatFld_diff, T1_flatFld_diff, alpha1_flatFld_diff, C_ph1_flatFld_diff, diffCellStencil[faceOwn], 0.1*MIN_ALPHA_DIFF, 1, alphap, betap, magt1p, magt2p, Yp, Y1p, Y2p, Tp, T1p, T2p, foundCell1p, foundCell2p, n, debug, os);
            //nei
            calcCellGradWeights(faceNei, -nf, Y1_flatFld_diff, T1_flatFld_diff, alpha1_flatFld_diff, C_ph1_flatFld_diff, diffCellStencil[faceNei], 0.1*MIN_ALPHA_DIFF, 1, alpham, betam, magt1m, magt2m, Ym, Y1m, Y2m, Tm, T1m, T2m, foundCell1m, foundCell2m, n, debug, os);

            calcFaceGradFromWeights(nf, Cf_ph1_faceI, Cp_ph1, alphap, betap, magt1p, magt2p, Yp, Y1p, Y2p, Tp, T1p, T2p, foundCell1p, foundCell2p, Cm_ph1, alpham, betam, magt1m, magt2m, Ym, Y1m, Y2m, Tm, T1m, T2m, foundCell1m, foundCell2m, gradf_Y1_faceI, gradf_T1_faceI, n, debug, os);
            //ph-0
            //own
            calcCellGradWeights(faceOwn, nf, Y0_flatFld_diff, T0_flatFld_diff, alpha1_flatFld_diff, C_ph0_flatFld_diff, diffCellStencil[faceOwn], 0.1*MIN_ALPHA_DIFF, 0, alphap, betap, magt1p, magt2p, Yp, Y1p, Y2p, Tp, T1p, T2p, foundCell1p, foundCell2p, n, debug, os);
            //nei
            calcCellGradWeights(faceNei, -nf, Y0_flatFld_diff, T0_flatFld_diff, alpha1_flatFld_diff, C_ph0_flatFld_diff, diffCellStencil[faceNei], 0.1*MIN_ALPHA_DIFF, 0, alpham, betam, magt1m, magt2m, Ym, Y1m, Y2m, Tm, T1m, T2m, foundCell1m, foundCell2m, n, debug, os);

            calcFaceGradFromWeights(nf, Cf_ph0_faceI, Cp_ph0, alphap, betap, magt1p, magt2p, Yp, Y1p, Y2p, Tp, T1p, T2p, foundCell1p, foundCell2p, Cm_ph0, alpham, betam, magt1m, magt2m, Ym, Y1m, Y2m, Tm, T1m, T2m, foundCell1m, foundCell2m, gradf_Y0_faceI, gradf_T0_faceI, n, debug, os);
        }//end if(curPhaseState == 2)

        for(i=0; i<n; i++)
        {
            gradf_Y1[i][faceI] = gradf_Y1_faceI[i];
            gradf_Y0[i][faceI] = gradf_Y0_faceI[i];
        }
        gradf_T1[faceI] = gradf_T1_faceI;
        gradf_T0[faceI] = gradf_T0_faceI;

        if(debug)
        {
            print_line(os, 80);
            os<< setw(7) << "Species" << "  " << setw(16) << "gradf_Y1" << "  " << setw(16) << "  " << "gradf_Y0" << endl;
            print_line(os, 80);
            for(i=0; i<n; i++)
            {
                os<< setw(7) << i << "  " << setw(16) << gradf_Y1_faceI[i] << "  " << setw(16) << "  " << gradf_Y0_faceI[i] << endl;
            }
            print_line(os, 80);
            os<< "gradf_T1 = " << gradf_T1_faceI << "  gradf_T0 = " << gradf_T0_faceI << endl;
            print_line(os, 80);
            os<< endl;
        }        
    }//end for(label faceI=0; faceI<mesh.nInternalFaces(); faceI++)

    //end internal faces

    if(debug)
    {
        os<< "Boundary faces" << nl
            << endl;
    }

    //Boundary faces    
    forAll(Y1[0].boundaryField(), patchI)
    {
        const polyPatch& pp = patches[patchI];
        const fvPatchScalarField& pY10 = Y1[0].boundaryField()[patchI];
        const fvsPatchVectorField& pSf = meshSf.boundaryField()[patchI];
        const fvsPatchScalarField& pMagSf = meshMagSf.boundaryField()[patchI];

        if(pp.coupled())
        {
            if(debug)
            {
                os<< "---------------------------------------------------------------------------------" << nl
                    << "Calculation of weights for own and nei cell gradient for coupled patch " << patchI << nl
                    << "---------------------------------------------------------------------------------" << nl
                    << endl;
            }

            faceI = pp.start();

            forAll(pY10, fcI)
            {
                bndFaceI = faceI - mesh.nInternalFaces();                

                nf = pSf[fcI]/pMagSf[fcI];
                faceOwn = own[faceI];
                alpha1Own = alpha1_flatFld_diff[faceOwn];
                const labelList& ownCells = diffCellStencil[faceOwn];
                curPhaseState = face_phaseState_diff[faceI];

                if(curPhaseState == 0)
                {                    
                    //ph-0
                    //own
                    calcCellGradWeights(faceOwn, nf, Y0_flatFld_diff, T0_flatFld_diff, alpha1_flatFld_diff, C_ph0_flatFld_diff, ownCells, 0.1*MIN_ALPHA_DIFF, 0, alphap, betap, magt1p, magt2p, Yp, Y1p, Y2p, Tp, T1p, T2p, foundCell1p, foundCell2p, n, debug, os);
                    
                    Cp_ph0_bnd[bndFaceI] = C_ph0_flatFld_diff[faceOwn];
                    alphap_ph0_bnd[bndFaceI] = alphap;
                    betap_ph0_bnd[bndFaceI] = betap;
                    magt1p_ph0_bnd[bndFaceI] = magt1p;
                    magt2p_ph0_bnd[bndFaceI] = magt2p;
                    Tp_ph0_bnd[bndFaceI] = Tp;
                    T1p_ph0_bnd[bndFaceI] = T1p;
                    T2p_ph0_bnd[bndFaceI] = T2p;
                    foundCell1p_ph0_bnd[bndFaceI] = foundCell1p;
                    foundCell2p_ph0_bnd[bndFaceI] = foundCell2p;
                    for(i=0; i<n; i++)
                    {
                        Yp_ph0_bnd[i][bndFaceI] = Yp[i];
                        Y1p_ph0_bnd[i][bndFaceI] = Y1p[i];
                        Y2p_ph0_bnd[i][bndFaceI] = Y2p[i];
                    }
                    //nei
                    Cm_ph0_bnd[bndFaceI] = C_ph0_flatFld_diff[faceOwn];
                    alpham_ph0_bnd[bndFaceI] = alphap;
                    betam_ph0_bnd[bndFaceI] = betap;
                    magt1m_ph0_bnd[bndFaceI] = magt1p;
                    magt2m_ph0_bnd[bndFaceI] = magt2p;
                    Tm_ph0_bnd[bndFaceI] = Tp;
                    T1m_ph0_bnd[bndFaceI] = T1p;
                    T2m_ph0_bnd[bndFaceI] = T2p;
                    foundCell1m_ph0_bnd[bndFaceI] = foundCell1p;
                    foundCell2m_ph0_bnd[bndFaceI] = foundCell2p;
                    for(i=0; i<n; i++)
                    {
                        Ym_ph0_bnd[i][bndFaceI] = Yp[i];
                        Y1m_ph0_bnd[i][bndFaceI] = Y1p[i];
                        Y2m_ph0_bnd[i][bndFaceI] = Y2p[i];
                    }

                    //ph-1
                    //own
                    Cp_ph1_bnd[bndFaceI] = C_ph1_flatFld_diff[faceOwn];
                    alphap_ph1_bnd[bndFaceI] = alphap;
                    betap_ph1_bnd[bndFaceI] = betap;
                    magt1p_ph1_bnd[bndFaceI] = magt1p;
                    magt2p_ph1_bnd[bndFaceI] = magt2p;
                    Tp_ph1_bnd[bndFaceI] = Tp;
                    T1p_ph1_bnd[bndFaceI] = T1p;
                    T2p_ph1_bnd[bndFaceI] = T2p;
                    foundCell1p_ph1_bnd[bndFaceI] = foundCell1p;
                    foundCell2p_ph1_bnd[bndFaceI] = foundCell2p;
                    for(i=0; i<n; i++)
                    {
                        Yp_ph1_bnd[i][bndFaceI] = Yp[i];
                        Y1p_ph1_bnd[i][bndFaceI] = Y1p[i];
                        Y2p_ph1_bnd[i][bndFaceI] = Y2p[i];
                    }
                    //nei
                    Cm_ph1_bnd[bndFaceI] = C_ph1_flatFld_diff[faceOwn];
                    alpham_ph1_bnd[bndFaceI] = alphap;
                    betam_ph1_bnd[bndFaceI] = betap;
                    magt1m_ph1_bnd[bndFaceI] = magt1p;
                    magt2m_ph1_bnd[bndFaceI] = magt2p;
                    Tm_ph1_bnd[bndFaceI] = Tp;
                    T1m_ph1_bnd[bndFaceI] = T1p;
                    T2m_ph1_bnd[bndFaceI] = T2p;
                    foundCell1m_ph1_bnd[bndFaceI] = foundCell1p;
                    foundCell2m_ph1_bnd[bndFaceI] = foundCell2p;
                    for(i=0; i<n; i++)
                    {
                        Ym_ph1_bnd[i][bndFaceI] = Yp[i];
                        Y1m_ph1_bnd[i][bndFaceI] = Y1p[i];
                        Y2m_ph1_bnd[i][bndFaceI] = Y2p[i];
                    }
                }//end if(curPhaseState == 0)
                else if(curPhaseState == 1)
                {
                    //ph-1
                    //own
                    calcCellGradWeights(faceOwn, nf, Y1_flatFld_diff, T1_flatFld_diff, alpha1_flatFld_diff, C_ph1_flatFld_diff, ownCells, 0.1*MIN_ALPHA_DIFF, 1, alphap, betap, magt1p, magt2p, Yp, Y1p, Y2p, Tp, T1p, T2p, foundCell1p, foundCell2p, n, debug, os);
                    
                    Cp_ph1_bnd[bndFaceI] = C_ph1_flatFld_diff[faceOwn];
                    alphap_ph1_bnd[bndFaceI] = alphap;
                    betap_ph1_bnd[bndFaceI] = betap;
                    magt1p_ph1_bnd[bndFaceI] = magt1p;
                    magt2p_ph1_bnd[bndFaceI] = magt2p;
                    Tp_ph1_bnd[bndFaceI] = Tp;
                    T1p_ph1_bnd[bndFaceI] = T1p;
                    T2p_ph1_bnd[bndFaceI] = T2p;
                    foundCell1p_ph1_bnd[bndFaceI] = foundCell1p;
                    foundCell2p_ph1_bnd[bndFaceI] = foundCell2p;
                    for(i=0; i<n; i++)
                    {
                        Yp_ph1_bnd[i][bndFaceI] = Yp[i];
                        Y1p_ph1_bnd[i][bndFaceI] = Y1p[i];
                        Y2p_ph1_bnd[i][bndFaceI] = Y2p[i];
                    }
                    //nei
                    Cm_ph1_bnd[bndFaceI] = C_ph1_flatFld_diff[faceOwn];
                    alpham_ph1_bnd[bndFaceI] = alphap;
                    betam_ph1_bnd[bndFaceI] = betap;
                    magt1m_ph1_bnd[bndFaceI] = magt1p;
                    magt2m_ph1_bnd[bndFaceI] = magt2p;
                    Tm_ph1_bnd[bndFaceI] = Tp;
                    T1m_ph1_bnd[bndFaceI] = T1p;
                    T2m_ph1_bnd[bndFaceI] = T2p;
                    foundCell1m_ph1_bnd[bndFaceI] = foundCell1p;
                    foundCell2m_ph1_bnd[bndFaceI] = foundCell2p;
                    for(i=0; i<n; i++)
                    {
                        Ym_ph1_bnd[i][bndFaceI] = Yp[i];
                        Y1m_ph1_bnd[i][bndFaceI] = Y1p[i];
                        Y2m_ph1_bnd[i][bndFaceI] = Y2p[i];
                    }

                    //ph-0
                    //own
                    Cp_ph0_bnd[bndFaceI] = C_ph0_flatFld_diff[faceOwn];
                    alphap_ph0_bnd[bndFaceI] = alphap;
                    betap_ph0_bnd[bndFaceI] = betap;
                    magt1p_ph0_bnd[bndFaceI] = magt1p;
                    magt2p_ph0_bnd[bndFaceI] = magt2p;
                    Tp_ph0_bnd[bndFaceI] = Tp;
                    T1p_ph0_bnd[bndFaceI] = T1p;
                    T2p_ph0_bnd[bndFaceI] = T2p;
                    foundCell1p_ph0_bnd[bndFaceI] = foundCell1p;
                    foundCell2p_ph0_bnd[bndFaceI] = foundCell2p;
                    for(i=0; i<n; i++)
                    {
                        Yp_ph0_bnd[i][bndFaceI] = Yp[i];
                        Y1p_ph0_bnd[i][bndFaceI] = Y1p[i];
                        Y2p_ph0_bnd[i][bndFaceI] = Y2p[i];
                    }
                    //nei
                    Cm_ph0_bnd[bndFaceI] = C_ph0_flatFld_diff[faceOwn];
                    alpham_ph0_bnd[bndFaceI] = alphap;
                    betam_ph0_bnd[bndFaceI] = betap;
                    magt1m_ph0_bnd[bndFaceI] = magt1p;
                    magt2m_ph0_bnd[bndFaceI] = magt2p;
                    Tm_ph0_bnd[bndFaceI] = Tp;
                    T1m_ph0_bnd[bndFaceI] = T1p;
                    T2m_ph0_bnd[bndFaceI] = T2p;
                    foundCell1m_ph0_bnd[bndFaceI] = foundCell1p;
                    foundCell2m_ph0_bnd[bndFaceI] = foundCell2p;
                    for(i=0; i<n; i++)
                    {
                        Ym_ph0_bnd[i][bndFaceI] = Yp[i];
                        Y1m_ph0_bnd[i][bndFaceI] = Y1p[i];
                        Y2m_ph0_bnd[i][bndFaceI] = Y2p[i];
                    }
                }//end else if(curPhaseState == 0)
                else if(curPhaseState == 2)
                {                    
                    //ph-1
                    //own
                    calcCellGradWeights(faceOwn, nf, Y1_flatFld_diff, T1_flatFld_diff, alpha1_flatFld_diff, C_ph1_flatFld_diff, ownCells, 0.1*MIN_ALPHA_DIFF, 1, alphap, betap, magt1p, magt2p, Yp, Y1p, Y2p, Tp, T1p, T2p, foundCell1p, foundCell2p, n, debug, os);
                    
                    Cp_ph1_bnd[bndFaceI] = C_ph1_flatFld_diff[faceOwn];
                    alphap_ph1_bnd[bndFaceI] = alphap;
                    betap_ph1_bnd[bndFaceI] = betap;
                    magt1p_ph1_bnd[bndFaceI] = magt1p;
                    magt2p_ph1_bnd[bndFaceI] = magt2p;
                    Tp_ph1_bnd[bndFaceI] = Tp;
                    T1p_ph1_bnd[bndFaceI] = T1p;
                    T2p_ph1_bnd[bndFaceI] = T2p;
                    foundCell1p_ph1_bnd[bndFaceI] = foundCell1p;
                    foundCell2p_ph1_bnd[bndFaceI] = foundCell2p;
                    for(i=0; i<n; i++)
                    {
                        Yp_ph1_bnd[i][bndFaceI] = Yp[i];
                        Y1p_ph1_bnd[i][bndFaceI] = Y1p[i];
                        Y2p_ph1_bnd[i][bndFaceI] = Y2p[i];
                    }
                    //nei
                    Cm_ph1_bnd[bndFaceI] = C_ph1_flatFld_diff[faceOwn];
                    alpham_ph1_bnd[bndFaceI] = alphap;
                    betam_ph1_bnd[bndFaceI] = betap;
                    magt1m_ph1_bnd[bndFaceI] = magt1p;
                    magt2m_ph1_bnd[bndFaceI] = magt2p;
                    Tm_ph1_bnd[bndFaceI] = Tp;
                    T1m_ph1_bnd[bndFaceI] = T1p;
                    T2m_ph1_bnd[bndFaceI] = T2p;
                    foundCell1m_ph1_bnd[bndFaceI] = foundCell1p;
                    foundCell2m_ph1_bnd[bndFaceI] = foundCell2p;
                    for(i=0; i<n; i++)
                    {
                        Ym_ph1_bnd[i][bndFaceI] = Yp[i];
                        Y1m_ph1_bnd[i][bndFaceI] = Y1p[i];
                        Y2m_ph1_bnd[i][bndFaceI] = Y2p[i];
                    }

                    //ph-0
                    //own
                    calcCellGradWeights(faceOwn, nf, Y0_flatFld_diff, T0_flatFld_diff, alpha1_flatFld_diff, C_ph0_flatFld_diff, ownCells, 0.1*MIN_ALPHA_DIFF, 0, alphap, betap, magt1p, magt2p, Yp, Y1p, Y2p, Tp, T1p, T2p, foundCell1p, foundCell2p, n, debug, os);
                    
                    Cp_ph0_bnd[bndFaceI] = C_ph0_flatFld_diff[faceOwn];
                    alphap_ph0_bnd[bndFaceI] = alphap;
                    betap_ph0_bnd[bndFaceI] = betap;
                    magt1p_ph0_bnd[bndFaceI] = magt1p;
                    magt2p_ph0_bnd[bndFaceI] = magt2p;
                    Tp_ph0_bnd[bndFaceI] = Tp;
                    T1p_ph0_bnd[bndFaceI] = T1p;
                    T2p_ph0_bnd[bndFaceI] = T2p;
                    foundCell1p_ph0_bnd[bndFaceI] = foundCell1p;
                    foundCell2p_ph0_bnd[bndFaceI] = foundCell2p;
                    for(i=0; i<n; i++)
                    {
                        Yp_ph0_bnd[i][bndFaceI] = Yp[i];
                        Y1p_ph0_bnd[i][bndFaceI] = Y1p[i];
                        Y2p_ph0_bnd[i][bndFaceI] = Y2p[i];
                    }
                    //nei
                    Cm_ph0_bnd[bndFaceI] = C_ph0_flatFld_diff[faceOwn];
                    alpham_ph0_bnd[bndFaceI] = alphap;
                    betam_ph0_bnd[bndFaceI] = betap;
                    magt1m_ph0_bnd[bndFaceI] = magt1p;
                    magt2m_ph0_bnd[bndFaceI] = magt2p;
                    Tm_ph0_bnd[bndFaceI] = Tp;
                    T1m_ph0_bnd[bndFaceI] = T1p;
                    T2m_ph0_bnd[bndFaceI] = T2p;
                    foundCell1m_ph0_bnd[bndFaceI] = foundCell1p;
                    foundCell2m_ph0_bnd[bndFaceI] = foundCell2p;
                    for(i=0; i<n; i++)
                    {
                        Ym_ph0_bnd[i][bndFaceI] = Yp[i];
                        Y1m_ph0_bnd[i][bndFaceI] = Y1p[i];
                        Y2m_ph0_bnd[i][bndFaceI] = Y2p[i];
                    }
                }//end if(curPhaseState == 2)
                else
                {
                    if(alpha1Own < MIN_ALPHA_DIFF)
                    {
                        //ph-0
                        //own
                        calcCellGradWeights(faceOwn, nf, Y0_flatFld_diff, T0_flatFld_diff, alpha1_flatFld_diff, C_ph0_flatFld_diff, ownCells, 0.1*MIN_ALPHA_DIFF, 0, alphap, betap, magt1p, magt2p, Yp, Y1p, Y2p, Tp, T1p, T2p, foundCell1p, foundCell2p, n, debug, os);
                    }
                    else
                    {
                        //ph-1
                        //own
                        calcCellGradWeights(faceOwn, nf, Y1_flatFld_diff, T1_flatFld_diff, alpha1_flatFld_diff, C_ph1_flatFld_diff, ownCells, 0.1*MIN_ALPHA_DIFF, 1, alphap, betap, magt1p, magt2p, Yp, Y1p, Y2p, Tp, T1p, T2p, foundCell1p, foundCell2p, n, debug, os);
                    }

                    //ph-1
                    //own
                    Cp_ph1_bnd[bndFaceI] = C_ph1_flatFld_diff[faceOwn];
                    alphap_ph1_bnd[bndFaceI] = alphap;
                    betap_ph1_bnd[bndFaceI] = betap;
                    magt1p_ph1_bnd[bndFaceI] = magt1p;
                    magt2p_ph1_bnd[bndFaceI] = magt2p;
                    Tp_ph1_bnd[bndFaceI] = Tp;
                    T1p_ph1_bnd[bndFaceI] = T1p;
                    T2p_ph1_bnd[bndFaceI] = T2p;
                    foundCell1p_ph1_bnd[bndFaceI] = foundCell1p;
                    foundCell2p_ph1_bnd[bndFaceI] = foundCell2p;
                    for(i=0; i<n; i++)
                    {
                        Yp_ph1_bnd[i][bndFaceI] = Yp[i];
                        Y1p_ph1_bnd[i][bndFaceI] = Y1p[i];
                        Y2p_ph1_bnd[i][bndFaceI] = Y2p[i];
                    }
                    //nei
                    Cm_ph1_bnd[bndFaceI] = C_ph1_flatFld_diff[faceOwn];
                    alpham_ph1_bnd[bndFaceI] = alphap;
                    betam_ph1_bnd[bndFaceI] = betap;
                    magt1m_ph1_bnd[bndFaceI] = magt1p;
                    magt2m_ph1_bnd[bndFaceI] = magt2p;
                    Tm_ph1_bnd[bndFaceI] = Tp;
                    T1m_ph1_bnd[bndFaceI] = T1p;
                    T2m_ph1_bnd[bndFaceI] = T2p;
                    foundCell1m_ph1_bnd[bndFaceI] = foundCell1p;
                    foundCell2m_ph1_bnd[bndFaceI] = foundCell2p;
                    for(i=0; i<n; i++)
                    {
                        Ym_ph1_bnd[i][bndFaceI] = Yp[i];
                        Y1m_ph1_bnd[i][bndFaceI] = Y1p[i];
                        Y2m_ph1_bnd[i][bndFaceI] = Y2p[i];
                    }

                    //ph-0
                    //own
                    Cp_ph0_bnd[bndFaceI] = C_ph0_flatFld_diff[faceOwn];
                    alphap_ph0_bnd[bndFaceI] = alphap;
                    betap_ph0_bnd[bndFaceI] = betap;
                    magt1p_ph0_bnd[bndFaceI] = magt1p;
                    magt2p_ph0_bnd[bndFaceI] = magt2p;
                    Tp_ph0_bnd[bndFaceI] = Tp;
                    T1p_ph0_bnd[bndFaceI] = T1p;
                    T2p_ph0_bnd[bndFaceI] = T2p;
                    foundCell1p_ph0_bnd[bndFaceI] = foundCell1p;
                    foundCell2p_ph0_bnd[bndFaceI] = foundCell2p;
                    for(i=0; i<n; i++)
                    {
                        Yp_ph0_bnd[i][bndFaceI] = Yp[i];
                        Y1p_ph0_bnd[i][bndFaceI] = Y1p[i];
                        Y2p_ph0_bnd[i][bndFaceI] = Y2p[i];
                    }
                    //nei
                    Cm_ph0_bnd[bndFaceI] = C_ph0_flatFld_diff[faceOwn];
                    alpham_ph0_bnd[bndFaceI] = alphap;
                    betam_ph0_bnd[bndFaceI] = betap;
                    magt1m_ph0_bnd[bndFaceI] = magt1p;
                    magt2m_ph0_bnd[bndFaceI] = magt2p;
                    Tm_ph0_bnd[bndFaceI] = Tp;
                    T1m_ph0_bnd[bndFaceI] = T1p;
                    T2m_ph0_bnd[bndFaceI] = T2p;
                    foundCell1m_ph0_bnd[bndFaceI] = foundCell1p;
                    foundCell2m_ph0_bnd[bndFaceI] = foundCell2p;
                    for(i=0; i<n; i++)
                    {
                        Ym_ph0_bnd[i][bndFaceI] = Yp[i];
                        Y1m_ph0_bnd[i][bndFaceI] = Y1p[i];
                        Y2m_ph0_bnd[i][bndFaceI] = Y2p[i];
                    }
                }//end if(curPhaseState == 3)

                faceI++;
            }//end forAll(pp, fcI)
        }//end if(pp.coupled())
    }//end forAll(patches,patchI)

    syncTools::swapBoundaryFaceList(mesh, Cm_ph0_bnd);
    syncTools::swapBoundaryFaceList(mesh, alpham_ph0_bnd);
    syncTools::swapBoundaryFaceList(mesh, betam_ph0_bnd);
    syncTools::swapBoundaryFaceList(mesh, magt1m_ph0_bnd);
    syncTools::swapBoundaryFaceList(mesh, magt2m_ph0_bnd);
    syncTools::swapBoundaryFaceList(mesh, Tm_ph0_bnd);
    syncTools::swapBoundaryFaceList(mesh, T1m_ph0_bnd);
    syncTools::swapBoundaryFaceList(mesh, T2m_ph0_bnd);
    syncTools::swapBoundaryFaceList(mesh, foundCell1m_ph0_bnd);
    syncTools::swapBoundaryFaceList(mesh, foundCell2m_ph0_bnd);
    for(i=0; i<n; i++)
    {
        syncTools::swapBoundaryFaceList(mesh, Ym_ph0_bnd[i]);
        syncTools::swapBoundaryFaceList(mesh, Y1m_ph0_bnd[i]);
        syncTools::swapBoundaryFaceList(mesh, Y2m_ph0_bnd[i]);
    }

    syncTools::swapBoundaryFaceList(mesh, Cm_ph1_bnd);
    syncTools::swapBoundaryFaceList(mesh, alpham_ph1_bnd);
    syncTools::swapBoundaryFaceList(mesh, betam_ph1_bnd);
    syncTools::swapBoundaryFaceList(mesh, magt1m_ph1_bnd);
    syncTools::swapBoundaryFaceList(mesh, magt2m_ph1_bnd);
    syncTools::swapBoundaryFaceList(mesh, Tm_ph1_bnd);
    syncTools::swapBoundaryFaceList(mesh, T1m_ph1_bnd);
    syncTools::swapBoundaryFaceList(mesh, T2m_ph1_bnd);
    syncTools::swapBoundaryFaceList(mesh, foundCell1m_ph1_bnd);
    syncTools::swapBoundaryFaceList(mesh, foundCell2m_ph1_bnd);
    for(i=0; i<n; i++)
    {
        syncTools::swapBoundaryFaceList(mesh, Ym_ph1_bnd[i]);
        syncTools::swapBoundaryFaceList(mesh, Y1m_ph1_bnd[i]);
        syncTools::swapBoundaryFaceList(mesh, Y2m_ph1_bnd[i]);
    }
    
    if(debug)
    {
        os<< "Done calculation of weights for own and nei cell gradient for coupled patches" << nl
            << "---------------------------------------------------------------------------------" << nl 
            << endl;
    }    

    forAll(Y1[0].boundaryField(), patchI)
    {
        const polyPatch& pp = patches[patchI];        
        const fvPatchScalarField& pY10 = Y1[0].boundaryField()[patchI];        
        const fvPatchScalarField& pT1 = T1.boundaryField()[patchI];        
        const fvsPatchVectorField& pSf = meshSf.boundaryField()[patchI];
        const fvsPatchScalarField& pMagSf = meshMagSf.boundaryField()[patchI];
        faceI = pp.start();

        if(pp.coupled())
        {
            forAll(pY10, fcI)
            {
                bndFaceI = faceI - mesh.nInternalFaces();
                curMagSf = pMagSf[fcI];        
                nf = pSf[fcI]/curMagSf;
                faceOwn = own[faceI];                
                curPhaseState = face_phaseState_diff[faceI];
                Cf_ph1_faceI = 0.5*(Cf_ph1_own[faceI] + Cf_ph1_nei[faceI]);
                Cf_ph0_faceI = 0.5*(Cf_ph0_own[faceI] + Cf_ph0_nei[faceI]);

                if(curPhaseState == 0)
                {
                    //ph-0
                    Cp_ph0 = Cp_ph0_bnd[bndFaceI];
                    alphap = alphap_ph0_bnd[bndFaceI];
                    betap = betap_ph0_bnd[bndFaceI];
                    magt1p = magt1p_ph0_bnd[bndFaceI];
                    magt2p = magt2p_ph0_bnd[bndFaceI];
                    Tp = Tp_ph0_bnd[bndFaceI];
                    T1p = T1p_ph0_bnd[bndFaceI];
                    T2p = T2p_ph0_bnd[bndFaceI];
                    foundCell1p = foundCell1p_ph0_bnd[bndFaceI];
                    foundCell2p = foundCell2p_ph0_bnd[bndFaceI];
                    Cm_ph0 = Cm_ph0_bnd[bndFaceI];
                    alpham = alpham_ph0_bnd[bndFaceI];
                    betam = betam_ph0_bnd[bndFaceI];
                    magt1m = magt1m_ph0_bnd[bndFaceI];
                    magt2m = magt2m_ph0_bnd[bndFaceI];
                    Tm = Tm_ph0_bnd[bndFaceI];
                    T1m = T1m_ph0_bnd[bndFaceI];
                    T2m = T2m_ph0_bnd[bndFaceI];
                    foundCell1m = foundCell1m_ph0_bnd[bndFaceI];
                    foundCell2m = foundCell2m_ph0_bnd[bndFaceI];
                    for(i=0; i<n; i++)
                    {
                        Yp[i] = Yp_ph0_bnd[i][bndFaceI];
                        Y1p[i] = Y1p_ph0_bnd[i][bndFaceI];
                        Y2p[i] = Y2p_ph0_bnd[i][bndFaceI];
                        Ym[i] = Ym_ph0_bnd[i][bndFaceI];
                        Y1m[i] = Y1m_ph0_bnd[i][bndFaceI];
                        Y2m[i] = Y2m_ph0_bnd[i][bndFaceI];
                    }

                    calcFaceGradFromWeights(nf, Cf_ph0_faceI, Cp_ph0, alphap, betap, magt1p, magt2p, Yp, Y1p, Y2p, Tp, T1p, T2p, foundCell1p, foundCell2p, Cm_ph0, alpham, betam, magt1m, magt2m, Ym, Y1m, Y2m, Tm, T1m, T2m, foundCell1m, foundCell2m, gradf_Y0_faceI, gradf_T0_faceI, n, debug, os);
                    //ph-1
                    gradf_T1_faceI = 0;
                    for(i=0; i<n; i++)
                    {
                        gradf_Y1_faceI[i] = 0;
                    }
                }//end if(curPhaseState == 0)
                else if(curPhaseState == 1)
                {
                    //ph-1
                    Cp_ph1 = Cp_ph1_bnd[bndFaceI];
                    alphap = alphap_ph1_bnd[bndFaceI];
                    betap = betap_ph1_bnd[bndFaceI];
                    magt1p = magt1p_ph1_bnd[bndFaceI];
                    magt2p = magt2p_ph1_bnd[bndFaceI];
                    Tp = Tp_ph1_bnd[bndFaceI];
                    T1p = T1p_ph1_bnd[bndFaceI];
                    T2p = T2p_ph1_bnd[bndFaceI];
                    foundCell1p = foundCell1p_ph1_bnd[bndFaceI];
                    foundCell2p = foundCell2p_ph1_bnd[bndFaceI];
                    Cm_ph1 = Cm_ph1_bnd[bndFaceI];
                    alpham = alpham_ph1_bnd[bndFaceI];
                    betam = betam_ph1_bnd[bndFaceI];
                    magt1m = magt1m_ph1_bnd[bndFaceI];
                    magt2m = magt2m_ph1_bnd[bndFaceI];
                    Tm = Tm_ph1_bnd[bndFaceI];
                    T1m = T1m_ph1_bnd[bndFaceI];
                    T2m = T2m_ph1_bnd[bndFaceI];
                    foundCell1m = foundCell1m_ph1_bnd[bndFaceI];
                    foundCell2m = foundCell2m_ph1_bnd[bndFaceI];
                    for(i=0; i<n; i++)
                    {
                        Yp[i] = Yp_ph1_bnd[i][bndFaceI];
                        Y1p[i] = Y1p_ph1_bnd[i][bndFaceI];
                        Y2p[i] = Y2p_ph1_bnd[i][bndFaceI];
                        Ym[i] = Ym_ph1_bnd[i][bndFaceI];
                        Y1m[i] = Y1m_ph1_bnd[i][bndFaceI];
                        Y2m[i] = Y2m_ph1_bnd[i][bndFaceI];
                    }

                    calcFaceGradFromWeights(nf, Cf_ph1_faceI, Cp_ph1, alphap, betap, magt1p, magt2p, Yp, Y1p, Y2p, Tp, T1p, T2p, foundCell1p, foundCell2p, Cm_ph1, alpham, betam, magt1m, magt2m, Ym, Y1m, Y2m, Tm, T1m, T2m, foundCell1m, foundCell2m, gradf_Y1_faceI, gradf_T1_faceI, n, debug, os);
                    //ph-0
                    gradf_T0_faceI = 0;
                    for(i=0; i<n; i++)
                    {
                        gradf_Y0_faceI[i] = 0;
                    }
                }//end if(curPhaseState == 1)
                else if(curPhaseState == 2)
                {
                    //ph-1
                    Cp_ph1 = Cp_ph1_bnd[bndFaceI];
                    alphap = alphap_ph1_bnd[bndFaceI];
                    betap = betap_ph1_bnd[bndFaceI];
                    magt1p = magt1p_ph1_bnd[bndFaceI];
                    magt2p = magt2p_ph1_bnd[bndFaceI];
                    Tp = Tp_ph1_bnd[bndFaceI];
                    T1p = T1p_ph1_bnd[bndFaceI];
                    T2p = T2p_ph1_bnd[bndFaceI];
                    foundCell1p = foundCell1p_ph1_bnd[bndFaceI];
                    foundCell2p = foundCell2p_ph1_bnd[bndFaceI];
                    Cm_ph1 = Cm_ph1_bnd[bndFaceI];
                    alpham = alpham_ph1_bnd[bndFaceI];
                    betam = betam_ph1_bnd[bndFaceI];
                    magt1m = magt1m_ph1_bnd[bndFaceI];
                    magt2m = magt2m_ph1_bnd[bndFaceI];
                    Tm = Tm_ph1_bnd[bndFaceI];
                    T1m = T1m_ph1_bnd[bndFaceI];
                    T2m = T2m_ph1_bnd[bndFaceI];
                    foundCell1m = foundCell1m_ph1_bnd[bndFaceI];
                    foundCell2m = foundCell2m_ph1_bnd[bndFaceI];
                    for(i=0; i<n; i++)
                    {
                        Yp[i] = Yp_ph1_bnd[i][bndFaceI];
                        Y1p[i] = Y1p_ph1_bnd[i][bndFaceI];
                        Y2p[i] = Y2p_ph1_bnd[i][bndFaceI];
                        Ym[i] = Ym_ph1_bnd[i][bndFaceI];
                        Y1m[i] = Y1m_ph1_bnd[i][bndFaceI];
                        Y2m[i] = Y2m_ph1_bnd[i][bndFaceI];
                    }

                    calcFaceGradFromWeights(nf, Cf_ph1_faceI, Cp_ph1, alphap, betap, magt1p, magt2p, Yp, Y1p, Y2p, Tp, T1p, T2p, foundCell1p, foundCell2p, Cm_ph1, alpham, betam, magt1m, magt2m, Ym, Y1m, Y2m, Tm, T1m, T2m, foundCell1m, foundCell2m, gradf_Y1_faceI, gradf_T1_faceI, n, debug, os);

                    //ph-0
                    Cp_ph0 = Cp_ph0_bnd[bndFaceI];
                    alphap = alphap_ph0_bnd[bndFaceI];
                    betap = betap_ph0_bnd[bndFaceI];
                    magt1p = magt1p_ph0_bnd[bndFaceI];
                    magt2p = magt2p_ph0_bnd[bndFaceI];
                    Tp = Tp_ph0_bnd[bndFaceI];
                    T1p = T1p_ph0_bnd[bndFaceI];
                    T2p = T2p_ph0_bnd[bndFaceI];
                    foundCell1p = foundCell1p_ph0_bnd[bndFaceI];
                    foundCell2p = foundCell2p_ph0_bnd[bndFaceI];
                    Cm_ph0 = Cm_ph0_bnd[bndFaceI];
                    alpham = alpham_ph0_bnd[bndFaceI];
                    betam = betam_ph0_bnd[bndFaceI];
                    magt1m = magt1m_ph0_bnd[bndFaceI];
                    magt2m = magt2m_ph0_bnd[bndFaceI];
                    Tm = Tm_ph0_bnd[bndFaceI];
                    T1m = T1m_ph0_bnd[bndFaceI];
                    T2m = T2m_ph0_bnd[bndFaceI];
                    foundCell1m = foundCell1m_ph0_bnd[bndFaceI];
                    foundCell2m = foundCell2m_ph0_bnd[bndFaceI];
                    for(i=0; i<n; i++)
                    {
                        Yp[i] = Yp_ph0_bnd[i][bndFaceI];
                        Y1p[i] = Y1p_ph0_bnd[i][bndFaceI];
                        Y2p[i] = Y2p_ph0_bnd[i][bndFaceI];
                        Ym[i] = Ym_ph0_bnd[i][bndFaceI];
                        Y1m[i] = Y1m_ph0_bnd[i][bndFaceI];
                        Y2m[i] = Y2m_ph0_bnd[i][bndFaceI];
                    }

                    calcFaceGradFromWeights(nf, Cf_ph0_faceI, Cp_ph0, alphap, betap, magt1p, magt2p, Yp, Y1p, Y2p, Tp, T1p, T2p, foundCell1p, foundCell2p, Cm_ph0, alpham, betam, magt1m, magt2m, Ym, Y1m, Y2m, Tm, T1m, T2m, foundCell1m, foundCell2m, gradf_Y0_faceI, gradf_T0_faceI, n, debug, os);
                }//end if(curPhaseState == 2)
                else
                {
                    //ph-1
                    gradf_T1_faceI = 0;
                    for(i=0; i<n; i++)
                    {
                        gradf_Y1_faceI[i] = 0;
                    }

                    //ph-0
                    gradf_T0_faceI = 0;
                    for(i=0; i<n; i++)
                    {
                        gradf_Y0_faceI[i] = 0;
                    }
                }//end if(curPhaseState == 3)
                
                for(i=0; i<n; i++)
                {
                    gradf_Y1[i].boundaryField()[patchI][fcI] = gradf_Y1_faceI[i];
                    gradf_Y0[i].boundaryField()[patchI][fcI] = gradf_Y0_faceI[i];
                }
                gradf_T1.boundaryField()[patchI][fcI] = gradf_T1_faceI;
                gradf_T0.boundaryField()[patchI][fcI] = gradf_T0_faceI;

                faceI++;
            }//end forAll(pY10, fcI)
        }//end if(pp.coupled())
        else if(isA<zeroGradientFvPatchScalarField>(pY10))
        {
            if(isA<zeroGradientFvPatchScalarField>(pT1))
            {
                forAll(pY10, fcI)
                {                
                    for(i=0; i<n; i++)
                    {
                        gradf_Y1[i].boundaryField()[patchI][fcI] = 0;
                        gradf_Y0[i].boundaryField()[patchI][fcI] = 0;
                    }
                    gradf_T1.boundaryField()[patchI][fcI] = 0;
                    gradf_T0.boundaryField()[patchI][fcI] = 0;                
                }
            }
            else
            {
                forAll(pY10, fcI)
                {
                    nf = pSf[fcI]/pMagSf[fcI];
                    faceOwn = own[faceI];
                    const labelList& ownCells = diffCellStencil[faceOwn];
                    curPhaseState = face_phaseState_diff[faceI];

                    if(curPhaseState == 0)
                    {   
                        //ph-0
                        calcCellGrad(faceOwn, nf, Y0_flatFld_diff, T0_flatFld_diff, alpha1_flatFld_diff, C_ph0_flatFld_diff, ownCells, 0.1*MIN_ALPHA_DIFF, 0, gradf_Y0_faceI, gradf_T0_faceI, n, debug, os);
                        //ph-1
                        for(i=0; i<n; i++)
                        {
                            gradf_Y1_faceI[i] = 0;
                        }
                        gradf_T1_faceI = 0;
                    }//end if(curPhaseState == 0)
                    else if(curPhaseState == 1)
                    {
                        //ph-1
                        calcCellGrad(faceOwn, nf, Y1_flatFld_diff, T1_flatFld_diff, alpha1_flatFld_diff, C_ph1_flatFld_diff, ownCells, 0.1*MIN_ALPHA_DIFF, 1, gradf_Y1_faceI, gradf_T1_faceI, n, debug, os);
                        //ph-0
                        for(i=0; i<n; i++)
                        {
                            gradf_Y0_faceI[i] = 0;
                        }
                        gradf_T0_faceI = 0;
                    }//end if(curPhaseState == 1)
                    else if(curPhaseState == 2)
                    {
                        //ph-1
                        calcCellGrad(faceOwn, nf, Y1_flatFld_diff, T1_flatFld_diff, alpha1_flatFld_diff, C_ph1_flatFld_diff, ownCells, 0.1*MIN_ALPHA_DIFF, 1, gradf_Y1_faceI, gradf_T1_faceI, n, debug, os);
                        //ph-0
                        calcCellGrad(faceOwn, nf, Y0_flatFld_diff, T0_flatFld_diff, alpha1_flatFld_diff, C_ph0_flatFld_diff, ownCells, 0.1*MIN_ALPHA_DIFF, 0, gradf_Y0_faceI, gradf_T0_faceI, n, debug, os);
                    }//end if(curPhaseState == 2)
                    else
                    {
                        for(i=0; i<n; i++)
                        {
                            gradf_Y1_faceI[i] = 0;
                            gradf_Y0_faceI[i] = 0;
                        }
                        gradf_T1_faceI = 0;
                        gradf_T0_faceI = 0;
                    }//end if(curPhaseState == 3)

                    for(i=0; i<n; i++)
                    {
                        gradf_Y1[i].boundaryField()[patchI][fcI] = 0;
                        gradf_Y0[i].boundaryField()[patchI][fcI] = 0;
                    }                                        
                    gradf_T1.boundaryField()[patchI][fcI] = gradf_T1_faceI;
                    gradf_T0.boundaryField()[patchI][fcI] = gradf_T0_faceI;

                    faceI++;
                }
            }
        }//end if(isA<zeroGradientFvPatchScalarField>(pY10))
        else if(isA<fixedValueFvPatchScalarField>(pY10))
        {            
            forAll(pY10, fcI)
            {
                nf = pSf[fcI]/pMagSf[fcI];
                faceOwn = own[faceI];
                const labelList& ownCells = diffCellStencil[faceOwn];
                curPhaseState = face_phaseState_diff[faceI];

                if(curPhaseState == 0)
                {   
                    //ph-0
                    calcCellGrad(faceOwn, nf, Y0_flatFld_diff, T0_flatFld_diff, alpha1_flatFld_diff, C_ph0_flatFld_diff, ownCells, 0.1*MIN_ALPHA_DIFF, 0, gradf_Y0_faceI, gradf_T0_faceI, n, debug, os);
                    //ph-1
                    for(i=0; i<n; i++)
                    {
                        gradf_Y1_faceI[i] = 0;
                    }
                    gradf_T1_faceI = 0;
                }//end if(curPhaseState == 0)
                else if(curPhaseState == 1)
                {
                    //ph-1
                    calcCellGrad(faceOwn, nf, Y1_flatFld_diff, T1_flatFld_diff, alpha1_flatFld_diff, C_ph1_flatFld_diff, ownCells, 0.1*MIN_ALPHA_DIFF, 1, gradf_Y1_faceI, gradf_T1_faceI, n, debug, os);
                    //ph-0
                    for(i=0; i<n; i++)
                    {
                        gradf_Y0_faceI[i] = 0;
                    }
                    gradf_T0_faceI = 0;
                }//end if(curPhaseState == 1)
                else if(curPhaseState == 2)
                {
                    //ph-1
                    calcCellGrad(faceOwn, nf, Y1_flatFld_diff, T1_flatFld_diff, alpha1_flatFld_diff, C_ph1_flatFld_diff, ownCells, 0.1*MIN_ALPHA_DIFF, 1, gradf_Y1_faceI, gradf_T1_faceI, n, debug, os);
                    //ph-0
                    calcCellGrad(faceOwn, nf, Y0_flatFld_diff, T0_flatFld_diff, alpha1_flatFld_diff, C_ph0_flatFld_diff, ownCells, 0.1*MIN_ALPHA_DIFF, 0, gradf_Y0_faceI, gradf_T0_faceI, n, debug, os);
                }//end if(curPhaseState == 2)
                else
                {
                    for(i=0; i<n; i++)
                    {
                        gradf_Y1_faceI[i] = 0;
                        gradf_Y0_faceI[i] = 0;
                    }
                    gradf_T1_faceI = 0;
                    gradf_T0_faceI = 0;
                }//end if(curPhaseState == 3)

                for(i=0; i<n; i++)
                {
                    gradf_Y1[i].boundaryField()[patchI][fcI] = gradf_Y1_faceI[i];
                    gradf_Y0[i].boundaryField()[patchI][fcI] = gradf_Y0_faceI[i];
                }

                if(isA<fixedValueFvPatchScalarField>(pT1))
                {                    
                    gradf_T1.boundaryField()[patchI][fcI] = gradf_T1_faceI;
                    gradf_T0.boundaryField()[patchI][fcI] = gradf_T0_faceI;
                }
                else
                {
                    gradf_T1.boundaryField()[patchI][fcI] = 0;
                    gradf_T0.boundaryField()[patchI][fcI] = 0;
                }

                faceI++;            
            }//end forAll(pY10, fcI)                        
        }//end if(isA<fixedValueFvPatchScalarField>(pY10))
        else
        {
            forAll(pY10, fcI)
            {                
                for(i=0; i<n; i++)
                {
                    gradf_Y1[i].boundaryField()[patchI][fcI] = 0;
                    gradf_Y0[i].boundaryField()[patchI][fcI] = 0;
                }
                gradf_T1.boundaryField()[patchI][fcI] = 0;
                gradf_T0.boundaryField()[patchI][fcI] = 0;                
            }
        }//end if(pp.coupled())        
    }//end forAll(Y1[0].boundaryField(), patchI)

    //end boundary faces
}


void calc_2ph_diffFluxes_Yi_Fick
(
    const fvMesh& mesh,
    const surfaceScalarField& rho1f,
    const surfaceScalarField& rho0f,
    const surfaceScalarField& D1fi,
    const surfaceScalarField& D0fi,
    const surfaceScalarField& gradf_Y1i,
    const surfaceScalarField& gradf_Y0i,
    const scalarField& magSf_ph1_own,
    const scalarField& magSf_ph0_own,
    const scalarField& magSf_ph1_nei,
    const scalarField& magSf_ph0_nei,
    const labelList& face_phaseState_diff,
    const volScalarField& rho1,
    const volScalarField& alpha1,
    const volScalarField& Y1i,
    const volScalarField& rho0,
    const volScalarField& alpha0,
    const volScalarField& Y0i,
    surfaceScalarField& diffFlux_Y1i,
    surfaceScalarField& diffFlux_Y0i,
    const scalar& Y1MIN,
    const scalar& Y1MAX,
    const scalar& Y0MIN,
    const scalar& Y0MAX,
    const label& i,
    const scalar& dt,
    const bool debug,
    OFstream& os
)
{
    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();    
    const surfaceScalarField& meshMagSf = mesh.magSf();    
    const scalarField& meshV = mesh.V();
    
    label faceI, faceOwn, faceNei, curPhaseState;
    scalar curMagSf_ph1, curMagSf_ph0, diffFlux_limiter;
    scalar VOwn, VNei, rho1Own, alpha1Own, Y1iOwn, rho1Nei, alpha1Nei, Y1iNei, rho0Own, alpha0Own, Y0iOwn, rho0Nei, alpha0Nei, Y0iNei;

    if(debug)
    {
        os<< "---------------------------------------------------------" << nl
            << "Diffusion flux calculation" << nl
            << "---------------------------------------------------------" << nl
            << nl
            << "Internal faces" << nl
            << endl;
    }

    //--------------------------------------------------------------//
    //Internal faces
    const scalarField& rho1Cells = rho1.internalField();
    const scalarField& alpha1Cells = alpha1.internalField();
    const scalarField& Y1iCells = Y1i.internalField();    

    const scalarField& rho0Cells = rho0.internalField();
    const scalarField& alpha0Cells = alpha0.internalField();
    const scalarField& Y0iCells = Y0i.internalField();    

    for(faceI=0; faceI<mesh.nInternalFaces(); faceI++)
    {
        curMagSf_ph1 = min(magSf_ph1_own[faceI], magSf_ph1_nei[faceI]);
        curMagSf_ph0 = min(magSf_ph0_own[faceI], magSf_ph0_nei[faceI]);
        curPhaseState = face_phaseState_diff[faceI];
        faceOwn = own[faceI];
        faceNei = nei[faceI];                        

        if(curPhaseState == 3)
        {            
            diffFlux_Y1i[faceI] = 0;
            diffFlux_Y0i[faceI] = 0;
        }
        else if(curPhaseState == 0)
        {
            diffFlux_Y1i[faceI] = 0;

            diffFlux_Y0i[faceI] = -rho0f[faceI]*D0fi[faceI]*curMagSf_ph0*gradf_Y0i[faceI];
            diffFlux_limiter = 1;
            calc_diffFlux_limiter2
            (
                rho0Cells[faceOwn],
                alpha0Cells[faceOwn],
                Y0iCells[faceOwn],
                meshV[faceOwn],
                rho0Cells[faceNei],
                alpha0Cells[faceNei],
                Y0iCells[faceNei],
                meshV[faceNei],
                dt,
                diffFlux_Y0i[faceI],
                diffFlux_limiter,
                Y0MIN,
                Y0MAX
            );
            diffFlux_Y0i[faceI] *= min(diffFlux_limiter, 1);

            if(diffFlux_limiter < 1)
            {
                if(debug)
                {
                    os<< "Face " << faceI << ": limiting diffFlux_Y0" << i 
                        << nl                    
                        << "Own " << faceOwn << "  Nei " << faceNei 
                        << nl
                        << "diffFlux_Y0" << i << ": " << diffFlux_Y0i[faceI] << "  diffFlux_limter: " << diffFlux_limiter 
                        << nl 
                        << "rho0Own: " << rho0Cells[faceOwn] << "  alpha0Own: " << alpha0Cells[faceOwn] << "  Y0iOwn: " <<  Y0iCells[faceOwn] << "  dt: " << dt << "  VOwn: " <<  meshV[faceOwn] 
                        << nl 
                        << "rho0Nei: " << rho0Cells[faceNei] << "  alpha0Nei: " << alpha0Cells[faceNei] << "  Y0iNei: " <<  Y0iCells[faceNei] << "  dt: " << dt << "  VNei: " <<  meshV[faceNei]
                        << nl
                        << endl;
                }
            }
        }
        else if(curPhaseState == 1)
        {            
            diffFlux_Y1i[faceI] = -rho1f[faceI]*D1fi[faceI]*curMagSf_ph1*gradf_Y1i[faceI];
            diffFlux_limiter = 1;
            calc_diffFlux_limiter2
            (
                rho1Cells[faceOwn],
                alpha1Cells[faceOwn],
                Y1iCells[faceOwn],
                meshV[faceOwn],
                rho1Cells[faceNei],
                alpha1Cells[faceNei],
                Y1iCells[faceNei],
                meshV[faceNei],
                dt,
                diffFlux_Y1i[faceI],
                diffFlux_limiter,
                Y1MIN,
                Y1MAX
            );
            diffFlux_Y1i[faceI] *= min(diffFlux_limiter, 1);

            if(diffFlux_limiter < 1)
            {
                if(debug)
                {
                    os<< "Face " << faceI << ": limiting diffFlux_Y1" << i 
                        << nl                    
                        << "Own " << faceOwn << "  Nei " << faceNei 
                        << nl
                        << "diffFlux_Y1" << i << ": " << diffFlux_Y1i[faceI] << "  diffFlux_limter: " << diffFlux_limiter 
                        << nl 
                        << "rho1Own: " << rho1Cells[faceOwn] << "  alpha1Own: " << alpha1Cells[faceOwn] << "  Y1iOwn: " <<  Y1iCells[faceOwn] << "  dt: " << dt << "  VOwn: " <<  meshV[faceOwn] 
                        << nl 
                        << "rho1Nei: " << rho1Cells[faceNei] << "  alpha1Nei: " << alpha1Cells[faceNei] << "  Y1iNei: " <<  Y1iCells[faceNei] << "  dt: " << dt << "  VNei: " <<  meshV[faceNei]
                        << nl
                        << endl;
                }            
            }

            diffFlux_Y0i[faceI] = 0;
        }
        else if(curPhaseState == 2)
        {            
            diffFlux_Y1i[faceI] = -rho1f[faceI]*D1fi[faceI]*curMagSf_ph1*gradf_Y1i[faceI];
            diffFlux_limiter = 1;
            calc_diffFlux_limiter2
            (
                rho1Cells[faceOwn],
                alpha1Cells[faceOwn],
                Y1iCells[faceOwn],
                meshV[faceOwn],
                rho1Cells[faceNei],
                alpha1Cells[faceNei],
                Y1iCells[faceNei],
                meshV[faceNei],
                dt,
                diffFlux_Y1i[faceI],
                diffFlux_limiter,
                Y1MIN,
                Y1MAX
            );
            diffFlux_Y1i[faceI] *= min(diffFlux_limiter, 1);

            if(diffFlux_limiter < 1)
            {
                if(debug)
                {
                    os<< "Face " << faceI << ": limiting diffFlux_Y1" << i 
                        << nl                    
                        << "Own " << faceOwn << "  Nei " << faceNei 
                        << nl
                        << "diffFlux_Y1" << i << ": " << diffFlux_Y1i[faceI] << "  diffFlux_limter: " << diffFlux_limiter 
                        << nl 
                        << "rho1Own: " << rho1Cells[faceOwn] << "  alpha1Own: " << alpha1Cells[faceOwn] << "  Y1iOwn: " <<  Y1iCells[faceOwn] << "  dt: " << dt << "  VOwn: " <<  meshV[faceOwn] 
                        << nl 
                        << "rho1Nei: " << rho1Cells[faceNei] << "  alpha1Nei: " << alpha1Cells[faceNei] << "  Y1iNei: " <<  Y1iCells[faceNei] << "  dt: " << dt << "  VNei: " <<  meshV[faceNei]
                        << nl
                        << endl;
                }            
            }

            diffFlux_Y0i[faceI] = -rho0f[faceI]*D0fi[faceI]*curMagSf_ph0*gradf_Y0i[faceI];
            diffFlux_limiter = 1;
            calc_diffFlux_limiter2
            (
                rho0Cells[faceOwn],
                alpha0Cells[faceOwn],
                Y0iCells[faceOwn],
                meshV[faceOwn],
                rho0Cells[faceNei],
                alpha0Cells[faceNei],
                Y0iCells[faceNei],
                meshV[faceNei],
                dt,
                diffFlux_Y0i[faceI],
                diffFlux_limiter,
                Y0MIN,
                Y0MAX
            );
            diffFlux_Y0i[faceI] *= min(diffFlux_limiter, 1);

            if(diffFlux_limiter < 1)
            {
                if(debug)
                {
                    os<< "Face " << faceI << ": limiting diffFlux_Y0" << i 
                        << nl                    
                        << "Own " << faceOwn << "  Nei " << faceNei 
                        << nl
                        << "diffFlux_Y0" << i << ": " << diffFlux_Y0i[faceI] << "  diffFlux_limter: " << diffFlux_limiter 
                        << nl 
                        << "rho0Own: " << rho0Cells[faceOwn] << "  alpha0Own: " << alpha0Cells[faceOwn] << "  Y0iOwn: " <<  Y0iCells[faceOwn] << "  dt: " << dt << "  VOwn: " <<  meshV[faceOwn] 
                        << nl 
                        << "rho0Nei: " << rho0Cells[faceNei] << "  alpha0Nei: " << alpha0Cells[faceNei] << "  Y0iNei: " <<  Y0iCells[faceNei] << "  dt: " << dt << "  VNei: " <<  meshV[faceNei]
                        << nl
                        << endl;
                }
            }
        }
        else
        {
            diffFlux_Y1i[faceI] = 0;
            diffFlux_Y0i[faceI] = 0;
        }

        if(debug)
        {
            os<< "Face: " << faceI << "  mag(Sf) = " << meshMagSf[faceI] << nl
                << "magSf_ph1_own = " << magSf_ph1_own[faceI] << "  alpha1f_own = " << magSf_ph1_own[faceI]/meshMagSf[faceI] 
                << "  magSf_ph1_nei = " << magSf_ph1_nei[faceI] << "  alpha1f_nei = " << magSf_ph1_nei[faceI]/meshMagSf[faceI] << nl
                << "magSf_ph0_own = " << magSf_ph0_own[faceI] << "  alpha0f_own = " << magSf_ph0_own[faceI]/meshMagSf[faceI] 
                << "  magSf_ph0_nei = " << magSf_ph0_nei[faceI] << "  alpha0f_nei = " << magSf_ph0_nei[faceI]/meshMagSf[faceI] << nl
                << "face phase state for diffusion flux calculation: " << curPhaseState << nl
                << "Own: " << faceOwn << "  Nei: " << faceNei << nl
                << "Own Y1" << i << ": " << Y1i[faceOwn] << "  Nei Y1" << i << ": " << Y1i[faceNei] << nl
                << "Own Y0" << i << ": " << Y0i[faceOwn] << "  Nei Y0" << i << ": " << Y0i[faceNei] << nl               
                << "diffusion flux ph1 for Y" << i << " = " << diffFlux_Y1i[faceI] << nl
                << "diffusion flux ph0 for Y" << i << " = " << diffFlux_Y0i[faceI] << nl
                << endl;
        }        
    }//end for(label faceI=0; faceI<mesh.nInternalFaces(); faceI++)

    //end internal faces    
    //--------------------------------------------------------------//

    //--------------------------------------------------------------//
    //Boundary faces
    if(debug)
    {
        os<< "Boundary faces" << nl
            << endl;
    }
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    //--------------------------------------------------------------//
    //Need volume of cell neighbour on neighbouring processor for 
    //each coupled patch face in order to calculate flux limiter
    //for that face
    const label nBnd = mesh.nFaces() - mesh.nInternalFaces();
    List<scalar> VNeiFld(nBnd);

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        if(pp.coupled())
        {
            faceI = pp.start();
            label bndFaceI = pp.start() - mesh.nInternalFaces();
            forAll(pp, fcI)
            {
                faceOwn = own[faceI];
                VNeiFld[bndFaceI] = meshV[faceOwn];
                faceI++;
                bndFaceI++;
            }
        }
    }

    syncTools::swapBoundaryFaceList(mesh, VNeiFld);
    //VNei now has cell volume of neighbouring cell for each coupled
    //patch face
    //--------------------------------------------------------------//
    
    forAll(Y1i.boundaryField(), patchI)
    {
        const polyPatch& pp = patches[patchI];
        const fvPatchScalarField& pY1i = Y1i.boundaryField()[patchI];
        const fvPatchScalarField& pY0i = Y0i.boundaryField()[patchI];
        fvsPatchScalarField& pdiffFlux_Y1i = diffFlux_Y1i.boundaryField()[patchI];
        const fvsPatchScalarField& pgradf_Y1i = gradf_Y1i.boundaryField()[patchI];
        const fvsPatchScalarField& prho1f = rho1f.boundaryField()[patchI];
        const fvsPatchScalarField& pD1fi = D1fi.boundaryField()[patchI];
        fvsPatchScalarField& pdiffFlux_Y0i = diffFlux_Y0i.boundaryField()[patchI];
        const fvsPatchScalarField& pgradf_Y0i = gradf_Y0i.boundaryField()[patchI];
        const fvsPatchScalarField& prho0f = rho0f.boundaryField()[patchI];
        const fvsPatchScalarField& pD0fi = D0fi.boundaryField()[patchI];
        const fvsPatchScalarField& pMagSf = meshMagSf.boundaryField()[patchI];        
        
        faceI = pp.start();

        if(pp.coupled())
        {
            label bndFaceI = pp.start() - mesh.nInternalFaces();

            const scalarField& rho1NeiFld = rho1.boundaryField()[patchI].patchNeighbourField();
            const scalarField& alpha1NeiFld = alpha1.boundaryField()[patchI].patchNeighbourField();
            const scalarField& Y1iNeiFld = Y1i.boundaryField()[patchI].patchNeighbourField();

            const scalarField& rho0NeiFld = rho0.boundaryField()[patchI].patchNeighbourField();
            const scalarField& alpha0NeiFld = alpha0.boundaryField()[patchI].patchNeighbourField();
            const scalarField& Y0iNeiFld = Y0i.boundaryField()[patchI].patchNeighbourField();
            
            forAll(pY1i, fcI)
            {                
                curMagSf_ph1 = min(magSf_ph1_own[faceI], magSf_ph1_nei[faceI]);
                curMagSf_ph0 = min(magSf_ph0_own[faceI], magSf_ph0_nei[faceI]);
                curPhaseState = face_phaseState_diff[faceI];
                faceOwn = own[faceI];
                VOwn = meshV[faceOwn];
                VNei = VNeiFld[bndFaceI];

                rho1Own = rho1Cells[faceOwn];
                alpha1Own = alpha1Cells[faceOwn];
                Y1iOwn = Y1iCells[faceOwn];
                rho1Nei = rho1NeiFld[fcI];
                alpha1Nei = alpha1NeiFld[fcI];
                Y1iNei = Y1iNeiFld[fcI];                

                rho0Own = rho0Cells[faceOwn];
                alpha0Own = alpha0Cells[faceOwn];
                Y0iOwn = Y0iCells[faceOwn];
                rho0Nei = rho0NeiFld[fcI];
                alpha0Nei = alpha0NeiFld[fcI];
                Y0iNei = Y0iNeiFld[fcI];

                if(curPhaseState == 3)
                {
                    pdiffFlux_Y1i[fcI] = 0;
                    pdiffFlux_Y0i[fcI] = 0;
                }
                if(curPhaseState == 0)
                {
                    pdiffFlux_Y1i[fcI] = 0;
                    pdiffFlux_Y0i[fcI] = -prho0f[fcI]*pD0fi[fcI]*curMagSf_ph0*pgradf_Y0i[fcI];
                    diffFlux_limiter = 1;
                    calc_diffFlux_limiter2
                    (
                        rho0Own,
                        alpha0Own,
                        Y0iOwn,
                        VOwn,
                        rho0Nei,
                        alpha0Nei,
                        Y0iNei,
                        VNei,
                        dt,
                        pdiffFlux_Y0i[fcI],
                        diffFlux_limiter,
                        Y0MIN,
                        Y0MAX
                    );
                    pdiffFlux_Y0i[fcI] *= min(diffFlux_limiter, 1);

                    if(diffFlux_limiter < 1)
                    {
                        if(debug)
                        {
                            os<< "Face " << faceI << ": limiting diffFlux_Y0" << i 
                                << nl                            
                                << "Own " << faceOwn 
                                << nl
                                << "diffFlux_Y0" << i << ": " << pdiffFlux_Y0i[fcI] << "  diffFlux_limter: " << diffFlux_limiter 
                                << nl 
                                << "rho0Own: " << rho0Own << "  alpha0Own: " << alpha0Own << "  Y0iOwn: " <<  Y0iOwn << "  dt: " << dt << "  VOwn: " <<  VOwn 
                                << nl 
                                << "rho0Nei: " << rho0Nei << "  alpha0Nei: " << alpha0Nei << "  Y0iNei: " <<  Y0iNei << "  dt: " << dt << "  VNei: " <<  VNei 
                                << nl
                                << endl;
                        }            
                    }
                }
                if(curPhaseState == 1)
                {
                    pdiffFlux_Y1i[fcI] = -prho1f[fcI]*pD1fi[fcI]*curMagSf_ph1*pgradf_Y1i[fcI];
                    diffFlux_limiter = 1;
                    calc_diffFlux_limiter2
                    (
                        rho1Own,
                        alpha1Own,
                        Y1iOwn,
                        VOwn,
                        rho1Nei,
                        alpha1Nei,
                        Y1iNei,
                        VNei,
                        dt,
                        pdiffFlux_Y1i[fcI],
                        diffFlux_limiter,
                        Y1MIN,
                        Y1MAX
                    );
                    pdiffFlux_Y1i[fcI] *= min(diffFlux_limiter, 1);
                    pdiffFlux_Y0i[fcI] = 0;

                    if(diffFlux_limiter < 1)
                    {
                        if(debug)
                        {
                            os<< "Face " << faceI << ": limiting diffFlux_Y1" << i 
                                << nl
                                << "Own " << faceOwn 
                                << nl
                                << "diffFlux_Y1" << i << ": " << pdiffFlux_Y1i[fcI] << "  diffFlux_limter: " << diffFlux_limiter 
                                << nl 
                                << "rho1Own: " << rho1Own << "  alpha1Own: " << alpha1Own << "  Y1iOwn: " <<  Y1iOwn << "  dt: " << dt << "  VOwn: " <<  VOwn 
                                << nl 
                                << "rho1Nei: " << rho1Nei << "  alpha1Nei: " << alpha1Nei << "  Y1iNei: " <<  Y1iNei << "  dt: " << dt << "  VNei: " <<  VNei 
                                << nl
                                << endl;
                        }            
                    }
                }
                if(curPhaseState == 2)
                {
                    pdiffFlux_Y1i[fcI] = -prho1f[fcI]*pD1fi[fcI]*curMagSf_ph1*pgradf_Y1i[fcI];                    
                    diffFlux_limiter = 1;
                    calc_diffFlux_limiter2
                    (
                        rho1Own,
                        alpha1Own,
                        Y1iOwn,
                        VOwn,
                        rho1Nei,
                        alpha1Nei,
                        Y1iNei,
                        VNei,
                        dt,
                        pdiffFlux_Y1i[fcI],
                        diffFlux_limiter,
                        Y1MIN,
                        Y1MAX
                    );
                    pdiffFlux_Y1i[fcI] *= min(diffFlux_limiter, 1);

                    if(diffFlux_limiter < 1)
                    {
                        if(debug)
                        {
                            os<< "Face " << faceI << ": limiting diffFlux_Y1" << i 
                                << nl
                                << "Own " << faceOwn 
                                << nl
                                << "diffFlux_Y1" << i << ": " << pdiffFlux_Y1i[fcI] << "  diffFlux_limter: " << diffFlux_limiter 
                                << nl 
                                << "rho1Own: " << rho1Own << "  alpha1Own: " << alpha1Own << "  Y1iOwn: " <<  Y1iOwn << "  dt: " << dt << "  VOwn: " <<  VOwn 
                                << nl 
                                << "rho1Nei: " << rho1Nei << "  alpha1Nei: " << alpha1Nei << "  Y1iNei: " <<  Y1iNei << "  dt: " << dt << "  VNei: " <<  VNei 
                                << nl
                                << endl;
                        }            
                    }

                    pdiffFlux_Y0i[fcI] = -prho0f[fcI]*pD0fi[fcI]*curMagSf_ph0*pgradf_Y0i[fcI];
                    diffFlux_limiter = 1;
                    calc_diffFlux_limiter2
                    (
                        rho0Own,
                        alpha0Own,
                        Y0iOwn,
                        VOwn,
                        rho0Nei,
                        alpha0Nei,
                        Y0iNei,
                        VNei,
                        dt,
                        pdiffFlux_Y0i[fcI],
                        diffFlux_limiter,
                        Y0MIN,
                        Y0MAX
                    );
                    pdiffFlux_Y0i[fcI] *= min(diffFlux_limiter, 1);

                    if(diffFlux_limiter < 1)
                    {
                        if(debug)
                        {
                            os<< "Face " << faceI << ": limiting diffFlux_Y0" << i 
                                << nl                            
                                << "Own " << faceOwn 
                                << nl
                                << "diffFlux_Y0" << i << ": " << pdiffFlux_Y0i[fcI] << "  diffFlux_limter: " << diffFlux_limiter 
                                << nl 
                                << "rho0Own: " << rho0Own << "  alpha0Own: " << alpha0Own << "  Y0iOwn: " <<  Y0iOwn << "  dt: " << dt << "  VOwn: " <<  VOwn 
                                << nl 
                                << "rho0Nei: " << rho0Nei << "  alpha0Nei: " << alpha0Nei << "  Y0iNei: " <<  Y0iNei << "  dt: " << dt << "  VNei: " <<  VNei 
                                << nl
                                << endl;
                        }            
                    }
                }
                else
                {
                    pdiffFlux_Y1i[fcI] = 0;
                    pdiffFlux_Y0i[fcI] = 0;
                }

                if(debug)
                {
                    os<< "Face: " << faceI << "  mag(Sf) = " << pMagSf[fcI] << nl
                        << "magSf_ph1_own = " << magSf_ph1_own[faceI] << "  alpha1f_own = " << magSf_ph1_own[faceI]/pMagSf[fcI] 
                        << "  magSf_ph1_nei = " << magSf_ph1_nei[faceI] << "  alpha1f_nei = " << magSf_ph1_nei[faceI]/pMagSf[fcI] << nl
                        << "magSf_ph0_own = " << magSf_ph0_own[faceI] << "  alpha0f_own = " << magSf_ph0_own[faceI]/pMagSf[fcI] 
                        << "  magSf_ph0_nei = " << magSf_ph0_nei[faceI] << "  alpha0f_nei = " << magSf_ph0_nei[faceI]/pMagSf[fcI] << nl
                        << "face phase state for diffusion flux calculation: " << curPhaseState << nl
                        << "Own: " << faceOwn << nl
                        << "Own Y1" << i << ": " << Y1iOwn << "  Nei Y1" << i << ": " << Y1iNei << nl
                        << "Own Y0" << i << ": " << Y0iOwn << "  Nei Y0" << i << ": " << Y0iNei << nl               
                        << "diffusion flux ph1 for Y" << i << " = " << pdiffFlux_Y1i[fcI] << nl
                        << "diffusion flux ph0 for Y" << i << " = " << pdiffFlux_Y0i[fcI] << nl
                        << endl;
                }
                
                faceI++;
                bndFaceI++;
            }//end forAll(pY1i, fcI)
        }//end if(pp.coupled())        
        else if(isA<fixedValueFvPatchScalarField>(pY1i))
        {
            forAll(pY1i, fcI)
            {                
                curMagSf_ph1 = magSf_ph1_own[faceI];
                curMagSf_ph0 = magSf_ph0_own[faceI];
                curPhaseState = face_phaseState_diff[faceI];
                faceOwn = own[faceI];                
                VOwn = meshV[faceOwn];

                rho1Own = rho1Cells[faceOwn];
                alpha1Own = alpha1Cells[faceOwn];
                Y1iOwn = Y1iCells[faceOwn];

                rho0Own = rho0Cells[faceOwn];
                alpha0Own = alpha0Cells[faceOwn];
                Y0iOwn = Y0iCells[faceOwn];                

                if(curPhaseState == 3)
                {
                    pdiffFlux_Y1i[fcI] = 0;
                    pdiffFlux_Y0i[fcI] = 0;
                }
                if(curPhaseState == 0)
                {
                    pdiffFlux_Y1i[fcI] = 0;
                    pdiffFlux_Y0i[fcI] = -prho0f[fcI]*pD0fi[fcI]*curMagSf_ph0*pgradf_Y0i[fcI];
                    diffFlux_limiter = 1;
                    calc_diffFlux_limiter
                    (
                        rho0Own,
                        alpha0Own,
                        Y0iOwn,
                        VOwn,
                        rho0Own,
                        alpha0Own,
                        Y0iOwn,
                        VOwn,
                        dt,
                        pdiffFlux_Y0i[fcI],
                        diffFlux_limiter,
                        Y0MIN,
                        Y0MAX
                    );
                    pdiffFlux_Y0i[fcI] *= min(diffFlux_limiter, 1);
                }
                if(curPhaseState == 1)
                {
                    pdiffFlux_Y1i[fcI] = -prho1f[fcI]*pD1fi[fcI]*curMagSf_ph1*pgradf_Y1i[fcI];
                    diffFlux_limiter = 1;
                    calc_diffFlux_limiter
                    (
                        rho1Own,
                        alpha1Own,
                        Y1iOwn,
                        VOwn,
                        rho1Own,
                        alpha1Own,
                        Y1iOwn,
                        VOwn,
                        dt,
                        pdiffFlux_Y1i[fcI],
                        diffFlux_limiter,
                        Y1MIN,
                        Y1MAX
                    );
                    pdiffFlux_Y1i[fcI] *= min(diffFlux_limiter, 1);
                    pdiffFlux_Y0i[fcI] = 0;
                }
                if(curPhaseState == 2)
                {
                    pdiffFlux_Y1i[fcI] = -prho1f[fcI]*pD1fi[fcI]*curMagSf_ph1*pgradf_Y1i[fcI];
                    pdiffFlux_Y0i[fcI] = -prho0f[fcI]*pD0fi[fcI]*curMagSf_ph0*pgradf_Y0i[fcI];
                }
                else
                {
                    pdiffFlux_Y1i[fcI] = 0;
                    pdiffFlux_Y0i[fcI] = 0;
                }

                if(debug)
                {
                    os<< "Face: " << faceI << "  mag(Sf) = " << pMagSf[faceI] << nl
                        << "magSf_ph1_own = " << magSf_ph1_own[faceI] << "  alpha1f_own = " << magSf_ph1_own[faceI]/pMagSf[fcI] 
                        << "  magSf_ph1_nei = " << magSf_ph1_nei[faceI] << "  alpha1f_nei = " << magSf_ph1_nei[faceI]/pMagSf[fcI] << nl
                        << "magSf_ph0_own = " << magSf_ph0_own[faceI] << "  alpha0f_own = " << magSf_ph0_own[faceI]/pMagSf[fcI] 
                        << "  magSf_ph0_nei = " << magSf_ph0_nei[faceI] << "  alpha0f_nei = " << magSf_ph0_nei[faceI]/pMagSf[fcI] << nl
                        << "face phase state for diffusion flux calculation: " << curPhaseState << nl
                        << "Own: " << faceOwn << nl
                        << "Own Y1" << i << ": " << pY1i[fcI] << nl
                        << "Own Y0" << i << ": " << pY0i[fcI] << nl               
                        << "diffusion flux ph1 for Y" << i << " = " << pdiffFlux_Y1i[fcI] << nl
                        << "diffusion flux ph0 for Y" << i << " = " << pdiffFlux_Y0i[fcI] << nl
                        << endl;
                }
                
                faceI++;
            }//end forAll(pY1i, fcI)
        }//end else if(isA<fixedValueFvPatchScalarField>(pY1i))
        else
        {
            forAll(pY1i, fcI)
            {
                pdiffFlux_Y1i[fcI] = 0;
                pdiffFlux_Y0i[fcI] = 0;
                faceI++;
            }
        }
    }//end forAll(patches,patchI)

    //end boundary faces

    if(debug)
    {
        os<<"---------------------------------------------" << nl
            << "Diffusive fluxes at faces for Y" << i << nl
            << "--------------------------------------------" << nl
            << endl;

        os<< "-----------------------------------------------------------------------------------------------------" << nl
            << " Face        face_phaseState            diffFlux_Y1i           diffFlux_Y0i" << nl
            << "-----------------------------------------------------------------------------------------------------" << nl
            << endl;
        
        for(faceI=0; faceI<mesh.nInternalFaces(); faceI++)
        {
            os<< "   " << faceI << "          " << face_phaseState_diff[faceI] << "                " << diffFlux_Y1i[faceI] << "                " << diffFlux_Y0i[faceI] << endl;
        }

        os<< "------------------------------------------------------------------------------------------------------" << nl
            << endl;
    }
}


void init_MS_flux
(
    label faceI,
    const PtrList<surfaceScalarField>& xf,
    const PtrList<surfaceScalarField>& gradf_x,    
    const PtrList<surfaceScalarField>& Dijf,
    int n,
    double *x,
    double *grad_x,
    double *Dij
)
{
    int idx;

    for(int i=0; i<n; i++)
    {
        x[i] = xf[i].internalField()[faceI];
        grad_x[i] = gradf_x[i].internalField()[faceI];

        for(int j=0; j<n; j++)
        {
            idx = i + j*n;
            Dij[idx] = Dijf[idx].internalField()[faceI];
        }
    }
}


void init_MS_flux
(
    label patchI,
    label faceI,
    const PtrList<surfaceScalarField>& xf,
    const PtrList<surfaceScalarField>& gradf_x,    
    const PtrList<surfaceScalarField>& Dijf,
    int n,
    double *x,
    double *grad_x,
    double *Dij
)
{
    int idx;

    for(int i=0; i<n; i++)
    {
        x[i] = xf[i].boundaryField()[patchI][faceI];
        grad_x[i] = gradf_x[i].boundaryField()[patchI][faceI];

        for(int j=0; j<n; j++)
        {
            idx = i + j*n;
            Dij[idx] = Dijf[idx].boundaryField()[patchI][faceI];
        }
    }
}


void calc_MS_flux
(
    double P,
    double T,
    double V,
    double *x,
    double R_gas,
    int n,
    double *Pc,
    double *Tc,
    double *w,
    double *MW,
    double *kij_T,
    double Ta_kij,
    double Tb_kij,
    int nT_kij,
    bool calcKij,
    double *kij,
    double *lnphi,
    double *dlnphi_dxj,
    double *grad_x,
    double *Dij,
    int n_flux_type,
    LPT_UMFPACK flux_umf,
    double *rhs_flux,
    double *flux_m
)
{
    int i, j, idx;
    double Z;

    Z = P*V/(R_gas*T);

    if(calcKij) calc_kij_from_table(T,n,Ta_kij,Tb_kij,nT_kij,kij_T,kij);
    
    fugacities_n_its_derivatives2_(&P,&T,&n,Pc,Tc,w,x,kij,lnphi,dlnphi_dxj);

    // rhs of Maxwell-Stefan equation
    for(i=0; i<n; i++)
    {
        rhs_flux[i] = -grad_x[i];

        if(n_flux_type > 0)
        {
            for(j=0; j<n; j++)
            {
                if(x[j]<1e-9) continue;

                idx = i*n + j;
                rhs_flux[i] = -x[i]*dlnphi_dxj[idx]*grad_x[j];
            }
        }

        rhs_flux[i] *= 1.0/V;    // 1/V[i] = c[i] mol/m^3
    }

    Maxwell_Stefan_mass_flux(Z,n,MW,x,Dij,rhs_flux,flux_m,flux_umf);
}


void calc_2ph_diffFluxes_Y_MS
(
    const fvMesh& mesh,
    const surfaceScalarField& Pf,
    const surfaceScalarField& T1f,
    const surfaceScalarField& T0f,
    const surfaceScalarField& V1f,
    const surfaceScalarField& V0f,
    double R_gas,
    int n,
    double *Pc,
    double *Tc,
    double *w,
    double *MW,
    double *kij_T,
    double Ta_kij,
    double Tb_kij,
    int nT_kij,
    bool calcKij,
    double *kij,
    const PtrList<surfaceScalarField>& Dij1f,
    const PtrList<surfaceScalarField>& Dij0f,
    const PtrList<surfaceScalarField>& x1f,
    const PtrList<surfaceScalarField>& x0f,
    const PtrList<surfaceScalarField>& gradf_x1,
    const PtrList<surfaceScalarField>& gradf_x0,    
    const scalarField& magSf_ph1_own,
    const scalarField& magSf_ph0_own,
    const scalarField& magSf_ph1_nei,
    const scalarField& magSf_ph0_nei,
    const labelList& face_phaseState_diff,
    const volScalarField& rho1,
    const volScalarField& alpha1,
    const PtrList<volScalarField>& Y1,
    const volScalarField& rho0,
    const volScalarField& alpha0,
    const PtrList<volScalarField>& Y0,
    PtrList<surfaceScalarField>& diffFlux_Y1,
    PtrList<surfaceScalarField>& diffFlux_Y0,
    LPT_UMFPACK flux_umf,
    int n_flux_type,
    scalar deltaT,
    const List<scalar>& Y1MIN,
    const List<scalar>& Y1MAX,
    const List<scalar>& Y0MIN,
    const List<scalar>& Y0MAX,
    const bool debug,
    OFstream& os
)
{
    int i;
    double P_tmp; double T1_tmp; double T0_tmp;
    double V1_tmp; double V0_tmp;

    double *x_tmp; double *grad_x_tmp; double *Dij_tmp;
    double *rhs_flux_tmp; double *flux_m_ph1; double *flux_m_ph0;
    double *lnphi_tmp; double *dlnphi_dxj_tmp;

    _NNEW_(x_tmp, double, n);
    _NNEW_(grad_x_tmp, double, n);
    _NNEW_(Dij_tmp, double, n*n);
    _NNEW_(rhs_flux_tmp, double, n);
    _NNEW_(flux_m_ph1, double, n);
    _NNEW_(flux_m_ph0, double, n);
    _NNEW_(lnphi_tmp, double, n);
    _NNEW_(dlnphi_dxj_tmp, double, n*n);

    label faceI, faceOwn, faceNei, curPhaseState, bndFaceI, nBnd;
    scalar curMagSf, curMagSf_ph1, curMagSf_ph0;
    scalar diffFlux_limiter_ph1, diffFlux_limiter_ph1i, diffFlux_limiter_ph0, diffFlux_limiter_ph0i;
    scalar VOwn, VNei, rho1Own, alpha1Own, rho1Nei, alpha1Nei, rho0Own, alpha0Own, rho0Nei, alpha0Nei;

    List<scalar> curDiffFlux_Y1(n);
    List<scalar> curDiffFlux_Y0(n);
    List<scalar> Y1Own(n);
    List<scalar> Y1Nei(n);
    List<scalar> Y0Own(n);
    List<scalar> Y0Nei(n);

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();    
    const surfaceScalarField& meshMagSf = mesh.magSf();    
    const scalarField& meshV = mesh.V();            

    if(debug)
    {
        os<< "---------------------------------------------------------" << nl
            << "Diffusion flux calculation" << nl
            << "---------------------------------------------------------" << nl
            << nl
            << "Internal faces" << nl
            << endl;
    }

    //--------------------------------------------------------------//
    //Internal faces
    const scalarField& rho1Cells = rho1.internalField();
    const scalarField& alpha1Cells = alpha1.internalField();

    const scalarField& rho0Cells = rho0.internalField();
    const scalarField& alpha0Cells = alpha0.internalField();

    List<scalarField> Y1Cells(n);
    List<scalarField> Y0Cells(n);
    for(i=0; i<n; i++)
    {
        Y1Cells[i] = Y1[i].internalField();
        Y0Cells[i] = Y0[i].internalField();
    }

    List<scalarField> Y1NeiFld(n);
    List<scalarField> Y0NeiFld(n);    

    //start loop over internal faces
    for(faceI=0; faceI<mesh.nInternalFaces(); faceI++)
    {
        curMagSf_ph1 = min(magSf_ph1_own[faceI], magSf_ph1_nei[faceI]);
        curMagSf_ph0 = min(magSf_ph0_own[faceI], magSf_ph0_nei[faceI]);
        curPhaseState = face_phaseState_diff[faceI];
        faceOwn = own[faceI];
        faceNei = nei[faceI];

        P_tmp = Pf[faceI];
        T1_tmp = T1f[faceI];
        T0_tmp = T0f[faceI];        
        V1_tmp = V1f[faceI];
        V0_tmp = V0f[faceI];

        if(curPhaseState == 3)
        {
            diffFlux_limiter_ph1 = 1;
            diffFlux_limiter_ph0 = 1;

            for(i=0; i<n; i++)
            {
                diffFlux_Y1[i][faceI] = 0;
                diffFlux_Y0[i][faceI] = 0;                
            }            
        }
        else if(curPhaseState == 0)
        {
            diffFlux_limiter_ph1 = 1;
            diffFlux_limiter_ph0 = 1;

            //ph-0
            init_MS_flux(faceI,x0f,gradf_x0,Dij0f,n,x_tmp,grad_x_tmp,Dij_tmp);
            calc_MS_flux(P_tmp,T0_tmp,V0_tmp,x_tmp,R_gas,n,Pc,Tc,w,MW,kij_T,Ta_kij,Tb_kij,nT_kij,calcKij,kij,lnphi_tmp,dlnphi_dxj_tmp,grad_x_tmp,Dij_tmp,n_flux_type,flux_umf,rhs_flux_tmp,flux_m_ph0);

            for(i=0; i<n; i++)
            {
                curDiffFlux_Y1[i] = 0;
                curDiffFlux_Y0[i] = curMagSf_ph0*flux_m_ph0[i];

                calc_diffFlux_limiter2
                (
                    rho0Cells[faceOwn],
                    alpha0Cells[faceOwn],
                    Y0Cells[i][faceOwn],
                    meshV[faceOwn],
                    rho0Cells[faceNei],
                    alpha0Cells[faceNei],
                    Y0Cells[i][faceNei],
                    meshV[faceNei],
                    deltaT,
                    curDiffFlux_Y0[i],
                    diffFlux_limiter_ph0i,
                    Y0MIN[i],
                    Y0MAX[i]
                );
                diffFlux_limiter_ph0 = min(diffFlux_limiter_ph0, diffFlux_limiter_ph0i);
            }
                        
            for(i=0; i<n; i++) 
            {
                diffFlux_Y1[i][faceI] = 0;
                diffFlux_Y0[i][faceI] = curDiffFlux_Y0[i]*min(diffFlux_limiter_ph0, 1);
            }
        }
        else if(curPhaseState == 1)
        {
            diffFlux_limiter_ph1 = 1;
            diffFlux_limiter_ph0 = 1;

            //ph-1
            init_MS_flux(faceI,x1f,gradf_x1,Dij1f,n,x_tmp,grad_x_tmp,Dij_tmp);
            calc_MS_flux(P_tmp,T1_tmp,V1_tmp,x_tmp,R_gas,n,Pc,Tc,w,MW,kij_T,Ta_kij,Tb_kij,nT_kij,calcKij,kij,lnphi_tmp,dlnphi_dxj_tmp,grad_x_tmp,Dij_tmp,n_flux_type,flux_umf,rhs_flux_tmp,flux_m_ph1);

            for(i=0; i<n; i++)
            {
                curDiffFlux_Y1[i] = curMagSf_ph1*flux_m_ph1[i];
                curDiffFlux_Y0[i] = 0;

                calc_diffFlux_limiter2
                (
                    rho1Cells[faceOwn],
                    alpha1Cells[faceOwn],
                    Y1Cells[i][faceOwn],
                    meshV[faceOwn],
                    rho1Cells[faceNei],
                    alpha1Cells[faceNei],
                    Y1Cells[i][faceNei],
                    meshV[faceNei],
                    deltaT,
                    curDiffFlux_Y1[i],
                    diffFlux_limiter_ph1i,
                    Y1MIN[i],
                    Y1MAX[i]
                );
                diffFlux_limiter_ph1 = min(diffFlux_limiter_ph1, diffFlux_limiter_ph1i);
            }
                        
            for(i=0; i<n; i++)
            {
                diffFlux_Y1[i][faceI] = curDiffFlux_Y1[i]*min(diffFlux_limiter_ph1, 1);
                diffFlux_Y0[i][faceI] = 0;
            }
        }
        else if(curPhaseState == 2)
        {
            diffFlux_limiter_ph1 = 1;
            diffFlux_limiter_ph0 = 1;

            //ph-1
            init_MS_flux(faceI,x1f,gradf_x1,Dij1f,n,x_tmp,grad_x_tmp,Dij_tmp);
            calc_MS_flux(P_tmp,T1_tmp,V1_tmp,x_tmp,R_gas,n,Pc,Tc,w,MW,kij_T,Ta_kij,Tb_kij,nT_kij,calcKij,kij,lnphi_tmp,dlnphi_dxj_tmp,grad_x_tmp,Dij_tmp,n_flux_type,flux_umf,rhs_flux_tmp,flux_m_ph1);

            //ph-0
            init_MS_flux(faceI,x0f,gradf_x0,Dij0f,n,x_tmp,grad_x_tmp,Dij_tmp);
            calc_MS_flux(P_tmp,T0_tmp,V0_tmp,x_tmp,R_gas,n,Pc,Tc,w,MW,kij_T,Ta_kij,Tb_kij,nT_kij,calcKij,kij,lnphi_tmp,dlnphi_dxj_tmp,grad_x_tmp,Dij_tmp,n_flux_type,flux_umf,rhs_flux_tmp,flux_m_ph0);

            for(i=0; i<n; i++)
            {
                curDiffFlux_Y1[i] = curMagSf_ph1*flux_m_ph1[i];
                curDiffFlux_Y0[i] = curMagSf_ph0*flux_m_ph0[i];

                calc_diffFlux_limiter2
                (
                    rho1Cells[faceOwn],
                    alpha1Cells[faceOwn],
                    Y1Cells[i][faceOwn],
                    meshV[faceOwn],
                    rho1Cells[faceNei],
                    alpha1Cells[faceNei],
                    Y1Cells[i][faceNei],
                    meshV[faceNei],
                    deltaT,
                    curDiffFlux_Y1[i],
                    diffFlux_limiter_ph1i,
                    Y1MIN[i],
                    Y1MAX[i]
                );
                diffFlux_limiter_ph1 = min(diffFlux_limiter_ph1, diffFlux_limiter_ph1i);

                calc_diffFlux_limiter2
                (
                    rho0Cells[faceOwn],
                    alpha0Cells[faceOwn],
                    Y0Cells[i][faceOwn],
                    meshV[faceOwn],
                    rho0Cells[faceNei],
                    alpha0Cells[faceNei],
                    Y0Cells[i][faceNei],
                    meshV[faceNei],
                    deltaT,
                    curDiffFlux_Y0[i],
                    diffFlux_limiter_ph0i,
                    Y0MIN[i],
                    Y0MAX[i]
                );
                diffFlux_limiter_ph0 = min(diffFlux_limiter_ph0, diffFlux_limiter_ph0i);
            }

            for(i=0; i<n; i++) 
            {
                diffFlux_Y1[i][faceI] = curDiffFlux_Y1[i]*min(diffFlux_limiter_ph1, 1);
                diffFlux_Y0[i][faceI] = curDiffFlux_Y0[i]*min(diffFlux_limiter_ph0, 1);
            }
        }
        else
        {
            diffFlux_limiter_ph1 = 1;
            diffFlux_limiter_ph0 = 1;

            for(i=0; i<n; i++)
            {
                diffFlux_Y1[i][faceI] = 0;
                diffFlux_Y0[i][faceI] = 0;                
            }
        }        

        if(debug)
        {
            os<< "Face: " << faceI << "  mag(Sf) = " << meshMagSf[faceI] << nl
                << "magSf_ph1_own = " << magSf_ph1_own[faceI] << "  alpha1f_own = " << magSf_ph1_own[faceI]/meshMagSf[faceI] 
                << "  magSf_ph1_nei = " << magSf_ph1_nei[faceI] << "  alpha1f_nei = " << magSf_ph1_nei[faceI]/meshMagSf[faceI] << nl
                << "magSf_ph0_own = " << magSf_ph0_own[faceI] << "  alpha0f_own = " << magSf_ph0_own[faceI]/meshMagSf[faceI] 
                << "  magSf_ph0_nei = " << magSf_ph0_nei[faceI] << "  alpha0f_nei = " << magSf_ph0_nei[faceI]/meshMagSf[faceI] << nl
                << "face phase state for diffusion flux calculation: " << curPhaseState << nl
                << "Vown = " << meshV[faceOwn] << "  VNei = " << meshV[faceNei] << "  dt = " << deltaT << nl
                << "--------------------------------------------------------------------------------------" << nl
                << "phase-1" << nl
                << "--------------------------------------------------------------------------------------" << nl
                << "Own: " << faceOwn << "  alpha1Own = " << alpha1Cells[faceOwn] << "  rho1Own = " << rho1Cells[faceOwn] 
                << endl;
            for(i=0; i<n; i++)
            {
                os<< "Y1[" << i << "] = " << Y1Cells[i][faceOwn] << "  ";
            }
            os<< nl
                << "--------------------------------------------------------------------------------------" << nl
                << "Nei: " << faceNei << "  alpha1Nei = " << alpha1Cells[faceNei] << "  rho1Nei = " << rho1Cells[faceNei]                
                << endl;
            for(i=0; i<n; i++)
            {
                os<< "Y1[" << i << "] = " << Y1Cells[i][faceNei] << "  ";
            }
            os<< nl
                << "--------------------------------------------------------------------------------------" 
                << endl;
            for(i=0; i<n; i++)
            {
                os<< "diffFlux_Y1[" << i << "] = " << diffFlux_Y1[i][faceI] << "  ";
            }
            os<< nl
                << "--------------------------------------------------------------------------------------" << nl
                << "diffFlux_limiter_ph1 = " << diffFlux_limiter_ph1 << nl
                << "--------------------------------------------------------------------------------------" << nl
                << "--------------------------------------------------------------------------------------" << nl
                << "phase-0" << nl
                << "--------------------------------------------------------------------------------------" << nl
                << "Own: " << faceOwn << "  alpha0Own = " << alpha0Cells[faceOwn] << "  rho0Own = " << rho0Cells[faceOwn] 
                << endl;
            for(i=0; i<n; i++)
            {
                os<< "Y0[" << i << "] = " << Y0Cells[i][faceOwn] << "  ";
            }
            os<< nl
                << "--------------------------------------------------------------------------------------" << nl
                << "Nei: " << faceNei << "  alpha0Nei = " << alpha0Cells[faceNei] << "  rho0Nei = " << rho0Cells[faceNei] 
                << endl;
            for(i=0; i<n; i++)
            {
                os<< "Y0[" << i << "] = " << Y0Cells[i][faceNei] << "  ";
            }
            os<< nl
                << "--------------------------------------------------------------------------------------" 
                << endl;
            for(i=0; i<n; i++)
            {
                os<< "diffFlux_Y0[" << i << "] = " << diffFlux_Y0[i][faceI] << "  ";
            }
            os<< nl
                << "--------------------------------------------------------------------------------------" << nl
                << "diffFlux_limiter_ph0 = " << diffFlux_limiter_ph0 << nl
                << "--------------------------------------------------------------------------------------" << nl
                << endl;
        }                
    }//end loop over internal faces 
     //for(label faceI=0; faceI<mesh.nInternalFaces(); faceI++)

    //end internal faces    
    //--------------------------------------------------------------//

    //--------------------------------------------------------------//
    //Boundary faces
    if(debug)
    {
        os<< "Boundary faces" << nl
            << endl;
    }
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    //--------------------------------------------------------------//
    //Need volume of cell neighbour on neighbouring processor for 
    //each coupled patch face in order to calculate flux limiter
    //for that face
    nBnd = mesh.nFaces() - mesh.nInternalFaces();
    List<scalar> VNeiFld(nBnd);

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        if(pp.coupled())
        {
            faceI = pp.start();
            bndFaceI = pp.start() - mesh.nInternalFaces();
            forAll(pp, fcI)
            {
                faceOwn = own[faceI];
                VNeiFld[bndFaceI] = meshV[faceOwn];
                faceI++;
                bndFaceI++;
            }
        }
    }

    syncTools::swapBoundaryFaceList(mesh, VNeiFld);
    //VNei now has cell volume of neighbouring cell for each coupled
    //patch face
    //--------------------------------------------------------------//
    
    forAll(alpha1.boundaryField(), patchI)
    {
        const polyPatch& pp = patches[patchI];
        const fvsPatchScalarField& pMagSf = meshMagSf.boundaryField()[patchI];

        const fvsPatchScalarField& pPf = Pf.boundaryField()[patchI];

        //ph-1        
        const fvsPatchScalarField& pT1f = T1f.boundaryField()[patchI];
        const fvsPatchScalarField& pV1f = V1f.boundaryField()[patchI];
        const fvPatchScalarField& pY10 = Y1[0].boundaryField()[patchI];
        //ph-0        
        const fvsPatchScalarField& pT0f = T0f.boundaryField()[patchI];
        const fvsPatchScalarField& pV0f = V0f.boundaryField()[patchI];        
        
        faceI = pp.start();

        if(pp.coupled())
        {
            bndFaceI = pp.start() - mesh.nInternalFaces();
            
            const scalarField& rho1NeiFld = rho1.boundaryField()[patchI].patchNeighbourField();
            const scalarField& alpha1NeiFld = alpha1.boundaryField()[patchI].patchNeighbourField();

            const scalarField& rho0NeiFld = rho0.boundaryField()[patchI].patchNeighbourField();
            const scalarField& alpha0NeiFld = alpha0.boundaryField()[patchI].patchNeighbourField();
            
            for(i=0; i<n; i++)
            {                
                Y1NeiFld[i] = Y1[i].boundaryField()[patchI].patchNeighbourField();
                Y0NeiFld[i] = Y0[i].boundaryField()[patchI].patchNeighbourField();
            }
            
            //start loop over all faces of coupled patch
            forAll(pY10, fcI)
            {                
                curMagSf_ph1 = min(magSf_ph1_own[faceI], magSf_ph1_nei[faceI]);
                curMagSf_ph0 = min(magSf_ph0_own[faceI], magSf_ph0_nei[faceI]);
                curPhaseState = face_phaseState_diff[faceI];
                faceOwn = own[faceI];
                VOwn = meshV[faceOwn];
                VNei = VNeiFld[bndFaceI];

                rho1Own = rho1Cells[faceOwn];
                alpha1Own = alpha1Cells[faceOwn];
                rho1Nei = rho1NeiFld[fcI];
                alpha1Nei = alpha1NeiFld[fcI];

                rho0Own = rho0Cells[faceOwn];
                alpha0Own = alpha0Cells[faceOwn];
                rho0Nei = rho0NeiFld[fcI];
                alpha0Nei = alpha0NeiFld[fcI];

                for(i=0; i<n; i++)
                {
                    Y1Own[i] = Y1Cells[i][faceOwn];
                    Y1Nei[i] = Y1NeiFld[i][fcI];
                    Y0Own[i] = Y0Cells[i][faceOwn];
                    Y0Nei[i] = Y0NeiFld[i][fcI];
                }

                P_tmp = pPf[fcI];
                T1_tmp = pT1f[fcI];
                T0_tmp = pT0f[fcI];        
                V1_tmp = pV1f[fcI];
                V0_tmp = pV0f[fcI];

                if(curPhaseState == 3)
                {
                    diffFlux_limiter_ph1 = 1;
                    diffFlux_limiter_ph0 = 1;

                    for(i=0; i<n; i++)
                    {
                        diffFlux_Y1[i].boundaryField()[patchI][fcI] = 0;
                        diffFlux_Y0[i].boundaryField()[patchI][fcI] = 0;                
                    }                    
                }
                if(curPhaseState == 0)
                {
                    diffFlux_limiter_ph1 = 1;
                    diffFlux_limiter_ph0 = 1;

                    //ph-0
                    init_MS_flux(patchI,fcI,x0f,gradf_x0,Dij0f,n,x_tmp,grad_x_tmp,Dij_tmp);
                    calc_MS_flux(P_tmp,T0_tmp,V0_tmp,x_tmp,R_gas,n,Pc,Tc,w,MW,kij_T,Ta_kij,Tb_kij,nT_kij,calcKij,kij,lnphi_tmp,dlnphi_dxj_tmp,grad_x_tmp,Dij_tmp,n_flux_type,flux_umf,rhs_flux_tmp,flux_m_ph0);

                    for(i=0; i<n; i++)
                    {
                        curDiffFlux_Y1[i] = 0;
                        curDiffFlux_Y0[i] = curMagSf_ph0*flux_m_ph0[i];

                        calc_diffFlux_limiter2
                        (
                            rho0Own,
                            alpha0Own,
                            Y0Own[i],
                            VOwn,
                            rho0Nei,
                            alpha0Nei,
                            Y0Nei[i],
                            VNei,
                            deltaT,
                            curDiffFlux_Y0[i],
                            diffFlux_limiter_ph0i,
                            Y0MIN[i],
                            Y0MAX[i]
                        );
                        diffFlux_limiter_ph0 = min(diffFlux_limiter_ph0, diffFlux_limiter_ph0i);
                    }
                        
                    for(i=0; i<n; i++) 
                    {
                        diffFlux_Y1[i].boundaryField()[patchI][fcI] = 0;
                        diffFlux_Y0[i].boundaryField()[patchI][fcI] = curDiffFlux_Y0[i]*min(diffFlux_limiter_ph0, 1);
                    }
                }
                if(curPhaseState == 1)
                {
                    diffFlux_limiter_ph1 = 1;
                    diffFlux_limiter_ph0 = 1;

                    //ph-1
                    init_MS_flux(patchI,fcI,x1f,gradf_x1,Dij1f,n,x_tmp,grad_x_tmp,Dij_tmp);
                    calc_MS_flux(P_tmp,T1_tmp,V1_tmp,x_tmp,R_gas,n,Pc,Tc,w,MW,kij_T,Ta_kij,Tb_kij,nT_kij,calcKij,kij,lnphi_tmp,dlnphi_dxj_tmp,grad_x_tmp,Dij_tmp,n_flux_type,flux_umf,rhs_flux_tmp,flux_m_ph1);

                    for(i=0; i<n; i++)
                    {
                        curDiffFlux_Y1[i] = curMagSf_ph1*flux_m_ph1[i];
                        curDiffFlux_Y0[i] = 0;                        

                        calc_diffFlux_limiter2
                        (
                            rho1Own,
                            alpha1Own,
                            Y1Own[i],
                            VOwn,
                            rho1Nei,
                            alpha1Nei,
                            Y1Nei[i],
                            VNei,
                            deltaT,
                            curDiffFlux_Y1[i],
                            diffFlux_limiter_ph1i,
                            Y1MIN[i],
                            Y1MAX[i]
                        );
                        diffFlux_limiter_ph1 = min(diffFlux_limiter_ph1, diffFlux_limiter_ph1i);
                    }
                        
                    for(i=0; i<n; i++) 
                    {                        
                        diffFlux_Y1[i].boundaryField()[patchI][fcI] = curDiffFlux_Y1[i]*min(diffFlux_limiter_ph1, 1);
                        diffFlux_Y0[i].boundaryField()[patchI][fcI] = 0;
                    }
                }
                if(curPhaseState == 2)
                {
                    diffFlux_limiter_ph1 = 1;
                    diffFlux_limiter_ph0 = 1;

                    //ph-1
                    init_MS_flux(patchI,fcI,x1f,gradf_x1,Dij1f,n,x_tmp,grad_x_tmp,Dij_tmp);
                    calc_MS_flux(P_tmp,T1_tmp,V1_tmp,x_tmp,R_gas,n,Pc,Tc,w,MW,kij_T,Ta_kij,Tb_kij,nT_kij,calcKij,kij,lnphi_tmp,dlnphi_dxj_tmp,grad_x_tmp,Dij_tmp,n_flux_type,flux_umf,rhs_flux_tmp,flux_m_ph1);

                    //ph-0
                    init_MS_flux(patchI,fcI,x0f,gradf_x0,Dij0f,n,x_tmp,grad_x_tmp,Dij_tmp);
                    calc_MS_flux(P_tmp,T0_tmp,V0_tmp,x_tmp,R_gas,n,Pc,Tc,w,MW,kij_T,Ta_kij,Tb_kij,nT_kij,calcKij,kij,lnphi_tmp,dlnphi_dxj_tmp,grad_x_tmp,Dij_tmp,n_flux_type,flux_umf,rhs_flux_tmp,flux_m_ph0);

                    for(i=0; i<n; i++)
                    {
                        curDiffFlux_Y1[i] = curMagSf_ph1*flux_m_ph1[i];
                        curDiffFlux_Y0[i] = curMagSf_ph0*flux_m_ph0[i];

                        calc_diffFlux_limiter2
                        (
                            rho1Own,
                            alpha1Own,
                            Y1Own[i],
                            VOwn,
                            rho1Nei,
                            alpha1Nei,
                            Y1Nei[i],
                            VNei,
                            deltaT,
                            curDiffFlux_Y1[i],
                            diffFlux_limiter_ph1i,
                            Y1MIN[i],
                            Y1MAX[i]
                        );
                        diffFlux_limiter_ph1 = min(diffFlux_limiter_ph1, diffFlux_limiter_ph1i);

                        calc_diffFlux_limiter2
                        (
                            rho0Own,
                            alpha0Own,
                            Y0Own[i],
                            VOwn,
                            rho0Nei,
                            alpha0Nei,
                            Y0Nei[i],
                            VNei,
                            deltaT,
                            curDiffFlux_Y0[i],
                            diffFlux_limiter_ph0i,
                            Y0MIN[i],
                            Y0MAX[i]
                        );
                        diffFlux_limiter_ph0 = min(diffFlux_limiter_ph0, diffFlux_limiter_ph0i);
                    }
                        
                    for(i=0; i<n; i++) 
                    {                        
                        diffFlux_Y1[i].boundaryField()[patchI][fcI] = curDiffFlux_Y1[i]*min(diffFlux_limiter_ph1, 1);
                        diffFlux_Y0[i].boundaryField()[patchI][fcI] = curDiffFlux_Y0[i]*min(diffFlux_limiter_ph0, 1);
                    }
                }
                else
                {
                    diffFlux_limiter_ph1 = 1;
                    diffFlux_limiter_ph0 = 1;

                    for(i=0; i<n; i++)
                    {
                        diffFlux_Y1[i].boundaryField()[patchI][fcI] = 0;
                        diffFlux_Y0[i].boundaryField()[patchI][fcI] = 0;                
                    }                    
                }
                
                if(debug)
                {
                    curMagSf = pMagSf[fcI];
                    os<< "Face: " << faceI << "  mag(Sf) = " << curMagSf << nl
                        << "magSf_ph1_own = " << magSf_ph1_own[faceI] << "  alpha1f_own = " << magSf_ph1_own[faceI]/curMagSf 
                        << "  magSf_ph1_nei = " << magSf_ph1_nei[faceI] << "  alpha1f_nei = " << magSf_ph1_nei[faceI]/curMagSf << nl
                        << "magSf_ph0_own = " << magSf_ph0_own[faceI] << "  alpha0f_own = " << magSf_ph0_own[faceI]/curMagSf 
                        << "  magSf_ph0_nei = " << magSf_ph0_nei[faceI] << "  alpha0f_nei = " << magSf_ph0_nei[faceI]/curMagSf << nl
                        << "face phase state for diffusion flux calculation: " << curPhaseState << nl
                        << "Vown = " << VOwn << "  VNei = " << VNei << "  dt = " << deltaT << nl
                        << "--------------------------------------------------------------------------------------" << nl
                        << "phase-1" << nl
                        << "--------------------------------------------------------------------------------------" << nl
                        << "Own: " << faceOwn << "  alpha1Own = " << alpha1Own << "  rho1Own = " << rho1Own 
                        << endl;
                    for(i=0; i<n; i++)
                    {
                        os<< "Y1[" << i << "] = " << Y1Own[i] << "  ";
                    }
                    os<< nl
                        << "--------------------------------------------------------------------------------------" << nl
                        << "Nei: " << "NA" << "  alpha1Nei = " << alpha1Nei << "  rho1Nei = " << rho1Nei                
                        << endl;
                    for(i=0; i<n; i++)
                    {
                        os<< "Y1[" << i << "] = " << Y1Nei[i] << "  ";
                    }
                    os<< nl
                        << "--------------------------------------------------------------------------------------" 
                        << endl;
                    for(i=0; i<n; i++)
                    {
                        os<< "diffFlux_Y1[" << i << "] = " << diffFlux_Y1[i].boundaryField()[patchI][fcI] << "  ";
                    }
                    os<< nl
                        << "--------------------------------------------------------------------------------------" << nl
                        << "diffFlux_limiter_ph1 = " << diffFlux_limiter_ph1 << nl
                        << "--------------------------------------------------------------------------------------" << nl
                        << "--------------------------------------------------------------------------------------" << nl
                        << "phase-0" << nl
                        << "--------------------------------------------------------------------------------------" << nl
                        << "Own: " << faceOwn << "  alpha0Own = " << alpha0Own << "  rho0Own = " << rho0Own 
                        << endl;
                    for(i=0; i<n; i++)
                    {
                        os<< "Y0[" << i << "] = " << Y0Own[i] << "  ";
                    }
                    os<< nl
                        << "--------------------------------------------------------------------------------------" << nl
                        << "Nei: " << "NA" << "  alpha0Nei = " << alpha0Nei << "  rho0Nei = " << rho0Nei                
                        << endl;
                    for(i=0; i<n; i++)
                    {
                        os<< "Y0[" << i << "] = " << Y0Nei[i] << "  ";
                    }
                    os<< nl
                        << "--------------------------------------------------------------------------------------" 
                        << endl;
                    for(i=0; i<n; i++)
                    {
                        os<< "diffFlux_Y0[" << i << "] = " << diffFlux_Y0[i].boundaryField()[patchI][fcI] << "  ";
                    }
                    os<< nl
                        << "--------------------------------------------------------------------------------------" << nl
                        << "diffFlux_limiter_ph0 = " << diffFlux_limiter_ph0 << nl
                        << "--------------------------------------------------------------------------------------" << nl
                        << "--------------------------------------------------------------------------------------" << nl
                        << endl;
                }

                faceI++;
                bndFaceI++;
            }//end loop over all faces for the coupled patch
             //forAll(alpha1Own, fcI)
        }//end if(pp.coupled())        
        else if(isA<fixedValueFvPatchScalarField>(pY10))
        {
            const fvPatchScalarField& rho1NeiFld = rho1.boundaryField()[patchI];
            const fvPatchScalarField& alpha1NeiFld = alpha1.boundaryField()[patchI];

            const fvPatchScalarField& rho0NeiFld = rho0.boundaryField()[patchI];
            const fvPatchScalarField& alpha0NeiFld = alpha0.boundaryField()[patchI];

            for(i=0; i<n; i++)
            {
                const scalarField& pY1i = Y1[i].boundaryField()[patchI];
                Y1NeiFld[i] = pY1i;
                const scalarField& pY0i = Y0[i].boundaryField()[patchI];
                Y0NeiFld[i] = pY0i;
            }

            forAll(pY10, fcI)
            {                
                curMagSf_ph1 = magSf_ph1_own[faceI];
                curMagSf_ph0 = magSf_ph0_own[faceI];
                curPhaseState = face_phaseState_diff[faceI];
                faceOwn = own[faceI];
                VOwn = meshV[faceOwn];                

                rho1Own = rho1Cells[faceOwn];
                alpha1Own = alpha1Cells[faceOwn];
                rho1Nei = rho1NeiFld[fcI];
                alpha1Nei = alpha1NeiFld[fcI];

                rho0Own = rho0Cells[faceOwn];
                alpha0Own = alpha0Cells[faceOwn];
                rho0Nei = rho0NeiFld[fcI];
                alpha0Nei = alpha0NeiFld[fcI];

                for(i=0; i<n; i++)
                {
                    Y1Own[i] = Y1Cells[i][faceOwn];
                    Y1Nei[i] = Y1NeiFld[i][fcI];
                    Y0Own[i] = Y0Cells[i][faceOwn];
                    Y0Nei[i] = Y0NeiFld[i][fcI];
                }

                P_tmp = pPf[fcI];
                T1_tmp = pT1f[fcI];
                T0_tmp = pT0f[fcI];        
                V1_tmp = pV1f[fcI];
                V0_tmp = pV0f[fcI];

                if(curPhaseState == 3)
                {
                    diffFlux_limiter_ph1 = 1;
                    diffFlux_limiter_ph0 = 1;

                    for(i=0; i<n; i++)
                    {
                        diffFlux_Y1[i].boundaryField()[patchI][fcI] = 0;
                        diffFlux_Y0[i].boundaryField()[patchI][fcI] = 0;                
                    }                    
                }
                if(curPhaseState == 0)
                {
                    diffFlux_limiter_ph1 = 1;
                    diffFlux_limiter_ph0 = 1;

                    //ph-0
                    init_MS_flux(patchI,fcI,x0f,gradf_x0,Dij0f,n,x_tmp,grad_x_tmp,Dij_tmp);
                    calc_MS_flux(P_tmp,T0_tmp,V0_tmp,x_tmp,R_gas,n,Pc,Tc,w,MW,kij_T,Ta_kij,Tb_kij,nT_kij,calcKij,kij,lnphi_tmp,dlnphi_dxj_tmp,grad_x_tmp,Dij_tmp,n_flux_type,flux_umf,rhs_flux_tmp,flux_m_ph0);

                    for(i=0; i<n; i++)
                    {
                        curDiffFlux_Y1[i] = 0;
                        curDiffFlux_Y0[i] = curMagSf_ph0*flux_m_ph0[i];

                        calc_diffFlux_limiter2
                        (
                            rho0Own,
                            alpha0Own,
                            Y0Own[i],
                            VOwn,
                            rho0Nei,
                            alpha0Nei,
                            Y0Nei[i],
                            VOwn,
                            deltaT,
                            curDiffFlux_Y0[i],
                            diffFlux_limiter_ph0i,
                            Y0MIN[i],
                            Y0MAX[i]
                        );
                        diffFlux_limiter_ph0 = min(diffFlux_limiter_ph0, diffFlux_limiter_ph0i);
                    }
                        
                    for(i=0; i<n; i++) 
                    {
                        diffFlux_Y1[i].boundaryField()[patchI][fcI] = 0;
                        diffFlux_Y0[i].boundaryField()[patchI][fcI] = curDiffFlux_Y0[i]*min(diffFlux_limiter_ph0, 1);
                    }                    
                }
                if(curPhaseState == 1)
                {
                    diffFlux_limiter_ph1 = 1;
                    diffFlux_limiter_ph0 = 1;

                    //ph-1
                    init_MS_flux(patchI,fcI,x1f,gradf_x1,Dij1f,n,x_tmp,grad_x_tmp,Dij_tmp);
                    calc_MS_flux(P_tmp,T1_tmp,V1_tmp,x_tmp,R_gas,n,Pc,Tc,w,MW,kij_T,Ta_kij,Tb_kij,nT_kij,calcKij,kij,lnphi_tmp,dlnphi_dxj_tmp,grad_x_tmp,Dij_tmp,n_flux_type,flux_umf,rhs_flux_tmp,flux_m_ph1);

                    for(i=0; i<n; i++)
                    {
                        curDiffFlux_Y1[i] = curMagSf_ph1*flux_m_ph1[i];
                        curDiffFlux_Y0[i] = 0;                        

                        calc_diffFlux_limiter2
                        (
                            rho1Own,
                            alpha1Own,
                            Y1Own[i],
                            VOwn,
                            rho1Nei,
                            alpha1Nei,
                            Y1Nei[i],
                            VOwn,
                            deltaT,
                            curDiffFlux_Y1[i],
                            diffFlux_limiter_ph1i,
                            Y1MIN[i],
                            Y1MAX[i]
                        );
                        diffFlux_limiter_ph1 = min(diffFlux_limiter_ph1, diffFlux_limiter_ph1i);
                    }
                        
                    for(i=0; i<n; i++) 
                    {                        
                        diffFlux_Y1[i].boundaryField()[patchI][fcI] = curDiffFlux_Y1[i]*min(diffFlux_limiter_ph1, 1);
                        diffFlux_Y0[i].boundaryField()[patchI][fcI] = 0;
                    }
                }
                if(curPhaseState == 2)
                {
                    diffFlux_limiter_ph1 = 1;
                    diffFlux_limiter_ph0 = 1;

                    //ph-1
                    init_MS_flux(patchI,fcI,x1f,gradf_x1,Dij1f,n,x_tmp,grad_x_tmp,Dij_tmp);
                    calc_MS_flux(P_tmp,T1_tmp,V1_tmp,x_tmp,R_gas,n,Pc,Tc,w,MW,kij_T,Ta_kij,Tb_kij,nT_kij,calcKij,kij,lnphi_tmp,dlnphi_dxj_tmp,grad_x_tmp,Dij_tmp,n_flux_type,flux_umf,rhs_flux_tmp,flux_m_ph1);

                    //ph-0
                    init_MS_flux(patchI,fcI,x0f,gradf_x0,Dij0f,n,x_tmp,grad_x_tmp,Dij_tmp);
                    calc_MS_flux(P_tmp,T0_tmp,V0_tmp,x_tmp,R_gas,n,Pc,Tc,w,MW,kij_T,Ta_kij,Tb_kij,nT_kij,calcKij,kij,lnphi_tmp,dlnphi_dxj_tmp,grad_x_tmp,Dij_tmp,n_flux_type,flux_umf,rhs_flux_tmp,flux_m_ph0);

                    for(i=0; i<n; i++)
                    {
                        curDiffFlux_Y1[i] = curMagSf_ph1*flux_m_ph1[i];
                        curDiffFlux_Y0[i] = curMagSf_ph0*flux_m_ph0[i];

                        calc_diffFlux_limiter2
                        (
                            rho1Own,
                            alpha1Own,
                            Y1Own[i],
                            VOwn,
                            rho1Nei,
                            alpha1Nei,
                            Y1Nei[i],
                            VOwn,
                            deltaT,
                            curDiffFlux_Y1[i],
                            diffFlux_limiter_ph1i,
                            Y1MIN[i],
                            Y1MAX[i]
                        );
                        diffFlux_limiter_ph1 = min(diffFlux_limiter_ph1, diffFlux_limiter_ph1i);

                        calc_diffFlux_limiter2
                        (
                            rho0Own,
                            alpha0Own,
                            Y0Own[i],
                            VOwn,
                            rho0Nei,
                            alpha0Nei,
                            Y0Nei[i],
                            VOwn,
                            deltaT,
                            curDiffFlux_Y0[i],
                            diffFlux_limiter_ph0i,
                            Y0MIN[i],
                            Y0MAX[i]
                        );
                        diffFlux_limiter_ph0 = min(diffFlux_limiter_ph0, diffFlux_limiter_ph0i);
                    }
                        
                    for(i=0; i<n; i++) 
                    {                        
                        diffFlux_Y1[i].boundaryField()[patchI][fcI] = curDiffFlux_Y1[i]*min(diffFlux_limiter_ph1, 1);
                        diffFlux_Y0[i].boundaryField()[patchI][fcI] = curDiffFlux_Y0[i]*min(diffFlux_limiter_ph0, 1);
                    }
                }
                else
                {
                    diffFlux_limiter_ph1 = 1;
                    diffFlux_limiter_ph0 = 1;

                    for(i=0; i<n; i++)
                    {
                        diffFlux_Y1[i].boundaryField()[patchI][fcI] = 0;
                        diffFlux_Y0[i].boundaryField()[patchI][fcI] = 0;                
                    }
                }

                if(debug)
                {
                    curMagSf = pMagSf[fcI];
                    os<< "Face: " << faceI << "  mag(Sf) = " << curMagSf << nl
                        << "magSf_ph1_own = " << magSf_ph1_own[faceI] << "  alpha1f_own = " << magSf_ph1_own[faceI]/curMagSf 
                        << "  magSf_ph1_nei = " << magSf_ph1_nei[faceI] << "  alpha1f_nei = " << magSf_ph1_nei[faceI]/curMagSf << nl
                        << "magSf_ph0_own = " << magSf_ph0_own[faceI] << "  alpha0f_own = " << magSf_ph0_own[faceI]/curMagSf 
                        << "  magSf_ph0_nei = " << magSf_ph0_nei[faceI] << "  alpha0f_nei = " << magSf_ph0_nei[faceI]/curMagSf << nl
                        << "face phase state for diffusion flux calculation: " << curPhaseState << nl
                        << "Vown = " << meshV[faceOwn] << "  VNei = " << meshV[faceOwn] << "  dt = " << deltaT << nl
                        << "--------------------------------------------------------------------------------------" << nl
                        << "phase-1" << nl
                        << "--------------------------------------------------------------------------------------" << nl
                        << "Own: " << faceOwn << "  alpha1Own = " << alpha1Own << "  rho1Own = " << rho1Own 
                        << endl;
                    for(i=0; i<n; i++)
                    {
                        os<< "Y1[" << i << "] = " << Y1Own[i] << "  ";
                    }
                    os<< nl
                        << "--------------------------------------------------------------------------------------" << nl
                        << "Nei: " << "NA" << "  alpha1Nei = " << alpha1Nei << "  rho1Nei = " << rho1Nei                
                        << endl;
                    for(i=0; i<n; i++)
                    {
                        os<< "Y1[" << i << "] = " << Y1Nei[i] << "  ";
                    }
                    os<< nl
                        << "--------------------------------------------------------------------------------------" 
                        << endl;
                    for(i=0; i<n; i++)
                    {
                        os<< "diffFlux_Y1[" << i << "] = " << diffFlux_Y1[i].boundaryField()[patchI][fcI] << "  ";
                    }
                    os<< nl
                        << "--------------------------------------------------------------------------------------" << nl
                        << "diffFlux_limiter_ph1 = " << diffFlux_limiter_ph1 << nl
                        << "--------------------------------------------------------------------------------------" << nl
                        << "--------------------------------------------------------------------------------------" << nl
                        << "phase-0" << nl
                        << "--------------------------------------------------------------------------------------" << nl
                        << "Own: " << faceOwn << "  alpha0Own = " << alpha0Own << "  rho0Own = " << rho0Own 
                        << endl;
                    for(i=0; i<n; i++)
                    {
                        os<< "Y0[" << i << "] = " << Y0Own[i] << "  ";
                    }
                    os<< nl
                        << "--------------------------------------------------------------------------------------" << nl
                        << "Nei: " << "NA" << "  alpha0Nei = " << alpha0Nei << "  rho0Nei = " << rho0Nei                
                        << endl;
                    for(i=0; i<n; i++)
                    {
                        os<< "Y0[" << i << "] = " << Y0Nei[i] << "  ";
                    }
                    os<< nl
                        << "--------------------------------------------------------------------------------------" 
                        << endl;
                    for(i=0; i<n; i++)
                    {
                        os<< "diffFlux_Y0[" << i << "] = " << diffFlux_Y0[i].boundaryField()[patchI][fcI] << "  ";
                    }
                    os<< nl
                        << "--------------------------------------------------------------------------------------" << nl
                        << "diffFlux_limiter_ph0 = " << diffFlux_limiter_ph0 << nl
                        << "--------------------------------------------------------------------------------------" << nl
                        << "--------------------------------------------------------------------------------------" << nl
                        << endl;
                }

                faceI++;
            }//end forAll(alpha1Own, fcI)
        }//end else if(isA<fixedValueFvPatchScalarField>(pY10))
        else
        {
            forAll(pY10, fcI)
            {
                for(i=0; i<n; i++)
                {
                    diffFlux_Y1[i].boundaryField()[patchI][fcI] = 0;
                    diffFlux_Y0[i].boundaryField()[patchI][fcI] = 0;                
                }
                faceI++;
            }
        }
    }//end loop over all boundary patches 
     //forAll(alpha1.boundaryField(),patchI)

    //end boundary faces

    if(debug)
    {
        os<<"------------------------------------------------------------------------------------------------------" << nl
            <<"------------------------------------------------------------------------------------------------------" << nl
            << "Species Diffusive fluxes at faces" << nl
            << "-----------------------------------------------------------------------------------------------------" << nl
            << endl;

        os<< "-----------------------------------------------------------------------------------------------------" << nl
            << " Face        face_phaseState            i        diffFlux_Y1[i]           diffFlux_Y0[i]" << nl
            << "-----------------------------------------------------------------------------------------------------" << nl
            << endl;
        
        for(faceI=0; faceI<mesh.nInternalFaces(); faceI++)
        {
            os<< "   " << faceI << "          " << face_phaseState_diff[faceI] << "                " << "0" << "        " << diffFlux_Y1[0][faceI] << "                " << diffFlux_Y0[0][faceI] << endl;
            for(i=0; i<n; i++)
            {
                os<< "   " << "  " << "          " << " " << "                " << i << "        " << diffFlux_Y1[i][faceI] << "                " << diffFlux_Y0[i][faceI] << endl;
            }
        }

        os<< "------------------------------------------------------------------------------------------------------" << nl
            << "------------------------------------------------------------------------------------------------------" << nl
            << endl;
    }

    _DDELETE_(x_tmp);
    _DDELETE_(grad_x_tmp);
    _DDELETE_(Dij_tmp);
    _DDELETE_(rhs_flux_tmp);
    _DDELETE_(flux_m_ph1);
    _DDELETE_(flux_m_ph0);
    _DDELETE_(lnphi_tmp);
    _DDELETE_(dlnphi_dxj_tmp);
}


void calc_2ph_diffFluxes_T
(
    const fvMesh& mesh,
    double R_gas,
    int n,
    double *Pc,
    double *Tc,
    double *w,
    double *MW,
    double *Tb,
    double *SG,
    double *H8,
    double *kij_T,
    double Ta_kij,
    double Tb_kij,
    int nT_kij,
    double *kij,
    double P,
    const surfaceScalarField& lambda1f,
    const surfaceScalarField& lambda0f,
    const surfaceScalarField& gradf_T1,
    const surfaceScalarField& gradf_T0,    
    const scalarField& magSf_ph1_own,
    const scalarField& magSf_ph0_own,
    const scalarField& magSf_ph1_nei,
    const scalarField& magSf_ph0_nei,
    const labelList& face_phaseState_diff,
    const volScalarField& rho1,
    const volScalarField& alpha1,
    const PtrList<volScalarField>& X1,
    const volScalarField& H1,
    const volScalarField& T1,
    const volScalarField& rho0,
    const volScalarField& alpha0,
    const PtrList<volScalarField>& X0,
    const volScalarField& H0,
    const volScalarField& T0,
    const scalar& deltaT,
    surfaceScalarField& diffFlux_T1,
    surfaceScalarField& diffFlux_T0,    
    const int& MAX_ITER_T,
    const double& T_TOL,
    const double& H_TOL,
    const scalar& T1MIN,
    const scalar& T1MAX,
    const scalar& T0MIN,
    const scalar& T0MAX,
    const bool debug,
    OFstream& os
)
{
    int i;
    label faceI, curPhaseState, faceOwn, faceNei, nBnd, bndFaceI;    
    scalar magSf_faceI, magSf_ph1_faceI, magSf_ph0_faceI, VOwn, VNei;
    scalar rho1Own, alpha1Own, T1Own, H1Own, rho0Own, alpha0Own, T0Own, H0Own;
    scalar rho1Nei, alpha1Nei, T1Nei, H1Nei, rho0Nei, alpha0Nei, T0Nei, H0Nei;
    scalar lambda1_faceI, lambda0_faceI, gradf_T1_faceI, gradf_T0_faceI;
    scalar diffFlux_T1_faceI, diffFlux_T0_faceI, diffFlux_limiter_ph1_faceI, diffFlux_limiter_ph0_faceI;
    List<scalar> x1Own(n); List<scalar> x0Own(n);
    List<scalar> x1Nei(n); List<scalar> x0Nei(n);

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();    
    const surfaceScalarField& meshMagSf = mesh.magSf();    
    const scalarField& meshV = mesh.V();

    const scalarField& rho1Cells = rho1.internalField();
    const scalarField& alpha1Cells = alpha1.internalField();
    const scalarField& T1Cells = T1.internalField();
    const scalarField& H1Cells = H1.internalField();

    const scalarField& rho0Cells = rho0.internalField();
    const scalarField& alpha0Cells = alpha0.internalField();
    const scalarField& T0Cells = T0.internalField();
    const scalarField& H0Cells = H0.internalField();

    List<scalarField> x1Cells(n);
    List<scalarField> x0Cells(n);
    for(i=0; i<n; i++)
    {
        x1Cells[i] = X1[i].internalField();
        x0Cells[i] = X0[i].internalField();
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    nBnd = mesh.nFaces() - mesh.nInternalFaces();
    List<scalar> VNei_bnd(nBnd);
    
    if(debug)
    {
        os<< "---------------------------------------------------------" << nl
            << "Diffusion flux calculation" << nl
            << "---------------------------------------------------------" << nl
            << nl
            << "Internal faces" << nl
            << endl;
    }

    //--------------------------------------------------------------//
    //Internal faces    
    for(faceI=0; faceI<mesh.nInternalFaces(); faceI++)
    {
        magSf_faceI = meshMagSf[faceI];
        magSf_ph1_faceI = min(magSf_ph1_own[faceI], magSf_ph1_nei[faceI]);
        magSf_ph0_faceI = min(magSf_ph0_own[faceI], magSf_ph0_nei[faceI]);
        curPhaseState = face_phaseState_diff[faceI];
        faceOwn = own[faceI];
        faceNei = nei[faceI];
        VOwn = meshV[faceOwn];
        VNei = meshV[faceNei];        
        rho1Own = rho1Cells[faceOwn];
        alpha1Own = alpha1Cells[faceOwn];
        T1Own = T1Cells[faceOwn];
        H1Own = H1Cells[faceOwn];
        rho0Own = rho0Cells[faceOwn];
        alpha0Own = alpha0Cells[faceOwn];
        T0Own = T0Cells[faceOwn];
        H0Own = H0Cells[faceOwn];
        rho1Nei = rho1Cells[faceNei];
        alpha1Nei = alpha1Cells[faceNei];
        T1Nei = T1Cells[faceNei];
        H1Nei = H1Cells[faceNei];
        rho0Nei = rho0Cells[faceNei];
        alpha0Nei = alpha0Cells[faceNei];
        T0Nei = T0Cells[faceNei];
        H0Nei = H0Cells[faceNei];
        for(i=0; i<n; i++)
        {
            x1Own[i] = x1Cells[i][faceOwn];
            x1Nei[i] = x1Cells[i][faceNei];
            x0Own[i] = x0Cells[i][faceOwn];
            x0Nei[i] = x0Cells[i][faceNei];
        }

        lambda1_faceI = lambda1f[faceI];
        lambda0_faceI = lambda0f[faceI];
        gradf_T1_faceI = gradf_T1[faceI];
        gradf_T0_faceI = gradf_T0[faceI];

        if(curPhaseState == 3)
        {
            diffFlux_T1_faceI = 0;
            diffFlux_T0_faceI = 0;
            diffFlux_limiter_ph1_faceI = 1;
            diffFlux_limiter_ph0_faceI = 1;
        }
        else if(curPhaseState == 0)
        {
            diffFlux_T0_faceI = -magSf_ph0_faceI*lambda0_faceI*gradf_T0_faceI;            
            calc_diffFlux_limiter_T(alpha0Own, rho0Own, x0Own, T0Own, H0Own, VOwn, alpha0Nei, rho0Nei, x0Nei, T0Nei, H0Nei, VNei, P, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, deltaT, diffFlux_T0_faceI, diffFlux_limiter_ph0_faceI, MAX_ITER_T, T_TOL, H_TOL, T0MIN, T0MAX);
            diffFlux_T0_faceI *= min(diffFlux_limiter_ph0_faceI, 1);

            diffFlux_limiter_ph1_faceI = 1;
            diffFlux_T1_faceI = 0;
        }
        else if(curPhaseState == 1)
        {
            diffFlux_T1_faceI = -magSf_ph1_faceI*lambda1_faceI*gradf_T1_faceI;
            calc_diffFlux_limiter_T(alpha1Own, rho1Own, x1Own, T1Own, H1Own, VOwn, alpha1Nei, rho1Nei, x1Nei, T1Nei, H1Nei, VNei, P, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, deltaT, diffFlux_T1_faceI, diffFlux_limiter_ph1_faceI, MAX_ITER_T, T_TOL, H_TOL, T1MIN, T1MAX);
            diffFlux_T1_faceI *= min(diffFlux_limiter_ph1_faceI, 1);

            diffFlux_limiter_ph0_faceI = 1;
            diffFlux_T0_faceI = 0;
        }
        else if(curPhaseState == 2)
        {
            diffFlux_T1_faceI = -magSf_ph1_faceI*lambda1_faceI*gradf_T1_faceI;
            calc_diffFlux_limiter_T(alpha1Own, rho1Own, x1Own, T1Own, H1Own, VOwn, alpha1Nei, rho1Nei, x1Nei, T1Nei, H1Nei, VNei, P, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, deltaT, diffFlux_T1_faceI, diffFlux_limiter_ph1_faceI, MAX_ITER_T, T_TOL, H_TOL, T1MIN, T1MAX);
            diffFlux_T1_faceI *= min(diffFlux_limiter_ph1_faceI, 1);

            diffFlux_T0_faceI = -magSf_ph0_faceI*lambda0_faceI*gradf_T0_faceI;            
            calc_diffFlux_limiter_T(alpha0Own, rho0Own, x0Own, T0Own, H0Own, VOwn, alpha0Nei, rho0Nei, x0Nei, T0Nei, H0Nei, VNei, P, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, deltaT, diffFlux_T0_faceI, diffFlux_limiter_ph0_faceI, MAX_ITER_T, T_TOL, H_TOL, T0MIN, T0MAX);
            diffFlux_T0_faceI *= min(diffFlux_limiter_ph0_faceI, 1);
        }
        else
        {
            diffFlux_T1_faceI = 0;
            diffFlux_T0_faceI = 0;
            diffFlux_limiter_ph1_faceI = 1;
            diffFlux_limiter_ph0_faceI = 1;
        }

        diffFlux_T1[faceI] = diffFlux_T1_faceI;
        diffFlux_T0[faceI] = diffFlux_T0_faceI;

        if(debug)
        {
            os<< "Face: " << faceI << "  mag(Sf) = " << magSf_faceI << nl
                << "magSf_ph1_own = " << magSf_ph1_own[faceI] << "  alpha1f_own = " << magSf_ph1_own[faceI]/magSf_faceI 
                << "  magSf_ph1_nei = " << magSf_ph1_nei[faceI] << "  alpha1f_nei = " << magSf_ph1_nei[faceI]/magSf_faceI << nl
                << "magSf_ph0_own = " << magSf_ph0_own[faceI] << "  alpha0f_own = " << magSf_ph0_own[faceI]/magSf_faceI 
                << "  magSf_ph0_nei = " << magSf_ph0_nei[faceI] << "  alpha0f_nei = " << magSf_ph0_nei[faceI]/magSf_faceI << nl
                << "face phase state for diffusion flux calculation: " << curPhaseState << nl
                << "Own: " << faceOwn << "  Nei: " << faceNei << nl
                << "T1Own = " << T1Own << "  T1Nei = " << T1Nei << nl
                << "T0Own = " << T0Own << "  T0Nei = " << T0Nei << nl
                << "diffFlux_T1 = " << diffFlux_T1_faceI << "  limiter_ph1 = " << diffFlux_limiter_ph1_faceI << nl
                << "diffFlux_T0 = " << diffFlux_T0_faceI << "  limiter_ph0 = " << diffFlux_limiter_ph0_faceI << nl
                << endl;
        }       
    }//end for(label faceI=0; faceI<mesh.nInternalFaces(); faceI++)

    //end internal faces    
    //--------------------------------------------------------------//

    //--------------------------------------------------------------//
    //Boundary faces
    if(debug)
    {
        os<< "Boundary faces" << nl
            << endl;
    }    
    //--------------------------------------------------------------//
    //Need volume of cell neighbour on neighbouring processor for 
    //each coupled patch face in order to calculate flux limiter
    //for that face
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        if(pp.coupled())
        {
            faceI = pp.start();
            bndFaceI = faceI - mesh.nInternalFaces();
            forAll(pp, fcI)
            {
                faceOwn = own[faceI];
                VNei_bnd[bndFaceI] = meshV[faceOwn];
                faceI++;
                bndFaceI++;
            }
        }
    }

    syncTools::swapBoundaryFaceList(mesh, VNei_bnd);
    //VNei now has cell volume of neighbouring cell for each coupled
    //patch face
    //--------------------------------------------------------------//
    
    forAll(T1.boundaryField(), patchI)
    {
        const polyPatch& pp = patches[patchI];
        const fvsPatchScalarField& pMagSf = meshMagSf.boundaryField()[patchI];
        const fvPatchScalarField& pT1 = T1.boundaryField()[patchI];
        const fvsPatchScalarField& pgradf_T1 = gradf_T1.boundaryField()[patchI];
        const fvsPatchScalarField& plambda1f = lambda1f.boundaryField()[patchI];
        fvsPatchScalarField& pdiffFlux_T1 = diffFlux_T1.boundaryField()[patchI];
        const fvsPatchScalarField& pgradf_T0 = gradf_T0.boundaryField()[patchI];
        const fvsPatchScalarField& plambda0f = lambda0f.boundaryField()[patchI];
        fvsPatchScalarField& pdiffFlux_T0 = diffFlux_T0.boundaryField()[patchI];
                
        faceI = pp.start();

        if(pp.coupled())
        {            
            const scalarField& rho1NeiFld = rho1.boundaryField()[patchI].patchNeighbourField();
            const scalarField& alpha1NeiFld = alpha1.boundaryField()[patchI].patchNeighbourField();
            const scalarField& T1NeiFld = T1.boundaryField()[patchI].patchNeighbourField();
            const scalarField& H1NeiFld = H1.boundaryField()[patchI].patchNeighbourField();

            const scalarField& rho0NeiFld = rho0.boundaryField()[patchI].patchNeighbourField();
            const scalarField& alpha0NeiFld = alpha0.boundaryField()[patchI].patchNeighbourField();
            const scalarField& T0NeiFld = T0.boundaryField()[patchI].patchNeighbourField();
            const scalarField& H0NeiFld = H0.boundaryField()[patchI].patchNeighbourField();

            List<scalarField> x1NeiFld(n);
            List<scalarField> x0NeiFld(n);
            for(i=0; i<n; i++)
            {               
                x1NeiFld[i] = X1[i].boundaryField()[patchI].patchNeighbourField();
                x0NeiFld[i] = X0[i].boundaryField()[patchI].patchNeighbourField();
            }
            
            forAll(pT1, fcI)
            {
                bndFaceI = faceI - mesh.nInternalFaces();
                magSf_faceI = pMagSf[fcI];
                magSf_ph1_faceI = min(magSf_ph1_own[faceI], magSf_ph1_nei[faceI]);
                magSf_ph0_faceI = min(magSf_ph0_own[faceI], magSf_ph0_nei[faceI]);
                curPhaseState = face_phaseState_diff[faceI];
                faceOwn = own[faceI];
                VOwn = meshV[faceOwn];
                VNei = VNei_bnd[bndFaceI];
                rho1Own = rho1Cells[faceOwn];
                alpha1Own = alpha1Cells[faceOwn];
                T1Own = T1Cells[faceOwn];
                H1Own = H1Cells[faceOwn];
                rho0Own = rho0Cells[faceOwn];
                alpha0Own = alpha0Cells[faceOwn];
                T0Own = T0Cells[faceOwn];
                H0Own = H0Cells[faceOwn];
                rho1Nei = rho1NeiFld[fcI];
                alpha1Nei = alpha1NeiFld[fcI];
                T1Nei = T1NeiFld[fcI];
                H1Nei = H1NeiFld[fcI];
                rho0Nei = rho0NeiFld[fcI];
                alpha0Nei = alpha0NeiFld[fcI];
                T0Nei = T0NeiFld[fcI];
                H0Nei = H0NeiFld[fcI];
                for(i=0; i<n; i++)
                {
                    x1Own[i] = x1Cells[i][faceOwn];
                    x1Nei[i] = x1NeiFld[i][fcI];
                    x0Own[i] = x1Cells[i][faceOwn];
                    x0Nei[i] = x0NeiFld[i][fcI];
                }

                lambda1_faceI = plambda1f[fcI];
                lambda0_faceI = plambda0f[fcI];
                gradf_T1_faceI = pgradf_T1[fcI];
                gradf_T0_faceI = pgradf_T0[fcI];

                if(curPhaseState == 3)
                {
                    diffFlux_T1_faceI = 0;
                    diffFlux_T0_faceI = 0;
                    diffFlux_limiter_ph1_faceI = 1;
                    diffFlux_limiter_ph0_faceI = 1;
                }
                if(curPhaseState == 0)
                {                    
                    diffFlux_T0_faceI = -magSf_ph0_faceI*lambda0_faceI*gradf_T0_faceI;
                    calc_diffFlux_limiter_T(alpha0Own, rho0Own, x0Own, T0Own, H0Own, VOwn, alpha0Nei, rho0Nei, x0Nei, T0Nei, H0Nei, VNei, P, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, deltaT, diffFlux_T0_faceI, diffFlux_limiter_ph0_faceI, MAX_ITER_T, T_TOL, H_TOL, T0MIN, T0MAX);
                    diffFlux_T0_faceI *= min(diffFlux_limiter_ph0_faceI, 1);

                    diffFlux_limiter_ph1_faceI = 1;
                    diffFlux_T1_faceI = 0;                    
                }
                if(curPhaseState == 1)
                {
                    diffFlux_T1_faceI = -magSf_ph1_faceI*lambda1_faceI*gradf_T1_faceI;
                    calc_diffFlux_limiter_T(alpha1Own, rho1Own, x1Own, T1Own, H1Own, VOwn, alpha1Nei, rho1Nei, x1Nei, T1Nei, H1Nei, VNei, P, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, deltaT, diffFlux_T1_faceI, diffFlux_limiter_ph1_faceI, MAX_ITER_T, T_TOL, H_TOL, T1MIN, T1MAX);            
                    diffFlux_T1_faceI *= min(diffFlux_limiter_ph1_faceI, 1);

                    diffFlux_limiter_ph0_faceI = 1;
                    diffFlux_T0_faceI = 0;
                }
                if(curPhaseState == 2)
                {
                    diffFlux_T1_faceI = -magSf_ph1_faceI*lambda1_faceI*gradf_T1_faceI;
                    calc_diffFlux_limiter_T(alpha1Own, rho1Own, x1Own, T1Own, H1Own, VOwn, alpha1Nei, rho1Nei, x1Nei, T1Nei, H1Nei, VNei, P, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, deltaT, diffFlux_T1_faceI, diffFlux_limiter_ph1_faceI, MAX_ITER_T, T_TOL, H_TOL, T1MIN, T1MAX);            
                    diffFlux_T1_faceI *= min(diffFlux_limiter_ph1_faceI, 1);

                    diffFlux_T0_faceI = -magSf_ph0_faceI*lambda0_faceI*gradf_T0_faceI;            
                    calc_diffFlux_limiter_T(alpha0Own, rho0Own, x0Own, T0Own, H0Own, VOwn, alpha0Nei, rho0Nei, x0Nei, T0Nei, H0Nei, VNei, P, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, deltaT, diffFlux_T0_faceI, diffFlux_limiter_ph0_faceI, MAX_ITER_T, T_TOL, H_TOL, T0MIN, T0MAX);            
                    diffFlux_T0_faceI *= min(diffFlux_limiter_ph0_faceI, 1);
                }
                else
                {
                    diffFlux_T1_faceI = 0;
                    diffFlux_T0_faceI = 0;
                    diffFlux_limiter_ph1_faceI = 1;
                    diffFlux_limiter_ph0_faceI = 1;
                }

                pdiffFlux_T1[fcI] = diffFlux_T1_faceI;
                pdiffFlux_T0[fcI] = diffFlux_T0_faceI;

                if(debug)
                {
                    os<< "Face: " << faceI << "  mag(Sf) = " << magSf_faceI << nl
                        << "magSf_ph1_own = " << magSf_ph1_own[faceI] << "  alpha1f_own = " << magSf_ph1_own[faceI]/magSf_faceI 
                        << "  magSf_ph1_nei = " << magSf_ph1_nei[faceI] << "  alpha1f_nei = " << magSf_ph1_nei[faceI]/magSf_faceI << nl
                        << "magSf_ph0_own = " << magSf_ph0_own[faceI] << "  alpha0f_own = " << magSf_ph0_own[faceI]/magSf_faceI 
                        << "  magSf_ph0_nei = " << magSf_ph0_nei[faceI] << "  alpha0f_nei = " << magSf_ph0_nei[faceI]/magSf_faceI << nl
                        << "face phase state for diffusion flux calculation: " << curPhaseState << nl
                        << "Own: " << faceOwn << nl
                        << "T1Own = " << T1Own << "  T1Nei = " << T1Nei << nl
                        << "T0Own = " << T0Own << "  T0Nei = " << T0Nei << nl
                        << "diffFlux_T1 = " << diffFlux_T1_faceI << "  limiter_ph1 = " << diffFlux_limiter_ph1_faceI << nl
                        << "diffFlux_T0 = " << diffFlux_T0_faceI << "  limiter_ph0 = " << diffFlux_limiter_ph0_faceI << nl
                        << endl;
                }
                
                faceI++;                
            }//end forAll(pT1, fcI)
        }//end if(pp.coupled())        
        else if(isA<fixedValueFvPatchScalarField>(pT1))
        {            
            forAll(pT1, fcI)
            {                
                magSf_faceI = pMagSf[fcI];
                magSf_ph1_faceI = magSf_ph1_own[faceI];
                magSf_ph0_faceI = magSf_ph0_own[faceI];
                curPhaseState = face_phaseState_diff[faceI];
                faceOwn = own[faceI];           
                VOwn = meshV[faceOwn];                
                rho1Own = rho1Cells[faceOwn];
                alpha1Own = alpha1Cells[faceOwn];
                T1Own = T1Cells[faceOwn];
                H1Own = H1Cells[faceOwn];
                rho0Own = rho0Cells[faceOwn];
                alpha0Own = alpha0Cells[faceOwn];
                T0Own = T0Cells[faceOwn];
                H0Own = H0Cells[faceOwn];                
                for(i=0; i<n; i++)
                {
                    x1Own[i] = x1Cells[i][faceOwn];
                    x0Own[i] = x0Cells[i][faceOwn];                    
                }

                lambda1_faceI = plambda1f[fcI];
                lambda0_faceI = plambda0f[fcI];
                gradf_T1_faceI = pgradf_T1[fcI];
                gradf_T0_faceI = pgradf_T0[fcI];

                if(curPhaseState == 3)
                {
                    diffFlux_T1_faceI = 0;
                    diffFlux_T0_faceI = 0;
                    diffFlux_limiter_ph1_faceI = 1;
                    diffFlux_limiter_ph0_faceI = 1;
                }
                if(curPhaseState == 0)
                {                    
                    diffFlux_T0_faceI = -magSf_ph0_faceI*lambda0_faceI*gradf_T0_faceI;
                    calc_diffFlux_limiter2_T(alpha0Own, rho0Own, x0Own, T0Own, H0Own, VOwn, P, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, deltaT, diffFlux_T0_faceI, diffFlux_limiter_ph0_faceI, T0MIN, T0MAX);            
                    diffFlux_T0_faceI *= min(diffFlux_limiter_ph0_faceI, 1);

                    diffFlux_limiter_ph1_faceI = 1;
                    diffFlux_T1_faceI = 0;                    
                }
                if(curPhaseState == 1)
                {
                    diffFlux_T1_faceI = -magSf_ph1_faceI*lambda1_faceI*gradf_T1_faceI;
                    calc_diffFlux_limiter2_T(alpha1Own, rho1Own, x1Own, T1Own, H1Own, VOwn, P, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, deltaT, diffFlux_T1_faceI, diffFlux_limiter_ph1_faceI, T1MIN, T1MAX);
                    diffFlux_T1_faceI *= min(diffFlux_limiter_ph1_faceI, 1);

                    diffFlux_limiter_ph0_faceI = 1;
                    diffFlux_T0_faceI = 0;
                }
                if(curPhaseState == 2)
                {
                    diffFlux_T1_faceI = -magSf_ph1_faceI*lambda1_faceI*gradf_T1_faceI;
                    calc_diffFlux_limiter2_T(alpha1Own, rho1Own, x1Own, T1Own, H1Own, VOwn, P, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, deltaT, diffFlux_T1_faceI, diffFlux_limiter_ph1_faceI, T1MIN, T1MAX);            
                    diffFlux_T1_faceI *= min(diffFlux_limiter_ph1_faceI, 1);

                    diffFlux_T0_faceI = -magSf_ph0_faceI*lambda0_faceI*gradf_T0_faceI;            
                    calc_diffFlux_limiter2_T(alpha0Own, rho0Own, x0Own, T0Own, H0Own, VOwn, P, n, Pc, Tc, w, MW, Tb, SG, H8, kij_T, Ta_kij, Tb_kij, nT_kij, kij, deltaT, diffFlux_T0_faceI, diffFlux_limiter_ph0_faceI, T0MIN, T0MAX);            
                    diffFlux_T0_faceI *= min(diffFlux_limiter_ph0_faceI, 1);
                }
                else
                {
                    diffFlux_T1_faceI = 0;
                    diffFlux_T0_faceI = 0;
                    diffFlux_limiter_ph1_faceI = 1;
                    diffFlux_limiter_ph0_faceI = 1;
                }

                pdiffFlux_T1[fcI] = diffFlux_T1_faceI;
                pdiffFlux_T0[fcI] = diffFlux_T0_faceI;

                if(debug)
                {
                    os<< "Face: " << faceI << "  mag(Sf) = " << magSf_faceI << nl
                        << "magSf_ph1_own = " << magSf_ph1_own[faceI] << "  alpha1f_own = " << magSf_ph1_own[faceI]/magSf_faceI 
                        << "  magSf_ph1_nei = " << magSf_ph1_nei[faceI] << "  alpha1f_nei = " << magSf_ph1_nei[faceI]/magSf_faceI << nl
                        << "magSf_ph0_own = " << magSf_ph0_own[faceI] << "  alpha0f_own = " << magSf_ph0_own[faceI]/magSf_faceI 
                        << "  magSf_ph0_nei = " << magSf_ph0_nei[faceI] << "  alpha0f_nei = " << magSf_ph0_nei[faceI]/magSf_faceI << nl
                        << "face phase state for diffusion flux calculation: " << curPhaseState << nl
                        << "Own: " << faceOwn << nl
                        << "T1Own = " << T1Own << "  T1Nei = " << T1Nei << nl
                        << "T0Own = " << T0Own << "  T0Nei = " << T0Nei << nl
                        << "diffFlux_T1 = " << diffFlux_T1_faceI << "  limiter_ph1 = " << diffFlux_limiter_ph1_faceI << nl
                        << "diffFlux_T0 = " << diffFlux_T0_faceI << "  limiter_ph0 = " << diffFlux_limiter_ph0_faceI << nl
                        << endl;
                }                

                faceI++;
            }//end forAll(pT1, fcI)
        }//end else if(isA<fixedValueFvPatchScalarField>(pT1))
        else
        {
            forAll(pT1, fcI)
            {
                pdiffFlux_T1[fcI] = 0;
                pdiffFlux_T0[fcI] = 0;
                faceI++;
            }
        }
    }//end forAll(T1.boundaryField(), patchI)
    //end boundary faces    
}


void calc_T_from_h
(
    double *x,
    double P,
    double h,
    int n,
    double *Pc,
    double *Tc,
    double *w,
    double *MW,
    double *Tb,
    double *SG,
    double *H8,
    double *kij_T,
    double Ta_kij,
    double Tb_kij,
    int nT_kij,
    double *kij,
    double& T,
    const int& MAX_ITER_T,
    const double& T_TOL,
    const double& H_TOL
)
{
    int i, iter;
    double hm, TOld, TNew, F, dF, T_err, F_err, MW_tmp;
    double v_tmp, Cp_tmp, h_tmp;
    
    MW_tmp = 0;
    for(i=0; i<n; i++)
    {
        MW_tmp += x[i]*MW[i];
    }
    MW_tmp *= 1e-3;

    hm = h*MW_tmp;

    TOld = T;
    TNew = T;

    for(iter=0; iter<MAX_ITER_T; iter++)
    {        
        calc_kij_from_table(TOld, n, Ta_kij, Tb_kij, nT_kij, kij_T, kij);
        calc_v_cp_h_(&P, &TOld, x, &n, Pc, Tc, w, MW, Tb, SG, H8, kij, &v_tmp, &Cp_tmp, &h_tmp);    

        F = h_tmp - hm;
        dF = Cp_tmp;

        TNew = TOld - F/dF;
        
        T_err = mag(TNew - TOld);
        F_err = mag(F);

        TOld = TNew;
        
        if(T_err < T_TOL || F_err < H_TOL) break;
    }

    T = TNew;
}


void calc_h_from_T
(
    double *x,
    double P,
    double T,
    int n,
    double *Pc,
    double *Tc,
    double *w,
    double *MW,
    double *Tb,
    double *SG,
    double *H8,
    double *kij_T,
    double Ta_kij,
    double Tb_kij,
    int nT_kij,
    double *kij,
    double& h
)
{
    int i;
    double MW_tmp;
    double v_tmp, h_tmp;
    
    MW_tmp = 0;
    for(i=0; i<n; i++)
    {
        MW_tmp += x[i]*MW[i];
    }
    MW_tmp *= 1e-3;
    
    calc_kij_from_table(T, n, Ta_kij, Tb_kij, nT_kij, kij_T, kij);
    calc_v_h_(&P, &T, x, &n, Pc, Tc, w, MW, Tb, SG, H8, kij, &v_tmp, &h_tmp);

    h = h_tmp/MW_tmp;
}


void calc_diffFlux_limiter_T
(
    const scalar& alphaOwn,
    const scalar& rhoOwn,
    const List<scalar>& xOwn,
    const scalar& TOwn,
    const scalar& HOwn,
    const scalar& VOwn,
    const scalar& alphaNei,
    const scalar& rhoNei,
    const List<scalar>& xNei,
    const scalar& TNei,
    const scalar& HNei,
    const scalar& VNei,
    double P,
    int n,
    double *Pc,
    double *Tc,
    double *w,
    double *MW,
    double *Tb,
    double *SG,
    double *H8,
    double *kij_T,
    double Ta_kij,
    double Tb_kij,
    int nT_kij,
    double *kij,
    const scalar& dt,
    const scalar& diffFlux,
    scalar& diffFlux_limiter,
    const int& MAX_ITER_T,
    const double& T_TOL,
    const double& H_TOL,
    double TMIN,
    double TMAX
)
{
    int i, iter;
    double TOld, TNew, F, dF, dFOwn, dFNei, T_err, F_err, MWOwn_tmp, MWNei_tmp;
    double vOwn_tmp, CpOwn_tmp, hOwn_tmp;
    double vNei_tmp, CpNei_tmp, hNei_tmp;
    double diffFlux_max, diffFlux_max_1, diffFlux_max_2, diffFlux_max_3;
    double *xOwn_tmp, *xNei_tmp;

    _NNEW_(xOwn_tmp, double, n);
    _NNEW_(xNei_tmp, double, n);

    for(i=0; i<n; i++)
    {
        xOwn_tmp[i] = xOwn[i];
        xNei_tmp[i] = xNei[i];
    }    

    if(mag(diffFlux) > SMALL)
    {
        MWOwn_tmp = 0; MWNei_tmp = 0;
        for(i=0; i<n; i++)
        {
            MWOwn_tmp += xOwn_tmp[i]*MW[i];
            MWNei_tmp += xNei_tmp[i]*MW[i];
        }
        MWOwn_tmp *= 1e-3;
        MWNei_tmp *= 1e-3;

        TOld = 0.5*(TOwn + TNei);

        for(iter=0; iter<MAX_ITER_T; iter++)
        {
            calc_kij_from_table(TOld, n, Ta_kij, Tb_kij, nT_kij, kij_T, kij);

            calc_v_cp_h_(&P, &TOld, xOwn_tmp, &n, Pc, Tc, w, MW, Tb, SG, H8, kij, &vOwn_tmp, &CpOwn_tmp, &hOwn_tmp);            

            calc_v_cp_h_(&P, &TOld, xNei_tmp, &n, Pc, Tc, w, MW, Tb, SG, H8, kij, &vNei_tmp, &CpNei_tmp, &hNei_tmp);

            F = VOwn*alphaOwn*rhoOwn*hOwn_tmp/MWOwn_tmp + VNei*alphaNei*rhoNei*hNei_tmp/MWNei_tmp - HOwn - HNei;
            dFOwn = VOwn*alphaOwn*rhoOwn*CpOwn_tmp/MWOwn_tmp;
            dFNei = VNei*alphaNei*rhoNei*CpNei_tmp/MWNei_tmp;
            dF = dFOwn + dFNei;

            TNew = TOld - F/dF;
        
            T_err = mag(TNew - TOld);
            F_err = mag(F);

            TOld = TNew;
        
            if(T_err < T_TOL || F_err < H_TOL) break;
        }

        if(diffFlux > 0)
        {
            if(TOwn > TNei)
            {
                diffFlux_max_1 = mag(0.25*VOwn*(HOwn - alphaOwn*rhoOwn*hOwn_tmp/MWOwn_tmp));
                calc_kij_from_table(TMIN, n, Ta_kij, Tb_kij, nT_kij, kij_T, kij);
                calc_v_h_(&P, &TMIN, xOwn_tmp, &n, Pc, Tc, w, MW, Tb, SG, H8, kij, &vOwn_tmp, &hOwn_tmp);
                calc_kij_from_table(TMAX, n, Ta_kij, Tb_kij, nT_kij, kij_T, kij);
                calc_v_h_(&P, &TMAX, xNei_tmp, &n, Pc, Tc, w, MW, Tb, SG, H8, kij, &vNei_tmp, &hNei_tmp);

                diffFlux_max_2 = mag(0.25*VOwn*(HOwn - alphaOwn*rhoOwn*hOwn_tmp/MWOwn_tmp));
                diffFlux_max_3 = mag(0.25*VNei*(alphaNei*rhoNei*hNei_tmp/MWNei_tmp - HNei));
                diffFlux_max = min(diffFlux_max_1, min(diffFlux_max_2, diffFlux_max_3));
            }
            else
            {
                diffFlux_max = 0;
            }
        }
        else
        {
            if(TOwn < TNei)
            {
                diffFlux_max_1 = mag(0.25*VNei*(HNei - alphaNei*rhoNei*hNei_tmp/MWNei_tmp));
                calc_kij_from_table(TMAX, n, Ta_kij, Tb_kij, nT_kij, kij_T, kij);
                calc_v_h_(&P, &TMAX, xOwn_tmp, &n, Pc, Tc, w, MW, Tb, SG, H8, kij, &vOwn_tmp, &hOwn_tmp);
                calc_kij_from_table(TMIN, n, Ta_kij, Tb_kij, nT_kij, kij_T, kij);
                calc_v_h_(&P, &TMIN, xNei_tmp, &n, Pc, Tc, w, MW, Tb, SG, H8, kij, &vNei_tmp, &hNei_tmp);
                diffFlux_max_2 = mag(0.25*VOwn*(alphaOwn*rhoOwn*hOwn_tmp/MWOwn_tmp - HOwn));
                diffFlux_max_3 = mag(0.25*VNei*(HNei - alphaNei*rhoNei*hNei_tmp/MWNei_tmp));
                diffFlux_max = min(diffFlux_max_1, min(diffFlux_max_2, diffFlux_max_3));
            }
            else
            {
                diffFlux_max = 0;
            }
        }        

        diffFlux_limiter = min(diffFlux_max/(mag(diffFlux)*dt), 1.0);
    }
    else
    {
        diffFlux_limiter = 0.0;
    }

    _DDELETE_(xOwn_tmp);
    _DDELETE_(xNei_tmp);
}


void calc_diffFlux_limiter2_T
(
    const scalar& alphaOwn,
    const scalar& rhoOwn,
    const List<scalar>& xOwn,
    const scalar& TOwn,
    const scalar& HOwn,
    const scalar& VOwn,
    double P,
    int n,
    double *Pc,
    double* Tc,
    double* w,
    double* MW,
    double *Tb,
    double *SG,
    double *H8,
    double *kij_T,
    double Ta_kij,
    double Tb_kij,
    int nT_kij,
    double *kij,
    const scalar& dt,
    const scalar& diffFlux,
    scalar& diffFlux_limiter,
    double TMIN,
    double TMAX
)
{
    int i;
    double vOwn_tmp, hOwn_tmp, MWOwn_tmp;
    double diffFlux_max;
    double *xOwn_tmp;

    _NNEW_(xOwn_tmp, double, n);

    for(i=0; i<n; i++)
    {
        xOwn_tmp[i] = xOwn[i];
    }

    if(mag(diffFlux) > SMALL)
    {
        MWOwn_tmp = 0;
        for(i=0; i<n; i++)
        {
            MWOwn_tmp += xOwn_tmp[i]*MW[i];
        }
        MWOwn_tmp *= 1e-3;

        if(diffFlux > 0)
        {
            calc_kij_from_table(TMIN, n, Ta_kij, Tb_kij, nT_kij, kij_T, kij);
            calc_v_h_(&P, &TMIN, xOwn_tmp, &n, Pc, Tc, w, MW, Tb, SG, H8, kij, &vOwn_tmp, &hOwn_tmp);                
            diffFlux_max = mag(0.25*VOwn*(HOwn - alphaOwn*rhoOwn*hOwn_tmp/MWOwn_tmp));            
        }
        else
        {
            calc_kij_from_table(TMAX, n, Ta_kij, Tb_kij, nT_kij, kij_T, kij);
            calc_v_h_(&P, &TMAX, xOwn_tmp, &n, Pc, Tc, w, MW, Tb, SG, H8, kij, &vOwn_tmp, &hOwn_tmp);
            diffFlux_max = mag(0.25*VOwn*(alphaOwn*rhoOwn*hOwn_tmp/MWOwn_tmp - HOwn));
        }        

        diffFlux_limiter = min(diffFlux_max/(mag(diffFlux)*dt), 1.0);
    }
    else
    {
        diffFlux_limiter = 0.0;
    }

    _DDELETE_(xOwn_tmp);
}


void calc_diffFlux_limiter
(
    const scalar& rhoOwn,
    const scalar& alphaOwn,
    const scalar& YOwn,
    const scalar& VOwn,
    const scalar& rhoNei,
    const scalar& alphaNei,
    const scalar& YNei,
    const scalar& VNei,
    const scalar& dt,
    const scalar& diffFlux,
    scalar& diffFlux_limiter,
    const scalar& YMIN,
    const scalar& YMAX
)
{
    scalar maxDiffFluxOwn;
    scalar maxDiffFluxNei;    
    scalar maxDiffFlux;
    scalar magDiffFlux = mag(diffFlux);
    diffFlux_limiter = 1;    

    if(diffFlux > 0)
    {
        maxDiffFluxOwn = 0.25*rhoOwn*alphaOwn*(YOwn - YMIN)*VOwn/dt;
        
        maxDiffFluxNei = 0.25*rhoNei*alphaNei*(YMAX - YNei)*VNei/dt;
        
        maxDiffFlux = min(maxDiffFluxOwn, maxDiffFluxNei);
        
        if(magDiffFlux > maxDiffFlux)
        {
            if(magDiffFlux < SMALL)
            {
                magDiffFlux += SMALL;
            }
            diffFlux_limiter = maxDiffFlux/magDiffFlux;
        }        
    }
    else
    {
        maxDiffFluxNei = 0.25*rhoNei*alphaNei*(YNei - YMIN)*VNei/dt;

        maxDiffFluxOwn = 0.25*rhoOwn*alphaOwn*(YMAX - YOwn)*VOwn/dt;

        maxDiffFlux = min(maxDiffFluxOwn, maxDiffFluxNei);

        if(magDiffFlux > maxDiffFlux)
        {
            if(magDiffFlux < SMALL)
            {
                magDiffFlux += SMALL;
            }
            diffFlux_limiter = maxDiffFlux/magDiffFlux;
        }
    }       
}


void calc_diffFlux_limiter2
(
    const scalar& rhoOwn,
    const scalar& alphaOwn,
    const scalar& YOwn,
    const scalar& VOwn,
    const scalar& rhoNei,
    const scalar& alphaNei,
    const scalar& YNei,
    const scalar& VNei,
    const scalar& dt,
    const scalar& diffFlux,
    scalar& diffFlux_limiter,
    const scalar& YMIN,
    const scalar& YMAX
)
{
    scalar maxDiffFluxOwn;
    scalar maxDiffFluxNei;
    scalar coeff;
    scalar den;
    scalar maxDiffFluxOwnNei;
    scalar maxDiffFlux;
    scalar magDiffFlux = mag(diffFlux);
    diffFlux_limiter = 1;    

    coeff = rhoOwn*rhoNei*alphaOwn*alphaNei*VOwn*VNei/dt;
    den = rhoOwn*alphaOwn*VOwn + rhoNei*alphaNei*VNei;
    if(den < SMALL)
    {
        den += SMALL;
    }
    coeff /= den;

    if(diffFlux > 0)
    {
        maxDiffFluxOwn = 0.25*rhoOwn*alphaOwn*(YOwn - YMIN)*VOwn/dt;
        
        maxDiffFluxNei = 0.25*rhoNei*alphaNei*(YMAX - YNei)*VNei/dt;

        if(YOwn > YNei)
        {
            maxDiffFluxOwnNei = coeff*(YOwn - YNei);        
        }
        else
        {
            maxDiffFluxOwnNei = 0;
        }

        maxDiffFlux = min(min(maxDiffFluxOwn, maxDiffFluxNei), maxDiffFluxOwnNei);
        
        if(magDiffFlux > maxDiffFlux)
        {
            if(magDiffFlux < SMALL)
            {
                magDiffFlux += SMALL;
            }
            diffFlux_limiter = maxDiffFlux/magDiffFlux;
        }        
    }
    else
    {
        maxDiffFluxNei = 0.25*rhoNei*alphaNei*(YNei - YMIN)*VNei/dt;

        maxDiffFluxOwn = 0.25*rhoOwn*alphaOwn*(YMAX - YOwn)*VOwn/dt;

        if(YNei > YOwn)
        {
            maxDiffFluxOwnNei = coeff*(YNei - YOwn);        
        }
        else
        {
            maxDiffFluxOwnNei = 0;
        }

        maxDiffFlux = min(min(maxDiffFluxOwn, maxDiffFluxNei), maxDiffFluxOwnNei);

        if(magDiffFlux > maxDiffFlux)
        {
            if(magDiffFlux < SMALL)
            {
                magDiffFlux += SMALL;
            }
            diffFlux_limiter = maxDiffFlux/magDiffFlux;
        }
    }       
}


template<class Type>
void linearInterpolate_2ph
(
    const GeometricField<Type, fvPatchField, volMesh>& Y,
    const fvMesh& mesh,
    const surfaceScalarField& weights,
    GeometricField<Type, fvsPatchField, surfaceMesh>& Yf
)
{
    label faceI, faceOwn, faceNei;
    scalar YOwn;

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();

    const Field<Type>& YCells = Y.internalField();

    for(faceI=0; faceI<mesh.nInternalFaces(); faceI++)
    {
        faceOwn = own[faceI];
        faceNei = nei[faceI];
        scalar wf = weights[faceI];

        Yf[faceI] = wf*YCells[faceOwn] + (1.0 - wf)*YCells[faceNei];
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(Yf.boundaryField(), patchI)
    {
        const polyPatch& pp = patches[patchI];
        fvsPatchField<Type>& pYf = Yf.boundaryField()[patchI];
        const fvPatchField<Type>& pY = Y.boundaryField()[patchI];
        const fvsPatchScalarField& pw = weights.boundaryField()[patchI];
        faceI = pp.start();

        if(pp.coupled())
        {
            //const Field<Type>& pYOwn = pY.patchInternalField();
            const Field<Type>& pYNei = pY.patchNeighbourField();
            forAll(pYf, fcI)
            {
                faceOwn = own[faceI];
                YOwn = YCells[faceOwn];
                scalar wf = pw[fcI];
                pYf[fcI] = wf*YOwn + (1.0 - wf)*pYNei[fcI];

                faceI++;
            }
        }
        else
        {
            forAll(pYf, fcI)
            {
                scalar wf = pw[fcI];
                pYf[fcI] = wf*pY[fcI] + (1.0 - wf)*pY[fcI];
            }
        }
    }
}


template<class Type>
void linearInterpolate_2ph
(
    const GeometricField<Type, fvPatchField, volMesh>& Y,
    const fvMesh& mesh,
    const surfaceScalarField& weights,
    GeometricField<Type, fvsPatchField, surfaceMesh>& Yf,
    const bool debug,
    OFstream& os
)
{
    label faceI, faceOwn, faceNei;
    scalar wf, YOwn;

    if(debug)
    {
        print_line(os, 100);
        os<< "Interpolating field " << Y.name() << endl; 
        print_line(os, 100);
        os<< endl;
    }

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();

    const Field<Type>& YCells = Y.internalField();

    for(faceI=0; faceI<mesh.nInternalFaces(); faceI++)
    {
        faceOwn = own[faceI];
        faceNei = nei[faceI];
        wf = weights[faceI];

        Yf[faceI] = wf*YCells[faceOwn] + (1.0 - wf)*YCells[faceNei];
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const wordList& patchNames = patches.names();

    forAll(Yf.boundaryField(), patchI)
    {
        const polyPatch& pp = patches[patchI];
        fvsPatchField<Type>& pYf = Yf.boundaryField()[patchI];
        const fvPatchField<Type>& pY = Y.boundaryField()[patchI];
        const fvsPatchScalarField& pw = weights.boundaryField()[patchI];
        faceI = pp.start();

        if(pp.coupled())
        {
            if(debug)
            {
                print_line(os, 100);
                os<< "Couples patch " << patchNames[patchI] << endl; 
                print_line(os, 100);
                os<< setw(6) << "faceI" << "  " << setw(8) << "own" << "  " << setw(8) << "own w" << setw(8) << "own val" << "  " << setw(8) << "pif val" << "  " << setw(8) << "pnf val" << "  " << setw(8) << "face val" << endl;
                print_line(os, 100);
            }

            const Field<Type>& pYOwn = pY.patchInternalField();
            const Field<Type>& pYNei = pY.patchNeighbourField();
            forAll(pYf, fcI)
            {
                faceOwn = own[faceI];
                YOwn = YCells[faceOwn];
                wf = pw[fcI];
                pYf[fcI] = wf*YOwn + (1.0 - wf)*pYNei[fcI];
                if(debug)
                {
                    os<< setw(6) << fcI << "  " << setw(8) << faceOwn << "  " << setw(8) << wf << "  " << setw(8) << YCells[faceOwn] << "  " << setw(8) << pYOwn[fcI] << "  " << setw(8) << pYNei[fcI] << "  " << setw(8) << pYf[fcI] << endl;
                }

                faceI++;
            }
            
            if(debug)
            {
                print_line(os, 100);
                os<< endl;
            }
        }
        else
        {
            forAll(pYf, fcI)
            {
                scalar wf = pw[fcI];
                pYf[fcI] = wf*pY[fcI] + (1.0 - wf)*pY[fcI];
            }
        }
    }
}


void calc_2ph_linearInterpolation_weights
(
    const fvMesh& mesh,
    const vectorField& C_ph,
    const vectorField& Cf_ph,
    surfaceScalarField& weights
)
{
    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();

    for(label faceI=0; faceI<mesh.nInternalFaces(); faceI++)
    {        
        scalar SfdOwn = mag(Sf[faceI] & (Cf_ph[faceI] - C_ph[own[faceI]]));
        scalar SfdNei = mag(Sf[faceI] & (C_ph[nei[faceI]] - Cf_ph[faceI]));
        weights[faceI] = SfdNei/(SfdOwn + SfdNei);        
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    //--------------------------------------------------------------//
    //Need phase centroid in cell neighbour on neighbouring 
    //processor for each coupled patch face in order to 
    //calculate weight for that face
    const label nBnd = mesh.nFaces() - mesh.nInternalFaces();
    List<vector> CNei(nBnd);

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        if(pp.coupled())
        {
            label faceI = pp.start();
            label bndFaceI = pp.start() - mesh.nInternalFaces();
            forAll(pp, fcI)
            {
                CNei[bndFaceI] = C_ph[own[faceI]];
                faceI++;
                bndFaceI++;
            }
        }
    }

    syncTools::swapBoundaryFaceList(mesh, CNei);
    //CNei now has phase centroid of neighbouring cell for each 
    //coupled patch face
    //--------------------------------------------------------------//

    forAll(weights.boundaryField(), patchI)
    {
        const polyPatch& pp = patches[patchI];
        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchI];
        fvsPatchScalarField& pw = weights.boundaryField()[patchI];
        if(pp.coupled())
        {
            label faceI = pp.start();
            label bndFaceI = pp.start() - mesh.nInternalFaces();
            forAll(pw, fcI)
            {                
                scalar SfdOwn = mag(pSf[fcI] & (Cf_ph[faceI] - C_ph[own[faceI]]));
                scalar SfdNei = mag(pSf[fcI] & (CNei[bndFaceI] - Cf_ph[faceI]));                
                pw[fcI] = SfdNei/(SfdOwn + SfdNei);
                faceI++;
                bndFaceI++;
            }
        }
        else
        {            
            forAll(pw, fcI)
            {                
                pw[fcI] = 1.0;                
            }
        }
    }
}


void calc_face_phaseState
(
    const scalar& curAlpha1f_own,
    const scalar& curAlpha0f_own,
    const scalar& curAlpha1f_nei,
    const scalar& curAlpha0f_nei,
    const scalar& MIN_ALPHA_DIFF,
    label& curPhaseState
)
{
    scalar MAX_ALPHA_DIFF = 1 - MIN_ALPHA_DIFF;

    if((curAlpha1f_own < MIN_ALPHA_DIFF && curAlpha0f_own < MIN_ALPHA_DIFF) || (curAlpha1f_nei < MIN_ALPHA_DIFF && curAlpha0f_nei < MIN_ALPHA_DIFF))
    {
        curPhaseState = 3;
    }
    else
    {
        if(curAlpha1f_own < MIN_ALPHA_DIFF || curAlpha1f_nei < MIN_ALPHA_DIFF)
        {
            curPhaseState = 0;
        }
        else
        {
            if(curAlpha1f_own > MAX_ALPHA_DIFF || curAlpha1f_nei > MAX_ALPHA_DIFF)
            {
                curPhaseState = 1;
            }
            else
            {
                curPhaseState = 2;
            }
        }
    }
}


void calc_face_phaseState_diff
(
    const fvMesh& mesh,
    const volScalarField& Y1i,
    const scalarField& magSf_ph1_own,
    const scalarField& magSf_ph0_own,
    const scalarField& magSf_ph1_nei,
    const scalarField& magSf_ph0_nei,
    const scalar& MIN_ALPHA_DIFF,
    labelList& face_phaseState,
    const bool debug,
    OFstream& os
)
{
    scalar curMagSf; 
    scalar curMagSf_ph1_own; 
    scalar curMagSf_ph0_own;
    scalar curAlpha1f_own; 
    scalar curAlpha0f_own; 
    scalar curMagSf_ph1_nei; 
    scalar curMagSf_ph0_nei;
    scalar curAlpha1f_nei; 
    scalar curAlpha0f_nei;
    label curPhaseState;    

    for(label faceI=0; faceI<mesh.nInternalFaces(); faceI++)
    {        
        curMagSf = mesh.magSf()[faceI];       
        curMagSf_ph1_own = magSf_ph1_own[faceI];
        curMagSf_ph0_own = magSf_ph0_own[faceI];
        curMagSf_ph1_nei = magSf_ph1_nei[faceI];
        curMagSf_ph0_nei = magSf_ph0_nei[faceI];

        curAlpha1f_own = curMagSf_ph1_own/curMagSf;
        curAlpha0f_own = curMagSf_ph0_own/curMagSf;                
        curAlpha1f_nei = curMagSf_ph1_nei/curMagSf;
        curAlpha0f_nei = curMagSf_ph0_nei/curMagSf;

        curPhaseState = 2;
        calc_face_phaseState(curAlpha1f_own, curAlpha0f_own, curAlpha1f_nei, curAlpha0f_nei, MIN_ALPHA_DIFF, curPhaseState);

        face_phaseState[faceI] = curPhaseState;

        if(debug)
        {
            os<< "Face: " << faceI << "  magSf = " << mesh.magSf()[faceI] << nl
                << "Af_ph1_own = " << magSf_ph1_own[faceI] << "  Af_ph1_nei = " << magSf_ph1_nei[faceI] << nl
                << "alpha1f_own = " << curAlpha1f_own << "  alpha1f_nei = " << curAlpha1f_nei << nl
                << "Af_ph0_own = " << magSf_ph0_own[faceI] << "  Af_ph0_nei = " << magSf_ph0_nei[faceI] << nl
                << "alpha0f_own = " << curAlpha0f_own << "  alpha0f_nei = " << curAlpha0f_nei << nl
                << "face phase state: " << curPhaseState << nl
                << endl;
        }
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(Y1i.boundaryField(), patchI)
    {
        const polyPatch& pp = patches[patchI];        
        const fvPatchScalarField& pY1i = Y1i.boundaryField()[patchI];        
        const fvsPatchScalarField& pMagSf = mesh.magSf().boundaryField()[patchI];
        label faceI = pp.start();        

        if(pp.coupled())
        {
            forAll(pY1i, fcI)
            {
                if(debug)
                {
                    os<< "Face: " << faceI << "  patch face index: " << fcI << "  magSf = " << pMagSf[fcI] << nl
                        << "Af_ph1_own = " << magSf_ph1_own[faceI] << "  Af_ph1_nei = " << magSf_ph1_nei[faceI] << nl                       
                        << "Af_ph0_own = " << magSf_ph0_own[faceI] << "  Af_ph0_nei = " << magSf_ph0_nei[faceI]                       
                        << endl;
                }

                curMagSf = pMagSf[fcI];       
                curMagSf_ph1_own = magSf_ph1_own[faceI];
                curMagSf_ph0_own = magSf_ph0_own[faceI];
                curMagSf_ph1_nei = magSf_ph1_nei[faceI];
                curMagSf_ph0_nei = magSf_ph0_nei[faceI];

                curAlpha1f_own = curMagSf_ph1_own/curMagSf;
                curAlpha0f_own = curMagSf_ph0_own/curMagSf;                
                curAlpha1f_nei = curMagSf_ph1_nei/curMagSf;
                curAlpha0f_nei = curMagSf_ph0_nei/curMagSf;

                curPhaseState = 2;
                calc_face_phaseState(curAlpha1f_own, curAlpha0f_own, curAlpha1f_nei, curAlpha0f_nei, MIN_ALPHA_DIFF, curPhaseState);

                face_phaseState[faceI] = curPhaseState;

                if(debug)
                {
                    os<< "face phase state: " << curPhaseState << nl
                        << endl;
                }

                faceI++;
            }
        }
        else
        {
            forAll(pY1i, fcI)
            {
                if(debug)
                {
                    os<< "Face: " << faceI << "  patch face index: " << fcI << "  magSf = " << pMagSf[fcI] << nl
                        << "Af_ph1_own = " << magSf_ph1_own[faceI] << "  Af_ph1_nei = " << magSf_ph1_nei[faceI] << nl                       
                        << "Af_ph0_own = " << magSf_ph0_own[faceI] << "  Af_ph0_nei = " << magSf_ph0_nei[faceI] 
                        << endl;
                }

                curMagSf = pMagSf[fcI];       
                curMagSf_ph1_own = magSf_ph1_own[faceI];
                curMagSf_ph0_own = magSf_ph0_own[faceI];                

                curAlpha1f_own = curMagSf_ph1_own/curMagSf;
                curAlpha0f_own = curMagSf_ph0_own/curMagSf;
                curAlpha1f_nei = curAlpha1f_own;
                curAlpha0f_nei = curAlpha0f_own;

                curPhaseState = 2;
                calc_face_phaseState(curAlpha1f_own, curAlpha0f_own, curAlpha1f_nei, curAlpha0f_nei, MIN_ALPHA_DIFF, curPhaseState);

                face_phaseState[faceI] = curPhaseState;             

                if(debug)
                {
                    os<< "face phase state: " << curPhaseState << nl
                        << endl;
                }

                faceI++;
            }
        }
    }
}


void calc_2ph_Cf
(
    const vectorField& Cf_ph1_own,
    const vectorField& Cf_ph0_own,
    const vectorField& Cf_ph1_nei,
    const vectorField& Cf_ph0_nei,
    const scalarField& Af_ph1_own,
    const scalarField& Af_ph0_own,
    const scalarField& Af_ph1_nei,
    const scalarField& Af_ph0_nei,
    vectorField& Cf_ph1,
    vectorField& Cf_ph0
)
{
    forAll(Cf_ph1_own, faceI)
    {
        if(Af_ph1_nei[faceI] < Af_ph1_own[faceI])
        {
            Cf_ph1[faceI] = Cf_ph1_nei[faceI];
        }
        else
        {
            Cf_ph1[faceI] = Cf_ph1_own[faceI];
        }

        if(Af_ph0_nei[faceI] < Af_ph0_own[faceI])
        {
            Cf_ph0[faceI] = Cf_ph0_nei[faceI];
        }
        else
        {
            Cf_ph0[faceI] = Cf_ph0_own[faceI];
        }
    }
}


void redistribute_alpha_field
(
    const fvMesh& mesh,
    scalarField& alphaCells,
    const labelListList& cellCellsAll,
    const scalar& ALPHA_MIN_BOUND,
    const label& ALPHA_BOUND_ITERS_MAX,
    label& minAlpha1Cell,
    label& maxAlpha1Cell,
    scalar& minAlpha1,
    scalar& maxAlpha1,
    const bool debug,
    OFstream& os
)
{
    label minAlphaCell, maxAlphaCell, nIters, curCell;
    int i, nCells;
    bool allNeiDone;
    scalar ALPHA_MAX_BOUND, minAlpha, maxAlpha, tAlpha;

    nCells = mesh.nCells();
    ALPHA_MAX_BOUND = 1 - ALPHA_MIN_BOUND;

    forAll(alphaCells, cellI)
    {
        if(alphaCells[cellI] < ALPHA_MIN_BOUND)
        {
            if(debug)
            {
                os<< "Correcting alpha1 in cell " << cellI
                    << "  alpha1 = " << alphaCells[cellI] << endl;
            }

            const labelList& curCellCells = cellCellsAll[cellI];
            minAlpha = 1;
            minAlphaCell = cellI;
            nIters = 0;
            allNeiDone = true;
            for(i=0; i<curCellCells.size(); i++)
            {
                curCell = curCellCells[i];
                if(!(curCell==cellI) && (curCell < nCells) && (alphaCells[curCell] > 0))
                {
                    allNeiDone = false;
                    break;
                }
            }

            if(!allNeiDone)
            {
                do
                {
                    if(debug)
                    {
                        os<< "alpha1 redistribution iteration no: " << nIters+1 << endl;
                    }
            
                    allNeiDone = true;
                    minAlpha = 1;
                    for(i=0; i<curCellCells.size(); i++)
                    {
                        curCell = curCellCells[i];
                        if(!(curCell==cellI) && (curCell < nCells) && (alphaCells[curCell] > 0))
                        {
                            allNeiDone = false;
                            if(alphaCells[curCell] < minAlpha)
                            {
                                minAlphaCell = curCell;
                                minAlpha = alphaCells[curCell];
                            }
                        }
                    }            

                    if(debug)
                    {
                        os<< "Nei cell with minimum non-zero alpha1: " << minAlphaCell
                            << "  Nei cell alpha1 = " << alphaCells[minAlphaCell] << endl;
                    }

                    tAlpha = alphaCells[minAlphaCell] + alphaCells[cellI];
                    alphaCells[minAlphaCell] = max(tAlpha, 0);
                    alphaCells[cellI] = min(tAlpha, 0);

                    if(debug)
                    {
                        os<< "New nei cell alpha1 = " << alphaCells[minAlphaCell]
                            << "  New cell alpha1 = " << alphaCells[cellI] << endl;
                    }
            
                    nIters++;
                }while(alphaCells[cellI] < ALPHA_MIN_BOUND && !allNeiDone && nIters < ALPHA_BOUND_ITERS_MAX);        
            }        
        }

        if(alphaCells[cellI] > ALPHA_MAX_BOUND)
        {
            if(debug)
            {
                os<< "Correcting alpha1 in cell " << cellI
                    << "  alpha1 = " << alphaCells[cellI] << endl;
            }

            const labelList& curCellCells = cellCellsAll[cellI];
            maxAlpha = 0;
            maxAlphaCell = cellI;
            nIters = 0;
            allNeiDone = true;
            for(i=0; i<curCellCells.size(); i++)
            {
                curCell = curCellCells[i];
                if(!(curCell==cellI) && (curCell < nCells) && (alphaCells[curCell] < 1))
                {
                    allNeiDone = false;
                    break;
                }
            }
        
            if(!allNeiDone)
            {
                do
                {
                    if(debug)
                    {
                        os<< "alpha1 correction iteration no: " << nIters+1 << endl;
                    }

                    maxAlpha = 0;
                    allNeiDone = true;
                    for(i=0; i<curCellCells.size(); i++)
                    {
                        curCell = curCellCells[i];
                        if(!(curCell==cellI) && (curCell < nCells) && (alphaCells[curCell] < 1))
                        {
                            allNeiDone = false;
                            if(alphaCells[curCell] > maxAlpha)
                            {
                                maxAlphaCell = curCell;
                                maxAlpha = alphaCells[curCell];
                            }
                        }
                    }

                    if(debug)
                    {
                        os<< "Nei cell with maximum alpha1 below 1: " << maxAlphaCell
                            << "  Nei cell alpha1 = " << alphaCells[maxAlphaCell] << endl;
                    }
            
                    tAlpha = alphaCells[maxAlphaCell] + alphaCells[cellI] - 1;
                    alphaCells[maxAlphaCell] = min(tAlpha, 1);
                    alphaCells[cellI] = max(tAlpha, 1);

                    if(debug)
                    {
                        os<< "New nei cell alpha1 = " << alphaCells[maxAlphaCell]
                            << "  New cell alpha1 = " << alphaCells[cellI] << nl << endl;
                    }

                    nIters++;
                }while(alphaCells[cellI] > ALPHA_MAX_BOUND && !allNeiDone && nIters < ALPHA_BOUND_ITERS_MAX);        
            }
        }

        if(alphaCells[cellI] < minAlpha1)
        {
            minAlpha1 = alphaCells[cellI];
            minAlpha1Cell = cellI;
        }

        if(alphaCells[cellI] > maxAlpha1)
        {
            maxAlpha1 = alphaCells[cellI];
            maxAlpha1Cell = cellI;
        }
    }
}


void redistribute_alpha_calc_dC_field
(
    const fvMesh& mesh,
    volScalarField& alpha1,
    const volScalarField& rho1,
    const PtrList<volScalarField>& Y1,
    const volScalarField& rho0,
    const PtrList<volScalarField>& Y0,
    List<scalarField>& dC1_redistCells,
    List<scalarField>& dC0_redistCells,
    const labelListList& cellCellsAll,
    int n,
    const scalar& ALPHA_MIN_BOUND,
    const label& ALPHA_BOUND_ITERS_MAX,    
    label& minAlpha1Cell,
    label& maxAlpha1Cell,
    scalar& minAlpha1,
    scalar& maxAlpha1,
    const bool debug,
    OFstream& os
)
{
    label minAlphaCell, maxAlphaCell, nIters, curCell, cellJ;
    int i, nCells;
    bool allNeiDone;
    scalar ALPHA_MAX_BOUND, minAlpha, maxAlpha, tAlpha;
    scalar alpha1_cellI, alpha1_curCell, alpha1_cellJ, dAlpha1_cellI, dAlpha1_cellJ, rho1_cellI, rho1_cellJ, dC1_cellI, dC1_cellJ, dAlpha0_cellI, dAlpha0_cellJ, rho0_cellI, rho0_cellJ, dC0_cellI, dC0_cellJ;

    List<scalarField> Y1Cells(n);
    List<scalarField> Y0Cells(n);
    for(i=0; i<n; i++)
    {
        const scalarField& Y1iCells = Y1[i].internalField();
        Y1Cells[i] = Y1iCells;
        const scalarField& Y0iCells = Y0[i].internalField();
        Y0Cells[i] = Y0iCells;
    }

    scalarField& alpha1Cells = alpha1.internalField();
    const scalarField& rho1Cells = rho1.internalField();
    const scalarField& rho0Cells = rho0.internalField();

    nCells = mesh.nCells();
    ALPHA_MAX_BOUND = 1 - ALPHA_MIN_BOUND;

    forAll(alpha1Cells, cellI)
    {
        for(i=0; i<(n-1); i++)
        {
            dC1_redistCells[i][cellI] = 0.0;
            dC0_redistCells[i][cellI] = 0.0;
        }
    }

    forAll(alpha1Cells, cellI)
    {
        alpha1_cellI = alpha1Cells[cellI];

        if(alpha1_cellI < ALPHA_MIN_BOUND)
        {
            if(debug)
            {
                os<< "Correcting alpha1 in cell " << cellI
                    << "  alpha1 = " << alpha1_cellI << endl;
            }

            const labelList& curCellCells = cellCellsAll[cellI];
            minAlpha = 1;
            minAlphaCell = cellI;
            nIters = 0;
            allNeiDone = true;
            for(cellJ=0; cellJ<curCellCells.size(); cellJ++)
            {
                curCell = curCellCells[cellJ];
                if(!(curCell==cellI) && (curCell < nCells) && (alpha1Cells[curCell] > 0))
                {
                    allNeiDone = false;
                    break;
                }
            }

            if(!allNeiDone)
            {
                do
                {
                    if(debug)
                    {
                        os<< "alpha1 redistribution iteration no: " << nIters+1 << endl;
                    }
            
                    allNeiDone = true;
                    minAlpha = 1;
                    for(cellJ=0; cellJ<curCellCells.size(); cellJ++)
                    {
                        curCell = curCellCells[cellJ];
                        alpha1_curCell = alpha1Cells[curCell];

                        if(!(curCell==cellI) && (curCell < nCells) && (alpha1_curCell > 0))
                        {
                            allNeiDone = false;
                            if(alpha1_curCell < minAlpha)
                            {
                                minAlphaCell = curCell;
                                minAlpha = alpha1_curCell;
                            }
                        }
                    }            

                    if(debug)
                    {
                        os<< "Nei cell with minimum non-zero alpha1: " << minAlphaCell
                            << "  Nei cell alpha1 = " << minAlpha << endl;
                    }
                                        
                    alpha1_cellJ = alpha1Cells[minAlphaCell];
                    tAlpha = alpha1_cellJ + alpha1_cellI;
                    alpha1Cells[minAlphaCell] = max(tAlpha, 0);
                    alpha1Cells[cellI] = min(tAlpha, 0);
                    dAlpha1_cellJ = alpha1Cells[minAlphaCell] - alpha1_cellJ;
                    dAlpha0_cellI = dAlpha1_cellJ;

                    rho1_cellJ = rho1Cells[minAlphaCell];
                    rho0_cellI = rho0Cells[cellI];
                    for(i=0; i<(n-1); i++)
                    {
                        dC1_cellJ = dAlpha1_cellJ*rho1_cellJ*Y1Cells[i][minAlphaCell];                        
                        dC1_redistCells[i][minAlphaCell] += dC1_cellJ;
                        dC1_redistCells[i][cellI] -= dC1_cellJ;

                        dC0_cellI = dAlpha0_cellI*rho0_cellI*Y0Cells[i][cellI];
                        dC0_redistCells[i][cellI] += dC0_cellI;
                        dC0_redistCells[i][minAlphaCell] -= dC0_cellI;
                    }

                    if(debug)
                    {
                        os<< "New nei cell alpha1 = " << alpha1Cells[minAlphaCell]
                            << "  New cell alpha1 = " << alpha1Cells[cellI] << endl;
                    }
            
                    nIters++;
                }while(alpha1Cells[cellI] < ALPHA_MIN_BOUND && !allNeiDone && nIters < ALPHA_BOUND_ITERS_MAX);        
            }        
        }

        if(alpha1_cellI > ALPHA_MAX_BOUND)
        {
            if(debug)
            {
                os<< "Correcting alpha1 in cell " << cellI
                    << "  alpha1 = " << alpha1_cellI << endl;
            }

            const labelList& curCellCells = cellCellsAll[cellI];
            maxAlpha = 0;
            maxAlphaCell = cellI;
            nIters = 0;
            allNeiDone = true;
            for(cellJ=0; cellJ<curCellCells.size(); cellJ++)
            {
                curCell = curCellCells[cellJ];
                if(!(curCell==cellI) && (curCell < nCells) && (alpha1Cells[curCell] < 1))
                {
                    allNeiDone = false;
                    break;
                }
            }
        
            if(!allNeiDone)
            {
                do
                {
                    if(debug)
                    {
                        os<< "alpha1 correction iteration no: " << nIters+1 << endl;
                    }

                    maxAlpha = 0;
                    allNeiDone = true;
                    for(cellJ=0; cellJ<curCellCells.size(); cellJ++)
                    {
                        curCell = curCellCells[cellJ];
                        alpha1_curCell = alpha1Cells[curCell];
                        if(!(curCell==cellI) && (curCell < nCells) && (alpha1_curCell < 1))
                        {
                            allNeiDone = false;
                            if(alpha1_curCell > maxAlpha)
                            {
                                maxAlphaCell = curCell;
                                maxAlpha = alpha1_curCell;
                            }
                        }
                    }

                    if(debug)
                    {
                        os<< "Nei cell with maximum alpha1 below 1: " << maxAlphaCell
                            << "  Nei cell alpha1 = " << maxAlpha << endl;
                    }
                    
                    alpha1_cellJ = alpha1Cells[maxAlphaCell];
                    tAlpha = alpha1_cellJ + alpha1_cellI - 1;
                    alpha1Cells[maxAlphaCell] = min(tAlpha, 1);
                    alpha1Cells[cellI] = max(tAlpha, 1);
                    dAlpha1_cellI = alpha1Cells[cellI] - alpha1_cellI;
                    dAlpha0_cellJ = dAlpha1_cellI;

                    rho1_cellI = rho1Cells[cellI];
                    rho0_cellJ = rho0Cells[cellJ];
                    for(i=0; i<(n-1); i++)
                    {
                        dC1_cellI = dAlpha1_cellI*rho1_cellI*Y1Cells[i][cellI];
                        dC1_redistCells[i][cellI] += dC1_cellI;
                        dC1_redistCells[i][maxAlphaCell] -= dC1_cellI;

                        dC0_cellJ = dAlpha0_cellJ*rho0_cellJ*Y0Cells[i][maxAlphaCell];
                        dC0_redistCells[i][maxAlphaCell] += dC0_cellJ;
                        dC0_redistCells[i][cellI] -= dC0_cellJ;
                    }

                    if(debug)
                    {
                        os<< "New nei cell alpha1 = " << alpha1Cells[maxAlphaCell]
                            << "  New cell alpha1 = " << alpha1Cells[cellI] << nl << endl;
                    }

                    nIters++;
                }while(alpha1Cells[cellI] > ALPHA_MAX_BOUND && !allNeiDone && nIters < ALPHA_BOUND_ITERS_MAX);        
            }
        }

        if(alpha1Cells[cellI] < minAlpha1)
        {
            minAlpha1 = alpha1Cells[cellI];
            minAlpha1Cell = cellI;
        }

        if(alpha1Cells[cellI] > maxAlpha1)
        {
            maxAlpha1 = alpha1Cells[cellI];
            maxAlpha1Cell = cellI;
        }
    }
}


void redistribute_Ci_field
(
    const fvMesh& mesh,
    scalarField& CCells,
    scalarField& YCells,
    const scalarField& alphaCells,
    const scalarField& rhoCells,
    const labelListList& cellStencil,
    const scalar& Y_MIN,
    const scalar& Y_MAX,
    const scalar& ALPHA_MIN,
    const label& Y_BOUND_ITERS_MAX,
    const bool debug,
    OFstream& os
)
{    
    forAll(YCells, cellI)
    {
        if(YCells[cellI] < Y_MIN && alphaCells[cellI] > ALPHA_MIN)
        {            
            if(debug)
            {
                os<< "---------------------------------------------------------------------" << nl
                    << "Cell: " << cellI << "  alpha = " << alphaCells[cellI] << "  Ci = " << CCells[cellI] << "  Yi = " << YCells[cellI] << endl; 
            }

            const labelList& curCellCells = cellStencil[cellI];
            scalar minY;
            label minYCell = cellI;
            label nIters = 0;
            bool allNeiDone = true;
            for(label i=0; i<curCellCells.size(); i++)
            {
                label curCell = curCellCells[i];

                if(!(curCell==cellI) && (curCell < mesh.nCells()) && (YCells[curCell] > Y_MIN) && (alphaCells[curCell] > ALPHA_MIN))
                {
                    allNeiDone = false;
                    break;
                }
            }

            if(!allNeiDone)
            {
                do
                {                     
                    allNeiDone = true;
                    minY = 1.1;
                    for(label i=0; i<curCellCells.size(); i++)
                    {
                        label curCell = curCellCells[i];

                        if(!(curCell==cellI) && (curCell < mesh.nCells()) && (YCells[curCell] > Y_MIN) && (alphaCells[curCell] > ALPHA_MIN))
                        {
                            allNeiDone = false;
                            if(YCells[curCell] < minY)
                            {
                                minYCell = curCell;
                                minY = YCells[curCell];
                            }
                        }
                    }

                    if(debug)
                    {
                        os<< "Redist iteration: " << nIters << nl
                            << "minYCell: " << minYCell << "  minY = " << minY << endl;
                    }

                    scalar tC1 = rhoCells[minYCell]*alphaCells[minYCell]*(YCells[minYCell] - Y_MIN);
                    scalar tC2 = rhoCells[cellI]*alphaCells[cellI]*(Y_MIN - YCells[cellI]);
                    scalar tC = min(tC1, tC2);

                    CCells[minYCell] -= tC;
                    CCells[cellI] += tC;

                    YCells[minYCell] = CCells[minYCell]/rhoCells[minYCell]/alphaCells[minYCell];
                    YCells[cellI] = CCells[cellI]/rhoCells[cellI]/alphaCells[cellI];
            
                    nIters++;
                }while(YCells[cellI] < Y_MIN && !allNeiDone && nIters < Y_BOUND_ITERS_MAX);        
            }

            if(debug)
            {
                os<< "---------------------------------------------------------------------" << endl;
            }
        }

        if(YCells[cellI] > Y_MAX && alphaCells[cellI] > ALPHA_MIN)
        {    
            if(debug)
            {
                os<< "---------------------------------------------------------------------" << nl
                    << "Cell: " << cellI << "  alpha = " << alphaCells[cellI] << "  Ci = " << CCells[cellI] << "  Yi = " << YCells[cellI] << endl; 
            }
        
            const labelList& curCellCells = cellStencil[cellI];
            scalar maxY;
            label maxYCell = cellI;
            label nIters = 0;
            bool allNeiDone = true;
            for(label i=0; i<curCellCells.size(); i++)
            {
                label curCell = curCellCells[i];

                if(!(curCell==cellI) && (curCell < mesh.nCells()) && (YCells[curCell] < Y_MAX) && (alphaCells[curCell] > ALPHA_MIN))
                {
                    allNeiDone = false;
                    break;
                }
            }

            if(!allNeiDone)
            {
                do
                {                                
                    allNeiDone = true;
                    maxY = -0.1;
                    for(label i=0; i<curCellCells.size(); i++)
                    {
                        label curCell = curCellCells[i];

                        if(!(curCell==cellI) && (curCell < mesh.nCells()) && (YCells[curCell] < Y_MAX) && (alphaCells[curCell] > ALPHA_MIN))
                        {
                            allNeiDone = false;
                            if(YCells[curCell] > maxY)
                            {
                                maxYCell = curCell;
                                maxY = YCells[curCell];
                            }
                        }
                    }

                    if(debug)
                    {
                        os<< "Redist iteration: " << nIters << nl
                            << "maxYCell: " << maxYCell << "  maxY = " << maxY << endl;
                    }

                    scalar tC1 = rhoCells[maxYCell]*alphaCells[maxYCell]*(Y_MAX - YCells[maxYCell]);
                    scalar tC2 = rhoCells[cellI]*alphaCells[cellI]*(YCells[cellI] - Y_MAX);
                    scalar tC = min(tC1, tC2);

                    CCells[maxYCell] += tC;
                    CCells[cellI] -= tC;

                    YCells[maxYCell] = CCells[maxYCell]/rhoCells[maxYCell]/alphaCells[maxYCell];
                    YCells[cellI] = CCells[cellI]/rhoCells[cellI]/alphaCells[cellI];
            
                    nIters++;
                }while(YCells[cellI] > Y_MAX && !allNeiDone && nIters < Y_BOUND_ITERS_MAX);        
            }
        
            if(debug)
            {
                os<< "---------------------------------------------------------------------" << endl;
            }
        }       
    }
}


label findCellInIntfcDir
(
    const fvMesh& mesh,   
    const List<scalar>& alpha1,
    const labelList& cells,
    const List<vector>& C,
    const vector& Cp,
    const vector& nf,
    const label& C_lbl,
    bool& foundCell,
    bool debug,
    OFstream& os
)
{
    scalar cosThetaMax = 0;
    label cell_lbl = 0;
    foundCell = false;
    vector meshCi;

    if(debug)
    {
        os<< "Finding cell in face normal direction" << nl 
          << "nf: " << nf
          << endl;
    }

    for(label cellI=0; cellI<cells.size(); cellI++)
    {
        label curCell = cells[cellI];
        vector Ci = C[curCell];
        vector CCi = Ci - Cp;

        if(debug)
        {
            if(cellI < mesh.nCells())
            {
                meshCi = mesh.C()[curCell];
            }
            else
            {
                meshCi = vector::one;
            }
            
            os<< "Cell " << cellI << " in stencil: " << curCell << nl
                << "C = " << Cp << "  Ci = " << Ci
                << "  alpha1i = " << alpha1[curCell] << "  meshCi = " << meshCi << nl
                << "CCi" << CCi
                << endl;
        }

        if(curCell != C_lbl)
        {            
            scalar magCCi = mag(CCi);
            if(magCCi < SMALL)
            {
                os<< "Time = " << mesh.time().timeName() << nl 
                    << "label findCellInIntfcDir" << nl
                    << "Cell " << C_lbl << "  CCi = " << CCi << "  mag(CCi) = " << magCCi << nl
                    << endl;
                magCCi += SMALL;
            }
            scalar magnf = mag(nf);
            if(magnf < SMALL)
            {
                magnf += SMALL;
            }

            scalar cosTheta = (nf & CCi)/magCCi/magnf;

            if(debug)
            {
                os<< "costheta = " << cosTheta 
                    << endl;
            }

            if(cosTheta > cosThetaMax)
            {                
                cell_lbl = curCell;
                cosThetaMax = cosTheta;
                foundCell = true;
            }
        }
    }

    if(debug)
    {
        os<< "Cell in nf direction: " << cell_lbl << nl
          << endl;
    }

    return cell_lbl;
}


label findCellInIntfcOrthDir
(
    const fvMesh& mesh,
    const List<scalar>& alpha1,
    const labelList& cells,
    const List<vector>& C,
    const vector& Cp,
    const vector& C1,
    const vector& nf,
    const label& C_lbl,
    const label& C1_lbl,
    bool& foundCell,
    bool debug,
    OFstream& os
)
{
    foundCell = false;
    scalar magCosThetaMin = 1;
    label cell_lbl = 0;

    if(debug)
    {
        os<< "Finding cell in face normal orthogonal direction" << nl 
            << "nf: " << nf
            << endl;
    }

    scalar CC1_CCi_dir;
    //vector Ci;
    vector CC1 = C1 - Cp;
    scalar magCC1 = mag(CC1);
    if(magCC1 < SMALL)
    {
        os<< "Time = " << mesh.time().timeName() << nl
            << "label findCellInIntfcOrthDir()" << nl
            << "Cell " << C_lbl << "  CC1 = " << CC1 << "  mag(CC1) = " << magCC1 << nl
            << endl;
        magCC1 += SMALL;
    }
    scalar magnf = mag(nf);
    if(magnf < SMALL)
    {
        os<< "Time = " << mesh.time().timeName() << nl
            << "label findCellInIntfcOrthDir()" << nl
            << "Cell " << C_lbl << "  nf = " << nf << "  mag(nf) = " << magnf << nl
            << endl;
        magnf += SMALL;
    }
    vector nfXCC1 = nf ^ CC1;
    nfXCC1 /= magCC1*magnf;
    vector CCi;
    vector meshCi;

    for(label cellI=0; cellI<cells.size(); cellI++)
    {
        label curCell = cells[cellI];
        CCi = C[curCell] - Cp;
        scalar magCCi = mag(CCi);
        if(magCCi < SMALL)
        {
            os<< "Time = " << mesh.time().timeName() << nl
                << "label findCellInIntfcOrthDir()" << nl
                << "Cell " << C_lbl << "  CCi = " << CCi << "  mag(CCi) = " << magCCi << nl
                << endl;
            magCCi += SMALL;
        }

        CC1_CCi_dir = nfXCC1 & (nf ^ CCi);
        CC1_CCi_dir /= magCCi*magnf;

        if(debug)
        {
            if(cellI < mesh.nCells())
            {
                meshCi = mesh.C()[curCell];
            }
            else
            {
                meshCi = vector::one;
            }

            os<< "Cell " << cellI << " in stencil: " << curCell << nl
                << "C = " << Cp << "  C1 = " << C1 << "  Ci = " << C[curCell] << nl
                << "alpha1i = " << alpha1[curCell] << "  meshCi = " << meshCi << nl
                << "CC1 = " << CC1 << "  CCi" << CCi << "  CC1_CCi_dir = " << CC1_CCi_dir
                << endl;
        }

        scalar cosTheta = (CC1 & CCi)/mag(CCi)/mag(CC1);

        if(curCell != C1_lbl && CC1_CCi_dir < 0 && cosTheta >= 0)
        {
            foundCell = true;
            //CCi /= mag(CCi);
            //scalar cosTheta = nf & CCi;            

            if(debug)
            {
                os<< "costheta = " << cosTheta << endl;
            }

            if(mag(cosTheta) < magCosThetaMin)
            {
                cell_lbl = curCell;
                magCosThetaMin = mag(cosTheta);
            }
        }
    }

    if(debug)
    {
        os<< "Cell in nf orthogonal direction: " << cell_lbl << nl
            << endl;
    }

    return cell_lbl;
}


void calc_cell_intfcGrad_coeffs
(
    const fvMesh& mesh,
    const label& curCell_lbl,
    const vector& nf,
    const vector& C_intfc,
    const List<List<scalar> >& Y,
    const List<scalar>& alpha1,
    const List<vector>& C,
    const labelList& curCellsAll,
    const label& nSpecies,
    const scalar& MIN_ALPHA_2PH,
    const label& phaseLbl,
    scalar& dn,
    List<scalar>& Yeff,    
    bool debug, 
    OFstream& os
)
{
    if(debug)
    {
        os<< nl
          << "Calculating cell intfc grad weights in cell " << curCell_lbl << nl
          << endl;
    }

    scalar MAX_ALPHA_2PH = 1 - MIN_ALPHA_2PH;
    
    labelList curCells(curCellsAll.size());
    label n_ph = 0;

    if(debug)
    {
        os<< "Reducing cell stencil for phase " << phaseLbl << nl
            << "Full cell stencil" << nl
            << curCellsAll << endl;
    }

    if(phaseLbl == 1)
    {
        for(label cellI=0; cellI<curCellsAll.size(); cellI++)
        {            
            label cellI_lbl = curCellsAll[cellI];

            if(debug)
            {
                os<< "cellI: " << cellI << "  cellI_lbl: " << cellI_lbl << nl
                    << "alpha1 = " << alpha1[cellI_lbl] << "  MIN_ALPHA_2PH = " << MIN_ALPHA_2PH
                    << endl;
            }

            if(alpha1[cellI_lbl] > MIN_ALPHA_2PH && cellI_lbl != curCell_lbl)
            {
                curCells[n_ph++] = cellI_lbl;
            }            
        }        
    }
    else
    {
        for(label cellI=0; cellI<curCellsAll.size(); cellI++)
        {
            label cellI_lbl = curCellsAll[cellI];

            if(debug)
            {
                os<< "cellI: " << cellI << "  cellI_lbl: " << cellI_lbl << nl
                    << "alpha1 = " << alpha1[cellI_lbl] << "  MAX_ALPHA_2PH = " << MAX_ALPHA_2PH
                    << endl;
            }

            if(alpha1[cellI_lbl] < MAX_ALPHA_2PH && cellI_lbl != curCell_lbl)
            {
                curCells[n_ph++] = curCellsAll[cellI];
            }            
        }
    }
    
    curCells.setSize(n_ph);

    if(debug)
    {
        os<< "Cell reduced stencil" << nl
            << curCells << nl
            << endl;
    }

    // suffix 1: direction closest to nf
    // suffix 2: direction closest orthogonal to nf in 2-D
    // further improvements needed for 3-D calculation
    bool foundCell1 = false;
    bool foundCell2 = false;
    vector Cp = C_intfc;
    label C1_lbl = findCellInIntfcDir(mesh,alpha1,curCells,C,Cp,nf,curCell_lbl,foundCell1,debug,os);
    if(debug)
    {
        os<< "C1_lbl: " << C1_lbl << "  found C1: " << foundCell1
            << endl;
    }
    if(foundCell1)
    {
        vector C1 = C[C1_lbl];
        vector t1 = C1 - Cp;
        scalar magt1;    
        magt1 = mag(t1);
        scalar magnf;
        magnf = mag(nf);
        if(magt1 < SMALL)
        {
            os<< "Time = " << mesh.time().timeName() << nl 
                << "Cell " << curCell_lbl << "  Cp = " << Cp << "  C1 = " << C1 << "  magt1 = " << magt1 << nl
                << endl;
            magt1 += SMALL;
        }
        if(magnf < SMALL)
        {
            os<< "Time = " << mesh.time().timeName() << nl 
                << "Cell " << curCell_lbl << "  nf = " << nf << "  mag(nf) = " << magnf << nl
                << endl;
            magnf += SMALL;
        }
        scalar costheta1 = (nf & t1)/magt1/magnf;
        if(costheta1 > 1)
        {        
            os<< "Time = " << mesh.time().timeName() << nl 
                << "Cell " << curCell_lbl << "  Cp = " << Cp << "C1 cell " << C1_lbl << "  C1 = " << C1 << nl
                << "t1 = " << t1 << "  mag(t1) = " << magt1 << "  nf = " << nf << "  mag(nf) = " << mag(nf) 
                << "  costheta1 = " << costheta1 << nl 
                << endl;
            costheta1 = 1;
        }
        scalar theta1 = acos(costheta1);
        if(debug)
        {
            os<< "Cp = " << Cp << "  C1 = " << C1 << "  t1 = " << t1 << "  mag(t1) = " << magt1
                << "  costheta1 = " << costheta1 << "  theta1 = " << theta1                      
                << endl;
        }

        if(theta1 > 1E-3)
        {
            label C2_lbl = findCellInIntfcOrthDir(mesh,alpha1,curCells,C,Cp,C1,nf,curCell_lbl,C1_lbl,foundCell2,debug,os);
            if(debug)
            {
                os<< "C2_lbl: " << C2_lbl << "  found C2: " << foundCell2
                    << endl;
            }   

            if(foundCell2)
            {                                        
                vector C2 = C[C2_lbl];    
                vector t2 = C2 - Cp;
                scalar magt2;
                magt2 = mag(t2);
                if(magt2 < SMALL)
                {
                    os<< "Time = " << mesh.time().timeName() << nl 
                        << "Cell " << curCell_lbl << "  Cp = " << Cp << "  C2 = " << C2 << "  magt2 = " << magt2 << nl
                        << endl;
                    magt2 += SMALL;
                }    
                scalar magt1t2 = magt1*magt2;
                if(magt1t2 < SMALL)
                {
                    os<< "Time = " << mesh.time().timeName() << nl 
                        << "Cell " << curCell_lbl << "  magt1 = " << magt1 << "  magt2 = " << magt2 << "  magt1t2 = " << magt1t2 << nl
                        << endl;
                    magt1t2 += SMALL;
                }        
        
                scalar costheta2 = (nf & t2)/magt2/magnf;
                if(costheta2 > 1)
                {                
                    os<< "Time = " << mesh.time().timeName() << nl 
                        << "Cell " << curCell_lbl << "  Cp = " << Cp << "C2 cell " << C2_lbl << "  C2 = " << C2 << nl
                        << "t2 = " << t2 << "  mag(t2) = " << magt2 << "  nf = " << nf << "  mag(nf) = " << mag(nf) 
                        << "  costheta2 = " << costheta2 << nl 
                        << endl;
                    costheta2 = 1;
                }
                scalar theta2 = acos(costheta2);
                if(debug)
                {
                    os<< "t2 = " << t2 << "  mag(t2) = " << magt2
                        << "  costheta2 = " << costheta2 << "  theta2 = " << theta2                      
                        << endl;
                }
                scalar theta2_sign = -((nf ^ t1) & (nf ^ t2))/magt1t2/magnf/magnf;
                if(debug)
                {
                    os<< "theta2_sign = " << theta2_sign            
                        << endl;
                }

                scalar sintheta12 = sin(theta1 + theta2);

                if(theta1 + theta2 > constant::mathematical::pi)
                {
                    os<< "Time = " << mesh.time().timeName() << nl 
                        << "Cell " << curCell_lbl << "  nf" << nf << "  Cp = " << Cp << "  phaseLabel: " << phaseLbl << "  cell alpha1 = " << alpha1[curCell_lbl] << "  cell meshC = " << mesh.C()[curCell_lbl] << nl 
                        << "C1_lbl: " << C1_lbl << "  C1 = " << C1 << "  C2_lbl: " << C2_lbl << "  C2 = " << C2 << nl
                        << "magt1 = " << magt1 << "  magt2 = " << magt2 << "  magt1t2 = " << magt1t2 << nl
                        << "theta1 = " << theta1 << "  theta2 = " << theta2 << "  (theta1+theta2) = " << theta1 + theta2 << " sin(theta1+theta2) = " << sintheta12 << nl
                        << endl;

                    C1_lbl = findCellInIntfcDir(mesh,alpha1,curCells,C,Cp,nf,curCell_lbl,foundCell1,1,os);
                    C2_lbl = findCellInIntfcOrthDir(mesh,alpha1,curCells,C,Cp,C1,nf,curCell_lbl,C1_lbl,foundCell2,1,os);

                    os<< endl;
                }

                if(sintheta12 < SMALL)
                {
                    os<< "Time = " << mesh.time().timeName() << nl 
                        << "Cell " << curCell_lbl << "  theta1 = " << theta1 << "  theta2 = " << theta2 << "  sintheta12 = " << sintheta12 << nl
                        << endl;
                    sintheta12 += SMALL;
                }

                scalar alpha = sin(theta2)/sintheta12;        
                scalar beta = sin(theta1)/sintheta12;
                scalar at2bt1 = alpha*magt2 + beta*magt1;
                if(mag(at2bt1) < SMALL)
                {
                    os<< "Time = " << mesh.time().timeName() << nl 
                        << "Cell " << curCell_lbl << "  theta1 = " << theta1 << "  theta2 = " << theta2 << "  alpha = " << alpha 
                        << "  beta = " << beta << "  at2bt1 = " << at2bt1 << nl
                        << endl;
                    at2bt1 += SMALL;
                }

                if(debug)
                {
                    os<< "alpha = " << alpha << "  beta = " << beta << nl
                        << "at2bt1 = " << at2bt1
                        << endl;
                }
    
                dn = magt1t2/at2bt1;

                if(debug)
                {
                    os<< "dn = " << dn
                        << endl;
                }

                scalar Y1;
                scalar Y2;
                for(label i=0; i<nSpecies; i++)
                {
                    Y1 = Y[i][C1_lbl];
                    Y2 = Y[i][C2_lbl];
                    Yeff[i] = (alpha*magt2*Y1 + beta*magt1*Y2)/at2bt1;

                    if(debug)
                    {
                        os<< "Species index: " << i << nl
                            << "Y1 = " << Y1 << "  Y2 = " << Y2 << "  Yeff = " << Yeff[i]
                            << endl;
                    }
                }             

                if(debug)
                {
                    os<< endl;        
                }
            }
            else
            {
                if(costheta1 < 1E-3)
                {        
                    os<< "Time = " << mesh.time().timeName() << nl 
                        << "Cell " << curCell_lbl << "  Cp = " << Cp << "  C1 cell " << C1_lbl << "  C1 = " << C1 << nl
                        << "t1 = " << t1 << "  mag(t1) = " << magt1 << "  nf = " << nf << "  mag(nf) = " << mag(nf) 
                        << "  costheta1 = " << costheta1 << nl 
                        << endl;
                    costheta1 = 1E-3;
                }
                dn = magt1/costheta1;
                if(debug)
                {
                    os<< "dn = " << dn
                        << endl;
                }
                for(label i=0; i<nSpecies; i++)
                {             
                    Yeff[i] = Y[i][C1_lbl];
                    if(debug)
                    {
                        os<< "Species index: " << i << nl
                            << "Y1 = " << Y[i][C1_lbl] << "  Yeff = " << Yeff[i]
                            << endl;
                    }
                }
            }
        }
        else
        {
            dn = magt1;
            if(debug)
            {
                os<< "dn = " << dn
                    << endl;
            }
            for(label i=0; i<nSpecies; i++)
            {             
                Yeff[i] = Y[i][C1_lbl];
                if(debug)
                {
                    os<< "Species index: " << i << nl
                        << "Y1 = " << Y[i][C1_lbl] << "  Yeff = " << Yeff[i]
                        << endl;
                }
            }
        }
    }
    else
    {
        os<< "Time = " << mesh.time().timeName() << nl 
            << "Cell " << curCell_lbl << ": Cell in interface normal direction not found! Using curCell centroid" << nl
            << endl;
        /*
        vector C1 = C[curCell_lbl];
        vector t1 = C1 - Cp;
        scalar magt1;    
        magt1 = mag(t1);
        scalar magnf;
        magnf = mag(nf);
        */
	scalar cellVol = mesh.V()[curCell_lbl];
        scalar min_dn = pow(cellVol, 1.0/3.0);
        /*
            if(magt1 < min_dn)
            {
            os<< "Time = " << mesh.time().timeName() << nl 
            << "Cell " << curCell_lbl << ": curCell centroid too close to interface!" << nl 
            << "Cp = " << Cp << "  C1 = " << C1 << "  magt1 = " << magt1 << nl
            << endl;
            if(magt1 < SMALL)
            {
            magt1 += SMALL;
            }
            t1 *= min_dn/magt1;
            magt1 = mag(t1);
            }
            if(magnf < SMALL)
            {
            os<< "Time = " << mesh.time().timeName() << nl 
            << "Cell " << curCell_lbl << "  nf = " << nf << "  mag(nf) = " << magnf << nl
            << endl;
            magnf += SMALL;
            }
            scalar costheta1 = (nf & t1)/magt1/magnf;
            if(costheta1 > 1)
            {        
            os<< "Time = " << mesh.time().timeName() << nl 
            << "Cell " << curCell_lbl << "  Cp = " << Cp << "C1 cell " << C1_lbl << "  C1 = " << C1 << nl
            << "t1 = " << t1 << "  mag(t1) = " << magt1 << "  nf = " << nf << "  mag(nf) = " << mag(nf) 
            << "  costheta1 = " << costheta1 << nl 
            << endl;
            costheta1 = 1;
            }
            */
        dn = min_dn;
        if(debug)
        {
            os<< "dn = " << dn
                << endl;
        }
        for(label i=0; i<nSpecies; i++)
        {             
            Yeff[i] = Y[i][curCell_lbl];
            if(debug)
            {
                os<< "Species index: " << i << nl
                    << "Y1 = " << Y[i][curCell_lbl] << "  Yeff = " << Yeff[i]
                    << endl;
            }
        }
    }
    
    if(dn < SMALL)
    {
        dn += SMALL;
    }
    if(debug)
    {
        os<< "dn stab = " << dn << endl;
    }
}


void calc_Js
(
    const fvMesh& mesh,
    const labelListList& cellStencil,    
    const List<List<scalar> >& Y1_flatFld,
    const List<List<scalar> >& Y0_flatFld,
    const List<scalar>& alpha1_flatFld,
    const volScalarField& alpha1,
    const volScalarField& rho1,
    const volScalarField& rho0,
    const PtrList<volScalarField>& D1,
    const PtrList<volScalarField>& D0,
    const List<vector>& C_ph1_flatFld,
    const List<vector>& C_ph0_flatFld,
    const volVectorField& C_intfc,
    const volScalarField& A_intfc,
    const volVectorField& nHat,
    const PtrList<volScalarField>& Ys1,
    const PtrList<volScalarField>& Ys0,
    PtrList<volScalarField>& Js1,
    PtrList<volScalarField>& Js0,    
    const label& nSpecies,
    const scalar& ALPHA_2PH_MIN,
    const bool debug,
    OFstream& os
)
{
    scalar ALPHA_2PH_MAX = 1 - ALPHA_2PH_MIN;

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();

    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();
    const surfaceVectorField& Cf = mesh.Cf();

    const scalarField& alpha1Cells = alpha1.internalField();    
    const vectorField& C_intfcCells = C_intfc.internalField();
    const scalarField& A_intfcCells = A_intfc.internalField();
    const vectorField& nHatCells = nHat.internalField();
    const scalarField& rho1Cells = rho1.internalField();
    const scalarField& rho0Cells = rho0.internalField();

    if(debug)
    {
        os<< "-------------------------------------------------------------------------" << nl
            << "Interfacial Species Flux Calculation" << nl
            << "-------------------------------------------------------------------------" << nl
            << nl
            << "-------------------------------------------------------------------------" << nl
            << "Internal cells" << nl
            << "-------------------------------------------------------------------------" << nl
            << endl;
    }

    //Js for all interface cells
    forAll(alpha1Cells, cellI)
    {
        scalar alpha1_cellI = alpha1Cells[cellI];        

        if(alpha1_cellI > ALPHA_2PH_MIN && alpha1_cellI < ALPHA_2PH_MAX)
        {
            vector nf = nHatCells[cellI];
            vector C_intfc_cellI = C_intfcCells[cellI];
            scalar A_intfc_cellI = A_intfcCells[cellI];
            labelList curCellsAll = cellStencil[cellI];

            if(debug)
            {
                os<< "Cell: " << cellI << "  alpha1 = " << alpha1_cellI << endl;
            }

            if(debug)
            {
                os<< "nf = " << nf << "  C_intfc = " << C_intfc_cellI << "  A_intfc = " << A_intfc_cellI << nl
                    << "C_ph1 = " << C_ph1_flatFld[cellI] << "  C_ph0 = " << C_ph0_flatFld[cellI] << nl
                    << "rho1 = " << rho1Cells[cellI] << "rho0 = " << rho0Cells[cellI] << endl;
            }

            scalar dn1;
            List<scalar> Yeff1(nSpecies);
            scalar dn0;
            List<scalar> Yeff0(nSpecies);
            scalar intfcGradi_cellI;

            //phase-1
            //ensure nf direction is into the phase
            calc_cell_intfcGrad_coeffs(mesh, cellI, nf, C_intfc_cellI, Y1_flatFld, alpha1_flatFld, C_ph1_flatFld, curCellsAll, nSpecies, ALPHA_2PH_MIN, 1, dn1, Yeff1, debug, os);
            if(debug)
            {
                os<< "dn1 = " << dn1 << endl;
            }
            if(dn1 < SMALL)
            {
                dn1 += SMALL;
            }
            if(debug)
            {
                os<< "dn1 stab = " << dn1 << endl;
            }

            for(label i=0; i<nSpecies; i++)
            {
                intfcGradi_cellI = (Yeff1[i] - Ys1[i].internalField()[cellI])/dn1;
                
                Js1[i].internalField()[cellI] = -A_intfc_cellI*rho1Cells[cellI]*D1[i].internalField()[cellI]*intfcGradi_cellI;

                if(debug)
                {
                    os<< "species: " << i << "  Ys1i = " << Ys1[i].internalField()[cellI] << "  Yeff1i = " << Yeff1[i] << "  intfcGrad1i = " << intfcGradi_cellI << "  D1i = " << D1[i].internalField()[cellI] << "  Js1i = " <<  Js1[i].internalField()[cellI] << endl;
                }
            }
            
            //phase-0
            //ensure nf direction is into the phase
            //then reverse nf again for Js0 calculation            
            calc_cell_intfcGrad_coeffs(mesh, cellI, -nf, C_intfc_cellI, Y0_flatFld, alpha1_flatFld, C_ph0_flatFld, curCellsAll, nSpecies, ALPHA_2PH_MIN, 0, dn0, Yeff0, debug, os);
            if(debug)
            {
                os<< "dn0 = " << dn0 << endl;
            }
            if(dn0 < SMALL)
            {
                dn0 += SMALL;
            }
            if(debug)
            {
                os<< "dn0 stab = " << dn0 << endl;
            }

            for(label i=0; i<nSpecies; i++)
            {
                intfcGradi_cellI = (Yeff0[i] - Ys0[i].internalField()[cellI])/dn0;
                
                Js0[i].internalField()[cellI] = A_intfc_cellI*rho0Cells[cellI]*D0[i].internalField()[cellI]*intfcGradi_cellI;

                if(debug)
                {
                    os<< "species: " << i << "  Ys0i = " << Ys0[i].internalField()[cellI] << "  Yeff0i = " << Yeff0[i] << "  intfcGrad0i = " << intfcGradi_cellI << "  D0i = " << D0[i].internalField()[cellI] << "  Js0i = " <<  Js0[i].internalField()[cellI] << endl;
                }
            }

            if(debug)
            {
                os<< "-------------------------------------------------------------------------" << endl;

                for(label i=0; i<nSpecies; i++)
                {
                    os<< "Js1[" << i << "] = " << Js1[i].internalField()[cellI] << "  Js0[" << i << "] = " << Js0[i].internalField()[cellI]
                        << endl;
                }

                os<< "-------------------------------------------------------------------------" << nl
                    << endl;
            }
        }
        else
        {
            for(label i=0; i<nSpecies; i++)
            { 
                Js0[i].internalField()[cellI] = 0;
                Js1[i].internalField()[cellI] = 0;
            }
        }        
    }

    if(debug)
    {
        os<< "-------------------------------------------------------------------------" << nl
            << "Internal faces" << nl
            << "-------------------------------------------------------------------------" << nl
            << endl;
    }

    //Js for cell faces that are almost a phase inteface
    //internal faces
    for(label faceI=0; faceI<mesh.nInternalFaces(); faceI++)
    {
        label faceOwn = own[faceI];
        label faceNei = nei[faceI];
        scalar alpha1Own = alpha1Cells[faceOwn];
        scalar alpha1Nei = alpha1Cells[faceNei];        

        if(alpha1Own >= ALPHA_2PH_MAX && alpha1Nei <= ALPHA_2PH_MIN)
        {                        
            scalar A_intfc_cellI = magSf[faceI];
            vector nf = -Sf[faceI]/A_intfc_cellI;
            vector C_intfc_cellI = Cf[faceI];
            labelList curCellsAll = cellStencil[faceOwn];

            if(debug)
            {
                os<< "Face: " << faceI << nl
                    << "own: " << faceOwn << "  alpha1Own = " << alpha1Own << "  nei: " << faceNei << "  alpha1Nei = " << alpha1Nei << endl;
            }

            if(debug)
            {
                os<< "ph-1 cell: " << faceOwn << "  ph-0 cell: " << faceNei << nl
                    << "nfOwn = " << nHatCells[faceOwn] << "  nfNei = " << nHatCells[faceNei] << nl
                    << "C_intfcOwn = " << C_intfcCells[faceOwn] << "  C_intfcNei = " << C_intfcCells[faceNei] << nl
                    << "A_intfcOwn = " << A_intfcCells[faceOwn] << "  A_intfcNei = " << A_intfcCells[faceNei] << nl
                    << "nf = " << nf << "  C_intfc = " << C_intfc_cellI << "  A_intfc = " << A_intfc_cellI << nl
                    << "C_ph1Own = " << C_ph1_flatFld[faceOwn] << "  C_ph0Own = " << C_ph0_flatFld[faceOwn] << nl
                    << "rho1Own = " << rho1Cells[faceOwn] << "  rho0Own = " << rho0Cells[faceOwn] << nl
                    << "C_ph1Nei = " << C_ph1_flatFld[faceNei] << "  C_ph0Nei = " << C_ph0_flatFld[faceNei] << nl
                    << "rho1Nei = " << rho1Cells[faceNei] << "  rho0Nei = " << rho0Cells[faceNei] << endl;
            }

            scalar dn1;
            List<scalar> Yeff1(nSpecies);
            scalar dn0;
            List<scalar> Yeff0(nSpecies);
            scalar intfcGradi_cellI;

            //phase-1
            //ensure nf direction is into the phase
            calc_cell_intfcGrad_coeffs(mesh, faceOwn, nf, C_intfc_cellI, Y1_flatFld, alpha1_flatFld, C_ph1_flatFld, curCellsAll, nSpecies, ALPHA_2PH_MIN, 1, dn1, Yeff1, debug, os);
            if(debug)
            {
                os<< "dn1 = " << dn1 << endl;
            }
            if(dn1 < SMALL)
            {
                dn1 += SMALL;
            }
            if(debug)
            {
                os<< "dn1 stab = " << dn1 << endl;
            }

            for(label i=0; i<nSpecies; i++)
            {
                intfcGradi_cellI = (Yeff1[i] - Ys1[i].internalField()[faceOwn])/dn1;
                
                Js1[i].internalField()[faceOwn] = -A_intfc_cellI*rho1Cells[faceOwn]*D1[i].internalField()[faceOwn]*intfcGradi_cellI;
                Js1[i].internalField()[faceNei] = Js1[i].internalField()[faceOwn];

                if(debug)
                {
                    os<< "species: " << i << "  Ys1iOwn = " << Ys1[i].internalField()[faceOwn] << "  Yeff1i = " << Yeff1[i] << "  intfcGrad1iOwn = " << intfcGradi_cellI << "  D1iOwn = " << D1[i].internalField()[faceOwn] << "  Js1iOwn = " <<  Js1[i].internalField()[faceOwn] << "  Js1iNei = " <<  Js1[i].internalField()[faceNei] << endl;
                }
            }
            
            //phase-0
            //ensure nf direction is into the phase
            //then reverse nf again for Js0 calculation
            curCellsAll = cellStencil[faceNei];
            calc_cell_intfcGrad_coeffs(mesh, faceNei, -nf, C_intfc_cellI, Y0_flatFld, alpha1_flatFld, C_ph0_flatFld, curCellsAll, nSpecies, ALPHA_2PH_MIN, 0, dn0, Yeff0, debug, os);
            if(debug)
            {
                os<< "dn0 = " << dn0 << endl;
            }
            if(dn0 < SMALL)
            {
                dn0 += SMALL;
            }
            if(debug)
            {
                os<< "dn0 stab = " << dn0 << endl;
            }

            for(label i=0; i<nSpecies; i++)
            {
                intfcGradi_cellI = (Yeff0[i] - Ys0[i].internalField()[faceNei])/dn0;
                
                Js0[i].internalField()[faceNei] = A_intfc_cellI*rho0Cells[faceNei]*D0[i].internalField()[faceNei]*intfcGradi_cellI;
                Js0[i].internalField()[faceOwn] = Js0[i].internalField()[faceNei];

                if(debug)
                {
                    os<< "species: " << i << "  Ys0iNei = " << Ys0[i].internalField()[faceNei] << "  Yeff0i = " << Yeff0[i] << "  intfcGrad0iOwn = " << intfcGradi_cellI << "  D0iNei = " << D0[i].internalField()[faceNei] << "  Js0iNei = " <<  Js0[i].internalField()[faceNei] << "  Js0iOwn = " <<  Js0[i].internalField()[faceOwn] << endl;
                }
            }
        } 

        if(alpha1Own <= ALPHA_2PH_MIN && alpha1Nei >= ALPHA_2PH_MAX)
        {
            scalar A_intfc_cellI = magSf[faceI];
            vector nf = Sf[faceI]/A_intfc_cellI;
            vector C_intfc_cellI = Cf[faceI];
            labelList curCellsAll = cellStencil[faceNei];

            if(debug)
            {
                os<< "Face: " << faceI << nl
                    << "own: " << faceOwn << "  alpha1Own = " << alpha1Own << "  nei: " << faceNei << "  alpha1Nei = " << alpha1Nei << endl;
            }

            if(debug)
            {
                os<< "ph-1 cell: " << faceNei << "  ph-0 cell: " << faceOwn << nl
                    << "nfOwn = " << nHatCells[faceOwn] << "  nfNei = " << nHatCells[faceNei] << nl
                    << "C_intfcOwn = " << C_intfcCells[faceOwn] << "  C_intfcNei = " << C_intfcCells[faceNei] << nl
                    << "A_intfcOwn = " << A_intfcCells[faceOwn] << "  A_intfcNei = " << A_intfcCells[faceNei] << nl
                    << "nf = " << nf << "  C_intfc = " << C_intfc_cellI << "  A_intfc = " << A_intfc_cellI << nl
                    << "C_ph1Own = " << C_ph1_flatFld[faceOwn] << "  C_ph0Own = " << C_ph0_flatFld[faceOwn] << nl
                    << "rho1Own = " << rho1Cells[faceOwn] << "  rho0Own = " << rho0Cells[faceOwn] << nl
                    << "C_ph1Nei = " << C_ph1_flatFld[faceNei] << "  C_ph0Nei = " << C_ph0_flatFld[faceNei] << nl
                    << "rho1Nei = " << rho1Cells[faceNei] << "  rho0Nei = " << rho0Cells[faceNei] << endl;
            }

            scalar dn1;
            List<scalar> Yeff1(nSpecies);
            scalar dn0;
            List<scalar> Yeff0(nSpecies);
            scalar intfcGradi_cellI;

            //phase-1
            //ensure nf direction is into the phase
            calc_cell_intfcGrad_coeffs(mesh, faceNei, nf, C_intfc_cellI, Y1_flatFld, alpha1_flatFld, C_ph1_flatFld, curCellsAll, nSpecies, ALPHA_2PH_MIN, 1, dn1, Yeff1, debug, os);
            if(debug)
            {
                os<< "dn1 = " << dn1 << endl;
            }
            if(dn1 < SMALL)
            {
                dn1 += SMALL;
            }
            if(debug)
            {
                os<< "dn1 stab = " << dn1 << endl;
            }            

            for(label i=0; i<nSpecies; i++)
            {
                intfcGradi_cellI = (Yeff1[i] - Ys1[i].internalField()[faceNei])/dn1;
                
                Js1[i].internalField()[faceNei] = -A_intfc_cellI*rho1Cells[faceNei]*D1[i].internalField()[faceNei]*intfcGradi_cellI;
                Js1[i].internalField()[faceOwn] = Js1[i].internalField()[faceNei];

                if(debug)
                {
                    os<< "species: " << i << "  Ys1iNei = " << Ys1[i].internalField()[faceNei] << "  Yeff1i = " << Yeff1[i] << "  intfcGrad1iOwn = " << intfcGradi_cellI << "  D1iNei = " << D1[i].internalField()[faceNei] << "  Js1iNei = " <<  Js1[i].internalField()[faceNei] << "  Js1iOwn = " <<  Js1[i].internalField()[faceOwn] << endl;
                }
            }
            
            //phase-0
            //ensure nf direction is into the phase
            //then reverse nf again for Js0 calculation
            curCellsAll = cellStencil[faceOwn];
            calc_cell_intfcGrad_coeffs(mesh, faceOwn, -nf, C_intfc_cellI, Y0_flatFld, alpha1_flatFld, C_ph0_flatFld, curCellsAll, nSpecies, ALPHA_2PH_MIN, 0, dn0, Yeff0, debug, os);
            if(debug)
            {
                os<< "dn0 = " << dn0 << endl;
            }
            if(dn0 < SMALL)
            {
                dn0 += SMALL;
            }
            if(debug)
            {
                os<< "dn0 stab = " << dn0 << endl;
            }
            
            for(label i=0; i<nSpecies; i++)
            {
                intfcGradi_cellI = (Yeff0[i] - Ys0[i].internalField()[faceOwn])/dn0;
                
                Js0[i].internalField()[faceOwn] = A_intfc_cellI*rho0Cells[faceOwn]*D0[i].internalField()[faceOwn]*intfcGradi_cellI;
                Js0[i].internalField()[faceNei] = Js0[i].internalField()[faceOwn];

                if(debug)
                {
                    os<< "species: " << i << "  Ys0iOwn = " << Ys0[i].internalField()[faceOwn] << "  Yeff0i = " << Yeff0[i] << "  intfcGrad0iOwn = " << intfcGradi_cellI << "  D0iOwn = " << D0[i].internalField()[faceOwn] << "  Js0iOwn = " <<  Js0[i].internalField()[faceOwn] << "  Js0iNei = " <<  Js0[i].internalField()[faceNei] << endl;
                }
            }
        } 
    }    

    //boundary faces
    if(debug)
    {
        os<< "-------------------------------------------------------------------------" << nl
            << "Boundary coupled faces step 1" << nl
            << "-------------------------------------------------------------------------" << nl
            << endl;
    }
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const wordList& patchNames = patches.names();
    const label nBnd = mesh.nFaces() - mesh.nInternalFaces();

    List<List<scalar> > JsOwn(nSpecies);
    List<label> phOwn(nBnd);
    List<List<scalar> > JsNei(nSpecies);
    List<label> phNei(nBnd);
    for(label i=0; i<nSpecies; i++)
    {
        JsOwn[i].setSize(nBnd);
        JsNei[i].setSize(nBnd);
    }

    forAll(Js1[0].boundaryField(), patchI)
    {
        const polyPatch& pp = patches[patchI];                

        if(pp.coupled())
        {
            if(debug)
            {
                os<< "-------------------------------------------------------------------------" << nl
                    << "Patch: " << patchNames[patchI] << nl
                    << "-------------------------------------------------------------------------" << nl
                    << endl;
            }

            label faceI = pp.start();
            label bndFaceI = faceI - mesh.nInternalFaces();

            const scalarField& alpha1NeiFld = alpha1.boundaryField()[patchI].patchNeighbourField();
            const vectorField& nHatNeiFld = nHat.boundaryField()[patchI].patchNeighbourField();
            const vectorField& C_intfcNeiFld = C_intfc.boundaryField()[patchI].patchNeighbourField();
            const scalarField& A_intfcNeiFld = A_intfc.boundaryField()[patchI].patchNeighbourField();
            const scalarField& rho1NeiFld = rho1.boundaryField()[patchI].patchNeighbourField();            
            const scalarField& rho0NeiFld = rho0.boundaryField()[patchI].patchNeighbourField();
            const fvsPatchVectorField& pSf = Sf.boundaryField()[patchI];
            const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchI];
            const fvsPatchVectorField& pCf = Cf.boundaryField()[patchI];

            forAll(Js1[0].boundaryField()[patchI], fcI)
            {
                label faceOwn = own[faceI];
                scalar alpha1Own = alpha1Cells[faceOwn];
                scalar alpha1Nei = alpha1NeiFld[fcI];                

                if(alpha1Own >= ALPHA_2PH_MAX && alpha1Nei <= ALPHA_2PH_MIN)
                {
                    if(debug)
                    {
                        os<< "Face: " << faceI << "  patch face index: " << fcI << "  bnd face index: " << bndFaceI << nl
                            << "own: " << faceOwn << "  alpha1Own = " << alpha1Own << "  alpha1Nei = " << alpha1Nei << endl;
                    }

                    phOwn[bndFaceI] = 1;
                    phNei[bndFaceI] = 1;
                    
                    scalar A_intfc_cellI = pMagSf[fcI];
                    vector nf = -pSf[fcI]/A_intfc_cellI;
                    vector C_intfc_cellI = pCf[fcI];
                    const labelList& curCellsAll = cellStencil[faceOwn];

                    if(debug)
                    {
                        os<< "ph-1 cell: " << faceOwn << nl
                            << "nfOwn = " << nHatCells[faceOwn] << "  nfNei = " << nHatNeiFld[fcI] << nl
                            << "C_intfcOwn = " << C_intfcCells[faceOwn] << "  C_intfcNei = " << C_intfcNeiFld[fcI] << nl
                            << "A_intfcOwn = " << A_intfcCells[faceOwn] << "  A_intfcNei = " << A_intfcNeiFld[fcI] << nl
                            << "nf = " << nf << "  C_intfc = " << C_intfc_cellI << "  A_intfc = " << A_intfc_cellI << nl
                            << "C_ph1Own = " << C_ph1_flatFld[faceOwn] << "  C_ph0Own = " << C_ph0_flatFld[faceOwn] << nl
                            << "rho1Own = " << rho1Cells[faceOwn] << "  rho0Own = " << rho0Cells[faceOwn] << nl                            
                            << "rho1Nei = " << rho1NeiFld[fcI] << "  rho0Nei = " << rho0NeiFld[fcI] << endl;
                    }

                    scalar dn1;
                    List<scalar> Yeff1(nSpecies);                    
                    scalar intfcGradi_cellI;

                    //phase-1
                    //ensure nf direction is into the phase
                    calc_cell_intfcGrad_coeffs(mesh, faceOwn, nf, C_intfc_cellI, Y1_flatFld, alpha1_flatFld, C_ph1_flatFld, curCellsAll, nSpecies, ALPHA_2PH_MIN, 1, dn1, Yeff1, debug, os);
                    if(debug)
                    {
                        os<< "dn1 = " << dn1 << endl;
                    }
                    if(dn1 < SMALL)
                    {
                        dn1 += SMALL;
                    }
                    if(debug)
                    {
                        os<< "dn1 stab = " << dn1 << endl;
                    }
             
                    for(label i=0; i<nSpecies; i++)
                    {
                        intfcGradi_cellI = (Yeff1[i] - Ys1[i].internalField()[faceOwn])/dn1;
                
                        JsOwn[i][bndFaceI] = -A_intfc_cellI*rho1Cells[faceOwn]*D1[i].internalField()[faceOwn]*intfcGradi_cellI;
                        JsNei[i][bndFaceI] = JsOwn[i][bndFaceI];

                        if(debug)
                        {
                            os<< "species: " << i << "  Ys1iOwn = " << Ys1[i].internalField()[faceOwn] << "  Yeff1i = " << Yeff1[i] << "  intfcGrad1iOwn = " << intfcGradi_cellI << "  D1iOwn = " << D1[i].internalField()[faceOwn] << "  Js1iOwn = " <<  JsOwn[i][bndFaceI] << "  Js1iNei = " <<  JsNei[i][bndFaceI] << endl;
                        }
                    }
                }

                if(alpha1Own <= ALPHA_2PH_MIN && alpha1Nei >= ALPHA_2PH_MAX)
                {
                    if(debug)
                    {
                        os<< "Face: " << faceI << "  patch face index: " << fcI << "  bnd face index: " << bndFaceI << nl
                            << "own: " << faceOwn << "  alpha1Own = " << alpha1Own << "  alpha1Nei = " << alpha1Nei << endl;
                    }

                    phOwn[bndFaceI] = 0;
                    phNei[bndFaceI] = 0;
  
                    scalar A_intfc_cellI = pMagSf[fcI];
                    vector nf = pSf[fcI]/A_intfc_cellI;
                    vector C_intfc_cellI = pCf[fcI];
                    const labelList& curCellsAll = cellStencil[faceOwn];

                    if(debug)
                    {
                        os<< "ph-0 cell: " << faceOwn << nl
                            << "nfOwn = " << nHatCells[faceOwn] << "  nfNei = " << nHatNeiFld[fcI] << nl
                            << "C_intfcOwn = " << C_intfcCells[faceOwn] << "  C_intfcNei = " << C_intfcNeiFld[fcI] << nl
                            << "A_intfcOwn = " << A_intfcCells[faceOwn] << "  A_intfcNei = " << A_intfcNeiFld[fcI] << nl
                            << "nf = " << nf << "  C_intfc = " << C_intfc_cellI << "  A_intfc = " << A_intfc_cellI << nl
                            << "C_ph1Own = " << C_ph1_flatFld[faceOwn] << "  C_ph0Own = " << C_ph0_flatFld[faceOwn] << nl
                            << "rho1Own = " << rho1Cells[faceOwn] << "  rho0Own = " << rho0Cells[faceOwn] << nl                            
                            << "rho1Nei = " << rho1NeiFld[fcI] << "  rho0Nei = " << rho0NeiFld[fcI] << endl;
                    }

                    scalar dn0;
                    List<scalar> Yeff0(nSpecies);                    
                    scalar intfcGradi_cellI;

                    //phase-0
                    //ensure nf direction is into the phase
                    //then reverse nf again for Js0 calculation
                    calc_cell_intfcGrad_coeffs(mesh, faceOwn, -nf, C_intfc_cellI, Y0_flatFld, alpha1_flatFld, C_ph0_flatFld, curCellsAll, nSpecies, ALPHA_2PH_MIN, 0, dn0, Yeff0, debug, os);
                    if(debug)
                    {
                        os<< "dn0 = " << dn0 << endl;
                    }
                    if(dn0 < SMALL)
                    {
                        dn0 += SMALL;
                    }
                    if(debug)
                    {
                        os<< "dn0 stab = " << dn0 << endl;
                    }
                    
                    for(label i=0; i<nSpecies; i++)
                    {
                        intfcGradi_cellI = (Yeff0[i] - Ys0[i].internalField()[faceOwn])/dn0;
                
                        JsOwn[i][bndFaceI] = A_intfc_cellI*rho0Cells[faceOwn]*D0[i].internalField()[faceOwn]*intfcGradi_cellI;
                        JsNei[i][bndFaceI] = JsOwn[i][bndFaceI];

                        if(debug)
                        {
                            os<< "species: " << i << "  Ys0iOwn = " << Ys0[i].internalField()[faceOwn] << "  Yeff0i = " << Yeff0[i] << "  intfcGrad0iOwn = " << intfcGradi_cellI << "  D0iOwn = " << D0[i].internalField()[faceOwn] << "  Js0iOwn = " <<  JsOwn[i][bndFaceI] << "  Js0iNei = " <<  JsNei[i][bndFaceI] << endl;
                        }
                    }
                }
                
                faceI++;
                bndFaceI++;
            }
        }
    }

    syncTools::swapBoundaryFaceList(mesh, phNei);
    for(label i=0; i<nSpecies; i++)
    {
        syncTools::swapBoundaryFaceList(mesh, JsNei[i]);
    }

    if(debug)
    {
        os<< "-------------------------------------------------------------------------" << nl
            << "Boundary coupled faces step 2" << nl
            << "-------------------------------------------------------------------------" << nl
            << endl;
    }

    forAll(Js1[0].boundaryField(), patchI)
    {
        const polyPatch& pp = patches[patchI];        

        if(pp.coupled())
        {
            if(debug)
            {
                os<< "-------------------------------------------------------------------------" << nl
                    << "Patch: " << patchNames[patchI] << nl
                    << "-------------------------------------------------------------------------" << nl
                    << endl;
            }

            label faceI = pp.start();
            label bndFaceI = faceI - mesh.nInternalFaces();

            const scalarField& alpha1NeiFld = alpha1.boundaryField()[patchI].patchNeighbourField();

            forAll(Js1[0].boundaryField()[patchI], fcI)
            {
                label faceOwn = own[faceI];
                scalar alpha1Own = alpha1Cells[faceOwn];
                scalar alpha1Nei = alpha1NeiFld[fcI];                

                if(alpha1Own >= ALPHA_2PH_MAX && alpha1Nei <= ALPHA_2PH_MIN)
                {
                    if(debug)
                    {
                        os<< "Face: " << faceI << "  patch face index: " << fcI << "  bnd face index: " << bndFaceI << nl
                            << "own: " << faceOwn << "  alpha1Own = " << alpha1Own << "  alpha1Nei = " << alpha1Nei << endl;
                    }

                    if(debug)
                    {
                        os<< "-------------------------------------------------------------------------" << endl;                            
                    }
                    for(label i=0; i<nSpecies; i++)
                    {
                        Js1[i].internalField()[faceOwn] = JsOwn[i][bndFaceI];
                        Js0[i].internalField()[faceOwn] = JsNei[i][bndFaceI];
                        
                        if(debug)
                        {
                            os<< "species: " << i << "  JsOwn = " << JsOwn[i][bndFaceI] << "  JsNei = " << JsNei[i][bndFaceI] << "  Js1 = " << Js1[i].internalField()[faceOwn] << "  Js0 = " << Js0[i].internalField()[faceOwn] << endl;
                        }
                    }
                    if(debug)
                    {
                        os<< "-------------------------------------------------------------------------" << nl
                            << endl;
                    }
                }

                if(alpha1Own <= ALPHA_2PH_MIN && alpha1Nei >= ALPHA_2PH_MAX)
                {
                    if(debug)
                    {
                        os<< "Face: " << faceI << "  patch face index: " << fcI << "  bnd face index: " << bndFaceI << nl
                            << "own: " << faceOwn << "  alpha1Own = " << alpha1Own << "  alpha1Nei = " << alpha1Nei << endl;
                    }

                    if(debug)
                    {
                        os<< "-------------------------------------------------------------------------" << endl;                            
                    }
                    for(label i=0; i<nSpecies; i++)
                    {
                        Js0[i].internalField()[faceOwn] = JsOwn[i][bndFaceI];
                        Js1[i].internalField()[faceOwn] = JsNei[i][bndFaceI];

                        if(debug)
                        {
                            os<< "species: " << i << "  JsOwn = " << JsOwn[i][bndFaceI] << "  JsNei = " << JsNei[i][bndFaceI] << "  Js1 = " << Js1[i].internalField()[faceOwn] << "  Js0 = " << Js0[i].internalField()[faceOwn] << endl;
                        }
                    }
                    if(debug)
                    {
                        os<< "-------------------------------------------------------------------------" << nl
                            << endl;
                    }
                }                

                faceI++;
                bndFaceI++;
            }
        }
    }

    if(debug)
    {
        os<< nl
            << " Done Interfacial Species Flux Calculation" << nl
            << "-------------------------------------------------------------------------" << nl
            << endl;
    }
}


void x2y(int n, double *MW, double *x, double *y)
{
    int i;
    double sum_MWx = 0.0, r1_sum;
    for (i=0; i<n; i++) sum_MWx += MW[i]*x[i];

    if(sum_MWx < SMALL) sum_MWx += SMALL;

    r1_sum = 1.0/sum_MWx;

    for (i=0; i<n; i++) y[i] = MW[i]*x[i]*r1_sum;
}


void y2x(int n, double *MW, double *y, double *x)
{
    int i;
    double sum_MWy = 0.0, r1_sum;
    for (i=0; i<n; i++) sum_MWy += y[i]/MW[i];

    if(sum_MWy < SMALL) sum_MWy += SMALL;

    r1_sum = 1.0/sum_MWy;

    for (i=0; i<n; i++) x[i] = y[i]*r1_sum/MW[i];
}


double calc_MW_from_x
(
    int n,
    double *MW,
    double *x
)
{
    int i;
    double MWm = 0;
    
    for(i=0; i<n; i++) MWm += MW[i]*x[i];

    return MWm;
}


double calc_rho_from_V(int n, double *x, double *MW, double V)
{
    int i;
    double sum = 0.0;

    for (i=0;i<n;i++) 
    {                             
        sum += MW[i]*x[i];
    }

    if(V < SMALL) V += SMALL;

    return 1e-3*sum/V;
}


double calc_CvIG_from_CpIG(int n, double *x, double *CpIG)
{
    int i;
    double sum = 0.0;

    for (i=0;i<n;i++)
    {
        sum += CpIG[i]*x[i];
    }

    return sum - 8.3144621;
}


void calc_kij_from_table
(
    double T,
    int n,
    double Ta,
    double Tb,
    int nT,
    double *kij_T,
    double *kij
)
{
    int i, j, idT, idx;
    double dT, T1;

    dT = (Tb - Ta)/(nT - 1);

    if(T > Ta && T < Tb)
    {
        idT = floor((T - Ta)/dT);
    }
    else if(T <= Ta)
    {
        idT = 0;
    }
    else
    {
        idT = nT - 1;
    }

    T1 = Ta + idT*dT;

    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
        {
            idx = i*n + j;
            kij[idx] = kij_T[idT*n*n + idx] + (kij_T[(idT+1)*n*n + idx] - kij_T[idT*n*n + idx])*(T - T1)/dT;
        }
    }
}


void calc_rho
(
    double P,    
    double T,
    double *x,
    int n,
    double *Pc,
    double *Tc,
    double *w,
    double *MW,    
    double *kij,    
    double& rho
)
{
    double rho_tmp;

    density_pr_eos2_(&P,&T,x,&n,Pc,Tc,w,MW,kij,&rho_tmp);    
    rho = rho_tmp;
}


void calc_mu
(
    double P,
    double T,
    double *y,
    int n,
    double *Pc,
    double *Tc,
    double *Vc,
    double *w,
    double *MW,    
    double *Tb,
    double *SG,
    double *H8,    
    double *k,
    double *dm,    
    double *kij,    
    double& mu
)
{
    int i, j;
    double vis, cond, V, CvIG;
    double *x_tmp;

    _NNEW_(x_tmp, double, n);

    mu = 0;

    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
        {
            if(j==i) x_tmp[j] = 1;
            else x_tmp[j] = 0;
        }
        calc_v_cvig_(&P,&T,x_tmp,&n,Pc,Tc,w,MW,Tb,SG,H8,kij,&V,&CvIG);
        vis_n_cond_(&P,&T,&n,Pc,Tc,Vc,w,MW,k,dm,x_tmp,&CvIG,&V,&cond,&vis);
        mu += y[i]*vis;
    }

    _DDELETE_(x_tmp);
}


void calc_mu_lambda
(
    double P,
    double T,
    double *y,
    int n,
    double *Pc,
    double *Tc,
    double *Vc,
    double *w,
    double *MW,    
    double *Tb,
    double *SG,
    double *H8,    
    double *k,
    double *dm,
    double *kij,    
    double& mu,
    double& lambda
)
{
    int i, j;
    double vis, cond, V, CvIG;
    double *x_tmp;

    _NNEW_(x_tmp, double, n);

    mu = 0;
    lambda = 0;

    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
        {
            if(j==i) x_tmp[j] = 1;
            else x_tmp[j] = 0;
        }
        calc_v_cvig_(&P,&T,x_tmp,&n,Pc,Tc,w,MW,Tb,SG,H8,kij,&V,&CvIG);
        vis_n_cond_(&P,&T,&n,Pc,Tc,Vc,w,MW,k,dm,x_tmp,&CvIG,&V,&cond,&vis);
        mu += y[i]*vis;
        lambda += y[i]*cond;
    }

    _DDELETE_(x_tmp);
}


void correct_rho
(
    double P,    
    const volScalarField& T,
    const PtrList<volScalarField>& X,
    int n,
    double *Pc,
    double *Tc,
    double *w,
    double *MW,
    double *kij_T,
    double Ta_kij,
    double Tb_kij,
    int nT_kij,
    bool calcKij,
    double *kij,
    volScalarField& rho
)
{
    int i;
    double T_tmp, rho_tmp;
    double *x_tmp;

    _NNEW_(x_tmp, double, n);

    const scalarField& TCells = T.internalField();
    scalarField& rhoCells = rho.internalField();

    forAll(TCells, cellI)
    {
        T_tmp = TCells[cellI];
        for(i=0; i<n; i++)
        {
            x_tmp[i] = X[i].internalField()[cellI];
        }
        
        if(calcKij) calc_kij_from_table(T_tmp, n, Ta_kij, Tb_kij, nT_kij, kij_T, kij);

        calc_rho(P,T_tmp,x_tmp,n,Pc,Tc,w,MW,kij,rho_tmp);
     
        rhoCells[cellI] = rho_tmp;
    }

    forAll(T.boundaryField(), patchI)
    {
        const fvPatchScalarField& pT = T.boundaryField()[patchI];
        fvPatchScalarField& prho = rho.boundaryField()[patchI];
        
        forAll(pT, fcI)
        {
            T_tmp = pT[fcI];
            for(i=0; i<n; i++)
            {
                x_tmp[i] = X[i].boundaryField()[patchI][fcI];
            }
            
            if(calcKij) calc_kij_from_table(T_tmp, n, Ta_kij, Tb_kij, nT_kij, kij_T, kij);

            calc_rho(P,T_tmp,x_tmp,n,Pc,Tc,w,MW,kij,rho_tmp);

            prho[fcI] = rho_tmp;
        }
    }

    _DDELETE_(x_tmp);
}


void correct_h
(
    double P,
    const volScalarField& T,
    const PtrList<volScalarField>& X,
    int n,
    double *Pc,
    double *Tc,
    double *w,
    double *MW,
    double *Tb,
    double *SG,
    double *H8,
    double *kij_T,
    double Ta_kij,
    double Tb_kij,
    int nT_kij,
    double *kij,
    volScalarField& h
)
{
    int i;
    double T_tmp, h_tmp, v_tmp, MW_tmp;
    double *x_tmp;

    _NNEW_(x_tmp, double, n);

    const scalarField& TCells = T.internalField();
    scalarField& hCells = h.internalField();

    forAll(TCells, cellI)
    {
        T_tmp = TCells[cellI];
        for(i=0; i<n; i++)
        {
            x_tmp[i] = X[i].internalField()[cellI];
        }
        MW_tmp = 0;
        for(i=0; i<n; i++)
        {
            MW_tmp += x_tmp[i]*MW[i];
        }
        MW_tmp *= 1e-3;

        calc_kij_from_table(T_tmp,n,Ta_kij,Tb_kij,nT_kij,kij_T,kij);
        calc_v_h_(&P, &T_tmp, x_tmp, &n, Pc, Tc, w, MW, Tb, SG, H8, kij, &v_tmp, &h_tmp);
        hCells[cellI] = h_tmp/MW_tmp;
    }

    forAll(T.boundaryField(), patchI)
    {
        const fvPatchScalarField& pT = T.boundaryField()[patchI];
        fvPatchScalarField& ph = h.boundaryField()[patchI];
        
        forAll(pT, fcI)
        {
            T_tmp = pT[fcI];
            for(i=0; i<n; i++)
            {
                x_tmp[i] = X[i].boundaryField()[patchI][fcI];
            }
            MW_tmp = 0;
            for(i=0; i<n; i++)
            {
                MW_tmp += x_tmp[i]*MW[i];
            }
            MW_tmp *= 1e-3;

            calc_kij_from_table(T_tmp,n,Ta_kij,Tb_kij,nT_kij,kij_T,kij);
            calc_v_h_(&P, &T_tmp, x_tmp, &n, Pc, Tc, w, MW, Tb, SG, H8, kij, &v_tmp, &h_tmp);
            ph[fcI] = h_tmp/MW_tmp;
        }
    }

    _DDELETE_(x_tmp);
}


void correct_boundaryField_C
(    
    volScalarField& Yi,
    const volScalarField& rho,
    const volScalarField& alpha,
    volScalarField& ci,
    volScalarField& Ci
)
{
    forAll(Yi.boundaryField(), patchI)
    {
        fvPatchScalarField& pYi = Yi.boundaryField()[patchI];
        const fvPatchScalarField& prho = rho.boundaryField()[patchI];
        const fvPatchScalarField& palpha = alpha.boundaryField()[patchI];
        fvPatchScalarField& pci = ci.boundaryField()[patchI];
        fvPatchScalarField& pCi = Ci.boundaryField()[patchI];

        forAll(pYi, fcI)
        {
            pci[fcI] = prho[fcI]*pYi[fcI];
            pCi[fcI] = palpha[fcI]*pci[fcI];
        }
    }
}


void correct_x_from_Y
(
    int n,
    double *MW,
    const PtrList<volScalarField>& Y,
    PtrList<volScalarField>& x
)
{
    int i;
    double *y_tmp, *x_tmp;

    _NNEW_(y_tmp, double, n);
    _NNEW_(x_tmp, double, n);

    forAll(Y[0].internalField(), cellI)
    {
        for(i=0; i<n; i++)
        {
            y_tmp[i] = Y[i].internalField()[cellI];
        }

        plicFuncs::y2x(n, MW, y_tmp, x_tmp);

        for(i=0; i<n; i++)
        {
            x[i].internalField()[cellI] = x_tmp[i];
        }
    }

    forAll(Y[0].boundaryField(), patchI)
    {
        forAll(Y[0].boundaryField()[patchI], fcI)
        {
            for(i=0; i<n; i++)
            {
                y_tmp[i] = Y[i].boundaryField()[patchI][fcI];
            }

            plicFuncs::y2x(n, MW, y_tmp, x_tmp);

            for(i=0; i<n; i++)
            {
                x[i].boundaryField()[patchI][fcI] = x_tmp[i];
            }
        }
    }

    _DDELETE_(y_tmp);
    _DDELETE_(x_tmp);
}


void correct_boundaryField_h_rhoh_H
(
    double P,
    const volScalarField& T,
    const PtrList<volScalarField>& X,
    const volScalarField& rho,
    const volScalarField& alpha,
    int n,
    double *Pc,
    double *Tc,
    double *w,
    double *MW,
    double *Tb,
    double *SG,
    double *H8,
    double *kij_T,
    double Ta_kij,
    double Tb_kij,
    int nT_kij,
    double *kij,    
    volScalarField& h,
    volScalarField& rhoh,
    volScalarField& H
)
{
    int i;
    double T_tmp, h_tmp, v_tmp, MW_tmp;
    double *x_tmp;

    _NNEW_(x_tmp, double, n);

    forAll(T.boundaryField(), patchI)
    {
        const fvPatchScalarField& pT = T.boundaryField()[patchI];
        const fvPatchScalarField& prho = rho.boundaryField()[patchI];
        const fvPatchScalarField& palpha = alpha.boundaryField()[patchI];
        fvPatchScalarField& ph = h.boundaryField()[patchI];
        fvPatchScalarField& prhoh = rhoh.boundaryField()[patchI];
        fvPatchScalarField& pH = H.boundaryField()[patchI];
        
        forAll(pT, fcI)
        {
            T_tmp = pT[fcI];
            MW_tmp = 0;
            for(i=0; i<n; i++)
            {
                x_tmp[i] = X[i].boundaryField()[patchI][fcI];
                MW_tmp += x_tmp[i]*MW[i];
            }

            calc_kij_from_table(T_tmp,n,Ta_kij,Tb_kij,nT_kij,kij_T,kij);
            calc_v_h_(&P, &T_tmp, x_tmp, &n, Pc, Tc, w, MW, Tb, SG, H8, kij, &v_tmp, &h_tmp);
            h_tmp /= MW_tmp;

            ph[fcI] = h_tmp;
            prhoh[fcI] = prho[fcI]*h_tmp;
            pH[fcI] = palpha[fcI]*prho[fcI];
        }
    }

    _DDELETE_(x_tmp);
}

void correct_thermo_trans_prop
(
    const fvMesh& mesh,
    double P,
    const volScalarField& T,
    const PtrList<volScalarField>& X,
    const PtrList<volScalarField>& Y,
    int n,
    double *Pc,
    double *Tc,
    double *Vc,
    double *w,
    double *MW,    
    double *Tb,
    double *SG,
    double *H8,
    double *k,
    double *dm,
    double *kij_T,
    double Ta_kij,
    double Tb_kij,
    int nT_kij,
    double *kij,
    volScalarField& v,
    PtrList<volScalarField>& hpar,
    volScalarField& lambda,
    volScalarField& Cp,
    volScalarField& mu,
    PtrList<volScalarField>& D
)
{
    int i, j, idx;
    double T_tmp, V, MW_tmp, CvIG, Cp_tmp, cond, vis;     
    double *x_tmp, *y_tmp, *Dij, *h_tmp;

    _NNEW_(x_tmp, double, n);
    _NNEW_(y_tmp, double, n);
    _NNEW_(Dij, double, n*n);
    _NNEW_(h_tmp, double, n);

    const scalarField& TCells = T.internalField();
    scalarField& vCells = v.internalField();
    scalarField& lambdaCells = lambda.internalField();
    scalarField& CpCells = Cp.internalField();
    scalarField& muCells = mu.internalField();

    forAll(TCells, cellI)
    {
        T_tmp = TCells[cellI];
        for(i=0; i<n; i++)
        {
            x_tmp[i] = X[i].internalField()[cellI];
            y_tmp[i] = Y[i].internalField()[cellI];
        }
        MW_tmp = 0;
        for(i=0; i<n; i++)
        {
            MW_tmp += x_tmp[i]*MW[i];
        }
        MW_tmp *= 1e-3;
       
        calc_kij_from_table(T_tmp,n,Ta_kij,Tb_kij,nT_kij,kij_T,kij);

        calc_v_cvig_cp_hpar_(&P,&T_tmp,x_tmp,&n,Pc,Tc,w,MW,Tb,SG,H8,kij,&V,&CvIG,&Cp_tmp,h_tmp);

        calc_mu_lambda(P,T_tmp,y_tmp,n,Pc,Tc,Vc,w,MW,Tb,SG,H8,k,dm,kij,vis,cond);

        new_tlsm_diffusion_krishna_model_(&P,&T_tmp,&n,Pc,Tc,Vc,w,MW,kij,x_tmp,Dij);

        vCells[cellI] = V;
        lambdaCells[cellI] = cond;
        CpCells[cellI] = Cp_tmp/MW_tmp;
        muCells[cellI] = vis;
        for(i=0; i<n; i++)
        {
            hpar[i].internalField()[cellI] = h_tmp[i]/(MW[i]*1e-3);
            
            for(j=0; j<n; j++)
            {
                idx = i*n + j;
                D[idx].internalField()[cellI] = Dij[idx];
            }
        }
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(T.boundaryField(), patchI)
    {
        const polyPatch& pp = patches[patchI];
        const fvPatchScalarField& pT = T.boundaryField()[patchI];
        fvPatchScalarField& pv = v.boundaryField()[patchI];
        fvPatchScalarField& plambda = lambda.boundaryField()[patchI];
        fvPatchScalarField& pCp = Cp.boundaryField()[patchI];
        fvPatchScalarField& pmu = mu.boundaryField()[patchI];

        if((!pp.coupled()))
        {
            forAll(pT, fcI)
            {
                T_tmp = pT[fcI];

                for(i=0; i<n; i++)
                {
                    x_tmp[i] = X[i].boundaryField()[patchI][fcI];
                    y_tmp[i] = Y[i].boundaryField()[patchI][fcI];
                }
                MW_tmp = 0;
                for(i=0; i<n; i++)
                {
                    MW_tmp += x_tmp[i]*MW[i];
                }
                MW_tmp *= 1e-3;
        
                calc_kij_from_table(T_tmp,n,Ta_kij,Tb_kij,nT_kij,kij_T,kij);

                calc_v_cvig_cp_hpar_(&P,&T_tmp,x_tmp,&n,Pc,Tc,w,MW,Tb,SG,H8,kij,&V,&CvIG,&Cp_tmp,h_tmp);

                calc_mu_lambda(P,T_tmp,y_tmp,n,Pc,Tc,Vc,w,MW,Tb,SG,H8,k,dm,kij,vis,cond);

                new_tlsm_diffusion_krishna_model_(&P,&T_tmp,&n,Pc,Tc,Vc,w,MW,kij,x_tmp,Dij);

                pv[fcI] = V;
                plambda[fcI] = cond;
                pCp[fcI] = Cp_tmp/MW_tmp;
                pmu[fcI] = vis;
                for(i=0; i<n; i++)
                {
                    hpar[i].boundaryField()[patchI][fcI] = h_tmp[i]/(MW[i]*1e-3);
            
                    for(j=0; j<n; j++)
                    {
                        idx = i*n + j;
                        D[idx].boundaryField()[patchI][fcI] = Dij[idx];
                    }
                }
            }
        }
    }

    _DDELETE_(x_tmp);
    _DDELETE_(y_tmp);
    _DDELETE_(Dij);
    _DDELETE_(h_tmp);
}


void correct_thermo_trans_prop
(
    const fvMesh& mesh,
    double P,
    const volScalarField& T,
    const PtrList<volScalarField>& X,
    const PtrList<volScalarField>& Y,    
    int n,
    double *Pc,
    double *Tc,
    double *Vc,
    double *w,
    double *MW,    
    double *Tb,
    double *SG,
    double *H8,
    double *k,
    double *dm,
    double *kij,    
    volScalarField& v,
    volScalarField& mu,
    PtrList<volScalarField>& D,
    const bool debug,
    OFstream& os
)
{
    int i, j, idx;
    double T_tmp, V, MW_tmp; 
    double CvIG;    
    double cond, vis;
    double *x_tmp, *y_tmp, *Dij;

    _NNEW_(x_tmp, double, n);
    _NNEW_(y_tmp, double, n);
    _NNEW_(Dij, double, n*n);

    const scalarField& TCells = T.internalField();
    scalarField& vCells = v.internalField();    
    scalarField& muCells = mu.internalField();

    forAll(TCells, cellI)
    {
        T_tmp = TCells[cellI];

        for(i=0; i<n; i++)
        {
            x_tmp[i] = X[i].internalField()[cellI];
            y_tmp[i] = Y[i].internalField()[cellI];            
        }
        
        MW_tmp = 0;
        for(i=0; i<n; i++)
        {
            MW_tmp += x_tmp[i]*MW[i];
        }
        MW_tmp *= 1e-3;                

        calc_v_cvig_(&P,&T_tmp,x_tmp,&n,Pc,Tc,w,MW,Tb,SG,H8,kij,&V,&CvIG);

        calc_mu_lambda(P,T_tmp,y_tmp,n,Pc,Tc,Vc,w,MW,Tb,SG,H8,k,dm,kij,vis,cond);

        new_tlsm_diffusion_krishna_model_(&P,&T_tmp,&n,Pc,Tc,Vc,w,MW,kij,x_tmp,Dij);

        vCells[cellI] = V;
        muCells[cellI] = vis;
        for(i=0; i<n; i++)
        {            
            for(j=0; j<n; j++)
            {
                idx = i*n + j;
                D[idx].internalField()[cellI] = Dij[idx];
            }
        }
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    //const wordList& patchNames = patches.names();
    
    forAll(T.boundaryField(), patchI)
    {
        const polyPatch& pp = patches[patchI];
        const fvPatchScalarField& pT = T.boundaryField()[patchI];
        fvPatchScalarField& pv = v.boundaryField()[patchI];        
        fvPatchScalarField& pmu = mu.boundaryField()[patchI];

        if(!(pp.coupled()))
        {
            forAll(pT, fcI)
            {
                T_tmp = pT[fcI];
                
                for(i=0; i<n; i++)
                {
                    x_tmp[i] = X[i].boundaryField()[patchI][fcI];
                    y_tmp[i] = Y[i].boundaryField()[patchI][fcI];
                }

                MW_tmp = 0;
                for(i=0; i<n; i++)
                {
                    MW_tmp += x_tmp[i]*MW[i];
                }
                MW_tmp *= 1e-3;

                calc_v_cvig_(&P,&T_tmp,x_tmp,&n,Pc,Tc,w,MW,Tb,SG,H8,kij,&V,&CvIG);

                calc_mu_lambda(P,T_tmp,y_tmp,n,Pc,Tc,Vc,w,MW,Tb,SG,H8,k,dm,kij,vis,cond);

                new_tlsm_diffusion_krishna_model_(&P,&T_tmp,&n,Pc,Tc,Vc,w,MW,kij,x_tmp,Dij);

                pv[fcI] = V;            
                pmu[fcI] = vis;
                for(i=0; i<n; i++)
                {
                    for(j=0; j<n; j++)
                    {
                        idx = i*n + j;
                        D[idx].boundaryField()[patchI][fcI] = Dij[idx];
                    }
                }
            }
        }
    }

    _DDELETE_(x_tmp);
    _DDELETE_(y_tmp);
    _DDELETE_(Dij);    
}


void calc_D_from_Dij
(
    int n,
    double *x,
    double *Dij,
    double *D
)
{
    int i, j, idx;
    double sum_xbyD;

    if(n==2)
    {
        for(i=0; i<n; i++)
        {
            D[i] = Dij[1];
        }
    }
    else
    {
        for(i=0; i<n; i++)
        {
            sum_xbyD = 0.0;
            for(j=0; j<n; j++)
            {
                if(j!=i)
                {
                    idx = j + i*n;
                    sum_xbyD += x[j]/Dij[idx];
                }
            }
            if(sum_xbyD < SMALL) 
            {
                D[i] = 0.0;
                for(j=0; j<n; j++)
                {
                    if(j!=i)
                    {
                        idx = j + i*n;
                        D[i] += Dij[idx];
                    }
                }
                D[i] /= (n-1);
            }
            else
            {
                D[i] = 1.0/sum_xbyD;
            }
        }
    }
}


void correct_D_from_Dij
(
    int n,
    const PtrList<volScalarField>& x,
    const PtrList<volScalarField>& Dij,
    PtrList<volScalarField>& D
)
{
    int i, j, idx;
    double *x_tmp;
    double *Dij_tmp; 
    double *D_tmp;    

    _NNEW_(x_tmp, double, n);
    _NNEW_(Dij_tmp, double, n*n);
    _NNEW_(D_tmp, double, n);

    forAll(x[0].internalField(), cellI)
    {
        for(i=0; i<n; i++)
        {
            x_tmp[i] = x[i].internalField()[cellI];
            for(j=0; j<n; j++)
            {
                idx = j + i*n;
                Dij_tmp[idx] = Dij[idx].internalField()[cellI];
            }
        }

        calc_D_from_Dij(n, x_tmp, Dij_tmp, D_tmp);
        
        for(i=0; i<n; i++)
        {
            D[i].internalField()[cellI] = D_tmp[i];
        }
    }

    forAll(x[0].boundaryField(), patchI)
    {
        const fvPatchScalarField& px0 = x[0].boundaryField()[patchI];

        forAll(px0, fcI)
        {
            for(i=0; i<n; i++)
            {
                x_tmp[i] = x[i].boundaryField()[patchI][fcI];
                for(j=0; j<n; j++)
                {
                    idx = j + i*n;
                    Dij_tmp[idx] = Dij[idx].boundaryField()[patchI][fcI];
                }
            }

            calc_D_from_Dij(n, x_tmp, Dij_tmp, D_tmp);            
        
            for(i=0; i<n; i++)
            {
                D[i].boundaryField()[patchI][fcI] = D_tmp[i];
            }
        }
    }

    _DDELETE_(x_tmp);
    _DDELETE_(Dij_tmp);
    _DDELETE_(D_tmp);
}

void calc_intfc_transLLE
(
    double P,
    double& Ts,
    int n,
    double *Pc,
    double *Tc,
    double* Vc,
    double *w,
    double *MW,    
    double *Tb,
    double *SG,
    double *H8,
    double *k,
    double *dm,
    double *kij_T,
    double Ta_kij,
    double Tb_kij,
    int nT_kij,
    double *kij,
    double dn1,
    double dn0,
    double *xeff1,
    double *xeff0,
    double Teff1,
    double Teff0,
    double *xs1,
    double *xs0,
    double *ys1,
    double *ys0,
    double& JsTot,
    double *flux_m_1,
    double *flux_m_0,
    double& cond1,
    double *h1,
    int& iLLE,
    int& iTs,
    int n_flux_type,
    LPT_UMFPACK flux_umf,
    double Ts_TOL,
    int MAX_ITERS_Ts,
    double MASS_FRAC_TOL,
    const bool debug,
    OFstream& os
)
{
    int i, num;    
    double Ts_tmp;
    double cond0, vis;
    double J1, J0, qs, TNum, TDen, Ts_err;
    double *lnphi1; double *lnphi0; double *Dij1; double *Dij0; 
    double *h0;

    _NNEW_(lnphi1, double, n);
    _NNEW_(lnphi0, double, n);    
    _NNEW_(Dij1, double, n*n);
    _NNEW_(Dij0, double, n*n);    
    _NNEW_(h0, double, n);    

    Ts_tmp = Ts;

    if(debug)
    {
        print_line(os, 100);
        os<< "Interface Tansport-LLE Calculation" << endl;
        print_line(os, 100);
        os<< endl;
    }

    for(iTs=0; iTs<MAX_ITERS_Ts; iTs++)
    {
        if(debug)
        {
            print_line(os, 100);
            os<< "iTs = " << iTs+1 << "  Ts = " << Ts_tmp << endl;
            print_line(os, 100);
            os<< endl;
        }
        //calculate BIPs for current iteration interface temperature Ts_tmp
        calc_kij_from_table(Ts_tmp,n,Ta_kij,Tb_kij,nT_kij,kij_T,kij);

        if(debug)
        {            
            os<< "BIP matrix:  ";
            for(i=0; i<n*n; i++) os<< kij[i] << "  ";
            os<< nl<< endl;
            print_line(os, 100);
            os<< "GSL optimization function" << endl;
            print_line(os, 100);
            os<< "dn1 = " << dn1 << "  dn0 = " << dn0 << endl;
            print_line(os, 100);
            os<< setw(8) << "Species" << setw(14) << "xeff1" << setw(14) << "xeff0" << setw(14) << "xs1" << setw(14) << "xs0" << endl;
            print_line(os, 100);
            for(i=0; i<n; i++)
            {
                os<< setw(8) << i << setw(14) << xeff1[i] << setw(14) << xeff0[i] << setw(14) << xs1[i] << setw(14) << xs0[i] << endl;
            }
            print_line(os, 100);
        }

        iLLE = my2_gsl_find_transport_LLE(n,Pc,Tc,Vc,w,MW,Tb,SG,H8,kij,n_flux_type,dn1,dn0,P,Ts_tmp,xeff1,xeff0,xs1,xs0,lnphi1,lnphi0,Dij1,Dij0,h1,h0,flux_m_1,flux_m_0,flux_umf);

        if(debug)
        {
            print_line(os, 100);
            os<< "Done GSL optimization function" << endl;
            os<< "iLLE = " << iLLE << endl;
            print_line(os, 100);
            os<< setw(8) << "Species" << setw(14) << "xs1" << setw(14) << "xs0" << setw(14) << "flux_m_1" << setw(14) << "flux_m_0" << endl;
            print_line(os, 100);
            for(i=0; i<n; i++)
            {
                os<< setw(8) << i << setw(14) << xs1[i] << setw(14) << xs0[i] << setw(14) << flux_m_1[i] << setw(14) << flux_m_0[i] << endl;
            }
            print_line(os, 100);
            os<< nl << endl;
        }

        //calculate fluid properties on both sides of interface
        //convert the mole fractions to mass fractions
        x2y(n,MW,xs1,ys1);
        x2y(n,MW,xs0,ys0);

        //phase 1 
        //thermal conductivity
        calc_mu_lambda(P,Ts_tmp,ys1,n,Pc,Tc,Vc,w,MW,Tb,SG,H8,k,dm,kij,cond1,vis);
        //partial molar enthalpies
        calc_hpar_(&P,&Ts_tmp,xs1,&n,Pc,Tc,w,MW,Tb,SG,H8,kij,h1);

        //phase 0 
        //thermal conductivity
        calc_mu_lambda(P,Ts_tmp,ys0,n,Pc,Tc,Vc,w,MW,Tb,SG,H8,k,dm,kij,cond0,vis);
        //partial molar enthalpies
        calc_hpar_(&P,&Ts_tmp,xs0,&n,Pc,Tc,w,MW,Tb,SG,H8,kij,h0);

        //calculate the total ph-1 to ph-0 phase change mass flux 
        //from ph-1 and ph-0 diffusive fluxes
        JsTot = 0;
        num = 0;
        for(i=0; i<n; i++)
        {
            if(mag(ys1[i]-ys0[i]) > MASS_FRAC_TOL)
            {
                JsTot += (flux_m_0[i] - flux_m_1[i])/(ys1[i] - ys0[i]);
                num++;
            }
        }

        if(num > 0)
        {
            JsTot = JsTot/(static_cast<double>(num));
        }

        //calculate the enthalpy flux due to species fluxes
        qs=0;
        for(i=0; i<n; i++)
        {
            J1 = JsTot*ys1[i] + flux_m_1[i];
            J0 = JsTot*ys0[i] + flux_m_0[i];
            qs += (J0*h0[i] - J1*h1[i])/(MW[i]*1e-3);
        }

        //calculate the interface temperature consistent with
        //energy balance across the interface
        TNum = cond1*Teff1/dn1 + cond0*Teff0/dn0 - qs;
        TDen = cond1/dn1 + cond0/dn0;
        if(TDen < SMALL) TDen += SMALL;
        Ts_tmp = TNum/TDen;

        //calculate error in Ts and break if converged
        Ts_err = mag(Ts_tmp - Ts);
        if(Ts_err < Ts_TOL) break;

        //Ts under-relaxation
        if (n_flux_type>0) 
        {
            Ts_tmp = Ts_tmp*0.5 + Ts*0.5;	
        }
        else
        {
            if (iTs<3)
                Ts_tmp = Ts_tmp*0.5 + Ts*0.5;
            else if (iTs<20)
                Ts_tmp = Ts_tmp*0.1 + Ts*0.9;
            else
                Ts_tmp = Ts_tmp*0.01+ Ts*0.99;
        }

        //update prev iteration value of Ts
        Ts = Ts_tmp;
    }    

    _DDELETE_(lnphi1);
    _DDELETE_(lnphi0);
    _DDELETE_(Dij1);
    _DDELETE_(Dij0);    
    _DDELETE_(h0);
}


void calc_limiter_mS1
(
    double& mS1Tot_cellI, 
    double *mS1_cellI, 
    double *C1_cellI, 
    double *C0_cellI,
    double& alpha1_cellI,
    double& rho1_cellI,
    double& rho0_cellI,
    double& V_cellI, 
    double& dt,
    int n,
    double& limiter_mS1Tot, 
    double *limiter_mS1, 
    double& limiter_min, 
    const bool debug, 
    OFstream& os
)
{
    int i;
    double mS1Tot_tmp, max_mS, mS1i_cellI, max_mSi;

    limiter_min = 1.0;
    mS1Tot_tmp = mS1Tot_cellI/V_cellI;

    if(mS1Tot_tmp > 0)
    {
        max_mS = rho0_cellI*(1 - alpha1_cellI)/dt;
        limiter_mS1Tot = min(1.0, max_mS/mS1Tot_tmp);
    }
    else
    {
        max_mS = rho1_cellI*alpha1_cellI/dt;
        limiter_mS1Tot = min(1.0, max_mS/mag(mS1Tot_tmp));
    }
    limiter_min = min(limiter_min, limiter_mS1Tot);

    for(i=0; i<n; i++)
    {
        mS1i_cellI = mS1_cellI[i]/V_cellI;
        if(mS1i_cellI > 0)
        {
            max_mSi = C0_cellI[i]/dt;
            limiter_mS1[i] = min(1.0, max_mSi/mS1i_cellI);            
        }
        else if(mS1i_cellI < 0)
        {
            max_mSi = C1_cellI[i]/dt;
            limiter_mS1[i] = min(1.0, max_mSi/mag(mS1i_cellI));
        }
        else
        {
            limiter_mS1[i] = 1.0;
        }
        limiter_min = min(limiter_mS1[i], limiter_min);
    }
}


void print_line(OFstream& os, int n)
{
    for(int i=0; i<n; i++)
    {
        os<< "-";
    }
    os<< endl;
}


void calc_cell_intfcGrad_coeffs
(
    const fvMesh& mesh,
    const label& curCell_lbl,
    const vector& nf,
    const vector& C_intfc,
    const List<List<scalar> >& Y,
    const List<scalar>& T,
    const List<scalar>& alpha1,
    const List<vector>& C,
    const labelList& curCellsAll,
    const label& nSpecies,
    const scalar& MIN_ALPHA_2PH,
    const label& phaseLbl,
    double& dn,
    double *Yeff,
    double& Teff,
    bool debug, 
    OFstream& os
)
{
    if(debug)
    {
        os<< nl
            << "Calculating cell intfc grad weights in cell " << curCell_lbl << nl
            << endl;
    }

    scalar MAX_ALPHA_2PH = 1 - MIN_ALPHA_2PH;
    
    labelList curCells(curCellsAll.size());
    label n_ph = 0;

    if(debug)
    {
        os<< "Reducing cell stencil for phase " << phaseLbl << nl
            << "Full cell stencil" << nl
            << curCellsAll << endl;
    }

    if(phaseLbl == 1)
    {
        for(label cellI=0; cellI<curCellsAll.size(); cellI++)
        {            
            label cellI_lbl = curCellsAll[cellI];

            if(debug)
            {
                os<< "cellI: " << cellI << "  cellI_lbl: " << cellI_lbl << nl
                    << "alpha1 = " << alpha1[cellI_lbl] << "  MIN_ALPHA_2PH = " << MIN_ALPHA_2PH
                    << endl;
            }

            if(alpha1[cellI_lbl] > MIN_ALPHA_2PH && cellI_lbl != curCell_lbl)
            {
                curCells[n_ph++] = cellI_lbl;
            }            
        }        
    }
    else
    {
        for(label cellI=0; cellI<curCellsAll.size(); cellI++)
        {
            label cellI_lbl = curCellsAll[cellI];

            if(debug)
            {
                os<< "cellI: " << cellI << "  cellI_lbl: " << cellI_lbl << nl
                    << "alpha1 = " << alpha1[cellI_lbl] << "  MAX_ALPHA_2PH = " << MAX_ALPHA_2PH
                    << endl;
            }

            if(alpha1[cellI_lbl] < MAX_ALPHA_2PH && cellI_lbl != curCell_lbl)
            {
                curCells[n_ph++] = curCellsAll[cellI];
            }            
        }
    }
    
    curCells.setSize(n_ph);

    if(debug)
    {
        os<< "Cell reduced stencil" << nl
            << curCells << nl
            << endl;
    }

    // suffix 1: direction closest to nf
    // suffix 2: direction closest orthogonal to nf in 2-D
    // further improvements needed for 3-D calculation
    bool foundCell1 = false;
    bool foundCell2 = false;
    vector Cp = C_intfc;
    label C1_lbl = findCellInIntfcDir(mesh,alpha1,curCells,C,Cp,nf,curCell_lbl,foundCell1,debug,os);
    if(debug)
    {
        os<< "C1_lbl: " << C1_lbl << "  found C1: " << foundCell1
            << endl;
    }
    if(foundCell1)
    {
        vector C1 = C[C1_lbl];
        vector t1 = C1 - Cp;
        scalar magt1;    
        magt1 = mag(t1);
        scalar magnf;
        magnf = mag(nf);
        if(magt1 < SMALL)
        {
            os<< "Time = " << mesh.time().timeName() << nl 
                << "Cell " << curCell_lbl << "  Cp = " << Cp << "  C1 = " << C1 << "  magt1 = " << magt1 << nl
                << endl;
            magt1 += SMALL;
        }
        if(magnf < SMALL)
        {
            os<< "Time = " << mesh.time().timeName() << nl 
                << "Cell " << curCell_lbl << "  nf = " << nf << "  mag(nf) = " << magnf << nl
                << endl;
            magnf += SMALL;
        }
        scalar costheta1 = (nf & t1)/magt1/magnf;
        if(costheta1 > 1)
        {        
            os<< "Time = " << mesh.time().timeName() << nl 
                << "Cell " << curCell_lbl << "  Cp = " << Cp << "C1 cell " << C1_lbl << "  C1 = " << C1 << nl
                << "t1 = " << t1 << "  mag(t1) = " << magt1 << "  nf = " << nf << "  mag(nf) = " << mag(nf) 
                << "  costheta1 = " << costheta1 << nl 
                << endl;
            costheta1 = 1;
        }
        scalar theta1 = acos(costheta1);
        if(debug)
        {
            os<< "Cp = " << Cp << "  C1 = " << C1 << "  t1 = " << t1 << "  mag(t1) = " << magt1
                << "  costheta1 = " << costheta1 << "  theta1 = " << theta1                      
                << endl;
        }

        if(theta1 > 1E-3)
        {
            label C2_lbl = findCellInIntfcOrthDir(mesh,alpha1,curCells,C,Cp,C1,nf,curCell_lbl,C1_lbl,foundCell2,debug,os);
            if(debug)
            {
                os<< "C2_lbl: " << C2_lbl << "  found C2: " << foundCell2
                    << endl;
            }   

            if(foundCell2)
            {                                        
                vector C2 = C[C2_lbl];    
                vector t2 = C2 - Cp;
                scalar magt2;
                magt2 = mag(t2);
                if(magt2 < SMALL)
                {
                    os<< "Time = " << mesh.time().timeName() << nl 
                        << "Cell " << curCell_lbl << "  Cp = " << Cp << "  C2 = " << C2 << "  magt2 = " << magt2 << nl
                        << endl;
                    magt2 += SMALL;
                }    
                scalar magt1t2 = magt1*magt2;
                if(magt1t2 < SMALL)
                {
                    os<< "Time = " << mesh.time().timeName() << nl 
                        << "Cell " << curCell_lbl << "  magt1 = " << magt1 << "  magt2 = " << magt2 << "  magt1t2 = " << magt1t2 << nl
                        << endl;
                    magt1t2 += SMALL;
                }        
        
                scalar costheta2 = (nf & t2)/magt2/magnf;
                if(costheta2 > 1)
                {                
                    os<< "Time = " << mesh.time().timeName() << nl 
                        << "Cell " << curCell_lbl << "  Cp = " << Cp << "C2 cell " << C2_lbl << "  C2 = " << C2 << nl
                        << "t2 = " << t2 << "  mag(t2) = " << magt2 << "  nf = " << nf << "  mag(nf) = " << mag(nf) 
                        << "  costheta2 = " << costheta2 << nl 
                        << endl;
                    costheta2 = 1;
                }
                scalar theta2 = acos(costheta2);
                if(debug)
                {
                    os<< "t2 = " << t2 << "  mag(t2) = " << magt2
                        << "  costheta2 = " << costheta2 << "  theta2 = " << theta2                      
                        << endl;
                }
                scalar theta2_sign = -((nf ^ t1) & (nf ^ t2))/magt1t2/magnf/magnf;
                if(debug)
                {
                    os<< "theta2_sign = " << theta2_sign            
                        << endl;
                }

                scalar sintheta12 = sin(theta1 + theta2);

                if(theta1 + theta2 > constant::mathematical::pi)
                {
                    os<< "Time = " << mesh.time().timeName() << nl 
                        << "Cell " << curCell_lbl << "  nf" << nf << "  Cp = " << Cp << "  phaseLabel: " << phaseLbl << "  cell alpha1 = " << alpha1[curCell_lbl] << "  cell meshC = " << mesh.C()[curCell_lbl] << nl 
                        << "C1_lbl: " << C1_lbl << "  C1 = " << C1 << "  C2_lbl: " << C2_lbl << "  C2 = " << C2 << nl
                        << "magt1 = " << magt1 << "  magt2 = " << magt2 << "  magt1t2 = " << magt1t2 << nl
                        << "theta1 = " << theta1 << "  theta2 = " << theta2 << "  (theta1+theta2) = " << theta1 + theta2 << " sin(theta1+theta2) = " << sintheta12 << nl
                        << endl;

                    C1_lbl = findCellInIntfcDir(mesh,alpha1,curCells,C,Cp,nf,curCell_lbl,foundCell1,1,os);
                    C2_lbl = findCellInIntfcOrthDir(mesh,alpha1,curCells,C,Cp,C1,nf,curCell_lbl,C1_lbl,foundCell2,1,os);

                    os<< endl;
                }

                if(sintheta12 < SMALL)
                {
                    os<< "Time = " << mesh.time().timeName() << nl 
                        << "Cell " << curCell_lbl << "  theta1 = " << theta1 << "  theta2 = " << theta2 << "  sintheta12 = " << sintheta12 << nl
                        << endl;
                    sintheta12 += SMALL;
                }

                scalar alpha = sin(theta2)/sintheta12;        
                scalar beta = sin(theta1)/sintheta12;
                scalar at2bt1 = alpha*magt2 + beta*magt1;
                if(mag(at2bt1) < SMALL)
                {
                    os<< "Time = " << mesh.time().timeName() << nl 
                        << "Cell " << curCell_lbl << "  theta1 = " << theta1 << "  theta2 = " << theta2 << "  alpha = " << alpha 
                        << "  beta = " << beta << "  at2bt1 = " << at2bt1 << nl
                        << endl;
                    at2bt1 += SMALL;
                }

                if(debug)
                {
                    os<< "alpha = " << alpha << "  beta = " << beta << nl
                        << "at2bt1 = " << at2bt1
                        << endl;
                }
    
                dn = magt1t2/at2bt1;
                scalar T1; scalar T2;
                T1 = T[C1_lbl];
                T2 = T[C2_lbl];
                Teff = (alpha*magt2*T1 + beta*magt1*T2)/at2bt1;
                if(debug)
                {
                    os<< "dn = " << dn << "  Teff = " << Teff
                        << endl;
                }

                scalar Y1; scalar Y2;                
                if(debug)
                {
                    os<< "Teff = " << Teff
                        << endl;
                }
                for(label i=0; i<nSpecies; i++)
                {
                    Y1 = Y[i][C1_lbl];
                    Y2 = Y[i][C2_lbl];
                    Yeff[i] = (alpha*magt2*Y1 + beta*magt1*Y2)/at2bt1;

                    if(debug)
                    {
                        os<< "Species index: " << i << nl
                            << "Y1 = " << Y1 << "  Y2 = " << Y2 << "  Yeff = " << Yeff[i]
                            << endl;
                    }
                }             

                if(debug)
                {
                    os<< endl;        
                }
            }
            else
            {
                if(costheta1 < 1E-3)
                {        
                    os<< "Time = " << mesh.time().timeName() << nl 
                        << "Cell " << curCell_lbl << "  Cp = " << Cp << "  C1 cell " << C1_lbl << "  C1 = " << C1 << nl
                        << "t1 = " << t1 << "  mag(t1) = " << magt1 << "  nf = " << nf << "  mag(nf) = " << mag(nf) 
                        << "  costheta1 = " << costheta1 << nl 
                        << endl;
                    costheta1 = 1E-3;
                }
                dn = magt1/costheta1;
                Teff = T[C1_lbl];
                if(debug)
                {
                    os<< "dn = " << dn << "  Teff = " << Teff
                        << endl;
                }
                
                for(label i=0; i<nSpecies; i++)
                {             
                    Yeff[i] = Y[i][C1_lbl];
                    if(debug)
                    {
                        os<< "Species index: " << i << nl
                            << "Y1 = " << Y[i][C1_lbl] << "  Yeff = " << Yeff[i]
                            << endl;
                    }
                }
            }
        }
        else
        {
            dn = magt1;
            Teff = T[C1_lbl];
            if(debug)
            {
                os<< "dn = " << dn << "  Teff = " << Teff
                    << endl;
            }
            
            for(label i=0; i<nSpecies; i++)
            {             
                Yeff[i] = Y[i][C1_lbl];
                if(debug)
                {
                    os<< "Species index: " << i << nl
                        << "Y1 = " << Y[i][C1_lbl] << "  Yeff = " << Yeff[i]
                        << endl;
                }
            }
        }
    }
    else
    {
        os<< "Time = " << mesh.time().timeName() << nl 
            << "Cell " << curCell_lbl << ": Cell in interface normal direction not found! Using curCell centroid" << nl
            << endl;
        /*
            vector C1 = C[curCell_lbl];
            vector t1 = C1 - Cp;
            scalar magt1;    
            magt1 = mag(t1);
            scalar magnf;
            magnf = mag(nf);
            */
	scalar cellVol = mesh.V()[curCell_lbl];
        scalar min_dn = pow(cellVol, 1.0/3.0);
        /*
            if(magt1 < min_dn)
            {
            os<< "Time = " << mesh.time().timeName() << nl 
            << "Cell " << curCell_lbl << ": curCell centroid too close to interface!" << nl 
            << "Cp = " << Cp << "  C1 = " << C1 << "  magt1 = " << magt1 << nl
            << endl;
            if(magt1 < SMALL)
            {
            magt1 += SMALL;
            }
            t1 *= min_dn/magt1;
            magt1 = mag(t1);
            }
            if(magnf < SMALL)
            {
            os<< "Time = " << mesh.time().timeName() << nl 
            << "Cell " << curCell_lbl << "  nf = " << nf << "  mag(nf) = " << magnf << nl
            << endl;
            magnf += SMALL;
            }
            scalar costheta1 = (nf & t1)/magt1/magnf;
            if(costheta1 > 1)
            {        
            os<< "Time = " << mesh.time().timeName() << nl 
            << "Cell " << curCell_lbl << "  Cp = " << Cp << "C1 cell " << C1_lbl << "  C1 = " << C1 << nl
            << "t1 = " << t1 << "  mag(t1) = " << magt1 << "  nf = " << nf << "  mag(nf) = " << mag(nf) 
            << "  costheta1 = " << costheta1 << nl 
            << endl;
            costheta1 = 1;
            }
            */
        dn = min_dn;
        Teff = T[curCell_lbl];
        if(debug)
        {
            os<< "dn = " << dn << "  Teff = " << Teff
                << endl;
        }
        
        for(label i=0; i<nSpecies; i++)
        {             
            Yeff[i] = Y[i][curCell_lbl];
            if(debug)
            {
                os<< "Species index: " << i << nl
                    << "Y1 = " << Y[i][curCell_lbl] << "  Yeff = " << Yeff[i]
                    << endl;
            }
        }
    }
    
    if(dn < SMALL)
    {
        dn += SMALL;
    }
    if(debug)
    {
        os<< "dn stab = " << dn << endl;
    }
}


void calc_Xs_Ys_Js_mS_alphaS
(
    const fvMesh& mesh,
    double *Pc,
    double *Tc,
    double* Vc,
    double *w,
    double *MW,
    double *Tb,
    double *SG,
    double *H8,
    double *k,
    double *dm,
    double *kij_T,
    double Ta_kij,
    double Tb_kij,
    int nT_kij,
    double *kij,
    const labelListList& cellStencil,
    const List<List<scalar> >& x1_flatFld,
    const List<List<scalar> >& x0_flatFld,
    const List<scalar>& T1_flatFld,
    const List<scalar>& T0_flatFld,
    const List<scalar>& alpha1_flatFld,
    const volScalarField& alpha1,
    const volScalarField& rho1,
    const volScalarField& rho0,
    const volScalarField& T1,
    const volScalarField& T0,
    volScalarField& Ts,
    const volScalarField& P,
    const PtrList<volScalarField>& C1,
    const PtrList<volScalarField>& C0,
    const PtrList<volScalarField>& X1,
    const PtrList<volScalarField>& X0,
    const PtrList<volScalarField>& Y1,
    const PtrList<volScalarField>& Y0,
    const List<vector>& C_ph1_flatFld,
    const List<vector>& C_ph0_flatFld,
    const volVectorField& C_intfc,
    const volScalarField& A_intfc,
    const volVectorField& nHat,
    double dt,
    PtrList<volScalarField>& Ys1,
    PtrList<volScalarField>& Ys0,
    PtrList<volScalarField>& Xs1,
    PtrList<volScalarField>& Xs0,
    PtrList<volScalarField>& Js1,
    PtrList<volScalarField>& Js0,
    PtrList<volScalarField>& mS1,
    PtrList<volScalarField>& mS0,
    volScalarField& JsTot,
    volScalarField& mS1Tot,
    volScalarField& mS0Tot,
    volScalarField& alphaS1,
    volScalarField& alphaS0,
    volScalarField& Qs,
    labelList& n_iters_Ts,
    labelList& status_transLLE,
    labelList& cell_had_intfc,
    int n_flux_type,
    LPT_UMFPACK flux_umf,
    const label& nSpecies,
    const scalar& ALPHA_2PH_MIN,
    const scalar& A_INTFC_2PH_MIN,
    double Ts_TOL,
    int MAX_ITERS_Ts,
    double MASS_FRAC_TOL,
    const List<scalar>& xs1_0,
    const List<scalar>& xs0_0,
    const bool debug,
    OFstream& os
)
{
    int n, i, iLLE, iTs;
    label curCell_had_intfc;
    double alpha1_cellI, A_intfc_cellI, rho1_cellI, rho0_cellI, V_cellI; 
    double Ts_cellI, Ts_cellI_old, T1_cellI, T0_cellI, P_cellI;
    double dn1, dn0, Teff1, Teff0, JsTot_cellI, mS1Tot_cellI, Qs_cellI, conds1;
    scalar mS1Tot_cellI_tmp;//, mS1i_cellI, max_mSi;
    scalar limiterTot, limiter_min;
    vector nf, C_intfc_cellI;
    labelList curCellsAll = cellStencil[0];    
    double *xeff1, *xeff0;//, *C1_cellI, *C0_cellI;
    double *x1_cellI, *x0_cellI;//, *y1_cellI, *y0_cellI;
    double *xs1, *xs0, *ys1, *ys0, *flux_m_1, *flux_m_0, *Js1_cellI, *Js0_cellI, *hs1;
    //double *mS1_cellI, *limiter_mS1;

    n = nSpecies;

    _NNEW_(xeff1, double, n);
    _NNEW_(xeff0, double, n);
    //_NNEW_(C1_cellI, double, n);
    //_NNEW_(C0_cellI, double, n);
    _NNEW_(x1_cellI, double, n);
    _NNEW_(x0_cellI, double, n);
    //_NNEW_(y1_cellI, double, n);
    //_NNEW_(y0_cellI, double, n);
    _NNEW_(xs1, double, n);
    _NNEW_(xs0, double, n);
    _NNEW_(ys1, double, n);
    _NNEW_(ys0, double, n);
    _NNEW_(flux_m_1, double, n);
    _NNEW_(flux_m_0, double, n);
    _NNEW_(Js1_cellI, double, n);
    _NNEW_(Js0_cellI, double, n);
    //_NNEW_(mS1_cellI, double, n);
    //_NNEW_(limiter_mS1, double, n);
    _NNEW_(hs1, double, n);

    List<scalar> limiterY(n);
    List<scalar> C1_cellI(n);
    List<scalar> C0_cellI(n);
    List<scalar> Y1_cellI(n);
    List<scalar> Y0_cellI(n);
    List<scalar> mS1_cellI(n);

    scalar ALPHA_2PH_MAX = 1 - ALPHA_2PH_MIN;

    /*
        const labelList& own = mesh.owner();
        const labelList& nei = mesh.neighbour();

        const surfaceVectorField& Sf = mesh.Sf();
        const surfaceScalarField& magSf = mesh.magSf();
        const surfaceVectorField& Cf = mesh.Cf();
        */

    const scalarField& V = mesh.V();

    const scalarField& alpha1Cells = alpha1.internalField();    
    const vectorField& C_intfcCells = C_intfc.internalField();
    const scalarField& A_intfcCells = A_intfc.internalField();
    const vectorField& nHatCells = nHat.internalField();
    const scalarField& rho1Cells = rho1.internalField();
    const scalarField& rho0Cells = rho0.internalField();
    const scalarField& T1Cells = T1.internalField();
    const scalarField& T0Cells = T0.internalField();
    scalarField& TsCells = Ts.internalField();
    const scalarField& PCells = P.internalField();
    scalarField& JsTotCells = JsTot.internalField();
    scalarField& mS1TotCells = mS1Tot.internalField();
    scalarField& mS0TotCells = mS0Tot.internalField();
    scalarField& alphaS1Cells = alphaS1.internalField();
    scalarField& alphaS0Cells = alphaS0.internalField();
    scalarField& QsCells = Qs.internalField();

    if(debug)
    {
        print_line(os, 80);
        print_line(os, 80);
        os<< "Interfacial Species Fluxes Calculation with Transport-LLE Constraints" << endl;
        print_line(os, 80);
        os<< endl;
    }

    //Js, mS, alphaS for all interface cells
    forAll(alpha1Cells, cellI)
    {
        alpha1_cellI = alpha1Cells[cellI];
        A_intfc_cellI = A_intfcCells[cellI];

        if(alpha1_cellI > ALPHA_2PH_MIN && alpha1_cellI < ALPHA_2PH_MAX)
        {
            V_cellI = V[cellI];
            nf = nHatCells[cellI];
            C_intfc_cellI = C_intfcCells[cellI];
            curCellsAll = cellStencil[cellI];
            rho1_cellI = rho1Cells[cellI];
            rho0_cellI = rho0Cells[cellI];
            T1_cellI = T1Cells[cellI];
            T0_cellI = T0Cells[cellI];
            Ts_cellI_old = TsCells[cellI];
            P_cellI = PCells[cellI];
            for(i=0; i<n; i++)
            {
                C1_cellI[i] = C1[i].internalField()[cellI];
                C0_cellI[i] = C0[i].internalField()[cellI];
                x1_cellI[i] = X1[i].internalField()[cellI];
                x0_cellI[i] = X0[i].internalField()[cellI];
                //y1_cellI[i] = Y1[i].internalField()[cellI];
                //y0_cellI[i] = Y0[i].internalField()[cellI];
                Y1_cellI[i] = Y1[i].internalField()[cellI];
                Y0_cellI[i] = Y0[i].internalField()[cellI];
            }
            curCell_had_intfc = cell_had_intfc[cellI];

            //phase-1
            //ensure nf direction is into the phase
            calc_cell_intfcGrad_coeffs(mesh, cellI, nf, C_intfc_cellI, x1_flatFld, T1_flatFld, alpha1_flatFld, C_ph1_flatFld, curCellsAll, nSpecies, ALPHA_2PH_MIN, 1, dn1, xeff1, Teff1, debug, os);            

            //phase-0
            //ensure nf direction is into the phase
            //then reverse nf again for Js0 calculation            
            calc_cell_intfcGrad_coeffs(mesh, cellI, -nf, C_intfc_cellI, x0_flatFld, T0_flatFld, alpha1_flatFld, C_ph0_flatFld, curCellsAll, nSpecies, ALPHA_2PH_MIN, 0, dn0, xeff0, Teff0, debug, os);

            if(curCell_had_intfc == 1)
            {
                Ts_cellI = Ts_cellI_old;

                for(i=0; i<n; i++)
                {
                    xs1[i] = x1_cellI[i];
                    xs0[i] = x0_cellI[i];
                }
            }
            else
            {
                Ts_cellI = 0.5*(T1_cellI + T0_cellI);

                for(i=0; i<n; i++)
                {
                    xs1[i] = Xs1[i].internalField()[cellI];
                    xs0[i] = Xs0[i].internalField()[cellI];
                }
            }

            x2y(n, MW, xs1, ys1);
            x2y(n, MW, xs0, ys0);

            calc_intfc_transLLE(P_cellI, Ts_cellI, n, Pc, Tc, Vc, w, MW, Tb, SG, H8, k, dm, kij_T, Ta_kij, Tb_kij, nT_kij, kij, dn1, dn0, xeff1, xeff0, Teff1, Teff0, xs1, xs0, ys1, ys0, JsTot_cellI, flux_m_1, flux_m_0, conds1, hs1, iLLE, iTs, n_flux_type, flux_umf, Ts_TOL, MAX_ITERS_Ts, MASS_FRAC_TOL, debug, os);

            if(iLLE <= 0)
            {
                os<< "Time = " << mesh.time().timeName() << nl
                    << "GSL calculation failed in cell " << cellI << "  iLLE = " << iLLE << endl;

                JsTotCells[cellI] = 0;
                mS1TotCells[cellI] = 0;
                mS0TotCells[cellI] = 0;
                for(i=0; i<n; i++)
                {                
                    mS1[i].internalField()[cellI] = 0;
                    mS0[i].internalField()[cellI] = 0;
                    Js1[i].internalField()[cellI] = 0;
                    Js0[i].internalField()[cellI] = 0;                
                }
                alphaS1Cells[cellI] = 0;
                alphaS0Cells[cellI] = 0;
                QsCells[cellI] = 0;

                n_iters_Ts[cellI] = 0;
                status_transLLE[cellI] = iLLE;
                cell_had_intfc[cellI] = 0;
            }
            else
            {            
                mS1Tot_cellI = 0; mS1Tot_cellI_tmp = 0;

                for(i=0; i<n; i++)
                {
                    Js1_cellI[i] = -A_intfc_cellI*flux_m_1[i]/V_cellI;
                    Js0_cellI[i] = -A_intfc_cellI*flux_m_0[i]/V_cellI;                
                }
                mS1Tot_cellI_tmp = -A_intfc_cellI*JsTot_cellI/V_cellI;

                if(mS1Tot_cellI_tmp > 0)
                {
                    for(i=0; i<n; i++)
                    {
                        mS1_cellI[i] = mS1Tot_cellI_tmp*Y0_cellI[i] + Js0_cellI[i];
                    }
                }
                else
                {
                    for(i=0; i<n; i++)
                    {
                        mS1_cellI[i] = mS1Tot_cellI_tmp*Y1_cellI[i] + Js1_cellI[i];
                    }
                }

                limiterTot = 1;
                limiter_min = 1;
                for(i=0; i<n; i++) limiterY[i] = 1;

                calc_mS_limiter(C1_cellI, C0_cellI, Y1_cellI, Y0_cellI, mS1_cellI, mS1Tot_cellI_tmp, dt, n, limiterTot, limiterY, limiter_min, debug, os);

                for(i=0; i<n; i++)
                {                                
                    mS1_cellI[i] = limiter_min*mS1_cellI[i];                
                }
            
                mS1Tot_cellI = limiter_min*mS1Tot_cellI_tmp;            

                //calculate the interfacial enthalpy transfer
                Qs_cellI = -A_intfc_cellI*conds1*(Teff1 - Ts_cellI)/dn1/V_cellI;
                for(i=0; i<n; i++)
                {
                    Qs_cellI += mS1_cellI[i]*hs1[i]/(MW[i]*1e-3);
                }
            
                //assign the calculated values to corresponding fields
                //TsCells[cellI] = Ts_cellI;
                mS1TotCells[cellI] = mS1Tot_cellI;
                mS0TotCells[cellI] = -mS1Tot_cellI;
                JsTotCells[cellI] = mS1Tot_cellI*V_cellI;
                for(i=0; i<n; i++)
                {                
                    mS1[i].internalField()[cellI] = mS1_cellI[i];
                    mS0[i].internalField()[cellI] = -mS1_cellI[i];
                    Js1[i].internalField()[cellI] = Js1_cellI[i];
                    Js0[i].internalField()[cellI] = Js0_cellI[i];
                    Xs1[i].internalField()[cellI] = xs1[i];
                    Xs0[i].internalField()[cellI] = xs0[i];
                    Ys1[i].internalField()[cellI] = ys1[i];
                    Ys0[i].internalField()[cellI] = ys0[i];
                }
                alphaS1Cells[cellI] = mS1Tot_cellI/rho1_cellI;
                alphaS0Cells[cellI] = -mS1Tot_cellI/rho0_cellI;
                QsCells[cellI] = Qs_cellI;            

                n_iters_Ts[cellI] = iTs;
                status_transLLE[cellI] = iLLE;            
                cell_had_intfc[cellI] = 1;

                if(debug)
                {
                    print_line(os, 100);
                    print_line(os, 100);                
                    os<< "Cell: " << cellI << "  alpha1 = " << alpha1_cellI << nl
                        << "nf = " << nf << "  C_intfc = " << C_intfc_cellI << "  A_intfc = " << A_intfc_cellI << nl
                        << "C_ph1 = " << C_ph1_flatFld[cellI] << "  C_ph0 = " << C_ph0_flatFld[cellI] << nl
                        << "T1 = " << T1Cells[cellI] << "  T0 = " << T0Cells[cellI] << "  rho1 = " << rho1Cells[cellI] << "  rho0 = " << rho0Cells[cellI] << endl;
                    print_line(os, 100);
                    os<< setw(8) << "Species" << "  " << setw(14) << "x1" << "  " << setw(14) << "x0" << "  " << setw(14) << "y1" << "  " << setw(14) << "y0" << "  " << setw(14) << "C1" << "  " << setw(14) << "C0" << endl;
                    print_line(os, 100);
                    for(i=0; i<n; i++)
                    {
                        os<< setw(8) << i << "  " << setw(14) << x1_cellI[i] << "  " << setw(14) << x0_cellI[i] << "  " << setw(14) << Y1_cellI[i] << "  " << setw(14) << Y0_cellI[i] << "  " << setw(14) << C1_cellI[i] << "  " << setw(14) << C0_cellI[i] << endl;
                    }
                    print_line(os, 100);
                    os<< "dn1 = " << dn1 << "  dn0 = " << dn0 << "  Teff1 = " << Teff1 << "  Teff0 = " << Teff0 << endl;
                    print_line(os, 100);
                    os<< setw(8) << "Species" << "  " << setw(14) << "xeff1" << "  " << setw(14) << "xeff0" << endl;
                    print_line(os, 100);
                    for(i=0; i<n; i++)
                    {
                        os<< setw(8) << i << "  " << setw(14) << xeff1[i] << "  " << setw(14) << xeff0[i] << endl;                
                    }
                    print_line(os, 100);
                    os<< "Ts = " << Ts_cellI << "  mS1Tot = " << mS1Tot_cellI << "  alphaS1 = " << alphaS1Cells[cellI] << "  alphaS0 = " << alphaS0Cells[cellI] << nl
                        //<< "limiter_mS1Tot = " << limiter_mS1Tot << "  limiter_min = " << limiter_min << nl
                        << "iLLE = " << iLLE << "  iTs = " << iTs << endl;
                    print_line(os, 100);
                    os<< setw(8) << "Species" << "  " << setw(14) << "xs1" << "  " << setw(14) << "xs0" << "  " << setw(14) << "ys1" << "  " << setw(14) << "ys0" << endl;
                    print_line(os, 100);
                    for(i=0; i<n; i++)
                    {
                        os<< setw(8) << i << "  " << setw(14) << xs1[i] << "  " << setw(14) << xs0[i] << "  " << setw(14) << ys1[i] << "  " << setw(14) << ys0[i] << endl;
                    }
                    print_line(os, 100);
                    os<< setw(8) << "Species" << "  " << setw(14) << "flux_m_1" << "  " << setw(14) << "flux_m_0" << "  " << setw(14) << "Js1" << "  " << setw(14) << "Js0" << "  " << setw(14) << "mS1" << endl;
                    print_line(os, 100);
                    for(i=0; i<n; i++)
                    {
                        os<< setw(8) << i << "  " << setw(14) << flux_m_1[i] << "  " << setw(14) << flux_m_0[i] << "  " << setw(14) << Js1_cellI[i] << "  " << setw(14) << Js0_cellI[i] << "  " << setw(14) << mS1_cellI[i] << endl;
                    }
                    print_line(os, 100);
                    print_line(os, 100);
                }
            }
        }
        else
        {
            //assign the calculated values to corresponding fields
            JsTotCells[cellI] = 0;
            mS1TotCells[cellI] = 0;
            mS0TotCells[cellI] = 0;
            for(i=0; i<n; i++)
            {                
                mS1[i].internalField()[cellI] = 0;
                mS0[i].internalField()[cellI] = 0;
                Js1[i].internalField()[cellI] = 0;
                Js0[i].internalField()[cellI] = 0;                
            }
            alphaS1Cells[cellI] = 0;
            alphaS0Cells[cellI] = 0;
            QsCells[cellI] = 0;

            n_iters_Ts[cellI] = 0;
            status_transLLE[cellI] = 0;            
            cell_had_intfc[cellI] = 0;
        }        
    }

    /*
    //Js, mS, alphaS for all internal faces very close to interface
    for(faceI=0; faceI<mesh.nInternalFaces(); faceI++)
    {
        faceOwn = own[faceI];
        faceNei = nei[faceI];
        alpha1Own = alpha1Cells[faceOwn];
        alpha1Nei = alpha1Cells[faceNei];        

        if(alpha1Own >= ALPHA_2PH_MAX && alpha1Nei <= ALPHA_2PH_MIN)
        {
            A_intfc_cellI = magSf[faceI];
            nf = -Sf[faceI]/A_intfc_cellI;
            C_intfc_cellI = Cf[faceI];

            VOwn = V[faceOwn];
            VNei = V[faceNei];
            rho1Own  = rho1Cells[faceOwn];
            rho0Own  = rho0Cells[faceOwn];
            rho1Nei = rho1Cells[faceNei];
            rho0Nei = rho0Cells[faceNei];
            TsOwn_old = TsCells[faceOwn];
            TsNei_old = TsCells[faceNei];
            T1Own = T1Cells[faceOwn];
            T0Nei = T0Cells[faceNei];            
            P_cellI = PCells[faceOwn];
            for(i=0; i<n; i++)
            {                
                C1Own[i] = C1[i].internalField()[faceOwn];
                C0Own[i] = C0[i].internalField()[faceOwn];
                C1Nei[i] = C1[i].internalField()[faceNei];
                C0Nei[i] = C0[i].internalField()[faceNei];
            }
            own_had_intfc = cell_had_intfc[faceOwn];
            nei_had_intfc = cell_had_intfc[faceNei];

            //phase-1
            //ensure nf direction is into the phase
            curCellsAll = cellStencil[faceOwn];
            calc_cell_intfcGrad_coeffs(mesh, cellI, nf, C_intfc_cellI, x1_flatFld, T1_flatFld, alpha1_flatFld, C_ph1_flatFld, curCellsAll, nSpecies, ALPHA_2PH_MIN, 1, dn1, xeff1, Teff1, debug, os);            

            //phase-0
            //ensure nf direction is into the phase
            //then reverse nf again for Js0 calculation
            curCellsAll = cellStencil[faceNei];
            calc_cell_intfcGrad_coeffs(mesh, cellI, -nf, C_intfc_cellI, x0_flatFld, T0_flatFld, alpha1_flatFld, C_ph0_flatFld, curCellsAll, nSpecies, ALPHA_2PH_MIN, 0, dn0, xeff0, Teff0, debug, os);

            if(own_had_intfc == 1)
            {
                Ts_cellI = TsOwn_old;
            }
            else
            {
                if(nei_had_intfc == 1)
                {
                    Ts_cellI = TsNei_old;
                }
                else
                {
                    Ts_cellI = 0.5*(T1Own + T0Nei);
                }
            }

            calc_intfc_transLLE(P_cellI, Ts_cellI, n, Pc, Tc, Vc, w, MW, tk, coef_ab, Tb, SG, H8, k, dm, dn1, dn0, xeff1, xeff0, Teff1, Teff0, xs1, xs0, ys1, ys0, JsTot_cellI, flux_m_1, flux_m_0, conds1, hs1, iLLE, iTs, n_flux_type, flux_umf, Ts_TOL, MAX_ITERS_Ts, MASS_FRAC_TOL, debug, os);

            for(i=0; i<n; i++)
            {
                Js1_cellI[i] = -A_intfc_cellI*flux_m_1[i];
                Js0_cellI[i] = -A_intfc_cellI*flux_m_0[i];
            }
            mS1Tot_cellI_tmp = -A_intfc_cellI*JsTot_cellI;

            mS1TotOwn = 0; mS1TotNei = 0;            

            for(i=0; i<n; i++)
            {
                mS1i_cellI = mS1Tot_cellI_tmp*ys1[i] + Js1_cellI[i];
                if(mS1i_cellI > 0)
                {
                    max_mSi = C0Own[i]*VOwn/dt;
                    mS1Own[i] = min(max_mSi, mS1i_cellI);
                    max_mSi = C0Nei[i]*VNei/dt;
                    mS1Nei[i] = min(max_mSi, (mS1i_cellI - mS1Own[i]));
                }
                else
                {
                    max_mSi = C1Nei[i]*VNei/dt;
                    mS1Nei[i] = -min(max_mSi, -mS1i_cellI);
                    max_mSi = C1Own[i]*VOwn/dt;
                    mS1Own[i] = -min(max_mSi, (-mS1i_cellI + mS1Nei[i]));
                }

                mS1TotOwn += mS1Own[i];
                mS1TotNei += mS1Nei[i];

                mS1[i].internalField()[faceOwn] = mS1Own[i]/VOwn;
                mS1[i].internalField()[faceNei] = mS1Nei[i]/VNei;
            }
            
            mS1TotCells[faceOwn] = mS1TotOwn/VOwn;
            mS1TotCells[faceNei] = mS1TotNei/VNei;

            alphaS1Cells[faceOwn] = mS1TotCells[faceOwn]/rho1Own;
            alphaS0Cells[faceOwn] = -mS1TotCells[faceOwn]/rho0Own;
            alphaS1Cells[faceNei] = mS1TotCells[faceNei]/rho1Nei;
            alphaS0Cells[faceNei] = -mS1TotCells[faceNei]/rho0Nei;

            //calculate the interfacial enthalpy transfer
            Qs_cellI = -A_intfc_cellI*conds1*(Teff1 - Ts_cellI)/dn1/V_cellI;
            for(i=0; i<n; i++)
            {
                Qs_cellI += mS1_cellI[i]*hs1[i]/(MW[i]*1e-3);
            }
            
            //assign the calculated values to corresponding fields
            TsCells[cellI] = Ts_cellI;
            mS1TotCells[cellI] = mS1Tot_cellI;
            for(i=0; i<n; i++)
            {                
                mS1[i].internalField()[cellI] = mS1_cellI[i];
                Js1[i].internalField()[cellI] = Js1_cellI[i];
                Js0[i].internalField()[cellI] = Js0_cellI[i];
                Xs1[i].internalField()[cellI] = xs1[i];
                Xs0[i].internalField()[cellI] = xs0[i];
                Ys1[i].internalField()[cellI] = ys1[i];
                Ys0[i].internalField()[cellI] = ys0[i];
            }
            alphaS1Cells[cellI] = mS1Tot_cellI/rho1_cellI;
            alphaS0Cells[cellI] = -mS1Tot_cellI/rho0_cellI;
            QsCells[cellI] = Qs_cellI;            

            n_iters_Ts[cellI] = iTs;
            status_transLLE[cellI] = iLLE;            
            cell_had_intfc[cellI] = 1;
        }

        if(alpha1Own <= ALPHA_2PH_MIN && alpha1Nei >= ALPHA_2PH_MAX)
        {
            VOwn = V[faceOwn];
            VNei = V[faceNei];
            rho1Own  = rho1Cells[faceOwn];
            rho0Own  = rho0Cells[faceOwn];
            rho1Nei = rho1Cells[faceNei];
            rho0Nei = rho0Cells[faceNei];
            for(i=0; i<n; i++)
            {
                Js1_cellI[i] = Js1[i].internalField()[faceOwn];
                Js0_cellI[i] = Js0[i].internalField()[faceOwn];
                Ys1_cellI[i] = Ys1[i].internalField()[faceOwn];
                Ys0_cellI[i] = Ys0[i].internalField()[faceOwn];
                C1Own[i] = C1[i].internalField()[faceOwn];
                C0Own[i] = C0[i].internalField()[faceOwn];
                C1Nei[i] = C1[i].internalField()[faceNei];
                C0Nei[i] = C0[i].internalField()[faceNei];
            }

            mS1TotOwn = 0; mS1TotNei = 0;

            mS1Tot_cellI_tmp = Js0_cellI[0] - Js1_cellI[0];
            mS1Tot_cellI_tmp /= (Ys1_cellI[0] - Ys0_cellI[0]);

            for(i=0; i<n; i++)
            {
                mS1i_cellI = mS1Tot_cellI_tmp*Ys1_cellI[i] + Js1_cellI[i];

                if(mS1i_cellI > 0)
                {
                    max_mSi = C0Nei[i]*VNei/dt;
                    mS1Nei[i] = min(max_mSi, mS1i_cellI);
                    max_mSi = C0Own[i]*VOwn/dt;
                    mS1Own[i] = min(max_mSi, (mS1i_cellI - mS1Nei[i]));
                }
                else
                {
                    max_mSi = C1Own[i]*VOwn/dt;
                    mS1Own[i] = -min(max_mSi, -mS1i_cellI);
                    max_mSi = C1Nei[i]*VNei/dt;
                    mS1Nei[i] = -min(max_mSi, (-mS1i_cellI + mS1Own[i]));
                }

                mS1TotOwn += mS1Own[i];
                mS1TotNei += mS1Nei[i];

                mS1[i].internalField()[faceOwn] = mS1Own[i]/VOwn;
                mS0[i].internalField()[faceOwn] = -mS1Own[i]/VOwn;

                mS1[i].internalField()[faceNei] = mS1Nei[i]/VNei;
                mS0[i].internalField()[faceNei] = -mS1Nei[i]/VNei;
            }
            
            mS1TotCells[faceOwn] = mS1TotOwn/VOwn;
            mS0TotCells[faceOwn] = -mS1TotOwn/VOwn;
            mS1TotCells[faceNei] = mS1TotNei/VNei;
            mS0TotCells[faceNei] = -mS1TotNei/VNei;

            alphaS1Cells[faceOwn] = mS1TotCells[faceOwn]/rho1Own;
            alphaS0Cells[faceOwn] = mS0TotCells[faceOwn]/rho0Own;
            alphaS1Cells[faceNei] = mS1TotCells[faceNei]/rho1Nei;
            alphaS0Cells[faceNei] = mS0TotCells[faceNei]/rho0Nei;
        }
    }*/

    _DDELETE_(xeff1);
    _DDELETE_(xeff0);
    //_DDELETE_(C1_cellI);
    //_DDELETE_(C0_cellI);
    _DDELETE_(x1_cellI);
    _DDELETE_(x0_cellI);
    //_DDELETE_(y1_cellI);
    //_DDELETE_(y0_cellI);
    _DDELETE_(xs1);
    _DDELETE_(xs0);
    _DDELETE_(ys1);
    _DDELETE_(ys0);
    _DDELETE_(flux_m_1);
    _DDELETE_(flux_m_0);
    _DDELETE_(Js1_cellI);
    _DDELETE_(Js0_cellI);
    //_DDELETE_(mS1_cellI);
    //_DDELETE_(limiter_mS1);
    _DDELETE_(hs1);

    if(debug)
    {
        os<< nl
            << "Done Interfacial Species Fluxes Calculation with Transport-LLE Constraints" << endl;
        print_line(os, 80);
        os<< endl;
    }
}


void calc_mS_limiter_eq
(
    const List<scalar>& C1_0,
    const List<scalar>& C0_0,
    const List<scalar>& Y1_0,
    const List<scalar>& Y0_0,
    const List<scalar>& Ys1,
    const List<scalar>& Ys0,
    const List<scalar>& mS1,
    const scalar& mS1Tot,
    const scalar& dt,    
    int n,
    scalar& limiter_min,
    const bool debug,
    OFstream& os
)
{
    int i;
    scalar m1_0, m0_0, mTot_0, m1_eq, mS1Tot_eq, limiterTot;
    List<scalar> mS1_eq(n);
    List<scalar> limiter(n);

    m1_0 = 0; m0_0 = 0;  
    for(i=0; i<n; i++)
    {
        m1_0 += C1_0[i];
        m0_0 += C0_0[i];
    }
    mTot_0 = m1_0 + m0_0;

    m1_eq = C1_0[0] + C0_0[0] - mTot_0*Ys0[0];
    m1_eq /= (Ys1[0] - Ys0[0]);

    if(m1_eq < 0)
    { 
        m1_eq = 0;
        mS1Tot_eq = (m1_eq - m1_0)/dt;
        for(i=0; i<n; i++) mS1_eq[i] = -C1_0[i]/dt;
    }
    else if(m1_eq > mTot_0)
    {
        m1_eq = mTot_0;
        mS1Tot_eq = (m1_eq - m1_0)/dt;
        for(i=0; i<n; i++) mS1_eq[i] = C0_0[i]/dt;
    }
    else
    {
        mS1Tot_eq = (m1_eq - m1_0)/dt;
        for(i=0; i<n; i++) mS1_eq[i] = (m1_eq*Ys1[i] - C1_0[i])/dt;
    }

    limiterTot = mag(mS1Tot_eq)/(max(mag(mS1Tot), SMALL));
    limiter_min = limiterTot;
    for(i=0; i<n; i++)
    {
        limiter[i] = mag(mS1_eq[i])/(max(mag(mS1[i]), SMALL));
        limiter_min = min(limiter[i], limiter_min);
    }
    limiter_min = min(limiter_min, 1);

    if(debug)
    {
        os<< "mTot_0 = " << mTot_0 << "  m1_0 = " << m1_0 << "  m0_0 = " << m0_0 << "  dt = " << dt << endl;
        print_line(os,100);
        os<< setw(8) << "Species" << " " << setw(14) << "C1_0" << " " << setw(14) << "C0_0" << " " << setw(14) << "Y1_0" << " " << setw(14) << "Y0_0" << endl;
        print_line(os,100);
        for(i=0; i<n; i++)
        {
            os<< setw(14) << i << " " << setw(14) << C1_0[i] << " " << setw(14) << C0_0[i] << " " << setw(14) << Y1_0[i] << " " << setw(14) << Y0_0[i] << endl;
        }
        print_line(os,100);
        os<< "m1_eq = " << m1_eq << "  mS1Tot_eq*dt = " << mS1Tot_eq*dt << "  mS1Tot*dt = " << mS1Tot*dt << "  limiterTot = " << limiterTot << endl;
        print_line(os,100);
        os<< setw(8) << "Species" << " " << setw(14) << "mS1_eq*dt" << " " << setw(14) << "mS1*dt" << " " << setw(14) << "limiter" << endl;
        print_line(os,100);
        for(i=0; i<n; i++)
        {
            os<< setw(14) << i << " " << setw(14) << mS1_eq[i]*dt << " " << setw(14) << mS1[i]*dt << " " << setw(14) << limiter[i] << endl;
        }
        print_line(os,100);
        os<< "m1_eq = " << m1_eq << "  mS1Tot_eq = " << mS1Tot_eq << "  mS1Tot = " << mS1Tot << "  limiterTot = " << limiterTot << endl;
        print_line(os,100);
        os<< setw(8) << "Species" << " " << setw(14) << "mS1_eq" << " " << setw(14) << "mS1" << " " << setw(14) << "limiter" << endl;
        print_line(os,100);
        for(i=0; i<n; i++)
        {
            os<< setw(14) << i << " " << setw(14) << mS1_eq[i] << " " << setw(14) << mS1[i] << " " << setw(14) << limiter[i] << endl;
        }
        print_line(os,100);
        os<< "limiter_min = " << limiter_min << endl;
        print_line(os,100);
        os<< endl;
    }
}


void calc_mS_limiter
(
    const List<scalar>& C1_0,
    const List<scalar>& C0_0,
    const List<scalar>& Y1_0,
    const List<scalar>& Y0_0,    
    const List<scalar>& mS1,
    const scalar& mS1Tot,
    const scalar& dt,    
    int n,
    scalar& limiterTot,
    List<scalar>& limiterY,
    scalar& limiter_min,
    const bool debug,
    OFstream& os
)
{
    int i;
    scalar m1_0, m0_0;    

    m1_0 = 0; m0_0 = 0;  
    for(i=0; i<n; i++)
    {
        m1_0 += C1_0[i];
        m0_0 += C0_0[i];
    }    

    if(mS1Tot > 0)
    {
        limiterTot = m0_0/dt/max(mS1Tot, SMALL);
    }
    else if(mS1Tot < 0)
    {
        limiterTot = m1_0/dt/max(-mS1Tot, SMALL);
    }
    else
    {
        limiterTot = 1;
    }

    limiter_min = limiterTot;

    for(i=0; i<n; i++)
    {
        if(mS1[i] > 0)
        {
            limiterY[i] = C0_0[i]/dt/max(mS1[i], SMALL);
        }
        else if(mS1[i] < 0)
        {
            limiterY[i] = C1_0[i]/dt/max(-mS1[i], SMALL);
        }
        else
        {
            limiterY[i] = 1;
        }

        limiter_min = min(limiterY[i], limiter_min);
    }

    limiter_min = min(limiter_min, 1);

    if(debug)
    {
        os<< "m1_0 = " << m1_0 << "  m0_0 = " << m0_0 << "  dt = " << dt << endl;
        print_line(os,100);
        os<< setw(8) << "Species" << " " << setw(14) << "C1_0" << " " << setw(14) << "C0_0" << " " << setw(14) << "Y1_0" << " " << setw(14) << "Y0_0" << endl;
        print_line(os,100);
        for(i=0; i<n; i++)
        {
            os<< setw(14) << i << " " << setw(14) << C1_0[i] << " " << setw(14) << C0_0[i] << " " << setw(14) << Y1_0[i] << " " << setw(14) << Y0_0[i] << endl;
        }
        print_line(os,100);
        os<< "mS1Tot*dt = " << mS1Tot*dt << "  limiterTot = " << limiterTot << endl;
        print_line(os,100);
        os<< setw(8) << "Species" << " " << setw(14) << "mS1*dt" << " " << setw(14) << "limiter" << endl;
        print_line(os,100);
        for(i=0; i<n; i++)
        {
            os<< setw(14) << i << " " << setw(14) << mS1[i]*dt << " " << setw(14) << limiterY[i] << endl;
        }
        print_line(os,100);
        os<< "mS1Tot = " << mS1Tot << "  limiterTot = " << limiterTot << endl;
        print_line(os,100);
        os<< setw(8) << "Species" << " " << setw(14) << "mS1" << " " << setw(14) << "limiter" << endl;
        print_line(os,100);
        for(i=0; i<n; i++)
        {
            os<< setw(14) << i << " " << setw(14) << mS1[i] << " " << setw(14) << limiterY[i] << endl;
        }
        print_line(os,100);
        os<< "limiter_min = " << limiter_min << endl;
        print_line(os,100);
        os<< endl;
    }
}


void calc_mS_alphaS
(
    const fvMesh& mesh,
    const PtrList<volScalarField>& C1,
    const PtrList<volScalarField>& C0,
    const PtrList<volScalarField>& Y1,
    const PtrList<volScalarField>& Y0,
    const volScalarField& alpha1,
    const volScalarField& rho1,
    const volScalarField& rho0,
    const PtrList<volScalarField>& Ys1,
    const PtrList<volScalarField>& Ys0,
    const PtrList<volScalarField>& Js1,
    const PtrList<volScalarField>& Js0,
    const label& nSpecies,
    const scalar& ALPHA_2PH_MIN,
    PtrList<volScalarField>& mS1,
    PtrList<volScalarField>& mS0,
    volScalarField& mS1Tot,
    volScalarField& mS0Tot,
    volScalarField& alphaS1,
    volScalarField& alphaS0,
    const scalar& dt,
    const bool debug,
    OFstream& os
)
{
    int n, i;
    //label faceI, bndFaceI, nBnd, faceOwn, faceNei;
    n = nSpecies;
    scalar alpha1_cellI, rho1_cellI, rho0_cellI, V_cellI, mS1Tot_cellI, mS1Tot_cellI_tmp, limiterTot, limiter_min;
    //scalar alpha1Own, alpha1Nei, VOwn, VNei, rho1Own, rho0Own, rho1Nei, rho0Nei, mS1TotOwn, mS1TotNei;
    List<scalar> mS1_cellI(n);
    List<scalar> Js1_cellI(n);
    List<scalar> Js0_cellI(n);
    List<scalar> Ys1_cellI(n);
    List<scalar> Ys0_cellI(n);
    List<scalar> Y1_cellI(n);
    List<scalar> Y0_cellI(n);
    List<scalar> C1_cellI(n);
    List<scalar> C0_cellI(n);
    List<scalar> limiterY(n);
    /*
    List<scalar> Y1Own(n);
    List<scalar> Y0Own(n);
    List<scalar> Y1Nei(n);
    List<scalar> Y0Nei(n);
    List<scalar> C1Own(n);
    List<scalar> C0Own(n);
    List<scalar> C1Nei(n);
    List<scalar> C0Nei(n);
    List<scalar> mS1Own(n);
    List<scalar> mS1Nei(n);
    List<scalar> mS0Own(n);
    List<scalar> mS0Nei(n);
        */    

    //nBnd = mesh.nFaces() - mesh.nInternalFaces();
    //List<scalar> VNeiFld(nBnd);

    
    scalar ALPHA_2PH_MAX = 1 - ALPHA_2PH_MIN;
    const scalarField& V = mesh.V();

    //const labelList& own = mesh.owner();
    //const labelList& nei = mesh.neighbour();

    const scalarField& alpha1Cells = alpha1.internalField();    
    const scalarField& rho1Cells = rho1.internalField();
    const scalarField& rho0Cells = rho0.internalField();
    scalarField& mS1TotCells = mS1Tot.internalField();
    scalarField& mS0TotCells = mS0Tot.internalField();
    scalarField& alphaS1Cells = alphaS1.internalField();
    scalarField& alphaS0Cells = alphaS0.internalField();

    if(debug)
    {
        print_line(os, 100);
        print_line(os, 100);
        os<< "Interfacial Species Transfer Source Terms Calculation" << endl;
        print_line(os, 100);
        print_line(os, 100);
        os<< "Internal cells" << endl;
        print_line(os, 100);
    }

    //mS for all interface cells
    forAll(alpha1Cells, cellI)
    {
        alpha1_cellI = alpha1Cells[cellI];        

        if(alpha1_cellI > ALPHA_2PH_MIN && alpha1_cellI < ALPHA_2PH_MAX)
        {
            rho1_cellI = rho1Cells[cellI];        
            rho0_cellI = rho0Cells[cellI];        
            V_cellI = V[cellI];        
            for(i=0; i<n; i++)
            {
                Js1_cellI[i] = Js1[i].internalField()[cellI]/V_cellI;
                Js0_cellI[i] = Js0[i].internalField()[cellI]/V_cellI;
                Ys1_cellI[i] = Ys1[i].internalField()[cellI];
                Ys0_cellI[i] = Ys0[i].internalField()[cellI];
                Y1_cellI[i] = Y1[i].internalField()[cellI];
                Y0_cellI[i] = Y0[i].internalField()[cellI];
                C1_cellI[i] = C1[i].internalField()[cellI];
                C0_cellI[i] = C0[i].internalField()[cellI];
            }
            mS1Tot_cellI = 0; mS1Tot_cellI_tmp = 0;
            
            mS1Tot_cellI_tmp = Js0_cellI[0] - Js1_cellI[0];
            if(mS1Tot_cellI_tmp > 0)
            {
                if((Ys1_cellI[0] - Y0_cellI[0]) > 0)
                {
                    mS1Tot_cellI_tmp /= (Ys1_cellI[0] - Y0_cellI[0]);
                }
                else
                {
                    mS1Tot_cellI_tmp /= (Ys1_cellI[0] - Ys0_cellI[0]);
                }

                for(i=0; i<n; i++)
                {                                
                    mS1_cellI[i] = mS1Tot_cellI_tmp*Ys1_cellI[i] + Js1_cellI[i];
                }
            }
            else
            {
                if((Y1_cellI[0] - Ys0_cellI[0]) > 0)
                {
                    mS1Tot_cellI_tmp /= (Y1_cellI[0] - Ys0_cellI[0]);
                }
                else
                {
                    mS1Tot_cellI_tmp /= (Ys1_cellI[0] - Ys0_cellI[0]);
                }

                for(i=0; i<n; i++)
                {                                
                    mS1_cellI[i] = mS1Tot_cellI_tmp*Ys0_cellI[i] + Js0_cellI[i];
                }
            }

            if(debug)
            {
                print_line(os, 100);
                print_line(os, 100);
                os<< "Cell: " << cellI << "  alpha1 = " << alpha1_cellI << endl;
                print_line(os, 100);
                os<< "mS limiter calculation" << endl;
                print_line(os, 100);
            }

            limiterTot = 1;
            limiter_min = 1;
            for(i=0; i<n; i++) limiterY[i] = 1;

            calc_mS_limiter(C1_cellI, C0_cellI, Y1_cellI, Y0_cellI, mS1_cellI, mS1Tot_cellI_tmp, dt, n, limiterTot, limiterY, limiter_min, debug, os);     
     
            for(i=0; i<n; i++)
            {                                
                mS1_cellI[i] = limiter_min*mS1_cellI[i];                

                mS1[i].internalField()[cellI] = mS1_cellI[i];
                mS0[i].internalField()[cellI] = -mS1_cellI[i];
            }
            
            mS1Tot_cellI = limiter_min*mS1Tot_cellI_tmp;

            mS1TotCells[cellI] = mS1Tot_cellI;
            mS0TotCells[cellI] = -mS1Tot_cellI;

            alphaS1Cells[cellI] = mS1Tot_cellI/rho1_cellI;
            alphaS0Cells[cellI] = -mS1Tot_cellI/rho0_cellI;

            if(debug)
            {
                print_line(os, 100);
                print_line(os, 100);
                os<< "Cell: " << cellI << "  alpha1 = " << alpha1_cellI << endl;
                print_line(os, 100);
                os<< setw(8) << "Species" << "  " << setw(14) << "C1" << "  " << setw(14) << "C0" << "  " << setw(14) << "Y1" << "  " << setw(14) << "Y0" << endl;
                print_line(os, 100);
                for(i=0; i<n; i++)
                {
                    os<< setw(7) << i << "  " << setw(14) << C1[i].internalField()[cellI] << "  " << setw(14) << C0[i].internalField()[cellI] << "  " << setw(14) << Y1[i].internalField()[cellI] << "  " << setw(14) << Y0[i].internalField()[cellI] << endl;
                }
                print_line(os, 100);
                os<< setw(8) << "Species" << "  " << setw(14) << "Js1" << "  " << setw(14) << "Js0" << "  " << setw(14) << "mS1" << "  " << setw(14) << "mS0" << endl;
                print_line(os, 100);
                for(i=0; i<n; i++)
                {
                    os<< setw(7) << i << "  " << setw(14) << Js1[i].internalField()[cellI] << "  " << setw(14) << Js0[i].internalField()[cellI] << "  " << setw(14) << mS1[i].internalField()[cellI] << "  " << setw(14) << mS0[i].internalField()[cellI] << endl;
                }
                print_line(os, 100);
                os<< "limiter_min = " << limiter_min << nl
                    << "mS1Tot = " << mS1TotCells[cellI] << "  mS0Tot = " << mS0TotCells[cellI] << nl
                    << "rho1 = " << rho1_cellI << "  rho0 = " << rho0_cellI << nl
                    << "alphaS1 = " << alphaS1Cells[cellI] << "  alphaS0 = " << alphaS0Cells[cellI] << endl;
                print_line(os, 100);
                print_line(os, 100);
            }
        }
        else
        {
            for(i=0; i<n; i++)
            {
                mS1[i].internalField()[cellI] = 0;
                mS0[i].internalField()[cellI] = 0;
            }

            mS1TotCells[cellI] = 0;
            mS0TotCells[cellI] = 0;

            alphaS1Cells[cellI] = 0;
            alphaS0Cells[cellI] = 0;
        }        
    }
    /*
    if(debug)
    {
        print_line(os, 100);
        os<< "Internal faces" << endl;
        print_line(os, 100);
    }
    
    for(faceI=0; faceI<mesh.nInternalFaces(); faceI++)
    {
        faceOwn = own[faceI];
        faceNei = nei[faceI];
        alpha1Own = alpha1Cells[faceOwn];
        alpha1Nei = alpha1Cells[faceNei];        

        if(alpha1Own >= ALPHA_2PH_MAX && alpha1Nei <= ALPHA_2PH_MIN)
        {
            VOwn = V[faceOwn];
            VNei = V[faceNei];
            rho1Own  = rho1Cells[faceOwn];
            rho0Own  = rho0Cells[faceOwn];
            rho1Nei = rho1Cells[faceNei];
            rho0Nei = rho0Cells[faceNei];
            for(i=0; i<n; i++)
            {
                Js1_cellI[i] = Js1[i].internalField()[faceOwn];
                Js0_cellI[i] = Js0[i].internalField()[faceOwn];
                Ys1_cellI[i] = Ys1[i].internalField()[faceOwn];
                Ys0_cellI[i] = Ys0[i].internalField()[faceOwn];
                C1Own[i] = C1[i].internalField()[faceOwn];
                C0Own[i] = C0[i].internalField()[faceOwn];
                C1Nei[i] = C1[i].internalField()[faceNei];
                C0Nei[i] = C0[i].internalField()[faceNei];
                Y1Own[i] = Y1[i].internalField()[faceOwn];
                Y0Own[i] = Y0[i].internalField()[faceOwn];
                Y1Nei[i] = Y1[i].internalField()[faceNei];
                Y0Nei[i] = Y0[i].internalField()[faceNei];
            }

            mS1TotOwn = 0; mS1TotNei = 0;

            mS1Tot_cellI_tmp = Js0_cellI[0] - Js1_cellI[0];
            if(mS1Tot_cellI_tmp > 0)
            {
                if((Ys1_cellI[0] - Y0Nei[0]) > 0)
                {
                    mS1Tot_cellI_tmp /= (Ys1_cellI[0] - Y0Nei[0]);
                }
                else
                {
                    mS1Tot_cellI_tmp /= (Ys1_cellI[0] - Ys0_cellI[0]);
                }

                for(i=0; i<n; i++)
                {                                
                    mS1_cellI[i] = mS1Tot_cellI_tmp*Ys1_cellI[i] + Js1_cellI[i];
                }
            }
            else
            {
                if((Y1Own[0] - Ys0_cellI[0]) > 0)
                {
                    mS1Tot_cellI_tmp /= (Y1Own[0] - Ys0_cellI[0]);
                }
                else
                {
                    mS1Tot_cellI_tmp /= (Ys1_cellI[0] - Ys0_cellI[0]);
                }

                for(i=0; i<n; i++)
                {                                
                    mS1_cellI[i] = mS1Tot_cellI_tmp*Ys0_cellI[i] + Js0_cellI[i];
                }
            }

            for(i=0; i<n; i++)
            {
                C1_cellI[i] = VOwn*C1Own[i] + VNei*C1Nei[i];
                C0_cellI[i] = VOwn*C0Own[i] + VNei*C0Nei[i];
            }

            if(debug)
            {
                print_line(os, 100);
                print_line(os, 100);
                os<< "Face: " << faceI
                    << "  own: " << faceOwn << "  alpha1Own = " << alpha1Own << "  nei: " << faceNei << "  alpha1Nei = " << alpha1Nei << endl;
                print_line(os, 100);
                os<< "mS limiter_min calculation" << endl;
                print_line(os, 100);
            }

            calc_mS_limiter_min(C1_cellI, C0_cellI, Y1Own, Y0Nei, mS1_cellI, mS1Tot_cellI_tmp, dt, n, limiter_min, debug, os);
            
            mS1Tot_cellI = limiter_min*mS1Tot_cellI_tmp;

            if(mS1Tot_cellI_tmp > 0)
            {
                for(i=0; i<n; i++)
                {
                    mS1Own[i] = limiter_min*Js1_cellI[i];
                    mS1Nei[i] = mS1Tot_cellI*Ys1_cellI[i];

                    mS1[i].internalField()[faceOwn] = mS1Own[i]/VOwn;
                    mS0[i].internalField()[faceOwn] = 0;

                    mS1[i].internalField()[faceNei] = mS1Nei[i]/VNei;
                    mS0[i].internalField()[faceNei] = (-mS1Nei[i] - mS1Own[i])/VNei;
                }

                mS1TotOwn = 0;
                mS1TotNei = mS1Tot_cellI/VNei;
            }
            else
            {
                for(i=0; i<n; i++)
                {
                    mS0Nei[i] = -limiter_min*Js0_cellI[i];
                    mS0Own[i] = -mS1Tot_cellI*Ys0_cellI[i];

                    mS1[i].internalField()[faceOwn] = (-mS0Own[i] - mS0Nei[i])/VOwn;
                    mS0[i].internalField()[faceOwn] = mS0Own[i]/VOwn;

                    mS1[i].internalField()[faceNei] = 0;
                    mS0[i].internalField()[faceNei] = mS0Nei[i]/VNei;
                }

                mS1TotOwn = mS1Tot_cellI/VOwn;
                mS1TotNei = 0;
            }

            mS1TotCells[faceOwn] = mS1TotOwn;
            mS0TotCells[faceOwn] = -mS1TotOwn;
            mS1TotCells[faceNei] = mS1TotNei;
            mS0TotCells[faceNei] = -mS1TotNei;

            alphaS1Cells[faceOwn] = mS1TotCells[faceOwn]/rho1Own;
            alphaS0Cells[faceOwn] = mS0TotCells[faceOwn]/rho0Own;
            alphaS1Cells[faceNei] = mS1TotCells[faceNei]/rho1Nei;
            alphaS0Cells[faceNei] = mS0TotCells[faceNei]/rho0Nei;

            if(debug)
            {
                print_line(os, 100);
                print_line(os, 100);
                os<< "Face: " << faceI
                    << "  own: " << faceOwn << "  alpha1Own = " << alpha1Own << "  nei: " << faceNei << "  alpha1Nei = " << alpha1Nei << endl;
                print_line(os, 100);
                os<< "Own" << endl;
                print_line(os, 100);
                os<< setw(8) << "Species" << "  " << setw(14) << "C1" << "  " << setw(14) << "C0" << "  " << setw(14) << "Y1" << "  " << setw(14) << "Y0" << endl;
                print_line(os, 100);
                for(i=0; i<n; i++)
                {
                    os<< setw(8) << i << "  " << setw(14) << C1[i].internalField()[faceOwn] << "  " << setw(14) << C0[i].internalField()[faceOwn] << "  " << setw(14) << Y1[i].internalField()[faceOwn] << "  " << setw(14) << Y0[i].internalField()[faceOwn] << endl;
                }
                print_line(os, 100);
                os<< setw(8) << "Species" << "  " << setw(14) << "Js1" << "  " << setw(14) << "Js0" << "  " << setw(14) << "mS1" << "  " << setw(14) << "mS0" << endl;
                print_line(os, 100);
                for(i=0; i<n; i++)
                {
                    os<< setw(8) << i << "  " << setw(14) << Js1[i].internalField()[faceOwn] << "  " << setw(14) << Js0[i].internalField()[faceOwn] << "  " << setw(14) << mS1[i].internalField()[faceOwn] << "  " << setw(14) << mS0[i].internalField()[faceOwn] << endl;
                }
                print_line(os, 100);
                os<< "limiter_min = " << limiter_min << nl
                    << "mS1Tot = " << mS1TotCells[faceOwn] << "  mS0Tot = " << mS0TotCells[faceOwn] << nl
                    << "rho1 = " << rho1Own << "  rho0 = " << rho0Own << nl
                    << "alphaS1 = " << alphaS1Cells[faceOwn] << "  alphaS0 = " << alphaS0Cells[faceOwn] << endl;
                print_line(os, 100);
                os<< "Nei" << endl;
                print_line(os, 100);
                os<< setw(8) << "Species" << "  " << setw(14) << "C1" << "  " << setw(14) << "C0" << "  " << setw(14) << "Y1" << "  " << setw(14) << "Y0" << endl;
                print_line(os, 100);
                for(i=0; i<n; i++)
                {
                    os<< setw(8) << i << "  " << setw(14) << C1[i].internalField()[faceNei] << "  " << setw(14) << C0[i].internalField()[faceNei] << "  " << setw(14) << Y1[i].internalField()[faceNei] << "  " << setw(14) << Y0[i].internalField()[faceNei] << endl;
                }
                print_line(os, 100);
                os<< setw(8) << "Species" << "  " << setw(14) << "Js1" << "  " << setw(14) << "Js0" << "  " << setw(14) << "mS1" << "  " << setw(14) << "mS0" << endl;
                print_line(os, 100);
                for(i=0; i<n; i++)
                {
                    os<< setw(8) << i << "  " << setw(14) << Js1[i].internalField()[faceNei] << "  " << setw(14) << Js0[i].internalField()[faceNei] << "  " << setw(14) << mS1[i].internalField()[faceNei] << "  " << setw(14) << mS0[i].internalField()[faceNei] << endl;
                }
                print_line(os, 100);
                os<< "limiter_min = " << limiter_min << nl
                    << "mS1Tot = " << mS1TotCells[faceNei] << "  mS0Tot = " << mS0TotCells[faceNei] << nl
                    << "rho1 = " << rho1Nei << "  rho0 = " << rho0Nei << nl
                    << "alphaS1 = " << alphaS1Cells[faceNei] << "  alphaS0 = " << alphaS0Cells[faceNei] << endl;
                print_line(os, 100);
                print_line(os, 100);
            }
        }

        if(alpha1Own <= ALPHA_2PH_MIN && alpha1Nei >= ALPHA_2PH_MAX)
        {
            VOwn = V[faceOwn];
            VNei = V[faceNei];
            rho1Own  = rho1Cells[faceOwn];
            rho0Own  = rho0Cells[faceOwn];
            rho1Nei = rho1Cells[faceNei];
            rho0Nei = rho0Cells[faceNei];
            for(i=0; i<n; i++)
            {
                Js1_cellI[i] = Js1[i].internalField()[faceOwn];
                Js0_cellI[i] = Js0[i].internalField()[faceOwn];
                Ys1_cellI[i] = Ys1[i].internalField()[faceOwn];
                Ys0_cellI[i] = Ys0[i].internalField()[faceOwn];
                C1Own[i] = C1[i].internalField()[faceOwn];
                C0Own[i] = C0[i].internalField()[faceOwn];
                C1Nei[i] = C1[i].internalField()[faceNei];
                C0Nei[i] = C0[i].internalField()[faceNei];
                Y1Own[i] = Y1[i].internalField()[faceOwn];
                Y0Own[i] = Y0[i].internalField()[faceOwn];
                Y1Nei[i] = Y1[i].internalField()[faceNei];
                Y0Nei[i] = Y0[i].internalField()[faceNei];
            }

            mS1TotOwn = 0; mS1TotNei = 0;

            mS1Tot_cellI_tmp = Js0_cellI[0] - Js1_cellI[0];
            if(mS1Tot_cellI_tmp > 0)
            {
                if((Ys1_cellI[0] - Y0Own[0]) > 0)
                {
                    mS1Tot_cellI_tmp /= (Ys1_cellI[0] - Y0Own[0]);
                }
                else
                {
                    mS1Tot_cellI_tmp /= (Ys1_cellI[0] - Ys0_cellI[0]);
                }

                for(i=0; i<n; i++)
                {                                
                    mS1_cellI[i] = mS1Tot_cellI_tmp*Ys1_cellI[i] + Js1_cellI[i];
                }
            }
            else
            {
                if((Y1Nei[0] - Ys0_cellI[0]) > 0)
                {
                    mS1Tot_cellI_tmp /= (Y1Nei[0] - Ys0_cellI[0]);
                }
                else
                {
                    mS1Tot_cellI_tmp /= (Ys1_cellI[0] - Ys0_cellI[0]);
                }

                for(i=0; i<n; i++)
                {                                
                    mS1_cellI[i] = mS1Tot_cellI_tmp*Ys0_cellI[i] + Js0_cellI[i];
                }
            }

            for(i=0; i<n; i++)
            {
                C1_cellI[i] = VOwn*C1Own[i] + VNei*C1Nei[i];
                C0_cellI[i] = VOwn*C0Own[i] + VNei*C0Nei[i];
            }

            if(debug)
            {
                print_line(os, 100);
                print_line(os, 100);
                os<< "Face: " << faceI
                    << "  own: " << faceOwn << "  alpha1Own = " << alpha1Own << "  nei: " << faceNei << "  alpha1Nei = " << alpha1Nei << endl;
                print_line(os, 100);
                os<< "mS limiter_min calculation" << endl;
                print_line(os, 100);
            }

            calc_mS_limiter_min(C1_cellI, C0_cellI, Y1Nei, Y0Own, mS1_cellI, mS1Tot_cellI_tmp, dt, n, limiter_min, debug, os);

            mS1Tot_cellI = limiter_min*mS1Tot_cellI_tmp;

            if(mS1Tot_cellI_tmp > 0)
            {
                for(i=0; i<n; i++)
                {
                    mS1Nei[i] = limiter_min*Js1_cellI[i];
                    mS1Own[i] = mS1Tot_cellI*Ys1_cellI[i];

                    mS1[i].internalField()[faceNei] = mS1Nei[i]/VNei;
                    mS0[i].internalField()[faceNei] = 0;

                    mS1[i].internalField()[faceOwn] = mS1Own[i]/VOwn;
                    mS0[i].internalField()[faceOwn] = (-mS1Nei[i] - mS1Own[i])/VOwn;
                }

                mS1TotNei = 0;
                mS1TotOwn = mS1Tot_cellI/VOwn;
            }
            else
            {
                for(i=0; i<n; i++)
                {
                    mS0Own[i] = -limiter_min*Js0_cellI[i];
                    mS0Nei[i] = -mS1Tot_cellI*Ys0_cellI[i];

                    mS1[i].internalField()[faceNei] = (-mS0Own[i] - mS0Nei[i])/VNei;
                    mS0[i].internalField()[faceNei] = mS0Nei[i]/VNei;

                    mS1[i].internalField()[faceOwn] = 0;
                    mS0[i].internalField()[faceOwn] = mS0Own[i]/VOwn;
                }

                mS1TotNei = mS1Tot_cellI/VNei;
                mS1TotOwn = 0;
            }

            mS1TotCells[faceOwn] = mS1TotOwn;
            mS0TotCells[faceOwn] = -mS1TotOwn;
            mS1TotCells[faceNei] = mS1TotNei;
            mS0TotCells[faceNei] = -mS1TotNei;

            alphaS1Cells[faceOwn] = mS1TotCells[faceOwn]/rho1Own;
            alphaS0Cells[faceOwn] = mS0TotCells[faceOwn]/rho0Own;
            alphaS1Cells[faceNei] = mS1TotCells[faceNei]/rho1Nei;
            alphaS0Cells[faceNei] = mS0TotCells[faceNei]/rho0Nei;

            if(debug)
            {
                print_line(os, 100);
                print_line(os, 100);
                os<< "Face: " << faceI
                    << "  own: " << faceOwn << "  alpha1Own = " << alpha1Own << "  nei: " << faceNei << "  alpha1Nei = " << alpha1Nei << endl;
                print_line(os, 100);
                os<< "Own" << endl;
                print_line(os, 100);
                os<< setw(8) << "Species" << "  " << setw(14) << "C1" << "  " << setw(14) << "C0" << "  " << setw(14) << "Y1" << "  " << setw(14) << "Y0" << endl;
                print_line(os, 100);
                for(i=0; i<n; i++)
                {
                    os<< setw(8) << i << "  " << setw(14) << C1[i].internalField()[faceOwn] << "  " << setw(14) << C0[i].internalField()[faceOwn] << "  " << setw(14) << Y1[i].internalField()[faceOwn] << "  " << setw(14) << Y0[i].internalField()[faceOwn] << endl;
                }
                print_line(os, 100);
                os<< setw(8) << "Species" << "  " << setw(14) << "Js1" << "  " << setw(14) << "Js0" << "  " << setw(14) << "mS1" << "  " << setw(14) << "mS0" << endl;
                print_line(os, 100);
                for(i=0; i<n; i++)
                {
                    os<< setw(8) << i << "  " << setw(14) << Js1[i].internalField()[faceOwn] << "  " << setw(14) << Js0[i].internalField()[faceOwn] << "  " << setw(14) << mS1[i].internalField()[faceOwn] << "  " << setw(14) << mS0[i].internalField()[faceOwn] << endl;
                }
                print_line(os, 100);
                os<< "limiter_min = " << limiter_min << nl
                    << "mS1Tot = " << mS1TotCells[faceOwn] << "  mS0Tot = " << mS0TotCells[faceOwn] << nl
                    << "rho1 = " << rho1Own << "  rho0 = " << rho0Own << nl
                    << "alphaS1 = " << alphaS1Cells[faceOwn] << "  alphaS0 = " << alphaS0Cells[faceOwn] << endl;
                print_line(os, 100);
                os<< "Nei" << endl;
                print_line(os, 100);
                os<< setw(8) << "Species" << "  " << setw(14) << "C1" << "  " << setw(14) << "C0" << "  " << setw(14) << "Y1" << "  " << setw(14) << "Y0" << endl;
                print_line(os, 100);
                for(i=0; i<n; i++)
                {
                    os<< setw(8) << i << "  " << setw(14) << C1[i].internalField()[faceNei] << "  " << setw(14) << C0[i].internalField()[faceNei] << "  " << setw(14) << Y1[i].internalField()[faceNei] << "  " << setw(14) << Y0[i].internalField()[faceNei] << endl;
                }
                print_line(os, 100);
                os<< setw(8) << "Species" << "  " << setw(14) << "Js1" << "  " << setw(14) << "Js0" << "  " << setw(14) << "mS1" << "  " << setw(14) << "mS0" << endl;
                print_line(os, 100);
                for(i=0; i<n; i++)
                {
                    os<< setw(8) << i << "  " << setw(14) << Js1[i].internalField()[faceNei] << "  " << setw(14) << Js0[i].internalField()[faceNei] << "  " << setw(14) << mS1[i].internalField()[faceNei] << "  " << setw(14) << mS0[i].internalField()[faceNei] << endl;
                }
                print_line(os, 100);
                os<< "limiter_min = " << limiter_min << nl
                    << "mS1Tot = " << mS1TotCells[faceNei] << "  mS0Tot = " << mS0TotCells[faceNei] << nl
                    << "rho1 = " << rho1Nei << "  rho0 = " << rho0Nei << nl
                    << "alphaS1 = " << alphaS1Cells[faceNei] << "  alphaS0 = " << alphaS0Cells[faceNei] << endl;
                print_line(os, 100);
                print_line(os, 100);
            }
        }
    }

    //boundary faces
    if(debug)
    {
        print_line(os, 100);
        os<< "Boundary coupled faces" << endl;
        print_line(os, 100);
    }
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const wordList& patchNames = patches.names();

    //--------------------------------------------------------------//
    //Need volume of cell neighbour on neighbouring processor for 
    //each coupled patch face in order to calculate flux limiter_min
    //for that face    
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        if(pp.coupled())
        {
            faceI = pp.start();
            bndFaceI = pp.start() - mesh.nInternalFaces();
            forAll(pp, fcI)
            {
                faceOwn = own[faceI];
                VNeiFld[bndFaceI] = V[faceOwn];
                faceI++;
                bndFaceI++;
            }
        }
    }
    syncTools::swapBoundaryFaceList(mesh, VNeiFld);
    //VNei now has cell volume of neighbouring cell for each coupled
    //patch face
    //--------------------------------------------------------------//

    forAll(mS1Tot.boundaryField(), patchI)
    {
        const polyPatch& pp = patches[patchI];
        const fvPatchScalarField& pmS1Tot = mS1Tot.boundaryField()[patchI];        

        if(pp.coupled())
        {
            if(debug)
            {
                print_line(os, 100);
                os<< "Patch: " << patchNames[patchI] << endl;
                print_line(os, 100);
            }

            faceI = pp.start();
            bndFaceI = pp.start() - mesh.nInternalFaces();
            const scalarField& alpha1NeiFld = alpha1.boundaryField()[patchI].patchNeighbourField();
            const scalarField& rho1NeiFld = rho1.boundaryField()[patchI].patchNeighbourField();
            const scalarField& rho0NeiFld = rho0.boundaryField()[patchI].patchNeighbourField();            

            forAll(pmS1Tot, fcI)
            {
                faceOwn = own[faceI];
                alpha1Own = alpha1Cells[faceOwn];
                alpha1Nei = alpha1NeiFld[fcI];

                if(alpha1Own >= ALPHA_2PH_MAX && alpha1Nei <= ALPHA_2PH_MIN)
                {
                    VOwn = V[faceOwn];
                    VNei = VNeiFld[bndFaceI];
                    rho1Own  = rho1Cells[faceOwn];
                    rho0Own  = rho0Cells[faceOwn];
                    rho1Nei = rho1NeiFld[fcI];
                    rho0Nei = rho0NeiFld[fcI];
                    for(i=0; i<n; i++)
                    {
                        Js1_cellI[i] = Js1[i].internalField()[faceOwn];
                        Js0_cellI[i] = Js0[i].internalField()[faceOwn];
                        Ys1_cellI[i] = Ys1[i].internalField()[faceOwn];
                        Ys0_cellI[i] = Ys0[i].internalField()[faceOwn];
                        C1Own[i] = C1[i].internalField()[faceOwn];
                        C0Own[i] = C0[i].internalField()[faceOwn];
                        const scalarField& C1iNeiFld = C1[i].boundaryField()[patchI].patchNeighbourField();
                        const scalarField& C0iNeiFld = C0[i].boundaryField()[patchI].patchNeighbourField();
                        C1Nei[i] = C1iNeiFld[fcI];
                        C0Nei[i] = C0iNeiFld[fcI];
                        Y1Own[i] = Y1[i].internalField()[faceOwn];
                        Y0Own[i] = Y0[i].internalField()[faceOwn];
                        const scalarField& Y1iNeiFld = Y1[i].boundaryField()[patchI].patchNeighbourField();
                        const scalarField& Y0iNeiFld = Y0[i].boundaryField()[patchI].patchNeighbourField();
                        Y1Nei[i] = Y1iNeiFld[fcI];
                        Y0Nei[i] = Y0iNeiFld[fcI];
                    }

                    mS1TotOwn = 0; mS1TotNei = 0;

                    mS1Tot_cellI_tmp = Js0_cellI[0] - Js1_cellI[0];
                    if(mS1Tot_cellI_tmp > 0)
                    {
                        if((Ys1_cellI[0] - Y0Nei[0]) > 0)
                        {
                            mS1Tot_cellI_tmp /= (Ys1_cellI[0] - Y0Nei[0]);
                        }
                        else
                        {
                            mS1Tot_cellI_tmp /= (Ys1_cellI[0] - Ys0_cellI[0]);
                        }

                        for(i=0; i<n; i++)
                        {                                
                            mS1_cellI[i] = mS1Tot_cellI_tmp*Ys1_cellI[i] + Js1_cellI[i];
                        }
                    }
                    else
                    {
                        if((Y1Own[0] - Ys0_cellI[0]) > 0)
                        {
                            mS1Tot_cellI_tmp /= (Y1Own[0] - Ys0_cellI[0]);
                        }
                        else
                        {
                            mS1Tot_cellI_tmp /= (Ys1_cellI[0] - Ys0_cellI[0]);
                        }

                        for(i=0; i<n; i++)
                        {                                
                            mS1_cellI[i] = mS1Tot_cellI_tmp*Ys0_cellI[i] + Js0_cellI[i];
                        }
                    }

                    for(i=0; i<n; i++)
                    {
                        C1_cellI[i] = VOwn*C1Own[i] + VNei*C1Nei[i];
                        C0_cellI[i] = VOwn*C0Own[i] + VNei*C0Nei[i];
                    }

                    if(debug)
                    {
                        print_line(os, 100);
                        print_line(os, 100);
                        os<< "Face: " << faceI
                            << "  own: " << faceOwn << "  alpha1Own = " << alpha1Own << "  nei: " << -1 << "  alpha1Nei = " << alpha1Nei << endl;
                        print_line(os, 100);
                        os<< "mS limiter_min calculation" << endl;
                        print_line(os, 100);
                    }

                    calc_mS_limiter_min(C1_cellI, C0_cellI, Y1Own, Y0Nei, mS1_cellI, mS1Tot_cellI_tmp, dt, n, limiter_min, debug, os);
            
                    mS1Tot_cellI = limiter_min*mS1Tot_cellI_tmp;

                    if(mS1Tot_cellI_tmp > 0)
                    {
                        for(i=0; i<n; i++)
                        {
                            mS1Own[i] = limiter_min*Js1_cellI[i];
                            mS1Nei[i] = mS1Tot_cellI*Ys1_cellI[i];

                            mS1[i].internalField()[faceOwn] = mS1Own[i]/VOwn;
                            mS0[i].internalField()[faceOwn] = 0;                            
                        }

                        mS1TotOwn = 0;
                        mS1TotNei = mS1Tot_cellI/VNei;
                    }
                    else
                    {
                        for(i=0; i<n; i++)
                        {
                            mS0Nei[i] = -limiter_min*Js0_cellI[i];
                            mS0Own[i] = -mS1Tot_cellI*Ys0_cellI[i];

                            mS1[i].internalField()[faceOwn] = (-mS0Own[i] - mS0Nei[i])/VOwn;
                            mS0[i].internalField()[faceOwn] = mS0Own[i]/VOwn;                            
                        }

                        mS1TotOwn = mS1Tot_cellI/VOwn;
                        mS1TotNei = 0;
                    }

                    mS1TotCells[faceOwn] = mS1TotOwn;
                    mS0TotCells[faceOwn] = -mS1TotOwn;                    

                    alphaS1Cells[faceOwn] = mS1TotCells[faceOwn]/rho1Own;
                    alphaS0Cells[faceOwn] = mS0TotCells[faceOwn]/rho0Own;                    

                    if(debug)
                    {
                        print_line(os, 100);
                        print_line(os, 100);
                        os<< "Face: " << faceI
                            << "  own: " << faceOwn << "  alpha1Own = " << alpha1Own << "  nei: " << -1 << "  alpha1Nei = " << alpha1Nei << endl;
                        print_line(os, 100);
                        os<< "Own" << endl;
                        print_line(os, 100);
                        os<< setw(8) << "Species" << "  " << setw(14) << "C1" << "  " << setw(14) << "C0" << "  " << setw(14) << "Y1" << "  " << setw(14) << "Y0" << endl;
                        print_line(os, 100);
                        for(i=0; i<n; i++)
                        {
                            os<< setw(8) << i << "  " << setw(14) << C1[i].internalField()[faceOwn] << "  " << setw(14) << C0[i].internalField()[faceOwn] << "  " << setw(14) << Y1[i].internalField()[faceOwn] << "  " << setw(14) << Y0[i].internalField()[faceOwn] << endl;
                        }
                        print_line(os, 100);
                        os<< setw(8) << "Species" << "  " << setw(14) << "Js1" << "  " << setw(14) << "Js0" << "  " << setw(14) << "mS1" << "  " << setw(14) << "mS0" << endl;
                        print_line(os, 100);
                        for(i=0; i<n; i++)
                        {
                            os<< setw(8) << i << "  " << setw(14) << Js1[i].internalField()[faceOwn] << "  " << setw(14) << Js0[i].internalField()[faceOwn] << "  " << setw(14) << mS1[i].internalField()[faceOwn] << "  " << setw(14) << mS0[i].internalField()[faceOwn] << endl;
                        }
                        print_line(os, 100);
                        os<< "mS1Tot = " << mS1TotCells[faceOwn] << "  mS0Tot = " << mS0TotCells[faceOwn] << nl
                            << "rho1 = " << rho1Own << "  rho0 = " << rho0Own << nl
                            << "alphaS1 = " << alphaS1Cells[faceOwn] << "  alphaS0 = " << alphaS0Cells[faceOwn] << endl;
                        print_line(os, 100);                        
                        print_line(os, 100);
                    }
                }

                if(alpha1Own <= ALPHA_2PH_MIN && alpha1Nei >= ALPHA_2PH_MAX)
                {
                    VOwn = V[faceOwn];
                    VNei = VNeiFld[bndFaceI];
                    rho1Own  = rho1Cells[faceOwn];
                    rho0Own  = rho0Cells[faceOwn];
                    rho1Nei = rho1NeiFld[fcI];
                    rho0Nei = rho0NeiFld[fcI];
                    for(i=0; i<n; i++)
                    {
                        Js1_cellI[i] = Js1[i].internalField()[faceOwn];
                        Js0_cellI[i] = Js0[i].internalField()[faceOwn];
                        Ys1_cellI[i] = Ys1[i].internalField()[faceOwn];
                        Ys0_cellI[i] = Ys0[i].internalField()[faceOwn];
                        C1Own[i] = C1[i].internalField()[faceOwn];
                        C0Own[i] = C0[i].internalField()[faceOwn];
                        const scalarField& C1iNeiFld = C1[i].boundaryField()[patchI].patchNeighbourField();
                        const scalarField& C0iNeiFld = C0[i].boundaryField()[patchI].patchNeighbourField();
                        C1Nei[i] = C1iNeiFld[fcI];
                        C0Nei[i] = C0iNeiFld[fcI];
                        Y1Own[i] = Y1[i].internalField()[faceOwn];
                        Y0Own[i] = Y0[i].internalField()[faceOwn];
                        const scalarField& Y1iNeiFld = Y1[i].boundaryField()[patchI].patchNeighbourField();
                        const scalarField& Y0iNeiFld = Y0[i].boundaryField()[patchI].patchNeighbourField();
                        Y1Nei[i] = Y1iNeiFld[fcI];
                        Y0Nei[i] = Y0iNeiFld[fcI];
                    }

                    mS1TotOwn = 0; mS1TotNei = 0;

                    mS1Tot_cellI_tmp = Js0_cellI[0] - Js1_cellI[0];
                    if(mS1Tot_cellI_tmp > 0)
                    {
                        if((Ys1_cellI[0] - Y0Own[0]) > 0)
                        {
                            mS1Tot_cellI_tmp /= (Ys1_cellI[0] - Y0Own[0]);
                        }
                        else
                        {
                            mS1Tot_cellI_tmp /= (Ys1_cellI[0] - Ys0_cellI[0]);
                        }

                        for(i=0; i<n; i++)
                        {                                
                            mS1_cellI[i] = mS1Tot_cellI_tmp*Ys1_cellI[i] + Js1_cellI[i];
                        }
                    }
                    else
                    {
                        if((Y1Nei[0] - Ys0_cellI[0]) > 0)
                        {
                            mS1Tot_cellI_tmp /= (Y1Nei[0] - Ys0_cellI[0]);
                        }
                        else
                        {
                            mS1Tot_cellI_tmp /= (Ys1_cellI[0] - Ys0_cellI[0]);
                        }

                        for(i=0; i<n; i++)
                        {                                
                            mS1_cellI[i] = mS1Tot_cellI_tmp*Ys0_cellI[i] + Js0_cellI[i];
                        }
                    }

                    for(i=0; i<n; i++)
                    {
                        C1_cellI[i] = VOwn*C1Own[i] + VNei*C1Nei[i];
                        C0_cellI[i] = VOwn*C0Own[i] + VNei*C0Nei[i];
                    }

                    if(debug)
                    {
                        print_line(os, 100);
                        print_line(os, 100);
                        os<< "Face: " << faceI
                            << "  own: " << faceOwn << "  alpha1Own = " << alpha1Own << "  nei: " << -1 << "  alpha1Nei = " << alpha1Nei << endl;
                        print_line(os, 100);
                        os<< "mS limiter_min calculation" << endl;
                        print_line(os, 100);
                    }

                    calc_mS_limiter_min(C1_cellI, C0_cellI, Y1Nei, Y0Own, mS1_cellI, mS1Tot_cellI_tmp, dt, n, limiter_min, debug, os);

                    mS1Tot_cellI = limiter_min*mS1Tot_cellI_tmp;

                    if(mS1Tot_cellI_tmp > 0)
                    {
                        for(i=0; i<n; i++)
                        {
                            mS1Nei[i] = limiter_min*Js1_cellI[i];
                            mS1Own[i] = mS1Tot_cellI*Ys1_cellI[i];                  

                            mS1[i].internalField()[faceOwn] = mS1Own[i]/VOwn;
                            mS0[i].internalField()[faceOwn] = (-mS1Nei[i] - mS1Own[i])/VOwn;
                        }

                        mS1TotNei = 0;
                        mS1TotOwn = mS1Tot_cellI/VOwn;
                    }
                    else
                    {
                        for(i=0; i<n; i++)
                        {
                            mS0Own[i] = -limiter_min*Js0_cellI[i];
                            mS0Nei[i] = -mS1Tot_cellI*Ys0_cellI[i];                            

                            mS1[i].internalField()[faceOwn] = 0;
                            mS0[i].internalField()[faceOwn] = mS0Own[i]/VOwn;
                        }

                        mS1TotNei = mS1Tot_cellI/VNei;
                        mS1TotOwn = 0;
                    }

                    mS1TotCells[faceOwn] = mS1TotOwn;
                    mS0TotCells[faceOwn] = -mS1TotOwn;                    

                    alphaS1Cells[faceOwn] = mS1TotCells[faceOwn]/rho1Own;
                    alphaS0Cells[faceOwn] = mS0TotCells[faceOwn]/rho0Own;                    

                    if(debug)
                    {
                        print_line(os, 100);
                        print_line(os, 100);
                        os<< "Face: " << faceI
                            << "  own: " << faceOwn << "  alpha1Own = " << alpha1Own << "  nei: " << -1 << "  alpha1Nei = " << alpha1Nei << endl;
                        print_line(os, 100);
                        os<< "Own" << endl;
                        print_line(os, 100);
                        os<< setw(8) << "Species" << "  " << setw(14) << "C1" << "  " << setw(14) << "C0" << "  " << setw(14) << "Y1" << "  " << setw(14) << "Y0" << endl;
                        print_line(os, 100);
                        for(i=0; i<n; i++)
                        {
                            os<< setw(8) << i << "  " << setw(14) << C1[i].internalField()[faceOwn] << "  " << setw(14) << C0[i].internalField()[faceOwn] << "  " << setw(14) << Y1[i].internalField()[faceOwn] << "  " << setw(14) << Y0[i].internalField()[faceOwn] << endl;
                        }
                        print_line(os, 100);
                        os<< setw(8) << "Species" << "  " << setw(14) << "Js1" << "  " << setw(14) << "Js0" << "  " << setw(14) << "mS1" << "  " << setw(14) << "mS0" << endl;
                        print_line(os, 100);
                        for(i=0; i<n; i++)
                        {
                            os<< setw(8) << i << "  " << setw(14) << Js1[i].internalField()[faceOwn] << "  " << setw(14) << Js0[i].internalField()[faceOwn] << "  " << setw(14) << mS1[i].internalField()[faceOwn] << "  " << setw(14) << mS0[i].internalField()[faceOwn] << endl;
                        }
                        print_line(os, 100);
                        os<< "mS1Tot = " << mS1TotCells[faceOwn] << "  mS0Tot = " << mS0TotCells[faceOwn] << nl
                            << "rho1 = " << rho1Own << "  rho0 = " << rho0Own << nl
                            << "alphaS1 = " << alphaS1Cells[faceOwn] << "  alphaS0 = " << alphaS0Cells[faceOwn] << endl;
                        print_line(os, 100);                        
                        print_line(os, 100);
                    }
                }
                
                faceI++;
                bndFaceI++;
            }
        }
    }
        */
    if(debug)
    {
        print_line(os, 100);
        print_line(os, 100);
        os<< "Done Interfacial Species Transfer Source Terms Calculation" << endl;
        print_line(os, 100);
        print_line(os, 100);
    }
}


void calc_mS_He
(
    const fvMesh& mesh,
    const labelListList& cellStencil,    
    const List<List<scalar> >& Y1_flatFld,
    const List<List<scalar> >& Y0_flatFld,
    const List<scalar>& alpha1_flatFld,
    const PtrList<volScalarField>& C1,
    const PtrList<volScalarField>& C0,
    const PtrList<volScalarField>& Y1,
    const PtrList<volScalarField>& Y0,
    const volScalarField& alpha1,
    const volScalarField& rho1,
    const volScalarField& rho0,
    const PtrList<volScalarField>& D1,
    const PtrList<volScalarField>& D0,
    const List<vector>& C_ph1_flatFld,
    const List<vector>& C_ph0_flatFld,
    const volVectorField& C_intfc,
    const volScalarField& A_intfc,
    const volVectorField& nHat,
    const List<scalar>& He,
    PtrList<volScalarField>& Ys1,
    PtrList<volScalarField>& Ys0,
    PtrList<volScalarField>& Js1,
    PtrList<volScalarField>& Js0,
    const label& nSpecies,
    const scalar& ALPHA_2PH_MIN,
    const scalar& A_INTFC_2PH_MIN,
    PtrList<volScalarField>& mS1,
    PtrList<volScalarField>& mS0,
    const bool debug,
    OFstream& os
)
{
    scalar ALPHA_2PH_MAX = 1 - ALPHA_2PH_MIN;
    scalar dt = mesh.time().deltaTValue();
    const scalarField& V = mesh.V();

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();

    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();
    const surfaceVectorField& Cf = mesh.Cf();

    const scalarField& alpha1Cells = alpha1.internalField();    
    const scalarField& rho1Cells = rho1.internalField();
    const scalarField& rho0Cells = rho0.internalField();
    const vectorField& C_intfcCells = C_intfc.internalField();
    const scalarField& A_intfcCells = A_intfc.internalField();
    const vectorField& nHatCells = nHat.internalField();

    if(debug)
    {
        os<< "-------------------------------------------------------------------------" << nl
          << "Interfacial Species Transfer Source Terms Calculation (Henry's Law Model)" << nl
          << "-------------------------------------------------------------------------" << nl
          << nl
          << "-------------------------------------------------------------------------" << nl
          << "Internal cells" << nl
          << "-------------------------------------------------------------------------" << nl
          << endl;
    }

    //Js for all interface cells
    forAll(alpha1Cells, cellI)
    {
        scalar alpha1_cellI = alpha1Cells[cellI];        
        scalar rho1_cellI = rho1Cells[cellI];        
        scalar rho0_cellI = rho0Cells[cellI];        
        scalar V_cellI = V[cellI];
        if(debug)
        {
            os<< "Cell: " << cellI << "  alpha1 = " << alpha1_cellI << nl
              << "-------------------------------------------------------------------------" << endl;
        }

        scalar mS1i_cellI = 0;
        //scalar mS0Tot_cellI = 0;
        //List<scalar> mS1_cellI(nSpecies);
        //List<scalar> mS0_cellI(nSpecies);
        //scalar limiter = 1;
        //scalar limiter_min = 1;
        scalar max_mSi;

        if(alpha1_cellI > ALPHA_2PH_MIN && alpha1_cellI < ALPHA_2PH_MAX && A_intfcCells[cellI] > A_INTFC_2PH_MIN)
        {                         
            vector nf = nHatCells[cellI];
            vector C_intfc_cellI = C_intfcCells[cellI];
            scalar A_intfc_cellI = A_intfcCells[cellI];
            labelList curCellsAll = cellStencil[cellI];

            if(debug)
            {
                os<< "nf = " << nf << "  C_intfc = " << C_intfc_cellI << "  A_intfc = " << A_intfc_cellI << nl
                    << "C_ph1 = " << C_ph1_flatFld[cellI] << "  C_ph0 = " << C_ph0_flatFld[cellI] << nl
                    << "rho1 = " << rho1_cellI << "  rho0 = " << rho0_cellI << nl << endl;
            }

            scalar dn1;
            List<scalar> Yeff1(nSpecies);
            scalar dn0;
            List<scalar> Yeff0(nSpecies);
            scalar intfcGradi_cellI;

            //phase-1
            //ensure nf direction is into the phase
            calc_cell_intfcGrad_coeffs(mesh, cellI, nf, C_intfc_cellI, Y1_flatFld, alpha1_flatFld, C_ph1_flatFld, curCellsAll, nSpecies, ALPHA_2PH_MIN, 1, dn1, Yeff1, debug, os);

            if(debug)
            {
                os<< "dn1 = " << dn1 << endl;
            }
            if(dn1 < SMALL)
            {
                dn1 += SMALL;
            }
            if(debug)
            {
                os<< "dn1 stab = " << dn1 << endl;
            }

            //phase-0
            //ensure nf direction is into the phase
            //then reverse nf again for Js0 calculation            
            calc_cell_intfcGrad_coeffs(mesh, cellI, -nf, C_intfc_cellI, Y0_flatFld, alpha1_flatFld, C_ph0_flatFld, curCellsAll, nSpecies, ALPHA_2PH_MIN, 0, dn0, Yeff0, debug, os);
            if(debug)
            {
                os<< "dn0 = " << dn0 << endl;
            }
            if(dn0 < SMALL)
            {
                dn0 += SMALL;
            }
            if(debug)
            {
                os<< "dn0 stab = " << dn0 << endl;
            }            

            scalar mS1Tot = 0;
            scalar mS0Tot = 0;

            for(label i=0; i<(nSpecies-1); i++)
            {     
                Ys0[i].internalField()[cellI] = rho0_cellI*D0[i].internalField()[cellI]*Yeff0[i]*dn1 + rho1_cellI*D1[i].internalField()[cellI]*Yeff1[i]*dn0;
                Ys0[i].internalField()[cellI] /= (rho0_cellI*D0[i].internalField()[cellI]*dn1 + He[i]*rho1_cellI*D1[i].internalField()[cellI]*dn0);

                Ys1[i].internalField()[cellI] = He[i]*Ys0[i].internalField()[cellI];

                intfcGradi_cellI = (Yeff1[i] - Ys1[i].internalField()[cellI])/dn1;
                Js1[i].internalField()[cellI] = -A_intfc_cellI*rho1_cellI*D1[i].internalField()[cellI]*intfcGradi_cellI;

                Js0[i].internalField()[cellI] = Js1[i].internalField()[cellI];

                mS1i_cellI = Js1[i].internalField()[cellI]/V_cellI;
                
                if(mS1i_cellI > 0)
                {
                    max_mSi = C0[i].internalField()[cellI]/dt;

                    mS1[i].internalField()[cellI] = min(mS1i_cellI, max_mSi);                                        
                }
                else if(mS1i_cellI < 0)
                {
                    max_mSi = C1[i].internalField()[cellI]/dt;

                    mS1[i].internalField()[cellI] = -min(-mS1i_cellI, max_mSi);
                }
                else
                {
                    mS1[i].internalField()[cellI] = 0;
                }

                mS0[i].internalField()[cellI] = -mS1[i].internalField()[cellI];

                mS1Tot += mS1[i].internalField()[cellI];
                mS0Tot += mS0[i].internalField()[cellI];
            }            

            mS1[nSpecies-1].internalField()[cellI] = -mS1Tot;
            mS0[nSpecies-1].internalField()[cellI] = -mS0Tot;
        }
        else
        {
            for(label i=0; i<nSpecies; i++)
            {
                mS1[i].internalField()[cellI] = 0;
                mS0[i].internalField()[cellI] = 0;
            }
        }

        if(debug)
        {
            os<< "-----------------------------------------------------------------------------------" << endl;
            for(label i=0; i<nSpecies; i++)
            {
                os<< "C1[" << i << "] = " << C1[i].internalField()[cellI] << "  C0[" << i << "] = " << C0[i].internalField()[cellI] << nl
                    << "Ys1[" << i << "] = " << Ys1[i].internalField()[cellI] << "  Ys0[" << i << "] = " << Ys0[i].internalField()[cellI] << nl 
                    << "Js1[" << i << "] = " << Js1[i].internalField()[cellI] << "  Js0[" << i << "] = " << Js0[i].internalField()[cellI] << nl 
                    << "mS1[" << i << "] = " << mS1[i].internalField()[cellI] << "  mS0[" << i << "] = " << mS0[i].internalField()[cellI] 
                    << endl;
            }
            os<< "rho1 = " << rho1_cellI << "  rho0 = " << rho0_cellI << nl                
                << "-----------------------------------------------------------------------------------" << nl
                << endl;
        }        
    }

    if(debug)
    {
        os<< "-------------------------------------------------------------------------" << nl
            << "Internal faces" << nl
            << "-------------------------------------------------------------------------" << nl
            << endl;
    }

    for(label faceI=0; faceI<mesh.nInternalFaces(); faceI++)
    {
        label faceOwn = own[faceI];
        label faceNei = nei[faceI];
        scalar alpha1Own = alpha1Cells[faceOwn];
        scalar alpha1Nei = alpha1Cells[faceNei];
        scalar VOwn = V[faceOwn];
        scalar VNei = V[faceNei];
        scalar rho1Own  = rho1Cells[faceOwn];
        scalar rho0Own  = rho0Cells[faceOwn];
        scalar rho1Nei = rho1Cells[faceNei];
        scalar rho0Nei = rho0Cells[faceNei];

        if(debug)
        {
            os<< "Face: " << faceI << nl
                << "own: " << faceOwn << "  alpha1Own = " << alpha1Own << "  nei: " << faceNei << "  alpha1Nei = " << alpha1Nei
                << endl;
        }

        scalar mS1Tot_own = 0;
        scalar mS1Tot_nei = 0;
        scalar mS0Tot_own = 0;
        scalar mS0Tot_nei = 0;
        scalar mS1i_cellI = 0;
        scalar max_mSi = 0;

        if(alpha1Own >= ALPHA_2PH_MAX && alpha1Nei <= ALPHA_2PH_MIN)
        {
            scalar A_intfc_cellI = magSf[faceI];
            vector nf = -Sf[faceI]/A_intfc_cellI;
            vector C_intfc_cellI = Cf[faceI];
            labelList curCellsAll = cellStencil[faceOwn];

            if(debug)
            {
                os<< "ph-1 cell: " << faceOwn << "  ph-0 cell: " << faceNei << nl
                    << "nfOwn = " << nHatCells[faceOwn] << "  nfNei = " << nHatCells[faceNei] << nl
                    << "C_intfcOwn = " << C_intfcCells[faceOwn] << "  C_intfcNei = " << C_intfcCells[faceNei] << nl
                    << "A_intfcOwn = " << A_intfcCells[faceOwn] << "  A_intfcNei = " << A_intfcCells[faceNei] << nl
                    << "nf = " << nf << "  C_intfc = " << C_intfc_cellI << "  A_intfc = " << A_intfc_cellI << nl
                    << "C_ph1Own = " << C_ph1_flatFld[faceOwn] << "  C_ph0Own = " << C_ph0_flatFld[faceOwn] << nl
                    << "rho1Own = " << rho1Cells[faceOwn] << "  rho0Own = " << rho0Cells[faceOwn] << nl
                    << "C_ph1Nei = " << C_ph1_flatFld[faceNei] << "  C_ph0Nei = " << C_ph0_flatFld[faceNei] << nl
                    << "rho1Nei = " << rho1Cells[faceNei] << "  rho0Nei = " << rho0Cells[faceNei] << endl;
            }

            scalar dn1;
            List<scalar> Yeff1(nSpecies);
            scalar dn0;
            List<scalar> Yeff0(nSpecies);
            scalar intfcGradi_cellI;

            //phase-1
            //ensure nf direction is into the phase
            calc_cell_intfcGrad_coeffs(mesh, faceOwn, nf, C_intfc_cellI, Y1_flatFld, alpha1_flatFld, C_ph1_flatFld, curCellsAll, nSpecies, ALPHA_2PH_MIN, 1, dn1, Yeff1, debug, os);

            if(debug)
            {
                os<< "dn1 = " << dn1 << endl;
            }
            if(dn1 < SMALL)
            {
                dn1 += SMALL;
            }
            if(debug)
            {
                os<< "dn1 stab = " << dn1 << endl;
            }

            //phase-0
            //ensure nf direction is into the phase
            //then reverse nf again for Js0 calculation
            curCellsAll = cellStencil[faceNei];
            calc_cell_intfcGrad_coeffs(mesh, faceNei, -nf, C_intfc_cellI, Y0_flatFld, alpha1_flatFld, C_ph0_flatFld, curCellsAll, nSpecies, ALPHA_2PH_MIN, 0, dn0, Yeff0, debug, os);
            if(debug)
            {
                os<< "dn0 = " << dn0 << endl;
            }
            if(dn0 < SMALL)
            {
                dn0 += SMALL;
            }
            if(debug)
            {
                os<< "dn0 stab = " << dn0 << endl;
            }

            for(label i=0; i<(nSpecies-1); i++)
            {     
                Ys0[i].internalField()[faceOwn] = rho0Nei*D0[i].internalField()[faceNei]*Yeff0[i]*dn1 + rho1Own*D1[i].internalField()[faceOwn]*Yeff1[i]*dn0;
                Ys0[i].internalField()[faceOwn] /= (rho0Nei*D0[i].internalField()[faceNei]*dn1 + He[i]*rho1Own*D1[i].internalField()[faceOwn]*dn0);                

                Ys1[i].internalField()[faceOwn] = He[i]*Ys0[i].internalField()[faceOwn];

                Ys0[i].internalField()[faceNei] = Ys0[i].internalField()[faceOwn];
                Ys1[i].internalField()[faceNei] = Ys1[i].internalField()[faceOwn];

                intfcGradi_cellI = (Yeff1[i] - Ys1[i].internalField()[faceOwn])/dn1;
                Js1[i].internalField()[faceOwn] = -A_intfc_cellI*rho1Own*D1[i].internalField()[faceOwn]*intfcGradi_cellI;

                Js0[i].internalField()[faceOwn] = Js1[i].internalField()[faceOwn];

                Js0[i].internalField()[faceNei] = Js0[i].internalField()[faceOwn];
                Js1[i].internalField()[faceNei] = Js1[i].internalField()[faceOwn];

                mS1i_cellI = Js1[i].internalField()[faceOwn];
                
                if(mS1i_cellI > 0)
                {
                    max_mSi = C0[i].internalField()[faceNei]*VNei/dt;

                    mS1[i].internalField()[faceOwn] = min(mS1i_cellI, max_mSi);                    
                }
                else if(mS1i_cellI < 0)
                {
                    max_mSi = C1[i].internalField()[faceOwn]*VOwn/dt;

                    mS1[i].internalField()[faceOwn] = -min(-mS1i_cellI, max_mSi);
                }
                else
                {
                    mS1[i].internalField()[faceOwn] = 0;
                }

                mS0[i].internalField()[faceNei] = -mS1[i].internalField()[faceOwn];
                mS1[i].internalField()[faceOwn] /= VOwn;
                mS0[i].internalField()[faceNei] /= VNei;
                mS0[i].internalField()[faceOwn] = 0;
                mS1[i].internalField()[faceNei] = 0;

                mS1Tot_own += mS1[i].internalField()[faceOwn];
                mS0Tot_own += mS0[i].internalField()[faceOwn];

                mS1Tot_nei += mS1[i].internalField()[faceNei];
                mS0Tot_nei += mS0[i].internalField()[faceNei];
            }            

            mS1[nSpecies-1].internalField()[faceOwn] = -mS1Tot_own;
            mS0[nSpecies-1].internalField()[faceOwn] = -mS0Tot_own;	   

            mS1[nSpecies-1].internalField()[faceNei] = -mS1Tot_nei;
            mS0[nSpecies-1].internalField()[faceNei] = -mS0Tot_nei;	   
        }

        if(alpha1Own <= ALPHA_2PH_MIN && alpha1Nei >= ALPHA_2PH_MAX)
        {            
            scalar A_intfc_cellI = magSf[faceI];
            vector nf = Sf[faceI]/A_intfc_cellI;
            vector C_intfc_cellI = Cf[faceI];
            labelList curCellsAll = cellStencil[faceNei];

            if(debug)
            {
                os<< "ph-1 cell: " << faceNei << "  ph-0 cell: " << faceOwn << nl
                    << "nfOwn = " << nHatCells[faceOwn] << "  nfNei = " << nHatCells[faceNei] << nl
                    << "C_intfcOwn = " << C_intfcCells[faceOwn] << "  C_intfcNei = " << C_intfcCells[faceNei] << nl
                    << "A_intfcOwn = " << A_intfcCells[faceOwn] << "  A_intfcNei = " << A_intfcCells[faceNei] << nl
                    << "nf = " << nf << "  C_intfc = " << C_intfc_cellI << "  A_intfc = " << A_intfc_cellI << nl
                    << "C_ph1Own = " << C_ph1_flatFld[faceOwn] << "  C_ph0Own = " << C_ph0_flatFld[faceOwn] << nl
                    << "rho1Own = " << rho1Cells[faceOwn] << "  rho0Own = " << rho0Cells[faceOwn] << nl
                    << "C_ph1Nei = " << C_ph1_flatFld[faceNei] << "  C_ph0Nei = " << C_ph0_flatFld[faceNei] << nl
                    << "rho1Nei = " << rho1Cells[faceNei] << "  rho0Nei = " << rho0Cells[faceNei] << endl;
            }

            scalar dn1;
            List<scalar> Yeff1(nSpecies);
            scalar dn0;
            List<scalar> Yeff0(nSpecies);
            scalar intfcGradi_cellI;

            //phase-1
            //ensure nf direction is into the phase
            calc_cell_intfcGrad_coeffs(mesh, faceNei, nf, C_intfc_cellI, Y1_flatFld, alpha1_flatFld, C_ph1_flatFld, curCellsAll, nSpecies, ALPHA_2PH_MIN, 1, dn1, Yeff1, debug, os);

            if(debug)
            {
                os<< "dn1 = " << dn1 << endl;
            }
            if(dn1 < SMALL)
            {
                dn1 += SMALL;
            }
            if(debug)
            {
                os<< "dn1 stab = " << dn1 << endl;
            }

            //phase-0
            //ensure nf direction is into the phase
            //then reverse nf again for Js0 calculation
            curCellsAll = cellStencil[faceOwn];
            calc_cell_intfcGrad_coeffs(mesh, faceOwn, -nf, C_intfc_cellI, Y0_flatFld, alpha1_flatFld, C_ph0_flatFld, curCellsAll, nSpecies, ALPHA_2PH_MIN, 0, dn0, Yeff0, debug, os);
            if(debug)
            {
                os<< "dn0 = " << dn0 << endl;
            }
            if(dn0 < SMALL)
            {
                dn0 += SMALL;
            }
            if(debug)
            {
                os<< "dn0 stab = " << dn0 << endl;
            }

            for(label i=0; i<(nSpecies-1); i++)
            {     
                Ys0[i].internalField()[faceNei] = rho0Own*D0[i].internalField()[faceOwn]*Yeff0[i]*dn1 + rho1Nei*D1[i].internalField()[faceNei]*Yeff1[i]*dn0;
                Ys0[i].internalField()[faceNei] /= (rho0Own*D0[i].internalField()[faceOwn]*dn1 + He[i]*rho1Nei*D1[i].internalField()[faceNei]*dn0);                

                Ys1[i].internalField()[faceNei] = He[i]*Ys0[i].internalField()[faceNei];

                Ys0[i].internalField()[faceOwn] = Ys0[i].internalField()[faceNei];
                Ys1[i].internalField()[faceOwn] = Ys1[i].internalField()[faceNei];

                intfcGradi_cellI = (Yeff1[i] - Ys1[i].internalField()[faceNei])/dn1;
                Js1[i].internalField()[faceNei] = -A_intfc_cellI*rho1Nei*D1[i].internalField()[faceNei]*intfcGradi_cellI;

                Js0[i].internalField()[faceNei] = Js1[i].internalField()[faceNei];

                Js0[i].internalField()[faceOwn] = Js0[i].internalField()[faceNei];
                Js1[i].internalField()[faceOwn] = Js1[i].internalField()[faceNei];

                mS1i_cellI = Js1[i].internalField()[faceNei];
                
                if(mS1i_cellI > 0)
                {
                    max_mSi = C0[i].internalField()[faceOwn]*VOwn/dt;

                    mS1[i].internalField()[faceNei] = min(mS1i_cellI, max_mSi);                    
                }
                else if(mS1i_cellI < 0)
                {
                    max_mSi = C1[i].internalField()[faceNei]*VNei/dt;

                    mS1[i].internalField()[faceNei] = -min(-mS1i_cellI, max_mSi);
                }
                else
                {
                    mS1[i].internalField()[faceNei] = 0;
                }

                mS0[i].internalField()[faceOwn] = -mS1[i].internalField()[faceNei];
                mS1[i].internalField()[faceNei] /= VNei;
                mS0[i].internalField()[faceOwn] /= VOwn;
                mS0[i].internalField()[faceNei] = 0;
                mS1[i].internalField()[faceOwn] = 0;

                mS1Tot_nei += mS1[i].internalField()[faceNei];
                mS0Tot_nei += mS0[i].internalField()[faceNei];

                mS1Tot_own += mS1[i].internalField()[faceOwn];
                mS0Tot_own += mS0[i].internalField()[faceOwn];
            }            

            mS1[nSpecies-1].internalField()[faceNei] = -mS1Tot_nei;
            mS0[nSpecies-1].internalField()[faceNei] = -mS0Tot_nei;	   

            mS1[nSpecies-1].internalField()[faceOwn] = -mS1Tot_own;
            mS0[nSpecies-1].internalField()[faceOwn] = -mS0Tot_own;	   
        }

        if(debug)
        {
            os<< "-----------------------------------------------------------------------------------" << endl;
            for(label i=0; i<nSpecies; i++)
            {
                os<< "face owner: " << nl
                    << "C1[" << i << "] = " << C1[i].internalField()[faceOwn] << "  C0[" << i << "] = " << C0[i].internalField()[faceOwn] << nl
                    << "Js1[" << i << "] = " << Js1[i].internalField()[faceOwn] << "  Js0[" << i << "] = " << Js0[i].internalField()[faceOwn] << nl 
                    << "mS1[" << i << "] = " << mS1[i].internalField()[faceOwn] << "  mS0[" << i << "] = " << mS0[i].internalField()[faceOwn] 
                    << endl;
            }
            for(label i=0; i<nSpecies; i++)
            {
                os<< "face neighbour: " << nl
                    << "C1[" << i << "] = " << C1[i].internalField()[faceNei] << "  C0[" << i << "] = " << C0[i].internalField()[faceNei] << nl
                    << "Js1[" << i << "] = " << Js1[i].internalField()[faceNei] << "  Js0[" << i << "] = " << Js0[i].internalField()[faceNei] << nl
                    << "mS1[" << i << "] = " << mS1[i].internalField()[faceNei] << "  mS0[" << i << "] = " << mS0[i].internalField()[faceNei] 
                    << endl;
            }            
            os<< "-----------------------------------------------------------------------------------" << nl
                << endl;
        }
    }

    //boundary faces
    if(debug)
    {
        os<< "-------------------------------------------------------------------------" << nl
          << "Boundary coupled faces step 1" << nl
          << "-------------------------------------------------------------------------" << nl
          << endl;
    }
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const wordList& patchNames = patches.names();
    const label nBnd = mesh.nFaces() - mesh.nInternalFaces();

    List<List<scalar> > YeffOwn(nSpecies);
    List<List<scalar> > YeffNei(nSpecies);
    List<scalar> dnOwn(nBnd);
    List<scalar> dnNei(nBnd);
    List<scalar> VNeiFld(nBnd);
    for(label i=0; i<nSpecies; i++)
    {
        YeffOwn[i].setSize(nBnd);
        YeffNei[i].setSize(nBnd);
    }

    forAll(Js1[0].boundaryField(), patchI)
    {
        const polyPatch& pp = patches[patchI];
        
        if(debug)
        {
            os<< "-------------------------------------------------------------------------" << nl
                << "Patch: " << patchNames[patchI] << nl
                << "-------------------------------------------------------------------------" << nl
                << endl;
        }

        if(pp.coupled())
        {
            label faceI = pp.start();
            label bndFaceI = faceI - mesh.nInternalFaces();

            const scalarField& alpha1NeiFld = alpha1.boundaryField()[patchI].patchNeighbourField();
            const vectorField& nHatNeiFld = nHat.boundaryField()[patchI].patchNeighbourField();
            const vectorField& C_intfcNeiFld = C_intfc.boundaryField()[patchI].patchNeighbourField();
            const scalarField& A_intfcNeiFld = A_intfc.boundaryField()[patchI].patchNeighbourField();
            const scalarField& rho1NeiFld = rho1.boundaryField()[patchI].patchNeighbourField();            
            const scalarField& rho0NeiFld = rho0.boundaryField()[patchI].patchNeighbourField();
            const fvsPatchVectorField& pSf = Sf.boundaryField()[patchI];
            const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchI];
            const fvsPatchVectorField& pCf = Cf.boundaryField()[patchI];

            forAll(Js1[0].boundaryField()[patchI], fcI)
            {
                label faceOwn = own[faceI];
                scalar alpha1Own = alpha1Cells[faceOwn];
                scalar alpha1Nei = alpha1NeiFld[fcI];

                if(debug)
                {
                    os<< "Face: " << faceI << "  patch face index: " << fcI << "  bnd face index: " << bndFaceI << nl
                        << "own: " << faceOwn << "  alpha1Own = " << alpha1Own << "  alpha1Nei = " << alpha1Nei << endl;
                }

                if(alpha1Own >= ALPHA_2PH_MAX && alpha1Nei <= ALPHA_2PH_MIN)
                {                                
                    scalar A_intfc_cellI = pMagSf[fcI];
                    vector nf = -pSf[fcI]/pMagSf[fcI];
                    vector C_intfc_cellI = pCf[fcI];
                    const labelList& curCellsAll = cellStencil[faceOwn];

                    if(debug)
                    {
                        os<< "ph-1 cell: " << faceOwn << nl
                            << "nfOwn = " << nHatCells[faceOwn] << "  nfNei = " << nHatNeiFld[fcI] << nl
                            << "C_intfcOwn = " << C_intfcCells[faceOwn] << "  C_intfcNei = " << C_intfcNeiFld[fcI] << nl
                            << "A_intfcOwn = " << A_intfcCells[faceOwn] << "  A_intfcNei = " << A_intfcNeiFld[fcI] << nl
                            << "nf = " << nf << "  C_intfc = " << C_intfc_cellI << "  A_intfc = " << A_intfc_cellI << nl
                            << "C_ph1Own = " << C_ph1_flatFld[faceOwn] << "  C_ph0Own = " << C_ph0_flatFld[faceOwn] << nl
                            << "rho1Own = " << rho1Cells[faceOwn] << "  rho0Own = " << rho0Cells[faceOwn] << nl                            
                            << "rho1Nei = " << rho1NeiFld[fcI] << "  rho0Nei = " << rho0NeiFld[fcI] << endl;
                    }

                    scalar dn1;
                    List<scalar> Yeff1(nSpecies);                    

                    //phase-1
                    //ensure nf direction is into the phase
                    calc_cell_intfcGrad_coeffs(mesh, faceOwn, nf, C_intfc_cellI, Y1_flatFld, alpha1_flatFld, C_ph1_flatFld, curCellsAll, nSpecies, ALPHA_2PH_MIN, 1, dn1, Yeff1, debug, os);
                    if(debug)
                    {
                        os<< "dn1 = " << dn1 << endl;
                    }
                    if(dn1 < SMALL)
                    {
                        dn1 += SMALL;
                    }
                    if(debug)
                    {
                        os<< "dn1 stab = " << dn1 << endl;
                    }
             
                    for(label i=0; i<nSpecies; i++)
                    {
                        YeffOwn[i][bndFaceI] = Yeff1[i];
                        YeffNei[i][bndFaceI] = Yeff1[i];
                        dnOwn[bndFaceI] = dn1;
                        dnNei[bndFaceI] = dn1;
                    }
                }// end if ownCell is ph-1: if(alpha1Own >= ALPHA_2PH_MAX && alpha1Nei <= ALPHA_2PH_MIN)

                if(alpha1Own <= ALPHA_2PH_MIN && alpha1Nei >= ALPHA_2PH_MAX)
                {
                    scalar A_intfc_cellI = pMagSf[fcI];
                    vector nf = pSf[fcI]/pMagSf[fcI];
                    vector C_intfc_cellI = pCf[fcI];
                    const labelList& curCellsAll = cellStencil[faceOwn];

                    if(debug)
                    {
                        os<< "ph-0 cell: " << faceOwn << nl
                            << "nfOwn = " << nHatCells[faceOwn] << "  nfNei = " << nHatNeiFld[fcI] << nl
                            << "C_intfcOwn = " << C_intfcCells[faceOwn] << "  C_intfcNei = " << C_intfcNeiFld[fcI] << nl
                            << "A_intfcOwn = " << A_intfcCells[faceOwn] << "  A_intfcNei = " << A_intfcNeiFld[fcI] << nl
                            << "nf = " << nf << "  C_intfc = " << C_intfc_cellI << "  A_intfc = " << A_intfc_cellI << nl
                            << "C_ph1Own = " << C_ph1_flatFld[faceOwn] << "  C_ph0Own = " << C_ph0_flatFld[faceOwn] << nl
                            << "rho1Own = " << rho1Cells[faceOwn] << "  rho0Own = " << rho0Cells[faceOwn] << nl                            
                            << "rho1Nei = " << rho1NeiFld[fcI] << "  rho0Nei = " << rho0NeiFld[fcI] << endl;
                    }

                    scalar dn0;
                    List<scalar> Yeff0(nSpecies);                    

                    //phase-0
                    //ensure nf direction is into the phase
                    //then reverse nf again for Js0 calculation
                    calc_cell_intfcGrad_coeffs(mesh, faceOwn, -nf, C_intfc_cellI, Y0_flatFld, alpha1_flatFld, C_ph0_flatFld, curCellsAll, nSpecies, ALPHA_2PH_MIN, 0, dn0, Yeff0, debug, os);
                    if(debug)
                    {
                        os<< "dn0 = " << dn0 << endl;
                    }
                    if(dn0 < SMALL)
                    {
                        dn0 += SMALL;
                    }
                    if(debug)
                    {
                        os<< "dn0 stab = " << dn0 << endl;
                    }
                    
                    for(label i=0; i<nSpecies; i++)
                    {
                        YeffOwn[i][bndFaceI] = Yeff0[i];
                        YeffNei[i][bndFaceI] = Yeff0[i];
                        dnOwn[bndFaceI] = dn0;
                        dnNei[bndFaceI] = dn0;
                    }                    
                }// end if ownCell is ph-0: if(alpha1Own <= ALPHA_2PH_MIN && alpha1Nei >= ALPHA_2PH_MAX)

                VNeiFld[bndFaceI] = V[faceOwn];
                
                faceI++;
                bndFaceI++;
            }// end loop over fcI for patchI: forAll(Js1[0].boundaryField()[patchI], fcI)
        }// end if(pp.coupled())
    }// end loop over patches: forAll(Js1[0].boundaryField(), patchI)

    //----------------------------------------------------------------------------------//
    // Swap the dnNei, YeffNei[i] and VNeiFld lists which contain
    // values of owner cells of patch faces with corresponding lists
    // on neighbouring processor or on coupled cyclic patch.  After
    // swapping, these lists will now contain the values of the
    // neighbour cells of coupled patch faces.  If owner cell is ph-1,
    // these lists will contain ph-0 side values and if owner cell is
    // ph-0, these lists will contain ph-1 side values
    //----------------------------------------------------------------------------------//
    syncTools::swapBoundaryFaceList(mesh, dnNei);
    for(label i=0; i<nSpecies; i++)
    {
        syncTools::swapBoundaryFaceList(mesh, YeffNei[i]);
    }
    syncTools::swapBoundaryFaceList(mesh, VNeiFld);
    //----------------------------------------------------------------------------------//
    //----------------------------------------------------------------------------------//

    if(debug)
    {
        os<< "-------------------------------------------------------------------------" << nl
          << "Boundary coupled faces step 2" << nl
          << "-------------------------------------------------------------------------" << nl
          << endl;
    }

    // start loop over all boundary patches
    forAll(Js1[0].boundaryField(), patchI)
    {
        const polyPatch& pp = patches[patchI];

        if(debug)
        {
            os<< "-------------------------------------------------------------------------" << nl
                << "Patch: " << patchNames[patchI] << nl
                << "-------------------------------------------------------------------------" << nl
                << endl;
        }

        // do only for coupled patches
        if(pp.coupled())
        {
            // faceI is face index in list of all faces
            label faceI = pp.start();
            // bndFaceI is face index in list of all boundary faces
            label bndFaceI = faceI - mesh.nInternalFaces();

            // Reference to neighbour cell fields 
            const scalarField& alpha1NeiFld = alpha1.boundaryField()[patchI].patchNeighbourField();            
            const vectorField& nHatNeiFld = nHat.boundaryField()[patchI].patchNeighbourField();
            const vectorField& C_intfcNeiFld = C_intfc.boundaryField()[patchI].patchNeighbourField();
            const scalarField& A_intfcNeiFld = A_intfc.boundaryField()[patchI].patchNeighbourField();
            const scalarField& rho1NeiFld = rho1.boundaryField()[patchI].patchNeighbourField();            
            const scalarField& rho0NeiFld = rho0.boundaryField()[patchI].patchNeighbourField();
            const fvsPatchVectorField& pSf = Sf.boundaryField()[patchI];
            const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchI];
            const fvsPatchVectorField& pCf = Cf.boundaryField()[patchI];
            List<scalarField> C1NeiFld(nSpecies);
            List<scalarField> C0NeiFld(nSpecies);
            List<scalarField> D1NeiFld(nSpecies);
            List<scalarField> D0NeiFld(nSpecies);
            for(label i=0; i<nSpecies; i++)
            {
                C1NeiFld[i] = C1[i].boundaryField()[patchI].patchNeighbourField();
                C0NeiFld[i] = C0[i].boundaryField()[patchI].patchNeighbourField();
                D1NeiFld[i] = D1[i].boundaryField()[patchI].patchNeighbourField();
                D0NeiFld[i] = D0[i].boundaryField()[patchI].patchNeighbourField();
            }

            // start loop over all faces fcI on patchI (only if coupled patch)
            forAll(Js1[0].boundaryField()[patchI], fcI)
            {
                label faceOwn = own[faceI];
                scalar alpha1Own = alpha1Cells[faceOwn];
                scalar alpha1Nei = alpha1NeiFld[fcI];
                scalar VOwn = V[faceOwn];
                scalar VNei = VNeiFld[bndFaceI];
                scalar rho1Own  = rho1Cells[faceOwn];
                scalar rho0Own  = rho0Cells[faceOwn];
                scalar rho1Nei = rho1NeiFld[fcI];
                scalar rho0Nei = rho0NeiFld[fcI];
                List<scalar> C1Own(nSpecies);
                List<scalar> C0Own(nSpecies);
                List<scalar> C1Nei(nSpecies);
                List<scalar> C0Nei(nSpecies);
                List<scalar> D1Own(nSpecies);
                List<scalar> D0Own(nSpecies);
                List<scalar> D1Nei(nSpecies);
                List<scalar> D0Nei(nSpecies);
                for(label i=0; i<nSpecies; i++)
                {
                    C1Own[i] = C1[i].internalField()[faceOwn];
                    C0Own[i] = C0[i].internalField()[faceOwn];
                    C1Nei[i] = C1NeiFld[i][fcI];
                    C0Nei[i] = C0NeiFld[i][fcI];
                    D1Own[i] = D1[i].internalField()[faceOwn];
                    D0Own[i] = D0[i].internalField()[faceOwn];
                    D1Nei[i] = D1NeiFld[i][fcI];
                    D0Nei[i] = D0NeiFld[i][fcI];
                }
                
                scalar A_intfc_cellI = pMagSf[fcI];                
                vector C_intfc_cellI = pCf[fcI];

                if(debug)
                {
                    os<< "Face: " << faceI << "  patch face index: " << fcI << "  bnd face index: " << bndFaceI << nl
                        << "own: " << faceOwn << "  alpha1Own = " << alpha1Own << "  alpha1Nei = " << alpha1Nei 
                        << endl;
                }

                scalar mS1Tot_own = 0;
                scalar mS1i_cellI = 0;
                scalar mS0Tot_own = 0;
                scalar mS0i_cellI = 0;
                scalar max_mSi = 0;

                // if own cell is ph-1
                if(alpha1Own >= ALPHA_2PH_MAX && alpha1Nei <= ALPHA_2PH_MIN)
                {
                    vector nf = -pSf[fcI]/pMagSf[fcI];

                    scalar dn1 = dnOwn[bndFaceI];
                    scalar dn0 = dnNei[bndFaceI];
                    List<scalar> Yeff1(nSpecies);                    
                    List<scalar> Yeff0(nSpecies);
                    for(label i=0; i<nSpecies; i++)
                    {
                        Yeff1[i] = YeffOwn[i][bndFaceI];
                        Yeff0[i] = YeffNei[i][bndFaceI];
                    }
                    scalar intfcGradi_cellI;

                    if(debug)
                    {
                        os<< "ph-1 cell: " << faceOwn << nl
                            << "nfOwn = " << nHatCells[faceOwn] << "  nfNei = " << nHatNeiFld[fcI] << nl
                            << "C_intfcOwn = " << C_intfcCells[faceOwn] << "  C_intfcNei = " << C_intfcNeiFld[fcI] << nl
                            << "A_intfcOwn = " << A_intfcCells[faceOwn] << "  A_intfcNei = " << A_intfcNeiFld[fcI] << nl
                            << "nf = " << nf << "  C_intfc = " << C_intfc_cellI << "  A_intfc = " << A_intfc_cellI << nl
                            << "C_ph1Own = " << C_ph1_flatFld[faceOwn] << "  C_ph0Own = " << C_ph0_flatFld[faceOwn] << nl
                            << "rho1Own = " << rho1Cells[faceOwn] << "  rho0Own = " << rho0Cells[faceOwn] << nl                            
                            << "rho1Nei = " << rho1NeiFld[fcI] << "  rho0Nei = " << rho0NeiFld[fcI] << endl;
                    }

                    for(label i=0; i<(nSpecies-1); i++)
                    {     
                        Ys0[i].internalField()[faceOwn] = rho0Nei*D0Nei[i]*Yeff0[i]*dn1 + rho1Own*D1Own[i]*Yeff1[i]*dn0;
                        Ys0[i].internalField()[faceOwn] /= (rho0Nei*D0Nei[i]*dn1 + He[i]*rho1Own*D1Own[i]*dn0);                

                        Ys1[i].internalField()[faceOwn] = He[i]*Ys0[i].internalField()[faceOwn];                        

                        intfcGradi_cellI = (Yeff1[i] - Ys1[i].internalField()[faceOwn])/dn1;
                        Js1[i].internalField()[faceOwn] = -A_intfc_cellI*rho1Own*D1Own[i]*intfcGradi_cellI;

                        Js0[i].internalField()[faceOwn] = Js1[i].internalField()[faceOwn];

                        mS1i_cellI = Js1[i].internalField()[faceOwn];
                
                        if(mS1i_cellI > 0)
                        {
                            max_mSi = C0Nei[i]*VNei/dt;

                            mS1[i].internalField()[faceOwn] = min(mS1i_cellI, max_mSi);                    
                        }
                        else if(mS1i_cellI < 0)
                        {
                            max_mSi = C1Own[i]*VOwn/dt;

                            mS1[i].internalField()[faceOwn] = -min(-mS1i_cellI, max_mSi);
                        }
                        else
                        {
                            mS1[i].internalField()[faceOwn] = 0;
                        }
                        
                        mS1[i].internalField()[faceOwn] /= VOwn;                        
                        mS0[i].internalField()[faceOwn] = 0;                        

                        mS1Tot_own += mS1[i].internalField()[faceOwn];
                        mS0Tot_own += mS0[i].internalField()[faceOwn];
                    }            

                    mS1[nSpecies-1].internalField()[faceOwn] = -mS1Tot_own;
                    mS0[nSpecies-1].internalField()[faceOwn] = -mS0Tot_own;	   
                }// end if own cell is ph-1:
                 // if(alpha1Own >= ALPHA_2PH_MAX && alpha1Nei <= ALPHA_2PH_MIN)

                // if own cell is ph-0
                if(alpha1Own <= ALPHA_2PH_MIN && alpha1Nei >= ALPHA_2PH_MAX)
                {
                    vector nf = pSf[fcI]/pMagSf[fcI];
                    
                    scalar dn1 = dnNei[bndFaceI];
                    scalar dn0 = dnOwn[bndFaceI];
                    List<scalar> Yeff1(nSpecies);                    
                    List<scalar> Yeff0(nSpecies);
                    for(label i=0; i<nSpecies; i++)
                    {
                        Yeff1[i] = YeffNei[i][bndFaceI];
                        Yeff0[i] = YeffOwn[i][bndFaceI];
                    }
                    scalar intfcGradi_cellI;

                    if(debug)
                    {
                        os<< "ph-0 cell: " << faceOwn << nl
                            << "nfOwn = " << nHatCells[faceOwn] << "  nfNei = " << nHatNeiFld[fcI] << nl
                            << "C_intfcOwn = " << C_intfcCells[faceOwn] << "  C_intfcNei = " << C_intfcNeiFld[fcI] << nl
                            << "A_intfcOwn = " << A_intfcCells[faceOwn] << "  A_intfcNei = " << A_intfcNeiFld[fcI] << nl
                            << "nf = " << nf << "  C_intfc = " << C_intfc_cellI << "  A_intfc = " << A_intfc_cellI << nl
                            << "C_ph1Own = " << C_ph1_flatFld[faceOwn] << "  C_ph0Own = " << C_ph0_flatFld[faceOwn] << nl
                            << "rho1Own = " << rho1Cells[faceOwn] << "  rho0Own = " << rho0Cells[faceOwn] << nl                            
                            << "rho1Nei = " << rho1NeiFld[fcI] << "  rho0Nei = " << rho0NeiFld[fcI] << endl;
                    }

                    for(label i=0; i<(nSpecies-1); i++)
                    {     
                        Ys0[i].internalField()[faceOwn] = rho0Own*D0Own[i]*Yeff0[i]*dn1 + rho1Nei*D1Nei[i]*Yeff1[i]*dn0;
                        Ys0[i].internalField()[faceOwn] /= (rho0Own*D0Own[i]*dn1 + He[i]*rho1Nei*D1Nei[i]*dn0);                

                        Ys1[i].internalField()[faceOwn] = He[i]*Ys0[i].internalField()[faceOwn];                        

                        intfcGradi_cellI = (Yeff0[i] - Ys0[i].internalField()[faceOwn])/dn1;
                        Js0[i].internalField()[faceOwn] = A_intfc_cellI*rho0Own*D0Own[i]*intfcGradi_cellI;

                        Js1[i].internalField()[faceOwn] = Js0[i].internalField()[faceOwn];

                        mS0i_cellI = -Js0[i].internalField()[faceOwn];
                
                        if(mS0i_cellI > 0)
                        {
                            max_mSi = C1Nei[i]*VNei/dt;

                            mS0[i].internalField()[faceOwn] = min(mS0i_cellI, max_mSi);                    
                        }
                        else if(mS0i_cellI < 0)
                        {
                            max_mSi = C0Own[i]*VOwn/dt;

                            mS0[i].internalField()[faceOwn] = -min(-mS0i_cellI, max_mSi);
                        }
                        else
                        {
                            mS0[i].internalField()[faceOwn] = 0;
                        }
                        
                        mS0[i].internalField()[faceOwn] /= VOwn;                        
                        mS1[i].internalField()[faceOwn] = 0;                        

                        mS1Tot_own += mS1[i].internalField()[faceOwn];
                        mS0Tot_own += mS0[i].internalField()[faceOwn];
                    }            

                    mS1[nSpecies-1].internalField()[faceOwn] = -mS1Tot_own;
                    mS0[nSpecies-1].internalField()[faceOwn] = -mS0Tot_own;	   
                }// end if own cell is ph-0:
                 // if(alpha1Own <= ALPHA_2PH_MIN && alpha1Nei >= ALPHA_2PH_MAX)

                faceI++;
                bndFaceI++;

                if(debug)
                {
                    os<< "-----------------------------------------------------------------------------------" << endl;
                    for(label i=0; i<nSpecies; i++)
                    {
                        os<< "face owner: " << nl
                            << "C1[" << i << "] = " << C1[i].internalField()[faceOwn] << "  C0[" << i << "] = " << C0[i].internalField()[faceOwn] << nl
                            << "Js1[" << i << "] = " << Js1[i].internalField()[faceOwn] << "  Js0[" << i << "] = " << Js0[i].internalField()[faceOwn] << nl 
                            << "mS1[" << i << "] = " << mS1[i].internalField()[faceOwn] << "  mS0[" << i << "] = " << mS0[i].internalField()[faceOwn] 
                            << endl;
                    }                                
                    os<< "-----------------------------------------------------------------------------------" << nl
                        << endl;
                }
            }// end loop over all faces fcI on patchI (only if coupled patch):
             // forAll(Js1[0].boundaryField()[patchI], fcI)
        }// end do only for coupled patches: if(pp.coupled())
    }// end loop over all boundary patches: forAll(Js1[0].boundaryField(), patchI)

    if(debug)
    {
        os<< nl
            << " Done Interfacial Species Transfer Source Terms Calculation (Henry's Law Model)" << nl
            << "---------------------------------------------------------------------------------" << nl
            << endl;
    }    
}


void print_advFluxFld
(
    const fvMesh& mesh,
    const surfaceScalarField& advFlux,
    const volScalarField& alpha,
    const volScalarField& C,
    const volScalarField& c,
    const volScalarField& Y,
    const word& cName,
    const word& YName,
    OFstream& os
)
{
    label faceI, faceOwn, faceNei;
    word ownYName = "own" + YName;
    word neiYName = "nei" + YName;

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();
    const scalarField& YCells = Y.internalField();

    print_line(os, 80);
    print_line(os, 80);
    os<< "Advective fluxes " << cName << endl;
    print_line(os, 80);
    os<< setw(5) << "faceI" << "  " << setw(12) << "advFlux" << "  " << setw(12) << ownYName << "  " << setw(12) << neiYName << endl;
    print_line(os, 80);
    for(faceI=0; faceI<mesh.nInternalFaces(); faceI++)
    {
        faceOwn = own[faceI];
        faceNei = nei[faceI];

        os<< setw(5) << faceI << "  " << setw(12) << advFlux[faceI] << "  " << setw(12) << YCells[faceOwn] << "  " << setw(12) << YCells[faceNei] << endl;        
    }
    print_line(os, 80);
    
    forAll(advFlux.boundaryField(), patchI)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchI];
        const fvsPatchScalarField& padvFlux = advFlux.boundaryField()[patchI];
        const fvPatchScalarField& pY = Y.boundaryField()[patchI];

        faceI = pp.start();

        print_line(os, 80);
        os<< "Patch " << mesh.boundaryMesh().names()[patchI] << endl;
        print_line(os, 80);
        os<< setw(5) << "faceI" << "  " << setw(12) << "advFlux" << "  " << setw(12) << ownYName << "  " << setw(12) << neiYName << endl;
        print_line(os, 80);
        forAll(padvFlux, fcI)
        {
            faceOwn = own[faceI];
            os<< setw(5) << faceI << "  " << setw(12) << padvFlux[fcI] << "  " << setw(12) << YCells[faceOwn] << "  " << setw(12) << pY[fcI] << endl;            
            faceI++;
        }
    }
    print_line(os, 80);
    print_line(os, 80);
}


void print_advFluxFld
(
    const fvMesh& mesh,
    const surfaceScalarField& advFlux,
    const volScalarField& alpha,
    const word& alphaName,
    OFstream& os
)
{
    label faceI, faceOwn, faceNei;
    word ownAlphaName = "own" + alphaName;
    word neiAlphaName = "nei" + alphaName;

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();
    const scalarField& alphaCells = alpha.internalField();

    print_line(os, 80);
    print_line(os, 80);
    os<< "Advective fluxes " << alphaName << endl;
    print_line(os, 80);
    os<< setw(5) << "faceI" << "  " << setw(12) << "advFlux" << "  " << setw(12) << ownAlphaName << "  " << setw(12) << neiAlphaName << endl;
    print_line(os, 80);
    for(faceI=0; faceI<mesh.nInternalFaces(); faceI++)
    {
        faceOwn = own[faceI];
        faceNei = nei[faceI];

        os<< setw(5) << faceI << "  " << setw(12) << advFlux[faceI] << "  " << setw(12) << alphaCells[faceOwn] << "  " << setw(12) << alphaCells[faceNei] << endl;        
    }
    print_line(os, 80);
    
    forAll(advFlux.boundaryField(), patchI)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchI];
        const fvsPatchScalarField& padvFlux = advFlux.boundaryField()[patchI];
        const fvPatchScalarField& palpha = alpha.boundaryField()[patchI];

        faceI = pp.start();

        print_line(os, 80);
        os<< "Patch " << mesh.boundaryMesh().names()[patchI] << endl;
        print_line(os, 80);
        os<< setw(5) << "faceI" << "  " << setw(12) << "advFlux" << "  " << setw(12) << ownAlphaName << "  " << setw(12) << neiAlphaName << endl;
        print_line(os, 80);
        forAll(padvFlux, fcI)
        {
            faceOwn = own[faceI];
            os<< setw(5) << faceI << "  " << setw(12) << padvFlux[fcI] << "  " << setw(12) << alphaCells[faceOwn] << "  " << setw(12) << palpha[fcI] << endl;            
            faceI++;
        }
    }
    print_line(os, 80);
    print_line(os, 80);
}


void print_advFluxIntData
(
    const fvMesh& mesh,
    const surfaceScalarField& advFlux,
    const scalarField& surfInt_advFlux,
    const scalar& dt,
    const volScalarField& alpha,
    const volScalarField& C,
    const volScalarField& c,
    const volScalarField& Y,
    const scalarField& Af_own,
    const scalarField& Af_nei,
    const word& phaseName,
    const word& alphaName,
    const word& CName,
    const word& cName,
    const word& YName,
    OFstream& os
)
{
    label curFaceLbl, curOwn, curNei;    
    word alphaOldName = alphaName + "Old";
    word advFluxName = "advFlux_" + cName;
    word AfOwnName = "Af_own_" + phaseName;
    word AfNeiName = "Af_nei_" + phaseName;

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();
    const cellList& meshCells = mesh.cells();
    const scalarField& meshV = mesh.V();
    const scalarField& alphaCells = alpha.internalField();
    const scalarField& CCells = C.internalField();
    const scalarField& cCells = c.internalField();
    const scalarField& YCells = Y.internalField();

    print_line(os, 80);
    print_line(os, 80);
    os<< "Advective flux integration in cells " << cName << " " << phaseName << endl;
    print_line(os, 80);

    forAll(YCells,cellI)
    {        
        const cell& curCell = meshCells[cellI];
        print_line(os, 80);
        os<< "Cell " << cellI << endl;
        print_line(os, 80);
        os<< setw(10) << CName << "  " << setw(10) << cName << "  " << setw(10) << YName << "  " << setw(10) << alphaName << "  " << setw(10) << alphaOldName << endl;
        os<< setw(10) << CCells[cellI] << "  " << setw(10) << cCells[cellI] << "  " << setw(10) << YCells[cellI] << "  " << setw(10) << alphaCells[cellI] << "  " << setw(10) << alpha.oldTime().internalField()[cellI] << endl;
        print_line(os, 80);

        forAll(curCell, faceI)
        {
            curFaceLbl = curCell[faceI];
            curOwn = own[curFaceLbl];

            os<< "face " << curFaceLbl << endl; 
            print_line(os, 80);
            os<< setw(6) << "Own" << "  " << setw(10) << CName << "  " << setw(10) << cName << "  " << setw(10) << YName << "  " << setw(10) << alphaName << "  " << setw(10) << alphaOldName << endl;
            os<< setw(6) << curOwn << "  " << setw(10) << CCells[curOwn] << "  " << setw(10) << cCells[curOwn] << "  " << setw(10) << YCells[curOwn] << "  " << setw(10) << alphaCells[curOwn] << "  " << setw(10) << alpha.oldTime().internalField()[curOwn] << endl;
              
            if(curFaceLbl < mesh.nInternalFaces())
            {
                curNei = nei[curFaceLbl];
                os<< setw(6) << "Nei" << "  " << setw(10) << CName << "  " << setw(10) << cName << "  " << setw(10) << YName << "  " << setw(10) << alphaName << "  " << setw(10) << alphaOldName << endl;
                os<< setw(6) << curNei << "  " << setw(10) << CCells[curNei] << "  " << setw(10) << cCells[curNei] << "  " << setw(10) << YCells[curNei] << "  " << setw(10) << alphaCells[curNei] << "  " << setw(10) << alpha.oldTime().internalField()[curNei] << endl;
                os<< setw(10) << AfOwnName << "  " << setw(10) << AfNeiName << "  " << advFluxName << endl;
                os<< setw(10) << Af_own[curFaceLbl] << "  " << setw(10) << Af_nei[curFaceLbl] << "  " <<  advFlux[curFaceLbl] << endl;
            }
            else
            {
                os<< setw(10) << AfOwnName << "  " << setw(10) << AfNeiName << endl;
                os<< setw(10) << Af_own[curFaceLbl] << "  " << setw(10) << Af_nei[curFaceLbl] << endl;
            }

            print_line(os, 80);
        }

        os<< setw(16) << "Cell Vol" << "  " << setw(16) << "div(advFlux)" << "  " << setw(10) << "div(advFlux)*dt" << endl;
        os<< setw(16) << meshV[cellI] << "  " << setw(16) << surfInt_advFlux[cellI] << "  " << setw(10) << surfInt_advFlux[cellI]*dt << endl;
        print_line(os, 80);        
    }
    print_line(os, 80);
    print_line(os, 80);
}


void print_advFluxIntData
(
    const fvMesh& mesh,
    const surfaceScalarField& advFlux,
    const scalarField& surfInt_advFlux,
    const scalar& dt,
    const volScalarField& alpha,    
    const scalarField& Af_own,
    const scalarField& Af_nei,
    const word& phaseName,
    const word& alphaName,    
    OFstream& os
)
{
    label curFaceLbl, curOwn, curNei;    
    word alphaOldName = alphaName + "Old";
    word advFluxName = "advFlux_" + alphaName;
    word AfOwnName = "Af_own_" + phaseName;
    word AfNeiName = "Af_nei_" + phaseName;

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();
    const cellList& meshCells = mesh.cells();
    const scalarField& meshV = mesh.V();
    const scalarField& alphaCells = alpha.internalField();    

    print_line(os, 80);
    print_line(os, 80);
    os<< "Advective flux integration in cells " << alphaName << " " << phaseName << endl;
    print_line(os, 80);

    forAll(alphaCells,cellI)
    {        
        const cell& curCell = meshCells[cellI];
        print_line(os, 80);
        os<< "Cell " << cellI << endl;
        print_line(os, 80);
        os<< setw(10) << alphaName << "  " << setw(10) << alphaOldName << endl;
        os<< setw(10) << alphaCells[cellI] << "  " << setw(10) << alpha.oldTime().internalField()[cellI] << endl;
        print_line(os, 80);

        forAll(curCell, faceI)
        {
            curFaceLbl = curCell[faceI];
            curOwn = own[curFaceLbl];

            os<< "face " << curFaceLbl << endl; 
            print_line(os, 80);
            os<< setw(6) << "Own" << "  " << setw(10) << alphaName << "  " << setw(10) << alphaOldName << endl;
            os<< setw(6) << curOwn << "  " << setw(10) << alphaCells[curOwn] << "  " << setw(10) << alpha.oldTime().internalField()[curOwn] << endl;              
            if(curFaceLbl < mesh.nInternalFaces())
            {
                curNei = nei[curFaceLbl];
                os<< setw(6) << "Nei" << "  " << setw(10) << alphaName << "  " << setw(10) << alphaOldName << endl;
                os<< setw(6) << curNei << "  " << setw(10) << alphaCells[curNei] << "  " << setw(10) << alpha.oldTime().internalField()[curNei] << endl;
                os<< setw(10) << AfOwnName << "  " << setw(10) << AfNeiName << "  " << advFluxName << endl;
                os<< setw(10) << Af_own[curFaceLbl] << "  " << setw(10) << Af_nei[curFaceLbl] << "  " <<  advFlux[curFaceLbl] << endl;
            }
            else
            {
                os<< setw(10) << AfOwnName << "  " << setw(10) << AfNeiName << endl;
                os<< setw(10) << Af_own[curFaceLbl] << "  " << setw(10) << Af_nei[curFaceLbl] << endl;
            }

            print_line(os, 80);
        }

        os<< setw(16) << "Cell Vol" << "  " << setw(16) << "div(advFlux)" << "  " << setw(10) << "div(advFlux)*dt" << endl;
        os<< setw(16) << meshV[cellI] << "  " << setw(16) << surfInt_advFlux[cellI] << "  " << setw(10) << surfInt_advFlux[cellI]*dt << endl;
        print_line(os, 80);        
    }
    print_line(os, 80);
    print_line(os, 80);
}


void print_diffFluxFld
(
    const fvMesh& mesh,
    const surfaceScalarField& diffFlux,
    const volScalarField& alpha,
    const volScalarField& C,
    const volScalarField& c,
    const volScalarField& Y,
    const word& cName,
    const word& YName,
    OFstream& os
)
{
    label faceI, faceOwn, faceNei;
    word ownYName = "own" + YName;
    word neiYName = "nei" + YName;

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();
    const scalarField& YCells = Y.internalField();

    print_line(os, 80);
    print_line(os, 80);
    os<< "Diffusive fluxes " << cName << endl;
    print_line(os, 80);
    os<< setw(5) << "faceI" << "  " << setw(12) << "diffFlux" << "  " << setw(12) << ownYName << "  " << setw(12) << neiYName << endl;
    print_line(os, 80);
    for(faceI=0; faceI<mesh.nInternalFaces(); faceI++)
    {
        faceOwn = own[faceI];
        faceNei = nei[faceI];

        os<< setw(5) << faceI << "  " << setw(12) << diffFlux[faceI] << "  " << setw(12) << YCells[faceOwn] << "  " << setw(12) << YCells[faceNei] << endl;        
    }
    print_line(os, 80);
    
    forAll(diffFlux.boundaryField(), patchI)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchI];
        const fvsPatchScalarField& pdiffFlux = diffFlux.boundaryField()[patchI];
        const fvPatchScalarField& pY = Y.boundaryField()[patchI];

        faceI = pp.start();

        print_line(os, 80);
        os<< "Patch " << mesh.boundaryMesh().names()[patchI] << endl;
        print_line(os, 80);
        os<< setw(5) << "faceI" << "  " << setw(12) << "diffFlux" << "  " << setw(12) << ownYName << "  " << setw(12) << neiYName << endl;
        print_line(os, 80);
        forAll(pdiffFlux, fcI)
        {
            faceOwn = own[faceI];
            os<< setw(5) << faceI << "  " << setw(12) << pdiffFlux[fcI] << "  " << setw(12) << YCells[faceOwn] << "  " << setw(12) << pY[fcI] << endl;            
            faceI++;
        }
    }
    print_line(os, 80);
    print_line(os, 80);
}


void print_diffFluxIntData
(
    const fvMesh& mesh,
    const surfaceScalarField& diffFlux,
    const scalarField& surfInt_diffFlux,
    const scalar& dt,
    const volScalarField& alpha,
    const volScalarField& C,
    const volScalarField& c,
    const volScalarField& Y,
    const scalarField& Af_own,
    const scalarField& Af_nei,
    const word& phaseName,
    const word& alphaName,
    const word& CName,
    const word& cName,
    const word& YName,
    OFstream& os
)
{
    label curFaceLbl, curOwn, curNei;    
    word alphaOldName = alphaName + "Old";
    word diffFluxName = "diffFlux_" + cName;
    word AfOwnName = "Af_own_" + phaseName;
    word AfNeiName = "Af_nei_" + phaseName;

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();
    const cellList& meshCells = mesh.cells();
    const scalarField& meshV = mesh.V();
    const scalarField& alphaCells = alpha.internalField();
    const scalarField& CCells = C.internalField();
    const scalarField& cCells = c.internalField();
    const scalarField& YCells = Y.internalField();

    print_line(os, 80);
    print_line(os, 80);
    os<< "Diffusive flux integration in cells " << cName << " " << phaseName << endl;
    print_line(os, 80);

    forAll(YCells,cellI)
    {        
        const cell& curCell = meshCells[cellI];
        print_line(os, 80);
        os<< "Cell " << cellI << endl;
        print_line(os, 80);
        os<< setw(10) << CName << "  " << setw(10) << cName << "  " << setw(10) << YName << "  " << setw(10) << alphaName << "  " << setw(10) << alphaOldName << endl;
        os<< setw(10) << CCells[cellI] << "  " << setw(10) << cCells[cellI] << "  " << setw(10) << YCells[cellI] << "  " << setw(10) << alphaCells[cellI] << "  " << setw(10) << alpha.oldTime().internalField()[cellI] << endl;
        print_line(os, 80);

        forAll(curCell, faceI)
        {
            curFaceLbl = curCell[faceI];
            curOwn = own[curFaceLbl];

            os<< "face " << curFaceLbl << endl; 
            print_line(os, 80);
            os<< setw(6) << "Own" << "  " << setw(10) << CName << "  " << setw(10) << cName << "  " << setw(10) << YName << "  " << setw(10) << alphaName << "  " << setw(10) << alphaOldName << endl;
            os<< setw(6) << curOwn << "  " << setw(10) << CCells[curOwn] << "  " << setw(10) << cCells[curOwn] << "  " << setw(10) << YCells[curOwn] << "  " << setw(10) << alphaCells[curOwn] << "  " << setw(10) << alpha.oldTime().internalField()[curOwn] << endl;
              
            if(curFaceLbl < mesh.nInternalFaces())
            {
                curNei = nei[curFaceLbl];
                os<< setw(6) << "Nei" << "  " << setw(10) << CName << "  " << setw(10) << cName << "  " << setw(10) << YName << "  " << setw(10) << alphaName << "  " << setw(10) << alphaOldName << endl;
                os<< setw(6) << curNei << "  " << setw(10) << CCells[curNei] << "  " << setw(10) << cCells[curNei] << "  " << setw(10) << YCells[curNei] << "  " << setw(10) << alphaCells[curNei] << "  " << setw(10) << alpha.oldTime().internalField()[curNei] << endl;
                os<< setw(10) << AfOwnName << "  " << setw(10) << AfNeiName << "  " << diffFluxName << endl;
                os<< setw(10) << Af_own[curFaceLbl] << "  " << setw(10) << Af_nei[curFaceLbl] << "  " <<  diffFlux[curFaceLbl] << endl;
            }
            else
            {
                os<< setw(10) << AfOwnName << "  " << setw(10) << AfNeiName << endl;
                os<< setw(10) << Af_own[curFaceLbl] << "  " << setw(10) << Af_nei[curFaceLbl] << endl;
            }

            print_line(os, 80);
        }

        os<< setw(16) << "Cell Vol" << "  " << setw(16) << "div(diffFlux)" << "  " << setw(10) << "div(diffFlux)*dt" << endl;
        os<< setw(16) << meshV[cellI] << "  " << setw(16) << surfInt_diffFlux[cellI] << "  " << setw(10) << surfInt_diffFlux[cellI]*dt << endl;
        print_line(os, 80);        
    }
    print_line(os, 80);
    print_line(os, 80);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace plicFuncs

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
