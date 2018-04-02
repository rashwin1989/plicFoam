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
    const bool debug,
    OFstream& os
)
{
    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();    
    const surfaceScalarField& meshMagSf = mesh.magSf();    
    const scalarField& meshV = mesh.V();
        
    scalar curMagSf_ph1;
    scalar curMagSf_ph0;
    label curPhaseState;    
    label faceOwn; 
    label faceNei;    
    scalar diffFlux_limiter;

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

    for(label faceI=0; faceI<mesh.nInternalFaces(); faceI++)
    {
        curMagSf_ph1 = min(magSf_ph1_own[faceI], magSf_ph1_nei[faceI]);
        curMagSf_ph0 = min(magSf_ph0_own[faceI], magSf_ph0_nei[faceI]);
        curPhaseState = face_phaseState_diff[faceI];
        faceOwn = own[faceI];
        faceNei = nei[faceI];                        

        if(curPhaseState == 3)
        {
            diffFlux_Y0i[faceI] = 0;
            diffFlux_Y1i[faceI] = 0;
        }
        else if(curPhaseState == 0)
        {
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
                mesh.time().deltaTValue(),
                diffFlux_Y0i[faceI],
                diffFlux_limiter,
                Y0MIN,
                Y0MAX
            );
            diffFlux_Y0i[faceI] *= min(diffFlux_limiter, 1);
            diffFlux_Y1i[faceI] = 0;
        }
        else if(curPhaseState == 1)
        {
            diffFlux_Y0i[faceI] = 0;
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
                mesh.time().deltaTValue(),
                diffFlux_Y1i[faceI],
                diffFlux_limiter,
                Y1MIN,
                Y1MAX
            );
        }
        else if(curPhaseState == 2)
        {
            diffFlux_Y0i[faceI] = -rho0f[faceI]*D0fi[faceI]*curMagSf_ph0*gradf_Y0i[faceI];
            diffFlux_Y1i[faceI] = -rho1f[faceI]*D1fi[faceI]*curMagSf_ph1*gradf_Y1i[faceI];
        }
        else
        {
            diffFlux_Y0i[faceI] = 0;
            diffFlux_Y1i[faceI] = 0;
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

        //----------------------------------------------------------//
        // check diffFlux > maxDiffFlux causing unboundedness and
        // limit diffFlux if that is the case 
        diffFlux_limiter = 1;
        //ph1        
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
            mesh.time().deltaTValue(),
            diffFlux_Y1i[faceI],
            diffFlux_limiter,
            Y1MIN,
            Y1MAX
        );

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
                    << "rho1Own: " << rho1Cells[faceOwn] << "  alpha1Own: " << alpha1Cells[faceOwn] << "  Y1iOwn: " <<  Y1iCells[faceOwn] << "  dt: " << mesh.time().deltaTValue() << "  VOwn: " <<  meshV[faceOwn] 
                    << nl 
                    << "rho1Nei: " << rho1Cells[faceNei] << "  alpha1Nei: " << alpha1Cells[faceNei] << "  Y1iNei: " <<  Y1iCells[faceNei] << "  dt: " << mesh.time().deltaTValue() << "  VNei: " <<  meshV[faceNei]
                    << nl
                    << endl;
            }            
        }

        diffFlux_Y1i[faceI] *= min(diffFlux_limiter, 1);

        //ph0
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
            mesh.time().deltaTValue(),
            diffFlux_Y0i[faceI],
            diffFlux_limiter,
            Y0MIN,
            Y0MAX
        );        

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
                    << "rho0Own: " << rho0Cells[faceOwn] << "  alpha0Own: " << alpha0Cells[faceOwn] << "  Y0iOwn: " <<  Y0iCells[faceOwn] << "  dt: " << mesh.time().deltaTValue() << "  VOwn: " <<  meshV[faceOwn] 
                    << nl 
                    << "rho0Nei: " << rho0Cells[faceNei] << "  alpha0Nei: " << alpha0Cells[faceNei] << "  Y0iNei: " <<  Y0iCells[faceNei] << "  dt: " << mesh.time().deltaTValue() << "  VNei: " <<  meshV[faceNei]
                    << nl
                    << endl;
            }
        }

        diffFlux_Y0i[faceI] *= min(diffFlux_limiter, 1);

        //end check diffFlux > maxDiffFlux
        //----------------------------------------------------------//

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
    List<scalar> VNei(nBnd);

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        if(pp.coupled())
        {
            label faceI = pp.start();
            label bndFaceI = pp.start() - mesh.nInternalFaces();
            forAll(pp, fcI)
            {
                const label& faceOwn = own[faceI];
                VNei[bndFaceI] = meshV[faceOwn];
                faceI++;
                bndFaceI++;
            }
        }
    }

    syncTools::swapBoundaryFaceList(mesh, VNei);
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
        
        label faceI = pp.start();

        if(pp.coupled())
        {
            label bndFaceI = pp.start() - mesh.nInternalFaces();

            const scalarField& rho1Own = rho1.boundaryField()[patchI].patchInternalField();
            const scalarField& alpha1Own = alpha1.boundaryField()[patchI].patchInternalField();
            const scalarField& Y1iOwn = Y1i.boundaryField()[patchI].patchInternalField();
            const scalarField& rho1Nei = rho1.boundaryField()[patchI].patchNeighbourField();
            const scalarField& alpha1Nei = alpha1.boundaryField()[patchI].patchNeighbourField();
            const scalarField& Y1iNei = Y1i.boundaryField()[patchI].patchNeighbourField();

            const scalarField& rho0Own = rho0.boundaryField()[patchI].patchInternalField();
            const scalarField& alpha0Own = alpha0.boundaryField()[patchI].patchInternalField();
            const scalarField& Y0iOwn = Y0i.boundaryField()[patchI].patchInternalField();
            const scalarField& rho0Nei = rho0.boundaryField()[patchI].patchNeighbourField();
            const scalarField& alpha0Nei = alpha0.boundaryField()[patchI].patchNeighbourField();
            const scalarField& Y0iNei = Y0i.boundaryField()[patchI].patchNeighbourField();
            
            forAll(pY1i, fcI)
            {                
                curMagSf_ph1 = min(magSf_ph1_own[faceI], magSf_ph1_nei[faceI]);
                curMagSf_ph0 = min(magSf_ph0_own[faceI], magSf_ph0_nei[faceI]);
                curPhaseState = face_phaseState_diff[faceI];
                faceOwn = own[faceI];                

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
                        rho0Own[fcI],
                        alpha0Own[fcI],
                        Y0iOwn[fcI],
                        meshV[faceOwn],
                        rho0Nei[fcI],
                        alpha0Nei[fcI],
                        Y0iNei[fcI],
                        VNei[bndFaceI],
                        mesh.time().deltaTValue(),
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
                    calc_diffFlux_limiter2
                    (
                        rho1Own[fcI],
                        alpha1Own[fcI],
                        Y1iOwn[fcI],
                        meshV[faceOwn],
                        rho1Nei[fcI],
                        alpha1Nei[fcI],
                        Y1iNei[fcI],
                        VNei[bndFaceI],
                        mesh.time().deltaTValue(),
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
                    os<< "Face: " << faceI << "  mag(Sf) = " << pMagSf[fcI] << nl
                        << "magSf_ph1_own = " << magSf_ph1_own[faceI] << "  alpha1f_own = " << magSf_ph1_own[faceI]/pMagSf[fcI] 
                        << "  magSf_ph1_nei = " << magSf_ph1_nei[faceI] << "  alpha1f_nei = " << magSf_ph1_nei[faceI]/pMagSf[fcI] << nl
                        << "magSf_ph0_own = " << magSf_ph0_own[faceI] << "  alpha0f_own = " << magSf_ph0_own[faceI]/pMagSf[fcI] 
                        << "  magSf_ph0_nei = " << magSf_ph0_nei[faceI] << "  alpha0f_nei = " << magSf_ph0_nei[faceI]/pMagSf[fcI] << nl
                        << "face phase state for diffusion flux calculation: " << curPhaseState << nl
                        << "Own: " << faceOwn << nl
                        << "Own Y1" << i << ": " << Y1iOwn[fcI] << "  Nei Y1" << i << ": " << Y1iNei[fcI] << nl
                        << "Own Y0" << i << ": " << Y0iOwn[fcI] << "  Nei Y0" << i << ": " << Y0iNei[fcI] << nl               
                        << "diffusion flux ph1 for Y" << i << " = " << pdiffFlux_Y1i[fcI] << nl
                        << "diffusion flux ph0 for Y" << i << " = " << pdiffFlux_Y0i[fcI] << nl
                        << endl;
                }

                // check diffFlux > maxDiffFlux causing unboundedness and
                // limit diffFlux if that is the case 
                
                //ph1        
                diffFlux_limiter = 1;

                calc_diffFlux_limiter2
                (
                    rho1Own[fcI],
                    alpha1Own[fcI],
                    Y1iOwn[fcI],
                    meshV[faceOwn],
                    rho1Nei[fcI],
                    alpha1Nei[fcI],
                    Y1iNei[fcI],
                    VNei[bndFaceI],
                    mesh.time().deltaTValue(),
                    pdiffFlux_Y1i[fcI],
                    diffFlux_limiter,
                    Y1MIN,
                    Y1MAX
                );

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
                            << "rho1Own: " << rho1Own[fcI] << "  alpha1Own: " << alpha1Own[fcI] << "  Y1iOwn: " <<  Y1iOwn[fcI] << "  dt: " << mesh.time().deltaTValue() << "  VOwn: " <<  meshV[faceOwn] 
                            << nl 
                            << "rho1Nei: " << rho1Nei[fcI] << "  alpha1Nei: " << alpha1Nei[fcI] << "  Y1iNei: " <<  Y1iNei[fcI] << "  dt: " << mesh.time().deltaTValue() << "  VNei: " <<  VNei[fcI] 
                            << nl
                            << endl;
                    }            
                }

                pdiffFlux_Y1i[fcI] *= min(diffFlux_limiter, 1);
                
                //ph0        
                diffFlux_limiter = 1;

                calc_diffFlux_limiter2
                (
                    rho0Own[fcI],
                    alpha0Own[fcI],
                    Y0iOwn[fcI],
                    meshV[faceOwn],
                    rho0Nei[fcI],
                    alpha0Nei[fcI],
                    Y0iNei[fcI],
                    VNei[bndFaceI],
                    mesh.time().deltaTValue(),
                    pdiffFlux_Y0i[fcI],
                    diffFlux_limiter,
                    Y0MIN,
                    Y0MAX
                );

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
                            << "rho0Own: " << rho0Own[fcI] << "  alpha0Own: " << alpha0Own[fcI] << "  Y0iOwn: " <<  Y0iOwn[fcI] << "  dt: " << mesh.time().deltaTValue() << "  VOwn: " <<  meshV[faceOwn] 
                            << nl 
                            << "rho0Nei: " << rho0Nei[fcI] << "  alpha0Nei: " << alpha0Nei[fcI] << "  Y0iNei: " <<  Y0iNei[fcI] << "  dt: " << mesh.time().deltaTValue() << "  VNei: " <<  VNei[fcI] 
                            << nl
                            << endl;
                    }            
                }

                pdiffFlux_Y0i[fcI] *= min(diffFlux_limiter, 1);
                
                faceI++;
                bndFaceI++;
            }//end forAll(pY1i, fcI)
        }//end if(pp.coupled())        
        else if(isA<fixedValueFvPatchScalarField>(pY1i))
        {
            const scalarField& rho1Own = rho1.boundaryField()[patchI].patchInternalField();
            const scalarField& alpha1Own = alpha1.boundaryField()[patchI].patchInternalField();
            const scalarField& Y1iOwn = Y1i.boundaryField()[patchI].patchInternalField();            

            const scalarField& rho0Own = rho0.boundaryField()[patchI].patchInternalField();
            const scalarField& alpha0Own = alpha0.boundaryField()[patchI].patchInternalField();
            const scalarField& Y0iOwn = Y0i.boundaryField()[patchI].patchInternalField();

            forAll(pY1i, fcI)
            {                
                curMagSf_ph1 = magSf_ph1_own[faceI];
                curMagSf_ph0 = magSf_ph0_own[faceI];
                curPhaseState = face_phaseState_diff[faceI];
                faceOwn = own[faceI];                

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
                        rho0Own[fcI],
                        alpha0Own[fcI],
                        Y0iOwn[fcI],
                        meshV[faceOwn],
                        rho0Own[fcI],
                        alpha0Own[fcI],
                        Y0iOwn[fcI],
                        meshV[faceOwn],
                        mesh.time().deltaTValue(),
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
                        rho1Own[fcI],
                        alpha1Own[fcI],
                        Y1iOwn[fcI],
                        meshV[faceOwn],
                        rho1Own[fcI],
                        alpha1Own[fcI],
                        Y1iOwn[fcI],
                        meshV[faceOwn],
                        mesh.time().deltaTValue(),
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

                // check diffFlux > maxDiffFlux causing unboundedness and
                // limit diffFlux if that is the case 
                
                //ph1        
                diffFlux_limiter = 1;
                
                calc_diffFlux_limiter
                (
                    rho1Own[fcI],
                    alpha1Own[fcI],
                    Y1iOwn[fcI],
                    meshV[faceOwn],
                    rho1Own[fcI],
                    alpha1Own[fcI],
                    Y1iOwn[fcI],
                    meshV[faceOwn],
                    mesh.time().deltaTValue(),
                    pdiffFlux_Y1i[fcI],
                    diffFlux_limiter,
                    Y1MIN,
                    Y1MAX
                );

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
                            << "rho1Own: " << rho1Own[fcI] << "  alpha1Own: " << alpha1Own[fcI] << "  Y1iOwn: " <<  Y1iOwn[fcI] << "  dt: " << mesh.time().deltaTValue() << "  VOwn: " <<  meshV[faceOwn] 
                            << nl                             
                            << endl;
                    }            
                }

                pdiffFlux_Y1i[fcI] *= min(diffFlux_limiter, 1);
                
                //ph0        
                diffFlux_limiter = 1;

                calc_diffFlux_limiter
                (
                    rho0Own[fcI],
                    alpha0Own[fcI],
                    Y0iOwn[fcI],
                    meshV[faceOwn],
                    rho0Own[fcI],
                    alpha0Own[fcI],
                    Y0iOwn[fcI],
                    meshV[faceOwn],
                    mesh.time().deltaTValue(),
                    pdiffFlux_Y0i[fcI],
                    diffFlux_limiter,
                    Y0MIN,
                    Y0MAX
                );

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
                            << "rho0Own: " << rho0Own[fcI] << "  alpha0Own: " << alpha0Own[fcI] << "  Y0iOwn: " <<  Y0iOwn[fcI] << "  dt: " << mesh.time().deltaTValue() << "  VOwn: " <<  meshV[faceOwn] 
                            << nl                            
                            << endl;
                    }            
                }

                pdiffFlux_Y0i[fcI] *= min(diffFlux_limiter, 1);

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
        
        for(label faceI=0; faceI<mesh.nInternalFaces(); faceI++)
        {
            os<< "   " << faceI << "          " << face_phaseState_diff[faceI] << "                " << diffFlux_Y1i[faceI] << "                " << diffFlux_Y0i[faceI] << endl;
        }

        os<< "------------------------------------------------------------------------------------------------------" << nl
            << endl;
    }
}


void init_MS_flux
(
    const label& faceI,
    PtrList<surfaceScalarField>& xf,
    PtrList<surfaceScalarField>& gradf_x,    
    PtrList<surfaceScalarField>& Dijf,
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
    const label& patchI,
    const label& faceI,
    PtrList<surfaceScalarField>& xf,
    PtrList<surfaceScalarField>& gradf_x,    
    PtrList<surfaceScalarField>& Dijf,
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
    int *tk,
    double *coef_ab,
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
    int idx;
    double Z; double V_tmp;

    Z = P*V/(R_gas*T);

    fugacities_n_its_derivatives_(&P,&T,&n,Pc,Tc,w,x,tk,coef_ab,lnphi,dlnphi_dxj,&V_tmp);

    // rhs of Maxwell-Stefan equation
    for(int i=0; i<n; i++)
    {
        rhs_flux[i] = -grad_x[i];

        if(n_flux_type > 0)
        {
            for(int j=0; j<n; j++)
            {
                if(x[j]<1e-9) continue;

                idx = i + j*n;
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
    int *tk,
    double *coef_ab,    
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

    _NEW_(x_tmp, double, n);
    _NEW_(grad_x_tmp, double, n);
    _NEW_(Dij_tmp, double, n*n);
    _NEW_(rhs_flux_tmp, double, n);
    _NEW_(flux_m_ph1, double, n);
    _NEW_(flux_m_ph0, double, n);
    _NEW_(lnphi_tmp, double, n);
    _NEW_(dlnphi_dxj_tmp, double, n*n);

    scalar curMagSf_ph1;
    scalar curMagSf_ph0;
    label curPhaseState;    
    label faceOwn; 
    label faceNei;
    scalar diffFlux_limiter_ph1; scalar diffFlux_limiter_ph1i; 
    scalar diffFlux_limiter_ph0; scalar diffFlux_limiter_ph0i;
    List<scalar> curDiffFlux_Y1(n);
    List<scalar> curDiffFlux_Y0(n);
    scalar curMagSf;

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

    //start loop over internal faces
    for(label faceI=0; faceI<mesh.nInternalFaces(); faceI++)
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
            calc_MS_flux(P_tmp,T0_tmp,V0_tmp,x_tmp,R_gas,n,Pc,Tc,w,MW,tk,coef_ab,lnphi_tmp,dlnphi_dxj_tmp,grad_x_tmp,Dij_tmp,n_flux_type,flux_umf,rhs_flux_tmp,flux_m_ph0);

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
            calc_MS_flux(P_tmp,T1_tmp,V1_tmp,x_tmp,R_gas,n,Pc,Tc,w,MW,tk,coef_ab,lnphi_tmp,dlnphi_dxj_tmp,grad_x_tmp,Dij_tmp,n_flux_type,flux_umf,rhs_flux_tmp,flux_m_ph1);

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
            calc_MS_flux(P_tmp,T1_tmp,V1_tmp,x_tmp,R_gas,n,Pc,Tc,w,MW,tk,coef_ab,lnphi_tmp,dlnphi_dxj_tmp,grad_x_tmp,Dij_tmp,n_flux_type,flux_umf,rhs_flux_tmp,flux_m_ph1);

            //ph-0
            init_MS_flux(faceI,x0f,gradf_x0,Dij0f,n,x_tmp,grad_x_tmp,Dij_tmp);
            calc_MS_flux(P_tmp,T0_tmp,V0_tmp,x_tmp,R_gas,n,Pc,Tc,w,MW,tk,coef_ab,lnphi_tmp,dlnphi_dxj_tmp,grad_x_tmp,Dij_tmp,n_flux_type,flux_umf,rhs_flux_tmp,flux_m_ph0);

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
    const label nBnd = mesh.nFaces() - mesh.nInternalFaces();
    List<scalar> VNei(nBnd);

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        if(pp.coupled())
        {
            label faceI = pp.start();
            label bndFaceI = pp.start() - mesh.nInternalFaces();
            forAll(pp, fcI)
            {
                const label& faceOwn = own[faceI];
                VNei[bndFaceI] = meshV[faceOwn];
                faceI++;
                bndFaceI++;
            }
        }
    }

    syncTools::swapBoundaryFaceList(mesh, VNei);
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
        
        label faceI = pp.start();

        if(pp.coupled())
        {
            label bndFaceI = pp.start() - mesh.nInternalFaces();

            const scalarField& rho1Own = rho1.boundaryField()[patchI].patchInternalField();
            const scalarField& alpha1Own = alpha1.boundaryField()[patchI].patchInternalField();            
            const scalarField& rho1Nei = rho1.boundaryField()[patchI].patchNeighbourField();
            const scalarField& alpha1Nei = alpha1.boundaryField()[patchI].patchNeighbourField();

            const scalarField& rho0Own = rho0.boundaryField()[patchI].patchInternalField();
            const scalarField& alpha0Own = alpha0.boundaryField()[patchI].patchInternalField();
            const scalarField& rho0Nei = rho0.boundaryField()[patchI].patchNeighbourField();
            const scalarField& alpha0Nei = alpha0.boundaryField()[patchI].patchNeighbourField();

            List<scalarField> Y1Own(n); List<scalarField> Y1Nei(n);
            List<scalarField> Y0Own(n); List<scalarField> Y0Nei(n);
            for(i=0; i<n; i++)
            {
                Y1Own[i] = Y1[i].boundaryField()[patchI].patchInternalField();
                Y1Nei[i] = Y1[i].boundaryField()[patchI].patchNeighbourField();
                Y0Own[i] = Y0[i].boundaryField()[patchI].patchInternalField();
                Y0Nei[i] = Y0[i].boundaryField()[patchI].patchNeighbourField();
            }
            
            //start loop over all faces of coupled patch
            forAll(alpha1Own, fcI)
            {                
                curMagSf_ph1 = min(magSf_ph1_own[faceI], magSf_ph1_nei[faceI]);
                curMagSf_ph0 = min(magSf_ph0_own[faceI], magSf_ph0_nei[faceI]);
                curPhaseState = face_phaseState_diff[faceI];
                faceOwn = own[faceI];

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
                    calc_MS_flux(P_tmp,T0_tmp,V0_tmp,x_tmp,R_gas,n,Pc,Tc,w,MW,tk,coef_ab,lnphi_tmp,dlnphi_dxj_tmp,grad_x_tmp,Dij_tmp,n_flux_type,flux_umf,rhs_flux_tmp,flux_m_ph0);

                    for(i=0; i<n; i++)
                    {
                        curDiffFlux_Y1[i] = 0;
                        curDiffFlux_Y0[i] = curMagSf_ph0*flux_m_ph0[i];

                        calc_diffFlux_limiter2
                        (
                            rho0Own[fcI],
                            alpha0Own[fcI],
                            Y0Own[i][fcI],
                            meshV[faceOwn],
                            rho0Nei[fcI],
                            alpha0Nei[fcI],
                            Y0Nei[i][fcI],
                            VNei[bndFaceI],
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
                    calc_MS_flux(P_tmp,T1_tmp,V1_tmp,x_tmp,R_gas,n,Pc,Tc,w,MW,tk,coef_ab,lnphi_tmp,dlnphi_dxj_tmp,grad_x_tmp,Dij_tmp,n_flux_type,flux_umf,rhs_flux_tmp,flux_m_ph1);

                    for(i=0; i<n; i++)
                    {
                        curDiffFlux_Y1[i] = curMagSf_ph1*flux_m_ph1[i];
                        curDiffFlux_Y0[i] = 0;                        

                        calc_diffFlux_limiter2
                        (
                            rho1Own[fcI],
                            alpha1Own[fcI],
                            Y1Own[i][fcI],
                            meshV[faceOwn],
                            rho1Nei[fcI],
                            alpha1Nei[fcI],
                            Y1Nei[i][fcI],
                            VNei[bndFaceI],
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
                    calc_MS_flux(P_tmp,T1_tmp,V1_tmp,x_tmp,R_gas,n,Pc,Tc,w,MW,tk,coef_ab,lnphi_tmp,dlnphi_dxj_tmp,grad_x_tmp,Dij_tmp,n_flux_type,flux_umf,rhs_flux_tmp,flux_m_ph1);

                    //ph-0
                    init_MS_flux(patchI,fcI,x0f,gradf_x0,Dij0f,n,x_tmp,grad_x_tmp,Dij_tmp);
                    calc_MS_flux(P_tmp,T0_tmp,V0_tmp,x_tmp,R_gas,n,Pc,Tc,w,MW,tk,coef_ab,lnphi_tmp,dlnphi_dxj_tmp,grad_x_tmp,Dij_tmp,n_flux_type,flux_umf,rhs_flux_tmp,flux_m_ph0);

                    for(i=0; i<n; i++)
                    {
                        curDiffFlux_Y1[i] = curMagSf_ph1*flux_m_ph1[i];
                        curDiffFlux_Y0[i] = curMagSf_ph0*flux_m_ph0[i];

                        calc_diffFlux_limiter2
                        (
                            rho1Own[fcI],
                            alpha1Own[fcI],
                            Y1Own[i][fcI],
                            meshV[faceOwn],
                            rho1Nei[fcI],
                            alpha1Nei[fcI],
                            Y1Nei[i][fcI],
                            VNei[bndFaceI],
                            deltaT,
                            curDiffFlux_Y1[i],
                            diffFlux_limiter_ph1i,
                            Y1MIN[i],
                            Y1MAX[i]
                        );
                        diffFlux_limiter_ph1 = min(diffFlux_limiter_ph1, diffFlux_limiter_ph1i);

                        calc_diffFlux_limiter2
                        (
                            rho0Own[fcI],
                            alpha0Own[fcI],
                            Y0Own[i][fcI],
                            meshV[faceOwn],
                            rho0Nei[fcI],
                            alpha0Nei[fcI],
                            Y0Nei[i][fcI],
                            VNei[bndFaceI],
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
                        << "Vown = " << meshV[faceOwn] << "  VNei = " << VNei[bndFaceI] << "  dt = " << deltaT << nl
                        << "--------------------------------------------------------------------------------------" << nl
                        << "phase-1" << nl
                        << "--------------------------------------------------------------------------------------" << nl
                        << "Own: " << faceOwn << "  alpha1Own = " << alpha1Own[fcI] << "  rho1Own = " << rho1Own[fcI] 
                        << endl;
                    for(i=0; i<n; i++)
                    {
                        os<< "Y1[" << i << "] = " << Y1Own[i][fcI] << "  ";
                    }
                    os<< nl
                        << "--------------------------------------------------------------------------------------" << nl
                        << "Nei: " << "NA" << "  alpha1Nei = " << alpha1Nei[fcI] << "  rho1Nei = " << rho1Nei[fcI]                
                        << endl;
                    for(i=0; i<n; i++)
                    {
                        os<< "Y1[" << i << "] = " << Y1Nei[i][fcI] << "  ";
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
                        << "Own: " << faceOwn << "  alpha0Own = " << alpha0Own[fcI] << "  rho0Own = " << rho0Own[fcI] 
                        << endl;
                    for(i=0; i<n; i++)
                    {
                        os<< "Y0[" << i << "] = " << Y0Own[i][fcI] << "  ";
                    }
                    os<< nl
                        << "--------------------------------------------------------------------------------------" << nl
                        << "Nei: " << "NA" << "  alpha0Nei = " << alpha0Nei[fcI] << "  rho0Nei = " << rho0Nei[fcI]                
                        << endl;
                    for(i=0; i<n; i++)
                    {
                        os<< "Y0[" << i << "] = " << Y0Nei[i][fcI] << "  ";
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
            const scalarField& rho1Own = rho1.boundaryField()[patchI].patchInternalField();
            const scalarField& alpha1Own = alpha1.boundaryField()[patchI].patchInternalField();
            const fvPatchScalarField& rho1Nei = rho1.boundaryField()[patchI];
            const fvPatchScalarField& alpha1Nei = alpha1.boundaryField()[patchI];

            const scalarField& rho0Own = rho0.boundaryField()[patchI].patchInternalField();
            const scalarField& alpha0Own = alpha0.boundaryField()[patchI].patchInternalField();
            const fvPatchScalarField& rho0Nei = rho0.boundaryField()[patchI];
            const fvPatchScalarField& alpha0Nei = alpha0.boundaryField()[patchI];

            List<scalarField> Y1Own(n); List<scalarField> Y1Nei(n);
            List<scalarField> Y0Own(n); List<scalarField> Y0Nei(n);
            for(i=0; i<n; i++)
            {
                Y1Own[i] = Y1[i].boundaryField()[patchI].patchInternalField();
                const scalarField& pY1i = Y1[i].boundaryField()[patchI];
                Y1Nei[i] = pY1i;
                Y0Own[i] = Y0[i].boundaryField()[patchI].patchInternalField();
                const scalarField& pY0i = Y0[i].boundaryField()[patchI];
                Y0Nei[i] = pY0i;
            }

            forAll(alpha1Own, fcI)
            {                
                curMagSf_ph1 = magSf_ph1_own[faceI];
                curMagSf_ph0 = magSf_ph0_own[faceI];
                curPhaseState = face_phaseState_diff[faceI];
                faceOwn = own[faceI];

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
                    calc_MS_flux(P_tmp,T0_tmp,V0_tmp,x_tmp,R_gas,n,Pc,Tc,w,MW,tk,coef_ab,lnphi_tmp,dlnphi_dxj_tmp,grad_x_tmp,Dij_tmp,n_flux_type,flux_umf,rhs_flux_tmp,flux_m_ph0);

                    for(i=0; i<n; i++)
                    {
                        curDiffFlux_Y1[i] = 0;
                        curDiffFlux_Y0[i] = curMagSf_ph0*flux_m_ph0[i];

                        calc_diffFlux_limiter2
                        (
                            rho0Own[fcI],
                            alpha0Own[fcI],
                            Y0Own[i][fcI],
                            meshV[faceOwn],
                            rho0Nei[fcI],
                            alpha0Nei[fcI],
                            Y0Nei[i][fcI],
                            meshV[faceOwn],
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
                    calc_MS_flux(P_tmp,T1_tmp,V1_tmp,x_tmp,R_gas,n,Pc,Tc,w,MW,tk,coef_ab,lnphi_tmp,dlnphi_dxj_tmp,grad_x_tmp,Dij_tmp,n_flux_type,flux_umf,rhs_flux_tmp,flux_m_ph1);

                    for(i=0; i<n; i++)
                    {
                        curDiffFlux_Y1[i] = curMagSf_ph1*flux_m_ph1[i];
                        curDiffFlux_Y0[i] = 0;                        

                        calc_diffFlux_limiter2
                        (
                            rho1Own[fcI],
                            alpha1Own[fcI],
                            Y1Own[i][fcI],
                            meshV[faceOwn],
                            rho1Nei[fcI],
                            alpha1Nei[fcI],
                            Y1Nei[i][fcI],
                            meshV[faceOwn],
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
                    calc_MS_flux(P_tmp,T1_tmp,V1_tmp,x_tmp,R_gas,n,Pc,Tc,w,MW,tk,coef_ab,lnphi_tmp,dlnphi_dxj_tmp,grad_x_tmp,Dij_tmp,n_flux_type,flux_umf,rhs_flux_tmp,flux_m_ph1);

                    //ph-0
                    init_MS_flux(patchI,fcI,x0f,gradf_x0,Dij0f,n,x_tmp,grad_x_tmp,Dij_tmp);
                    calc_MS_flux(P_tmp,T0_tmp,V0_tmp,x_tmp,R_gas,n,Pc,Tc,w,MW,tk,coef_ab,lnphi_tmp,dlnphi_dxj_tmp,grad_x_tmp,Dij_tmp,n_flux_type,flux_umf,rhs_flux_tmp,flux_m_ph0);

                    for(i=0; i<n; i++)
                    {
                        curDiffFlux_Y1[i] = curMagSf_ph1*flux_m_ph1[i];
                        curDiffFlux_Y0[i] = curMagSf_ph0*flux_m_ph0[i];

                        calc_diffFlux_limiter2
                        (
                            rho1Own[fcI],
                            alpha1Own[fcI],
                            Y1Own[i][fcI],
                            meshV[faceOwn],
                            rho1Nei[fcI],
                            alpha1Nei[fcI],
                            Y1Nei[i][fcI],
                            meshV[faceOwn],
                            deltaT,
                            curDiffFlux_Y1[i],
                            diffFlux_limiter_ph1i,
                            Y1MIN[i],
                            Y1MAX[i]
                        );
                        diffFlux_limiter_ph1 = min(diffFlux_limiter_ph1, diffFlux_limiter_ph1i);

                        calc_diffFlux_limiter2
                        (
                            rho0Own[fcI],
                            alpha0Own[fcI],
                            Y0Own[i][fcI],
                            meshV[faceOwn],
                            rho0Nei[fcI],
                            alpha0Nei[fcI],
                            Y0Nei[i][fcI],
                            meshV[faceOwn],
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
                        << "Own: " << faceOwn << "  alpha1Own = " << alpha1Own[fcI] << "  rho1Own = " << rho1Own[fcI] 
                        << endl;
                    for(i=0; i<n; i++)
                    {
                        os<< "Y1[" << i << "] = " << Y1Own[i][fcI] << "  ";
                    }
                    os<< nl
                        << "--------------------------------------------------------------------------------------" << nl
                        << "Nei: " << "NA" << "  alpha1Nei = " << alpha1Nei[fcI] << "  rho1Nei = " << rho1Nei[fcI]                
                        << endl;
                    for(i=0; i<n; i++)
                    {
                        os<< "Y1[" << i << "] = " << Y1Nei[i][fcI] << "  ";
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
                        << "Own: " << faceOwn << "  alpha0Own = " << alpha0Own[fcI] << "  rho0Own = " << rho0Own[fcI] 
                        << endl;
                    for(i=0; i<n; i++)
                    {
                        os<< "Y0[" << i << "] = " << Y0Own[i][fcI] << "  ";
                    }
                    os<< nl
                        << "--------------------------------------------------------------------------------------" << nl
                        << "Nei: " << "NA" << "  alpha0Nei = " << alpha0Nei[fcI] << "  rho0Nei = " << rho0Nei[fcI]                
                        << endl;
                    for(i=0; i<n; i++)
                    {
                        os<< "Y0[" << i << "] = " << Y0Nei[i][fcI] << "  ";
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
        
        for(label faceI=0; faceI<mesh.nInternalFaces(); faceI++)
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

    _DELETE_(x_tmp);
    _DELETE_(grad_x_tmp);
    _DELETE_(Dij_tmp);
    _DELETE_(rhs_flux_tmp);
    _DELETE_(flux_m_ph1);
    _DELETE_(flux_m_ph0);
    _DELETE_(lnphi_tmp);
    _DELETE_(dlnphi_dxj_tmp);
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
    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();

    const Field<Type>& YCells = Y.internalField();

    for(label faceI=0; faceI<mesh.nInternalFaces(); faceI++)
    {
        label faceOwn = own[faceI];
        label faceNei = nei[faceI];
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
        if(pp.coupled())
        {
            const Field<Type>& pYOwn = pY.patchInternalField();
            const Field<Type>& pYNei = pY.patchNeighbourField();
            forAll(pYf, fcI)
            {
                scalar wf = pw[fcI];
                pYf[fcI] = wf*pYOwn[fcI] + (1.0 - wf)*pYNei[fcI];
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
        scalar min_dn = mag(mesh.C()[curCell_lbl] - mesh.C()[curCellsAll[1]]);
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
        if(debug)
        {
            os<< "Cell: " << cellI << "  alpha1 = " << alpha1_cellI << endl;
        }

        if(alpha1_cellI > ALPHA_2PH_MIN && alpha1_cellI < ALPHA_2PH_MAX)
        {
            vector nf = nHatCells[cellI];
            vector C_intfc_cellI = C_intfcCells[cellI];
            scalar A_intfc_cellI = A_intfcCells[cellI];
            labelList curCellsAll = cellStencil[cellI];

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
        }
        else
        {
            for(label i=0; i<nSpecies; i++)
            { 
                Js0[i].internalField()[cellI] = 0;
                Js1[i].internalField()[cellI] = 0;
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

        if(debug)
        {
            os<< "Face: " << faceI << nl
                << "own: " << faceOwn << "  alpha1Own = " << alpha1Own << "  nei: " << faceNei << "  alpha1Nei = " << alpha1Nei << endl;
        }

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
    int *tk,
    double *coef_ab,
    double *Tb,
    double *SG,
    double *H8,
    double *k,
    double *dm,
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
    int& iLLE,
    int& iTs,
    LPT_UMFPACK flux_umf,
    double Ts_TOL,
    int MAX_ITERS_Ts,
    const bool debug,
    OFstream& os
)
{
    int i, bKijSet, num;
    double Ts_tmp, V1, V0; 
    double Cp1, Cp0, Cv, dVdT, G, am, bm, G_only, V_tmp, CvIG;
    double rho1, rho0, cond1, cond0, vis1, vis0;
    double J1, J0, qs, TNum, TDen, Ts_err;
    double *kij_tmp; 
    double *lnphi1; double *lnphi0; double *dlnphi1; double *dlnphi0; double *Dij1; double *Dij0; 
    double *h1; double *h0; 
    double *H1; double *H0; double *CpIG; double *Hdep1; double *Hdep0; double *Vpar;

    _NEW_(kij_tmp, double, n*n);
    _NEW_(lnphi1, double, n);
    _NEW_(lnphi0, double, n);
    _NEW_(dlnphi1, double, n*n);
    _NEW_(dlnphi0, double, n*n);
    _NEW_(Dij1, double, n*n);
    _NEW_(Dij0, double, n*n);
    _NEW_(h1, double, n);
    _NEW_(h0, double, n);
    _NEW_(H1, double, n);
    _NEW_(H0, double, n);
    _NEW_(Hdep1, double, n);
    _NEW_(Hdep0, double, n);
    _NEW_(CpIG, double, n);
    _NEW_(Vpar, double, n);

    Ts_tmp = Ts;

    for(iTs=0; iTs<MAX_ITERS_Ts; iTs++)
    {
        //calculate BIPs for current iteration interface temperature Ts_tmp
        bKijSet = 1;
        calculate_kij_(&bKijSet,&Ts_tmp,&n,Pc,Tc,w,tk,kij_tmp);

        iLLE = my2_gsl_find_transport_LLE(n,Pc,Tc,Vc,w,MW,tk,Tb,SG,H8,n_flux_type,dn1,dn0,P,Ts_tmp,xeff1,xeff0,xs1,xs0,lnphi1,lnphi0,Dij1,Dij0,H1,H0,flux_m_1,flux_m_0,flux_umf);

        //calculate fluid properties on both sides of interface
        //convert the mole fractions to mass fractions
        x2y(n,MW,xs1,ys1);
        x2y(n,MW,xs0,ys0);

        //phase 1
        thermo_properties_(&P,&Ts_tmp,&n,Pc,Tc,w,MW,xs1,Tb,SG,H8,tk,&V1,&Cp1,&Cv,CpIG,h1,Hdep1,Vpar,&dVdT,&G,lnphi1,&am,&bm,&G_only);

        fugacities_n_its_derivatives_(&P,&Ts_tmp,&n,Pc,Tc,w,xs1,tk,coef_ab,lnphi1,dlnphi1,&V_tmp);

        CvIG = calc_CvIG_from_CpIG(n,xs1,CpIG);

        vis_n_cond_(&P,&Ts_tmp,&n,Pc,Tc,Vc,w,MW,k,dm,xs1,&CvIG,&V1,&cond1,&vis1);

        rho1 = calc_rho_from_V(n,xs1,MW,V1);

        new_tlsm_diffusion_krishna_model_(&P,&Ts_tmp,&n,Pc,Tc,Vc,w,tk,coef_ab,MW,xs1,Dij1);

        //phase 0
        thermo_properties_(&P,&Ts_tmp,&n,Pc,Tc,w,MW,xs0,Tb,SG,H8,tk,&V0,&Cp0,&Cv,CpIG,h0,Hdep0,Vpar,&dVdT,&G,lnphi0,&am,&bm,&G_only);

        fugacities_n_its_derivatives_(&P,&Ts_tmp,&n,Pc,Tc,w,xs0,tk,coef_ab,lnphi0,dlnphi0,&V_tmp);

        CvIG = calc_CvIG_from_CpIG(n,xs0,CpIG);

        vis_n_cond_(&P,&Ts_tmp,&n,Pc,Tc,Vc,w,MW,k,dm,xs0,&CvIG,&V0,&cond0,&vis0);

        rho0 = calc_rho_from_V(n,xs0,MW,V0);

        new_tlsm_diffusion_krishna_model_(&P,&Ts_tmp,&n,Pc,Tc,Vc,w,tk,coef_ab,MW,xs0,Dij0);

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
            JsTot = JsTot/(double)(num);
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
        if(Tden < SMALL) TDen += SMALL;
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

    _DELETE_(kij_tmp);
    _DELETE_(lnphi1);
    _DELETE_(lnphi0);
    _DELETE_(dlnphi1);
    _DELETE_(dlnphi0);
    _DELETE_(Dij1);
    _DELETE_(Dij0);
    _DELETE_(h1);
    _DELETE_(h0);
    _DELETE_(H1);
    _DELETE_(H0);
    _DELETE_(Hdep1);
    _DELETE_(Hdep0);
    _DELETE_(CpIG);
    _DELETE_(Vpar);
}


void calc_Js_Ys
(
    const fvMesh& mesh,
    const labelListList& cellStencil,
    const List<List<scalar> >& x1_flatFld,
    const List<List<scalar> >& x0_flatFld,
    const List<scalar>& alpha1_flatFld,
    const volScalarField& alpha1,
    const volScalarField& rho1,
    const volScalarField& rho0,
    const List<vector>& C_ph1_flatFld,
    const List<vector>& C_ph0_flatFld,
    const volVectorField& C_intfc,
    const volScalarField& A_intfc,
    const volVectorField& nHat,
    PtrList<volScalarField>& Ys1,
    PtrList<volScalarField>& Ys0,
    PtrList<volScalarField>& xs1,
    PtrList<volScalarField>& xs0,
    PtrList<volScalarField>& Js1,
    PtrList<volScalarField>& Js0,    
    const label& nSpecies,
    const scalar& ALPHA_2PH_MIN,
    const scalar& A_INTFC_2PH_MIN,
    const bool debug,
    OFstream& os
)
{
    int n = nSpecies;
    scalar alpha1_cellI;
    vector nf; vector C_intfc_cellI;
    labelList curCellsAll = cellStencil[0];
    scalar dn1; scalar dn0;
    List<scalar> xeff1(n); List<scalar> xeff0(n);    
    scalar intfcGradi_cellI;

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
            << "-------------------------------------------------------------------------" << nl
            << "Interfacial Species Fluxes Calculation with Transport-LLE Constraints" << nl
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
        alpha1_cellI = alpha1Cells[cellI];
        A_intfc_cellI = A_intfcCells[cellI];
        if(debug)
        {
            os<< "-------------------------------------------------------------------------" << nl
                << "Cell: " << cellI << "  alpha1 = " << alpha1_cellI << nl
                << "-------------------------------------------------------------------------" << endl;
        }

        if(alpha1_cellI > ALPHA_2PH_MIN && alpha1_cellI < ALPHA_2PH_MAX && A_intfc_cellI > A_INTFC_2PH_MIN)
        {
            nf = nHatCells[cellI];
            C_intfc_cellI = C_intfcCells[cellI];            
            curCellsAll = cellStencil[cellI];

            if(debug)
            {
                os<< "nf = " << nf << "  C_intfc = " << C_intfc_cellI << "  A_intfc = " << A_intfc_cellI << nl
                    << "C_ph1 = " << C_ph1_flatFld[cellI] << "  C_ph0 = " << C_ph0_flatFld[cellI] << nl
                    << "rho1 = " << rho1Cells[cellI] << "rho0 = " << rho0Cells[cellI] << endl;
            }

            //phase-1
            //ensure nf direction is into the phase
            calc_cell_intfcGrad_coeffs(mesh, cellI, nf, C_intfc_cellI, x1_flatFld, T1_flatFld, alpha1_flatFld, C_ph1_flatFld, curCellsAll, nSpecies, ALPHA_2PH_MIN, 1, dn1, xeff1, Teff1, debug, os);            

            //phase-0
            //ensure nf direction is into the phase
            //then reverse nf again for Js0 calculation            
            calc_cell_intfcGrad_coeffs(mesh, cellI, -nf, C_intfc_cellI, x0_flatFld, T0_flatFld, alpha1_flatFld, C_ph0_flatFld, curCellsAll, nSpecies, ALPHA_2PH_MIN, 0, dn0, xeff0, Teff0, debug, os);

            if(curCell_had_intfc)
            {
                curTs = curTs_old;
            }
            else
            {
                curTs = 0.5*(curT1 + curT0);
            }

            calc_intfc_transLLE(curP, curTs, n, Pc, Tc, Vc, w, MW, tk, coef_ab, Tb, SG, H8, k, dm, dn1, dn0, xeff1, xeff0, Teff1, Teff0, xs1, xs0, ys1, ys0, JsTot, flux_m_1, flux_m_0, iLLE, iTs, flux_umf, Ts_TOL, MAX_ITERS_Ts, debug, os); 
        }
        else
        {
            for(label i=0; i<nSpecies; i++)
            { 
                Js0[i].internalField()[cellI] = 0;
                Js1[i].internalField()[cellI] = 0;
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

        if(debug)
        {
            os<< "Face: " << faceI << nl
                << "own: " << faceOwn << "  alpha1Own = " << alpha1Own << "  nei: " << faceNei << "  alpha1Nei = " << alpha1Nei << endl;
        }

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
    const bool debug,
    OFstream& os
)
{
    scalar ALPHA_2PH_MAX = 1 - ALPHA_2PH_MIN;
    scalar dt = mesh.time().deltaTValue();
    const scalarField& V = mesh.V();

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();

    const scalarField& alpha1Cells = alpha1.internalField();    
    const scalarField& rho1Cells = rho1.internalField();
    const scalarField& rho0Cells = rho0.internalField();

    if(debug)
    {
        os<< "-------------------------------------------------------------------------" << nl
          << "Interfacial Species Transfer Source Terms Calculation" << nl
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
            os<< "Cell: " << cellI << "  alpha1 = " << alpha1_cellI << endl;
        }

        scalar mS1Tot_cellI = 0;
        scalar mS1i_cellI = 0;
        //scalar mS0Tot_cellI = 0;
        //List<scalar> mS1_cellI(nSpecies);
        //List<scalar> mS0_cellI(nSpecies);
        scalar limiter = 1;
        scalar limiter_min = 1;
        scalar max_mS;
        scalar max_mSi;

        if(alpha1_cellI > ALPHA_2PH_MIN && alpha1_cellI < ALPHA_2PH_MAX)
        {
            mS1Tot_cellI = Js0[0].internalField()[cellI] - Js1[0].internalField()[cellI];
            mS1Tot_cellI /= (Ys1[0].internalField()[cellI] - Ys0[0].internalField()[cellI]);
            mS1Tot_cellI /= V_cellI;

            if(mS1Tot_cellI > 0)
            {
                max_mS = rho0_cellI*(1 - alpha1_cellI)/dt;
                mS1Tot_cellI = min(max_mS, mS1Tot_cellI);
            }
            else
            {
                max_mS = rho1_cellI*alpha1_cellI/dt;
                mS1Tot_cellI = -min(max_mS, -mS1Tot_cellI);
            }

            
            for(label i=0; i<nSpecies; i++)
            {                                
                mS1i_cellI = mS1Tot_cellI*Ys1[i].internalField()[cellI] + Js1[i].internalField()[cellI]/V_cellI;
                
                if(mS1i_cellI > 0)
                {
                    max_mSi = C0[i].internalField()[cellI]/dt;

                    if(mag(mS1i_cellI) > max_mSi)
                    {
                        limiter = max_mSi/mag(mS1i_cellI);
                        limiter_min = min(limiter, limiter_min);                        
                    }                    
                }
                else if(mS1i_cellI < 0)
                {
                    max_mSi = C1[i].internalField()[cellI]/dt;

                    if(mag(mS1i_cellI) > max_mSi)
                    {
                        limiter = max_mSi/mag(mS1i_cellI);
                        limiter_min = min(limiter, limiter_min);
                    }
                }
                else
                {
                    limiter = 1;
                    limiter_min = min(limiter, limiter_min);
                }
            }

            for(label i=0; i<nSpecies; i++)
            {
                mS1[i].internalField()[cellI] = limiter_min*(mS1Tot_cellI*Ys1[i].internalField()[cellI] + Js1[i].internalField()[cellI]/V_cellI);
                mS0[i].internalField()[cellI] = -mS1[i].internalField()[cellI];
            }

            mS1Tot.internalField()[cellI] = limiter_min*mS1Tot_cellI;
            mS0Tot.internalField()[cellI] = -limiter_min*mS1Tot_cellI;

            alphaS1.internalField()[cellI] = limiter_min*mS1Tot_cellI/rho1Cells[cellI];
            alphaS0.internalField()[cellI] = -limiter_min*mS1Tot_cellI/rho0Cells[cellI];
        }
        else
        {
            for(label i=0; i<nSpecies; i++)
            {
                mS1[i].internalField()[cellI] = 0;
                mS0[i].internalField()[cellI] = 0;
            }

            mS1Tot.internalField()[cellI] = 0;
            mS0Tot.internalField()[cellI] = 0;

            alphaS1.internalField()[cellI] = 0;
            alphaS0.internalField()[cellI] = 0;
        }

        if(debug)
        {
            os<< "-----------------------------------------------------------------------------------" << endl;
            for(label i=0; i<nSpecies; i++)
            {
                os<< "C1[" << i << "] = " << C1[i].internalField()[cellI] << "  C0[" << i << "] = " << C0[i].internalField()[cellI] << nl
                    << "Js1[" << i << "] = " << Js1[i].internalField()[cellI] << "  Js0[" << i << "] = " << Js0[i].internalField()[cellI] << nl 
                    << "mS1[" << i << "] = " << mS1[i].internalField()[cellI] << "  mS0[" << i << "] = " << mS0[i].internalField()[cellI] 
                    << endl;
            }
            os<< "mS1Tot = " << mS1Tot.internalField()[cellI] << "  mS0Tot = " << mS0Tot.internalField()[cellI] << nl
                << "rho1 = " << rho1Cells[cellI] << "  rho0 = " << rho0Cells[cellI] << nl
                << "alphaS1 = " << alphaS1.internalField()[cellI] << "  alphaS0 = " << alphaS0.internalField()[cellI] << nl
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
        //scalar rho1Own  = rho1Cells[faceOwn];
        scalar rho0Own  = rho0Cells[faceOwn];
        //scalar rho1Nei = rho1Cells[faceNei];
        scalar rho0Nei = rho0Cells[faceNei];

        if(debug)
        {
            os<< "Face: " << faceI << nl
                << "own: " << faceOwn << "  alpha1Own = " << alpha1Own << "  nei: " << faceNei << "  alpha1Nei = " << alpha1Nei
                << endl;
        }

        scalar mS1Tot_cellI = 0;
        scalar mS1Tot_own = 0;
        scalar mS1Tot_nei = 0;                       
        scalar max_mS;

        if(alpha1Own >= ALPHA_2PH_MAX && alpha1Nei <= ALPHA_2PH_MIN)
        {
            mS1Tot_cellI = Js0[0].internalField()[faceOwn] - Js1[0].internalField()[faceOwn];
            mS1Tot_cellI /= (Ys1[0].internalField()[faceOwn] - Ys0[0].internalField()[faceOwn]);

            if(mS1Tot_cellI > 0)
            {
                max_mS = rho0Own*(1 - alpha1Own)*VOwn/dt;
                if(mS1Tot_cellI > max_mS)
                {
                    mS1Tot_own = max_mS;                    
                }
                mS1Tot_nei = mS1Tot_cellI - mS1Tot_own;            
                
                for(label i=0; i<nSpecies; i++)
                {
                    mS1[i].internalField()[faceOwn] = mS1Tot_own*Ys1[i].internalField()[faceOwn] + Js1[i].internalField()[faceOwn];
                    mS1[i].internalField()[faceNei] = mS1Tot_nei*Ys1[i].internalField()[faceOwn];
                    mS1[i].internalField()[faceOwn] /= VOwn;
                    mS1[i].internalField()[faceNei] /= VNei;

                    mS0[i].internalField()[faceOwn] = -mS1Tot_own*Y0[i].internalField()[faceOwn];
                    mS0[i].internalField()[faceNei] = -mS1Tot_nei*Ys0[i].internalField()[faceOwn] - Js0[i].internalField()[faceOwn] - mS0[i].internalField()[faceOwn];
                    mS0[i].internalField()[faceOwn] /= VOwn;
                    mS0[i].internalField()[faceNei] /= VNei;
                }
            }
            else
            {
                mS1Tot_own = mS1Tot_cellI;
                mS1Tot_nei = mS1Tot_cellI - mS1Tot_own;            

                for(label i=0; i<nSpecies; i++)
                {
                    mS1[i].internalField()[faceOwn] = mS1Tot_own*Ys1[i].internalField()[faceOwn] + Js1[i].internalField()[faceOwn];
                    mS1[i].internalField()[faceNei] = 0;
                    mS1[i].internalField()[faceOwn] /= VOwn;
                    mS1[i].internalField()[faceNei] /= VNei;

                    mS0[i].internalField()[faceOwn] = -mS1Tot_own*Ys0[i].internalField()[faceOwn];
                    mS0[i].internalField()[faceNei] = -Js0[i].internalField()[faceOwn];
                    mS0[i].internalField()[faceOwn] /= VOwn;
                    mS0[i].internalField()[faceNei] /= VNei;
                }
            }            

            mS1Tot.internalField()[faceOwn] = mS1Tot_own/VOwn;
            mS0Tot.internalField()[faceOwn] = -mS1Tot_own/VOwn;

            mS1Tot.internalField()[faceNei] = mS1Tot_nei/VNei;
            mS0Tot.internalField()[faceNei] = -mS1Tot_nei/VNei;      
        }

        if(alpha1Own <= ALPHA_2PH_MIN && alpha1Nei >= ALPHA_2PH_MAX)
        {            
            mS1Tot_cellI = Js0[0].internalField()[faceNei] - Js1[0].internalField()[faceNei];
            mS1Tot_cellI /= (Ys1[0].internalField()[faceNei] - Ys0[0].internalField()[faceNei]);

            if(mS1Tot_cellI > 0)
            {
                max_mS = rho0Nei*(1 - alpha1Nei)*VNei/dt;
                if(mS1Tot_cellI > max_mS)
                {
                    mS1Tot_nei = max_mS;                    
                }                
                mS1Tot_own = mS1Tot_cellI - mS1Tot_nei;

                for(label i=0; i<nSpecies; i++)
                {
                    mS1[i].internalField()[faceNei] = mS1Tot_nei*Ys1[i].internalField()[faceNei] + Js1[i].internalField()[faceNei];
                    mS1[i].internalField()[faceOwn] = mS1Tot_own*Ys1[i].internalField()[faceNei];
                    mS1[i].internalField()[faceOwn] /= VOwn;
                    mS1[i].internalField()[faceNei] /= VNei;

                    mS0[i].internalField()[faceNei] = -mS1Tot_nei*Y0[i].internalField()[faceNei];
                    mS0[i].internalField()[faceOwn] = -mS1Tot_own*Ys0[i].internalField()[faceNei] - Js0[i].internalField()[faceNei] - mS0[i].internalField()[faceNei];
                    mS0[i].internalField()[faceOwn] /= VOwn;
                    mS0[i].internalField()[faceNei] /= VNei;
                }
            }
            else
            {
                mS1Tot_nei = mS1Tot_cellI;
                mS1Tot_own = mS1Tot_cellI - mS1Tot_nei;

                for(label i=0; i<nSpecies; i++)
                {
                    mS1[i].internalField()[faceNei] = mS1Tot_nei*Ys1[i].internalField()[faceNei] + Js1[i].internalField()[faceNei];
                    mS1[i].internalField()[faceOwn] = 0;
                    mS1[i].internalField()[faceOwn] /= VOwn;
                    mS1[i].internalField()[faceNei] /= VNei;

                    mS0[i].internalField()[faceNei] = -mS1Tot_nei*Ys0[i].internalField()[faceNei];
                    mS0[i].internalField()[faceOwn] = -Js0[i].internalField()[faceNei];
                    mS0[i].internalField()[faceOwn] /= VOwn;
                    mS0[i].internalField()[faceNei] /= VNei;
                }
            }                        

            mS1Tot.internalField()[faceOwn] = mS1Tot_own/VOwn;
            mS0Tot.internalField()[faceOwn] = -mS1Tot_own/VOwn;

            mS1Tot.internalField()[faceNei] = mS1Tot_nei/VNei;
            mS0Tot.internalField()[faceNei] = -mS1Tot_nei/VNei;      
        }

        alphaS1.internalField()[faceOwn] = mS1Tot.internalField()[faceOwn]/rho1Cells[faceOwn];
        alphaS0.internalField()[faceOwn] = mS0Tot.internalField()[faceOwn]/rho0Cells[faceOwn];

        alphaS1.internalField()[faceNei] = mS1Tot.internalField()[faceNei]/rho1Cells[faceNei];
        alphaS0.internalField()[faceNei] = mS0Tot.internalField()[faceNei]/rho0Cells[faceNei];

        if(debug)
        {
            os<< "-----------------------------------------------------------------------------------" << endl;
            for(label i=0; i<nSpecies; i++)
            {
                os<< "C1[" << i << "] = " << C1[i].internalField()[faceOwn] << "  C0[" << i << "] = " << C0[i].internalField()[faceOwn] << nl
                    << "Js1[" << i << "] = " << Js1[i].internalField()[faceOwn] << "  Js0[" << i << "] = " << Js0[i].internalField()[faceOwn] << nl 
                    << "mS1[" << i << "] = " << mS1[i].internalField()[faceOwn] << "  mS0[" << i << "] = " << mS0[i].internalField()[faceOwn] 
                    << endl;
            }
            os<< "mS1Tot = " << mS1Tot.internalField()[faceOwn] << "  mS0Tot = " << mS0Tot.internalField()[faceOwn] << nl
                << "rho1 = " << rho1Cells[faceOwn] << "  rho0 = " << rho0Cells[faceOwn] << nl
                << "alphaS1 = " << alphaS1.internalField()[faceOwn] << "  alphaS0 = " << alphaS0.internalField()[faceOwn] << endl;
            for(label i=0; i<nSpecies; i++)
            {
                os<< "C1[" << i << "] = " << C1[i].internalField()[faceNei] << "  C0[" << i << "] = " << C0[i].internalField()[faceNei] << nl
                    << "Js1[" << i << "] = " << Js1[i].internalField()[faceNei] << "  Js0[" << i << "] = " << Js0[i].internalField()[faceNei] << nl
                    << "mS1[" << i << "] = " << mS1[i].internalField()[faceNei] << "  mS0[" << i << "] = " << mS0[i].internalField()[faceNei] 
                    << endl;
            }
            os<< "mS1Tot = " << mS1Tot.internalField()[faceNei] << "  mS0Tot = " << mS0Tot.internalField()[faceNei] << nl
                << "rho1 = " << rho1Cells[faceNei] << "  rho0 = " << rho0Cells[faceNei] << nl
                << "alphaS1 = " << alphaS1.internalField()[faceNei] << "  alphaS0 = " << alphaS0.internalField()[faceNei] << endl;
            os<< "-----------------------------------------------------------------------------------" << nl
                << endl;
        }
    }

    //boundary faces
    if(debug)
    {
        os<< "-------------------------------------------------------------------------" << nl
          << "Boundary coupled faces" << nl
          << "-------------------------------------------------------------------------" << nl
          << endl;
    }
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const wordList& patchNames = patches.names();

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
            label faceI = pp.start();
            label bndFaceI = pp.start() - mesh.nInternalFaces();
            forAll(pp, fcI)
            {
                const label& faceOwn = own[faceI];
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
            label bndFaceI = pp.start() - mesh.nInternalFaces();
            const scalarField& alpha1NeiFld = alpha1.boundaryField()[patchI].patchNeighbourField();
            //const scalarField& rho1NeiFld = rho1.boundaryField()[patchI].patchNeighbourField();            
            const scalarField& rho0NeiFld = rho0.boundaryField()[patchI].patchNeighbourField();            

            forAll(pmS1Tot, fcI)
            {
                label faceOwn = own[faceI];                
                scalar alpha1Own = alpha1Cells[faceOwn];
                scalar alpha1Nei = alpha1NeiFld[fcI];
                scalar VOwn = V[faceOwn];
                scalar VNei = VNeiFld[bndFaceI];
                //scalar rho1Own  = rho1Cells[faceOwn];
                scalar rho0Own  = rho0Cells[faceOwn];
                //scalar rho1Nei = rho1NeiFld[fcI];
                scalar rho0Nei = rho0NeiFld[fcI];

                if(debug)
                {
                    os<< "Face: " << faceI << "  patch face index: " << fcI << nl
                        << "own: " << faceOwn << "  alpha1Own = " << alpha1Own << "  alpha1Nei = " << alpha1Nei << endl;
                }

                scalar mS1Tot_cellI = 0;
                scalar mS1Tot_own = 0;
                scalar mS1Tot_nei = 0;                       
                scalar max_mS;

                if(alpha1Own >= ALPHA_2PH_MAX && alpha1Nei <= ALPHA_2PH_MIN)
                {
                    mS1Tot_cellI = Js0[0].internalField()[faceOwn] - Js1[0].internalField()[faceOwn];
                    mS1Tot_cellI /= (Ys1[0].internalField()[faceOwn] - Ys0[0].internalField()[faceOwn]);

                    if(mS1Tot_cellI > 0)
                    {
                        max_mS = rho0Own*(1 - alpha1Own)*VOwn/dt;
                        if(mS1Tot_cellI > max_mS)
                        {
                            mS1Tot_own = max_mS;                    
                        }
                        mS1Tot_nei = mS1Tot_cellI - mS1Tot_own;

                        for(label i=0; i<nSpecies; i++)
                        {
                            mS1[i].internalField()[faceOwn] = mS1Tot_own*Ys1[i].internalField()[faceOwn] + Js1[i].internalField()[faceOwn];
                            mS1[i].internalField()[faceOwn] /= VOwn;

                            mS0[i].internalField()[faceOwn] = -mS1Tot_own*Y0[i].internalField()[faceOwn];
                            mS0[i].internalField()[faceOwn] /= VOwn;
                        }
                    }
                    else
                    {
                        mS1Tot_own = mS1Tot_cellI;
                        mS1Tot_nei = mS1Tot_cellI - mS1Tot_own;

                        for(label i=0; i<nSpecies; i++)
                        {
                            mS1[i].internalField()[faceOwn] = mS1Tot_own*Ys1[i].internalField()[faceOwn] + Js1[i].internalField()[faceOwn];
                            mS1[i].internalField()[faceOwn] /= VOwn;

                            mS0[i].internalField()[faceOwn] = -mS1Tot_own*Ys0[i].internalField()[faceOwn];
                            mS0[i].internalField()[faceOwn] /= VOwn;
                        }
                    }                    
                    
                    mS1Tot.internalField()[faceOwn] = mS1Tot_own/VOwn;
                    mS0Tot.internalField()[faceOwn] = -mS1Tot_own/VOwn;
                }

                if(alpha1Own <= ALPHA_2PH_MIN && alpha1Nei >= ALPHA_2PH_MAX)
                {            
                    mS1Tot_cellI = Js0[0].internalField()[faceOwn] - Js1[0].internalField()[faceOwn];
                    mS1Tot_cellI /= (Ys1[0].internalField()[faceOwn] - Ys0[0].internalField()[faceOwn]);

                    if(mS1Tot_cellI > 0)
                    {
                        max_mS = rho0Nei*(1 - alpha1Nei)*VNei/dt;
                        if(mS1Tot_cellI > max_mS)
                        {
                            mS1Tot_nei = max_mS;                    
                        }                
                        mS1Tot_own = mS1Tot_cellI - mS1Tot_nei;

                        for(label i=0; i<nSpecies; i++)
                        {
                            const scalarField& Y0iNeiFld = Y0[i].boundaryField()[patchI].patchNeighbourField();

                            mS1[i].internalField()[faceOwn] = mS1Tot_own*Ys1[i].internalField()[faceOwn];
                            mS1[i].internalField()[faceOwn] /= VOwn;
                            
                            mS0[i].internalField()[faceOwn] = -mS1Tot_own*Ys0[i].internalField()[faceOwn] - Js0[i].internalField()[faceOwn] + mS1Tot_nei*Y0iNeiFld[fcI];
                            mS0[i].internalField()[faceOwn] /= VOwn;
                        }
                    }
                    else
                    {
                        mS1Tot_nei = mS1Tot_cellI;
                        mS1Tot_own = mS1Tot_cellI - mS1Tot_nei;

                        for(label i=0; i<nSpecies; i++)
                        {
                            mS1[i].internalField()[faceOwn] = 0;
                            mS1[i].internalField()[faceOwn] /= VOwn;

                            mS0[i].internalField()[faceOwn] = -Js0[0].internalField()[faceOwn];
                            mS0[i].internalField()[faceOwn] /= VOwn;
                        }
                    }                                        

                    mS1Tot.internalField()[faceOwn] = mS1Tot_own/VOwn;
                    mS0Tot.internalField()[faceOwn] = -mS1Tot_own/VOwn;
                }

                alphaS1.internalField()[faceOwn] = mS1Tot.internalField()[faceOwn]/rho1Cells[faceOwn];
                alphaS0.internalField()[faceOwn] = mS0Tot.internalField()[faceOwn]/rho0Cells[faceOwn];

                if(debug)
                {
                    os<< "-----------------------------------------------------------------------------------" << endl;
                    for(label i=0; i<nSpecies; i++)
                    {
                        os<< "C1[" << i << "] = " << C1[i].internalField()[faceOwn] << "  C0[" << i << "] = " << C0[i].internalField()[faceOwn] << nl
                            << "Js1[" << i << "] = " << Js1[i].internalField()[faceOwn] << "  Js0[" << i << "] = " << Js0[i].internalField()[faceOwn] << nl 
                            << "mS1[" << i << "] = " << mS1[i].internalField()[faceOwn] << "  mS0[" << i << "] = " << mS0[i].internalField()[faceOwn] 
                            << endl;
                    }
                    os<< "mS1Tot = " << mS1Tot.internalField()[faceOwn] << "  mS0Tot = " << mS0Tot.internalField()[faceOwn] << nl
                        << "rho1 = " << rho1Cells[faceOwn] << "  rho0 = " << rho0Cells[faceOwn] << nl
                        << "alphaS1 = " << alphaS1.internalField()[faceOwn] << "  alphaS0 = " << alphaS0.internalField()[faceOwn] << endl;
                    /*
                    for(label i=0; i<nSpecies; i++)
                    {
                        os<< "C1[" << i << "] = " << C1[i].internalField()[faceNei] << "  C0[" << i << "] = " << C0[i].internalField()[faceNei] << nl
                            << "Js1[" << i << "] = " << Js1[i].internalField()[faceNei] << "  Js0[" << i << "] = " << Js0[i].internalField()[faceNei] << nl
                            << "mS1[" << i << "] = " << mS1[i].internalField()[faceNei] << "  mS0[" << i << "] = " << mS0[i].internalField()[faceNei] 
                            << endl;
                    }
                    os<< "mS1Tot = " << mS1Tot.internalField()[faceNei] << "  mS0Tot = " << mS0Tot.internalField()[faceNei] << nl
                        << "rho1 = " << rho1Cells[faceNei] << "  rho0 = " << rho0Cells[faceNei] << nl
                        << "alphaS1 = " << alphaS1.internalField()[faceNei] << "  alphaS0 = " << alphaS0.internalField()[faceNei] << endl;
                        */
                    os<< "-----------------------------------------------------------------------------------" << nl
                        << endl;
                }

                faceI++;
                bndFaceI++;
            }
        }
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace plicFuncs

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
