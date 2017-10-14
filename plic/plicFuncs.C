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
    bool debug,
    OFstream& os
)
{
    scalar cosThetaMax = -1;
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
    bool debug,
    OFstream& os
)
{
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
    label C1_lbl = findCellInFaceDir(curCells,C,Cp,nf,-1,debug,os);
    vector C1 = C[C1_lbl];
    label C2_lbl = findCellInFaceOrthDir(curCells,C,Cp,C1,nf,C1_lbl,debug,os);

    if(debug)
    {
        os<< "Cell " << curCell_lbl << nl
            << "Cell phase centroid: " << Cp << "  nf: " << nf << nl
            << "C1_lbl: " << C1_lbl << "  C2_lbl: " << C2_lbl
            << endl;
    }
    
    scalar Yp = Y[curCell_lbl];
    scalar Y1 = Y[C1_lbl];
    vector C2 = C[C2_lbl];
    scalar Y2 = Y[C2_lbl];
    scalar Y2_1 = Y2;

    vector t1 = C1 - Cp;
    magt1 = mag(t1);
    if(magt1 < SMALL)
    {
        magt1 += SMALL;
    }
    vector t2 = C2 - Cp;
    magt2 = mag(t2);
    if(magt2 < SMALL)
    {
        magt2 += SMALL;
    }    
    scalar magt1t2 = magt1*magt2;
    if(magt1t2 < SMALL)
    {
        magt1t2 += SMALL;
    }    

    scalar theta2_sign = -((nf ^ t1) & (nf ^ t2))/magt1t2;

    scalar theta1 = acos((nf & t1)/magt1);
    scalar theta2;
    scalar theta2_1;
    if(theta2_sign >= 0)
    {
        theta2 = acos((nf & t2)/magt2);
        theta2_1 = theta2;
    }
    else
    {
        theta2 = constant::mathematical::pi - acos((nf & t2)/magt2);
        theta2_1 = acos((nf & t2)/magt2);
        Y2 = Yp + (Yp - Y2_1);
    }

    if(debug)
    {
        os<< "t1 = " << t1 << "  t2 = " << t2 << nl
            << "mag(t1) = " << magt1 << "  mag(t2) = " << magt2 << nl
            << "theta1 = " << theta1 << "  theta2_1 = " << theta2_1 << nl
            << "theta2_sign = " << theta2_sign << "  theta2 = " << theta2 << nl
            << "Yp = " << Yp << "  Y1 = " << Y1 << "  Y2 = " << Y2 << "  Y2_1 = " << Y2_1 
            << endl;
    }

    scalar sintheta12 = sin(theta1 + theta2);
    if(sintheta12 < SMALL)
    {
        sintheta12 += SMALL;
    }

    alpha = sin(theta2)/sintheta12;
        
    beta = sin(theta1)/sintheta12;            
    //beta = -sin(theta1)/sin(theta1 + theta2);
    
    //d = alpha*Y1/magt1 + beta*Y2/magt2;
    d = beta*Y2/magt2;

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
    scalar Yp = Y[curCell_lbl];
    label C1_lbl = findCellInFaceDir(curCells,C,Cp,nf,-1,debug,os);
    vector C1 = C[C1_lbl];
    scalar Y1 = Y[C1_lbl];
    label C2_lbl = findCellInFaceOrthDir(curCells,C,Cp,C1,nf,C1_lbl,debug,os);

    if(debug)
    {
        os<< "Cell " << curCell_lbl << nl
            << "Cell phase centroid: " << Cp << "  nf: " << nf << nl
            << "C1_lbl: " << C1_lbl << "  C2_lbl: " << C2_lbl
            << endl;
    }

    vector C2 = C[C2_lbl];
    scalar Y2 = Y[C2_lbl];

    vector t1 = C1 - Cp;
    scalar magt1 = mag(t1);
    if(magt1 < SMALL)
    {
        magt1 += SMALL;
    }
    vector t2 = C2 - Cp;
    scalar magt2 = mag(t2);
    magt2 = mag(t2);
    if(magt2 < SMALL)
    {
        magt2 += SMALL;
    }    
    scalar magt1t2 = magt1*magt2;
    if(magt1t2 < SMALL)
    {
        magt1t2 += SMALL;
    }

    scalar theta2_sign = -((nf ^ t1) & (nf ^ t2))/magt1t2;

    scalar theta1 = acos((nf & t1)/magt1);
    scalar theta2;
    scalar theta2_1;
    if(theta2_sign >= 0)
    {
        theta2 = acos((nf & t2)/magt2);
        theta2_1 = theta2;
    }
    else
    {
        theta2 = constant::mathematical::pi - acos((nf & t2)/magt2);
        theta2_1 = acos((nf & t2)/magt2);
    }

    if(debug)
    {
        os<< "t1 = " << t1 << "  t2 = " << t2 << nl
            << "mag(t1) = " << magt1 << "  mag(t2) = " << magt2 << nl
            << "theta1 = " << theta1 << "  theta2_1 = " << theta2_1 << nl
            << "theta2_sign = " << theta2_sign << "  theta2 = " << theta2 << nl            
            << "Y1 = " << Y1 << "  Y2 = " << Y2 
            << endl;
    }

    scalar sintheta12 = sin(theta1 + theta2);
    if(sintheta12 < SMALL)
    {
        sintheta12 += SMALL;
    }

    scalar alpha = sin(theta2)/sintheta12;
    scalar beta;
    if(theta2_sign >= 0)
    {
        beta = sin(theta1)/sintheta12;
    }
    else
    {
        beta = -sin(theta1)/sintheta12;
    }

    cellGrad = alpha*(Y1 - Yp)/magt1 + beta*(Y2 - Yp)/magt2;

    if(debug)
    {
        os<< "alpha = " << alpha << "  beta = " << beta << nl
            << "cell grad = " << cellGrad << nl
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
            diffFlux_Y1i[faceI] = 0;
        }
        else if(curPhaseState == 1)
        {
            diffFlux_Y0i[faceI] = 0;
            diffFlux_Y1i[faceI] = -rho1f[faceI]*D1fi[faceI]*curMagSf_ph1*gradf_Y1i[faceI];
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
        calc_diffFlux_limiter
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
            diffFlux_limiter
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

        calc_diffFlux_limiter
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
            diffFlux_limiter
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
                }
                if(curPhaseState == 1)
                {
                    pdiffFlux_Y1i[fcI] = -prho1f[fcI]*pD1fi[fcI]*curMagSf_ph1*pgradf_Y1i[fcI];
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

                calc_diffFlux_limiter
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
                    diffFlux_limiter
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

                calc_diffFlux_limiter
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
                    diffFlux_limiter
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
                }
                if(curPhaseState == 1)
                {
                    pdiffFlux_Y1i[fcI] = -prho1f[fcI]*pD1fi[fcI]*curMagSf_ph1*pgradf_Y1i[fcI];
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
                    diffFlux_limiter
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
                    diffFlux_limiter
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
    scalar& diffFlux_limiter
)
{
    scalar maxDiffFluxOwn;
    scalar maxDiffFluxNei;
    scalar maxDiffFlux;
    diffFlux_limiter = 1;    

    if(diffFlux > 0)
    {
        maxDiffFluxOwn = 0.25*rhoOwn*alphaOwn*YOwn*VOwn/dt;
        
        maxDiffFluxNei = 0.25*rhoNei*alphaNei*(1 - YNei)*VNei/dt;

        maxDiffFlux = min(maxDiffFluxOwn, maxDiffFluxNei);
        
        if(mag(diffFlux) > maxDiffFlux)
        {
            diffFlux_limiter = maxDiffFlux/(mag(diffFlux) + SMALL);
        }        
    }
    else
    {
        maxDiffFluxNei = 0.25*rhoNei*alphaNei*YNei*VNei/dt;

        maxDiffFluxOwn = 0.25*rhoOwn*alphaOwn*(1 - YOwn)*VOwn/dt;

        maxDiffFlux = min(maxDiffFluxOwn, maxDiffFluxNei);

        if(mag(diffFlux) > maxDiffFlux)
        {
            diffFlux_limiter = maxDiffFlux/(mag(diffFlux) + SMALL);
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
    const label& Y_BOUND_ITERS_MAX
)
{
    forAll(YCells, cellI)
    {
        if(YCells[cellI] < Y_MIN)
        {            
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
        }

        if(YCells[cellI] > Y_MAX)
        {            
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
        }       
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace plicFuncs

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
