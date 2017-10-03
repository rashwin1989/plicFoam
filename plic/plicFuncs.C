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
    bool debug
)
{
    scalar cosThetaMax = -1;
    label cell_lbl = 0;

    if(debug)
    {
        Info<< "Finding cell in face normal direction" << nl 
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
            Info<< "Cell " << cellI << " in stencil: " << curCell << nl
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
                Info<< "costheta = " << cosTheta << endl;
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
        Info<< "Cell in nf direction: " << cell_lbl << nl
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
    bool debug
)
{
    scalar magCosThetaMin = 1;
    label cell_lbl = 0;

    if(debug)
    {
        Info<< "Finding cell in face normal orthogonal direction" << nl 
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
            Info<< "Cell " << cellI << " in stencil: " << curCell << nl
                << "C = " << Cp << "  C1 = " << C1 << "  Ci = " << C[curCell] << nl
                << "CC1 = " << CC1 << "  CCi" << CCi << "  CC1_CCi_dir = " << CC1_CCi_dir << nl
                << endl;
        }

        if(curCell != C1_lbl)
        {
            //CCi /= mag(CCi);
            //scalar cosTheta = nf & CCi;
            scalar cosTheta = (CC1 & CCi)/mag(CCi)/mag(CC1);

            if(debug)
            {
                Info<< "costheta = " << cosTheta << endl;
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
        Info<< "Cell in nf orthogonal direction: " << cell_lbl << nl
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
    bool debug
)
{        
    if(debug)
    {
        Info<< "Calculating cell grad weights in cell " << curCell_lbl << nl
            << endl;
    }

    scalar MAX_ALPHA_DIFF = 1 - MIN_ALPHA_DIFF;
    
    labelList curCells(curCellsAll.size());
    label n_ph = 0;

    if(debug)
    {
        Info<< "Reducing cell stencil for phase " << phaseLbl << nl
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
                Info<< "cellI: " << cellI << "  cellI_lbl: " << cellI_lbl << nl
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
                Info<< "cellI: " << cellI << "  cellI_lbl: " << cellI_lbl << nl
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
        Info<< "Cell reduced stencil" << nl
            << curCells << nl
            << endl;
    }

    // suffix 1: direction closest to nf
    // suffix 2: direction closest orthogonal to nf in 2-D
    // further improvements needed for 3-D calculation
    vector Cp = C[curCell_lbl];
    label C1_lbl = findCellInFaceDir(curCells,C,Cp,nf,-1,debug);
    vector C1 = C[C1_lbl];
    label C2_lbl = findCellInFaceOrthDir(curCells,C,Cp,C1,nf,C1_lbl,debug);

    if(debug)
    {
        Info<< "Cell " << curCell_lbl << nl
            << "Cell phase centroid: " << Cp << "  nf: " << nf << nl
            << "C1_lbl: " << C1_lbl << "  C2_lbl: " << C2_lbl
            << endl;
    }
    
    scalar Y1 = Y[C1_lbl];
    vector C2 = C[C2_lbl];
    scalar Y2 = Y[C2_lbl];

    vector t1 = C1 - Cp;
    magt1 = mag(t1);
    vector t2 = C2 - Cp;
    magt2 = mag(t2);
    scalar theta2_sign = -((nf ^ t1) & (nf ^ t2))/mag(t1)/mag(t2);

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
        Info<< "t1 = " << t1 << "  t2 = " << t2 << nl
            << "mag(t1) = " << magt1 << "  mag(t2) = " << magt2 << nl
            << "theta1 = " << theta1 << "  theta2_1 = " << theta2_1 << nl
            << "theta2_sign = " << theta2_sign << "  theta2 = " << theta2 << nl
            << "Y1 = " << Y1 << "  Y2 = " << Y2 
            << endl;
    }

    alpha = sin(theta2)/sin(theta1 + theta2);
    if(theta2_sign >= 0)
    {
        beta = sin(theta1)/sin(theta1 + theta2);
    }
    else
    {
        beta = -sin(theta1)/sin(theta1 + theta2);
    }
    //d = alpha*Y1/magt1 + beta*Y2/magt2;
    d = beta*Y2/magt2;

    if(debug)
    {
        Info<< "alpha = " << alpha << "  beta = " << beta << "  d = " << d << nl
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
    bool debug
)
{    
    if(debug)
    {
        Info<< "Calculating cell grad in cell " << curCell_lbl << nl
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
        Info<< "Cell reduced stencil" << nl
            << curCells << endl;    
    }

    // suffix 1: direction closest to nf
    // suffix 2: direction closest orthogonal to nf in 2-D
    // further improvements needed for 3-D calculation
    vector Cp = C[curCell_lbl];
    scalar Yp = Y[curCell_lbl];
    label C1_lbl = findCellInFaceDir(curCells,C,Cp,nf,-1,debug);
    vector C1 = C[C1_lbl];
    scalar Y1 = Y[C1_lbl];
    label C2_lbl = findCellInFaceOrthDir(curCells,C,Cp,C1,nf,C1_lbl,debug);

    if(debug)
    {
        Info<< "Cell " << curCell_lbl << nl
            << "Cell phase centroid: " << Cp << "  nf: " << nf << nl
            << "C1_lbl: " << C1_lbl << "  C2_lbl: " << C2_lbl
            << endl;
    }

    vector C2 = C[C2_lbl];
    scalar Y2 = Y[C2_lbl];

    vector t1 = C1 - Cp;
    scalar magt1 = mag(t1);
    vector t2 = C2 - Cp;
    scalar magt2 = mag(t2);    
    scalar theta2_sign = -((nf ^ t1) & (nf ^ t2))/mag(t1)/mag(t2);

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
        Info<< "t1 = " << t1 << "  t2 = " << t2 << nl
            << "mag(t1) = " << magt1 << "  mag(t2) = " << magt2 << nl
            << "theta1 = " << theta1 << "  theta2_1 = " << theta2_1 << nl
            << "theta2_sign = " << theta2_sign << "  theta2 = " << theta2 << nl            
            << "Y1 = " << Y1 << "  Y2 = " << Y2 
            << endl;
    }

    scalar alpha = sin(theta2)/sin(theta1 + theta2);
    scalar beta;
    if(theta2_sign >= 0)
    {
        beta = sin(theta1)/sin(theta1 + theta2);
    }
    else
    {
        beta = -sin(theta1)/sin(theta1 + theta2);
    }

    cellGrad = alpha*(Y1 - Yp)/magt1 + beta*(Y2 - Yp)/magt2;

    if(debug)
    {
        Info<< "alpha = " << alpha << "  beta = " << beta << nl
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
    bool debug
)
{
    if(debug)
    {
        Info<< "Calculating two-sided grad weights for face " << faceI << "  Phase: " << phaseLbl << nl
            << "Own: " << own << "  Nei: " << nei << nl
            << endl;
    }

    if(debug)
    {
        Info<< "Own side calculation" << nl 
            << endl;
    }
    scalar alphap = 1;
    scalar betap = 1;
    scalar magt1p = 1;
    scalar magt2p = 1;
    scalar dp = 1;
    const labelList& ownCells = diffCellStencil[own];
    calcCellGradWeights(own,nf,Y,alpha1,C,ownCells,MIN_ALPHA_DIFF,phaseLbl,alphap,betap,magt1p,magt2p,dp,debug);

    if(debug)
    {
        Info<< "Nei side calculation" << nl
            << endl;
    }
    scalar alpham = 1;
    scalar betam = 1;
    scalar magt1m = 1;
    scalar magt2m = 1;
    scalar dm = 1;
    const labelList& neiCells = diffCellStencil[nei];
    calcCellGradWeights(nei,-nf,Y,alpha1,C,neiCells,MIN_ALPHA_DIFF,phaseLbl,alpham,betam,magt1m,magt2m,dm,debug);

    scalar mup;
    scalar mum;
    if(mag(dm) < SMALL && mag(dp) < SMALL)
    {
        mup = 0.5;
    }
    else
    {
        mup = dm/(dm + dp);
    }
    mum = 1 - mup;

    wOwn = -(mup*(alphap/magt1p + betap/magt2p) + mum*alpham/magt1m);
    wNei = mum*(alpham/magt1m + betam/magt2m) + mup*alphap/magt1p;

    if(debug)
    {
        Info<< nl
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
    bool debug
)
{
    if(debug)
    {
        Info<< "Calculating face grad from weights" << nl            
            << endl;
    }

    scalar mup;
    scalar mum;
    if(mag(dm) < SMALL && mag(dp) < SMALL)
    {
        mup = 0.5;
    }
    else
    {
        mup = dm/(dm + dp);
    }    
    mum = 1 - mup;

    scalar wOwn = -(mup*(alphap/magt1p + betap/magt2p) + mum*alpham/magt1m);
    scalar wNei = mum*(alpham/magt1m + betam/magt2m) + mup*alphap/magt1p;

    faceGrad = wOwn*Yp + wNei*Ym;

    if(debug)
    {
        Info<< nl
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace plicFuncs

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
