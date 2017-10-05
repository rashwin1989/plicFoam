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

InNamespace
    Foam::plicFuncs

Description
    Functions to be used by plic class and other solver code

SourceFiles
    plicFuncs.C

\*---------------------------------------------------------------------------*/


#ifndef plicFuncs_H
#define plicFuncs_H

#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "cell.H"
#include "cellInfo.H"
#include "faceList.H"
#include "pointField.H"
#include "Plane.H"
#include "labelList.H"
#include "pointFields.H"
#include "fvMesh.H"
#include "word.H"
#include "label.H"
#include "plic.H"
#include "dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace plicFuncs functions Declaration
\*---------------------------------------------------------------------------*/

namespace plicFuncs
{
    template <class Type>
    void display_field
    (
        const GeometricField<Type, fvPatchField, volMesh>&
    );

    template <class Type>
    void write_field
    (
        const GeometricField<Type, fvPatchField, volMesh>&
    );

    template <class Type>
    void display_point_field
    (
        const GeometricField<Type, pointPatchField, pointMesh>&
    );

    template <class Type>
    void write_point_field
    (
        const GeometricField<Type, pointPatchField, pointMesh>&,
        const fvMesh&
    );

    void display_cell
    (
        const cell&,
        const faceList&,
        const pointField&
    );

    void display_cellInfo
    (
        const cellInfo&        
    );
 
    void display_face
    (
        const face&,
        const pointField&
    );
 
    void display_pointField
    (
        const pointField&
    );

    void display_Plane
    (
        const Plane&
    );

    void display_labelList
    (
        const labelList&
    );

    void write_labelList
    (
        const labelList&,
        const fvMesh&,
        const word&
    );

    void write_boolList
    (
        const boolList&,
        const fvMesh&,
        const word&
    );

    void write_stencil
    (
        const labelListList&,
        const fvMesh&,
        const word&
    );

    template<class Type>
    void write_flatFld
    (
        const List<Type>&,
        const GeometricField<Type, fvPatchField, volMesh>&
    );

    template<class Type>
    void write_surfaceField
    (
        const GeometricField<Type, fvsPatchField, surfaceMesh>&,
        const fvMesh&
    );

    template<class Type>
    void display_surfaceField
    (
        const GeometricField<Type, fvsPatchField, surfaceMesh>&,
        const fvMesh&
    );

    point centre
    (
        const face&,
        const pointField&
    );

    label findCellInFaceDir
    (
        const labelList&,
        const List<vector>&,
        const vector&,
        const vector&,
        const label&,
        bool
    );

    label findCellInFaceOrthDir
    (
        const labelList&,
        const List<vector>&,
        const vector&,
        const vector&,
        const label&,
        bool
    );    

    void calcCellGradWeights
    (
        const label&,
        const vector&,
        const List<scalar>&,
        const List<scalar>&,
        const List<vector>&,
        const labelList&,
        const scalar&,
        const label&,
        scalar&,
        scalar&,
        scalar&,
        scalar&,
        scalar&,
        bool
    );

    void calcCellGrad
    (
        const label&,
        const vector&,
        const List<scalar>&,
        const List<scalar>&,
        const List<vector>&,
        const labelList&,
        const scalar&,
        const label&,
        scalar&,
        bool
    );

    void calcTwoSidedFaceGradWeights
    (
        const label&,
        const label&,
        const label&,
        const vector&,
        const List<scalar>&,
        const List<scalar>&,
        const List<vector>&,
        const labelListList&,
        const scalar&,
        const label&,
        scalar&,
        scalar&,        
        bool
    );

    void calcFaceGradFromWeights
    (
        const scalar&,
        const scalar&,
        const scalar&,
        const scalar&,
        const scalar&,
        const scalar&,
        const scalar&,
        const scalar&,
        const scalar&,
        const scalar&,
        const scalar&,
        const scalar&,    
        scalar&, 
        bool
    );

    template<class Type>
    void makeFlatFld
    (
        const GeometricField<Type, fvPatchField, volMesh>&,
        const mapDistribute&,
        List<Type>&
    );

    void calc_2ph_gradf
    (
        const label& nSpecies,
        const plic& interface,
        const fvMesh& mesh,
        const labelListList& diffCellStencil,
        const volScalarField& Y1i,
        const volScalarField& Y0i,
        const List<scalar>& Y1i_flatFld_diff,
        const List<scalar>& Y0i_flatFld_diff,
        const List<scalar>& alpha1_flatFld_diff,
        const List<vector>& C_ph1_flatFld_diff,
        const List<vector>& C_ph0_flatFld_diff,
        surfaceScalarField& gradf_Y1i,
        surfaceScalarField& gradf_Y0i,
        const label& i,
        const scalar& MIN_ALPHA_DIFF,
        const bool debug
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "plicFuncs.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
