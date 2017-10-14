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

#include "fvCFD.H"
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
#include "syncTools.H"
#include "OFstream.H"

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
        bool, 
        OFstream&
    );

    label findCellInFaceOrthDir
    (
        const labelList&,
        const List<vector>&,
        const vector&,
        const vector&,
        const label&,
        bool, 
        OFstream&
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
        bool, 
        OFstream&
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
        bool, 
        OFstream&
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
        bool, 
        OFstream&
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
        bool, 
        OFstream&
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
    );

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
    );

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
    );

    template<class Type>
    void linearInterpolate_2ph
    (
        const GeometricField<Type, fvPatchField, volMesh>& Y,
        const fvMesh& mesh,
        const surfaceScalarField& weights,
        GeometricField<Type, fvsPatchField, surfaceMesh>& Yf
    );

    void calc_2ph_linearInterpolation_weights
    (
        const fvMesh& mesh,
        const vectorField& C_ph,
        const vectorField& Cf_ph,
        surfaceScalarField& weights
    );

    void calc_face_phaseState
    (
        const scalar& curAlpha1f_own,
        const scalar& curAlpha0f_own,
        const scalar& curAlpha1f_nei,
        const scalar& curAlpha0f_nei,
        const scalar& MIN_ALPHA_DIFF,
        label& curPhaseState
    );

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
    );

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
    );

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
