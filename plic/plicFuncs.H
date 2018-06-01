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
#include "IOmanip.H"

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include "MACROS.H"
#include "PR_EoS.h"
#include "myUmfpack.h"
#include "vis_n_therm.h"
#include "transport_LLE2.h"
#include "gsl_optimizations2.h"
#include "Maxwell_Stefan_flux.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace plicFuncs functions Declaration
\*---------------------------------------------------------------------------*/

namespace plicFuncs
{

#define    _NNEW_(var,type,size)    var = new type[size];
#define    _DDELETE_(var)           delete[] var;
    
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
    void print_field
    (
        const GeometricField<Type, fvPatchField, volMesh>&,
        OFstream& os
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
    void print_surfaceField
    (
        const GeometricField<Type, fvsPatchField, surfaceMesh>&,
        const fvMesh&,
        OFstream& os
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
        bool&,
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
        bool&,
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
        const scalar& Y1MIN,
        const scalar& Y1MAX,
        const scalar& Y0MIN,
        const scalar& Y0MAX,
        const label& i,
        const scalar& dt,
        const bool debug, 
        OFstream& os
    );

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
    );

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
    );

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
    );

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
    );

    void calc_2ph_diffFluxes_T
    (
        const fvMesh& mesh,
        double R_gas,
        int n,
        double *Pc,
        double *Tc,
        double *w,
        double *MW,
        int *tk,
        double *Tb,
        double *SG,
        double *H8,
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
    );

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
        double& T,
        const int& MAX_ITER_T,
        const double& T_TOL,
        const double& H_TOL
    );

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
        double& h
    );

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
        const scalar& dt,
        const scalar& diffFlux,
        scalar& diffFlux_limiter,
        const int& MAX_ITER_T,
        const double& T_TOL,
        const double& H_TOL,
        double TMIN,
        double TMAX
    );

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
        const scalar& dt,
        const scalar& diffFlux,
        scalar& diffFlux_limiter,
        double TMIN,
        double TMAX
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
        scalar& diffFlux_limiter,
        const scalar& YMIN,
        const scalar& YMAX
    );

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
        const label& Y_BOUND_ITERS_MAX,
        const bool debug,
        OFstream& os
    );

    label findCellInIntfcDir
    (
        const fvMesh& mesh,
        const List<scalar>& alpha1,
        const labelList& cells,
        const List<vector>& C,
        const vector& Cp,
        const vector& nf,
        const label& C1_lbl,
        bool& foundCell,
        bool debug,
        OFstream& os
    );

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
    );

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
    );

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
    );

    void x2y(int n, double *MW, double *x, double *y);
    
    void y2x(int n, double *MW, double *y, double *x);

    double calc_rho_from_V(int n, double *x, double *MW, double V);

    double calc_CvIG_from_CpIG(int n, double *x, double *CpIG);

    void correct_rho
    (
        int n,
        double *Pc,
        double *Tc,
        double *w,
        double *MW,
        int *tk,
        double *coef_ab,
        double P,
        const PtrList<volScalarField>& X,
        const volScalarField& T,
        volScalarField& rho
    );

    void correct_h
    (
        int n,
        double *Pc,
        double *Tc,
        double *w,
        double *MW,
        double *Tb,
        double *SG,
        double *H8,
        double P,
        const PtrList<volScalarField>& X,
        const volScalarField& T,
        volScalarField& h
    );

    void correct_boundaryField_C
    (    
        volScalarField& Yi,
        const volScalarField& rho,
        const volScalarField& alpha,
        volScalarField& ci,
        volScalarField& Ci
    );

    void correct_x_from_Y
    (
        int n,
        double *MW,
        const PtrList<volScalarField>& Y,
        PtrList<volScalarField>& x
    );

    void correct_boundaryField_h_rhoh_H
    (
        int n,
        double *Pc,
        double *Tc,
        double *w,
        double *MW,
        double *Tb,
        double *SG,
        double *H8,
        double P,
        const PtrList<volScalarField>& X,
        const volScalarField& T,
        const volScalarField& rho,
        const volScalarField& alpha,
        volScalarField& h,
        volScalarField& rhoh,
        volScalarField& H
    );

    void correct_thermo_trans_prop
    (
        const fvMesh& mesh,    
        int n,
        double *Pc,
        double *Tc,
        double *Vc,
        double *w,
        double *MW,
        int *tk,
        double *coef_ab,
        double *Tb,
        double *SG,
        double *H8,
        double *k,
        double *dm,
        double P,
        const PtrList<volScalarField>& X,
        const volScalarField& T,
        volScalarField& v,
        PtrList<volScalarField>& hpar,
        volScalarField& lambda,
        volScalarField& Cp,
        volScalarField& mu,
        PtrList<volScalarField>& D
    );

    void correct_thermo_trans_prop
    (
        const fvMesh& mesh,    
        int n,
        double *Pc,
        double *Tc,
        double *Vc,
        double *w,
        double *MW,
        int *tk,
        double *coef_ab,
        double *Tb,
        double *SG,
        double *H8,
        double *k,
        double *dm,
        double P,
        const PtrList<volScalarField>& X,
        const volScalarField& T,
        volScalarField& v,
        volScalarField& mu,
        PtrList<volScalarField>& D
    );

    void correct_D_from_Dij
    (
        const PtrList<volScalarField>& x,
        const PtrList<volScalarField>& Dij,
        PtrList<volScalarField>& D,
        int n
    );

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
    );

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
    );

    void print_line(OFstream& os, int n);

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
    );

    void calc_Xs_Ys_Js_mS_alphaS
    (
        const fvMesh& mesh,
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
        volScalarField& JsTot,
        volScalarField& mS1Tot,
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
    );

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
    );

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
    );

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
    );

    void print_advFluxFld
    (
        const fvMesh& mesh,
        const surfaceScalarField& advFlux,
        const volScalarField& alpha,
        const word& alphaName,
        OFstream& os
    );

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
    );

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
    );

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
    );

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
