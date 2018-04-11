/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Class    
    Foam::plic

Description
    Piecewise-linear interface construction class definition

SourceFiles
    plicI.H
    plic.C

\*---------------------------------------------------------------------------*/

#ifndef plic_H
#define plic_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "scalarField.H"
#include "vectorField.H"
#include "faceList.H"
#include "labelList.H"
#include "boolList.H"
#include "cell.H"
#include "face.H"
#include "edge.H"
#include "point.H"
#include "Plane.H"
#include "centredFPCCellToFaceStencil.H"
#include "tetrahedron.H"
#include "tetPoints.H"
#include "cellInfo.H"
#include "volPointInterpolation.H"
#include "IOdictionary.H"
#include "plicLeastSquaresGrad.H"
#include "plicFuncs.H"
#include "wordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class plic Declaration
\*---------------------------------------------------------------------------*/

class plic
{
    // Private data

    //- Reference to mesh
    const fvMesh& mesh_;

    //- extended face stencil
    centredFPCCellToFaceStencil faceStencil_;

    //- least squares grad
    plicLeastSquaresGrad lsGrad_;

    //- Reference to volume fraction field
    const volScalarField& alpha_ph1_;

    //- Input dictionary
    IOdictionary transPropDict_;

    //- Gradient of volume fraction field
    //volVectorField& gradAlpha_;

    //- Reference to velocity field
    const volVectorField& U_;

    //- Cell volumes to cell points interpolation object
    volPointInterpolation ptInterp_;
    
    //- velocity field interpolated to mesh points
    pointVectorField U_pts_;

    //- volume fraction field interpolated to mesh points
    pointScalarField alpha_ph1_pts_;

    //- volume fraction field on faces
    surfaceScalarField alpha1f_;

    //- Reference to density field
    const volScalarField& rho1_;
    const volScalarField& rho0_;
    const volScalarField& rho_;

    //- Reference to flux field
    const surfaceScalarField& phi_;

    //- alpha1 flux at faces
    surfaceScalarField phiAlpha1_;

    surfaceScalarField phiAlpha0_;
  
    //- Number of mixed phase cells
    label nMixCells_;
    //- Number of mixed phase faces
    label nMixFcs_;

    //- Label to indicate phase state of cells
    labelList cell_phaseState_;

    //- Label to indicate phase state of faces
    //- own side intfc reconstruction
    labelList face_phaseState_own_;
    //- nei side intfc reconstruction
    labelList face_phaseState_nei_;

    //- Bool list to indicate if cell is a
    //- near interface cell
    boolList cell_near_intfc_;

    //- Bool list to indicate if face needs 
    //- 2-phase flux calculation
    boolList fc_2ph_flux_needed_;

    //- Gradient of volume fraction field
    volVectorField gradAlpha1_;

    //- Interface unit normal vector field
    volVectorField nHat_;        
    surfaceScalarField nHatf_;

    //- Interface face center point field
    volVectorField C_intfc_;
    //- Interface area field
    volScalarField A_intfc_;    

    //- Phase centroid fields
    volVectorField C_ph1_;
    volVectorField C_ph0_;

    //- Phase centroid fields on faces
    //- own side intfc reconstruction
    vectorField Cf_ph1_own_;
    vectorField Cf_ph0_own_;
    //- nei side intfc reconstruction
    vectorField Cf_ph1_nei_;
    vectorField Cf_ph0_nei_;

    //- Phase areas fields on faces
    //- own side intfc reconstruction
    scalarField Af_ph1_own_;
    scalarField Af_ph0_own_;
    //- nei side intfc reconstruction
    scalarField Af_ph1_nei_;
    scalarField Af_ph0_nei_;

    //- Cell center field including boundary cells
    volVectorField meshC_;

    wordList species_;
    label nSpecies_;

    surfaceScalarField fAlpha1_;
    surfaceScalarField frho1_;
    surfaceScalarField frho0_;
    PtrList<surfaceScalarField> fY1_;
    PtrList<surfaceScalarField> fY0_;
    PtrList<surfaceScalarField> fc1_;
    PtrList<surfaceScalarField> fc0_;
    surfaceScalarField fh1_;
    surfaceScalarField fh0_;

    //- Flat fields for extended stencil use
    //- Use mapDistribute to obtain remote values
    List<cellInfo> cells_flatFld_;
    List<List<Plane> > cellPlns_flatFld_;
    List<scalar> alpha_ph1_flatFld_;    
    List<vector> nHat_flatFld_;
    List<vector> gradAlpha1_flatFld_;
    List<vector> C_intfc_flatFld_;
    
    List<scalar> rho1_flatFld_;
    List<scalar> rho0_flatFld_;
    List<List<scalar> > Y1_flatFld_;
    List<List<scalar> > Y0_flatFld_;
    List<List<scalar> > c1_flatFld_;
    List<List<scalar> > c0_flatFld_;
    List<scalar> h1_flatFld_;
    List<scalar> h0_flatFld_;

    //-mesh faces info
    vectorField Cf_;
    vectorField Sf_;
    scalarField magSf_;

    scalar ALPHA_2PH_MIN;
    scalar ALPHA_2PH_MAX;
    scalar GRADALPHA_MIN;
    scalar ALPHA_TOL;
    scalar DIST_TOL;
    scalar PT_DIST_TOL;
    scalar SMALLEST_TET_VOL;
    scalar SMALLEST_PHI;
    label MAX_BRENT_ITERS;
    
    bool debug_;
    bool debug2_;
    bool debugIR_;
    bool debugIR2_;
    bool debugF_;
    bool debugF2_;

    label brent_iters_max_;
    scalar brent_err_max_;
    label brent_max_cell_;
    vector brent_max_cellC_;
    scalar brent_max_cellAlpha_;
    vector brent_max_cell_nHat_;
    vector brent_max_cell_gradAlpha1_;
    label brent_iters_tmp_;
    scalar brent_err_tmp_;
    
    // Member Functions

    void intfc_cell_reconstruct
    (
        const vector&,
        const scalar&,
        const cell&,
        const faceList&,
        const pointField&,
        cell&,
        cell&,
        faceList&,
        pointField&,
        face&,
        labelList&,
        labelList&
    );
  
    void plane_cell_intersect
    (
        const Plane&, 
        const cell&,
        const faceList&,
        const pointField&,
        cell&,
        cell&,
        faceList&,
        pointField&,
        face&,
        labelList&,
        labelList&
    ); 

    scalar calc_signed_dist(const Plane&, const point&);

    point plane_edge_intersect
    (
        const point&, 
        const point&,
        const scalar&,
        const scalar&
    );

    label append_point_if_new(pointField&, const point&);

    void remove_elems
    (
        labelList&,
        const label,
        const label
    );

    void remove_duplicateLabels(labelList&);

    void subCells_reset
    (
        cell&, 
        cell&, 
        faceList&, 
        pointField&, 
        face&
    );

    void swap_scalar(scalar&, scalar&);

    void calcFaceFluxTets(const label, const scalar&, List<tetPoints>&, List<scalar>&);

    void calcFaceFluxTets
    (
        const label curFcI,
        const scalar& curPhi,
        List<tetPoints>& tets,
        List<scalar>& tetVols,
        const scalar& deltaT
    );

    face face_from_lbls(const label, const label, const label, const label);

    void makePatchFaceCellInfo(const fvPatch&, const label, cellInfo&);

    void cellInfoAndPlanes_collectData();

    void intfcInfo_collectData();

    void Y_collectData
    (
        const PtrList<volScalarField>&,
        const PtrList<volScalarField>&
    );

    void c_h_collectData
    (
        const PtrList<volScalarField>& c1,
        const PtrList<volScalarField>& c0,
        const volScalarField& h1,
        const volScalarField& h0
    );

    scalar tet_cell_intersect(tetPoints&, const label&);

    void tet_cell_intersect(tetPoints&, const label&, scalar&, scalar&);
    
    void tet_slice_with_plns(tetPoints&, const List<Plane>&, scalar&, const label&); 

    void tet_slice_with_plns(tetPoints&, const List<Plane>&, scalar&, scalar&, const label&); 

  public:

    //- Runtime type information
    TypeName("plic");

    // Constructors

    //- Construct from mesh, volume fraction field, velocity field,
    //- volumetric flux field and density field
    plic
    (
        const fvMesh&,
        const volScalarField&,
        const volVectorField&,
        const surfaceScalarField&,
        const volScalarField&,
        const volScalarField&,
        const volScalarField&     
    );


    // Destructors
    virtual ~plic();

    // Public member functions

    void intfc_normal_correct_lsq();
    
    void intfc_normal_correct();
    
    void intfc_correct();

    void calc_face_phaseFluxes();

    void calc_2ph_advFluxes
    (
        const PtrList<volScalarField>&,
        const PtrList<volScalarField>&,
        PtrList<surfaceScalarField>&,
        PtrList<surfaceScalarField>&
    );

    void calc_2ph_advFluxes
    (
        const PtrList<volScalarField>& c1,
        const PtrList<volScalarField>& c0,    
        const volScalarField& h1,
        const volScalarField& h0,        
        const scalar& deltaT,
        surfaceScalarField& advFlux_rho1,
        surfaceScalarField& advFlux_rho0,
        PtrList<surfaceScalarField>& advFlux_c1,
        PtrList<surfaceScalarField>& advFlux_c0,
        surfaceScalarField& advFlux_h1,
        surfaceScalarField& advFlux_h0,
        const bool debug,
        const bool debug2,
        OFstream& os
    );

    void plic_correct();

    inline const fvMesh& mesh() const;

    inline const centredFPCCellToFaceStencil& faceStencil() const;    

    inline surfaceScalarField& phiAlpha1();

    inline const surfaceScalarField& phiAlpha1() const;

    inline const label& nMixCells() const;

    inline const label& nMixFaces() const;
    
    inline const labelList& cell_phaseState() const;

    inline const labelList& face_phaseState_own() const;

    inline const labelList& face_phaseState_nei() const;

    inline const boolList& cell_near_intfc() const;

    inline const boolList& face_2ph_flux_needed() const;

    inline const volVectorField& gradAlpha1() const;

    inline const volVectorField& nHat() const;    

    const surfaceScalarField& nHatf() const
    { 
        return nHatf_;
    }

    inline const volVectorField& C_intfc() const;

    inline const volScalarField& A_intfc() const;    

    inline const volVectorField& C_ph1() const;

    inline const volVectorField& C_ph0() const;

    inline const vectorField& Cf_ph1_own() const;

    inline const vectorField& Cf_ph0_own() const;

    inline const scalarField& Af_ph1_own() const;

    inline const scalarField& Af_ph0_own() const;

    inline const vectorField& Cf_ph1_nei() const;

    inline const vectorField& Cf_ph0_nei() const;

    inline const scalarField& Af_ph1_nei() const;

    inline const scalarField& Af_ph0_nei() const;

    inline const List<scalar>& alpha_ph1_flatFld() const;

    inline const List<vector>& nHat_flatFld() const;

    inline const List<vector>& C_intfc_flatFld() const;
    
    bool isEq_pt(const point&, const point&);

    const vectorField& meshCf() const
    { 
        return Cf_;
    }

    const vectorField& meshSf() const
    { 
        return Sf_;
    }

    const scalarField& meshMagSf() const
    { 
        return magSf_;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "plicI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
