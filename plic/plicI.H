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

Description


    Calculate the centre by breaking the face into triangles using the face
    centre and area-weighted averaging their centres.  This method copes with
    small face-concavity.

\*---------------------------------------------------------------------------*/

#include "plic.H"


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * //

inline const Foam::fvMesh& Foam::plic::mesh() const
{
    return mesh_;
}


inline const Foam::centredFPCCellToFaceStencil& Foam::plic::faceStencil() const
{
    return faceStencil_;
}


inline Foam::surfaceScalarField& Foam::plic::phiAlpha1()
{
    return phiAlpha1_;
}


inline const Foam::surfaceScalarField& Foam::plic::phiAlpha1() const
{
    return phiAlpha1_;
}


inline const Foam::label& Foam::plic::nMixCells() const
{
    return nMixCells_;
}


inline const Foam::label& Foam::plic::nMixFaces() const
{
    return nMixFcs_;
}


inline const Foam::labelList& Foam::plic::cell_phaseState() const
{
    return cell_phaseState_;
}


inline const Foam::labelList& Foam::plic::face_phaseState() const
{
    return face_phaseState_;
}


inline const Foam::boolList& Foam::plic::cell_near_intfc() const
{
    return cell_near_intfc_;
}


inline const Foam::boolList& Foam::plic::face_2ph_flux_needed() const
{
    return fc_2ph_flux_needed_;
}


inline const Foam::volVectorField& Foam::plic::gradAlpha1() const
{
    return gradAlpha1_;
}


inline const Foam::volVectorField& Foam::plic::nHat() const
{
    return nHat_;
}


inline const Foam::volVectorField& Foam::plic::C_intfc() const
{
    return C_intfc_;
}


inline const Foam::volScalarField& Foam::plic::A_intfc() const
{
    return A_intfc_;
}


inline const Foam::volVectorField& Foam::plic::C_ph1() const
{
    return C_ph1_;
}


inline const Foam::volVectorField& Foam::plic::C_ph0() const
{
    return C_ph0_;
}


inline const Foam::vectorField& Foam::plic::Cf_ph1() const
{
    return Cf_ph1_;
}


inline const Foam::vectorField& Foam::plic::Cf_ph0() const
{
    return Cf_ph0_;
}


inline const Foam::scalarField& Foam::plic::Af_ph1() const
{
    return Af_ph1_;
}


inline const Foam::scalarField& Foam::plic::Af_ph0() const
{
    return Af_ph0_;
}


inline const Foam::List<Foam::scalar>& Foam::plic::alpha_ph1_flatFld() const
{
    return alpha_ph1_flatFld_;
}


inline const Foam::List<Foam::vector>& Foam::plic::nHat_flatFld() const
{
    return nHat_flatFld_;
}


inline const Foam::List<Foam::vector>& Foam::plic::C_intfc_flatFld() const
{
    return C_intfc_flatFld_;
}


// ************************************************************************* //
