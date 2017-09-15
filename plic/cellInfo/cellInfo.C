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

\*---------------------------------------------------------------------------*/

#include "cellInfo.H"
#include "dictionaryEntry.H"
#include "pyramidPointFaceRef.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(cellInfo, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::cellInfo::append_point_if_new
(
    pointField& pts,
    const point& pt
)
{    
    for(label i=0; i<pts.size(); i++)
    {
        if(pt == pts[i])
        { 
            return i;
        }
    }
    
    pts.append(pt);
    return (pts.size() - 1);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellInfo::cellInfo()
:
    nPoints_(8),
    nFaces_(6),
    points_(nPoints_),
    faces_(nFaces_)
{
    centre_ = vector::zero;
    vol_ = 0;
}


Foam::cellInfo::cellInfo
(
    const pointField& pts,
    const faceList& fcs,
    const labelList& own,
    const cell& curCell,
    const label& cellI
)
:
    nPoints_(0),
    nFaces_(curCell.size()),
    points_(nPoints_),
    faces_(nFaces_)
{    
    for(label faceI=0; faceI<curCell.size(); faceI++)
    {
        face curFc = fcs[curCell[faceI]];
        face curFace(curFc.size());

        for(label pointI=0; pointI<curFc.size(); pointI++)
        {
            point curPt = pts[curFc[pointI]];
            curFace[pointI] = append_point_if_new(points_,curPt);
        }
        
        if(own[faceI] == cellI)
        {
            curFace.flip();
        }

        faces_[faceI] = curFace;
    }

    calc_centreAndVol();
    correct_faceDir();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellInfo::~cellInfo()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::cellInfo::calc_centreAndVol()
{
    vector cEst = vector::zero;
    scalar sumArea = 0;

    forAll(faces_, faceI)
    {
        scalar a = faces_[faceI].mag(points_);
        cEst += faces_[faceI].centre(points_)*a;
        sumArea += a;
    }

    cEst /= sumArea + VSMALL;

    // Calculate the centre by breaking the cell into pyramids and
    // volume-weighted averaging their centres
    vector sumVc = vector::zero;

    scalar sumV = 0;

    forAll(faces_, faceI)
    {
        // calculate pyramid volume. If it is greater than zero, OK.
        // If not, the pyramid is inside-out. Create a face with the opposite
        // order and recalculate pyramid centre!
        scalar pyrVol = pyramidPointFaceRef(faces_[faceI], cEst).mag(points_);
        vector pyrCentre = pyramidPointFaceRef(faces_[faceI], cEst).centre(points_);

        // if pyramid inside-out because face points inwards invert
        // N.B. pyramid remains unchanged
        if (pyrVol < 0)
        {
            pyrVol = -pyrVol;
        }

        sumVc += pyrVol*pyrCentre;
        sumV += pyrVol;
    }

    centre_ =  sumVc/(sumV + VSMALL);
    vol_ = sumV;
}


void Foam::cellInfo::correct_faceDir()
{
    forAll(faces_,faceI)
    {
        face& curFace = faces_[faceI];
        const point& curFaceCtr = curFace.centre(points_);
        const vector& curFaceNormal = curFace.normal(points_);
        scalar cDotN = (centre_ - curFaceCtr) & curFaceNormal;

        if(cDotN < 0)
        {
            curFace.flip();
        }
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

bool Foam::operator==(const cellInfo& a, const cellInfo& b)
{
    // Trivial reject: cells are different size
    if (a.nPoints_ != b.nPoints_)
    {
        return false;
    }

    if (a.nFaces_ != b.nFaces_)
    {
        return false;
    }    
    
    if (a.points_ != b.points_)
    {
        return false;
    }

    if (a.faces_ != b.faces_)
    {
        return false;
    }

    return true;
}


bool Foam::operator!=(const cellInfo& a, const cellInfo& b)
{
    return !(a == b);
}


Foam::Istream& Foam::operator>>(Istream& is, cellInfo& a)
{
    dictionaryEntry entry(dictionary::null, is);
    entry.lookup("nPoints") >> a.nPoints_;
    entry.lookup("nFaces") >> a.nFaces_;
    entry.lookup("points") >> a.points_;
    entry.lookup("faces") >> a.faces_;

    return is;
}

Foam::Ostream& Foam::operator<<(Ostream& os, const cellInfo& a)
{
    os  << "nPoints" << tab << a.nPoints_ << tab
        << "nFaces" << tab << a.nFaces_ << tab
        << "points" << tab << a.points_ << tab
        << "faces" << tab << a.faces_ << tab;

    return os;
}


// ************************************************************************* //
