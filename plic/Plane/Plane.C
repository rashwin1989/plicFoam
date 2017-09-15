/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "Plane.H"
#include "tensor.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Calculate base point and unit normal vector from Plane equation
void Foam::Plane::calcPntAndVec(const scalarList& C)
{
    if (mag(C[0]) > VSMALL)
    {
        basePoint_ = vector((-C[3]/C[0]), 0, 0);
    }
    else
    {
        if (mag(C[1]) > VSMALL)
        {
            basePoint_ = vector(0, (-C[3]/C[1]), 0);
        }
        else
        {
            if (mag(C[2]) > VSMALL)
            {
                basePoint_ = vector(0, 0, (-C[3]/C[2]));
            }
            else
            {
                FatalErrorIn("void Plane::calcPntAndVec(const scalarList&)")
                    << "At least one Plane coefficient must have a value"
                    << abort(FatalError);
            }
        }
    }

    unitVector_ = vector(C[0], C[1], C[2]);
    scalar magUnitVector(mag(unitVector_));

    if (magUnitVector < VSMALL)
    {
        FatalErrorIn("void Plane::calcPntAndVec(const scalarList&)")
            << "Plane normal defined with zero length"
            << abort(FatalError);
    }

    unitVector_ /= magUnitVector;
}


void Foam::Plane::calcPntAndVec
(
    const point& point1,
    const point& point2,
    const point& point3
)
{
    basePoint_ = (point1 + point2 + point3)/3;
    vector line12 = point1 - point2;
    vector line23 = point2 - point3;

    if
    (
        mag(line12) < VSMALL
     || mag(line23) < VSMALL
     || mag(point3-point1) < VSMALL
    )
    {
        FatalErrorIn
        (
            "void Plane::calcPntAndVec\n"
            "(\n"
            "    const point&,\n"
            "    const point&,\n"
            "    const point&\n"
            ")\n"
        )   << "Bad points:" << point1 << ' ' << point2 << ' ' << point3
            << abort(FatalError);
    }

    unitVector_ = line12 ^ line23;
    scalar magUnitVector(mag(unitVector_));

    if (magUnitVector < VSMALL)
    {
        FatalErrorIn
        (
            "void Plane::calcPntAndVec\n"
            "(\n"
            "    const point&,\n"
            "    const point&,\n"
            "    const point&\n"
            ")\n"
        )   << "Plane normal defined with zero length" << nl
            << "Bad points:" << point1 << ' ' << point2 << ' ' << point3
            << abort(FatalError);
    }

    unitVector_ /= magUnitVector;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct null
Foam::Plane::Plane()
{}

// Construct from normal vector through the origin
Foam::Plane::Plane(const vector& normalVector)
:
    unitVector_(normalVector),
    basePoint_(vector::zero)
{
    scalar magUnitVector(mag(unitVector_));

    if (magUnitVector > VSMALL)
    {
        unitVector_ /= magUnitVector;
    }
    else
    {
        FatalErrorIn("Plane::Plane(const vector&)")
            << "Plane normal has zero length. basePoint:" << basePoint_
            << abort(FatalError);
    }
}


// Construct from point and normal vector
Foam::Plane::Plane(const point& basePoint, const vector& normalVector)
:
    unitVector_(normalVector),
    basePoint_(basePoint)
{
    scalar magUnitVector(mag(unitVector_));

    if (magUnitVector > VSMALL)
    {
        unitVector_ /= magUnitVector;
    }
    else
    {
        FatalErrorIn("Plane::Plane(const point&, const vector&)")
            << "Plane normal has zero length. basePoint:" << basePoint_
            << abort(FatalError);
    }
}


// Construct from Plane equation
Foam::Plane::Plane(const scalarList& C)
{
    calcPntAndVec(C);
}


// Construct from three points
Foam::Plane::Plane
(
    const point& a,
    const point& b,
    const point& c
)
{
    calcPntAndVec(a, b, c);
}


// Construct from dictionary
Foam::Plane::Plane(const dictionary& dict)
:
    unitVector_(vector::zero),
    basePoint_(point::zero)
{
    const word PlaneType(dict.lookup("PlaneType"));

    if (PlaneType == "PlaneEquation")
    {
        const dictionary& subDict = dict.subDict("PlaneEquationDict");
        scalarList C(4);

        C[0] = readScalar(subDict.lookup("a"));
        C[1] = readScalar(subDict.lookup("b"));
        C[2] = readScalar(subDict.lookup("c"));
        C[3] = readScalar(subDict.lookup("d"));

        calcPntAndVec(C);

    }
    else if (PlaneType == "embeddedPoints")
    {
        const dictionary& subDict = dict.subDict("embeddedPointsDict");

        point point1(subDict.lookup("point1"));
        point point2(subDict.lookup("point2"));
        point point3(subDict.lookup("point3"));

        calcPntAndVec(point1, point2, point3);
    }
    else if (PlaneType == "pointAndNormal")
    {
        const dictionary& subDict = dict.subDict("pointAndNormalDict");

        basePoint_ = subDict.lookup("basePoint");
        unitVector_ = subDict.lookup("normalVector");
        unitVector_ /= mag(unitVector_);
    }
    else
    {
        FatalIOErrorIn("Plane::Plane(const dictionary&)", dict)
            << "Invalid Plane type: " << PlaneType << nl
            << "Valid options include: PlaneEquation, embeddedPoints and "
            << "pointAndNormal"
            << abort(FatalIOError);
    }
}


// Construct from Istream. Assumes point and normal vector.
Foam::Plane::Plane(Istream& is)
:
    unitVector_(is),
    basePoint_(is)
{
    scalar magUnitVector(mag(unitVector_));

    if (magUnitVector > VSMALL)
    {
        unitVector_ /= magUnitVector;
    }
    else
    {
        FatalErrorIn("Plane::Plane(Istream& is)")
            << "Plane normal has zero length. basePoint:" << basePoint_
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return Plane normal vector
const Foam::vector& Foam::Plane::normal() const
{
    return unitVector_;
}


// Return Plane base point
const Foam::point& Foam::Plane::refPoint() const
{
    return basePoint_;
}


// Return coefficients for Plane equation: ax + by + cz + d = 0
Foam::FixedList<Foam::scalar, 4> Foam::Plane::PlaneCoeffs() const
{
    FixedList<scalar, 4> C(4);

    scalar magX = mag(unitVector_.x());
    scalar magY = mag(unitVector_.y());
    scalar magZ = mag(unitVector_.z());

    if (magX > magY)
    {
        if (magX > magZ)
        {
            C[0] = 1;
            C[1] = unitVector_.y()/unitVector_.x();
            C[2] = unitVector_.z()/unitVector_.x();
        }
        else
        {
            C[0] = unitVector_.x()/unitVector_.z();
            C[1] = unitVector_.y()/unitVector_.z();
            C[2] = 1;
        }
    }
    else
    {
        if (magY > magZ)
        {
            C[0] = unitVector_.x()/unitVector_.y();
            C[1] = 1;
            C[2] = unitVector_.z()/unitVector_.y();
        }
        else
        {
            C[0] = unitVector_.x()/unitVector_.z();
            C[1] = unitVector_.y()/unitVector_.z();
            C[2] = 1;
        }
    }

    C[3] = - C[0] * basePoint_.x()
           - C[1] * basePoint_.y()
           - C[2] * basePoint_.z();

    return C;
}


// Return nearest point in the Plane for the given point
Foam::point Foam::Plane::nearestPoint(const point& p) const
{
    return p - unitVector_*((p - basePoint_) & unitVector_);
}


// Return distance from the given point to the Plane
Foam::scalar Foam::Plane::distance(const point& p) const
{
    return mag((p - basePoint_) & unitVector_);
}


// Cutting point for Plane and line defined by origin and direction
Foam::scalar Foam::Plane::normalIntersect
(
    const point& pnt0,
    const vector& dir
) const
{
    scalar denom = stabilise((dir & unitVector_), VSMALL);

    return ((basePoint_ - pnt0) & unitVector_)/denom;
}


// Cutting line of two Planes
Foam::Plane::ray Foam::Plane::PlaneIntersect(const Plane& Plane2) const
{
    // Mathworld Plane-Plane intersection. Assume there is a point on the
    // intersection line with z=0 and solve the two Plane equations
    // for that (now 2x2 equation in x and y)
    // Better: use either z=0 or x=0 or y=0.

    const vector& n1 = normal();
    const vector& n2 = Plane2.normal();

    const point& p1 = refPoint();
    const point& p2 = Plane2.refPoint();

    scalar n1p1 = n1&p1;
    scalar n2p2 = n2&p2;

    vector dir = n1 ^ n2;

    // Determine zeroed out direction (can be x,y or z) by looking at which
    // has the largest component in dir.
    scalar magX = mag(dir.x());
    scalar magY = mag(dir.y());
    scalar magZ = mag(dir.z());

    direction iZero, i1, i2;

    if (magX > magY)
    {
        if (magX > magZ)
        {
            iZero = 0;
            i1 = 1;
            i2 = 2;
        }
        else
        {
            iZero = 2;
            i1 = 0;
            i2 = 1;
        }
    }
    else
    {
        if (magY > magZ)
        {
            iZero = 1;
            i1 = 2;
            i2 = 0;
        }
        else
        {
            iZero = 2;
            i1 = 0;
            i2 = 1;
        }
    }

    vector pt;

    pt[iZero] = 0;
    pt[i1] = (n2[i2]*n1p1 - n1[i2]*n2p2) / (n1[i1]*n2[i2] - n2[i1]*n1[i2]);
    pt[i2] = (n2[i1]*n1p1 - n1[i1]*n2p2) / (n1[i2]*n2[i1] - n1[i1]*n2[i2]);

    return ray(pt, dir);
}


// Cutting point of three Planes
Foam::point Foam::Plane::PlanePlaneIntersect
(
    const Plane& Plane2,
    const Plane& Plane3
) const
{
    FixedList<scalar, 4> coeffs1(PlaneCoeffs());
    FixedList<scalar, 4> coeffs2(Plane2.PlaneCoeffs());
    FixedList<scalar, 4> coeffs3(Plane3.PlaneCoeffs());

    tensor a
    (
        coeffs1[0],coeffs1[1],coeffs1[2],
        coeffs2[0],coeffs2[1],coeffs2[2],
        coeffs3[0],coeffs3[1],coeffs3[2]
    );

    vector b(coeffs1[3],coeffs2[3],coeffs3[3]);

    return (inv(a) & (-b));
}


Foam::Plane::side Foam::Plane::sideOfPlane(const point& p) const
{
    const scalar angle((p - basePoint_) & unitVector_);

    return (angle < 0 ? FLIP : NORMAL);
}


Foam::point Foam::Plane::mirror(const point& p) const
{
    const vector mirroredPtDir = p - nearestPoint(p);

    if ((normal() & mirroredPtDir) > 0)
    {
        return p - 2.0*distance(p)*normal();
    }
    else
    {
        return p + 2.0*distance(p)*normal();
    }
}


void Foam::Plane::writeDict(Ostream& os) const
{
    os.writeKeyword("PlaneType") << "pointAndNormal"
        << token::END_STATEMENT << nl;
    os  << indent << "pointAndNormalDict" << nl
        << indent << token::BEGIN_BLOCK << incrIndent << nl;
    os.writeKeyword("basePoint") << basePoint_ << token::END_STATEMENT << nl;
    os.writeKeyword("normalVector") << unitVector_ << token::END_STATEMENT
        << nl;
    os << decrIndent << indent << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

bool Foam::operator==(const Plane& a, const Plane& b)
{
    if (a.basePoint_ == b.basePoint_ && a.unitVector_ == b.unitVector_)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool Foam::operator!=(const Plane& a, const Plane& b)
{
    return !(a == b);
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, Plane& a)
{
    is.readBegin("Plane");
    is  >> a.unitVector_ >> a.basePoint_;
    is.readEnd("Plane");

    is.check("Istream& operator>>(Istream&, Plane&)");
    
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const Plane& a)
{
    os  << nl
        << token::BEGIN_LIST
        << a.unitVector_ << token::SPACE 
        << a.basePoint_
        << token::END_LIST;

    return os;
}


// ************************************************************************* //
