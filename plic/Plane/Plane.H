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

Class
    Foam::Plane

Description
    Geometric class that creates a 2D Plane and can return the intersection
    point between a line and the Plane.

SourceFiles
    Plane.C

\*---------------------------------------------------------------------------*/

#ifndef Plane_H
#define Plane_H

#include "point.H"
#include "scalarList.H"
#include "dictionary.H"
#include "line.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class Plane;
bool operator==(const Plane&, const Plane&);
bool operator!=(const Plane&, const Plane&);
Istream& operator>>(Istream&, Plane&);
Ostream& operator<<(Ostream&, const Plane&);


/*---------------------------------------------------------------------------*\
                           Class Plane Declaration
\*---------------------------------------------------------------------------*/

class Plane
{
public:

    //- Side of the Plane
    enum side
    {
        NORMAL,
        FLIP
    };


        //- A direction and a reference point
        class ray
        {
            point refPoint_;

            vector dir_;

        public:

            ray(const point& refPoint, const vector& dir)
            :
                refPoint_(refPoint),
                dir_(dir)
            {}

            const point& refPoint() const
            {
                return refPoint_;
            }

            const vector& dir() const
            {
                return dir_;
            }
        };


private:

    // Private data

        //- Plane normal
        vector unitVector_;

        //- Base point
        point basePoint_;


    // Private Member Functions

        //- Calculates basePoint and normal vector given Plane coefficients
        void calcPntAndVec(const scalarList& C);

        //- Calculates basePoint and normal vector given three points
        //- Normal vector determined using right hand rule
        void calcPntAndVec
        (
            const point& point1,
            const point& point2,
            const point& point3
        );


public:

    // Constructors

        //- Construct null
        Plane();

        //- Construct from normal vector through the origin
        Plane(const vector& normalVector);

        //- Construct from normal vector and point in Plane
        Plane(const point& basePoint, const vector& normalVector);

        //- Construct from three points in Plane
        Plane(const point& point1, const point& point2, const point& point3);

        //- Construct from coefficients for the
        //  Plane equation: ax + by + cz + d = 0
        Plane(const scalarList& C);

        //- Construct from dictionary
        Plane(const dictionary& PlaneDict);

        //- Construct from Istream. Assumes the base + normal notation.
        Plane(Istream& is);


    // Member Functions

        //- Return Plane normal
        const vector& normal() const;

        //- Return or return Plane base point
        const point& refPoint() const;

        //- Return coefficients for the
        //  Plane equation: ax + by + cz + d = 0
        FixedList<scalar, 4> PlaneCoeffs() const;

        //- Return nearest point in the Plane for the given point
        point nearestPoint(const point& p) const;

        //- Return distance from the given point to the Plane
        scalar distance(const point& p) const;

        //- Return cut coefficient for Plane and line defined by
        //  origin and direction
        scalar normalIntersect(const point& pnt0, const vector& dir) const;

        //- Return cut coefficient for Plane and ray
        scalar normalIntersect(const ray& r) const
        {
            return normalIntersect(r.refPoint(), r.dir());
        }

        //- Return the cutting point between the Plane and
        // a line passing through the supplied points
        template<class Point, class PointRef>
        scalar lineIntersect(const line<Point, PointRef>& l) const
        {
            return normalIntersect(l.start(), l.vec());
        }

        //- Return the cutting line between this Plane and another.
        //  Returned as direction vector and point line goes through.
        ray PlaneIntersect(const Plane&) const;

        //- Return the cutting point between this Plane and two other Planes
        point PlanePlaneIntersect(const Plane&, const Plane&) const;

        //- Return the side of the Plane that the point is on.
        //  If the point is on the Plane, then returns NORMAL.
        side sideOfPlane(const point& p) const;

        //- Mirror the supplied point in the Plane. Return the mirrored point.
        point mirror(const point& p) const;

        //- Write to dictionary
        void writeDict(Ostream&) const;


    // friend Operators
 
        friend bool operator==(const Plane&, const Plane&);
        friend bool operator!=(const Plane&, const Plane&);


    // IOstream Operators
        
        friend Istream& operator>>(Istream&, Plane&);
        //- Write Plane properties
        friend Ostream& operator<<(Ostream&, const Plane&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
