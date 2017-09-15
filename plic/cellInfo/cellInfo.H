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
    Foam::cellInfo

Description
    Cell object with all faces and points of cell info stored
    within

SourceFiles
    cellInfoI.H
    cellInfo.C

\*---------------------------------------------------------------------------*/

#ifndef cellInfo_H
#define cellInfo_H

#include "pointField.H"
#include "faceList.H"
#include "labelList.H"
#include "cell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class cellInfo;
bool operator==(const cellInfo&, const cellInfo&);
bool operator!=(const cellInfo&, const cellInfo&);
Ostream& operator<<(Ostream&, const cellInfo&);
Istream& operator>>(Istream&, cellInfo&);


/*---------------------------------------------------------------------------*\
                             Class cellInfo Declaration
\*---------------------------------------------------------------------------*/

class cellInfo
{
    // Private data
    
    //- no. of points
    label nPoints_;

    //- no. of faces
    label nFaces_;
    
    //- cell points
    pointField points_;

    //- cell faces
    faceList faces_;

    //- cell centre
    point centre_;

    //- cell volume
    scalar vol_;

    // Private member functions
    label append_point_if_new(pointField&, const point&);

public:

    //- Runtime type information
    TypeName("cellInfo");

    // Constructors
    
    //- Construct hexahedral 
    cellInfo();
    
    //- Construct from cell, master pointField and master faceList
    cellInfo
    (
        const pointField&,
        const faceList&,
        const labelList&,
        const cell&,
        const label&
    );
 
    //- Destructor
    virtual ~cellInfo();


    // Member functions

    //- Calculate cell-centre and volume
    void calc_centreAndVol();

    //- Correct face orientations
    void correct_faceDir();

    //- Access to cell centre
    inline const point& centre() const;

    //- Access to volume
    inline const scalar& mag() const;

    //- Access to the cell points
    inline pointField& points();

    //- Access to the cell points
    inline const pointField& points() const;

    //- Access to cell faces
    inline faceList& faces(); 

    //- Access to cell faces
    inline const faceList& faces() const;


    // friend Operators

    friend bool operator==(const cellInfo&, const cellInfo&);
    friend bool operator!=(const cellInfo&, const cellInfo&);


    // IOstream Operators

    //- Write cellInfo properties
    friend Ostream& operator<<(Ostream&, const cellInfo&);
    friend Istream& operator>>(Istream&, cellInfo&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellInfoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
