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
    Foam::incompTwoPhaseMixture

Description
    A two-phase incomp transportModel

SourceFiles
    incompTwoPhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef incompTwoPhaseMixture_H
#define incompTwoPhaseMixture_H

#include "IOdictionary.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class incompTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class incompTwoPhaseMixture
:
    public IOdictionary    
{
protected:

    // Protected data        

        const volScalarField& alpha1_;                

        dimensionedScalar rho1_;
        dimensionedScalar rho0_;
    
        dimensionedScalar mu1_;
        dimensionedScalar mu0_;

        volScalarField mu_;        


    // Private Member Functions

        //- Calculate and return the laminar viscosity
        void calcMu();


public:

    // Constructors

        //- Construct from components
        incompTwoPhaseMixture
        (
            const volScalarField&            
        );


    //- Destructor
    virtual ~incompTwoPhaseMixture()
    {}


    // Member Functions

        //- Return const-access to phase1 viscosity
        const dimensionedScalar& mu1() const
        {
            return mu1_;
        }

        //- Return const-access to phase0 viscosity
        const dimensionedScalar& mu0() const
        {
            return mu0_;
        };

        //- Return const-access to phase1 density
        const dimensionedScalar& rho1() const
        {
            return rho1_;
        }

        //- Return const-access to phase0 density
        const dimensionedScalar& rho0() const
        {
            return rho0_;
        };

        //- Return the dynamic laminar viscosity
        const volScalarField& mu() const;        

        //- Correct the laminar viscosity
        virtual void correct()
        {
            calcMu();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
