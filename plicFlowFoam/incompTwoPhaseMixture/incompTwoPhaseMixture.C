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

#include "incompTwoPhaseMixture.H"
#include "surfaceFields.H"
#include "fvc.H"


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

//- Calculate and return the laminar viscosity
void Foam::incompTwoPhaseMixture::calcMu()
{
    const volScalarField limitedAlpha1
    (
        "limitedAlpha1",
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    // Average dynamic viscosity 
    mu_ = limitedAlpha1*mu1_ + (scalar(1) - limitedAlpha1)*mu0_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompTwoPhaseMixture::incompTwoPhaseMixture
(
    const volScalarField& alpha1    
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            alpha1.time().constant(),
            alpha1.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),         
    alpha1_(alpha1),
    rho1_("rho1", dimDensity, this->lookup("rho1")),
    rho0_("rho0", dimDensity, this->lookup("rho0")),
    mu1_("mu1", dimensionSet(1, -1, -1, 0, 0), this->lookup("mu1")),
    mu0_("mu0", dimensionSet(1, -1, -1, 0, 0), this->lookup("mu0")),
    mu_
    (
        IOobject
        (
            "mu",
            alpha1_.time().timeName(),
            alpha1_.db()
        ),
        alpha1_.mesh(),
        dimensionedScalar("mu", dimensionSet(1, -1, -1, 0, 0), 0),
        calculatedFvPatchScalarField::typeName
    )
{
    bool debug = true;    

    if(debug)
    {
        Info<< "Reading transport properties from file" << nl
            << "Done reading transport properties from file" << nl << endl;
    }

    calcMu();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::volScalarField&
Foam::incompTwoPhaseMixture::mu() const
{
    return mu_;
}


// ************************************************************************* //
