/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "Vreman.H"
#include "fvModels.H"
#include "fvConstraints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<volScalarField> Vreman<BasicMomentumTransportModel>::k
(
    const tmp<volTensorField>& gradU
) const
{
    volSymmTensorField D(symm(gradU));

    volScalarField a(this->Ce_/this->delta());
    volScalarField b((2.0/3.0)*tr(D));
    volScalarField c(2*this->Ck_*this->delta()*(dev(D) && D));

    return volScalarField::New
    (
        this->groupName("k"),
        sqr((-b + sqrt(sqr(b) + 4*a*c))/(2*a))
    );
}


template<class BasicMomentumTransportModel>
void Vreman<BasicMomentumTransportModel>::correctNut()
{
    volScalarField k(this->k(fvc::grad(this->U_)));

    const volTensorField alp(fvc::grad(this->U_));

    const volScalarField alpalp(alp && alp); //double inner product of two tensors
    const volTensorField beta(sqr(this->delta())*alp.T() & alp); //inner product of two tensors

    volScalarField betaxx(beta.component(tensor::XX));
    volScalarField betayy(beta.component(tensor::YY));
    volScalarField betazz(beta.component(tensor::ZZ));
    volScalarField betaxy(beta.component(tensor::XY));
    volScalarField betaxz(beta.component(tensor::XZ));
    volScalarField betayz(beta.component(tensor::YZ));
    volScalarField Bbeta(betaxx*betayy - sqr(betaxy)
                        + betaxx*betazz - sqr(betaxz)
                        + betayy*betazz - sqr(betayz));

    const dimensionedScalar alpalp_small(
                                        alpalp.dimensions(),
                                        small
                                      );

    const dimensionedScalar Bbeta_small(
                                        Bbeta.dimensions(),
                                        small
                                      );

    this->nut_ = 2.5*sqr(this->Cs_)*sqrt(max(Bbeta, Bbeta_small)/max(alpalp, alpalp_small));
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
Vreman<BasicMomentumTransportModel>::Vreman
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity,
    const word& type
)
:
    LESeddyViscosity<BasicMomentumTransportModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),

    Cs_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cs",
            this->coeffDict_,
            0.094
        )
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool Vreman<BasicMomentumTransportModel>::read()
{
    return LESeddyViscosity<BasicMomentumTransportModel>::read();
}


template<class BasicMomentumTransportModel>
void Vreman<BasicMomentumTransportModel>::correct()
{
    LESeddyViscosity<BasicMomentumTransportModel>::correct();
    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
