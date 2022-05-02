/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "PANSkOmegaSST.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> PANSkOmegaSST<BasicTurbulenceModel>::PANSkOmegaSST::F1
(
    const volScalarField& CDkOmega
) const
{
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
    );

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/this->betaStar_)*sqrt(kU_)/(omegaU_*this->y_),
                scalar(500)*(this->mu()/this->rho_)/(sqr(this->y_)*omegaU_)
            ),
            (4*this->alphaOmega2_*(fK_/fOmega_))*kU_
            /(CDkOmegaPlus*sqr(this->y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
PANSkOmegaSST<BasicTurbulenceModel>::PANSkOmegaSST::F2() const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/this->betaStar_)*sqrt(kU_)/(omegaU_*this->y_),
            scalar(500)*(this->mu()/this->rho_)/(sqr(this->y_)*omegaU_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
PANSkOmegaSST<BasicTurbulenceModel>::PANSkOmegaSST::F3() const
{
    tmp<volScalarField> arg3 = min
    (
        150*(this->mu()/this->rho_)/(omegaU_*sqr(this->y_)),
        scalar(10)
    );

    return 1 - tanh(pow4(arg3));
}

template<class BasicEddyViscosityModel>
tmp<volScalarField> PANSkOmegaSST<BasicEddyViscosityModel>::F23() const
{
    tmp<volScalarField> f23(F2());

    if (this->F3_)
    {
        f23.ref() *= this->F3();
    }

    return f23;
}


template<class BasicTurbulenceModel>
void PANSkOmegaSST<BasicTurbulenceModel>::correctNut
(
    const volScalarField& S2
    // const volScalarField& F2
)
{
    this->nut_ = this->a1_*kU_/max(this->a1_*omegaU_, this->b1_*F23()*sqrt(S2));
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void PANSkOmegaSST<BasicTurbulenceModel>::correctNut()
{
    // correctNut(2*magSqr(symm(fvc::grad(this->U_))), this->F23());
    correctNut(2*magSqr(symm(fvc::grad(this->U_))));
}

template<class BasicEddyViscosityModel>
tmp<volScalarField::Internal> PANSkOmegaSST<BasicEddyViscosityModel>::GbyNu
(
    const volScalarField::Internal& GbyNu0,
    const volScalarField::Internal& F2,
    const volScalarField::Internal& S2
) const
{
    return min
    (
        GbyNu0,
        (this->c1_/this->a1_)*this->betaStar_*omegaU_()
       *max(this->a1_*omegaU_(), this->b1_*F2*sqrt(S2))
    );
}

template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> PANSkOmegaSST<BasicTurbulenceModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            kU_,
            dimVolume*this->rho_.dimensions()*kU_.dimensions()/dimTime
        )
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> PANSkOmegaSST<BasicTurbulenceModel>::omegaSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omegaU_,
            dimVolume*this->rho_.dimensions()*omegaU_.dimensions()/dimTime
        )
    );
}

template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> PANSkOmegaSST<BasicTurbulenceModel>::Qsas
(
    const volScalarField::Internal& S2,
    const volScalarField::Internal& gamma,
    const volScalarField::Internal& beta
) const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omegaU_,
            dimVolume*this->rho_.dimensions()*omegaU_.dimensions()/dimTime
        )
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
PANSkOmegaSST<BasicTurbulenceModel>::PANSkOmegaSST
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    kOmegaSST<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
    ),

    fEpsilon_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "fEpsilon",
            this->coeffDict_,
            1.0
        )
    ),

    fK_
    (
        IOobject
        (
            IOobject::groupName("fK", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    fOmega_
    (
        IOobject
        (
            "fOmega",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fEpsilon_/fK_
    ),

    kU_
    (
        IOobject
        (
            IOobject::groupName("kU", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    omegaU_
    (
        IOobject
        (
            IOobject::groupName("omegaU", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    bound(kU_, min(fK_)*this->kMin_);
    bound(omegaU_, min(fOmega_)*this->omegaMin_);

    if (type == typeName)
    {
        correctNut();
        this->printCoeffs(type);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool PANSkOmegaSST<BasicTurbulenceModel>::read()
{
    if (kOmegaSST<BasicTurbulenceModel>::read())
    {
        fEpsilon_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicTurbulenceModel>
void PANSkOmegaSST<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    BasicTurbulenceModel::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(symm(tgradU())));
    volScalarField::Internal GbyNu0
    (
        this->type() + ":GbyNu",
        (tgradU() && dev(twoSymm(tgradU())))
    );
    volScalarField::Internal G(this->GName(), nut*GbyNu0);

    // Update omegaU and G at the wall
    omegaU_.boundaryFieldRef().updateCoeffs();

    volScalarField CDkOmega
    (
        (2*this->alphaOmega2_*(fOmega_/fK_))*
        (fvc::grad(kU_) & fvc::grad(omegaU_))/omegaU_
    );

    volScalarField F1(this->F1(CDkOmega));
    volScalarField F23(this->F23());

    {
        volScalarField::Internal gamma(this->gamma(F1));
        volScalarField::Internal beta(this->beta(F1));
        volScalarField::Internal betaL
        (
            gamma*this->betaStar_ - (gamma *this->betaStar_/fOmega_)
            + (beta/fOmega_)
        );
        GbyNu0 = GbyNu(GbyNu0, F23(), S2());


        // Unresolved Turbulent frequency equation
        tmp<fvScalarMatrix> omegaUEqn
        (
            fvm::ddt(alpha, rho, omegaU_)
          + fvm::div(alphaRhoPhi, omegaU_)
          - fvm::laplacian(alpha*rho*DomegaUEff(F1), omegaU_)
         ==
            alpha()*rho()*gamma*GbyNu0
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, omegaU_)
          - fvm::Sp(alpha()*rho()*betaL*omegaU_(), omegaU_)
          - fvm::SuSp
            (
                alpha()*rho()*(F1() - scalar(1))*CDkOmega()/omegaU_(),
                omegaU_
            )
          + Qsas(S2(), gamma, beta)
          // + omegaSource()
          + fvOptions(alpha, rho, omegaU_)
        );

        omegaUEqn.ref().relax();
        fvOptions.constrain(omegaUEqn.ref());
        omegaUEqn.ref().boundaryManipulate(omegaU_.boundaryFieldRef());
        solve(omegaUEqn);
        fvOptions.correct(omegaU_);
        bound(omegaU_, min(fOmega_)*this->omegaMin_);
;    }

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kUEqn
    (
        fvm::ddt(alpha, rho, kU_)
      + fvm::div(alphaRhoPhi, kU_)
      - fvm::laplacian(alpha*rho*DkUEff(F1), kU_)
     ==
        alpha()*rho()*min(G, (this->c1_*this->betaStar_)*kU_()*omegaU_())
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, kU_)
      - fvm::Sp(alpha()*rho()*this->betaStar_*omegaU_, kU_)
      // + kSource()
      + fvOptions(alpha, rho, kU_)
    );

    tgradU.clear();

    kUEqn.ref().relax();
    fvOptions.constrain(kUEqn.ref());
    solve(kUEqn);
    fvOptions.correct(kU_);
    bound(kU_, min(fK_)*this->kMin_);


    // Calculation of total Turbulent kinetic energy and Frequency
    this->k_ = kU_/fK_;
    this->omega_ = omegaU_/fOmega_;

    correctNut(S2);

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //

// ************************************************************************* //
