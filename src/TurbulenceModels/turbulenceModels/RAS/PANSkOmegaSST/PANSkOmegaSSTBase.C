/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
#include "PANSkOmegaSSTBase.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"
// #include <iostream>
// SMALL = 1e-15 (predefined in C++)
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicEddyViscosityModel>
tmp<volScalarField> PANSkOmegaSSTBase<BasicEddyViscosityModel>::F1
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
                (scalar(1)/betaStar_)*sqrt(kU_)/(omegaU_*y_),
                scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omegaU_)
            ),
            (4*alphaOmega2_*(fK_/fOmega_))*kU_/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}


template<class BasicEddyViscosityModel>
tmp<volScalarField> PANSkOmegaSSTBase<BasicEddyViscosityModel>::F2() const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/this->betaStar_)*sqrt(kU_)/(omegaU_*y_),
            scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omegaU_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}


template<class BasicEddyViscosityModel>
tmp<volScalarField> PANSkOmegaSSTBase<BasicEddyViscosityModel>::F3() const
{
    tmp<volScalarField> arg3 = min
    (
        150*(this->mu()/this->rho_)/(omegaU_*sqr(y_)),
        scalar(10)
    );

    return 1 - tanh(pow4(arg3));
}

template<class BasicEddyViscosityModel>
tmp<volScalarField> PANSkOmegaSSTBase<BasicEddyViscosityModel>::F23() const
{
    tmp<volScalarField> f23(F2());

    if (F3_)
    {
        f23.ref() *= F3();
    }

    return f23;
}


template<class BasicEddyViscosityModel>
void PANSkOmegaSSTBase<BasicEddyViscosityModel>::correctNut
(
    const volScalarField& S2
)
{
    this->nut_ = a1_*kU_/max(a1_*omegaU_, b1_*F23()*sqrt(S2));
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

}

template<class BasicEddyViscosityModel>
void PANSkOmegaSSTBase<BasicEddyViscosityModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))));
}

template<class BasicEddyViscosityModel>
tmp<volScalarField::Internal> PANSkOmegaSSTBase<BasicEddyViscosityModel>::Pk
(
    const volScalarField::Internal& G
) const
{
    return min(G, (c1_*betaStar_)*this->kU_()*this->omegaU_());
}


template<class BasicEddyViscosityModel>
tmp<volScalarField::Internal>
PANSkOmegaSSTBase<BasicEddyViscosityModel>::epsilonByk
(
    const volScalarField& F1,
    const volTensorField& gradU
) const
{
    return betaStar_*omegaU_();
}


template<class BasicEddyViscosityModel>
tmp<volScalarField::Internal> PANSkOmegaSSTBase<BasicEddyViscosityModel>::GbyNu
(
    const volScalarField::Internal& GbyNu0,
    const volScalarField::Internal& F2,
    const volScalarField::Internal& S2
) const
{
    return min
    (
        GbyNu0,
        (c1_/a1_)*betaStar_*omegaU_()
       *max(a1_*omegaU_(), b1_*F2*sqrt(S2))
    );
}


template<class BasicEddyViscosityModel>
tmp<fvScalarMatrix> PANSkOmegaSSTBase<BasicEddyViscosityModel>::kSource() const
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


template<class BasicEddyViscosityModel>
tmp<fvScalarMatrix> PANSkOmegaSSTBase<BasicEddyViscosityModel>::omegaSource() const
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

template<class BasicEddyViscosityModel>
tmp<fvScalarMatrix> PANSkOmegaSSTBase<BasicEddyViscosityModel>::Qsas
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

template<class BasicEddyViscosityModel>
PANSkOmegaSSTBase<BasicEddyViscosityModel>::PANSkOmegaSSTBase
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    BasicEddyViscosityModel
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    alphaK1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            0.85
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaK2",
            this->coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmega1",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            0.856
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "gamma1",
            this->coeffDict_,
            5.0/9.0
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "gamma2",
            this->coeffDict_,
            0.44
        )
    ),
    beta1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "beta1",
            this->coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "beta2",
            this->coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "a1",
            this->coeffDict_,
            0.31
        )
    ),
    b1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "b1",
            this->coeffDict_,
            1.0
        )
    ),
    c1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "c1",
            this->coeffDict_,
            10.0
        )
    ),
    F3_
    (
        Switch::getOrAddToDict
        (
            "F3",
            this->coeffDict_,
            false
        )
    ),

    y_(wallDist::New(this->mesh_).y()),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    decayControl_
    (
        Switch::getOrAddToDict
        (
            "decayControl",
            this->coeffDict_,
            false
        )
    ),
    kInf_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "kInf",
            this->coeffDict_,
            k_.dimensions(),
            0
        )
    ),
    omegaInf_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "omegaInf",
            this->coeffDict_,
            omega_.dimensions(),
            0
        )
    )


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
            IOobject::NO_WRITE
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
            //IOobject::MUST_READ,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->k_*fK_,
        this->k_.boundaryField().types()
    ),

    omegaU_
    (
        IOobject
        (
            IOobject::groupName("omegaU", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            //IOobject::MUST_READ,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->omega_*fOmega_,
        this->omega_.boundaryField().types()
    )
{
    bound(kU_, min(fK_)*this->kMin_);    
    bound(omegaU_, max(fOmega_)*this->omegaMin_);

    setDecayControl(this->coeffDict_);

    // if (type == typeName)
    // {
    //     correctNut();
    //     this->printCoeffs(type);
    // }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicEddyViscosityModel>
void PANSkOmegaSSTBase<BasicEddyViscosityModel>::setDecayControl
(
    const dictionary& dict
)
{
    decayControl_.readIfPresent("decayControl", dict);

    if (decayControl_)
    {
        kInf_.read(dict);
        omegaInf_.read(dict);

        Info<< "    Employing decay control with kInf:" << kInf_
            << " and omegaInf:" << omegaInf_ << endl;
    }
    else
    {
        kInf_.value() = 0;
        omegaInf_.value() = 0;
    }
}

template<class BasicEddyViscosityModel>
bool PANSkOmegaSSTBase<BasicEddyViscosityModel>::read()
{
    if (BasicEddyViscosityModel::read())
    {
        fEpsilon_.readIfPresent(this->coeffDict());
        alphaK1_.readIfPresent(this->coeffDict());
        alphaK2_.readIfPresent(this->coeffDict());
        alphaOmega1_.readIfPresent(this->coeffDict());
        alphaOmega2_.readIfPresent(this->coeffDict());
        gamma1_.readIfPresent(this->coeffDict());
        gamma2_.readIfPresent(this->coeffDict());
        beta1_.readIfPresent(this->coeffDict());
        beta2_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());
        a1_.readIfPresent(this->coeffDict());
        b1_.readIfPresent(this->coeffDict());
        c1_.readIfPresent(this->coeffDict());
        F3_.readIfPresent("F3", this->coeffDict());

        setDecayControl(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicEddyViscosityModel>
void PANSkOmegaSSTBase<BasicEddyViscosityModel>::correct()
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

    BasicEddyViscosityModel::correct();
    
    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(symm(tgradU())));
    volScalarField::Internal GbyNu0
    (
        this->type() + ":GbyNu",
        (tgradU() && dev(twoSymm(tgradU())))
    );
    volScalarField::Internal G(this->GName(), nut*GbyNu);
    // tgradU.clear();

    // Update omegaU and G at the wall
    omegaU_.boundaryFieldRef().updateCoeffs();

    volScalarField CDkOmega
    (
        (2*alphaOmega2_*(fK_/fOmega_))*(fvc::grad(kU_) & fvc::grad(omegaU_))/omegaU_
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
          + alpha()*rho()*betaL*sqr(omegaInf_)
          + Qsas(S2, gamma, betaL)
          + omegaSource()
          + fvOptions(alpha, rho, omegaU_)
        );
        // std::cout << "Pass 5" << "\n";

        omegaUEqn.ref().relax();
        fvOptions.constrain(omegaUEqn.ref());
        omegaUEqn.ref().boundaryManipulate(omegaU_.boundaryFieldRef());
        solve(omegaUEqn);
        fvOptions.correct(omegaU_);
        bound(omegaU_, max(fOmega_)*this->omegaMin_);
        // bound(omegaU_, fOmega_*this->omegaMin_);
    }

    // Turbulent kinetic energy equation
    // std::cout << "Begin 6" << "\n";
    tmp<fvScalarMatrix> kUEqn
    (
        fvm::ddt(alpha, rho, kU_)
      + fvm::div(alphaRhoPhi, kU_)
      - fvm::laplacian(alpha*rho*DkUEff(F1), kU_)
     ==
        // min(alpha*rho*G, (this->c1_*this->betaStar_)*alpha()*rho()*kU_()*omegaU_())
        alpha()*rho()*Pk(G)
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, kU_)
      - fvm::Sp(alpha()*rho()*epsilonByk(F1, tgradU()), kU_)
      + alpha()*rho()*betaStar_*omegaInf_*kInf_
      + kSource()
      + fvOptions(alpha, rho, kU_)
    );
    // std::cout << "Pass 6" << "\n";
    tgradU.clear();

    kUEqn.ref().relax();
    fvOptions.constrain(kUEqn.ref());
    solve(kUEqn);
    fvOptions.correct(kU_);
    bound(kU_, min(fK_)*this->kMin_);
    // bound(kU_, fK_*this->kMin_);


    // Calculation of Turbulent kinetic energy and Frequency
    k_ = kU_/fK_;
    // this->k_.correctBoundaryConditions();

    omega_ = omegaU_/fOmega_;
    // this->omega_.correctBoundaryConditions(); //does not seem to affect the result

    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    correctNut(S2);

// FOR VARIABLE fK
    // // Recalculate fK with new kU and epsilonU
    
    // // Calculate the unresolved turbulence length scale
    // volScalarField lu_
    // (
    //     sqrt(kU_)/(this->betaStar_*omegaU_)
    // );

    // // update fK
    // fK_.primitiveFieldRef() = min
    // (
    //     max(sqrt(1.0/this->betaStar_)*pow(delta()/lu_,2.0/3.0),loLim_),
    //     uLim_
    // );

    // // update fOmega
    // fOmega_ = fEpsilon_/fK_;
}
// {
//     if (type == typeName)
//     {
//         this->printCoeffs(type);
//     }
// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
