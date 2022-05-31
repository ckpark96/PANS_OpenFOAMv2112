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
#include <vector>
#include <iostream>
#include "frozenInterpPANSkOmegaSST.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> frozenInterpPANSkOmegaSST<BasicTurbulenceModel>::frozenInterpPANSkOmegaSST::F1
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
                (scalar(1)/betaStar_)*sqrt(kU_LES_)/(omegaU_*y_),
                scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omegaU_)
            ),
            (4*alphaOmega2_*(fK_/fOmega_))*kU_LES_
            /(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
frozenInterpPANSkOmegaSST<BasicTurbulenceModel>::frozenInterpPANSkOmegaSST::F2() const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(kU_LES_)/(omegaU_*y_),
            scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omegaU_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
frozenInterpPANSkOmegaSST<BasicTurbulenceModel>::frozenInterpPANSkOmegaSST::F3() const
{
    tmp<volScalarField> arg3 = min
    (
        150*(this->mu()/this->rho_)/(omegaU_*sqr(y_)),
        scalar(10)
    );

    return 1 - tanh(pow4(arg3));
}

template<class BasicEddyViscosityModel>
tmp<volScalarField> frozenInterpPANSkOmegaSST<BasicEddyViscosityModel>::F23() const
{
    tmp<volScalarField> f23(F2());

    if (F3_)
    {
        f23.ref() *= F3();
    }

    return f23;
}


template<class BasicTurbulenceModel>
void frozenInterpPANSkOmegaSST<BasicTurbulenceModel>::correctNut
(
    const volScalarField& S2
    // const volScalarField& F2
)
{
    this->nut_ = a1_*kU_LES_/max(a1_*omegaU_, b1_*F23()*sqrt(S2));
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();

}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void frozenInterpPANSkOmegaSST<BasicTurbulenceModel>::correctNut()
{
    // correctNut(2*magSqr(symm(fvc::grad(this->U_))), this->F23());
    correctNut(2*magSqr(symm(fvc::grad(this->U_))));
}

template<class BasicEddyViscosityModel>
tmp<volScalarField::Internal> frozenInterpPANSkOmegaSST<BasicEddyViscosityModel>::GbyNu
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

template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> frozenInterpPANSkOmegaSST<BasicTurbulenceModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            kU_LES_,
            dimVolume*this->rho_.dimensions()*kU_LES_.dimensions()/dimTime
        )
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> frozenInterpPANSkOmegaSST<BasicTurbulenceModel>::omegaSource() const
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
tmp<fvScalarMatrix> frozenInterpPANSkOmegaSST<BasicTurbulenceModel>::Qsas
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

/////////////// Custom function to aid in interpolation ///////////////

class customClass
{
public:
  const float period_;
  const float timeStep_;
  std::vector<double> times_hifi_;

  customClass();
  ~customClass();

};

customClass::customClass()
:
period_(0.00825617),
timeStep_(1e-4),
times_hifi_(arange<double>(0, round_up(period_,4), round_up(timeStep_,4)))
{
    Info << "Declaring custom class" << endl;
}


customClass::~customClass()
{}

customClass myClass;

template<class BasicTurbulenceModel>
frozenInterpPANSkOmegaSST<BasicTurbulenceModel>::frozenInterpPANSkOmegaSST
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
    eddyViscosity<RASModel<BasicTurbulenceModel>>
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
            IOobject::MUST_READ, //MUST_READ,
            IOobject::NO_WRITE
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

    //========================== LES fields ==============================
    k_LES_
    (
        IOobject
        (
            IOobject::groupName("k_LES", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_
    ),

    kU_LES_
    (
        IOobject
        (
            IOobject::groupName("kU_LES", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        k_LES_ * fK_
    ),

    tauij_LES_
    (
        IOobject
        (
            "tauij_LES",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        // (tauij_LES_post_ - tauij_LES_pre_) / (postTime_ - preTime_) * (currentTime_ - preTime_) + tauij_LES_pre_
        this->mesh_
    ),
    tauijU_LES_
    (
        IOobject
        (
            "tauijU_LES",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tauij_LES_ * fK_ * fK_
    ),
    aijU_LES_
    (
        IOobject
        (
            "aijU_LES",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tauijU_LES_ - ((2.0/3.0)*I)*kU_LES_
    ),
    bijU_LES_
    (
        IOobject
        (
            "bijU_LES",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aijU_LES_ / 2.0 / (kU_LES_ + this->kMin_)
    ),
    PkULES_
    (
        IOobject
        (
            "PkULES",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("PkULES", dimensionSet(0,2,-3,0,0,0,0), 0.0)
    ),

    //========================== Unknown fields - MUST be written ============================
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
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
    ),

    kUDeficit_
    (
        IOobject(
	          "kUDeficit",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("kUDeficit", dimensionSet(0,2,-3,0,0,0,0), 0.0)
    ),

    bijUDelta_
    (
        IOobject
        (
            "bijDelta",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0.0*symm(fvc::grad(this->U_))/omegaU_
    ),

    //========================== Misc ============================
    // y_(
    //    IOobject
    //    (
    //         "walldist",
    //         this->runTime_.timeName(),
    //         this->mesh_,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //    ),
    //    wallDist::New(this->mesh_).y()
    // ),
    gradU_
    (
        IOobject
        (
            "gradU",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::grad(this->U_)
    ),
    gradkU_LES_
    (
        IOobject
        (
            "gradkU_LES",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::grad(kU_LES_)
    ),
    gradomegaU_
    (
       IOobject
       (
           "gradomegaU",
           this->runTime_.timeName(),
           this->mesh_,
           IOobject::NO_READ,
           IOobject::NO_WRITE
       ),
       this->mesh_,
       dimensionedVector("gradomegaU", dimensionSet(0,-1,-1,0,0,0,0), Zero)
    )
{
    if (type == typeName)
    {
        correctNut();
        this->printCoeffs(type);
    }
    bound(kU_LES_, min(fK_)*this->kMin_);
    bound(omegaU_, min(fOmega_)*this->omegaMin_);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool frozenInterpPANSkOmegaSST<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
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
        fEpsilon_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicTurbulenceModel>
void frozenInterpPANSkOmegaSST<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    Info << "current run TIME: " << this->runTime_.value() << endl;

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    volScalarField& omegaU_ = this->omegaU_;
    volScalarField& k_LES_ = this->k_LES_;
    volScalarField& kU_LES_ = this->kU_LES_;

    BasicTurbulenceModel::correct();

    ////////////////// RE-READ FIELDS FOR NEW TIMESTEP ////////////////////

    double currentTime_(this->runTime_.value());
    int lowerIndex_(std::get<0>(searchBounds(myClass.period_, currentTime_, myClass.times_hifi_)));
    int upperIndex_(std::get<1>(searchBounds(myClass.period_, currentTime_, myClass.times_hifi_)));
    double preTime_(myClass.times_hifi_[lowerIndex_]);
    double postTime_(myClass.times_hifi_[upperIndex_]);
    double remainTime_(fmod(currentTime_,myClass.period_) - preTime_);
    Info << "preTime_: " << preTime_ <<endl;
    Info << "postTime_: " << postTime_ <<endl;
    Info << "Reading k_LES_pre" << endl;
    volScalarField k_LES_pre_
    (
        IOobject
        (
            IOobject::groupName("k_LES", U.group()),
            name(preTime_),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_
    );

    Info << "Reading k_LES_post" << endl;
    volScalarField k_LES_post_
    (
        IOobject
        (
            IOobject::groupName("k_LES", U.group()),
            name(postTime_),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_
    );

    Info << "Calculating k_LES and kU_LES" << endl;
    // k_LES_ = (k_LES_post_ - k_LES_pre_) / (postTime_ - preTime_) * (currentTime_ - preTime_) + k_LES_pre_;
    k_LES_ = (k_LES_post_ - k_LES_pre_) / myClass.timeStep_ * remainTime_ + k_LES_pre_;
    kU_LES_ = k_LES_ * fK_;

    Info << "Reading tauij_LES_pre" << endl;
    volSymmTensorField tauij_LES_pre_
    (
        IOobject
        (
            "tauij_LES",
            name(preTime_),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_
    );

    Info << "Reading tauij_LES_post" << endl;
    volSymmTensorField tauij_LES_post_
    (
        IOobject
        (
            "tauij_LES",
            name(postTime_),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_
    );

    Info << "Calculating tauij_LES, tauijU_LES, aijU_LES, bijU_LES and gradkU_LES" << endl;
    tauij_LES_ =(tauij_LES_post_ - tauij_LES_pre_) / myClass.timeStep_ * remainTime_ + tauij_LES_pre_;
    tauijU_LES_ = tauij_LES_ * fK_ * fK_;
    aijU_LES_ = tauijU_LES_ - ((2.0/3.0)*I)*kU_LES_;
    bijU_LES_ = aijU_LES_ / 2.0 / (kU_LES_ + this->kMin_);
    gradkU_LES_ = fvc::grad(kU_LES_);

    ///////////////// END OF RE-READ /////////////////////

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();
    const dimensionedScalar time = this->runTime_;
    Info << "Confirming runTime: " << this->runTime_.value() << endl;

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(symm(tgradU())));
    volScalarField::Internal GbyNu0
    (
        this->type() + ":GbyNu",
        (tgradU() && dev(twoSymm(tgradU())))
    );
    volScalarField::Internal G(this->GName(), nut*GbyNu0);

    // Info << "Finding ERROR 0 " << endl;
    // Production term from HiFi dataset
    PkULES_ = -tauijU_LES_ && tgradU();

    // "Free" temporary variable
    tgradU.clear();
    // Info << "Finding ERROR 1" << endl;

    // Update omegaU and G at the wall
    omegaU_.boundaryFieldRef().updateCoeffs();

    // Info << "Finding ERROR 2" << endl;

    volScalarField CDkOmega
    (
        (2*this->alphaOmega2_*(fOmega_/fK_))*
        (fvc::grad(kU_LES_) & fvc::grad(omegaU_))/omegaU_
    );

    volScalarField F1(this->F1(CDkOmega));
    volScalarField F23(this->F23());

    {
        volScalarField::Internal gamma(this->gamma(F1));
        volScalarField::Internal beta(this->beta(F1));
        volScalarField::Internal betaL
        (
            gamma*betaStar_ - (gamma *betaStar_/fOmega_)
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
            alpha()*rho()*gamma*
            (
        	   // G/nut from HiFi data
        	   PkULES_ *omegaU_()/kU_LES_() // omega/k = 1/nut
          	)
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
        Info << "Finding ERROR 3" << endl;
        omegaUEqn.ref().relax();
        fvOptions.constrain(omegaUEqn.ref());
        omegaUEqn.ref().boundaryManipulate(omegaU_.boundaryFieldRef());
        solve(omegaUEqn);
        fvOptions.correct(omegaU_);
        bound(omegaU_, min(fOmega_)*this->omegaMin_);
;    }

    tgradU.clear();

    bound(kU_LES_, min(fK_)*this->kMin_);

    // Notes:
    //   - Incompressible only due to the change from fvm to fvc (why?)
    //   - Is ddt necessary here?
    //   - kSource() and fvOptions() not included

    // kUDeficit_ refers to R term (correction term in kU-equation)

    kUDeficit_ =  fvc::ddt(alpha, rho*kU_LES_)
                + fvc::div(alphaRhoPhi, kU_LES_)
                - fvc::laplacian(alpha*rho*DkUEff(F1), kU_LES_)
                - alpha()*rho()*PkULES_
                //+ (2.0/3.0)*alpha*rho*divU*k_LES_ // Incompressible: divU = 0
                + fvc::Sp(betaStar_*alpha*rho*omegaU_, kU_LES_); // betastar = C_mu


    // Calculation of Turbulent kinetic energy and Frequency
    omega_ = omegaU_/fOmega_;

    // Re-calculate eddy viscosity (k/omega)
    correctNut(S2);

    // Calculate bijUDelta, the model correction term for RST equation
    bijUDelta_ = bijU_LES_ + nut / kU_LES_ * symm(fvc::grad(this->U_));

    Info << "Last check: current time: " << this->runTime_.value() << endl;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //

// ************************************************************************* //
