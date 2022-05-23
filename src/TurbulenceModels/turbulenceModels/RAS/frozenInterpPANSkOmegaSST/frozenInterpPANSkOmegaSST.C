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
// #include "temporalInterpolate.H"

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
                (scalar(1)/this->betaStar_)*sqrt(kU_LES_)/(omegaU_*this->y_),
                scalar(500)*(this->mu()/this->rho_)/(sqr(this->y_)*omegaU_)
            ),
            (4*this->alphaOmega2_*(fK_/fOmega_))*kU_LES_
            /(CDkOmegaPlus*sqr(this->y_))
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
            (scalar(2)/this->betaStar_)*sqrt(kU_LES_)/(omegaU_*this->y_),
            scalar(500)*(this->mu()/this->rho_)/(sqr(this->y_)*omegaU_)
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
        150*(this->mu()/this->rho_)/(omegaU_*sqr(this->y_)),
        scalar(10)
    );

    return 1 - tanh(pow4(arg3));
}

template<class BasicEddyViscosityModel>
tmp<volScalarField> frozenInterpPANSkOmegaSST<BasicEddyViscosityModel>::F23() const
{
    tmp<volScalarField> f23(F2());

    if (this->F3_)
    {
        f23.ref() *= this->F3();
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
    this->nut_ = this->a1_*kU_LES_/max(this->a1_*omegaU_, this->b1_*F23()*sqrt(S2));
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

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
        (this->c1_/this->a1_)*this->betaStar_*omegaU_()
       *max(this->a1_*omegaU_(), this->b1_*F2*sqrt(S2))
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

// vectorField  readkU_LES;
// fileName caseDir = "./0";
// IFstream dataStream(caseDir/"k_LES_test");
// dataStream >> readkU_LES;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

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
times_hifi_(arange<double>(round_up(timeStep_,4), round_up(period_,4), round_up(timeStep_,4)))
// times_hifi_()
// times_hifi_(0.0 0.000000001 0.0001 0.0002 0.0003 0.0004 0.0005 0.0006 0.0007 0.0008 0.0009 0.001 0.0011 0.0012 0.0013 0.0014 0.0015 0.0016 0.0017 0.0018 0.0019 0.002 0.0021 0.0022 0.0023 0.0024 0.0025 0.0026 0.0027 0.0028 0.0029 0.003 0.0031 0.0032 0.0033 0.0034 0.0035 0.0036 0.0037 0.0038 0.0039 0.004 0.0041 0.0042 0.0043 0.0044 0.0045 0.0046 0.0047 0.0048 0.0049 0.005 0.0051 0.0052 0.0053 0.0054 0.0055 0.0056 0.0057 0.0058 0.0059 0.006 0.0061 0.0062 0.0063 0.0064 0.0065 0.0066 0.0067 0.0068 0.0069 0.007 0.0071 0.0072 0.0073 0.0074 0.0075 0.0076 0.0077 0.0078 0.0079 0.008 0.0081 0.0082 )
{
    times_hifi_.insert(times_hifi_.begin(),1e-9);
    times_hifi_.insert(times_hifi_.begin(),0);
    Info << "checking TIMES: " << times_hifi_ << endl;
    // std::cout << "dt = " << timeStep_ << endl;
}

// auto [value1, value2] = searchBounds(pe, 1.125, times);

customClass::~customClass()
{}

customClass myClass;
// auto times = arange<double>(0, round_up(myClass.period_,4), round_up(myClass.timeStep_,4));
// auto times_hifi_ = arange<double>(round_up(myClass.timeStep_,4), round_up(myClass.period_,4), round_up(myClass.timeStep_,4));
// times_hifi_.insert(times_hifi.begin(),1e-9);
// times_hifi_.insert(times_hifi.begin(),0);
// myClass.times_hifi_ = { 0.0 0.000000001 0.0001 0.0002 0.0003 0.0004 0.0005 0.0006 0.0007 0.0008 0.0009 0.001 0.0011 0.0012 0.0013 0.0014 0.0015 0.0016 0.0017 0.0018 0.0019 0.002 0.0021 0.0022 0.0023 0.0024 0.0025 0.0026 0.0027 0.0028 0.0029 0.003 0.0031 0.0032 0.0033 0.0034 0.0035 0.0036 0.0037 0.0038 0.0039 0.004 0.0041 0.0042 0.0043 0.0044 0.0045 0.0046 0.0047 0.0048 0.0049 0.005 0.0051 0.0052 0.0053 0.0054 0.0055 0.0056 0.0057 0.0058 0.0059 0.006 0.0061 0.0062 0.0063 0.0064 0.0065 0.0066 0.0067 0.0068 0.0069 0.007 0.0071 0.0072 0.0073 0.0074 0.0075 0.0076 0.0077 0.0078 0.0079 0.008 0.0081 0.0082 };

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
    // period_(0.00825617),
    // timeStep_(1e-4),
    // times_hifi_(arange<double>(round_up(timeStep_,4), round_up(period_,4), round_up(timeStep_,4))),
    currentTime_(this->runTime_.value()),
    lowerIndex_(searchLowerBound(myClass.period_, currentTime_, myClass.times_hifi_)),
    preTime_(myClass.times_hifi_[lowerIndex_]),
    postTime_(myClass.times_hifi_[lowerIndex_+1]),

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
    k_LES_pre_
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
    ),
    k_LES_post_
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
    ),
    k_LES_
    (
        IOobject
        (
            IOobject::groupName("k_LES", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (k_LES_post_ - k_LES_pre_) / (postTime_ - preTime_) * (currentTime_ - preTime_) + k_LES_pre_
    ),


    kU_LES_
    (
        IOobject
        (
            IOobject::groupName("kU_LES", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        k_LES_ * fK_
    ),
    tauij_LES_pre_
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
    ),
    tauij_LES_post_
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
    ),

    tauij_LES_
    (
        IOobject
        (
            "tauij_LES",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (tauij_LES_post_ - tauij_LES_pre_) / (postTime_ - preTime_) * (currentTime_ - preTime_) + tauij_LES_pre_
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
    y_(
       IOobject
       (
            "walldist",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
       ),
       wallDist::New(this->mesh_).y()
    ),
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
        fvc::grad(this->kU_LES_)  // Only compute once - k=k_LES is fixed in frozen
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
    bound(kU_LES_, min(fK_)*this->kMin_);
    bound(omegaU_, min(fOmega_)*this->omegaMin_);

    Info<< "\n name(postTime_),: " << name(postTime_) << endl;
    Info<< "\n name(preTime_),: " << name(preTime_) << endl;

    if (type == typeName)
    {
        correctNut();
        this->printCoeffs(type);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool frozenInterpPANSkOmegaSST<BasicTurbulenceModel>::read()
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
void frozenInterpPANSkOmegaSST<BasicTurbulenceModel>::correct()
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

    // k_LES_pre_
    // (
    //     IOobject
    //     (
    //         IOobject::groupName("k_LES", U.group()),
    //         name(preTime_),
    //         this->mesh_,
    //         IOobject::MUST_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     this->mesh_
    // );
    // k_LES_post_
    // (
    //     IOobject
    //     (
    //         IOobject::groupName("k_LES", U.group()),
    //         name(postTime_),
    //         this->mesh_,
    //         IOobject::MUST_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     this->mesh_
    // );
    //
    // Info << "Calculating k_LES" << endl;
    // k_LES_
    // (
    //     IOobject
    //     (
    //         IOobject::groupName("k_LES", U.group()),
    //         this->runTime_.timeName(),
    //         this->mesh_,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     (k_LES_post_ - k_LES_pre_) / (postTime_ - preTime_) * (currentTime_ - preTime_) + k_LES_pre_
    // );

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();
    // This is the current iteration, 1000, etc.
    const dimensionedScalar time = this->runTime_;
    Info << "runTime: " << this->runTime_.value() << endl;


    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(symm(tgradU())));
    volScalarField::Internal GbyNu0
    (
        this->type() + ":GbyNu",
        (tgradU() && dev(twoSymm(tgradU())))
    );
    volScalarField::Internal G(this->GName(), nut*GbyNu0);


    // Production term from HiFi dataset
    PkULES_ = -tauijU_LES_ && tgradU();

    // "Free" temporary variable
    tgradU.clear();

    // Update omegaU and G at the wall
    omegaU_.boundaryFieldRef().updateCoeffs();

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

    kUDeficit_ = //  fvc::ddt(alpha, rho, k_LES_)
                  fvc::div(alphaRhoPhi, kU_LES_)
                - fvc::laplacian(alpha*rho*DkUEff(F1), kU_LES_)
                - alpha()*rho()*PkULES_
                //+ (2.0/3.0)*alpha*rho*divU*k_LES_ // Incompressible: divU = 0
                + fvc::Sp(this->betaStar_*alpha*rho*omegaU_, kU_LES_); // betastar = C_mu


    // Calculation of Turbulent kinetic energy and Frequency
    this->omega_ = omegaU_/fOmega_;

    // Re-calculate eddy viscosity (k/omega)
    correctNut(S2);

    // Calculate bijUDelta, the model correction term for RST equation
    bijUDelta_ = bijU_LES_ + nut / kU_LES_ * symm(fvc::grad(this->U_));

    // word tim=this->runTime_.timeName();
    // float tim2 = this->runTime_.value();
    // word tim3 = name(tim2);
    //
    // // Foam::word tim=this->runTime_.timeName();
    // temporalInterpolate(tim2);
    // temporalInterpolate(tim3);

    // customClass myClass;
    // float pe = myClass.period_;
    // printStuff(pe);
    // auto times = arange<double>(0, round_up(myClass.period_,4), round_up(myClass.timeStep_,4));
    // printStuff(times);
    // auto [value1, value2] = searchBounds(pe, 1.125, times);
    // printStuff(value1);
    // printStuff(value2);

    //
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //

// ************************************************************************* //
