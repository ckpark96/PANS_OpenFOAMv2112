/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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
    Foam::kOmegaSSTBase

Description
    Base class implementation of the k-omega-SST turbulence model for
    incompressible and compressible flows.

    Turbulence model described in:
    \verbatim
        Menter, F. R. & Esch, T. (2001).
        Elements of Industrial Heat Transfer Prediction.
        16th Brazilian Congress of Mechanical Engineering (COBEM).
    \endverbatim

    with updated coefficients from
    \verbatim
        Menter, F. R., Kuntz, M., and Langtry, R. (2003).
        Ten Years of Industrial Experience with the SST Turbulence Model.
        Turbulence, Heat and Mass Transfer 4, ed: K. Hanjalic, Y. Nagano,
        & M. Tummers, Begell House, Inc., 625 - 632.
    \endverbatim

    but with the consistent production terms from the 2001 paper as form in the
    2003 paper is a typo, see
    \verbatim
        http://turbmodels.larc.nasa.gov/sst.html
    \endverbatim

    and the addition of the optional F3 term for rough walls from
    \verbatim
        Hellsten, A. (1998).
        Some Improvements in Menter's k-omega-SST turbulence model
        29th AIAA Fluid Dynamics Conference, AIAA-98-2554.
    \endverbatim

    and the optional decay control from:
    \verbatim
        Spalart, P. R. and Rumsey, C. L. (2007).
        Effective Inflow Conditions for Turbulence Models in Aerodynamic
        Calculations
        AIAA Journal, 45(10), 2544 - 2553.
    \endverbatim

    Note that this implementation is written in terms of alpha diffusion
    coefficients rather than the more traditional sigma (alpha = 1/sigma) so
    that the blending can be applied to all coefficients in a consistent
    manner.  The paper suggests that sigma is blended but this would not be
    consistent with the blending of the k-epsilon and k-omega models.

    Also note that the error in the last term of equation (2) relating to
    sigma has been corrected.

    Wall-functions are applied in this implementation by using equations (14)
    to specify the near-wall omega as appropriate.

    The blending functions (15) and (16) are not currently used because of the
    uncertainty in their origin, range of applicability and that if y+ becomes
    sufficiently small blending u_tau in this manner clearly becomes nonsense.

    The default model coefficients are
    \verbatim
        kOmegaSSTBaseCoeffs
        {
            alphaK1         0.85;
            alphaK2         1.0;
            alphaOmega1     0.5;
            alphaOmega2     0.856;
            beta1           0.075;
            beta2           0.0828;
            betaStar        0.09;
            gamma1          5/9;
            gamma2          0.44;
            a1              0.31;
            b1              1.0;
            c1              10.0;
            F3              no;

            // Optional decay control
            decayControl    yes;
            kInf            \<far-field k value\>;
            omegaInf        \<far-field omega value\>;
        }
    \endverbatim

SourceFiles
    PANSkOmegaSST.C

SeeAlso
    kOmegaSSTBase.H

\*---------------------------------------------------------------------------*/

#ifndef frozenPANSkOmegaSST_H
#define frozenPANSkOmegaSST_H

#include "kOmegaSST.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                     Class PANSkOmegaSST Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class frozenPANSkOmegaSST
:
    public kOmegaSST<BasicTurbulenceModel>
{
    // Private Member Functions

        //- No copy construct
        frozenPANSkOmegaSST(const frozenPANSkOmegaSST&) = delete;

        //- No copy assignment
        void operator=(const frozenPANSkOmegaSST&) = delete;


protected:

    // Protected data


        // PANS coefficients

            dimensionedScalar fEpsilon_;
            volScalarField fK_;
            volScalarField fOmega_;

          // LES fields

              volScalarField k_LES_;
              volScalarField kU_LES_;
              volSymmTensorField tauij_LES_;
              volSymmTensorField tauijU_LES_;
              volSymmTensorField aijU_LES_;
              volSymmTensorField bijU_LES_;
              volScalarField PkULES_;

          // Fields to solve for

              volScalarField omegaU_;
              volScalarField kUDeficit_;
              volSymmTensorField bijUDelta_;

          // Misc fields

              //- Wall distance
              //  Note: different to wall distance in parent RASModel
              //  which is for near-wall cells only (TODO: check this)
              volScalarField y_;
              volTensorField gradU_;  // Needed for output (to compute S, W etc.)
              volVectorField gradkU_LES_;
              volVectorField gradomegaU_;

    // Private Member Functions

        virtual tmp<volScalarField> F1(const volScalarField& CDkOmega) const;
        virtual tmp<volScalarField> F2() const;
        virtual tmp<volScalarField> F3() const;
        virtual tmp<volScalarField> F23() const;

        // void correctNut(const volScalarField& S2, const volScalarField& F2);
        virtual void correctNut(const volScalarField& S2);

    // Protected Member Functions

        virtual void correctNut();

        //- Return G/nu
        virtual tmp<volScalarField::Internal> GbyNu
        (
            const volScalarField::Internal& GbyNu0,
            const volScalarField::Internal& F2,
            const volScalarField::Internal& S2
        ) const;

        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> omegaSource() const;

        //- PANS omega source
        virtual tmp<fvScalarMatrix> Qsas
        (
            const volScalarField::Internal& S2,
            const volScalarField::Internal& gamma,
            const volScalarField::Internal& beta
        ) const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("frozenPANSkOmegaSST");


    // Constructors

        //- Construct from components
        frozenPANSkOmegaSST
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~frozenPANSkOmegaSST() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for unresolved k
        tmp<volScalarField> DkUEff(const volScalarField& F1) const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkUEff",
                    (fOmega_/fK_)*this->alphaK(F1)*this->nut_ + this->nu()
                )
            );
        }

        //- Return the effective diffusivity for unresolved omega
        tmp<volScalarField> DomegaUEff(const volScalarField& F1) const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DomegaUEff",
                    (fOmega_/fK_)*this->alphaOmega(F1)*this->nut_ + this->nu()
                )
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return this->k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "epsilon",
                        this->mesh_.time().timeName(),
                        this->mesh_
                    ),
                    this->betaStar_*this->k_*this->omega_,
                    this->omega_.boundaryField().types()
                )
            );
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return this->omega_;
        }


        //- Return the unresolved turbulence kinetic energy
        virtual tmp<volScalarField> kU_LES() const
        {
            return kU_LES_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> omegaU() const
        {
            return omegaU_;
        }

        //- Return the unresolved turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilonU() const
        {
            return volScalarField::New
            (
                "epsilon",
                this->betaStar_*kU_LES_*omegaU_,
                omegaU_.boundaryField().types()
            );
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "frozenPANSkOmegaSST.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
