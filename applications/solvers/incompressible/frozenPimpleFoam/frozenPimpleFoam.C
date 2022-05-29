/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

Application
    pimpleFoam.C

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, turbulent flow of Newtonian fluids
    on a moving mesh.

    \heading Solver details
    The solver uses the PIMPLE (merged PISO-SIMPLE) algorithm to solve the
    continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \ddt{\vec{U}} + \div \left( \vec{U} \vec{U} \right) - \div \gvec{R}
          = - \grad p + \vec{S}_U
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
        \vec{R} | Stress tensor
        \vec{S}_U | Momentum source
    \endvartable

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
        \<turbulence fields\> | As required by user selection
    \endplaintable

Note
   The motion frequency of this solver can be influenced by the presence
   of "updateControl" and "updateInterval" in the dynamicMeshDict.

\*---------------------------------------------------------------------------*/

// #include "vector"

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

//--------------------------OWN STUFF------------------------//

// Similar to numpy.arange()
template<typename T>
std::vector<T> arange(T start, T stop, T step)
{
    std::vector<T> values;
    for (T value = start; value < stop; value += step)
        values.push_back(value);
    return values;
}

// Round double to certain decimal places
double round_up(double value, int decimal_places)
{
    const double multiplier = std::pow(10.0, decimal_places);
    return std::ceil(value * multiplier) / multiplier;
}

// Search for the smallest time in that is larger than the instantaneous time
int searchLowerBound(double per, double val, std::vector<double> vec)
{
  double remain;
  int lowerBoundIndex = 0;
  for(std::size_t i = 0; i < vec.size(); ++i)
    {
      // remain = remainder(val, per); // not absolute remainder but scaled with respect to the divider
      remain = fmod(val, per);
      if (vec[i] > remain)
      {
        Info << "vec[i]: " << vec[i] << endl;
        Info << "remain: " << remain << endl;
        lowerBoundIndex = i-1;
        break;
      }
    }
  return lowerBoundIndex;
}

class customClass
{
public:
  const float period_;
  const float timeStep_;
  double currentTime_;
  int lowerIndex_;
  double preTime_;
  double postTime_;

  customClass();
  ~customClass();

};

customClass::customClass()
:
period_(0.00825617),
timeStep_(1e-4),
currentTime_(0.0),
lowerIndex_(0)
{
    std::cout << "period = " << period_ << endl;
    std::cout << "\ndt = " << timeStep_ << endl;
    Info << "\nextra digit '1' is printed behind the floats but actual values remain" << endl;
}

customClass::~customClass()
{}

customClass myClass;
auto times_hifi = arange<double>(0, round_up(myClass.period_,4), round_up(myClass.timeStep_,4));

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // std::cout << typeid(times_hifi).name() << '\n';
    Info << "\nTIMES: "<< times_hifi << endl;
    // Initialise the preTime and postTime
    Info << "\nThis should only be printed once\n" << endl;
    myClass.preTime_ = times_hifi[0];
    myClass.postTime_ = times_hifi[1];

    argList::addNote
    (
        "Transient solver for incompressible, turbulent flow"
        " of Newtonian fluids on a moving mesh."
    );

    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H" // runTime initialised
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        myClass.currentTime_ = runTime.value();
        myClass.lowerIndex_ = searchLowerBound(myClass.period_, myClass.currentTime_, times_hifi);
        myClass.preTime_ = times_hifi[myClass.lowerIndex_];
        myClass.postTime_ = times_hifi[myClass.lowerIndex_+1];

        Info<< "myClass.currentTime = " << myClass.currentTime_ << endl;
        Info<< "myClass.preTime_: " << myClass.preTime_ << endl;
        Info<< "myClass.postTime_: " << myClass.postTime_ << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                // Do any mesh changes
                mesh.controlledUpdate();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            // #include "UEqn.H"
            Info<< "Reading field U_LES_pre" << endl;

            volVectorField U_LES_pre
            (
                IOobject
                (
                    "U_LES",
                    name(myClass.preTime_),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );
            Info<< "Reading field U_LES_post" << endl;

            volVectorField U_LES_post
            (
                IOobject
                (
                    "U_LES",
                    name(myClass.postTime_),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            Info<< "Calculating field U_LES" << endl;
            U_LES = (U_LES_post - U_LES_pre) / (myClass.postTime_ - myClass.preTime_) * (myClass.currentTime_ - myClass.preTime_) + U_LES_pre;
            Info<< "DONE Calculating field U_LES\n" << endl;

            // --- Pressure corrector loop
            // while (pimple.correct())
            // {
            //     //#include "pEqn.H"
            // }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
