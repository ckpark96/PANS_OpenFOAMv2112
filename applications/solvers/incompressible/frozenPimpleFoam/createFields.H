#include "createRDeltaT.H"
// Info<< "createFields: myClass.postTime_: " << myClass.postTime_ << endl;
// Info<< "createFields: myClass.preTime_: " << myClass.preTime_ << endl;
Info << "Check that this goes through only once for t=0" << endl;
Info << "createFields: myClass.currentTime_: " << myClass.currentTime_ << endl;

//////////////// COPY OF FUNCTIONS FROM FROZEN MODEL HEADER FILE

Info << "createFields: Reading field p" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);


// Info<< "createFields: Reading field U_LES_pre" << endl;
// volVectorField U_LES_pre
// (
//     IOobject
//     (
//         "U_LES",
//         name(myClass.preTime_),
//         mesh,
//         IOobject::MUST_READ,
//         IOobject::NO_WRITE
//     ),
//     mesh
// );
// Info<< "createFields: Reading field U_LES_post" << endl;
//
// volVectorField U_LES_post
// (
//     IOobject
//     (
//         "U_LES",
//         name(myClass.postTime_),
//         mesh,
//         IOobject::MUST_READ,
//         IOobject::NO_WRITE
//     ),
//     mesh
// );

Info<< "createFields: Initialising U_LES" << endl;
volVectorField U_LES
(
    IOobject
    (
        "U_LES",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
  );

/*
For t=0 boundary conditions are required and READ BUT these conditions do not allow for math operations
Thus, a small t value of 1e-8 with same field values as t=0 is included
For non-zero time, boundary conditions are NOT READ
*/
// if (runTime.value() == 0)
// {
//   Info<< "createFields: Reading field U_LES for t=0" << endl;
//   volVectorField U_LES
//   (
//       IOobject
//       (
//           "U_LES",
//           runTime.timeName(),
//           mesh,
//           IOobject::MUST_READ,
//           IOobject::NO_WRITE
//       ),
//       mesh
//     );
// } else
// {
//   Info<< "createFields: Calculating field U_LES" << endl;
//
//   volVectorField U_LES
//   (
//       IOobject
//       (
//           "U_LES",
//           runTime.timeName(),
//           mesh,
//           IOobject::NO_READ,
//           IOobject::NO_WRITE
//       ),
//       (U_LES_post - U_LES_pre) / (myClass.postTime_ - myClass.preTime_) * (myClass.currentTime_ - myClass.preTime_) + U_LES_pre
//   );
//   Info<< "createFields: DONE Calculating field U_LES\n" << endl;
// }

// volScalarField &p = p_LES;  // unused variable
volVectorField &U = U_LES;  // As "U" is used in other #includes

#include "createPhi.H"


// label pRefCell = 0;
// scalar pRefValue = 0.0;
// setRefCell(p, pimple.dict(), pRefCell, pRefValue);
// mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi); // Calls: Selecting incompressible transport model Newtonian

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"
