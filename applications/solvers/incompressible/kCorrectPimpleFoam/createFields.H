#include "createRDeltaT.H"

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p_LES",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "Reading field U_LES\n" << endl;
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
// volScalarField &p = p_LES;  // As "U" is used in other #includes
volVectorField &U = U_LES;  // As "U" is used in other #includes
#include "createPhi.H"


// label pRefCell = 0;
// scalar pRefValue = 0.0;
// setRefCell(p, pimple.dict(), pRefCell, pRefValue);
// mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"
