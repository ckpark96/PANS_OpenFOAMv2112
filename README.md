# PANS_OpenFOAM
MSc graduation thesis done from Oct 2021 till Aug 2022
Published thesis can be found here: https://www.researchgate.net/publication/364709155_k-corrective_frozen_PANS_A_data-driven_stochastic_turbulence_closure_model
PANS = Partially averaged Navier-Stokes

Built for OpenFOAM v2112

This PANS solver(s) is for f_k fixed in both space and time which is backed in many literature to perform better than f_k varying in space and/or time.

The primary focus is on the k-omega SST model which has been validated whereas the k-epsilon model isn't. 

Derivation on the PANS equations used is included in the link above.

Secondly, k-corrective-frozen-PANS approach which is based on Discovery of Algebraic Reynolds-Stress Models Using Sparse Symbolic Regression (2020) by Schmelzer, Dwight and Cinnella (10.1007/s10494-019-00089-x) is implemented with some modification/correction.

Instructions on how to implement your own turbulence model can be found in a course by Nilsson .H. The link: http://www.tfd.chalmers.se/~hani/kurser/OS_CFD_2020/lectureNotes/implementTurbulenceModel.pdf

All the developed turbulence models can be found under src/TurbulenceModels/turbulenceModels/RAS/
- PANSkOmegaSST: is the original PANS developed purely from the PANS equations without any modification
- frozen PANSkOmegaSST: is the frozen method for steadystate flows solved using PANS
- forzenInterpPANSkOmegaSST: is the main model for the thesis whereby frozen method is implemented for unsteady flow in which HiFi data is implemented for 82 time steps in a single period and interpolation is executed for other times.

The forzenInterpPANSkOmegaSST works hand-in-hand with frozenPimpleFoam which can be found under applications/solvers/incompressible
- It also has the interpolation implemented in the same way
