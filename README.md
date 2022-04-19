# PANS_OpenFOAM
MSc graduation thesis

Built for OpenFOAM v2112

This PANS solver(s) is for f_k fixed in both space and time which is backed in many literature to perform better than f_k varying in space and/or time.

The primary focus is on the k-omega SST model which has been validated whereas the k-epsilon model isn't. 

Derivation on the PANS equations used will be published as part of thesis in the near future.

Secondly, k-corrective-frozen-PANS approach which is based on Discovery of Algebraic Reynolds-Stress Models Using Sparse Symbolic Regression (2020) by Schmelzer, Dwight and Cinnella (10.1007/s10494-019-00089-x) is implemented with slight modification.

Instructions on how to implement your own turbulence model can be found in a course by Nilsson .H. The link: http://www.tfd.chalmers.se/~hani/kurser/OS_CFD_2020/lectureNotes/implementTurbulenceModel.pdf
