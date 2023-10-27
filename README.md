# R-EpiNetPerco
A package in R for modeling infectious diease transmission on population with social and biological heterogenity, based on SIR-compartmental complex epidemiology random networks with percolation process method and its modifications.

The prediction functions are developed based on following works withmodifications:

1. Typical bond percolation method original developed by M.E.J Newman.
2. Multiple type network percolation method by S. Bansal & L.A. Meyers.
3. Modified network percolation process with disease dynamics by J.C. Miller, A.C. Slim & E.M.Volz.

Also include different random network generator with specified degree sequence, based on following algrorithms with modifications:
1. Orignal algorithms of developers.
2. Sequential importance sampling algorithm by J. Blitzstein & P. Diaconis.
3. Asymptotic sequential algorithm by M. Bayati, J.H. Kim & A. Saberi, with modification.

The disease simulation function is developed based on Doob–Gillespie algorithm (Stochastic Simulation Algorithm, SSA) 

The package is developed under R version 4.2.2

Developing, optimizing and testing in progress now! A realted paper for systematic interperation of theoritical background, methods and modifications will be submitted for review soon. 
Looking forward to comments, issues, bug reports, using opinions and discussions from all users. 
Please contact with primary developer (Richard Zhao, Ph.D. candidate, Dept. Math & Stat, Queen's University, Kingston, ON, Canada) by email 20sz11@queensu.ca

## Main. R
Including basic functions for data sorting & formatting, PGFs & modeling tool functions and all prediction functions for percolation methods.

Required Packages and Environment: pracma gsl deSolve

Data sorting & formating functions: DegreeDF DDistPK

PGFs & Modeling tools functions: PGFG0 PGFD1G0 PGFD2G0 Kn PGFG1 PGFD1G1 ueqn  

Typical bond percolation method prediction: Tconst_Neman Tconst_exp TypProc

Multiple type network percolation method prediction (with uncertain correction of the original method): P2u_polar P2s_polar Pl2_polar TypProc_Polar P2_leaky Q2_leaky Tmatrix PGFdf1_leaky PGFdf2_leaky s_leaky fsys_soln ES_leaky TypProc_leaky

Modified network percolation method prediction: ModProc_CM
