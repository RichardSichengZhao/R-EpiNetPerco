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

## Required Packages and Environment: pracma gsl deSolve

## Main. R
Including basic functions for data sorting & formatting, PGFs & modeling tool functions and all prediction functions for percolation methods.
Except the uncertain corrections in Multiple type network percolation method, the functions here are all tested and verified with original article results.
Data sorting & formating functions: DegreeDF DDistPK

PGFs & Modeling tools functions: PGFG0 PGFD1G0 PGFD2G0 Kn PGFG1 PGFD1G1 ueqn  

Typical bond percolation method prediction: Tconst_Neman Tconst_exp TypProc

Multiple type network percolation method prediction (with some uncertain corrections of the original method): P2u_polar P2s_polar Pl2_polar TypProc_Polar P2_leaky Q2_leaky Tmatrix PGFdf1_leaky PGFdf2_leaky s_leaky fsys_soln ES_leaky TypProc_leaky

Modified network percolation method prediction: ModProc_CM

## Random Network Generator with Specified Degree Sequence
Currently, there is not a well reconized optimal/best alogrithm to realize networks generally under all specified degree sequence.
There exist many different algorithms under several different theories and ideas.
Based on the core mathematical theories and idea, these methods can be generally classified into the following 4 family with each of its own pros and cons: Configuration model algorithms, Sequential algorithm, Switching-based algorithms and Markov chain algorithms.
There are multiple algorithms in each class optimized to fit different scenario and there are also algorithms using mixed ideas from multiple families.
So it depends on the reseacher/user to consider the trade offs of pros and cons case by case and determine which algorithm is optimal for the given 

As for now, the pacakge originaly targets on providing better interpreation of the percolation process methods in epidemiology networks and trying to extend their model applications to a wider range of problems with different assumptions. 
So we decided to start with the family of Sequential algorithms, considering its following advatages:

1.   Its interperation is simpler and relatively close to the intuitive understanding for percolation process: both assuming "asymptotical independence" and considering a degree related random paring process when generating edges.
2.   The assumptions for our epidemiology network models also fits better with its optimal condition: large sparse networks.
3.   It provides a relatively high calculation coeffiency among all families, especially for large sparse networks in general.
4.   It guaranteed that every possible graph has a postive possiblity to appear as the generation result.

However, the disadvantage of sequential algorithms is that it is harder to prove that its results is uniformly distributed among all possilbe random graphs, which is required by application and interpretion of the percolation method.
It might require extra sampling techniques on simulation results or extra conditions on original degree sequece 

In the package, we plan to provide 5 different sequential algorithms, considering the method illustration purpose and optimal performances under several possible application scenrario.
In each algorithms, provides with and Erdos-Gallai condition checking function modified to fit each algorithm.

These algorithm are temporary, we are still developing these generators while comparing their performance (and with other alorithms) under testing, to figure out if there is an optimal choice for our models, or if there is some systematical principle to choose algorithms under given condition.

### Algorithm 1 : Original algorithm developed by us, inspired by Newman's algorithm. 
Its is fastest in all scenario, but has random failure, lack of rigorous analysis of distribution result.
It is used to quickly generate relatively larger size networks for testing and verification, and also to illustrate Newman's idea for comparision.

### Algorithm 2.0 : Sequential importance algorithm with importance sampling technique, developed with some modification, based on J. Blitzstein & P. Diaconis.
Its mainly for illustration the idea of Blitzstein and Diaconis who originally develope the first squential algorithm for general degree sequences.
It come with an extra sampling procedure to adjust the simulation result to mathematically close to uniform distributed, while increasing some calculation cost.
An optimized version with several modification are provided by Algorithm 2.1

### Algorithm 2.1 (In Developing): Optimized version of Algorithm 2.0
With several modification on Algorithm 2.0, it proiveds better calculation efficiency than 2.0.
From testing, it performs best when generating network not so "sparse"---for degree sequences with higher mean or with longer tail (like powerlaw distribtion).
The corresponding version of importance sampling technique are in developing now.

### Algorithm 3 : Mixture algorithm of 1 and 2


The generation and epidemic simulation are built based on sparce matrix at the beginning.
While it is easiear to understand and modify, the matrix data form in R might not be a good choice to store data, considering the limitation of comupting resourses, especially RAM, so the network size is highly limited.
It might also slows down the calculation speed when loops and iterations are involved.
We are currently working on, together with collaborating data scientists and programmers, is to optimize the data structure in these algorithms to enlarge the capacity of network size for same given RAM limitation and perhaps accelate the calculation.
Another ongoing improvment will be optimize these generators and the whole package for parallel computing in multple servers and computing clusters, so simulation and verification involving large amount of simulations could be finished faster and on larger network size. 



















