# R-EpiNetPerco
A package in R for modeling infectious disease transmission on population with social and biological heterogeneity, based on SIR-compartmental complex epidemiology random networks with percolation process method and its modifications.

## Developer & Contributor
Sicheng (Richard) Zhao \[aut, cre\], M.Sc., Ph.D. Candidate, Dept. of Mathematics and Statistics, Queen's University, Kingston, ON, Canada (20sz11@queensu.ca)

Felicia M. G. Magpantay \[aut\], Ph.D., Associate Professor, Dept. of Mathematics and Statistics, Queen's University, Kingston, ON, Canada (Felicia.Magpantay@queensu.ca)

## Introduction
The prediction functions are developed based on following works with modifications:

1. Typical bond percolation method original developed by M.E.J Newman.
  [Newman, M. E. J. (2002). Spread of epidemic disease on networks, Phys. Rev. E 66, 016128.](https://journals.aps.org/pre/pdf/10.1103/PhysRevE.66.016128)

2. Multiple type network percolation method by S. Bansal & L.A. Meyers.
  [Bansal, S. & Meyers, L. A. (2012). The impact of past epidemics on future disease dynamics, Journal of Theoretical Biology 309, 176–184.](https://www.bio.utexas.edu/research/meyers/_docs/publications/BansalJTB12.pdf)

3. Modified network percolation process with disease dynamics by J.C. Miller, A.C. Slim & E.M.Volz.
  [Miller, J. C., Slim, A. C. & Volz, E. M. (2012). Edge-based compartmental modelling for infectious disease spread, Journal of The Royal Society Interface 9, 890–906.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3306633/pdf/rsif20110403.pdf)

Also include different random network generator with specified degree sequence, based on following algorithms with modifications:
1. Original test algorithms by developers.

2. Sequential importance sampling algorithm by J. Blitzstein & P. Diaconis.
[Blitzstein, J. & Diaconis, P. (2011). A sequential importance sampling algorithm for generating random graphs with prescribed degrees, Internet mathematics 6, 489–522.](https://statweb.stanford.edu/~cgates/PERSI/papers/GraphAlgorithm.pdf)

3. Asymptotic sequential algorithm by M. Bayati, J.H. Kim & A. Saberi, with modification.
[Bayati, M., Kim, J. H. & Saberi, A. (2010). A sequential algorithm for generating random graphs, Algorithmica 58, 860–910.](https://arxiv.org/abs/cs/0702124)

The disease simulation function is developed based on Doob–Gillespie algorithm (Stochastic Simulation Algorithm, SSA) 


## Dependent R Packages and Environment:
The package is developed under R version 4.2.2 and the dependent packages are:

1. [pracma](https://cran.r-project.org/web/packages/pracma/index.html)

2. [gsl](https://cran.r-project.org/web/packages/gsl/index.html)

3. [deSolve](https://cran.r-project.org/web/packages/deSolve/)

# Looking for comments and collaboration! 
Developing, optimizing and testing are still in progress now!

Looking forward to comments, issues, bug reports, using opinions and discussions from all users. 
Please contact with primary developer by email 20sz11@queensu.ca

## Developing roadmap

### 1st phase: Basic theories and functions (2022 Sep - 2023 May) \[finished\]: 

\[Finished\] Basic theoretical structure, prediction functions, network generator and outbreak simulator and first round of function testing

### 2nd phase: Verification and network generator algorithm optimization (2023 Jun - 2023 Dec) \[Ongoing\]:

\[Finished\] Verify testing results with existing literature, testing on common degree distribution

\[Ongoing\] Enrich options of network generator for different cases

### 3rd phase: Data Structure and parallel computing optimization (2023 Nov - 2024 Feb) \[initiated\]

\[Ongoing\]Finding better data structure to enlarge network size capacity under resourse limitation

Optimize codes for parallel computing and calculation efficiency

### 4th phase: Publication, Documentation and Promotion (2023 Jun - 2024 Aug) \[Ongoing\]

\[Ongoing\] GitHub Repository creation

\[Ongoing\] User Manual and Documentation

CRAN package submission

\[Ongoing\] Package promotion within social media, academic conferences and communications

### 5th phase: Test and optimization with real application (2023 Sep - 2024 Jun) \[initiated\] 

Collaboration with Knowldege Management Team of KFL&A Public Health for application with real community data

Optimization based on feed back

## Related Paper
A related paper for systematic interpretation of theoretical background, methodology and modifications will be submitted for review soon. 

# Details

## Prediction Functions based on Percolation Process

### [PercolationProcess.R](PercolationProcess.R)
Including basic functions for data sorting & formatting, PGFs & modeling tool functions and all prediction functions for percolation methods.
Except the uncertain corrections in Multiple type network percolation method, the functions here are all tested and verified with original article results.

1. Data sorting & formatting functions: DegreeDF DDistPK
2. PGFs & Modeling tools functions: PGFG0 PGFD1G0 PGFD2G0 Kn PGFG1 PGFD1G1 ueqn
3. Typical bond percolation method prediction: Tconst_Neman Tconst_exp TypProc
4. Multiple type network percolation method prediction (with some uncertain corrections of the original method): P2u_polar P2s_polar Pl2_polar TypProc_Polar P2_leaky Q2_leaky Tmatrix PGFdf1_leaky PGFdf2_leaky s_leaky fsys_soln ES_leaky TypProc_leaky
5. Modified network percolation method prediction: ModProc_CM

## Random Network Generator with Specified Degree Sequence 
Currently, there is not a well recognized optimal/best algorithm to realize networks generally under all specified degree sequence.
There exist many different algorithms under several different theories and ideas.
Based on the core mathematical theories and idea, these methods can be generally classified into the following 4 family with each of its own pros and cons: Configuration model algorithms, Sequential algorithm, Switching-based algorithms and Markov chain algorithms.
There are multiple algorithms in each class optimized to fit different scenario and there are also algorithms using mixed ideas from multiple families.
So it depends on the researcher/user to consider the tradeoffs of pros and cons case by case and determine which algorithm is optimal for the given 

As for now, the package originally targets on providing better interpretation of the percolation process methods in epidemiology networks and trying to extend their model applications to a wider range of problems with different assumptions. 
So we decided to start with the family of Sequential algorithms, considering its following advantages:

1. Its interpretation is simpler and relatively close to the intuitive understanding for percolation process: both assuming "asymptotical independence" and considering a degree related random paring process when generating edges.
2. The assumptions for our epidemiology network models also fits better with its optimal condition: large sparse networks.
3. It provides a relatively high calculation efficiency among all families, especially for large sparse networks in general.
4. It guaranteed that every possible graph has a positive possibility to appear as the generation result.

However, the disadvantage of sequential algorithms is that it is harder to prove that its results is uniformly distributed among all possible random graphs, which is required by application and interpretation of the percolation method.
It might require extra sampling techniques on simulation results or extra conditions on original degree sequences 

In the package, we plan to provide 4 different sequential algorithms, considering the method illustration purpose and optimal performances under several possible application scenrario.
In each algorithms, provides with and Erdos-Gallai condition checking function modified to fit each algorithm.

These algorithms are temporary, we are still developing these generators while comparing their performance (and with other algorithms) under testing, to figure out if there is an optimal choice for our models, or if there is some systematical principle to choose algorithms under given condition.

### [RandomNetworkGenerator_Algorithm1.R](RandomNetworkGenerator_Algorithm1.R): Original algorithm developed by developers, inspired by Newman 
It is fastest in all scenario, but has random failure, lack of rigorous analysis of distribution result.
It is used to quickly generate relatively larger size networks for testing and verification, and also to illustrate Newman's idea for comparison.

### [RandomNetworkGenerator_Algorithm2.0.R](RandomNetworkGenerator_Algorithm2.0.R): Sequential importance algorithm with importance sampling technique, developed with some modification, based on J. Blitzstein & P. Diaconis
Its mainly for illustration the idea of Blitzstein and Diaconis who originally develop the first sequential algorithm for general degree sequences.
It come with an extra sampling procedure to adjust the simulation result to mathematically close to uniform distributed, while increasing some calculation cost.
An optimized version with several modification are provided by Algorithm 2.1

### [RandomNetworkGenerator_Algorithm2.1.R](RandomNetworkGenerator_Algorithm2.1.R) (Developing in Progress): Optimized version of Algorithm 2.0
With several modification on Algorithm 2.0, it provides better calculation efficiency than 2.0.
Currently this include the Durfee number optimization on each round of Erdős–Gallai condition checking and lower down the candidate list length within a single host vertex by reapply list from last round for each round of generation.
These modifications significantly increase the efficiency for large networks.
We are still trying to include the improvement inspired by Moseman in the next step, but it might require more mathematical evidence for equivalence in result.

From testing, it performs best when generating network not so "sparse"---for degree sequences with higher mean or with longer tail (like powerlaw distribution).
The corresponding version of importance sampling technique are also in developing now.

### [RandomNetworkGenerator_Algorithm3.R](RandomNetworkGenerator_Algorithm3.R): Mixture algorithm of 1 and 2
The algorithm is inspired by both Algorithm 1 and 2.1, which intuitively should be equivalent to Algorithm 2, but might need rigorous proof.
The idea is generate candidate list for each host vertex all at once and then check the EG condition, while guarantee success in each round of generation.
This significantly reduce the frequency of sampling and EG checking, however, leads to longer restarting duration if the sampling is failed.

From testing, it performs best for relatively smaller network with concentrated distribution like Mean Field and Poisson, and with higher mean degree.
The problem is the restarting time is much longer for large degree vertices at the beginning.
The sampling technique should be similar with Algorithm 2.1.

### RandomNetworkGenerator_Algorithm4.R (Developing in Progress, not uploaded yet): Asymptotic Uniform Algorithm based on M. Bayati, J.H. Kim & A. Saberi
Theoretically, this algorithm seems to be the "best" in general large enough networks.
It does not require sampling techniques for the simulation results, thus guarantee asymptotically uniform, but require a relatively strict condition on network size, based on maximal degree of the sequence.
So it might still be not ideal for smaller size networks or networks with long tail degree distribution, thus other algorithm are still useful in application.

This algorithm requires better understanding on random network theories, thus we are trying to development it soon after the submission of the theoretical paper.

### General Optimization of data structure and parallel computation (Ongoing, looking for collaboration)
The generation and epidemic simulation are built based on sparce matrix at the beginning.
While it is easier to understand and modify, the matrix data form in R might not be a good choice to store data, considering the limitation of computing resources, especially RAM, so the network size is highly limited.
It might also slows down the calculation speed when loops and iterations are involved.
We are currently working on, together with collaborating data scientists and programmers, is to optimize the data structure in these algorithms to enlarge the capacity of network size for same given RAM limitation and perhaps accelerate the calculation.
Another ongoing improvement will be optimizing these generators and the whole package for parallel computing in multiple servers and computing clusters, so simulation and verification involving large amount of simulations could be finished faster and on larger network size. 

## Epidemic Outbreak Simulator on Network
### [EpidemicDynamicSimulator.R](EpidemicDynamicSimulator.R)
Include the epidemic dynamic simulator function GilAlgo, which is developed based on the well-known widely used Doob–Gillespie algorithm (Stochastic Simulation Algorithm, SSA).
This simulator provides dynamic information of a randomly initiated disease on a given network, represented now by a sparse matrix, with the SIR disease parameters.
Current version is tested but might need to be optimized after the data structure optimization on network generators.

















