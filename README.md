This folder contains all files necessary to reproduce the numerical results in the manuscript, 
"Bayesian estimation of stochastic gene expression using multifidelity models" 
Huy D. Vo, Zachary R. Fox, Ania Baetica, Brian E. Munsky.

Before running the examples, please run the script 'setup.m' provided in the 
main folder to add the necessary functions and objects to MATLAB's search path.

The following third-party softwares are prequisites for generating some of the figures and tables in 
the manuscript:
- multiESS (https://github.com/lacerbi/multiESS/blob/master/multiESS.m).
- MCMCSTAT (https://github.com/mjlaine/mcmcstat).

Information about the subfolders:

- DataSynthesizer contains functions to simulate smFISH experimental observations at specific time points using stochastic simulation.

- FSP contains the data structure and methods to store and evaluate the full parameter-dependent FSP model.

- ROM contains the data structure and methods to store and evaluate the low-fidelity models of the full FSP.

- MCMC contains the MCMC routines for sampling from the posterior distribution. 

- GA contains the custom mutation used for the genetic algorithm in two examples in the paper.

- MatrixExponential contains a modified version of the Expokit (https://www.maths.uq.edu.au/expokit/) to output the solutions of a linear dynamical system at multiple time points.

- Examples contains three subfolders with codes to reproduce the examples in the manuscript and a 'Toy_problem' subfolder. Note that the 'spatial_gene_expression' example will output very large matrices (up to 37 GB). Please make sure your machine has sufficient memory before running the codes in the 'spatial_gene_expression' subfolder.

Questions should be addressed to Brian Munsky (munsky@colostate.edu) and Huy Vo (huydvo@colostate.edu).
