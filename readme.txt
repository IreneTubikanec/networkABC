
Author: Irene Tubikanec
Date:   2025-03-28

Description: This code provides an implementation of the nSMC-ABC method for network inference and parameter estimation in the stochastic multi-population JRNMM, proposed in the paper:
         
Network inference via approximate Bayesian computation. Illustration on a stochastic multi-population neural mass model, by S. Ditlevsen, M. Tamborrino and I. Tubikanec.

In particular, it reproduces the estimation results of Setting 2 "Partially connected network" considered in the paper (and can be easily adapted to other settings by storing a different set of reference data into the folder Reference_Data).

How does the code work?

1. Install relevant packages (see the R-file required_packages)

2. Install the package SplittingJRNMM (splitting-simulation of paths of the JRNMM using C++ Code) via:

install.packages("SplittingJRNMM_1.0.tar.gz")

3. Run the R-file main_SMC_ABC_JRNMM (nSMC-ABC inference for the JRNMM)
After termination of the algorithm, the kept ABC posterior samples are stored into the folder ABC_Results.

4. Run the R-file plot_results (visualization of estimation results)
It will create figures reporting the marginal ABC posterior densities/histograms of the continuuos model parameters and discrete (coupling direction) parameters, respectively. The figures also report the underlying true parameters values. The figure for the continuous parameters also shows the respective uniform prior distributions and the weighted ABC posterior means. 

Code description:
For a detailed description of the code, please consider the respective files

Licence information:
Please consider the txt-files LICENCE and COPYING.GPL.v3.0




