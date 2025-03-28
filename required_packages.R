
#-----------------------------------------------------------------------------------
# Author: Irene Tubikanec
# Date:  2025-03-28
#-----------------------------------------------------------------------------------
# Required packages for nSMC-ABC inference in the JRNMM
#-----------------------------------------------------------------------------------

#Rcpp
library(Rcpp)
library(RcppNumerical)
library(devtools)
find_rtools(T)

#Parallelization
library(foreach)
library(doSNOW)
library(doParallel)

#Simulation of a path of the JRNMM (Cpp-code in the SplittingJRNMM-package)
library(SplittingJRNMM)

#Further packages
library(rexpokit) #matrix exponentiation
library(mvnfast) #multivariate normal distribution

