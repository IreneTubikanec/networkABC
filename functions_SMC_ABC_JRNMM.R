
#-----------------------------------------------------------------------------------
# Author: Irene Tubikanec
# Date:  2025-03-28
#-----------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
# Functions required for nSMC-ABC inference in the stochastic multi-population JRNMM
#-----------------------------------------------------------------------------------

#-------------------------------------------------------------------------
# Multivariate normal perturbation kernel 

#Input:
#x          vector at which the multivariate normal density should be evaluated
#mu_vec     mean vector of multivariate normal distribution
#sigma_mat  covariance matrix of multivariate normal distribution

#Output:    value of multivariate normal density evaluated in x
#-------------------------------------------------------------------------

continuous_Kernel<-function(x,mu_vec,sigma_mat)
{
  kernelval<-dmvn(x,mu_vec,sigma_mat)
  return(kernelval)
}

#-------------------------------------------------------------------------
# Perturb a (continuous) sample

#Input: 
#theta_c_sampled  sampled continuous parameter vector 
#sigma_kernel     covariance matrix of multivariate perturbation kernel
#Pr_cont          Prior distribution for continuous model parameters

#Output: perturbed continuous sample, which still lies within the prior region
#-------------------------------------------------------------------------

perturb_continuous<-function(theta_c_sampled,sigma_kernel,Pr_cont){ 
  
  c<-length(theta_c_sampled)
  
  help_while<-TRUE
  help_while_vec<-rep(TRUE,N)
  
  while(help_while){
    #perturb 
    theta_c_perturbed<-rmvn(1,mu=theta_c_sampled,sigma=sigma_kernel) 
    #check if perturbed values lie inside the prior region
    for(i in 1:c){
      help_while_vec[i]<-(theta_c_perturbed[i]<=Pr_cont[i,1]||theta_c_perturbed[i]>=Pr_cont[i,2])
    }
    help_while<-sum(help_while_vec)
  }  
  
  #return perturbed sample
  return(theta_c_perturbed)
}

#-------------------------------------------------------------------------
# Generate (discrete) sample according to success-prob. from previous population

#Input: 
#N          number of populations in the JRNMM
#hat_p_vec  Bernoulli success probabilities estimated from previous population

#Output:    sample from the discrete Bernoulli parameters (using hat_p_vec) stored in a matrix
#-------------------------------------------------------------------------

sample_discrete<-function(N,hat_p_vec){
  
  #sample according to hat_p_vec
  Rho<-matrix(Inf,nrow=N,ncol=N)
  counter<-0
  for(j in 1:N){
    for(k in 1:N){
      if(k!=j){
        counter<-counter+1
        Rho[j,k]<-sample(c(1,0),1,prob=c(hat_p_vec[counter],1-hat_p_vec[counter]))
      }
    }
  }
  ret<-Rho
  
  #return sample
  return(ret)
}

#-------------------------------------------------------------------------
# Determination of parameters needed to determine ABC summaries of a dataset

#Input: 
#X       observed dataset/path of the JRNMM
#T       time horizon for path simulation
#h       step size for path simulation

#Output: vector of required summary parameters
#-------------------------------------------------------------------------

ABC_summaries_parameters<-function(X,T,h){
  
  #parameters for density
  startSupp<--40
  endSupp<-40
  Lsupport<-1001
  stepD<-(endSupp-startSupp)/(Lsupport-1)
  
  #parameters for spectral densities 
  span_val<-5*T
  lag_val<-100
  specDens<-spectrum(X[1,],log="no",span=span_val,plot=FALSE)
  spx<-specDens$freq 
  stepP<-diff(spx)[1]
  Lspec<-length(spx)
  
  #parameters for cross-correlation function
  Lccf<-201
  stepCcf<-h
  
  #compute the vector of summary parameters and return it
  summaries_parameters<-c(startSupp,endSupp,Lsupport,stepD,span_val,Lspec,stepP,lag_val,Lccf,stepCcf)
  return(summaries_parameters)
}

#-------------------------------------------------------------------------
# Computation of ABC summaries of a dataset

#Input:
#X       observed dataset/path of the JRNMM
#T       time horizon for path simulation
#h       step size for path simulation

#Output: summaries (densities, spectral densities, ccfs) of a dataset X
#-------------------------------------------------------------------------

ABC_summaries<-function(in_X,in_T,in_hsim,in_summaries_parameters)
{
  #Input general
  X<-in_X #data
  T<-in_T #time horizon
  hsim<-in_hsim #step size
  
  N<-length(X[,1]) #number of populations
  
  #Input summary parameters
  summaries_parameters<-in_summaries_parameters
  startSupp<-summaries_parameters[1]
  endSupp<-summaries_parameters[2]
  Lsupport<-summaries_parameters[3] 
  span_val<-summaries_parameters[5]
  Lspec<-summaries_parameters[6]
  lag_val<-summaries_parameters[8]
  Lccf<-summaries_parameters[9]
  
  #-----------------------------------------------------
  # Densties, spectral densities and ccf's of X 
  #-----------------------------------------------------
  
  dens_mat<-matrix(0,nrow=N,ncol=Lsupport)
  spec_mat<-matrix(0,nrow=N,ncol=Lspec)
  ccf_mat<-matrix(0,nrow=(N*(N-1)),ncol=Lccf)
  
  count_ccf<-0
  for(j in 1:(N-1)){
    for(k in (j+1):N){
      #ccf
      Xj<-ts(X[j,], frequency = 1/hsim)
      Xk<-ts(X[k,], frequency = 1/hsim)
      count_ccf<-count_ccf+1
      ccfXjk<-ccf(Xj,Xk,lag=lag_val,plot=FALSE)
      ccf_mat[count_ccf,]<-ccfXjk$acf
      count_ccf<-count_ccf+1
      ccfXkj<-ccf(Xk,Xj,lag=lag_val,plot=FALSE)
      ccf_mat[count_ccf,]<-ccfXkj$acf
    }
    #spec
    spec_mat[j,]<-spectrum(X[j,],log="no",span=span_val,plot=FALSE)$spec #specXjk$spec1
    #dens
    dens_mat[j,]<-density(X[j,],n=Lsupport,from=startSupp,to=endSupp)$y
  }
  spec_mat[N,]<-spectrum(X[N,],log="no",span=span_val,plot=FALSE)$spec #specXjk$spec2
  dens_mat[N,]<-density(X[N,],n=Lsupport,from=startSupp,to=endSupp)$y
  
  ret<-list()
  ret[[1]]<-dens_mat
  ret[[2]]<-spec_mat
  ret[[3]]<-ccf_mat
  
  return(ret)
}


#-------------------------------------------------------------------------
# Computation of weights (for distance measure among two sets of summaries) from observed data

#Input:
#summaries             summaries of observed dataset
#summaries_parameters  corresponding summary parameters
#N                     number of populations in the JRNMM

#Output:               vector of summary weights
#-------------------------------------------------------------------------

ABC_summaries_weights<-function(summaries,summaries_parameters,N){
  
  #access specific summary parameters
  stepP<-summaries_parameters[7]
  stepCcf<-summaries_parameters[10]
  
  #access different summaries
  dens_mat<-summaries[[1]] #dens
  spec_mat<-summaries[[2]] #spec
  ccf_mat<-summaries[[3]] #ccf
  
  specA<-rep(0,N)
  for(j in 1:N){
    specA[j]<-sum(stepP*abs(spec_mat[j,]))
  }
  mean_specA<-mean(specA)
  
  ccfA<-rep(0,N)
  for(j in 1:(N*(N-1))){
    ccfA[j]<-sum(stepCcf*abs(ccf_mat[j,]))
  }
  mean_ccfA<-mean(ccfA)
  
  #weight density
  we_dens<-mean_specA
  we_dens
  
  #weight ccf
  we_ccf<-mean_specA/mean_ccfA
  we_ccf
  
  #determine vector of summaries and return it
  summaries_weights<-c(we_dens,we_ccf)
  return(summaries_weights)
}

#-------------------------------------------------------------------------
# Computation of the ABC distance

#Input:
#summaries_X            summaries of observed dataset X
#summaries_Y            summaries of synthetic dataset Y
#summaries_parameters   corresponding summary parameters
#summaries_weights      weights for distance computation derived from X

#Output: distance value among two sets of summaries
#-------------------------------------------------------------------------

ABC_distance<-function(summaries_X,summaries_Y,summaries_parameters,summaries_weights)
{
  
  #access the different summaries
  #densities
  dens_matX<-summaries_X[[1]]
  dens_matY<-summaries_Y[[1]]
  #spectral densities
  spec_matX<-summaries_X[[2]]
  spec_matY<-summaries_Y[[2]]
  #ccf's
  ccf_matX<-summaries_X[[3]]
  ccf_matY<-summaries_Y[[3]]
  
  #access summary weights
  we_dens<-summaries_weights[1]
  we_ccf<-summaries_weights[2]
  
  #access specific summary parameters
  stepD<-summaries_parameters[4]
  stepP<-summaries_parameters[7]
  stepCcf<-summaries_parameters[10]

  #Compute distance
  dist_dens<-0
  dist_spec<-0
  dist_ccf<-0
  
  for(j in 1:(N*(N-1))){
    if(j<=N){
    dist_dens<-dist_dens+sum(abs(stepD*(dens_matX[j,]-dens_matY[j,])))
    dist_spec<-dist_spec+sum(abs(stepP*(spec_matX[j,]-spec_matY[j,])))
    }
    dist_ccf<-dist_ccf+sum(abs(stepCcf*(ccf_matX[j,]-ccf_matY[j,])))
  }
  
  ret<-dist_spec/N+(dist_ccf/(N*(N-1)))*we_ccf+(dist_dens/N)*we_dens
  
  #return distance
  return(ret)
}

#-------------------------------------------------------------------------
# ABC: Pilot simulation
# Procedure happening within the for-loop of the reference table ABC algorithm

#Input:
#N       number of populations in the JRNMM
#T       time horizon for path simulation
#h       step size for path simulation
#grid    time grid for path simulation
#startv  starting value X0 for path simulation
#Theta   continuous model parameters of the JRNMM
#Rho     discrete coupling direction parameters of the JRNMM   
#K       continuous coupling strength parameters of the JRNMM
#dm      exponential matrix used for path simulation
#meanVec zero vector used for path simulation
#cm      covariance matrix (after Cholesky decomposition) used for path simulation
#summaries_X            summaries of observed dataset X
#summaries_parameters   corresponding summary parameters
#summaries_weights      weights for distance computation derived from X
#Pr_cont                prior distribution for continuous model parameters

#Output: Distance value between the summaries of the observed dataset and a simulated synthetic dataset
#-------------------------------------------------------------------------

ABC_pilot<-function(N,T,h,grid,startv,Theta,Rho,K,dm,meanVec,cm,summaries_X,summaries_parameters,summaries_weights,Pr_cont){
  
  #-----------------------------------------------------------
  #sample theta according to the prior
  
  #update Theta - all A's are unknown
  for(j in 1:N){
    Aj<-runif(1,Pr_cont[j,1],Pr_cont[j,2])
    Theta[j,1]<-Aj
  }
  
  #update Rho - all Rho's are unknown
  Rho_vec<-sample(c(0,1),N^2,prob=c(1/2,1/2),replace=TRUE) #Bernoulli prior with p=1/2
  Rho<-matrix(Rho_vec,byrow=TRUE,nrow=N)
  diag(Rho)<-Inf
  
  #Update K - all K's (L and c) are unknown
  L<-runif(1,Pr_cont[N+1,1],Pr_cont[N+1,2])
  c<-runif(1,Pr_cont[N+2,1],Pr_cont[N+2,2])
  for(j in 1:(N-1)){
    for(k in (j+1):N){
      K_value<-c^(abs(k-j)-1)*L
      K[j,k]<-K_value
      K[k,j]<-K_value
    }
  }
  
  #-----------------------------------------------------------
  #conditioned on theta, simulate a synthetic dataset Y
  
  Y<-matrix(0,nrow=N,ncol=length(grid))
  sol<-JRNMM_Splitting_Cpp(N,grid,h,startv,dm,meanVec,cm,Theta,Rho,K)
  for(j in 1:N){
    sol_l<-3*(j-1)+2
    sol_r<-3*(j-1)+3
    Y[j,]<-sol[sol_l,]-sol[sol_r,]
  }
  
  #-----------------------------------------------------------
  #compute the summaries of the synthetic dataset Y  

  summaries_Y<-ABC_summaries(Y,T,h,summaries_parameters)

  #-----------------------------------------------------------
  #calculate the distance to the observed reference dataset X
  
  Dist<-ABC_distance(summaries_X,summaries_Y,summaries_parameters,summaries_weights)

  #-----------------------------------------------------------
  #Return the calculated distance
  
  dim_arr<-1
  ret<-array(0,dim=c(1,dim_arr,1))
  ret[,,1]<-c(Dist)  
  
  return(ret)
}

#-------------------------------------------------------------------------
#SMC-ABC: Iteration 1 (r=1)
# Procedure happening within the first for-loop of the SMC-ABC algorithm

#Input:
#N        number of populations in the JRNMM
#T        time horizon for path simulation
#h        step size for path simulation
#grid     time grid for path simulation
#startv   starting value X0 for path simulation
#Theta    continuous model parameters of the JRNMM
#Rho      discrete coupling direction parameters of the JRNMM   
#K        continuous coupling strength parameters of the JRNMM
#dm       exponential matrix used for path simulation
#meanVec  zero vector used for path simulation
#cm       covariance matrix (after Cholesky decomposition) used for path simulation
#summaries_X            summaries of observed dataset X
#summaries_parameters   corresponding summary parameters
#summaries_weights      weights for distance computation derived from X
#Pr_cont                prior distribution for continuous model parameters
#delta_1  initial threshold for distance calculation obtained from ABC pilot simulation

#Output: Distance value between the summaries of the observed dataset and a simulated synthetic dataset
#        and sampled (kept) parameter values 
#-------------------------------------------------------------------------

ABC_SMC_r1<-function(N,T,h,grid,startv,Theta,Rho,K,dm,meanVec,cm,summaries_X,summaries_parameters,summaries_weights,Pr_cont,delta_1){

  #-----------------------------------------------------------
  #set initial Dist>=delta_1 to enter while-loop
  Dist<-delta_1+1.0 

  #while-loop
  while(Dist>=delta_1){
    
    #-----------------------------------------------------------
    #sample theta according to the prior
    
    #Update Theta (all A's are unknown)
    for(j in 1:N){
      Aj<-runif(1,Pr_cont[j,1],Pr_cont[j,2])
      Theta[j,1]<-Aj
    }
    
    #Update Rho (all rho's are unknown)
    Rho_vec<-sample(c(0,1),N^2,prob=c(1/2,1/2),replace=TRUE) #Bernoulli prior with p=1/2
    Rho<-matrix(Rho_vec,byrow=TRUE,nrow=N)
    diag(Rho)<-Inf
    
    #Update K - all K's (L and c) are unknown
    L<-runif(1,Pr_cont[N+1,1],Pr_cont[N+1,2])
    c<-runif(1,Pr_cont[N+2,1],Pr_cont[N+2,2])
    for(j in 1:(N-1)){
      for(k in (j+1):N){
        K_value<-c^(abs(k-j)-1)*L
        K[j,k]<-K_value
        K[k,j]<-K_value
      }
    }
    
    #-----------------------------------------------------------
    #conditioned on theta, simulate a synthetic dataset Y
    
    Y<-matrix(0,nrow=N,ncol=length(grid))
    sol<-JRNMM_Splitting_Cpp(N,grid,h,startv,dm,meanVec,cm,Theta,Rho,K)
    for(j in 1:N){
      sol_l<-3*(j-1)+2
      sol_r<-3*(j-1)+3
      Y[j,]<-sol[sol_l,]-sol[sol_r,]
    }
    
    #-----------------------------------------------------------
    #compute the summaries of the synthetic dataset Y  
    
    summaries_Y<-ABC_summaries(Y,T,h,summaries_parameters)
    
    #-----------------------------------------------------------
    #calculate the distance to the observed reference dataset X
    
    Dist<-ABC_distance(summaries_X,summaries_Y,summaries_parameters,summaries_weights)
    
  } #end while-loop
  
  #-----------------------------------------------------------
  #Return the calculated distance and sampled values
  
  dim_arr<-N*(N-1)+N+3
  
  #A's
  A_vector<-Theta[,1]
  
  #Rho
  Rho_vec<-as.vector(t(Rho))
  Rho_vec2<-Rho_vec[-which(Rho_vec==Inf)]
  
  ret<-array(0,dim=c(1,dim_arr,1))
  ret[,,1]<-c(Dist,A_vector,L,c,Rho_vec2)  
  
  return(ret)
}

#-------------------------------------------------------------------------
#SMC-ABC: Iterations r>1
# Procedure happening within the second for-loop of the SMC-ABC algorithm

#Input:
#N                number of populations in the JRNMM
#T, h, grid       time horizon, step size and time grid for path simulation
#startv           starting value X0 for path simulation
#Theta, Rho, K    model parameters of the JRNMM
#dm, meanVec, cm  exponential matrix, zero vector, and covariance matrix (after Cholesky decomposition) used for path simulation
#summaries_X                          summaries of observed dataset X
#summaries_parameters                 corresponding summary parameters
#summaries_weights                    weights for distance computation derived from X
#kept_A, kept_L, kept_c               kept samples of previous iteration
#kept_weights_c                       corresponding weights of previous iteration
#sigma_kernel, hat_p_vec, stay_prop   estimates (derived from previous population) used for perturbation
#Pr_cont          prior distribution for continuous model parameters
#delta_r          current threshold for distance calculation obtained from previous population

#Output: Distance value between the summaries of the observed dataset and a simulated synthetic dataset,
#        sampled (kept) parameter values, and a counter for the acceptance rate of particles 
#-------------------------------------------------------------------------

ABC_SMC<-function(N,T,h,grid,startv,Theta,Rho,K,dm,meanVec,cm,summaries_X,summaries_parameters,summaries_weights,kept_A,kept_L,kept_c,kept_weights_c,sigma_kernel,hat_p_vec,stay_prob,Pr_cont,delta_r){

  #-----------------------------------------------------------
  #set initial Dist>=delta_r to enter while-loop
  Dist<-delta_r+0.1 
  
  #-----------------------------------------------------------
  #Counter for acceptance rate
  count_acc<-0
  
  #-----------------------------------------------------------
  #while-loop
  while(Dist>=delta_r){
    
    #-----------------------------------------------------------------------
    #sample theta from the weighted set and perturb it
    
    #sample: continuous
    index_c<-sample(length(kept_weights_c),1,prob=kept_weights_c) 
    As<-kept_A[index_c,] 
    L<-kept_L[index_c] 
    c<-kept_c[index_c] 
    theta_c_sampled<-c(As,L,c)
    
    #perturb: continuous
    theta_c_perturbed<-perturb_continuous(theta_c_sampled,sigma_kernel,Pr_cont)
    As<-theta_c_perturbed[1:N]
    L<-theta_c_perturbed[N+1]
    c<-theta_c_perturbed[N+2]
    
    #sample: discrete
    Rho<-sample_discrete(N,hat_p_vec)
    
    #perturb: discrete
    counter<-0
    for(j in 1:N){
      for(k in 1:N){
        if(k!=j){
          counter<-counter+1
          #0: not perturb (stay), 1: perturb (move=1-stay)
          perturb<-sample(c(0,1),1,prob=c(stay_prob,1-stay_prob)) 
          if(perturb==TRUE){
            Rho[j,k]<-(Rho[j,k]+1)%%2
          }
        }
      }
    }
    
    #Update Theta
    Theta[,1]<-As
    
    #Update K
    for(j in 1:(N-1)){
      for(k in (j+1):N){
        K_value<-c^(abs(k-j)-1)*L
        K[j,k]<-K_value
        K[k,j]<-K_value
      }
    }
    
    #-----------------------------------------------------------------------
    #conditioned on theta, simulate a synthetic dataset Y
    
    Y<-matrix(0,nrow=N,ncol=length(grid))
    sol<-JRNMM_Splitting_Cpp(N,grid,h,startv,dm,meanVec,cm,Theta,Rho,K)
    for(j in 1:N){
      sol_l<-3*(j-1)+2
      sol_r<-3*(j-1)+3
      Y[j,]<-sol[sol_l,]-sol[sol_r,]
    }
    
    #-----------------------------------------------------------
    #compute the summaries of the synthetic dataset Y  
    
    summaries_Y<-ABC_summaries(Y,T,h,summaries_parameters)
    
    #-----------------------------------------------------------------------
    #calculate the distance to the observed reference dataset X
    
    Dist<-ABC_distance(summaries_X,summaries_Y,summaries_parameters,summaries_weights)    
  
    #-----------------------------------------------------------------------
    #increase counter for acceptance rate
    count_acc<-count_acc+1
    
  } #end while-loop
  
  #-----------------------------------------------------------
  #Return the calculated distance, the sampled values and the counter for the acceptance rate
  dim_arr<-N*(N-1)+N+4
  
  Rho_vec<-as.vector(t(Rho))
  Rho_vec2<-Rho_vec[-which(Rho_vec==Inf)]
  
  ret<-array(0,dim=c(1,dim_arr,1))
  ret[,,1]<-c(Dist,As,L,c,Rho_vec2,count_acc)  
  
  return(ret)
}









