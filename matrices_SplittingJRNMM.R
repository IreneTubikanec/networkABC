
#-----------------------------------------------------------------------------------
# Author: Irene Tubikanec
# Date:  2023-05-23
#-----------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
# Exponential and covariance matrices for splitting-simulation of the JRNMM
#-----------------------------------------------------------------------------------

#Input:
#N        number of populations in the JRNMM
#Gamma    diagonal matrix containing the model parameters a and b of the JRNMM
#h        step size for path simulation

#Output:  exponential matrix appearing in the exact solution of the first subequation in the splitting procedure

exponential_matrix_JRNMM<-function(N,Gamma,h)
{
  #determine matrix
  dim<-3*N
  exp_Gamma<-expm(-Gamma*h)
  diag_mat<-diag(dim)
  
  mat1<-exp_Gamma%*%(diag_mat+Gamma*h)
  mat2<-exp_Gamma*h
  mat3<--Gamma%*%Gamma%*%exp_Gamma*h
  mat4<-exp_Gamma%*%(diag_mat-Gamma*h)
  
  mat_top<-cbind(mat1,mat2)
  mat_bottom<-cbind(mat3,mat4)
  
  ret<-rbind(mat_top,mat_bottom)
  
  #return matrix
  return(ret)
}

#-----------------------------------------------------------------------------------

#Input:
#N        number of populations in the JRNMM
#Gamma    diagonal matrix containing the model parameters a and b of the JRNMM
#Sigma    diagonal matrix containing the noise parameters epsilon and sigma of the JRNMM
#h        step size used for path simulation

#Output:  covariance matrix appearing in the exact solution of the first subequation in the splitting procedure

covariance_matrix_JRNMM<-function(N,Gamma,Sigma,h)
{
  #determine matrix
  dim<-3*N
  exp_Gamma<-expm(-Gamma*h)
  diag_mat<-diag(dim)
  
  v1<-exp_Gamma%*%(diag_mat+Gamma*h)
  v2<-exp_Gamma*h
  v3<--Gamma%*%Gamma%*%exp_Gamma*h
  v4<-exp_Gamma%*%(diag_mat-Gamma*h)
  
  inv_Gamma<-diag(diag(Gamma)^-1,dim,dim)
  squared_Sigma<-Sigma%*%Sigma
  
  mat1<-(1.0/4.0)*(inv_Gamma%*%inv_Gamma%*%inv_Gamma)%*%squared_Sigma%*%(diag_mat+v2%*%v3-(v1%*%v1))
  mat2<-(1.0/2.0)*(squared_Sigma)%*%(v2%*%v2)
  mat3<-mat2
  mat4<-(1.0/4.0)*(inv_Gamma)%*%squared_Sigma%*%(diag_mat+v2%*%v3-(v4%*%v4))
  
  mat_top<-cbind(mat1,mat2)
  mat_bottom<-cbind(mat3,mat4)
  
  ret<-rbind(mat_top,mat_bottom)
  
  #return matrix
  return(ret)
}

