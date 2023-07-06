#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

//------------------------------------------------------------------------------
// Author: Irene Tubikanec
// Date: 2023-05-23
//
// Description: Strang splitting method for path simulation of the stochastic
//              multi-population JRNMM, proposed in the paper:
//
//              Network inference in a stochastic multi-population neural mass 
//              model, by S. Ditlevsen, M. Tamborrino and I. Tubikanec
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//Matrix-Vector Multiplication
//
// Input: matrix mat, vector vec
// Output: product of mat and vec
//------------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector mv_mult_JRNMM_(NumericMatrix mat, NumericVector vec)
{
  NumericVector ret(mat.nrow());
  double temp=0;

  for(int i = 0; i < mat.nrow(); i++)
  {
    for(int j = 0; j < vec.size(); j++)
    {
      temp = temp + mat(i,j) * vec[j];
    }
    ret[i]=temp;
    temp=0;
  }
  return ret;
};

//------------------------------------------------------------------------------
//Sigmoid function of the JRNMM
//
// Input:
// x             value at which the sigmoid function should be evaluated
// vmax, v0, r   model parameters of the JRNMM
// 
// Output:       value of sigmoid function evaluated in x
//------------------------------------------------------------------------------
// [[Rcpp::export]]
double sigmoid_JRNMM_Cpp_(double x, double vmax, double v0, double r)
{
  double ret=vmax/(1+exp(r*(v0-x)));
  return ret;
};

//------------------------------------------------------------------------------
//linear SDE: exact solution of 1st subequation of splitting procedure for the JRNMM
//
// Input:
// vec       current value X^[1](ti) of the solution of the 1st subequation
// dm        exponential matrix appearing in the exact solution of the 1st subequation
// xi        normally distributed random vector xi appearing in the exact solution of the 1st subequation
//
// Output:   next value X^[1](ti+1) of the solution of the 1st subequation
//------------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector linear_JRNMM_Cpp_(NumericVector vec, NumericMatrix dm, NumericVector xi)
{
  NumericVector ret=mv_mult_JRNMM_(dm,vec)+xi;
  return ret;
};

//------------------------------------------------------------------------------
//nonlinear ODE: exact solution of 2nd subequation of splitting procedure for the JRNMM
//
// Input:
// N        number of populations in the JRNMM
// vec      current value X^[2](ti) of the solution of the 2nd subequation
// h        step size for path simulation
// Theta    matrix containing continuous model parameters of the JRNMM
// Rho      matrix containing the coupling direction parameters of the JRNMM
// K        matrix containing the coupling strength parameters of the JRNMM
//
// Output:  next value X^[2](ti+1) of the solution of the 2nd subequation
//------------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector nonlinear_JRNMM_Cpp_(int N, NumericVector vec, double h, NumericMatrix Theta, NumericMatrix Rho, NumericMatrix K)
{
  NumericVector help(vec.size());
  double sum=0;
  int l=(6*N)/2;

  for(int k = 0; k < N; k++)
  {
    help(3*k+l)=Theta(k,0)*Theta(k,2)*sigmoid_JRNMM_Cpp_(vec(3*k+1)-vec(3*k+2),Theta(k,8),Theta(k,6),Theta(k,7));
    help(3*k+(l+2))=Theta(k,1)*Theta(k,3)*0.25*Theta(k,4)*sigmoid_JRNMM_Cpp_(0.25*Theta(k,4)*vec(3*k),Theta(k,8),Theta(k,6),Theta(k,7));

    sum=0;
    for(int j = 0; j < N; j++)
    {
      if(j!=k){
        sum=sum+Rho(j,k)*K(j,k)*vec(3*j);
      }
    }
    help(3*k+(l+1))=Theta(k,0)*Theta(k,2)*(Theta(k,5)+0.8*Theta(k,4)*sigmoid_JRNMM_Cpp_(Theta(k,4)*vec(3*k),Theta(k,8),Theta(k,6),Theta(k,7))+sum);

  }

  NumericVector ret=vec+h*help;
  return ret;
};

//------------------------------------------------------------------------------
//Strang splitting method for the JRNMM
//
// Input:
// N        number of populations in the JRNMM
// grid     time grid for path simulation
// h        step size for path simulation
// startv   starting value X0 for path simulation
// dm       exponential matrix appearing in the exact solution of the 1st subequation
// meanVec  zero vector used to generate the random vector xi appearing in the exact solution of the 1st subequation
// covMat   covariance matrix (after Cholesky decomposition) appearing in the exact solution of the 1st subequation
// Theta    matrix containing continuous model parameters of the JRNMM
// Rho      matrix containing the coupling direction parameters of the JRNMM
// K        matrix containing the coupling strength parameters of the JRNMM
// 
// Output:  path of the stochastic N-population JRNMM
//------------------------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix JRNMM_Splitting_Cpp_(int N_i, NumericVector grid_i, double h_i, NumericVector startv, NumericMatrix dm_i, NumericVector meanVec_i, NumericMatrix covMat_i, NumericMatrix Theta_i, NumericMatrix Rho_i, NumericMatrix K_i)
{
  int N=N_i;
  double h=h_i;
  NumericVector start=startv;
  NumericVector grid=grid_i;
  int iter=grid.size();
  int dim=6*N;

  NumericMatrix dm=dm_i;
  NumericVector meanVec=meanVec_i;
  NumericMatrix covMat=covMat_i;
  NumericMatrix Theta=Theta_i;
  NumericMatrix Rho=Rho_i;
  NumericMatrix K=K_i;

  NumericMatrix sol(dim,iter);
  sol(_, 0)=start;
  NumericVector newv=start;

  Function f_rmvn("rmvn");
  NumericMatrix xi=f_rmvn(iter,meanVec,covMat,Named("isChol", true));
  NumericVector xi_rel;

  for(int i=1;i<iter;i++)
  {
    xi_rel=xi(i,_);
    newv=nonlinear_JRNMM_Cpp_(N,newv,h/2,Theta,Rho,K);
    newv=linear_JRNMM_Cpp_(newv,dm,xi_rel);
    newv=nonlinear_JRNMM_Cpp_(N,newv,h/2,Theta,Rho,K);
    sol(_,i)=newv;
  }

  NumericMatrix ret=sol;

  return sol;
};


