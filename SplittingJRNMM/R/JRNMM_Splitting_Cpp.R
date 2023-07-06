#'@rdname JRNMM_Splitting_Cpp
#'@title JRNMM_Splitting_Cpp
#'@description Strang splitting method for path simulation of the 
#' stochastic multi-population JRNMM
#'@return path of the stochastic N-population JRNMM
#'@export
JRNMM_Splitting_Cpp <- function(N, grid, h, startv, dm, meanVec, covMat, Theta, Rho, K){
  return(JRNMM_Splitting_Cpp_(N, grid, h, startv, dm, meanVec, covMat, Theta, Rho, K))
}
