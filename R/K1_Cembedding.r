#' Fit mixed model for K=1 (one mobility class)
#' 
#' @param timeall numeric vector containing the time points at which measurements are available
#' @param concall numeric vector containing the measurements at time points specified in argument time
#' @param nobs numeric vector indicating how many measurements have been done for each nucleus
#' @param ttotal total number of time points at which measurements are available (for all nuclei)
#' @param J number of cell nuclei
#' @param sv numeric vector containing the parameters of the IG distributions of the random effects
#' @param size number of iterations
#' @param burnin length of burnin
#' @param thin thinning parameter (e.g., if thin=3, only every third iteration is used)
#' @param tune at which iterations tuning should be done 
#' @param nobs_cum numeric vector containing the cumulated number of time points at which measurements are available for all cells, starting with 0
#' @param sigma2b1 proposal variance for the MH step for b1off (may be different after the fit because of tuning during the fit)
#' @param sigma2beta1 proposal variance for the MH step for beta1joff (may be different after the fit because of tuning during the fit)
#' @param a0start starting value for a0
#' @param a1start starting value for a1
#' @param b1start starting value for log(b1off)
#' 
#' @return timeall, concall, nobs, ttotal, J, sv, size, burnin, thin, tune, nobs_cum, sigma2b1, sigma2beta1, a0start, a1start, b1start: see description above
#' @return tau2_alpha0, tau2_alpha1, tau2_beta1: variances of the random effects
#' @return sigma2: estimated variance of the error
#' @return a0, a1, b1: numeric vectors containing the saved results for the fixed effects a0, a1, log(b1off)
#' @return alpha0, alpha1, beta1: numeric vectors containing the saved results for the random effects
#' @return acc_b1, acc_beta1: number of accepted proposals in MH step for parameters b1off and beta1joff
#' 
#' @export fitk1
#' 
fitk1 <- function(timeall,concall,nobs,ttotal,J,sv,size,burnin,thin,tune,nobs_cum,sigma2b1,sigma2beta1,a0start,a1start,b1start){
  k1mhgibbs <- .C("k1_mhgibbs",
      as.double(timeall),      
      as.double(concall),
      as.integer(nobs),
      as.integer(ttotal),
      as.integer(J),
      as.double(sv), 
      as.integer(size),
      as.integer(burnin),
      as.integer(thin),
      as.integer(tune),
      as.integer(nobs_cum),
      as.double(sigma2b1),
      as.double(sigma2beta1),  
      as.double(a0start),
      as.double(a1start),
      as.double(b1start),
                  
      # create empty vectors for saving the results
      as.double(rep(0,((size-burnin)/thin))),   # tau2_alpha0          
      as.double(rep(0,((size-burnin)/thin))),   # tau2_alpha1
      as.double(rep(0,((size-burnin)/thin))),   # tau2_beta1
      as.double(rep(0,((size-burnin)/thin))),   # sigma2
                  
      as.double(rep(0,((size-burnin)/thin))),   # a0
      as.double(rep(0,((size-burnin)/thin))),   # a1
      as.double(rep(0,((size-burnin)/thin))),   # b1
      
      as.double(rep(0,((size-burnin)/thin)*J)), # alpha0
      as.double(rep(0,((size-burnin)/thin)*J)), # alpha1
      as.double(rep(0,((size-burnin)/thin)*J)), # beta1
      
      as.integer(0), # acc_b1
      as.integer(0), # acc_beta1
      
      PACKAGE="frapmm")

      return(list(
      "timeall"=k1mhgibbs[[1]],
      "concall"=k1mhgibbs[[2]],
      "nobs"=k1mhgibbs[[3]],
      "ttotal"=k1mhgibbs[[4]],
      "J"=k1mhgibbs[[5]],
      "sv"=k1mhgibbs[[6]],
      "size"=k1mhgibbs[[7]],
      "burnin"=k1mhgibbs[[8]],
      "thin"=k1mhgibbs[[9]],
      "tune"=k1mhgibbs[[10]],
      "nobs_cum"=k1mhgibbs[[11]],
      "sigma2b1"=k1mhgibbs[[12]],
      "sigma2beta1"=k1mhgibbs[[13]],
      "a0start"=k1mhgibbs[[14]],
      "a1start"=k1mhgibbs[[15]],
      "b1start"=k1mhgibbs[[16]],
      "tau2_alpha0"=k1mhgibbs[[17]],
      "tau2_alpha1"=k1mhgibbs[[18]],
      "tau2_beta1"=k1mhgibbs[[19]],
      "sigma2"=k1mhgibbs[[20]],
      "a0"=k1mhgibbs[[21]],
      "a1"=k1mhgibbs[[22]],
      "b1"=k1mhgibbs[[23]],
      "alpha0"=k1mhgibbs[[24]],
      "alpha1"=k1mhgibbs[[25]],
      "beta1"=k1mhgibbs[[26]],
      "acc_b1"=k1mhgibbs[[27]],
      "acc_beta1"=k1mhgibbs[[28]]
      ))
}