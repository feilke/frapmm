#' Fit mixed model for K=2 (two mobility classes)
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
#' @param sigma2b2 proposal variance for the MH step for b2off (may be different after the fit because of tuning during the fit)
#' @param sigma2beta2 proposal variance for the MH step for beta2joff (may be different after the fit because of tuning during the fit)
#' @param a0start starting value for a0
#' @param a1start starting value for a1
#' @param b1start starting value for log(b1off)
#' @param a2start starting value for a2
#' @param b2start starting value for log(b2off)
#' 
#' @return timeall, concall, nobs, ttotal, J, sv, size, burnin, thin, tune, nobs_cum, sigma2b1, sigma2beta1, sigma2b2, sigma2beta2, a0start, a1start, b1start, a2start, b2start: see description above
#' @return tau2_alpha0, tau2_alpha1, tau2_beta1, tau2_alpha2, tau2_beta2: variances of the random effects
#' @return sigma2: estimated variance of the error
#' @return a0, a1, b1, a2, b2: numeric vectors containing the saved results for the fixed effects a0, a1, log(b1off)
#' @return alpha0, alpha1, beta1, alpha2, beta2: numeric vectors containing the saved results for the random effects
#' @return acc_b1, acc_beta1, acc_b2, acc_beta2: number of accepted proposals in MH step for parameters b1off, b2off, beta1joff and beta2joff
#' 
#' @export fitk2
#' 
fitk2 <- function(timeall,concall,nobs,ttotal,J,sv,size,burnin,thin,tune,nobs_cum,sigma2b1,sigma2beta1,sigma2b2,sigma2beta2,a0start,a1start,b1start,a2start,b2start){
    k2mhgibbs <- .C("k2_mhgibbs",
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
      as.double(sigma2b2),
      as.double(sigma2beta2),  
      as.double(a0start),
      as.double(a1start),
      as.double(b1start),
      as.double(a2start),
      as.double(b2start),
                  
      # create empty vectors for saving the results
      as.double(rep(0,((size-burnin)/thin))),   # tau2_alpha0
      as.double(rep(0,((size-burnin)/thin))),   # tau2_alpha1
      as.double(rep(0,((size-burnin)/thin))),   # tau2_beta1
      as.double(rep(0,((size-burnin)/thin))),   # tau2_alpha2
      as.double(rep(0,((size-burnin)/thin))),   # tau2_beta2
      as.double(rep(0,((size-burnin)/thin))),   # sigma2
      
      as.double(rep(0,((size-burnin)/thin))),   # a0
      as.double(rep(0,((size-burnin)/thin))),   # a1      
      as.double(rep(0,((size-burnin)/thin))),   # b1
      as.double(rep(0,((size-burnin)/thin))),   # a2
      as.double(rep(0,((size-burnin)/thin))),   # b2

      as.double(rep(0,((size-burnin)/thin)*J)), # alpha0
      as.double(rep(0,((size-burnin)/thin)*J)), # alpha1
      as.double(rep(0,((size-burnin)/thin)*J)), # beta1
      as.double(rep(0,((size-burnin)/thin)*J)), # alpha2
      as.double(rep(0,((size-burnin)/thin)*J)), # beta2
      
      as.integer(0), # acc_b1
      as.integer(0), # acc_beta1
      as.integer(0), # acc_b2
      as.integer(0), # acc_beta2
      
      PACKAGE="frapmm")

      return(list(
      "timeall"=k2mhgibbs[[1]],
      "concall"=k2mhgibbs[[2]],
      "nobs"=k2mhgibbs[[3]],
      "ttotal"=k2mhgibbs[[4]],
      "J"=k2mhgibbs[[5]],
      "sv"=k2mhgibbs[[6]],
      "size"=k2mhgibbs[[7]],
      "burnin"=k2mhgibbs[[8]],
      "thin"=k2mhgibbs[[9]],
      "tune"=k2mhgibbs[[10]],
      "nobs_cum"=k2mhgibbs[[11]],
      "sigma2b1"=k2mhgibbs[[12]],
      "sigma2beta1"=k2mhgibbs[[13]],
      "sigma2b2"=k2mhgibbs[[14]],
      "sigma2beta2"=k2mhgibbs[[15]],
      "a0start"=k2mhgibbs[[16]],
      "a1start"=k2mhgibbs[[17]],
      "b1start"=k2mhgibbs[[18]],
      "a2start"=k2mhgibbs[[19]],
      "b2start"=k2mhgibbs[[20]],
      "tau2_alpha0"=k2mhgibbs[[21]],
      "tau2_alpha1"=k2mhgibbs[[22]],
      "tau2_beta1"=k2mhgibbs[[23]],
      "tau2_alpha2"=k2mhgibbs[[24]],
      "tau2_beta2"=k2mhgibbs[[25]],
      "sigma2"=k2mhgibbs[[26]],
      "a0"=k2mhgibbs[[27]],
      "a1"=k2mhgibbs[[28]],
      "b1"=k2mhgibbs[[29]],
      "a2"=k2mhgibbs[[30]],
      "b2"=k2mhgibbs[[31]],
      "alpha0"=k2mhgibbs[[32]],
      "alpha1"=k2mhgibbs[[33]],
      "beta1"=k2mhgibbs[[34]],
      "alpha2"=k2mhgibbs[[35]],
      "beta2"=k2mhgibbs[[36]],
      "acc_b1"=k2mhgibbs[[37]],
      "acc_beta1"=k2mhgibbs[[38]],
      "acc_b2"=k2mhgibbs[[39]],
      "acc_beta2"=k2mhgibbs[[40]]
      ))
}