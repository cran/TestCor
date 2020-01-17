
#' Simulates Gaussian data with a given correlation matrix and applies a FWER controlling procedure on the correlations.
#'
#' @useDynLib TestCor
#' @importFrom Rcpp sourceCpp
#'
#' @param corr_theo     the correlation matrix of Gaussien data simulated
#' @param n             sample size
#' @param Nsimu         number of simulations
#' @param alpha         level of multiple testing  
#' @param stat_test  
#' \describe{
#'   \item{'empirical'}{\eqn{\sqrt{n}*abs(corr)}}
#'   \item{'fisher'}{\eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
#'   \item{'student'}{\eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
#'   \item{'gaussian'}{\eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}}
#' }
#' @param method        choice between 'Bonferroni', 'Sidak', 'BootRW', 'MaxTinfty'
#' @param Nboot         number of iterations for Monte-Carlo of bootstrap quantile evaluation
#' @param stepdown      logical, if TRUE a stepdown procedure is applied
#' @param seed          seed for the Gaussian simulations
#'
#' @return Returns a line vector containing estimated fwer, estimated fdr, estimated power, estimated true discovery rate.
#'
#' @importFrom MASS mvrnorm
#' @export
#'
#' @references  Bonferroni, C. E. (1935). Il calcolo delle assicurazioni su gruppi di teste. Studi in onore del professore salvatore ortu carboni, 13-60.
#' @references  Drton, M., & Perlman, M. D. (2007). Multiple testing and error control in Gaussian graphical model selection. Statistical Science, 22(3), 430-449.
#' @references  Romano, J. P., & Wolf, M. (2005). Exact and approximate stepdown methods for multiple hypothesis testing. Journal of the American Statistical Association, 100(469), 94-108.
#' @references  Roux, M. (2018). Graph inference by multiple testing with application to Neuroimaging, Ph.D., Université Grenoble Alpes, France, https://tel.archives-ouvertes.fr/tel-01971574v1.
#' @references  Westfall, P.H. & Young, S. (1993) Resampling-based multiple testing: Examples and methods for p-value adjustment, John Wiley & Sons, vol. 279.
#' @seealso ApplyFwerCor, SimuFwer_oracle, SimuFdr
#'
#' @examples
#' Nsimu <- 1000 
#' n <- 100
#' p <- 10
#' corr_theo <- diag(1,p)
#' alpha <- 0.05
#' res <- SimuFwer(corr_theo,n,Nsimu,alpha,stat_test='empirical',method='Bonferroni',stepdown=FALSE)
SimuFwer <- function(corr_theo,n=100,Nsimu=1,alpha=0.05,stat_test='empirical',method='MaxTinfty',Nboot=1000,stepdown=TRUE,seed=NULL){

  if((stat_test=='gaussian')&(method=='MaxTinfty')){
        stop('MaxTinfty procedure is not implemented for Gaussian type statistics.\n')
  }


  if(!is.null(seed)){ set.seed(seed) }
  p <- nrow(corr_theo)

  corr_theo_vect <- vectorize(corr_theo)
  H1 <- (corr_theo_vect!=0)

  # we simulate all data outside the loop to improve time
  data_all <- mvrnorm(n*Nsimu,rep(0,p),corr_theo)

	false_positive <- array(0,dim=c(Nsimu))
	true_positive <- array(0,dim=c(Nsimu))
    nb_rejet <- array(0,dim=c(Nsimu))
	for(nsimu in 1:Nsimu){	
	   data <- data_all[(1:n)+n*(nsimu-1), ]
	   rejet <- ApplyFwerCor(data,alpha,stat_test,method,Nboot,stepdown,vect=TRUE)
	   true_positive[nsimu] <- ifelse(sum(H1)==0,0,sum(H1[rejet])/sum(H1))
       false_positive[nsimu] <- ifelse(sum(rejet)==0,0,(sum(rejet)-sum(H1[rejet]))/sum(rejet))
       nb_rejet[nsimu] <- sum(rejet)
	}

  # return fwer, fdr, power, tdr
  res  <- c(mean(false_positive>0),mean(false_positive),mean(true_positive),mean(ifelse(nb_rejet==0,0,true_positive/nb_rejet*sum(H1))))
  names(res) <- c('fwer','fdr','power','tdr')
  

  return(res)


}

#' Simulates Gaussian data with a given correlation matrix and applies oracle MaxTinfty on the correlations.
#' @description Simulates Gaussian data with a given correlation matrix and applies oracle MaxTinfty (i.e. Drton & Perlman (2007)'s procedure with the true correlation matrix) on the correlations.
#'
#' @param corr_theo     the correlation matrix of Gaussien data simulated
#' @param n             sample size
#' @param Nsimu         number of simulations
#' @param alpha         level of multiple testing  
#' @param stat_test  
#' \describe{
#'   \item{'empirical'}{\eqn{\sqrt{n}*abs(corr)}}
#'   \item{'fisher'}{\eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
#'   \item{'student'}{\eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
#'   \item{'gaussian'}{\eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}}
#' }
#' @param method        only 'MaxTinfty' available
#' @param Nboot         number of iterations for Monte-Carlo of bootstrap quantile evaluation
#' @param stepdown      logical, if TRUE a stepdown procedure is applied
#' @param seed          seed for the Gaussian simulations
#'
#' @return Returns a line vector containing estimated fwer, estimated fdr, estimated power, estimated true discovery rate.
#'
#' @references  Drton, M., & Perlman, M. D. (2007). Multiple testing and error control in Gaussian graphical model selection. Statistical Science, 22(3), 430-449.
#' @references  Roux, M. (2018). Graph inference by multiple testing with application to Neuroimaging, Ph.D., Université Grenoble Alpes, France, https://tel.archives-ouvertes.fr/tel-01971574v1.
#' @export
#' @seealso ApplyFwerCor_Oracle, SimuFwer
#'
#' @examples
#' Nsimu <- 1000
#' n <- 50
#' p <- 10
#' corr_theo <- diag(1,p)
#' alpha <- 0.05
#' res <- SimuFwer_oracle(corr_theo,n,Nsimu,alpha,stat_test='empirical',stepdown=FALSE,Nboot=100)
SimuFwer_oracle <- function(corr_theo,n=100,Nsimu=1,alpha=0.05,stat_test='empirical',method='MaxTinfty',Nboot=1000,stepdown=TRUE,seed=NULL){

  if(method!='MaxTinfty'){
       stop('The only oracle procedure available is MaxTinfty.\n')
  }
  if(stat_test=='gaussian'){
       stop('MaxTinfty procedure is not implemented for Gaussian type statistics.\n')
  }


  if(!is.null(seed)){ set.seed(seed) }
  p <- nrow(corr_theo)

  corr_theo_vect <- vectorize(corr_theo)
  H1 <- (corr_theo_vect!=0)

	# we simulate all data outside the loop to improve time
  data_all <- mvrnorm(n*Nsimu,rep(0,p),corr_theo)

	false_positive <- array(0,dim=c(Nsimu))
	true_positive <- array(0,dim=c(Nsimu))
    nb_rejet <- array(0,dim=c(Nsimu))
	for(nsimu in 1:Nsimu){	
	   data <- data_all[(1:n)+n*(nsimu-1), ]
	   rejet <- ApplyFwerCor_oracle(data,corr_theo,alpha,stat_test,method,Nboot,stepdown,vect=TRUE)
	   true_positive[nsimu] <- ifelse(sum(rejet)==0,0,sum(H1[rejet])/sum(rejet))
       false_positive[nsimu] <- ifelse(sum(rejet)==0,0,(sum(rejet)-sum(H1[rejet]))/sum(rejet))
       nb_rejet[nsimu] <- sum(rejet)
	}

  # return fwer, fdr, power, tdr
res  <- c(mean(false_positive>0),mean(false_positive),mean(true_positive),mean(ifelse(nb_rejet==0,0,true_positive/nb_rejet*sum(H1))))
  names(res) <- c('fwer','fdr','power','tdr')

	
  return(res)
}




#' Simulates Gaussian data with a given correlation matrix and applies a FDR controlling procedure on the correlations.
#'
#' @param corr_theo   the correlation matrix of Gaussien data simulated
#' @param n           sample size
#' @param Nsimu       number of simulations
#' @param alpha       level of multiple testing  
#' @param stat_test  
#' \describe{
#'   \item{'empirical'}{\eqn{\sqrt{n}*abs(corr)}}
#'   \item{'fisher'}{\eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
#'   \item{'student'}{\eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
#'   \item{'gaussian'}{\eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}}
#' }
#' @param method      choice between 'LCTnorm' and 'LCTboot', developped by Cai & Liu (2016), 
#'                    'BH', traditional Benjamini-Hochberg (1995)'s procedure,
#'                    and 'BHboot', Benjamini-Hochberg (1995)'s procedure with bootstrap evaluation of pvalues
#' @param Nboot       number of iterations for Monte-Carlo of bootstrap quantile evaluation
#' @param seed        seed for the Gaussian simulations
#'
#' @return Returns a line vector containing estimated fwer, estimated fdr, estimated power, estimated true discovery rate.
#' @importFrom MASS mvrnorm
#' @export
#'
#' @references Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the royal statistical society. Series B (Methodological), 289-300.
#' @references  Cai, T. T., & Liu, W. (2016). Large-scale multiple testing of correlations. Journal of the American Statistical Association, 111(513), 229-240.
#' @references  Roux, M. (2018). Graph inference by multiple testing with application to Neuroimaging, Ph.D., Université Grenoble Alpes, France, \url{https://tel.archives-ouvertes.fr/tel-01971574v1}.
#' @seealso ApplyFdrCor, SimuFwer
#'
#' @examples
#' Nsimu <- 1000
#' n <- 100
#' p <- 10
#' corr_theo <- diag(1,p)
#' alpha <- 0.05
#' res <- SimuFdr(corr_theo,n,Nsimu,alpha,stat_test='empirical',method='LCTnorm')
SimuFdr <- function(corr_theo,n=100,Nsimu=1,alpha=0.05,stat_test='empirical',method='LCTnorm',Nboot=1000,seed=NULL){

  if(!is.null(seed)){ set.seed(seed) }
  p <- nrow(corr_theo)

  corr_theo_vect <- vectorize(corr_theo)
  H1 <- (corr_theo_vect!=0)

  # we simulate all data outside the loop to improve time
  data_all <- mvrnorm(n*Nsimu,rep(0,p),corr_theo)

	false_positive <- array(0,dim=c(Nsimu))
	true_positive <- array(0,dim=c(Nsimu))
    nb_rejet <- array(0,dim=c(Nsimu))
	for(nsimu in 1:Nsimu){	
	   data <- data_all[(1:n)+n*(nsimu-1), ]
	   rejet <- ApplyFdrCor(data,alpha,stat_test,method,Nboot,vect=TRUE)
	   true_positive[nsimu] <- ifelse(sum(H1)==0,0,sum(H1[rejet])/sum(H1))
       false_positive[nsimu] <- ifelse(sum(rejet)==0,0,(sum(rejet)-sum(H1[rejet]))/sum(rejet))
       nb_rejet[nsimu] <- sum(rejet)
	}

  # return fwer, fdr, power, tdr
  res  <- c(mean(false_positive>0),mean(false_positive),mean(true_positive),mean(ifelse(nb_rejet==0,0,true_positive/nb_rejet*sum(H1))))
  names(res) <- c('fwer','fdr','power','tdr')

  return(res)
}

