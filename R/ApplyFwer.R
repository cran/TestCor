#' Applies multiple testing procedures controlling (asymptotically) the FWER
#' for tests on a correlation matrix.
#' @description Applies multiple testing procedures controlling (asymptotically) the FWER
#' for tests on a correlation matrix.
#' Methods are described in Chapter 5 of \cite{roux}.
#'
#'
#' @return Returns the list of significative correlations according to the multiple testing procedure applied
#
#' @param data        matrix of observations
#' @param alpha       level of multiple testing
#' @param stat_test  
#' \describe{
#'   \item{'empirical'}{\eqn{\sqrt{n}*abs(corr)}}
#'   \item{'fisher'}{   \eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
#'   \item{'student'}{  \eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
#'   \item{'gaussian'}{ \eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}}
#' }
#' @param method      choice between 'Bonferroni', 'Sidak', 'BootRW', 'MaxTinfty'
#' @param Nboot       number of iterations for Monte-Carlo of bootstrap quantile evaluation
#' @param stepdown    logical, if TRUE a stepdown procedure is applied
#' @param vect        if TRUE returns a vector of TRUE/FALSE values, corresponding to \code{vectorize(cor(data))};
#'                    if FALSE, returns an array containing rows and columns of significative correlations 
#'
#' @importFrom stats cor
#' @importFrom MASS mvrnorm
#' @export
#'
#' @references Bonferroni, C. E. (1935). Il calcolo delle assicurazioni su gruppi di teste. Studi in onore del professore salvatore ortu carboni, 13-60.
#' @references  Drton, M., & Perlman, M. D. (2007). Multiple testing and error control in Gaussian graphical model selection. Statistical Science, 22(3), 430-449.
#' @references  Romano, J. P., & Wolf, M. (2005). Exact and approximate stepdown methods for multiple hypothesis testing. Journal of the American Statistical Association, 100(469), 94-108.
#' @references  Roux, M. (2018). Graph inference by multiple testing with application to Neuroimaging, Ph.D., Université Grenoble Alpes, France, https://tel.archives-ouvertes.fr/tel-01971574v1.
#' @references  Šidák, Z. (1967). Rectangular confidence regions for the means of multivariate normal distributions. Journal of the American Statistical Association, 62(318), 626-633.
#' @seealso ApplyFwerCor_SD, ApplyFdrCor
#' @seealso BonferroniCor, SidakCor, BootRWCor, maxTinftyCor
#' @seealso BonferroniCor_SD, SidakCor_SD, BootRWCor_SD, maxTinftyCor_SD 
#'
#' @examples
#' n <- 100
#' p <- 10
#' corr_theo <- diag(1,p)
#' data <- MASS::mvrnorm(n,rep(0,p),corr_theo)
#' alpha <- 0.05
#' res <- ApplyFwerCor(data,alpha,stat_test='empirical',method='Bonferroni',stepdown=FALSE)
ApplyFwerCor <- function(data,alpha=0.05,stat_test='empirical',method='MaxTinfty',Nboot=1000,stepdown=TRUE, vect=FALSE){

 if(method=='Bonferroni'){
	if(stepdown==FALSE){ res <- BonferroniCor(data,alpha,stat_test,vect=TRUE) }
	else{ res <- BonferroniCor_SD(data,alpha,stat_test,vect=TRUE) }
 } 
 if(method=='Sidak'){
	if(stepdown==FALSE){ res <- SidakCor(data,alpha,stat_test,vect=TRUE) }
	else{ res <- SidakCor_SD(data,alpha,stat_test,vect=TRUE) }
 }
 if(method=='BootRW'){
	if(stepdown==FALSE){ res <- BootRWCor(data,alpha,stat_test,Nboot,vect=TRUE) }
	else{ res <- BootRWCor_SD(data,alpha,stat_test,Nboot,vect=TRUE) }
 }
 if(method=='MaxTinfty'){
    if(stat_test=='gaussian'){
        stop('MaxTinfty procedure is not implemented for Gaussian type statistics.\n')
    }else{
	    if(stepdown==FALSE){ res <- maxTinftyCor(data,alpha,stat_test,Nboot,vect=TRUE) }
	    else{ res <- maxTinftyCor_SD(data,alpha,stat_test,Nboot,vect=TRUE) }
    }
 }
 res <- as.logical(res)


 if(vect==TRUE){
   return(res)
 }else{
   p <- ncol(data)
   rows <- vectorize(matrix(1:p,nrow=p,ncol=p))
   columns <- vectorize(t(matrix(1:p,nrow=p,ncol=p)))
   return(cbind(rows[which(res)],columns[which(res)]))
 }

}
 
#' Applies an oracle version of MaxTinfty procedure described in Drton & Perlman (2007) for correlation testing.             
#' @description Applies oracle MaxTinfty procedure described in Drton & Perlman (2007) which controls asymptotically the FWER
#' for tests on a correlation matrix. It needs the true correlation matrix. 
#'
#' @return Returns the list of significative correlations according to the multiple testing procedure applied
#' Oracle estimation of the quantile is used, based on the true correlation matrix
#'
#' @param data          matrix of observations
#' @param corr_theo     true matrix of corrlations
#' @param alpha         level of multiple testing
#' @param stat_test     
#' \describe{
#'   \item{'empirical'}{\eqn{\sqrt{n}*abs(corr)}}
#'   \item{'fisher'}{\eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
#'   \item{'student'}{\eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
#'   \item{'gaussian'}{\eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}}
#' }
#'                      
#' @param method        only 'MaxTinfty' implemented
#' @param Nboot         number of iterations for Monte-Carlo of bootstrap quantile evaluation
#' @param vect          optional, logical, if TRUE returns a vector of TRUE/FALSE values, corresponding to \code{vectorize(cor(data))};
#'                     if FALSE, returns an array containing rows and columns of significative correlations 
#' @param stepdown      logical, if TRUE a stepdown procedure is applied
#'                      if FALSE, returns an array containing rows and columns of significative correlations   
#'                      
#' @importFrom stats cor
#' @importFrom MASS mvrnorm
#' 
#' @references  Drton, M., & Perlman, M. D. (2007). Multiple testing and error control in Gaussian graphical model selection. Statistical Science, 22(3), 430-449.
#' @references  Roux, M. (2018). Graph inference by multiple testing with application to Neuroimaging, Ph.D., Université Grenoble Alpes, France, https://tel.archives-ouvertes.fr/tel-01971574v1.
#' @seealso ApplyFwerCor 
#' @seealso maxTinftyCor, maxTinftyCor_SD 
#'
#' @export
#' @examples
#' n <- 100
#' p <- 10
#' corr_theo <- diag(1,p)
#' data <- MASS::mvrnorm(n,rep(0,p),corr_theo)
#' alpha <- 0.05
#' res <- ApplyFwerCor_oracle(data,corr_theo,alpha,stat_test='empirical',Nboot=1000,stepdown=FALSE)
ApplyFwerCor_oracle <- function(data,corr_theo,alpha=0.05,stat_test='empirical',method='MaxTinfty',Nboot=1000,stepdown=TRUE,vect=FALSE){

 Gtheo <- covDcorNorm(corr_theo,stat_test)
 if(method=='MaxTinfty'){
    if(stat_test=='gaussian'){
        stop('MaxTinfty procedure is not implemented for Gaussian type statistics.\n')
    }else{
    	if(stepdown==FALSE){ res <- maxTinftyCor(data,alpha,stat_test,Nboot,Gtheo,vect=TRUE) }
	    else{ res <- maxTinftyCor_SD(data,alpha,stat_test,Nboot,Gtheo,vect=TRUE) }
    }
 }
 else{ stop('The method is not implemented with oracle version.\n') }
 res <- as.logical(res)

 if(vect==TRUE){
   return(res)
 }else{
   p <- ncol(data)
   rows <- vectorize(matrix(1:p,nrow=p,ncol=p))
   columns <- vectorize(t(matrix(1:p,nrow=p,ncol=p)))
   return(cbind(rows[which(res)],columns[which(res)]))
 }

}

#' Applies multiple testing procedures built to control (asymptotically) the FDR for correlation testing.
#' @description Applies multiple testing procedures built to control (asymptotically) the FDR for correlation testing.
#' Some have no theoretical proofs for tests on a correlation matrix.
#'
#'
#' @return Returns the list of significative correlations according to the multiple testing procedure applied.
#
#' @param data      matrix of observations
#' @param alpha     level of multiple testing
#' @param stat_test  
#' \describe{
#'   \item{'empirical'}{\eqn{\sqrt{n}*abs(corr)}}
#'   \item{'fisher'}{\eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
#'   \item{'student'}{\eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
#'   \item{'gaussian'}{\eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}}
#' }
#' @param method    choice between 'LCTnorm' and 'LCTboot' developped by Cai & Liu (2016), 
#'                  'BH', traditional Benjamini-Hochberg's procedure Benjamini & Hochberg (1995)'s
#'                  and 'BHboot', Benjamini-Hochberg (1995)'s procedure with bootstrap evaluation of p-values
#' @param Nboot     number of iterations for bootstrap p-values evaluation
#' @param vect      if TRUE returns a vector of TRUE/FALSE values, corresponding to \code{vectorize(cor(data))};
#'                  if FALSE, returns an array containing rows and columns of significative correlations 
#'
#' @importFrom stats cor
#' @importFrom MASS mvrnorm
#' @export
#'
#' @references Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the royal statistical society. Series B (Methodological), 289-300.
#' @references  Cai, T. T., & Liu, W. (2016). Large-scale multiple testing of correlations. Journal of the American Statistical Association, 111(513), 229-240.
#' @references  Roux, M. (2018). Graph inference by multiple testing with application to Neuroimaging, Ph.D., Université Grenoble Alpes, France, https://tel.archives-ouvertes.fr/tel-01971574v1.
#' @seealso ApplyFwerCor 
#' @seealso LCTnorm, LCTboot, BHCor, BHBootCor
#'
#' @examples
#' n <- 100
#' p <- 10
#' corr_theo <- diag(1,p)
#' data <- MASS::mvrnorm(n,rep(0,p),corr_theo)
#' alpha <- 0.05
#' res <- ApplyFdrCor(data,alpha,stat_test='empirical',method='LCTnorm')
ApplyFdrCor <- function(data,alpha=0.05,stat_test='empirical',method='LCTnorm',Nboot=1000,vect=FALSE){

 if(method=='LCTnorm'){
        res <- LCTnorm(data,alpha=alpha,stat_test=stat_test,vect=TRUE)
 } 
 if(method=='LCTboot'){
        res <- LCTboot(data,alpha=alpha,stat_test=stat_test,Nboot=Nboot,vect=TRUE)
 } 
 if(method=='BH'){
       res <- BHCor(data,alpha=alpha,stat_test=stat_test,vect=TRUE)
 }
if(method=='BHboot'){
       res <- BHBootCor(data,alpha=alpha,stat_test=stat_test,Nboot=Nboot,vect=TRUE)
 }
 res <- as.logical(res)


 if(vect==TRUE){
   return(res)
 }else{
   p <- ncol(data)
   rows <- vectorize(matrix(1:p,nrow=p,ncol=p))
   columns <- vectorize(t(matrix(1:p,nrow=p,ncol=p)))
   return(cbind(rows[which(res)],columns[which(res)]))
 }

}


