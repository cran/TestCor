#--------------------------------method 1 : Bonferroni--------------------------------------
#--------------------------------------------------------------------------------------------

#' Bonferroni multiple testing procedure for correlations.
#'
#' @param data         matrix of observations
#' @param alpha        level of multiple testing
#' @param stat_test  
#' \describe{
#'   \item{'empirical'}{\eqn{\sqrt{n}*abs(corr)}}
#'   \item{'fisher'}{   \eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
#'   \item{'student'}{  \eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
#'   \item{'gaussian'}{ \eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}}
#' }
#' @param vect         if TRUE returns a vector of TRUE/FALSE values, corresponding to \code{vectorize(cor(data))};
#'                     if FALSE, returns an array containing rows and columns of significative correlations 
#'
#' @return Returns \itemize{\item{a vector of logicals, equal to TRUE if the corresponding element of the statistic vector is rejected, if \code{vect=TRUE},} \item{a vector containing indexes \eqn{\lbrace(i,j),\,i<j\rbrace} for which correlation between variables \eqn{i} and \eqn{j} is significative, if \code{vect=FALSE}.}}
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats cor
#' @export
#'
#' @references  Bonferroni, C. E. (1935). Il calcolo delle assicurazioni su gruppi di teste. Studi in onore del professore salvatore ortu carboni, 13-60.
#' @references  Roux, M. (2018). Graph inference by multiple testing with application to Neuroimaging, Ph.D., Université Grenoble Alpes, France, https://tel.archives-ouvertes.fr/tel-01971574v1.
#' @seealso ApplyFwerCor, BonferroniCor_SD
#'
#' @examples 
#' n <- 100
#' p <- 10
#' corr_theo <- diag(1,p)
#' data <- MASS::mvrnorm(n,rep(0,p),corr_theo)
#' corr_mat <- cor(data)
#' corr_vect <- corr_mat[upper.tri(corr_mat)]
#' alpha <- 0.05
#' res <- BonferroniCor(data,alpha,stat_test='empirical')
BonferroniCor <- function(data,alpha=0.05,stat_test='empirical',vect=FALSE){

	stat <- abs(eval_stat(data,stat_test))
    m <- length(stat)
    t_bonf <- qnorm(1-alpha/(2*m))
    result <- ( stat > t_bonf )

  if(vect==TRUE){
    return(result)
  }else{
    n <- nrow(data)
    rows <- vectorize(matrix(1:n,nrow=n,ncol=n))
    columns <- vectorize(t(matrix(1:n,nrow=n,ncol=n)))
    return(cbind(rows[which(result)],columns[which(result)]))
  }
}

#--------------------------------method 2 : Sidak-------------------------------------------
#--------------------------------------------------------------------------------------------

#' Sidak multiple testing procedure for correlations.
#'
#' @param data         matrix of observations
#' @param alpha        level of multiple testing
#' @param stat_test  
#' \describe{
#'   \item{'empirical'}{\eqn{\sqrt{n}*abs(corr)}}
#'   \item{'fisher'}{   \eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
#'   \item{'student'}{  \eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
#'   \item{'gaussian'}{ \eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}}
#' }
#' @param vect         if TRUE returns a vector of TRUE/FALSE values, corresponding to \code{vectorize(cor(data))};
#'                     if FALSE, returns an array containing rows and columns of significative correlations 
#'
#' @return Returns \itemize{\item{a vector of logicals, equal to TRUE if the corresponding element of the statistic vector is rejected, if \code{vect=TRUE},} \item{a vector containing indexes \eqn{\lbrace(i,j),\,i<j\rbrace} for which correlation between variables \eqn{i} and \eqn{j} is significative, if \code{vect=FALSE}.}}
#'
#' @importFrom stats qnorm cor
#' @importFrom MASS mvrnorm
#' @export
#'
#' @references  Roux, M. (2018). Graph inference by multiple testing with application to Neuroimaging, Ph.D., Université Grenoble Alpes, France, https://tel.archives-ouvertes.fr/tel-01971574v1.
#' @references  Šidák, Z. (1967). Rectangular confidence regions for the means of multivariate normal distributions. Journal of the American Statistical Association, 62(318), 626-633.
#' @seealso ApplyFwerCor, SidakCor_SD
#'
#' @examples  
#' n <- 100
#' p <- 10
#' corr_theo <- diag(1,p)
#' data <- MASS::mvrnorm(n,rep(0,p),corr_theo)
#' corr_mat <- cor(data)
#' corr_vect <- corr_mat[upper.tri(corr_mat)]
#' alpha <- 0.05
#' res <- SidakCor(data,alpha,stat_test='empirical')
SidakCor <- function(data,alpha=0.05,stat_test='empirical',vect=FALSE){

	stat <- abs(eval_stat(data,stat_test))
    m <- length(stat)
    t_sidak <- qnorm(0.5*((1-alpha)^(1/m))+0.5)

    result <- (stat > t_sidak)

  if(vect==TRUE){
    return(result)
  }else{
    n <- nrow(data)
    rows <- vectorize(matrix(1:n,nrow=n,ncol=n))
    columns <- vectorize(t(matrix(1:n,nrow=n,ncol=n)))
    return(cbind(rows[which(result)],columns[which(result)]))
  }
}

#--------------------------------method 3 : Bootstrap Romano/Wolf---------------------------
#--------------------------------------------------------------------------------------------

#' Bootstrap multiple testing method of Romano & Wolf (2005) for correlations.
#' @description Multiple testing method based on the evaluation of quantile by bootstrap 
#' in the initial dataset (Romano & Wolf (2005)).
#'
#' @param data       matrix of observations
#' @param alpha      level of multiple testing
#' @param stat_test  
#' \describe{
#'   \item{'empirical'}{\eqn{\sqrt{n}*abs(corr)}}
#'   \item{'fisher'}{   \eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
#'   \item{'student'}{  \eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
#'   \item{'gaussian'}{ \eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}}
#' }
#' @param Nboot      number of iterations for Bootstrap quantile evaluation
#' @param vect       if TRUE returns a vector of TRUE/FALSE values, corresponding to \code{vectorize(cor(data))};
#'                   if FALSE, returns an array containing rows and columns of significative correlations 
#'
#' @return Returns \itemize{\item{a vector of logicals, equal to TRUE if the corresponding element of the statistic vector is rejected, if \code{vect=TRUE},} \item{a vector containing indexes \eqn{\lbrace(i,j),\,i<j\rbrace} for which correlation between variables \eqn{i} and \eqn{j} is significative, if \code{vect=FALSE}.}}
#'
#' @importFrom stats cor quantile
#' @importFrom MASS mvrnorm
#' 
#' @export
#'
#' @references  Romano, J. P., & Wolf, M. (2005). Exact and approximate stepdown methods for multiple hypothesis testing. Journal of the American Statistical Association, 100(469), 94-108.
#' @references  Roux, M. (2018). Graph inference by multiple testing with application to Neuroimaging, Ph.D., Université Grenoble Alpes, France, https://tel.archives-ouvertes.fr/tel-01971574v1.
#' @seealso ApplyFwerCor, BootRWCor_SD
#'
#' @examples  
#' n <- 100
#' p <- 10
#' corr_theo <- diag(1,p)
#' data <- MASS::mvrnorm(n,rep(0,p),corr_theo)
#' alpha <- 0.05
#' res <- BootRWCor(data,alpha,stat_test='empirical',Nboot=1000)
BootRWCor <- function(data,alpha=0.05,stat_test='empirical',Nboot=1000,vect=FALSE){

    # test statistic	
    n <- nrow(data)
    p <- ncol(data)
    stat <- eval_stat(data,stat_test)

    # evaluation of quantile by bootstrap
    max_boot <- rep(0,Nboot)
	
 	for (nboot in 1:Nboot){
		
		indb <- sample(seq(1,n,1),replace=TRUE)
		datab <- data[indb,] 		
		stat_boot <- eval_stat(datab,stat_test)
        stat_boot <- abs(stat_boot-stat)
	
	  max_boot[nboot] <- max(stat_boot)
  }
	
    t_boot <- quantile(max_boot,1-alpha,names=FALSE)
	
    result <- (abs(stat) > t_boot)

  if(vect==TRUE){
    return(result)
  }else{
    n <- nrow(data)
    rows <- vectorize(matrix(1:n,nrow=n,ncol=n))
    columns <- vectorize(t(matrix(1:n,nrow=n,ncol=n)))
    return(cbind(rows[which(result)],columns[which(result)]))
  }

}

#--------------------------------method 4 : maxT_infty--------------------------------------
#--------------------------------------------------------------------------------------------

#' Multiple testing method of Drton & Perlman (2007) for correlations.
#' @description Multiple testing method based on the evaluation of quantile by simulation of observations 
#' from the asymptotic distribution (Drton & Perlman (2007)).
#'
#' @param data         matrix of observations
#' @param alpha        level of multiple testing
#' @param stat_test  
#' \describe{
#'   \item{'empirical'}{\eqn{\sqrt{n}*abs(corr)}}
#'   \item{'fisher'}{   \eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
#'   \item{'student'}{  \eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
#' }
#' Notice that 'gaussian' is not available.
#' @param Nboot        number of iterations for Monte-Carlo quantile evaluation
#' @param OmegaChap    matrix of covariance of empirical correlations used for quantile evaluation;
#'                     optional, useful for oracle estimation and step-down
#' @param vect         if TRUE returns a vector of TRUE/FALSE values, corresponding to \code{vectorize(cor(data))};
#'                     if FALSE, returns an array containing rows and columns of significative correlations  
#'
#' @return Returns a vector of logicals, equal to TRUE if the corresponding element of stat is rejected.
#'
#' @importFrom stats quantile cor
#' @importFrom MASS mvrnorm
#' @export
#'
#' @references  Drton, M., & Perlman, M. D. (2007). Multiple testing and error control in Gaussian graphical model selection. Statistical Science, 22(3), 430-449.
#' @references  Roux, M. (2018). Graph inference by multiple testing with application to Neuroimaging, Ph.D., Université Grenoble Alpes, France, https://tel.archives-ouvertes.fr/tel-01971574v1.
#' @seealso ApplyFwerCor, maxTinftyCor_SD
#'
#' @examples  
#' n <- 100
#' p <- 10
#' corr_theo <- diag(1,p)
#' data <- MASS::mvrnorm(n,rep(0,p),corr_theo)
#' alpha <- 0.05
#' res <- maxTinftyCor(data,alpha,stat_test='empirical',Nboot=1000)
maxTinftyCor <- function(data,alpha=0.05,stat_test='empirical',Nboot=1000,OmegaChap = covDcorNorm(cor(data),stat_test),vect=FALSE){


    if(stat_test=='gaussian'){
        stop('MaxTinfty procedure is not implemented for Gaussian type statistics.\n')
    }

    stat <- abs(eval_stat(data,stat_test))
    OmegaChap <- as.matrix(OmegaChap)
 
    # evaluation of the (1-alpha/2)-quantile of a N(0,OmegaChap) by simulation
    dataq <- mvrnorm(Nboot,rep(0,nrow(OmegaChap)),OmegaChap)
    maxq <- apply(as.matrix(dataq),1,function(x){max(abs(x))})
    t_maxTinfty <- quantile(maxq,1-alpha,names=FALSE)
	
    result <- (stat > t_maxTinfty)

   if(vect==TRUE){
     return(result)
   }else{
     n <- nrow(data)
     rows <- vectorize(matrix(1:n,nrow=n,ncol=n))
     columns <- vectorize(t(matrix(1:n,nrow=n,ncol=n)))
     return(cbind(rows[which(result)],columns[which(result)]))
   }

}



