# Multiple testing procedures controlling (asymptotically) the FWER
# for tests on a correlation matrix ** with Stepdown **
###############################################################################

#--------------------------------method 1 : Bonferroni--------------------------------------
#--------------------------------------------------------------------------------------------

#' Bonferroni multiple testing method for correlations 
#' with stepdown procedure.
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
#' @param vect      if TRUE returns a vector of TRUE/FALSE values, corresponding to \code{vectorize(cor(data))};
#'                    if FALSE, returns an array containing rows and columns of significative correlations 
#'
#' @return Returns \itemize{\item{a vector of logicals, equal to TRUE if the corresponding element of the statistic vector is rejected, if \code{vect=TRUE},} \item{a vector containing indexes \eqn{\lbrace(i,j),\,i<j\rbrace} for which correlation between variables \eqn{i} and \eqn{j} is significative, if \code{vect=FALSE}.}}
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats cor
#' @export
#'
#' @references Bonferroni, C. E. (1935). Il calcolo delle assicurazioni su gruppi di teste. Studi in onore del professore salvatore ortu carboni, 13-60.
#' @references  Roux, M. (2018). Graph inference by multiple testing with application to Neuroimaging, Ph.D., Université Grenoble Alpes, France, https://tel.archives-ouvertes.fr/tel-01971574v1.
#' @seealso ApplyFwerCor, BonferroniCor
#'
#' @examples  
#' n <- 100
#' p <- 10
#' corr_theo <- diag(1,p)
#' data <- MASS::mvrnorm(n,rep(0,p),corr_theo)
#' alpha <- 0.05
#' res <- BonferroniCor_SD(data,alpha,stat_test='empirical')
BonferroniCor_SD <- function(data,alpha,stat_test='empirical',vect=FALSE){

	stat <- abs(eval_stat(data, stat_test))
    m <- length(stat)
    t_bonf <- qnorm(1 - alpha/(2 * m))
	vBonf <- (stat > t_bonf)

    res_SD <- rep(0, length(vBonf))
    indNR <- which(vBonf == 0)
    indR <- which(vBonf != 0)
    res_SD[indR] <- 1
	
  while((sum(vBonf)!=0)&&(sum(indNR)!=0)){
      stat_SD <- stat[indNR]
      m <- length(stat_SD)
      t_bonf <- qnorm(1 - alpha/(2 * m))
	  vBonf <- (stat_SD > t_bonf)

	  indR <- which(vBonf !=0)            
	  res_SD[indNR[indR]] <- 1
	  indNR <- indNR[-indR]	
  }
	
 res_SD <- as.logical(res_SD)
 if(vect==TRUE){
   return(res_SD)
 }else{
   p <- ncol(data)
   rows <- vectorize(matrix(1:p,nrow=p,ncol=p))
   columns <- vectorize(t(matrix(1:p,nrow=p,ncol=p)))
   return(cbind(rows[which(res_SD)],columns[which(res_SD)]))
 }
}


#--------------------------------method 2 : Sidak-------------------------------------------
#--------------------------------------------------------------------------------------------

#' Sidak multiple testing method for correlations
#' with stepdown procedure.
#'
#' @param data         matrix of observations
#' @param alpha        level of multiple testing
#' @param stat_test  
#' \describe{
#'   \item{'empirical'}{\eqn{\sqrt{n}*abs(corr)}}
#'   \item{'fisher'}{\eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
#'   \item{'student'}{\eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
#'   \item{'gaussian'}{\eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}}
#' }
#' @param vect        if TRUE returns a vector of TRUE/FALSE values, corresponding to \code{vectorize(cor(data))};
#'                    if FALSE, returns an array containing rows and columns of significative correlations 
#'
#' @return Returns \itemize{\item{a vector of logicals, equal to TRUE if the corresponding element of the statistic vector is rejected, if \code{vect=TRUE},} \item{a vector containing indexes \eqn{\lbrace(i,j),\,i<j\rbrace} for which correlation between variables \eqn{i} and \eqn{j} is significative, if \code{vect=FALSE}.}}
#'
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom stats cor
#'
#' @references  Roux, M. (2018). Graph inference by multiple testing with application to Neuroimaging, Ph.D., Université Grenoble Alpes, France, https://tel.archives-ouvertes.fr/tel-01971574v1.
#' @references Šidák, Z. (1967). Rectangular confidence regions for the means of multivariate normal distributions. Journal of the American Statistical Association, 62(318), 626-633.
#' @seealso ApplyFwerCor, SidakCor
#'
#' @examples   
#' n <- 100
#' p <- 10
#' corr_theo <- diag(1,p)
#' data <- MASS::mvrnorm(n,rep(0,p),corr_theo)
#' alpha <- 0.05
#' res <- SidakCor_SD(data,alpha,stat_test='empirical')
SidakCor_SD <- function(data,alpha,stat_test='empirical',vect=FALSE){

	stat <- abs(eval_stat(data,stat_test))
    m <- length(stat)
    t_sidak <- qnorm(0.5*((1-alpha)^(1/m))+0.5)
    vSidak  <- (stat > t_sidak)

  res_SD <- rep(0,length(vSidak))
  indNR <- which(vSidak==0)   # indexes of non rejected hypothesis 
  indR <- which(vSidak!=0)               
  res_SD[indR] <- 1
	
  while((sum(vSidak)!=0)&&(sum(indNR)!=0)){
      stat_SD <- stat[indNR]
      m <- length(stat_SD)
      t_sidak <- qnorm(0.5*((1-alpha)^(1/m))+0.5)
      vSidak  <- (stat_SD > t_sidak)


	  indR <- which(vSidak !=0)            
	  res_SD[indNR[indR]] <- 1
	  indNR =indNR[-indR]
  }
	
 res_SD <- as.logical(res_SD)
 if(vect==TRUE){
   return(res_SD)
 }else{
   p <- ncol(data)
   rows <- vectorize(matrix(1:p,nrow=p,ncol=p))
   columns <- vectorize(t(matrix(1:p,nrow=p,ncol=p)))
   return(cbind(rows[which(res_SD)],columns[which(res_SD)]))
 }
}



#--------------------------------method 3 : Bootstrap Romano/Wolf---------------------------
#--------------------------------------------------------------------------------------------

#' Boootstrap multiple testing method of Romano & Wolf (2005) for correlations, with stepdown procedure.
#'
#' @description  Multiple testing method based on the evaluation of quantile by bootstrap 
#' in the initial dataset (Romano & Wolf (2005)),
#' with stepdown procedure.
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
#' @param Nboot        number of iterations for Bootstrap quantile evaluation
#' @param vect         if TRUE returns a vector of TRUE/FALSE values, corresponding to \code{vectorize(cor(data))};
#'                     if FALSE, returns an array containing rows and columns of significative correlations 
#'
#' @return Returns \itemize{\item{a vector of logicals, equal to TRUE if the corresponding element of the statistic vector is rejected, if \code{vect=TRUE},} \item{a vector containing indexes \eqn{\lbrace(i,j),\,i<j\rbrace} for which correlation between variables \eqn{i} and \eqn{j} is significative, if \code{vect=FALSE}.}}
#' 
#' @importFrom stats cor
#' @importFrom MASS mvrnorm
#' @export
#'
#' @references Romano, J. P., & Wolf, M. (2005). Exact and approximate stepdown methods for multiple hypothesis testing. Journal of the American Statistical Association, 100(469), 94-108.
#' @references  Roux, M. (2018). Graph inference by multiple testing with application to Neuroimaging, Ph.D., Université Grenoble Alpes, France, https://tel.archives-ouvertes.fr/tel-01971574v1.
#' @seealso ApplyFwerCor, BootRWCor
#' 
#' @examples  
#' n <- 100
#' p <- 10
#' corr_theo <- diag(1,p)
#' data <- MASS::mvrnorm(n,rep(0,p),corr_theo)
#' alpha <- 0.05
#' res <- BootRWCor_SD(data,alpha,stat_test='empirical',Nboot=1000)
BootRWCor_SD <- function(data,alpha,stat_test='empirical',Nboot=1000,vect=FALSE){

  vBootRW <- BootRWCor(data,alpha,stat_test,Nboot,vect=TRUE)

  res_SD <- rep(0,length(vBootRW))
  indNR <- which(vBootRW==0)   # indexes of non rejected hypothesis 
  indR <- which(vBootRW!=0)               
  res_SD[indR] <- 1
	
  n <- ncol(data)

  stat <- eval_stat(data,stat_test)
   
  while((sum(vBootRW)!=0)&&(sum(indNR)!=0)){
		
    # test statistic
    stat_SD <- stat[indNR]

    # evaluation of quantile by bootstrap
    max_boot <- rep(0,Nboot)

    for (nboot in 1:Nboot){

       indb <- sample(seq(1,n,1),replace=TRUE)
	   data_boot <- data[indb,]
	   stat_boot <- eval_stat(data_boot,stat_test)

      stat_boot <- abs(stat_boot[indNR]-stat_SD)
      max_boot[nboot] <- max(stat_boot)
  }

  t_boot <- quantile(max_boot,1-alpha,names=FALSE)
  vBootRW <- (abs(stat_SD) > t_boot)

  indR <- which(vBootRW !=0)            
  res_SD[indNR[indR]] <- 1
  indNR <- indNR[-indR]
	
  }
	
 res_SD <- as.logical(res_SD)
 if(vect==TRUE){
   return(res_SD)
 }else{
   p <- ncol(data)
   rows <- vectorize(matrix(1:p,nrow=p,ncol=p))
   columns <- vectorize(t(matrix(1:p,nrow=p,ncol=p)))
   return(cbind(rows[which(res_SD)],columns[which(res_SD)]))
 }
	
}



#--------------------------------method 4 : maxT_infty--------------------------------------
#--------------------------------------------------------------------------------------------
#' Multiple testing method of Drton & Perlman (2007) for correlations, with stepdown procedure.
#'
#' @description Multiple testing method based on the evaluation of quantile by simulation of observations 
#' from the asymptotic distribution (Drton & Perlman (2007)), 
#' with stepdown procedure.
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
#' @param OmegaChap    matrix of covariance of test statistics;
#'                     optional, useful for oracle estimation and step-down
#' @param vect         if TRUE returns a vector of TRUE/FALSE values, corresponding to \code{vectorize(cor(data))};
#'                     if FALSE, returns an array containing rows and columns of significative correlations 
#' 
#' @return Returns \itemize{\item{a vector of logicals, equal to TRUE if the corresponding element of the statistic vector is rejected, if \code{vect=TRUE},} \item{a vector containing indexes \eqn{\lbrace(i,j),\,i<j\rbrace} for which correlation between variables \eqn{i} and \eqn{j} is significative, if \code{vect=FALSE}.}}
#'
#' @export                     
#' @importFrom MASS mvrnorm
#' @importFrom stats cor
#'
#' @references Drton, M., & Perlman, M. D. (2007). Multiple testing and error control in Gaussian graphical model selection. Statistical Science, 22(3), 430-449.
#' @references  Roux, M. (2018). Graph inference by multiple testing with application to Neuroimaging, Ph.D., Université Grenoble Alpes, France, https://tel.archives-ouvertes.fr/tel-01971574v1.
#' @seealso ApplyFwerCor, maxTinftyCor
#'
#' @examples  
#' n <- 100
#' p <- 10
#' corr_theo <- diag(1,p)
#' data <- MASS::mvrnorm(n,rep(0,p),corr_theo)
#' alpha <- 0.05
#' res <- maxTinftyCor_SD(data,alpha,stat_test='empirical',Nboot=1000)
maxTinftyCor_SD <- function(data,alpha=0.05,stat_test='empirical',Nboot=1000,OmegaChap=covDcorNorm(cor(data),stat_test),vect=FALSE){ 

  if(stat_test=='gaussian'){
        stop('MaxTinfty procedure is not implemented for Gaussian type statistics.\n')
  }

  vmaxTinf <- maxTinftyCor(data,alpha,stat_test,Nboot,OmegaChap,vect=TRUE)

  stat <- abs(eval_stat(data,stat_test))
  res_SD <- rep(0,length(vmaxTinf))
  indNR <- which(vmaxTinf==0)   # indexes of non rejected hypothesis 
  indR <- which(vmaxTinf!=0)               
  res_SD[indR] <- 1

	
  while((sum(vmaxTinf)!=0)&&(sum(indNR)!=0)){
    stat_SD <- stat[indNR]
    OmegaChap_SD <- OmegaChap[indNR,indNR]
 
    # evaluation of the (1-alpha/2)-quantile of a N(0,OmegaChap) by simulation
    dataq <- mvrnorm(Nboot,rep(0,nrow(OmegaChap_SD)),OmegaChap_SD)
    maxq <- apply(as.matrix(dataq),1,function(x){max(abs(x))})
    t_maxTinfty <- quantile(maxq,1-alpha,names=FALSE)
    vmaxTinf <- (stat_SD > t_maxTinfty)
    
      indR <- which(vmaxTinf !=0)            
      res_SD[indNR[indR]] <- 1
      indNR <- indNR[-indR]
  }

 res_SD <- as.logical(res_SD)	
 if(vect==TRUE){
   return(res_SD)
 }else{
   p <- ncol(data)
   rows <- vectorize(matrix(1:p,nrow=p,ncol=p))
   columns <- vectorize(t(matrix(1:p,nrow=p,ncol=p)))
   return(cbind(rows[which(res_SD)],columns[which(res_SD)]))
 }
	
}



