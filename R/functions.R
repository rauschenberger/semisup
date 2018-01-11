
################################################################################
### User functions #############################################################
################################################################################

#-------------------------------------------------------------------------------
#--- User: testing -------------------------------------------------------------
#-------------------------------------------------------------------------------

#' @title
#' Hypothesis testing
#' 
#' @export
#' @keywords methods
#' 
#' @description
#' This function tests whether the unlabelled observations
#' come from a mixture of two distributions.
#' 
#' @usage
#' scrutor(Y, Z, dist = "norm",
#'         phi = NULL, pi = NULL, gamma = NULL,
#'         test = "perm", iter = NULL, kind = NULL,
#'         debug = TRUE, ...)
#' 
#' @inheritParams arguments
#' 
#' @param phi
#' dispersion parameter(s)\strong{:}
#' numeric vector of length \code{q},
#' or \code{NULL} (\code{norm:} none, \code{nbinom:} MLE)
#' 
#' @param pi
#' zero-inflation parameter(s)\strong{:}
#' numeric vector of length \code{q},
#' or \code{NULL} (\code{norm:} none,\code{nbinom:} MLE)
#' 
#' @details
#' By default, \code{phi} and \code{pi}
#' are estimated by the maximum likelihood method,
#' and \code{gamma} is replaced by a vector of ones.
#' 
#' @return
#' This function tests a one-component (\code{H0})
#' against a two-component mixture model (\code{H1}).
#' 
#' \item{y}{index observations}
#' \item{z}{index class labels}
#' \item{lrts}{test statistic}
#' \item{p.value}{\code{p}-value}
#' 
#' @section Reference:
#' A Rauschenberger, RX Menezes, MA van de Wiel,
#' NM van Schoor, and MA Jonker (2017).
#' "Detecting SNPs with interactive effects on a quantitative trait",
#' \emph{Manuscript in preparation}.
#' 
#' @seealso
#' Use \code{\link{mixtura}} for model fitting. 
#' All other functions are \code{\link{internal}}.
#' 
#' @examples
#' # data simulation
#' n <- 100
#' z <- rep(0:1,each=n/2)
#' y <- rnorm(n=n,mean=2*z,sd=1)
#' z[(n/4):n] <- NA
#' 
#' # hypothesis testing
#' scrutor(y,z,dist="norm")
#' 
scrutor <- function(Y,Z,dist="norm",phi=NULL,pi=NULL,gamma=NULL,test="perm",iter=NULL,kind=NULL,debug=TRUE,...){
    
    # initialisation
    if(class(Y) %in% c("RangedSummarizedExperiment","SummarizedExperiment")){
        Y <- SummarizedExperiment::assays(Y)$counts
    }
    Y <- as.matrix(Y)
    Z <- as.matrix(Z)
    
    # combinations
    py <- ncol(Y)
    pz <- ncol(Z)
    if(py==pz){
        comb <- data.frame(y=seq_len(py),z=seq_len(pz))
        p <- py
    } else {
        comb <- expand.grid(y=seq_len(py),z=seq_len(pz))
        comb <- comb[order(comb$y),]
        p <- py*pz
        rownames(comb) <- seq_len(p)
    }
    
    # default values
    if(is.null(iter)){iter <- p/0.05}
    if(is.null(kind)){kind <- 0.05/p}
    if(is.null(test)){iter <- kind <- 1}
    
    # debugging
    if(debug){
        debug(y=Y,z=Z,dist=dist,phi=phi,pi=pi,gamma=gamma,
              test=test,iter=iter,kind=kind,...)
        n <- nrow(Y)
        ties <- apply(Y,2,function(y) length(unique(y))!=n)
        if(any(ties)){
            noise <- stats::runif(n=sum(ties)*n,min=0,max=0.001)
            Y[,ties] <- Y[,ties,drop=FALSE] + noise
            warning("Argument \"Y\" contains ties.",call.=FALSE)
        }
    }
    
    # testing
    fit <- list()
    pb <- utils::txtProgressBar(min=0,max=p,width=0.5*getOption("width"),style=3)
    for(i in seq_len(p)){
        utils::setTxtProgressBar(pb=pb,value=i)
        fit[[i]] <- semisup::mixtura(y=Y[,comb$y[i]],
                            z=Z[,comb$z[i]],
                            dist=dist,
                            phi=phi[comb$y[i]],
                            pi=pi[comb$y[i]],
                            gamma=gamma,
                            kind=kind,iter=iter,
                            debug=FALSE,test=test,...)
    }
    close(pb)
    
    # output
    list <- list()
    list$y <- comb$y
    list$z <- comb$z
    list$tau <- sapply(fit,function(x) x$estim1$p1)
    if(dist=="norm"){
        list$delta <- sapply(fit,function(x) 
            abs(x$estim1$mean0-x$estim1$mean1)/abs(x$estim0$mean0))
    }
    if(dist %in% c("nbinom","zinb")){
        list$delta <- sapply(fit,function(x) 
            abs(x$estim1$mu0-x$estim1$mu1)/abs(x$estim0$mu0))
        list$phi <- sapply(fit,function(x) x$estim0$phi)
        if(dist=="zinb"){
            list$pi <- sapply(fit,function(x) x$estim0$pi)
        }
    }
    list$lrts <- sapply(fit,function(x) x$lrts)
    if(!is.null(test)){
        list$p.value <- sapply(fit,function(x) x$p.value)
    }
    as.data.frame(list)
}

#-------------------------------------------------------------------------------
#--- User: fitting -------------------------------------------------------------
#-------------------------------------------------------------------------------

#' @title
#' Model fitting
#' 
#' @export
#' @keywords methods
#' 
#' @description
#' This function fits a semi-supervised mixture model.
#' It simultaneously estimates two mixture components,
#' and assigns the unlabelled observations to these.
#' 
#' @usage
#' mixtura(y, z, dist = "norm",
#'         phi = NULL, pi = NULL, gamma = NULL,
#'         test = NULL, iter = 100, kind = 0.05,
#'         debug = TRUE, ...)
#'         
#' @inheritParams arguments
#' 
#' @details
#' By default, \code{phi} and \code{pi}
#' are estimated by the maximum likelihood method,
#' and \code{gamma} is replaced by a vector of ones.
#' 
#' @return
#' This function fits and compares a one-component (\code{H0})
#' and a two-component (\code{H1}) mixture model.
#' 
#' \item{posterior}{probability of belonging to class 1\strong{:}
#' numeric vector of length \code{n}}
#' \item{converge}{path of the log-likelihood\strong{:}
#' numeric vector with maximum length
#' \code{it.em}}
#' \item{estim0}{parameter estimates under \code{H0}\strong{:}
#' data frame}
#' \item{estim1}{parameter estimates under \code{H1}\strong{:}
#' data frame}
#' \item{loglik0}{log-likelihood under \code{H0}\strong{:}
#' numeric}
#' \item{loglik1}{log-likelihood under \code{H1}\strong{:}
#' numeric}
#' \item{lrts}{likelihood-ratio test statistic\strong{:}
#' positive numeric}
#' \item{p.value}{\code{H0} versus \code{H1}\strong{:}
#' numeric between \code{0} and \code{1}, or \code{NULL}}
#' 
#' @seealso
#' Use \code{\link{scrutor}} for hypothesis testing.
#' All other functions are \code{\link{internal}}.
#' 
#' @section Reference:
#' A Rauschenberger, RX Menezes, MA van de Wiel,
#' NM van Schoor, and MA Jonker (2017).
#' "Detecting SNPs with interactive effects on a quantitative trait",
#' \emph{Manuscript in preparation}.
#' 
#' @examples
#' # data simulation
#' n <- 100
#' z <- rep(0:1,each=n/2)
#' y <- rnorm(n=n,mean=2,sd=1)
#' z[(n/4):n] <- NA
#' 
#' # model fitting
#' mixtura(y,z,dist="norm",test="perm")
#' 
mixtura <- function(y,z,dist="norm",phi=NULL,pi=NULL,gamma=NULL,test=NULL,iter=100,kind=0.05,debug=TRUE,...){
    
    if(debug){
        debug(y=y,z=z,dist=dist,phi=phi,pi=pi,gamma=gamma,
                 test=test,iter=iter,kind=kind,...)
        if(length(unique(y))!=length(y)){
            noise <- stats::runif(n=length(y),min=0,max=0.001)
            y <- y + noise
            warning("Argument \"y\" contains ties.",call.=FALSE)
        }
    }
    
    # --- population parameters ------------------------------------------------
    
    pass <- list()
    pass$n <- n <- length(y)
    pass$keep <- which(z==1)
    
    if(dist=="norm"){
            pass$mean <- sum(y)/n
            pass$sd <- sqrt(sum((y-pass$mean)^2)/(n-1))
    } else if(dist=="nbinom"){
            if(is.null(gamma)){
                gamma <- rep(1,times=n)
            }
            pass$mu <- sum(y)/n
            if(is.null(phi)){
                phi <- semisup::estim.nbinom(y=y,z=z,gamma=gamma)$phi
            }
    } else if(dist=="zinb"){
        if(is.null(gamma)){
            gamma <- rep(1,times=n)
        }
        if(is.null(phi) & is.null(pi)){
            temp <- semisup::estim.zinb(y=round(y),z=z,gamma=gamma)
            phi <- temp$phi
            pi <- temp$pi
        }
        pass$mu <- sum(y)/((1-pi)*n)
    }
    
    #--- fitting ---------------------------------------------------------------
    
    fit <- semisup::fit.wrap(y=y,z=z,dist=dist,phi=phi,pi=pi,gamma=gamma,...)
    
    #--- testing ---------------------------------------------------------------
    
    if(!is.null(test)){
        if(kind==1){
            
            #- - testing without interruption- - - - - - - - - - - - - - - - - - 
    
            lrts <- rep(NA,times=iter)
            lrts[1] <- fit$lrts
        
            for(i in seq_len(iter)[-1]){
               lrts[i] <- semisup::resam.lrts(y=y,z=z,dist=dist,
                               phi=phi,pi=pi,gamma=gamma,
                               test=test,pass=pass,...)
            }
            fit$p.value <- sum(lrts >= lrts[1])/iter

            
        } else if(kind > 0 & kind < 1){
            
            #- - testing with interruption - - - - - - - - - - - - - - - - - - -
            
            lrts <- rep(NA,times=iter)
            lrts[1] <- fit$lrts
            target <- ceiling(kind * iter)
            i <- -Inf
            counts <- 1
            pos <- c(2^(seq_len(floor(log(iter, base = 2)))),iter+1)
            for (j in seq_len(length(pos) - 1)) {
                if (counts <= target & i <= iter) {
                    for (i in seq(from=pos[j],to=pos[j+1]-1,by=1)) {
                        lrts[i] <- semisup::resam.lrts(y=y,z=z,dist=dist,
                                        phi=phi,pi=pi,gamma=gamma,
                                        test=test,pass=pass,...)
                    }
                    counts <- counts + 
                                sum(lrts[pos[j]:(pos[j + 1] - 1)] >= lrts[1])
                } else {
                    i <- max(i,1)
                    break
                }
            }
            fit$p.value <- counts/i
            
        } else if(kind==0){
            
            #- - approximate distribution- - - - - - - - - - - - - - - - - - - -
            
            quan <- semisup::table$quantile[[dist]]
            pval <- semisup::table$pvalue[[dist]]
            iter <- semisup::table$iter[[dist]]
            temp <- pval[pmax(1,sum(fit$lrts > quan))]
            fit$p.value <- (temp*iter+1)/(iter+1)
            
        }
    }
    
    fit
}

################################################################################
### Internal functions #########################################################
################################################################################

#-------------------------------------------------------------------------------
#--- Internal: wrapper ---------------------------------------------------------
#-------------------------------------------------------------------------------

#' @title
#' Internal function
#' 
#' @export
#' @keywords internal
#'
#' @description
#' This function fits the semi-supervised mixture model multiple times.
#' It is called by \code{\link{mixtura}} and \code{\link{scrutor}}.
#' 
#' @inheritParams arguments
#' 
#' @details
#' The distributions are parametrised as follows:
#' \itemize{
#' \item Gaussian \cr
#' \code{y ~ N(mean,sd^2)} \cr
#' \code{E[y]=mean}  \cr
#' \code{Var[y]=sd^2}
#' \item Negative binomial \cr
#' \code{y ~ NB(mu,phi)} \cr
#' \code{E[y]=mu}  \cr
#' \code{Var[y]=mu+phi*mu^2}
#' \item Zero-inflated negative binomial \cr
#' \code{y ~ ZINB(mu,phi,pi)} \cr
#' \code{E[y]=(1-pi)*mu}
#' }
#' 
#' @return
#' This function returns the parameter estimates, the posterior probabilities,
#' and the likelihood.
#' 
#' \item{posterior}{probability of belonging to class 1\strong{:}
#' numeric vector of length \code{n}}
#' \item{converge}{path of the log-likelihood\strong{:}
#' numeric vector with maximum length \code{it.em}}
#' \item{estim0}{parameter estimates under \code{H0}\strong{:}
#' data frame}
#' \item{estim1}{parameter estimates under \code{H1}\strong{:}
#' data frame}
#' \item{loglik0}{log-likelihood under \code{H0}\strong{:}
#' numeric}
#' \item{loglik1}{log-likelihood under \code{H1}\strong{:}
#' numeric}
#' \item{lrts}{likelihood-ratio test statistic\strong{:}
#' positive numeric}
#' 
#' @seealso
#' This is an \code{\link{internal}} function.
#' The user functions are \code{\link{mixtura}} and \code{\link{scrutor}}.
#' 
#' @examples
#' # data simulation
#' n <- 100
#' z <- rep(0:1,each=n/2)
#' y <- rnorm(n=n,mean=2*z,sd=1)
#' z[(n/4):n] <- NA
#' 
#' # model fitting
#' fit.wrap(y,z,dist="norm")
#' 
fit.wrap <- function(y,z,dist,phi,pi,gamma,starts=1,it.em=100,epsilon=1e-04){
    if(starts==1){
        
    #--- single starting point -------------------------------------------------
        
        if(dist=="norm"){
            semisup::fit.norm(y=y,z=z,
                        it.em=it.em,epsilon=epsilon)
        } else if(dist=="nbinom"){
            semisup::fit.nbinom(y=y,z=z,phi=phi,gamma=gamma,
                        it.em=it.em,epsilon=epsilon)
        } else if(dist=="zinb"){
            semisup::fit.zinb(y=y,z=z,phi=phi,pi=pi,gamma=gamma,
                        it.em=it.em,epsilon=epsilon)
        }
    } else {
        
    #--- multiple starting points ----------------------------------------------
        
        if(dist=="norm"){
            trial <- lapply(seq_len(starts),FUN=function(x)
                semisup::fit.norm(y=y,z=z,it.em=it.em,epsilon=epsilon))
        } else if(dist=="nbinom"){
            trial <- lapply(seq_len(starts),FUN=function(x)
                semisup::fit.nbinom(y=y,z=z,phi=phi,gamma=gamma,
                            it.em=it.em,epsilon=epsilon))
        } else if(dist=="zinb"){
            trial <- lapply(seq_len(starts),FUN=function(x)
                semisup::fit.zinb(y=y,z=z,phi=phi,pi=pi,gamma=gamma,
                            it.em=it.em,epsilon=epsilon))
        }
        likelihood <- sapply(trial,function(x) x$loglik1)
        id <- order(likelihood,decreasing=TRUE)[1]
        trial[[id]]
    }
}

#' @title
#' Internal function
#' 
#' @export
#' @keywords internal
#'
#' @description
#' This function resamples the data,
#' fits the semi-supervised mixture model,
#' and returns the likelihood ratio test statistic.
#' It is called by \code{\link{mixtura}}.
#' 
#' @inheritParams arguments
#' 
#' @return
#' This function returns a numeric.
#' 
#' @seealso
#' This is an \code{\link{internal}} function.
#' The user functions are \code{\link{mixtura}} and \code{\link{scrutor}}.
#' 
#' @examples
#' # data simulation
#' n <- 100
#' z <- rep(0:1,each=n/2)
#' y <- rnorm(n=n,mean=2*z,sd=1)
#' z[(n/4):n] <- NA
#' 
#' # observed test statistic
#' fit.wrap(y=y,z=z,dist="norm")$lrts
#' 
#' # simulated test statistic
#' resam.lrts(y=y,z=z,dist="norm",
#'            phi=NULL,pi=NULL,gamma=NULL,
#'            test="perm",pass=NULL)
#' 
resam.lrts <- function(y,z,dist,phi,pi,gamma,test,pass,...){
    
    if(test=="perm"){
        
    #--- permutation -----------------------------------------------------------
        
        order <- sample(seq_along(y))
        y_sim <- y[order]
        gamma_sim <- gamma[order]
        semisup::fit.wrap(y=y_sim,z=z,dist=dist,
                phi=phi,pi=pi,gamma=gamma_sim,...)$lrts
        
    } else if(test=="boot") {
        
    #--- parametric bootstrap --------------------------------------------------
     
        if(dist=="norm"){
            y_sim <- stats::rnorm(pass$n,mean=pass$mean,sd=pass$sd)
        } else if(dist=="nbinom"){
            y_sim <- stats::rnbinom(pass$n,mu=gamma*pass$mu,size=1/phi)
        } else if(dist=="zinb"){
            y_sim <- stats::rnbinom(pass$n,mu=gamma*pass$mu,size=1/phi)
            zero <- as.logical(stats::rbinom(pass$n,size=1,prob=pi))
            y_sim[zero] <- 0
        }
        y_sim[pass$keep] <- y[pass$keep]
        semisup::fit.wrap(y=y_sim,z=z,dist=dist,
                phi=phi,pi=pi,gamma=gamma,...)$lrts
    }
}

#-------------------------------------------------------------------------------
#--- Internal: Gaussian --------------------------------------------------------
#-------------------------------------------------------------------------------

#' @title
#' Internal function
#'
#' @description
#' This function fits the semi-supervised Gaussian mixture model.
#' It is called by \code{\link{fit.wrap}}.
#' 
#' @export
#' @keywords internal
#' @inheritParams arguments
#' 
#' @return
#' This function returns
#' the parameter estimates,
#' the posterior probabilities,
#' and the likelihood.
#' 
#' @seealso
#' This is an \code{\link{internal}} function.
#' The user functions are \code{\link{mixtura}} and \code{\link{scrutor}}.
#'
#' @examples
#' # data simulation
#' n <- 100
#' z <- rep(0:1,each=n/2)
#' y <- rnorm(n=n,mean=2*z,sd=1)
#' z[(n/4):n] <- NA
#' 
#' # model fitting
#' fit.norm(y,z,it.em=100,epsilon=1e-04)
#' 
fit.norm <- function(y,z,it.em,epsilon){
    x <- 1*is.na(z)
    n <- length(y)
    
    # variance penalty
    z1 <- z
    z1[x==1] <- 0.5
    z0 <- 1 - z1
    s0 <- sqrt(sum(z0*(y-sum(z0*y)/sum(z0))^2)/sum(z0))
    s1 <- sqrt(sum(z1*(y-sum(z1*y)/sum(z1))^2)/sum(z1))
    
    # null likelihood
    z1[x==1] <- 0
    z0 <- 1 - z1
    
    loglik.alt <- rep(NA,times=it.em)
    for(i in c(0,seq_len(it.em))){
        
        #--- maximisation step -------------------------------------------------
        
        mean0 <- sum(z0*y)/sum(z0)
        mean1 <- sum(z1*y)/sum(z1)
        
        sd0 <- sqrt((2*s0^2+sum(z0*(y-mean0)^2))/(2+sum(z0)))
        sd1 <- sqrt((2*s1^2+sum(z1*(y-mean1)^2))/(2+sum(z1)))
        
        p0 <- sum(z0*x)/sum(x)
        if(is.na(p0)){p0 <- 0.5} # check whether necessary
        p1 <- 1 - p0
        
        #--- expectation step --------------------------------------------------
        
        d0 <- 1/sqrt(2*pi*sd0^2)*exp(-(y-mean0)^2/(2*sd0^2))
        d1 <- 1/sqrt(2*pi*sd1^2)*exp(-(y-mean1)^2/(2*sd1^2))
        
        logd0 <- log(d0)
        logd1 <- log(d1)
        
        #--- convergence -------------------------------------------------------
        
        if(i==0){
            # null estimates
            estim.null <- data.frame(p0=p0,mean0=mean0,sd0=sd0,
                                    p1=p1,mean1=mean1,sd1=sd1)
            # null likelihood
            logd0[!is.finite(logd0)] <- -99e99
            logd1[!is.finite(logd1)] <- -99e99
            pen0 <- 1 - log(sd0^2/s0^2) - s0^2/sd0^2
            pen1 <- 1 - log(sd1^2/s1^2) - s1^2/sd1^2
            loglik.null <- sum(z0*logd0+z1*logd1)+sum(c(pen0,pen1),na.rm=TRUE)
            # initialise membership probabilities
            z1[x==1] <- stats::runif(n=n,min=0,max=1)[x==1]
            z0 <- 1 - z1
        } else {
            # calculate membership probabilities
            z0 <- (1-x)*z0 + x*(p0*d0)/(p0*d0+p1*d1)
            z1 <- 1 - z0
            # alternative likelihood
            logd0[logd0==-Inf] <- -99e99
            logd1[logd1==-Inf] <- -99e99
            loglik.alt[i] <- 
                sum((1-x)*z0*logd0 + (1-x)*z1*logd1 + x*log(p0*d0+p1*d1)) + 
                2 - log(sd0^2/s0^2) - s0^2/sd0^2 - log(sd1^2/s1^2) - s1^2/sd1^2
            # convergence
            if(i>1){if(loglik.alt[i]<loglik.alt[i-1]+epsilon){break}}
        }
    }
    
    # auxiliary estimates: s0=s0, s1=s1
    estim.alt <- data.frame(p0=p0,mean0=mean0,sd0=sd0,p1=p1,mean1=mean1,sd1=sd1)
    list(posterior=z1,converge=loglik.alt[seq_len(i)],
        estim0=estim.null,estim1=estim.alt,
        loglik0=loglik.null,loglik1=loglik.alt[i],
        lrts=2*(loglik.alt[i]-loglik.null))
}


#-------------------------------------------------------------------------------
#--- Internal: negative binomial -----------------------------------------------
#-------------------------------------------------------------------------------

#' @title
#' Internal function
#'
#' @description
#' This function fits the semi-supervised negative binomial mixture model.
#' It is called by \code{\link{fit.wrap}}.
#' 
#' @export
#' @keywords internal
#' @inheritParams arguments
#' 
#' @return
#' This function returns
#' the parameter estimates,
#' the posterior probabilities,
#' and the likelihood.
#' 
#' @seealso
#' This is an \code{\link{internal}} function.
#' The user functions are \code{\link{mixtura}} and \code{\link{scrutor}}.
#'
#' @examples
#' # data simulation
#' n <- 100
#' z <- rep(0:1,each=n/2)
#' gamma <- runif(n=n,min=0,max=2)
#' y <- rnbinom(n=n,mu=gamma*(5+2*z),size=1/0.05)
#' z[(n/4):n] <- NA
#' 
#' # model fitting
#' fit.nbinom(y,z,phi=0.05,gamma=gamma,
#' it.em=100,epsilon=1e-04)
#' 
fit.nbinom <- function(y,z,phi,gamma,it.em,epsilon){
    
    x <- 1*is.na(z)
    n <- length(y)
    
    # null likelihood
    z1 <- z
    z1[x==1] <- 0
    z0 <- 1 - z1
    
    loglik.alt <- rep(NA,times=it.em)
    for(i in c(0,seq_len(it.em))){
        
        #--- maximisation step -------------------------------------------------
        
        mu0 <- sum(z0*y/gamma)/sum(z0)
        mu1 <- sum(z1*y/gamma)/sum(z1)
        
        p0 <- sum(z0*x)/sum(x)
        if(is.na(p0)){p0 <- 0.5}
        p1 <- 1 - p0
        
        #--- expectation step --------------------------------------------------
        
        d0 <- exp(lgamma(y+1/phi) - lgamma(1/phi) - lgamma(y+1) +
                    (1/phi)*log(1/(1+gamma*mu0*phi)) + 
                    y*log(gamma*mu0/(1/phi+gamma*mu0)))
        d1 <- exp(lgamma(y+1/phi) - lgamma(1/phi) - lgamma(y+1) +
                    (1/phi)*log(1/(1+gamma*mu1*phi)) + 
                    y*log(gamma*mu1/(1/phi+gamma*mu1)))
        
        if(i==0){
            if(!is.na(mu0)){if(mu0==0){d0[y==0] <- 1}}
            if(!is.na(mu1)){if(mu1==0){d1[y==0] <- 1}}
        } else {
            if(mu0==0){d0[y==0] <- 1}
            if(mu1==0){d1[y==0] <- 1}
        }
        
        logd0 <- log(d0)
        logd1 <- log(d1)
        
        #--- convergence -------------------------------------------------------
        
        if(i==0){
            # null estimates
            estim.null <- data.frame(p0=p0,mu0=mu0,
                                    p1=p1,mu1=mu1,phi=phi)
            # null likelihood
            logd0[!is.finite(logd0)] <- -99e99
            logd1[!is.finite(logd1)] <- -99e99
            loglik.null <- sum(z0*logd0+z1*logd1)
            # initialise membership probabilities
            z1[x==1] <- stats::runif(n=n,min=0,max=1)[x==1]
            z0 <- 1 - z1
        } else {
            # calculate membership probabilities
            d0[d0==0] <- 1e-99
            d1[d1==0] <- 1e-99
            z0 <- (1-x)*z0 + x*(p0*d0)/(p0*d0+p1*d1)
            z1 <- 1 - z0
            if(sum(z0)==0|sum(z1)==0){
                z1[x==1] <- stats::runif(n=n,min=0,max=1)[x==1]
                z0 <- 1 - z1
            }
            # alternative likelihood
            logd0[logd0==-Inf] <- -99e99
            logd1[logd1==-Inf] <- -99e99
            loglik.alt[i] <- 
                sum((1-x)*z0*logd0 + (1-x)*z1*logd1 + x*log(p0*d0+p1*d1))
            if(i>1){if(loglik.alt[i]<loglik.alt[i-1]+epsilon){break}}
        }
    }
    
    estim.alt <- data.frame(p0=p0,mu0=mu0,p1=p1,mu1=mu1,phi=phi)
    list(posterior=z1,converge=loglik.alt[seq_len(i)],estim0=estim.null,
        estim1=estim.alt,loglik0=loglik.null,loglik1=loglik.alt[i],
        lrts=2*(loglik.alt[i]-loglik.null))
}

#' @title
#' Internal function
#'
#' @description
#' This function fits the semi-supervised zero-inflated
#' negative binomial mixture model.
#' It is called by \code{\link{fit.wrap}}.
#' 
#' @export
#' @keywords internal
#' @inheritParams arguments
#' 
#' @return
#' This function returns
#' the parameter estimates,
#' the posterior probabilities,
#' and the likelihood.
#' 
#' @seealso
#' This is an \code{\link{internal}} function.
#' The user functions are \code{\link{mixtura}} and \code{\link{scrutor}}.
#'
#' @examples
#' # data simulation
#' n <- 100
#' z <- rep(0:1,each=n/2)
#' gamma <- runif(n=n,min=0,max=2)
#' y <- rnbinom(n=n,mu=gamma*(5+2*z),size=1/0.05)
#' y[sample(1:n,size=0.2*n)] <- 0
#' z[(n/4):n] <- NA
#' 
#' # model fitting
#' fit.zinb(y,z,phi=0.05,pi=0.2,gamma=gamma,
#' it.em=100,epsilon=1e-04)
#' 
fit.zinb <- function(y,z,phi,pi,gamma,it.em,epsilon){
    
    x <- 1*is.na(z)
    n <- length(y)
    
    # null likelihood
    z1 <- z
    z1[x==1] <- 0
    z0 <- 1 - z1
    
    loglik1 <- rep(NA,times=it.em)
    for(i in c(0,seq_len(it.em))){
        
        #--- maximisation step -------------------------------------------------
        
        mu0 <- sum(z0*y/((1-pi)*gamma))/sum(z0)
        mu1 <- sum(z1*y/((1-pi)*gamma))/sum(z1)
        
        p0 <- sum(z0*x)/sum(x)
        if(is.na(p0)){p0 <- 0.5}
        p1 <- 1 - p0
        
        #--- expectation step --------------------------------------------------
        
        d0 <- pi*(y<=0.001) +
                (1-pi)*exp(lgamma(y+1/phi) - lgamma(1/phi) - lgamma(y+1) + 
                (1/phi)*log(1/(1+gamma*mu0*phi)) + 
                y*log(gamma*mu0/(1/phi+gamma*mu0)))
        d1 <- pi*(y<=0.001) +
                (1-pi)*exp(lgamma(y+1/phi) - lgamma(1/phi) - lgamma(y+1) + 
                (1/phi)*log(1/(1+gamma*mu1*phi)) + 
                y*log(gamma*mu1/(1/phi+gamma*mu1)))
        
        if(!is.na(mu0)){if(mu0==0){d0[y==0] <- 1}}
        if(!is.na(mu1)){if(mu1==0){d1[y==0] <- 1}}
        
        logd0 <- log(d0)
        logd1 <- log(d1)

        #--- convergence -------------------------------------------------------
        
        if(i==0){
            # null estimates
            estim0 <- data.frame(p0=p0,mu0=mu0,p1=p1,mu1=mu1,phi=phi,pi=pi)
            # null likelihood
            logd0[!is.finite(logd0)] <- -99e99
            logd1[!is.finite(logd1)] <- -99e99
            loglik0 <- sum(z0*logd0+z1*logd1)
            # initialise membership probabilities
            z1[x==1] <- stats::runif(n=n,min=0,max=1)[x==1]
            z0 <- 1 - z1
        } else {
            # calculate membership probabilities
            d0[d0==0] <- 1e-99
            d1[d1==0] <- 1e-99
            z0 <- (1-x)*z0 + x*(p0*d0)/(p0*d0+p1*d1)
            z1 <- 1 - z0
            if(sum(z0)==0|sum(z1)==0){
                z1[x==1] <- stats::runif(n=n,min=0,max=1)[x==1]
                z0 <- 1 - z1
            }
            # alternative likelihood
            logd0[logd0==-Inf] <- -99e99
            logd1[logd1==-Inf] <- -99e99
            loglik1[i] <- 
                sum((1-x)*z0*logd0 + (1-x)*z1*logd1 + x*log(p0*d0+p1*d1))
            if(i>1){if(loglik1[i]<loglik1[i-1]+epsilon){break}}
        }
    }
    
    estim1 <- data.frame(p0=p0,mu0=mu0,p1=p1,mu1=mu1,phi=phi,pi=pi)
    
    list(posterior=z1,converge=loglik1[seq_len(i)],estim0=estim0,
         estim1=estim1,loglik0=loglik0,loglik1=loglik1[i],
         lrts=2*(loglik1[i]-loglik0))
}


#-------------------------------------------------------------------------------
#--- Internal: Estimation of fixed parameters ----------------------------------
#-------------------------------------------------------------------------------

#' @title
#' Internal function
#'
#' @description
#' These functions estimate the parameters of the
#' (zero-inflated) negative binomial distribution
#' by applying the maximum likelihood method
#' to the labelled observations in class 0.
#' 
#' @export
#' @keywords internal
#' @inheritParams arguments
#' 
#' @return
#' These functions return a list of numerics.
#' 
#' @seealso
#' These are \code{\link{internal}} functions.
#' The user functions are \code{\link{mixtura}} and \code{\link{scrutor}}.
#'
#' @examples
#' # data simulation
#' n <- 100
#' y <- stats::rnbinom(n=n,mu=5,size=1/0.05)
#' y[sample(1:n,size=0.2*n)] <- 0
#' z <- rep(0,times=n)
#' gamma <- rep(1,times=n)
#' 
#' # parameter estimation
#' estim.nbinom(y,z,gamma)
#' estim.zinb(y,z,gamma)
#' 
estim.nbinom <- function(y,z,gamma){
    y <- y[which(z==0)]
    gamma <- gamma[which(z==0)]
    mu <- sum(y)/length(y)
    loglik <- function(phi) sum(lgamma(y + 1/phi) - 
                                lgamma(1/phi) - 
                                lgamma(y + 1) - 
                                1/phi * log(1 + gamma * mu * phi) + 
                                y * log(gamma * mu) - 
                                y * log(1/phi + gamma * mu))
    phi <- suppressWarnings(stats::optimize(loglik,interval=c(0,1000),
                                            tol=10^{-10},maximum=TRUE)$maximum)
    list(mu=mu,phi=phi)
}
#'@rdname estim.nbinom
#'@export
#'@keywords internal
estim.zinb <- function(y,z,gamma){
    y <- y[which(z==0)]
    gamma <- gamma[which(z==0)]
    
    # log-likelihood
    fn <- function(par){
        -sum(VGAM::dzinegbin(x=y,munb=gamma*par[1],size=1/par[2],
                             pstr0=par[3],log=TRUE))
    }
    
    # initial values
    par <- function(y){
        mu <- max(0.01,mean(y[y!=0]),na.rm=TRUE)
        sd <- max(0.01,sd(y[y!=0]),na.rm=TRUE)
        phi <- max(0.01,(sd^2-mu)/mu^2)
        pi <- min(0.99,max(0.01,mean(y==0)))
        c(mu,phi,pi)
    }
    
    # minimisation
    epsilon <- 1e-06
    par <- stats::optim(par = par(y),fn = fn,
        lower = c(epsilon,epsilon,epsilon),
        upper = c(99e+99,99e+99,1-epsilon),
        method="L-BFGS-B"
        )$par  
    
    list(mu = par[1], phi = par[2],pi = par[3])
}

#-------------------------------------------------------------------------------
#--- Internal: debugging -------------------------------------------------------
#-------------------------------------------------------------------------------


#' @title
#' Internal function
#'
#' @description
#' This function verifies whether the arguments fulfill some formal requirements.
#'
#' @keywords internal
#' @inheritParams arguments
#'
#' @details
#' If one or more entries of \code{z} are equal to \code{1},
#' the mixture model can be fitted but not tested.
#' Accordingly, \code{kind} is replaced by \code{NULL}.
#'
#' Resampling-based testing cannot reach \code{p}-values below \code{1/iter}.
#' If \code{kind} is smaller than \code{1/iter}, it is replaced by \code{0}.
#'
#' @return
#' This function returns warnings and errors.
#' It also returns \code{kind} (see details).
#'
#' @seealso
#' This is an \code{\link{internal}} function.
#' The user functions are \code{\link{mixtura}} and \code{\link{scrutor}}.
#'
#' @examples
#' NULL
#'
debug <- function(y,z,dist,phi,pi,gamma,test,iter,kind,...){
    
    dots <- list(...)
    starts <- dots$starts
    it.em <- dots$it.em
    epsilon <- dots$epsilon
    
    Y <- as.matrix(y)
    Z <- as.matrix(z)
    
    #--- observations ----------------------------------------------------------
    
    if(!is.numeric(Y) & !is.integer(Y)){
        stop("Argument \"y\" must be of type \"numeric\".",call.=FALSE)
    }
    
    if(any(!is.finite(Y))){
        stop("Argument \"y\" must be finite.",call.=FALSE)
    }
    
    if(any((Y < 0) & (dist=="nbinom"))){
        stop("Argument \"y\" must be non-negative (\"nbinom\").",call.=FALSE)
    }
    
    #--- labels ----------------------------------------------------------------
    
    if(!is.numeric(Z) & !is.integer(Z) & !is.logical(Z)){
        stop("Argument \"z\" must be of type \"numeric\".",call.=FALSE)
    }
    
    if(nrow(Z)!=nrow(Y)){
        stop("Arguments \"y\" and \"z\" must have the same sample size.",
             call.=FALSE)
    }
    
    if(!all(Z %in% c(0,1,NA))){
        stop("Argument \"z\" must equal 0 or NA.",call.=FALSE)
    }
    
    if((1 %in% Z) & !is.null(test)){
        stop("Set \"test\" to NULL (because \"z\" contains 1).",call.=FALSE)
    }
    
    zeros <- apply(Z,2,function(z) sum(z==0,na.rm=TRUE))
    ones <- apply(Z,2,function(z) sum(z==1,na.rm=TRUE))
    nolab <- apply(Z,2,function(z) sum(is.na(z)))
    if(any(zeros+nolab < 2) | any(ones+nolab < 2)){
        stop("Argument \"Z\" must allow two observations per class.",
             call.=FALSE)
    }
    
    #--- distribution ----------------------------------------------------------
    
    if(!is.character(dist) & !is.factor(dist)){
        stop("Argument \"dist\" must of type \"character\".",call.=FALSE)
    }
    
    if(length(dist)!=1){
        stop("Argument \"dist\" must have length 1.",call.=FALSE)
    }
    
    if(!dist %in% c("norm","nbinom","zinb")){
        stop("Argument \"dist\" must equal \"norm\", \"nbinom\", or \"zinb\".",
             call.=FALSE)
    }
    
    if(dist=="zinb" & (is.null(phi) != is.null(pi))){
        stop("Provide both or none of \"phi\" and \"pi\".",
             call.=FALSE)
    }
    
    #--- dispersion ------------------------------------------------------------
    
    if(!is.null(phi)){
        if(dist=="norm"){
            warning("Ignoring \"phi\" (because \"dist\" is \"norm\").",
                    call.=FALSE)
        } else {
            if(!is.numeric(phi)){
                stop("Argument \"phi\" must be of type \"numeric\".",
                     call.=FALSE)
            }
            if(length(phi)!=ncol(Y)){
                stop("Argument \"phi\" must have length q.",call.=FALSE)
            }
            if(any(!is.finite(phi))){
                stop("Argument \"phi\" must be finite.",call.=FALSE)
            }
            if(any(phi<=0)){
                stop("Argument \"phi\" must be positive.",call.=FALSE)
            }
        }
    }
    
    #--- zero-inflation --------------------------------------------------------
    
    if(!is.null(pi)){
        if(dist!="zinb"){
            warning("Ignoring \"pi\" (because \"dist\" is not \"zinb\").",
                call.=FALSE)
        } else {
            if(!is.numeric(pi)){
                stop("Argument \"pi\" must be of type \"numeric\".",
                     call.=FALSE)
            }
            if(length(pi)!=ncol(Y)){
                stop("Argument \"pi\" must have length p.",call.=FALSE)
            }
            if(any(!is.finite(pi))){
                stop("Argument \"pi\" must be finite.",call.=FALSE)
            }
            if(any(pi < 0 | pi > 1)){
                stop("Argument \"pi\" must be between 0 and 1.",call.=FALSE)
            }
        }
    }
    
    #--- offset ----------------------------------------------------------------
    
    if(!is.null(gamma)){
        if(dist=="norm"){
            warning("Ignoring \"gamma\" (because \"dist\" is \"norm\").",
                call.=FALSE)
        } else {
            if(!is.numeric(gamma)){
                stop("Argument \"gamma\" must be of type \"numeric\".",
                     call.=FALSE)
            }
            if(length(gamma)!=nrow(Y)){
                stop("Argument \"gamma\" must have length n.",call.=FALSE)
            }
            if(any(!is.finite(gamma))){
                stop("Argument \"gamma\" must be finite.",call.=FALSE)
            }
            if(any(gamma <= 0)){
                stop("Argument \"gamma\" must be positive.",call.=FALSE)
            }
        }
    }
    
    #--- test ------------------------------------------------------------------
    
    if(!is.null(test)){
        if(!is.character(test) & !is.factor(test)){
            stop("Argument \"test\" must of type \"character\".",call.=FALSE)
        }
        
        if(length(test)!=1){
            stop("Argument \"test\" must have length 1.",call.=FALSE)
        }
        
        if(!test %in% c("boot","perm")){
            stop("Argument \"test\" must equal \"boot\" or \"perm\".",
                 call.=FALSE)
        }
    }
    
    
    #--- iter ------------------------------------------------------------------
    
    if(!is.null(iter)){
        if(!is.integer(iter) & !is.numeric(iter)){
            stop("Argument \"iter\" must be of type \"integer\".",call.=FALSE)
        }
        if(length(iter)!=1){
            stop("Argument \"iter\" must have length 1.",call.=FALSE)
        }
        if(!is.finite(iter)){
            stop("Argument \"iter\" must be finite.",call.=FALSE)
        }
        if(iter!=round(iter)){
            stop("Argument \"iter\" must be an integer.",call.=FALSE)
        }
        if(iter<1){
            stop("Argument \"iter\" must be positive.",call.=FALSE)
        }
    }
    
    #--- kind ------------------------------------------------------------------
    
    if(!is.null(kind)){
        if(!is.numeric(kind)){
            stop("Argument \"kind\" must be of type \"numeric\".",call.=FALSE)
        }
        if(length(kind)!=1){
            stop("Argument \"kind\" must have length 1.",call.=FALSE)
        }
        if(!is.finite(kind)){
            stop("Argument \"kind\" must be finite.",call.=FALSE)
        }
        if((kind < 0) | (kind > 1)){
            stop("Argument \"kind\" must be between 0 and 1.",call.=FALSE)
        }
        if(kind==0){
            stop("Not yet implemented: \"kind=0\".",call.=FALSE)
        } else {
            if(iter < (1/kind)){
                stop("Set \"iter\" to 1/kind or larger.",call.=FALSE)
            }
        }
    }
    
    #--- it.em -----------------------------------------------------------------
    
    if(!is.null(it.em)){
        if(!is.integer(it.em) & !is.numeric(it.em)){
            stop("Argument \"it.em\" must be of type \"integer\".",call.=FALSE)
        }
        if(length(it.em)!=1){
            stop("Argument \"it.em\" must have length 1.",call.=FALSE)
        }
        if(!is.finite(it.em)){
            stop("Argument \"it.em\" must be finite.",call.=FALSE)
        }
        if(it.em!=round(it.em)){
            stop("Argument \"it.em\" must be an integer.",call.=FALSE)
        }
        if(it.em<1){
            stop("Argument \"it.em\" must be positive.",call.=FALSE)
        }
    }
    
    #--- epsilon ---------------------------------------------------------------
    
    if(!is.null(epsilon)){
        if(!is.numeric(epsilon)){
            stop("Argument \"epsilon\" must be of type \"numeric\".",call.=FALSE)
        }
        if(length(it.em)!=1){
            stop("Argument \"epsilon\" must have length 1.",call.=FALSE)
        }
        if(!is.finite(it.em)){
            stop("Argument \"epsilon\" must be finite.",call.=FALSE)
        }
        if(it.em<0){
            stop("Argument \"epsilon\" must be non-negative.",call.=FALSE)
        }
    }
    
    #--- starts ----------------------------------------------------------------
    
    if(!is.null(starts)){
        if(!is.integer(starts) & !is.numeric(starts)){
            stop("Argument \"starts\" must be of type \"integer\".",call.=FALSE)
        }
        if(length(starts)!=1){
            stop("Argument \"starts\" must have length 1.",call.=FALSE)
        }
        if(!is.finite(starts)){
            stop("Argument \"starts\" must be finite.",call.=FALSE)
        }
        if(starts!=round(starts)){
            stop("Argument \"starts\" must be an integer.",call.=FALSE)
        }
        if(starts<1){
            stop("Argument \"starts\" must be positive.",call.=FALSE)
        }
    }
}

