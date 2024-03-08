##
## Notes
## 1) When doing intercept, slope, and bivar based sampling, we are assuming we are sampling based on the variables that determine the random effects (Z)
##
  

#library(MASS, lib.loc="/usr/local/biostat_r/lib/R")
#library(mvtnorm, lib.loc="/usr/local/biostat_r/lib/R")
library(MASS)
library(mvtnorm)
#library(lme4, lib.loc="/usr/local/biostat_r/lib/R")

vi.calc <- function(zi, sigma0, sigma1, rho, sigmae){
    zi %*% matrix(c(sigma0^2, rho*sigma0*sigma1,rho*sigma0*sigma1,sigma1^2), nrow=2) %*% t(zi) + 
    sigmae*sigmae*diag(length(zi[,2]))
}
## Ascertainment correction piece for univariate sampling
lci <- function(cutpoints, SampProb, mu_q, sigma_q){ 
	CDFs <- pnorm(c(-Inf, cutpoints, Inf), mu_q, sigma_q)
        sum( SampProb*(CDFs[2:length(CDFs)] - CDFs[1:(length(CDFs)-1)]) )
}
## Ascertainment correction piece for bivariate sampling
lci.bivar <- function(cutpoints, SampProb, mu_q, sigma_q){ 
    (SampProb[1]-SampProb[2])*pmvnorm(lower=c(cutpoints[c(1,3)]), upper=c(cutpoints[c(2,4)]), mean=mu_q, sigma=sigma_q)[[1]] + SampProb[2]
    #print("blah2")
}
subject.ll.lme <- function(yi, xi, beta, vi){ 
	resid <- yi - xi %*% beta
    -(1/2) * (length(xi[,1])*log(2*pi) + log(det(vi)) + t(resid) %*% solve(vi) %*% resid )[1,1]
}

## Calculate log of the ascertainment correction for the univariate sampling case
ascertainment.correction <- function(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb){
    vi      <- vi.calc(zi, sigma0, sigma1, rho, sigmae)
    mu      <- xi %*% beta
    mu_q    <- (wi %*% mu)[,1]  
    sigma_q <- sqrt((wi %*% vi %*% t(wi))[1,1])
    log(lci(cutpoints, SampProb, mu_q, sigma_q))
}
## Calculate log of the ascertainment correction for the bivariate sampling case
ascertainment.correction.bivar <- function(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb){
    vi      <- vi.calc(zi, sigma0, sigma1, rho, sigmae)
    mu      <- xi %*% beta
    mu_q    <- as.vector(wi %*% mu)
    sigma_q <- wi %*% vi %*% t(wi)
    log((SampProb[1]-SampProb[2])*pmvnorm(lower=c(cutpoints[c(1,3)]), upper=c(cutpoints[c(2,4)]),
                                          mean=mu_q, sigma=sigma_q)[[1]] + SampProb[2])
}

## Calculate conditional likelihood for the univariate and bivariate sampling cases
total.nll.lme <- function(y, x, z, w.function, id, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb, SampProbi){
    total <- 0
    for(i in unique(id)){ 
        yi <- y[id==i]
        ni <- length(yi)
        xi <- x[id==i,]
        zi <- z[id==i,]
        if (w.function != "bivar"){
            if (w.function=="mean")      wi <- t(rep(1/ni, ni))
            if (w.function=="intercept") wi<- (solve(t(zi[,1:2])%*% zi[,1:2]) %*% t(zi[,1:2]))[1,]
            if (w.function=="slope")     wi<- (solve(t(zi[,1:2])%*% zi[,1:2]) %*% t(zi[,1:2]))[2,]
            wi    <- matrix(wi, 1, ni)
            IPWi  <- 1/ unique(SampProbi[id==i])
            vi    <- vi.calc(zi, sigma0, sigma1, rho, sigmae)
            total <- total + subject.ll.lme(yi, xi, beta, vi)*IPWi -
                             ascertainment.correction(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb)           
        }else{  
            wi    <- (solve(t(zi[,1:2])%*% zi[,1:2]) %*% t(zi[,1:2]))
            IPWi  <- 1/ unique(SampProbi[id==i])
            vi    <- vi.calc(zi, sigma0, sigma1, rho, sigmae)
            total <- total + subject.ll.lme(yi, xi, beta, vi)*IPWi -
                             ascertainment.correction.bivar(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb)
        }
    }    
    -total    
}

ascertainment.gradient.correction.bivar <- function(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb){
    eps     <- 1e-6
    param   <- c(beta, sigma0, sigma1, rho, sigmae)
    npar    <- length(param)
    vi      <- vi.calc(zi, sigma0, sigma1, rho, sigmae)
    mu      <- xi %*% beta
    mu_q    <- as.vector(wi %*% mu)
    sigma_q <- wi %*% vi %*% t(wi)
    start   <- pmvnorm(lower=c(cutpoints[c(1,3)]), upper=c(cutpoints[c(2,4)]), mean=mu_q, sigma=sigma_q)[[1]]

    ## for this bivariate sampling case, right now we calculate gradients numerically
    new.area <- NULL
    eps.mtx <- diag(c(rep(eps,npar))) 
    for (rr in 1:npar){
        par.new      <- param+eps.mtx[rr,]  
        vi.tmp       <- vi.calc(zi, par.new[(npar-3)], par.new[(npar-2)], par.new[(npar-1)], par.new[npar])
        mu.tmp       <- xi %*% par.new[1:(npar-4)]
        mu_q.tmp     <- as.vector(wi %*% mu.tmp)
        sigma_q.tmp  <- wi %*% vi.tmp %*% t(wi)
         new.area <- c(new.area, pmvnorm(lower=c(cutpoints[c(1,3)]), upper=c(cutpoints[c(2,4)]), mean=mu_q.tmp, sigma=sigma_q.tmp)[[1]])
    }
    Deriv <- (new.area-start)/eps
    out <- (SampProb[1]-SampProb[2])*Deriv / lci.bivar(cutpoints, SampProb, mu_q, sigma_q)
    out
}

ascertainment.gradient.correction <- function(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb){
    vi      <- vi.calc(zi, sigma0, sigma1, rho, sigmae)
    mu      <- xi %*% beta
    mu_q    <- (wi %*% mu)[,1]
    sigma_q <- sqrt((wi %*% vi %*% t(wi))[1,1])

    l <- lci(cutpoints, SampProb, mu_q, sigma_q)
    p <- SampProb[1:(length(SampProb)-1)] - SampProb[2:(length(SampProb))]
    f <- dnorm(cutpoints, mu_q, sigma_q)
  
    d_li_beta <- (wi %*% xi) * sum(p*f) / l

    f_alpha_k <- sum(p*f*(cutpoints - mu_q)) / (l * sigma_q *2 * sqrt(wi %*% vi %*% t(wi))[1,1])
    a5        <- (wi %*% zi %*% matrix(c(2*sigma0, rho*sigma1,    rho*sigma1,     0),        nrow=2) %*% t(zi) %*% t(wi))[1,1]
    a6        <- (wi %*% zi %*% matrix(c(0,        rho*sigma0,    rho*sigma0,     2*sigma1), nrow=2) %*% t(zi) %*% t(wi))[1,1]
    a7        <- (wi %*% zi %*% matrix(c(0,        sigma0*sigma1, sigma0*sigma1,  0),        nrow=2) %*% t(zi) %*% t(wi))[1,1]
    a8        <- (wi %*% (2 * sigmae * diag(length(yi))) %*% t(wi))[1,1]
    c(d_li_beta, c(f_alpha_k * c(a5, a6, a7, a8)))
}

subject.gradient.ll.lme <- function(yi, xi, zi, beta, sigma0, sigma1, rho, sigmae){
    resid <- yi - xi %*% beta
    vi    <- vi.calc(zi, sigma0, sigma1, rho, sigmae)
    inv.v <- solve(vi)
  
    dvk5 <- zi %*% matrix(c(2*sigma0, rho*sigma1, rho*sigma1, 0),nrow=2) %*% t(zi)
    dvk6 <- zi %*% matrix(c(0, rho*sigma0, rho*sigma0, 2*sigma1),nrow=2) %*% t(zi)
    dvk7 <- zi %*% matrix(c(0, sigma0*sigma1, sigma0*sigma1, 0),nrow=2) %*% t(zi)
    dvk8 <- 2 * sigmae  * diag(length(yi))
    
    l14 <- t(xi) %*% inv.v %*% resid # for Beta
    l5  <- -0.5*(sum(diag(inv.v %*% dvk5)) - t(resid) %*% inv.v %*% dvk5 %*% inv.v %*% resid)[1,1]
    l6  <- -0.5*(sum(diag(inv.v %*% dvk6)) - t(resid) %*% inv.v %*% dvk6 %*% inv.v %*% resid)[1,1]
    l7  <- -0.5*(sum(diag(inv.v %*% dvk7)) - t(resid) %*% inv.v %*% dvk7 %*% inv.v %*% resid)[1,1]
    l8  <- -0.5*(sum(diag(inv.v %*% dvk8)) - t(resid) %*% inv.v %*% dvk8 %*% inv.v %*% resid)[1,1]
    list(gr=append(l14, c(l5,l6,l7,l8)),
         vi=vi)
}

# 

gradient.nll.lme <- function(y, x, z, w.function, id, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb, SampProbi, CheeseCalc=FALSE){
    total <- 0
    param.vec <- c(beta, log(sigma0),log(sigma1),log((1+rho)/(1-rho)),log(sigmae))
    n.par <- length(param.vec)
    cheese <- matrix(0,n.par,n.par)
    for(i in unique(id)){ 
        yi <- y[id==i]
        ni <- length(yi)
        xi <- x[id==i,]
        zi <- z[id==i,]
        IPWi <- 1/ unique(SampProbi[id==i])
        if (w.function != "bivar"){
            if (w.function=="mean")      wi <- t(rep(1/ni, ni))
            if (w.function=="intercept") wi<- (solve(t(xi[,1:2])%*% xi[,1:2]) %*% t(xi[,1:2]))[1,]
            if (w.function=="slope")     wi<- (solve(t(xi[,1:2])%*% xi[,1:2]) %*% t(xi[,1:2]))[2,]
            wi   <- matrix(wi, 1, ni)
            subject <- subject.gradient.ll.lme(yi, xi, zi, beta, sigma0, sigma1, rho, sigmae)
            correct <- ascertainment.gradient.correction(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb)
            Gradi  <- subject[['gr']]*IPWi  + correct ## Gradient for ith subject
            total  <- total + Gradi
        }else{
            wi   <- (solve(t(zi[,1:2])%*% zi[,1:2]) %*% t(zi[,1:2]))
            subject <- subject.gradient.ll.lme(yi, xi, zi, beta, sigma0, sigma1, rho, sigmae)
            correct <- ascertainment.gradient.correction.bivar(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb)
            Gradi  <- subject[['gr']]*IPWi  - correct ## Gradient for ith subject: Notice the minus here versus the plus in the univariate case
            total  <- total + Gradi
        }

 #       if(CheeseCalc==TRUE) 
 #       {
 #           print(paste("subject",i, "total",total[1]))
 #       }

        if (CheeseCalc==TRUE){
           ## Need to use the chain rule: note that param,vec is on the unconstrained scale but Gradi was calculated on the constrained parameters
           Gradi[(n.par-3)] <- Gradi[(n.par-3)]*exp(param.vec[(n.par-3)])
           Gradi[(n.par-2)] <- Gradi[(n.par-2)]*exp(param.vec[(n.par-2)])
           Gradi[(n.par-1)] <- Gradi[(n.par-1)]*2*exp(param.vec[(n.par-1)])/((exp(param.vec[(n.par-1)])+1)^2)
           Gradi[n.par]     <- Gradi[n.par]  *exp(param.vec[n.par])
           cheese  <- cheese + outer(Gradi, Gradi)
       }
    }
    if (CheeseCalc==TRUE) total <- -cheese
    -total
}
   
                               
LogLikeAndScore <- function(params, y, x, z, id, w.function, cutpoints, SampProb, SampProbi, ProfileCol=NA){
     npar   <- length(params)
     beta   <- params[1:(npar-4)]
     sigma0 <- exp(params[(npar-3)])
     sigma1 <- exp(params[(npar-2)])
     rho    <- (exp(params[(npar-1)])-1) / (exp(params[(npar-1)])+1)
     sigmae <- exp(params[npar])

     out     <- total.nll.lme(   y, x, z, w.function, id, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb, SampProbi)
     GRAD    <- gradient.nll.lme(y, x, z, w.function, id, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb, SampProbi)

     ## Need to use the chain rule: note that params is on the unconstrained scale but GRAD was calculated on the constrained parameters
     GRAD[(npar-3)] <- GRAD[(npar-3)]*exp(params[(npar-3)])
     GRAD[(npar-2)] <- GRAD[(npar-2)]*exp(params[(npar-2)])
     GRAD[(npar-1)] <- GRAD[(npar-1)]*2*exp(params[(npar-1)])/((exp(params[(npar-1)])+1)^2)
     GRAD[npar]     <- GRAD[npar]*exp(params[npar])
     ## Force the gradient of the fixed parameter to be zero, so that it does not move
     #if (!is.na(ProfileCol)) GRAD[ProfileCol] <- 0
     attr(out,"gradient") <- GRAD
     
     out
  }

## If you do not want to use the ascertainment correction term in the conditional likelihood
## set all SampProb values equal to each other.  This would be the case if you were doing 
## straightforward maximum likelihood (albeit inefficient) or weighted likelihood.
ACML.LME <- function(y,    ## response vector
                     x,    ## fixed effects design matrix
                     z,    ## random effects design matrix (right now this should be an intercept and a time-varying covariate (?)
                     id,   ## subject id variable
                     w.function="mean",           ## Function upon which sampling is based. Choices are the univariate "mean", "intercept", "slope", and the bivariate "bivar"
                     InitVals,                    ## Starting values 
                     cutpoints = c(0,5),          ## the cutpoints: when w.function="bivar", this is a vector of length 4 that define a central, rectangular region with vertices (x_lower, x_upper, y_lower, y_upper).
                     SampProb = c(1, 1, 1),       ## Sampling probabilities within the sampling strata to be used for ACML
                     SampProbi=rep(1, length(y)), ## Subject specific sampling probabilities to only be used if doing IPWL.  Note if doing IPWL, only use robcov (robust variances) and not covar
                     ProfileCol=NA){              ## Columns to be held fixed while doing profile likelihood.  It is fixed at its initial value.
   
    out <- nlm(LogLikeAndScore, InitVals, y=y, x=x, z=z, id=id, w.function=w.function, cutpoints=cutpoints, SampProb=SampProb,
               SampProbi=SampProbi, ProfileCol=ProfileCol, gradtol=1e-12,
               stepmax=4, iterlim=250, check.analyticals = TRUE) #, print.level=1)

    ## Calculate the observed information and then invert to get the covariance matrix       
    npar <- length(out$estimate)
    Hessian.eps <- 1e-7
    eps.mtx     <- diag(rep(Hessian.eps, npar))
    grad.at.max <- out$gradient
    ObsInfo     <- matrix(NA, npar, npar)
  
    ## Observed Information
    for (j in 1:npar){
    	temp        <- LogLikeAndScore(out$estimate+eps.mtx[j,], y=y, x=x, z=z, id=id, 
                                       w.function=w.function, cutpoints=cutpoints, 
                                       SampProb=SampProb,SampProbi=SampProbi, ProfileCol=ProfileCol)
  	    ObsInfo[j,] <- (attr(temp,"gradient")-grad.at.max)/(Hessian.eps)
  	}
    ## Cheese part of the sandwich estimator
    Cheese <- gradient.nll.lme(y=y, x=x, z=z, w.function=w.function, 
                             id=id, beta=out$estimate[c(1:(npar-4))], 
                             sigma0=exp(out$estimate[(npar-3)]), 
                             sigma1=exp(out$estimate[(npar-2)]), 
                             rho=   (exp(out$estimate[(npar-1)])-1) / (exp(out$estimate[(npar-1)])+1), 
                             sigmae=exp(out$estimate[npar]), 
                             cutpoints=cutpoints, 
                             SampProb=SampProb, 
                             SampProbi=SampProbi,
                             CheeseCalc=TRUE)
    
    if (!is.na(ProfileCol)){ 
        out$estimate <- out$estimate[-ProfileCol]
        ObsInfo <- ObsInfo[-ProfileCol, -ProfileCol]
        Cheese  <- Cheese[-ProfileCol, -ProfileCol]
    }

    list(Ests=out$estimate, covar=solve(ObsInfo), LogL= -out$minimum,  Code=out$code, robcov=solve(ObsInfo)%*%Cheese%*%solve(ObsInfo))
}

