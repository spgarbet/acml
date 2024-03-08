

IndividualScores.C <- function(params, y, x, z, id, w.function, cutpoints, SampProb, SampProbi, ProfileCol=NA)
{
    .Call("lldetail",
                params,y, x, z, id*1.0, w.function, 
                cutpoints, SampProb, SampProbi)
}

LogLikeAndScore.C <- function(params, y, x, z, id, w.function, cutpoints, SampProb, SampProbi, ProfileCol=NA)
{
    ans <- .Call("llscore",
                    params,y, x, z, id*1.0, w.function, 
                    cutpoints, SampProb, SampProbi)

     ## Force the gradient of the fixed parameter to be zero, so that it does not move
     if (!is.na(ProfileCol))
     {
        attr(ans, "gradient")[ProfileCol] <- 0
     }
     
     ans
}

## If you do not want to use the ascertainment correction term in the conditional likelihood
## set all SampProb values equal to each other.  This would be the case if you were doing 
## straightforward maximum likelihood (albeit inefficient) or weighted likelihood.
ACML.LME.C <- function( y,                            ## response vector
                        x,                            ## fixed effects design matrix
                        z,                            ## random effects design matrix (right now this should be an intercept and a time-varying covariate (?)
                        id,                           ## subject id variable
                        w.function="mean",            ## Function upon which sampling is based. Choices are the univariate "mean", "intercept", "slope", and the bivariate "bivar"
                        InitVals,                     ## Starting values 
                        cutpoints = c(0,5),           ## the cutpoints: when w.function="bivar", this is a vector of length 4 that define a central, rectangular region with vertices (x_lower, x_upper, y_lower, y_upper).
                        SampProb  = c(1, 1, 1),       ## Sampling probabilities within the sampling strata to be used for ACML
                        SampProbi = rep(1, length(y)),## Subject specific sampling probabilities to only be used if doing IPWL.  Note if doing IPWL, only use robcov (robust variances) and not covar
                        ProfileCol= NA,               ## Columns to be held fixed while doing profile likelihood.  It is fixed at its initial value.
                        print.level=0)
{

    out <- nlm(LogLikeAndScore.C, InitVals, y=y, x=x, z=z, id=id, w.function=w.function, cutpoints=cutpoints, SampProb=SampProb,
               SampProbi=SampProbi, ProfileCol=ProfileCol, gradtol=1e-12,
               stepmax=4, iterlim=250, check.analyticals = TRUE, print.level=print.level)

    ## Calculate the observed information and then invert to get the covariance matrix       
    npar <- length(out$estimate)
    Hessian.eps <- 1e-7
    eps.mtx     <- diag(rep(Hessian.eps, npar))
    grad.at.max <- out$gradient
    ObsInfo     <- matrix(NA, npar, npar)
  
    ## Observed Information
    for (j in 1:npar){
        temp        <- LogLikeAndScore.C(
                            out$estimate+eps.mtx[j,], y=y, x=x, z=z, id=id, 
                            w.function=w.function, cutpoints=cutpoints,
                            SampProb=SampProb,SampProbi=SampProbi, ProfileCol=ProfileCol)
        ObsInfo[j,] <- (attr(temp,"gradient")-grad.at.max)/(Hessian.eps)
    }
    ## Cheese part of the sandwich estimator
    Cheese <- gradient.nll.lme( y=y, x=x, z=z,
                                w.function  = w.function, 
                                id          = id,
                                beta        = out$estimate[c(1:(npar-4))], 
                                sigma0      = exp(out$estimate[(npar-3)]), 
                                sigma1      = exp(out$estimate[(npar-2)]), 
                                rho         = (exp(out$estimate[(npar-1)])-1) / (exp(out$estimate[(npar-1)])+1), 
                                sigmae      = exp(out$estimate[npar]), 
                                cutpoints   = cutpoints, 
                                SampProb    = SampProb, 
                                SampProbi   = SampProbi,
                                CheeseCalc  = TRUE)

    if (!is.na(ProfileCol))
    {
        out$estimate    <- out$estimate[-ProfileCol]
        ObsInfo         <- ObsInfo[-ProfileCol, -ProfileCol]
        Cheese          <- Cheese[-ProfileCol, -ProfileCol]
    }

    list(Ests   = out$estimate,
         covar  = solve(ObsInfo),
         LogL   = -out$minimum,
         Code   = out$code,
         robcov = solve(ObsInfo)%*%Cheese%*%solve(ObsInfo) )
}
