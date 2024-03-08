#library(MASS, lib.loc="/usr/local/biostat_r/lib/R")
library(MASS)
GenPopnData <- function(N = 1000, n = 11, beta = c(1, 0.25, 0, 0.25),
                     sig.b0 = 0.25, sig.b1 = 0.25, rho = 0, sig.e = 0.5, prob.grp = 0.1, n.low=11, n.high=11){
    id      <- rep(1:N, each=n)
    grp     <- rep(rbinom(N, 1, prob.grp), each=n)
    time    <- rep(seq(-1,1, length.out=n), N)
    cov.mat <- matrix(c(sig.b0^2, rho*sig.b0*sig.b1, rho*sig.b0*sig.b1, sig.b1^2),2,2)
    bi      <- mvrnorm(N, mu=c(0,0), Sigma=cov.mat)
    b       <- cbind(rep(bi[,1], each=n),rep(bi[,2], each=n))
    error   <- rnorm(N*n, 0, sig.e)

   
    X <- as.matrix(cbind(1, time, grp, time*grp))
    Z <- X[,c(1:2)]
    Y <- X %*% beta + Z[,1]*b[,1] + Z[,2]*b[,2] + error

    ## Induce MCAR dropout
    n.obs <- rep(sample(rep(c(n.low:n.high),2), N, replace=TRUE), each=n) 
    obs.num <- rep(c(1:n), N)
    X <- X[obs.num<= n.obs,]
    Z <- Z[obs.num<= n.obs,]
    Y <- Y[obs.num<= n.obs]
    id <- id[obs.num<= n.obs]
    list(id=id, X=X, Y=Y, Z=Z, ni=c(unlist(tapply(Y,id,length))),
         N=N, n=n, beta=beta,sig.b0=sig.b0, sig.b1=sig.b1, rho=rho, sig.e=sig.e, prob.grp=prob.grp)
}

LinRegFn <- function(data){ X <- cbind(1, data[,3])
	                          Y <- data$Y
	                          solve(t(X)%*%X) %*% t(X) %*% Y}
	                          
est.quants <- function(N, n, beta, sig.b0, sig.b1, rho, sig.e, prob.grp, quant, n.low=11, n.high=11){
    d        <- GenPopnData(N, n, beta, sig.b0, sig.b1, rho, sig.e, prob.grp)
    data.tmp <- data.frame(id=d$id, Y=d$Y, X.time=d$X[,2])
    out      <- matrix(unlist(lapply(split(data.tmp, data.tmp$id), LinRegFn)), byrow=TRUE, ncol=2)
    out      <- cbind(c(tapply(data.tmp$Y, data.tmp$id, mean)), out)
    r        <- rbind( quantile(out[,1], probs=quant), 
                       quantile(out[,2], probs=quant),
                       quantile(out[,3], probs=quant))
    r}

## Do a search to find the quantiles that correspond to the central rectangle that contains 60 and 80 percent of 
## the subject specific intercepts and slopes.  This is not necessary if slope and intercepts are independent
## but with unequal followup they were positiviely correlated.  Searches for the smallest 'rectangle' defined
## by quantiles that contains 60 and 80 percent of the data
est.bivar.lims <- function(N, n, beta, sig.b0, sig.b1, rho, sig.e, prob.grp, quants, n.low=11, n.high=11){
        d        <- GenPopnData(N, n, beta, sig.b0, sig.b1, rho, sig.e, prob.grp, n.low=n.low, n.high=n.high)
        data.tmp <- data.frame(id=d$id, Y=d$Y, X.time=d$X[,2])
        out      <- matrix(unlist(lapply(split(data.tmp, data.tmp$id), LinRegFn)), byrow=TRUE, ncol=2)
        out      <- cbind(c(tapply(data.tmp$Y, data.tmp$id, mean)), out)
        print(cor(out[,2], out[,3]))
        
        q1 <- .99
        Del <- 1
        while (Del>0.001){ q1 <- q1-.00025
                           Del <- abs(mean(out[,2] > quantile(out[,2], probs=1-q1) & 
                                           out[,2] < quantile(out[,2], probs=q1) &
                                           out[,3] > quantile(out[,3], probs=1-q1) & 
                                           out[,3] < quantile(out[,3], probs=q1)) - quants[1])
                           #print(c(q1,Del))
                           q1}
        q2 <- .99
        Del <- 1
        while (Del>0.001){ q2 <- q2-.00025
                           Del <- abs(mean(out[,2] > quantile(out[,2], probs=1-q2) & 
                                           out[,2] < quantile(out[,2], probs=q2) &
                                           out[,3] > quantile(out[,3], probs=1-q2) & 
                                           out[,3] < quantile(out[,3], probs=q2)) - quants[2])
                           q2}
       
        rbind( c( quantile(out[,2], probs=1-q1), quantile(out[,2], probs=q1), quantile(out[,3], probs=1-q1), quantile(out[,3], probs=q1)),
               c( quantile(out[,2], probs=1-q2), quantile(out[,2], probs=q2), quantile(out[,3], probs=1-q2), quantile(out[,3], probs=q2)))

}
#  est.bivar.lims <- function(N, n, beta, sig.b0, sig.b1, rho, sig.e, prob.grp, quants, n.low=11, n.high=11){
#         d        <- GenPopnData(N, n, beta, sig.b0, sig.b1, rho, sig.e, prob.grp)
#         data.tmp <- data.frame(id=d$id, Y=d$Y, X.time=d$X[,2])
#         out      <- matrix(unlist(lapply(split(data.tmp, data.tmp$id), LinRegFn)), byrow=TRUE, ncol=2)
#         out      <- cbind(c(tapply(data.tmp$Y, data.tmp$id, mean)), out)
#     
#         ## Do not forget that we are assuming the intercept and slope are independent here.  If we alter that assumption
#         ## we need to rethink this.  pr(a<Int<b, c<Slope<d)=pr(a<Int<b)*pr(c<Slope<d)
#         ReturnLims <- rbind( c( quantile(out[,2], probs=c( (1-sqrt(quants[1]))/2,  (1+sqrt(quants[1]))/2)), quantile(out[,3], probs= c( (1-sqrt(quants[1]))/2,  (1+sqrt(quants[1]))/2)  )),
#                       c( quantile(out[,2], probs=c( (1-sqrt(quants[2]))/2,  (1+sqrt(quants[2]))/2)), quantile(out[,3], probs= c( (1-sqrt(quants[2]))/2,  (1+sqrt(quants[2]))/2)  )))

        #mean((out[,2]>tmp[1,1] & out[,2]<tmp[1,2] & out[,3]>tmp[1,3] & out[,3]<tmp[1,4])) 
        #mean((out[,2]>tmp[2,1] & out[,2]<tmp[2,2] & out[,3]>tmp[2,3] & out[,3]<tmp[2,4])) 
        #a1       <- qnorm((sqrt(quants[1])+1)/2)
        #a2       <- qnorm((sqrt(quants[2])+1)/2) 
        #ReturnLims <- rbind( c(mean(out[,2])-a1*sqrt(var(out[,2])), mean(out[,2])+a1*sqrt(var(out[,2])),
        #                       mean(out[,3])-a1*sqrt(var(out[,3])), mean(out[,3])+a1*sqrt(var(out[,3]))),
        #                     c(mean(out[,2])-a2*sqrt(var(out[,2])), mean(out[,2])+a2*sqrt(var(out[,2])),
        #                       mean(out[,3])-a2*sqrt(var(out[,3])), mean(out[,3])+a2*sqrt(var(out[,3]))))
        #mean((out[,2]>ReturnLims[1,1] & out[,2]<ReturnLims[1,2] & out[,3]>ReturnLims[1,3] & out[,3]<ReturnLims[1,4])) 
        #mean((out[,2]>ReturnLims[2,1] & out[,2]<ReturnLims[2,2] & out[,3]>ReturnLims[2,3] & out[,3]<ReturnLims[2,4])) 
# 

#dat     <- GenPopnData(N, n=nhigh[count], beta[count,], sig.b0[count], sig.b1[count], rho, sig.e, prob.grp, n.low=nlow[count], n.high=nhigh[count])
#PopnQuantsBivar <- est.bivar.lims(N=20000, n=nhigh[count], beta=dat$beta, sig.b0=dat$sig.b0, sig.b1=dat$sig.b1, 
#                                 rho=dat$rho, sig.e=dat$sig.e, prob.grp=dat$prob.grp, 
#                                 quants=c(.6, .8), n.low=nhigh[count], n.high=nhigh[count])

#dat     <- GenPopnData(N, n=nhigh[count], beta[count,], sig.b0[count], sig.b1[count], rho, sig.e, prob.grp, n.low=nlow[count], n.high=nhigh[count])
#tmp <- ODS.Sampling.Bivar(dat=dat,                      ## a list generated from the GenPopnData() function
#                   PopnQuantsBivar=PopnQuantsBivar,          ## a matrix from est.bivar.lims function
#                   PopnPropInRectangle=0.8,          ## Proportion of subjects in the central rectangle
#                   TargetNSampledPerStratum=c(100,150), ## Theoretical (and possibly observed) number sampled per stratum
#                   SamplingStrategy="IndepODS")        ## Options are "IndepODS" and "DepODS"
#table(tmp$SampStratum)   
#cbind(tmp$SampStratum, tmp$Qi, tmp$Lims[1], tmp$Lims[2],tmp$Lims[3],tmp$Lims[4])

ODS.Sampling.Bivar <- function(dat,                      ## a list generated from the GenPopnData() function
                               PopnQuantsBivar,          ## a matrix from est.bivar.lims function
                               PopnPropInRectangle,          ## Proportion of subjects in the central rectangle
                               TargetNSampledPerStratum, ## Theoretical (and possibly observed) number sampled per stratum
                               SamplingStrategy){        ## Options are "IndepODS" and "DepODS"
         
    #######################################################################################################
   # dat=dat           ## a list generated from the GenPopnData() function
    #               PopnQuantsBivar=PopnQuantsBivar    ## a matrix from est.bivar.lims function
     #              PopnPropInRectangle=0.8            ## Proportion of subjects in the central rectangle
      #             TargetNSampledPerStratum=c(50,200) ## Theoretical (and possibly observed) number sampled per stratum
       #            SamplingStrategy="IndepODS"
    ########################################################################################################

    NCohort                    <- length(unique(dat$id))
    NStratumThry               <- round(NCohort*c(PopnPropInRectangle, (1-PopnPropInRectangle)))
    SampProbThry               <- TargetNSampledPerStratum / NStratumThry

    Lims <- (PopnPropInRectangle==.6)*PopnQuantsBivar[1,] + (PopnPropInRectangle==.8)*PopnQuantsBivar[2,] 
   
    uid <- unique(dat$id)
    ni  <- c(unlist(tapply(dat$id, dat$id, length)))  
    SampVar <- NULL  
    for(i in uid){ yi      <- dat$Y[dat$id==i]
                   xi      <- dat$X[dat$id==i,1:2]
                   SampVar <- rbind(SampVar, c(mean(yi), solve(t(xi)%*%xi) %*% t(xi) %*% yi))
    }
#   print(sum(SampVar[,2]>Lims[1] & SampVar[,2]<Lims[2]))
#   print(sum(SampVar[,3]>Lims[3] & SampVar[,3]<Lims[4]))
 
    SampStratum  <- ifelse(SampVar[,2]>Lims[1] & SampVar[,2]<Lims[2] & 
                           SampVar[,3]>Lims[3] & SampVar[,3]<Lims[4], 1,2) 
#    print(table(SampStratum))

    NperStratum  <- unlist(tapply(uid, SampStratum, length))
  
    SampProbiThry <- ifelse(SampStratum==1, SampProbThry[1],
                     ifelse(SampStratum==2, SampProbThry[2], NA))

    Sampled     <- rbinom(length(SampProbiThry), 1, SampProbiThry)
    SampProbObs <- c(tapply(Sampled, SampStratum, mean))

    SampProbiObs  <- ifelse(SampStratum==1, SampProbObs[1],
                     ifelse(SampStratum==2, SampProbObs[2], NA))

    ## Independent Sampling
    if (SamplingStrategy=="IndepODS") InODS <- uid[ Sampled==1] 
        
    TheSample <- dat$id %in% InODS                     
    X.ods     <- dat$X[TheSample,]
    Y.ods     <- dat$Y[TheSample]
    Z.ods     <- dat$Z[TheSample,]
    id.ods    <- dat$id[TheSample]
    if (SamplingStrategy=="IndepODS"){
        SampProbi.ods <- rep(SampProbiThry, ni)[TheSample]
        SampProb.ods  <- SampProbThry
    }
    #if (SamplingStrategy=="DepODS"){
    #	   SampProbi.ods <- rep(SampProbiObs, ni)[TheSample]
    #	   SampProb.ods  <- SampProbObs
    #}
    SampStrat.ods <- SampStratum[InODS]
    Qi <- SampVar[InODS,2:3]
    dup.id <- duplicated(dat$id)
    dat.univariate <- dat$X[!dup.id,]
    dat.univ.ods   <- dat.univariate[InODS,]

    list(X=X.ods, Y=Y.ods, Z=Z.ods, id=id.ods,
      SampProb=SampProb.ods, SampProbi=SampProbi.ods,
      N=dat$N, n=dat$n, beta=dat$beta, sig.b0=dat$sig.b0, 
      sig.b1=dat$sig.b1, rho=dat$rho, sig.e=dat$sig.e,
      prob.grp=dat$prob.grp,
      #cutpoint=c(C1,C2), 
      SampStratum=SampStrat.ods, Qi=Qi, dat.univ.ods=dat.univ.ods,
      cutpoint=Lims)
}

ODS.Sampling <- function(dat,                      ## a list generated from the GenPopnData() function
                         PopnQuants,               ## a matrix from est.quants function
                         w.function,               ## Response summary to sample on ("mean","intercept", or "slope")
                         quants,                   ## Population quantiles to define the theoretical sampling strata
                         TargetNSampledPerStratum, ## Theoretical (and possibly observed) number sampled per stratum
                         SamplingStrategy,         ## Options are "IndepODS" and "DepODS"
                         Univariate=FALSE){
                         	
    NCohort                    <- length(unique(dat$id))
    NStratumThry               <- round(NCohort*c(quants[1], quants[2]-quants[1], 1-quants[2]))
#    print(TargetNSampledPerStratum)
#    print(NStratumThry)
    SampProbThry               <- TargetNSampledPerStratum / NStratumThry

    C1 <- ifelse(quants[1]==.1   & w.function=="mean",       PopnQuants[1,1],
          ifelse(quants[1]==.1   & w.function=="intercept",  PopnQuants[2,1],
          ifelse(quants[1]==.1   & w.function=="slope",      PopnQuants[3,1],
          ifelse(quants[1]==.125   & w.function=="mean",       PopnQuants[1,2],
          ifelse(quants[1]==.125   & w.function=="intercept",  PopnQuants[2,2],
          ifelse(quants[1]==.125   & w.function=="slope",      PopnQuants[3,2],
          ifelse(quants[1]==.2   & w.function=="mean",       PopnQuants[1,3],
          ifelse(quants[1]==.2   & w.function=="intercept",  PopnQuants[2,3],
          ifelse(quants[1]==.2   & w.function=="slope",      PopnQuants[3,3])))))))))
    C2 <- ifelse(quants[2]==.8   & w.function=="mean",       PopnQuants[1,4],
          ifelse(quants[2]==.8   & w.function=="intercept",  PopnQuants[2,4],
          ifelse(quants[2]==.8   & w.function=="slope",      PopnQuants[3,4],
          ifelse(quants[2]==.875   & w.function=="mean",       PopnQuants[1,5],
          ifelse(quants[2]==.875   & w.function=="intercept",  PopnQuants[2,5],
          ifelse(quants[2]==.875   & w.function=="slope",      PopnQuants[3,5],
          ifelse(quants[2]==.9   & w.function=="mean",       PopnQuants[1,6],
          ifelse(quants[2]==.9   & w.function=="intercept",  PopnQuants[2,6],
          ifelse(quants[2]==.9   & w.function=="slope",      PopnQuants[3,6])))))))))   
 

    uid <- unique(dat$id)
    ni  <- c(unlist(tapply(dat$id, dat$id, length)))  
    SampVar <- NULL  
    for(i in uid){ yi      <- dat$Y[dat$id==i]
                   xi      <- dat$X[dat$id==i,1:2]
                   SampVar <- rbind(SampVar, c(mean(yi), solve(t(xi)%*%xi) %*% t(xi) %*% yi))
    }
    SampVar <- (w.function=="mean")*SampVar[,1] + 
               (w.function=="intercept")*SampVar[,2] + 
               (w.function=="slope")*SampVar[,3]

    SampStratum  <- ifelse(SampVar<C1, 1, 
                    ifelse(SampVar<C2, 2, 3))
    NperStratum  <- unlist(tapply(uid, SampStratum, length))
  
    SampProbiThry <- ifelse(SampVar<C1, SampProbThry[1],
                     ifelse(SampVar<C2, SampProbThry[2], SampProbThry[3]))

    Sampled     <- rbinom(length(SampProbiThry), 1, SampProbiThry)
    SampProbObs <- c(tapply(Sampled, SampStratum, mean))

    SampProbiObs  <- ifelse(SampVar<C1, SampProbObs[1],
                     ifelse(SampVar<C2, SampProbObs[2], SampProbObs[3]))
    #print(rbind(SampProbThry, SampProbObs))
    #print(cbind(SampProbiThry,SampProbiObs))
    

    ## Independent Sampling
    if (SamplingStrategy=="IndepODS") InODS <- uid[ Sampled==1] 
    #if (SamplingStrategy=="IndepODS"){
    #    InODS         <- uid[ rbinom(length(SampProbiThry), 1, SampProbiThry)==1] 
    #}
    ## Fixed sample size per stratum with dependent sampling
    if (SamplingStrategy=="DepODS"){
        InODS         <- c(sample(uid[SampStratum==1], TargetNSampledPerStratum[1], replace=FALSE),   
                           sample(uid[SampStratum==2], TargetNSampledPerStratum[2], replace=FALSE),
                           sample(uid[SampStratum==3], TargetNSampledPerStratum[3], replace=FALSE))
    }
    
    TheSample <- dat$id %in% InODS                     
    X.ods     <- dat$X[TheSample,]
    Y.ods     <- dat$Y[TheSample]
    Z.ods     <- dat$Z[TheSample,]
    id.ods    <- dat$id[TheSample]
    if (SamplingStrategy=="IndepODS"){
        SampProbi.ods <- rep(SampProbiThry, ni)[TheSample]
        SampProb.ods  <- SampProbThry
    }
    #if (SamplingStrategy=="IndepODS"){
    #    SampProbi.ods <- rep(SampProbiObs, ni)[TheSample]
    #    SampProb.ods  <- SampProbObs
    #} 
    if (SamplingStrategy=="DepODS"){
    	   SampProbi.ods <- rep(SampProbiObs, ni)[TheSample]
    	   SampProb.ods  <- SampProbObs
    }
    SampStrat.ods <- SampStratum[InODS]
    Qi <- SampVar[InODS]
    dup.id <- duplicated(dat$id)
    dat.univariate <- dat$X[!dup.id,]
    dat.univ.ods   <- dat.univariate[InODS,]

    list(X=X.ods, Y=Y.ods, Z=Z.ods, id=id.ods,
      SampProb=SampProb.ods, SampProbi=SampProbi.ods,
      N=dat$N, n=dat$n, beta=dat$beta, sig.b0=dat$sig.b0, 
      sig.b1=dat$sig.b1, rho=dat$rho, sig.e=dat$sig.e,
      prob.grp=dat$prob.grp,
      cutpoint=c(C1,C2), SampStratum=SampStrat.ods, Qi=Qi, dat.univ.ods=dat.univ.ods)
}

Random.Sampling <- function(d, n=225){
    s <- sample(unique(d$id), n)     
    TheSample <- d$id %in% s
    X.rand     <- d$X[TheSample,]
    Y.rand     <- d$Y[TheSample]
    Z.rand     <- d$Z[TheSample,]
    id.rand    <- d$id[TheSample]

    list(X=X.rand, Y=Y.rand, Z=Z.rand, id=id.rand,
         N=d$N, n=d$n, beta=d$beta, sig.b0=d$sig.b0, 
         sig.b1=d$sig.b1, rho=d$rho, sig.e=d$sig.e,
         prob.grp=d$prob.grp)
}
