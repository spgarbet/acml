# This file comprises unit test for the estfit.c code
# At present it is does not include a test for each function, ideally it would
# These "tests" are simple assertions that a single set of parameters return an expected value
# They do not indicate that the code as a whole produces the correct result, just that each
# function tested in isolation produces a valid result
#
# Note: It is difficult to test the memory frame code used by BLAS here. If these unit tests pass, 
#       and the integration tests in 2.R fail, it is quite possible to have misspecified lda, ldb or ldc when
#       class BLAS or LAPACK.

test.vi_calc <- function()
{
    sigma0    <- sqrt(3.87)
    sigma1    <- sqrt(0.25)
    sigmae    <- sqrt(0.49)
    rho       <- 0.2
    sigmai    <- as.matrix(c(sigma0*sigma0, rho*sigma0*sigma1,  rho*sigma0*sigma1, sigma1*sigma1), nrow=2)
    ni        <- 11
    n         <- ni
    zi        <- as.matrix(data.frame(intercept=rep(1,ni), time=(1:11/5 - 1.2)))
    vi        <- rep(0, ni*ni)
    vi_tmp    <- rep(0, ni*ni)

    correct   <- vi.calc(zi, sigma0, sigma1, rho, sigmae)

    result    <- .C("vi_calc", 
                    zi      = as.double(zi),
                    n       = as.integer(n),
                    nz      = as.integer(ni),
                    sigma_e = as.double(sigmae),
                    Sigma_i = as.double(sigmai),
                    vi      = as.double(vi),
                    vi_tmp  = as.double(vi_tmp))

    checkTrue( sum(abs(correct-result$vi)) < tolerance)
}

test.subject_ll_lme <- function()
{
    beta      <- c(8.5,-0.27, -0.77, 0.48)
    sigma0    <- sqrt(4.04)
    sigma1    <- sqrt(0.255)
    sigmae    <- sqrt(0.489)
    rho       <- 0.18

    yi        <- c(6.216108,6.275344,5.740897,5.011603,5.564250,5.703305,4.940070,5.913187,5.434388,5.805954,6.564057)
    
    ni        <- length(yi)
    xi        <- as.matrix(data.frame(intercept=rep(1,ni), time=(1:11/5 - 1.2), grp=rep(1, ni), snp=rep(1,ni) ))
    zi        <- xi[,1:2]
    
    resid     <- rep(0, ni)
    vi        <- vi.calc(zi, sigma0, sigma1, rho, sigmae)
    vi_inv    <- rep(0, ni*ni)

    n         <- ni
    nbeta     <- length(beta)

    vi_tmp    <- matrix(rep(0, ni*ni), ncol=ni)
    vi2_tmp   <- matrix(rep(0, ni*ni), ncol=ni)
    ll_lme    <- 0.0
    ipvt      <- rep(0, ni)

    correct <- subject.ll.lme(yi, xi, beta, vi)
    
    result  <- .C("subject_ll_lme",
                   yi       = as.double(yi),
                   xi       = as.double(xi),
                   beta     = as.double(beta),
                   vi       = as.double(vi),
                   n        = as.integer(n),
                   ni       = as.integer(ni),
                   nbeta    = as.integer(nbeta),
                   ll_lme   = as.double(ll_lme),
                   vi_inv   = as.double(vi_inv),
                   resid    = as.double(resid),
                   vi_tmp   = as.double(vi_tmp),
                   vi2_tmp  = as.double(vi2_tmp),
                   ipvt     = as.integer(ipvt))
                   
    checkTrue( abs(correct-result$ll_lme) < tolerance)
    checkTrue( sum(abs(result$vi_inv-solve(vi))) < tolerance)
    checkTrue( sum(abs((yi - xi %*% beta)-result$resid)) < tolerance)
}

test.lci <- function()
{
    cutpoints  <- c(-1, 1)
    ncutpoints <- length(cutpoints)
    SampProb   <- c(0.8, 0.4, 0.8)
    mu_q       <- 1
    sigma_q    <- 1

    correct <- lci(cutpoints, SampProb, mu_q, sigma_q)
    
    answer        <- 0

    result  <- .C("lci", cutpoints  = as.double(cutpoints),
                         ncutpoints = as.integer(ncutpoints),
                         sampprob   = as.double(SampProb),
                         mu_q       = as.double(mu_q),
                         Sigma_q    = as.double(sigma_q),
                         lci        = as.double(answer))

    checkTrue( abs(correct-result$lci) < tolerance)
}


test.dgtrpr <- function()
{
    m     <- 2
    alpha <- 3.0
    lda   <- 7
    A     <- matrix(1.0*(1:(lda*2)), ncol=2)
    ldb   <- 11
    B     <- matrix(1.0*((lda*2+1):(lda*2+ldb*2)), ncol=2)
    y     <- 0
    result <- .C("dgtrpr", m     = as.integer(m), 
                           alpha = as.double(alpha),
                           A     = as.double(A),
                           lda   = as.integer(lda),
                           B     = as.double(B),
                           ldb   = as.integer(ldb),
                           y     = as.double(y))

    checkTrue( (round(alpha*sum(diag((A[1:m,1:m] %*% B[1:m, 1:m])))) - round(result$y)) == 0)
}

test.slope <- function()
{
    grad        <- 0
    ni          <- 4
    n           <- ni*2
    inv_v       <- matrix(1.0*(1:16), ncol=ni)
    resid       <- (1:ni) - 0.5*(ni+1)
    vi_tmp      <- matrix(1.0*(1:16), ncol=ni) # dvk
    vi2_tmp     <- matrix(rep(0, 16), ncol=ni)
    
    result      <- .C("slope", grad     = as.double(grad), 
                               inv_v    = as.double(inv_v),
                               resid    = as.double(resid),
                               n        = as.integer(n),
                               ni       = as.integer(ni), 
                               vi_tmp   = as.double(vi_tmp),
                               vi2_tmp  = as.double(vi2_tmp)
                     )
                     
    checkTrue(result$grad == 6142)
}

test.subject_gradient_ll_lme <- function()
{
    beta      <- c(10.1,-0.24, -0.8, 0.45)
    sigma0    <- sqrt(4.1)
    sigma1    <- sqrt(0.24)
    sigmae    <- sqrt(0.51)
    rho       <- 0.2
    
    yi        <- c(6.216108,6.275344,5.740897,5.011603,5.564250,5.703305,4.940070,5.913187,5.434388,5.805954,6.564057)
    ni        <- length(yi)
    xi        <- as.matrix(data.frame(intercept=rep(1,ni), time=(1:11/5 - 1.2), grp=rep(1, ni), snp=rep(1,ni) ))
    zi        <- xi[,1:2]
    grad      <- rep(0, length(beta)+4)
    
    resid     <- yi - xi %*% beta
    vi        <- vi.calc(zi, sigma0, sigma1, rho, sigmae)
    inv.v     <- solve(vi)
    
    n         <- ni
    npar      <- length(beta)
    
    vi_tmp      <- matrix(rep(0, ni*ni), ncol=ni)
    vi2_tmp     <- matrix(rep(0, ni*ni), ncol=ni)
    
    correct <- subject.gradient.ll.lme(yi, xi, zi, beta, sigma0, sigma1, rho, sigmae)
    
    result      <- .C("subject_gradient_ll_lme",
                      sub_grad  = as.double(grad),
                      resid     = as.double(resid),
                      vi        = as.double(vi),
                      vi_inv    = as.double(inv.v),
                      xi        = as.double(xi),
                      yi        = as.double(yi),
                      zi        = as.double(zi),
                      n         = as.integer(n),
                      ni        = as.integer(ni),
                      npar      = as.integer(npar),
                      sigma0    = as.double(sigma0),
                      sigma1    = as.double(sigma1),
                      rho       = as.double(rho),
                      sigmae    = as.double(sigmae),
                      vi_tmp    = as.double(vi_tmp),
                      vi2_tmp   = as.double(vi2_tmp) )
                      
    checkTrue(sum(abs(result$sub_grad - correct$gr)) < tolerance)
}

test.comp_ac_dvk <- function()
{
    ni   <- 11
    n    <- ni
    dvk  <- rep(0, ni*ni)
    zi   <- as.matrix(data.frame( intercept=rep(1,ni), time=(1:11/5 - 1.2) ))
    sigma_p <- matrix(c(1,2,2,3), ncol=2)
    vi2_tmp <- rep(0, ni*ni)
    
    correct <- zi %*% sigma_p %*% t(zi)

    result      <- .C("comp_ac_dvk",
                      dvk = as.double(dvk),
                      zi = as.double(zi),
                      n = as.integer(n),
                      ni = as.integer(ni),
                      sigma_p = as.double(sigma_p),
                      vi2_tmp = as.double(vi2_tmp)
                     )

    checkTrue(sum(abs(result[['dvk']] - correct)) < tolerance)
}


test.ac_slope <- function()
{
    ni   <- 11
    n    <- ni
    dvk  <- matrix(rnorm(ni*ni),ncol=ni)
    correction <- 0.0
    f_alpha_k <- rnorm(1)
    vi2_tmp <- rep(0, ni*ni)
    wi <- t(rnorm(ni))
    
    reference <- f_alpha_k *  wi %*% dvk %*% t(wi)

    result      <- .C("ac_slope",
                      correction = as.double(correction),
                      wi = as.double(wi),
                      dvk = as.double(dvk),
                      ni = as.integer(ni),
                      f_alpha_k = as.double(f_alpha_k),
                      vi2_tmp = as.double(vi2_tmp)
                     )

    checkTrue(abs((result[['correction']] - reference[1,1])/reference[1,1]) < tolerance)
}

test.lci_bivar <- function()
{
    ni        <- 11
    wi        <- rbind(t(rep(1/ni, ni)), 1:ni / (2*(ni-1)) - (0.55/2) )
    beta      <- rnorm(4, sd=4)
    cutpoints <- c(-1, 1, -2, 2)
    ncutpoints<- length(cutpoints)
    sampprob  <- c(0.4, 0.8)
    sigma0    <- sqrt(4.1)
    sigma1    <- sqrt(0.24)
    sigmae    <- sqrt(0.51)
    rho       <- 0.2


    xi        <- as.matrix(data.frame(intercept=rep(1,ni), time=(1:11/5 - 1.2), grp=rep(1, ni), snp=rep(1,ni) ))
    zi        <- xi[,1:2]
    vi        <- vi.calc(zi, sigma0, sigma1, rho, sigmae)

    mu          <- xi %*% beta
    mu_q        <- as.vector(wi %*% mu)
    sigma_q     <- wi %*% vi %*% t(wi)
                
    ac          <- 0

    result <- .C("lci_bivar",
                 cutpoints  = as.double(cutpoints),
                 ncutpoints = as.integer(ncutpoints),
                 sampprob   = as.double(sampprob),
                 mu_q       = as.double(mu_q),
                 sigma_q    = as.double(sigma_q),
                 ac         = as.double(ac))
    reference <- lci.bivar(cutpoints, sampprob, mu_q, sigma_q)
    checkTrue(abs((result[['ac']] - reference)/reference) < tolerance)
}

ascertainment.correction <- function(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb){
    vi      <- vi.calc(zi, sigma0, sigma1, rho, sigmae)
    mu      <- xi %*% beta
    mu_q    <- (wi %*% mu)[,1]  
    sigma_q <- sqrt((wi %*% vi %*% t(wi))[1,1])
    log(lci(cutpoints, SampProb, mu_q, sigma_q))
}

test.ascertainment_corr_univar <- function()
{
   beta      <- rnorm(4, sd=4)
   sigma0    <- sqrt(runif(1, 0.1, 4))
   sigma1    <- sqrt(runif(1, 0.1, 4))
   sigmae    <- sqrt(runif(1, 0.1, 4))
   rho       <- runif(1, -0.2, 0.2)
   
   cutpoints <- c(-1, 1)
   ncutpoints<- length(cutpoints)
   sampprob <- c(0.8, 0.4, 0.8)


   ni        <- 11
   xi        <- as.matrix(data.frame(intercept=rep(1,ni), time=(1:11/5 - 1.2), grp=rep(1, ni), snp=rep(1,ni) ))
   zi        <- xi[,1:2]
   grad      <- rep(0, length(beta)+4)

   vi        <- vi.calc(zi, sigma0, sigma1, rho, sigmae)

   n         <- ni
   npar      <- length(beta)

   vi_tmp    <- matrix(rep(0, ni*ni), ncol=ni)

   wi <- t(rnorm(ni))

   ac        <- 23

   result <- .C("exp_ascertainment_correction_univar",
                 n          = as.integer(n),
                 ni         = as.integer(ni),
                 npar       = as.integer(npar),
                 xi         = as.double(xi),
                 vi         = as.double(vi),
                 wi         = as.double(wi),
                 beta       = as.double(beta),
                 cutpoints  = as.double(cutpoints),
                 ncutpoints = as.integer(ncutpoints),
                 sampprob   = as.double(sampprob),
                 vi_tmp     = as.double(vi_tmp),
                 ac         = as.double(ac))
                 
    reference <- ascertainment.correction(NA, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, sampprob)
    
    checkTrue(abs((log(result[['ac']]) - reference)/reference) < tolerance)
}

test.ascertainment_corr_bivar <- function()
{
   beta      <- rnorm(4, sd=4)
   sigma0    <- sqrt(runif(1, 0.1, 4))
   sigma1    <- sqrt(runif(1, 0.1, 4))
   sigmae    <- sqrt(runif(1, 0.1, 4))
   rho       <- runif(1, -0.2, 0.2)
   
   cutpoints <- c(-1, 1, -2, 2)
   ncutpoints<- length(cutpoints)
   sampprob  <- c(0.4, 0.8)

   ni        <- 11
   xi        <- as.matrix(data.frame(intercept=rep(1,ni), time=(1:11/5 - 1.2), grp=rep(1, ni), snp=rep(1,ni) ))
   zi        <- xi[,1:2]
   grad      <- rep(0, length(beta)+4)

   vi        <- vi.calc(zi, sigma0, sigma1, rho, sigmae)

   n         <- ni
   npar      <- length(beta)

   vi_tmp    <- matrix(rep(0, ni*ni), ncol=ni)

   wi        <- rbind(t(rep(1/ni, ni)), 1:ni / (2*(ni-1)) - (0.55/2) )

   ac        <- 0

   result <- .C("exp_ascertainment_correction_bivar",
                 n          = as.integer(n),
                 ni         = as.integer(ni),
                 npar       = as.integer(npar),
                 xi         = as.double(xi),
                 vi         = as.double(vi),
                 wi         = as.double(wi),
                 beta       = as.double(beta),
                 cutpoints  = as.double(cutpoints),
                 ncutpoints = as.integer(ncutpoints),
                 sampprob   = as.double(sampprob),
                 vi_tmp     = as.double(vi_tmp),
                 ac         = as.double(ac))
    reference <- ascertainment.correction.bivar(NA, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, sampprob)

    checkTrue(abs((log(result[['ac']]) - reference)/reference) < tolerance)
}

test.ascertainment_gradient_corr_univar <- function()
{
    beta      <- c(10.1,-0.24, -0.8, 0.45)
    sigma0    <- sqrt(4.1)
    sigma1    <- sqrt(0.24)
    sigmae    <- sqrt(0.51)
    rho       <- 0.2
    cutpoints <- c(-1, 1)
    ncutpoints<- length(cutpoints)
    samp_prob <- c(0.8, 0.4, 0.8)

    ni        <- 11
    yi        <- sort(rnorm(ni, 5))
    xi        <- as.matrix(data.frame(intercept=rep(1,ni), time=(1:11/5 - 1.2), grp=rep(1, ni), snp=rep(1,ni) ))
    zi        <- xi[,1:2]
    grad      <- rep(0, length(beta)+4)

    vi        <- vi.calc(zi, sigma0, sigma1, rho, sigmae)

    n         <- ni
    npar      <- length(beta)

    vi_tmp    <- matrix(rep(0, ni*ni), ncol=ni)
    vi2_tmp   <- matrix(rep(0, ni*ni), ncol=ni)

    wi        <- t(rep(1/ni, ni))

    corr      <- rep(0, npar+4)

    correct <- ascertainment.gradient.correction(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, samp_prob)

    result      <- .C("ascertainment_gradient_corr_univar",
                      corr      = as.double(corr),
                      yi        = as.double(yi),
                      xi        = as.double(xi),
                      zi        = as.double(zi),
                      wi        = as.double(wi),
                      vi        = as.double(vi),
                      beta      = as.double(beta),
                      n         = as.integer(n),
                      ni        = as.integer(ni),
                      npar      = as.integer(npar),
                      sigma0    = as.double(sigma0),
                      sigma1    = as.double(sigma1),
                      rho       = as.double(rho),
                      sigmae    = as.double(sigmae),
                      cutpoints = as.double(cutpoints),
                      ncutpoints= as.integer(ncutpoints),
                      samp_prob = as.double(samp_prob),
                      vi_tmp    = as.double(vi_tmp),
                      vi2_tmp   = as.double(vi2_tmp) )
                      
    checkTrue(sum(abs(result[['corr']] - correct)) < tolerance)
}

test.ascertainment_gradient_corr_bivar <- function()
{
    beta      <- c(10.1,-0.24, -0.8, 0.45)
    sigma0    <- sqrt(4.15)
    sigma1    <- sqrt(0.19)
    sigmae    <- sqrt(0.534)
    rho       <- 0.2
    cutpoints <- c(-1, 1, -1, 1)
    ncutpoints<- length(cutpoints)
    samp_prob <- c(0.8, 0.4, 0.8)

    ni        <- 2
    yi        <- sort(rnorm(ni, 5))
    xi        <- as.matrix(data.frame(intercept=rep(1,ni), time=(1:ni/5 - 1.2), grp=rep(1, ni), snp=rep(1,ni) ))
    zi        <- xi[,1:2]
    grad      <- rep(0, length(beta)+4)

    vi        <- vi.calc(zi, sigma0, sigma1, rho, sigmae)

    n         <- ni
    npar      <- length(beta)

    vi_tmp    <- matrix(rep(0, ni*ni), ncol=ni)
    vi2_tmp   <- matrix(rep(0, ni*ni), ncol=ni)

    wi        <- (solve(t(zi[,1:2])%*% zi[,1:2]) %*% t(zi[,1:2]))

    corr      <- rep(0, 2*(npar+4)) # Needs double for working tmp memory

    correct   <- ascertainment.gradient.correction.bivar(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, samp_prob)
    result    <- .C("ascertainment_gradient_corr_bivar",
                      corr      = as.double(corr),
                      yi        = as.double(yi),
                      xi        = as.double(xi),
                      zi        = as.double(zi),
                      wi        = as.double(wi),
                      vi        = as.double(vi),
                      beta      = as.double(beta),
                      n         = as.integer(n),
                      ni        = as.integer(ni),
                      npar      = as.integer(npar),
                      sigma0    = as.double(sigma0),
                      sigma1    = as.double(sigma1),
                      rho       = as.double(rho),
                      sigmae    = as.double(sigmae),
                      cutpoints = as.double(cutpoints),
                      ncutpoints= as.integer(ncutpoints),
                      samp_prob = as.double(samp_prob),
                      vi_tmp    = as.double(vi_tmp),
                      vi2_tmp   = as.double(vi2_tmp) )

    for(i in 1:(npar+4))
    {
        checkTrue(abs( (result[['corr']][i] - correct[i])/correct[i] ) < 1e-3)
    }
}

test.llscore <- function()
{
    beta        <- rnorm(4, sd=4)
    sigma0      <- sqrt(runif(1, 0.1, 4))
    sigma1      <- sqrt(runif(1, 0.1, 4))
    sigmae      <- sqrt(runif(1, 0.1, 4))
    rho         <- runif(1, -0.2, 0.2)

    rho         <- 0.1
    nlow        <- 11
    nhigh       <- 11
    prob.grp    <- 0.10

    count       <- 0
    N           <- 1000

    answer      <- append(beta, c(log(sigma0), log(sigma1), rho, log(sigmae)))

    inits       <- c(beta, log(sigma0), log(sigma1), log((1+rho)/(1-rho)), log(sigmae))
    dat         <- GenPopnData(N, n=nhigh, beta, sigma0, sigma1, rho, sigmae, prob.grp, n.low=nlow, n.high=nhigh)

    PopnQuants <- est.quants(2500, nhigh, beta, sigma0, sigma1, rho, sigmae, prob.grp, c(0.1,.125, .2,.8, .875, .9), n.low=nlow, n.high=nhigh)

    for(method in c("mean", "intercept", "slope"))
    {
        ods <- ODS.Sampling(dat=dat,
                            PopnQuants,  
                            method,
                            c(0.2,0.8),
                            TargetNSampledPerStratum=c(20,50,20), 
                            "IndepODS")
        reference <- LogLikeAndScore(inits, ods$Y, ods$X, ods$Z, ods$id*1.0, method,
                                            ods$cutpoint,
                                            ods$SampProb,
                                            ods$SampProbi)
        result    <- LogLikeAndScore.C(inits, ods$Y, ods$X, ods$Z, ods$id*1.0, method,
                                            ods$cutpoint,
                                            ods$SampProb,
                                            ods$SampProbi)
    }

    checkTrue(abs( (reference - result)/reference ) < tolerance )
    for(i in 1:8)
    {
        checkTrue(abs( (attr(reference, "gradient")[i] - attr(result, "gradient")[i])/attr(reference, "gradient")[i] ) < tolerance)
    }
}
