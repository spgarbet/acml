# aclm methods and utilities

# Most of this code was lifted from the lme4 package by Douglas Bates and adapted as needed 
# for working with ascertainment corrected models.

setClass("acmlObj",
    representation(cov = "numeric", logl="numeric", robcov="numeric", data="data.frame"),
    contains = "numeric")


  #####################################################################
 ## Create the nonzero pattern for the sparse matrix Cm from A.
##  ncol(A) is s * ncol(Cm).  The s groups of ncol(Cm) consecutive
##  columns in A are overlaid to produce Cm.
##
#createCm <- function(A, s)
#{
#    stopifnot(is(A, "dgCMatrix"))
#    s <- as.integer(s)[1]
#    if (s == 1L) return(A)
#    if ((nc <- ncol(A)) %% s)
#        stop(gettextf("ncol(A) = %d is not a multiple of s = %d",
#                      nc, s))
#    ncC <- as.integer(nc / s)
#    TA <- as(A, "TsparseMatrix")
#    as(new("dgTMatrix", Dim = c(nrow(A), ncC),
#           i = TA@i, j = as.integer(TA@j %% ncC), x = TA@x),
#       "CsparseMatrix")
#}

  #####################################################################
 ## Create model frame
##
## Create the model frame, X, Y, offset and terms
##
## mc        - matched call of calling function
## formula   - two-sided formula
## contrasts - contrasts argument
##
aclmeFrames <- function(mc, formula, contrasts)
{
    mf <- mc
    m <- match(c("data", "subset", "na.action", "offset"),
               names(mf), 0)
    mf <- mf[c(1, m)]

    ## The model formula for evaluation of the model frame.  It looks
    ## like a linear model formula but includes any random effects
    ## terms and any names of parameters used in a nonlinear mixed model.
    frame.form <- subbars(formula)      # substitute `+' for `|'

    ## The model formula for the fixed-effects terms only.
    fixed.form <- nobars(formula)       # remove any terms with `|'
    if (!inherits(fixed.form, "formula"))
      ## RHS is empty - use `y ~ 1'
      fixed.form <- as.formula(substitute(foo ~ 1, list(foo = fixed.form)))

    ## attach the correct environment
    environment(fixed.form) <- environment(frame.form) <- environment(formula)

    ## evaluate a model frame
    mf$formula <- frame.form
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    fe <- mf                            # save a copy of the call
    mf <- eval(mf, parent.frame(2))

    ## evaluate the terms for the fixed-effects only (used in anova)
    fe$formula <- fixed.form
    fe <- eval(fe, parent.frame(2)) # allow model.frame to update them

    ## response vector
    Y <- model.response(mf, "any")
    ## avoid problems with 1D arrays, but keep names
    if(length(dim(Y)) == 1) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if(!is.null(nm)) names(Y) <- nm
    }
    mt <- attr(fe, "terms")

    ## Extract X checking for a null model. This check shouldn't be
    ## needed because an empty formula is changed to ~ 1 but it can't hurt.
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf, contrasts) else matrix(,NROW(Y),0)
    storage.mode(X) <- "double"      # when ncol(X) == 0, X is logical
    fixef <- numeric(ncol(X))
    names(fixef) <- colnames(X)
    dimnames(X) <- NULL

    ## Extract the weights and offset.  For S4 classes we want the
    ## `not used' condition to be numeric(0) instead of NULL
    off <- model.offset(mf); if (is.null(off)) off <- numeric(0)

    ## check weights and offset

    if(length(off) && length(off) != NROW(Y))
        stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                  length(off), NROW(Y)))

    ## remove the terms attribute from mf
    attr(mf, "terms") <- mt

    bars        <- expandSlash(findbars(formula[[3]]))
    ids <-      names(bars) <- unlist(lapply(bars, function(x) deparse(x[[3]])))

    if(length(bars) != 1)
        stop("Only one random effect is supported for acml method at present")
    Z <- lapply(bars,
                 function(x)
                 {
                     ff <- eval(substitute(as.factor(fac)[,drop = TRUE],
                                           list(fac = x[[3]])), mf)
                     im <- as(ff, "sparseMatrix") # transpose of indicators
                     ## Could well be that we should rather check earlier .. :
                     if(!isTRUE(validObject(im, test=TRUE)))
                         stop("invalid conditioning factor in random effect: ", format(x[[3]]))

                     mm <- model.matrix(eval(substitute(~ expr, # model matrix
                                                        list(expr = x[[2]]))),
                                        mf)
# This removes the intercept
#                     if (rmInt) {
#                         if (is.na(icol <- match("(Intercept)", colnames(mm)))) break
#                         if (ncol(mm) < 2)
#                             stop("lhs of a random-effects term cannot be an intercept only")
#                         mm <- mm[ , -icol , drop = FALSE]
#                     }
                     mm
                 })[[1]]

    list(Y=Y, X=X, Z=Z, ids=mf[,ids], off = as.double(off), mf = mf, fixef = fixef,
        ranef = unlist(lapply(bars, function(x) deparse(x[[2]])))[[1]]
    )
}

  #####################################################################
 ## Create model matrices
## 
## Create the list of model matrices from the random-effects terms in
## the formula and the model frame.
##
## formula : model formula
## fr      : list with '$mf': model frame; '$X': .. matrix
## rmInt   : logical scalar - should the `(Intercept)` column be removed before creating Zt
## drop    : logical scalar indicating if elements with numeric
##           value 0 should be dropped from the sparse model matrices
##
## returns : a list with components named trms, fl and dims
##   ba
##aclmeFactorList <- function(formula, fr, rmInt, drop)
#{
#    mf <- fr$mf
#    ## record dimensions and algorithm settings
#
#    ## create factor list for the random effects
#    bars <- expandSlash(findbars(formula[[3]]))
#    if (!length(bars)) stop("No random effects terms specified in formula")
#    names(bars) <- unlist(lapply(bars, function(x) deparse(x[[3]])))
#    fl <- lapply(bars,
#                 function(x)
#             {
#                 ff <- eval(substitute(as.factor(fac)[,drop = TRUE],
#                                       list(fac = x[[3]])), mf)
#                 im <- as(ff, "sparseMatrix") # transpose of indicators
#                 ## Could well be that we should rather check earlier .. :
#                 if(!isTRUE(validObject(im, test=TRUE)))
#                     stop("invalid conditioning factor in random effect: ", format(x[[3]]))
#
#                 mm <- model.matrix(eval(substitute(~ expr, # model matrix
#                                                    list(expr = x[[2]]))),
#                                    mf)
#                 if (rmInt) {
#                     if (is.na(icol <- match("(Intercept)", colnames(mm)))) break
#                     if (ncol(mm) < 2)
#                         stop("lhs of a random-effects term cannot be an intercept only")
#                     mm <- mm[ , -icol , drop = FALSE]
#                 }
#                 ans <- list(f = ff,
#                             A = do.call(rBind,
#                             lapply(seq_l(ncol(mm)), function(j) im)),
#                             Zt = do.call(rBind,
#                             lapply(seq_l(ncol(mm)),
#                                    function(j) {im@x <- mm[,j]; im})),
#                             ST = matrix(0, ncol(mm), ncol(mm),
#                             dimnames = list(colnames(mm), colnames(mm))))
#                 if (drop) {
#                     ## This is only used for naclme models.
#                     ## Need to do something more complicated for A
#                     ## here.  Essentially you need to create a copy
#                     ## of im for each column of mm, im@x <- mm[,j],
#                     ## create the appropriate number of copies,
#                     ## prepend matrices of zeros, then rBind and drop0.
#                     ans$A@x <- rep(0, length(ans$A@x))
#                     ans$Zt <- drop0(ans$Zt)
#                 }
#                 ans
#             })
#    dd <-
#        VecFromNames(dimsNames, "integer",
#                     c(list(n = nrow(mf), p = ncol(fr$X), nt = length(fl),
#                            q = sum(sapply(fl, function(el) nrow(el$Zt)))),
#                       dimsDefault))
#    ## order terms by decreasing number of levels in the factor but don't
#    ## change the order if this is already true
#    nlev <- sapply(fl, function(el) length(levels(el$f)))
#    ## determine the number of random effects at this point
#    if (any(diff(nlev)) > 0) fl <- fl[rev(order(nlev))]
#    ## separate the terms from the factor list
#    trms <- lapply(fl, "[", -1)
#    names(trms) <- NULL
#    fl <- lapply(fl, "[[", "f")
#    attr(fl, "assign") <- seq_along(fl)
#    ## check for repeated factors
#    fnms <- names(fl)
#    if (length(fnms) > length(ufn <- unique(fnms))) {
#        ## check that the lengths of the number of levels coincide
#        fl <- fl[match(ufn, fnms)]
#        attr(fl, "assign") <- match(fnms, ufn)
#    }
#    names(fl) <- ufn
#    ## check for nesting of factors
#    dd["nest"] <- all(sapply(seq_along(fl)[-1],
#                             function(i) isNested(fl[[i-1]], fl[[i]])))
#
#    list(trms = trms, fl = fl, dims = dd)
#}


  #####################################################################
 ## dimsNames and devNames are in the package's namespace rather than
## in the function aclmeFactorList because the function sparseRasch
## needs to access them.
#dimsNames <- c("nt", "n", "p", "q", "s", "np", "LMM",
#               "fTyp", "lTyp", "vTyp", "nest", "useSc", "nAGQ",
#               "verb", "mxit", "mxfn", "cvg")
#
#dimsDefault <- list(s    = 1L,   # identity mechanistic model
#                    mxit = 300L, # maximum number of iterations
#                    mxfn = 900L, # maximum number of function evaluations
#                    verb = 0L,   # no verbose output
#                    np   = 0L,   # number of parameters in ST
#                    LMM  = 0L,   # not a linear mixed model
#                    fTyp = 2L,   # default family is "gaussian"
#                    lTyp = 5L,   # default link is "identity"
#                    vTyp = 1L,   # default variance function is "constant"
#                    useSc= 1L,   # default is to use the scale parameter
#                    nAGQ = 1L,   # default is Laplace
#                    cvg  = 0L)   # no optimization yet attempted
#
#devNames <- c("ML", "ldL2", "ldRX2", "sigmaML",
#              "pwrss", "disc", "usqr", "wrss",
#              "dev", "llik", "NULLdev")

  #####################################################################
 ## Generate a named vector of the given mode
##
VecFromNames <- function(nms, mode = "numeric", defaults = list())
{
    ans <- vector(mode = mode, length = length(nms))
    names(ans) <- nms
    ans[] <- NA
    if ((nd <- length(defaults <- as.list(defaults))) > 0) {
        if (length(dnms <- names(defaults)) < nd)
            stop("defaults must be a named list")
        stopifnot(all(dnms %in% nms))
        ans[dnms] <- as(unlist(defaults), mode)
    }
    ans
}

  #####################################################################
 ## Control parameters for aclme search
##
aclmeControl <- function(msVerbose = getOption("verbose"),
                         maxIter   = 300L,
                         maxFN     = 900L)
{
    stopifnot(maxIter >= 0, maxFN >= 0)
    
    list(maxIter = as.integer(maxIter),
         maxFN = as.integer(maxFN),
         msVerbose = as.integer(msVerbose)) # "integer" on purpose
}

  #####################################################################
 ## Check that the 'STnew' argument matches the form of ST.
##
#checkSTform <- function(ST, STnew)
#
#{
#    stopifnot(is.list(STnew), length(STnew) == length(ST),
#              all.equal(names(ST), names(STnew)))
#    lapply(seq_along(STnew), function (i)
#           stopifnot(class(STnew[[i]]) == class(ST[[i]]),
#                     all.equal(dim(STnew[[i]]), dim(ST[[i]]))))
#    all(unlist(lapply(STnew, function(m) all(diag(m) > 0))))
#}

  #####################################################################
 ## Create the standard versions of flist, Zt, Gp, ST, A, Cm, Cx, and L.
##  Update dd.
##
#mkZt <- function(FL, start, s = 1L)
#{
#    dd <- FL$dims
#    fl <- FL$fl
#    asgn <- attr(fl, "assign")
#    trms <- FL$trms
#    ST <- lapply(trms, `[[`, "ST")
#    Ztl <- lapply(trms, `[[`, "Zt")
#    Zt <- do.call(rBind, Ztl)
#    Zt@Dimnames <- vector("list", 2)
#    Gp <- unname(c(0L, cumsum(sapply(Ztl, nrow))))
#
##FIXMEFIXMEFIXMEFIXMEFIXMEFIXME
### OUCH!!!
### This zeros all the allocated matrices and slots
####    .Call(mer_ST_initialize, ST, Gp, Zt)
#
#    A <- do.call(rBind, lapply(trms, `[[`, "A"))
#    rm(Ztl, FL)                         # because they could be large
#    nc <- sapply(ST, ncol)         # of columns in els of ST
#    Cm <- createCm(A, s)
#
##FIXMEFIXMEFIXMEFIXMEFIXMEFIXME
### OUCH!!!
####   L <- .Call(mer_create_L, Cm)
#    L <- chol(Cm)
#    if (s < 2) Cm <- new("dgCMatrix")
#    if (!is.null(start) && checkSTform(ST, start)) ST <- start
#
#    nvc <- sapply(nc, function (qi) (qi * (qi + 1))/2) # no. of var. comp.
#    
#### FIXME: Check number of variance components versus number of
#### levels in the factor for each term. Warn or stop as appropriate
#
#    dd["np"] <- as.integer(sum(nvc))    # number of parameters in optimization
#    dev <- VecFromNames(devNames, "numeric")
#    fl <- do.call(data.frame, c(fl, check.names = FALSE))
#    attr(fl, "assign") <- asgn
#
#    list(Gp = Gp, ST = ST, A = A, Cm = Cm, L = L, Zt = Zt,
#         dd = dd, dev = dev, flist = fl)
#}

  #####################################################################
 ## Modifications to aclme often involve modifying model matrices before
## creating and optimizing the mer object.  Everything past the model
## matrices is encapsulated in this function
#acml_finalize <- function(fr, FL, start, verbose)
#{
#    Y <- as.double(fr$Y)
#    if (is.list(start) && all(sort(names(start)) == sort(names(FL))))
#        start <- list(ST = start)
#    if (is.numeric(start)) start <- list(STpars = start)
#    dm <- mkZt(FL, start[["ST"]])
#    
#    ### This checks that the number of levels in a grouping factor < n
#    ### Only need to check the first factor because it is the one with
#    ### the most levels.
#    if (!(length(levels(dm$flist[[1]])) < length(Y)))
#        stop(paste("Number of levels of a grouping factor for the random effects",
#                   "must be less than the number of observations", sep = "\n"))
#
#    dm$dd["verb"] <- as.integer(verbose)
#    swts <- sqrt(unname(fr$wts))
#    p <- dm$dd[["p"]]
#    n <- length(Y)
#
#    result <- new(Class = "mer",
#               env = new.env(),
#               nlmodel = (~I(x))[[2]],
#               frame = fr$mf,
#               call = call("foo"),      # later overwritten
#               flist = dm$flist,
#               X = fr$X,
#               Zt = dm$Zt,
#               pWt = unname(fr$wts),
#               offset = unname(fr$off),
#               y = unname(Y),
#               Gp = unname(dm$Gp),
#               dims = dm$dd,
#               ST = dm$ST,
#               A = dm$A,
#               Cm = dm$Cm,
#               Cx = if (length(swts)) (dm$A)@x else numeric(0),
#               L = dm$L,
#               deviance = dm$dev,
#               fixef = fr$fixef,
#               ranef = numeric(dm$dd[["q"]]),
#               u = numeric(dm$dd[["q"]]),
#               eta = numeric(n),
#               mu = numeric(n),
#               resid = numeric(n),
#               sqrtrWt = swts,
#               sqrtXWt = as.matrix(swts),
#               RZX = matrix(0, dm$dd[["q"]], p),
#               RX = matrix(0, p, p))
#
#               
##    if (!is.null(stp <- start$STpars) && is.numeric(stp)) {
##        STp <- .Call(mer_ST_getPars, ans)
##        if (length(STp) == length(stp))
##            .Call(mer_ST_setPars, ans, stp)
##    }
##    
##    mer_finalize(result)
#
#result
#}

LogLikeAndScore.C <- function(
  params,
  y, x, z, id,
  w.function,
  cutpoints,
  SampProb,
  SampProbi,
  ProfileCol=NA)
{
  ans <- .Call("llscore",
                params,y, x, z, id*1.0, w.function, 
                cutpoints, SampProb, SampProbi)

  if (!is.na(ProfileCol)) ans@grad[ProfileCol] <- 0
  ans
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
                     cutpoints,                   ## the cutpoints: when w.function="bivar", this is a vector of length 4 that define a central, rectangular region with vertices (x_lower, x_upper, y_lower, y_upper).
                     SampProb,                    ## Sampling probabilities within the sampling strata to be used for ACML
                     SampProbi,                   ## Subject specific sampling probabilities to only be used if doing IPWL.  Note if doing IPWL, only use robcov (robust variances) and not covar
                     ProfileCol=NA){              ## Columns to be held fixed while doing profile likelihood.  It is fixed at its initial value.
 
    out <- nlm(LogLikeAndScore.C, InitVals, y=y, x=x, z=z, id=id, w.function=w.function, cutpoints=cutpoints, SampProb=SampProb,
               SampProbi=SampProbi, ProfileCol=ProfileCol, 
               stepmax=4, iterlim=250, check.analyticals = TRUE) #, print.level=2)

    ## Calculate the observed information and then invert to get the covariance matrix       
    npar        <- length(out$estimate)
    Hessian.eps <- 1e-7
    eps.mtx     <- diag(rep(Hessian.eps, npar))
    grad.at.max <- out$gradient
    ObsInfo     <- matrix(NA, npar, npar)
 
    ## Observed Information
    for (j in 1:npar){
    	temp        <- LogLikeAndScore.C(out$estimate+eps.mtx[j,], y=y, x=x, z=z, id=id, 
                                       w.function=w.function, cutpoints=cutpoints, 
                                       SampProb=SampProb,SampProbi=SampProbi, ProfileCol=ProfileCol)
  	    ObsInfo[j,] <- (attr(temp,"gradient")-grad.at.max)/(Hessian.eps)
  	}
  	cheese <- rep(0, npar*npar)
  	cheeseresult <- .C(     "total_nll_lme",
            beta        = as.double(out$estimate[c(1:(npar-4))]),
            nbeta       = as.integer(npar-4),
            y           = as.double(y),
            x           = as.double(x),
            z           = as.double(z),
            id          = as.double(id),
            n           = as.integer(length(y)),
            wfunction   = as.character(w.function),
            cutpoints   = as.double(cutpoints),
            ncutpoints  = as.integer(length(cutpoints)),
            sampprob    = as.double(SampProb),
            sampprobi   = as.double(SampProbi),
            sigma0      = as.double(exp(out$estimate[npar-3])),
            sigma1      = as.double(exp(out$estimate[npar-2])),
            rho         = as.double((exp(out$estimate[npar-1])-1) / (exp(out$estimate[npar-1])+1)),
            sigmae      = as.double(exp(out$estimate[npar])),
            grad        = as.double(out$gradient),
            cheese      = as.double(cheese) )

    if (!is.na(ProfileCol)){ 
        out$estimate <- out$estimate[-ProfileCol]
        ObsInfo <- ObsInfo[-ProfileCol, -ProfileCol]
        Cheese  <- Cheese[-ProfileCol, -ProfileCol]
    }

    list(Ests=out$estimate, covar=solve(ObsInfo), LogL= -out$minimum,  Code=out$code, robcov=solve(ObsInfo)%*%matrix(cheeseresult$cheese, nrow=npar, byrow=FALSE)%*%solve(ObsInfo))
}

  #####################################################################
 ## The main event
## Linear Mixed-Effects in R with Ascertainment Correction
##
acml <- function(   formula,
                    data,
                    sample.function,
                    sample.cutpoints,
                    sample.prob,
                    sample.probi = NULL,
                    control      = list(),
                    start        = NULL,
                    verbose      = FALSE,
                    subset,
                    na.action,
                    offset,
                    contrasts = NULL,
                    ...)
{
    mc <- match.call()

    stopifnot(length(formula <- as.formula(formula)) == 3)

    fr <- aclmeFrames(mc, formula, contrasts)  # model frame, X, etc.

    # Warn on unused arguments
    largs <- list(...)
    if(length(largs))
        warning("the following '...' arguments have  *not* been used: ",
                sub("^list", "", deparse(largs, control=NULL)))

    # Setup search control
    cv <- do.call(aclmeControl, control)

    if (missing(verbose)) verbose <- cv$msVerbose
    
    # Need an initial guess

    naive <- lmer(formula, data=data)
    beta  <- fixef(naive)
    #beta <- rep(0, 4)
 
    result  <- ACML.LME(y=fr$Y, x=fr$X, z=fr$Z, id=fr$id,
                        w.function=sample.function,  InitVals=append(beta, c(1, 1, 0, 1)),
                        SampProb=sample.prob,        cutpoints=sample.cutpoints, SampProbi=sample.probi)
    
    #FL$dims["mxit"] <- cv$maxIter
    #FL$dims["mxfn"] <- cv$maxFN
    
    # Create the output variable
    #result <- list(fr = fr, FL = FL, start = start, verbose = verbose)
    #
    #if (doFit)
    #{
    #    result <- do.call(acml_finalize, result)
    #    result@call <- mc
    #}
    obj <- new("acmlObj", result$Ests )
    
    attr(obj, "cov")    <- result$cov
    attr(obj, "logl")   <- result$LogL
    attr(obj, "robcov") <- result$robcov
    attr(obj, "data")   <- data
    
    obj
}

setMethod("show", "acmlObj",
    function(object)
    {
        sigma <- sqrt(diag(object@robcov))
        
        upper <- object + qnorm(0.025, lower.tail=FALSE)*sigma
        lower <- object - qnorm(0.025, lower.tail=FALSE)*sigma
        
        l <- length(object)
        cat("Estimate for\n");
        cat("  beta      ",object[1:(l-4)], "\n")
        cat("    lower95 ",lower[1:(l-4)],  "\n")
        cat("    upper95 ",upper[1:(l-4)],  "\n")
        cat("  sigma[0]  ",exp(object[l-3])," (", exp(lower[l-3]), ",", exp(upper[l-3]), ")\n")
        cat("  sigma[1]  ",exp(object[l-2])," (", exp(lower[l-2]), ",", exp(upper[l-2]), ")\n")
        cat("  rho       ",(exp(object[l-1])-1)/(exp(object[l-1])+1),
            " (", (exp(lower[l-1])-1)/(exp(lower[l-1])+1), ",", (exp(upper[l-1])-1)/(exp(upper[l-1])+1),
           ")\n")
        cat("  sigma[e]  ",exp(object[l])," (", exp(lower[l]), ",", exp(upper[l]), ")\n")
    }
)

setMethod("summary", "acmlObj",
    show
)