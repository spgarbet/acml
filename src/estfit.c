//
// To compile for R use: `R CMD SHLIB estfit.c mvt.f`
// On linux, the link fails, so this follows (Thanks Dexter!)
// For Linux: `gcc -shared -o estfit.so estfit.o mvt.o -lgfortran -lm -L/usr/lib64/R/lib -lR -llapack`
//
//  This now fails on 10.7.4 Lion on a Mac. Arrrgh.
//
// Contains all C fns used in generating and fitting marginalized longitudinal models.
//
// Helpful sites: http://www.netlib.org/lapack/lug/node145.html
// http://www.netlib.org/blas/dgemm.f
// http://www.netlib.org/lapack/double/dgetri.f

// Include the R headers 
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>
#include <R_ext/Parse.h>
#include <R_ext/Utils.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

static const double one      = 1.0;   // coefficients for the matrix multiplication using the blas
static const double zero     = 0.0;
static const double minusone = -1.0;
static const int    izero    = 0;
static const int    ione     = 1;
static const int    itwo     = 2;
static const double neg_half = -0.5;
static const double half     = 0.5;

#define NINVSQRT2  (-0.707106781186547524400844362104849039)
#define LOG2PI     (1.8378770664093454835606594728112352797)
#define LN10       (2.3025850929940456840179914546843642076)
#define INVSQRT2PI (0.3989422804014326779399460599343818685)

double F77_SUB(mvtdst)(const int*,const int*,double*,double*,int*,double*,double*,int*,double*,double*,double*,double*,int*);

// Put here to make linker happy (??)
void   F77_SUB(rndstart)(void) { GetRNGstate();      }
void   F77_SUB(rndend)  (void) { PutRNGstate();      }
double F77_SUB(unifrnd) (void) { return unif_rand(); }

// Inner-most computation, called from multiple places
/* Original R Code
vi.calc <- function(zi, sigma0, sigma1, rho, sigmae)
{
    zi %*% matrix(c(sigma0^2, rho*sigma0*sigma1,rho*sigma0*sigma1,sigma1^2), nrow=2) %*% t(zi) + 
    sigmae*sigmae*diag(length(zi[,2]))
}
*/
void vi_calc(   double* zi,      // Zi matrix (nz x 2)
                int*    n,       // length of dataset
                int*    nz,      // rows in Z
                double* sigma_e, // value of sigma_e
                double* Sigma_i, // 2 x 2 covariance matrix
                double* vi,      // Place to put answer, nz x nz
                double* vi_tmp)  // Temp computation space at least nz x 2
{
    double var_e = (*sigma_e)*(*sigma_e);
    int    i;

    // multiply zi times sigma-i to result in temporary
    F77_CALL(dgemm)("N","N",nz,&itwo,&itwo,&one,zi,n,Sigma_i,&itwo,&zero,vi_tmp,nz);

    // multiply temporary result times t(zi)
    F77_CALL(dgemm)("N","T",nz,nz,&itwo,&one,vi_tmp,nz,zi,n,&zero,vi,nz);

    // Add variance to diagonal components
    for(i=0; i<(*nz); ++i) { vi[i*(*nz)+i] += var_e; }
}

/* Original R Code
subject.ll.lme <- function(yi, xi, beta, vi){ 
	resid <- yi - xi %*% beta
    -(1/2) * (length(xi[,1])*log(2*pi) + log(det(vi)) + t(resid) %*% solve(vi) %*% resid )[1,1]
}
*/
void subject_ll_lme(    double* yi,
                        double* xi,
                        double* beta, 
                        double* vi,
                        int*    n,
                        int*    ni,
                        int*    nbeta,
                        
                        double* ll_lme,  // OUTPUT value
                        double* vi_inv,  // OUTPUT Inverse of vi
                        double* resid,   // OUTPUT: Residual for this id
                        
                        // Required temporary working space
                        double* vi_tmp,  // Same size as vi
                        double* vi2_tmp, // Same size as vi
                        int*    ipvt     // ni long working space
                     )
{ 
    int    job=3;  // assumes binary spec, or should this be 11?
    int    info;
    int    i;
    double det=1.0; // for lapack
    double ldet;


    // Duplicate yi into resid, due to destructive nature of next call
    memcpy(resid, yi, (*ni)*sizeof(double));

    //resid <- yi - xi %*% beta (resid is now resid)
    F77_CALL(dgemv)("N", ni, nbeta, &minusone, xi, n, beta, &ione, &one, resid, &ione);

    // Duplicate vi into vi_inv, due to destructive nature of next call
    memcpy(vi_inv, vi, (*ni)*(*ni)*sizeof(double));

    // LU Decomposition of vi
    F77_CALL(dgetrf)(ni, ni, vi_inv, ni, ipvt, &info);
    if(info < 0) {error("dgetrf illegal value");}
    if(info > 0) {error("dgetrf singular matrix");}

    // Compute determinate of vi from LU decomposition
    for(i=0; i<(*ni); ++i) { det *= vi_inv[ i*(*ni)+i ]; } // printf("vi_inv[%d] = %lf, det=%lf\n", i*(*ni)+i, vi_inv[ i*(*ni)+i ], det);}
    ldet = log(fabs(det));
    //printf("det = %lf\n", det);
    //printf("ldet = %lf\n\n", ldet);

    // vi_inv <- solve(vi)
    F77_CALL(dgetri)(ni, vi_inv, ni, ipvt, vi_tmp, ni, &info);
    if(info < 0) {error("dgetri illegal value");}
    if(info > 0) {error("dgetri singular matrix");}

    // vi_tmp now contains ( t(solve(vi)) * resid)
    F77_CALL(dgemv)("T", ni, ni, &one, vi_inv, ni, resid, &ione, &zero, vi_tmp, &ione);

    // vi2_tmp now contains t(resid) %*% solve(vi) %*% resid 
    F77_CALL(dgemv)("T", ni, &ione, &one, resid, ni, vi_tmp, &ione, &zero, vi2_tmp, &ione);

    // NOTE: vi_inv contains the inverse of vi after this call as an output

    // Final computation on subject log-likelihood
    // -(1/2) * (length(xi[,1])*log(2*pi) +t(resid) %*% solve(vi) %*% resid )[1,1]+log(det)
    *ll_lme = -0.5*(LOG2PI*(double)(*ni) + vi2_tmp[0] + ldet);
}

double dpnorm(double x, double mu, double sigma)
{
    return( 0.5*erfc((x-mu) * NINVSQRT2/sigma) );
}

/* Original R routine
## Ascertainment correction piece for univariate sampling
lci <- function(cutpoints, SampProb, mu_q, sigma_q){ 
	CDFs <- pnorm(c(-Inf, cutpoints, Inf), mu_q, sigma_q)
        sum( SampProb*(CDFs[2:length(CDFs)] - CDFs[1:(length(CDFs)-1)]) )
}
*/
void lci(   double* cutpoints,
            int*    ncutpoints,
            double* sampprob,
            double* mu_q,
            double* sigma_q,
            double* lci)
{
    double  left=0.0;
    double  right;
    int     i;
    
    (*lci) = 0.0;

    for(i=0; i< (*ncutpoints); ++i)
    {
        right = dpnorm(cutpoints[i], *mu_q, *sigma_q);
        (*lci) += sampprob[i]*(right - left);
        left = right;
    }
    right = 1.0;
    (*lci)  += sampprob[(*ncutpoints)]*(right - left);
}


/* Original R routine
## Calculate log of the ascertainment correction for the univariate sampling case
ascertainment.correction <- function(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb){
    vi      <- vi.calc(zi, sigma0, sigma1, rho, sigmae)
    mu      <- xi %*% beta
    mu_q    <- (wi %*% mu)[,1]  
    sigma_q <- sqrt((wi %*% vi %*% t(wi))[1,1])
    log(lci(cutpoints, SampProb, mu_q, sigma_q))
}

*/
void exp_ascertainment_correction_univar(
    int*    n,
    int*    ni,
    int*    npar,
    double* xi,
    double* vi,
    double* wi,
    double* beta,
    double* cutpoints,
    int*    ncutpoints,
    double* sampprob,
    double* vi_tmp,
    double* ac)
{
    double  expac;
    double  mu_q;
    double  sigma_q;

// FIXME, This is duplicated code from the nll case
    // Compute mu    vi_tmp <- xi %*% beta
    F77_CALL(dgemv)("N", ni, npar, &one, xi, n, beta, &ione, &zero, vi_tmp, &ione);

    // Compute mu_q <- (wi %*% mu)[,1]
    mu_q = F77_CALL(ddot)(ni, wi, &ione, vi_tmp, &ione);

    // Compute sigma_q <- sqrt((wi %*% vi %*% t(wi))[1,1])
    // vi_tmp <- vi %*% t(wi)  (vi_tmp is now a vector)
    F77_CALL(dgemv)("N", ni, ni, &one, vi, ni, wi, &ione, &zero, vi_tmp, &ione);
    sigma_q = sqrt(F77_CALL(ddot)(ni,wi,&ione,vi_tmp,&ione));
    
    lci(cutpoints, ncutpoints, sampprob, &mu_q, &sigma_q, &expac);
// End of FIXME
    
    *ac = expac;
}

/* Original R routine
## Ascertainment correction piece for bivariate sampling
lci.bivar <- function(cutpoints, SampProb, mu_q, sigma_q){ 
       (SampProb[1]-SampProb[2])*pmvnorm(
                 lower=c(cutpoints[c(1,3)]), 
                 upper=c(cutpoints[c(2,4)]),
                 mean=mu_q,
                 sigma=sigma_q)[[1]] 
       + SampProb[2]
       #print("blah2")
}
*/
void lci_bivar(
    double*  cutpoints,
    int*     ncutpoints,
    double*  sampprob,
    double*  mu_q,
    double*  sigma_q,
    double*  lci)
{
    double   lower[2];
    double   upper[2];
    int      infin[2];
    double   err;
    int      maxpts=25000;
    double   pnorm;
    int      inform;
    double   abseps = 0.001;
    double   releps = 0.0;
    double   corr;
    double   delta[2];
    
    delta[0] = 0.0;
    delta[1] = 0.0;
    
    lower[0] = (cutpoints[0]-mu_q[0])/sqrt(sigma_q[0]);
    lower[1] = (cutpoints[2]-mu_q[1])/sqrt(sigma_q[3]);
    upper[0] = (cutpoints[1]-mu_q[0])/sqrt(sigma_q[0]);
    upper[1] = (cutpoints[3]-mu_q[1])/sqrt(sigma_q[3]);
    infin[0] = 2; // limits are [LOWER(I), UPPER(I)]
    infin[1] = 2; // limits are [LOWER(I), UPPER(I)]
    corr=sigma_q[1]/(sqrt(sigma_q[0])*sqrt(sigma_q[3]));
    
    F77_CALL(mvtdst)(
        &itwo,
        &izero,
        lower,
        upper,
        infin,
        &corr,
        delta,
        &maxpts,
        &abseps,
        &releps,
        &err,
        &pnorm,
        &inform);

    if(inform !=0)
    {
        printf("%d",inform);
        error("Unable to compute mvtnorm");
    }

    //printf("mvtdst %lg\n",pnorm);
    *lci = (sampprob[0]-sampprob[1])*pnorm+sampprob[1];
}

/* Original R routine
## Calculate log of the ascertainment correction for the bivariate sampling case
ascertainment.correction.bivar <- function(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb){
    vi      <- vi.calc(zi, sigma0, sigma1, rho, sigmae)
    mu      <- xi %*% beta
    mu_q    <- as.vector(wi %*% mu)
    sigma_q <- wi %*% vi %*% t(wi)
    log((SampProb[1]-SampProb[2])*pmvnorm(lower=c(cutpoints[c(1,3)]), upper=c(cutpoints[c(2,4)]),
                                          mean=mu_q, sigma=sigma_q)[[1]] + SampProb[2])
}
*/
void exp_ascertainment_correction_bivar(
    int*    n,
    int*    ni,
    int*    npar,
    double* xi,
    double* vi,
    double* wi,
    double* beta,
    double* cutpoints,
    int*    ncutpoints,
    double* sampprob,
    double* vi_tmp,
    double* ac)
{
    double  mu_q[2];
    double  sigma_q[4];
    double  expac;
    int     i; 

    // Compute mu    vi_tmp <- xi %*% beta
    F77_CALL(dgemv)("N", ni, npar, &one, xi, n, beta, &ione, &zero, vi_tmp, &ione);

    // Compute mu_q <- (wi %*% mu)
    F77_CALL(dgemv)("N", &itwo, ni, &one, wi, &itwo, vi_tmp, &ione, &zero, mu_q, &ione);

    // Compute sigma_q <- wi %*% vi %*% t(wi)
    // vi_tmp <- vi %*% t(wi)  (vi_tmp is now a vector)
    F77_CALL(dgemm)("N", "T", ni, &itwo, ni, &one, vi, ni, wi, &itwo, &zero, vi_tmp, ni);
    F77_CALL(dgemm)("N", "N", &itwo, &itwo, ni, &one, wi, &itwo, vi_tmp, ni, &zero, sigma_q, &itwo);

    lci_bivar(cutpoints, ncutpoints, sampprob, mu_q, sigma_q, &expac);

//    printf("expac %lg\n",expac);

    *ac = expac;
}

/* Routine to compute the maximum length of any subject's data */
int max_i_length(double *id, R_len_t n)
{
    int ni_max=0;
    int current_id;
    int i;
    int ni;

    i=0;
    do
    {
        // Initialize
        ni = i;
        current_id = id[i];

        // Find next id group
        do {++i;} while (i < n && current_id == id[i]);
        ni = i - ni;
        if(ni > ni_max) {ni_max = ni;}
    } while(i < n);

    return(ni_max);
}

void compute_wi(double*      zi,
                int*         ni,
                int*         n,
                const char** wfunction,
                double*      wi)
{
    int    i;
    double x;
    int    info;
    double wi_tmp[4];
    double work[16];
    int    lwork=16;
    int    ipvt[4];

    if(strcasecmp(*wfunction, "mean")==0)
    {
        // wi <- t(rep(1/ni, ni))
        x = 1.0/((double)(*ni));
        for(i=0; i<(*ni); ++i) { wi[i] = x; }
        return;
    }

    // Must solve linear model for subject variable of sampling

    // wi_tmp <- t(zi[,1:2])%*% zi[,1:2]
    F77_CALL(dgemm)("T","N",&itwo, &itwo, ni, &one, zi, n, zi, n, &zero, wi_tmp, &itwo);

    // LU Decomposition of contents of wi_tmp
    F77_CALL(dgetrf)(&itwo, &itwo, wi_tmp, &itwo, ipvt, &info);
    if(info < 0) {error("dgetrf illegal value (wi)");}
    if(info > 0) {error("dgetrf singular matrix (wi)");}

    // wi_tmp <- solve(wi_tmp)
    F77_CALL(dgetri)(&itwo, wi_tmp, &itwo, ipvt, work, &lwork, &info);
    if(info < 0) {error("dgetri illegal value (wi)");}
    if(info > 0) {error("dgetri singular matrix (wi)");}

    // wi_tmp %*% t(zi)
    F77_CALL(dgemm)("N", "T", &itwo, ni, &itwo, &one, wi_tmp, &itwo, zi, n, &zero, wi, &itwo);

    if(strcasecmp(*wfunction, "bivar")==0)
    {
        // wi <- (solve(t(zi[,1:2])%*% zi[,1:2]) %*% t(zi[,1:2]))
        return;
    }
    else if(strcasecmp(*wfunction, "intercept")==0)
    {
        // wi<- (solve(t(zi[,1:2])%*% zi[,1:2]) %*% t(zi[,1:2]))[1,]
        F77_CALL(dcopy)(ni, wi, &itwo, wi, &ione); // Copy vector out of wi matrix
        return;
    }
    else if(strcasecmp(*wfunction, "slope")==0)
    {
        // wi<- (solve(t(zi[,1:2])%*% zi[,1:2]) %*% t(zi[,1:2]))[2,]
        F77_CALL(dcopy)(ni, &wi[1], &itwo, wi, &ione); // Copy vector out of wi matrix
        return;
    }
    else
    {
        error("Invalid Wi function specified");
    }
}

// Double, General Matrix, trace product of two matrices
// Uses an accumulator of the Haramard product, very efficient computationally
//  y := alpha*trace(A %*% B)
// Using conventions of BLAS for naming of function/variables
void dgtrpr(int*    m,              // INPUT: Number of rows and columns of matrix A and B (required to be same size and square)
            const double* alpha,    // INPUT: scalar alpha
            double* A,              // INPUT: Matrix A
            int*    lda,            // INPUT: the first dimension of A, must be at least max(1,m)
            double* B,              // INPUT: Matrix B
            int*    ldb,            // INPUT: the first dimension of B, must be at least max(1,m)
            double* y               // OUTPUT: a single double of the result
            )
{
    int i;
    int j;
    
    int bump_a = (*lda) - (*m);
    int bump_b = -((*ldb)*(*m))+1;
    
    *y = 0.0;
    for(i=0; i<(*m); ++i)
    {
        for(j=0; j<(*m); ++j)
        {
            (*y) += (*A) * (*B);  // This is the core accumulation multiplication
            A += 1;               // Bump pointer to next row in A
            B += (*ldb);          // Bump pointer to next col in B
        }
        A += bump_a; // Next col in A, reset row
        B += bump_b; // Reset to 1st col in B, first row
    }
    (*y) *= (*alpha); // Final scalar
}

void comp_dvk(double* dvk,        // OUTPUT: dvk
              double* zi,         // INPUT: zi
              double* sig_prime,  // INPUT: S' matrix
              int*    n,          // INPUT
              int*    ni,         // INPUT
              double* vi2_tmp     // Working memory
             )
{
    // Compute dkv in vi_tmp
    // vi2_tmp <- zi %*% sig_prime    (ni x 2) 
    F77_CALL(dgemm)("N","N",ni,&itwo,&itwo,&one,zi,n,sig_prime,&itwo,&zero,vi2_tmp,ni);
    // dvk  <- vi2_tmp %*% t(zi)   (ni x ni) i.e. dvk
    F77_CALL(dgemm)("N","T",ni,ni,&itwo,&one,vi2_tmp,ni,zi,n,&zero,dvk,ni);
}

// dvk5 <- zi %*% matrix(c(2*sigma0, rho*sigma1, rho*sigma1, 0),nrow=2) %*% t(zi)
// l5  <- -0.5*(sum(diag(inv.v %*% dvk5)) - t(resid) %*% inv.v %*% dvk5 %*% inv.v %*% resid)[1,1]
void slope(double* grad,       // OUTPUT: Gradient
           double* inv_v,      // INPUT: ni x ni
           double* resid,      // INPUT: ni x 1 (continguous vector)
           int*    n,          // INPUT
           int*    ni,         // INPUT
           double* vi_tmp,     // INPUT: dvk, Destructive as Working memory
           double* vi2_tmp)    // Working memory
{
    // Compute first product, writes to grad
    // grad <- -0.5*(sum(diag(inv_v %*% vi_tmp)))
    dgtrpr(ni, &neg_half, inv_v, ni, vi_tmp, ni, grad);

    // 0.5 * t(resid) %*% inv.v %*% dvk5 %*% inv.v %*% resid)
    // vi2_tmp <- vi_tmp %*% inv_v
    F77_CALL(dgemm)("N","N",ni,ni,ni,&one,vi_tmp,ni,inv_v,ni,&zero,vi2_tmp,ni);

    // vi_tmp <- vi2_tmp %*% resid
    F77_CALL(dgemv)("N",ni,ni,&one,vi2_tmp,ni,resid,&ione,&zero,vi_tmp,&ione);

    // vi2_tmp <- inv_v %*% vi_tmp
    F77_CALL(dgemv)("N",ni,ni,&one,inv_v,ni,vi_tmp,&ione,&zero,vi2_tmp,&ione);

    // grad += 0.5 * t(resid) %*% vi2_tmp
    (*grad) += 0.5*F77_CALL(ddot)(ni, resid, &ione, vi2_tmp, &ione);
}

double ddnorm(double x, double mu, double sigma)
{
    return(INVSQRT2PI * exp(-0.5 * pow((x-mu)/sigma, 2) ) / sigma);
}

// zi %*% sigma_p %*% t(zi)
void comp_ac_dvk(double* dvk,  // vi_tmp
                 double* zi,
                 int*    n,
                 int*    ni, 
                 double* sigma_p,
                 double* vi2_tmp)
{
    // sigma_p %*% t(zi)
    F77_CALL(dgemm)("N","T",&itwo,ni,&itwo,&one,sigma_p,&itwo,zi,n,&zero,vi2_tmp,&itwo);

    // zi %*% sigma_p %*% t(zi)
    F77_CALL(dgemm)("N","N",ni,ni,&itwo,&one,zi,n,vi2_tmp,&itwo,&zero,dvk,ni);
}

// wi %*% dvk %*% t(wi)
void ac_slope(double* correction,
              double* wi,
              double* dvk,  // From vi_tmp
              int*    ni,
              double* f_alpha_k,
              double* vi2_tmp)
{
    F77_CALL(dgemv)("N", ni, ni, &one, dvk, ni, wi, &ione, &zero, vi2_tmp, &ione);

    *correction = (*f_alpha_k) * F77_CALL(ddot)(ni, wi, &ione, vi2_tmp, &ione);
}

/* Original R Code
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
*/
void ascertainment_gradient_corr_univar(    double* correction,  // OUTPUT: gradient vector correction
                                            double* yi,          // INPUT: ni length vector
                                            double* xi,          // INPUT: ni x npar matrix
                                            double* zi,          // INPUT: ni x 2 matrix
                                            double* wi,          // INPUT: 2 x ni
                                            double* vi,          // INPUT: ni x ni matrix
                                            double* beta,        // INPUT: npar length vector
                                            int*    n,           // INPUT: total single dimension length in memory
                                            int*    ni,          // INPUT: number of individual points 
                                            int*    npar,        // INPUT: number of parameters in beta
                                            double* sigma0,      // INPUT: sigma_0
                                            double* sigma1,      // INPUT: sigma_1
                                            double* rho,         // INPUT: rho
                                            double* sigmae,      // INPUT: sigma_e
                                            double* cutpoints,   // INPUT: ncutpoints length vector
                                            int*    ncutpoints,  // INPUT: Number of cutpoints
                                            double* samp_prob,   // INPUT: Sampling probilities (ncutpoints+1) length vector
                                            
                                            double* vi_tmp,      // Working memory
                                            double* vi2_tmp      // Working memory
                                        )
{
    double  mu_q;
    double  sigma_q;
    double  expac;
    double  f_alpha_k;
    double  sigma_prime[4];
    double  spf;
    double  tmp;
    int     i;
    double  sigma_p[4];

// FIXME, This is duplicated code from the nll case
    // Compute mu    vi_tmp <- xi %*% beta
    F77_CALL(dgemv)("N", ni, npar, &one, xi, n, beta, &ione, &zero, vi_tmp, &ione);

    // Compute mu_q <- (wi %*% mu)[,1]
    mu_q = F77_CALL(ddot)(ni, wi, &ione, vi_tmp, &ione);

    // Compute sigma_q <- sqrt((wi %*% vi %*% t(wi))[1,1])
    // vi_tmp <- vi %*% t(wi)  (vi_tmp is now a vector)
    F77_CALL(dgemv)("N", ni, ni, &one, vi, ni, wi, &ione, &zero, vi_tmp, &ione);
    sigma_q = sqrt(F77_CALL(ddot)(ni,wi,&ione,vi_tmp,&ione));

    lci(cutpoints, ncutpoints, samp_prob, &mu_q, &sigma_q, &expac);
// END OF FIXME

    for(i=0, spf=0.0, f_alpha_k=0.0; i<(*ncutpoints); ++i)
    {
        tmp = (samp_prob[i] - samp_prob[i+1])*ddnorm(cutpoints[i], mu_q, sigma_q);
        spf += tmp;
        f_alpha_k += tmp*(cutpoints[i]-mu_q);
    }
    spf /= expac; // sum(p*f) / l

    f_alpha_k /= 2.0 * expac * sigma_q * sigma_q;

    // First part of correction: d_li_beta <- (wi %*% xi) * spf
    F77_CALL(dgemv)("T", ni, npar, &spf, xi, n, wi, &ione, &zero, correction, &ione);
    
    /*
    a5        <- (wi %*% zi %*% matrix(c(2*sigma0, rho*sigma1,    rho*sigma1,     0),        nrow=2) %*% t(zi) %*% t(wi))[1,1]
    a6        <- (wi %*% zi %*% matrix(c(0,        rho*sigma0,    rho*sigma0,     2*sigma1), nrow=2) %*% t(zi) %*% t(wi))[1,1]
    a7        <- (wi %*% zi %*% matrix(c(0,        sigma0*sigma1, sigma0*sigma1,  0),        nrow=2) %*% t(zi) %*% t(wi))[1,1]
    a8        <- (wi %*% (2 * sigmae * diag(length(yi))) %*% t(wi))[1,1]
    */
 
    sigma_p[0] = 2*(*sigma0);
    sigma_p[1] = (*rho)*(*sigma1);
    sigma_p[2] = sigma_p[1];
    sigma_p[3] = 0.0;
    comp_ac_dvk(vi_tmp,zi,n,ni, sigma_p, vi2_tmp);
    ac_slope(&correction[*npar],wi,vi_tmp,ni,&f_alpha_k,vi2_tmp);

    sigma_p[0] = 0.0;
    sigma_p[1] = (*rho)*(*sigma0);
    sigma_p[2] = sigma_p[1];
    sigma_p[3] = 2.0*(*sigma1);
    comp_ac_dvk(vi_tmp,zi,n,ni, sigma_p, vi2_tmp);
    ac_slope(&correction[*npar+1],wi,vi_tmp,ni,&f_alpha_k,vi2_tmp);

    sigma_p[0] = 0.0;
    sigma_p[1] = (*sigma0)*(*sigma1);
    sigma_p[2] = sigma_p[1];
    sigma_p[3] = 0.0;
    comp_ac_dvk(vi_tmp,zi,n,ni, sigma_p, vi2_tmp);
    ac_slope(&correction[*npar+2],wi,vi_tmp,ni,&f_alpha_k,vi2_tmp);

    memset(vi_tmp, 0, sizeof(double)*(*ni)*(*ni));
    for(i=0; i<(*ni); ++i) {vi_tmp[i+i*(*ni)] = 2.0*(*sigmae);}
    ac_slope(&correction[*npar+3],wi,vi_tmp,ni,&f_alpha_k,vi2_tmp);
}

/*
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
*/
void ascertainment_gradient_corr_bivar( double* correction,  // OUTPUT: gradient vector correction
                                        double* yi,          // INPUT: ni length vector
                                        double* xi,          // INPUT: ni x npar matrix
                                        double* zi,          // INPUT: ni x 2 matrix
                                        double* wi,          // INPUT: 2 x ni
                                        double* vi,          // INPUT: ni x ni matrix
                                        double* beta,        // INPUT: npar length vector
                                        int*    n,           // INPUT: total single dimension length in memory
                                        int*    ni,          // INPUT: number of individual points 
                                        int*    npar,        // INPUT: number of parameters in beta
                                        double* sigma0,      // INPUT: sigma_0
                                        double* sigma1,      // INPUT: sigma_1
                                        double* rho,         // INPUT: rho
                                        double* sigmae,      // INPUT: sigma_e
                                        double* cutpoints,   // INPUT: ncutpoints length vector
                                        int*    ncutpoints,  // INPUT: Number of cutpoints
                                        double* samp_prob,   // INPUT: Sampling probilities (ncutpoints+1) length vector
                                    
                                        double* vi_tmp,      // Working memory
                                        double* vi2_tmp      // Working memory
                                        )
{
    double   eps = 1e-6;        // stepsize for numerical estimation
    double   lower[2];
    double   upper[2];
    int      infin[2];
    double   err;
    int      maxpts=25000;
    double   start;
    int      inform;
    double   abseps = 1e-6;
    double   releps = 0.0;
    double   corr;
    double   delta[2];
    double   lci;

    double   mu_q[2];
    double   sigma_q[4];
    double   mu_q_tmp[2];
    double   sigma_q_tmp[4];
    int      i,j;
    double   Sigma_perturbed[4];
    double   param_perturbed[100];    // Spare working space, FIXME limits to 100 fixed effects
    
    // FIXME, This is duplicated code from the bivar case
    // Compute mu    vi_tmp <- xi %*% beta
    F77_CALL(dgemv)("N", ni, npar, &one, xi, n, beta, &ione, &zero, vi_tmp, &ione);

    // Compute mu_q <- (wi %*% mu)
    F77_CALL(dgemv)("N", &itwo, ni, &one, wi, &itwo, vi_tmp, &ione, &zero, mu_q, &ione);

    // Compute sigma_q <- wi %*% vi %*% t(wi)
    // vi_tmp <- vi %*% t(wi)  (vi_tmp is now a vector)
    F77_CALL(dgemm)("N", "T", ni, &itwo, ni, &one, vi, ni, wi, &itwo, &zero, vi_tmp, ni);
    // sigma_q <- wi %*% vi_tmp
    F77_CALL(dgemm)("N", "N", &itwo, &itwo, ni, &one, wi, &itwo, vi_tmp, ni, &zero, sigma_q, &itwo);
    // END OF FIXME
    
    //start   <- pmvnorm(lower=c(cutpoints[c(1,3)]), upper=c(cutpoints[c(2,4)]), mean=mu_q, sigma=sigma_q)[[1]]
    
    delta[0] = 0.0;
    delta[1] = 0.0;
    
    lower[0] = (cutpoints[0]-mu_q[0])/sqrt(sigma_q[0]);
    lower[1] = (cutpoints[2]-mu_q[1])/sqrt(sigma_q[3]);
    upper[0] = (cutpoints[1]-mu_q[0])/sqrt(sigma_q[0]);
    upper[1] = (cutpoints[3]-mu_q[1])/sqrt(sigma_q[3]);
    infin[0] = 2; // limits are [LOWER(I), UPPER(I)]
    infin[1] = 2; // limits are [LOWER(I), UPPER(I)]
    corr=sigma_q[1]/(sqrt(sigma_q[0])*sqrt(sigma_q[3]));
    
    F77_CALL(mvtdst)(
        &itwo,
        &izero,
        lower,
        upper,
        infin,
        &corr,
        delta,
        &maxpts,
        &abseps,
        &releps,
        &err,
        &start,
        &inform);

    if(inform !=0)
    {
        printf("%d",inform);
        error("Unable to compute mvtnorm (2)");
    }

    /* 
    // for this bivariate sampling case, right now we calculate gradients numerically
    eps.mtx <- diag(c(rep(eps,npar))) 
    */

    /*
    for (rr in 1:npar){ */
    for(i=0; i<(*npar+4); ++i)
    {
        /* par.new      <- param+eps.mtx[rr,] */
        for(j=0; j<(*npar); ++j) {param_perturbed[j] = beta[j];}
        param_perturbed[(*npar)]   = *sigma0;
        param_perturbed[(*npar)+1] = *sigma1;
        param_perturbed[(*npar)+2] = *rho;
        param_perturbed[(*npar)+3] = *sigmae;
        param_perturbed[i] += eps;
        
        /* vi.tmp       <- vi.calc(zi, par.new[(npar-3)], par.new[(npar-2)], par.new[(npar-1)], par.new[npar]) */
        Sigma_perturbed[0] = (param_perturbed[(*npar)])*(param_perturbed[(*npar)]);
        Sigma_perturbed[1] = (param_perturbed[(*npar)+2])*(param_perturbed[(*npar)])*(param_perturbed[(*npar)+1]);
        Sigma_perturbed[2] = Sigma_perturbed[1];
        Sigma_perturbed[3] = (param_perturbed[(*npar)+1])*(param_perturbed[(*npar)+1]);
        vi_calc(zi, n, ni, &param_perturbed[(*npar)+3], Sigma_perturbed, vi_tmp, vi2_tmp);
        
        /* mu.tmp       <- xi %*% par.new[1:(npar-4)] */
        // vi2_tmp <- xi %*% beta
        F77_CALL(dgemv)("N", ni, npar, &one, xi, n, param_perturbed, &ione, &zero, vi2_tmp, &ione);
        
        /*  mu_q.tmp     <- as.vector(wi %*% mu.tmp) */
        // Compute mu_q_tmp <- (wi %*% vi2_tmp)
        F77_CALL(dgemv)("N", &itwo, ni, &one, wi, &itwo, vi2_tmp, &ione, &zero, mu_q_tmp, &ione);
        
        /*  sigma_q.tmp  <- wi %*% vi.tmp %*% t(wi) */
        // vi2_tmp <- vi_tmp %*% t(wi)  (vi_tmp is now a vector)
        F77_CALL(dgemm)("N", "T", ni, &itwo, ni, &one, vi_tmp, ni, wi, &itwo, &zero, vi2_tmp, ni);
        // sigma_q <- wi %*% vi2_tmp
        F77_CALL(dgemm)("N", "N", &itwo, &itwo, ni, &one, wi, &itwo, vi2_tmp, ni, &zero, sigma_q_tmp, &itwo);
        
        /* new.area <- c(new.area, pmvnorm(lower=c(cutpoints[c(1,3)]), upper=c(cutpoints[c(2,4)]), mean=mu_q.tmp, sigma=sigma_q.tmp)[[1]]) */
        lower[0] = (cutpoints[0]-mu_q_tmp[0])/sqrt(sigma_q_tmp[0]);
        lower[1] = (cutpoints[2]-mu_q_tmp[1])/sqrt(sigma_q_tmp[3]);
        upper[0] = (cutpoints[1]-mu_q_tmp[0])/sqrt(sigma_q_tmp[0]);
        upper[1] = (cutpoints[3]-mu_q_tmp[1])/sqrt(sigma_q_tmp[3]);
        corr=sigma_q_tmp[1]/(sqrt(sigma_q_tmp[0])*sqrt(sigma_q_tmp[3]));

        F77_CALL(mvtdst)(
            &itwo,
            &izero,
            lower,
            upper,
            infin,
            &corr,
            delta,
            &maxpts,
            &abseps,
            &releps,
            &err,
            &correction[i],
            &inform);

        if(inform !=0)
        {
            printf("%d",inform);
            error("Unable to compute mvtnorm (3)");
        }
    }
    
    // Deriv <- (new.area-start)/eps
    // out <- (SampProb[1]-SampProb[2])*Deriv / lci.bivar(cutpoints, SampProb, mu_q, sigma_q)
    
    lci_bivar(cutpoints, ncutpoints, samp_prob, mu_q, sigma_q, &lci);


    for(i=0; i<(*npar+4); ++i)
    {
        correction[i] -= start;
        correction[i] /= eps;
        correction[i] *= (samp_prob[0] - samp_prob[1]);
        correction[i] /= lci;
    }
}


/* Original R Code
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
*/
void subject_gradient_ll_lme(   double* sub_grad,// OUTPUT: Gradient for this individual 
                                double* resid,   // INPUT: Residual for this id
                                double* vi,      // INPUT: vi
                                double* vi_inv,  // INPUT: solve(vi)
                                double* xi,      // INPUT: xi
                                double* yi,      // INPUT: yi
                                double* zi,      // INPUT: zi
                                int*    n,       // INPUT: length of dataset
                                int*    ni,      // INPUT: rows in dataset
                                int*    npar,    // INPUT: number of parameters
                                double* sigma0,  // INPUT
                                double* sigma1,  // INPUT
                                double* rho,     // INPUT
                                double* sigmae,  // INPUT
                                
                                // Required temporary working space
                                double* vi_tmp,  // Same size as vi
                                double* vi2_tmp  // Same size as vi
                            )
{
    int    i;
    double sig_prime[4];
    
    // l14 <- t(xi) %*% inv.v %*% resid # for Beta
    // multiply t(xi) times vi_inv and store result in vi_tmp
    // vi_tmp <- inv.v %*% resid
    F77_CALL(dgemv)("N", ni, ni, &one, vi_inv, ni, resid, &ione, &zero, vi_tmp, &ione);
    // sub_grad <- t(xi) %*% vi_tmp 
    F77_CALL(dgemv)("T", ni, npar, &one, xi, n, vi_tmp, &ione, &zero, sub_grad, &ione);
    
    // dvk5 <- zi %*% matrix(c(2*sigma0, rho*sigma1, rho*sigma1, 0),nrow=2) %*% t(zi)
    // l5  <- -0.5*(sum(diag(inv.v %*% dvk5)) - t(resid) %*% inv.v %*% dvk5 %*% inv.v %*% resid)[1,1]
    sig_prime[0] = 2.0 * (*sigma0);
    sig_prime[1] = (*rho) * (*sigma1);
    sig_prime[2] = sig_prime[1];
    sig_prime[3] = 0.0;
    comp_dvk(vi_tmp, zi, sig_prime, n, ni, vi2_tmp);
    slope(&sub_grad[*npar], vi_inv, resid, n, ni, vi_tmp, vi2_tmp);

    // dvk6 <- zi %*% matrix(c(0, rho*sigma0, rho*sigma0, 2*sigma1),nrow=2) %*% t(zi)
    // l6  <- -0.5*(sum(diag(inv.v %*% dvk6)) - t(resid) %*% inv.v %*% dvk6 %*% inv.v %*% resid)[1,1]
    sig_prime[0] = 0.0;
    sig_prime[1] = (*rho) * (*sigma0);
    sig_prime[2] = sig_prime[1];
    sig_prime[3] =  2.0 * (*sigma1);
    comp_dvk(vi_tmp, zi, sig_prime, n, ni, vi2_tmp);
    slope(&sub_grad[*npar+1], vi_inv, resid, n, ni, vi_tmp, vi2_tmp);

    // dvk7 <- zi %*% matrix(c(0, sigma0*sigma1, sigma0*sigma1, 0),nrow=2) %*% t(zi)
    // l7  <- -0.5*(sum(diag(inv.v %*% dvk7)) - t(resid) %*% inv.v %*% dvk7 %*% inv.v %*% resid)[1,1]
    sig_prime[0] = 0.0;
    sig_prime[1] = (*sigma0) * (*sigma1);
    sig_prime[2] = sig_prime[1];
    sig_prime[3] = 0.0;
    comp_dvk(vi_tmp, zi, sig_prime, n, ni, vi2_tmp);
    slope(&sub_grad[*npar+2], vi_inv, resid, n, ni, vi_tmp, vi2_tmp);

    // dvk8 <- 2 * sigmae  * diag(length(yi))
    // l8  <- -0.5*(sum(diag(inv.v %*% dvk8)) - t(resid) %*% inv.v %*% dvk8 %*% inv.v %*% resid)[1,1]
    memset(vi_tmp, 0, sizeof(double)*(*ni)*(*ni));
    for(i=0; i<(*ni); ++i) {vi_tmp[i+i*(*ni)] = 2.0*(*sigmae);}
    slope(&sub_grad[*npar+3], vi_inv, resid, n, ni, vi_tmp, vi2_tmp);
}


/* Original R Code
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
*/
double total_nll_lme(double*      beta,
                     int*         nbeta,
                     double*      y, 
                     double*      x, 
                     double*      z,
                     double*      id,
                     int*         n,
                     const char** wfunction, 
                     double*      cutpoints,
                     int*         ncutpoints,
                     double*      sampprob,
                     double*      sampprobi,
                     double*      sigma0,
                     double*      sigma1,
                     double*      rho,
                     double*      sigmae,
                     double*      grad,
                     double*      cheese)
{
    double   total=0.0;
    double   current_id;
    double   iprwi;
    int      i;
    int      ni;
    int      ni_max;
    int      ni_max2;
    double*  xi;
    double*  zi;
    double*  vi;
    double*  vi_inv;
    double*  vi_tmp;
    double*  resid;
    double*  grad_tmp;
    double*  grad_tmp2;
    double*  vi2_tmp;
    double*  wi;
    double*  yi;
    int*     ipvt;
    double   Sigma_i[4];
    int      j;
    double   ac; //acertainment correction1
    double   ll_lme;
    int      npar;

    npar = (*nbeta)+4;
    
    // Set gradient to zero
    memset(grad, 0, npar * sizeof(double));

    // Create Sigma_i on call stack
    Sigma_i[0] = (*sigma0)*(*sigma0);
    Sigma_i[1] = (*rho)*(*sigma0)*(*sigma1);
    Sigma_i[2] = (*rho)*(*sigma0)*(*sigma1);
    Sigma_i[3] = (*sigma1)*(*sigma1);

    // Find Longest set of measurements for an individual
    ni_max  = max_i_length(id, *n);
    ni_max2 = ni_max*ni_max;

    // ONE Malloc to rule them all -- malloc is a VERY expensive operation, minimizing it's usage by pointer tricks
    // resid   => ni_max x 1
    // vi2_tmp  => ni_max x ni_max
    // vi       => ni_max x ni_max
    // vi_inv   => ni_max x ni_max
    // vi_tmp   => ni_max x ni_max
    // wi       => ni_max x 2
    // ipvt     => ni_max x 1
    // grad_tmp => 2*npar
    // grad_tmp2 => npar
    // TOTAL    => (4*ni_max^2+3*ni_max)  doubles + ni_max ints
    resid  = (double *)malloc((4*ni_max2+3*ni_max + 3*npar)*sizeof(double)+ni_max*sizeof(int));
    vi2_tmp = &resid[ni_max];
    vi      = &vi2_tmp[ni_max2];
    vi_inv  = &vi[ni_max2];
    vi_tmp  = &vi_inv[ni_max2];
    wi      = &vi_tmp[ni_max2];
    grad_tmp= &wi[ni_max*2];
    grad_tmp2=&grad_tmp[npar]; // For local cheese temporary
    ipvt    = (int *)&grad_tmp2[npar];

    // Zero the cheese
    if(cheese != 0)
    {
        memset(cheese, 0, npar * npar * sizeof(double));
    }

    // Loop through individuals
    i = 0; // i holds global offset in data frame, NOT individual number
    do
    {
        // Use Pointers to find them
        current_id = id[i];     // Get current id, i is the offset in the overall data
        xi         = &x[i];     // xi is offset in x
        zi         = &z[i];     // zi is offset in z
        
        //printf("current_id = %lf\n", current_id);

        // Find beginning of next id group
        ni = i;                                         // Save position to compute size
        do {++i;} while (i < (*n) && current_id == id[i]); // Find next id or end of data
        ni = i - ni;                                    // How long is this group?

        // Functions to bring them all
        compute_wi(zi, &ni, n, wfunction, wi);
        vi_calc(zi, n, &ni, sigmae, Sigma_i, vi, vi_tmp);

        // In the darkness total them
        yi = &y[i-ni];
        subject_ll_lme(yi, xi, beta, vi, n, &ni, nbeta, &ll_lme, vi_inv, resid, vi_tmp, vi2_tmp, ipvt);
        total += ll_lme / sampprobi[i-ni];
        
        subject_gradient_ll_lme(grad_tmp,
                                resid, vi, vi_inv,
                                xi, yi, zi,
                                n, &ni, nbeta, 
                                sigma0, sigma1, rho, sigmae,
                                vi_tmp, vi2_tmp);

        for(j=0; j<npar; ++j) { grad_tmp2[j] = -1.0*grad_tmp[j] / sampprobi[i-ni]; }

        if(strcasecmp(*wfunction, "bivar")==0)
        {
            exp_ascertainment_correction_bivar(n, &ni, nbeta, xi, vi, wi, beta, cutpoints, ncutpoints, sampprob,vi_tmp, &ac);
            total -= log(ac);
            ascertainment_gradient_corr_bivar(grad_tmp,
                                              yi, xi, zi, wi, vi,
                                              beta, n, &ni, nbeta,
                                              sigma0, sigma1, rho, sigmae,
                                              cutpoints, ncutpoints, sampprob, vi_tmp, vi2_tmp);
            for(j=0; j<npar; ++j)
            {
                grad_tmp2[j] += grad_tmp[j];  // Apply Correction into grad2
                grad[j]      += grad_tmp2[j]; // Sum gradient
            }
        }
        else
        {
            // Adjust Log-Likelihood
            exp_ascertainment_correction_univar(n, &ni, nbeta, xi, vi, wi, beta, cutpoints, ncutpoints, sampprob, vi_tmp, &ac);
            total -= log(ac);

            // Adjust gradient of log-likelihood
            ascertainment_gradient_corr_univar(grad_tmp, yi, xi, zi, wi, vi, beta,
                                               n, &ni, nbeta,
                                               sigma0, sigma1, rho, sigmae,
                                               cutpoints, ncutpoints, sampprob, 
                                               vi_tmp, vi2_tmp);
            for(j=0; j<npar; ++j)
            {
                grad_tmp2[j] -= grad_tmp[j];  // Apply Correction into grad2
                grad[j]      += grad_tmp2[j]; // Sum gradient
             }
        }
/*        if(cheese != 0 ) {
            printf("gradi[id=%lf] = %lf\n", current_id, grad[0]);
        } */
        if(cheese != 0)
        {
            /* Gradi[(n.par-3)] <- Gradi[(n.par-3)]*exp(param.vec[(n.par-3)])
               Gradi[(n.par-2)] <- Gradi[(n.par-2)]*exp(param.vec[(n.par-2)])
               Gradi[(n.par-1)] <- Gradi[(n.par-1)]*2*exp(param.vec[(n.par-1)])/((exp(param.vec[(n.par-1)])+1)^2)
               Gradi[n.par]     <- Gradi[n.par]  *exp(param.vec[n.par])
            */
            grad_tmp2[(*nbeta)  ] *= (*sigma0);
            grad_tmp2[(*nbeta)+1] *= (*sigma1);
            grad_tmp2[(*nbeta)+2] *= 2.0*(1+(*rho))/(1-(*rho))/pow((1+(*rho))/(1-(*rho)) +1.0,2.0);
            grad_tmp2[(*nbeta)+3] *= (*sigmae);
//            //cheese  <- cheese + outer(Gradi, Gradi)
            F77_CALL(dger)(&npar, &npar, &one, grad_tmp2, &ione, grad_tmp2, &ione, cheese, &npar);

        }
    } while(i < (*n)); // Loop while there is more groups to process (in Middle Earth)

    free(resid);   // Throw malloc in the crack of doom

    return(-1.0*total);  // The total!!! The hobbits live
}

double detail_nll_lme(double*     beta,
                      int         npar,
                      double*     y, 
                      double*     x, 
                      double*     z,
                      double*     id,
                      int         n,
                      const char* wfunction, 
                      double*     cutpoints,
                      int         ncutpoints,
                      double*     sampprob,
                      double*     sampprobi,
                      double      sigma0,
                      double      sigma1,
                      double      rho,
                      double      sigmae,
                      double*     scores,
                      double*     corrections)
{
    double   total=0.0;
    double   current_id;
    double   iprwi;
    int      i;
    int      ni;
    int      ni_max;
    int      ni_max2;
    int      n_id;
    double*  xi;
    double*  zi;
    double*  vi;
    double*  vi_inv;
    double*  vi_tmp;
    double*  resid;
    double*  grad_tmp;
    double*  vi2_tmp;
    double*  wi;
    double*  yi;
    int*     ipvt;
    double   Sigma_i[4];
    int      j;
    double   ac; //acertainment correction1
    double   ll_lme;

    // Create Sigma_i on call stack
    Sigma_i[0] = sigma0*sigma0;
    Sigma_i[1] = rho*sigma0*sigma1;
    Sigma_i[2] = rho*sigma0*sigma1;
    Sigma_i[3] = sigma1*sigma1;

    // Find Longest set of measurements for an individual
    ni_max  = max_i_length(id, n);
    ni_max2 = ni_max*ni_max;

    // ONE Malloc to rule them all -- malloc is a VERY expensive operation, minimizing it's usage by pointer tricks
    // resid   => ni_max x 1
    // vi2_tmp  => ni_max x ni_max
    // vi       => ni_max x ni_max
    // vi_inv   => ni_max x ni_max
    // vi_tmp   => ni_max x ni_max
    // wi       => ni_max x 2
    // ipvt     => ni_max x 1
    // grad_tmp => 2*(npar+4)
    // TOTAL    => (4*ni_max^2+3*ni_max)  doubles + ni_max ints
    resid  = (double *)malloc((4*ni_max2+3*ni_max + 2*npar + 8)*sizeof(double)+ni_max*sizeof(int));
    vi2_tmp = &resid[ni_max];
    vi      = &vi2_tmp[ni_max2];
    vi_inv  = &vi[ni_max2];
    vi_tmp  = &vi_inv[ni_max2];
    wi      = &vi_tmp[ni_max2];
    grad_tmp= &wi[ni_max*2];
    ipvt    = (int *)&grad_tmp[2*(npar+4)];

    // Loop through individuals
    i = 0; // i holds global offset in data frame, NOT individual number
    n_id = 0; // count of id
    do
    {
        // Use Pointers to find them
        current_id = id[i];     // Get current id
        xi         = &x[i];     // xi is offset in x
        zi         = &z[i];     // zi is offset in z
        
        // Find beginning of next id group
        ni = i;                                         // Save position to compute size
        do {++i;} while (i < n && current_id == id[i]); // Find next id or end of data
        ni = i - ni;                                    // How long is this group?

        // Functions to bring them all
        compute_wi(zi, &ni, &n, &wfunction, wi);
        vi_calc(zi, &n, &ni, &sigmae, Sigma_i, vi, vi_tmp);

        // In the darkness total them
        yi = &y[i-ni];
        subject_ll_lme(yi, xi, beta, vi, &n, &ni, &npar, &ll_lme, vi_inv, resid, vi_tmp, vi2_tmp, ipvt);
        scores[n_id] = ll_lme / sampprobi[i-ni];

        if(strcasecmp(wfunction, "bivar")==0)
        {
            exp_ascertainment_correction_bivar(&n, &ni, &npar, xi, vi, wi, beta, cutpoints, &ncutpoints, sampprob,vi_tmp, &ac);
        }
        else
        {
            exp_ascertainment_correction_univar(&n, &ni, &npar, xi, vi, wi, beta, cutpoints, &ncutpoints, sampprob, vi_tmp, &ac);
        }
        corrections[n_id] =  ac;
        
        n_id++;
    } while(i < n); // Loop while there is more groups to process (in Middle Earth)

    free(resid);   // Throw malloc in the crack of doom
}

SEXP llscore(SEXP s_params,
             SEXP s_y, 
             SEXP s_x, 
             SEXP s_z,
             SEXP s_id,
             SEXP s_wfunction, 
             SEXP s_cutpoints,
             SEXP s_sampprob,
             SEXP s_sampprobi)
{
    // Unconstrained Transform of Vector
    double*     params     = REAL(s_params);
    R_len_t     npar       = length(s_params);
    R_len_t     n          = length(s_y);
    double*     beta       = params;
    double      sigma0     = exp(params[npar - 4]);
    double      sigma1     = exp(params[npar - 3]);
    double      rho        = (exp(params[npar-2])-1) / (exp(params[npar-2])+1);
    double      sigmae     = exp(params[npar - 1]);
    double*     y          = REAL(s_y);
    double*     x          = REAL(s_x);
    double*     z          = REAL(s_z);
    double*     id         = REAL(s_id);
    const char* wfunction  = translateChar(STRING_ELT(s_wfunction, 0));
    double*     cutpoints  = REAL(s_cutpoints);
    R_len_t     ncutpoints = length(s_cutpoints);
    double*     sampprob   = REAL(s_sampprob);
    double*     sampprobi  = REAL(s_sampprobi);
    double*     grad       = malloc(npar * sizeof(double)); // Unavoidable malloc
    int         nbeta      = npar-4;
    
    // Set gradient to zero
    memset(grad, 0, npar * sizeof(double));

    // The actual fit values
    double      out     = total_nll_lme(beta, &nbeta, y, x, z, id, &n, &wfunction,
                                        cutpoints, &ncutpoints, 
                                        sampprob, sampprobi, 
                                        &sigma0, &sigma1, &rho, &sigmae,
                                        grad, 0); 

    // R interface variables
    SEXP        SEXP_retval;
    SEXP        SEXP_grad;

    // Some internal looping
    int         i;

    // Additional chain rule transform of gradient for unconstrained
    grad[npar-4] *= sigma0;
    grad[npar-3] *= sigma1;
    //grad[npar-2] *= 2.0*exp(params[npar-2])/pow(exp(params[(npar-2)])+1.0,2.0);
    grad[npar-2] *= 2.0*(1+rho)/(1-rho)/pow((1+rho)/(1-rho) +1.0,2.0);
    grad[npar-1] *= sigmae;

    // Prepare Return to R
    PROTECT(SEXP_retval = NEW_NUMERIC(1));
    REAL(SEXP_retval)[0] = out;

    PROTECT(SEXP_grad = NEW_NUMERIC(npar));

    // Copy gradient into output
    for(i=0; i < npar; ++i) {REAL(SEXP_grad)[i] = grad[i];}
    free(grad);  // Free given gradient
    setAttrib(SEXP_retval, install("gradient"), SEXP_grad);

    UNPROTECT(2);   // unprotects SEXP_retval and gradient
    return SEXP_retval;
}

SEXP lldetail(SEXP s_params,
              SEXP s_y, 
              SEXP s_x, 
              SEXP s_z,
              SEXP s_id,
              SEXP s_wfunction, 
              SEXP s_cutpoints,
              SEXP s_sampprob,
              SEXP s_sampprobi)
{
    // Unconstrained Transform of Vector
    double*     params     = REAL(s_params);
    R_len_t     npar       = length(s_params);
    R_len_t     n          = length(s_y);
    double*     beta       = params;
    double      sigma0     = exp(params[npar - 4]);
    double      sigma1     = exp(params[npar - 3]);
    double      rho        = (exp(params[npar-2])-1) / (exp(params[npar-2])+1);
    double      sigmae     = exp(params[npar - 1]);
    double*     y          = REAL(s_y);
    double*     x          = REAL(s_x);
    double*     z          = REAL(s_z);
    double*     id         = REAL(s_id);
    const char* wfunction  = translateChar(STRING_ELT(s_wfunction, 0));
    double*     cutpoints  = REAL(s_cutpoints);
    R_len_t     ncutpoints = length(s_cutpoints);
    double*     sampprob   = REAL(s_sampprob);
    double*     sampprobi  = REAL(s_sampprobi);
    double*     scores; 
    double*     corrections;
    int         i;
    int         current_id;
    int         nid;
    
    // Count individuals
    i = 0; // i holds global offset in data frame, NOT individual number
    nid = 0;
    do
    {
        // Bump count
        ++nid;
        
        for(current_id = id[i];  i < n && current_id == id[i]; ++i);
    } while(i < n); // Loop while there is more groups to process

    scores      = malloc(nid * sizeof(double));
    corrections = malloc(nid * sizeof(double));

    // Set output to zero
    memset(scores,      0, nid * sizeof(double));
    memset(corrections, 0, nid * sizeof(double));

    // The actual fit values
    detail_nll_lme(beta, npar-4, y, x, z, id, n, wfunction,
                   cutpoints, ncutpoints, 
                   sampprob, sampprobi, 
                   sigma0, sigma1, rho, sigmae,
                   scores, corrections); 

    // R interface variables
    SEXP        SEXP_retval;
    SEXP        SEXP_corrections;

    // Prepare Return to R
    PROTECT(SEXP_retval = NEW_NUMERIC(nid));
    for(i=0; i<nid; ++i) REAL(SEXP_retval)[i] = scores[i];

    PROTECT(SEXP_corrections = NEW_NUMERIC(nid));

    // Copy gradient into output
    for(i=0; i<nid; ++i) {REAL(SEXP_corrections)[i] = corrections[i];}
    free(scores);  // Free given gradient
    free(corrections);
    setAttrib(SEXP_retval, install("corrections"), SEXP_corrections);

    UNPROTECT(2);   // unprotects SEXP_retval and corrections
    return SEXP_retval;
}

static R_CallMethodDef callMethods[] =
{ 
  {"llscore",  (DL_FUNC)&llscore,  9},
  {"lldetail", (DL_FUNC)&lldetail, 9},
  {NULL, NULL, 0}
};

void R_init_MTLVM(DllInfo *info)
{
    R_useDynamicSymbols(info, FALSE);
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

