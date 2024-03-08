\name{GroupByTimeInteraction}
\alias{GroupByTimeInteraction}
\title{Simulated Group by Time Interaction dataset}

\description{
This is a simulated data set with known ``truth''. The dataset consists of 1000 patients and a measured response to a hypothetical drug over the course of 5 months. Further, some ethnicity information is available weather they are of Hispanic origin or not. An expensive genotype analysis exists, and to minimize exposure a biased sample is performed. Patients that are in the outer 20\% of responses for intercept and slope are sampled for the genotyping with a probability of 1. Patients that are in the inner 80\% of responses for intercept and slope are sampled with a probability of 0.25 for genotyping. Thus on average 400 patients of the 1000 would be genotyped.

The model was generated using the following,
\deqn{Y_{ij} = 10 + 0.5 month - 0.5 genotype  + 0.25 month \cdot genotype + 0.5 race + b_{0i} + b_{1i} month + \epsilon_{ij}}
     {Y_ij = 10 + 0.5 month - 0.5 genotype  + 0.25 month * genotype + 0.5 race + b_0i + b_1i * month + e_ij}
where
\deqn{b_{0i} \sim \mathcal{N}(\mu=0, \sigma^2_0=4)}
     {b_0i ~ N(mu=0, sigma^2_0=4)}
\deqn{b_{1i} \sim \mathcal{N}(\mu=0, \sigma^2_1=0.25)}
     {b_1i ~ N(mu=0, sigma^2_1=0.25)}
with a 0.1 correlation them. Further,
\deqn{pr(genotype_i = 1 | race_i) = 0.1 + 0.15*I(race_i=1)}
     {pr(genotype_i = 1 | race_i) = 0.1 + 0.15*I(race_i=1)}

FIXMEFIXMEFIXME Description of bivariate sampling done for genotype. FIXMEFIXMEFIXME
}

\usage{
  data(GroupByTimeInteraction)
}
\value{
   Patient responses, with stratified sampling of genotype in a biased manner based on outcome.
}
\author{ Shawn Garbett }

\examples{
  data(GroupByTimeInteraction)
}