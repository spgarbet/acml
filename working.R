setwd("~/Projects/schildcrout")
setwd("acml")
setwd("R")
source("parse.R")
source("acml.R")
setwd("../data")
load("GroupByTimeInteraction.rda")



# naive example
setwd("../../validated")
dyn.load("estfit.so")
setwd("../acml")
source("R/parse.R")
source("R/acml.R")

load("data/GroupByTimeInteraction.rda")

# sample
data <- GroupByTimeInteraction
sample <- subset(GroupByTimeInteraction, !is.na(genotype))

library(lme4)
lmer(response ~ month*genotype + race+ (month|patient), data=sample)

target <- 1-sqrt(0.8)
probs  <- c(target/2, 1-target/2)
quantile.intercept <- quantile(data[,2], probs=probs)
quantile.slope     <- quantile(data[,3], probs=probs)
acml(response ~ month*genotype + race + (month|patient), data=sample,
    sample.function="bivar",
    sample.cutpoints=c(quantile.intercept, quantile.slope), 
    sample.prob = c(1.0, 0.25),
    sample.probi = sample$sample.prob
    )

formula <- response ~ month*genotype + race + (month|patient)
sample.function <- "bivar"
sample.cutpoints <- c(quantile.intercept, quantile.slope)
sample.prob <- c(1.0, 0.25)
sample.probi <- sample$sample.prob

mc.stub <- function(   formula,
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
    match.call()
}

mc <- mc.stub(response ~ month*genotype + race + (month|patient), data=sample,
    sample.function="bivar",
    sample.cutpoints=c(quantile.intercept, quantile.slope), 
    sample.prob = c(1.0, 0.25),
    sample.probi = sample$sample.prob)


best <- c(10.495883676,0.363930810,-1.053602789,-0.212373113,0.245794628,1.693826586,-1.134179693,-0.007173192,0.020267984)


######
library(acml)
data(GroupByTimeInteraction)
d <- GroupByTimeInteraction
# Data that was sampled for genotyping 
sample <- subset(GroupByTimeInteraction, !is.na(genotype))
sample.prob <- c(1.0, 0.25)
sample.probi <- sample$sample.prob
sample.function <- "bivar"

# Estimate cutpoints in bivariate space
target <- 1-sqrt(0.8)
probs  <- c(target/2, 1-target/2)
quantile.intercept <- quantile(d[,2], probs=probs)
quantile.slope     <- quantile(d[,3], probs=probs)
sample.cutpoints <- c(quantile.intercept, quantile.slope)


result <- acml(response ~ month*genotype + race + (month|patient),
               data=sample,
               sample.function="bivar",
               sample.cutpoints=sample.cutpoints,
               sample.prob=sample.prob,
               sample.probi=sample$sample.prob
)
