# Experimental design and power analysis for feeding trials. 
# Author: Nathan Spindel

# Set working directory
wd <- "C:/Users/nates/OneDrive/Documents/abalone_feeding_trials"
setwd(wd)

# Load functions
functions <- "/functions.R"
source(paste0(wd, functions))

# Load packages
packages <- c("tidyverse", "pwr", "pwr2", "ez", "BayesFactor", "snowfall", "reshape2")
ipak(packages)

# Expected size and temperature-dependent scaling relationship for whole-organism metabolic rate (Gillooly et al 2001)
# I = M 3/4e -Ei/kT 
# Where:
# I = whole-organism metabolic rate
# M = body mass
# Ei = energy of activation (eV) = approximately 0.6 eV
# k = Boltzmann's constant = 8.6173303e-05 eV/K
# T = temperature in Kelvin

# Store temperature values in Kelvin for relevant °C temperatures:
T.low <- 281.15 # 8 degrees Celsius (trough upwelling)
T.med <- 284.15 # 11 degrees Celsius (medium)
T.high <- 288.15 # 15 degrees Celsius (peak summer)

# Store sets of hypothetical body mass, temperature, and calculated whole-organism metabolic rate values:
mass <- c(400, 800, 1600, 3200) # grams
temperature <- c(8, 11, 15) # °C      ##NOTE: calculation done using Kelvin
metabolic.rate.8.C <- mass^(3/4) * exp(-0.6/(8.6173303e-05 * T.low))
metabolic.rate.11.C <- mass^(3/4) * exp(-0.6/(8.6173303e-05 * T.med))
metabolic.rate.15.C <- mass^(3/4) * exp(-0.6/(8.6173303e-05 * T.high))

# Reformat data for graphing 
dat <- c(metabolic.rate.8.C, metabolic.rate.11.C, metabolic.rate.15.C)
dat.matrix <- matrix(data = dat, nrow = 4, ncol = 3, dimnames = list(mass, temperature))
df <- as.data.frame(dat.matrix)
df.long <- reshape(data = df, idvar = "body mass", v.names = "metabolic rate",ids = row.names(df),
                   times = names(df), timevar = "temp",
                   varying = list(names(df)), direction = "long")
df.long$`body mass` <- as.numeric(df.long$`body mass`)
df.long$`temp` <- as.numeric(df.long$`temp`)
df.long$`temp` <- as.factor(df.long$`temp`)

# Visualize body mass versus whole-organism metabolic rate relationship normalized to 3 relevant temperatures
metabolic.rate.vs.body.size <- ggplot(data = df.long, aes(x = `body mass`, y = `metabolic rate`, shape = `temp`, color = `temp`)) +
  geom_point(size = 3) +
  geom_smooth(method = lm) +
  xlab("Body mass (g)") +
  ylab("Whole-organism metabolic rate") +
  ggtitle("Expected metabolic scaling at three temperatures") 

metabolic.rate.vs.body.size

# Code Gillooly's scaling formula:
# exponent <- -0.6/(8.6173303e-05 * T.low)

# whole.organism.metabolic.rate <- I ~ M^(3/4) * exp(exponent)

# Frequentist Power Analysis. 
# See https://cran.r-project.org/web/packages/pwr/pwr.pdf and https://cran.r-project.org/web/packages/pwr2/pwr2.pdf
# Power calculations for a two factor ANOVA:
a	<- 3 # Number of groups in Factor A (i.e. temperature)
b	<- 4 # Number of groups in Factor B (i.e. body size)
alpha	<- 0.05 # Significant level (Type I error probability)
size.A <- 15 # Sample size per group in Factor A
size.B <- 13 # Sample size per group in Factor B
f.A	<- 0.25 # Effect size of Factor A 
f.B	<- 0.25 # Effect size of Factor B
delta.A	<- NULL # The smallest difference among a groups in Factor A
delta.B	<- NULL # The smallest difference among b groups in Factor B
sigma.A	<- NULL # Standard deviation, i.e. square root of variance in Factor A
sigma.B	<- NULL # Standard deviation, i.e. square root of variance in Factor B

pwr.2way(a=a, b=b, alpha=alpha, size.A=size.A, size.B=size.B, f.A=f.A, f.B=f.B, 
         delta.A=delta.A, delta.B=delta.B, sigma.A=sigma.A, sigma.B=sigma.B)
# Balanced two-way analysis of variance power calculation 
# 
# a = 3
# b = 4
# n.A = 15
# n.B = 15
# sig.level = 0.05
# power.A = 0.8545298
# power.B = 0.8038759
# power = 0.8038759
# 
# NOTE: power is the minimum power among two factors

power.anova.test(groups=2,n=NULL,between.var=1,within.var=2,
                 + sig.level=0.05,power=0.80)

# Power calculations for the general linear model:
# # Large effect size:
# # Assume: p = 2 and R2 = 0.8 , so that f2 = 4
# deg.freedom.numerator <- 2
# deg.freedom.denominator <- NULL
# effect.size <- 4
# type.1.error.probability <- 0.05
# one.minus.type.2.error.probability <- 0.80
# 
# pwr.f2.test(u=deg.freedom.numerator,v=deg.freedom.denominator,f2=effect.size,sig.level=type.1.error.probability, power = one.minus.type.2.error.probability)
# # Multiple regression power calculation 
# # 
# # u = 2
# # v = 3.478442
# # f2 = 4
# # sig.level = 0.05
# # power = 0.8
# # Need n = v + (p + 1) = 4 + (2 + 1) = 7 samples

# Bayesian Power Analysis
# SOURCE1: http://daniellakens.blogspot.com/2016/01/power-analysis-for-default-bayesian-t.html
# SOURCE2: https://datashenanigan.wordpress.com/2016/01/15/speeding-bayesian-power-analysis-t-test-up-with-snowfall/
# The function takes the sample size n, the true effect size D, as well as the effect size of the alternative hypothesis rscaleBF as arguments.
simFun <- function(n, D, rscaleBF){
  library(BayesFactor)
  x <- rnorm(n = n, mean = 0, sd = 1)
  y <- rnorm(n = n, mean = D, sd = 1)
  
  return(exp((ttestBF(x, y, rscale = rscaleBF))@bayesFactor$bf))
}

D<-0.0 #Set the true effect size
n<-180 #Set sample size of your study (number in each group)
nSim<-100000 #Set number of simulations (it takes a while, be patient)
rscaleBF<-sqrt(2)/2 #Set effect size of alternative hypothesis (default = sqrt(2)/2, or 0.707)
threshold<-3 #Set threshold for 'support' - e.g., 3, 10, or 30

# for time-keeping
t0 <- Sys.time()

# initiate a parallel cluster with 12 cpus
sfInit(parallel = T, cpus = 12)

# export the function to the clusters
sfExport("simFun", "n", "D", "rscaleBF")

# execute the code
bf <- sfClusterApplyLB(1:nSim, function(i) simFun(n = n, D = D, rscaleBF = rscaleBF))

# stop the clusters
sfStop()

# print the time it took for the calculation
Sys.time() - t0

# and finally the result
supportH0 <- sum(bf<(1/threshold))/nSim
supportH1 <- sum(bf>threshold)/nSim

hist(log(as.numeric(bf)), breaks=20, main = "Histogram of Bayes Factor", xlab = "Bayes Factor")