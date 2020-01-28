# Data is stored in table "antidepressants_suicideattempt.xls" and "antidepressants_suicide.xls".
# 
# 
# R Code adopted from Efthimiou(1) :
  
### Start by installing necessary libraries
install.packages("meta")
install.packages("mmeta")# for fitting the beta-binomial model
### Call the libraries
library("meta")
library("mmeta")
library(readxl)

#for suicide_attempts execute the following 7 lines don't execute lines for suicide below

suicide_attempts <- read_excel("antidepressants_suicideattempt.xls")

y1<-suicide_attempts$ADsuicideattempts
n1<-suicide_attempts$ADparticipants
y2<-suicide_attempts$PLCsuicideattempts
n2<-suicide_attempts$PLCparticipants
y1<-as.numeric(y1)
y2<-as.numeric(y2)

#for suicide execute the following 4 lines
suicide <- read_excel("antidepressants_suicide.xls")

y1<-suicide$ADsuicides
n1<-suicide$ADparticipants
y2<-suicide$PLCsuicides
n2<-suicide$PLCparticipants

### Perform the analyses
# inverse-variance odds-ratio with 0.5 continuity correction
OR.IV <- metabin(y1, n1, y2, n2, sm="OR", method = "Inverse", incr=0.5)
print(summary(OR.IV), digits=2)

# inverse-variance odds-ratio with “treatment-arm” continuity correction
OR.IV2 <- metabin(y1, n1, y2, n2, sm="OR", method = "Inverse", incr="TACC")
print(summary(OR.IV2), digits=2)

# Peto odds-ratio fixed effects
MH.Peto <- metabin(y1, n1, y2, n2, sm="OR", method = "Peto")
print(summary(MH.Peto), digits=2)

# Mantel-Haenszel odds-ratio with no continuity correction fixed effects 
MH.OR <- metabin(y1, n1, y2, n2, sm="OR", MH.exact=TRUE)
print(summary(MH.OR), digits=2)

# Mantel-Haenszel odds-ratio with “treatment-arm” continuity correction for forest plot
MH.OR1 <- metabin(y1, n1, y2, n2, sm="OR", MH.exact=F, incr="TACC")
print(summary(MH.OR1), digits=2)

# Mantel-Haenszel odds-ratio with 0.5 continuity correction not in table 
MH.OR2 <- metabin(y1, n1, y2, n2, sm="OR", MH.exact=F, incr=0.5)
print(summary(MH.OR2), digits=2)

# Mantel-Haenszel risk-difference with no continuity correction not in table 
MH.RD<- metabin(y1, n1, y2, n2, sm="RD", MH.exact=TRUE)
print(summary(MH.RD), digits=5)

# Beta-binomial with correlated responses
B=data.frame(y1=y2,y2=y1,n1=n2,n2=n1)
B$studynames=suicide_attempts$study #or
#B$studynames=suicide$study
str(B)
B<-B[c(-3,-6),] #remove NAs for suicide attempts
Beta.Bin<-multipletables(data=B, measure="OR", model="Sarmanov", method="sampling", nsam=1000, alpha = 0.05)
summary(Beta.Bin)
#calculate p-value https://www.bmj.com/content/343/bmj.d2304

#1calculate the standard error:

log(Beta.Bin$overall$CI[2])
log(Beta.Bin$overall$overall)
SE=(Beta.Bin$overall$CI[2]-Beta.Bin$overall$CI[1])/(2*1.96)

SE=(log(Beta.Bin$overall$CI[2])-log(Beta.Bin$overall$CI[1]))/(2*1.96)


#2calculate the test statistic: 

z = Beta.Bin$overall$overall/SE
z = log(Beta.Bin$overall$overall)/SE
#3calculate the P value:
P = exp(-0.717*abs(z)-0.416*abs(z)^2)

# Arcsine difference fixed effect model from Efthimiou
ASD <- metabin(y1, n1, y2, n2, sm="ASD",comb.random=T,method.tau = "DL")
print(summary(ASD), digits=3)

# Arcsine difference for random effects from Hengartner
s.es = escalc(measure="AS", ai=y1, n1i=n1, ci=y2, n2i=n2, drop00=TRUE, slab=Study,data=data.s)
m.s  <- rma(yi = yi, vi = vi, data = s.es, method = "DL", slab = Study)
summary(m.s)

##plots


funnel(OR.IV)
forest(OR.IV,  comb.random = F, text.random = NULL, text.random.w = NULL, label.right = "higher risk in antidepressant", label.left = "higher risk in placebo",
       overall=T,squaresize = 0, studlab = study, overall.hetstat=F)




# Manual Bayesian Meta-Analysis from Hengartner with jags
# https://rpubs.com/mbounthavong/272658

library(rjags,coda)

#### Alternatively, you can enter the model directly into R.
cat("model
    {
    ## Define the likelihood (logit) model:
    for( i in 1:N ) {
    rD[i] ~ dbin( pD[i], nD[i] )
    rP[i] ~ dbin( pP[i], nP[i] )
    
    logit( pP[i] ) <- mu[i]
    logit( pD[i] ) <- mu[i] + delta[i]
    
    mu[i] ~ dnorm( 0.0, 1.0E-5 )
    delta[i] ~ dnorm( d, prec )
    
    }
    
    ## Priors
    d ~ dnorm( 0.0, 0.01) # Variance = 100, therfore use 1/100
    tau <- sqrt( tau.sq )
    tau.sq <- 1 / ( prec )
    prec ~ dgamma( 0.0001, 0.0001 )
    
    ## Outcomes of interest
    OR <- exp( d )
    # prob.OR1 <- step( d )
    
    }", file="model-random.txt")


##########################################

nnn = 5000000 # resamples
n.ch= 8 # number of chains
burnin = 1000000  #burn in
start = burnin +1



N = length(y1) # N Studies
dat <- list("N" = N, "rD" = y1, "rP" =y2  , "nD" =n1, "nP" = n2)  # names list of numbers

inits <- list( "d" = 0, "prec"=1, "delta" = rep(0,N), "mu" =  rep(0,N))
jags.m <- jags.model( file = "model-random.txt", data=dat, inits=inits, n.chains=n.ch, n.adapt=burnin )
params <- c("d","OR")
samps <- coda.samples( jags.m, params, n.iter=nnn )
summary(window(samps, start=start))

plot(window(samps, start=start))

####MetaStan from Günhan et al. https://doi.org/10.1002/jrsm.1370
##needs installation of RStan
##https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
##for general info https://twiecki.io/blog/2015/1
##for further reading implemented  metropolis hastings in stata https://blog.stata.com/2016/11/15/introduction-to-bayesian-statistics-part-2-mcmc-and-the-metropolis-hastings-algorithm/
##
library(rstan)
library(MetaStan)
library(tidyverse)

suicide_attempts<-suicide_attempts[c(-3,-6),] #exclude studies with no info available 
bnhmstan<-suicide_attempts%>%
  meta_stan(nctrl=PLCparticipants,
            ntrt=ADparticipants,
            rctrl=PLCsuicideattempts,
            rtrt=ADsuicideattempts,
            data=.,
            tau_prior_dist="half-normal",
            tau_prior = 0.5, delta=250, chains=4, iter=2000, warmup=1000,adapt_delta=0.95)

print(bnhmstan)
bnhmstan$fit
#https://mc-stan.org/users/documentation/case-studies/rstan_workflow.html

check_hmc_diagnostics(bnhmstan$fit)
#calculate exp() on values for obtaining OR
check_div(bnhmstan$fit)

library("bayesplot")
library("ggplot2")
#diagnostic plots
diagn_trace<-mcmc_trace(bnhmstan$fit, pars = "theta", transformations = "exp") + 
  xlab("Iteration")+ylab("OR")
diagn_dens<-mcmc_dens_overlay(bnhmstan$fit, pars = "theta",transformations = "exp")+xlab("OR")
gridExtra::grid.arrange(diagn_trace,diagn_dens)



#alternetive further diagnostics
library(shinystan)
diagnplot<-as.shinystan(bnhmstan$fit)
launch_shinystan(diagnplot)

bnhmstan_suicide<-suicide%>%meta_stan(nctrl=PLCparticipants,
                                      ntrt=ADparticipants,
                                      rctrl=PLCsuicides,
                                      rtrt=ADsuicides,
                                      data=.,
                                      tau_prior_dist="half-normal",
                                      tau_prior = 0.5, delta=250, chains=4, iter=2000, warmup=1000, adapt_delta=0.95)
print(bnhmstan_suicide)
#calculate exp() on values for OR
bnhmstan_suicide$fit
#diagnostic plots
diagn_trace<-mcmc_trace(bnhmstan_suicide$fit, pars = "theta", transformations = "exp") + 
  xlab("Iteration")+ylab("OR")
diagn_dens<-mcmc_dens_overlay(bnhmstan_suicide$fit, pars = "theta",transformations = "exp")+xlab("OR")
gridExtra::grid.arrange(diagn_trace,diagn_dens)

#https://mc-stan.org/users/documentation/case-studies/rstan_workflow.html
check_hmc_diagnostics(bnhmstan_suicide$fit)


#alternetive further diagnostics
library(shinystan)
diagnplot<-as.shinystan(bnhmstan$fit)
launch_shinystan(diagnplot)



#calculate pooled mean and sd of age from Khan et al. according to Cochrane Handbook
m1=45.2
m2=42.9
n1=17991+4136
n2=13790+4940
meanage=(n1*m1+n2*m2)/(n1+n2)

sd1=3
sd2=1.7

sdage=sqrt( (((n1-1)*sd1^2)+((n2-1)*sd2^2)+((n1*n2)
                                            
                                            /(n1+n2))*(m1^2+m2^2-2*m1*m2))/(n1+n2-1) )



# References:
# 1. 	Efthimiou O. Practical guide to the meta-analysis of rare events. Evid Based Ment Health. 2018;21(2):72–6. 
# 2. 	Bradburn MJ, Deeks JJ, Berlin JA, Localio AR. Much ado about nothing: A comparison of the performance of meta-analytical methods with rare events. Stat Med. 2007;26(1):53–77. 
# 3. 	Wasserstein RL, Lazar NA. The ASA Statement on p -Values: Context, Process, and Purpose. Am Stat [Internet]. 2016 Apr 2;70(2):129–33. Available from: http://www.tandfonline.com/doi/full/10.1080/00031305.2016.1154108
# 4. 	Windish DM, Huot SJ, Green ML. Medicine residents’ understanding of the biostatistics and results in the medical literature. J Am Med Assoc. 2007;298(9):1010–22. 
# 5. 	Kruschke JK. Doing Bayesiabn Data Analysis. 2nd ed. Academic Press; 2014. 776 p. 

