#KKICv demo

# Taper, M.L., S.R. Lele, J.M. Ponciano, B. Dennis and C.L. Jerde in 
# "Assessing the global and local uncertainty in scientific evidence in the 
# presence of model misspecification" (submitted) use the R function KKICv to 
# construct confidence intervals for evidential comparisons. This file 
# demonstrate the uses of the KKICv function for several problems.

# The KKICv use the reduced variance bootstrap bias correction to general 
# functional estimates developed by S. Konishi and G. Kitagawa in a long series  
# of articles. :

# Kitagawa, G., and S. Konishi. 2010. 
# Bias and variance reduction techniques for bootstrap information criteria. 
# Annals of the Institute of Statistical Mathematics 62:209-234.
#
# Konishi, S., and G. Kitagawa. 1996. 
# Generalised information criteria in model selection. Biometrika 83:875-890.

# Konishi, S., and G. Kitagawa. 2003. Asymptotic theory for information criteria 
# in model selection - functional approach. 
# Journal of Statistical Planning and Inference 114:45-61.

# Konishi, S., and G. Kitagawa. 2008. 
# Information Criteria and Statistical Modeling. Springer, New York.

# The KKICv requires packages:  stringr, foreach, tictoc, plyr, bbmle, lava
# rlist, doFuture, and kde1d. These packages should be installed before proceeding.


source(file = "KKICv.R") #sources KKICv and support functions
KKICv.lib() # loads required libraries

shazam("multi") #registers a DoFuture parallel processing back end

#readline(prompt = "Pause. Press <Enter> to continue...")

GK.d <- readRDS(file = "GK.d.rds") # example data from Grace and Keeley 2006

#To use KKICv one needs to construct two wrapper functions for the statistic to 
#bias corrected. A fes function that estimates the necessary parameters given a 
# data set and then evaluates the desired statistic and a fev function that 
# evaluates the statistic given a data set and set of parameters.

#In Taper et al.  The principal statistic of interest is the difference of 
#Schwarz information criterion (SIC, a.k.a. BIC) values for two different estimated
#models. We will compare as the reference model:

GKBM <- lvm()#Grace and Keeley 2006 best model.  This is a package lava SEM
{
  regression(GKBM) <- list(
    distance ~ L, # observed ~ latent
    abiotic ~ C,
    age ~ A,
    hetero ~ H,
    firesev ~ F,
    cover ~ P,
    rich ~ R,
    A ~ L,
    C ~ L,
    H ~ L,
    F ~ A,
    P ~ F,
    R ~ C + H + L + P)
  latent(GKBM) <- ~ L+C+A+H+F+P+R 
  regression(GKBM,distance ~ L) <- 1
  regression(GKBM,abiotic ~ C) <- 1
  regression(GKBM,age ~ A) <- 1
  regression(GKBM,hetero ~ H) <- 1
  regression(GKBM,firesev ~ F) <- 1
  regression(GKBM,cover ~ P) <- 1
  regression(GKBM,rich ~ R) <- 1
  covariance(GKBM, ~ distance + abiotic + age + cover + rich) <- 0 # fixed indicator coef
  covariance(GKBM, ~ hetero + firesev) <- .1 # Known measurement error variances
  covariance(GKBM) <- P ~ H # error correlations between latents as per G&K
}

#to an alternative model
GKBM_m_F.A <- GKBM #copy GKBM to GKBM_m_F.a
cancel(GKBM_m_F.A) <- ~ F + A #remove path from m1

#First we need to construct wrappers fes and fev for statistic estimation and 
#evaluation.

fes.sem.DelBIC <- function(X,FixedP){
  # estimates 2 sem models given data using the lava package
  # and returns Delta BIC and fitted models.
  # FixedP is list of the models with names "refmd" and "altmd"
  # e.g. FP=list(refmd=m0,altmd=m1) #example FixedP
  # return value list(l=l, theta=theta).
  # l is your statistic of interest and theta is an estimated parameter vector.
  # YOU MUST USE THE NAMES "l" AND "theta"! 
  # code by Mark L. Taper <MarkLTaper@gmail.com>
  refmd <- FixedP$refmd
  altmd <- FixedP$altmd
  erefmd <- lava::estimate(x = refmd,data = X) 
  ealtmd <- lava::estimate(x = altmd,data = X)
  l <- BIC(ealtmd)-BIC(erefmd) 
  theta <- list(erefmd=erefmd,ealtmd=ealtmd)
  out <- list(l=l,theta=theta)
  return(out)
}

fev.sem.DelBIC <- function(X,FixedP,theta){
  #evaluates Delta(BIC) of 2 fitted sem models with new data
  #FixedP is of the form FP=list(refmd=m0,altmd=m1)
  #theta is the theta value from the output of corresponding fes
  #code by Mark L. Taper <MarkLTaper@gmail.com>
  erefmd <- theta$erefmd
  ealtmd <- theta$ealtmd 
  l <- BIC(logLik(object = ealtmd,data=X))-BIC(logLik(object = erefmd,data=X))
  #Note the logLik method for objects of class lvmfit supports new data
  out=list(l=l,nuD = X,nunobs = T )
  return(out)
}

#FixedP is a list of whatever information besides the data needed to make 
# the estimate.  The structure of FixedP must be matched to fes

FixedP.sem.DelBIC = list(refmd=GKBM,altmd=GKBM_m_F.A)#example FixedP

#if you haven't registered the DoFuture back end using our function shazam,
#do so now!

shazam("multi")

#Now test how long per bootstrap the KKICv takes on your machine.
#tic()and toc()are convenient stopwatch start and stop functions from the 
#package tictoc. Start with B a smaller number
tic()
tst <- KKICv(X = GK.d,fes = fes.sem.DelBIC,fev = fev.sem.DelBIC,B = 128,
             M=0,FixedP = FixedP.sem.DelBIC,cnf_lvls = c(.95,.9,.5,0),svRngOnly = T,lBout = F)
toc()

# 16.53 sec elapsed

#now time it again.  
tic()
tst <- KKICv(X = GK.d,fes = fes.sem.DelBIC,fev = fev.sem.DelBIC,B = 128,
             M=0,FixedP = FixedP.sem.DelBIC,cnf_lvls = c(.95,.9,.5,0),svRngOnly = T,lBout = F)
toc()
# 4.82 sec elapsed
# The primary difference in time is due to initial checks the DoFuture backend 
# runs the first time a function requiring parallel processing is run.
# We now know that on this machine, if we set B=2000, computation should take
# a little over a minute.

tic()
tst <- KKICv(X = GK.d,fes = fes.sem.DelBIC,fev = fev.sem.DelBIC,B = 2000,
             M=0,FixedP = FixedP.sem.DelBIC,cnf_lvls = c(.95,.9,.5,0),svRngOnly = T,lBout = F)
toc()

#75.11 sec elapsed

#Let's extract some useful information from our object tst.
#Let's look at the mean and confidence bounds for the global & local evidence.

mvskkde(tst$bndKD.objs$uncon.kde.obj)
tst$cnf_tables$uncon.bnds

mvskkde(tst$bndKD.objs$con.l01.kde.obj)
tst$cnf_tables$con.l01.bnds

#We can plot the two density estimate
Ylim=c(0,max(tst$bndKD.objs$con.l01.kde.obj$values))
plot(tst$bndKD.objs$uncon.kde.obj,ylim=Ylim)
points.kde(tst$bndKD.objs$con.l01.kde.obj,type="l",lty=2)
