# Functions supporting KKICv

KKICv.lib <-function(){
  library(stringr) #rational string functions
  library(foreach)  # iteration and looping compatible with parallel processing
  library(tictoc) # timing functions
  library(plyr)  #rational wrapper for apply like functions
  library(bbmle) # B.Bolker likelihood tools used for ICtables
  library(lava)  # Latent variable analysis
#  library(rlist)
  library(doFuture) #parallel processing back end
  library(kde1d) # spline smoothed kernel density estimation
}

KKICv.lib()

balsmpl2 <- {  function( #generates indicies for a balanced bootstrap
  B,           # number of bootstraps desired
  sz=NULL,     # number of observations in dataset to be bootstrapped
  strata=NULL, # vector of group idicators
  incl_bt0=T,  # include original data as a bootstrap
  rwname="bt",  # prefix for rownames for bootstraps
  k=1           # number of data-clone replications
) # end of function arguments
{ #function body
  {# Header documentation
    #Generates indicies for a balanced bootstrap (unstratified of stratified)
    #All indicies occur = number of times in matrix  
    #Balanced bootstraps are more efficient than regular bootstraps
    #Use each row of the matrix to select a bootstrapped sample
    #Each row (i.e. bootstrap) is given an identifier in the rowname e.g. "bt324"
    #strata does not need to be sorted
    #One and only one of sz or strata must be given
    #k allows for data-cloning like inflation of data set size
    #Returns a BX(k*sz) matrix of indices for a balanced bootstrap if incl_bt0=F.  
    #Returns a (B)X(k*sz) matrix of indices for a balanced bootstrap if incl_bt0=T.  
    #If incl_b0=T then first row is 1:sz (i.e original estimate)
    #This function replaces balsmpl, and behaves the same if sz is given and k=1.
    #Version 2 adds data cloning facility and changes the behavior when incl_bt0=T
    #  num of output rows now = B regardless of value of incl_bt0 version 1 had B+1 rows
    # Version 1: "2019-05-03 19:44:05 MDT"
    # Version 2: "2019-12-10 10:56:11 MST"
    #code by Mark L. Taper <MarkLTaper@gmail.com>
  } # end Header documentation
  
  { # utility functions
    balsmpli <- function(B,idvec,incl_bt0=T,rwname="bt"){
      #"2019-04-20" by Mark L. Taper
      #Modified balsmpl.  Samples from a vector of indexes rather than assuming 1:sz
      #specifying idvec=1:sz gives same behavior as balsmpl
      #Returns a B X length(idvec) matrix of indices for a balanced bootstrap if bt0=F.  
      #Returns a (B) X length(idvec) matrix of indices for a balanced bootstrap if bt0=T.  
      #All indicies occur = number of times in matrix  
      #Balanced bootstraps are more efficient than regular bootstraps
      #Use each row of the matrix to select a bootstrapped sample
      #Each row (i.e. bootstrap) is given an identifier in the rowname e.g. "bt324"
      #If incl_b0=T then first row is 1:sz if sz is scalar or sz if sz is vector
      if (incl_bt0) {
        if (B==1){
          out <- matrix(idvec,nrow=1)
        } else {
          idxes=sample(rep(idvec,times = B-1))
          out=matrix(idxes,nrow=B-1)
          out <-rbind(1:length(idvec),out)
        }
        rownames(out) <- stringr::str_c(rwname,0:(B-1),sep = "")
      } else {
        idxes=sample(rep(idvec,times = B))
        out=matrix(idxes,nrow=B)
        rownames(out) <- stringr::str_c(rwname,1:B,sep = "")
        
      }
      return(out)
    } #makes indecies for a balanced bootstrap
    
  } # end of utility functions
  
  { # checks and initialization
    if (is.null(sz) & is.null(strata)) stop("one of sz or strata must be given")
    if (!is.null(sz) & !is.null(strata)) stop("only one of sz or strata can be given")
    if (incl_bt0 & (B==1)) warning("(incl_bt0 & (B==1)) This may not be what you want")
    if (!is.null(sz)) {strata <- rep(1,sz)
    } else sz <- length(strata) #data set number of rows
    univl <- unique(strata) #unique values in strata
    num_strt <- length(univl) # number of strata
    idx_strt <- cbind(1:sz,strata)
    colnames(idx_strt) <- c("idx", "strt")
    out <- NULL
  } # end of checks and initialization
  
  for (s in univl){ # build bootstrap indicies by strata
    idvec <- idx_strt[idx_strt[,"strt"]==s,][,"idx"]
    btmat_s <- balsmpli(B,idvec=idvec,incl_bt0 = incl_bt0,rwname = rwname)
    out <- cbind(out,btmat_s)
  } # end of build bootstrap indicies by strata
  
  if (k>1){
    for (i in 1:(k-1)){ # adds data-clone replicate to matrix
      for (s in univl){ # build bootstrap indicies by strata
        idvec <- idx_strt[idx_strt[,"strt"]==s,][,"idx"]
        btmat_s <- balsmpli(B,idvec=idvec,incl_bt0 = incl_bt0,rwname = rwname)
        out <- cbind(out,btmat_s)
      } # end of build bootstrap indicies by strata
    }
    
  }
  return(out)
} 
}# end of function balsmpl2
#function balsmpl2 generates indicies for a balanced bootstrap

BIC.lm.nu <- {function(ftd_lm,nuD=NULL,nunobs=T){
  # modified to single model  "2020-04-25 11:14:05 MDT"
  # revised 2020-05-18 17:08:26 MDT. 
  # requires function logLik.lm.nu
  # default bases the penalty off of dimension of nuD
  #code by Mark L. Taper <MarkLTaper@gmail.com>
  out <- BIC(logLik.lm.nu(ftd_lm = ftd_lm,nuD=nuD,nunobs = nunobs)) 
}} #end of function BIC.lm.nu

DeltaBIC.lm.nu <- {function(erefmd,ealtmd,nuD,nunobs=T) {
  #code by Mark L. Taper <MarkLTaper@gmail.com>
  BIC.ref.nu <- BIC.lm.nu(ftd_lm = erefmd,nuD = nuD,nunobs = nunobs)
  BIC.alt.nu <- BIC.lm.nu(ftd_lm = ealtmd,nuD = nuD,nunobs = nunobs)
  out <- BIC.alt.nu - BIC.ref.nu
  return(out)
}} #end of function DeltaBIC.lm.nu

GKBM <- lvm()#Grace and Keeley 2006 best model
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

GKBM_m_F.A <- GKBM #copy GKBM to GKBM_m_F.a
cancel(GKBM_m_F.A) <- ~ F + A #remove path from m1

kdeBndd <- function( #bounded 1d kernel density estimation using kde1d::kde1d
  #code by Mark L. Taper <MarkLTaper@gmail.com>
  obvec, # vector of observations to be converted to a density
  bndType="TB", # one of "TB","NB","LB","UB" "LUB" for Test, No, Lower,Upper or Upper&Lower bound 
  xmin=NaN, # prespecified lower bound
  xmax=NaN, # prespecified upper bound
  frac=0.01, # fraction of approximate bandwidth bounds are displaced from range of obvec
  svRangeNot_x=F,#if TRUE, reduces size of output object by not saving input data
  ...) # pass through arguments to kde1d::kde1d
{
  {#header documentation
    #wrapper function for kde1d::kde1d
    #Returns a object of class kde1d
    #Estimation can be unbounded i.e. bndType="NB"
    #Estimation can be lower bounded i.e. bndType="LB"
    #Estimation can be upper bounded i.e. bndType="UB"
    #Estimation can be tested for best bounding i.e. bndType="TB" (default)
    #For bounded estimation if a bound is prespecified, it will be used otherwise 
    #   xmin <- min(obvec)-espln(obvec,frac) or xmax <- max(obvec)+espln(obvec,frac)
    #ML Taper "2020-08-20 16:06:36 MDT"
    #ML Taper "2020-09-02 17:18:24 MDT" added bndType "LUB"
  } # end header documentation
  {# functions
    epsln <- function(d,frac=0.001) {
      #Returns a fraction of a generalized bandwidth estimator.
      #Based on eq 3.31 from Silverman 1986. Density estimation for statistical data analysis
      n=length(d)
      IQL <- function(d){
        IQR <- quantile(d,probs = c(.25,.75))
        IQL <- IQR[2]-IQR[1]
        return(IQL)
      }
      A <- min(sd(d),IQL(d)/1.34)
      epsln <- frac*.9*A/(n^.2)
      return(epsln)
    }
    
  }# end functions  
  {
    mnm <- min(obvec)
    mx <- max(obvec)
    if (is.nan(xmin)) {xmin <- mnm-epsln(obvec,frac = frac)}#calculate bound if needed
    if (is.nan(xmax)) {xmax <- mx+epsln(obvec,frac=frac)}#calculate bound if needed
    
    switch (bndType,
            "TB" = {
              kdeNB <- kde1d::kde1d(x = obvec,xmin=NaN,xmax = NaN,...)
              kdeLB <- kde1d::kde1d(x = obvec,xmin = xmin,...)
              kdeUB <- kde1d::kde1d(x = obvec,xmax=xmax,...)
              #            kdeLUB <- kde1d::kde1d(x = obvec,xmin=xmin, xmax=xmax,...)
              #            out <- get(attr(bbmle::AICtab(kdeNB,kdeLB,kdeUB,kdeLUB),which = "row.names")[1])
              out <- get(attr(bbmle::AICtab(kdeNB,kdeLB,kdeUB),which = "row.names")[1])
            },
            "NB" = {out <- kde1d::kde1d(x = obvec,xmin=NaN,xmax = NaN,...)},
            "LB" = {out <- kde1d::kde1d(x = obvec,xmin = xmin,...)},
            "UB" = {out <- kde1d::kde1d(x = obvec,xmax=xmax,...)}#,
            #         "LUB" = {out <- kde1d::kde1d(x = obvec,xmin=xmin,xmax=xmax,...)}
    )
    if (svRangeNot_x){out$x <- range(out$x)} #reduce object size
  }
  return(out)
} # end kdeBndd

KKICv <- function(X,fes,fev,B=2000,M=0,FixedP,cnf_lvls=c(0.95),svRngOnly=T,lBout=F,note="")
{
  #Arguments: FixedP, a object containing additional arguments.  It will be
  #           passed to both fes and fev 
  #  X, data set to be bootstrapped.  either matrix or dataframe observations must be 1 per row.
  # FixedP is a list of fixed parameter for fes probably should contain names of models to be compared
  #  fes, function for EStimating statistic arguments (X, FixedP) fes, returns list with at least 
  #           "l"=statistic and "theta" parameter estimates
  #  fev, function for EValuating statistic. arguments (X,FixedP, theta).
  #  fev must return a list of named objects containing at least "l" the estimated statistic
  #  B, number of bootstraps used
  #  M, number of third level subbootstraps. Used to calculate 3rd order bias correction
  # If svRngOnly=T (default) the bootstrapped vector of values is replace with a range vector in the kde1d object
  # this reduces the object size considerably if B is large.
  # lBout boolean control on whether the raw (potentially large) bootstrap matrix  is included in output
  #End of Arguments:
  #Algorithm from Kitagawa and Konishi 2010. 
  #The goal of K&K 2010 was a reduced variance bias correction for the maximum likelihood.
  #Variable names attempt to mirror that paper's notation so the statistic is named "l".
  #However, THE BIAS CORRECTION CAN BE APPLIED TO ANY STATISTIC THAT IS A FUNCTIONAL OF AN EMPIRICAL CDF.
  #Throughout the code (and K&K 2010) you will see variable names of the form "l" followed by 2 integers.
  #i.e. c("l00","l01","l10","l11","l12","l21","l22"). The left integer indicates the bootstrap level of
  #the data used to evaluate the statistic and the right integer indicates the bootstrap level of the 
  #parameters estimates used to evaluate the statistic. For example l00 is the statistic evaluated
  #with the original data and parameters estimated from the original data.  l11 is the statistic evaluated
  #using bootstrapped data and parameter estimates from bootstrapped data. l01 is the statistic using original
  #data and bootstrapped parameter. l10 is the statistic using bootstrapped data and original parameters.
  #Value: list(est.bc, cnf_tables, run_info, bndKD.objs)
  #est.bc: list(theta00 (estimated parameter vector), l00 (estimated statistic), 
  #    b1bse (standard error of the 1st level bias correction), b1b (1st level bias correction),
  #    b2b (2nd level bias correction), l1bc (1st level bias corrected statistic ),
  #    l2bc (2nd level bias corrected statistic)
  #     )
  #cnf_tables: list(
    # uncon.bnds (unconditional bounds at the levels specified by cnf_lvls)
    # con.l01.bnds (bounds conditional on original data)
    # con.l10.bnds (bounds conditional on original parameter estimates)
  #) confidence tables contain both quantiles of the bootstrap distribution and
  # bounds that are from kde1d smoothing of the bootstrap distribution.
  # run_info list(all inital arguments and run-time information)
  # bndKD.objs: list(uncon.kde.obj, con.l01.kde.obj, con.l10.kde.obj)
  # In bndKD.objs, the kde.objs are the kde1d smoothed representations of the 
  # bootstrap distributions.
  # lB: if(lBout==T) the (B x 8) matrix of bootstrapped statistics
  #Coding Mark L. Taper: <MarkLTaper@gmail.com>
  #Warning: The KKIC's  variance correction doesn't kick in until moderate sample size ~100
  #2018-12-28 08:18:26 MST corrected apparent error in b2b see eq 54 K&K 2010
  #"2020-06-17 11:12:51 MDT": Version 2 calibrated intervals eliminated
  #2020-06-18 20:07:30 MDT: Foreach parallel processing implemented
  #2020-10-12 16:41:42 MDT, Version 4: conditional bounds based on l(x;theta*) renamed con.l01
  #new conditional kde object based on l(x*;theta^) added.  This is named con.l10.
  #"2021-02-24 17:06:45 MST": Version 5: code & documentation cleanup
  version <- 5
  {run <- vector(mode = "list", length = 13)
    names(run) <- c("note","Version","start","stop","elapsed",
                    "B","M","smpl_sz","data0","FixedP","fes","fev","cnf_lvls")
    run$note <- note
    run$start <- Sys.time()
    run$Version <- version
    run$B <- B
    run$M <- M
    run$smpl_sz <- dim(X[1])
    run$data0 <- X
    run$FixedP <- FixedP
    run$fes <- fes
    run$fev <- fev
    run$cnf_lvls <- cnf_lvls
  }
  
  cnf_pnts <-c((1-cnf_lvls)/2,1-rev((1-cnf_lvls)/2))
  bnd_type=c(rep("L",length(cnf_lvls)),rep("U",length(cnf_lvls)))
  if(is.vector(X)) {X=as.matrix(X)}
  Xlen = dim(X)[1]
  est0=fes(X,FixedP)
  theta0=est0$theta
  l00=est0$l
  lB=matrix(NA,ncol=7,nrow=B,dimnames = list(NULL,c("l00","l01","l10","l11","l12","l21","l22")))
  l22.mat <- l21.mat <- l12.mat <- l02.mat <- matrix(NA,nrow = M,ncol = B) #initializing second level bootstrap matrices
  Bidx=balsmpl2(B=B,sz=Xlen,incl_bt0 = F) #incl_bt0=F because theta0 already estimated
  wrkrs <- getDoParWorkers() 
  blksz <- B/wrkrs
  lB <- foreach (i = 1:B,.combine = rbind,.inorder = F,.options.future = list(chunk_size = blksz))%dopar%{
    Xi=as.data.frame(X[Bidx[i,],])
    esti=fes(Xi,FixedP)
    thetai=esti$theta
    l11=esti$l
    l10=fev(Xi,FixedP,theta0)$l
    l01=fev(X,FixedP,thetai)$l
    if (M == 0) {l12=NA; l21 = NA; l22 = NA}
    else {
      if (M == 1) Midx=matrix(sample(1:Xlen,size = Xlen,replace = T),nrow=1)
      else Midx = balsmpl2(M,sz=Xlen)
      l12 <- l21 <- l22 <- 0
      for (j in 1:M){
        Xib=as.data.frame(Xi[Midx[j, ], ])
        estib=fes(Xib,FixedP)
        thetaib=estib$theta
        l22.mat[j,i]=fev(Xib,FixedP,thetaib)$l 
        l21.mat[j,i]=fev(Xib,FixedP,thetai)$l 
        l12.mat[j,i]=fev(Xi,FixedP,thetaib)$l
        l02.mat[j,i] <- fev(X,FixedP,thetaib)$l
      }
      l22 = mean(l22.mat[,i])
      l21 = mean(l22.mat[,i])
      l12 = mean(l22.mat[,i])
    }
    rowi=c(l00,l01,l10,l11,l12,l21,l22)
    return(rowi)
  }
  colnames(lB) <- c("l00","l01","l10","l11","l12","l21","l22")
  lB <- as.data.frame(lB)
  b1bse <- with(lB,sd((l11-l10+l00-l01),na.rm = T)/sqrt(B))
  b1b <- with(lB,mean(l11-l10+l00-l01))
  b2b <- b1b-with(lB,mean(l22-l21+l11-l12,na.rm=T))
  l00 <- est0$l
  n <- dim(X)[1]
  lb1c <- l00 - b1b # non parametric bias corrections
  if (M>1) # 3rd order corrections.  Warning high variance unless both ss & B are large
  {
    lb2c <- l00 -b2b
  } else {
    lb2c <- NULL
  }
  l00.vec <- lB[,"l00"] #vector repeating original estimate
  l11.vec <- lB[,"l11"] #vector of unconditional bootstrap estimates.  
  l01.vec <- lB[,"l01"] #vector of conditional bootstrap estimates. fixed x bootstrap theta
  l10.vec <- lB[,"l10"] #vector of conditional bootstrap estimates. fixed theta, bootstrap x
  lb1c.vec <- l00.vec - (l11.vec - l10.vec +l00.vec - l01.vec)
  
  uncon.bnds <- t(as.matrix(quantile(x = l11.vec,probs = cnf_pnts)))
  rownames(uncon.bnds) <-c("qnt.bnds")
  uncon.kde.obj <- kdeBndd(obvec = l11.vec,bndType = "TB")
  if (svRngOnly) {uncon.kde.obj$x <- range(uncon.kde.obj$x)}
  uncon.bnds <- rbind(uncon.bnds,kde.bnds=qkde1d(p = cnf_pnts,obj = uncon.kde.obj))
  
  con.l01.bnds <- t(as.matrix(quantile(x = l01.vec,probs = cnf_pnts)))
  rownames(con.l01.bnds) <-c("qnt.bnds")
  con.l01.kde.obj <- kdeBndd(obvec = l01.vec,bndType = "TB")
  if (svRngOnly) {con.l01.kde.obj$x <- range(con.l01.kde.obj$x)}
  con.l01.bnds <- rbind(con.l01.bnds,kde.bnds=qkde1d(p = cnf_pnts,obj = con.l01.kde.obj))
  
  con.l10.bnds <- t(as.matrix(quantile(x = l10.vec,probs = cnf_pnts)))
  rownames(con.l10.bnds) <-c("qnt.bnds")
  con.l10.kde.obj <- kdeBndd(obvec = l10.vec,bndType = "TB")
  if (svRngOnly) {con.l10.kde.obj$x <- range(con.l10.kde.obj$x)}
  con.l10.bnds <- rbind(con.l10.bnds,kde.bnds=qkde1d(p = cnf_pnts,obj = con.l10.kde.obj))
  
  # lb1c.bnds <- t(as.matrix(quantile(x = lb1c.vec,probs = cnf_pnts)))
  # rownames(lb1c.bnds) <-c("qnt.bnds")
  # lb1c.kde.obj <- kdeBndd(obvec = lb1c.vec,bndType = "TB")
  # if (svRngOnly) {lb1c.kde.obj$x <- range(lb1c.kde.obj$x)}
  # lb1c.bnds <- rbind(lb1c.bnds,kde.bnds=qkde1d(p = cnf_pnts,obj = con.l10.kde.obj))
  
  est.bc <- list( theta00=est0$theta,l00=l00,b1bse=b1bse,b1b=b1b,b2b=b2b,lb1c=lb1c,lb2c=lb2c)
  cnf_tables <- list(uncon.bnds=uncon.bnds,con.l01.bnds=con.l01.bnds,con.l10.bnds=con.l10.bnds)
  run$stop <- Sys.time()
  run$elapsed <- run$stop - run$start
  bndKD.objs <- list(uncon.kde.obj=uncon.kde.obj,con.l01.kde.obj=con.l01.kde.obj,
                     con.l10.kde.obj=con.l10.kde.obj)
  if (lBout == F) {lB <- NA} # suppresses saving of raw bootstrap matrix
  out <- list(est.bc=est.bc,cnf_tables=cnf_tables,run_info=run,bndKD.objs=bndKD.objs,lB=lB)
  return(out)
} #end of function KKICv version 5

fes.DeltaBIC.lm <- {# fes function for linear models
  #code by Mark L. Taper <MarkLTaper@gmail.com>
  function(X,FixedP){
    refmd <- FixedP$refmd
    altmd <- FixedP$altmd
    erefmd <- lm(data=as.data.frame(X),as.formula(refmd))
    ealtmd <- lm(data=as.data.frame(X),as.formula(altmd))
    l <- DeltaBIC.lm.nu(erefmd = erefmd,ealtmd = ealtmd,nuD = NULL,nunobs = F)
    theta <- list(erefmd=erefmd,ealtmd=ealtmd)
    out <- list(l=l,theta=theta)
    return(out)
  }} #end of function fes.DeltaBIC.lm

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

fes.sem.TIC <- function(X,FixedP){
  # estimates a sem models given data and returns -2*logLik and fitted models
  #code by Mark L. Taper <MarkLTaper@gmail.com>
  md <- FixedP$md
  emd <- lava::estimate(x = md,data = X) 
  l <- -2*logLik(emd)
  theta <- list(emd=emd)
  out <- list(l=l,theta=theta)
  return(out)
}

fev.DeltaBIC.lm <- {# fev for a linear models
  #code by Mark L. Taper <MarkLTaper@gmail.com>
  function(X,FixedP,theta){
    erefmd <- theta$erefmd
    ealtmd <- theta$ealtmd 
    out=list(l=DeltaBIC.lm.nu(erefmd =erefmd,ealtmd = ealtmd,nuD = X,nunobs = T ))
    return(out)
  }} #end of function fev.DeltaBIC.lm

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

fev.sem.TIC <- function(X,FixedP,theta){
  #evaluates 2*logLik of a fitted sem models with new data
  #code by Mark L. Taper <MarkLTaper@gmail.com>
  emd <- theta$emd
  l <- -2*logLik(object = emd,data=X)
  #Note the logLik method for objects of class lvmfit supports new data
  out=list(l=l,nuD = X,nunobs = T )
  return(out)
}

#FixedP.DeltaBIC.lm <- list(refmd=,altmd=) # example FixedP.DeltaBIC.lm

FixedP.sem.DelBIC = list(refmd=GKBM,altmd=GKBM_m_F.A)#example FixedP

FixedP.sem.TIC=list(md=GKBM) #example FixedP


logLik.lm.nu <- {function(ftd_lm,nuD=NULL,nunobs=T){
  #logLik of new data for a fitted lm
  #currently only for a single dependent variable
  #revised 2020-05-18 17:08:26 MDT
  #default now bases nobs off the dim(nuD)
  #code by Mark L. Taper <MarkLTaper@gmail.com>
  lLf <- logLik(ftd_lm)
  if (is.null(nuD)){
    out <- lLf
    attributes(out)$nuD <- FALSE
  } else {
    dv.nm <- colnames(model.frame(ftd_lm)[1]) #extracting the dependent variable name
    dvs.mle=sd(ftd_lm$residuals)*sqrt((nobs(ftd_lm)-1)/nobs(ftd_lm))
    out=sum(dnorm(x = (nuD[,dv.nm]-predict.lm(ftd_lm,newdata = nuD)),mean = 0,sd = dvs.mle,log = T))
    attributes(out) <- attributes(lLf)
    if (nunobs) {attributes(out)$nobs <- dim(nuD)[1]}
    attributes(out)$nuD <- TRUE
  }
  return(out)
}} # end of function logLik.lm.nu

mvskkde <- function(kde.obj,n=100000,reps=10) {
  #returns a vector of estimates of the mean, variance, skew, and extra kurtosis
  # of a kde1d object
  #moments are calculated by integration.
  #if there is an integration error, moments are calculated by simulation. 
  #Only informative digits are shown in simulated results 
  #n and reps influence simulation if integration has failed
  # total number of observation in simulation = n*reps*(number of execution workers)
  # Mark L. Taper "2021-02-24 14:01:32 MST"
  #code by Mark L. Taper <MarkLTaper@gmail.com>
  if ( !(class(kde.obj)=="kde1d")) stop("object class not kde1d")
  mvskkdeSB <- function(kde.obj,n=n,reps=reps){
    # mean, variance, skew and kurtosis of a kde1d object by simulation
    # requires packages foreach and dofuture
    # total sample size involved in estimation is n*reps*wrkrs
    # only approximately significant digits retained
    se <- function(x){sd(x)/sqrt(length(x))}
    wrkrs <- getDoParWorkers()
    statmat <- foreach(i=1:(reps*wrkrs),.combine = rbind,.inorder = F)%dopar%{
      vec <- rkde1d(n = n,obj = kde.obj)
      stats <- c(mean(vec),var(vec),SkewKurtosis_se(vec)[1,1],SkewKurtosis_se(vec)[2,1])   
      return(stats)
    }
    statmat=matrix(statmat,ncol=4) #makes sure satmat is a matrix even if reps=1
    out <- colMeans(statmat) #approximate resolution with 100000 points
    ses <- plyr::aaply(.data = statmat,.margins = 2,.fun = se)
    for (i in 1:4) {out[i] <- round(out[i],digits=round(log10(1/ses[i])))}
    names(out) <- c("mean","var","skew","xkurtosis")
    return(out)
    
  } # end of mvskkdeSB
  
  wd <- 2*(max(kde.obj$x)- min(kde.obj$x))
  lo <- min(kde.obj$x)-wd #lower integration limit
  hi <- max(kde.obj$x)+wd #upper integration limit
  xd <- function(x,kde.obj) {x*dkde1d(x,kde.obj)}
  out <- rep(NA,4)
  names(out) <- c("mean","var","skew","xkurtosis")
  
  mn <- tryCatch(error = function(cnd){"is error"},integrate(f = xd,lower = lo,upper = hi,subdivisions = 1000,kde.obj=kde.obj)$value)
  if(mn=="is error") {out <- mvskkdeSB(kde.obj)
  } else {
    out["mean"] <- mn
    devpd <- function(x,mn,p,kde.obj){((x-mn)^p)*dkde1d(x,kde.obj)  }
  }
  if (is.na(out["var"])) {
    mn <- out["mean"]
    vr <- tryCatch(error = function(cnd){"is error"},integrate(f = devpd,lower = lo,upper = hi,subdivisions = 1000,mn=mn,p=2,kde.obj=kde.obj)$value)
    if (vr == "is error") {out <- mvskkdeSB(kde.obj)
    } else {out["var"] <- vr}
  }
  if (is.na(out["skew"])) {
    mn <- out["mean"]
    vr <- out["var"]
    sd <- sqrt(vr)
    skw <- tryCatch(error = function(cnd){"is error"}, integrate(f = devpd,lower = lo,upper = hi,subdivisions = 1000,mn=mn,p=3,kde.obj=kde.obj)$value/(sd^3))
    if (skw == "is error") {out <- mvskkdeSB(kde.obj)
    } else {out["skew"] <- skw}
  }
  if (is.na(out["xkurtosis"])) {
    mn <- out["mean"]
    vr <- out["var"]
    sd <- sqrt(vr)
    exKrt <- tryCatch(error = function(cnd){"is error"}, (integrate(f = devpd,lower = lo,upper = hi,subdivisions = 1000,mn=mn,p=4,kde.obj=kde.obj)$value/(sd^4))-3)
    if (exKrt == "is error") {out <- mvskkdeSB(kde.obj)
    } else {out["xkurtosis"] <- exKrt}
  }
  return(out)
} #end of mvskkde

points.kde <- function(kde.obj,...){
  # adds a plot of kde1d density to an existing plot
  #code by Mark L. Taper <MarkLTaper@gmail.com>
  if ( !(class(kde.obj)=="kde1d")) stop("object class not kde1d")
  x <- kde.obj$grid_points
  y <- kde.obj$values
  xy <- list(x=x,y=y)
  points(xy,...)
}

shazam <- function(pln="multi"){#partial matching of sequential or multiprocess 
  # registers a doFuture parallel processing back sets processing to either
  # multiprocess or sequential
  #code by Mark L. Taper <MarkLTaper@gmail.com>
  registerDoFuture()
  switch(  pmatch(pln,table = c("sequential","multiprocess")),
           sequential=plan(sequential),
           multiprocess=plan(multiprocess)
  )
  cat("DoFuture parallel processing backend with",getDoParWorkers(),"workers\n")
}

SkewKurtosis_se=function(x) # skewnes & Kurtosis with s.e.
{
  {
    # Skewness and kurtosis and their standard errors as implement by SPSS
    #
    # Reference: pp 451-452 of
    # http://support.spss.com/ProductsExt/SPSS/Documentation/Manuals/16.0/SPSS 16.0 Algorithms.pdf
    # ported to R by Howard Seltman http://www.stat.cmu.edu/~hseltman/files/spssSkewKurtosis.R
    # See also: Suggestion for Using Powerful and Informative Tests of Normality,
    # Ralph B. D'Agostino, Albert Belanger, Ralph B. D'Agostino, Jr.,
    # The American Statistician, Vol. 44, No. 4 (Nov., 1990), pp. 316-321
    # current wrapper by Mark L. Taper <MarkLTaper@gmail.com>
  } # end of header documentation
  w=length(x)
  m1=mean(x)
  m2=sum((x-m1)^2)
  m3=sum((x-m1)^3)
  m4=sum((x-m1)^4)
  s1=sd(x)
  skew=w*m3/(w-1)/(w-2)/s1^3
  sdskew=sqrt( 6*w*(w-1) / ((w-2)*(w+1)*(w+3)) )
  kurtosis=(w*(w+1)*m4 - 3*m2^2*(w-1)) / ((w-1)*(w-2)*(w-3)*s1^4)
  sdkurtosis=sqrt( 4*(w^2-1) * sdskew^2 / ((w-3)*(w+5)) )
  mat=matrix(c(skew,kurtosis, sdskew,sdkurtosis), 2,
             dimnames=list(c("skew","kurtosis"), c("estimate","se")))
  return(mat)
}

STIC <- function(lvm.obj){#a consistent TIC for sem
  TIC.lvmfit(lvm.obj,fnob=log)}

TIC <- function(lvm.obj){TIC.lvmfit(lvm.obj)}

TIC.lvmfit <- function #TIC based information criterion for an lvmfit class object
( lvm.obj, # lvmfit object
  fnob=NULL # consistency multiplier function of nobs 
) # end arguments for TIC.lvmfit
{ # TIC.lvmfit function body
  {# TIC.lvmfit header documentation
    # Calculates a TIC based information criterion
    # TIC = -2*lL + fnob(nobs)*tr(I%*%J^-1)
    # by default fnob(nobs)== 2 and the information criterion of Takeuchi is returned
    # Mark L. Taper "Mon Nov 18 17:00:38 2019" <MarkLTaper@gmail.com>
  }
  nobs.lvmfit <- function(ob){return(ob$data$n)}
  if (is.null(fnob)) { #sets consitency multiplier to default
    fnob <- function(n){2}
  }
  out <- -2*logLik(lvm.obj)+fnob(nobs(lvm.obj))*tr(information(lvm.obj,type="outer")%*%solve(information(lvm.obj,type="hessian")))
  out <- as.numeric(out)
  return(out)
} #end TIC.lvmfit function






