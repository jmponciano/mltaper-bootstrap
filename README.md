mltaper-bootstrap:  repository for the programs in the manuscript

Taper, M.L., S.R. Lele, J.M. Ponciano, B. Dennis and C.L. Jerde in 
"Assessing the global and local uncertainty in scientific evidence in the 
presence of model misspecification" (submitted to Frontiers in Ecology and evolution) 
use the R function KKICv to construct confidence intervals for evidential comparisons. 
The file "KKICv.demo.R" demonstrate the uses of the KKICv function for several problems.

The KKICv function uses the reduced variance bootstrap bias correction to general 
functional estimates developed by S. Konishi and G. Kitagawa in this long series  
of articles:

Kitagawa, G., and S. Konishi. 2010. 
Bias and variance reduction techniques for bootstrap information criteria. 
Annals of the Institute of Statistical Mathematics 62:209-234.

Konishi, S., and G. Kitagawa. 1996. 
Generalised information criteria in model selection. Biometrika 83:875-890.

Konishi, S., and G. Kitagawa. 2003. Asymptotic theory for information criteria 
in model selection - functional approach. 
Journal of Statistical Planning and Inference 114:45-61.

Konishi, S., and G. Kitagawa. 2008. 
Information Criteria and Statistical Modeling. Springer, New York.

The KKICv function requires packages:  stringr, foreach, tictoc, plyr, bbmle, lava
rlist, doFuture, and kde1d. These packages should be installed before proceeding.

Please contact Mark L. Taper for questions regarding the implementation of these functions
e-mail: markltaper@gmail.com.  Thanks!   


