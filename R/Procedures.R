### Main procedures

#' Global SSR minimizers procedure
#'
#' A helper function to identify if the estimated break model is i) pure change or ii)
#' partial change model. The procedure then calls appropriate functions \link{} to estimate
#' the pure change model and \link{} to estimate the partial change model. This helper function
#' is required for supF, UDMax, WDMax and supF(l+1|l) test functions invoked via \link{}
#'
#'@aliases doglob
#'@param y_name dependent variables in matrix form
#'@param z_name matrix of independent variables with coefficients are allowed to change across
#'regimes
#'@param x_name  matrix of independent variables with coefficients constant across regimes
#'@param eps1 trimming level
#'@param eps convergence criterion for iterative recursive computation
#'@param maxi maximum number of iterations
#'@param fixb option to use fixed initial input \eqn{\beta}. If \code{1},
#'the model will use values given in \code{betaini}. If \code{0}, betaini is skipped
#'@param betaini Initial \eqn{beta_0} to use in estimation
#'@param printd option to print results of iterations for partial change model
#'@return A list containing the following components:
#'\itemize{
#'\item{glb} {Minimum global SSR}
#'\item{datevec} {Vector of dates (optimal minimizers)}
#'\item{bigvec} {Associated SSRs with possible break dates combination}}
#'@export
#'
doglob = function (y,z,x,m,eps,eps1,maxi,fixb,betaini,printd){

#check if model is pure or partial change
if (is.null(x)) {p = 0}
else {p = dim(x)[2]}

q = dim(z)[2]
bigT = dim(y)[1]
h = round(eps1*bigT)


if(p == 0) {
  if (printd == 1){
  cat('This is a pure structural change model with the following specifications:\n')
  cat(paste(q,'regressors z with allowed to change coefficients\n'))
  cat(paste('maximum number of breaks:',m),'\n') }
  out = dating(y,z,h,m,q,bigT)
  }
else{
  if (printd == 1){
  cat('This is a partial structural change model with the following specifications:\n')
  cat(paste(q,'regressors z with allowed to change coefficients\n'))
  cat(paste(p,'regressors x with fixed coefficients\n'))
  cat(paste('maximum number of breaks:',m),'\n')
  if(fixb == 1) {cat('initial regime-wise coefficients: \n')
    cat(betaini)}
  cat(paste('convergence criterion:',eps),'\n')
  cat(paste('print iteration option (1)=TRUE/(0)=FALSE',printd),'\n')}
  out = nldat(y,z,x,h,m,p,q,bigT,fixb,eps,maxi,betaini,printd)
}
glb = out$glb
datevec = out$datevec
bigvec = out$bigvec
#printing results
if (printd == 1) {
for (i in 1:m){
  cat(paste('Model with',i,'breaks has SSR:',glb[i,1]),'\n')
  cat('The dates of breaks are:\n')
  cat(datevec[1:i,i])
  }
}
  return (out)
}

#' SupF, UDMax & WDMax testing procedure
#'
#' @description
#' The procedure calculate the test statistics and print results of the 2 main tests:
#' \itemize{
#' \item{SupF test} {F test of 0 vs m breaks}
#' \item{Double Max test} {UDMax: the unweighted version
#' and WDMax: the weighted version}
#'}
#'
#'@param y_name name of dependent variable in the data set
#'@param z_name name of independent variables in the data set which coefficients are allowed to change
#'across regimes. \code{default} is vector of 1 (Mean-shift model)
#'@param x_name name of independent variables in the data set which coefficients are constant across
#' regimes. \code{default} is NULL
#'@param data the data set for estimation
#'@param m maximum number of breaks
#'@param eps1 trimming level
#'@param eps convergence criterion for iterative recursive computation
#'@param maxi maximum number of iterations
#'@param fixb option to use fixed initial input \eqn{\beta}. If \code{1},
#'the model will use values given in \code{betaini}. If \code{0}, betaini is skipped
#'@param betaini Initial \eqn{beta_0} to use in estimation
#'@param printd option to print results of iterations for partial change model
#'@param prewhit option to use AR(1) for prewhitening
#'@param robust,hetdat,hetvar options on error terms assumptions
#'@return A list that contains following:
#'\itemize{
#'\item{ftest: SupF test of 0 vs m (1 to maximum) breaks statistics}
#'\item{cv_supF: Critical values for Sup F test }
#'\item{cv_Dmax: Critical values for Double Max test}
#'\item{supF1: table summarizing the SupF test (for viewing purposes)}
#'\item{UDMax: table summarizing the Double Max test (including UDMax statistics and CVs)}
#'}
#'@export


dotest = function(y_name,z_name=NULL,x_name=NULL,data,
                  m=5,eps=0.00001,eps1=0.15,maxi=10,fixb=0,betaini=0,printd=0,prewhit=1,robust=1,
                  hetdat=1,hetvar=1,hetq=1,hetomega=1,const=1){
  eps1=check_trimming(eps1)
  siglev=matrix(c(10,5,2.5,1),4,1)
  df = process_data(y_name = y_name,z_name = z_name,x_name = x_name,data=data,const)
  y = df$y
  z = df$z
  x = df$x


  if(m==0){
    warning('The test is undefined for no breaks model')
    out = list()
    out$mbreak = 0;
    class(out) = 'sbtests'
    return(out)
  }else{

    if (is.null(x)) {p = 0}
    else {p = dim(x)[2]}

    q = dim(z)[2]
    bigT = dim(y)[1]
    h = round(eps1*bigT)

    m = check_m(bigT,eps1,m)


  out = doglob(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd)
  datevec = out$datevec

  #procedure for F test
  #print('a) supF tests against a fixed number of breaks')

  ftest = matrix(0L, nrow = m, ncol = 1)
  wftest = matrix(0L, nrow = m, ncol = 1)

  for (i in 1:m){
    ftest[i,1] = pftest(y,z,i,q,bigT,datevec,prewhit,robust,x,p,hetdat,hetvar)
    if(printd==1){
    print(paste('supF test for 0 versus',i,'breaks (scaled by q):',ftest[i,1]))}
  }

  cv_supF = matrix(0L,4,m)
  for (sl in 1:4){
    #critical values for supF test
    cv = getcv1(sl,eps1)
    cv_supF[sl,] = cv[q,1:m,drop=FALSE]
    if (printd==1){
    print(paste('The critical values at the',siglev[sl,1],'% level are (for k = 1 to',m,'):'))
    print(cv[q,1:m,drop=FALSE])}
  }

  #procedure for Dmax and UDmax test

  #print('b) Dmax test against an unknown number of breaks')
  #print(paste('The UDmax test is:',max(ftest)))
  cv_Dmax = matrix(0L,4,1)
  for (sl in 1:4) {
    #critical values for Dmax test
    cvm = getdmax(sl,eps1)
    cv_Dmax[sl,1] = cvm[q,1]
    if(printd==1){
    print(paste('The critical values at the',siglev[sl,1],'% level is:',
                cvm[q,1]))}
  }


  for (sl in 1:4){
    #computation of WDmax test
    cv = getcv1(sl,eps1)
    for( i in 1:m){
      wftest[i,1] = cv[q,1] * ftest[i,1] / cv[q,i]
    }
    if (printd==1){
    print(paste('WDmax test at the',siglev[sl,1],'% level is:',max(wftest)))}
  }
  rownames(cv_supF) = siglev
  rownames(cv_Dmax) = siglev


  out = list('ftest' = ftest, 'cv_supF' = cv_supF,
             'cv_Dmax' = cv_Dmax, 'UDmax' = max(ftest))
  out$mbreak = m
  class(out) = 'sbtests'

  out = compile_sbtests(out)
  return(out)}
}


#' SupF(l+1|l) test
#'
#' Function computes the procedure of SupF(l+1|l) test. The function returns
#' the test statistics of supF(l+1|l) test
#' with null hypothesis is maximum number of break is l
#' and alternative hypothesis is l+1.
#' The l breaks under the null hypothesis are taken from the
#' global minimization. Also, new date (if available) and critical values based on
#' significant levels are returned for plotting and inference
#'
#'@param y dependent variables in matrix form
#'@param z matrix of independent variables with coefficients are allowed to change across
#'regimes
#'@param x matrix of independent variables with coefficients constant across regimes
#'@param m maximum number of breaks
#'@param eps1 trimming level
#'@param eps convergence criterion for iterative recursive computation
#'@param maxi maximum number of iterations
#'@param fixb option to use fixed initial input \eqn{\beta}. If \code{1},
#'the model will use values given in \code{betaini}. If \code{0}, betaini is skipped
#'@param betaini Initial \eqn{beta_0} to use in estimation
#'@param printd option to print results of iterations for partial change model
#'@param prewhit option to use AR(1) for prewhitening
#'@param robust,hetdat,hetvar options on error terms assumptions
#' @return A list that contains following:
#' \itemize{
#'\item{supfl: SupF(l+1|l) test statistics}
#'\item{cv: Critical values for SupF(l+1|l) test}
#'\item{ndat: New date (if available)} }
#'@export

doseqtests = function(y_name,z_name=NULL,x_name=NULL,data,
                    m=5,eps=0.00001,eps1=0.15,maxi=10,fixb=0,betaini=0,printd=0,
                    prewhit=1,
                    robust=1,hetdat=1,hetvar=1,hetq=1,hetomega=1,const=1) {
  eps1 = check_trimming(eps1)

  siglev=matrix(c(10,5,2.5,1),4,1)
  df = process_data(y_name = y_name,z_name = z_name,x_name = x_name,data=data,const)
  y = df$y
  z = df$z
  x = df$x


  if (is.null(x)) {p = 0}
  else {p = dim(x)[2]}

  q = dim(z)[2]
  bigT = dim(y)[1]
  h = round(eps1*bigT)

  upper_m = floor(bigT/h)-1
  if(m>upper_m){
    warning(paste('\nNot enough observations for',m+1,'segments with minimum length per segment =',
              h,'.The total required observations for such',m,'breaks would be ',(m+1)*h,'>T=',bigT,'\n'))
    message(paste('Set m to',upper_m,'\n'))
    m=upper_m
  }

  if(m<=0){
    warning('\nMaximum number of breaks cannot be non-positive\n')
    message(paste('Set m to',upper_m,'\n'))
    m=upper_m
  }

  if(m<=1){
    out=list()
    out$mbreak = m
    class(out) = 'seqtests'
    out = compile_seqtests(out)
    return(out)
  }
  else{
    out = doglob(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd)
    datevec = out$datevec
    bigvec = out$bigvec
  supfl = matrix (0L,m,1)
  ndat = matrix (0L,m,1)
  for (i in seq(1,m-1,1)){
    out1 = spflp1(bigvec,datevec[1:i,i,drop=FALSE],i+1,y,z,h,q,prewhit,robust,x,p,hetdat,hetvar)
    supfl[i+1,1] = out1$maxf
    #print(paste('The supF(',i+1,'|',i,') test is',supfl[i,1]))
    #print(paste('It corresponds to a new break at:',ndat[i,1]))
  }
  supfl[1,1] = pftest(y,z,1,q,bigT,datevec,prewhit,robust,x,p,hetdat,hetvar)
  cv_supFl = matrix(0L,4,m)

  for (c in 1:4){
    cv = getcv2(c,eps1)
    cv_supFl[c,] = cv[q,1:m,drop=FALSE]
  }

  rownames(cv_supFl) = siglev


  out = list('supfl' = supfl, 'cv' = cv_supFl)
  out$mbreak = m
  class(out) = 'seqtests'
  out = compile_seqtests(out)
  return(out)
  }

}

#' Order estimation procedure
#'
#' The function carry out the procedure to estimate order
#'  using BIC and the criterion of Liu, Wu and Zidek
#'
#'@param y dependent variables in matrix form
#'@param z matrix of independent variables with coefficients are allowed to change across
#'regimes
#'@param x matrix of independent variables with coefficients constant across regimes
#'@param m maximum number of breaks
#'@param eps1 trimming level
#'@param eps convergence criterion for iterative recursive computation
#'@param maxi maximum number of iterations
#'@param fixb option to use fixed initial input \eqn{\beta}. If \code{1},
#'the model will use values given in \code{betaini}. If \code{0}, betaini is skipped
#'@param betaini Initial \eqn{beta_0} to use in estimation
#'@param printd option to print results of iterations for partial change model
#'@param bic indicator which criterion is used in selecting number of breaks
#'@return A list that contains following:
#'\item{mBIC}{number of breaks selected by BIC}
#'\item{mLWZ}{number of breaks selected by LWZ}
#'@export
#'@references

doorder = function(y_name,z_name = NULL,x_name = NULL,data,
                   m=5,eps=0.00001,eps1=0.15,maxi=10,fixb=0,
                   betaini=0,printd=0,opt='BIC',const=1) {

  #need to use new checker's functions:
  #check eps/check m/ check regressors for collinearity



  df = process_data(y_name = y_name,z_name = z_name,x_name = x_name,data=data,const)
  y = df$y
  z = df$z
  x = df$x


  bigT = dim(y)[1]
  #check trimming value
  eps1 = check_trimming(eps1)
  #check input number of breaks
  m = check_m(bigT,eps1,m)
  h = round(eps1*bigT)

  if (is.null(x)) {p = 0}
  else {p = dim(x)[2]}
  q = dim(z)[2]
  if (p == 0){zz = z}
  else{zz = cbind(z,x)}
  temp = doglob(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd)
  glb = temp$glb
  bigT = dim(y)[1]
  ssr0 = nssr(y,zz)
  delta0 = 0.1 #optimal parameters in LWZ paper
  c0 = 0.299
  glob= matrix(0L, nrow = m+1, ncol=1)
  glob[1,1] = ssr0
  glob[seq(2,m+1),1] = glb
  datevec = temp$datevec

  bic = matrix(0L,nrow = m+1, ncol = 1)
  lwz = matrix(0L,nrow = m+1, ncol = 1)
  kt  = matrix(0L,nrow = m+1, ncol = 1)
  for (i in seq(1,m+1)){
    #BIC criterion
    bic [i,1] = log(glob[i,1]/bigT) + log(bigT)*(i-1)*(q+1)/bigT
    #LWZ criterion
    lwz[i,1] = log(glob[i,1]/(bigT-i*q-i+1)) +
      ((i-1)*(q+1)*c0*(log(bigT))^(2+delta0))/bigT
    #Kurozumi and Tuvaandori (2011)
    if (i==1){
      bd=c(0,bigT)}
    else{
      bd=c(0,datevec[1:i-1,i-1],bigT)}
   for (l in seq(1,i)){
    segy   = y[seq(bd[l]+1,bd[l+1],1),,drop=FALSE];
    segz   = z[seq(bd[l]+1,bd[l+1],1),,drop=FALSE];
    segres = segy-segz%*%solve(t(segz)%*%segz)%*%t(segz)%*%segy;
    dt     = bd[l+1]-bd[l];
    kt[i,1]= kt[i,1]+(dt*log(t(segres)%*%segres/dt)+q*log(dt));
    }

    kt[i,1]    = kt[i,1]+2*i*log(bigT);
  }

  mBIC = which.min(bic) - 1
  mLWZ = which.min(lwz) - 1
  mKT = which.min(kt) - 1


  if (opt == 'BIC'){
    mSEL=mBIC
    p_name = 'BIC'
  }else if(opt == 'LWZ') {
    mSEL=mLWZ
    p_name = 'LWZ'
  }else if(opt == 'KT') {
    mSEL=mKT
    p_name = 'KT'
  }else{
    stop('No such criterion found. Please select either BIC, LWZ or KT')
  }

  if (mSEL == 0){
    message('\nThere are no breaks selected by ',p_name,'and estimation is skipped\n')
    out = list()
    out$p_name = p_name
    out$nbreak = mSEL
    class(out) = 'model'
    return(out)
  }
  else{
    date = temp$datevec[seq(1,mSEL,1),mSEL,drop=FALSE]
    hetq=1
    hetomega=1
    hetdat=1
    hetomega=1
    hetvar=1
    robust=1
    prewhit=1
    out = estim(mSEL,q,z,y,date,robust,prewhit,hetomega,hetq,x,p,hetdat,hetvar)
    out$p_name = p_name
    out$nbreak = mSEL
    class(out) = 'model'
    out$numz = q
    out$numx = p
    out$const = const
    out$y_name = y_name
    out$z_name = z_name
    out$x_name = x_name
    out$y = y
    out$x = x
    out$z = z
    out = compile_model(out)
    return(out)
  }

}




#' Sequential procedure
#'
#'function to apply sequential procedure to obtain number of breaks and break
#'dates. Current version only allows pure structural changes. This will be
#'generalized
#'
#'@param y dependent variables in matrix form
#'@param z matrix of independent variables with coefficients are allowed to change across
#'regimes
#'@param x matrix of independent variables with coefficients constant across regimes
#'@param m maximum number of breaks
#'@param eps1 trimming level
#'@param eps convergence criterion for iterative recursive computation
#'@param maxi maximum number of iterations
#'@param fixb option to use fixed initial input \eqn{\beta}. If \code{1},
#'the model will use values given in \code{betaini}. If \code{0}, betaini is skipped
#'@param betaini Initial \eqn{beta_0} to use in estimation
#'@param printd option to print results of iterations for partial change model
#'@param prewhit option to use AR(1) for prewhitening process
#'@param robust,hetdat,hetvar options on error terms assumptions
#' @return A list that contains following:
#' \itemize{
#'\item{nbreak}{Number of breaks}
#'\item{dateseq}{Sequence of break dates}}
#'@export

dosequa = function(y_name,z_name=NULL,x_name=NULL,data,
                   m=5,eps=0.00001,eps1=0.15,maxi=10,fixb=0,betaini=0,printd=0,
                   prewhit=1,robust=1,hetdat=1,hetvar=1,const=1,signif=2) {


  df = process_data(y_name = y_name,z_name = z_name,x_name = x_name,data=data,const)
  y = df$y
  z = df$z
  x = df$x

  if (is.null(x)) {p = 0}
  else {p = dim(x)[2]}

  q = dim(z)[2]
  bigT = dim(y)[1]
  #check trimming value
  eps1 = check_trimming(eps1)
  #check input number of breaks
  m = check_m(bigT,eps1,m)
  h = round(eps1*bigT)
  nbreak = 0
  siglev=matrix(c(10,5,2.5,1),4,1)



   # print(paste('Output from the sequential procedure at significance level',
   #              siglev[j,1],'%'))
    out_seq = sequa(m,signif,q,h,bigT,robust,prewhit,z,y,x,p,hetdat,hetvar,eps1)
    nbr = out_seq$nbreak


    #print(paste('The sequential procedure estimated the number of breaks at:',nbr))
    if (nbr > 0) {datese = as.matrix(out_seq$dv0)}
    else{cat("\nThere are no breaks selected by sequential procedure and estimation is skipped\n")
      out = list()
      out$p_name = 'dosequa'
      out$nbreak = nbr
      class(out) = 'model'
      return(out)
    }

    nbreak = nbr
    if (nbr!=0){
      dateseq = t(datese)
    }

  mSEL = nbreak

  if (mSEL == 0){
    message('\nThere are no breaks selected by sequential and estimation is skipped\n')
    out = list()
    out$p_name = 'dosequa'
    out$nbreak = mSEL
    class(out) = 'model'
    return(out)
    }
  else{
    date = dateseq
    date = t(date)
    hetq=1
    hetomega=1
    out = estim(mSEL,q,z,y,date,robust,prewhit,hetomega,hetq,x,p,hetdat,hetvar)
    out$p_name = 'dosequa'
    out$nbreak = mSEL
    class(out) = 'model'
    out$numz = q
    out$numx = p
    out$const = const
    out$y_name = y_name
    out$z_name = z_name
    out$x_name = x_name
    out$y = y
    out$x = x
    out$z = z
    out = compile_model(out)
    return(out)
  }


}



#'Repartition procedure
#'
#'The following procedure constructs the so-called repartition
#'estimates of the breaks obtained by the sequential method (see Bai
#'(1995), Estimating Breaks one at a time, Econometric Theory, 13,
#'315-352. It alows estimates that have the same asymptotic
#'distribution as those obtained by global minimization. Otherwise, the
#'output from the procedure "estim" below do not deliver asymptotically
#'correct confidence intervals for the break dates.
#'
#'@param y dependent variables in matrix form
#'@param z matrix of independent variables with coefficients are allowed to change across
#'regimes
#'@param x matrix of independent variables with coefficients constant across regimes
#'@param m maximum number of breaks
#'@param eps1 trimming level
#'@param eps convergence criterion for iterative recursive computation
#'@param maxi maximum number of iterations
#'@param fixb option to use fixed initial input \eqn{\beta}. If \code{1},
#'the model will use values given in \code{betaini}. If \code{0}, betaini is skipped
#'@param betaini Initial \eqn{beta_0} to use in estimation
#'@param printd option to print results of iterations for partial change model
#'@param prewhit option to use AR(1) for prewhitening process
#'@param robust,hetdat,hetvar options on error terms assumptions
#'@return reparv Repartition method estimation of break dates
#'
#'@export
dorepart = function(y_name,z_name = NULL,x_name = NULL,data,
                    m=5,eps=0.00001,eps1=0.15,maxi=10,fixb=0,betaini=0,printd=0,
                    prewhit=1,robust=1,hetdat=1,hetvar=1,const=1,signif=2){

  if(eps1 <0.05 || eps1 >0.5){
    warning('Invalid trimming level, set trimming level to 15%')
    eps1 = 0.15
  }

  df = process_data(y_name = y_name,z_name = z_name,x_name = x_name,data=data,const)
  y = df$y
  z = df$z
  x = df$x

  if (is.null(x)) {p = 0}
  else {p = dim(x)[2]}

  q = dim(z)[2]
  bigT = dim(y)[1]
  #check trimming value
  eps1 = check_trimming(eps1)
  #check input number of breaks
  m = check_m(bigT,eps1,m)
  h = round(eps1*bigT)

  reparv = matrix (0L,4,m)
  siglev=matrix(c(10,5,2.5,1),4,1)


  temp = sequa(m,signif,q,h,bigT,robust,prewhit,z,y,x,p,hetdat,hetvar,eps1)
  nbreak = temp$nbreak
    if (temp$nbreak == 0){
      message('\nThere are no breaks selected by sequential procedure and the repartition procedure is skipped.\n')
      out=list()
      out$p_name = 'dorepart'
      out$nbreak = 0
      class(out) = 'model'
      return(out)
    }
    else {
      repartda = preparti(y,z,temp$nbreak,
                          temp$dv0,
                          h,x,p)
      reparv = repartda
    }


  #estimate the date at 5% significant level
  mSEL = nbreak


  if (mSEL == 0){
    message('\nThere are no breaks selected by sequential procedure and estimation is skipped\n')
    out = list()
    out$p_name = 'dorepart'
    out$nbreak = mSEL
    class(out) = 'model'
    return(out)
  }
  else{
    date = reparv
    hetq=1
    hetomega=1
    out = estim(mSEL,q,z,y,date,robust,prewhit,hetomega,hetq,x,p,hetdat,hetvar)
    out$p_name = 'dorepart'
    out$nbreak = mSEL
    class(out) = 'model'
    out$numz = q
    out$numx = p
    out$const = const
    out$y_name = y_name
    out$z_name = z_name
    out$x_name = x_name
    out$y = y
    out$x = x
    out$z = z
    out = compile_model(out)
    return(out)
  }


  return (reparv)
}

#'Estimate a model with pre-specified number of breaks
#'
#'The following procedure constructs the so-called repartition
#'estimates of the breaks obtained by the sequential method (see Bai
#'(1995), Estimating Breaks one at a time, Econometric Theory, 13,
#'315-352. It alows estimates that have the same asymptotic
#'distribution as those obtained by global minimization. Otherwise, the
#'output from the procedure "estim" below do not deliver asymptotically
#'correct confidence intervals for the break dates.
#'
#'@param y dependent variables in matrix form
#'@param z matrix of independent variables with coefficients are allowed to change across
#'regimes
#'@param x matrix of independent variables with coefficients constant across regimes
#'@param m maximum number of breaks
#'@param eps1 trimming level
#'@param eps convergence criterion for iterative recursive computation
#'@param maxi maximum number of iterations
#'@param fixb option to use fixed initial input \eqn{\beta}. If \code{1},
#'the model will use values given in \code{betaini}. If \code{0}, betaini is skipped
#'@param betaini Initial \eqn{beta_0} to use in estimation
#'@param printd option to print results of iterations for partial change model
#'@param prewhit option to use AR(1) for prewhitening process
#'@param robust,hetdat,hetvar options on error terms assumptions
#'@return reparv Repartition method estimation of break dates
#'
#'@export
#'
#Estimate a model with pre-specified number of breaks
dofix = function(y_name,z_name = NULL,x_name=NULL,data,
                    fixn=5,eps=0.00001,eps1=0.15,maxi=10,fixb=0,betaini=0,printd=0,
                    prewhit=1,robust=1,hetdat=1,hetvar=1,hetq=1,hetomega=1,const=1){



  if(fixn<0){
    warning('\nThe maximum number of breaks cannot be negative, set prespecified breaks = 2\n')
    fixn = 2
  }

  df = process_data(y_name = y_name,z_name = z_name,x_name = x_name,data=data,const)
  y = df$y
  z = df$z
  x = df$x


  if (is.null(x)) {p = 0}
  else {p = dim(x)[2]}

  q = dim(z)[2]
  bigT = dim(y)[1]
  #check trimming value
  eps1 = check_trimming(eps1)
  #check input number of breaks
  fixn = check_m(bigT,eps1,fixn)
  h = round(eps1*bigT)
  t_out = doglob(y,z,x,fixn,eps,eps1,maxi,fixb,betaini,printd)
  datevec = t_out$datevec
if(length(datevec) == 0){
  message('\nThere are no breaks selected by the procedure\n')
  out = list()
  out$p_name = 'dofix'
  out$nbreak = length(datevec)
  class(out) = 'model'
  return(out)
}else{
  date = datevec[,fixn,drop=FALSE]
  hetq=1
  hetomega=1
  fix_mdl = estim(fixn,q,z,y,date,robust,prewhit,hetomega,hetq,x,p,hetdat,hetvar)
out = fix_mdl
out$p_name = 'fix'
out$nbreak = fixn
class(out) = 'model'
out$numz = q
out$numx = p
out$const = const
out$y_name =y_name
out$x_name =x_name
out$z_name =z_name
out$y = y
out$z = z
out$x = x
out = compile_model(out)
return(out)}
}


#Format output of n break model
#'
#'@noRd
compile_model = function (x,digits = -1,...){
  if (digits==-1){digits = 3}
  if (x$nbreak == 0){return(NULL)}
  else {

    #format date estimation
    CI_95 = c()
    CI_90 = c()
    coln = c()
    for (i in 1:x$nbreak){
      coln = c(coln, paste('Break',i,sep=''))
      CI_95 = c(CI_95,paste('(',x$CI[i,1],',',x$CI[i,2],')',sep=''))
      CI_90 = c(CI_90,paste('(',x$CI[i,3],',',x$CI[i,4],')',sep=''))
    }
    date_tab = data.frame(date = t(x$date),stringsAsFactors = FALSE)
    date_tab = rbind(date_tab,CI_95)
    date_tab = rbind(date_tab,CI_90)
    colnames(date_tab) = coln
    rownames(date_tab) = c('Date','95% CI','90% CI')


    #format full-sample coefficients
    if (x$const == 1){
      rnameRS = 'Const'
      for (i in 1:x$numz-1){
        rnameRS = cbind(rnameRS,x$z_name[i])
      }
    }else{
      rnameRS = c()
      for (i in 1:x$numz){
        rnameRS = cbind(rnameRS,x$z_name[i])
      }
    }

    cnameRS = c()
    coefRS = c()
    for (i in 1:(x$nbreak+1)){
    cnameRS = cbind(cnameRS,paste('Regime',i))
    }
    for (j in 1:x$numz){
      coefRSj = c()
    for (i in 1:(x$nbreak+1)){
      coefRSj = cbind(coefRSj,paste(format(round(x$beta[(j-1)*(x$nbreak+1)+i,1],digits),nsmall=digits),
                              paste('(',format(round(x$SE[(j-1)*(x$nbreak+1)+i,1],digits),nsmall=digits),')',sep='')))
      }
      coefRS = rbind(coefRS,coefRSj)
    }

  }
    rnameRSf = paste(rnameRS,'(SE)',sep=' ')
    coef_tabRS=data.frame(co = coefRS)


    rownames(coef_tabRS) = rnameRSf
    colnames(coef_tabRS) = cnameRS


    #format full sample coefficients
    if(x$numx == 0){coef_tabRW = NULL}
    else{
      rnameRW = c()
      for (i in 1:x$numx){
        rnameRW = cbind(rnameRW,x$x_name[i])
      }

    cnameRW = 'Full sample'
    coefRW = c()
    for (j in 1:x$numx){
      coefRW = cbind(coefRW,paste(format(round(x$beta[x$numz*(x$nbreak+1)+j,1],digits),nsmall=digits),
                                  paste('(',format(round(x$SE[x$numz*x$nbreak+j,1],digits),nsmall=digits),')',sep='')))}

    rnameRWf = paste(rnameRW,'(SE)',sep=' ')

    coef_tabRW=data.frame(co = coefRW)
    colnames(coef_tabRW) = rnameRWf
    rownames(coef_tabRW) = cnameRW}

    table = list('date_tab' = date_tab,'RS_tab' = coef_tabRS, 'FS_tab' = coef_tabRW)
    x$tab = table
    return(x)
  }



#'Summary output of a n breaks model
#'
#'Function to format the output of the n-break model
#'@export

print.model <- function(x,...)
{
  #print procedure used to select number of breaks
  proc = switch(x$p_name,'dosequa'='sequential procedure', 'BIC' = 'BIC', 'LWZ' = 'LWZ',
                'KT'='KT',
                'dorepart' = 'repartition procedure', 'fix' = 'specified number of breaks')
  digits = max(3L, getOption("digits") - 3L)

  if (x$nbreak == 0){
    cat(paste('\nNo breaks were found using',proc),'\n')
  }else{
    cat(paste('\nThe number of breaks is estimated by',proc,'\n'))
  if(x$numx == 0){
    cat(paste('Pure change model with',x$nbreak,'estimated breaks.',sep=' '))
  }else if(x$numx > 0){
    cat(paste('Partial change model with', x$nbreak,'estimated breaks.',sep=' '))
  }
  cat('\nMinimum SSR =',
              format(round(x$SSR,3),nsmall=3),'\n')

  cat('\nEstimated date:\n')
  print(x$tab$date_tab,quote=FALSE)

  cat('\nEstimated regime-specific coefficients:\n')
  print(x$tab$RS_tab,quote=FALSE)


  if(x$numx == 0) {cat('\nNo full sample regressors\n')}
  else{
    cat('\nEstimated full-sample coefficients:\n\n')
    print(x$tab$RW_tab,quote='FALSE')}}


  invisible(x)
}


#'Summary output of Sup Wald test
#'
#'Function to format the output of the Sup F test
#'@export


compile_sbtests <- function(x,digits = -1,...)
{
  if(x$mbreak == 0){
    return(x)
  }
  else{
  cnames1 = c()
  for (i in 1:x$mbreak){
    if (i == 1){
    cnames1 = cbind(cnames1, paste(i,'break'))}
    else
    {
      cnames1 = cbind(cnames1,paste(i,'breaks'))
    }
  }
  ftest = t(x$ftest)


  supF1 = data.frame(ftest = format(round(ftest,3),nsmall=3))
  cv_supF = format(round(x$cv_supF,3),nsmall = 3)
  colnames(cv_supF) = colnames(supF1)
  supF1 = rbind(supF1,cv_supF)

  rownames(supF1) = c('Sup F','10% CV','5% CV','2.5% CV','1% CV')
  colnames(supF1) = cnames1

  UDmax = data.frame(UDmax = format(round(x$UDmax,3),nsmall = 3),
                     cv = format(round(t(x$cv_Dmax),3),nsmall = 3))

  colnames(UDmax) = c('UDMax','10% CV','5% CV','2.5% CV','1% CV')


  x$supF1 = supF1
  x$UDmax = UDmax
  return(x)}
}


#'S3 print function for sup F tests
#'#'Function to format the output of the Sup F test
#'@export

print.sbtests <- function(x,...)
{ if(x$mbreak == 0){
  warning('\nThe test is undefined for no break model\n')
}else{
  cat('\na) SupF tests against a fixed number of breaks\n\n')
  print(x$supF1,quote=FALSE)
  cat('\nb) UDmax tests against an unknown number of breaks\n\n')
  print(x$UDmax,quote=FALSE)}
  invisible(x)
}


compile_seqtests = function(x){
  if(x$mbreak==1){
    #message('\nThe test is exactly 0 versus 1 break, hence the sequential test is not repeated\n')
    #x$sfl = NULL
    return(x)
  }else if (x$mbreak==0){
    warning('\nThe test is undefined for maximum break = 0\n')
    x$sfl = NULL
    return(x)
  }
  else{
  nbreak = x$mbreak
  cnames = c()
  for (i in seq(0,nbreak-1)){
    cnames = cbind(cnames, paste('supF(',i+1,'|',i,')',sep=''))
  }
  temp_supfl = format(round(x$supfl,3),nsmall = 3)
  #temp_ndat = format(x$ndat,nsmall=0)
  temp_cv = format(round(x$cv,3),nsmall = 3)

  sfl =data.frame(supfl = t(temp_supfl[seq(1,nbreak,1),1,drop=FALSE]))
  #ndat = t(temp_ndat[seq(1,nbreak,1),1,drop=FALSE])
  cv = temp_cv[,seq(1,nbreak,1),drop=FALSE]
  #colnames(ndat) = colnames(sfl)
  colnames(cv) = colnames(sfl)
  sfl = rbind(sfl,cv)
  colnames(sfl) = cnames
  rownames(sfl) = c('Seq supF','10% CV','5% CV', '2.5% CV', '1% CV')
  x$sfl = sfl
  return (x)}
}

#'S3 function to print sequential tests
#'@export
print.seqtests = function(x,...){
  if(x$mbreak==1){
      cat('\nThe test is exactly 0 versus 1 break, hence the sequential test is not repeated\n')
    }else if (x$mbreak==0){
      cat('\nThe test is undefined for maximum break = 0\n')
    }
  else{
  cat('\nsupF(l+1|l) tests using global optimizers under the null\n\n')
  print(x$sfl,quote=FALSE)}
  invisible(x)
}



#' new plot function for class model
#' @importFrom ggplot2 ggplot aes annotate geom_segment geom_line ggtitle .data coord_cartesian scale_color_manual geom_ribbon
#' @export
plot_model = function(model,CI=0.95,title=NULL){
  m = model$nbreak
  if(m==0){
    warning('The model has no break. Visualization for comparison between structural breaks
        versus no breaks is skipped')
    return(NULL)
  }


  zreg = model$z
  xreg = model$x
  p = model$numx
  q = model$numz
  date = model$date
  beta = model$beta

  #comparison between structural break vs no break
  y = model$y
  ypred_break = model$fitted
  fixreg = cbind(xreg,zreg)
  fixbeta = OLS(y,fixreg)
  ypred_fix = fixreg%*%fixbeta
  x_t = seq(1,dim(y)[1],1)
  tbl = data.frame(x_t,y,ypred_break,ypred_fix,stringsAsFactors = TRUE)
  colnames(tbl) = c('time','y','ypred_break','ypred_fix')


  #labels and annotations for date's CIs
  date_lab = c()
  for (i in 1:m){
    date_lab = c(date_lab,paste('Date',i))
  }
  vline_seg = data.frame(date,date,rep(Inf,m),rep(-Inf,m))
  colnames(vline_seg) = c('x','xend','y','yend')
  y_pos=c()
  for (i in 1:m){
    y_pos = rbind(y_pos,(10.2+i/2)/10*min(y))
  }
  model$CI[model$CI<0] = 0
  model$CI[model$CI>dim(y)[1]] = dim(y)[1]
  CI_seg95 = data.frame(model$CI[,1],model$CI[,2],y_pos)
  CI_seg90 = data.frame(model$CI[,3],model$CI[,4],y_pos)

  #compute CIs for estimation y
  zbar = diag_par(zreg,m,date)
  if (p == 0){
    reg = zbar
  }
  else{
    reg = cbind(xreg,zbar)
  }

  if(!CI==0.95&&!CI==0.90){
    warning('Not available CI level, set to 95%')
    CI = 0.95
  }
  if(CI == 0.95){
   CI_seg = CI_seg95
   sd = 1.960*model$SE
   beta_lb = beta-sd
   beta_ub = beta+sd
  }else if (CI==0.90){
    CI_seg = CI_seg90
    sd = 1.645*model$SE
    beta_lb = beta-sd
    beta_ub = beta+sd
  }
  tbl$ypred_ub = reg%*%beta_ub
  tbl$ypred_lb = reg%*%beta_lb
  colnames(CI_seg) = c('lb','ub','y')

  if(is.null(title)){
  if(model$numx==0){
    if (m>1){
    title = paste('Pure change with',m,'breaks')}
    else{
      title = paste('Pure change with',m,'break')
    }
  }else{
    if (m>1){
      title = paste('Partial change with',m,'breaks')}
    else{
      title = paste('Partial change with',m,'break')
    }
  }}

  grph = ggplot2::ggplot(data = tbl, ggplot2::aes(x=.data$time))+
    ggplot2::coord_cartesian(xlim=c(0,max(x_t)+1))+
    ggplot2::geom_line(ggplot2::aes(y=.data$ypred_break,color='y_break'),size=0.3)+
    ggplot2::geom_line(ggplot2::aes(y=.data$ypred_fix,color='y_fix'),size=0.3)+
    ggplot2::geom_line(ggplot2::aes(y=.data$y,color='y'),size=0.2)+
    ggplot2::geom_ribbon(ggplot2::aes(ymax=.data$ypred_ub, ymin=.data$ypred_lb), fill="gray", alpha=.35)+
    #ggplot2::geom_ribbon(ggplot2::aes(ymax=.datax.upper, ymin=x.lower), fill="pink", alpha=.5)
    ggplot2::scale_color_manual(name=paste('Legends'),
                                breaks = c('y','y_break','y_fix'),
                                values = c('y'='black','y_break' = 'blue', 'y_fix' = 'red'),
      labels = c(expression(y),expression(hat(y)[m]), expression(hat(y)[0]))
      )+
    ggplot2::annotate('text',x = model$date[,1]-3*max(x_t)/100, y=10/11*max(y),label= date_lab,angle = 90)+
    ggplot2::geom_segment(data=vline_seg,
                          ggplot2::aes(x=.data$x,y=.data$y,
                                       xend=.data$xend,yend=.data$yend),
                          alpha =0.85,colour='purple',linetype='dashed')+
    ggplot2::geom_segment(data=CI_seg,
                          ggplot2::aes(x=.data$lb,xend=.data$ub,y=.data$y,yend=.data$y),
                          colour='red',alpha=0.6)+
  ggplot2::ggtitle(title)+ggplot2::ylab(model$y_name)
  grph
}


