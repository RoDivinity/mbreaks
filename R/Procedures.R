### Main procedures

#' Global SSR minimizers procedure
#'
#' Function obtains global minimizers using the recursive algorithm to compute
#' and minimize SSR over all possible segments. The procedure is required to
#' conduct supF, UDMax, WDMax and supF(l+1|l) test
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
#'\item{bigvec} {Associated SSRs}}
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
#'
#'
dotest = function(y_name,z_name=NULL,x_name=NULL,data,
                  m=5,eps=0.00001,eps1=0.15,maxi=10,fixb=0,betaini=0,printd=0,prewhit=1,robust=1,
                  hetdat=1,hetvar=1){
  siglev=matrix(c(10,5,2.5,1),4,1)
  df = process_data(y_name = y_name,z_name = z_name,x_name = x_name,data=data)
  y = df$y
  z = df$z
  x = df$x
  if(m<0){
    cat('\nThe maximum number of breaks cannot be negative, set m = 5\n')
    m = 5
  }
  
  if(m==0){
    cat('The test is undefined for no breaks model')
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
    
    upper_m = floor(bigT/h)-1
    if(m>upper_m){
      cat(paste('\nNot enough observations for',m+1,'segments with minimum length per segment =',
                h,'.The total required observations for such',m,'breaks would be ',(m+1)*h,'>T=',bigT,'\n'))
      cat(paste('Set m to 5\n'))
      m=5
    }
  
  
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
  for (c in 1:4){
    #critical values for supF test
    cv = getcv1(c,eps1)
    cv_supF[c,] = cv[q,1:m,drop=FALSE]
    if (printd==1){
    print(paste('The critical values at the',siglev[c,1],'% level are (for k = 1 to',m,'):'))
    print(cv[q,1:m,drop=FALSE])}
  }

  #procedure for Dmax and UDmax test

  #print('b) Dmax test against an unknown number of breaks')
  #print(paste('The UDmax test is:',max(ftest)))
  cv_Dmax = matrix(0L,4,1)
  for (c in 1:4) {
    #critical values for Dmax test
    cvm = getdmax(c,eps1)
    cv_Dmax[c,1] = cvm[q,1]
    if(printd==1){
    print(paste('The critical values at the',siglev[c,1],'% level is:',
                cvm[q,1]))}
  }


  for (c in 1:4){
    #computation of WDmax test
    cv = getcv1(c,eps1)
    for( i in 1:m){
      wftest[i,1] = cv[q,1] * ftest[i,1] / cv[q,i]
    }
    if (printd==1){
    print(paste('WDmax test at the',siglev[c,1],'% level is:',max(wftest)))}
  }
  rownames(cv_supF) = siglev
  rownames(cv_Dmax) = siglev
  
 
  out = list('ftest' = ftest, 'cv_supF' = cv_supF,
             'cv_Dmax' = cv_Dmax, 'UDmax' = max(ftest))
  out$mbreak = m
  class(out) = 'sbtests'
  
  out = compile.sbtests(out)
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
#'
#'@export
doseqtests = function(y_name,z_name=NULL,x_name=NULL,data,
                    m=5,eps=0.00001,eps1=0.15,maxi=10,fixb=0,betaini=0,printd=0,
                    prewhit=1,
                    robust=1,hetdat=1,hetvar=1) {
  
  siglev=matrix(c(10,5,2.5,1),4,1)
  df = process_data(y_name = y_name,z_name = z_name,x_name = x_name,data=data)
  y = df$y
  z = df$z
  x = df$x
  
  
  if(m<0){
    cat('\nThe maximum number of breaks cannot be negative, set m = 5\n')
    m = 5
  }

  if (is.null(x)) {p = 0}
  else {p = dim(x)[2]}

  q = dim(z)[2]
  bigT = dim(y)[1]
  h = round(eps1*bigT)
  
  upper_m = floor(bigT/h)-1
  if(m>upper_m){
    cat(paste('\nNot enough observations for',m+1,'segments with minimum length per segment =',
              h,'.The total required observations for such',m,'breaks would be ',(m+1)*h,'>T=',bigT,'\n'))
    cat(paste('Set m to 5\n'))
    m=5
  }
  
  if(m<=1){
    out=list()
    out$mbreak = m
    class(out) = 'seqtests'
    out = compile.seqtests(out)
    return(out)
  }
  else{
    out = doglob(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd)
    datevec = out$datevec
    bigvec = out$bigvec
  supfl = matrix (0L,m-1,1)
  ndat = matrix (0L,m-1,1)
  for (i in seq(1,m-1,1)){
    out1 = spflp1(bigvec,datevec[1:i,i,drop=FALSE],i+1,y,z,h,q,prewhit,robust,x,p,hetdat,hetvar)
    supfl[i,1] = out1$maxf
    ndat[i,1] = out1$newd
    #print(paste('The supF(',i+1,'|',i,') test is',supfl[i,1]))
    #print(paste('It corresponds to a new break at:',ndat[i,1]))
  }

  cv_supFl = matrix(0L,4,m)

  for (c in 1:4){
    cv = getcv2(c,eps1)
    cv_supFl[c,] = cv[q,1:m,drop=FALSE]
  }
  rownames(cv_supFl) = siglev

  
  out = list('supfl' = supfl, 'ndat' = ndat, 'cv' = cv_supFl)
  out$mbreak = m
  class(out) = 'seqtests'
  out = compile.seqtests(out)
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
#


doorder = function(y_name,z_name = NULL,x_name = NULL,data,
                   m=5,eps=0.00001,eps1=0.15,maxi=10,fixb=0,
                   betaini=0,printd=0,bic_opt=1) {

  df = process_data(y_name = y_name,z_name = z_name,x_name = x_name,data=data)
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
    cat(paste('\nNot enough observations for',m+1,'segments with minimum length per segment =',
              h,'.The total required observations for such',m,'breaks would be ',(m+1)*h,'>T=',bigT,'\n'))
    cat(paste('Set m to 5\n'))
    m=5
  }
  
  if (m<=0){
    cat('\nMaximum number of breaks cannot be less than 1, m is set to 5\n')
    m=5
  }
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

  bic = matrix(0L,nrow = m+1, ncol = 1)
  lwz = matrix(0L,nrow = m+1, ncol = 1)

  for (i in seq(1,m+1)){
    bic [i,1] = log(glob[i,1]/bigT) + log(bigT)*(i-1)*(q+1)/bigT
    lwz[i,1] = log(glob[i,1]/(bigT-i*q-i+1)) +
      ((i-1)*(q+1)*c0*(log(bigT))^(2+delta0))/bigT

   
  }

  mBIC = which.min(bic) - 1
  mLWZ = which.min(lwz) - 1
  
  if (bic_opt == 1){
    mSEL=mBIC
    p_name = 'BIC'
  }else{
    mSEL=mLWZ
    p_name = 'LWZ'
  }
  
  if (mSEL == 0){
    cat('\nThere are no breaks selected by ',p_name,'and estimation is skipped\n')
    out = list()
    out$p_name = p_name
    out$nbreak = mSEL
    class(out) = 'model'
    return(out)
  }
  else{
    date = temp$datevec[seq(1,mSEL,1),mSEL,drop=FALSE]
    hetq=0
    hetomega=0
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
    out$y_name = y_name
    out$z_name = z_name
    out$x_name = x_name
    out$y = y
    out$x = x
    out$z = z
    out = compile.model(out)
    return(out)
  }
  
}


#sequential procedure
sequa = function(m,signif,q,h,bigT,robust,prewhit,z,y,x,p,hetdat,hetvar,eps1){

  dv = matrix(0L, nrow = m+2, ncol = 1)
  dv2 = matrix(0L, nrow = m+2, ncol = 1)
  ftestv = matrix(0L, nrow = m+1,ncol = 1)

  cv = getcv2(signif,eps1)
  dv[1,1] = 0

  if (p == 0){
    y_rev = rot90(rot90(y))
    z_rev = rot90(rot90(z))
    vssrev = ssr(1,y_rev,z_rev,h,bigT)
    vssr = ssr(1,y,z,h,bigT)
    out = partione(h,bigT-h,bigT,vssr,vssrev)
    datx = out$dx
    ssrmin = out$ssrmin
  }
  else{
    out = onebp(y,z,x,h,1,bigT)
    datx = out$bd
    ssrmin = out$ssrind
  }

  dv[2,1] = datx

  ftest=pftest(y,z,1,q,bigT,dv[2,1,drop=FALSE],prewhit,robust,x,p,hetdat,hetvar)

  if (ftest < cv[q,1]) {
    nbreak = 0
    dv[2,1] = 0
    #dv0 = dv[seq(2,nbreak+1,1),1]
    nseg = 1
  }
  else{
    #print(paste('First break found at:',datx))
    nbreak = 1
    nseg = 2
    dv[nseg+1,1] = bigT
  }

  while(nseg <= m){
    ds = matrix(0L,nseg+1,1)
    ftestv = matrix(0L,nseg+1,1)

    i_s = 1

    while(i_s <= nseg){
      length = dv[i_s+1,1] - dv[i_s,1]

      if(length >= 2*h){
        if(p==0){
          y_temp = y[seq(dv[i_s,1]+1,dv[i_s+1,1]),1,drop=FALSE]
          z_temp = z[seq(dv[i_s,1]+1,dv[i_s+1,1]),,drop=FALSE]
          vssr = ssr(1,y_temp,z_temp,h,length)
          y_temp_rev = rot90(rot90(y_temp))
          z_temp_rev = rot90(rot90(z_temp))
          vssrev = ssr(1,y_temp_rev,z_temp_rev,h,length)
          out = partione(h,length-h,length,vssr,vssrev)
          ds[i_s,1] = out$dx
          ftestv[i_s,1] = pftest(y_temp,z_temp,1,q,length,ds[i_s,1,drop=FALSE],prewhit,
                                 robust,0,p,hetdat,hetvar)
        }
        else{
          y_temp = y[seq(dv[i_s,1]+1,dv[i_s+1,1]),1,drop=FALSE]
          z_temp = z[seq(dv[i_s,1]+1,dv[i_s+1,1]),,drop=FALSE]
          x_temp = x[seq(dv[i_s,1]+1,dv[i_s+1,1]),,drop=FALSE]
          out = onebp(y,z,x,h,dv[i_s,1]+1,dv[i_s+1,1])
          ds[i_s,1] = out$bd - dv[i_s,1]
          ftestv[i_s,1] = pftest(y_temp,z_temp,1,q,length,ds[i_s,1],
                                 prewhit,robust,x_temp,p,hetdat,hetvar)
        }
      }
      else{
        ftestv[i_s,1] = 0.0
      }
      i_s = i_s + 1
    }

    maxf = max(ftestv[seq(1,nseg,1),1])

    if (maxf < cv[q,nseg]){
      #print(nbreak)
      #dv0 = dv[seq(2,nbreak+1,1),1]
    }
    else {
      newseg = which.max(ftestv[seq(1,nseg),1])
      dv[nseg+2,1] = ds[newseg,1] + dv[newseg,1]
      nbreak = nbreak + 1
      #check this sort
      dv2 = sort(dv[seq(2,nseg+2,1),1])
      dv2 = matrix(dv2, ncol = 1)
      dv[1,1] = 0
      dv[seq(2,nseg+2,1),1] = dv2

    }
    nseg = nseg + 1
  }

  #print('The sequential procedure has reached the upper limit')
  if (nbreak < 1) {dv0 = c()}
  else{
    dv0 = dv[seq(2,nbreak+1,1),1]}
  out = list('nbreak' = nbreak, 'dv0' = dv0)
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
                   prewhit=1,robust=1,hetdat=1,hetvar=1) {
  if (m<=0){
    cat('\nThe maximum number of breaks cannot be negative, set m = 5\n')
    m = 5;
  }
  
  
  df = process_data(y_name = y_name,z_name = z_name,x_name = x_name,data=data)
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
    cat(paste('\nNot enough observations for',m+1,'segments with minimum length per segment =',
              h,'.The total required observations for such',m,'breaks would be ',(m+1)*h,'>T=',bigT,'\n'))
    cat(paste('Set m to 5\n'))
    m=5
  }
  
  nbreak = matrix(0L, nrow = 4, ncol = 1)
  dateseq = matrix(0L,nrow = 4, ncol = m)
  siglev=matrix(c(10,5,2.5,1),4,1)

  
  for (j in 1:4){
   # print(paste('Output from the sequential procedure at significance level',
   #              siglev[j,1],'%'))
    out_seq = sequa(m,j,q,h,bigT,robust,prewhit,z,y,x,p,hetdat,hetvar,eps1)
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
    nbreak[j,1] =nbr

    if (nbr!=0){
      dateseq[j,seq(1,nbreak[j,1])] = t(datese)
    }
  }
  mSEL = nbreak[2,1]
  
  if (mSEL == 0){
    cat('\nThere are no breaks selected by sequential and estimation is skipped\n')
    out = list()
    out$p_name = 'dosequa'
    out$nbreak = mSEL
    class(out) = 'model'
    return(out)
    }
  else{
    date = as.matrix(dateseq[2,seq(1,nbreak[2,1],1),drop=FALSE])
    date = t(date)
    hetq=0
    hetomega=0
    out = estim(mSEL,q,z,y,date,robust,prewhit,hetomega,hetq,x,p,hetdat,hetvar)
    out$p_name = 'dosequa'
    out$nbreak = mSEL
    class(out) = 'model'
    out$numz = q
    out$numx = p
    out$y_name = y_name
    out$z_name = z_name
    out$x_name = x_name
    out$y = y
    out$x = x
    out$z = z
    out = compile.model(out)
    return(out)
  }
  
  
}

#preparti
preparti = function(y,z,nbreak,dateseq,h,x,p) {
  bigT = dim(z)[1]
  q = dim(z)[2]

  #care if nbreak is matrix or scalar
  dv = matrix(0L, nrow = nbreak+2, ncol = 1)
  dv[1,1] = 0
  dv[seq(2,nbreak+1,1),1] = dateseq

  dv[nbreak+2,1] = bigT
  ds = matrix(0L,nrow = nbreak, ncol = 1)
  dr = matrix(0L,nrow = nbreak, ncol = 1)

  for (is in 1:nbreak){
    length = dv[is+2,1] - dv[is,1]
    if (p == 0){
      index = seq(dv[is,1]+1,dv[is+2,1],1)
      y_temp = y[index,1,drop=FALSE]
      z_temp = z[index,,drop=FALSE]
      vssr = ssr(1,y_temp,z_temp,h,length)
      y_temp_rev = rot90(rot90(y_temp))
      z_temp_rev = rot90(rot90(z_temp))
      vssrev = ssr(1,y_temp_rev,z_temp_rev,h,length)
      out = partione(h,length-h,length,vssr,vssrev)
      ds[is,1] = out$dx
      dr[is,1] = ds[is,1] + dv[is,1]
    }
    else{
      out = onebp(y,z,x,h,dv[is,1]+1,dv[is+2,1])
      ds[is,1] = out$bd
      dr[is,1] = ds[is,1]
    }

  }
  return(dr)
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
dorepart = function(y_name,z_name,x_name,data,
                    m=5,eps=0.00001,eps1=0.15,maxi=10,fixb=0,betaini=0,printd=0,
                    prewhit=1,robust=1,hetdat=1,hetvar=1){
  
  if(m<0){
    cat('\nThe maximum number of breaks cannot be negative, set m = 5\n')
    m = 5
  }
  
  df = process_data(y_name = y_name,z_name = z_name,x_name = x_name,data=data)
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
    cat(paste('\nNot enough observations for',m+1,'segments with minimum length per segment =',
              h,'.The total required observations for such',m,'breaks would be ',(m+1)*h,'>T=',bigT,'\n'))
    cat(paste('Set m to 5\n'))
    m=5
  }
  
  reparv = matrix (0L,4,m)
  siglev=matrix(c(10,5,2.5,1),4,1)
  
  nbreak = matrix(0L,4,1)
  
  
  for (j in 1:4){
    temp = sequa(m,j,q,h,bigT,robust,prewhit,z,y,x,p,hetdat,hetvar,eps1)
    nbreak[j,1] = temp$nbreak
    if (temp$nbreak == 0){
      cat('\nThere are no breaks selected by sequential procedure and the repartition procedure is skipped.\n')
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
      reparv[j,seq(1:temp$nbreak)] = repartda
    }
  }
  
  #estimate the date at 5% significant level
  mSEL = nbreak[2,1]
  
  
  if (mSEL == 0){
    cat('\nThere are no breaks selected by sequential procedure and estimation is skipped\n')
    out = list()
    out$p_name = 'dorepart'
    out$nbreak = mSEL
    class(out) = 'model'
    return(out)
  }
  else{
    date = as.matrix(reparv[2,seq(1,nbreak[2,1],1),drop=FALSE])
    date = t(date)
    hetq=0
    hetomega=0
    out = estim(mSEL,q,z,y,date,robust,prewhit,hetomega,hetq,x,p,hetdat,hetvar)
    out$p_name = 'dorepart'
    out$nbreak = mSEL
    class(out) = 'model'
    out$numz = q
    out$numx = p
    out$y_name = y_name
    out$z_name = z_name
    out$x_name = x_name
    out$y = y
    out$x = x
    out$z = z
    out = compile.model(out)
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
dofix = function(y_name,z_name,x_name,data,
                    fixn=5,eps=0.00001,eps1=0.15,maxi=10,fixb=0,betaini=0,printd=0,
                    prewhit=1,robust=1,hetdat=1,hetvar=1){
  
  if(fixn<0){
    cat('\nThe maximum number of breaks cannot be negative, set prespecified breaks = 2\n')
    fixn = 2
  }
  
  df = process_data(y_name = y_name,z_name = z_name,x_name = x_name,data=data)
  y = df$y
  z = df$z
  x = df$x
  
  if (is.null(x)) {p = 0}
  else {p = dim(x)[2]}
  
  q = dim(z)[2]
  bigT = dim(y)[1]
  h = round(eps1*bigT)
  
  upper_m = floor(bigT/h)-1
  if(fixn>upper_m){
    cat(paste('\nNot enough observations for',fixn+1,'segments with minimum length per segment =',
              h,'.The total required observations for such',fixn,'breaks would be ',(fixn+1)*h,'>T=',bigT,'\n'))
    cat(paste('Set m to 5\n'))
    fixn=5
  }
  
t_out = doglob(y,z,x,fixn,eps,eps1,maxi,fixb,betaini,printd)
datevec = t_out$datevec
if(length(datevec) == 0){
  cat('\nThere are no breaks selected by the procedure\n')
  out = list()
  out$p_name = 'dofix'
  out$nbreak = length(datevec)
  class(out) = 'model'
  return(out)
}else{
  date = datevec[,fixn,drop=FALSE]
  hetq=0
  hetomega=0
  fix_mdl = estim(fixn,q,z,y,date,robust,prewhit,hetomega,hetq,x,p,hetdat,hetvar)
out = fix_mdl
out$p_name = 'fix'
out$nbreak = fixn
class(out) = 'model'
out$numz = q
out$numx = p
out$y_name =y_name
out$x_name =x_name
out$z_name =z_name
out$y = y
out$z = z
out$x = x
out = compile.model(out)
return(out)}
}







#Format output of n break model
#'
#'@noRD
compile.model = function (x,digits = -1,...){
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
    rnameRS = c()
      for (i in 1:x$numz){
        rnameRS = cbind(rnameRS,x$z_name[i])
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
    rnameRS = paste(rnameRS,'(SE)',sep=' ')
    coef_tabRS=data.frame(co = coefRS)
    
    
    rownames(coef_tabRS) = rnameRS
    colnames(coef_tabRS) = cnameRS
    
    
    #format regime-wise coefficients
    
    if(x$numx == 0){coef_tabRW = NULL}
    else{ 
      rnameRW = c()
      for (i in 1:x$numx){
        rnameRW = cbind(rnameRW,x$x_name[i])
      }
    
    cnameRW = 'All regimes'
    coefRW = c()
    for (j in 1:x$numx){
      coefRW = cbind(coefRW,paste(format(round(x$beta[x$numz*(x$nbreak+1)+j,1],digits),nsmall=digits),
                                  paste('(',format(round(x$SE[x$numz*x$nbreak+j,1],digits),nsmall=digits),')',sep='')))}
    
    rnameRW = paste(rnameRW,'(SE)',sep=' ')
    coef_tabRW=data.frame(co = coefRW)
    
    rownames(coef_tabRW) = rnameRW
    colnames(coef_tabRW) = cnameRW}
    
    table = list('date_tab' = date_tab,'RS_tab' = coef_tabRS, 'RW_tab' = coef_tabRW)
    x$tab = table
    return(x)
  }
  
  




#'Summary output of a n breaks model
#'
#'Function to format the output of the n-break model
#'@param nbreak number of breaks in the model
#'@param glb global minimum SSR of the model
#'@param datevec estimated break date
#'@param y dependent variables
#'@param z independent variables with regime-specific coefficients
#'@param x independent variables with regime-wise coefficients
#'@return tbl Tables containings:
#'i) global minimum SSR
#'ii) estimated date
#'iii) estimated coefficients

print.model <- function(x,digits = -1,...)
{
  #print procedure used to select number of breaks
  proc = switch(x$p_name,'dosequa'='sequential procedure', 'BIC' = 'BIC', 'LWZ' = 'LWZ',
                'dorepart' = 'repartition procedure', 'fix' = 'specified number of breaks')

  
  if (digits==-1){digits = 3}
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
#'@param nbreak number of breaks in the model
#'@param glb global minimum SSR of the model
#'@param datevec estimated break date
#'@param y dependent variables
#'@param z independent variables with regime-specific coefficients
#'@param x independent variables with regime-wise coefficients
#'@return tbl Tables containings:
#'i) global minimum SSR
#'ii) estimated date
#'iii) estimated coefficients
#'


compile.sbtests <- function(x,digits = -1,...)
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

#'Summary output of Sup(l+1|l) test
#'
#'Function to format the output of the Sup F test
#'@param nbreak number of breaks in the model
#'@param glb global minimum SSR of the model
#'@param datevec estimated break date
#'@param y dependent variables
#'@param z independent variables with regime-specific coefficients
#'@param x independent variables with regime-wise coefficients
#'@return tbl Tables containings:
#'i) global minimum SSR
#'ii) estimated date
#'iii) estimated coefficients
#'

#'@export
print.sbtests <- function(x,...)
{ if(x$mbreak == 0){
  cat('\nThe test is undefined for no break model\n')
}else{
  cat('\na) SupF tests against a fixed number of breaks\n\n')
  print(x$supF1,quote=FALSE)
  cat('\nb) UDmax tests against an unknown number of breaks\n\n')
  print(x$UDmax,quote=FALSE)}
  invisible(x)
}


compile.seqtests = function(x){
  if(x$mbreak==1){
    cat('\nThe test is exactly 0 versus 1 break, hence the sequential test is not repeated\n')
    x$sfl = NULL
    return(x)
  }else if (x$mbreak==0){
    cat('\nThe test is undefined for maximum break = 0\n')
    x$sfl = NULL
    return(x)
  }
  else{
  nbreak = x$mbreak-1
  cnames = c()
  for (i in 1:nbreak){
    cnames = cbind(cnames, paste('supF(',i+1,'|',i,')',sep=''))
  }
  x$supfl = format(round(x$supfl,3),nsmall = 3)
  x$ndat = format(x$ndat,nsmall=0)
  x$cv = format(round(x$cv,3),nsmall = 3)
  
  sfl =data.frame(supfl = t(x$supfl[seq(1,nbreak,1),1,drop=FALSE]))
  ndat = t(x$ndat[seq(1,nbreak,1),1,drop=FALSE])
  cv = x$cv[,seq(1,nbreak,1),drop=FALSE]
  colnames(ndat) = colnames(sfl)
  colnames(cv) = colnames(sfl)
  sfl = rbind(sfl,ndat)
  sfl = rbind(sfl,cv)
  colnames(sfl) = cnames
  rownames(sfl) = c('Seq supF','Estimated date','10% CV','5% CV', '2.5% CV', '1% CV')
  x$sfl = sfl
  return (x)}
}

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


#' #' Plot model of estimated n breaks
#' plot.model = function(x,...){
#'   #get data from the model
#'   m = x$nbreak
#'   y = x$y
#'   zreg = x$z
#'   xreg = x$x
#'   p = x$numx
#'   q = x$numz
#'   date = x$date
#'   beta = x$beta
#'   zbar = diag_par(zreg,m,date)
#'   T = length(y)
#'   if (p == 0){
#'     reg = zbar
#'   }
#'   else{
#'     reg = cbind(xreg,zbar)
#'   }
#'   
#'   #compute model with no breaks
#'   fixreg = cbind(xreg,zreg)
#'   fixbeta = OLS(y,fixreg)
#'   fity_fix = fixreg%*%fixbeta
#'   
#'   fity = reg%*%beta
#'   tx = seq(1,T,1)
#'   
#'   range_y = max(y)-min(y);
#'   
#'   #plot original series
#'   graphics::plot(tx,y,type='l',col="black", xlab='time',ylab="y", 
#'        ylim=c(min(y)-range_y/10,max(y)+range_y/10),lty=1)
#'   
#'   #plot fitted values series
#'   
#'   graphics::lines(tx, fity,type='l', col="blue",lty=2)
#'   
#'   
#'   #plot fitted values series
#'   
#'   graphics::lines(tx, fity_fix,type='l', col="dark red",lty=2)
#'   
#'   for (i in 1:m){
#'     graphics::abline(v=date[i,1],lty=2)
#'     if (x$CI[i,1] < 0 || x$CI[i,2]>T){}else{
#'    graphics::segments(x$CI[i,1],min(y)*(12+i/m)/10,x$CI[i,2],min(y)*(12+i/m)/10,lty=1,col='red')}
#'   }
#'   
#'   legend(0,max(y)+range_y/10,legend=c("observed y",paste(m,'break y'),"0 break y"),
#'         lty=c(1,2,2), col=c("black","blue","red"), ncol=1)
#'   
#'   
#' }
#' 
#' 

#' new plot function for class model
#' @importFrom ggplot2 ggplot aes annotate geom_segment geom_line ggtitle .data coord_cartesian scale_color_manual geom_ribbon
#' @export
plot_model = function(model,CI=0.95,title=NULL){
  m = model$nbreak
  if(m==0){
    cat('The model has no break. Visualization for comparison between structural breaks 
        versus no breaks is skipped')
    return(null)
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
  colnames(tbl) = c('time',model$y_name,'ypred_break','ypred_fix')
  
  
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
    cat('Not available CI level, set to 95%')
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
    ggplot2::geom_line(ggplot2::aes(y=.data$ypred_break,colour='y_break'),size=0.3)+
    ggplot2::geom_line(ggplot2::aes(y=.data$ypred_fix,colour='y_fix'),size=0.3)+
    ggplot2::geom_line(ggplot2::aes(y=.data$y,colour='y'),size=0.2)+
    ggplot2::geom_ribbon(ggplot2::aes(ymax=.data$ypred_ub, ymin=.data$ypred_lb), fill="gray", alpha=.35)+
    #ggplot2::geom_ribbon(ggplot2::aes(ymax=.datax.upper, ymin=x.lower), fill="pink", alpha=.5)
    ggplot2::scale_color_manual(name=paste('Legends'),
                                values = c("y_break" = "blue", "y_fix" = "red",'y'='black'),
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


