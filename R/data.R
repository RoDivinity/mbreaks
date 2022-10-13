# Functions and sub-functions to carry out the procedure to find global
# minimizers of the model

#' Calculate optimal point for segment partition
#'
#' Function calculates the optimal one break partition for a segment starting
#' from start to last. The possible range of the break is within
#' \eqn{[b_{last},b_{end}]} based on the vector of recursive SSR of the model
#'
#' @param start start date index of the segment
#' @param last end date index of the segment
#' @param b_start first possible break date
#' @param b_end last possible breakdate
#' @param bigT sample period T
#' @return A list containing the following components:
#' \item{ssrmin}{associated SSR of optimal break}
#' \item{dx}{optimal date (global minimizer)}
#' @references
#' @export
parti = function (start,b_start,b_end,last,bigvec,bigT) {
  #initialize array to store betas and value of index for procedure
  dvec = matrix(0L , nrow = bigT, ncol = 1)
  ini = (start-1)*bigT - (start-2)*(start-1)/2 + 1
  j = b_start

  #start the loop on the segment to find optimal break
  while (j <= b_end){
    jj = j - start
    k = j*bigT - (j-1)*j/2+last-j
    dvec[j,1] = bigvec[ini+jj,1] + bigvec[k,1]
    j = j+1
  }
  #get min SSR
  ssrmin = min( dvec[seq(b_start,b_end,1),1] )
  #get the date with lowest SSR
  dx = (b_start - 1) + which.min(dvec[seq(b_start,b_end,1)])
  out = list('ssrmin' = ssrmin, 'dx' = dx)
  return(out)
}

#' Optimal one break partitions with sequential method
#'
#' Function calculates an optimal one break partitions for a segment that
#' starts at date start and ends at date last. It returns the optimal break
#' and the associated SSR. Procedure used with the sequential method when
#' the T*(T+1)/2 vector of SSR is not computed.
#'@param b_start first possible break date
#'@param b_end last possible break date
#'@param last end of segment considered
#'@param vssr minimum SSRs of associated break date
#'@param vssrev recursive SSRs of the model
#'@return A list containing the following components:
#' \item{ssrmin}{associated SSR of optimal break}
#' \item{dx}{optimal date (global minimizer)}
#' @references
#' @export
partione = function(b1,b2,last,vssr,vssrev){
  dvec = matrix(0L,nrow = last, ncol = 1)
  j = b1;
  while (j <= b2) {
    dvec[j,1]=vssr[j,1]+vssrev[last-j,1];
    j=j+1;
  }
  ssrmin=min(dvec[seq(b1,b2),1]);
  dx=(b1-1)+which.min(dvec[seq(b1,b2),1]);
  out = list('ssrmin' = ssrmin, 'dx' = dx)
  return(out)
}


#' Optimal one break partition in partial structural change model
#'
#' Function computes the optimal one break partition in partial structural
#' change model by searching over all possible breaks
#'
#' @param y dependent variables
#' @param z independent variables with coefficients
#' allowed to change across regimes
#' @param x independent variables with constant coefficients across regimes
#' @param h minimal length of segment
#' @param start initial date to search
#' @param last last date to search
#' @return A list containing the following components
#' \item{ssrmin}{associated SSR of optimal break}
#' \item{dx}{optimal date (global minimizer)}
####
onebp = function(y,z,x,h,start,last) {
  ssrind = 999999999999999
  i = matrix(h,ncol=1)

  while(i <= last - start + 1 - h){

    zb = diag_par(z[start:last,,drop=FALSE],1,i)
    y_reg = y[start:last,1]
    x_reg = cbind(x[start:last,],zb)
    bb = OLS(y_reg,x_reg)
    resid = y_reg - x_reg %*% bb
    ssrn = t(resid) %*% resid

    if (ssrn < ssrind){
      ssrind = ssrn
      bdat = i
    }

    i = i + 1
  }
  bd = bdat + start - 1
  out = list('ssrind' = ssrind,'bd' = bd) #convert code use ssrind to ssrmin, bd to dx
  return(out)
}



#'Calculate optimal break points for linear dating procedure
#'
#'Function computes break points that globally minimizes SSR
#'
#'@param y dependent variable
#'@param z matrix of independent variables
#'@param h minimum length of segment
#'@param m maximum number of breaks
#'@param q number of regressors z
#'@param bigT sample period T
#'@return A list containing the following components:
#'\item{glb}{minimum global SSR}
#'\item{datevec}{Vector of dates (optimal minimizers)}
#'\item{bigvec}{Associated SSRs}
#'@export
dating = function(y,z,h,m,q,bigT){
  ####
  #
  # ************
  # glb = min global SSR
  # datevec = vector of optimal dates
  # bigvec =  associated SSRs
  ###

  #initialize arrays to store results
  datevec = matrix(0L, nrow = m, ncol = m)
  optdat = matrix(0L, nrow = bigT, ncol = m)
  optssr = matrix(0L, nrow = bigT, ncol = m)
  dvec = matrix(0L, nrow = bigT, ncol = 1)
  glb = matrix(0L,nrow = m, ncol = 1)
  bigvec = matrix(0L,nrow = bigT*(bigT+1)/2,ncol = 1)

  #calculate all possible SSR and store
  for (i in 1:(bigT-h+1)) {
    vecssr = ssr(i,y,z,h,bigT)
    bigvec[seq((i-1)*bigT+i-(i-1)*i/2, i*bigT - (i-1)*i/2 ,1),1] = vecssr[seq(i,bigT,1),1]
  }

  #base case: 1 break
  if (m == 1) {
    out = parti(1,h,bigT-h,bigT,bigvec,bigT)
    datevec[1,1] = out$dx
    glb[1,1] = out$ssrmin
  }
  #
  #more than 1 break
  else {
    #change the end point from smallest to full sample T, with case m = 1
    for (j1 in seq(2*h,bigT,1)){
      out = parti(1,h,j1-h,j1,bigvec,bigT)
      optssr[j1,1] = out$ssrmin
      optdat[j1,1] = out$dx
    }
    glb[1,1] = optssr[bigT,1]
    datevec[1,1] = optdat[bigT,1]

    #with case m >= 2
    for (ib in 2:m){
      if (ib == m) {
        jlast = bigT
        for (jb in seq(ib*h,jlast-h,1)) {
          dvec[jb,1] = optssr[jb,ib-1] + bigvec[(jb+1)*bigT - jb*(jb+1)/2,1]
        }
        optssr[jlast,ib] = matrix(t(min(dvec[seq(ib*h,jlast-h,1)])))
        optdat[jlast,ib] = ib*h-1 + which.min(dvec[seq(ib*h,jlast-h,1)])
      }

      else {
        for (jlast in seq((ib+1)*h,bigT,1)){
          for (jb in seq(ib*h,jlast-h,1)){
            dvec[jb,1] = optssr[jb,ib-1] + bigvec[jb*bigT - jb*(jb-1)/2 + jlast -jb,1]
          }
          optssr[jlast,ib] = min(dvec[seq(ib*h,jlast-h,1),1])
          optdat[jlast,ib] = ib*h-1 + which.min(dvec[seq(ib*h,jlast-h,1),1])
        }
      }

      datevec[ib,ib] = optdat[bigT,ib]

      for (i in seq(1,ib-1,1)){
        xx = ib-i
        datevec[xx,ib] = optdat[datevec[xx+1,ib],xx]
      }
      glb[ib,1] = optssr[bigT,ib]
    }
  }

  out = list('glb' = glb, 'datevec' = datevec, 'bigvec' = bigvec)
  return (out)
}

#'Non-linear dating procedure
#'
#'Function computes the break dates of a partial structural change model
#'with pre-specified number of breaks m. The procedure is iterations of
#'recursive computations of SSR and finding the m optimal global minimizers
#'
#'
#'@param y time-series observations of dependent variable
#'@param z matrix of regressors which coefficients are allowed to change across regimes
#'@param x matrix of regressors which coefficients are constant across regime
#'@param p number of \code{z} regressors
#'@param q number of \code{x} regressors
#'@param bigT the sample size T
#'@param fixb option to use initial \eqn{\beta} If \code{1}, procedure requires \code{betaini}.
#'If \code{0}, procedure will not use initial beta values
#'@param betaini initial beta values. Required when use with option \code{fixb}
#'@param eps Convergence criterion
#'@param printd option to print output from iterated estimations. If \code{1}, the results
#'for each iteration will be printed in console log. If \code{0}, no output will be printed
#'@return A list containing the following components:
#'\item{glb}{minimum global SSR}
#'\item{datevec}{Vector of dates (optimal minimizers)}
#'\item{bigvec}{Associated SSRs}
#'@references
#'@export
nldat = function(y,z,x,h,m,p,q,bigT,fixb,eps,maxi,betaini,printd){

  #initialize storage
  glb = matrix(0L , nrow = m, ncol = 1)
  globnl = matrix(0L , nrow = m, ncol = 1)
  datevec = matrix(0L, nrow = m, ncol = m)
  datenl = matrix(0L, nrow = m, ncol = m)

  #initialize current max break
  mi = 1
  while(mi <= m){
    if (printd == 1){
      print(paste('Breaks of model',mi))
    }
    if (fixb == 0){
      qq = p+q
      zz = cbind(x,z)

      #initial estimate of the model
      out = dating(y,zz,h,mi,qq,bigT)
      date = out$datevec[1:mi,mi,drop=FALSE]
      bigvec = out$bigvec

      #partition regressors with initial estimated date
      xbar = diag_par(x,mi,date)
      zbar = diag_par(z,mi,date)

      #calculate initial values of estimate
      teta = OLS(y,cbind(zbar,xbar))
      delta1 = teta[seq(1,q*(mi+1),1),1,drop=FALSE]
      beta1 = OLS(y - zbar %*% delta1, x)

      #calculate initial SSR of the model
      resid1 = y - x%*%beta1 - zbar%*%delta1
      ssr1 = t(resid1) %*% resid1

      if(printd==1) {
        print('The iterations are initialized with')
        prmatrix(delta1)
        prmatrix(beta1)
        print('With break date')
        print(date)}
    }
    else {
      beta1 = betaini
      ssr1 = -5
    }

    #start the iterations

    length = 999999999999
    i = 1

    while (length > eps) {

      out = dating(y-x%*%beta1,z,h,mi,q,bigT)
      #store the date vector for current max
      date = out$datevec[1:mi,mi,drop=FALSE]
      bigvec = out$bigvec
      zbar = diag_par(z,mi,date)

      #update estimates based on new partition
      teta1 = OLS(y,cbind(x,zbar))
      beta1 = teta1[seq(1,p,1),1,drop=FALSE]
      delta1 = teta1[seq(p+1,p+q*(mi+1),1),1,drop=FALSE]


      #check convergence condition
      resid_n = y - cbind(x,zbar) %*% teta1
      ssrn = t(resid_n) %*% resid_n
      length = abs(ssrn - ssr1)

      if(printd==1){
        print(paste('Iteration',i))
      }

      #check upper bound of iterations and update
      if (i >= maxi){
        print('Number of iterations has reached MAX ITER')
      }
      else {
        i = i+1
        ssr1 = ssrn
        glb[mi,1]=ssrn
        datevec[1:mi,mi] = date
      }

    }
    #finished current max breaks & update
    mi = mi + 1

  }


  out = list('glb' = glb, 'datevec' = datevec, 'bigvec' = bigvec)
  return(out)
}
