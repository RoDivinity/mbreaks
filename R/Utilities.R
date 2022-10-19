# Utilities functions

###### Matrix auxiliary functions (only available in MATLAB) ######
#' Function to calculate power of a matrix
#'  @noRd
mpower = function(A,t){
  ####
  # A = input matrix
  # t = exponent
  # ****
  # out = A^t
  ####
  B = A
  if (t == 1) {return(B)}
  else{for (i in 2:t) {B = B %*% A}}
  return (B)
}


#' Rotate a matrix
#'
#'Function rotate sa matrix 90 degree counterclockwise
#'@author Linh Nguyen
#'@param A matrix to be rotated 90 degree counterclockwise
#'@return B: rotated matrix A
#'@examples
#'A = matrix(c(1:10),5,2)
#'rot90(A)
#' @noRd

rot90 = function(A){
  rA = dim(A)[1]
  cA = dim(A)[2]
  B = matrix (0L, nrow = cA, ncol = rA)
  for (i in 1:rA){
    tempA = A[i,,drop=FALSE]
    tempA = tempA[,seq(dim(A)[2],1,-1)] #reverse the row elements
    tempB = t(tempA) #transpose the row
    B[,i] = tempB
  }
  return(B)
}

#'Kronecker Tensor Product
#'
#'Function computes Kronecker Tensor Product between matrix A and matrix B
#'
#'@author Linh Nguyen
#'@param A Matrix A
#'@param B Matrix B
#'@return out: Kronecker Tensor Product of A (X) B
#'@examples
#'A = matrix(c(1,2),1,2)
#'B = matrix(c(1:4),2,2)
#'kron(A,B)
#'@noRd
kron = function(A,B){
  ###
  # A,B = m x n matrix A & p x q matrix B
  # ****
  # out = A (x) B
  ###
  rA = dim(A)[1]
  cA = dim(A)[2]
  rB = dim(B)[1]
  cB = dim(B)[2]

  out = matrix(0L, nrow = rA * rB, ncol = cA * cB)

  for (i in 1:rA){
    for (j in 1:cA) {
      r_index = seq((i-1)*rB+1,i*rB)
      c_index = seq((j-1)*cB+1,j*cB)
      out[r_index,c_index] = A[i,j] * B
    }
  }
  return(out)
}

#'Diagonal partition given break dates
#'
#' Function constructs the matrix of regressors z which coefficients are changed
#' on the estimated break dates
#'
#'
#'@param input matrix of independent variables z with coefficients allowed to
#'change overtime
#'@param m number of breaks in the series
#'@param date vector of estimated break dates
#'@examples
#' z = matrix(c(1:100),50,2)
#' m = 2  #2 breaks
#' date = matrix(c(15,30),2,1) #first break at t = 15; second break at t = 30
#' diag_par(z,m,date)
#'
#'@return output: matrix of partitioned variables corresponds to break dates
#'@export
###
diag_par = function(input,m,date) {
  date =as.matrix(date)
  nt = dim(input)[1]
  q1 = dim(input)[2]
  #create output matrix of zeros with size: nt x (break+1)*q1
  output = matrix(0L,nrow = nt, ncol = (m+1)*q1)
  #copy values of 1st segment of z to the output matrix
  output[c(1:date[1,1]),c(1:q1)] = input[c(1:date[1,1]),,drop=FALSE]
  i = 2
  while (i <= m){
    #copy values of i-th segment of z to output matrix corresponding to date vector
    r_i = seq(date[i-1,1]+1,date[i,1],1) #indices of rows to copy input values
    c_i = seq((i-1)*q1+1,i*q1,1)     #indices of cols to copy input values
    output[r_i,c_i]=input[r_i,,drop=FALSE]
    i = i+1
  }
  rl_i = seq(date[m,1]+1,nt,1)
  c_i = seq(m*q1+1,(m+1)*q1,1)
  output[rl_i,c_i] = input[rl_i,,drop=FALSE]
  return (output)
}


######Computation auxiliary functions ############
#' OLS regression in matrix form
#'
#' Function computes OLS estimates of the regression y on x
#' in matrix form, \eqn{\beta = (X'X)^{-1} X'Y}
#'
#'@param y matrix of dependent variables
#'@param x matrix of independent variables
#'@return b: OLS estimates of the regression
#'
#'@noRd
OLS = function(y,x) {
  b = solve((t(x) %*% x)) %*% t(x) %*% y
  b = matrix(b)
  return (b)
}


#' SSR computation
#'
#' Function computes Sum of Squared Residuals based on OLS estimates, \eqn{\beta}
#'
#'@param y matrix of dependent variables
#'@param zz matrix of independent variables
#'@return ssr: Sum of Squared Residuals
#'@noRd
nssr = function(y,zz) {
  delta = OLS(y,zz)
  resid = y - zz %*% delta
  ssr = t(resid) %*% resid
  return(ssr)
}

#'Recursive SSR computation of a segment
#'
#'Function computes the recursive residuals of period from 'start' to 'last'
#'
#'
#'@param start starting index of the sample used
#'@param last ending index of the last segment considered
#'@param y dependent vars
#'@param z matrix of time-variant regressors (dimension q)
#'@param h minimal length of segment
#'@return out The vector of recursive SSR of length (last-start+1)
#'@noRd
ssr = function(start,y,z,h,last) {
  out = matrix(0L,nrow = last, ncol = 1)

  #initialize the recursion
  index_0 = seq(start,start+h-1,1)
  z_0 = z[index_0,,drop=FALSE]
  y_0 = y[index_0,,drop=FALSE]
  inv1 = solve(t(z_0) %*% z_0)
  delta1 = OLS(y_0,z_0) #calculate initial beta
  res = y_0 - z_0 %*% delta1 #calculate initial residuals
  out[start+h-1,1] = t(res) %*% res # store initial SSR

  #loop to construct the recursive residuals and update SSR
  r = start+h
  while (r <= last) {
    v = y[r,1] - z[r,]%*%delta1  #residual of next obs
    v = drop(v)
    invz = inv1 %*% matrix(t(z[r,]))
    f = 1 + z[r,,drop=FALSE] %*% invz
    f = drop(f)
    delta2 = delta1 + (invz * v) / f
    inv2 = inv1 - (invz %*% t(invz)) / f
    inv1 = inv2
    delta1 = delta2
    out_t = out[r-1,1] + v*v/f
    out[r,1] = out_t
    r=r+1
  }
  return(out)
}


#'Bandwidth calculation based on AR(1) approximation
#'
#'Function computes automatic bandwidth based on AR(1) approximation of
#'estimated vector of vhat with equal weight 1
#'
#'@param vhat Estimated variance
#'@return st Bandwidth-corrected variance
#'@references
#'@noRd
bandw = function(vhat) {
  nt = dim(vhat)[1]
  d = dim(vhat)[2]

  a2n = 0
  a2d = 0

  for (i in seq(1,d,1)){
    y = vhat[seq(2,nt,1),i,drop=FALSE]
    dy = vhat[seq(1,nt-1,1),i,drop=FALSE]
    b = OLS(y,dy)
    sig = nssr(y,dy) #calculate SSR based on AR(1) approximation
    sig = sig/(nt-1) #in-sample correction
    a2n = a2n + 4 * b * b * sig * sig / (1-b)^8
    a2d = a2d + sig * sig / (1-b)^4
  }

  a2 = a2n/a2d
  st = 1.3221 * (a2 * nt) ^.2
  return (st)
}

#'Evaluate quadratic kernel of x
#'
#'Function computes the quadratic kernel by transforming
#'value of x according to \eqn{k(\delta) = 3\frac{sin(\delta)}{\delta} -
#'\frac{cos(\delta)}{\delta^2}} where \eqn{\delta(x) = \frac{x6\pi}{5}} is
#'the transformation of x.
#'
#'@param x Value to evaluate kernel
#'@return ker Kernel value of x
#'@noRd
kern = function(x) {
  #####
  # x = value of some x
  ###
  # ker = quadratic kernel
  ####
  del=6*pi*x/5;
  ker=3*(sin(del)/del-cos(del))/(del*del);
  return (ker)
}


####### Auxiliary Variance-Covariance matrix constructors ###########
#'Construct diagonal matrix according to break date
#'
#'Function constructs a diagonal matrix of dimension (m+1) by (m+1) matrix
#'with i-th entry \eqn{\frac{T_{i} - T_{i-1}}{ T}}
#'
#'@param b Estimated date of changes
#'@param m Number of breaks
#'@param bigT The sample size T
#'@return lambda (m+1) by (m+1) diagonal matrix
#'with i-th entry \eqn{\frac{T_{i} - T_{i-1}}{ T}}
#

plambda = function(b,m,bigT) {

  lambda = matrix(0L , nrow = m+1, ncol = m+1)
  lambda[1,1] = b[1,1]/bigT

  k = 2
  while (k<= m){
    lambda[k,k] = ( b[k,1] - b[k-1,1] ) / bigT
    k = k+1
  }
  lambda[m+1,m+1] = ( bigT - b[m,1] ) / bigT

  return (lambda)
}

#' Construct diagonal matrix of estimated variance
#'
#'Function computes a diagonal matrix of dimension m+1 by m+1
#'with i-th entry is the estimated variance of residuals of segment i
#'
#'@param res big residual vector of the model
#'@param b Estimated date of changes
#'@param m Number of breaks
#'@param nt The sample size
#'@return sigmat (i+1)x(i+1) diagonal matrix with i-th entry
#'equal to estimated variance of regime i
#'
psigmq = function (res,b,q,m,nt) {

  sigmat = matrix(0L, nrow = m+1, ncol = m+1)
  resid_t = res[seq(1,b[1,1],1),1,drop=FALSE] #check this line if has problem
  sigmat[1,1] = t(resid_t) %*% resid_t / b[1,1]

  kk = 2
  while (kk <= m) {
    #check this loop for problem if errors
    bf = b[kk-1,1]
    bl = b[kk,1]
    resid_temp = res[seq(bf,bl,1),1]
    sigmat[kk,kk] = t(resid_temp) %*% resid_temp / (bl - bf)
    kk = kk+1
  }

  res_f = res[seq(b[m,1]+1,nt,1),1,drop = FALSE]
  sigmat[m+1,m+1] = t(res_f) %*% res_f / (nt - b[m,1])

  return(sigmat)
}

######## Auxiliary Computation Functions for Statistical Inference #######
#'Calculate p-value
#'
#'Function computes the p-value of the test
#'@importFrom stats pnorm
#'@param x
#'@param bet
#'@param alph
#'@param b
#'@param deld
#'@param gam
#'@return g The p-value
funcg = function (x,bet,alph,b,deld,gam){
  #g = pval

  if(x <= 0){
    xb = bet * sqrt(abs(x))
    if (abs(xb) <= 30){
      g = -sqrt(-x/(2*pi)) * exp(x/8) - (bet/alph)*exp(-alph*x)*stats::pnorm(-bet*sqrt(abs(x)))+
        ((2*bet*bet/alph)-2-x/2)*pnorm(-sqrt(abs(x))/2)
    }
    else{
      aa=log(bet/alph)-alph*x-xb^2/2-log(sqrt(2*pi))-log(xb)
      g=-sqrt(-x/(2*pi))*exp(x/8)-exp(aa)*pnorm(-sqrt(abs(x))/2)+
        ((2*bet*bet/alph)-2-x/2)*pnorm(-sqrt(abs(x))/2)
    }
  }
  else{
    xb=deld*sqrt(x)

    if (abs(xb) <= 30){

      g=1+(b/sqrt(2*pi))*sqrt(x)*exp(-b*b*x/8)+
        (b*deld/gam)*exp(gam*x)*pnorm(-deld*sqrt(x))+
        (2-b*b*x/2-2*deld*deld/gam)*pnorm(-b*sqrt(x)/2);
    }
    else{

      aa=log((b*deld/gam))+gam*x-xb^2/2-log(sqrt(2*pi))-log(xb);
      g=1+(b/sqrt(2*pi))*sqrt(x)*exp(-b*b*x/8)+
        exp(aa)+(2-b*b*x/2-2*deld*deld/gam)*pnorm(-b*sqrt(x)/2);
    }
  }
  return(g)
}

#'Critical values computation
#'
#'Function calculates the critical values for break date
#'
#'@param eta
#'@param phi1s
#'@param phi2s
#'@return cvec Critical values of break dates

cvg = function(eta,phi1s,phi2s){
  cvec = matrix(0L,nrow = 4, ncol = 1)
  a=phi1s/phi2s
  gam=((phi2s/phi1s)+1)*eta/2
  b=sqrt(phi1s*eta/phi2s)
  deld=sqrt(phi2s*eta/phi1s)+b/2
  alph=a*(1+a)/2
  bet=(1+2*a)/2
  sig = c(0.025,0.05,0.95,0.975)

  isig = 1
  while(isig <= 4){
    #initialize upper, lower bound and critical value
    upb = 2000
    lwb = -2000
    crit = 999999
    cct = 1

    while(abs(crit) >= 0.000001) {
      cct = cct + 1
      if (cct > 100){
        cat('the procedure to get critical values for the break dates has reached the upper bound on the number of iterations. This may happens in the procedure cvg. The resulting confidence interval for this break date is incorect')
        break
      }
      else{
        xx = lwb + (upb-lwb)/2
        pval=funcg(xx,bet,alph,b,deld,gam)
        crit = pval - sig[isig]
        if (crit <= 0) {
          lwb = xx}
        else {
          upb = xx}
      }
    }
    cvec[isig,1] = xx
    isig = isig+ 1
  }

  return(cvec)
}

#'SupF test Critical Values
#'
#'Function to retrieve critical values of supF test stored in
#'/SysData/SupF/cv_{x}.csv where \code{x} corresponds to the trimming level:
#'\itemize{
#'\item{1: \code{eps1} = 5%}
#'\item{2: \code{eps1} = 10%}
#'\item{3: \code{eps1} = 15%}
#'\item{4: \code{eps1} = 20%}
#'\item{5: \code{eps1} = 25%}}
#'The critical values are tabulated from @references
#'@param signif significant level
#'@param eps1 trimming level
#'@return cv Critical value of SupF test
getcv1 = function(signif,eps1){
  if(eps1 == .05){
    out = supFcv1
    cv = out[seq((signif-1)*10+1,signif*10,1),]
  }

  if(eps1 == .10) {
    out = supFcv2
    cv = out[seq((signif-1)*10+1,signif*10,1),]
  }

  if(eps1 == .15) {
    out = supFcv3
    cv = out[seq((signif-1)*10+1,signif*10,1),]
  }

  if(eps1 == .20) {
    out = supFcv4
    cv = out[seq((signif-1)*10+1,signif*10,1),]
  }

  if(eps1 == .25) {
    out = supFcv5
    cv = out[seq((signif-1)*10+1,signif*10,1),]
  }
  cv = data.matrix(cv)
  colnames(cv) = NULL
  rownames(cv) = NULL
  return(cv)
}

#' SupF(l+1|l) test Critical Values
#'
#'Function to retrieve critical values of SupF(l+1|l) test stored in
#'/SysData/SupF_next/cv_{x}.csv where \code{x} corresponds to the trimming level:
#'\itemize{
#'\item{1: \code{eps1} = 5%}
#'\item{2: \code{eps1} = 10%}
#'\item{3: \code{eps1} = 15%}
#'\item{4: \code{eps1} = 20%}
#'\item{5: \code{eps1} = 25%}}
#'The critical values are tabulated from @references
#'@param signif significant level
#'@param eps1 trimming level
#'@return cv Critical value of SupF(l+1|l) test
getcv2 = function(signif,eps1){
  if(eps1 == .05){
    out = supF_next_cv1
    cv = out[seq((signif-1)*10+1,signif*10,1),]
  }

  if(eps1 == .10) {
    out = supF_next_cv2
    cv = out[seq((signif-1)*10+1,signif*10,1),]
  }

  if(eps1 == .15) {
    out = supF_next_cv3
    cv = out[seq((signif-1)*10+1,signif*10,1),]
  }

  if(eps1 == .20) {
    out = supF_next_cv4
    cv = out[seq((signif-1)*10+1,signif*10,1),]
  }

  if(eps1 == .25) {
    out = supF_next_cv5
    cv = out[seq((signif-1)*10+1,signif*10,1),]
  }
  cv = data.matrix(cv)
  colnames(cv) = NULL
  rownames(cv) = NULL
  return(cv)
}

#'Dmax test Critical Values
#'
#'Function to retrieve critical values of supF test stored in
#'/SysData/Dmax/cv_{x}.csv where \code{x} corresponds to the trimming level:
#'\itemize{
#'\item{1: \code{eps1} = 5%}
#'\item{2: \code{eps1} = 10%}
#'\item{3: \code{eps1} = 15%}
#'\item{4: \code{eps1} = 20%}
#'\item{5: \code{eps1} = 25%}}
#'The critical values are tabulated from @references
#'
#'@param signif significant level
#'@param eps1 trimming level
#'@return cv Critical value of SupF test
#'@noRd
getdmax = function(signif,eps1){

  if(eps1 == .05){
    out = Dmax_cv1
    cv = out[seq((signif-1)*10+1,signif*10,1),]
  }

  if(eps1 == .10) {
    out = Dmax_cv2
    cv = out[seq((signif-1)*10+1,signif*10,1),]
  }

  if(eps1 == .15) {
    out = Dmax_cv3
    cv = out[seq((signif-1)*10+1,signif*10,1),]
  }

  if(eps1 == .20) {
    out = Dmax_cv4
    cv = out[seq((signif-1)*10+1,signif*10,1),]
  }

  if(eps1 == .25) {
    out = Dmax_cv5
    cv = out[seq((signif-1)*10+1,signif*10,1),]
  }
  cv = data.matrix(cv)
  colnames(cv) = NULL
  rownames(cv) = NULL
  return(cv)
}


#' function to process data before invoking procedures
#' @param y_name
#' @param z_name
#' @param x_name
#' @noRd
process_data = function(y_name,z_name = NULL,x_name = NULL,data,const){
  #handle data
  y_ind = match(y_name,colnames(data))
  x_ind = match(x_name,colnames(data))
  z_ind = match(z_name,colnames(data))

  if(is.na(y_ind)){
    stop('No matching dependent variable y. Please try again')
    return(NULL)}
  else{
  y = data[,y_ind]
  y = data.matrix(y)
  T = dim(y)[1]


  if (is.null(x_name)) {x = c()}
  else{
    if(anyNA(x_ind)){stop('No x regressors found. Please try again')
     }
    else{x = data.matrix(data[,x_ind,drop=FALSE])}}
  if (is.null(z_name)) {
    if (const == 1) {z = matrix(1L,T,1)}
    else {stop('No regressors in the model. Please specify at least 1 regressors')}
    }
  else{
    if(anyNA(z_ind)){stop('No z regressors found. Please try again')}
    else{
      if (const == 1){
        z = data.matrix(data[,z_ind,drop=FALSE])
        z = cbind(matrix(1L,T,1),z)}
      else{
        z = data.matrix(data[,z_ind,drop=FALSE])
      }}}
  out = list();
  #store y,x,z name and corresponding x,y and z data
  out$y_name = y_name
  out$z_name = z_name
  out$x_name = x_name
  out$y = y
  out$z = z
  out$x = x
  }

  return(out)
}


