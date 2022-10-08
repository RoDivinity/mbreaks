#' Estimation procedure
#'
#' Function to estimate the model by OLS given the obtained break dates
#' It also computes and reports confidence intervals for the break dates and
#' corrected standard errors of coefficients estimates
#' The method used depends on specification for robustness
#'
#'@references
#'
#'@param m number of break
#'@param z matrix of regressor z with coefficients are allowed to change across
#'regimes
#'@param x matrix of regressor x with coefficients are constant across regimes
#'@param q number of regressors z
#'@param p number of regressors x
#'@param b break dates
#'@param hetomega,hetq,hetdat,hetvar options for assumptions on the error terms
#'@param prewhit option to use prewhitening process based on AR(1) approximation
#'@return A list containing the following components:
#'\itemize{
#'\item{date} {List of estimated breaks}
#'\item{CI} {List of Confidence Intervals for each corresponding break}
#'\item{beta} {Estimated coefficients of the regression. The first
#'(\code{m}+1)*\code{q} are coefficients of \code{q} variables \code{z} that change across regimes.
#'The last \code{p} are coefficients of \code{p} variables \code{x}
#'that are constant across regimes}
#'\item{SE} {Corrected standard errors of the corresponding coefficients}}
#'@export

estim = function(m,q,z,y,b,robust,prewhit,hetomega,hetq,x,p,hetdat,hetvar){
  if (m == 0){
    print('There are no breaks in this model and estimation is skipped')
    return (NULL)}
  else{
    bigT = dim(z)[1]
    d = (m+1)*q + p
    vdel = matrix(0L,nrow = d, ncol = d)
    #construct zbar matrix. Diagonal partition of Z
    #at the estimated break date
    zbar = diag_par(z,m,b)

    #estimation and printing
    if (p == 0){
      reg = zbar
    }
    else{
      reg = cbind(x,zbar)
    }

    #estimation of β and δ in pure/partial model
    beta = OLS(y,reg)
    vdel = pvdel(y,z,m,q,bigT,b,prewhit,robust,x,p,1,hetdat,hetvar)
    colnames(beta) = 'coefficients'

    SE = matrix(0L,d,1)
    colnames(SE) = 'corrected SE'
    for (i in 1:d){
      #print(paste('Corrected SE for coefficient',i,'is',sqrt(vdel[i,i])))
      SE[i,1] = sqrt(vdel[i,i])
    }

    if (robust == 0 && hetdat == 1 && hetvar == 0){
      #print('In this case robust=0, hetdat=1 and hetvar=0, the "corrected" are the same as that of the printout except for a different small sample correction.')
      SE[i,1] = sqrt(vdel[i,i])
    }

    #confidence interval for break date
    bound = interval(y,z,zbar,b,q,m,robust,prewhit,hetomega,hetq,x,p)
    CI_95 = bound[,c(1,2)]
    CI_90 = bound[,c(3,4)]

    for (i in 1:m){
      
      #print(paste('The 95% C.I for the',i,'th break is:',bound[i,1],' ',bound[i,2]))
      #print(paste('The 90% C.I for the',i,'th break is:',bound[i,3],' ',bound[i,4]))
    }
    CI = cbind(bound[,1],bound[,2],bound[,3],bound[,4])
    colnames(CI) = c('lower 95% CI','upper 95% CI','lower 90% CI','upper 90% CI')
    rownames(CI) = c(1:m)
    fitted = as.matrix(reg%*%beta)
    resid = as.matrix(y - fitted)
    
    colnames(resid) = 'residuals'
    colnames(fitted) = 'fitted.values'
    SSR = t(resid)%*%resid
    out = list('SE' = SE, 'CI' = CI, 'beta' = beta, 'date' = b, 
               'SSR' = SSR, 'resid' = resid,'fitted.values' = fitted)
    return (out)
  }

}

