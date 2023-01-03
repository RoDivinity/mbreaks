## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(mbreaks)

## -----------------------------------------------------------------------------
#load US ex-post exchange rate data 
USrate = data(real)
#carry out all testing and estimating procedures via a single call
rate = mdl('rate',data=real)
#display the results from the procedures
rate

## ---- graph_USr, fig.width=7--------------------------------------------------
#estimating the mean-shift model with BIC (the default option bic_opt=1, which use BIC as criterion)
rate_BIC = doorder('rate',data=real)
#NOTE: equivalent to rate$BIC; type rate$BIC to compare with new result


# visualization of estimated model with BIC (in the argument, we can replace rate$BIC with rate_BIC for exact same graph; recall that `$` is the operator to refer to field BIC in return list from mdl())
plot_model(rate$BIC, title = 'US Exchange rate', CI=0.90)

## ---- reproduce_graph_USr, fig.width=7----------------------------------------
  #collect model information
  m = rate_BIC$nbreak           #number of breaks
  y = rate_BIC$y                #vector of dependent var
  zreg = rate_BIC$z             #matrix of regressors with changing coefs
  date = rate_BIC$date          #estimated date
  fity = rate_BIC$fitted.values #fitted values of model
  bigT = length(y)
  #compute the null model
  fixb = solve((t(zreg) %*% zreg)) %*% t(zreg) %*% y
  fity_fix = zreg%*%fixb    #fitted values of null model
  
  
  #plots the model
  tx = seq(1,bigT,1)
  range_y = max(y)-min(y);
  plot(tx,y,type='l',col="black", xlab='time',ylab="y", 
       ylim=c(min(y)-range_y/10,max(y)+range_y/10),lty=1)
  #plot fitted values series for break model
  lines(tx, fity,type='l', col="blue",lty=2)
  #plot fitted values series for null model
  lines(tx, fity_fix,type='l', col="dark red",lty=2)
  
  #plot estimated dates + CIs
  for (i in 1:m){
    abline(v=date[i,1],lty=2)
    if (rate_BIC$CI[i,1] < 0){rate_BIC$CI[i,1] = 0}
    if(rate_BIC$CI[i,2]>bigT){ rate_BIC$CI[i,2]=bigT}
    segments(rate_BIC$CI[i,1],min(y)*(12+i/m)/10,rate_BIC$CI[i,2],min(y)*(12+i/m)/10,lty=1,col='red')
  }
  
  legend(0,max(y)+range_y/10,legend=c("observed y",paste(m,'break y'),"0 break y"),
        lty=c(1,2,2), col=c("black","blue","red"), ncol=1)

## ---- reproduce_table_PY------------------------------------------------------
#x_t is GDP gap
  z_name = c('inflag','ygap','inffut')
  #we can invoke each test separately by using dotest() and doseqtests()
  supF_ygap = dotest('inf',z_name,data=nkpc,prewhit = 0)
  #z regressors' names are passed in the argument as an array, which equivalent to above argument call with z_name
  seqF_ygap = doseqtests('inf',c('inflag','ygap','inffut'),data=nkpc,prewhit = 0)
  #see test results
  supF_ygap
  seqF_ygap
  
  
#x_t is labor income share 
  #or invoke all tests using mdl() 
  nkpc_lbs = mdl('inf',c('inflag','lbs','inffut'),data=nkpc,prewhit = 0)
  nkpc_lbs$sbtests
  nkpc_lbs$seqtests
  

