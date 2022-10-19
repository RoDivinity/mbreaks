## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "")
library(mbreaks)

## ----estimate US rate---------------------------------------------------------
#the data set for the example is real.Rda
data(real)
#carry out all testing and estimating procedures via a single call
rate = mdl('y',data=real)
#display the results from the procedures
rate

## ----graph_USr, fig.width=7---------------------------------------------------
#estimating the mean-shift model with BIC (the default option bic_opt=1, which use BIC as criterion)
rate_BIC = doorder('y',data=real)
#NOTE: equivalent to rate$BIC; type rate$BIC to compare with new result


# visualization of estimated model with BIC (in the argument, we can replace rate$BIC with rate_BIC for exact same graph; recall that `$` is the operator to refer to field BIC in return list from mdl())
plot_model(rate$BIC, title = 'US Exchange rate')

## ----reproduce_graph_USr, fig.width=7-----------------------------------------
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
  

## ----reproduce_table_PY-------------------------------------------------------
data(nkpc)
#x_t is GDP gap
  z_name = c('inflag','ygap','inffut')
  #we can invoke each test separately by using dotest() and doseqtests()
  supF_ygap = dotest('inf',z_name,data=nkpc,prewhit = 0, eps1 = 0.1,m=1)
  #z regressors' names are passed in the argument as an array, which equivalent to above argument call with z_name
  seqF_ygap = doseqtests('inf',c('inflag','ygap','inffut'),data=nkpc,prewhit = 0, eps1=0.1)
  #see test results
  supF_ygap
  seqF_ygap
  
  
#x_t is labor income share 
  #or invoke all tests using mdl() 
  nkpc_lbs = mdl('inf',c('inflag','lbs','inffut'),data=nkpc,prewhit = 0, eps1=0.1, m=5)
  nkpc_lbs$sbtests
  nkpc_lbs$seqtests
  

## ----re_estimate model given known breaks, echo = FALSE-----------------------
#only need to re-estimate model with output gap since we use mdl() for income share, we can obtain the estimated sequential model from SEQ (which is returned from mdl() as a list element)
# It is recommended to store desirable options as variables and set arguments = variables to avoid mistakes and save time
eps1 = 0.1
prewhit = 0
ygap_fixn = dofix('inf',z_name,data=nkpc,fixn=1,prewhit=prewhit,eps1=eps1)
#or use data-dependent sequential approach
ygap_SEQ = dosequa('inf',z_name,data=nkpc,prewhit=prewhit,eps1=eps1)


## ----index_date, include=FALSE------------------------------------------------
ygap_date = ygap_fixn$date

## ----reproduce sub-sample IV estimates, include=FALSE-------------------------
k=4
#list of instruments
instruments = c('inflag','lbslag','ygaplag','spreadlag','dwlag','dcplag')
#list of endogenous
regressors = c('inffut','inflag','lbs')

bigT = dim(nkpc)[1]
#independent variable
Y = as.matrix(nkpc[,'inf',drop=FALSE])

#form matrix of instruments
Z = as.matrix(nkpc[,instruments])
Z = cbind(rep(1,151),Z)

#endogenous variable
X_e = as.matrix(nkpc$inffut,drop=FALSE)

#estimate breaks for first stage with 2 breaks
qz = dim(Z)[2]
mdl_z = dating(X_e,Z,round(bigT*0.1),2,qz,bigT)
bdz = mdl_z$datevec[,2]

Z1   = Z[seq(1,bdz[1],1),];
Z2   = Z[seq(bdz[1]+1,bdz[2],1),];
Z3   = Z[seq(bdz[2]+1,bigT,1),];
X_e1   = X_e[seq(1,bdz[1],1),];
X_e2   = X_e[seq(bdz[1]+1,bdz[2],1),];
X_e3   = X_e[seq(bdz[2]+1,bigT,1),];

Xh1  = Z1%*%solve(t(Z1)%*%Z1)%*%t(Z1)%*%X_e1; 
Xh2  = Z2%*%solve(t(Z2)%*%Z2)%*%t(Z2)%*%X_e2; 
Xh3  = Z3%*%solve(t(Z3)%*%Z3)%*%t(Z3)%*%X_e3; 
Xpred   = rbind(Xh1,Xh2,Xh3)


#2nd stage regressors
X = as.matrix(nkpc[,regressors])
X = cbind(rep(1,151),X)
Xh = as.matrix(nkpc[,c('inflag','lbs')])
Xh = cbind(rep(1,151),Xh,Xpred)

#partition the regressors
T1 = seq(1,ygap_date)
T2 = seq(ygap_date+1,bigT)
Xh1 = Xh[T1,]
Xh2 = Xh[T2,]
Y1 = Y[T1,1,drop=FALSE]
Y2 = Y[T2,1,drop=FALSE]

#full sample estimate:
beta = solve(t(Xh)%*%Xh)%*%t(Xh)%*%Y


#subsample
beta1 = solve(t(Xh1)%*%Xh1)%*%t(Xh1)%*%Y1
beta2 = solve(t(Xh2)%*%Xh2)%*%t(Xh2)%*%Y2

#compute variance
res1 = Y1-Xh1%*%beta1
res2 = Y2-Xh2%*%beta2
#no prewhitening to match paper
hac1 = correct(Xh1,res1,0)
hac2 = correct(Xh2,res2,0)
vhac1 = solve(t(Xh1)%*%Xh1)%*%hac1%*%solve(t(Xh1)%*%Xh1) #4 regressors
vhac1 = vhac1*(125-k)
vhac2 = solve(t(Xh2)%*%Xh2)%*%hac2%*%solve(t(Xh2)%*%Xh2)
vhac2 = vhac2*(bigT-125-k)
stdhac1 = sqrt(diag(vhac1))
stdhac2 = sqrt(diag(vhac2))



## ----display results, echo=FALSE----------------------------------------------
colnames = c('$\\mu$','$\\gamma(\\pi_{t-1})$','$\\kappa(x_t)$','$\\beta(E_t\\pi_{t+1})$')
rownames = c('1960:Q1-1991:Q1','$SE_1$','1991:Q2-1997:Q4','$SE_2$')
IV_estimates = data.frame(round(t(beta1),3))
IV_estimates = rbind(IV_estimates,round(stdhac1,3))
IV_estimates = rbind(IV_estimates,round(t(beta2),3))
IV_estimates = rbind(IV_estimates,round(stdhac2,3))
colnames(IV_estimates) = colnames
rownames(IV_estimates) = rownames
knitr::kable(IV_estimates)

## -----------------------------------------------------------------------------
#create a separate data set to add constant regressor
nkpc_c = nkpc
#create constant regressor
nkpc_c$const = rep(1,151)
# if included a constant explicitly, users need to turn off option const, otherwise matrix of regressors is linearly dependent (2 const regressors) so non-invertible
# use sequential approach to estimate number of breaks for
# i) only expected inflation change
nkpc_inffut = dosequa('inf','inffut',c('const','inflag','ygap'),data=nkpc_c,signif = 1,const=0)

# ii) only expected inflation change
nkpc_ygap = dosequa('inf','ygap',c('const','inflag','inffut'),data=nkpc_c,signif = 1,const=0)

# iii) only expected inflation change
nkpc_inflag = dosequa('inf','inflag',c('const','inffut','ygap'),data=nkpc_c,signif = 1,const=0)

nkpc_inflag

## ---- fig.width=7, echo=FALSE-------------------------------------------------
plot_model(nkpc_inflag,CI=0.9)

