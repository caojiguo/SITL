#================ Simulation II ===============

rm(list = ls())


#============ Similarity-Informed Transfer Learning for Multivariate Functional Censored Quantile Regression =========

library(survival)
library(quantreg)
library(fda)
library(pracma)
library(MASS)
library(xtable)
library(abind)
library(purrr)
library(glmnet)

#============ Generator functions ============

generator_fun_part1 <- function(t, v, nk,S1){
  
  kk = pi*(seq(1, nk, 1) -1) 
  phik = sqrt(2)*cos(t%*%t(kk))  
  phik[, 1] = rep(1, length(t))  
  if(S1==0){
    uk = runif(nk, -sqrt(3), sqrt(3))  
  }else{
    uk = rnorm(nk, 0, 1)  
  }
  xik = ((-1)^(seq(2, (nk+1), 1)))*(seq(1, nk, 1)^(-v/2))
  ZK = xik*uk*t(phik)  
  Z = abs(colSums(ZK))  
  
  psik = 4*((-1)^(seq(1, nk, 1)))*(seq(1, nk, 1)^(-2))  
  psi1 = colSums(psik*t(sqrt(2)*cos(t%*%t(kk))))  
  fun1 = Z*psi1
  
  # for source cohort
  psi2 = exp(t)
  fun2 = Z*psi2 
  
  return(list(Z = Z, fun1 = fun1, psi1 = psi1,fun2=fun2,psi2=psi2)) 
}

generator_fun_part2 <- function(t, nk, sigmaxi){
  
  basis = create.bspline.basis(c(0,1),nbasis = nk, norder = 4)
  Basist = eval.basis(t, basis)
  xi = rnorm(nk,0,sigmaxi)  
  Z = Basist%*%xi 
  
  psi1 = 4*cos(3*pi*t)+4*sin(3*pi*t) 
  fun1 = Z*psi1 
  
  psi2 = 10*exp(t)
  fun2 = Z*psi2 
  
  psi3 = 12*cos(3*pi*t)+12*sin(3*pi*t)
  fun3 = Z*psi3 
  
  return(list(Z = Z, fun1 = fun1, psi1 = psi1,fun2=fun2,psi2=psi2,fun3=fun3,psi3=psi3)) 
}

generator_dataM <- function(t, v, n, sdx1,sdx2,rhox, sde, sxi, b0, cenrate, seed,S1,S2,S3,S4){
  # set.seed(seed) 
  
  Z1 =  matrix(NA, n, m) 
  fun11 = matrix(NA, n, m) 
  fun12 = matrix(NA, n, m)
  
  Z2 =  matrix(NA, n, m) 
  fun21 = matrix(NA, n, m) 
  fun22 = matrix(NA, n, m)
  fun23 = matrix(NA, n, m)
  
  Timedata = matrix(NA, n, 2)  # (survival time, censor time)
  
  # error term
  epsilon = rnorm(n, 0, sde)
  
  # functional covariates part
  for(i in 1:n){
    # set.seed(seed * 100 + i) 
    res_fun1 = generator_fun_part1(t[,1], v, nk,S1) 
    Z1[i, ] = res_fun1$Z 
    fun11[i, ] = res_fun1$fun1
    psi11 = res_fun1$psi1
    fun12[i, ] = res_fun1$fun2
    psi12 = res_fun1$psi2
    
    res_fun2 = generator_fun_part2(t[,2],nk,sxi) 
    Z2[i, ] = res_fun2$Z 
    fun21[i, ] = res_fun2$fun1
    psi21 = res_fun2$psi1
    fun22[i, ] = res_fun2$fun2
    psi22 = res_fun2$psi2
    fun23[i, ] = res_fun2$fun3
    psi23 = res_fun2$psi3
    
  }
  
  iz1 = apply(Z1, 1, sum)*deltat  # int Z_1i(s) 
  iza11 = apply(fun11, 1, sum)*deltat  # int Z_1i(s)psi_1(s)
  iza12 = apply(fun12, 1, sum)*deltat  # int Z_1i(s)psi_2(s)
  
  iz2 = apply(Z2, 1, sum)*deltat  # int Z_2i(s) 
  iza21 = apply(fun21, 1, sum)*deltat  # int Z_2i(s)psi_2(s)
  iza22 = apply(fun22, 1, sum)*deltat  # int Z_2i(s)psi_2(s)
  iza23 = apply(fun23, 1, sum)*deltat  # int Z_2i(s)psi_2(s)
  
  
  if(S2==0){
    mT0 = b0[1,1] + iza11 + iza21 + iz1*epsilon 
  }
  if(S2==1){
    mT0 = b0[1,1] + iza11 + iza21 + iza22 +  iz1*epsilon 
  }
  if(S2==2){
    mT0 = b0[1,1] + iza11 + iza23 +  iz1*epsilon 
  }
  if(S2==3){
    mT0 = b0[1,1] + 3*iza11 + iza21 + iza22 +  iz1*epsilon 
  }
  
  C = exp(mT0)
  C[which(mT0>quantile(mT0,1-cenrate))] = exp(quantile(mT0,1-cenrate))
  
  Timedata = cbind(exp(mT0), C)
  
  
  Z = array(c(Z1,Z2),dim = c(n,dim(tgrid)[1],2))
  # return: Y, X 
  return(list(Timedata = Timedata,  Z = Z, psi11=psi11, psi12=psi12,psi2=psi21,err = epsilon))
}

#===== Estimate function===========
FCLR <- function(surtime,dataZ, tgrid,taugrd,intercept,nbasis,norder){
  #---------------------------------------------------------#
  # event time and censor indicator
  # surtime = cbind(pmin(dataY[, 1], dataY[, 2]), (dataY[, 1] <= dataY[, 2]))
  
  q = dim(dataZ)[3]
  Basismat = list()
  ZBMat = list()
  for (l in 1:q){
    basis = create.bspline.basis(range(tgrid[,l]),nbasis = nbasis, norder = norder)
    Basismat[[l]] = eval.basis(tgrid[,l], basis) 
    ZBMat[[l]] = dataZ[,,l]%*%Basismat[[l]]*deltat[l]
  }
  
  # design matrix
  covdata = do.call(cbind,ZBMat)
  
  T = log(surtime[,1])
  
  # estimate the coefficients
  if(intercept==0){
    fit <- crq(survival::Surv(log(surtime[, 1]), surtime[,2])~.-1, 
               data = data.frame(covdata),  taus=taugrd, method= "PengHuang")
    # estimate the coefficients for each quantile level
    coefest = t(coef(fit,taugrd))
    gammahat1 = t(coefest[,1:nbasis])
    gammahat2 = t(coefest[,(1+nbasis):(2*nbasis)])
    
    # functional coefficient (each one is length(tgrid)*length(tau))
    alphahat = list()
    alphahat[[1]] = Basismat[[1]]%*%gammahat1
    alphahat[[2]] = Basismat[[2]]%*%gammahat2
  }
  
  return(list(covdata = covdata,
              ZBMat = do.call(cbind,ZBMat),
              alphahat = alphahat,
              gammahat = rbind(gammahat1,gammahat2),
              coefest = t(coefest))) 
}

FCLR_boot <- function(res,surtime,dataZ, tgrid,taugrd,intercept,nbasis,norder,times=100){
  #---------------------------------------------------------#
  # event time and censor indicator
  # surtime = cbind(pmin(dataY[, 1], dataY[, 2]), (dataY[, 1] <= dataY[, 2]))
  
  n = length(tgrid[,1])
  q = dim(dataZ)[3]
  Basismat = list() 
  ZBMat = list() 
  for (l in 1:q){
    basis = create.bspline.basis(range(tgrid[,l]),nbasis = nbasis, norder = norder)
    Basismat[[l]] = eval.basis(tgrid[,l], basis) 
    ZBMat[[l]] = dataZ[,,l]%*%Basismat[[l]]*deltat[l]
  }
  m = dim(ZBMat[[1]])[2] + dim(ZBMat[[2]])[2]
  # design matrix
  covdata = do.call(cbind,ZBMat)
  
  G_boot = array(NA,dim=c(times,m,length(taugrd)))
  A1_boot = array(NA,dim=c(times,n,length(taugrd)))
  A2_boot = array(NA,dim=c(times,n,length(taugrd)))
  for(boot in 1:times){
    zeta_boot = rexp(n,1)
    T = log(surtime[,1]) + zeta_boot * rowMeans(res)
    
    # estimate the coefficients
    if(intercept==0){
      fit <- crq(survival::Surv(T, surtime[,2])~.-1, 
                 data = data.frame(covdata),  taus=taugrd, method= "PengHuang")
      coefest = t(coef(fit,taugrd))
      
      gammahat1 = t(coefest[,1:nbasis])
      gammahat2 = t(coefest[,(1+nbasis):(2*nbasis)])
      
      # functional coefficient (each one is length(tgrid)*length(tau))
      alphahat = list()
      alphahat[[1]] = Basismat[[1]]%*%gammahat1
      alphahat[[2]] = Basismat[[2]]%*%gammahat2
    }
    G_boot[boot,,] = rbind(gammahat1,gammahat2)
    A1_boot[boot,,] = Basismat[[1]]%*%gammahat1
    A2_boot[boot,,] = Basismat[[2]]%*%gammahat2
  }
  
  return(list(G_boot = G_boot,
              A1_boot = A1_boot,
              A2_boot = A2_boot,
              BMat = Basismat
  ))
}

#====== source weight =======
#==== Loss function ====
Loss_mFCQ <- function(dataY, dataX, delta, coefhat, taugrd){
  # dataY: min(T,C)
  # dataX: scalar+bspline (do not include intercept)
  # delta: censoring indicator
  # coefhat: length(tau)*(1+nbasis) or length(tau)*(nbasis)
  # taugrd: given quantile levels
  
  # the dimension of dataX and coefhat should be matched
  # dataX = cbind(rep(1,n),dataX) # if the dataX don't include the intercept
  
  r = dim(coefhat)[1]
  taugrd0 = c(0,taugrd)
  coefhat =cbind(rep(0,r),coefhat)
  dh = -log(1 - taugrd0[-1]) + log(1 - taugrd0[-(length(taugrd0))])
  Loss = rep(0,length(taugrd))
  for (j in 2:length(taugrd0)){
    # calculate u_j
    u = (log(dataY)>= dataX%*%coefhat[,1:(j-1)])%*%dh[1:j-1]
    e = log(dataY)- dataX%*%coefhat[,j]
    Loss[j-1] = t(e) %*% (u - (e <= 0 & delta == 1))
  }
  Loss_Total = sum(Loss)
  return(Loss_Total)
}

#==== Loss function2 ====
Loss2_mFCQ <- function(dataY, delta, alpha1_res, alpha2_res, taugrd){
  # dataY: min(T,C)
  # dataX: scalar+bspline (do not include intercept)
  # delta: censoring indicator
  # coefhat: length(tau)*(1+nbasis) or length(tau)*(nbasis)
  # taugrd: given quantile levels
  
  # the dimension of dataX and coefhat should be matched
  # dataX = cbind(rep(1,n),dataX) # if the dataX don't include the intercept
  
  r2 = dim(alpha1_res)[1]
  taugrd0 = c(0,taugrd)
  alpha1_res = cbind(rep(0,r2),alpha1_res)
  alpha2_res = cbind(rep(0,r2),alpha2_res)
  dh = -log(1 - taugrd0[-1]) + log(1 - taugrd0[-(length(taugrd0))])
  Loss = rep(0,length(taugrd))
  for (j in 2:length(taugrd0)){
    # calculate u_j
    u = as.matrix((log(dataY)>= (alpha1_res[,1:(j-1)] +alpha2_res[,1:(j-1)])))%*%dh[1:j-1]
    e = log(dataY) - alpha1_res[,j] - alpha2_res[,j]
    Loss[j-1] = t(e) %*% (u - (e <= 0 & delta == 1))
  }
  Loss_Total = sum(Loss)
  return(Loss_Total)
}

#==== Loss function Bootstrap ====
Loss_Boot <- function(dataY, delta, alpha1_res, alpha2_res, taugrd){
  # dataY: min(T,C)
  # dataX: scalar+bspline (do not include intercept)
  # delta: censoring indicator
  # coefhat: length(tau)*(1+nbasis) or length(tau)*(nbasis)
  # taugrd: given quantile levels
  
  # the dimension of dataX and coefhat should be matched
  # dataX = cbind(rep(1,n),dataX) # if the dataX don't include the intercept
  n = length(dataY)
  r2 = dim(alpha1_res)[1]
  taugrd0 = c(0,taugrd)
  
  alpha1_res = cbind(rep(0,r2),alpha1_res)
  alpha2_res = cbind(rep(0,r2),alpha2_res)
  dh = -log(1 - taugrd0[-1]) + log(1 - taugrd0[-(length(taugrd0))])
  Loss = rep(0,length(taugrd))
  for (j in 2:length(taugrd0)){
    # calculate u_j
    zeta = rexp(n,rate=1)
    u = as.matrix((log(dataY)>= (alpha1_res[,1:(j-1)] +alpha2_res[,1:(j-1)])))%*%dh[1:j-1]
    e = log(dataY) - alpha1_res[,j] - alpha2_res[,j]
    Loss[j-1] = zeta * t(e) %*% (u - (e <= 0 & delta == 1))
  }
  Loss_Total = sum(Loss)
  return(Loss_Total)
}

#==== SITL-threshold ====
Weight_Gauss <- function(dataY, dataX, delta, coefhat_Tar, coefhat_Sour, taugrd,h){
  # dataY: min(T,C) from the target cohort
  # dataX: scalar+bspline (do not include intercept), from the target cohort
  # delta: censoring indicator, from the target cohort
  # coefhat: length(tau)*(1+nbasis) or length(tau)*(nbasis)
  # taugrd: given quantile levels
  n0 = dim(dataX)[1]
  Loss_Tar = Loss_mFCQ(dataY, dataX, delta, coefhat_Tar, taugrd)
  Loss_Sour = Loss_mFCQ(dataY, dataX, delta, coefhat_Sour, taugrd)
  weight = (1/h)*dnorm((Loss_Tar-Loss_Sour)/(n0*h))
  return(weight)
}

#==== Hard-threshold ====
Weight_Unif <- function(dataY, dataX, delta, coefhat_Tar,  coefhat_Sour, taugrd,h,c0){
  # dataY: min(T,C) from the target cohort
  # dataX: scalar+bspline (do not include intercept), from the target cohort
  # delta: censoring indicator, from the target cohort
  # coefhat: length(tau)*(1+nbasis) or length(tau)*(nbasis)
  # taugrd: given quantile levels
  
  Loss_Tar = Loss_mFCQ(dataY, dataX, delta, coefhat_Tar, taugrd)
  Loss_Sour = Loss_mFCQ(dataY, dataX, delta, coefhat_Sour, taugrd)
  weight = (abs(Loss_Tar-Loss_Sour)<h*max(Loss_Tar,c0))*1
  return(weight)
}


#===== Simulation I =============
#parameters for generating data
v= 2
nk = 20
sdx1 = 2
sdx2 = 3
rhox = 1
sde = sqrt(0.2) 
sxi = 2
b0=rbind(0,0)

n0 = 100 # sample size for target cohort
cenrate = 0.3 # censoring rate
c0 = 6 # constant in hard-threshold algorithm
m = 100 # the number of the grid of functional covariate
deltat = c(1/m,1/m)
tgrid = matrix((1:100)/100,100,2)
taugrd = seq(0.1, 0.9, 0.1/10)
# eleq function find the position of specified elements.(which function doesn't  work)
eleq <- Vectorize(function(x, y) {isTRUE(all.equal(x, y))})

# select quantile level 0.3, 0.4, 0.5, 0.6, 0.7
qv <- c(which(eleq(0.3, taugrd)), which(eleq(0.4, taugrd)), which(eleq(0.5, taugrd)), 
        which(eleq(0.6, taugrd)), which(eleq(0.7, taugrd)))

t = tgrid
K = 4
n1 = c(500,1000,500,1000)
S2 = c(0,1,2,3)
S4 = c(0,0,0,0)
intercept = 0
h = 2*log(n0*5)
nbasis = c(6,4)
norder = c(3,3)


#===== simulation function =====
debias2_loss_boot <- function(eta,res,dataX,delta,u,zeta){
  e = res - dataX%*%eta
  gamma_eta = eta[1:m]
  Rho = t(zeta*e) %*% (u - (e <= 0 & delta == 1))
  return(Rho)
}

debias2_gradient_boot <- function(eta,res,dataX,delta,u,zeta){
  e = res - dataX%*%eta
  beta_eta = eta[1:m]
  rho = t(zeta*dataX) %*% (u - (e <= 0 & delta == 1))
  return(rho)
}

coordinate_descent2 <- function(Y, W,taugrd, Tar_TOnly,surtime, converge_tol=1e-8, max_iter=20, G0=NULL){
  # function to conduct coordinate descent and optimize E and H (and return A_est and H_est)
  n <- dim(Y)[1]
  m <- dim(W)[2]
  if(is.null(G0)){
    G0 <- matrix(0,m,length(taugrd)) # initial values of gamma
  }
  
  G <- matrix(0,m,length(taugrd))
  for(i in 1:length(taugrd)){
    fitG <- crq(survival::Surv((Y[,i]), surtime[,2])~.-1, 
                data = data.frame(W),  taus=taugrd, method= "PengHuang")
    G[,i] = coef(fitG,taugrd[i])
  }
  result <- list(G_est=G)
}

coordinate_descent_boot <- function(Y, W,taugrd, Tar_TOnly,surtime,Sour_gamma_weight, converge_tol=1e-6, max_iter=20, G0=NULL,times=100){
  # function to conduct coordinate descent and optimize E and H (and return A_est and H_est)
  n <- dim(Y)[1]
  m <- dim(W)[2]
  if(is.null(G0)){
    G0 <- matrix(0,m,length(taugrd)) # initial values of gamma
  }
  
  G_boot = array(NA,dim=c(m,length(taugrd),times))
  # optimization starts
  for(boot in 1:times){
    G <- matrix(0,m,length(taugrd))
    zeta_boot = rexp(n,1)
    Ytilde = cbind(rep(0,n0),Y)
    taugrd0 = c(0,taugrd)
    G_eta = matrix(0,dim(W)[2],length(taugrd0))
    dh = -log(1 - taugrd0[-1]) + log(1 - taugrd0[-(length(taugrd0))])
    
    for (j in 2:length(taugrd0)){
      # calculate u_j
      u = (Ytilde[,1:(j-1)]>= W%*%G_eta[,1:(j-1)])%*%dh[1:j-1]
      result = optim(0.8*(rbind(Tar_TOnly$gammahat)[,j-1]-rbind(Sour_gamma_weight)[,j-1])[1:8],debias2_loss_boot,debias2_gradient_boot,res=Ytilde[,j],
                     dataX=W,
                     delta=surtime[,2],
                     zeta = zeta_boot,
                     u=u,control = list(maxit = 50000))
      G_eta[,j] = result$par
    }
    G =G_eta[1:m,-1] # p*length(tau)
    G_boot[,,boot] = G
  }
  
  opt_result <- list(G_est=G_boot)
}

Simu_fun2 <- function(tgrid,v,n0,n1,K,sdx1,sdx2,rhox,sde,sxi,b0,cenrate,seed,nbasis,norder,taugrd,intercept,S2, S4,h){
  
  # generate target cohort
  Tar_cohort = generator_dataM(tgrid, v, n0, sdx1,sdx2,rhox, sde,sxi, b0, cenrate, seed,0,0,0,0)
  Tar_dataY = Tar_cohort$Timedata
  Tar_dataZ = Tar_cohort$Z
  Tar_surtime = cbind(pmin(Tar_dataY[, 1], Tar_dataY[, 2]), (Tar_dataY[, 1] <= Tar_dataY[, 2]))
  q = dim(Tar_dataZ)[3]
  
  alpha1true = matrix(Tar_cohort$psi11,length(tgrid[,1]),length(taugrd)) + t(matrix(quantile(Tar_cohort$err,taugrd),length(taugrd),length(tgrid[,1])))  
  alpha2true = matrix(Tar_cohort$psi2,length(tgrid[,2]),length(taugrd))
  
  # result using mFCQR-Tonly
  Tar_TOnly = FCLR(Tar_surtime, Tar_dataZ, tgrid,taugrd,intercept,nbasis[1],norder[1])
  mFCQR_Tonly_alpha1_RMSE = apply((alpha1true-Tar_TOnly$alphahat[[1]])^2, 2, FUN= function (x) sqrt(mean(x)))
  mFCQR_Tonly_alpha2_RMSE = apply((alpha2true-Tar_TOnly$alphahat[[2]])^2, 2, FUN= function (x) sqrt(mean(x)))
  
  #==== Source Cohorts ===============
  # generate source cohort and estimation result 
  Sour_cohort_comb = list()
  Sour_cohort_hard = list()
  iter = 1
  Sour_est = list()
  loss_weight = rep(0,K)
  hard_weight = rep(0,K)
  
  for (k in 1:K){
    Sour_cohort_k = generator_dataM(tgrid, v, n1[k], sdx1,sdx2,rhox, sde, sxi, b0, cenrate, seed*k,0,S2[k],0,S4[k])
    Sour_k_dataY = Sour_cohort_k$Timedata
    Sour_k_dataZ = Sour_cohort_k$Z
    Sour_k_surtime = cbind(pmin(Sour_k_dataY[, 1], Sour_k_dataY[, 2]), (Sour_k_dataY[, 1] <= Sour_k_dataY[, 2]))
    Sour_cohort_comb[[k]] = Sour_cohort_k
    
    Sour_k_est = FCLR(Sour_k_surtime, Sour_k_dataZ, tgrid,taugrd,intercept,nbasis[1],norder[1])
    Sour_est[[k]] = Sour_k_est
    
    # weight of source cohort
    loss_weight[k] = Weight_Gauss(pmin(Tar_dataY[, 1], Tar_dataY[, 2]), Tar_TOnly$covdata, (Tar_dataY[, 1] <= Tar_dataY[, 2]), Tar_TOnly$coefest, Sour_k_est$coefest, taugrd,h)
    hard_weight[k] = Weight_Unif(pmin(Tar_dataY[, 1], Tar_dataY[, 2]), Tar_TOnly$covdata, (Tar_dataY[, 1] <= Tar_dataY[, 2]), Tar_TOnly$coefest, Sour_k_est$coefest, taugrd,h,c0)
    if(hard_weight[k]==1){
      Sour_cohort_hard[[iter]] = Sour_cohort_k
      iter = iter + 1
    }
  } 
  weight = n1*loss_weight/sum(n1*loss_weight)
  
  Sour_alpha1_weight = Sour_est[[1]]$alphahat[[1]]*weight[1] + Sour_est[[2]]$alphahat[[1]]*weight[2] + Sour_est[[3]]$alphahat[[1]]*weight[3]  + Sour_est[[4]]$alphahat[[1]]*weight[4] 
  Sour_alpha2_weight = Sour_est[[1]]$alphahat[[2]]*weight[1] + Sour_est[[2]]$alphahat[[2]]*weight[2] + Sour_est[[3]]$alphahat[[2]]*weight[3] + Sour_est[[4]]$alphahat[[2]]*weight[4]
  Sour_coef_weight = Sour_est[[1]]$coefest*weight[1] + Sour_est[[2]]$coefest*weight[2] + Sour_est[[3]]$coefest*weight[3]  + Sour_est[[4]]$coefest*weight[4] 
  Sour_gamma_weight = Sour_est[[1]]$gammahat*weight[1] + Sour_est[[2]]$gammahat*weight[2] + Sour_est[[3]]$gammahat*weight[3] + Sour_est[[4]]$gammahat*weight[4] 
  
  # data combining target and cohort
  Sour_cohort_comb_re = do.call(Map, c(f = abind,along=1, Sour_cohort_comb))
  Sour_comb_dataY = Sour_cohort_comb_re[[1]]
  Sour_comb_dataZ = Sour_cohort_comb_re[[2]]
  
  Sour_cohort_hard_re = do.call(Map, c(f = abind,along=1, Sour_cohort_hard))
  Sour_hard_dataY = Sour_cohort_hard_re[[1]]
  Sour_hard_dataZ = Sour_cohort_hard_re[[2]]
  
  Comb_dataY = rbind(Tar_dataY,Sour_comb_dataY)
  Comb_dataZ = abind(Tar_dataZ,Sour_comb_dataZ,along = 1)
  Comb_surtime = cbind(pmin(Comb_dataY[, 1], Comb_dataY[, 2]), (Comb_dataY[, 1] <= Comb_dataY[, 2]))
  
  # result using combined data
  mFCQR_Comb = FCLR(Comb_surtime, Comb_dataZ, tgrid, taugrd,intercept,nbasis[1],norder[1])
  
  # debias step
  res = cbind(log(matrix(Tar_surtime[,1],n0,length(taugrd)))-Tar_TOnly$covdata%*%Sour_coef_weight)
  basis2 = create.bspline.basis(range(tgrid[,1]),nbasis = nbasis[2], norder = norder[2])
  Basismat2 = eval.basis(tgrid[,1], basis2) 
  ZBMat1 = Tar_dataZ[,,1]%*%Basismat2*deltat[1]
  ZBMat2 = Tar_dataZ[,,2]%*%Basismat2*deltat[1]
  
  # coordinate descent
  eta = coordinate_descent2(res,cbind(ZBMat1,ZBMat2),taugrd,Tar_TOnly,Tar_surtime)
  eta_gamma = eta$G_est
  
  Trans_alpha1 = Sour_alpha1_weight+Basismat2 %*% eta_gamma[1:nbasis[2],]
  Trans_alpha2 = Sour_alpha2_weight+Basismat2 %*% eta_gamma[(1+nbasis[2]):(2*nbasis[2]),]
  
  # Bootstrap
  # Trans
  eta_boot = coordinate_descent_boot(res,cbind(ZBMat1,ZBMat2),taugrd,Tar_TOnly,Tar_surtime,Sour_gamma_weight,G0=NULL)
  eta_gamma_boot = eta_boot$G_est
  A1_boot = array(NA,dim=c(length(tgrid[,1]),length(taugrd),100))
  A2_boot = array(NA,dim=c(length(tgrid[,1]),length(taugrd),100))
  for(boot in 1:100){
    A1_boot[,,boot] = Basismat2 %*% eta_gamma_boot[1:nbasis[2],,boot]
    A2_boot[,,boot] = Basismat2 %*% eta_gamma_boot[(1+nbasis[2]):(2*nbasis[2]),,boot]
  }
  
  # Comb & Tonly
  Comb_res = log(matrix(rep(Comb_surtime[,1],length(taugrd)),dim(Comb_dataY)[1],length(taugrd)))  - mFCQR_Comb$ZBMat %*% mFCQR_Comb$coefest
  Boot_Comb = FCLR_boot(Comb_res,Comb_surtime,Comb_dataZ, tgrid, taugrd,intercept,nbasis[1],norder[1])
  
  # Confidence Interval
  sgm_a1_A = matrix(NA,length(tgrid[,1]),length(taugrd))
  sgm_a2_A = matrix(NA,length(tgrid[,2]),length(taugrd))
  sgm_a1_C = matrix(NA,length(tgrid[,1]),length(taugrd))
  sgm_a2_C = matrix(NA,length(tgrid[,2]),length(taugrd))
  for(t in 1:length(taugrd)){
    sgm_a1_A[,t] = sqrt(apply(A1_boot[,,],c(1,2),var)[,t])
    sgm_a2_A[,t] = sqrt(apply(A2_boot[,,],c(1,2),var)[,t])
    sgm_a1_C[,t] = sqrt(apply(Boot_Comb$A1_boot[,,],c(2,3),var,na.rm=TRUE)[,t])
    sgm_a2_C[,t] = sqrt(apply(Boot_Comb$A2_boot[,,],c(2,3),var,na.rm=TRUE)[,t])
  }
  
  CI_Trans_alpha1_low = Trans_alpha1 - qnorm(0.975) * sgm_a1_A
  CI_Trans_alpha1_up = Trans_alpha1 + qnorm(0.975) * sgm_a1_A
  CI_Trans_alpha2_low = Trans_alpha2 - qnorm(0.975) * sgm_a2_A
  CI_Trans_alpha2_up = Trans_alpha2 + qnorm(0.975) * sgm_a2_A
  CI_Comb_alpha1_low = mFCQR_Comb$alphahat[[1]] - qnorm(0.975) * sgm_a1_C
  CI_Comb_alpha1_up = mFCQR_Comb$alphahat[[1]] + qnorm(0.975) * sgm_a1_C
  CI_Comb_alpha2_low = mFCQR_Comb$alphahat[[2]] - qnorm(0.975) * sgm_a2_C
  CI_Comb_alpha2_up = mFCQR_Comb$alphahat[[2]] + qnorm(0.975) * sgm_a2_C
  
  result = list(Trans_alpha1 = Trans_alpha1,
                Trans_alpha2 = Trans_alpha2,
                Comb = mFCQR_Comb,
                Tonly = Tar_TOnly,
                CI_Trans_alpha1_low=CI_Trans_alpha1_low,
                CI_Trans_alpha1_up=CI_Trans_alpha1_up,
                CI_Trans_alpha2_low=CI_Trans_alpha2_low,
                CI_Trans_alpha2_up=CI_Trans_alpha2_up,
                CI_Comb_alpha1_low=CI_Comb_alpha1_low,
                CI_Comb_alpha1_up=CI_Comb_alpha1_up,
                CI_Comb_alpha2_low=CI_Comb_alpha2_low,
                CI_Comb_alpha2_up=CI_Comb_alpha2_up)
  return(result)
}

Simu_fun_Tonly <- function(tgrid,v,n0,n1,K,sdx1,sdx2,rhox,sde,sxi,b0,cenrate,seed,nbasis,norder,taugrd,intercept,S2, S4,h){
  
  # generate target cohort
  
  Tar_cohort = generator_dataM(tgrid, v, n0, sdx1,sdx2,rhox, sde,sxi, b0, cenrate, seed,0,0,0,0)
  Tar_dataY = Tar_cohort$Timedata
  Tar_dataZ = Tar_cohort$Z
  Tar_surtime = cbind(pmin(Tar_dataY[, 1], Tar_dataY[, 2]), (Tar_dataY[, 1] <= Tar_dataY[, 2]))
  q = dim(Tar_dataZ)[3]
  
  alpha1true = matrix(Tar_cohort$psi11,length(tgrid[,1]),length(taugrd)) + t(matrix(quantile(Tar_cohort$err,taugrd),length(taugrd),length(tgrid[,1])))  
  alpha2true = matrix(Tar_cohort$psi2,length(tgrid[,2]),length(taugrd))
  
  # result using mFCQR-Tonly
  Tar_TOnly = FCLR(Tar_surtime, Tar_dataZ, tgrid,taugrd,intercept,nbasis[1],norder[1])
  mFCQR_Tonly_alpha1_RMSE = apply((alpha1true-Tar_TOnly$alphahat[[1]])^2, 2, FUN= function (x) sqrt(mean(x)))
  mFCQR_Tonly_alpha2_RMSE = apply((alpha2true-Tar_TOnly$alphahat[[2]])^2, 2, FUN= function (x) sqrt(mean(x)))
  
  #==== Source Cohorts ===============
  # Bootstrap
  Tonly_res = log(matrix(rep(Tar_surtime[,1],length(taugrd)),n0,length(taugrd))) - Tar_TOnly$ZBMat %*% Tar_TOnly$coefest
  Boot_TOnly = FCLR_boot(Tonly_res,Tar_surtime, Tar_dataZ, tgrid,taugrd,intercept,nbasis[1],norder[1])
  
  # Confidence Interval
  sgm_a1_T = matrix(NA,length(tgrid[,1]),length(taugrd))
  sgm_a2_T = matrix(NA,length(tgrid[,2]),length(taugrd))
  for(t in 1:length(taugrd)){
    sgm_a1_T[,t] = sqrt(apply(Boot_TOnly$A1_boot[,,],c(2,3),var)[,t])
    sgm_a2_T[,t] = sqrt(apply(Boot_TOnly$A2_boot[,,],c(2,3),var)[,t])
  }
  
  CI_T_alpha1_low = Tar_TOnly$alphahat[[1]] - qnorm(0.975) * sgm_a1_T
  CI_T_alpha1_up = Tar_TOnly$alphahat[[1]] + qnorm(0.975) * sgm_a1_T
  CI_T_alpha2_low = Tar_TOnly$alphahat[[2]] - qnorm(0.975) * sgm_a2_T
  CI_T_alpha2_up = Tar_TOnly$alphahat[[2]] + qnorm(0.975) * sgm_a2_T
  
  result = list(Tonly = Tar_TOnly,
                CI_T_alpha1_low=CI_T_alpha1_low,
                CI_T_alpha1_up=CI_T_alpha1_up,
                CI_T_alpha2_low=CI_T_alpha2_low,
                CI_T_alpha2_up=CI_T_alpha2_up)
  return(result)
}


S = 100
#------------------------------------------
Comb_alpha1 = array(0,dim=c(m,length(taugrd),S))
Comb_alpha1_low = array(0,dim=c(m,length(taugrd),S))
Comb_alpha1_up = array(0,dim=c(m,length(taugrd),S))

Comb_alpha2 = array(0,dim=c(m,length(taugrd),S))
Comb_alpha2_low = array(0,dim=c(m,length(taugrd),S))
Comb_alpha2_up = array(0,dim=c(m,length(taugrd),S))

Tonly_alpha1 = array(0,dim=c(m,length(taugrd),S))
Tonly_alpha1_low = array(0,dim=c(m,length(taugrd),S))
Tonly_alpha1_up = array(0,dim=c(m,length(taugrd),S))

Tonly_alpha2 = array(0,dim=c(m,length(taugrd),S))
Tonly_alpha2_low = array(0,dim=c(m,length(taugrd),S))
Tonly_alpha2_up = array(0,dim=c(m,length(taugrd),S))

Trans_alpha1 = array(0,dim=c(m,length(taugrd),S))
Trans_alpha1_low = array(0,dim=c(m,length(taugrd),S))
Trans_alpha1_up = array(0,dim=c(m,length(taugrd),S))

Trans_alpha2 = array(0,dim=c(m,length(taugrd),S))
Trans_alpha2_low = array(0,dim=c(m,length(taugrd),S))
Trans_alpha2_up = array(0,dim=c(m,length(taugrd),S))

result = Simu_fun2(tgrid,v,n0,n1,K,sdx1,sdx2,rhox,sde,sxi,b0,cenrate,seed,nbasis,norder,taugrd,intercept,S2, S4,h)
for (i in 1:S){
  try({
    result = Simu_fun2(tgrid,v,n0,n1,K,sdx1,sdx2,rhox,sde,sxi,b0,cenrate,seed,nbasis,norder,taugrd,intercept,S2, S4,h)
  }, silent=TRUE)
  
  Comb_alpha1[,,i] = result$Comb$alphahat[[1]]
  Comb_alpha1_low[,,i] = result$CI_Comb_alpha1_low
  Comb_alpha1_up[,,i] = result$CI_Comb_alpha1_up
  
  Comb_alpha2[,,i] = result$Comb$alphahat[[2]]
  Comb_alpha2_low[,,i] = result$CI_Comb_alpha2_low
  Comb_alpha2_up[,,i] = result$CI_Comb_alpha2_up
  
  Trans_alpha1[,,i] = result$Trans_alpha1
  Trans_alpha1_low[,,i] = result$CI_Trans_alpha1_low
  Trans_alpha1_up[,,i] = result$CI_Trans_alpha1_up
  
  Trans_alpha2[,,i] = result$Trans_alpha2
  Trans_alpha2_low[,,i] = result$CI_Trans_alpha2_low
  Trans_alpha2_up[,,i] = result$CI_Trans_alpha2_up
  
  print(i)
}

it1_l=0
it1_u=0
it2_l=0
it2_u=0
for (i in 1:S){
  try({
    result = Simu_fun_Tonly(tgrid,v,n0,n1,K,sdx1,sdx2,rhox,sde,sxi,b0,cenrate,seed,nbasis,norder,taugrd,intercept,S2, S4,h)
  }, silent=TRUE)
  
  Tonly_alpha1[,,i] = result$Tonly$alphahat[[1]]
  if(any(is.na(result$CI_T_alpha1_low))){
    Tonly_alpha1_low[,,i] = Tonly_alpha1_low[,,i-1]
    it1_l = it1_l+1
  } else{
    Tonly_alpha1_low[,,i] = result$CI_T_alpha1_low
  }
  
  if(any(is.na(result$CI_T_alpha1_up))){
    Tonly_alpha1_up[,,i] = Tonly_alpha1_up[,,i-1]
    it1_u = it1_u+1
  }else{
    Tonly_alpha1_up[,,i] = result$CI_T_alpha1_up
  }
  
  Tonly_alpha2[,,i] = result$Tonly$alphahat[[2]]
  if(any(is.na(result$CI_T_alpha2_low))){
    Tonly_alpha2_low[,,i] = Tonly_alpha2_low[,,i-1]
    it2_l = it2_l+1
  }else{
    Tonly_alpha2_low[,,i] = result$CI_T_alpha2_low
  }
  
  if(any(is.na(result$CI_T_alpha2_up))){
    Tonly_alpha2_up[,,i] = Tonly_alpha2_up[,,i-1]
    it2_u = it2_u+1
  }else{
    Tonly_alpha2_up[,,i] = result$CI_T_alpha2_up
  }  
  print(i)
}







