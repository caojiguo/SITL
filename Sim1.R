#================ Simulation I ================

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

#===== Estimate function ===========
# surtime = cbind(pmin(dataY[, 1], dataY[, 2]), (dataY[, 1] <= dataY[, 2]))
FCLR <- function(surtime,dataZ, tgrid,taugrd,intercept,nbasis,norder){
  #---------------------------------------------------------#
  # event time and censor indicator
  # surtime = cbind(pmin(dataY[, 1], dataY[, 2]), (dataY[, 1] <= dataY[, 2]))
  
  
  q = dim(dataZ)[3]
  Basismat = list() # q matrix, each one is length(t_l)*nbasis_l
  ZBMat = list() # q matrix, each one is n*nbasis_l
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


#====== self-defined functions used in debias step ======

debias_loss <- function(eta,res,dataX,delta,u,lambda1){
  e = res - dataX%*%eta
  beta_eta = eta[1:p]
  Rho = t(e) %*% (u - (e <= 0 & delta == 1)) + lambda1 * sum(abs(beta_eta))
  return(Rho)
}

debias_gradient <- function(eta,res,dataX,delta,u,lambda1){
  e = res - dataX%*%eta
  beta_eta = eta[1:p]
  rho = t(dataX) %*% (u - (e <= 0 & delta == 1)) + lambda1 * sum(sign(beta_eta))
  return(rho)
}


#====== source weight =======
#==== Loss function ====
Loss_mFCQ <- function(dataY, dataX, delta, coefhat, taugrd){
  # dataY: min(T,C)
  # dataX: bspline (do not include intercept)
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


#===== Set up=============

#parameters for generating data
v= 2
nk = 20
sdx1 = 2
sdx2 = 3
rhox = 1
sde = sqrt(0.2) 
sxi = 2
b0 = rbind(0,0)


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
coordinate_descent2 <- function(Y, W,taugrd, Tar_TOnly,Sour_beta_weight,surtime, converge_tol=1e-8, max_iter=20, G0=NULL){
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
    hard_weight[k] = Weight_Unif(pmin(Tar_dataY[, 1], Tar_dataY[, 2]), Tar_TOnly$covdata, (Tar_dataY[, 1] <= Tar_dataY[, 2]), Tar_TOnly$coefest, Sour_k_est$coefest, taugrd,12,c0)
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
  
  Hard_dataY = rbind(Tar_dataY,Sour_hard_dataY)
  Hard_dataZ = abind(Tar_dataZ,Sour_hard_dataZ,along = 1)
  Hard_surtime = cbind(pmin(Hard_dataY[, 1], Hard_dataY[, 2]), (Hard_dataY[, 1] <= Hard_dataY[, 2]))
  
  # result using combined data
  mFCQR_Comb = FCLR(Comb_surtime, Comb_dataZ, tgrid, taugrd,intercept,nbasis[1],norder[1])
  mFCQR_Hard = FCLR(Hard_surtime, Hard_dataZ, tgrid, taugrd,intercept,nbasis[1],norder[1])
  
  # debias step
  res = cbind(log(matrix(Tar_surtime[,1],n0,length(taugrd)))-Tar_TOnly$covdata%*%Sour_coef_weight)
  basis2 = create.bspline.basis(range(tgrid[,1]),nbasis = nbasis[2], norder = norder[2])
  Basismat2 = eval.basis(tgrid[,1], basis2) 
  ZBMat1 = Tar_dataZ[,,1]%*%Basismat2*deltat[1]
  ZBMat2 = Tar_dataZ[,,2]%*%Basismat2*deltat[1]
  
  # coordinate descent
  eta = coordinate_descent2(res,cbind(ZBMat1,ZBMat2),taugrd,Tar_TOnly,Sour_beta_weight,Tar_surtime)
  eta_gamma = eta$G_est
  
  Trans_alpha1 = Sour_alpha1_weight+Basismat2 %*% eta_gamma[1:nbasis[2],]
  Trans_alpha2 = Sour_alpha2_weight+Basismat2 %*% eta_gamma[(1+nbasis[2]):(2*nbasis[2]),]
  
  # RMSE of alpha
  mFCQR_Comb_alpha1_RMSE = rbind(apply((alpha1true-mFCQR_Comb$alphahat[[1]])^2, 2, FUN= function (x) sqrt(mean(x))))
  mFCQR_Hard_alpha1_RMSE = rbind(apply((alpha1true-mFCQR_Hard$alphahat[[1]])^2, 2, FUN= function (x) sqrt(mean(x))))
  mFCQR_Sourave_alpha1_RMSE = rbind(apply((alpha1true-Sour_alpha1_weight)^2, 2, FUN= function (x) sqrt(mean(x))))
  Trans_mFCQR_alpha1_RMSE = rbind(apply((alpha1true-Trans_alpha1)^2, 2, FUN= function (x) sqrt(mean(x))))
  
  mFCQR_Comb_alpha2_RMSE = rbind(apply((alpha2true-mFCQR_Comb$alphahat[[2]])^2, 2, FUN= function (x) sqrt(mean(x))))
  mFCQR_Hard_alpha2_RMSE = rbind(apply((alpha2true-mFCQR_Hard$alphahat[[2]])^2, 2, FUN= function (x) sqrt(mean(x))))
  mFCQR_Sourave_alpha2_RMSE = rbind(apply((alpha2true-Sour_alpha2_weight)^2, 2, FUN= function (x) sqrt(mean(x))))
  Trans_mFCQR_alpha2_RMSE = rbind(apply((alpha2true-Trans_alpha2)^2, 2, FUN= function (x) sqrt(mean(x))))
  
  # Residual
  Test_cohort = generator_dataM(tgrid, v, 100, sdx1,sdx2,rhox, sde,sxi, b0, cenrate, seed,0,0,0,0)
  Test_dataY = Test_cohort$Timedata
  Test_dataZ = Test_cohort$Z
  Test_surtime = cbind(pmin(Test_dataY[, 1], Test_dataY[, 2]), (Test_dataY[, 1] <= Test_dataY[, 2]))
  
  basis3 = create.bspline.basis(range(tgrid[,1]),nbasis = nbasis[1], norder = norder[1])
  Basismat3 = eval.basis(tgrid[,1], basis3) 
  ZBMat3 = Test_dataZ[,,1]%*%Basismat3*deltat[1]
  ZBMat4 = Test_dataZ[,,2]%*%Basismat3*deltat[1]
  
  Trans_alpha1_res = matrix(NA,nrow = 100,length(taugrd))
  Trans_alpha2_res = matrix(NA,nrow = 100,length(taugrd))
  TO_alpha1_res = matrix(NA,nrow = 100,length(taugrd))
  TO_alpha2_res = matrix(NA,nrow = 100,length(taugrd))
  Comb_alpha1_res = matrix(NA,nrow = 100,length(taugrd))
  Comb_alpha2_res = matrix(NA,nrow = 100,length(taugrd))
  Hard_alpha1_res = matrix(NA,nrow = 100,length(taugrd))
  Hard_alpha2_res = matrix(NA,nrow = 100,length(taugrd))
  
  for(i in 1:length(taugrd)){
    Trans_alpha1_res[,i] = apply(Test_dataZ[,,1]%*%Trans_alpha1[,i], 1, sum)*deltat
    Trans_alpha2_res[,i] = apply(Test_dataZ[,,2]%*%Trans_alpha2[,i], 1, sum)*deltat
    TO_alpha1_res[,i] = apply(Test_dataZ[,,1]%*%Tar_TOnly$alphahat[[1]][,i], 1, sum)*deltat
    TO_alpha2_res[,i] = apply(Test_dataZ[,,2]%*%Tar_TOnly$alphahat[[2]][,i], 1, sum)*deltat
    Comb_alpha1_res[,i] = apply(Test_dataZ[,,1]%*%mFCQR_Comb$alphahat[[1]][,i], 1, sum)*deltat
    Comb_alpha2_res[,i] = apply(Test_dataZ[,,2]%*%mFCQR_Comb$alphahat[[2]][,i], 1, sum)*deltat
    Hard_alpha1_res[,i] = apply(Test_dataZ[,,1]%*%mFCQR_Hard$alphahat[[1]][,i], 1, sum)*deltat
    Hard_alpha2_res[,i] = apply(Test_dataZ[,,2]%*%mFCQR_Hard$alphahat[[2]][,i], 1, sum)*deltat
    
  }
  
  mFCQR_Comb_loss = Loss2_mFCQ(Test_surtime[,1],Test_surtime[,2],Comb_alpha1_res,Comb_alpha2_res,taugrd)
  mFCQR_Hard_loss = Loss2_mFCQ(Test_surtime[,1],Test_surtime[,2],Hard_alpha1_res,Hard_alpha2_res,taugrd)
  mFCQR_Tonly_loss = Loss2_mFCQ(Test_surtime[,1],Test_surtime[,2],TO_alpha1_res,TO_alpha2_res,taugrd)
  Trans_mFCQR_loss = Loss2_mFCQ(Test_surtime[,1],Test_surtime[,2],Trans_alpha1_res,Trans_alpha2_res,taugrd)
  
  result = list(alpha1_RMSE=rbind(mFCQR_Tonly_alpha1_RMSE,mFCQR_Comb_alpha1_RMSE,mFCQR_Hard_alpha1_RMSE,Trans_mFCQR_alpha1_RMSE),
                alpha2_RMSE=rbind(mFCQR_Tonly_alpha2_RMSE,mFCQR_Comb_alpha2_RMSE,mFCQR_Hard_alpha2_RMSE,Trans_mFCQR_alpha2_RMSE),
                alphahat = list(Tar_TOnly$alphahat,mFCQR_Comb$alphahat,mFCQR_Hard$alphahat,Trans_alpha1,Trans_alpha2),
                mse_res = rbind(mFCQR_Tonly_loss,mFCQR_Comb_loss,mFCQR_Hard_loss,Trans_mFCQR_loss),
                Trans_alpha1 = Trans_alpha1,
                Trans_alpha2 = Trans_alpha2,
                Comb = mFCQR_Comb,
                Tonly = Tar_TOnly,
                Hard = mFCQR_Hard)
  return(result)
}


S = 100
alpha1_RMSE_result  = array(0, dim = c(4,length(taugrd),S))
alpha2_RMSE_result  = array(0, dim = c(4,length(taugrd),S))
mse_result = array(0,dim=(c(4,length(taugrd),S)))
Comb_alpha1 = matrix(0,m,length(taugrd))
Comb_alpha2 = matrix(0,m,length(taugrd))
Hard_alpha1 = matrix(0,m,length(taugrd))
Hard_alpha2 = matrix(0,m,length(taugrd))
Tonly_alpha1 = matrix(0,m,length(taugrd))
Tonly_alpha2 = matrix(0,m,length(taugrd))
Trans_alpha1 = matrix(0,m,length(taugrd))
Trans_alpha2 = matrix(0,m,length(taugrd))

result = Simu_fun2(tgrid,v,n0,n1,K,sdx1,sdx2,rhox,sde,sxi,b0,cenrate,seed,nbasis,norder,taugrd,intercept,S2, S4,h)
for (i in 1:S){
  try({
    result = Simu_fun2(tgrid,v,n0,n1,K,sdx1,sdx2,rhox,sde,sxi,b0,cenrate,seed,nbasis,norder,taugrd,intercept,S2, S4,h)
  }, silent=TRUE)
  alpha1_RMSE_result[,,i] = result$alpha1_RMSE
  alpha2_RMSE_result[,,i] = result$alpha2_RMSE
  mse_result[,,i] = result$mse_res
  print(i)
}

alpha1_RMSE_mean = apply(abs(alpha1_RMSE_result),c(1,2),mean,na.rm=TRUE)
alpha1_RMSE_SE = apply(alpha1_RMSE_result,c(1,2),sd,na.rm=TRUE)
alpha2_RMSE_mean = apply(abs(alpha2_RMSE_result),c(1,2),mean,na.rm=TRUE)
alpha2_RMSE_SE = apply(alpha2_RMSE_result,c(1,2),sd,na.rm=TRUE)
mse_mean = apply(mse_result/(n0*length(taugrd)),c(1,2),mean,na.rm=TRUE)
mse_se = apply(mse_result/(n0*length(taugrd)),c(1,2),sd,na.rm=TRUE)

alpha1_RMSE_mean[,qv]
alpha1_RMSE_SE[,qv]
alpha2_RMSE_mean[,qv]
alpha2_RMSE_SE[,qv]
mse_mean[,qv]
mse_se[,qv]

