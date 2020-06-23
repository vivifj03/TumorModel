#Following code runs for all patch definitions

# Packages ----------------------------------------------------------
library(lme4)
library(nlme)
library(LearnBayes)
library(plyr)
library(pracma)
library(ks)
library(matrixcalc)
library(MASS)
library(mvtnorm)
library(loo)
library(scorer)

setwd("")
load("") #RData file with patch definition

filelklhd<-"dataset_for_lklhd200425.RData"
filestimates<-"estimates200425.RData"
filepred<-"predictions200425.RData"


# Initial values ----------------------------------------------------
dataset<-dataset[order(dataset$patch,dataset$date),]

total.inicio<-Sys.time()
###### Mixed effects logistic regression -  gets initial values 
#This models is using time in years
inicio.naive<-Sys.time()
md1 <- glmer(nonZeroHD ~ Dose_med + years + years2  + (1 | patch), 
             data = dataset, family = binomial)
temp<-dataset[is.na(dataset$Dense.lt)==F | is.na(dataset$Hazy.lt)==F,]
###Dense
md2d <- lme(Dense.lt ~ Dose_med + years + years2, random = ~ 1| patch, 
            data = temp)
md2d2 <- lmer(Dense.lt ~ Dose_med + years + years2 + (1| patch), 
              data = temp, control = lmerControl(optCtrl=
                                                   list(maxfun=2e10)))
md2d2.r<-ranef(md2d2)
save.image(filelklhd)
###Hazy
md2h <- lme(Hazy.lt ~ Dose_med + years + years2, random= ~ 1| patch, 
            data = temp)
md2h2 <- lmer(Hazy.lt ~ Dose_med + years + years2  + (1| patch), 
              data = temp, control = lmerControl(optCtrl=
                                                   list(maxfun=2e10)))
save.image(filelklhd)

#Initial values
#Ocurrence variable
Beta<- as.matrix(md1@beta,ncol=1)
if (((md1@theta)^2)<0.1) {
  Psi_cc<-matrix(0.1,1,1)
} else {
  Psi_cc<-matrix(md1@theta)^2
}

#Intensity variable
Gamma<-c(md2d2@beta,md2h2@beta)
Psi_dd<-matrix(c(as.numeric(VarCorr(md2d)[1]),
                 sqrt(as.numeric(VarCorr(md2d)[1]))*
                   sqrt(as.numeric(VarCorr(md2h)[1]))*0.1,
                 sqrt(as.numeric(VarCorr(md2d)[1]))*
                   sqrt(as.numeric(VarCorr(md2h)[1]))*0.1,
                 as.numeric(VarCorr(md2h)[1])),
               nrow = length(md2d2@theta)+length(md2h2@theta),
               ncol = length(md2d2@theta)+length(md2h2@theta))
corr<-cor(temp$Dense.lt, temp$Hazy.lt)
Var.V1<-attr(VarCorr(md2d2), "sc")^2
Var.V2<-attr(VarCorr(md2h2), "sc")^2  

#Corr between random efects between models
Psi_dc<-matrix(c(sqrt(diag(Psi_dd))%*%sqrt(Psi_cc)*0.2),
               nrow=ncol(Psi_dd),ncol=ncol(Psi_cc)) 

m<-0.55 #set probability for normal mixture distribution for 
#random effects 
mu_1<-c(-0.05,0.05) #Set mean first normal dristribution for
#random effects

par.l<-list(Beta=Beta,Gamma=Gamma, Psi_cc=Psi_cc,Psi_dd=Psi_dd,
            Var.V1=Var.V1,Var.V2=Var.V2,corr=corr,
            Psi_dc=Psi_dc,m=m,mu_1=mu_1)
par<-c(Beta=Beta,Gamma=Gamma, Psi_cc=Psi_cc,Psi_dd=Psi_dd,
       Var.V1=Var.V1,Var.V2=Var.V2,corr=corr,
       Psi_dc=Psi_dc,m=m,mu_1=mu_1)

patch<-dataset$patch
patch<-unique(patch)

patch.ch<-dataset[dataset$chageID==1,]$patch
patch.ch<-unique(patch.ch)

patch.nch<-dataset[dataset$chageID==0,]$patch
patch.nch<-unique(patch.nch)

model.data<-list(datafile=dataset, id.var = patch,
                 X.var=c("Dose_med","years","years2"),Z.var=NULL,
                 X_a.var=c("Dose_med","years","years2"),Z_a.var=NULL,
                 U.var="nonZeroHD",V1="Dense.lt0",V2="Hazy.lt0")

save.image(filelklhd)


# Likelihood - Model ------------------------------------------------
source("tpmemmzl.R")

# Estimation of fixed effects ======================================
inicio.coef<-Sys.time()
maxit=1000
step=1
criterion=0
estimate<- vector("list", maxit+1) 
estimate[[1]]<-par.l
dif.theta<- vector("list", maxit)
rel.dif.theta<- vector("list", maxit)
dev<-vector("numeric", maxit)
reldev<-vector("numeric", maxit)
eu.dis<- vector("numeric", maxit) 

theta<-c(Beta=par.l$Beta,
         Gamma=par.l$Gamma,
         mu=par.l$mu_1,
         Sigma.vech=vech(matrix(c(par.l$Var.V1,rep(par.l$corr*
                                                     sqrt(par.l$Var.V1)*sqrt(par.l$Var.V2),2),
                                  par.l$Var.V2),nrow = 2,ncol = 2)),
         Psi_cc.vech=vech(par.l$Psi_cc),
         Psi_dc.i=par.l$Psi_dc %*% ginv(par.l$Psi_cc),
         H.vech=vech(par.l$Psi_dd-(par.l$Psi_dc %*% 
                                     ginv(par.l$Psi_cc) %*% t(par.l$Psi_dc)) ))

in.dev<-lapply(patch,function(x) lklhd.dev(id=x,datafile=dataset,
                                           initial=par.l,
                                           X.var=c("Dose_med","years","years2"),Z.var=NULL,
                                           X_a.var=c("Dose_med","years","years2"),Z_a.var=NULL,
                                           U.var="nonZeroHD",V1="Dense.lt0",V2="Hazy.lt0"))
t.dev<-matrix(0,nrow=length(in.dev), ncol = 1)
for (i in 1:length(in.dev)) {
  t.dev[i,1]<-in.dev[[i]]
}
dev[[1]]<- -2*colSums(t.dev)
print(dev[[1]])

while (criterion<=0) {
  print(step)
  inicio2<-Sys.time()
  
  par.l0<-par.l
  theta<-c(Beta=par.l$Beta,
           Gamma=par.l$Gamma,
           mu=par.l$mu_1,
           Sigma.vech=vech(matrix(c(par.l$Var.V1,rep(par.l$corr*
                                                       sqrt(par.l$Var.V1)*sqrt(par.l$Var.V2),2),
                                    par.l$Var.V2),nrow = 2,ncol = 2)),
           Psi_cc.vech=vech(par.l$Psi_cc),
           Psi_dc.i=par.l$Psi_dc %*% ginv(par.l$Psi_cc),
           H.vech=vech(par.l$Psi_dd-(par.l$Psi_dc %*% 
                                       ginv(par.l$Psi_cc) %*% t(par.l$Psi_dc)) ))
  theta0<-theta
  
  t<-lapply(patch,function(x) lklhd.em(id=x,datafile=dataset,
                                       initial=par.l,
                                       X.var=c("Dose_med","years","years2"),Z.var=NULL,
                                       X_a.var=c("Dose_med","years","years2"),Z_a.var=NULL,
                                       U.var="nonZeroHD",V1="Dense.lt0",V2="Hazy.lt0"))
  
  print("Done with lklhd.em")
  subm<-matrix(0,nrow=length(t), ncol = 1)
  for (i in 1:length(t)) {
    subm[i,1]<-t[[i]]$m  
  }
  score<-matrix(0,nrow=length(t), ncol = length(t[[1]]$scorei))
  for (i in 1:length(t)) {
    score[i,]<-t[[i]]$scorei  
  }
  
  trunk<- c(18:23)
  for (i in trunk) {
    cuts<-0+3*sd(score[,i])
    for (j in 1:nrow(score)) {
      score[j,i]<-ifelse(abs(score[j,i])>cuts,0,score[j,i])   
    }
  }
  
  score<-colSums(score,na.rm = T)
  theta<-theta+score
  
  #Check values of Sigma, Psi_cc
  while (theta[15]<=0)
  {
    test0<-theta0[15]
    #test1<-theta[15]
    test1<-0.001
    theta[15]<-(test0+test1)/2
  }
  while (theta[17]<=0)
  {
    test0<-theta0[17]
    #test1<-theta[17]
    test1<-0.001
    theta[17]<-(test0+test1)/2
  }
  corr.test<-theta[16]/(sqrt(theta[15])*sqrt(theta[17]))
  while (corr.test<= (-0.95)) {
    corr.test0<- theta0[16]/(sqrt(theta0[15])*sqrt(theta0[17]))
    #corr.test1<- theta[16]/(sqrt(theta[15])*sqrt(theta[17]))
    corr.test1<- -0.90
    corr.new<- (corr.test0+corr.test1)/2
    theta[16]<-corr.new*sqrt(theta[15])*sqrt(theta[17])
    corr.test<-theta[16]/(sqrt(theta[15])*sqrt(theta[17]))
  }
  while (corr.test >= 0.95) {
    corr.test0<- theta0[16]/(sqrt(theta0[15])*sqrt(theta0[17]))
    #corr.test1<- theta[16]/(sqrt(theta[15])*sqrt(theta[17]))
    corr.test1<- 0.90
    corr.new<- (corr.test0+corr.test1)/2
    theta[16]<-corr.new*sqrt(theta[15])*sqrt(theta[17])
    corr.test<-theta[16]/(sqrt(theta[15])*sqrt(theta[17]))
  }
  while (theta[18]<=0)
  {
    test0<-theta0[18]
    #test1<-theta[18]
    test1<-0.001
    theta[18]<-(test0+test1)/2
  }
  
  #New estimates
  m<-mean(subm)
  m<-ifelse(m<0.5,(1-m),m)
  while (m>=0.95)
  {
    m0<-par.l$m 
    m<-(0.94+m0)/2
  }
  Beta=matrix(theta[1:4],ncol=1)
  Gamma=matrix(theta[5:12],ncol=1) 
  mu_1=matrix(theta[13:14],ncol=1)
  Var.V1=theta[15]
  Var.V2=theta[17]
  corr=theta[16]/(sqrt(theta[15])*sqrt(theta[17]))
  Psi_cc=matrix(theta[18],ncol=1)
  Psi_dc.i=matrix(theta[19:20],ncol=1)
  Psi_dc=Psi_dc.i%*%Psi_cc
  Psi_cd=t(Psi_dc)
  H<-matrix(c(theta[21],rep(theta[22],2),theta[23]),ncol=2)
  
  while (det(H) <= 0 )
  {
    theta[21] <- (theta0[21]+theta[21])/2
    theta[22] <- (theta0[22]+theta[22])/2
    theta[23] <- (theta0[23]+theta[23])/2
    H<-matrix(c(theta[21],rep(theta[22],2),theta[23]),ncol=2)
  }
  
  while (theta[21]<=0)
  {
    test0<-theta0[21]
    test1<-0.001
    theta[21]<-(test0+test1)/2
  }
  
  while (theta[23]<=0)
  {
    test0<-theta0[23]
    test1<-0.001
    theta[23]<-(test0+test1)/2
  }
  H<-matrix(c(theta[21],theta[22],theta[22],theta[23]),ncol=2)
  
  while (det(H) <= 0 )
  {
    theta[21] <- (theta0[21]+theta[21])/2
    theta[22] <- (theta0[22]+theta[22])/2
    theta[23] <- (theta0[23]+theta[23])/2
    H<-matrix(c(theta[21],rep(theta[22],2),theta[23]),ncol=2)
  }
  Psi_dd=H+(Psi_dc.i %*% Psi_cd)
  
  par.l<-list(Beta=Beta,Gamma=Gamma, Psi_cc=Psi_cc,Psi_dd=Psi_dd,
              Var.V1=Var.V1,Var.V2=Var.V2,corr=corr,
              Psi_dc=Psi_dc,m=m,mu_1=mu_1)
  
  
  step=step+1
  
  #Deviance
  in.dev<-lapply(patch,function(x) lklhd.dev(id=x,datafile=dataset,
                                             initial=par.l,
                                             X.var=c("Dose_med","years","years2"),Z.var=NULL,
                                             X_a.var=c("Dose_med","years","years2"),Z_a.var=NULL,
                                             U.var="nonZeroHD",V1="Dense.lt0",V2="Hazy.lt0"))
  
  t.dev<-matrix(0,nrow=length(in.dev), ncol = 1)
  for (i in 1:length(in.dev)) {
    t.dev[i,1]<-in.dev[[i]]
  }
  dev[[step]]<- -2*colSums(t.dev)
  
  print( paste("new deviance",dev[[step]]))
  #Step-halving in case deviance is larger than previous iteration
  reldev[[step]]<-((dev[[step]] - dev[[step-1]])/(0.1 + 
                                                    abs(dev[[step]])) )
  loop3<-0
  if (reldev[[step]]>=1e-8){
    ii<-1
    while (reldev[[step]] > -1e-8) {
      if (ii>50) break
      ii<-ii+1
      loop3<-ifelse(ii==50,1,loop3)
      par.l$Beta<-(par.l$Beta+par.l0$Beta)/2
      par.l$Gamma<-(par.l$Gamma+par.l0$Gamma)/2
      par.l$Psi_cc<-(par.l$Psi_cc+par.l0$Psi_cc)/2
      par.l$Psi_dd<-(par.l$Psi_dd+par.l0$Psi_dd)/2
      par.l$Var.V1<-(par.l$Var.V1+par.l0$Var.V1)/2
      par.l$Var.V2<-(par.l$Var.V2+par.l0$Var.V2)/2
      par.l$corr<-(par.l$corr+par.l0$corr)/2
      par.l$Psi_dc<-(par.l$Psi_dc+par.l0$Psi_dc)/2
      par.l$m<-(par.l$m+par.l0$m)/2
      par.l$mu_1<-(par.l$mu_1+par.l0$mu_1)/2
      
      in.dev<-lapply(patch,function(x) lklhd.dev(id=x,
                                                 datafile=dataset,initial=par.l,
                                                 X.var=c("Dose_med","years","years2"),
                                                 Z.var=NULL,
                                                 X_a.var=c("Dose_med","years","years2"),
                                                 Z_a.var=NULL,
                                                 U.var="nonZeroHD",V1="Dense.lt0",
                                                 V2="Hazy.lt0"))
      t.dev<-matrix(0,nrow=length(in.dev), ncol = 1)
      for (i in 1:length(in.dev)) {
        t.dev[i,1]<-in.dev[[i]]
      }
      dev[[step]]<- -2*colSums(t.dev)
      reldev[[step]]<-((dev[[step]] - dev[[step-1]])/(0.1 + 
                                                        abs(dev[[step]])) )
      print( paste("base deviance",dev[[step-1]], "new deviance",
                   dev[[step]]))
    }
  }
  
  reldev[[step]]<-abs((dev[[step]] - dev[[step-1]])/(0.1 + 
                                                       abs(dev[[step]])) )
  theta<-c(Beta=par.l$Beta,
           Gamma=par.l$Gamma,
           mu=par.l$mu_1,
           Sigma.vech=vech(matrix(c(par.l$Var.V1,rep(par.l$corr*
                                                       sqrt(par.l$Var.V1)*sqrt(par.l$Var.V2),2),
                                    par.l$Var.V2),nrow = 2,ncol = 2)),
           Psi_cc.vech=vech(par.l$Psi_cc),
           Psi_dc.i=par.l$Psi_dc %*% solve(par.l$Psi_cc),
           H.vech=vech(par.l$Psi_dd-(par.l$Psi_dc %*% 
                                       solve(par.l$Psi_cc) %*% t(par.l$Psi_dc)) ))
  
  estimate[[step]]<-par.l
  eu.dis[step] <-(dist(rbind(theta0,theta), method = "euclidean")^2)/
    length(theta)
  dif.theta[[step]] <- (theta0-theta)
  rel.dif.theta[[step]] <- (theta-theta0)/theta
  absdev<-ifelse(step==2,1,abs(dev[[step]] - dev[[step-1]]))
  
  
  criterion.step<-ifelse(step<=maxit,0,1)
  criterion.distance<-ifelse(max(abs(dif.theta[[step]]))>1e-8,0,1)
  criterion.max<-ifelse(max(abs(rel.dif.theta[[step]]))>1e-8,0,1)
  criterion.dev<-ifelse(absdev>1e-8,0,1)
  criterion.reldev<-ifelse(reldev[[step]]>1e-8,0,1)
  
  criterion<-criterion.step+criterion.distance+criterion.max+
    criterion.dev+criterion.reldev
  
  criterion<-ifelse(loop3==1 & step<5,0,criterion)
  
  print(c(criterion.step=criterion.step,criterion.distance=
            criterion.distance,criterion.max=criterion.max,
          criterion=criterion))
  print(c(eu.dis=eu.dis[step],max.abs.dif=max(abs(dif.theta[[step]])),
          max.rel=max(abs(rel.dif.theta[[step]])),
          dev=dev[step],absdev.dif=absdev,reldev=reldev[step]))
  
  end2<-Sys.time()
  timeit<- end2-inicio2
  print(theta)
  print(timeit)
  
}
end.coef<-Sys.time()
end.coef-inicio.coef

# For figures of coefficientes at each iteration
Beta<-matrix(NA,step,4)
for(i in 1:step){
  Beta[i,]<-(estimate[[i]]$Beta)
}

Gamma<-matrix(NA,step,8)
for(i in 1:step){
  Gamma[i,]<-(estimate[[i]]$Gamma)
}

Psi_cc<-matrix(NA,step,1)
for(i in 1:step){
  Psi_cc[i,]<-(estimate[[i]]$Psi_cc)
}

Psi_dd<-matrix(NA,step,4)
for(i in 1:step){
  Psi_dd[i,]<-(estimate[[i]]$Psi_dd)
}

Var.V1<-matrix(NA,step,1)
for(i in 1:step){
  Var.V1[i,]<-(estimate[[i]]$Var.V1)
}
Var.V2<-matrix(NA,step,1)
for(i in 1:step){
  Var.V2[i,]<-(estimate[[i]]$Var.V2)
}
corr<-matrix(NA,step,1)
for(i in 1:step){
  corr[i,]<-(estimate[[i]]$corr)
}

m<-matrix(NA,step,1)
for(i in 1:step){
  m[i,]<-(estimate[[i]]$m)
}
mu_1<-matrix(NA,step,2)
for(i in 1:step){
  mu_1[i,]<-(estimate[[i]]$mu_1)
}

save.image(filestimates)

### Standard errors ###############################################
se<-lapply(patch,function(x) lklhd.var(id=x,datafile=dataset,
                                       est=par.l,
                                       X.var=c("Dose_med","years","years2"),Z.var=NULL,
                                       X_a.var=c("Dose_med","years","years2"),Z_a.var=NULL,
                                       U.var="nonZeroHD",V1="Dense.lt0",V2="Hazy.lt0"))

S_0b.total<-matrix(0,nrow=4, ncol = 4)
S_0g.total<-matrix(0,nrow=8, ncol = 8)
S_0mu.total<-matrix(0,nrow=2, ncol = 2)

for (i in 1:length(se)) {
  S_0b.total<-S_0b.total+se[[i]]$S_0b
  S_0g.total<-S_0g.total+se[[i]]$S_0g
  S_0mu.total<-S_0mu.total+se[[i]]$S_0mu
}
S_0b.total.i<-ginv(S_0b.total)
S_0g.total.i<-ginv(S_0g.total)
S_0mu.total.i<-ginv(S_0mu.total)

sqrt(diag(S_0b.total.i))
sqrt(diag(S_0g.total.i))

par.final<-c(unlist(par.l))
se.final1<-c(sqrt(diag(S_0b.total.i)),sqrt(diag(S_0g.total.i)))
lci1<-par.final[1:12]-(1.96*se.final1)
uci1<-par.final[1:12]+(1.96*se.final1)
save.image(filestimates)

### Random effects ################################################
try(if(criterion.step==1) stop("estimation did not converge"))

iniciop<-Sys.time()
random.est <- lapply(patch,function(x) lklhd.ram(id=x,
                                                 datafile=dataset,est=par.l,
                                                 X.var=c("Dose_med","years","years2"),
                                                 Z.var=NULL,
                                                 X_a.var=c("Dose_med","years","years2"),
                                                 Z_a.var=NULL,
                                                 U.var="nonZeroHD",V1="Dense.lt0",
                                                 V2="Hazy.lt0"))

endp<-Sys.time()
endp-iniciop
save.image(filepred)

### Predicted values ################################################
ran.est <- as.data.frame(matrix(NA,length(patch),7))
ran.est[,1] <- patch
for(i in 1:length(patch)){
  ran.est[i,2]<-(random.est[[i]]$ram.est$c.est)
  ran.est[i,3]<-(random.est[[i]]$ram.est$cc.t.est)
  ran.est[i,4]<-(random.est[[i]]$ram.est$d1.est)
  ran.est[i,5]<-(random.est[[i]]$ram.est$d2.est)
  ran.est[i,6]<-(random.est[[i]]$ram.est$ess1)
  ran.est[i,7]<-(random.est[[i]]$ram.est$ess2)
}
ran.est <- as.data.frame(ran.est)
colnames(ran.est) <- c("patch","c.est","cc.t.est","d1.est","d2.est",
                       "ess","ess2")
dataset<-dataset[,1:38]
dataset<-merge(dataset,ran.est,by="patch",all.x = T)
dataset$logp.1 <-  par.l$Beta[1]  +  par.l$Beta[2]*dataset$Dose_med + 
  par.l$Beta[3] * dataset$years + par.l$Beta[4] * dataset$years2 +
  dataset$c.est
dataset$odds.1 <- exp(dataset$logp.1)
dataset$p.1 <- dataset$odds.1 / (1 + dataset$odds.1)
dataset$subpred1 <-  par.l$Gamma[1]  +  par.l$Gamma[2]*dataset$Dose_med + 
  par.l$Gamma[3] * dataset$years + par.l$Gamma[4] * dataset$years2 +
  dataset$d1.est
dataset$subpred2 <-  par.l$Gamma[5]+par.l$Gamma[6]*dataset$Dose_med + 
  par.l$Gamma[7] * dataset$years + par.l$Gamma[8] * dataset$years2 +
  dataset$d2.est
dataset$subres1 <-   dataset$Dense.lt0-dataset$subpred1
dataset$subres2 <-   dataset$Hazy.lt0-dataset$subpred2  
dataset$pred1<- dataset$p.1*dataset$subpred1
dataset$pred2<- dataset$p.1*dataset$subpred2
dataset$residuals1<-(dataset$Dense.lt0-dataset$pred1)
dataset$residuals2<-(dataset$Hazy.lt0-dataset$pred2)
dataset$sresiduals1<-dataset$residuals1^2
dataset$sresiduals2<-dataset$residuals2^2

duan1<- mean(exp(dataset[dataset$nonZeroHD==1,"subres1"]))
duan2<- mean(exp(dataset[dataset$nonZeroHD==1,"subres2"]))

dataset$pred1.e<- dataset$p.1*exp(dataset$subpred1)*duan1
dataset$pred2.e<- dataset$p.1*exp(dataset$subpred2)*duan2
dataset$residuals1.e<-(dataset$Dense.t0-dataset$pred1.e)
dataset$residuals2.e<-(dataset$Hazy.t0-dataset$pred2.e)

ram.m_Vu.V1<-mean(dataset$Dense.t0)
ram.TSS.V1<-sum((dataset$Dense.t0-ram.m_Vu.V1)^2)
ram.sse.V.V1<-sum((dataset$Dense.t0-dataset$pred1.e)^2)
ram.ssr.V.V1<-sum((dataset$pred1.e-ram.m_Vu.V1)^2)
ram.tss.V1<-ram.sse.V.V1+ram.ssr.V.V1
ram.corr2.V1<-cor(dataset$Dense.t0,dataset$pred1.e)^2     
ram.R2.V1<-1-(var(dataset$residuals1.e)/var(dataset$Dense.t0))
ram.mse.V1<-mean_squared_error(dataset$Dense.t0,dataset$pred1.e)
ram.mae.V1<-mean_absolute_error(dataset$Dense.t0,dataset$pred1.e)
ram.mase.V1<-mean_absolute_scaled_error(dataset$Dense.t0,dataset$pred1.e)
ram.rmse.V1<-rmse(dataset$Dense.t0,dataset$pred1.e)
ram.tvs.V1<-total_variance_score(dataset$Dense.t0,dataset$pred1.e)
ram.uvs.V1<-unexplained_variance_score(dataset$Dense.t0,dataset$pred1.e)
ram.evs.V1<-explained_variance_score(dataset$Dense.t0,dataset$pred1.e) 

ram.m_Vu.V2<-mean(dataset$Hazy.t0)
ram.TSS.V2<-sum((dataset$Hazy.t0-ram.m_Vu.V2)^2)
ram.sse.V.V2<-sum((dataset$Hazy.t0-dataset$pred2.e)^2)
ram.ssr.V.V2<-sum((dataset$pred2.e-ram.m_Vu.V2)^2)
ram.tss.V2<-ram.sse.V.V2+ram.ssr.V.V2
ram.corr2.V2<-cor(dataset$Hazy.t0,dataset$pred2.e)^2     
ram.R2.V2<-1-(var(dataset$residuals2.e)/var(dataset$Hazy.t0))
ram.mse.V2<-mean_squared_error(dataset$Hazy.t0,dataset$pred2.e)
ram.mae.V2<-mean_absolute_error(dataset$Hazy.t0,dataset$pred2.e)
ram.mase.V2<-mean_absolute_scaled_error(dataset$Hazy.t0,dataset$pred2.e)
ram.rmse.V2<-rmse(dataset$Hazy.t0,dataset$pred2.e)
ram.tvs.V2<-total_variance_score(dataset$Hazy.t0,dataset$pred2.e)
ram.uvs.V2<-unexplained_variance_score(dataset$Hazy.t0,dataset$pred2.e)
ram.evs.V2<-explained_variance_score(dataset$Hazy.t0,dataset$pred2.e) 

dataset.t<-dataset[,1:38]
ran.est <- as.data.frame(matrix(NA,length(patch),7))
ran.est[,1] <- patch
for(i in 1:length(patch)){
  ran.est[i,2]<-(random.est[[i]]$ram.est.t$c.est)
  ran.est[i,3]<-(random.est[[i]]$ram.est.t$cc.t.est)
  ran.est[i,4]<-(random.est[[i]]$ram.est.t$d1.est)
  ran.est[i,5]<-(random.est[[i]]$ram.est.t$d2.est)
  ran.est[i,6]<-(random.est[[i]]$ram.est.t$ess1)
  ran.est[i,7]<-(random.est[[i]]$ram.est.t$ess2)
}
ran.est <- as.data.frame(ran.est)
colnames(ran.est) <- c("patch","c.est","cc.t.est","d1.est","d2.est",
                       "ess","ess2")
dataset.t<-merge(dataset.t,ran.est,by="patch",all.x = T)
dataset.t$logp.1 <-  par.l$Beta[1]+ par.l$Beta[2]*dataset.t$Dose_med + 
  par.l$Beta[3] * dataset.t$years + par.l$Beta[4] * dataset.t$years2 +
  dataset.t$c.est
dataset.t$odds.1 <- exp(dataset.t$logp.1)
dataset.t$p.1 <- dataset.t$odds.1 / (1 + dataset.t$odds.1)
dataset.t$subpred1 <-  par.l$Gamma[1]+par.l$Gamma[2]*dataset.t$Dose_med + 
  par.l$Gamma[3] * dataset.t$years + par.l$Gamma[4] * dataset.t$years2 +
  dataset.t$d1.est
dataset.t$subpred2 <-  par.l$Gamma[5]+par.l$Gamma[6]*dataset.t$Dose_med + 
  par.l$Gamma[7] * dataset.t$years + par.l$Gamma[8] * dataset.t$years2 +
  dataset.t$d2.est

dataset.t$subres1 <-   dataset.t$Dense.lt0-dataset.t$subpred1
dataset.t$subres2 <-   dataset.t$Hazy.lt0-dataset.t$subpred2  
dataset.t$pred1<- dataset.t$p.1*dataset.t$subpred1
dataset.t$pred2<- dataset.t$p.1*dataset.t$subpred2
dataset.t$residuals1<-(dataset.t$Dense.lt0-dataset.t$pred1)
dataset.t$residuals2<-(dataset.t$Hazy.lt0-dataset.t$pred2)
dataset.t$sresiduals1<-dataset.t$residuals1^2
dataset.t$sresiduals2<-dataset.t$residuals2^2

duan1<- mean(exp(dataset.t[dataset.t$nonZeroHD==1,"subres1"]))
duan2<- mean(exp(dataset.t[dataset.t$nonZeroHD==1,"subres2"]))
dataset.t$pred1.e<- dataset.t$p.1*exp(dataset.t$subpred1)*duan1
dataset.t$pred2.e<- dataset.t$p.1*exp(dataset.t$subpred2)*duan2
dataset.t$residuals1.e<-(dataset.t$Dense.t0-dataset.t$pred1.e)
dataset.t$residuals2.e<-(dataset.t$Hazy.t0-dataset.t$pred2.e)

ram.t.m_Vu.V1<-mean(dataset.t$Dense.t0)
ram.t.TSS.V1<-sum((dataset.t$Dense.t0-ram.t.m_Vu.V1)^2)
ram.t.sse.V.V1<-sum((dataset.t$Dense.t0-dataset.t$pred1.e)^2)
ram.t.ssr.V.V1<-sum((dataset.t$pred1.e-ram.t.m_Vu.V1)^2)
ram.t.tss.V1<-ram.t.sse.V.V1+ram.t.ssr.V.V1
ram.t.corr2.V1<-cor(dataset.t$Dense.t0,dataset.t$pred1.e)^2 
ram.t.R2.V1<-1-(var(dataset.t$residuals1.e)/var(dataset.t$Dense.t0))
ram.t.mse.V1<-mean_squared_error(dataset.t$Dense.t0,dataset.t$pred1.e)
ram.t.mae.V1<-mean_absolute_error(dataset.t$Dense.t0,dataset.t$pred1.e)
ram.t.mase.V1<-mean_absolute_scaled_error(dataset.t$Dense.t0,
                                          dataset.t$pred1.e)
ram.t.rmse.V1<-rmse(dataset.t$Dense.t0,dataset.t$pred1.e)
ram.t.tvs.V1<-total_variance_score(dataset.t$Dense.t0,dataset.t$pred1.e)
ram.t.uvs.V1<-unexplained_variance_score(dataset.t$Dense.t0,
                                         dataset.t$pred1.e)
ram.t.evs.V1<-explained_variance_score(dataset.t$Dense.t0,
                                       dataset.t$pred1.e) 

ram.t.m_Vu.V2<-mean(dataset.t$Hazy.t0)
ram.t.TSS.V2<-sum((dataset.t$Hazy.t0-ram.t.m_Vu.V2)^2)
ram.t.sse.V.V2<-sum((dataset.t$Hazy.t0-dataset.t$pred2.e)^2)
ram.t.ssr.V.V2<-sum((dataset.t$pred2.e-ram.t.m_Vu.V2)^2)
ram.t.tss.V2<-ram.t.sse.V.V2+ram.t.ssr.V.V2
ram.t.corr2.V2<-cor(dataset.t$Hazy.t0,dataset.t$pred2.e)^2
ram.t.R2.V2<-1-(var(dataset.t$residuals2.e)/var(dataset.t$Hazy.t0))
ram.t.mse.V2<-mean_squared_error(dataset.t$Hazy.t0,dataset.t$pred2.e)
ram.t.mae.V2<-mean_absolute_error(dataset.t$Hazy.t0,dataset.t$pred2.e)
ram.t.mase.V2<-mean_absolute_scaled_error(dataset.t$Hazy.t0,
                                          dataset.t$pred2.e)
ram.t.rmse.V2<-rmse(dataset.t$Hazy.t0,dataset.t$pred2.e)
ram.t.tvs.V2<-total_variance_score(dataset.t$Hazy.t0,dataset.t$pred2.e)
ram.t.uvs.V2<-unexplained_variance_score(dataset.t$Hazy.t0,
                                         dataset.t$pred2.e)
ram.t.evs.V2<-explained_variance_score(dataset.t$Hazy.t0,dataset.t$pred2.e)

dataset.p<-dataset[,1:38]
ran.est <- as.data.frame(matrix(NA,length(patch),7))
ran.est[,1] <- patch
for(i in 1:length(patch)){
  ran.est[i,2]<-(random.est[[i]]$ram.est.p$c.est)
  ran.est[i,3]<-(random.est[[i]]$ram.est.p$cc.t.est)
  ran.est[i,4]<-(random.est[[i]]$ram.est.p$d1.est)
  ran.est[i,5]<-(random.est[[i]]$ram.est.p$d2.est)
  ran.est[i,6]<-(random.est[[i]]$ram.est.p$ess1)
  ran.est[i,7]<-(random.est[[i]]$ram.est.p$ess2)
}
ran.est <- as.data.frame(ran.est)
colnames(ran.est) <- c("patch","c.est","cc.t.est","d1.est","d2.est",
                       "ess","ess2")
dataset.p<-merge(dataset.p,ran.est,by="patch",all.x = T)

dataset.p$logp.1 <-  par.l$Beta[1]+par.l$Beta[2]*dataset.p$Dose_med + 
  par.l$Beta[3] * dataset.p$years + par.l$Beta[4] * dataset.p$years2 +
  dataset.p$c.est
dataset.p$odds.1 <- exp(dataset.p$logp.1)
dataset.p$p.1 <- dataset.p$odds.1 / (1 + dataset.p$odds.1)

dataset.p$subpred1 <- par.l$Gamma[1]+par.l$Gamma[2]*dataset.p$Dose_med + 
  par.l$Gamma[3] * dataset.p$years + par.l$Gamma[4] * dataset.p$years2 +
  dataset.p$d1.est
dataset.p$subpred2 <- par.l$Gamma[5]+par.l$Gamma[6]*dataset.p$Dose_med + 
  par.l$Gamma[7] * dataset.p$years + par.l$Gamma[8] * dataset.p$years2 +
  dataset.p$d2.est
dataset.p$subres1 <-   dataset.p$Dense.lt0-dataset.p$subpred1
dataset.p$subres2 <-   dataset.p$Hazy.lt0-dataset.p$subpred2  
dataset.p$pred1<- dataset.p$p.1*dataset.p$subpred1
dataset.p$pred2<- dataset.p$p.1*dataset.p$subpred2
dataset.p$residuals1<-(dataset.p$Dense.lt0-dataset.p$pred1)
dataset.p$residuals2<-(dataset.p$Hazy.lt0-dataset.p$pred2)
dataset.p$sresiduals1<-dataset.p$residuals1^2
dataset.p$sresiduals2<-dataset.p$residuals2^2

duan1<- mean(exp(dataset.p[dataset.p$nonZeroHD==1,"subres1"]))
duan2<- mean(exp(dataset.p[dataset.p$nonZeroHD==1,"subres2"]))

dataset.p$pred1.e<- dataset.p$p.1*exp(dataset.p$subpred1)*duan1
dataset.p$pred2.e<- dataset.p$p.1*exp(dataset.p$subpred2)*duan2
dataset.p$residuals1.e<-(dataset.p$Dense.t0-dataset.p$pred1.e)
dataset.p$residuals2.e<-(dataset.p$Hazy.t0-dataset.p$pred2.e)

ram.p.m_Vu.V1<-mean(dataset.p$Dense.t0)
ram.p.TSS.V1<-sum((dataset.p$Dense.t0-ram.p.m_Vu.V1)^2)
ram.p.sse.V.V1<-sum((dataset.p$Dense.t0-dataset.p$pred1.e)^2)
ram.p.ssr.V.V1<-sum((dataset.p$pred1.e-ram.p.m_Vu.V1)^2)
ram.p.tss.V1<-ram.p.sse.V.V1+ram.p.ssr.V.V1
ram.p.corr2.V1<-cor(dataset.p$Dense.t0,dataset.p$pred1.e)^2
ram.p.R2.V1<-1-(var(dataset.p$residuals1.e)/var(dataset.p$Dense.t0))
ram.p.mse.V1<-mean_squared_error(dataset.p$Dense.t0,dataset.p$pred1.e)
ram.p.mae.V1<-mean_absolute_error(dataset.p$Dense.t0,dataset.p$pred1.e)
ram.p.mase.V1<-mean_absolute_scaled_error(dataset.p$Dense.t0,
                                          dataset.p$pred1.e)
ram.p.rmse.V1<-rmse(dataset.p$Dense.t0,dataset.p$pred1.e)
ram.p.tvs.V1<-total_variance_score(dataset.p$Dense.t0,dataset.p$pred1.e)
ram.p.uvs.V1<-unexplained_variance_score(dataset.p$Dense.t0,
                                         dataset.p$pred1.e)
ram.p.evs.V1<-explained_variance_score(dataset.p$Dense.t0,
                                       dataset.p$pred1.e)

ram.p.m_Vu.V2<-mean(dataset.p$Hazy.t0)
ram.p.TSS.V2<-sum((dataset.p$Hazy.t0-ram.p.m_Vu.V2)^2)
ram.p.sse.V.V2<-sum((dataset.p$Hazy.t0-dataset.p$pred2.e)^2)
ram.p.ssr.V.V2<-sum((dataset.p$pred2.e-ram.p.m_Vu.V2)^2)
ram.p.tss.V2<-ram.p.sse.V.V2+ram.p.ssr.V.V2
ram.p.corr2.V2<-cor(dataset.p$Hazy.t0,dataset.p$pred2.e)^2
ram.p.R2.V2<-1-(var(dataset.p$residuals2.e)/var(dataset.p$Hazy.t0))
ram.p.mse.V2<-mean_squared_error(dataset.p$Hazy.t0,dataset.p$pred2.e)
ram.p.mae.V2<-mean_absolute_error(dataset.p$Hazy.t0,dataset.p$pred2.e)
ram.p.mase.V2<-mean_absolute_scaled_error(dataset.p$Hazy.t0,
                                          dataset.p$pred2.e)
ram.p.rmse.V2<-rmse(dataset.p$Hazy.t0,dataset.p$pred2.e)
ram.p.tvs.V2<-total_variance_score(dataset.p$Hazy.t0,dataset.p$pred2.e)
ram.p.uvs.V2<-unexplained_variance_score(dataset.p$Hazy.t0,
                                         dataset.p$pred2.e)
ram.p.evs.V2<-explained_variance_score(dataset.p$Hazy.t0,dataset.p$pred2.e) 

total.end<-Sys.time()
total.end-total.inicio
save.image(filepred)


# Results -----------------------------------------------------------
library(compositions)
library(ggtern)
library(png)
library(dplyr)
library(mvnTest)

### Scaled residuals for U ########################################
agg<-dataset.p[,c("nonZeroHD","Dense.lt","Hazy.lt","Dense.lt0",
                  "Hazy.lt0","p.1","subpred1","subpred2","pred1",
                  "pred2","patch")]

agg2<-  agg %>%
  group_by(patch) %>% 
  summarise_each(sum)
agg2$scaled.u.res<-(agg2$nonZeroHD-agg2$p.1)/sqrt(agg2$p.1 * 
                                                    (5-agg2$p.1) /5)
agg2$scaled.u.res.ct<-ifelse(abs(agg2$scaled.u.res)<2.5,0,1)
summary(agg2$scaled.u.res)
table(agg2$scaled.u.res.ct)
prop.table(table(agg2$scaled.u.res.ct))  
table(agg2$scaled.u.res.ct,agg2$nonZeroHD)

### Scaled residuals for v for those with at least one nonzero #####
agg3<-agg2[agg2$nonZeroHD!=0,]
agg3$scaled.v1.res<-(agg3$Dense.lt0 - agg3$subpred1)/
  sqrt(agg3$nonZeroHD*par.l$Var.V1)
agg3$scaled.v2.res<-(agg3$Hazy.lt0 - agg3$subpred2)/
  sqrt(agg3$nonZeroHD*par.l$Var.V2 )
agg3$scaled.v1.res.ct<-ifelse(abs(agg3$scaled.v1.res)<2.5,0,1)
agg3$scaled.v2.res.ct<-ifelse(abs(agg3$scaled.v2.res)<2.5,0,1)
table(agg3$scaled.v1.res.ct)
prop.table(table(agg3$scaled.v1.res.ct))
table(agg3$scaled.v1.res.ct,agg3$nonZeroHD)
table(agg3$scaled.v2.res.ct)
prop.table(table(agg3$scaled.v2.res.ct))
table(agg3$scaled.v2.res.ct,agg3$nonZeroHD)

### Residuals ######################################################
v <- dataset.p[order(dataset.p$patch, dataset.p$date),]
v <- v[!duplicated(v$patch),]
test<-qqnorm(v$c.est, pch = 1, frame = FALSE,col = "#926b8d")
x <- c(quantile(test$x,0.25),quantile(test$x,0.75))
y <- c(quantile(log(test$y+(abs(min(test$y))+1)),0.25),
       quantile(log(test$y+(abs(min(test$y))+1)),0.75))
slope <- diff(y)/diff(x)
intercept <- y[1]-slope*x[1]
plot(test$x,log(test$y+(abs(min(test$y))+1)),col = "#926b8d")
la<-c(0,1,2,3,4) #cubic
la<-c(0,1,2,3,4,5) #spherical
la<-c(0,1,2,3,4) #iso
la2<-round(exp(la)-1-abs(min(test$y)),0)
plot(test$x,log(test$y+(abs(min(test$y))+1)),col = "#926b8d",
     xlab="Theoretical Quantiles", ylab="Sample Quantiles", yaxt="n",
     main="Normal Q-Q Plot")
abline(intercept, slope, col="#B7B6BA", lwd = 2)
axis(2,at=la,labels=la2)

x <- c(quantile(test$x,0.25),quantile(test$x,0.75))
y <- c(quantile(sqrt(test$y+(abs(min(test$y)))),0.25),
       quantile(sqrt(test$y+(abs(min(test$y)))),0.75))
slope <- diff(y)/diff(x)
intercept <- y[1]-slope*x[1]
plot(test$x,sqrt(test$y+(abs(min(test$y)))),col = "#926b8d",
     xlab="Theoretical Quantiles", ylab="Sample Quantiles",
     main="Normal Q-Q Plot")
la<-c(0,2,4,6,8)
la<-c(0,2,4,6,8,10,12,14)
la<-c(0,2,4,6,8)
la2<-round((la^2)-1-abs(min(test$y)),0)
plot(test$x,sqrt(test$y+(abs(min(test$y)))),col = "#926b8d",
     xlab="Theoretical Quantiles", ylab="Sample Quantiles", yaxt="n",
     main="Normal Q-Q Plot")
abline(intercept, slope, col="#B7B6BA", lwd = 2)
axis(2,at=la,labels=la2)

x <- c(quantile(test$x,0.25),quantile(test$x,0.75))
y <- c(quantile(sign(test$y)*(abs(test$y)^(1/3)),0.25),
       quantile(sign(test$y)*(abs(test$y)^(1/3)),0.75))
slope <- diff(y)/diff(x)
intercept <- y[1]-slope*x[1]
plot(test$x,sign(test$y)*(abs(test$y)^(1/3)),col = "#926b8d",
     xlab="Theoretical Quantiles", ylab="Sample Quantiles",
     main="Normal Q-Q Plot")
la<-c(-1,0,1,2,3,4)
la<-c(-4,-2,0,2,4)
la<- c(-4,-3,-2,-1,0,1)
la2<-round((la^3),0)
plot(test$x,sign(test$y)*(abs(test$y)^(1/3)),col = "#926b8d",
     xlab="Theoretical Quantiles", ylab="Sample Quantiles", yaxt="n",
     main="Normal Q-Q Plot")
abline(intercept, slope, col="#B7B6BA", lwd = 2)
axis(2,at=la,labels=la2)

x <- c(quantile(test$x,0.25),quantile(test$x,0.75))
y <- c(quantile(sign(test$y)*(abs(test$y)^(1/5)),0.25),
       quantile(sign(test$y)*(abs(test$y)^(1/5)),0.75))
slope <- diff(y)/diff(x)
intercept <- y[1]-slope*x[1]
plot(test$x,sign(test$y)*(abs(test$y)^(1/5)),col = "#926b8d",
     xlab="Theoretical Quantiles", ylab="Sample Quantiles",
     main="Normal Q-Q Plot")
la<-c(-1,0,1,2)
la<-c(-2,-1,0,1,2)
la<-c(-2,-1,0,1)
la2<-round((la^5),0)
plot(test$x,sign(test$y)*(abs(test$y)^(1/5)),col = "#926b8d",
     xlab="Theoretical Quantiles", ylab="Sample Quantiles", yaxt="n",
     main="Normal Q-Q Plot")
abline(intercept, slope, col="#B7B6BA", lwd = 2)
axis(2,at=la,labels=la2)

qqnorm(v[v$Com!=1,]$c.est, pch = 1, frame = FALSE,col = "#926b8d")
qqline(v[v$Com!=1,]$c.est, col ="#B7B6BA", lwd = 2)

qqnorm(v[v$Com==1,]$c.est, pch = 1, frame = FALSE,col = "#926b8d")
qqline(v[v$Com==1,]$c.est, col ="#B7B6BA", lwd = 2)
boxplot(c.est~nonZeroHD, data=v)

qqnorm(dataset.t$subres1, pch = 1, frame = FALSE)
qqline(dataset.t$subres1, col = "steelblue", lwd = 2)
qqnorm(dataset.t$subres2, pch = 1, frame = FALSE)
qqline(dataset.t$subres2, col = "steelblue", lwd = 2)

qqnorm(dataset.t$residuals1, pch = 1, frame = FALSE)
qqline(dataset.t$residuals1, col = "steelblue", lwd = 2)
qqnorm(dataset.t$residuals2, pch = 1, frame = FALSE)
qqline(dataset.t$residuals2, col = "steelblue", lwd = 2)

AD.test(dataset.t[,c("subres1","subres2")], qqplot=T)
CM.test(dataset.t[,c("subres1","subres2")], qqplot=F)
R.test(dataset.t[,c("subres1","subres2")], qqplot=F)

AD.test(dataset.t[,c("residuals1","residuals2")], qqplot=T)
CM.test(dataset.t[,c("residuals1","residuals2")], qqplot=F)
R.test(dataset.t[,c("residuals1","residuals2")], qqplot=F)

### Additive logistic transform ###################################
alrInv(c(par.l$Gamma[[1]],par.l$Gamma[[5]]))
ef.dose<-c(par.l$Gamma[[2]],par.l$Gamma[[6]])
alrInv(ef.dose)
sqrt(ef.dose %*% solve(diag(2) + c(1,1) %*% t(c(1,1))) %*% (ef.dose))
ef.time<-c(par.l$Gamma[[3]],par.l$Gamma[[7]])
alrInv(ef.time)
sqrt(ef.time %*% solve(diag(2) + c(1,1) %*% t(c(1,1))) %*% (ef.time))
ef.time2<-c(par.l$Gamma[[4]],par.l$Gamma[[8]])
alrInv(ef.time2)
sqrt(ef.time2 %*% solve(diag(2) + c(1,1) %*% t(c(1,1))) %*% (ef.time2))

### Goodness of fit ##############################################
v2<-dataset.p
v2$e.subpred.f1.v2<-exp(v2$subpred1)*exp(0.5*par.l$Var.V1)
v2$e.subpred.f2.v2<-exp(v2$subpred2)*exp(0.5*par.l$Var.V2)
v2$e.pred.f1<-v2$p.1*v2$e.subpred.f1
v2$e.pred.f2<-v2$p.1*v2$e.subpred.f2
v2$dense.pred.f<-v2$e.pred.f1/(1+v2$e.pred.f1+v2$e.pred.f2)
v2$hazy.pred.f<-v2$e.pred.f2/(1+v2$e.pred.f1+v2$e.pred.f2)
v2$com.pred.f<-1/(1+v2$e.pred.f1+v2$e.pred.f2)
v2[,c("dense.pred.f","hazy.pred.f","com.pred.f")]<-clo(v2[,
                                                          c("dense.pred.f","hazy.pred.f","com.pred.f")])

v2[,c("dense.res.f")]<-v2$Dense-v2$dense.pred.f
v2[,c("hazy.res.f")]<-v2$Hazy-v2$hazy.pred.f
v2[,c("com.res.f")]<-v2$Com-v2$com.pred.f

sum(v2$dense.res.f^2 )
sum(v2$hazy.res.f^2 )
sum(v2$com.res.f^2 )
sum(v2$dense.res.f^2 ) + sum(v2$hazy.res.f^2 ) + 
  sum(v2$com.res.f^2 )
sqrt(sum(v2$dense.res.f^2 ) + sum(v2$hazy.res.f^2 ) + 
       sum(v2$com.res.f^2 ))

v2[,c("dense.res.p1")]<-v2$Dense-v2$dense.pred
v2[,c("hazy.res.p1")]<-v2$Hazy-v2$hazy.pred
v2[,c("com.res.p1")]<-v2$Com-v2$com.pred

sum(v2$dense.res.p1^2 )
sum(v2$hazy.res.p1^2 )
sum(v2$com.res.p1^2 )
sum(v2$dense.res.p1^2 ) + sum(v2$hazy.res.p1^2 ) + 
  sum(v2$com.res.p1^2 )
sqrt(sum(v2$dense.res.p1^2 ) + sum(v2$hazy.res.p1^2 ) + 
       sum(v2$com.res.p1^2 ))

v2[,c("dense.res.p2")]<-v2$Dense-v2$dense.pred2
v2[,c("hazy.res.p2")]<-v2$Hazy-v2$hazy.pred2
v2[,c("com.res.p2")]<-v2$Com-v2$com.pred2

sum(v2$dense.res.p2^2 )
sum(v2$hazy.res.p2^2 )
sum(v2$com.res.p2^2 )
sum(v2$dense.res.p2^2 ) + sum(v2$hazy.res.p2^2 ) + 
  sum(v2$com.res.p2^2 )
sqrt(sum(v2$dense.res.p2^2 ) + sum(v2$hazy.res.p2^2 ) + 
       sum(v2$com.res.p2^2 ))

# Figures ===========================================================
#Figures 5.1 and 5.2 are in results/residuals

#Figures 5.3 
#Standardized residuals
dataset.p$subres1.st<-dataset.p$subres1/sqrt(par.l$Var.V1)
dataset.p$subres2.st<-dataset.p$subres2/sqrt(par.l$Var.V2)
dataset.p$residuals1.st<-dataset.p$residuals1/sqrt(par.l$Var.V1)
dataset.p$residuals2.st<-dataset.p$residuals2/sqrt(par.l$Var.V2)

ggplot(dataset.p, aes(x=Dose_med, y=subres1.st)) + 
  geom_point(alpha=0.7,col="#926b8d") +
  xlab("Median Dose") + ylab("Standardized residual - alr(Dense)")+
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = 
          element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = 
          element_line(colour = "black"),text = element_text(size=16),
        legend.position="bottom",legend.title = element_blank())

ggplot(dataset.p, aes(x=Dose_med, y=residuals1.st)) + 
  geom_point(alpha=0.7,col="#926b8d") +
  xlab("Median Dose") + ylab("Standardized residual - alr(Dense)")+
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = 
          element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = 
          element_line(colour = "black"),text = element_text(size=16),
        legend.position="bottom",legend.title = element_blank())

ggplot(dataset.p, aes(x=Dose_med, y=subres2.st)) + 
  geom_point(alpha=0.7,col="#926b8d") +
  xlab("Median Dose") + ylab("Standardized residual - alr(Hazy)")+
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = 
          element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = 
          element_line(colour = "black"),text = element_text(size=16),
        legend.position="bottom",legend.title = element_blank())

ggplot(dataset.p, aes(x=Dose_med, y=residuals2.st)) + 
  geom_point(alpha=0.7,col="#926b8d") +
  xlab("Median Dose") + ylab("Standardized residual - alr(Hazy)")+
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = 
          element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = 
          element_line(colour = "black"),text = element_text(size=16),
        legend.position="bottom",legend.title = element_blank())

ggplot(dataset.p, aes(x=factor(date2), y=subres1.st)) + 
  geom_point(alpha=0.7,col="#926b8d") +
  xlab("Months after RT") + ylab("Standardized residual - alr(Dense)")+
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = 
          element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = 
          element_line(colour = "black"),text = element_text(size=16),
        legend.position="bottom",legend.title = element_blank())

ggplot(dataset.p, aes(x=factor(date2), y=residuals1.st)) + 
  geom_point(alpha=0.7,col="#926b8d") +
  xlab("Months after RT") + ylab("Standardized residual - alr(Dense)")+
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = 
          element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = 
          element_line(colour = "black"),text = element_text(size=16),
        legend.position="bottom",legend.title = element_blank())

ggplot(dataset.p, aes(x=factor(date2), y=subres2.st)) + 
  geom_point(alpha=0.7,col="#926b8d") +
  xlab("Months after RT") + ylab("Standardized residual - alr(Hazy)")+
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = 
          element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = 
          element_line(colour = "black"),text = element_text(size=16),
        legend.position="bottom",legend.title = element_blank())

ggplot(dataset.p, aes(x=factor(date2), y=residuals2.st)) + 
  geom_point(alpha=0.7,col="#926b8d") +
  xlab("Months after RT") + ylab("Standardized residual - alr(Hazy)")+
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = 
          element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = 
          element_line(colour = "black"),text = element_text(size=16),
        legend.position="bottom",legend.title = element_blank())

#Figure 5.5
ggplot(dataset.p, aes(x=logp.1, y=subpred1)) + 
  geom_point(alpha=0.7,col="#926b8d") +
  xlab("Logit-probability of RILD") + ylab("Mean alr(Dense)")+
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = 
          element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = 
          element_line(colour = "black"),text = element_text(size=16),
        legend.position="bottom",legend.title = element_blank())

ggplot(dataset.p, aes(x=logp.1, y=subpred2)) + 
  geom_point(alpha=0.7,col="#926b8d") +
  xlab("Logit-probability of RILD") + ylab("Mean alr(Hazy)")+
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = 
          element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = 
          element_line(colour = "black"),text = element_text(size=16),
        legend.position="bottom",legend.title = element_blank())

#Figure 5.6
ggplot(agg2, aes(x=nonZeroHD, y=p.1)) + 
  geom_point(pch = 19, position = position_jitter(),alpha=0.3, 
             size=1,col="#B7B6BA") +
  geom_point(alpha=0.7,col="#926b8d") +
  xlab("Observed") + ylab("Predicted")+
  scale_x_continuous(breaks=c(0,1,2,3,4,5),labels=c(0,1,2,3,4,5)) +
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = 
          element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = 
          element_line(colour = "black"),text = element_text(size=16),
        legend.position="bottom",legend.title = element_blank())

#Figure 5.7
ggplot(dataset.p, aes(x=Dense.lt0, y=subpred1)) + 
  geom_point(alpha=0.7,col="#926b8d") +
  geom_abline(slope=1, intercept=0,col="#B7B6BA") +
  xlab("Observed alr(Dense)") + ylab("Sub-predicted alr(Dense)")+
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = 
          element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = 
          element_line(colour = "black"),text = element_text(size=16),
        legend.position="bottom",legend.title = element_blank())

ggplot(dataset.p, aes(x=Hazy.lt0, y=subpred2)) + 
  geom_point(alpha=0.7,col="#926b8d") +
  geom_abline(slope=1, intercept=0,col="#B7B6BA") +
  xlab("Observed alr(Hazy)") + ylab("Sub-predicted alr(Hazy)")+
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = 
          element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = 
          element_line(colour = "black"),text = element_text(size=16),
        legend.position="bottom",legend.title = element_blank())

ggplot(dataset.p, aes(x=Dense.lt0, y=pred1)) + 
  geom_point(alpha=0.7,col="#926b8d") +
  geom_abline(slope=1, intercept=0,col="#B7B6BA") +
  xlab("Observed alr(Dense)") + ylab("Predicted alr(Dense)")+
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = 
          element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = 
          element_line(colour = "black"),text = element_text(size=16),
        legend.position="bottom",legend.title = element_blank())

ggplot(dataset.p, aes(x=Hazy.lt0, y=pred2)) + 
  geom_point(alpha=0.7,col="#926b8d") +
  geom_abline(slope=1, intercept=0,col="#B7B6BA") +
  xlab("Observed alr(Hazy)") + ylab("Predicted alr(Hazy)")+
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = 
          element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = 
          element_line(colour = "black"),text = element_text(size=16),
        legend.position="bottom",legend.title = element_blank())

#Figure 5.8
v<-dataset.p
v$e.pred.f1<-exp(v$pred1)
v$e.pred.f2<-exp(v$pred2)
v$dense.pred.f<-v$e.pred.f1/(1+v$e.pred.f1+v$e.pred.f2)
v$hazy.pred.f<-v$e.pred.f2/(1+v$e.pred.f1+v$e.pred.f2)
v$com.pred.f<-1/(1+v$e.pred.f1+v$e.pred.f2)
v[,c("dense.pred.f","hazy.pred.f","com.pred.f")]<-clo(v[,
                                                        c("dense.pred.f","hazy.pred.f","com.pred.f")])

a<-c(3,6,12,18,24)
k<-a[1]
x<-acomp(v[v$date==k,c("dense.pred.f","hazy.pred.f","com.pred.f")])
xm<-(mean(x))
ggtern(data=v[v$date==k ,c("dense.pred.f","hazy.pred.f","com.pred.f")], 
       aes(x=dense.pred.f, y=hazy.pred.f, z=com.pred.f)) + 
  geom_point() +
  geom_point(aes(x=xm[1], y=xm[2], z=xm[3],colour="Mean")) +
  scale_colour_manual(name = element_blank(), labels = 
                        "Mean \ncomposition",values=c("#7aadb1"))+
  xlab("Dense") + ylab("Hazy") + zlab("Normal") +
  theme_showarrows() + ggtitle(paste(k, " Months", sep = "")) +
  theme(plot.title = element_text(face="bold", hjust = 0.5,vjust=8),
        plot.margin = unit(c(0,0,0,0), "cm"))

#Figure 5.9
ef.dose.data<-as.data.frame(rep(1:60,5))
colnames(ef.dose.data)<-c("Dose")
ef.dose.data$lDose<-log(ef.dose.data$Dose)
ef.dose.data$month<-rep(c(3,6,12,18,24),each=60)
ef.dose.data$time<-ef.dose.data$month/12
ef.dose.data$log<-par.l$Beta[[1]] + (ef.dose.data$lDose * 
                                       par.l$Beta[[2]]) + (ef.dose.data$time * par.l$Beta[[3]]) + 
  ((ef.dose.data$time^2) * par.l$Beta[[4]])
ef.dose.data$p.1<-exp(ef.dose.data$log)/(1+exp(ef.dose.data$log))
ef.dose.data$subpred1<-par.l$Gamma[[1]] + (ef.dose.data$lDose * 
                                             par.l$Gamma[[2]]) + (ef.dose.data$time * par.l$Gamma[[3]])+
  ((ef.dose.data$time^2) * par.l$Gamma[[4]])
ef.dose.data$subpred2<-par.l$Gamma[[5]] + (ef.dose.data$lDose * 
                                             par.l$Gamma[[6]]) + (ef.dose.data$time * par.l$Gamma[[7]]) +
  ((ef.dose.data$time^2) * par.l$Gamma[[8]])
ef.dose.data$e.subpred.f1<-exp(ef.dose.data$subpred1)*
  exp(0.5*par.l$Var.V1)
ef.dose.data$e.subpred.f2<-exp(ef.dose.data$subpred2)*
  exp(0.5*par.l$Var.V2)
ef.dose.data$e.pred.f1<-ef.dose.data$p.1*ef.dose.data$e.subpred.f1
ef.dose.data$e.pred.f2<-ef.dose.data$p.1*ef.dose.data$e.subpred.f2
ef.dose.data$dense.pred.f<-ef.dose.data$e.pred.f1/
  (1+ef.dose.data$e.pred.f1+ef.dose.data$e.pred.f2)
ef.dose.data$hazy.pred.f<-ef.dose.data$e.pred.f2/
  (1+ef.dose.data$e.pred.f1+ef.dose.data$e.pred.f2)
ef.dose.data$com.pred.f<-1/(1+ef.dose.data$e.pred.f1+
                              ef.dose.data$e.pred.f2)
ef.dose.data[,c("dense.pred.f","hazy.pred.f","com.pred.f")]<-
  clo(ef.dose.data[,c("dense.pred.f","hazy.pred.f","com.pred.f")])
ef.dose.data$t<-ef.dose.data$dense.pred.f+ef.dose.data$hazy.pred.f+
  ef.dose.data$com.pred.f

ggtern(data=ef.dose.data[ef.dose.data$month==24,], 
       aes(x=dense.pred.f, y=hazy.pred.f, z=com.pred.f, color=Dose)) + 
  geom_point()+
  geom_point(aes(x=b[1], y=b[2], z=b[3]),colour = "#7aadb1") +
  scale_color_gradient(breaks= c(10,20,30,40,50,60),labels= 
                         c(10,20,30,40,50,60),limits=c(0,60),
                       low = "#926b8d" , high = "#E8CA47") +
  labs(fill = "Median dose")  + xlab("Dense") + ylab("Hazy") + 
  zlab("Normal") +
  theme_showarrows() + ggtitle(paste(k, " Months", sep = "")) +
  theme(plot.title = element_text(face="bold", hjust = 0.5,vjust=8),
        plot.margin = unit(c(0,0,0,0), "cm")) 