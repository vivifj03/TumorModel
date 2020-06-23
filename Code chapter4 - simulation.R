
# Simulation ---------------------------------------------------------
library(MASS)

## Simulation base on spherical data
set.seed(124)
for (n in c(300,500,1000)) {
  for (ni in c(3,4,5)) {
    for (rp in 1:1000) {
      
      # Dose
      Dose<-rlnorm(n, meanlog = 0, sdlog =1)
      Dose<-Dose*(60/max(Dose))
      #plot(density(Dose))
      Dose <- rep(Dose,times=ni)
      
      id <- rep(1:n,times=ni)
      t <- c(3,6,12,18,24)/12
      t <- t[1:ni]
      time <- rep(t, each=n)
      data <- as.data.frame(cbind(id,Dose,time))
      data$time2 <- (data$time)^2
      data$int <- 1
      data <- data[,c("id","int","Dose","time","time2")]
      #dataU <- data[,c("id","int","Dose","time")]
      dataU<-data
      
      # Parameters
      m <- 0.6
      mu_1 <- c(-0.1,0.1)
      mu_2 <- -(m/(1-m))*mu_1
      Psi_diag<-diag(sqrt(c(1,0.5,1)))
      R <- matrix(c(1,0.1,0.2,0.1,1,0.1,0.2,0.1,1),nrow = 3)
      Psi<-Psi_diag %*% R %*% Psi_diag
      Psi_cc <- as.matrix(Psi[1,1])
      Psi_dd <- as.matrix(Psi[2:3,2:3])
      Psi_dc <- as.matrix(Psi[2:3,1])
      Psi_cd <- t(Psi_dc)
      #is.positive.definite(Psi_cc)
      #is.positive.definite(Psi_dd)
      #is.positive.definite(Psi)
      mu1 <- matrix(c(0,mu_1),ncol = 1)
      mu2 <- matrix(c(0,mu_2),ncol = 1)
      
      corr <- 0.2
      Var.V1 <- 1
      Var.V2 <- 0.5
      Sigma <- matrix(c(Var.V1,rep(corr*sqrt(Var.V1)*sqrt(Var.V2),2),
                        Var.V2),nrow = 2,ncol = 2)
      Beta <- as.matrix(c(-10,1.5,-2,-3),ncol=1)
      Gamma <- as.matrix(c(-7,1,0.5,-1,
                           -3,2,0.5,-1),ncol=1)
      b <- (m*mvrnorm(n = n, mu1, Psi))+((1-m)*mvrnorm(n = n, mu2,
            Psi))
      d <- data.frame(b[,2:3],i=rep(1:ni,ea=n))
      d <- d[,1:2]
      c <- data.frame(b[,1],i=rep(1:ni,ea=n))
      c <- c[,1]
      e <- mvrnorm(n=ni*n,c(0,0),Sigma) 
      
      X <- as.matrix(dataU[,-1])
      Z <- as.matrix(dataU[,c("int")])
      Zc <- (Z * c) 
      logit_U <- X %*% Beta + Zc
      pU <- exp(logit_U)/(1+exp(logit_U))
      U <-rbinom(length(pU),1,pU) 
      
      ps<-prop.table(table(U))[2]
      
      Zd1 <- Z * d[,1]
      Zd2 <- Z * d[,2]
      X <- as.matrix(data[,-1])
      Z <- as.matrix(data[,c("int")])
      V1 <- X %*% Gamma[1:4,] + Zd1 + e[,1]
      V2 <- X %*% Gamma[5:8,] + Zd2 + e[,2]
      
      data <- cbind(data,pU,U,V1,V2)
      data$V1u <- ifelse(data$U==0,0,data$V1) 
      data$V2u <- ifelse(data$U==0,0,data$V2) 
      data$nonZeroHD <- data$U
      
      
      sim.val<-c(Beta=Beta,Gamma=Gamma, Psi_cc=Psi_cc,Psi_dd=Psi_dd,
                 Var.V1=Var.V1,Var.V2=Var.V2,corr=corr,
                 Psi_dc=Psi_dc,m=m,mu_1=mu_1)
      
      #plot(data$Dose,data$V1u)
      #plot(data$Dose,data$V2u)
      #boxplot(Dose~nonZeroHD, data=data)
      
      #########################################################
      rm(list=setdiff(ls(), c("data","maxit","n","ni","rp","ps",
                              "sim.val")))
      
      dataset<-data[,-2]
      
      rm(data)
      names(dataset)<-c("patch","Dose","time","time2","pU","U","V1",
                        "V2","V1u","V2u","nonZeroHD")
      filename<-paste("s1/",n,"/",ni,"/data/","ds1","_",n,"_",ni,"_",
                      rp,".RData",sep = "")     
      save.image(file = filename) 
     
    }
  }
}


# Estimates for each data set ---------------------------------------
# Packages
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

n<- NUMBER #300, 500, 1000
rp<-REP #1-1000
prop<-PROP # 3,4,5

filename<-paste("s1/",n,"/",ni,"/data/","ds1","_",n,"_",prop,"_",
                rp,".RData",sep = "")
load(filename)

total.inicio<-Sys.time()
###### Mixed effects logistic regression -  gets initial values 
#This models is using time in years
data<-dataset
inicio.naive<-Sys.time()
md1 <- glmer(nonZeroHD ~ Dose + time + time2  + (1 | patch), 
             data = data, family = binomial,
             control = glmerControl(optimizer = "bobyqa",
             optCtrl=list(maxfun=2e8)),
             nAGQ = 10)

#summary(md1)

temp<-data
temp$V1u<- ifelse(temp$V1u==0,NA,temp$V1u)
temp$V2u<- ifelse(temp$V2u==0,NA,temp$V2u)
temp<- temp[is.na(temp$V1u)==F | is.na(temp$V1u)==F,]
md2d2 <- lmer(V1u ~ Dose + time + time2  + (1| patch), 
              data = temp, control = lmerControl(optCtrl=
              list(maxfun=2e10)))
md2d2.r<-ranef(md2d2)

md2h2 <- lmer(V2u ~ Dose + time + time2  + (1| patch), 
              data = temp, control = lmerControl(optCtrl=
              list(maxfun=2e10)))
md2h2.r<-ranef(md2h2)

# Predictions in naive model =====================================
coefp <- md1@beta

temp2<-unlist(unique(temp$patch))
temp2<-cbind(temp2,md2d2.r$patch[,1])
temp2<-cbind(temp2,md2h2.r$patch[,1])
temp2<-as.data.frame(temp2)
colnames(temp2)<-c("patch","d1","d2")
temp3<-setdiff(seq(1:n),temp2$patch)
temp3<-as.data.frame(temp3)
colnames(temp3)<-c("patch")
temp3$d1 <- mean(temp2$d1)
temp3$d2 <- mean(temp2$d2)
temp2<-rbind(temp2,temp3)
temp2<-temp2[order(temp2$patch),]
rm(temp3)

data<-merge(data,temp2,by="patch", all.x = T)
coefd <- md2d2@beta
coefh <- md2h2@beta


data$naive.p1 <- predict(md1, type = "response")
data$naive.d<- coefd[1] +  coefd[2]*data$Dose + coefd[3]*data$time + 
               coefd[4]*data$time2 + data$d1
data$naive.h<- coefh[1] +  coefh[2]*data$Dose + coefh[3]*data$time +
               coefh[4]*data$time2 + data$d2
data$naive.dr<- data$V1u -  data$naive.d
data$naive.hr<- data$V2u -  data$naive.h
data$V1u.e<-ifelse(data$V1u==0,0,exp(data$V1u))
data$V2u.e<-ifelse(data$V2u==0,0,exp(data$V2u))

data$naive.de<- data$naive.p1 * exp(data$naive.d) * 
                mean(exp(md2d2@resp$wtres))
data$naive.he<- data$naive.p1 * exp(data$naive.h) * 
                mean(exp(md2h2@resp$wtres))

data$naive.res.de<- data$V1u.e - data$naive.de
data$naive.res.he<- data$V2u.e - data$naive.he

naive.m_Vu.V1<-mean(data$V1u.e)
naive.TSS.V1<-sum((data$V1u.e-naive.m_Vu.V1)^2)
naive.sse.V.V1<-sum((data$V1u.e-data$naive.de)^2)
naive.ssr.V.V1<-sum((data$naive.de-naive.m_Vu.V1)^2)
naive.tss.V1<-naive.sse.V.V1+naive.ssr.V.V1
naive.corr2.V1<-cor(data$V1u.e,data$naive.de)^2 
naive.R2.V1<-1-(var(data$naive.res.de)/var(data$V1u.e))
naive.mse.V1<-mean_squared_error(data$V1u.e,data$naive.de)
naive.mae.V1<-mean_absolute_error(data$V1u.e,data$naive.de)
naive.mase.V1<-mean_absolute_scaled_error(data$V1u.e,data$naive.de)
naive.rmse.V1<-rmse(data$V1u.e,data$naive.de)
naive.tvs.V1<-total_variance_score(data$V1u.e,data$naive.de)
naive.uvs.V1<-unexplained_variance_score(data$V1u.e,data$naive.de)
naive.evs.V1<-explained_variance_score(data$V1u.e,data$naive.de) 
               #Best possible score is 1.0

naive.m_Vu.V2<-mean(data$V2u.e)
naive.TSS.V2<-sum((data$V2u.e-naive.m_Vu.V2)^2)
naive.sse.V.V2<-sum((data$V2u.e-data$naive.he)^2)
naive.ssr.V.V2<-sum((data$naive.he-naive.m_Vu.V2)^2)
naive.tss.V2<-naive.sse.V.V2+naive.ssr.V.V2
naive.corr2.V2<-cor(data$V2u.e,data$naive.he)^2 
naive.R2.V2<-1-(var(data$naive.res.he)/var(data$V2u.e))
naive.mse.V2<-mean_squared_error(data$V2u.e,data$naive.he)
naive.mae.V2<-mean_absolute_error(data$V2u.e,data$naive.he)
naive.mase.V2<-mean_absolute_scaled_error(data$V2u.e,data$naive.he)
naive.rmse.V2<-rmse(data$V2u.e,data$naive.he)
naive.tvs.V2<-total_variance_score(data$V2u.e,data$naive.he)
naive.uvs.V2<-unexplained_variance_score(data$V2u.e,data$naive.he)
naive.evs.V2<-explained_variance_score(data$V2u.e,data$naive.he) 

end.naive<-Sys.time()
end.naive-inicio.naive  

# Likelihood (tpmemmzl) ===========================================

#Initial values
#Ocurrence variable
Beta<- as.matrix(md1@beta,ncol=1)
if (((md1@theta)^2)<0.1) {
 #In case the initial value of Psi_cc was zero then put a small value
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
corr<-cor(temp$V1u, temp$V2u)
Var.V1<-attr(VarCorr(md2d2), "sc")^2
Var.V2<-attr(VarCorr(md2h2), "sc")^2  

Psi_dc<-matrix(c(sqrt(diag(Psi_dd))%*%sqrt(Psi_cc)*0.2),
               nrow=ncol(Psi_dd),ncol=ncol(Psi_cc)) #Corr between 
        #random efects between models
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


patch<-data$patch
patch<-unique(patch)

model.data<-list(datafile=data, id.var = patch,
                 X.var=c("Dose","time","time2"),Z.var=NULL,
                 X_a.var=c("Dose","time","time2"),Z_a.var=NULL,
                 U.var="nonZeroHD",V1="V1u",V2="V2u")

# Model ============================================================

#Functions 
source("tpmemmzl.R") #code in Chapter 3

# Deviance with the sim values 
sim<-par.l 
sim$Beta<-as.matrix(sim.val[1:4],ncol=1)
sim$Gamma<-as.matrix(sim.val[5:12],ncol=1)
sim$Psi_cc<-as.matrix(sim.val[13],ncol=1)
sim$Psi_dd<-matrix(c(sim.val[14:17]),nrow=2,ncol=2)
sim$Var.V1<-sim.val[18]
sim$Var.V2<-sim.val[19]
sim$corr<-sim.val[20]
sim$Psi_dc<-matrix(c(sim.val[21:22]),nrow=2,ncol=1)
sim$m<-sim.val[23]
sim$mu_1<-matrix(c(sim.val[24:25]),nrow=2,ncol=1)
in.dev<-lapply(patch,function(x) lklhd.dev(id=x,
               datafile=dataset,initial=sim,
               X.var=c("Dose","time","time2"),Z.var=NULL,
               X_a.var=c("Dose","time","time2"),Z_a.var=NULL,
               U.var="nonZeroHD",V1="V1u",V2="V2u"))
t.dev<-matrix(0,nrow=length(in.dev), ncol = 1)
for (i in 1:length(in.dev)) {
  t.dev[i,1]<-in.dev[[i]]
}
dev0 <- -2*colSums(t.dev)

### Estimation of fixed effects ####################################

inicio.coef<-Sys.time()
maxit=1000  #For empirical bayes estimates
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
         Psi_dc.i=par.l$Psi_dc %*% solve(par.l$Psi_cc),
         H.vech=vech(par.l$Psi_dd-(par.l$Psi_dc %*% 
                solve(par.l$Psi_cc) %*% t(par.l$Psi_dc)) ))
in.dev<-lapply(patch,function(x) lklhd.dev(id=x,datafile=dataset,
               initial=par.l,
               X.var=c("Dose","time","time2"),Z.var=NULL,
               X_a.var=c("Dose","time","time2"),Z_a.var=NULL,
               U.var="nonZeroHD",V1="V1u",V2="V2u"))
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
                       sqrt(par.l$Var.V1)*sqrt(par.l$Var.V2),2),par.l$Var.V2),nrow = 2,ncol = 2)),
           Psi_cc.vech=vech(par.l$Psi_cc),
           Psi_dc.i=par.l$Psi_dc %*% solve(par.l$Psi_cc),
           H.vech=vech(par.l$Psi_dd-(par.l$Psi_dc %*% 
                  solve(par.l$Psi_cc) %*% t(par.l$Psi_dc)) ))
  theta0<-theta
  
  t<-lapply(patch,function(x) lklhd.em(id=x,datafile=dataset,
                   initial=par.l,
                   X.var=c("Dose","time","time2"),Z.var=NULL,
                   X_a.var=c("Dose","time","time2"),Z_a.var=NULL,
                   U.var="nonZeroHD",V1="V1u",V2="V2u"))
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
  
  # Check H is positive define
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
                       X.var=c("Dose","time","time2"),Z.var=NULL,
                       X_a.var=c("Dose","time","time2"),Z_a.var=NULL,
                       U.var="nonZeroHD",V1="V1u",V2="V2u"))
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
    while (reldev[[step]] >-1e-8) {
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
                      X.var=c("Dose","time","time2"),Z.var=NULL,
                      X_a.var=c("Dose","time","time2"),Z_a.var=NULL,
                      U.var="nonZeroHD",V1="V1u",V2="V2u"))
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
  eu.dis[step] <-(dist(rbind(theta0,theta), method = 
                         "euclidean")^2)/length(theta)
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

### Standard errors ###############################################
se<-lapply(patch,function(x) lklhd.var(id=x,datafile=dataset,
                    est=par.l,
                    X.var=c("Dose","time","time2"),Z.var=NULL,
                    X_a.var=c("Dose","time","time2"),Z_a.var=NULL,
                    U.var="nonZeroHD",V1="V1u",V2="V2u"))

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

par.final<-c(unlist(par.l)) #Coefficients
bias<-par.final-sim.val #Bias
rel.bias<-bias*100/sim.val #Relative Bias
sqrt(diag(S_0b.total.i)) #SE for occurence model
sqrt(diag(S_0g.total.i)) #SE for intensity model
se.final1<-c(sqrt(diag(S_0b.total.i)),sqrt(diag(S_0g.total.i)))
lci1<-par.final[1:12]-(1.96*se.final1)
uci1<-par.final[1:12]+(1.96*se.final1)
coverage1<-rep(0,12)
for(e in 1:12) {
  coverage1[e]<-ifelse(sim.val[e]>lci1[e] & sim.val[e]<uci1[e],
                       1,coverage1[e])
}

### Random effects #################################################

try(if(criterion.step==1) stop("estimation did not converge"))

iniciop<-Sys.time()
random.est <- lapply(patch,function(x) lklhd.ram(id=x,
                        datafile=dataset,est=par.l,
                        X.var=c("Dose","time","time2"),Z.var=NULL,
                        X_a.var=c("Dose","time","time2"),Z_a.var=NULL,
                        U.var="nonZeroHD",V1="V1u",V2="V2u"))
endp<-Sys.time()
endp-iniciop

### Predicted values ################################################
#Predicted using basic EB RE
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
colnames(ran.est) <- c("patch","c.est","cc.t.est","d1.est",
                       "d2.est","ess","ess2")
dataset<-dataset[,1:11]
dataset$V1u.e<-ifelse(dataset$nonZeroHD==0,0,exp(dataset$V1))
dataset$V2u.e<-ifelse(dataset$nonZeroHD==0,0,exp(dataset$V2))
dataset<-merge(dataset,ran.est,by="patch",all.x = T)

dataset$logp.1 <-  par.l$Beta[1]  +  par.l$Beta[2] * dataset$Dose + 
  par.l$Beta[3] * dataset$time + par.l$Beta[4] * dataset$time2  +
  dataset$c.est
dataset$odds.1 <- exp(dataset$logp.1)
dataset$p.1 <- dataset$odds.1 / (1 + dataset$odds.1)

dataset$subpred1 <-  par.l$Gamma[1]  +  par.l$Gamma[2] * dataset$Dose + 
  par.l$Gamma[3] * dataset$time + par.l$Gamma[4] * dataset$time2  +
  dataset$d1.est
dataset$subpred2 <-  par.l$Gamma[5] + par.l$Gamma[6] * dataset$Dose + 
  par.l$Gamma[7] * dataset$time + par.l$Gamma[8] * dataset$time2 +
  dataset$d2.est

dataset$subres1 <-   dataset$V1u-dataset$subpred1
dataset$subres2 <-   dataset$V2u-dataset$subpred2  

duan1<- mean(exp(dataset[dataset$nonZeroHD==1,"subres1"]))
duan2<- mean(exp(dataset[dataset$nonZeroHD==1,"subres2"]))

dataset$pred1<- dataset$p.1*exp(dataset$subpred1)*duan1
dataset$pred2<- dataset$p.1*exp(dataset$subpred2)*duan2
dataset$residuals1<-(dataset$V1u.e-dataset$pred1)
dataset$residuals2<-(dataset$V2u.e-dataset$pred2)

ram.m_Vu.V1<-mean(dataset$V1u.e)
ram.TSS.V1<-sum((dataset$V1u.e-ram.m_Vu.V1)^2)
ram.sse.V.V1<-sum((dataset$V1u.e-dataset$pred1)^2)
ram.ssr.V.V1<-sum((dataset$pred1-ram.m_Vu.V1)^2)
ram.tss.V1<-ram.sse.V.V1+ram.ssr.V.V1
ram.corr2.V1<-cor(dataset$V1u.e,dataset$pred1)^2    
ram.R2.V1<-1-(var(dataset$residuals1)/var(dataset$V1u.e))
ram.mse.V1<-mean_squared_error(dataset$V1u.e,dataset$pred1)
ram.mae.V1<-mean_absolute_error(dataset$V1u.e,dataset$pred1)
ram.mase.V1<-mean_absolute_scaled_error(dataset$V1u.e,dataset$pred1)
ram.rmse.V1<-rmse(dataset$V1u.e,dataset$pred1)
ram.tvs.V1<-total_variance_score(dataset$V1u.e,dataset$pred1)
ram.uvs.V1<-unexplained_variance_score(dataset$V1u.e,dataset$pred1)
ram.evs.V1<-explained_variance_score(dataset$V1u.e,dataset$pred1) 

ram.m_Vu.V2<-mean(dataset$V2u.e)
ram.TSS.V2<-sum((dataset$V2u.e-ram.m_Vu.V2)^2)
ram.sse.V.V2<-sum((dataset$V2u.e-dataset$pred2)^2)
ram.ssr.V.V2<-sum((dataset$pred2-ram.m_Vu.V2)^2)
ram.tss.V2<-ram.sse.V.V2+ram.ssr.V.V2
ram.corr2.V2<-cor(dataset$V2u.e,dataset$pred2)^2     #https://wiki.bcs.rochester.edu/HlpLab/StatsCourses?action=AttachFile&do=get&target=Groningen11.pdf
ram.R2.V2<-1-(var(dataset$residuals2)/var(dataset$V2u.e))
ram.mse.V2<-mean_squared_error(dataset$V2u.e,dataset$pred2)
ram.mae.V2<-mean_absolute_error(dataset$V2u.e,dataset$pred2)
ram.mase.V2<-mean_absolute_scaled_error(dataset$V2u.e,dataset$pred2)
ram.rmse.V2<-rmse(dataset$V2u.e,dataset$pred2)
ram.tvs.V2<-total_variance_score(dataset$V2u.e,dataset$pred2)
ram.uvs.V2<-unexplained_variance_score(dataset$V2u.e,dataset$pred2)
ram.evs.V2<-explained_variance_score(dataset$V2u.e,dataset$pred2) #Best possible score is 1.0

#Predicted using truncated EB RE
dataset.t<-dataset[,1:13]
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
dataset.t$logp.1 <-  par.l$Beta[1] +par.l$Beta[2] * dataset.t$Dose + 
  par.l$Beta[3] * dataset.t$time + par.l$Beta[4] * dataset.t$time2  +
  dataset.t$c.est
dataset.t$odds.1 <- exp(dataset.t$logp.1)
dataset.t$p.1 <- dataset.t$odds.1 / (1 + dataset.t$odds.1)

dataset.t$subpred1 <- par.l$Gamma[1]+par.l$Gamma[2]*dataset.t$Dose + 
  par.l$Gamma[3] * dataset.t$time + par.l$Gamma[4]*dataset.t$time2  +
  dataset.t$d1.est
dataset.t$subpred2 <-  par.l$Gamma[5] +par.l$Gamma[6]*dataset.t$Dose + 
  par.l$Gamma[7] * dataset.t$time + par.l$Gamma[8] * dataset.t$time2 +
  dataset.t$d2.est

dataset.t$subres1 <-   dataset.t$V1u-dataset.t$subpred1
dataset.t$subres2 <-   dataset.t$V2u-dataset.t$subpred2  

duan1<- mean(exp(dataset.t[dataset.t$nonZeroHD==1,"subres1"]))
duan2<- mean(exp(dataset.t[dataset.t$nonZeroHD==1,"subres2"]))

dataset.t$pred1<- dataset.t$p.1*exp(dataset.t$subpred1)*duan1
dataset.t$pred2<- dataset.t$p.1*exp(dataset.t$subpred2)*duan2
dataset.t$residuals1<-(dataset.t$V1u.e-dataset.t$pred1)
dataset.t$residuals2<-(dataset.t$V2u.e-dataset.t$pred2)

ram.t.m_Vu.V1<-mean(dataset.t$V1u.e)
ram.t.TSS.V1<-sum((dataset.t$V1u.e-ram.t.m_Vu.V1)^2)
ram.t.sse.V.V1<-sum((dataset.t$V1u.e-dataset.t$pred1)^2)
ram.t.ssr.V.V1<-sum((dataset.t$pred1-ram.t.m_Vu.V1)^2)
ram.t.tss.V1<-ram.t.sse.V.V1+ram.t.ssr.V.V1
ram.t.corr2.V1<-cor(dataset.t$V1u.e,dataset.t$pred1)^2     
ram.t.mse.V1<-mean_squared_error(dataset.t$V1u.e,dataset.t$pred1)
ram.t.mae.V1<-mean_absolute_error(dataset.t$V1u.e,dataset.t$pred1)
ram.t.mase.V1<-mean_absolute_scaled_error(dataset.t$V1u.e,dataset.t$pred1)
ram.t.rmse.V1<-rmse(dataset.t$V1u.e,dataset.t$pred1)
ram.t.tvs.V1<-total_variance_score(dataset.t$V1u.e,dataset.t$pred1)
ram.t.uvs.V1<-unexplained_variance_score(dataset.t$V1u.e,dataset.t$pred1)
ram.t.evs.V1<-explained_variance_score(dataset.t$V1u.e,dataset.t$pred1)

ram.t.m_Vu.V2<-mean(dataset.t$V2u.e)
ram.t.TSS.V2<-sum((dataset.t$V2u.e-ram.t.m_Vu.V2)^2)
ram.t.sse.V.V2<-sum((dataset.t$V2u.e-dataset.t$pred2)^2)
ram.t.ssr.V.V2<-sum((dataset.t$pred2-ram.t.m_Vu.V2)^2)
ram.t.tss.V2<-ram.t.sse.V.V2+ram.t.ssr.V.V2
ram.t.corr2.V2<-cor(dataset.t$V2u.e,dataset.t$pred2)^2 
ram.t.R2.V2<-1-(var(dataset.t$residuals2)/var(dataset.t$V2u.e))
ram.t.mse.V2<-mean_squared_error(dataset.t$V2u.e,dataset.t$pred2)
ram.t.t.mae.V2<-mean_absolute_error(dataset.t$V2u.e,dataset.t$pred2)
ram.t.mase.V2<-mean_absolute_scaled_error(dataset.t$V2u.e,dataset.t$pred2)
ram.t.rmse.V2<-rmse(dataset.t$V2u.e,dataset.t$pred2)
ram.t.tvs.V2<-total_variance_score(dataset.t$V2u.e,dataset.t$pred2)
ram.t.uvs.V2<-unexplained_variance_score(dataset.t$V2u.e,dataset.t$pred2)
ram.t.evs.V2<-explained_variance_score(dataset.t$V2u.e,dataset.t$pred2) 

#Predicted using pareto EB RE
dataset.p<-dataset[,1:13]
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
dataset.p<-merge(dataset.p,ran.est,by="patch",all.x = T)

dataset.p$logp.1 <-  par.l$Beta[1]+ par.l$Beta[2]*dataset.p$Dose + 
  par.l$Beta[3] * dataset.p$time + par.l$Beta[4]*dataset.p$time2  +
  dataset.p$c.est
dataset.p$odds.1 <- exp(dataset.p$logp.1)
dataset.p$p.1 <- dataset.p$odds.1 / (1 + dataset.p$odds.1)

dataset.p$subpred1 <-  par.l$Gamma[1]+par.l$Gamma[2]*dataset.p$Dose + 
  par.l$Gamma[3] * dataset.p$time + par.l$Gamma[4]*dataset.p$time2  +
  dataset.p$d1.est
dataset.p$subpred2 <-  par.l$Gamma[5]+ par.l$Gamma[6]*dataset.p$Dose + 
  par.l$Gamma[7] * dataset.p$time + par.l$Gamma[8] * dataset.p$time2 +
  dataset.p$d2.est

dataset.p$subres1 <-   dataset.p$V1u-dataset.p$subpred1
dataset.p$subres2 <-   dataset.p$V2u-dataset.p$subpred2  

duan1<- mean(exp(dataset.p[dataset.p$nonZeroHD==1,"subres1"]))
duan2<- mean(exp(dataset.p[dataset.p$nonZeroHD==1,"subres2"]))

dataset.p$pred1<- dataset.p$p.1*exp(dataset.p$subpred1)*duan1
dataset.p$pred2<- dataset.p$p.1*exp(dataset.p$subpred2)*duan2
dataset.p$residuals1<-(dataset.p$V1u.e-dataset.p$pred1)
dataset.p$residuals2<-(dataset.p$V2u.e-dataset.p$pred2)

ram.p.m_Vu.V1<-mean(dataset.p$V1u.e)
ram.p.pSS.V1<-sum((dataset.p$V1u.e-ram.p.m_Vu.V1)^2)
ram.p.sse.V.V1<-sum((dataset.p$V1u.e-dataset.p$pred1)^2)
ram.p.ssr.V.V1<-sum((dataset.p$pred1-ram.p.m_Vu.V1)^2)
ram.p.pss.V1<-ram.p.sse.V.V1+ram.p.ssr.V.V1
ram.p.corr2.V1<-cor(dataset.p$V1u.e,dataset.p$pred1)^2 
ram.p.R2.V1<-1-(var(dataset.p$residuals1)/var(dataset.p$V1u.e))
ram.p.mse.V1<-mean_squared_error(dataset.p$V1u.e,dataset.p$pred1)
ram.p.mae.V1<-mean_absolute_error(dataset.p$V1u.e,dataset.p$pred1)
ram.p.mase.V1<-mean_absolute_scaled_error(dataset.p$V1u.e,dataset.p$pred1)
ram.p.rmse.V1<-rmse(dataset.p$V1u.e,dataset.p$pred1)
ram.p.pvs.V1<-total_variance_score(dataset.p$V1u.e,dataset.p$pred1)
ram.p.uvs.V1<-unexplained_variance_score(dataset.p$V1u.e,dataset.p$pred1)
ram.p.evs.V1<-explained_variance_score(dataset.p$V1u.e,dataset.p$pred1)

ram.p.m_Vu.V2<-mean(dataset.p$V2u.e)
ram.p.pSS.V2<-sum((dataset.p$V2u.e-ram.p.m_Vu.V2)^2)
ram.p.sse.V.V2<-sum((dataset.p$V2u.e-dataset.p$pred2)^2)
ram.p.ssr.V.V2<-sum((dataset.p$pred2-ram.p.m_Vu.V2)^2)
ram.p.pss.V2<-ram.p.sse.V.V2+ram.p.ssr.V.V2
ram.p.corr2.V2<-cor(dataset.p$V2u.e,dataset.p$pred2)^2 
ram.p.R2.V2<-1-(var(dataset.p$residuals2)/var(dataset.p$V2u.e))
ram.p.mse.V2<-mean_squared_error(dataset.p$V2u.e,dataset.p$pred2)
ram.p.p.mae.V2<-mean_absolute_error(dataset.p$V2u.e,dataset.p$pred2)
ram.p.mase.V2<-mean_absolute_scaled_error(dataset.p$V2u.e,dataset.p$pred2)
ram.p.rmse.V2<-rmse(dataset.p$V2u.e,dataset.p$pred2)
ram.p.pvs.V2<-total_variance_score(dataset.p$V2u.e,dataset.p$pred2)
ram.p.uvs.V2<-unexplained_variance_score(dataset.p$V2u.e,dataset.p$pred2)
ram.p.evs.V2<-explained_variance_score(dataset.p$V2u.e,dataset.p$pred2)

save.image(file=filenamep)
total.end<-Sys.time()
total.end-total.inicio

# Combine Simulation Results =======================================  
files1<-list.files(path = "s1/300/3/pred/",pattern = "*.RData")
files2<-list.files(path = "s1/300/4/pred/",pattern = "*.RData")
files3<-list.files(path = "s1/300/5/pred/",pattern = "*.RData")
files4<-list.files(path = "s1/500/3/pred/",pattern = "*.RData")
files5<-list.files(path = "s1/500/4/pred/",pattern = "*.RData")
files6<-list.files(path = "s1/500/5/pred/",pattern = "*.RData")
files7<-list.files(path = "s1/1000/3/pred/",pattern = "*.RData")
files8<-list.files(path = "s1/1000/4/pred/",pattern = "*.RData")
files9<-list.files(path = "s1/1000/5/pred/",pattern = "*.RData")
files<-c(files1,files2,files3,files4,files5,files6,
         files7,files8,files9)
file_name<-"syymmdd_combine_sim.RData"

rp<-c(1:1000)
NN<- c(300,500,1000) 
prop <- c(3,4,5)

tot<-length(rp)*length(NN)*length(prop)
results<-as.data.frame(matrix(0,nrow =tot,ncol = 280))
line<-1

for (i in 1:3) {
  for (k in 1:1000) {
    for (g in 1:3) {
      results[line,1] <- NN[i]
      results[line,2] <- prop[g]
      results[line,3] <- rp[k]
      results[line,4]<-paste("s1/",NN[i],"/",prop[g],"/pred24c/",
                             "s1","_",NN[i],"_",prop[g],"_",rp[k],
                             ".RData",sep = "")
      results[line,5]<-paste("s1","_",NN[i],"_",prop[g],"_",rp[k],
                             ".RData",sep = "")
      line<-line+1
    }
  }  
}

save.image(file=file_name)

(tot2<-nrow(results))

results$rdata<-results[,5] %in% files
resultsold<-results
results<-results[results$rdata==TRUE,]

(tot2<-nrow(results))

naive <- as.data.frame(matrix(0,nrow =tot2  ,ncol = 57))
naive[,1:5] <- results[,1:5]

V1 <- as.data.frame(matrix(0,nrow =tot2  ,ncol = 19))
V1[,1:5] <- results[,1:5]
V2 <- V1
V1_t <- V1
V2_t <- V1
V1_p <- V1
V2_p <- V1

# Predicted 
for (l in 1:tot2) {
  load(file = results[l,4]) 
  
  naive[l,6:9] <- coefp
  naive[l,10:13] <- coefd
  naive[l,14:17] <- coefh
  naive[l,18:21] <- sqrt(diag(vcov(md1)))
  naive[l,22:25] <- sqrt(diag(vcov(md2d2)))
  naive[l,26:29] <- sqrt(diag(vcov(md2h2)))
  naive[l,30] <- naive.m_Vu.V1
  naive[l,31] <- naive.TSS.V1
  naive[l,32] <- naive.sse.V.V1
  naive[l,33] <- naive.ssr.V.V1
  naive[l,34] <- naive.tss.V1
  naive[l,35] <- naive.corr2.V1
  naive[l,36] <- naive.R2.V1
  naive[l,37] <- naive.mse.V1
  naive[l,38] <- naive.mae.V1
  naive[l,39] <- naive.mase.V1
  naive[l,40] <- naive.rmse.V1
  naive[l,41] <- naive.tvs.V1
  naive[l,42] <- naive.uvs.V1
  naive[l,43] <- naive.evs.V1
  naive[l,44] <- naive.m_Vu.V2
  naive[l,45] <- naive.TSS.V2
  naive[l,46] <- naive.sse.V.V2
  naive[l,47] <- naive.ssr.V.V2
  naive[l,48] <- naive.tss.V2
  naive[l,49] <- naive.corr2.V2
  naive[l,50] <- naive.R2.V2
  naive[l,51] <- naive.mse.V2
  naive[l,52] <- naive.mae.V2
  naive[l,53] <- naive.mase.V2
  naive[l,54] <- naive.rmse.V2
  naive[l,55] <- naive.tvs.V2
  naive[l,56] <- naive.uvs.V2
  naive[l,57] <- naive.evs.V2
  
  results[l,6:30]<-par.final
  results[l,31:42]<-se.final1
  #results[l,43:54]<-se.final2
  results[l,55:79]<-bias
  results[l,80:104]<-rel.bias
  results[l,105:116]<-coverage1
  #results[l,117:128]<-coverage2
  initial.val<-unlist(estimate[[1]])
  results[l,129:153] <- initial.val
  results[l,154:178]<-par.final-initial.val
  results[l,179:203]<-(par.final-initial.val)*100/initial.val
  coverage3<-rep(0,12)
  for(e in 1:12) {
    coverage3[e]<-ifelse(initial.val[e]>lci1[e] & initial.val[e]<uci1[e],1,coverage3[e])
  }    
  results[l,204:215]<-coverage3
  results[l,216]<-step
  results[l,217]<-criterion
  results[l,218]<-criterion.distance
  results[l,219]<-criterion.max
  results[l,220]<-criterion.step
  results[l,221]<-criterion.dev
  results[l,222]<-criterion.reldev
  results[l,223]<-ii
  results[l,224]<-end.coef-inicio.coef
  results[l,225]<-ps
  results[l,226]<-max(abs(dif.theta[[step]]))
  results[l,227]<-max(abs(rel.dif.theta[[step]]))
  results[l,228]<-absdev
  results[l,229]<-reldev[[step]]
  results[l,230:241]<-lci1
  results[l,242:253]<-uci1
  results[l,254:265]<-lci2
  results[l,266:277]<-uci2
  results[l,278] <- dev[[step]]
  results[l,279] <- dev0
  results[l,280] <- loop3
  
  V1[l,6]<-ram.m_Vu.V1
  V1[l,7]<-ram.TSS.V1
  V1[l,8]<-ram.sse.V.V1
  V1[l,9]<-ram.ssr.V.V1
  V1[l,10]<-ram.tss.V1
  V1[l,11]<-ram.corr2.V1
  V1[l,12]<-ram.R2.V1
  V1[l,13]<-ram.mse.V1
  V1[l,14]<-ram.mae.V1
  V1[l,15]<-ram.mase.V1
  V1[l,16]<-ram.rmse.V1
  V1[l,17]<-ram.tvs.V1
  V1[l,18]<-ram.uvs.V1
  V1[l,19]<-ram.evs.V1
  
  V2[l,6]<-ram.m_Vu.V2
  V2[l,7]<-ram.TSS.V2
  V2[l,8]<-ram.sse.V.V2
  V2[l,9]<-ram.ssr.V.V2
  V2[l,10]<-ram.tss.V2
  V2[l,11]<-ram.corr2.V2
  V2[l,12]<-ram.R2.V2
  V2[l,13]<-ram.mse.V2
  V2[l,14]<-ram.mae.V2
  V2[l,15]<-ram.mase.V2
  V2[l,16]<-ram.rmse.V2
  V2[l,17]<-ram.tvs.V2
  V2[l,18]<-ram.uvs.V2
  V2[l,19]<-ram.evs.V2  
  
  V1_t[l,6]<-ram.t.m_Vu.V1
  V1_t[l,7]<-ram.t.TSS.V1
  V1_t[l,8]<-ram.t.sse.V.V1
  V1_t[l,9]<-ram.t.ssr.V.V1
  V1_t[l,10]<-ram.t.tss.V1
  V1_t[l,11]<-ram.t.corr2.V1
  V1_t[l,12]<-ram.t.R2.V1
  V1_t[l,13]<-ram.t.mse.V1
  V1_t[l,14]<-ram.t.mae.V1
  V1_t[l,15]<-ram.t.mase.V1
  V1_t[l,16]<-ram.t.rmse.V1
  V1_t[l,17]<-ram.t.tvs.V1
  V1_t[l,18]<-ram.t.uvs.V1
  V1_t[l,19]<-ram.t.evs.V1
  
  V2_t[l,6]<-ram.t.m_Vu.V2
  V2_t[l,7]<-ram.t.TSS.V2
  V2_t[l,8]<-ram.t.sse.V.V2
  V2_t[l,9]<-ram.t.ssr.V.V2
  V2_t[l,10]<-ram.t.tss.V2
  V2_t[l,11]<-ram.t.corr2.V2
  V2_t[l,12]<-ram.t.R2.V2
  V2_t[l,13]<-ram.t.mse.V2
  V2_t[l,14]<-ram.t.t.mae.V2
  V2_t[l,15]<-ram.t.mase.V2
  V2_t[l,16]<-ram.t.rmse.V2
  V2_t[l,17]<-ram.t.tvs.V2
  V2_t[l,18]<-ram.t.uvs.V2
  V2_t[l,19]<-ram.t.evs.V2  
  
  V1_p[l,6]<-ram.p.m_Vu.V1
  V1_p[l,7]<-ram.p.pSS.V1
  V1_p[l,8]<-ram.p.sse.V.V1
  V1_p[l,9]<-ram.p.ssr.V.V1
  V1_p[l,10]<-ram.p.pss.V1
  V1_p[l,11]<-ram.p.corr2.V1
  V1_p[l,12]<-ram.p.R2.V1
  V1_p[l,13]<-ram.p.mse.V1
  V1_p[l,14]<-ram.p.mae.V1
  V1_p[l,15]<-ram.p.mase.V1
  V1_p[l,16]<-ram.p.rmse.V1
  V1_p[l,17]<-ram.p.pvs.V1
  V1_p[l,18]<-ram.p.uvs.V1
  V1_p[l,19]<-ram.p.evs.V1
  
  V2_p[l,6]<-ram.p.m_Vu.V2
  V2_p[l,7]<-ram.p.pSS.V2
  V2_p[l,8]<-ram.p.sse.V.V2
  V2_p[l,9]<-ram.p.ssr.V.V2
  V2_p[l,10]<-ram.p.pss.V2
  V2_p[l,11]<-ram.p.corr2.V2
  V2_p[l,12]<-ram.p.R2.V2
  V2_p[l,13]<-ram.p.mse.V2
  V2_p[l,14]<-ram.p.p.mae.V2
  V2_p[l,15]<-ram.p.mase.V2
  V2_p[l,16]<-ram.p.rmse.V2
  V2_p[l,17]<-ram.p.pvs.V2
  V2_p[l,18]<-ram.p.uvs.V2
  V2_p[l,19]<-ram.p.evs.V2  
}

colnames(naive) <- c("N","ni","rep","file","file2","naive.beta1",
 "naive.beta2","naive.beta3","naive.beta4", "naive.gamma1",
 "naive.gamma2","naive.gamma3","naive.gamma4","naive.gamma5",
 "naive.gamma6","naive.gamma7","naive.gamma8","naive.se.Beta1",
 "naive.se.Beta2","naive.se.Beta3","naive.se.Beta4",
 "naive.se.Gamma1","naive.se.Gamma2","naive.se.Gamma3",
 "naive.se.Gamma4","naive.se.Gamma5","naive.se.Gamma6",
 "naive.se.Gamma7","naive.se.Gamma8","naive.m_Vu.V1","naive.TSS.V1",
 "naive.sse.V.V1","naive.ssr.V.V1","naive.tss.V1","naive.corr2.V1",
 "naive.R2.V1","naive.mse.V1","naive.mae.V1","naive.mase.V1",
 "naive.rmse.V1","naive.tvs.V1","naive.uvs.V1","naive.evs.V1",
 "naive.m_Vu.V2","naive.TSS.V2","naive.sse.V.V2","naive.ssr.V.V2",
 "naive.tss.V2","naive.corr2.V2","naive.R2.V2","naive.mse.V2",
 "naive.mae.V2","naive.mase.V2","naive.rmse.V2","naive.tvs.V2",
 "naive.uvs.V2","naive.evs.V2")

colnames(resultsold)<-c("N","ni","rep","file","file2",
  "beta1","beta2","beta3","beta4","gamma1","gamma2","gamma3",
  "gamma4","gamma5","gamma6","gamma7","gamma8","Psi_cc", "Psi_dd1",
  "Psi_dd2","Psi_dd3","Psi_dd4","Var.V1","Var.V2","corr","Psi_dc1",
  "Psi_dc2","m","mu_1_1","mu_1_2","se.Beta1","se.Beta2","se.Beta3",
  "se.Beta4","se.Gamma1","se.Gamma2","se.Gamma3","se.Gamma4",
  "se.Gamma5","se.Gamma6","se.Gamma7","se.Gamma8","se2.Beta1",
  "se2.Beta2","se2.Beta3","se2.Beta4","se2.Gamma1","se2.Gamma2",
  "se2.Gamma3","se2.Gamma4","se2.Gamma5","se2.Gamma6","se2.Gamma7",
  "se2.Gamma8","b.Beta1","b.Beta2","b.Beta3","b.Beta4","b.Gamma1",
  "b.Gamma2","b.Gamma3","b.Gamma4","b.Gamma5","b.Gamma6","b.Gamma7",
  "b.Gamma8","b.Psi_cc","b.Psi_dd1","b.Psi_dd2","b.Psi_dd3",
  "b.Psi_dd4", "b.Var.V1","b.Var.V2","b.corr","b.Psi_dc1","b.Psi_dc2",
  "b.m","b.mu_11","b.mu_12","rb.Beta1","rb.Beta2","rb.Beta3",
  "rb.Beta4","rb.Gamma1","rb.Gamma2","rb.Gamma3","rb.Gamma4",
  "rb.Gamma5","rb.Gamma6","rb.Gamma7","rb.Gamma8","rb.Psi_cc",
  "rb.Psi_dd1","rb.Psi_dd2","rb.Psi_dd3","rb.Psi_dd4", "rb.Var.V1",
  "rb.Var.V2","rb.corr","rb.Psi_dc1","rb.Psi_dc2","rb.m","rb.mu_11",
  "rb.mu_12","cv.Beta1","cv.Beta2","cv.Beta3","cv.Beta4","cv.Gamma1",
  "cv.Gamma2","cv.Gamma3","cv.Gamma4","cv.Gamma5","cv.Gamma6",
  "cv.Gamma7","cv.Gamma8","cv2.Beta1","cv2.Beta2","cv2.Beta3",
  "cv2.Beta4","cv2.Gamma1","cv2.Gamma2","cv2.Gamma3","cv2.Gamma4",
  "cv2.Gamma5","cv2.Gamma6","cv2.Gamma7","cv2.Gamma8","beta1.i",
  "beta2.i","beta3.i","beta4.i","gamma1.i","gamma2.i","gamma3.i",
  "gamma4.i","gamma5.i","gamma6.i","gamma7.i","gamma8.i","Psi_cc.i",
  "Psi_dd1.i","Psi_dd2.i","Psi_dd3.i","Psi_dd4.i","Var.V1.i",
  "Var.V2.i","corr.i","Psi_dc1.i","Psi_dc2.i","m.i","mu_1_1.i",
  "mu_1_2.i","b.Beta1.i","b.Beta2.i","b.Beta3.i","b.Beta4.i",
  "b.Gamma1.i","b.Gamma2.i","b.Gamma3.i","b.Gamma4.i","b.Gamma5.i",
  "b.Gamma6.i","b.Gamma7.i","b.Gamma8.i","b.Psi_cc.i","b.Psi_dd1.i",
  "b.Psi_dd2.i","b.Psi_dd3.i","b.Psi_dd4.i", "b.Var.V1.i",
  "b.Var.V2.i","b.corr.i","b.Psi_dc1.i","b.Psi_dc2.i","b.m.i",
  "b.mu_11.i","b.mu_12.i","rb.Beta1.i","rb.Beta2.i","rb.Beta3.i",
  "rb.Beta4.i","rb.Gamma1.i","rb.Gamma2.i","rb.Gamma3.i",
  "rb.Gamma4.i","rb.Gamma5.i","rb.Gamma6.i","rb.Gamma7.i",
  "rb.Gamma8.i","rb.Psi_cc.i","rb.Psi_dd1.i","rb.Psi_dd2.i",
  "rb.Psi_dd3.i","rb.Psi_dd4.i","rb.Var.V1.i","rb.Var.V2.i",
  "rb.corr.i","rb.Psi_dc1.i","rb.Psi_dc2.i","rb.m.i","rb.mu_11.i",
  "rb.mu_12.i","cv.Beta1.i","cv.Beta2.i","cv.Beta3.i","cv.Beta4.i",
  "cv.Gamma1.i","cv.Gamma2.i","cv.Gamma3.i","cv.Gamma4.i",
  "cv.Gamma5.i","cv.Gamma6.i","cv.Gamma7.i","cv.Gamma8.i",
  "step","criterion","criterion.distance","criterion.max",
  "criterion.step","criterion.dev","criterion.reldev","ii",
  "time","ps","max.abs.dif.theta","max.abs.rel.dif.theta","absdev",
  "reldev","lci1.beta1","lci1.beta2","lci1.beta3","lci1.beta4",
  "lci1.gamma1","lci1.gamma2","lci1.gamma3","lci1.gamma4",
  "lci1.gamma5","lci1.gamma6","lci1.gamma7","lci1.gamma8",
  "uci1.beta1","uci1.beta2","uci1.beta3","uci1.beta4", 
  "uci1.gamma1","uci1.gamma2","uci1.gamma3","uci1.gamma4",
  "uci1.gamma5","uci1.gamma6","uci1.gamma7","uci1.gamma8",
  "lci2.beta1","lci2.beta2","lci2.beta3","lci2.beta4", 
  "lci2.gamma1","lci2.gamma2","lci2.gamma3","lci2.gamma4",
  "lci2.gamma5","lci2.gamma6","lci2.gamma7","lci2.gamma8",
  "uci2.beta1","uci2.beta2","uci2.beta3","uci2.beta4", 
  "uci2.gamma1","uci2.gamma2","uci2.gamma3","uci2.gamma4",
  "uci2.gamma5","uci2.gamma6","uci2.gamma7","uci2.gamma8",
  "dev","dev0","loop3") 
colnames(results)<-colnames(resultsold)
results$conv<-ifelse(results$criterion.distance==1 | 
                     results$criterion.max==1 |
                     results$criterion.dev==1 | 
                     results$criterion.reldev==1,1,0)
results$step<-results$step-1
results$time<-results$time/60

colnames(V1)<-c("N","ni","rep","file","file2","m_Vu","TSS","sse.V",
                "ssr.V","tss","corr2","R2","mse","mae","mase",
                "rmse","tvs","uvs","evs")

colnames(V2) <- colnames(V1)
colnames(V1_t) <- colnames(V1)
colnames(V2_t) <- colnames(V1)
colnames(V1_p) <- colnames(V1)
colnames(V2_p) <- colnames(V1)

save.image(file=file_name)


# Summary ======================================= 

# This code uses the compile simulation 
# (scenario|sample size|RILD proportion)
# and produce the summary

#Data
load(file=file_name)

#First, delete objects we do not need
rm(list=setdiff(ls(), c("results", "resultsold", "naive", "V1", "V2",
                       "V1.m", "V2.m", "V1_t", "V2_t", "V1.m_t", 
                       "V2.m_t", "V1_p", "V2_p", "V1.m_p", "V2.m_p")))

# We need to trimm the outliers. 
# We did 1% (0.5% at each end) 
# Next lines, find the cutoff for trimming
resultsorig<-results
results.p05<-aggregate(results[,c(1:2,6:30)], 
             by = list(results$N,results$ni), 
             function(x) quantile(x, c(.005), na.rm=T))
results.p995<-aggregate(results[,c(1:2,6:30)], 
              by = list(results$N,results$ni), 
              function(x) quantile(x, c(.995), na.rm=T))
# Next lines are very inefficient and take long time but they 
# check if the simulation values are part of the outliers (1%)
results$beta1trim<-0
results$beta2trim<-0
results$beta3trim<-0
results$beta4trim<-0
results$gamma1trim<-0
results$gamma2trim<-0
results$gamma3trim<-0
results$gamma4trim<-0
results$gamma5trim<-0
results$gamma6trim<-0
results$gamma7trim<-0
results$gamma8trim<-0
results$Psi_cctrim<-0
results$Psi_dd1trim<-0
results$Psi_dd2trim<-0
results$Psi_dd3trim<-0
results$Psi_dd4trim<-0
results$Var.V1trim<-0
results$Var.V2trim<-0
results$corrtrim<-0
results$Psi_dc1trim<-0
results$Psi_dc2trim<-0
results$mtrim<-0
results$mu_1_1trim<-0
results$mu_1_2trim<-0

for(j in 1:nrow(results)) {
  N<-results[j,]$N
  ni<-results[j,]$ni
  beta1a<-results.p05[results.p05$N==N & results.p05$ni==ni,"beta1"]
  beta2a<-results.p05[results.p05$N==N & results.p05$ni==ni,"beta2"]
  beta3a<-results.p05[results.p05$N==N & results.p05$ni==ni,"beta3"]
  beta4a<-results.p05[results.p05$N==N & results.p05$ni==ni,"beta4"]
  gamma1a<-results.p05[results.p05$N==N & results.p05$ni==ni,"gamma1"]
  gamma2a<-results.p05[results.p05$N==N & results.p05$ni==ni,"gamma2"]
  gamma3a<-results.p05[results.p05$N==N & results.p05$ni==ni,"gamma3"]
  gamma4a<-results.p05[results.p05$N==N & results.p05$ni==ni,"gamma4"]
  gamma5a<-results.p05[results.p05$N==N & results.p05$ni==ni,"gamma5"]
  gamma6a<-results.p05[results.p05$N==N & results.p05$ni==ni,"gamma6"]
  gamma7a<-results.p05[results.p05$N==N & results.p05$ni==ni,"gamma7"]
  gamma8a<-results.p05[results.p05$N==N & results.p05$ni==ni,"gamma8"]
  Psi_cca<-results.p05[results.p05$N==N & results.p05$ni==ni,"Psi_cc"]
  Psi_dd1a<-results.p05[results.p05$N==N & results.p05$ni==ni,"Psi_dd1"]
  Psi_dd2a<-results.p05[results.p05$N==N & results.p05$ni==ni,"Psi_dd2"]
  Psi_dd3a<-results.p05[results.p05$N==N & results.p05$ni==ni,"Psi_dd3"]
  Psi_dd4a<-results.p05[results.p05$N==N & results.p05$ni==ni,"Psi_dd4"]
  Var.V1a<-results.p05[results.p05$N==N & results.p05$ni==ni,"Var.V1"]
  Var.V2a<-results.p05[results.p05$N==N & results.p05$ni==ni,"Var.V2"]
  corra<-results.p05[results.p05$N==N & results.p05$ni==ni,"corr"]
  Psi_dc1a<-results.p05[results.p05$N==N & results.p05$ni==ni,"Psi_dc1"]
  Psi_dc2a<-results.p05[results.p05$N==N & results.p05$ni==ni,"Psi_dc2"]
  ma<-results.p05[results.p05$N==N & results.p05$ni==ni,"m"]
  mu_1_1a<-results.p05[results.p05$N==N & results.p05$ni==ni,"mu_1_1"]
  mu_1_2a<-results.p05[results.p05$N==N & results.p05$ni==ni,"mu_1_2"]
  
  beta1b<-results.p995[results.p995$N==N & results.p995$ni==ni,"beta1"]
  beta2b<-results.p995[results.p995$N==N & results.p995$ni==ni,"beta2"]
  beta3b<-results.p995[results.p995$N==N & results.p995$ni==ni,"beta3"]
  beta4b<-results.p995[results.p995$N==N & results.p995$ni==ni,"beta4"]
  gamma1b<-results.p995[results.p995$N==N & results.p995$ni==ni,"gamma1"]
  gamma2b<-results.p995[results.p995$N==N & results.p995$ni==ni,"gamma2"]
  gamma3b<-results.p995[results.p995$N==N & results.p995$ni==ni,"gamma3"]
  gamma4b<-results.p995[results.p995$N==N & results.p995$ni==ni,"gamma4"]
  gamma5b<-results.p995[results.p995$N==N & results.p995$ni==ni,"gamma5"]
  gamma6b<-results.p995[results.p995$N==N & results.p995$ni==ni,"gamma6"]
  gamma7b<-results.p995[results.p995$N==N & results.p995$ni==ni,"gamma7"]
  gamma8b<-results.p995[results.p995$N==N & results.p995$ni==ni,"gamma8"]
  Psi_ccb<-results.p995[results.p995$N==N & results.p995$ni==ni,"Psi_cc"]
  Psi_dd1b<-results.p995[results.p995$N==N & results.p995$ni==ni,"Psi_dd1"]
  Psi_dd2b<-results.p995[results.p995$N==N & results.p995$ni==ni,"Psi_dd2"]
  Psi_dd3b<-results.p995[results.p995$N==N & results.p995$ni==ni,"Psi_dd3"]
  Psi_dd4b<-results.p995[results.p995$N==N & results.p995$ni==ni,"Psi_dd4"]
  Var.V1b<-results.p995[results.p995$N==N & results.p995$ni==ni,"Var.V1"]
  Var.V2b<-results.p995[results.p995$N==N & results.p995$ni==ni,"Var.V2"]
  corrb<-results.p995[results.p995$N==N & results.p995$ni==ni,"corr"]
  Psi_dc1b<-results.p995[results.p995$N==N & results.p995$ni==ni,"Psi_dc1"]
  Psi_dc2b<-results.p995[results.p995$N==N & results.p995$ni==ni,"Psi_dc2"]
  mb<-results.p995[results.p995$N==N & results.p995$ni==ni,"m"]
  mu_1_1b<-results.p995[results.p995$N==N & results.p995$ni==ni,"mu_1_1"]
  mu_1_2b<-results.p995[results.p995$N==N & results.p995$ni==ni,"mu_1_2"]
  
  results[j,]$beta1trim<-ifelse(results[j,]$beta1<beta1a | results[j,]$beta1>beta1b,1,0)
  results[j,]$beta2trim<-ifelse(results[j,]$beta2<beta2a | results[j,]$beta2>beta2b,1,0)
  results[j,]$beta3trim<-ifelse(results[j,]$beta3<beta3a | results[j,]$beta3>beta3b,1,0)
  results[j,]$beta4trim<-ifelse(results[j,]$beta4<beta4a | results[j,]$beta4>beta4b,1,0)
  results[j,]$gamma1trim<-ifelse(results[j,]$gamma1<gamma1a | results[j,]$gamma1>gamma1b,1,0)
  results[j,]$gamma2trim<-ifelse(results[j,]$gamma2<gamma2a | results[j,]$gamma2>gamma2b,1,0)
  results[j,]$gamma3trim<-ifelse(results[j,]$gamma3<gamma3a | results[j,]$gamma3>gamma3b,1,0)
  results[j,]$gamma4trim<-ifelse(results[j,]$gamma4<gamma4a | results[j,]$gamma4>gamma4b,1,0)
  results[j,]$gamma5trim<-ifelse(results[j,]$gamma5<gamma5a | results[j,]$gamma5>gamma5b,1,0)
  results[j,]$gamma6trim<-ifelse(results[j,]$gamma6<gamma6a | results[j,]$gamma6>gamma6b,1,0)
  results[j,]$gamma7trim<-ifelse(results[j,]$gamma7<gamma7a | results[j,]$gamma7>gamma7b,1,0)
  results[j,]$gamma8trim<-ifelse(results[j,]$gamma8<gamma8a | results[j,]$gamma8>gamma8b,1,0)
  results[j,]$Psi_cctrim<-ifelse(results[j,]$Psi_cc<Psi_cca | results[j,]$Psi_cc>Psi_ccb,1,0)
  results[j,]$Psi_dd1trim<-ifelse(results[j,]$Psi_dd1<Psi_dd1a | results[j,]$Psi_dd1>Psi_dd1b,1,0)
  results[j,]$Psi_dd2trim<-ifelse(results[j,]$Psi_dd2<Psi_dd2a | results[j,]$Psi_dd2>Psi_dd2b,1,0)
  results[j,]$Psi_dd3trim<-ifelse(results[j,]$Psi_dd3<Psi_dd3a | results[j,]$Psi_dd3>Psi_dd3b,1,0)
  results[j,]$Psi_dd4trim<-ifelse(results[j,]$Psi_dd4<Psi_dd4a | results[j,]$Psi_dd4>Psi_dd4b,1,0)
  results[j,]$Var.V1trim<-ifelse(results[j,]$Var.V1<Var.V1a | results[j,]$Var.V1>Var.V1b,1,0)
  results[j,]$Var.V2trim<-ifelse(results[j,]$Var.V2<Var.V2a | results[j,]$Var.V2>Var.V2b,1,0)
  results[j,]$corrtrim<-ifelse(results[j,]$corr<corra | results[j,]$corr>corrb,1,0)
  results[j,]$Psi_dc1trim<-ifelse(results[j,]$Psi_dc1<Psi_dc1a | results[j,]$Psi_dc1>Psi_dc1b,1,0)
  results[j,]$Psi_dc2trim<-ifelse(results[j,]$Psi_dc2<Psi_dc2a | results[j,]$Psi_dc2>Psi_dc2b,1,0)
  results[j,]$mtrim<-ifelse(results[j,]$m<ma | results[j,]$m>mb,1,0)
  results[j,]$mu_1_1trim<-ifelse(results[j,]$mu_1_1<mu_1_1a | results[j,]$mu_1_1>mu_1_1b,1,0)
  results[j,]$mu_1_2trim<-ifelse(results[j,]$mu_1_2<mu_1_2a | results[j,]$mu_1_2>mu_1_2b,1,0)
}
results$trim<-ifelse(results$beta1trim==1 & results$beta2trim==1 & results$beta3trim==1 & results$beta4trim==1,1,0)
resultstrim<-results[results$trim==0,] #this data set deletes the simulations where the betas are outliers (all of them)

#summary
#Predictions - Estimates convergence
conv.summary<-aggregate(resultstrim$criterion.step,
              by = list(resultstrim$N,resultstrim$ni), length)
ps.summary<-aggregate(resultstrim$ps,
              by = list(resultstrim$N,resultstrim$ni), mean)
#Next line creates time for estimation of the coef in min
resultstrim$time.coef<-(resultstrim$end.coef-
                        resultstrim$inicio.coef)/60
#Next line creates time for estimation of the RE in min
resultstrim$time.RE<-(resultstrim$endp-resultstrim$iniciop)/60
time.summary<-aggregate(resultstrim$time.coef, 
              by = list(resultstrim$N,resultstrim$ni), mean)
timep.summary<-aggregate(resultstrim$time.RE, 
              by = list(resultstrim$N,resultstrim$ni), mean)
step.summary<-aggregate(resultstrim$step, 
              by = list(resultstrim$N,resultstrim$ni), mean)

mean_result<-aggregate(resultstrim, by = list(resultstrim$N,
             resultstrim$ni), function(x) mean(x, na.rm=TRUE))
write.csv(mean_result, "mean_result.csv")
sd_result<-aggregate(resultstrim, by = list(resultstrim$N,
             resultstrim$ni), function(x) sd(x, na.rm=TRUE))
write.csv(sd_result, "sd_result.csv")

#Subset other datasets based on resultstrim
resultstrim$reps<-paste(resultstrim$N,"_",resultstrim$ni,"_",
                        resultstrim$rep,sep = "")
V1$reps<-paste(V1$N,"_",V1$ni,"_",V1$rep,sep = "")
V1.trim<-subset(V1, reps %in% resultstrim$reps)
V1_p$reps<-paste(V1_p$N,"_",V1_p$ni,"_",V1_p$rep,sep = "")
V1_p.trim<-subset(V1_p, reps %in% resultstrim$reps)
V1_t$reps<-paste(V1_t$N,"_",V1_t$ni,"_",V1_t$rep,sep = "")
V1_t.trim<-subset(V1_t, reps %in% resultstrim$reps)

V2$reps<-paste(V2$N,"_",V2$ni,"_",V2$rep,sep = "")
V2.trim<-subset(V2, reps %in% resultstrim$reps)
V2_p$reps<-paste(V2_p$N,"_",V2_p$ni,"_",V2_p$rep,sep = "")
V2_p.trim<-subset(V2_p, reps %in% resultstrim$reps)
V2_t$reps<-paste(V2_t$N,"_",V2_t$ni,"_",V2_t$rep,sep = "")
V2_t.trim<-subset(V2_t, reps %in% resultstrim$reps)

t.V1<-aggregate(V1.trim, by = list(V1.trim$N,V1.trim$ni), 
                function(x) mean(x, na.rm=TRUE))
t.V1_p<-aggregate(V1_p.trim, by = list(V1_p.trim$N,V1_p.trim$ni), 
                  function(x) mean(x, na.rm=TRUE))
t.V1_t<-aggregate(V1_t.trim, by = list(V1_t.trim$N,V1_t.trim$ni), 
                  function(x) mean(x, na.rm=TRUE))
t.V2<-aggregate(V2.trim, by = list(V2.trim$N,V2.trim$ni), 
                function(x) mean(x, na.rm=TRUE))
t.V2_p<-aggregate(V2_p.trim, by = list(V2_p.trim$N,V2_p.trim$ni), 
                  function(x) mean(x, na.rm=TRUE))
t.V2_t<-aggregate(V2_t.trim, by = list(V2_t.trim$N,V2_t.trim$ni), 
                  function(x) mean(x, na.rm=TRUE))

s.V1<-aggregate(V1.trim, by = list(V1.trim$N,V1.trim$ni), 
                function(x) sd(x, na.rm=TRUE))
s.V1_p<-aggregate(V1_p.trim, by = list(V1_p.trim$N,V1_p.trim$ni), 
                  function(x) sd(x, na.rm=TRUE))
s.V1_t<-aggregate(V1_t.trim, by = list(V1_t.trim$N,V1_t.trim$ni), 
                  function(x) sd(x, na.rm=TRUE))
s.V2<-aggregate(V2.trim, by = list(V2.trim$N,V2.trim$ni), 
                function(x) sd(x, na.rm=TRUE))
s.V2_p<-aggregate(V2_p.trim, by = list(V2_p.trim$N,V2_p.trim$ni), 
                  function(x) sd(x, na.rm=TRUE))
s.V2_t<-aggregate(V2_t.trim, by = list(V2_t.trim$N,V2_t.trim$ni), 
                  function(x) sd(x, na.rm=TRUE))

#Subset other datasets based on resultstrim
naive$reps<-paste(naive$N,"_",naive$ni,"_",naive$rep,sep = "")
naive.trim<-subset(naive, reps %in% resultstrim$reps)
t.naive<-aggregate(naive.trim, by = list(naive.trim$N,naive.trim$ni),
                   function(x) mean(x, na.rm=TRUE))
s.naive<-aggregate(naive.trim, by = list(naive.trim$N,naive.trim$ni),
                   function(x) sd(x, na.rm=TRUE))

### Plots ############################# 
library(xlsx)
library(ggplot2)
# File Simulations 2020 05 10.xlsx has the organized results 
# of the simulations
pred<-read.xlsx("Simulations 2020 05 10.xlsx", sheetIndex=2, header=TRUE)
coef<-read.xlsx("Simulations 2020 05 10.xlsx", sheetIndex=6, header=TRUE)
bias<-read.xlsx("Simulations 2020 05 10.xlsx", sheetIndex=7, header=TRUE)
se_sim<-read.xlsx("Simulations 2020 05 10.xlsx", sheetIndex=8, header=TRUE)
sd_sim<-read.xlsx("Simulations 2020 05 10.xlsx", sheetIndex=9, header=TRUE)
cov_se<-read.xlsx("Simulations 2020 05 10.xlsx", sheetIndex=10, header=TRUE)

pred<-pred[complete.cases(pred),]   
coef<-coef[complete.cases(coef),]   
bias<-bias[complete.cases(bias),]   
cov_se<-cov_se[complete.cases(cov_se),]   
sd_sim<-sd_sim[complete.cases(sd_sim),]
se_sim<-se_sim[complete.cases(se_sim),]

pred<-pred[order(pred$Patches,pred$ni),]   
coef<-coef[order(coef$Patches,coef$ni),]   
bias<-bias[order(bias$Patches,bias$ni),]   
cov_se<-cov_se[order(cov_se$Patches,cov_se$ni),]   
sd_sim<-sd_sim[order(sd_sim$Patches,sd_sim$ni),]   
se_sim<-se_sim[order(se_sim$Patches,se_sim$ni),]   

#Convrgence
ggplot(data=pred,aes(x = ni,y = Estimates..Convergence,
                     col=factor(Patches))) +
  geom_jitter(width = 0.005,aes(shape=factor(Patches)),size = 2) + 
  scale_shape_manual(name="Number of\nPatches",values=c(16, 3, 17)) +
  scale_color_manual(name="Number of\nPatches",
                     values=c("#926b8d","#7aadb1","#E8CA47")) +
  scale_x_continuous(name="Time points", breaks=c(3,4,5)) + 
  scale_y_continuous(name="Convergence (%)", 
                     breaks=seq(98,100,0.5),limits = c(98,100)) + 

ggplot(data=pred,aes(x = ni,y = Average..time.RE..m.,
                     col=factor(Patches))) +
  geom_jitter(width = 0.005,aes(shape=factor(Patches)),size = 2) + 
  scale_shape_manual(name="Number of\nPatches",values=c(16, 3, 17)) +
  scale_color_manual(name="Number of\nPatches",
                     values=c("#926b8d","#7aadb1","#E8CA47")) +
  scale_x_continuous(name="Time points", breaks=c(3,4,5)) + 
  scale_y_continuous(name="Average Time RE (min)") + 
  
ggplot(data=pred,aes(x = ni,y = Average..time.coef..m.,
                     col=factor(Patches))) +
  geom_jitter(width = 0.005,aes(shape=factor(Patches)),size = 2) + 
  scale_shape_manual(name="Number of\nPatches",values=c(16, 3, 17)) +
  scale_color_manual(name="Number of\nPatches",
                     values=c("#926b8d","#7aadb1","#E8CA47")) +
  scale_x_continuous(name="Time points", breaks=c(3,4,5)) + 
  scale_y_continuous(name="Average Time Coefficients (min)") + 
  
ggplot(data=pred,aes(x = ni,y = Average..steps,col=factor(Patches)))+
  geom_jitter(width = 0.005,aes(shape=factor(Patches)),size = 2) + 
  scale_shape_manual(name="Number of\nPatches",values=c(16, 3, 17)) +
  scale_color_manual(name="Number of\nPatches",
                     values=c("#926b8d","#7aadb1","#E8CA47")) +
  scale_x_continuous(name="Time points", breaks=c(3,4,5)) + 
  scale_y_continuous(name="Average iterations before convergence",
                     limits = c(10,15)) + 
  
### Coefficients  ####################################
coeff<-cbind(coef,sd_sim)
coeff2<-cbind(coef,se_sim)
colnames(coeff)<-c("Patches","ni","Scenario","Beta1","Beta2","Beta3",
                   "Beta4","Gamma1","Gamma2","Gamma3","Gamma4",
                   "Gamma5","Gamma6","Gamma7","Gamma8","Patches.s",
                   "ni.s","Scenario.s","Beta1.s","Beta2.s","Beta3.s",   
                   "Beta4.s","Gamma1.s","Gamma2.s","Gamma3.s",
                   "Gamma4.s","Gamma5.s","Gamma6.s","Gamma7.s",
                   "Gamma8.s")
colnames(coeff2)<-c("Patches","ni","Scenario","Beta1","Beta2","Beta3",
                    "Beta4","Gamma1","Gamma2","Gamma3","Gamma4",
                    "Gamma5","Gamma6","Gamma7","Gamma8","Patches.s",
                    "ni.s","Scenario.s","Beta1.s","Beta2.s","Beta3.s",   
                    "Beta4.s","Gamma1.s","Gamma2.s","Gamma3.s",
                    "Gamma4.s","Gamma5.s","Gamma6.s","Gamma7.s",
                    "Gamma8.s")
#plot with no bars
ggplot(data=coeff,aes(x = ni,y = Beta1,col=factor(Patches))) +
  geom_jitter(width = 0.1,aes(shape=factor(Patches)),size = 2) + 
  geom_hline(yintercept=-10,color="#B7B6BA",linetype="dashed")+
  scale_shape_manual(name="Number of\nPatches",values=c(16, 3, 17))+
  scale_color_manual(name="Number of\nPatches",
                     values=c("#926b8d","#7aadb1","#E8CA47")) +
  scale_x_continuous(name="Time points", breaks=c(3,4,5)) + 
  scale_y_continuous(name="Intercept") + 

#plot with sd of the estimates
ggplot(coeff, aes(x=ni, y=Beta1,color=factor(Patches))) + 
  geom_point(position=position_dodge(0.5),aes(shape=factor(Patches)),
             size = 2)+
  geom_errorbar(aes(ymin=Beta1-(2*Beta1.s), ymax=Beta1+(2*Beta1.s)),
                width=.2,position=position_dodge(0.5))+
  geom_hline(yintercept=-10,color="#B7B6BA",linetype="dashed")+
  scale_shape_manual(name="Number of\nPatches",values=c(16, 3, 17)) +
  scale_color_manual(name="Number of\nPatches",
                     values=c("#926b8d","#7aadb1","#E8CA47")) +
  scale_x_continuous(name="Time points", breaks=c(3,4,5)) + 
  scale_y_continuous(name="Intercept") + 
 
#plot with average se
ggplot(coeff2, aes(x=ni, y=Beta1,color=factor(Patches))) + 
  geom_point(position=position_dodge(0.5),aes(shape=factor(Patches)),
             size = 2)+
  geom_errorbar(aes(ymin=Beta1-(1.96*Beta1.s), ymax=Beta1+
                      (1.96*Beta1.s)), width=.2,
                position=position_dodge(0.5))+
  geom_hline(yintercept=-10,color="#B7B6BA",linetype="dashed")+
  scale_shape_manual(name="Number of\nPatches",values=c(16, 3, 17)) +
  scale_color_manual(name="Number of\nPatches",
                     values=c("#926b8d","#7aadb1","#E8CA47")) +
  scale_x_continuous(name="Time points", breaks=c(3,4,5)) + 
  scale_y_continuous(name="Intercept") + 
 
# Same plots for Beta2, Beta3, Beta4, Gamma1, Gamma2, Gamma3, Gamma4,
# Gamma5, Gamma6, Gamma7, Gamma8  
  
### Bias        #############################################
ggplot(data=bias,aes(x = ni,y = Beta1,col=factor(Patches))) +
  geom_jitter(width = 0.1,aes(shape=factor(Patches)),size = 2) + 
  scale_shape_manual(name="Number of\nPatches",values=c(16, 3, 17)) +
  scale_color_manual(name="Number of\nPatches",values=c("#926b8d","#7aadb1","#E8CA47")) +
  scale_x_continuous(name="Time points", breaks=c(3,4,5)) + 
  scale_y_continuous(name="Bias in Intercept") + 

# Same plots for Beta2, Beta3, Beta4, Gamma1, Gamma2, Gamma3, Gamma4,
# Gamma5, Gamma6, Gamma7, Gamma8  
  
### Plots predictions  #############################################

t.V1<- merge(t.V1,s.V1,by=c("Group.1","Group.2"),all=T)
t.V1_p<-merge(t.V1_p,s.V1_p,by=c("Group.1","Group.2"),all=T)
t.V1_t<-merge(t.V1_t,s.V1_t,by=c("Group.1","Group.2"),all=T)
t.V2<- merge(t.V2,s.V2,by=c("Group.1","Group.2"),all=T)
t.V2_p<-merge(t.V2_p,s.V2_p,by=c("Group.1","Group.2"),all=T)
t.V2_t<-merge(t.V2_t,s.V2_t,by=c("Group.1","Group.2"),all=T)

colnames(t.V1)<- c("Group.1","Group.2","N","ni","rep","file.x",
                   "file2.x","m_Vu","TSS","sse.V","ssr.V","tss",
                   "corr2","R2","mse","mae","mase","rmse","tvs",
                   "uvs","evs","reps","N.s","ni.s","rep.s","file.s",
                   "file2.s","m_Vu.s","TSS.s","sse.V.s","ssr.V.s",
                   "tss.s","corr2.s","R2.s","mse.s","mae.s","mase.s",
                   "rmse.s","tvs.s","uvs.s","evs.s","reps.y")
colnames(t.V1_p)<-colnames(t.V1)
colnames(t.V1_t)<-colnames(t.V1)
colnames(t.V2)<-colnames(t.V1)
colnames(t.V2_p)<-colnames(t.V1)
colnames(t.V2_t)<-colnames(t.V1)

t.V1<-t.V1[order(t.V1$N,t.V1$ni),]   
t.V1_p<-t.V1_p[order(t.V1_p$N,t.V1_p$ni),] 
t.V1_t<-t.V1_t[order(t.V1_t$N,t.V1_t$ni),] 
t.V2<-t.V2[order(t.V2$N,t.V2$ni),]   
t.V2_p<-t.V2_p[order(t.V2_p$N,t.V2_p$ni),] 
t.V2_t<-t.V2_t[order(t.V2_t$N,t.V2_t$ni),] 

ggplot(t.V1_p, aes(x=ni, y=corr2,color=factor(N))) + 
  geom_point(position=position_dodge(0.5),aes(shape=factor(N)),
             size = 2)+
  geom_errorbar(aes(ymin=corr2-(1.96*corr2.s), 
                    ymax=min(1,corr2+(1.96*corr2.s))), width=.2,
                position=position_dodge(0.5))+
  scale_shape_manual(name="Number of\nPatches",values=c(16, 3, 17)) +
  scale_color_manual(name="Number of\nPatches",
                     values=c("#926b8d","#7aadb1","#E8CA47")) +
  scale_x_continuous(name="Time points", breaks=c(3,4,5)) + 
  scale_y_continuous(name=expression(paste("Corr"^2,
                     "observed vs. predicted")),limits = c(-1,1)) + 

ggplot(t.V1_p, aes(x=ni, y=mase,color=factor(N))) + 
  geom_point(position=position_dodge(0.5),aes(shape=factor(N)),
             size = 2)+
  geom_errorbar(aes(ymin=mase-(1.96*mase.s), ymax=mase+1.96*mase.s),
                width=.2,position=position_dodge(0.5))+
  geom_hline(yintercept=1,color="#B7B6BA",linetype="dashed")+
  scale_shape_manual(name="Number of\nPatches",values=c(16, 3, 17)) +
  scale_color_manual(name="Number of\nPatches",
                     values=c("#926b8d","#7aadb1","#E8CA47")) +
  scale_x_continuous(name="Time points", breaks=c(3,4,5)) + 
  scale_y_continuous(name="MASE observed vs. predicted",
                     limits = c(0,1.6)) + 

# Same plots for t.V2_p

### Plots predictions with naive ###################################
t.naive<- merge(t.naive,s.naive,by=c("Group.1","Group.2"),all=T)
t.naive<-t.naive[order(t.naive$N.x,t.naive$ni.x),]   

ggplot(t.naive, aes(x=ni.x, y=naive.corr2.V1.x,color=factor(N.x))) + 
  geom_point(position=position_dodge(0.5),aes(shape=factor(N.x)),
             size = 2)+
  geom_errorbar(aes(ymin=naive.corr2.V1.x-(1.96*naive.corr2.V1.y),
                    ymax=min(1,naive.corr2.V1.x+
                               (1.96*naive.corr2.V1.y))), width=.2,
                position=position_dodge(0.5))+
  scale_shape_manual(name="Number of\nPatches",values=c(16, 3, 17)) +
  scale_color_manual(name="Number of\nPatches",
                     values=c("#926b8d","#7aadb1","#E8CA47")) +
  scale_x_continuous(name="Time points", breaks=c(3,4,5)) + 
  scale_y_continuous(name=expression(paste("Corr"^2,
                    "observed vs. predicted")),limits = c(-1,1)) + 

ggplot(t.naive, aes(x=ni.x, y=naive.mase.V1.x,color=factor(N.x))) + 
  geom_point(position=position_dodge(0.5),aes(shape=factor(N.x)),
             size = 2)+
  geom_errorbar(aes(ymin=naive.mase.V1.x-(1.96*naive.mase.V1.y),
                    ymax=naive.mase.V1.x+(1.96*naive.mase.V1.y)),
                width=.2,position=position_dodge(0.5))+
  geom_hline(yintercept=1,color="#B7B6BA",linetype="dashed")+
  scale_shape_manual(name="Number of\nPatches",values=c(16, 3, 17)) +
  scale_color_manual(name="Number of\nPatches",
                     values=c("#926b8d","#7aadb1","#E8CA47")) +
  scale_x_continuous(name="Time points", breaks=c(3,4,5)) + 
  scale_y_continuous(name="MASE observed vs. predicted",
                     limits = c(0,1.6)) + 

# Same plots for naive.corr2.V2, naive.mase.V1.x
