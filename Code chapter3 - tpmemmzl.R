#Functions --------------------------
p.inf.l<-function(x){
  #Some propi<-exp(eta)/(1+exp(eta)) can be inf, too small or large, then we are truncating the values
  x<-ifelse(is.finite(x),x,0.999999)
  x<-ifelse(x<0.000001,0.000001,x)
  x<-ifelse(x>0.999999,0.999999,x)
}

p.inf<-function(x){
  #Some propi<-exp(eta)/(1+exp(eta)) can be inf, too small or large, then we are truncating the values
  x<-ifelse(is.finite(x),x,0.99)
  x<-ifelse(x<1e-6,1e-6,x)
  x<-ifelse(x>0.98,0.98,x)
}

laplace2<-function (logpost, mode, ...) 
{
  options(warn = -1)
  fit = optim(mode, logpost, gr = NULL, ..., hessian = FALSE, 
              control = list(fnscale = -1))
  options(warn = 0)
  mode = fit$par
  stuff = list(mode = mode, converge = fit$convergence == 0)
  return(stuff)
}

lap.for.c<-function(ram,other){ #find mode for random effects
  ram<-as.matrix(ram)
  U<-other$U
  X<-as.matrix(other$X)
  Z<-as.matrix(other$Z)
  Beta<-as.matrix(other$Beta)
  Psi_cc.i<-as.matrix(other$Psi_cc.i)
  D<-as.matrix(other$D)
  g.To<-as.matrix(other$g.To)
  g.o<-t(g.To)
  eta<- (X %*% Beta) + (Z %*% ram)
  propi<-exp(eta)/(1+exp(eta))
  propi<-apply(propi, 1, p.inf.l) 
  eta<- log(propi/(1-propi))
  lu<-sum((U * eta) + log(1-propi))
  mat<-Psi_cc.i+D
  mat.i<-solve(Psi_cc.i+D)
  qua<-t(ram-(mat.i %*% g.o)) %*% mat %*% (ram-(mat.i %*% g.o))
  fun<-lu-((1/2)*qua)
  return(fun)
}

lklhd.em<-function(id,datafile,initial,X.var=X.var,Z.var=Z.var,X_a.var=X_a.var,Z_a.var=Z_a.var,U.var=U.var,V1=V1,V2=V2){    
  #Estimation of mixed model for random effects (m) 
  Beta=initial$Beta
  Gamma=initial$Gamma
  Gamma<-matrix(Gamma,ncol=1)
  Psi_cc=initial$Psi_cc
  Psi_dd=initial$Psi_dd
  corr<-initial$corr
  Var.V1<-initial$Var.V1
  Var.V2<-initial$Var.V2
  Psi_dc=initial$Psi_dc
  m=initial$m
  mu_1=initial$mu_1
  
  p<-length(Beta)
  q<-ifelse(length(Psi_cc)>1,dim(Psi_cc)[1],length(Psi_cc))
  p_a<-length(Gamma)
  q_a<-ifelse(length(Psi_dd)>1,dim(Psi_dd)[1],length(Psi_dd))
  Sigma<-matrix(c(Var.V1,rep(corr*sqrt(Var.V1)*sqrt(Var.V2),2),Var.V2),nrow = 2,ncol = 2)
  
  #Data
  base<-datafile[datafile$patch==id,]
  base$int <- rep(1,nrow(base))
  U<-base[,U.var]
  X.var<-c("int",X.var)
  X<-base[,X.var]
  X<-as.matrix(X)
  Z.var<-c("int",Z.var)
  Z<-base[,Z.var]
  Z<-as.matrix(Z)
  ni<-nrow(Z)
  
  X_a.var<-c("int",X_a.var)
  X_a<-base[,X_a.var]
  X_a<-as.matrix(X_a)
  X_a<-as.matrix(U * X_a)
  
  Z_a.var<-c("int",Z_a.var)
  Z_a<-base[,Z_a.var]
  Z_a<-as.matrix(Z_a)
  Z_a<-as.matrix(U * Z_a)
  n_a<- sum(Z_a[,1]) 
  V_a<-c(base[,V1],base[,V2])  
  
  
  #Calculations
  c<-as.matrix(rep(0,length(Z.var)),ncol=1)
  mu_2<- -(m/(1-m))*mu_1
  Psi_cc.i<-solve(Psi_cc)
  Psi_cd<-t(Psi_dc)
  H<-Psi_dd-(Psi_dc %*% Psi_cc.i %*% Psi_cd )
  H.i<-ginv(H)
  X_at<-diag(2) %x% as.matrix(X_a)
  Z_at<-diag(2) %x% as.matrix(Z_a)
  Sigma.i<-ginv(Sigma)
  Sigma.t<-Sigma %x% diag(nrow(X_a))
  Sigma.it<-Sigma.i %x% diag(nrow(X_a))
  Sigma.vech<-vech(Sigma)
  
  Psi_cc.vech<-vech(Psi_cc)
  Psi_ddcci<-Psi_dc %*% solve(Psi_cc)  #Check this in the future (when more random effects, if still a vector)
  H.vech<-vech(H)
  
  A_a<-matrix(V_a-X_at %*% Gamma)
  
  tempcheck<-ifelse(sum(V_a==0)==length(V_a),0,1) #if sum is zero, the all observations are zero
  if(tempcheck==1){
    B<-(t(Z_at) %*% Sigma.it %*%  Z_at) + H.i
    B.i<-ginv(B)
  } else {
    B<- H.i
    B.i<-H
  }
  
  E.T<-t(A_a) %*% Sigma.it %*% Z_at %*% B.i %*% H.i %*% Psi_dc %*% Psi_cc.i 
  M_1.T<-m*t(mu_1)
  M_2.T<-(1-m)*t(mu_2)
  M.T<-M_1.T+M_2.T  
  D<- t(Psi_dc %*% Psi_cc.i) %*% H.i %*% B.i %*% t(Z_at) %*%  Sigma.it %*% Z_at %*% Psi_dc %*% Psi_cc.i
  g_1.T<-E.T+(M_1.T %*% H.i %*% B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% Psi_dc %*% Psi_cc.i) 
  g_2.T<-E.T+(M_2.T %*% H.i %*% B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% Psi_dc %*% Psi_cc.i) 
  g.T <- E.T+(M.T   %*% H.i %*% B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% Psi_dc %*% Psi_cc.i) 
  
  # Max part for m
  other<-list(U=U, X=X, Z=Z, Beta=Beta, Psi_cc.i = Psi_cc.i, D=D, g.To=g_1.T)
  c_1t<- laplace2(lap.for.c,c,other)$mode 
  other<-list(U=U, X=X, Z=Z, Beta=Beta, Psi_cc.i = Psi_cc.i, D=D, g.To=g_2.T)
  c_2t<- laplace2(lap.for.c,c,other)$mode 
  
  eta.1<- (X %*% Beta) + (Z %*% c_1t)
  propi.1<-exp(eta.1)/(1+exp(eta.1))
  propi.1<-apply(propi.1, 1, p.inf)
  eta.1<- log(propi.1/(1-propi.1))
  W_t.1<- diag(c(propi.1),ncol=length(propi.1)) %*% diag(c(1-propi.1),ncol=length(propi.1))
  G.1<-(t(Z) %*% W_t.1 %*% Z ) + (Psi_cc.i+D)
  G.1.i<-solve(G.1)
  
  eta.2<- (X %*% Beta) + (Z %*% c_2t)
  propi.2<-exp(eta.2)/(1+exp(eta.2))
  propi.2<-apply(propi.2, 1, p.inf)
  eta.2<- log(propi.2/(1-propi.2))
  W_t.2<- diag(c(propi.2),ncol=length(propi.2)) %*% diag(c(1-propi.2),ncol=length(propi.2))
  G.2<-(t(Z) %*% W_t.2 %*% Z ) + (Psi_cc.i+D)
  G.2.i<-solve(G.2)
  
  lu<-sum((U * eta.1) + log(1-propi.1))
  qua<- t(c_1t-(solve(Psi_cc.i+D)%*%t(g_1.T))) %*% (Psi_cc.i+D) %*% (c_1t-(solve(Psi_cc.i+D)%*%t(g_1.T)))
  f_1<-lu-((1/2)*qua)
  lu<-sum((U * eta.2) + log(1-propi.2))
  qua<-t(c_2t-(solve(Psi_cc.i+D)%*%t(g_2.T))) %*% (Psi_cc.i+D) %*% (c_2t-(solve(Psi_cc.i+D)%*%t(g_2.T)))
  f_2<-lu-((1/2)*qua)
  
  m3_t.1<- diag(W_t.1) * (1-2*propi.1)
  m4_t.1<- diag(W_t.1 %*% (diag(dim(W_t.1)[1])-6*W_t.1))
  m6_t.1<- m4_t.1*diag(diag(dim(W_t.1)[1])-12*W_t.1)-(12*(m3_t.1^2))      
  P1.1<- (1/8)* m4_t.1 %*% diag(Z %*% G.1.i %*% t(Z))^2
  P2.1<- (1/48)*t(m6_t.1) %*% diag(Z %*% G.1.i %*% t(Z))^3
  P3.1<- (15/72)* t(colSums(Z *  m3_t.1 * diag(Z %*% G.1.i %*% t(Z)))) %*% G.1.i %*% colSums(Z *  m3_t.1 * diag(Z %*% G.1.i %*% t(Z)))  #Check that this runs well with more random efects
  P.1<- 1 - P1.1 -P2.1 + P3.1
  P.1<- ifelse(P.1>=1e-25,P.1,0.5)
  lf1<- -(n_a+(q/2))*log(2*pi) - (n_a/2)*log(det(Sigma)) - (1/2)*log(det(H)) + (1/2)*log(det(B.i)) - (1/2)*log(det(Psi_cc)) -
    ((1/2)* t(A_a) %*% (Sigma.it - Sigma.it %*% Z_at %*% B.i %*% t(Z_at) %*% Sigma.it) %*% A_a) + 
    (t(A_a) %*% Sigma.it %*% Z_at %*% B.i %*% H.i %*% t(M_1.T)) -
    ((1/2)*m*t(mu_1) %*% H.i %*% B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% mu_1) + 
    ((1/2)* g_1.T %*% solve(Psi_cc.i+D) %*% t(g_1.T)) + f_1 + (1/2)*log(det(G.1))+
    log(P.1)
  f1<-exp(lf1)
  m3_t.2<- diag(W_t.2) * (1-2*propi.2)
  m4_t.2<- diag(W_t.2 %*% (diag(dim(W_t.2)[1])-6*W_t.2))
  m6_t.2<- m4_t.2*diag(diag(dim(W_t.2)[1])-12*W_t.2)-(12*m3_t.2^2)             
  P1.2<- (1/8)* m4_t.2 %*% diag(Z %*% G.2.i %*% t(Z))^2
  P2.2<- (1/48)*t(m6_t.2) %*% diag(Z %*% G.2.i %*% t(Z))^3
  P3.2<- (15/72)* t(colSums(Z *  m3_t.2 * diag(Z %*% G.2.i %*% t(Z)))) %*% G.2.i %*% colSums(Z *  m3_t.2 * diag(Z %*% G.2.i %*% t(Z)))  #Check that this runs well with more random efects
  P.2<- 1 - P1.2 -P2.2 + P3.2
  P.2<- ifelse(P.2>=1e-25,P.1,0.5)
  lf2<- -(n_a+(q/2))*log(2*pi) - (n_a/2)*log(det(Sigma)) - (1/2)*log(det(H)) + (1/2)*log(det(B.i)) - (1/2)*log(det(Psi_cc)) -
    ((1/2)* t(A_a) %*% (Sigma.it - Sigma.it %*% Z_at %*% B.i %*% t(Z_at) %*% Sigma.it) %*% A_a) + 
    (t(A_a) %*% Sigma.it %*% Z_at %*% B.i %*% H.i %*% t(M_2.T)) -
    ((1/2)*(1-m)*t(mu_2) %*% H.i %*% B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% mu_2) + 
    ((1/2)* g_2.T %*% solve(Psi_cc.i+D) %*% t(g_2.T)) + f_2 + (1/2)*log(det(G.2))+
    log(P.2)  
  f2<-exp(lf2)
  
  m.new<-f1/(f1+f2)
  m<-ifelse(is.na(m.new),0,m.new)
  m<-ifelse(m<0.01,0.01,m)
  m<-ifelse(m>0.99,0.99,m)
  
  #max part for other par
  other<-list(U=U, X=X, Z=Z, Beta=Beta, Psi_cc.i = Psi_cc.i, D=D, g.To=E.T)
  c_t<- laplace2(lap.for.c,c,other)$mode 
  
  eta<- (X %*% Beta) + (Z %*% c_t)
  propi<-exp(eta)/(1+exp(eta))
  propi<-apply(propi, 1, p.inf)
  eta<- log(propi/(1-propi))
  lu<-sum((U * eta) + log(1-propi))
  qua<-t(c_t-(solve(Psi_cc.i+D)%*%t(g.T))) %*% (Psi_cc.i+D) %*% (c_t-(solve(Psi_cc.i+D)%*%t(g.T)))
  f_v2<-lu-((1/2)*qua)
  W_t<- diag(c(propi),ncol=length(propi)) %*% diag(c(1-propi),ncol=length(propi))
  U_a<- (diag(1/diag(W_t)) %*% (U - propi)) + eta
  G<-(t(Z) %*% W_t %*% Z ) + (Psi_cc.i+D)
  G.i<-solve(G)
  
  m3_t<-diag(W_t) * (1-2*propi)
  m4_t<-diag(W_t %*% (diag(dim(W_t)[1])-6*W_t))
  m5_t<-m3_t*diag(diag(dim(W_t)[1])-12*W_t)
  m6_t<-m4_t*diag(diag(dim(W_t)[1])-12*W_t)-(12*m3_t^2)             
  m7_t<-m5_t*diag(diag(dim(W_t)[1])-12*W_t)-(36*m3_t*m4_t) 
  
  P1<- (1/8)* m4_t %*% diag(Z %*% G.i %*% t(Z))^2
  P2<- (1/48)*t(m6_t) %*% diag(Z %*% G.i %*% t(Z))^3
  P3<- (15/72) * t(colSums(Z *  m3_t * diag(Z %*% G.i %*% t(Z)))) %*% G.i %*% colSums(Z *  m3_t * diag(Z %*% G.i %*% t(Z)))  #Check that this runs well with more random efects
  P<- (1 - P1 - P2 + P3)
  P<- ifelse(P>1e-25,P,0.5)
  
  #Scoring calculations
  Ki<-matrix(0,nrow=ncol(Z) ,ncol=ncol(Z) )
  for (k in 1:nrow(Z)) {
    Ki.temp<-t(Z[k,]) %*% m3_t[k] %*% Z[k,] %*% G.i %*%  t(Z[k,])
    Ki=Ki+Ki.temp }
  hi<-matrix(0,nrow=ncol(Z),ncol =ncol(Z) )
  for (k in 1:nrow(Z)) {
    h.temp<-m3_t[k] %*% G.i %*%t(Z[k,]) %*% (Z[k,] %*% G.i) %*% (Z[k,] %*% G.i %*% Ki)
    hi=hi+h.temp }
  fi<-matrix(0,nrow=ncol(Z),ncol =ncol(Z) )
  for (k in 1:nrow(Z)) {
    fi.temp<-m6_t[k] %*% G.i %*%t(Z[k,]) %*% ((Z[k,] %*% G.i %*% t(Z[k,]))^2) %*% (Z[k,] %*% G.i)
    fi=fi+fi.temp }
  Fi<-matrix(0,nrow=ncol(Z),ncol =ncol(Z) )
  for (k in 1:nrow(Z)) {
    Fi.temp<-m4_t[k] %*% G.i %*% t(Z[k,]) %*% (Z[k,] %*% G.i %*% t(Z[k,])) %*% (Z[k,] %*% G.i)
    Fi=Fi+Fi.temp }
  Rij<-matrix(0,nrow=nrow(Z),ncol=1)
  for (k in 1:nrow(Z)) {
    Rij[k,1] <-  -((1/8)*m5_t[k] %*% (Z[k,] %*% G.i %*% t(Z[k,]))^2) + ((1/4)* m3_t[k] %*% Z[k,] %*% Fi %*% t(Z[k,])) - 
      ((1/48)*m7_t[k] %*% (Z[k,] %*% G.i %*% t(Z[k,]))^3) +  ((1/16)* m3_t[k] %*% Z[k,] %*% fi %*% t(Z[k,])) - 
      ((15/72)*m3_t[k] %*% (Z[k,] %*% G.i %*% Ki)^2) + ((15/36)* m4_t[k] %*% (t(Ki) %*% G.i %*% t(Z[k,])) %*% (Z[k,] %*% G.i %*% t(Z[k,]))) -
      ((15/36)*m3_t[k] %*% Z[k,] %*% hi %*% t(Z[k,]))
  }
  Fij<-matrix(0,nrow=nrow(Z),ncol=1)
  for (k in 1: nrow(Z)) {
    Fij[k,1]<- ((1/4)*m4_t[k] %*% (Z[k,] %*% G.i %*% t(Z[k,]))) + ((1/16)*m6_t[k] %*% (Z[k,] %*% G.i %*% t(Z[k,]))^2) -
      ((15/36)*m3_t[k] %*% (Z[k,] %*% G.i %*% Ki))      }
  
  Yi<- t(Z) %*% W_t %*% (U_a - (eta-(Z %*% c_t))) + t(E.T)
  Hi<- t(Z_at) %*% Sigma.it %*% Z_at %*% B.i %*% H.i %*% Psi_dc %*% Psi_cc.i
  
  Eij.int<-matrix(0,nrow=nrow(Z),ncol=1)
  for (k in 1:nrow(Z)) {
    Eij.int[k,1]<- t(Z[k,]) %*% (((-1/2)*m3_t[k] %*% (Z[k,] %*% G.i %*% t(Z[k,]))) + (U[k] - propi[k]) + (Rij[k]/P) )}
  Eij<- t(E.T)-(Psi_cc.i+D)%*%c_t +colSums(Eij.int)
  
  Si<-((t(A_a) %*% Sigma.it %*% Z_at %*% B.i %*% H.i) %x% G.i)-(G.i %x% t(Yi) %*% G.i %*% t(Hi))-(t(Yi) %*% G.i %x% G.i %*% t(Hi))
  ti<- -G.i %*% t(Psi_dc %*% Psi_cc.i) %*% H.i %*% B.i %*% t(Z_at) %*% Sigma.it %*% X_at
  Di<- -(t(Yi) %*% G.i %x% G.i)
  
  Mi<- -(t(Yi) %*% G.i %*% t(Psi_dc %*% Psi_cc.i) %*% t(Z_at) %*% Sigma.it %*% Z_at %*% B.i %x% (G.i %*% t(Psi_dc %*% Psi_cc.i)))+
    (t(A_a) %*% Sigma.it %*% Z_at %*% B.i %x% (G.i %*% t(Psi_dc %*% Psi_cc.i)))-
    (t(A_a) %*% Sigma.it %*% Z_at %*% B.i %x% (G.i %*% t(Psi_dc %*% Psi_cc.i) %*% H.i %*% B.i))+
    (t(Yi) %*% G.i %*% t(Psi_dc %*% Psi_cc.i) %*% t(Z_at) %*% Sigma.it %*% Z_at %*% B.i %x% (G.i %*% t(Psi_dc %*% Psi_cc.i) %*% H.i %*% B.i)) 
  
  Pi<- (t(Yi) %*% G.i %*% t(Psi_dc %*% Psi_cc.i) %*% t(Z_at) %*% Sigma.it %*% Z_at %*% B.i %*% t(Z_at) %x% (G.i %*% t(Psi_dc %*% Psi_cc.i) %*% H.i %*% B.i %*% t(Z_at)))-
    (t(Yi) %*% G.i %*% t(Psi_dc %*% Psi_cc.i) %*% t(Z_at) %x% (G.i %*% t(Psi_dc %*% Psi_cc.i) %*% H.i %*% B.i %*% t(Z_at)))-
    (t(A_a) %*% Sigma.it %*% Z_at %*% B.i %*% t(Z_at) %x% (G.i %*% t(Psi_dc %*% Psi_cc.i) %*% H.i %*% B.i %*% t(Z_at)))+
    (t(A_a) %x% (G.i %*% t(Psi_dc %*% Psi_cc.i)%*% H.i %*% B.i %*% t(Z_at)))
  
  #Score vectors
  S_b.int1<-matrix(0,nrow=nrow(Z),ncol=ncol(X))
  for (k in 1:nrow(Z)) {
    S_b.int1[k,]<-c(((-1/2)*m3_t[k] %*% (Z[k,] %*% G.i %*% t(Z[k,])))) * (t(X[k,]) - t(t(X) %*% W_t %*% Z %*% G.i %*% t(Z[k,])))}
  S_b.int1<-t(S_b.int1)
  S_b.int1<-rowSums(S_b.int1)
  S_b.int2<-matrix(0,nrow=nrow(Z),ncol=ncol(X))
  for (k in 1:nrow(Z)) {
    S_b.int2[k,]<-Rij[k]*(t(X[k,]) - t(t(X) %*% W_t %*% Z %*% G.i %*% t(Z[k,])))}
  S_b.int2<-t(S_b.int2)
  S_b.int2<-rowSums(S_b.int2)
  S_b<- t(X) %*% W_t %*% (U_a - eta) + S_b.int1 + (c((1/P)) * S_b.int2)
  rule1<-sum(abs(S_b)<1e-1)
  if(rule1>0){
    S_b<-matrix(0,nrow=ncol(X),1)
  }
  
  S_g<-(t(X_at) %*% Sigma.it %*% V_a) - (t(X_at) %*% Sigma.it %*%  Z_at %*% B.i %*% t(Z_at) %*%  Sigma.it %*% V_a) -
    (t(X_at) %*% Sigma.it %*% X_at %*% Gamma) + (t(X_at) %*% Sigma.it %*%  Z_at %*% B.i %*% t(Z_at) %*%  Sigma.it %*% X_at %*% Gamma) -
    (t(X_at) %*% Sigma.it %*%  Z_at %*% B.i %*% H.i %*% Psi_dc %*% Psi_cc.i %*% c_t) + (t(ti) %*% Eij)
  
  mu_1<-matrix(mu_1,ncol=1)
  S_mu.int1<- (H.i %*% B.i %*% t(Z_at) %*% Sigma.it %*% Z_at)
  S_mu.int2<- (t(Z_at) %*% Sigma.it %*% Z_at %*% B.i %*% H.i)
  S_mu.int3<- (S_mu.int1 + S_mu.int2) %*% mu_1
  S_mu.int4<- -(1/2)*(m+((m^2)/(1-m)))
  S_mu<- c(S_mu.int4) * S_mu.int3
  
  S_pcc.int1<-ifelse(q==1,t((Psi_cc.i %x% Psi_cc.i)),t((Psi_cc.i %x% Psi_cc.i) %*% duplication.matrix(n=q)))
  S_pcc.int2<-matrix(0,nrow=ncol(Z)*ncol(Z),ncol=ni)
  for (k in 1:nrow(Z)) {
    S_pcc.int2[,k]<- vec(G.i %*% t(Z[k,]) %*% Z[k,] %*% G.i) *  Fij[k]}
  S_pcc.int2<-rowSums(S_pcc.int2)
  S_pcc<- -S_pcc.int1 %*% ( (1/2)*vec(Psi_cc - G.i - c_t %*% t(c_t)) + t(Di) %*% Eij + 
                              (1/P) %*% ( S_pcc.int2 - (15/72)*vec(G.i %*% Ki %*% t(Ki) %*% G.i ) ) ) 
  S_pcc<- as.matrix(apply(S_pcc, 1, function(x) {ifelse(abs(x)<1e-1 & x!=0, 0, x)}),ncol=1)
  
  S_pdc.int1<-matrix(0,nrow=ncol(Z_at)*ncol(Z),ncol=ni)
  for (k in 1:nrow(Z)) {
    S_pdc.int1[,k]<- vec(Hi %*% G.i %*% t(Z[k,]) %*% Z[k,] %*% G.i) * Fij[k]}
  S_pdc.int1<-rowSums(S_pdc.int1)
  S_pdc<- -(diag(x = 1, ncol(Z)) %x% Hi) %*% (vec(G.i) + vec(c_t %*% t(c_t))) + 
    vec(H.i %*% B.i %*% t(Z_at) %*% Sigma.it %*% A_a %*% t(c_t)) + t(Si) %*% Eij +
    c(2/P) * (S_pdc.int1 - (15/36) * vec(Hi %*% G.i %*% Ki %*% t(Ki) %*% G.i))
  
  S_h.int1<-matrix(0,nrow=ncol(Z_at)*2,ncol=ni)
  for (k in 1:nrow(Z)) {
    S_h.int1[,k]<- vec((diag(x = 1, ncol(Z_at)) - B.i %*% H.i) %*% Psi_dc %*% Psi_cc.i %*% G.i %*% t(Z[k,]) %*% Z[k,] %*% G.i %*%
                         t(Psi_dc %*% Psi_cc.i) %*% t(Z_at) %*% Sigma.it %*% Z_at %*% B.i) *  Fij[k]}
  S_h.int1<-rowSums(S_h.int1)
  S_h<-t((-H.i %x% H.i) %*% duplication.matrix(n=ncol(Z_at))) %*%
    ( -(1/2)*((B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% Psi_dc %*%  Psi_cc.i %x% (Psi_dc %*% Psi_cc.i))-
                (B.i %*% H.i %*% Psi_dc %*% Psi_cc.i %x%(B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% Psi_dc %*% Psi_cc.i))) %*% vec(G.i) +
        (1/2)*vec(H - B.i - B.i %*% t(Z_at) %*% Sigma.it %*% A_a %*% t(A_a) %*% Sigma.it %*% Z_at %*% B.i +
                    c(m) * B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% mu_1 %*% t(mu_1) %*% H.i %*% B.i - c(m)*B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% mu_1 %*% t(mu_1)+
                    (1-c(m)) * B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% mu_2 %*% t(mu_2) %*% H.i %*% B.i - (1-c(m))*B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% mu_2%*% t(mu_2)- 
                    Psi_dc %*% Psi_cc.i %*% c_t %*% t(c_t) %*% t(Psi_dc %*% Psi_cc.i) %*% t(Z_at) %*% Sigma.it %*% Z_at %*% B.i +
                    2*B.i %*% t(Z_at) %*% Sigma.it %*% A_a %*% t(c_t) %*% t(Psi_dc %*% Psi_cc.i) +
                    (B.i %*% H.i %*% Psi_dc %*% Psi_cc.i %*% c_t %*% t(c_t) %*% t(Psi_dc %*% Psi_cc.i) %*% t(Z_at) %*% Sigma.it %*% Z_at %*% B.i) - 
                    2*B.i %*%  t(Z_at) %*% Sigma.it %*% A_a %*% t(c_t) %*% t(Psi_dc %*% Psi_cc.i) %*% H.i %*% B.i) +
        t(Mi)%*%Eij + 
        c(1/P)*(S_h.int1 - (15/72)*vec((diag(x = 1, ncol(Z_at))-B.i %*% H.i) %*% Psi_dc %*% Psi_cc.i %*% G.i %*% Ki %*% t(Ki) %*% G.i %*%
                                         t(Psi_dc %*% Psi_cc.i) %*% t(Z_at) %*% Sigma.it %*% Z_at %*% B.i) ) )
  
  S_sigma.int1<-matrix(0,nrow=(dim(Sigma.it)[1]*dim(Sigma.it)[2]),ncol=nrow(Z))
  for (k in 1:nrow(Z)) {
    S_sigma.int1[,k]<- vec(Z_at %*% B.i %*% H.i %*% (Psi_dc %*% Psi_cc.i) %*% G.i %*% t(Z[k,]) %*% Z[k,] %*% G.i %*% t(Psi_dc %*% Psi_cc.i) %*% 
                             (t(Z_at)-t(Z_at) %*% Sigma.it %*% Z_at %*% B.i %*% t(Z_at))) * Fij[k]}
  S_sigma.int1<-rowSums(S_sigma.int1)
  S_sigma.int2<- ((1/2)*vec(- Z_at %*% B.i %*% t(Z_at) - A_a %*% t(A_a) + 
                              (2*A_a %*% t(A_a) %*% Sigma.it %*% Z_at %*% B.i %*% t(Z_at)) -
                              Z_at %*% B.i %*% t(Z_at) %*% Sigma.it  %*% A_a %*% t(A_a) %*% Sigma.it %*% Z_at %*% B.i %*% t(Z_at))-
                    (c(m)/2)*vec(Z_at %*% B.i %*% H.i %*%  mu_1 %*% t(mu_1) %*% t(Z_at)- Z_at %*% B.i %*% H.i %*%  mu_1 %*% t(mu_1) %*% t(Z_at) %*% Sigma.it %*% Z_at %*% B.i %*% t(Z_at))-
                    ((1-c(m))/2)*vec(Z_at %*% B.i %*% H.i %*%  mu_2 %*% t(mu_2) %*% t(Z_at)- Z_at %*% B.i %*% H.i %*%  mu_2 %*% t(mu_2) %*% t(Z_at) %*% Sigma.it %*% Z_at %*% B.i %*% t(Z_at))+  
                    t(Pi) %*% Eij - 
                    (1/2)* (((Z_at %*% B.i  %*% H.i %*% Psi_dc %*% Psi_cc.i) %x% (Z_at %*% Psi_dc %*% Psi_cc.i))-((Z_at %*% B.i %*% H.i %*% Psi_dc %*% Psi_cc.i) %x% (Z_at %*% B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% Psi_dc %*% Psi_cc.i))) %*% vec(G.i)+
                    (1/2)*vec(-(Z_at %*% B.i %*% H.i %*% Psi_dc %*% Psi_cc.i %*% c_t %*% t(c_t) %*% t(Psi_dc %*% Psi_cc.i) %*% t(Z_at))+
                                (Z_at %*% B.i %*% H.i %*% Psi_dc %*% Psi_cc.i %*% c_t %*% t(c_t) %*% t(Psi_dc %*% Psi_cc.i) %*% t(Z_at) %*% Sigma.it %*% Z_at %*% B.i %*% t(Z_at))+
                                2*(A_a %*% t(c_t) %*% t(Psi_dc %*% Psi_cc.i) %*% H.i %*% B.i %*% t(Z_at))-
                                2*(Z_at %*% B.i %*% t(Z_at) %*% Sigma.it %*% A_a %*% t(c_t) %*% t(Psi_dc %*% Psi_cc.i) %*% H.i %*% B.i %*% t(Z_at)))+
                    (1/c(P))*(S_sigma.int1 - (15/72)*vec(Z_at %*% B.i %*% H.i %*% (Psi_dc %*% Psi_cc.i) %*% G.i %*% Ki %*% t(Ki) %*% G.i %*% t(Psi_dc %*% Psi_cc.i) %*% 
                                                           (t(Z_at)-t(Z_at) %*% Sigma.it %*% Z_at %*% B.i %*% t(Z_at)))) )
  Q_sigma<-(diag(2) %x% commutation.matrix(ni, 2) %x% diag(ni)) %*% (diag(2) %x%  diag(2) %x% vec(diag(ni)))
  S_sigma.int3<- t((-Sigma.i %x% Sigma.i) %*% duplication.matrix(n=ncol(Z_at))) %*% ((n_a/2)*vec(Sigma) )
  S_sigma<- S_sigma.int3 - t(duplication.matrix(n=2)) %*% t(Q_sigma) %*% (t(Sigma.it) %x% Sigma.it) %*% S_sigma.int2
  
  
  S<-as.matrix(c(S_b=S_b,S_g=S_g,S_mu=S_mu,S_sigma=S_sigma,S_pcc=S_pcc,S_pdc=S_pdc,S_h=S_h),nrow=23,ncol=1)
  S_0<-matrix(0,nrow=23,ncol=1) #the future change this to nrow related with S
  for(l in 1:23){
    S_0[l,1]<-ifelse(is.na(S[l,1])==T,0,S[l,1])
    #S_0[l,1]<-ifelse(abs(S[l,1])<1e-4,0,S[l,1])
  }

  S_0m<-S_0 %*% t(S_0)
  S_0m[is.infinite(S_0m)] <- 0
  scorei<-t(ginv(S_0m) %*% S_0)
  max.vec<-list(m=m,scorei=scorei)
  return(max.vec)
}  

lklhd.dev<-function(id,datafile,initial,X.var=X.var,Z.var=Z.var,X_a.var=X_a.var,Z_a.var=Z_a.var,U.var=U.var,V1=V1,V2=V2){    
  #Deviance
  Beta=initial$Beta
  Gamma=initial$Gamma
  Gamma<-matrix(Gamma,ncol=1)
  Psi_cc=initial$Psi_cc
  Psi_dd=initial$Psi_dd
  corr<-initial$corr
  Var.V1<-initial$Var.V1
  Var.V2<-initial$Var.V2
  Psi_dc=initial$Psi_dc
  m=initial$m
  mu_1=initial$mu_1
  
  p<-length(Beta)
  q<-ifelse(length(Psi_cc)>1,dim(Psi_cc)[1],length(Psi_cc))
  p_a<-length(Gamma)
  q_a<-ifelse(length(Psi_dd)>1,dim(Psi_dd)[1],length(Psi_dd))
  Sigma<-matrix(c(Var.V1,rep(corr*sqrt(Var.V1)*sqrt(Var.V2),2),Var.V2),nrow = 2,ncol = 2)
  
  #Data
  base<-datafile[datafile$patch==id,]
  base$int <- rep(1,nrow(base))
  U<-base[,U.var]
  X.var<-c("int",X.var)
  X<-base[,X.var]
  X<-as.matrix(X)
  Z.var<-c("int",Z.var)
  Z<-base[,Z.var]
  Z<-as.matrix(Z)
  ni<-nrow(Z)
  
  X_a.var<-c("int",X_a.var)
  X_a<-base[,X_a.var]
  X_a<-as.matrix(X_a)
  X_a<-as.matrix(U * X_a)
  
  Z_a.var<-c("int",Z_a.var)
  Z_a<-base[,Z_a.var]
  Z_a<-as.matrix(Z_a)
  Z_a<-as.matrix(U * Z_a)
  n_a<- sum(Z_a[,1]) 
  V_a<-c(base[,V1],base[,V2])  
  
  
  #Calculations
  c<-as.matrix(rep(0,length(Z.var)),ncol=1)
  mu_2<- -(m/(1-m))*mu_1
  Psi_cc.i<-solve(Psi_cc)
  Psi_cd<-t(Psi_dc)
  H<-Psi_dd-(Psi_dc %*% Psi_cc.i %*% Psi_cd )
  H.i<-ginv(H)
  X_at<-diag(2) %x% as.matrix(X_a)
  Z_at<-diag(2) %x% as.matrix(Z_a)
  Sigma.i<-ginv(Sigma)
  Sigma.t<-Sigma %x% diag(nrow(X_a))
  Sigma.it<-Sigma.i %x% diag(nrow(X_a))
  Sigma.vech<-vech(Sigma)
  
  Psi_cc.vech<-vech(Psi_cc)
  Psi_ddcci<-Psi_dc %*% solve(Psi_cc)  #Check this in the future (when more random effects, if still a vector)
  H.vech<-vech(H)
  
  A_a<-matrix(V_a-X_at %*% Gamma)
  
  tempcheck<-ifelse(sum(V_a==0)==length(V_a),0,1) #if sum is zero, the all observations are zero
  if(tempcheck==1){
    B<-(t(Z_at) %*% Sigma.it %*%  Z_at) + H.i
    B.i<-ginv(B)
  } else {
    B<- H.i
    B.i<-H
  }
  
  E.T<-t(A_a) %*% Sigma.it %*% Z_at %*% B.i %*% H.i %*% Psi_dc %*% Psi_cc.i 
  M_1.T<-m*t(mu_1)
  M_2.T<-(1-m)*t(mu_2)
  M.T<-M_1.T+M_2.T  
  D<- t(Psi_dc %*% Psi_cc.i) %*% H.i %*% B.i %*% t(Z_at) %*%  Sigma.it %*% Z_at %*% Psi_dc %*% Psi_cc.i
  g_1.T<-E.T+(M_1.T %*% H.i %*% B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% Psi_dc %*% Psi_cc.i) 
  g_2.T<-E.T+(M_2.T %*% H.i %*% B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% Psi_dc %*% Psi_cc.i) 
  g.T <- E.T+(M.T   %*% H.i %*% B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% Psi_dc %*% Psi_cc.i) 
  
  #max part for other par
  other<-list(U=U, X=X, Z=Z, Beta=Beta, Psi_cc.i = Psi_cc.i, D=D, g.To=E.T)
  c_t<- laplace2(lap.for.c,c,other)$mode 
  
  eta<- (X %*% Beta) + (Z %*% c_t)
  propi<-exp(eta)/(1+exp(eta))
  propi<-apply(propi, 1, p.inf)
  eta<- log(propi/(1-propi))
  lu<-sum((U * eta) + log(1-propi))
  qua<-t(c_t-(solve(Psi_cc.i+D)%*%t(g.T))) %*% (Psi_cc.i+D) %*% (c_t-(solve(Psi_cc.i+D)%*%t(g.T)))
  f_v2<-lu-((1/2)*qua)
  W_t<- diag(c(propi),ncol=length(propi)) %*% diag(c(1-propi),ncol=length(propi))
  U_a<- (diag(1/diag(W_t)) %*% (U - propi)) + eta
  G<-(t(Z) %*% W_t %*% Z ) + (Psi_cc.i+D)
  G.i<-solve(G)
  
  m3_t<-diag(W_t) * (1-2*propi)
  m4_t<-diag(W_t %*% (diag(dim(W_t)[1])-6*W_t))
  m5_t<-m3_t*diag(diag(dim(W_t)[1])-12*W_t)
  m6_t<-m4_t*diag(diag(dim(W_t)[1])-12*W_t)-(12*m3_t^2)             
  m7_t<-m5_t*diag(diag(dim(W_t)[1])-12*W_t)-(36*m3_t*m4_t) 
  
  P1<- (1/8)* m4_t %*% diag(Z %*% G.i %*% t(Z))^2
  P2<- (1/48)*t(m6_t) %*% diag(Z %*% G.i %*% t(Z))^3
  P3<- (15/72) * t(colSums(Z *  m3_t * diag(Z %*% G.i %*% t(Z)))) %*% G.i %*% colSums(Z *  m3_t * diag(Z %*% G.i %*% t(Z)))  #Check that this runs well with more random efects
  P<- (1 - P1 - P2 + P3)
  P<- ifelse(P>1e-25,P,0.5)

  #loglike
  lf<-  -(n_a+(q/2))*log(2*pi) - (n_a/2)*log(det(Sigma)) - (1/2)*log(det(H)) + (1/2)*log(det(B.i)) - (1/2)*log(det(Psi_cc)) -
    ((1/2)* t(A_a) %*% (Sigma.it - Sigma.it %*% Z_at %*% B.i %*% t(Z_at) %*% Sigma.it) %*% A_a) + 
    (t(A_a) %*% Sigma.it %*% Z_at %*% B.i %*% H.i %*% t(M.T)) -
    ((1/2)*(m)*t(mu_1) %*% H.i %*% B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% mu_1) - 
    ((1/2)*(1-m)*t(mu_2) %*% H.i %*% B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% mu_2) + 
    ((1/2)* g.T %*% solve(Psi_cc.i+D) %*% t(g.T)) + f_v2 + (1/2)*log(det(G))+log(P) 
  
  return(dev=lf)
}  

lklhd.var<-function(id,datafile,est,X.var=X.var,Z.var=Z.var,X_a.var=X_a.var,Z_a.var=Z_a.var,U.var=U.var,V1=V1,V2=V2){    
  #print(id)
  initial<-est
  Beta=initial$Beta
  Gamma=initial$Gamma
  Gamma<-matrix(Gamma,ncol=1)
  Psi_cc=initial$Psi_cc
  Psi_dd=initial$Psi_dd
  corr<-initial$corr
  Var.V1<-initial$Var.V1
  Var.V2<-initial$Var.V2
  Psi_dc=initial$Psi_dc
  m=initial$m
  mu_1=initial$mu_1
  
  p<-length(Beta)
  q<-ifelse(length(Psi_cc)>1,dim(Psi_cc)[1],length(Psi_cc))
  p_a<-length(Gamma)
  q_a<-ifelse(length(Psi_dd)>1,dim(Psi_dd)[1],length(Psi_dd))
  Sigma<-matrix(c(Var.V1,rep(corr*sqrt(Var.V1)*sqrt(Var.V2),2),Var.V2),nrow = 2,ncol = 2)
  
  #Data
  base<-datafile[datafile$patch==id,]
  base$int <- rep(1,nrow(base))
  U<-base[,U.var]
  X.var<-c("int",X.var)
  X<-base[,X.var]
  X<-as.matrix(X)
  Z.var<-c("int",Z.var)
  Z<-base[,Z.var]
  Z<-as.matrix(Z)
  ni<-nrow(Z)
  
  X_a.var<-c("int",X_a.var)
  X_a<-base[,X_a.var]
  X_a<-as.matrix(X_a)
  X_a<-as.matrix(U * X_a)
  
  Z_a.var<-c("int",Z_a.var)
  Z_a<-base[,Z_a.var]
  Z_a<-as.matrix(Z_a)
  Z_a<-as.matrix(U * Z_a)
  n_a<- sum(Z_a[,1]) 
  V_a<-c(base[,V1],base[,V2])  
  
  
  #Calculations
  c<-as.matrix(rep(0,length(Z.var)),ncol=1)
  mu_2<- -(m/(1-m))*mu_1
  Psi_cc.i<-solve(Psi_cc)
  Psi_cd<-t(Psi_dc)
  H<-Psi_dd-(Psi_dc %*% Psi_cc.i %*% Psi_cd )
  H.i<-ginv(H)
  X_at<-diag(2) %x% as.matrix(X_a)
  Z_at<-diag(2) %x% as.matrix(Z_a)
  Sigma.i<-ginv(Sigma)
  Sigma.t<-Sigma %x% diag(nrow(X_a))
  Sigma.it<-Sigma.i %x% diag(nrow(X_a))
  Sigma.vech<-vech(Sigma)
  
  Psi_cc.vech<-vech(Psi_cc)
  Psi_ddcci<-Psi_dc %*% solve(Psi_cc)  #Check this in the future (when more random effects, if still a vector)
  H.vech<-vech(H)
  
  A_a<-matrix(V_a-X_at %*% Gamma)
  
  tempcheck<-ifelse(sum(V_a==0)==length(V_a),0,1) #if sum is zero, the all observations are zero
  if(tempcheck==1){
    B<-(t(Z_at) %*% Sigma.it %*%  Z_at) + H.i
    B.i<-ginv(B)
  } else {
    B<- H.i
    B.i<-H
  }
  
  E.T<-t(A_a) %*% Sigma.it %*% Z_at %*% B.i %*% H.i %*% Psi_dc %*% Psi_cc.i 
  M_1.T<-m*t(mu_1)
  M_2.T<-(1-m)*t(mu_2)
  M.T<-M_1.T+M_2.T  
  D<- t(Psi_dc %*% Psi_cc.i) %*% H.i %*% B.i %*% t(Z_at) %*%  Sigma.it %*% Z_at %*% Psi_dc %*% Psi_cc.i
  g_1.T<-E.T+(M_1.T %*% H.i %*% B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% Psi_dc %*% Psi_cc.i) 
  g_2.T<-E.T+(M_2.T %*% H.i %*% B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% Psi_dc %*% Psi_cc.i) 
  g.T <- E.T+(M.T   %*% H.i %*% B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% Psi_dc %*% Psi_cc.i) 
  
  #max part for other par
  other<-list(U=U, X=X, Z=Z, Beta=Beta, Psi_cc.i = Psi_cc.i, D=D, g.To=E.T)
  c_t<- laplace2(lap.for.c,c,other)$mode 
  
  eta<- (X %*% Beta) + (Z %*% c_t)
  propi<-exp(eta)/(1+exp(eta))
  propi<-apply(propi, 1, p.inf)
  eta<- log(propi/(1-propi))
  lu<-sum((U * eta) + log(1-propi))
  qua<-t(c_t-(solve(Psi_cc.i+D)%*%t(g.T))) %*% (Psi_cc.i+D) %*% (c_t-(solve(Psi_cc.i+D)%*%t(g.T)))
  f_v2<-lu-((1/2)*qua)
  W_t<- diag(c(propi),ncol=length(propi)) %*% diag(c(1-propi),ncol=length(propi))
  U_a<- (diag(1/diag(W_t)) %*% (U - propi)) + eta
  G<-(t(Z) %*% W_t %*% Z ) + (Psi_cc.i+D)
  G.i<-solve(G)
  
  m3_t<-diag(W_t) * (1-2*propi)
  m4_t<-diag(W_t %*% (diag(dim(W_t)[1])-6*W_t))
  m5_t<-m3_t*diag(diag(dim(W_t)[1])-12*W_t)
  m6_t<-m4_t*diag(diag(dim(W_t)[1])-12*W_t)-(12*m3_t^2)             
  m7_t<-m5_t*diag(diag(dim(W_t)[1])-12*W_t)-(36*m3_t*m4_t) 
  
  P1<- (1/8)* m4_t %*% diag(Z %*% G.i %*% t(Z))^2
  P2<- (1/48)*t(m6_t) %*% diag(Z %*% G.i %*% t(Z))^3
  P3<- (15/72) * t(colSums(Z *  m3_t * diag(Z %*% G.i %*% t(Z)))) %*% G.i %*% colSums(Z *  m3_t * diag(Z %*% G.i %*% t(Z)))  #Check that this runs well with more random efects
  P<- (1 - P1 - P2 + P3)
  P<- ifelse(P>1e-25,P,0.5)
  
 
  #Scoring calculations
  Ki<-matrix(0,nrow=ncol(Z) ,ncol=ncol(Z) )
  for (k in 1:nrow(Z)) {
    Ki.temp<-t(Z[k,]) %*% m3_t[k] %*% Z[k,] %*% G.i %*%  t(Z[k,])
    Ki=Ki+Ki.temp }
  hi<-matrix(0,nrow=ncol(Z),ncol =ncol(Z) )
  for (k in 1:nrow(Z)) {
    h.temp<-m3_t[k] %*% G.i %*%t(Z[k,]) %*% (Z[k,] %*% G.i) %*% (Z[k,] %*% G.i %*% Ki)
    hi=hi+h.temp }
  fi<-matrix(0,nrow=ncol(Z),ncol =ncol(Z) )
  for (k in 1:nrow(Z)) {
    fi.temp<-m6_t[k] %*% G.i %*%t(Z[k,]) %*% ((Z[k,] %*% G.i %*% t(Z[k,]))^2) %*% (Z[k,] %*% G.i)
    fi=fi+fi.temp }
  Fi<-matrix(0,nrow=ncol(Z),ncol =ncol(Z) )
  for (k in 1:nrow(Z)) {
    Fi.temp<-m4_t[k] %*% G.i %*% t(Z[k,]) %*% (Z[k,] %*% G.i %*% t(Z[k,])) %*% (Z[k,] %*% G.i)
    Fi=Fi+Fi.temp }
  Rij<-matrix(0,nrow=nrow(Z),ncol=1)
  for (k in 1:nrow(Z)) {
    Rij[k,1] <-  -((1/8)*m5_t[k] %*% (Z[k,] %*% G.i %*% t(Z[k,]))^2) + ((1/4)* m3_t[k] %*% Z[k,] %*% Fi %*% t(Z[k,])) - 
      ((1/48)*m7_t[k] %*% (Z[k,] %*% G.i %*% t(Z[k,]))^3) +  ((1/16)* m3_t[k] %*% Z[k,] %*% fi %*% t(Z[k,])) - 
      ((15/72)*m3_t[k] %*% (Z[k,] %*% G.i %*% Ki)^2) + ((15/36)* m4_t[k] %*% (t(Ki) %*% G.i %*% t(Z[k,])) %*% (Z[k,] %*% G.i %*% t(Z[k,]))) -
      ((15/36)*m3_t[k] %*% Z[k,] %*% hi %*% t(Z[k,]))
  }
  Fij<-matrix(0,nrow=nrow(Z),ncol=1)
  for (k in 1: nrow(Z)) {
    Fij[k,1]<- ((1/4)*m4_t[k] %*% (Z[k,] %*% G.i %*% t(Z[k,]))) + ((1/16)*m6_t[k] %*% (Z[k,] %*% G.i %*% t(Z[k,]))^2) -
      ((15/36)*m3_t[k] %*% (Z[k,] %*% G.i %*% Ki))      }
  
  Yi<- t(Z) %*% W_t %*% (U_a - (eta-(Z %*% c_t))) + t(E.T)
  Hi<- t(Z_at) %*% Sigma.it %*% Z_at %*% B.i %*% H.i %*% Psi_dc %*% Psi_cc.i
  
  Eij.int<-matrix(0,nrow=nrow(Z),ncol=1)
  for (k in 1:nrow(Z)) {
    Eij.int[k,1]<- t(Z[k,]) %*% (((-1/2)*m3_t[k] %*% (Z[k,] %*% G.i %*% t(Z[k,]))) + (U[k] - propi[k]) + (Rij[k]/P) )}
  Eij<- t(E.T)-(Psi_cc.i+D)%*%c_t +colSums(Eij.int)
  
  Si<-((t(A_a) %*% Sigma.it %*% Z_at %*% B.i %*% H.i) %x% G.i)-(G.i %x% t(Yi) %*% G.i %*% t(Hi))-(t(Yi) %*% G.i %x% G.i %*% t(Hi))
  ti<- -G.i %*% t(Psi_dc %*% Psi_cc.i) %*% H.i %*% B.i %*% t(Z_at) %*% Sigma.it %*% X_at
  Di<- -(t(Yi) %*% G.i %x% G.i)
  
  Mi<- -(t(Yi) %*% G.i %*% t(Psi_dc %*% Psi_cc.i) %*% t(Z_at) %*% Sigma.it %*% Z_at %*% B.i %x% (G.i %*% t(Psi_dc %*% Psi_cc.i)))+
    (t(A_a) %*% Sigma.it %*% Z_at %*% B.i %x% (G.i %*% t(Psi_dc %*% Psi_cc.i)))-
    (t(A_a) %*% Sigma.it %*% Z_at %*% B.i %x% (G.i %*% t(Psi_dc %*% Psi_cc.i) %*% H.i %*% B.i))+
    (t(Yi) %*% G.i %*% t(Psi_dc %*% Psi_cc.i) %*% t(Z_at) %*% Sigma.it %*% Z_at %*% B.i %x% (G.i %*% t(Psi_dc %*% Psi_cc.i) %*% H.i %*% B.i)) 
  
  Pi<- (t(Yi) %*% G.i %*% t(Psi_dc %*% Psi_cc.i) %*% t(Z_at) %*% Sigma.it %*% Z_at %*% B.i %*% t(Z_at) %x% (G.i %*% t(Psi_dc %*% Psi_cc.i) %*% H.i %*% B.i %*% t(Z_at)))-
    (t(Yi) %*% G.i %*% t(Psi_dc %*% Psi_cc.i) %*% t(Z_at) %x% (G.i %*% t(Psi_dc %*% Psi_cc.i) %*% H.i %*% B.i %*% t(Z_at)))-
    (t(A_a) %*% Sigma.it %*% Z_at %*% B.i %*% t(Z_at) %x% (G.i %*% t(Psi_dc %*% Psi_cc.i) %*% H.i %*% B.i %*% t(Z_at)))+
    (t(A_a) %x% (G.i %*% t(Psi_dc %*% Psi_cc.i)%*% H.i %*% B.i %*% t(Z_at)))
  
  #Score vectors
  S_b.int1<-matrix(0,nrow=nrow(Z),ncol=ncol(X))
  for (k in 1:nrow(Z)) {
    S_b.int1[k,]<-c(((-1/2)*m3_t[k] %*% (Z[k,] %*% G.i %*% t(Z[k,])))) * (t(X[k,]) - t(t(X) %*% W_t %*% Z %*% G.i %*% t(Z[k,])))}
  S_b.int1<-t(S_b.int1)
  S_b.int1<-rowSums(S_b.int1)
  S_b.int2<-matrix(0,nrow=nrow(Z),ncol=ncol(X))
  for (k in 1:nrow(Z)) {
    S_b.int2[k,]<-Rij[k]*(t(X[k,]) - t(t(X) %*% W_t %*% Z %*% G.i %*% t(Z[k,])))}
  S_b.int2<-t(S_b.int2)
  S_b.int2<-rowSums(S_b.int2)
  S_b<- t(X) %*% W_t %*% (U_a - eta) + S_b.int1 + (c((1/P)) * S_b.int2)
  rule1<-sum(abs(S_b)<1e-1)
  if(rule1>0){
    S_b<-matrix(0,nrow=ncol(X),1)
  }
  
  S_g<-(t(X_at) %*% Sigma.it %*% V_a) - (t(X_at) %*% Sigma.it %*%  Z_at %*% B.i %*% t(Z_at) %*%  Sigma.it %*% V_a) -
    (t(X_at) %*% Sigma.it %*% X_at %*% Gamma) + (t(X_at) %*% Sigma.it %*%  Z_at %*% B.i %*% t(Z_at) %*%  Sigma.it %*% X_at %*% Gamma) -
    (t(X_at) %*% Sigma.it %*%  Z_at %*% B.i %*% H.i %*% Psi_dc %*% Psi_cc.i %*% c_t) + (t(ti) %*% Eij)
  
  mu_1<-matrix(mu_1,ncol=1)
  S_mu.int1<- (H.i %*% B.i %*% t(Z_at) %*% Sigma.it %*% Z_at)
  S_mu.int2<- (t(Z_at) %*% Sigma.it %*% Z_at %*% B.i %*% H.i)
  S_mu.int3<- (S_mu.int1 + S_mu.int2) %*% mu_1
  S_mu.int4<- -(1/2)*(m+((m^2)/(1-m)))
  S_mu<- c(S_mu.int4) * S_mu.int3
  
  S_pcc.int1<-ifelse(q==1,t((Psi_cc.i %x% Psi_cc.i)),t((Psi_cc.i %x% Psi_cc.i) %*% duplication.matrix(n=q)))
  S_pcc.int2<-matrix(0,nrow=ncol(Z)*ncol(Z),ncol=ni)
  for (k in 1:nrow(Z)) {
    S_pcc.int2[,k]<- vec(G.i %*% t(Z[k,]) %*% Z[k,] %*% G.i) *  Fij[k]}
  S_pcc.int2<-rowSums(S_pcc.int2)
  S_pcc<- -S_pcc.int1 %*% ( (1/2)*vec(Psi_cc - G.i - c_t %*% t(c_t)) + t(Di) %*% Eij + 
                              (1/P) %*% ( S_pcc.int2 - (15/72)*vec(G.i %*% Ki %*% t(Ki) %*% G.i ) ) ) 
  S_pcc<- as.matrix(apply(S_pcc, 1, function(x) {ifelse(abs(x)<1e-1 & x!=0, 0, x)}),ncol=1)
  
  S_pdc.int1<-matrix(0,nrow=ncol(Z_at)*ncol(Z),ncol=ni)
  for (k in 1:nrow(Z)) {
    S_pdc.int1[,k]<- vec(Hi %*% G.i %*% t(Z[k,]) %*% Z[k,] %*% G.i) * Fij[k]}
  S_pdc.int1<-rowSums(S_pdc.int1)
  S_pdc<- -(diag(x = 1, ncol(Z)) %x% Hi) %*% (vec(G.i) + vec(c_t %*% t(c_t))) + 
    vec(H.i %*% B.i %*% t(Z_at) %*% Sigma.it %*% A_a %*% t(c_t)) + t(Si) %*% Eij +
    c(2/P) * (S_pdc.int1 - (15/36) * vec(Hi %*% G.i %*% Ki %*% t(Ki) %*% G.i))
  
  S_h.int1<-matrix(0,nrow=ncol(Z_at)*2,ncol=ni)
  for (k in 1:nrow(Z)) {
    S_h.int1[,k]<- vec((diag(x = 1, ncol(Z_at)) - B.i %*% H.i) %*% Psi_dc %*% Psi_cc.i %*% G.i %*% t(Z[k,]) %*% Z[k,] %*% G.i %*%
                         t(Psi_dc %*% Psi_cc.i) %*% t(Z_at) %*% Sigma.it %*% Z_at %*% B.i) *  Fij[k]}
  S_h.int1<-rowSums(S_h.int1)
  S_h<-t((-H.i %x% H.i) %*% duplication.matrix(n=ncol(Z_at))) %*%
    ( -(1/2)*((B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% Psi_dc %*%  Psi_cc.i %x% (Psi_dc %*% Psi_cc.i))-
                (B.i %*% H.i %*% Psi_dc %*% Psi_cc.i %x%(B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% Psi_dc %*% Psi_cc.i))) %*% vec(G.i) +
        (1/2)*vec(H - B.i - B.i %*% t(Z_at) %*% Sigma.it %*% A_a %*% t(A_a) %*% Sigma.it %*% Z_at %*% B.i +
                    c(m) * B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% mu_1 %*% t(mu_1) %*% H.i %*% B.i - c(m)*B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% mu_1 %*% t(mu_1)+
                    (1-c(m)) * B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% mu_2 %*% t(mu_2) %*% H.i %*% B.i - (1-c(m))*B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% mu_2%*% t(mu_2)- 
                    Psi_dc %*% Psi_cc.i %*% c_t %*% t(c_t) %*% t(Psi_dc %*% Psi_cc.i) %*% t(Z_at) %*% Sigma.it %*% Z_at %*% B.i +
                    2*B.i %*% t(Z_at) %*% Sigma.it %*% A_a %*% t(c_t) %*% t(Psi_dc %*% Psi_cc.i) +
                    (B.i %*% H.i %*% Psi_dc %*% Psi_cc.i %*% c_t %*% t(c_t) %*% t(Psi_dc %*% Psi_cc.i) %*% t(Z_at) %*% Sigma.it %*% Z_at %*% B.i) - 
                    2*B.i %*%  t(Z_at) %*% Sigma.it %*% A_a %*% t(c_t) %*% t(Psi_dc %*% Psi_cc.i) %*% H.i %*% B.i) +
        t(Mi)%*%Eij + 
        c(1/P)*(S_h.int1 - (15/72)*vec((diag(x = 1, ncol(Z_at))-B.i %*% H.i) %*% Psi_dc %*% Psi_cc.i %*% G.i %*% Ki %*% t(Ki) %*% G.i %*%
                                         t(Psi_dc %*% Psi_cc.i) %*% t(Z_at) %*% Sigma.it %*% Z_at %*% B.i) ) )
  
  S_sigma.int1<-matrix(0,nrow=(dim(Sigma.it)[1]*dim(Sigma.it)[2]),ncol=nrow(Z))
  for (k in 1:nrow(Z)) {
    S_sigma.int1[,k]<- vec(Z_at %*% B.i %*% H.i %*% (Psi_dc %*% Psi_cc.i) %*% G.i %*% t(Z[k,]) %*% Z[k,] %*% G.i %*% t(Psi_dc %*% Psi_cc.i) %*% 
                             (t(Z_at)-t(Z_at) %*% Sigma.it %*% Z_at %*% B.i %*% t(Z_at))) * Fij[k]}
  S_sigma.int1<-rowSums(S_sigma.int1)
  S_sigma.int2<- ((1/2)*vec(- Z_at %*% B.i %*% t(Z_at) - A_a %*% t(A_a) + 
                              (2*A_a %*% t(A_a) %*% Sigma.it %*% Z_at %*% B.i %*% t(Z_at)) -
                              Z_at %*% B.i %*% t(Z_at) %*% Sigma.it  %*% A_a %*% t(A_a) %*% Sigma.it %*% Z_at %*% B.i %*% t(Z_at))-
                    (c(m)/2)*vec(Z_at %*% B.i %*% H.i %*%  mu_1 %*% t(mu_1) %*% t(Z_at)- Z_at %*% B.i %*% H.i %*%  mu_1 %*% t(mu_1) %*% t(Z_at) %*% Sigma.it %*% Z_at %*% B.i %*% t(Z_at))-
                    ((1-c(m))/2)*vec(Z_at %*% B.i %*% H.i %*%  mu_2 %*% t(mu_2) %*% t(Z_at)- Z_at %*% B.i %*% H.i %*%  mu_2 %*% t(mu_2) %*% t(Z_at) %*% Sigma.it %*% Z_at %*% B.i %*% t(Z_at))+  
                    t(Pi) %*% Eij - 
                    (1/2)* (((Z_at %*% B.i  %*% H.i %*% Psi_dc %*% Psi_cc.i) %x% (Z_at %*% Psi_dc %*% Psi_cc.i))-((Z_at %*% B.i %*% H.i %*% Psi_dc %*% Psi_cc.i) %x% (Z_at %*% B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% Psi_dc %*% Psi_cc.i))) %*% vec(G.i)+
                    (1/2)*vec(-(Z_at %*% B.i %*% H.i %*% Psi_dc %*% Psi_cc.i %*% c_t %*% t(c_t) %*% t(Psi_dc %*% Psi_cc.i) %*% t(Z_at))+
                                (Z_at %*% B.i %*% H.i %*% Psi_dc %*% Psi_cc.i %*% c_t %*% t(c_t) %*% t(Psi_dc %*% Psi_cc.i) %*% t(Z_at) %*% Sigma.it %*% Z_at %*% B.i %*% t(Z_at))+
                                2*(A_a %*% t(c_t) %*% t(Psi_dc %*% Psi_cc.i) %*% H.i %*% B.i %*% t(Z_at))-
                                2*(Z_at %*% B.i %*% t(Z_at) %*% Sigma.it %*% A_a %*% t(c_t) %*% t(Psi_dc %*% Psi_cc.i) %*% H.i %*% B.i %*% t(Z_at)))+
                    (1/c(P))*(S_sigma.int1 - (15/72)*vec(Z_at %*% B.i %*% H.i %*% (Psi_dc %*% Psi_cc.i) %*% G.i %*% Ki %*% t(Ki) %*% G.i %*% t(Psi_dc %*% Psi_cc.i) %*% 
                                                           (t(Z_at)-t(Z_at) %*% Sigma.it %*% Z_at %*% B.i %*% t(Z_at)))) )
  Q_sigma<-(diag(2) %x% commutation.matrix(ni, 2) %x% diag(ni)) %*% (diag(2) %x%  diag(2) %x% vec(diag(ni)))
  S_sigma.int3<- t((-Sigma.i %x% Sigma.i) %*% duplication.matrix(n=ncol(Z_at))) %*% ((n_a/2)*vec(Sigma) )
  S_sigma<- S_sigma.int3 - t(duplication.matrix(n=2)) %*% t(Q_sigma) %*% (t(Sigma.it) %x% Sigma.it) %*% S_sigma.int2
  
  
  for(l in 1:length(S_b)){
    S_b[l]<-ifelse(is.na(S_b[l])==T,0,S_b[l])
  }
  for(l in 1:length(S_g)){
    S_g[l]<-ifelse(is.na(S_g[l])==T,0,S_g[l])
  }
  for(l in 1:length(S_mu)){
    S_mu[l]<-ifelse(is.na(S_mu[l])==T,0,S_mu[l])
  }
  
  S_0b <- S_b %*% t(S_b)
  S_0g <- S_g %*% t(S_g)
  S_0mu <- S_mu %*% t(S_mu)
  
  S_0b[is.infinite(S_0b)] <- 0
  S_0g[is.infinite(S_0g)] <- 0
  S_0mu[is.infinite(S_0mu)] <- 0
  
  FI<- list(S_0b=S_0b,S_0g=S_0g,S_0mu=S_0mu)
  return(FI)
} 

lklhd.ram<-function(id,datafile,est,X.var=X.var,Z.var=Z.var,X_a.var=X_a.var,Z_a.var=Z_a.var,U.var=U.var,V1=V1,V2=V2,sim=10000){    
  #Estimation of random effects 
  print(id)
  initial<-est
  Beta=initial$Beta
  Gamma=initial$Gamma
  Gamma<-matrix(Gamma,ncol=1)
  Psi_cc=initial$Psi_cc
  Psi_dd=initial$Psi_dd
  corr<-initial$corr
  Var.V1<-initial$Var.V1
  Var.V2<-initial$Var.V2
  Psi_dc=initial$Psi_dc
  m=initial$m
  mu_1=initial$mu_1
  
  p<-length(Beta)
  q<-ifelse(length(Psi_cc)>1,dim(Psi_cc)[1],length(Psi_cc))
  p_a<-length(Gamma)
  q_a<-ifelse(length(Psi_dd)>1,dim(Psi_dd)[1],length(Psi_dd))
  Sigma<-matrix(c(Var.V1,rep(corr*sqrt(Var.V1)*sqrt(Var.V2),2),Var.V2),nrow = 2,ncol = 2)
  
  #Data
  base<-datafile[datafile$patch==id,]
  base$int <- rep(1,nrow(base))
  U<-base[,U.var]
  X.var<-c("int",X.var)
  X<-base[,X.var]
  X<-as.matrix(X)
  Z.var<-c("int",Z.var)
  Z<-base[,Z.var]
  Z<-as.matrix(Z)
  ni<-nrow(Z)
  
  X_a.var<-c("int",X_a.var)
  X_a<-base[,X_a.var]
  X_a<-as.matrix(X_a)
  X_a<-as.matrix(U * X_a)
  
  Z_a.var<-c("int",Z_a.var)
  Z_a<-base[,Z_a.var]
  Z_a<-as.matrix(Z_a)
  Z_a<-as.matrix(U * Z_a)
  n_a<- sum(Z_a[,1]) 
  V_a<-c(base[,V1],base[,V2])  
  
  #Calculations
  mu_2<- -(m/(1-m))*mu_1
  Psi_cc.i<-solve(Psi_cc)
  Psi_cd<-t(Psi_dc)
  H<-Psi_dd-Psi_dc %*% Psi_cc.i %*% Psi_cd  
  H.i<-ginv(H)
  X_at<-diag(2) %x% as.matrix(X_a)
  Z_at<-diag(2) %x% as.matrix(Z_a)
  Sigma.i<-ginv(Sigma)
  Sigma.t<-Sigma %x% diag(nrow(X_a))
  Sigma.it<-Sigma.i %x% diag(nrow(X_a))
  Sigma.vech<-vech(Sigma)
  
  Psi_cc.vech<-vech(Psi_cc)
  Psi_ddcci<-Psi_dc %*% solve(Psi_cc)  #Check this in the future (when more random effects, if still a vector)
  H.vech<-vech(H)
  
  A_a<-matrix(V_a-X_at %*% Gamma)
  
  tempcheck<-ifelse(sum(V_a==0)==length(V_a),0,1) #if sum is zero, the all observations are zero
  if(tempcheck==1){
    B<-t(Z_at) %*% Sigma.it %*%  Z_at + H.i
    B.i<-ginv(B)
  } else {
    B<- H.i
    B.i<-H
  }
  
  E.T<-t(A_a) %*% Sigma.it %*% Z_at %*% B.i %*% H.i %*% Psi_dc %*% Psi_cc.i 
  D<- Psi_cc.i %*% Psi_cd %*% H.i %*% B.i %*% t(Z_at) %*%  Sigma.it %*% Z_at %*% Psi_dc %*% Psi_cc.i
  g_1.T<-E.T+(t(mu_1) %*% H.i %*% B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% Psi_dc %*% Psi_cc.i) 
  g_2.T<-E.T+(t(mu_2) %*% H.i %*% B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% Psi_dc %*% Psi_cc.i) 
  
  c<-as.matrix(rep(0,length(Z.var)),ncol=1)
  other<-list(U=U, X=X, Z=Z, Beta=Beta, Psi_cc.i = Psi_cc.i, D=D, g.To=E.T)
  c_t<- laplace2(lap.for.c,c,other)$mode 
  
  #Here I had to truncate c to get better predictions but I will have to check later
  c_t <- ifelse(c_t>60,60,c_t)
  eta<- (X %*% Beta) + (Z %*% c_t)
  propi<-exp(eta)/(1+exp(eta))
  propi<-apply(propi, 1, p.inf)
  eta<- log(propi/(1-propi))
  W_t<- diag(c(propi),ncol=length(propi)) %*% diag(c(1-propi),ncol=length(propi))
  G<-(t(Z) %*% W_t %*% Z ) + (Psi_cc.i+D)
  G.i.ram<-solve(G)
  
  random<-matrix(data=NA,nrow=sim,ncol=15)
  colnames(random)<-c("c","cc^t","g","q","w","cw","cc^tw","d1.1","d1.2","d2.1","d2.2","d1.1w","d1.2w","d2.1w","d2.2w")
  #random[,1]<-rmvt(sim, sigma = G.i.ram, df = 4, delta = c_t,type = c("shifted"))
  random[,1]<-rep(101,sim)
  for(e in 1:sim){
    while (random[e,1]>100) {
      random[e,1]<-rt(1, df=4, ncp=c_t)  
    }
    #random[e,1]<-rt(1, df=4, ncp=c_t) 
    random[e,2]<-random[e,1]*random[e,1]
    random[e,8:9]<- (mu_1 + (Psi_dc %*% Psi_cc.i * random[e,1])) + (B.i %*% ((t(Z_at) %*% Sigma.it %*% A_a) - (t(Z_at) %*% Sigma.it %*% Z_at %*% Psi_dc %*% Psi_cc.i * random[e,1]) ))
    random[e,10:11]<- (mu_2 + (Psi_dc %*% Psi_cc.i * random[e,1])) + (B.i %*% ((t(Z_at) %*% Sigma.it %*% A_a) - (t(Z_at) %*% Sigma.it %*% Z_at %*% Psi_dc %*% Psi_cc.i * random[e,1]) ))
  }
  
  for(e in 1:sim){
    #random[e,3]<-dmvt(random[e,1],sigma = G.i.ram, df = 4, delta = c_t, type = c("shifted"))
    random[e,3]<-dt(random[e,1], df=4, ncp=c_t)
    random[e,4]<-dmvnorm(random[e,1], mean = rep(0, q), sigma = Psi_cc, log = FALSE)
    
    eta<- (X %*% Beta) + (Z %*% random[e,1])
    propi<-exp(eta)/(1+exp(eta))
    propi<-apply(propi, 1, p.inf)
    eta<- log(propi/(1-propi))
    lu<-sum((U * eta) + log(1-propi))
    
    lf1<- m*exp( (t(A_a) %*% Sigma.it %*% Z_at %*% B.i %*% H.i %*% mu_1) - (1/2)*(t(mu_1) %*% H.i %*% B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% mu_1))*
      exp((E.T %*% random[e,1]) - (t(mu_1) %*% H.i %*% B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% Psi_dc %*% Psi_cc.i %*% random[e,1]) - (1/2)*(t(random[e,1]) %*% D %*% random[e,1]) ) 
    
    lf2<- (1-m)*exp( (t(A_a) %*% Sigma.it %*% Z_at %*% B.i %*% H.i %*% mu_2) - (1/2)*(t(mu_2) %*% H.i %*% B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% mu_2))*
      exp((E.T %*% random[e,1]) - (t(mu_2) %*% H.i %*% B.i %*% t(Z_at) %*% Sigma.it %*% Z_at %*% Psi_dc %*% Psi_cc.i %*% random[e,1]) - (1/2)*(t(random[e,1]) %*% D %*% random[e,1]) ) 
    
    lf<- -(n_a+(q/2))*log(2*pi) - (n_a/2)*log(det(Sigma)) - (1/2)*log(det(H)) + (1/2)*log(det(B.i)) - 
      ((1/2)* t(A_a) %*% (Sigma.it - Sigma.it %*% Z_at %*% B.i %*% t(Z_at) %*% Sigma.it) %*% A_a) + 
      log(lf1 + lf2)
    qu<- exp(lu) * exp(lf)
    
    random[e,4]<- random[e,4]*qu
    random[e,4]<-ifelse(is.finite(random[e,4]),random[e,4],0)
  }
  
  if (sum(random[,4])!=0) {
    minn<-max(min(random[random[,4]!=0,4]),1e-10)
    for(v in 1:sim){
      random[v,4]<-ifelse(random[v,4]<1e-10,minn,random[v,4])
    }
  } else {
    for(v in 1:sim){
      random[v,4]<-dmvnorm(random[v,1], mean = rep(0, q), sigma = Psi_cc, log = FALSE)
      random[v,4]<-random[v,4]*0.01
    }
    minn<-max(min(random[random[,4]!=0,4]),1e-10)
    for(v in 1:sim){
      random[v,4]<-ifelse(random[v,4]<1e-10,minn,random[v,4])
    }
  }
  
  random[,5]<-random[,4]/random[,3]
  
  sumsum<-ifelse(is.infinite(sum(random[,5]*random[,5])) & sum(random[,5]*random[,5])<0 ,0,1)
  while (sumsum<1e-10) {
    random[,4]<-random[,4]*(1e+10)
    random[,5]<-random[,4]/random[,3]     
    sumsum<-ifelse(is.infinite(sum(random[,5]*random[,5])) & sum(random[,5]*random[,5])<0 ,0,1)
  }
  
  sumsum<-ifelse(is.infinite(sum(random[,5]*random[,5])) & sum(random[,5]*random[,5])>0 ,0,1)
  while (sumsum<1e-10) {
    random[,4]<-random[,4]/(1e+10)
    random[,5]<-random[,4]/random[,3] 
    sumsum<-ifelse(is.infinite(sum(random[,5]*random[,5])) & sum(random[,5]*random[,5])>0 ,0,1)
    sumsum
  }
  
  random.t<-random
  random.p<-random
  
  #### estimates without adjustments on weights
  random[,6]<-random[,1]*random[,5]
  random[,7]<-random[,2]*random[,5]
  random[,12]<-random[,8]*random[,5]
  random[,13]<-random[,9]*random[,5]
  random[,14]<-random[,10]*random[,5]
  random[,15]<-random[,11]*random[,5]
  
  random <- random[complete.cases(random), ] 
  c.est <- mean(random[,6])/ mean(random[,5])
  #c.est <- ifelse(c.est < -100,-100,c.est)
  #c.est <- ifelse(c.est > 100,100,c.est)
  cc.t.est <- mean(random[,7])/ mean(random[,5])
  d1.est <- (m*mean(random[,12])/ mean(random[,5])) + ((1-m)*mean(random[,14])/ mean(random[,5]))
  d2.est <- (m*mean(random[,13])/ mean(random[,5])) + ((1-m)*mean(random[,15])/ mean(random[,5]))
  
  mean.w<-mean(random[,5])
  ess1<-((random[,5]/mean.w)-1)^2
  ess1<-sum(ess1)
  ess1<- sqrt(ess1/nrow(random))
  ess2<-(sum(random[,5])^2)/sum(random[,5]*random[,5])
  ram.est<-list(c_t=c_t,c.est=c.est,cc.t.est=cc.t.est,d1.est=d1.est,d2.est=d2.est,ess1=ess1, ess2=ess2)
  
  ###### Truncated weights
  def<-sqrt(sim*mean(random.t[,5]))
  for(rt in 1:sim){
    random.t[rt,5]<-min(random.t[rt,5],def)
  }
  
  
  sumsum<-ifelse(is.infinite(sum(random.t[,5]*random.t[,5])) & sum(random.t[,5]*random.t[,5])<0 ,0,1)
  while (sumsum<1e-10) {
    random.t[,4]<-random.t[,4]*(1e+10)
    random.t[,5]<-random.t[,4]/random.t[,3]     
    sumsum<-ifelse(is.infinite(sum(random.t[,5]*random.t[,5])) & sum(random.t[,5]*random.t[,5])<0 ,0,1)
  }
  
  sumsum<-ifelse(is.infinite(sum(random.t[,5]*random.t[,5])) & sum(random.t[,5]*random.t[,5])>0 ,0,1)
  while (sumsum<1e-10) {
    random.t[,4]<-random.t[,4]/(1e+10)
    random.t[,5]<-random.t[,4]/random.t[,3] 
    sumsum<-ifelse(is.infinite(sum(random.t[,5]*random.t[,5])) & sum(random.t[,5]*random.t[,5])>0 ,0,1)
    sumsum
  }
  
  random.t[,6]<-random.t[,1]*random.t[,5]
  random.t[,7]<-random.t[,2]*random.t[,5]
  random.t[,12]<-random.t[,8]*random.t[,5]
  random.t[,13]<-random.t[,9]*random.t[,5]
  random.t[,14]<-random.t[,10]*random.t[,5]
  random.t[,15]<-random.t[,11]*random.t[,5]
  
  random.t <- random.t[complete.cases(random.t), ] 
  c.est.t <- mean(random.t[,6])/ mean(random.t[,5])
  cc.t.est.t <- mean(random.t[,7])/ mean(random.t[,5])
  d1.est.t <- (m*mean(random.t[,12])/ mean(random.t[,5])) + ((1-m)*mean(random.t[,14])/ mean(random.t[,5]))
  d2.est.t <- (m*mean(random.t[,13])/ mean(random.t[,5])) + ((1-m)*mean(random.t[,15])/ mean(random.t[,5]))
  
  mean.w<-mean(random.t[,5])
  ess1.t<-((random.t[,5]/mean.w)-1)^2
  ess1.t<-sum(ess1.t)
  ess1.t<- sqrt(ess1.t/nrow(random.t))
  ess2.t<-(sum(random.t[,5])^2)/sum(random.t[,5]*random.t[,5])
  ram.est.t<-list(c.est.t=c.est.t,cc.t.est.t=cc.t.est.t,d1.est.t=d1.est.t,d2.est.t=d2.est.t,ess1.t=ess1.t,ess2.t=ess2.t)
  
  
  ###### Pareto weights
  log_ratios <- log(random.p[,5])
  r_eff <- relative_eff(random.p[,5], chain_id = rep(1,sim))
  psis_result <- psis(log_ratios, r_eff = r_eff)
  random.p[,5]<-weights(psis_result, log=FALSE, normalize = FALSE)
  
  sumsum<-ifelse(is.infinite(sum(random.p[,5]*random.p[,5])) & sum(random.p[,5]*random.p[,5])<0 ,0,1)
  while (sumsum<1e-10) {
    random.p[,4]<-random.p[,4]*(1e+10)
    random.p[,5]<-random.p[,4]/random.p[,3]     
    sumsum<-ifelse(is.infinite(sum(random.p[,5]*random.p[,5])) & sum(random.p[,5]*random.p[,5])<0 ,0,1)
  }
  
  sumsum<-ifelse(is.infinite(sum(random.p[,5]*random.p[,5])) & sum(random.p[,5]*random.p[,5])>0 ,0,1)
  while (sumsum<1e-10) {
    random.p[,4]<-random.p[,4]/(1e+10)
    random.p[,5]<-random.p[,4]/random.p[,3] 
    sumsum<-ifelse(is.infinite(sum(random.p[,5]*random.p[,5])) & sum(random.p[,5]*random.p[,5])>0 ,0,1)
    sumsum
  }
  
  random.p[,6]<-random.p[,1]*random.p[,5]
  random.p[,7]<-random.p[,2]*random.p[,5]
  random.p[,12]<-random.p[,8]*random.p[,5]
  random.p[,13]<-random.p[,9]*random.p[,5]
  random.p[,14]<-random.p[,10]*random.p[,5]
  random.p[,15]<-random.p[,11]*random.p[,5]
  
  random.p <- random.p[complete.cases(random.p), ] 
  c.est.p <- mean(random.p[,6])/ mean(random.p[,5])
  cc.t.est.p <- mean(random.p[,7])/ mean(random.p[,5])
  d1.est.p <- (m*mean(random.p[,12])/ mean(random.p[,5])) + ((1-m)*mean(random.p[,14])/ mean(random.p[,5]))
  d2.est.p <- (m*mean(random.p[,13])/ mean(random.p[,5])) + ((1-m)*mean(random.p[,15])/ mean(random.p[,5]))
  
  mean.w<-mean(random.p[,5])
  ess1.p<-((random.p[,5]/mean.w)-1)^2
  ess1.p<-sum(ess1.p)
  ess1.p<- sqrt(ess1.p/nrow(random))
  ess2.p<-(sum(random.p[,5])^2)/sum(random.p[,5]*random.p[,5])
  
  ram.est.p<-list(c.est.p=c.est.p,cc.t.est.p=cc.t.est.p,d1.est.p=d1.est.p,d2.est.p=d2.est.p,ess1.p=ess1.p,ess2.p=ess2.p)
  ram.est.f<-list(ram.est=ram.est, ram.est.t=ram.est.t, ram.est.p=ram.est.p)
  return(ram.est.f)
} 

