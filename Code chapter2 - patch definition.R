#First, export CT scan data from Slicer to Matlab, and from Matlab 
#to R using Kingston Kang's code.

#This code is for data of one patient (c1)

# Packages --------------------------------------------------------------
library(plyr)
library(tidyr)

# Read data -------------------------------------------------------------

files<-list.files(pattern = "*.rda")
dta<- new.env()
lapply(files,load,dta)

#Have to change the names
dta$c7_011.rda<-dta$c7_01_3.rda
dta$c7_012.rda<-dta$c7_01_6.rda
dta$c7_013.rda<-dta$c7_01_12.rda
dta$c7_014.rda<-dta$c7_01_18.rda
dta$c7_015.rda<-dta$c7_01_24.rda
dta$c7_101.rda<-dta$c7_10_3.rda
dta$c7_102.rda<-dta$c7_10_6.rda
dta$c7_103.rda<-dta$c7_10_12.rda
dta$c7_104.rda<-dta$c7_10_18.rda
dta$c7_105.rda<-dta$c7_10_24.rda
dta$c7_111.rda<-dta$c7_11_0.rda
remove("c7_00_0.rda","c7_01_3.rda","c7_00_18.rda","c7_01_24.rda",
       "c7_00_12.rda","c7_10_24.rda","c7_00_6.rda","c7_00_3.rda",
       "c7_01_18.rda","c7_01_12.rda","c7_10_18.rda","c7_00_24.rda",
       "c7_10_12.rda","c7_10_6.rda","c7_10_3.rda","c7_01_6.rda",
       "c7_11_0.rda", envir = dta)

# There is overlap in the images. Hazy has dense data, so next lines 
# remove dense from hazy
dta$c7_011.rda$dvf<-ifelse(dta$c7_011.rda$dvf==dta$c7_101.rda$dvf,
                           0,dta$c7_011.rda$dvf)
dta$c7_012.rda$dvf<-ifelse(dta$c7_012.rda$dvf==dta$c7_102.rda$dvf, 
                           0,dta$c7_012.rda$dvf)
dta$c7_013.rda$dvf<-ifelse(dta$c7_013.rda$dvf==dta$c7_103.rda$dvf, 
                           0,dta$c7_013.rda$dvf)
dta$c7_014.rda$dvf<-ifelse(dta$c7_014.rda$dvf==dta$c7_104.rda$dvf, 
                           0,dta$c7_014.rda$dvf)
dta$c7_015.rda$dvf<-ifelse(dta$c7_015.rda$dvf==dta$c7_105.rda$dvf, 
                           0,dta$c7_015.rda$dvf)

#I have to fix dose to be positive. because there are some negative 
#values
#The min dose value in the original data (small size) is 0.0035
dta$c7_111.rda$dvf<-ifelse(dta$c7_111.rda$dvf<0,0,dta$c7_111.rda$dvf) 

#Make sure all images have the same size
# c7_111 id the dose file

#If there are more files in dat than the images we need, then delete
remove(c7_000.rda, c7_000ptv.rda, c7_001.rda, c7_002.rda, c7_003.rda,
       c7_004.rda, c7_005.rda, envir = dta)

# Cubic Patch ----------------------------------------------------------
in_avg_s<-function(x, s){
  #Calculates de average of the values in the cube
  d1<-dim(x$dvf)[1]
  d3<-dim(x$dvf)[3]
  l<-seq(from = 1, to = d1, by = s)
  w<-seq(from = 1, to = d3, by = s)
  l2<-seq(from = 1, to = length(l), by = 1)
  w2<-seq(from = 1, to = length(w), by = 1)
  res<-data.frame(matrix(nrow = length(l)*length(l)*length(w), 
                         ncol = 13))
  colnames(res)<-c("patch_w", "patch_l", "patch_d", "patch_w2", 
                   "patch_l2", "patch_d2", "patch", "avg", "min", "max", 
                   "sd", "med", "riq")
  res[,1]<-rep(l,times=length(l)*length(w))
  res[,2]<-rep(sapply(l, function(x) rep(x,length(l))), length(w))
  res[,3]<-unlist(lapply(w, function(x) rep(x,length(l)^2)))
  res[,4]<-rep(l2,times=length(l)*length(w))
  res[,5]<-rep(sapply(l2, function(x) rep(x,length(l))), length(w))
  res[,6]<-unlist(lapply(w2, function(x) rep(x,length(l)^2)))
  res[,7]<-paste(res[,4],res[,5], res[,6], sep="_" )
  for(i in 1:nrow(res)){
    a<-as.numeric(res[i,1])
    b<-as.numeric(min(a+s-1,d1))
    c<-as.numeric(res[i,2])
    d<-as.numeric(min(c+s-1,d1))
    e<-as.numeric(res[i,3])
    f<-as.numeric(min(e+s-1,d3))
    res[i,8]<-mean(x$dvf[a:b,c:d,e:f], na.rm=TRUE)
    res[i,9]<-min(x$dvf[a:b,c:d,e:f], na.rm = TRUE)
    res[i,10]<-max(x$dvf[a:b,c:d,e:f], na.rm = TRUE)
    res[i,11]<-sd(x$dvf[a:b,c:d,e:f], na.rm=TRUE)
    res[i,12]<-median(x$dvf[a:b,c:d,e:f], na.rm=TRUE)
    res[i,13]<-quantile(x$dvf[a:b,c:d,e:f], 0.75, na.rm=TRUE)-
      quantile(x$dvf[a:b,c:d,e:f], 0.25, na.rm=TRUE)
  }
  return(res)
}

numb<-function(x,s){
  d1<-dim(x$dvf)[1]
  d3<-dim(x$dvf)[3]
  
  matl<-array(NA,c(d1,d1,d3))
  matw<-array(NA,c(d1,d1,d3))
  matd<-array(NA,c(d1,d1,d3))
  for(l in 1:d1){
    matl[l, ,]<-ceiling(l/s)
  }
  for(w in 1:d1){
    matw[,w ,]<-ceiling(w/s)
  }
  for(d in 1:d3){
    matd[,,d]<-ceiling(d/s)
  }
  mat<-data.frame(matrix(NA, nrow=d1*d1*d3, ncol=4))
  colnames(mat)<-c("patch_l", "patch_w",  "patch_d", "patch")
  mat[,1]<-as.vector(matl)
  mat[,2]<-as.vector(matw)
  mat[,3]<-as.vector(matd)
  mat[,4]<-paste(mat[,1], mat[,2], mat[,3], sep = "_", collapse = NULL)
  return(mat)
}

in_avg_s2<-function(x, s){
  #this fuction is summarizing dose in the patch
  #(max, min ...)
  d1<-dim(x$dvf)[1]
  d3<-dim(x$dvf)[3]
  
  subx<-data.frame(matrix(NA, nrow=d1*d1*d3, ncol=5))
  colnames(subx)<-c("val", "patch_l", "patch_w",  "patch_d", "patch")
  subx[,1]<-as.vector(x$dvf[1:d1,1:d1,1:d3])
  subx[,2:5]<-mat[,1:4]
  n_t<-aggregate(list(n=subx$val), by=list(patch=subx$patch), 
                 function(x) length(which(!is.na(x))))
  mean_t<-aggregate(list(avg=subx$val), by=list(patch=subx$patch), 
                    FUN=mean, na.rm=TRUE)
  min_t<-aggregate(list(min=subx$val), by=list(patch=subx$patch), 
                   FUN=min, na.rm=TRUE)
  max_t<-aggregate(list(max=subx$val), by=list(patch=subx$patch), 
                   FUN=max, na.rm=TRUE)
  sd_t<-aggregate(list(sd=subx$val), by=list(patch=subx$patch), 
                  FUN=sd, na.rm=TRUE)
  median_t<-aggregate(list(med=subx$val), by=list(patch=subx$patch), 
                      FUN=median, na.rm=TRUE)
  iqr_t<-aggregate(list(riq=subx$val), by=list(patch=subx$patch), 
                   function(x) quantile(x,probs=0.75, na.rm=TRUE)-
                     quantile(x,probs=0.25, na.rm=TRUE))
  
  res<-Reduce(function(x, y) merge(x, y, by="patch"), 
              list(n_t, mean_t, min_t, max_t, sd_t, median_t, iqr_t))
  return(res)
}

P_num<-function(x, s){
  #Build the cubic patch
  d1<-dim(x$dvf)[1]
  d3<-dim(x$dvf)[3]
  dosel_111<-data.frame(matrix(NA, nrow=d1*d1*d3, ncol=8))
  colnames(dosel_111)<-c("dose", "patch_l", "patch_w",  "patch_d",
                         "patch","row", "col", "layer")
  dosel_111[,1]<-as.vector(x$dvf[1:d1,1:d1,1:d3])
  
  dosel_111[,2:5]<-mat[,1:4]
  dosel_111[,6]<-rep(1:d1, by=1, time=d1*d3)
  dosel_111[,7]<-rep(1:d1, by=1, each=d1, time=d3)
  dosel_111[,8]<-rep(1:d3, by=1, each=d1*d1)
  return(dosel_111)
}

k<-8 #This is patch size, if want different size then change 
mat<-numb(dta$c7_111.rda,s=k)
dosel<-P_num(dta$c7_111.rda,s=k)
Resumen<-eapply(dta, function(x)in_avg_s2(x,s=k), USE.NAMES = TRUE)

v <- in_avg_s2(dta$c7_011.rda,s=k)

Resumen_l<-ldply(Resumen)
Resumen_l$Id<-sub("_.*", "", Resumen_l$.id)
Resumen_l$file<-substr(sub(".*_", "", Resumen_l$.id) ,1,2) 
Resumen_l$date<-substr(sub(".*_", "", Resumen_l$.id) ,3,3) 
Resumen_l<-Resumen_l[,c(".id", "Id", "file", "date", "patch", "n", 
                        "avg", "min", "max", "sd", "med", "riq")]

Resumen_w<-reshape(Resumen_l[,-1], idvar = c("Id", "date", "patch"),
                   timevar = "file", direction =  "wide")
Resumen_w<-Resumen_w[,c("Id", "date", "patch", "n.01", "avg.01", 
                        "avg.10","avg.11", "min.11", "max.11", "sd.11",
                        "med.11", "riq.11")]
colnames(Resumen_w)<-c("Id", "date", "patch", "n" , "Hazy", "Dense", 
                       "Dose_avg", "Dose_min", "Dose_max", "Dose_sd", 
                       "Dose_med", "Dose_riq")

Resumen_w$Com<-(1-(Resumen_w$Hazy+Resumen_w$Dense))

Resumen_w2<-reshape(Resumen_w, idvar = c("Id","patch"),  timevar =
                      "date",direction =  "wide")
Resumen_w2<-Resumen_w2[,c("Id", "patch", "n.1", "Dose_avg.1",
                          "Dose_min.1","Dose_max.1", "Dose_sd.1", 
                          "Dose_med.1", "Dose_riq.1","Dense.1", "Hazy.1",
                          "Com.1", "Dense.2", "Hazy.2", "Com.2",
                          "Dense.3", "Hazy.3","Com.3", "Dense.4",
                          "Hazy.4", "Com.4", "Dense.5", "Hazy.5", 
                          "Com.5")]
Resumen_w$date2<-ifelse(Resumen_w$date==1,3,ifelse(Resumen_w$date==2,6,
                                                   ifelse(Resumen_w$date==3,12,
                                                          ifelse(Resumen_w$date==4,18,24))))
#ZeroChange are patches that were always 100% normal
list_ZeroChange <- Resumen_w2[which(Resumen_w2$Com.1==1 &
                                      Resumen_w2$Com.2==1 &
                                      Resumen_w2$Com.3==1 & 
                                      Resumen_w2$Com.4==1 &  
                                      Resumen_w2$Com.5==1), 2]
length(list_ZeroChange)
list_maxDoseZero<-Resumen_w[which(Resumen_w$Dose_max<=0),3]

Resumen_w$chageID<-ifelse(Resumen_w$patch %in% list_ZeroChange,0, 1)
Resumen_w$someDose<-ifelse(Resumen_w$patch %in% list_maxDoseZero,0, 1)
dosel$chageID<-ifelse(dosel$patch %in% list_ZeroChange,0,1)
dosel$someDose<-ifelse(dosel$patch %in% list_maxDoseZero,0, 1)
dosel<-dosel[!is.na(dosel$dose),]      

rm(dta,files, mat, matl, matd, matw)
rm(end, inicio, in_avg_s, in_avg_s2, numb, P_num)

Resumen_w$keep<-ifelse(Resumen_w$chageID==1,1,0)
Resumen_w$keep<-ifelse(Resumen_w$someDose==1,1,Resumen_w$keep)
dataset<-Resumen_w[which(Resumen_w$keep==1),]
dataset$ZeroHD<-ifelse(dataset$Com==1,1,0) #Both hazy and dense are zero
datadose<-dataset[dataset$date2==3,c("patch","Dose_avg","Dose_min",
                                     "Dose_max","Dose_sd","Dose_med","Dose_riq")]
dataset<-dataset[,c("Id","patch","date","date2","Hazy","Dense","Com",
                    "chageID", "someDose","keep","ZeroHD")]
dataset<-merge(dataset,datadose,by=c("patch"), all.x=T)
dataset_w<-Resumen_w2[Resumen_w2$patch %in% dataset$patch,]

save.image(file = paste("cubic_patch", k,".RData", sep = "")) 

# Spherical patch ------------------------------------------------------
max.val<-which(dta$c7_111.rda$dvf == max(dta$c7_111.rda$dvf, na.rm=T),
               arr.ind = TRUE)
numb.centered<-function(x,s,v){
  d1<-dim(x$dvf)[1]
  d3<-dim(x$dvf)[3]
  
  matl<-array(NA,c(d1,d1,d3))
  matw<-array(NA,c(d1,d1,d3))
  matd<-array(NA,c(d1,d1,d3))
  for(l in 1:d1){
    matl[l, ,]<-ceiling(l/s)
  }
  for(w in 1:d1){
    matw[,w ,]<-ceiling(w/s)
  }
  for(d in 1:d3){
    matd[,,d]<-ceiling(d/s)
  }
  mat<-data.frame(matrix(NA, nrow=d1*d1*d3, ncol=7))
  colnames(mat)<-c("patch_l","patch_w","patch_d","dif_l","dif_w",
                   "dif_d","patch")
  mat[,1]<-as.vector(matl)
  mat[,2]<-as.vector(matw)
  mat[,3]<-as.vector(matd)
  mat[,4]<-mat[,1]-v[1]
  mat[,5]<-mat[,2]-v[2]
  mat[,6]<-mat[,3]-v[3]
  mat[,7]<-floor(sqrt(mat[,4]^2 + mat[,5]^2 + mat[,6]^2))
  mat[,7]<-ifelse(mat[,7]==0,1,mat[,7])
  return(mat)
}

D_sum<-function(x, s,name="value"){
  d1<-dim(x$dvf)[1]
  d3<-dim(x$dvf)[3]
  unwrap<-data.frame(matrix(NA, nrow=d1*d1*d3, ncol=8))
  colnames(unwrap)<-c(name, "patch_l", "patch_w",  "patch_d",
                      "patch", "row", "col", "layer")
  unwrap[,1]<-as.vector(x$dvf[1:d1,1:d1,1:d3])
  
  unwrap[,2:5]<-mat[,c(1:3,7)]
  unwrap[,6]<-rep(1:d1, by=1, time=d1*d3)
  unwrap[,7]<-rep(1:d1, by=1, each=d1, time=d3)
  unwrap[,8]<-rep(1:d3, by=1, each=d1*d1)
  
  Resumen<-aggregate(list(n=unwrap[name]), by=list(patch = unwrap$patch),
                     FUN=length)
  colnames(Resumen)<-c("patch","n")
  
  Resumen<-merge(Resumen,aggregate(list(Dose_min=unwrap[name]),
                                   by=list(patch = unwrap$patch), 
                                   FUN=min),by="patch")
  colnames(Resumen)<-c("patch","n","Dose_min")
  Resumen<-merge(Resumen,aggregate(list(Dose_max=unwrap[name]),
                                   by=list(patch = unwrap$patch), 
                                   FUN=max),by="patch")
  colnames(Resumen)<-c("patch","n","Dose_min","Dose_max")
  Resumen<-merge(Resumen,aggregate(list(Dose_med=unwrap[name]), 
                                   by=list(patch = unwrap$patch), 
                                   FUN=median),by="patch")
  colnames(Resumen)<-c("patch","n","Dose_min","Dose_max","Dose_med")
  Resumen<-merge(Resumen,aggregate(list(Dose_avg=unwrap[name]), 
                                   by=list(patch = unwrap$patch), 
                                   FUN=mean),by="patch")
  colnames(Resumen)<-c("patch","n","Dose_min","Dose_max","Dose_med",
                       "Dose_avg")
  Resumen<-merge(Resumen,aggregate(list(Dose_sd=unwrap[name]), 
                                   by=list(patch = unwrap$patch), 
                                   FUN=sd),by="patch")
  colnames(Resumen)<-c("patch","n","Dose_min","Dose_max","Dose_med",
                       "Dose_avg","Dose_sd")
  
  dosel<-list(Resumen=Resumen, dose_raw=unwrap)
  return(dosel)
}

DH_sum<-function(x){
  d1<-dim(x$dvf)[1]
  d3<-dim(x$dvf)[3]
  unwrap<-data.frame(matrix(NA, nrow=d1*d1*d3, ncol=5))
  colnames(unwrap)<-c("value", "patch_l", "patch_w",  "patch_d", "patch")
  unwrap[,1]<-as.vector(x$dvf[1:d1,1:d1,1:d3])
  
  unwrap[,2:5]<-mat[,c(1:3,7)]
  
  Resumen<-aggregate(list(count=unwrap$value), by=list(patch = 
                                                         unwrap$patch), FUN=sum)
  return(Resumen)
}

max.d<-max.val[1,]
mat<-numb.centered(dta$c7_111.rda,s=1,v=max.d)
dosel.f<-D_sum(dta$c7_111.rda,s=1,name="dose")
dosel<-dosel.f$Resumen

Resumen<-eapply(dta, function(x) DH_sum(x), USE.NAMES = TRUE) 

Resumen_l<-ldply(Resumen)
Resumen_l<-separate(Resumen_l, .id, c("Id","file"), "_")
Resumen_l$date<-as.numeric(as.character(substr(Resumen_l$file, 3, 3)))
Resumen_l$file<-as.character(substr(Resumen_l$file, 1, 2))
Resumen_l<-Resumen_l[,c(1,2,5,3,4)]
Resumen_l<-Resumen_l[Resumen_l$date!=0,]
Resumen_l<-Resumen_l[order(Resumen_l$date),]

Resumen_w<-reshape(Resumen_l, idvar = c("Id","date", "patch"),  
                   timevar = "file", direction =  "wide")
Resumen_w<-Resumen_w[,c("Id","patch","date","count.01","count.10")]
names(Resumen_w)<- c("Id", "patch", "date",  "Hazy.n", "Dense.n")
Resumen_w<-merge(Resumen_w,dosel,by="patch",all.x =T)
Resumen_w$Hazy<-Resumen_w$Hazy.n/Resumen_w$n
Resumen_w$Dense<-Resumen_w$Dense.n/Resumen_w$n
Resumen_w$Com<-(1-(Resumen_w$Hazy+Resumen_w$Dense))

Resumen_w<-Resumen_w[order(Resumen_w$patch, Resumen_w$date),]
Resumen_w<-Resumen_w[Resumen_w$date!=0,]

Resumen_w2<-reshape(Resumen_w, idvar = c("Id","patch"),  
                    timevar = "date", direction =  "wide")

list_ZeroChange <- Resumen_w2[ which(Resumen_w2$Com.3==1 & 
                                       Resumen_w2$Com.6==1 &  
                                       Resumen_w2$Com.12==1 & 
                                       Resumen_w2$Com.18==1 & 
                                       Resumen_w2$Com.24==1), 1]
length(list_ZeroChange)

list_maxDoseZero<-Resumen_w[which(Resumen_w$Dose_max<=0),1]

Resumen_w$chageID<-ifelse(Resumen_w$patch %in% list_ZeroChange,0, 1)
Resumen_w$someDose<-ifelse(Resumen_w$patch %in% list_maxDoseZero,0, 1)
dosel<-dosel.f$dose_raw
save.image(file = paste("spheric_patch.RData", sep = "")) 

# Isodose patch --------------------------------------------------------
m<-ceiling(max(dta$c7_111.rda$dvf,na.rm = T))

#Categorize dose for isodose lines
s1<-c(-Inf,seq(0,m,by=1)) #0, 
s2<-c(-Inf,seq(0,m,by=0.5)) #0, ..
s3<-c(-Inf,seq(0,m,by=0.2)) #0, .
s4<-c(-Inf,seq(0,m,by=0.1)) #0, .
s5<-c(-Inf,seq(0,m,by=0.05)) #0, .

d1<-dim(dta$c7_111.rda$dvf)[1]
d3<-dim(dta$c7_111.rda$dvf)[3]
dosel<-data.frame(matrix(NA, nrow=d1*d1*d3, ncol=4))
colnames(dosel)<-c("dose", "row", "col", "layer")
dosel[,1]<-as.vector(dta$c7_111.rda$dvf[1:d1,1:d1,1:d3])
dosel[,2]<-rep(1:d1, by=1, time=d1*d3)
dosel[,3]<-rep(1:d1, by=1, each=d1, time=d3)
dosel[,4]<-rep(1:d3, by=1, each=d1*d1)
dosel$dosecat_s1<-cut(dosel$dose, breaks=s1)
dosel$dosecat_s2<-cut(dosel$dose, breaks=s2)
dosel$dosecat_s3<-cut(dosel$dose, breaks=s3)
dosel$dosecat_s4<-cut(dosel$dose, breaks=s4)
dosel$dosecat_s5<-cut(dosel$dose, breaks=s5)

length(unique(dosel$dosecat_s1))
length(unique(dosel$dosecat_s2))
length(unique(dosel$dosecat_s3))
length(unique(dosel$dosecat_s4))
length(unique(dosel$dosecat_s5))

max.val<-which(dta$c7_111.rda$dvf == max(dta$c7_111.rda$dvf, na.rm=T),
               arr.ind = TRUE)

for (i in 1:5) {
  nam <- paste("Resumen_s", i, sep = "")
  v<-aggregate(list(n=dosel$dose), by=list(patch = dosel[,4+i]),
               FUN=length)
  v<-merge(v,aggregate(list(Dose_min=dosel$dose), by=list(patch = 
                                                            dosel[,4+i]), FUN=min),by="patch")
  v<-merge(v,aggregate(list(Dose_max=dosel$dose), by=list(patch = 
                                                            dosel[,4+i]), FUN=max),by="patch")
  v<-merge(v,aggregate(list(Dose_med=dosel$dose), by=list(patch = 
                                                            dosel[,4+i]), FUN=median),by="patch")
  v<-merge(v,aggregate(list(Dose_avg=dosel$dose), by=list(patch = 
                                                            dosel[,4+i]), FUN=mean),by="patch")
  v<-merge(v,aggregate(list(Dose_sd=dosel$dose), by=list(patch =
                                                           dosel[,4+i]), FUN=sd),by="patch")
  assign(nam, v)
}
rm(v,nam)
rm(s1,s2,s3,s4,s5)

DH_sum<-function(x,dosel=dosel,cat=1){
  d1<-dim(x$dvf)[1]
  d3<-dim(x$dvf)[3]
  unwrap<-data.frame(matrix(NA, nrow=d1*d1*d3, ncol=5))
  colnames(unwrap)<-c("value", "patch_l", "patch_w",  "patch_d", "patch")
  unwrap[,1]<-as.vector(x$dvf[1:d1,1:d1,1:d3])
  unwrap[,2:5]<-dosel[,c(2:4,4+cat)]
  Resumen<-aggregate(list(count=unwrap$value), by=list(patch = 
                                                         unwrap$patch), FUN=sum)
  
  return(Resumen)
}

Resumen_1<-eapply(dta, function(x) DH_sum(x,dosel=dosel,cat=1), 
                  USE.NAMES = TRUE) 
Resumen_2<-eapply(dta, function(x) DH_sum(x,dosel=dosel,cat=2), 
                  USE.NAMES = TRUE) 
Resumen_3<-eapply(dta, function(x) DH_sum(x,dosel=dosel,cat=3), 
                  USE.NAMES = TRUE) 
Resumen_4<-eapply(dta, function(x) DH_sum(x,dosel=dosel,cat=4), 
                  USE.NAMES = TRUE) 
Resumen_5<-eapply(dta, function(x) DH_sum(x,dosel=dosel,cat=5), 
                  USE.NAMES = TRUE) 
# Code for isodose patch of 1Gy interval. 
# Repeat this code for other sized changing Resumen_1
Resumen_l<-ldply(Resumen_1)

Resumen_l<-separate(Resumen_l, .id, c("Id","file"), "_")
Resumen_l$date<-as.numeric(as.character(substr(Resumen_l$file, 3, 3)))
Resumen_l$file<-as.character(substr(Resumen_l$file, 1, 2))
Resumen_l<-Resumen_l[,c(1,2,5,3,4)]
Resumen_l<-Resumen_l[Resumen_l$date!=0,]
Resumen_l<-Resumen_l[order(Resumen_l$date),]

Resumen_w<-reshape(Resumen_l, idvar = c("Id","date", "patch"),  
                   timevar = "file", direction =  "wide")
Resumen_w<-Resumen_w[,c("Id","patch","date","count.01","count.10")]
names(Resumen_w)<- c("Id", "patch", "date",  "Hazy.n", "Dense.n")
Resumen_w<-Resumen_w[order(Resumen_w$date,Resumen_w$patch),]
Resumen_w<-merge(Resumen_w,Resumen_s1,by="patch",all.x =T)
Resumen_w<-Resumen_w[order(Resumen_w$date,Resumen_w$patch),]

Resumen_w$Hazy<-Resumen_w$Hazy.n/Resumen_w$n
Resumen_w$Dense<-Resumen_w$Dense.n/Resumen_w$n
Resumen_w$Com<-(1-(Resumen_w$Hazy+Resumen_w$Dense))

Resumen_w<-Resumen_w[order(Resumen_w$patch, Resumen_w$date),]
Resumen_w<-Resumen_w[Resumen_w$date!=0,]

Resumen_w2<-reshape(Resumen_w, idvar = c("Id","patch"),  
                    timevar = "date", direction =  "wide")

list_ZeroChange <- Resumen_w2[ which(Resumen_w2$Com.3==1 & 
                                       Resumen_w2$Com.6==1 & Resumen_w2$Com.12==1 & 
                                       Resumen_w2$Com.18==1 & Resumen_w2$Com.24==1), 1]
length(list_ZeroChange)
list_maxDoseZero<-Resumen_w[which(Resumen_w$Dose_max<=0),1]

Resumen_w$chageID<-ifelse(Resumen_w$patch %in% list_ZeroChange,0, 1)
Resumen_w$someDose<-ifelse(Resumen_w$patch %in% list_maxDoseZero,0, 1)
# After creating Resumen_1_w to Resumen_5_w then save data
save.image(file = paste("isodose_patch", k,".RData", sep = "")) 

# Compositional data ---------------------------------------------------
# After creating the data sets, some adjustments to the composition 
# were done (e.g. Fry zero transformation)

Resumen_w$keep<-ifelse(Resumen_w$chageID==1,1,0)
Resumen_w$keep<-ifelse(Resumen_w$someDose==1,1,Resumen_w$keep)
Resumen_w$ZeroHD<-ifelse(Resumen_w$Com==1,1,0) #Both hazy&dense are zero
dataset<-Resumen_w[which(Resumen_w$keep==1),]
dataset_w<-Resumen_w2[Resumen_w2$patch %in% dataset$patch,]

#Transform data
d<-0.003
N<-3
dataset$M<-apply(dataset[,c("Hazy","Dense","Com")],1,function(x)
  sum(x==0))
dataset$M<-ifelse(dataset$ZeroHD==1, NA,dataset$M)
dataset$ta<-ifelse(dataset$ZeroHD==1,
                   NA,d*(dataset$M+1)*(N-dataset$M)/(N*N))
dataset$ts<-dataset$ta*dataset$M
dataset$Dense.03<-ifelse(dataset$Dense==0,dataset$ta,ifelse(dataset$M==0,
                                                            dataset$Dense,dataset$Dense-(dataset$Dense*dataset$ts)))
dataset$Hazy.03<-ifelse(dataset$Hazy==0,dataset$ta,ifelse(dataset$M==0,
                                                          dataset$Hazy,dataset$Hazy-(dataset$Hazy*dataset$ts)))
dataset$Com.03<-ifelse(dataset$Com==0,dataset$ta,ifelse(dataset$M==0,
                                                        dataset$Com,dataset$Com-(dataset$Com*dataset$ts)))
dataset$Hazy.t<-dataset$Hazy.03/dataset$Com.03
dataset$Dense.t<-dataset$Dense.03/dataset$Com.03
dataset$Hazy.lt<-log(dataset$Hazy.t)
dataset$Dense.lt<-log(dataset$Dense.t)

dataset$Hazy.t0<-ifelse(is.na(dataset$Hazy.t)==T,0,dataset$Hazy.t)
dataset$Dense.t0<-ifelse(is.na(dataset$Dense.t)==T,0,dataset$Dense.t)
dataset$Hazy.Inon0<-ifelse(dataset$Hazy.t0==0,0,1)
dataset$Dense.Inon0<-ifelse(dataset$Dense.t0==0,0,1)
dataset$Hazy.I0<-1-dataset$Hazy.Inon0
dataset$Dense.I0<-1-dataset$Dense.Inon0

dataset$Hazy.lt0<-ifelse(is.na(dataset$Hazy.lt)==T,0,dataset$Hazy.lt)
dataset$Dense.lt0<-ifelse(is.na(dataset$Dense.lt)==T,0,dataset$Dense.lt)

dataset$date2<-ifelse(dataset$date==1,3,
                      ifelse(dataset$date==2,6,
                             ifelse(dataset$date==3,12,
                                    ifelse(dataset$date==4,18,24))))
dataset$years<-dataset$date2/12
dataset$years2<-as.numeric(as.character((dataset$years)))**2
dataset$years3<-as.numeric(as.character((dataset$years)))**3

dataset<-dataset[,c("Id","patch","n","date","date2","years",
                    "years2","years3","Dose_min",
                    "Dose_max","Dose_med","Dose_avg","Dose_sd",
                    "Hazy","Dense","Com","chageID","someDose","keep",
                    "ZeroHD","M","ta","ts","Dense.03","Hazy.03",
                    "Com.03","Hazy.t","Dense.t","Hazy.lt","Dense.lt",
                    "Hazy.t0","Dense.t0","Hazy.Inon0","Dense.Inon0",
                    "Hazy.I0","Dense.I0","Hazy.lt0","Dense.lt0")]
rm(list=setdiff(ls(), c("dataset", "dosel","list_maxDoseZero",
                        "list_ZeroChange")))

# Calinski-Harabasz (CH) index
library(nlme)
library(reghelper)
ctrl <- lmeControl(maxIter = 100, msMaxIter = 100, singular.ok=TRUE,
                   returnObject=TRUE, opt='optim')
#For cubic
mod <- lme( dose~ 1,  data=dosel, random = ~ 1|patch, control=ctrl)
#For spherical
mod <- lme( dose~ 1,  data=dosel, random = ~ 1|patch, control=ctrl)
#For isodose
mod <- lme( dose~ 1,  data=dosel, random = ~ 1|dosecat_s1, control=ctrl)

icc<-ICC(mod)
sse<-mod$sigma^2
ssm<-diag(sqrt(getVarCov(mod)))^2 
icc.h<-ssm/(ssm+sse)
ch<-(ssm/(mod$dims$ngrps[1]-1))/(sse/(mod$dims$N - mod$dims$ngrps[1]))

# Figures --------------------------------------------------------------
library(ggplot2)
library(reshape2)
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., 
                                                 draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - 
                                                 violinwidth * (x - xmin), xmaxv = x + 
                                                 violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x =
                                                                  if (grp %% 2 == 1) xminv else xmaxv), if 
                                                      (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, 
                                              newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, 
                                       nrow(newdata)), "x"] <- round(newdata[1, 
                                                                             "x"])
                             
                             if (length(draw_quantiles) > 0 & 
                                 !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), 
                                         all(draw_quantiles <=1))
                               quantiles <-
                                 gplot2:::create_quantile_segment_frame(data,
                                                                        draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), 
                                                  setdiff(names(data), c("x", "y")), drop = 
                                                    FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, 
                                                                    ...)
                               ggplot2:::ggname("geom_split_violin", 
                                                grid::grobTree(GeomPolygon$draw_panel(
                                                  newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", 
                                                GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = 
                                "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale =
                                "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = 
          GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = 
          inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = 
                        draw_quantiles, na.rm = na.rm, ...))
}

# Figure 2.12 (a)  [2.17, 2.23]
c<-round(prop.table(table(dataset$ZeroHD,dataset$date2),2)*100,2)
c<- as.data.frame(cbind(Months=unique(dataset$date2),P=c[1,]) )
len<-nrow(dataset)/length(unique(dataset$date2))
ggplot(data=c, aes(x=factor(Months), y=P)) +
  geom_bar(stat="identity",fill=alpha(c("#926B8D"),0.7)) +
  xlab("Months after RT \n") + 
  ylab(paste("% of patches with RILD (N=",len,")",sep = "")) +
  scale_y_continuous(sec.axis = sec_axis(~ . * len/100, 
                                         name="Number of patches with RILD"),limits=c(0,1))

# Figure 2.12 (b) [2.17, 2.23]
boxplot<-as.data.frame(dataset[dataset$date==1,])
ggplot(data = boxplot, aes(x=factor(chageID), y=Dose_med)) + 
  geom_boxplot(aes(color=factor(chageID), fill = factor(chageID)), 
               alpha=0.7) +
  xlab("Patch with RILD Presence") + ylab("Median Dose x Patch") +
  geom_point(pch = 19, position = position_jitterdodge(), 
             alpha=0.3, size=1,aes(colour = factor(chageID)))+
  scale_fill_manual(values = c("#B7B6BA", "#926B8D")) +
  scale_color_manual(values = c("#B7B6BA", "#926B8D")) +
  scale_x_discrete(labels=c("0"= "No", "1"="Yes")) +
  scale_y_continuous(breaks=seq(0,max(boxplot$Dose_med),by=10))

# Figure 2.12 (c) [2.17, 2.23]
my_dat2<-dataset[,c("date2","Dense.lt","Hazy.lt")]
my_dat2 <- melt(my_dat2, id.vars=c("date2"))

ggplot(data = my_dat2, aes(x=as.factor(date2), y=value, 
                           fill =as.factor(variable),
                           color=as.factor(variable))) + 
  geom_split_violin(trim=FALSE) +
  xlab("\nMonths after RT") + 
  ylab("alr transformation \n ln(Y/N)\n") + 
  scale_color_manual(name = "Type of RILD", 
                     values=c("#926B8D","#7aadb1"),labels=c("Dense", 
                                                            "Hazy")) +
  scale_fill_manual(name = "Type of RILD",  
                    values=alpha(c("#926B8D","#7aadb1"),0.7),
                    labels=c("Dense", "Hazy"))

#Figure 2.14 [2.18, 2.24]
library(compositions)
library(ggtern)
library(png)
a<-c(3,6,12,18,21)
k<-1 #Change k for different time points
ggtern(data=dataset[dataset$date2==a[[k]],c("Hazy","Dense","Com")], 
       aes(x=Dense, y=Hazy, z=Com)) + geom_point() + 
  theme_showarrows() + ggtitle(paste(a[[k]], " Months", sep = "")) +
  theme(plot.title = element_text(face="bold", hjust = 0.5,vjust=8),
        plot.margin = unit(c(0,0,0,0), "cm")) 

#Figure 2.15 [2.19, 2.25]
ggplot(dataset, aes(x=Dose_med, y=Dense.lt0)) +
  geom_point(alpha=0.7,col="#926b8d") +
  xlab("Median dose") + ylab("alr(Dense)")+
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major =
          element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour =
                                                                       "black"),text = element_text(size=16),
        legend.position="bottom",legend.title = element_blank())

ggplot(dataset, aes(x=Dose_med, y=Hazy.lt0)) + 
  geom_point(alpha=0.7,col="#926b8d") +
  xlab("Median dose") + ylab("alr(Hazy)")+
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = 
          element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour=
                                                                      "black"),text = element_text(size=16),
        legend.position="bottom",legend.title = element_blank())