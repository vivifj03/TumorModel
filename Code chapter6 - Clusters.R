# Packages ----------------------------------------------------------------
library(compositions)
library(cluster)
library(Rfast)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(pspline)

# Data --------------------------------------------------------------------
setwd("]")
load(".RData") #Data set with prediction (v2 in code Chapter 5)
#v2 have the final predicted compositions

dcluster<-v2[,c("patch","Id","date","Hazy","Dense","Com","chageID","someDose" ,"nonZeroHD",
                "Dose_avg","Dose_min","Dose_max","Dose_sd","Dose_med",
                "dense.pred.f","hazy.pred.f","com.pred.f")]

patch2<-v2[!duplicated(v2$patch),c("patch","chageID") ]
rm(list=setdiff(ls(),c("v2","dcluster","patch2")))

dcluster[c("dense.p","hazy.p","com.p")]<-clo(dcluster[c("dense.pred.f","hazy.pred.f","com.pred.f")])
dcluster[,c("dense.p.clr","hazy.p.clr","com.p.clr")]<-clr(dcluster[c("dense.p","hazy.p","com.p")])

# Functions ---------------------------------------------------------


# Clustering -------------------------------------------------------
ni<-nrow(dcluster)/length(patch2)
df_list <- split(dcluster, as.factor(dcluster$date))
df_list$`3`<-df_list$`3`[,c("patch","dense.p","hazy.p","com.p",
                            "dense.p.clr","hazy.p.clr","com.p.clr")]
df_list$`6`<-df_list$`6`[,c("patch","dense.p","hazy.p","com.p",
                            "dense.p.clr","hazy.p.clr","com.p.clr")]
df_list$`12`<-df_list$`12`[,c("patch","dense.p","hazy.p","com.p",
                              "dense.p.clr","hazy.p.clr","com.p.clr")]
df_list$`18`<-df_list$`18`[,c("patch","dense.p","hazy.p","com.p",
                              "dense.p.clr","hazy.p.clr","com.p.clr")]
df_list$`24`<-df_list$`24`[,c("patch","dense.p","hazy.p","com.p",
                              "dense.p.clr","hazy.p.clr","com.p.clr")]

# Euclidian distance  ++++++++++++++=================================
rownames(df_list$`3`)<-df_list$`3`$patch
rownames(df_list$`6`)<-df_list$`6`$patch
rownames(df_list$`12`)<-df_list$`12`$patch
rownames(df_list$`18`)<-df_list$`18`$patch
rownames(df_list$`24`)<-df_list$`24`$patch
d3<-dist(df_list$`3`[,5:7], method = "euclidean")
d6<-dist(df_list$`6`[,5:7], method = "euclidean")
d12<-dist(df_list$`12`[,5:7], method = "euclidean")
d18<-dist(df_list$`18`[,5:7], method = "euclidean")
d24<-dist(df_list$`24`[,5:7], method = "euclidean")
dis<-(d3+d6+d12+d18+d24)/5

c2<-pam(dis,2,diss=T)
c3<-pam(dis,3,diss=T)
c4<-pam(dis,4,diss=T)
c5<-pam(dis,5,diss=T)
c6<-pam(dis,6,diss=T)
c7<-pam(dis,7,diss=T)

patch_c<-cbind(patch2,c2$clustering)
patch_c<-cbind(patch_c,c3$clustering)
patch_c<-cbind(patch_c,c4$clustering)
patch_c<-cbind(patch_c,c5$clustering)
patch_c<-cbind(patch_c,c6$clustering)
patch_c<-cbind(patch_c,c7$clustering)
patch_c<-merge(patch_c,v2[v2$date==3,c("patch","Dose_med")],
               by="patch")
colnames(patch_c)<-c("patch","changeID" ,"c2","c3","c4","c5","c6",
                     "c7","Dose")
patch_c<-as.data.frame(patch_c)

# Table 6.1
table(patch_c$changeID,patch_c$c2)
table(patch_c$changeID,patch_c$c3)
table(patch_c$changeID,patch_c$c4)
table(patch_c$changeID,patch_c$c5)
prop.table(table(patch_c$changeID,patch_c$c5),1)
table(patch_c$changeID,patch_c$c6)
prop.table(table(patch_c$changeID,patch_c$c6),1)
table(patch_c$changeID,patch_c$c7)
prop.table(table(patch_c$changeID,patch_c$c7),1)
dcluster<-merge(dcluster,patch_c,by="patch",all.x = T)

#Average silhouette coefficient
plot(c2,main=" ") 
plot(c3) 
plot(c4) 
plot(c5) 
plot(c6) 
plot(c7) 

# Figure 6.1
#Peter J. Rousseeuw (1987). "Silhouettes: a Graphical Aid to the 
#Interpretation and Validation of Cluster Analysis". Computational 
#and Applied Mathematics. 20: 53-65. doi:10.1016/0377-0427(87)90125-7.
x<-c(2,3,4,5,6,7)
y<-c(0.75,0.58,0.52,0.57,0.52,0.51)
dxy<-as.data.frame(cbind(x,y))
ggplot(dxy,aes(x = x,y = y)) +
  geom_line(col="#926b8d") + geom_point(col="#926b8d")+
  xlab("Number of clusters") + ylab("Average Silhouette Score")+
 
#Figure 6.2 
#replace c2 for c3-c7
ggplot(data = patch_c, aes(x=factor(c2), y=Dose)) + 
  geom_boxplot(alpha=0.7) +
  geom_point(pch = 19, position = position_jitterdodge(), 
             alpha=0.7, size=1,aes(colour = factor(changeID)))+
  xlab("Cluster") + ylab("Median Dose") +
  scale_color_manual(values = c("#B7B6BA", "#926B8D"),
                     labels = c("No RILD", "RILD"), name="") 

#Figure 6.3
# for Figure 6.4 replace dense.pred.f, hazy.pred.f, com.pred.f, for 
# Dense, Hazy, Com
cluster.labs <- c("Cluster 1","Cluster 2")
names(cluster.labs) <- c("1", "2")
ggplot(dcluster, aes(x=factor(date), y=dense.pred.f)) + 
  geom_boxplot(alpha=0.7, color="#926B8D")+
  facet_wrap(~c2,labeller = labeller(c2 = cluster.labs)) +
  xlab("Months after RT") + ylab("Predicted Dense %") + ylim(0,1) 

ggplot(dcluster, aes(x=factor(date), y=hazy.pred.f)) + 
  geom_boxplot(alpha=0.7, color="#7aadb1")+
  facet_wrap(~c2,labeller = labeller(c2 = cluster.labs)) +
  xlab("Months after RT") + ylab("Predicted Hazy %") + ylim(0,1) 
  
ggplot(dcluster, aes(x=factor(date), y=com.pred.f)) + 
  geom_boxplot(alpha=0.7, color="#B7B6BA")+
  facet_wrap(~c2,labeller = labeller(c2 = cluster.labs)) +
  xlab("Months after RT") + ylab("Predicted Normal %") + ylim(0,1) 

cluster.labs <- c("Cluster 1","Cluster 2","Cluster 3")
names(cluster.labs) <- c("1", "2", "3")
ggplot(dcluster, aes(x=factor(date), y=dense.pred.f)) + 
  geom_boxplot(alpha=0.7, color="#926B8D")+
  facet_wrap(~c3,labeller = labeller(c3 = cluster.labs)) +
  xlab("Months after RT") + ylab("Predicted Dense %") + ylim(0,1) 

ggplot(dcluster, aes(x=factor(date), y=hazy.pred.f)) + 
  geom_boxplot(alpha=0.7, color="#7aadb1")+
  facet_wrap(~c3,labeller = labeller(c3 = cluster.labs)) +
  xlab("Months after RT") + ylab("Predicted Hazy %") + ylim(0,1) 

ggplot(dcluster, aes(x=factor(date), y=com.pred.f)) + 
  geom_boxplot(alpha=0.7, color="#B7B6BA")+
  facet_wrap(~c3,labeller = labeller(c3 = cluster.labs)) +
  xlab("Months after RT") + ylab("Predicted Normal %") + ylim(0,1)

cluster.labs <- c("Cluster 1","Cluster 2","Cluster 3","Cluster 4")
names(cluster.labs) <- c("1", "2", "3","4")
ggplot(dcluster, aes(x=factor(date), y=dense.pred.f)) + 
  geom_boxplot(alpha=0.7, color="#926B8D")+
  facet_wrap(~c4,labeller = labeller(c4 = cluster.labs),nrow = 1) +
  xlab("Months after RT") + ylab("Predicted Dense %") + ylim(0,1) 

ggplot(dcluster, aes(x=factor(date), y=hazy.pred.f)) + 
  geom_boxplot(alpha=0.7, color="#7aadb1")+
  facet_wrap(~c4,labeller = labeller(c4 = cluster.labs),nrow = 1) +
  xlab("Months after RT") + ylab("Predicted Hazy %") + ylim(0,1) 

ggplot(dcluster, aes(x=factor(date), y=com.pred.f)) + 
  geom_boxplot(alpha=0.7, color="#B7B6BA")+
  facet_wrap(~c4,labeller = labeller(c4 = cluster.labs),nrow = 1)+
  xlab("Months after RT") + ylab("Predicted Normal %") + ylim(0,1)

cluster.labs <- c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5")
names(cluster.labs) <- c("1", "2", "3","4","5")
ggplot(dcluster, aes(x=factor(date), y=dense.pred.f)) + 
  geom_boxplot(alpha=0.7, color="#926B8D")+
  facet_wrap(~c5,labeller = labeller(c5 = cluster.labs),nrow = 1) +
  xlab("Months after RT") + ylab("Predicted Dense %") + ylim(0,1) 

ggplot(dcluster, aes(x=factor(date), y=hazy.pred.f)) + 
  geom_boxplot(alpha=0.7, color="#7aadb1")+
  facet_wrap(~c5,labeller = labeller(c5 = cluster.labs),nrow = 1) +
  xlab("Months after RT") + ylab("Predicted Hazy %") + ylim(0,1) 

ggplot(dcluster, aes(x=factor(date), y=com.pred.f)) + 
  geom_boxplot(alpha=0.7, color="#B7B6BA")+
  facet_wrap(~c5,labeller = labeller(c5 = cluster.labs),nrow = 1)+
  xlab("Months after RT") + ylab("Predicted Normal %") + ylim(0,1) 