# Packages  ------------------------------------------------------------
library(xlsx)
library(ggplot2)
library(compositions)
library(e1071)
library(caret)
library(pROC)

# Data  ----------------------------------------------------------------
# 1. Run model with data set of each patient, save model's coefficients 
and transform them back to compositional data. 
# 2. Use the predictions for each patient and cluster them using the 
medoids of the reference patient. 

coef<-read.xlsx("File with coefficients for each patient's model")
coef_com<-read.xlsx("File with coefficients in compositional form for 
                     each patient's model")
cluster<-read.xlsx("File with distribution (in percentage) of patches in 
                     each clustering configuration for each patient")

# Figures  =============================================================
# Figure 7.1 
# replace beta2 for beta3-beta4, gamma2-gamma8
ggplot(data=coef, aes(x = beta2)) +
  geom_histogram(aes(color = factor(rp), fill = factor(rp)), 
                 position = "identity", bins = 30, alpha = 0.4) +
  xlab("Coefficient for dose") + ylab("Count") + 
  scale_color_manual(name="RP",values = c("#926b8d","#E8CA47"),
                     labels = c("No", "Yes")) +
  scale_fill_manual(name="RP",values = c("#926b8d","#E8CA47"),
                    labels = c("No", "Yes")) 

### Regression ###########################################################
cluster$rpf<-ifelse(cluster$rp==0,"No","Yes")

X = acomp( cluster[,13:17] )
X = zeroreplace(X,0.001)
mylogit <- glm(rp ~ ilr(X), data = cluster, family = "binomial")
(a = coef(mylogit)[1])
(b = ilrInv(coef(mylogit)[-1],orig=X))
summary(mylogit)
predict(mylogit)
predpr<-predict(mylogit, type = "response")
mylogit_pred = ifelse(predict(mylogit, type = "link") > 0, "Yes", "No")

calc_class_err = function(actual, predicted) {
  mean(actual != predicted)
}
calc_class_err(actual = cluster$rpf, predicted = mylogit_pred)

train_tab = table(predicted = mylogit_pred, actual = cluster$rpf)
train_con_mat = confusionMatrix(train_tab, positive = "Yes")
c(train_con_mat$overall["Accuracy"], 
  train_con_mat$byClass["Sensitivity"], 
  train_con_mat$byClass["Specificity"])

### ROC #################################################################
roccurve <- roc(cluster$rp ~ predpr)
roccurve 
plot(roccurve)