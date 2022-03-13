###################################################

######### SADD data - bootstrap analysis (LMM) for all cases

####################################################

#read data
library(readstata13)
sadd_main<-read.dta13("SADD_main.dta")
sadd<-read.dta13("SADD_pre_MI.dta")
sadd<-sadd[,c("newid1","eq5d0","eq5d13","eq5d39","arm","cost0","cost13","cost39","qaly")]
data_wide<-sadd
data_wide$trt <- as.numeric(data_wide$arm)
data_wide$subjects<-rep(1:length(data_wide$arm))
data_wide$tcost<-data_wide$cost13 + data_wide$cost39

#remove arm sertraline - (comparison only for placebo and mirtazapine)
data_wide <- data_wide[data_wide$arm %in% c("Placebo","Mirtazapine"),]
data_wide.v1<-data_wide
data_wide$trt <- ifelse(data_wide$trt == 1, 0, 1)

#wide format
data_wide$id <- rep(1:nrow(data_wide))
data_wide$eq5d_1 <- data_wide$eq5d0
data_wide$eq5d_2 <- data_wide$eq5d13
data_wide$eq5d_3 <- data_wide$eq5d39
data_wide$cost_1 <- data_wide$cost0
data_wide$cost_2 <- data_wide$cost13
data_wide$cost_3 <- data_wide$cost39
data_wide$trt_f<-factor(data_wide$trt)
data_wide<-data_wide[,c("id","eq5d0","eq5d13","eq5d39","trt","trt_f","cost0","cost13","cost39",
                        "eq5d_1","eq5d_2","eq5d_3","cost_1","cost_2","cost_3")]
data_wide$qaly<-(13/104)*data_wide$eq5d0+(13/104 + 26/104)*data_wide$eq5d13+(26/104)*data_wide$eq5d39
data_wide$tcost<-data_wide$cost13+data_wide$cost39


############################################ bootstrapping + LMM
library(emmeans)
library(nlme)
library(lme4)
library(sommer)
library(boot)

#bootstrapping utilities

#define boot function
boot_lmm_u <- function(data, i){
  data2 <- data[i,]
  data_long_boot<-reshape(data2[i,], varying = c("eq5d_1","eq5d_2","eq5d_3"),
                          direction = "long",idvar = "id",sep = "_", new.row.names=sequence(657))
  data_long_boot$time_f2<-ifelse(data_long_boot$time==2,1,0)
  data_long_boot$time_f3<-ifelse(data_long_boot$time==3,1,0)
  data_long_boot$trt_time_f2<-data_long_boot$trt*data_long_boot$time_f2
  data_long_boot$trt_time_f3<-data_long_boot$trt*data_long_boot$time_f3
  data_long_boot$time_f<-as.factor(data_long_boot$time)
  cgm3_u_ml<-lme(eq5d ~ -1 + time_f + trt_time_f2 + trt_time_f3, random = ~ 1 | id, 
                    data=data_long_boot, method = "ML",
                    na.action = na.omit)
  em3_u_ml.eq5d<-emmeans(cgm3_u_ml,~ -1 + time_f + trt_time_f2 + trt_time_f3)
  em3_u_ml.eq5d.delta<-contrast(em3_u_ml.eq5d,list(mu0 = c(13/104,13/104 + 26/104,26/104,0,0,0,0,0,0,0,0,0), 
                                                  mu1=c(13/104,0,0,0,13/104 + 26/104,0,0,0,26/104,0,0,0)))
  em3_u_ml.eq5d.delta_ci<-confint(em3_u_ml.eq5d.delta)
  cgm3_u_ml.qalys.delta<-contrast(em3_u_ml.eq5d,list(diff = c(0,-13/104 - 26/104,-26/104,0,13/104 + 26/104,0,0,0,26/104,0,0,0)))
  cgm3_u_ml.qalys.delta_ci<-confint(cgm3_u_ml.qalys.delta)
  return(c(summary(em3_u_ml.eq5d)$emmean[c(1,2,3,5,9)], 
           summary(em3_u_ml.eq5d)$emmean[5]- summary(em3_u_ml.eq5d)$emmean[2],
           summary(em3_u_ml.eq5d)$emmean[9]- summary(em3_u_ml.eq5d)$emmean[3],
           summary(em3_u_ml.eq5d.delta)$estimate, 
           summary(cgm3_u_ml.qalys.delta)$estimate))
}

#apply boot function to data
set.seed(4567) 
boot_est_lmm_u <- boot(data_wide, boot_lmm_u, R=10000)
delta_e_boot<-boot_est_lmm_u$t[,10]
mu_e1_boot<-boot_est_lmm_u$t[,8]
mu_e2_boot<-boot_est_lmm_u$t[,9]
mu_e_boot<-cbind(mu_e1_boot,mu_e2_boot)
  

#bootstrapping costs

#define boot function
boot_lmm_c <- function(data, i){
  data2 <- data[i,]
  data_long_boot<-reshape(data2[i,], varying = c("cost_1","cost_2","cost_3"),
                          direction = "long",idvar = "id",sep = "_", new.row.names=sequence(657))
  data_long_boot$time_f2<-ifelse(data_long_boot$time==2,1,0)
  data_long_boot$time_f3<-ifelse(data_long_boot$time==3,1,0)
  data_long_boot$trt_time_f2<-data_long_boot$trt*data_long_boot$time_f2
  data_long_boot$trt_time_f3<-data_long_boot$trt*data_long_boot$time_f3
  data_long_boot$time_f<- factor(data_long_boot$time)
  cgm3_c_ml<-lme(cost ~ -1 + time_f + trt_time_f2 + trt_time_f3, random = ~ 1 | id, 
                    data=data_long_boot, method = "ML",
                    na.action = na.omit)
  em3_c_ml.cost<-emmeans(cgm3_c_ml,~ -1 + time_f + trt_time_f2 + trt_time_f3)
  em3_c_ml.cost.delta<-contrast(em3_c_ml.cost,list(mu0=c(0,1,1,0,0,0,0,0,0,0,0,0), 
                                                   mu1=c(0,0,0,0,1,0,0,0,1,0,0,0)))
  em3_c_ml.cost.delta_ci<-confint(em3_c_ml.cost.delta)
  cgm3_c_ml.cc.tcosts<-contrast(em3_c_ml.cost,list(diff = c(0,-1,-1,0,1,0,0,0,1,0,0,0)))
  cgm3_c_ml.cc.tcosts_ci<-confint(cgm3_c_ml.cc.tcosts)
  return(c(summary(em3_c_ml.cost)$emmean[c(1,2,3,5,9)], 
           summary(em3_c_ml.cost)$emmean[5]- summary(em3_c_ml.cost)$emmean[2],
           summary(em3_c_ml.cost)$emmean[9]- summary(em3_c_ml.cost)$emmean[3],
           summary(em3_c_ml.cost.delta)$estimate, 
           summary(cgm3_c_ml.cc.tcosts)$estimate))
}

#apply boot function to data
set.seed(4567)
boot_est_lmm_c <- boot(data_wide, boot_lmm_c, R=10000)
delta_tc_boot<-boot_est_lmm_c$t[,10]
mu_tc1_boot<-boot_est_lmm_c$t[,8]
mu_tc2_boot<-boot_est_lmm_c$t[,9]
mu_tc_boot<-cbind(mu_tc1_boot,mu_tc2_boot)


###########################################################
#obtain graphical outputs

#cep and ceac
library(BCEA)
cea_res<-bcea(e=mu_e_boot,c=mu_tc_boot, ref=2, Kmax = 40000)
icer<-cea_res$ICER
paste_icer<-paste("ICER",round(icer,0),sep = "=")

#CEP
opar <- par(ps=18) 
plot(cea_res$delta.e,cea_res$delta.c,xlim = c(-0.1,0.12),ylim = c(-10000,10000),main ="Cost-Effectiveness Plane",
     xlab = "Expected QALY differential",ylab = "Expected cost differential",type="n")
points(cea_res$delta.e,cea_res$delta.c,col="#7cb9b9",pch=16,cex=.7)
points(mean(cea_res$delta.e),mean(cea_res$delta.c),col="#007C7C",pch=16,cex=1)
abline(h=0,lty="dashed")
abline(v=0,lty="dashed")
abline(a=0,b=25000)
text(-0.04,-5000,"k=25000",cex = 0.7)
legend(0.065,-8000,legend = paste_icer,col = c("#007C7C"),pch = c(16),cex = 0.7,bty="n",
       y.intersp = 0.5,x.intersp = 0.5)

#CEAC
opar <- par(ps=18) 
plot(cea_res$k,cea_res$ceac,xlim = c(0,40000),ylim = c(0,1),main ="Cost-Effectiveness Acceptability Curve",
     xlab = "Acceptance threshold",ylab = "Probability of Cost-Effectiveness",type="n",axes=F)
axis(1)
axis(2)
lines(cea_res$k,cea_res$ceac,col="#007C7C",lty=1,lwd=3)
