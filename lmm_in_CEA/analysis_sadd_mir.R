###################################################

######### SADD data - main analysis (CCA, MI, LMM) for complete and all cases

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

####QALYs/utilities

########################### complete case analysis (n=101; n1=54;n2=47)
library(emmeans)
data_wide.cc<-data_wide[complete.cases(data_wide),]
#lm on QALYs
lm.cc_qaly<-lm(qaly~trt_f+eq5d_1,data = data_wide.cc)
lm.cc_qaly_means<-emmeans(lm.cc_qaly,~trt_f)
lm.cc_qaly_diff<-contrast(lm.cc_qaly_means,list(diff=c(-1,1)))
lm.cc_qaly_diff_ci<-confint(lm.cc_qaly_diff)

#lmm on utilities (for comparison)
library(reshape)
library(nlme)
library(lme4)
data_long.eq5d.cc<-reshape(data_wide.cc, varying = c("eq5d_1","eq5d_2","eq5d_3"),
                        direction = "long",idvar = "id",sep = "_")

data_long.eq5d.cc$time_f2<-ifelse(data_long.eq5d.cc$time==2,1,0)
data_long.eq5d.cc$time_f3<-ifelse(data_long.eq5d.cc$time==3,1,0)
data_long.eq5d.cc$trt_time_f2<-data_long.eq5d.cc$trt*data_long.eq5d.cc$time_f2
data_long.eq5d.cc$trt_time_f3<-data_long.eq5d.cc$trt*data_long.eq5d.cc$time_f3
data_long.eq5d.cc$time_f<- factor(data_long.eq5d.cc$time)

lmm.cc_qaly<-lme(eq5d ~ -1 + time_f + trt_time_f2 + trt_time_f3, random = ~ 1 | id, 
                 data=data_long.eq5d.cc, method = "ML", 
                 correlation = corSymm(form =~ as.numeric(time_f)|id),
                 weights=varIdent(form=~1|as.numeric(time_f)),
                 na.action = na.omit)
summary(lmm.cc_qaly)
intervals(lmm.cc_qaly, level = 0.95, which = "fixed")

lmm.cc_u_means<-emmeans(lmm.cc_qaly,~ -1 + time_f + trt_time_f2 + trt_time_f3)
lmm.cc_qaly_means<-contrast(lmm.cc_u_means,list(mu0 = c(13/104,13/104 + 26/104,26/104,0,0,0,0,0,0,0,0,0), 
                                          mu1=c(13/104,0,0,0,13/104 + 26/104,0,0,0,26/104,0,0,0)))
lmm.cc_qaly_means_ci<-confint(lmm.cc_qaly_means)
lmm.cc_qaly_diff<-contrast(lmm.cc_u_means,list(diff = c(0,-13/104 - 26/104,-26/104,0,13/104 + 26/104,0,0,0,26/104,0,0,0)))
lmm.cc_qaly_diff_ci<-confint(lmm.cc_qaly_diff)


########################### all case analysis 

########################### multiple imputation under MAR
library(mice)
#MICE
data_wide.mi<-data_wide[,c("trt","eq5d0","eq5d13","eq5d39","qaly")]
mice.wide.qalys.ini<-mice(data_wide.mi,maxit = 0, method = "norm.nob", print = FALSE)
meth<- mice.wide.qalys.ini$meth
pred <-mice.wide.qalys.ini$pred
pred[c("trt"), c("qaly","eq5d0","eq5d13","eq5d39")] <- 0
pred[c("eq5d0"), c("qaly","eq5d13","eq5d39","trt")] <- 0
pred["eq5d13", c("qaly")] <- 0
#pred["eq5d13", c("qaly","eq5d39")] <- 0
pred["eq5d39", c("qaly")] <- 0
pred[c("qaly"), c("trt","eq5d0","eq5d13","eq5d39")] <- 0
meth["qaly"]<- "~ I(13/104*eq5d0 + 13/104*eq5d13 + 26/104*eq5d13 + 26/104*eq5d39)"
mice.wide.qalys<-mice(data_wide.mi, meth=meth, pred=pred, maxit=100,m=500, seed=345, print=F)
mice.fit.wide.qalys<-with(data = mice.wide.qalys, exp = lm(qaly ~ trt + eq5d0))
mice_qalys_means<-emmeans(mice.fit.wide.qalys, ~trt)
mice_qaly_diff<-contrast(mice_qalys_means,list(diff=c(-1,1)))
mice_qaly_diff_ci<-confint(mice_qaly_diff)

mice.fit.wide.u1<-with(data = mice.wide.qalys, exp = lm(eq5d0 ~ 1))
mice_u1_means<-emmeans(mice.fit.wide.u1, ~1)
mice.fit.wide.u2<-with(data = mice.wide.qalys, exp = lm(eq5d13 ~ trt + eq5d0))
mice_u2_means<-emmeans(mice.fit.wide.u2, ~trt)
mice_u2_diff<-contrast(mice_u2_means,list(diff=c(-1,1)))
mice_u2_diff_ci<-confint(mice_u2_diff)
mice.fit.wide.u3<-with(data = mice.wide.qalys, exp = lm(eq5d39 ~ trt + eq5d0))
mice_u3_means<-emmeans(mice.fit.wide.u3, ~trt)
mice_u3_diff<-contrast(mice_u3_means,list(diff=c(-1,1)))
mice_u3_diff_ci<-confint(mice_u3_diff)

########################### LMM
library(reshape)
library(nlme)
library(lme4)
data_long.eq5d<-reshape(data_wide, varying = c("eq5d_1","eq5d_2","eq5d_3"),
                        direction = "long",idvar = "id",sep = "_")
data_long.eq5d$time_f2<-ifelse(data_long.eq5d$time==2,1,0)
data_long.eq5d$time_f3<-ifelse(data_long.eq5d$time==3,1,0)
data_long.eq5d$trt_time_f2<-data_long.eq5d$trt*data_long.eq5d$time_f2
data_long.eq5d$trt_time_f3<-data_long.eq5d$trt*data_long.eq5d$time_f3
data_long.eq5d$time_f<- factor(data_long.eq5d$time)

lmm_qaly<-lme(eq5d ~ -1 + time_f + trt_time_f2 + trt_time_f3, random = ~ 1 | id, 
                 data=data_long.eq5d, method = "ML", 
                 correlation = corSymm(form =~ as.numeric(time_f)|id),
                 weights=varIdent(form=~1|as.numeric(time_f)),
                 na.action = na.omit)
summary(lmm_qaly)
intervals(lmm_qaly, level = 0.95)
lmm_u_means<-emmeans(lmm_qaly,~ -1 + time_f + trt_time_f2 + trt_time_f3)
lmm_qaly_means<-contrast(lmm_u_means,list(mu0 = c(13/104,13/104 + 26/104,26/104,0,0,0,0,0,0,0,0,0), 
                                          mu1=c(13/104,0,0,0,13/104 + 26/104,0,0,0,26/104,0,0,0)))
lmm_qaly_means_ci<-confint(lmm_qaly_means)
lmm_qaly_diff<-contrast(lmm_u_means,list(diff = c(0,-13/104 - 26/104,-26/104,0,13/104 + 26/104,0,0,0,26/104,0,0,0)))
lmm_qaly_diff_ci<-confint(lmm_qaly_diff)





####TCs/costs

########################### complete case analysis
data_wide.cc<-data_wide[complete.cases(data_wide),]
#lm on total costs
library(emmeans)
lm.cc_totc<-lm(tcost~trt_f+cost0,data = data_wide.cc)
lm.cc_totc_means<-emmeans(lm.cc_totc,~trt_f)
lm.cc_totc_diff<-contrast(lm.cc_totc_means,list(diff=c(-1,1)))
lm.cc_totc_diff_ci<-confint(lm.cc_totc_diff)

#lmm (for comparison)
library(reshape)
library(nlme)
library(lme4)
data_long.cost.cc<-reshape(data_wide.cc, varying = c("cost_1","cost_2","cost_3"),
                           direction = "long",idvar = "id",sep = "_")
data_long.cost.cc$time_f2<-ifelse(data_long.cost.cc$time==2,1,0)
data_long.cost.cc$time_f3<-ifelse(data_long.cost.cc$time==3,1,0)
data_long.cost.cc$trt_time_f2<-data_long.cost.cc$trt*data_long.cost.cc$time_f2
data_long.cost.cc$trt_time_f3<-data_long.cost.cc$trt*data_long.cost.cc$time_f3
data_long.cost.cc$time_f<- factor(data_long.cost.cc$time)

lmm.cc_tcost<-lme(cost ~ -1 + time_f + trt_time_f2 + trt_time_f3, random = ~ 1 | id, 
               data=data_long.cost.cc, method = "ML", 
               correlation = corSymm(form =~ as.numeric(time_f)|id),
               weights=varIdent(form=~1|as.numeric(time_f)),
               na.action = na.omit)
summary(lmm.cc_tcost)
intervals(lmm.cc_tcost, level = 0.95, which = "fixed")

lmm.cc_c_means<-emmeans(lmm.cc_tcost,~ -1 + time_f + trt_time_f2 + trt_time_f3)
lmm.cc_tcost_means<-contrast(lmm.cc_c_means,list(mu0=c(0,1,1,0,0,0,0,0,0,0,0,0), 
                                           mu1=c(0,0,0,0,1,0,0,0,1,0,0,0)))
lmm.cc_tcost_means_ci<-confint(lmm.cc_tcost_means)
lmm.cc_tcost_diff<-contrast(lmm.cc_c_means,list(diff = c(0,-1,-1,0,1,0,0,0,1,0,0,0)))
lmm.cc_tcost_diff_ci<-confint(lmm.cc_tcost_diff)

########################### all case analysis

########################### multiple imputation 
library(mice)
#MICE
data_wide.mi<-data_wide[,c("trt","cost0","cost13","cost39","tcost","id")]
mice.wide.totcs.ini<-mice(data_wide.mi,maxit = 0, method = "norm.nob", print = FALSE)
meth<- mice.wide.totcs.ini$meth
pred <-mice.wide.totcs.ini$pred
pred[,"id"] <- 0
pred["id",] <- 0
pred[c("trt"), c("tcost","cost0","cost13","cost39")] <- 0
pred[c("cost0"), c("tcost","cost13","cost39","trt")] <- 0
pred["cost13", c("tcost")] <- 0
#pred["cost13", c("tcost","cost39")] <- 0
pred["cost39", c("tcost")] <- 0
pred[c("tcost"), c("trt","cost0","cost13","cost39")] <- 0
meth["tcost"]<- "~ I(cost13 + cost39)"
mice.wide.totcs<-mice(data_wide.mi, meth=meth, pred=pred, maxit=100,m=500, seed=345, print=F)
mice.fit.wide.totcs<-with(data = mice.wide.totcs, exp = lm(tcost ~ trt + cost0))
mice_totcs_means<-emmeans(mice.fit.wide.totcs, ~trt)
mice_totc_diff<-contrast(mice_totcs_means,list(diff=c(-1,1)))
mice_totc_diff_ci<-confint(mice_totc_diff)

########################### LMM
data_wide$cost_3<-data_wide$cost39
data_wide$cost_2<-data_wide$cost13
data_wide$cost_1<-data_wide$cost0
data_long.cost<-reshape(data_wide, varying = c("cost_1","cost_2","cost_3"),
                        direction = "long",idvar = "id",sep = "_")

data_long.cost$time_f2<-ifelse(data_long.cost$time==2,1,0)
data_long.cost$time_f3<-ifelse(data_long.cost$time==3,1,0)
data_long.cost$trt_time_f2<-data_long.cost$trt*data_long.cost$time_f2
data_long.cost$trt_time_f3<-data_long.cost$trt*data_long.cost$time_f3
data_long.cost$time_f<- factor(data_long.cost$time)

lmm_tcost<-lme(cost ~ -1 + time_f + trt_time_f2 + trt_time_f3, random = ~ 1 | id, 
              data=data_long.cost, method = "ML", 
              correlation = corSymm(form =~ as.numeric(time_f)|id),
              weights=varIdent(form=~1|as.numeric(time_f)),
              na.action = na.omit)
summary(lmm_tcost)
intervals(lmm_tcost, level = 0.95, which = "fixed")

lmm_c_means<-emmeans(lmm_tcost,~ -1 + time_f + trt_time_f2 + trt_time_f3)
lmm_tcost_means<-contrast(lmm_c_means,list(mu0=c(0,1,1,0,0,0,0,0,0,0,0,0), 
                                          mu1=c(0,0,0,0,1,0,0,0,1,0,0,0)))
lmm_tcost_means_ci<-confint(lmm_tcost_means)
lmm_tcost_diff<-contrast(lmm_c_means,list(diff = c(0,-1,-1,0,1,0,0,0,1,0,0,0)))
lmm_tcost_diff_ci<-confint(lmm_tcost_diff)
