###################################################

######### SADD data

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
#data_wide$sex <- factor(as.numeric(data_wide$sex2)-1)
#data_wide$age <- factor(data_wide$agecat3)
#data_wide$site <- factor(data_wide$site2)
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
#lm
lm.cc_qaly<-lm(qaly~trt_f+eq5d_1,data = data_wide.cc)
lm.cc_qaly_means<-emmeans(lm.cc_qaly,~trt_f)
lm.cc_qaly_diff<-contrast(lm.cc_qaly_means,list(diff=c(-1,1)))
lm.cc_qaly_diff_ci<-confint(lm.cc_qaly_diff)

#lmm
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

########################### multiple imputation (without costs as auxiliary variables)
library(mice)
#MICE
#mice-lm
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

#mice-lmm
data_wide_long.eq5d.mi<-data_wide[,c("eq5d_1","eq5d_2","eq5d_3","trt","id")]
data_long.eq5d.mi<-reshape(data_wide_long.eq5d.mi, varying = c("eq5d_1","eq5d_2","eq5d_3"),
                           direction = "long",idvar = "id",sep = "_")
data_long.eq5d.mi$time_f2<-ifelse(data_long.eq5d.mi$time==2,1,0)
data_long.eq5d.mi$time_f3<-ifelse(data_long.eq5d.mi$time==3,1,0)
data_long.eq5d.mi$trt_time_f2<-data_long.eq5d.mi$trt*data_long.eq5d.mi$time_f2
data_long.eq5d.mi$trt_time_f3<-data_long.eq5d.mi$trt*data_long.eq5d.mi$time_f3
data_long.eq5d.mi$time_f<- factor(data_long.eq5d.mi$time)
data_long.eq5d.mi$time<-NULL
data_long.eq5d.mi$time_f2<-NULL
data_long.eq5d.mi$time_f3<-NULL
mice.long.qalys.ini<-mice(data_long.eq5d.mi,maxit = 0, method = "norm.nob", print = FALSE)
meth<- mice.long.qalys.ini$meth
pred <-mice.long.qalys.ini$pred
pred[,c("id")] <- 0
pred["id",] <- 0
pred[c("trt"),] <- 0
pred[c("time_f","time_f"),] <- 0
mice.long.qalys<-mice(data_long.eq5d.mi, meth=meth, pred=pred, maxit=100,m=100, seed=345, print=F)
mice.fit.long.qalys<-with(data = mice.long.qalys, exp = lme(eq5d ~ -1 + time_f + trt_time_f2 + trt_time_f3, random = ~ 1 | id, 
                                                            method = "ML", correlation = corSymm(form =~ as.numeric(time_f)|id),
                                                            weights=varIdent(form=~1|as.numeric(time_f)),
                                                            na.action = na.omit))
mice.long_u_means<-emmeans(mice.fit.long.qalys,~-1 + time_f + trt_time_f2 + trt_time_f3)
mice.long_qaly_means<-contrast(mice.long_u_means,list(mu0 = c(13/104,13/104 + 26/104,26/104,0,0,0,0,0,0,0,0,0), 
                                                      mu1=c(13/104,0,0,0,13/104 + 26/104,0,0,0,26/104,0,0,0)))
mice.long_qaly_means_ci<-confint(mice.long_qaly_means)
mice.long_qaly_diff<-contrast(mice.long_u_means,list(diff = c(0,-13/104 - 26/104,-26/104,0,13/104 + 26/104,0,0,0,26/104,0,0,0)))
mice.long_qaly_diff_ci<-confint(mice.long_qaly_diff)

########################### multiple imputation (with costs as auxiliary variables)
library(mice)
#MICE
#mice-lm
data_wide.mi<-data_wide[,c("trt","eq5d0","eq5d13","eq5d39","qaly","cost0","cost13","cost39")]
mice.wide.qalys.ini<-mice(data_wide.mi,maxit = 0, method = "norm.nob", print = FALSE)
meth<- mice.wide.qalys.ini$meth
pred <-mice.wide.qalys.ini$pred
pred[c("trt"), c("qaly","eq5d0","eq5d13","eq5d39","cost0","cost13","cost39")] <- 0
pred[c("eq5d0"), c("qaly","eq5d13","eq5d39","trt","cost13","cost39")] <- 0
pred[c("cost0"), c("qaly","eq5d13","eq5d39","trt","cost13","cost39")] <- 0
pred["eq5d13", c("qaly")] <- 0
pred["cost13", c("qaly")] <- 0
#pred["eq5d13", c("qaly","eq5d39","cost39")] <- 0
#pred["cost13", c("qaly","eq5d39","cost39")] <- 0
pred["eq5d39", c("qaly")] <- 0
pred["cost39", c("qaly")] <- 0
pred[c("qaly"), c("trt","eq5d0","eq5d13","eq5d39","cost0","cost13","cost39")] <- 0
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
#lmm
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
#lm
library(emmeans)
lm.cc_totc<-lm(tcost~trt_f+cost0,data = data_wide.cc)
lm.cc_totc_means<-emmeans(lm.cc_totc,~trt_f)
lm.cc_totc_diff<-contrast(lm.cc_totc_means,list(diff=c(-1,1)))
lm.cc_totc_diff_ci<-confint(lm.cc_totc_diff)

#lmm
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

########################### multiple imputation (without utilities as auxiliary variables)
library(mice)
#MICE
#lm
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


########################### multiple imputation (with utilities as auxiliary variables)
library(mice)
#MICE
#lm
data_wide.mi<-data_wide[,c("trt","cost0","cost13","cost39","tcost","id","eq5d0","eq5d13","eq5d39")]
mice.wide.totcs.ini<-mice(data_wide.mi,maxit = 0, method = "norm.nob", print = FALSE)
meth<- mice.wide.totcs.ini$meth
pred <-mice.wide.totcs.ini$pred
pred[,"id"] <- 0
pred["id",] <- 0
pred[c("trt"), c("tcost","cost0","cost13","cost39","eq5d0","eq5d13","eq5d39")] <- 0
pred[c("cost0"), c("tcost","cost13","cost39","trt","eq5d13","eq5d39")] <- 0
pred[c("eq5d0"), c("tcost","cost13","cost39","trt","eq5d13","eq5d39")] <- 0
pred["cost13", c("tcost")] <- 0
#pred["cost13", c("tcost","cost39","eq5d39)] <- 0
pred["cost39", c("tcost")] <- 0
pred[c("tcost"), c("trt","cost0","cost13","cost39","eq5d0","eq5d13","eq5d39")] <- 0
meth["tcost"]<- "~ I(cost13 + cost39)"
mice.wide.totcs<-mice(data_wide.mi, meth=meth, pred=pred, maxit=100,m=500, seed=345, print=F)
mice.fit.wide.totcs<-with(data = mice.wide.totcs, exp = lm(tcost ~ trt + cost0))
mice_totcs_means<-emmeans(mice.fit.wide.totcs, ~trt)
mice_totc_diff<-contrast(mice_totcs_means,list(diff=c(-1,1)))
mice_totc_diff_ci<-confint(mice_totc_diff)

mice.fit.wide.c1<-with(data = mice.wide.totcs, exp = lm(cost0 ~ 1))
mice_c1_means<-emmeans(mice.fit.wide.c1, ~1)
mice.fit.wide.c2<-with(data = mice.wide.totcs, exp = lm(cost13 ~ trt + cost0))
mice_c2_means<-emmeans(mice.fit.wide.c2, ~trt)
mice_c2_diff<-contrast(mice_c2_means,list(diff=c(-1,1)))
mice_c2_diff_ci<-confint(mice_c2_diff)
mice.fit.wide.c3<-with(data = mice.wide.totcs, exp = lm(cost39 ~ trt + cost0))
mice_c3_means<-emmeans(mice.fit.wide.c3, ~trt)
mice_c3_diff<-contrast(mice_c3_means,list(diff=c(-1,1)))
mice_c3_diff_ci<-confint(mice_c3_diff)


########################### LMM
#lmm
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


################################### SAVE RESULTS FOR PLOTTING

ggdf_mean <- data.frame(
  arm = c("Control","Intervention","Incremental","Control","Intervention","Incremental",
          "Control","Intervention","Incremental","Control","Intervention","Incremental","Control","Intervention","Incremental",
          "Control","Intervention","Incremental"),
  outcome = c("QALYs","QALYs","QALYs","QALYs","QALYs","QALYs","QALYs","QALYs","QALYs",
              "TCs","TCs","TCs","TCs","TCs","TCs",
              "TCs","TCs","TCs"),
  method = c("CCA","CCA","CCA","MI","MI","MI",
             "LMM","LMM","LMM","CCA","CCA","CCA","MI","MI","MI",
             "LMM","LMM","LMM"),
  value = c(0.560,0.577,0.0166,0.540,0.563,0.0223,0.546,0.568,0.022,
            3464,3792,328,3440,3984,545,3348,3898,550),
  lower = c(0.529,0.543,-0.03,0.513,0.534,-0.0177,0.527,0.548,mu_e_diff_mi_ci_low,
            1863,2069,-2052,2267,2692,-1192,2167,2663,-1164),
  upper = c(0.592,0.611,0.0631,0.567,0.592,0.0624,0.565,0.588,mu_e_diff_mi_ci_up,
            5064,5516,2709,4612,5277,2281,4529,5133,2263)
)

ggdf_mean$arm <- factor(ggdf_mean$arm, levels = c("Control","Intervention","Incremental"))
ggdf_mean$outcome <- factor(ggdf_mean$outcome, levels = c("QALYs","TCs"))
ggdf_mean$method <- factor(ggdf_mean$method, levels = c("CCA","MI","LMM"))

ggdf_mean_marginal_e<-ggdf_mean[ggdf_mean$arm%in%c("Control","Intervention") & ggdf_mean$outcome%in%c("QALYs"),]
ggdf_mean_incremental_e<-ggdf_mean[ggdf_mean$arm%in%c("Incremental")& ggdf_mean$outcome%in%c("QALYs"),]

ggdf_mean_marginal_c<-ggdf_mean[ggdf_mean$arm%in%c("Control","Intervention") & ggdf_mean$outcome%in%c("TCs"),]
ggdf_mean_incremental_c<-ggdf_mean[ggdf_mean$arm%in%c("Incremental")& ggdf_mean$outcome%in%c("TCs"),]

#ggplot
library(ggplot2)

ggplot(ggdf_mean_marginal_e, aes(x = method, y = value) ) +
  geom_point(position = position_dodge(width = .75) ) +
  facet_grid(~ arm) +
  geom_errorbar( aes(ymin = lower, ymax = upper, width = 0.1 ),
                 position = position_dodge(width = .75) ) +
  theme_classic() + labs(y = "Mean QALYs", x = "Method") 

ggplot(ggdf_mean_marginal_c, aes(x = method, y = value) ) +
  geom_point(position = position_dodge(width = .75) ) +
  facet_grid(~ arm) +
  geom_errorbar( aes(ymin = lower, ymax = upper, width = 0.1 ),
                 position = position_dodge(width = .75) ) +
  theme_classic() + labs(y = "Mean TCs", x = "Method")


ggplot(ggdf_mean_incremental_e, aes(x = method, y = value) ) +
  geom_point(position = position_dodge(width = .75) ) +
  geom_errorbar( aes(ymin = lower, ymax = upper, width = 0.1 ),
                 position = position_dodge(width = .75) ) +
  theme_classic() + labs(y = "Mean QALYs", x = "Method")

ggplot(ggdf_mean_incremental_c, aes(x = method, y = value) ) +
  geom_point(position = position_dodge(width = .75) ) +
  geom_errorbar( aes(ymin = lower, ymax = upper, width = 0.1 ),
                 position = position_dodge(width = .75) ) +
  theme_classic() + labs(y = "Mean TCs", x = "Method")


















#save for plot

ggdf_ue_mean <- data.frame(
  arm = c("Control","Intervention","Control","Intervention","Control","Intervention","Control","Intervention","Control","Intervention",
          "Control","Intervention","Control","Intervention","Control","Intervention","Control","Intervention","Control","Intervention"),
  outcome = c("QALYs","QALYs","u0","u0","u1","u1","u2","u2","QALYs","QALYs",
              "QALYs","QALYs","u0","u0","u1","u1","u2","u2","QALYs","QALYs"),
  method = c("OLS","OLS","LMM","LMM","LMM","LMM","LMM","LMM","LMM","LMM",
             "OLS","OLS","LMM","LMM","LMM","LMM","LMM","LMM","LMM","LMM"),
  type = c("CC","CC","CC","CC","CC","CC","CC","CC","CC","CC",
           "ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL"),
  value = c(OLS_CC_mean_qalys_control[1],OLS_CC_mean_qalys_intervention[1],LMM_CC_mean_eq5d0_control[1],LMM_CC_mean_eq5d0_intervention[1],LMM_CC_mean_eq5d1_control[1],LMM_CC_mean_eq5d1_intervention[1],
            LMM_CC_mean_eq5d2_control[1],LMM_CC_mean_eq5d2_intervention[1],LMM_CC_mean_qalys_control[1],LMM_CC_mean_qalys_intervention[1],
            OLS_ALL_mean_qalys_control[1],OLS_ALL_mean_qalys_intervention[1],LMM_ALL_mean_eq5d0_control[1],LMM_ALL_mean_eq5d0_intervention[1],LMM_ALL_mean_eq5d1_control[1],LMM_ALL_mean_eq5d1_intervention[1],
            LMM_ALL_mean_eq5d2_control[1],LMM_ALL_mean_eq5d2_intervention[1],LMM_ALL_mean_qalys_control[1],LMM_ALL_mean_qalys_intervention[1]),
  lower = c(OLS_CC_mean_qalys_control[2],OLS_CC_mean_qalys_intervention[2],LMM_CC_mean_eq5d0_control[2],LMM_CC_mean_eq5d0_intervention[2],LMM_CC_mean_eq5d1_control[2],LMM_CC_mean_eq5d1_intervention[2],
            LMM_CC_mean_eq5d2_control[2],LMM_CC_mean_eq5d2_intervention[2],LMM_CC_mean_qalys_control[2],LMM_CC_mean_qalys_intervention[2],
            OLS_ALL_mean_qalys_control[2],OLS_ALL_mean_qalys_intervention[2],LMM_ALL_mean_eq5d0_control[2],LMM_ALL_mean_eq5d0_intervention[2],LMM_ALL_mean_eq5d1_control[2],LMM_ALL_mean_eq5d1_intervention[2],
            LMM_ALL_mean_eq5d2_control[2],LMM_ALL_mean_eq5d2_intervention[2],LMM_ALL_mean_qalys_control[2],LMM_ALL_mean_qalys_intervention[2]),
  upper = c(OLS_CC_mean_qalys_control[3],OLS_CC_mean_qalys_intervention[3],LMM_CC_mean_eq5d0_control[3],LMM_CC_mean_eq5d0_intervention[3],LMM_CC_mean_eq5d1_control[3],LMM_CC_mean_eq5d1_intervention[3],
            LMM_CC_mean_eq5d2_control[3],LMM_CC_mean_eq5d2_intervention[3],LMM_CC_mean_qalys_control[3],LMM_CC_mean_qalys_intervention[3],
            OLS_ALL_mean_qalys_control[3],OLS_ALL_mean_qalys_intervention[3],LMM_ALL_mean_eq5d0_control[3],LMM_ALL_mean_eq5d0_intervention[3],LMM_ALL_mean_eq5d1_control[3],LMM_ALL_mean_eq5d1_intervention[3],
            LMM_ALL_mean_eq5d2_control[3],LMM_ALL_mean_eq5d2_intervention[3],LMM_ALL_mean_qalys_control[3],LMM_ALL_mean_qalys_intervention[3])
)


ggdf_ue_diff <- data.frame(
  outcome = c("QALYs","u0","u1","u2","QALYs",
              "QALYs","u0","u1","u2","QALYs"),
  method = c("OLS","LMM","LMM","LMM","LMM",
             "OLS","LMM","LMM","LMM","LMM"),
  type = c("CC","CC","CC","CC","CC",
           "ALL","ALL","ALL","ALL","ALL"),
  value = c(OLS_CC_mean_qalys_diff[1],LMM_CC_mean_eq5d0_diff[1],LMM_CC_mean_eq5d1_diff[1],
            LMM_CC_mean_eq5d2_diff[1],LMM_CC_mean_qalys_diff[1],
            OLS_ALL_mean_qalys_diff[1],LMM_ALL_mean_eq5d0_diff[1],LMM_ALL_mean_eq5d1_diff[1],
            LMM_ALL_mean_eq5d2_diff[1],LMM_ALL_mean_qalys_diff[1]),
  lower = c(OLS_CC_mean_qalys_diff[2],LMM_CC_mean_eq5d0_diff[2],LMM_CC_mean_eq5d1_diff[2],
            LMM_CC_mean_eq5d2_diff[2],LMM_CC_mean_qalys_diff[2],
            OLS_ALL_mean_qalys_diff[2],LMM_ALL_mean_eq5d0_diff[2],LMM_ALL_mean_eq5d1_diff[2],
            LMM_ALL_mean_eq5d2_diff[2],LMM_ALL_mean_qalys_diff[2]),
  upper = c(OLS_CC_mean_qalys_diff[3],LMM_CC_mean_eq5d0_diff[3],LMM_CC_mean_eq5d1_diff[3],
            LMM_CC_mean_eq5d2_diff[3],LMM_CC_mean_qalys_diff[3],
            OLS_ALL_mean_qalys_diff[3],LMM_ALL_mean_eq5d0_diff[3],LMM_ALL_mean_eq5d1_diff[3],
            LMM_ALL_mean_eq5d2_diff[3],LMM_ALL_mean_qalys_diff[3])
)


ggdf_ctc_mean <- data.frame(
  arm = c("Control","Intervention","Control","Intervention","Control","Intervention","Control","Intervention",
          "Control","Intervention","Control","Intervention","Control","Intervention","Control","Intervention"),
  outcome = c("Total costs","Total costs","c1","c1","c2","c2","Total costs","Total costs",
              "Total costs","Total costs","c1","c1","c2","c2","Total costs","Total costs"),
  method = c("OLS","OLS","LMM","LMM","LMM","LMM","LMM","LMM",
             "OLS","OLS","LMM","LMM","LMM","LMM","LMM","LMM"),
  type = c("CC","CC","CC","CC","CC","CC","CC","CC",
           "ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL"),
  value = c(OLS_CC_mean_tcosts_control[1],OLS_CC_mean_tcosts_intervention[1],LMM_CC_mean_cost1_control[1],LMM_CC_mean_cost1_intervention[1],
            LMM_CC_mean_cost2_control[1],LMM_CC_mean_cost2_intervention[1],LMM_CC_mean_tcosts_control[1],LMM_CC_mean_tcosts_intervention[1],
            OLS_ALL_mean_tcosts_control[1],OLS_ALL_mean_tcosts_intervention[1],LMM_ALL_mean_cost1_control[1],LMM_ALL_mean_cost1_intervention[1],
            LMM_ALL_mean_cost2_control[1],LMM_ALL_mean_cost2_intervention[1],LMM_ALL_mean_tcosts_control[1],LMM_ALL_mean_tcosts_intervention[1]),
  lower = c(OLS_CC_mean_tcosts_control[2],OLS_CC_mean_tcosts_intervention[2],LMM_CC_mean_cost1_control[2],LMM_CC_mean_cost1_intervention[2],
            LMM_CC_mean_cost2_control[2],LMM_CC_mean_cost2_intervention[2],LMM_CC_mean_tcosts_control[2],LMM_CC_mean_tcosts_intervention[2],
            OLS_ALL_mean_tcosts_control[2],OLS_ALL_mean_tcosts_intervention[2],LMM_ALL_mean_cost1_control[2],LMM_ALL_mean_cost1_intervention[2],
            LMM_ALL_mean_cost2_control[2],LMM_ALL_mean_cost2_intervention[2],LMM_ALL_mean_tcosts_control[2],LMM_ALL_mean_tcosts_intervention[2]),
  upper = c(OLS_CC_mean_tcosts_control[3],OLS_CC_mean_tcosts_intervention[3],LMM_CC_mean_cost1_control[3],LMM_CC_mean_cost1_intervention[3],
            LMM_CC_mean_cost2_control[3],LMM_CC_mean_cost2_intervention[3],LMM_CC_mean_tcosts_control[3],LMM_CC_mean_tcosts_intervention[3],
            OLS_ALL_mean_tcosts_control[3],OLS_ALL_mean_tcosts_intervention[3],LMM_ALL_mean_cost1_control[3],LMM_ALL_mean_cost1_intervention[3],
            LMM_ALL_mean_cost2_control[3],LMM_ALL_mean_cost2_intervention[3],LMM_ALL_mean_tcosts_control[3],LMM_ALL_mean_tcosts_intervention[3])
)


ggdf_ctc_diff <- data.frame(
  outcome = c("Total costs","c1","c2","Total costs",
              "Total costs","c1","c2","Total costs"),
  method = c("OLS","LMM","LMM","LMM",
             "OLS","LMM","LMM","LMM"),
  type = c("CC","CC","CC","CC",
           "ALL","ALL","ALL","ALL"),
  value = c(OLS_CC_mean_tcosts_diff[1],LMM_CC_mean_cost1_diff[1],
            LMM_CC_mean_cost2_diff[1],LMM_CC_mean_tcosts_diff[1],
            OLS_ALL_mean_tcosts_diff[1],LMM_ALL_mean_cost1_diff[1],
            LMM_ALL_mean_cost2_diff[1],LMM_ALL_mean_tcosts_diff[1]),
  lower = c(OLS_CC_mean_tcosts_diff[2],LMM_CC_mean_cost1_diff[2],
            LMM_CC_mean_cost2_diff[2],LMM_CC_mean_tcosts_diff[2],
            OLS_ALL_mean_tcosts_diff[2],LMM_ALL_mean_cost1_diff[2],
            LMM_ALL_mean_cost2_diff[2],LMM_ALL_mean_tcosts_diff[2]),
  upper = c(OLS_CC_mean_tcosts_diff[3],LMM_CC_mean_cost1_diff[3],
            LMM_CC_mean_cost2_diff[3],LMM_CC_mean_tcosts_diff[3],
            OLS_ALL_mean_tcosts_diff[3],LMM_ALL_mean_cost1_diff[3],
            LMM_ALL_mean_cost2_diff[3],LMM_ALL_mean_tcosts_diff[3])
)

#separate plots
ggdf_u_mean <- ggdf_ue_mean[!ggdf_ue_mean$outcome %in% c("QALYs"),]
ggdf_u_mean$arm <- factor(ggdf_u_mean$arm, levels = c("Control","Intervention"))
ggdf_u_mean$outcome <- factor(ggdf_u_mean$outcome, levels = c("u0","u1","u2"))
ggdf_u_mean$method <- factor(ggdf_u_mean$method, levels = c("OLS","LMM"))
ggdf_u_mean$type <- factor(ggdf_u_mean$type, levels = c("CC","ALL"))

ggdf_e_mean <- ggdf_ue_mean[ggdf_ue_mean$outcome %in% c("QALYs"),]
ggdf_e_mean$arm <- factor(ggdf_e_mean$arm, levels = c("Control","Intervention"))
ggdf_e_mean$method <- factor(ggdf_e_mean$method, levels = c("OLS","LMM"))
ggdf_e_mean$type <- factor(ggdf_e_mean$type, levels = c("CC","ALL"))

ggdf_u_diff <- ggdf_ue_diff[!ggdf_ue_diff$outcome %in% c("QALYs"),]
ggdf_u_diff$outcome <- factor(ggdf_u_diff$outcome, levels = c("u0","u1","u2"))
ggdf_u_diff$method <- factor(ggdf_u_diff$method, levels = c("OLS","LMM"))
ggdf_u_diff$type <- factor(ggdf_u_diff$type, levels = c("CC","ALL"))

ggdf_e_diff <- ggdf_ue_diff[ggdf_ue_diff$outcome %in% c("QALYs"),]
ggdf_e_diff$method <- factor(ggdf_e_diff$method, levels = c("OLS","LMM"))
ggdf_e_diff$type <- factor(ggdf_e_diff$type, levels = c("CC","ALL"))

ggdf_c_mean <- ggdf_ctc_mean[!ggdf_ctc_mean$outcome %in% c("Total costs"),]
ggdf_c_mean$arm <- factor(ggdf_c_mean$arm, levels = c("Control","Intervention"))
ggdf_c_mean$outcome <- factor(ggdf_c_mean$outcome, levels = c("c1","c2"))
ggdf_c_mean$method <- factor(ggdf_c_mean$method, levels = c("OLS","LMM"))
ggdf_c_mean$type <- factor(ggdf_c_mean$type, levels = c("CC","ALL"))

ggdf_tc_mean <- ggdf_ctc_mean[ggdf_ctc_mean$outcome %in% c("Total costs"),]
ggdf_tc_mean$arm <- factor(ggdf_tc_mean$arm, levels = c("Control","Intervention"))
ggdf_tc_mean$method <- factor(ggdf_tc_mean$method, levels = c("OLS","LMM"))
ggdf_tc_mean$type <- factor(ggdf_tc_mean$type, levels = c("CC","ALL"))

ggdf_c_diff <- ggdf_ctc_diff[!ggdf_ctc_diff$outcome %in% c("Total costs"),]
ggdf_c_diff$outcome <- factor(ggdf_c_diff$outcome, levels = c("c1","c2"))
ggdf_c_diff$method <- factor(ggdf_c_diff$method, levels = c("OLS","LMM"))
ggdf_c_diff$type <- factor(ggdf_c_diff$type, levels = c("CC","ALL"))

ggdf_tc_diff <- ggdf_ctc_diff[ggdf_ctc_diff$outcome %in% c("Total costs"),]
ggdf_tc_diff$method <- factor(ggdf_tc_diff$method, levels = c("OLS","LMM"))
ggdf_tc_diff$type <- factor(ggdf_tc_diff$type, levels = c("CC","ALL"))

#ggplot
library(ggplot2)

ggplot(ggdf_u_mean, aes(x = arm, y = value) ) +
  geom_point(position = position_dodge(width = .75) ) +
  facet_grid(type ~ outcome) +
  geom_errorbar( aes(ymin = lower, ymax = upper, width = 0.1 ),
                 position = position_dodge(width = .75) ) +
  theme_classic() + labs(y = "Mean utility", x = "arm") 

ggplot(ggdf_u_diff, aes(x = type, y = value) ) +
  geom_point(position = position_dodge(width = .75) ) +
  facet_grid(. ~ outcome) +
  geom_errorbar( aes(ymin = lower, ymax = upper, width = 0.1 ),
                 position = position_dodge(width = .75) ) +
  theme_classic() + labs(y = "Mean utility difference", x = "type") +
  geom_hline(yintercept=0, linetype="dashed")
  
ggplot(ggdf_e_mean, aes(x = arm, y = value) ) +
  geom_point(position = position_dodge(width = .75) ) +
  facet_grid(type ~ method) +
  geom_errorbar( aes(ymin = lower, ymax = upper, width = 0.1 ),
                 position = position_dodge(width = .75) ) +
  theme_classic() + labs(y = "Mean QALYs", x = "arm") 

ggplot(ggdf_e_diff, aes(x = method, y = value) ) +
  geom_point(position = position_dodge(width = .75) ) +
  facet_grid(. ~ type) +
  geom_errorbar( aes(ymin = lower, ymax = upper, width = 0.1 ),
                 position = position_dodge(width = .75) ) +
  theme_classic() + labs(y = "Mean QALYs difference", x = "method") +
  geom_hline(yintercept=0, linetype="dashed")


ggplot(ggdf_c_mean, aes(x = arm, y = value) ) +
  geom_point(position = position_dodge(width = .75) ) +
  facet_grid(type ~ outcome) +
  geom_errorbar( aes(ymin = lower, ymax = upper, width = 0.1 ),
                 position = position_dodge(width = .75) ) +
  theme_classic() + labs(y = "Mean cost", x = "arm") 

ggplot(ggdf_c_diff, aes(x = type, y = value) ) +
  geom_point(position = position_dodge(width = .75) ) +
  facet_grid(. ~ outcome) +
  geom_errorbar( aes(ymin = lower, ymax = upper, width = 0.1 ),
                 position = position_dodge(width = .75) ) +
  theme_classic() + labs(y = "Mean cost difference", x = "type") +
  geom_hline(yintercept=0, linetype="dashed")

ggplot(ggdf_tc_mean, aes(x = arm, y = value) ) +
  geom_point(position = position_dodge(width = .75) ) +
  facet_grid(type ~ method) +
  geom_errorbar( aes(ymin = lower, ymax = upper, width = 0.1 ),
                 position = position_dodge(width = .75) ) +
  theme_classic() + labs(y = "Mean total costs", x = "arm") 

ggplot(ggdf_tc_diff, aes(x = method, y = value) ) +
  geom_point(position = position_dodge(width = .75) ) +
  facet_grid(. ~ type) +
  geom_errorbar( aes(ymin = lower, ymax = upper, width = 0.1 ),
                 position = position_dodge(width = .75) ) +
  theme_classic() + labs(y = "Mean total costs difference", x = "method") +
  geom_hline(yintercept=0, linetype="dashed")

res_univ_lmm <- list("plot"= list("u"=ggdf_u_diff,"e"=ggdf_e_diff, 
                                  "c"=ggdf_c_diff, "tc"=ggdf_tc_diff))
#saveRDS(res_univ_lmm,"res_univ_lmm.rds")


#table of results
table.res_u.mle<-matrix(NA,9,6)
colnames(table.res_u.mle)<-c("u0 (CC)","u0 (ALL)","u1 (CC)","u1 (ALL)","u2 (CC)","u2 (ALL)")
rownames(table.res_u.mle)<-c("con - est","con - LB","con - UB", "int - est","int - LB","int - UB", "diff - est","diff - LB","diff - UB")
table.res_u.mle[1,]<-c(LMM_CC_mean_eq5d0_control[1], LMM_ALL_mean_eq5d0_control[1], LMM_CC_mean_eq5d1_control[1], LMM_ALL_mean_eq5d1_control[1], LMM_CC_mean_eq5d2_control[1], LMM_ALL_mean_eq5d2_control[1])
table.res_u.mle[2,]<-c(LMM_CC_mean_eq5d0_control[2], LMM_ALL_mean_eq5d0_control[2], LMM_CC_mean_eq5d1_control[2], LMM_ALL_mean_eq5d1_control[2], LMM_CC_mean_eq5d2_control[2], LMM_ALL_mean_eq5d2_control[2])
table.res_u.mle[3,]<-c(LMM_CC_mean_eq5d0_control[3], LMM_ALL_mean_eq5d0_control[3], LMM_CC_mean_eq5d1_control[3], LMM_ALL_mean_eq5d1_control[3], LMM_CC_mean_eq5d2_control[3], LMM_ALL_mean_eq5d2_control[3])
table.res_u.mle[4,]<-c(LMM_CC_mean_eq5d0_intervention[1], LMM_ALL_mean_eq5d0_intervention[1], LMM_CC_mean_eq5d1_intervention[1], LMM_ALL_mean_eq5d1_intervention[1], LMM_CC_mean_eq5d2_intervention[1], LMM_ALL_mean_eq5d2_intervention[1])
table.res_u.mle[5,]<-c(LMM_CC_mean_eq5d0_intervention[2], LMM_ALL_mean_eq5d0_intervention[2], LMM_CC_mean_eq5d1_intervention[2], LMM_ALL_mean_eq5d1_intervention[2], LMM_CC_mean_eq5d2_intervention[2], LMM_ALL_mean_eq5d2_intervention[2])
table.res_u.mle[6,]<-c(LMM_CC_mean_eq5d0_intervention[3], LMM_ALL_mean_eq5d0_intervention[3], LMM_CC_mean_eq5d1_intervention[3], LMM_ALL_mean_eq5d1_intervention[3], LMM_CC_mean_eq5d2_intervention[3], LMM_ALL_mean_eq5d2_intervention[3])
table.res_u.mle[7,]<-c(LMM_CC_mean_eq5d0_diff[1], LMM_ALL_mean_eq5d0_diff[1], LMM_CC_mean_eq5d1_diff[1], LMM_ALL_mean_eq5d1_diff[1], LMM_CC_mean_eq5d2_diff[1], LMM_ALL_mean_eq5d2_diff[1])
table.res_u.mle[8,]<-c(LMM_CC_mean_eq5d0_diff[2], LMM_ALL_mean_eq5d0_diff[2], LMM_CC_mean_eq5d1_diff[2], LMM_ALL_mean_eq5d1_diff[2], LMM_CC_mean_eq5d2_diff[2], LMM_ALL_mean_eq5d2_diff[2])
table.res_u.mle[9,]<-c(LMM_CC_mean_eq5d0_diff[3], LMM_ALL_mean_eq5d0_diff[3], LMM_CC_mean_eq5d1_diff[3], LMM_ALL_mean_eq5d1_diff[3], LMM_CC_mean_eq5d2_diff[3], LMM_ALL_mean_eq5d2_diff[3])
table.res_u.mle<-round(table.res_u.mle, digits = 3)


table.res_e.mle<-matrix(NA,9,4)
colnames(table.res_e.mle)<-c("OLS (CC)","LMM (CC)","OLS (ALL)","LMM (ALL)")
rownames(table.res_e.mle)<-c("con - est","con - LB","con - UB", "int - est","int - LB","int - UB", "diff - est","diff - LB","diff - UB")
table.res_e.mle[1,]<-c(OLS_CC_mean_qalys_control[1], LMM_CC_mean_qalys_control[1],OLS_ALL_mean_qalys_control[1], LMM_ALL_mean_qalys_control[1])
table.res_e.mle[2,]<-c(OLS_CC_mean_qalys_control[2], LMM_CC_mean_qalys_control[2],OLS_ALL_mean_qalys_control[2], LMM_ALL_mean_qalys_control[2])
table.res_e.mle[3,]<-c(OLS_CC_mean_qalys_control[3], LMM_CC_mean_qalys_control[3],OLS_ALL_mean_qalys_control[3], LMM_ALL_mean_qalys_control[3])
table.res_e.mle[4,]<-c(OLS_CC_mean_qalys_intervention[1], LMM_CC_mean_qalys_intervention[1],OLS_ALL_mean_qalys_intervention[1], LMM_ALL_mean_qalys_intervention[1])
table.res_e.mle[5,]<-c(OLS_CC_mean_qalys_intervention[2], LMM_CC_mean_qalys_intervention[2],OLS_ALL_mean_qalys_intervention[2], LMM_ALL_mean_qalys_intervention[2])
table.res_e.mle[6,]<-c(OLS_CC_mean_qalys_intervention[3], LMM_CC_mean_qalys_intervention[3],OLS_ALL_mean_qalys_intervention[3], LMM_ALL_mean_qalys_intervention[3])
table.res_e.mle[7,]<-c(OLS_CC_mean_qalys_diff[1], LMM_CC_mean_qalys_diff[1],OLS_ALL_mean_qalys_diff[1], LMM_ALL_mean_qalys_diff[1])
table.res_e.mle[8,]<-c(OLS_CC_mean_qalys_diff[2], LMM_CC_mean_qalys_diff[2],OLS_ALL_mean_qalys_diff[2], LMM_ALL_mean_qalys_diff[2])
table.res_e.mle[9,]<-c(OLS_CC_mean_qalys_diff[3], LMM_CC_mean_qalys_diff[3],OLS_ALL_mean_qalys_diff[3], LMM_ALL_mean_qalys_diff[3])
table.res_e.mle<-round(table.res_e.mle, digits = 3)


table.res_c.mle<-matrix(NA,9,4)
colnames(table.res_c.mle)<-c("c1 (CC)","c1 (ALL)","c2 (CC)","c2 (ALL)")
rownames(table.res_c.mle)<-c("con - est","con - LB","con - UB", "int - est","int - LB","int - UB", "diff - est","diff - LB","diff - UB")
table.res_c.mle[1,]<-c(LMM_CC_mean_cost1_control[1], LMM_ALL_mean_cost1_control[1], LMM_CC_mean_cost2_control[1], LMM_ALL_mean_cost2_control[1])
table.res_c.mle[2,]<-c(LMM_CC_mean_cost1_control[2], LMM_ALL_mean_cost1_control[2], LMM_CC_mean_cost2_control[2], LMM_ALL_mean_cost2_control[2])
table.res_c.mle[3,]<-c(LMM_CC_mean_cost1_control[3], LMM_ALL_mean_cost1_control[3], LMM_CC_mean_cost2_control[3], LMM_ALL_mean_cost2_control[3])
table.res_c.mle[4,]<-c(LMM_CC_mean_cost1_intervention[1], LMM_ALL_mean_cost1_intervention[1], LMM_CC_mean_cost2_intervention[1], LMM_ALL_mean_cost2_intervention[1])
table.res_c.mle[5,]<-c(LMM_CC_mean_cost1_intervention[2], LMM_ALL_mean_cost1_intervention[2], LMM_CC_mean_cost2_intervention[2], LMM_ALL_mean_cost2_intervention[2])
table.res_c.mle[6,]<-c(LMM_CC_mean_cost1_intervention[3], LMM_ALL_mean_cost1_intervention[3], LMM_CC_mean_cost2_intervention[3], LMM_ALL_mean_cost2_intervention[3])
table.res_c.mle[7,]<-c(LMM_CC_mean_cost1_diff[1], LMM_ALL_mean_cost1_diff[1], LMM_CC_mean_cost2_diff[1], LMM_ALL_mean_cost2_diff[1])
table.res_c.mle[8,]<-c(LMM_CC_mean_cost1_diff[2], LMM_ALL_mean_cost1_diff[2], LMM_CC_mean_cost2_diff[2], LMM_ALL_mean_cost2_diff[2])
table.res_c.mle[9,]<-c(LMM_CC_mean_cost1_diff[3], LMM_ALL_mean_cost1_diff[3], LMM_CC_mean_cost2_diff[3], LMM_ALL_mean_cost2_diff[3])
table.res_c.mle<-round(table.res_c.mle, digits = 0)


table.res_tc.mle<-matrix(NA,9,4)
colnames(table.res_tc.mle)<-c("OLS (CC)","LMM (CC)","OLS (ALL)","LMM (ALL)")
rownames(table.res_tc.mle)<-c("con - est","con - LB","con - UB", "int - est","int - LB","int - UB", "diff - est","diff - LB","diff - UB")
table.res_tc.mle[1,]<-c(OLS_CC_mean_tcosts_control[1], LMM_CC_mean_tcosts_control[1],OLS_ALL_mean_tcosts_control[1], LMM_ALL_mean_tcosts_control[1])
table.res_tc.mle[2,]<-c(OLS_CC_mean_tcosts_control[2], LMM_CC_mean_tcosts_control[2],OLS_ALL_mean_tcosts_control[2], LMM_ALL_mean_tcosts_control[2])
table.res_tc.mle[3,]<-c(OLS_CC_mean_tcosts_control[3], LMM_CC_mean_tcosts_control[3],OLS_ALL_mean_tcosts_control[3], LMM_ALL_mean_tcosts_control[3])
table.res_tc.mle[4,]<-c(OLS_CC_mean_tcosts_intervention[1], LMM_CC_mean_tcosts_intervention[1],OLS_ALL_mean_tcosts_intervention[1], LMM_ALL_mean_tcosts_intervention[1])
table.res_tc.mle[5,]<-c(OLS_CC_mean_tcosts_intervention[2], LMM_CC_mean_tcosts_intervention[2],OLS_ALL_mean_tcosts_intervention[2], LMM_ALL_mean_tcosts_intervention[2])
table.res_tc.mle[6,]<-c(OLS_CC_mean_tcosts_intervention[3], LMM_CC_mean_tcosts_intervention[3],OLS_ALL_mean_tcosts_intervention[3], LMM_ALL_mean_tcosts_intervention[3])
table.res_tc.mle[7,]<-c(OLS_CC_mean_tcosts_diff[1], LMM_CC_mean_tcosts_diff[1],OLS_ALL_mean_tcosts_diff[1], LMM_ALL_mean_tcosts_diff[1])
table.res_tc.mle[8,]<-c(OLS_CC_mean_tcosts_diff[2], LMM_CC_mean_tcosts_diff[2],OLS_ALL_mean_tcosts_diff[2], LMM_ALL_mean_tcosts_diff[2])
table.res_tc.mle[9,]<-c(OLS_CC_mean_tcosts_diff[3], LMM_CC_mean_tcosts_diff[3],OLS_ALL_mean_tcosts_diff[3], LMM_ALL_mean_tcosts_diff[3])
table.res_tc.mle<-round(table.res_tc.mle, digits = 0)


library(xtable)
xtable(table.res_u.mle)
xtable(table.res_e.mle)

xtable(table.res_c.mle, digits = 0)
xtable(table.res_tc.mle, digits = 0)


######################## plot only mean and delta for QALYs and total costs

ggdf_e_mean_plot<-ggdf_e_mean[ggdf_e_mean$type=="ALL" & ggdf_e_mean$method=="LMM",]
ggdf_e_mean_plot$arm<- as.character(ggdf_e_mean_plot$arm)
ggdf_e_mean_plot$arm[ggdf_e_mean_plot$arm=="Control"]<-"Placebo"
ggdf_e_mean_plot$arm[ggdf_e_mean_plot$arm=="Intervention"]<-"Mirtazapine"
ggdf_e_mean_plot$arm<-factor(ggdf_e_mean_plot$arm, levels = c("Placebo","Mirtazapine"))

ggplot(ggdf_e_mean_plot, aes(x = arm, y = value, color=arm) ) +
  geom_point(position = position_dodge(width = .75) ) +
  geom_errorbar( aes(ymin = lower, ymax = upper, width = 0.1 ),
                 position = position_dodge(width = .75) ) +
  scale_color_manual(values=c('red','blue')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title = element_blank(), legend.key = element_rect(fill = "white"),
        legend.position = c(0.9,0.1), legend.text = element_text(size = 9)) +
  guides(color = guide_legend(override.aes = list(size = 0.7))) +
  labs(y = "Mean QALYs", x = "") 

ggdf_e_diff_plot <-ggdf_e_diff[ggdf_e_diff$method=="LMM" & ggdf_e_diff$type=="ALL",]

ggplot(ggdf_e_diff_plot, aes(x = method, y = value) ) +
  geom_point(position = position_dodge(width = .75) ) +
  geom_errorbar( aes(ymin = lower, ymax = upper, width = 0.1 ),
                 position = position_dodge(width = .75) ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title = element_blank(), legend.key = element_rect(fill = "white"),
        axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  labs(y = "Mean QALYs difference", x = "") +
  geom_hline(yintercept=0, linetype="dashed")







ggdf_tc_mean_plot<-ggdf_tc_mean[ggdf_tc_mean$type=="ALL" & ggdf_tc_mean$method=="LMM",]
ggdf_tc_mean_plot$arm<- as.character(ggdf_tc_mean_plot$arm)
ggdf_tc_mean_plot$arm[ggdf_tc_mean_plot$arm=="Control"]<-"Placebo"
ggdf_tc_mean_plot$arm[ggdf_tc_mean_plot$arm=="Intervention"]<-"Mirtazapine"
ggdf_tc_mean_plot$arm<-factor(ggdf_tc_mean_plot$arm, levels = c("Placebo","Mirtazapine"))

ggplot(ggdf_tc_mean_plot, aes(x = arm, y = value, color=arm) ) +
  geom_point(position = position_dodge(width = .75) ) +
  geom_errorbar( aes(ymin = lower, ymax = upper, width = 0.1 ),
                 position = position_dodge(width = .75) ) +
  scale_color_manual(values=c('red','blue')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title = element_blank(), legend.key = element_rect(fill = "white"),
        legend.position = c(0.9,0.1), legend.text = element_text(size = 9)) +
  guides(color = guide_legend(override.aes = list(size = 0.7))) +
  labs(y = "Mean total costs", x = "") 

ggdf_tc_diff_plot <-ggdf_tc_diff[ggdf_tc_diff$method=="LMM" & ggdf_tc_diff$type=="ALL",]

ggplot(ggdf_tc_diff_plot, aes(x = method, y = value) ) +
  geom_point(position = position_dodge(width = .75) ) +
  geom_errorbar( aes(ymin = lower, ymax = upper, width = 0.1 ),
                 position = position_dodge(width = .75) ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title = element_blank(), legend.key = element_rect(fill = "white"),
        axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  labs(y = "Mean total costs difference", x = "") +
  geom_hline(yintercept=0, linetype="dashed")






