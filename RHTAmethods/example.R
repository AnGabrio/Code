#this file contains the code for generating artificial data for an hypothetical two-arm
#trial-based CUA presenting all the typical data features for these types of analyses 
#NOTE: the data are simulated in a simplified context to put more emphasis on the objective
#of learning how to use R functions to fit the desired type of statistical methods

# data are simulated with following characteristics:
# one year follow-up with baseline utility and cost values
# outcomes are QALY and TC at 1 year
# data are clustered and treatment allocation depends on clusters
# utility/cost and QALY/TC data are correlated and skewed
# this version of the data have NO missing values

#set parameter values for generating data
n <- 50 #sample size by cluster
J <- 16 #n of clusters
######generate TC and QALY data
beta0_c <- 100
beta1_c <- 20
tau_c <- 10*2
trt_j <- c(rep(0,J/2),rep(1,J/2))
#simulate cluster-level means
set.seed(2345)
phi_cj <- rnorm(J,beta0_c+beta1_c*trt_j, tau_c) #TC cluster means
c_ij <- matrix(NA, nrow = n, ncol = J)
sigma_c <- 15*2
shape_cj <- phi_cj^2/sigma_c^2
scale_cj <- sigma_c^2/phi_cj
#simulate individual-level TC values
for(j in 1:J){
  c_ij[,j] <- rgamma(n,shape = shape_cj[j], scale = scale_cj[j]) #cost data
}
beta0_e <- 0.5
beta1_e <- 0.1
tau_e <- 0.1*2
gamma <- tau_e/tau_c  #parameter capturing correlation between cluster QALY and cost means
#simulate cluster-level means
phi_ej <- rnorm(J,beta0_e+beta1_e+gamma*(phi_cj - (beta0_c+beta1_c*trt_j)), tau_e) #QALY cluster means
e_ij <- matrix(NA, nrow = n, ncol = J)
sigma_e <- 0.1
rho_ec <- 0.5
theta <- rho_ec*(sigma_e/sigma_c) #parameter capturing correlation between QALY and cost
#simulate individual-level values
phi_ej <- phi_ej[j]+theta*(c_ij[,j]-phi_cj[j])
tau_ej <- (phi_ej*(1-phi_ej)/sigma_e^2-1)
shape1_ej <- phi_ej*tau_ej
shape2_ej <- (1-phi_ej)*tau_ej
for(j in 1:J){
  e_ij[,j] <- 1- rbeta(n,shape1 = shape1_ej[j],shape2 = shape2_ej[j]) #QALY data
}
#compute ICC by outcome (cluster variance/total variance)
icc_c <- tau_c^2/(tau_c^2+sigma_c^2)
icc_e <- tau_e^2/(tau_e^2+sigma_e^2)
#generate baseline cost and utility data separately
mu_c0 <- 50 + 0.5*c_ij 
sigma_c0 <- 17
shape_c0 <- mu_c0^2/sigma_c0^2
scale_c0 <- sigma_c0^2/mu_c0
c0_ij <- matrix(NA, nrow = n, ncol = J)
for(j in 1:J){
  c0_ij[,j] <- rgamma(n, shape = shape_c0[,j], scale = scale_c0[,j])
}
mu_u0 <- 0.4 + 0.45*e_ij
sigma_u0 <- 0.1
tau_u0 <- (mu_u0*(1-mu_u0)/sigma_u0^2-1)
shape1_u0 <- mu_u0*tau_u0
shape2_u0 <- (1-mu_u0)*tau_u0
u0_ij <- matrix(NA, nrow = n, ncol = J)
for(j in 1:J){
  u0_ij[,j] <- rbeta(n, shape1 = shape1_u0[,j], shape2 = shape2_u0[,j])
}
#generate dataset
cluster <- rep(1:J, each=n)
TC <- c(c_ij)
QALY <- c(e_ij)
c0 <- c(c0_ij)
u0 <- c(u0_ij)
trt <- ifelse(cluster<=J/2,"old","new")
id <- rep(1:n*J) #individual id number
data.clus.df <- data.frame(id,QALY,TC,u0,c0,trt,cluster)
data.clus.df$trt <- factor(data.clus.df$trt, levels = c("old","new"))
#randomly shuffle rows of the dataset
dataset <- data.clus.df[sample(1:nrow(data.clus.df)), ]

############################ end simulation data

#try different techniques to handle data features

#VERSION 1: adjustment for baseline values
#fit OLS models with baseline values included as covariates
lm_e <- lm(QALY ~ trt + u0, data = dataset)
lm_c <- lm(TC ~ trt + c0, data = dataset)
#summarise output
summary(lm_e)
summary(lm_c)
#get mean QALY/TC by trt group
library(emmeans)
lm_e_em <- emmeans(lm_e, ~ trt)
#see output
lm_e_em
lm_c_em <- emmeans(lm_c, ~ trt)
#see output
lm_c_em
#get difference as New - Old group
diff_new_vs_old <- list("New vs Old" = c(-1, 1))
#compute linear combination
lm_em_delta_e <- contrast(lm_e_em, diff_new_vs_old) 
lm_em_delta_c <- contrast(lm_c_em, diff_new_vs_old) 
#obtain results in terms of confidence intervals
confint(lm_em_delta_e, level = 0.95)
confint(lm_em_delta_c, level = 0.95)

#VERSION 2: accounting for correlation
library(systemfit) 
#fit SUR to QALY and TC
sur_ec <- systemfit(list(QALYreg = QALY~trt + u0, TCreg = TC~trt + c0), 
                    method="SUR", data=dataset)
#summarise output
summary(sur_ec)
#currently not automatic to obtain mean outcome by group using SUR in R since
#function emmeans not compatible (need to compute them manually) - here for simplicity 
#just focus on mean differences (i.e. coefficient of trt)
#obtain results in terms of confidence intervals
confint(sur_ec, level = 0.95)

#combine OLS/SUR with bootstrapping to obtain distributions of statistics of interest 
#such as mean QALY/TC by arm and mean differences between groups
source("code_functions.R") #this loads all pre-built functions stored in the R file called code_functions.R
#get bootstrap results (200 iterations)
set.seed(2345)
boot_res <- boot_ec(dataset, QALYreg = QALY ~ trt + u0,
                    TCreg = TC ~ trt + c0, method = "SUR", B=200)
#summarise output
summary(boot_res$QALY_boot$mu_e_ctr) #mean old arm
summary(boot_res$QALY_boot$mu_e_int) #mean new arm
summary(boot_res$QALY_boot$Delta_e) #mean difference
summary(boot_res$TC_boot$mu_c_ctr) #mean old arm
summary(boot_res$TC_boot$mu_c_int) #mean new arm
summary(boot_res$TC_boot$Delta_c) #mean difference
#compute percentile or BCa CIs
boot_ci_bca <- boot_ci(x = boot_res, method = "BCa")
#summarise output
boot_ci_bca

#VERSION 3: accounting for skewness
#fit GLM models assuming specific link functions and distributions
glm_e <- betareg(QALY ~ trt + u0, data = dataset,
                 link = "logit")
glm_c <- glm(TC ~ trt + c0, data = dataset, 
             family = Gamma(link = "identity"))
#summarise output
summary(glm_e)
summary(glm_c)
#get mean QALY/TC by trt group
glm_e_em <- emmeans(glm_e, ~ trt, type="response")
#see output
glm_e_em
glm_c_em <- emmeans(glm_c, ~ trt, type="response")
#see output
glm_c_em
#get difference as New - Old group
diff_new_vs_old <- list("New vs Old" = c(-1, 1))
#compute linear combination
glm_em_delta_e <- contrast(glm_e_em, diff_new_vs_old) 
glm_em_delta_c <- contrast(glm_c_em, diff_new_vs_old) 
#obtain results in terms of confidence intervals
confint(glm_em_delta_e, level = 0.95)
confint(glm_em_delta_c, level = 0.95)

#combine GLM with bootstrapping to obtain distributions of statistics of interest 
#such as mean QALY/TC by arm and mean differences between groups
#get bootstrap results (200 iterations)
set.seed(2345)
boot_res_glm <- boot_ec_glm(dataset, QALYreg = QALY ~ trt + u0,
                    TCreg = TC ~ trt + c0, QALY_dist = "Beta", 
                    TC_dist = "Gamma", QALY_link = "logit", TC_link = "identity",
                    B=200)
#summarise output
summary(boot_res_glm$QALY_boot$mu_e_ctr) #mean old arm
summary(boot_res_glm$QALY_boot$mu_e_int) #mean new arm
summary(boot_res_glm$QALY_boot$Delta_e) #mean difference
summary(boot_res_glm$TC_boot$mu_c_ctr) #mean old arm
summary(boot_res_glm$TC_boot$mu_c_int) #mean new arm
summary(boot_res_glm$TC_boot$Delta_c) #mean difference
#compute percentile or BCa CIs
boot_ci_bca_glm <- boot_ci_glm(x = boot_res_glm, method = "BCa")
#summarise output
boot_ci_bca_glm

#VERSION 4: accounting for clustering
#fit MLMs 
library(lme4)
library(nlme)
mlm_e <- lme(QALY ~ trt + u0, random = ~1|cluster, 
             data = dataset)
mlm_c <- lme(TC ~ trt + c0, random = ~1|cluster, 
             data = dataset)
#summarise output
summary(mlm_e)
summary(mlm_c)
#get mean QALY/TC by trt group
mlm_e_em <- emmeans(mlm_e, ~ trt)
#see output
mlm_e_em
mlm_c_em <- emmeans(mlm_c, ~ trt)
#see output
mlm_c_em
#get difference as New - Old group
diff_new_vs_old <- list("New vs Old" = c(-1, 1))
#compute linear combination
mlm_em_delta_e <- contrast(mlm_e_em, diff_new_vs_old) 
mlm_em_delta_c <- contrast(mlm_c_em, diff_new_vs_old) 
#obtain results in terms of confidence intervals
confint(mlm_em_delta_e, level = 0.95)
confint(mlm_em_delta_c, level = 0.95)

#combine OLS/SUR with Two-stage bootstrapping to obtain distributions of statistics of interest 
#such as mean QALY/TC by arm and mean differences between groups
#get TSB results (200 iterations)
set.seed(2345)
tsboot_res <- tsboot_ec(data = dataset, QALYreg = QALY ~ trt + u0, 
                        TCreg = TC ~ trt + c0, cluster = "cluster", B=200)
summary(tsboot_res$QALY_boot$Delta_e)
summary(tsboot_res$QALY_boot$mu_e_ctr)
summary(tsboot_res$QALY_boot$mu_e_int)
summary(tsboot_res$TC_boot$Delta_c)
summary(tsboot_res$TC_boot$mu_c_ctr)
summary(tsboot_res$TC_boot$mu_c_int)
#compute percentile CIs (note: BCa CI not appropriate for tsboot results)
tsboot_ci_perc <- boot_ci(x = tsboot_res, method = "perc")
tsboot_ci_perc

#combine MLM with bootstrapping to obtain distributions of statistics of interest 
#such as mean QALY/TC by arm and mean differences between groups
#get bootstrap results (200 iterations)
set.seed(2345)
boot_res_mlm <- boot_ec_mlm(dataset, QALYreg = QALY ~ trt + u0,
                            TCreg = TC ~ trt + c0, 
                            QALYrandom = ~1|cluster, 
                            TCrandom = ~1|cluster, B=200)
#summarise output
summary(boot_res_mlm$QALY_boot$mu_e_ctr) #mean old arm
summary(boot_res_mlm$QALY_boot$mu_e_int) #mean new arm
summary(boot_res_mlm$QALY_boot$Delta_e) #mean difference
summary(boot_res_mlm$TC_boot$mu_c_ctr) #mean old arm
summary(boot_res_mlm$TC_boot$mu_c_int) #mean new arm
summary(boot_res_mlm$TC_boot$Delta_c) #mean difference
#compute percentile or BCa CIs
boot_ci_bca_mlm <- boot_ci_mlm(x = boot_res_mlm, method = "BCa")
#summarise output
boot_ci_bca_mlm


#VERSION 5: dealing with missing data

#here need to introduce in dataset some missing values in QALY/TC
#the following code does this assuming that missingness in outcomes only 
#depend on their baseline values (which are left fully observed) - MAR
#introduce MAR missingness in QALY/TC given u0/c0
library(boot)
#set parameter values for mechanisms
eta0_mar <- -9
eta_1_mar <- 10.5
p_mar_e <- inv.logit(eta0_mar + eta_1_mar*dataset$u0) #about 0.3 missing on avg
iota0_mar <- -4
iota_1_mar <- 3
p_mar_c <- inv.logit(iota0_mar + iota_1_mar*dataset$c0/100) #about 0.3 missing on avg
#generate missing data indicators (0=observed,1=missing)
set.seed(2345)
m_mar_e <- rbinom(dim(dataset)[1], p_mar_e, size = 1)
m_mar_c <- rbinom(dim(dataset)[1], p_mar_c, size = 1)
dataset.mar <- dataset
#introduce MAR missing QALY data and distinguish between a variable showing only the observed QALY values (QALY_obs) and another showing only the missing QALY values (QALY_mis)
dataset.mar$m_QALY <- m_mar_e
dataset.mar$m_TC <- m_mar_c
dataset.mar$QALY_obs <- ifelse(dataset.mar$m_QALY==0,dataset.mar$QALY,NA)
dataset.mar$TC_obs <- ifelse(dataset.mar$m_TC==0,dataset.mar$TC,NA)
#randomly shuffle rows of the MAR dataset
dataset.mar <- dataset.mar[sample(1:nrow(dataset.mar)), ]
#keep only relevant variables and rename them
dataset.mis <- dataset.mar[,c("QALY_obs","TC_obs","trt","u0","c0")]
names(dataset.mis) <- c("QALY","TC","trt","u0","c0")
#end simulation of missigness in outcomes

#impute data using MICE with pmm
library(mice)
#split data by treatment group (Old and New)
dataset.mis_Old <- dataset.mis[dataset.mis$trt=="old",]
dataset.mis_New <- dataset.mis[dataset.mis$trt=="new",]
#set up MICE inputs for old group
mice_Old <- mice(dataset.mis_Old, print = FALSE, method = 'pmm', maxit = 0)
pM_Old <- mice_Old$predictorMatrix #extract default predictor matrix
#customise pred matrix to your needs (row=variable to be imputed, column=variable used as predictor)
#in the matrix an entry of 1 means that the corresponding column variable is used as predictor for the corresponding row variable
pM_Old["QALY",c("TC","u0")] <- 1 #require that QALY always imputed using TC and u0
pM_Old["QALY",c("c0")] <- 0 #require that QALY is never imputed using c0
pM_Old["TC",c("u0")] <- 0 #require that TC is never imputed using u0
pM_Old[,c("trt")] <- 0 #require that any variable is never imputed using trt
meth_Old <- mice_Old$method #extract default imputation methods
#by default PMM assumed as imputation methods for numeric variables
#set up MICE inputs for new group
mice_New <- mice(dataset.mis_New, print = FALSE, method = 'pmm', maxit = 0)
pM_New <- mice_New$predictorMatrix #extract default predictor matrix
pM_New["QALY",c("TC","u0")] <- 1 #require that QALY always imputed using TC and u0
pM_New["QALY",c("c0")] <- 0 #require that QALY is never imputed using c0
pM_New["TC",c("u0")] <- 0 #require that TC is never imputed using u0
pM_New[,c("trt")] <- 0 #require that any variable is never imputed using trt
meth_New <- mice_New$method #extract default imputation methods
M <- 30 #number of imputations
set.seed(2345)
#implement MICE to old and new group
mice_Old_fit <- mice(dataset.mis_Old, predictorMatrix = pM_Old, method=meth_Old, m = M, print = FALSE)
mice_New_fit <- mice(dataset.mis_New, predictorMatrix = pM_New, method=meth_New, m = M, print = FALSE)
#combine the imputed datasets across groups
mice_fit <- rbind(mice_Old_fit, mice_New_fit)
#analyse imputed data with OLS for QALY and TC
lme_mice_data_e <- with(mice_fit, lm(QALY ~ trt + u0))
lme_mice_data_c <- with(mice_fit, lm(TC ~ trt + c0))
#get pooled mean QALY per group
em_lme_mu_data.mi_e <- emmeans(lme_mice_data_e, ~ trt)
em_lme_mu_data.mi_c <- emmeans(lme_mice_data_c, ~ trt)
#get pooled mean difference between groups
diff_new_vs_old <- list("New vs Old" = c(-1, 1))
confint(contrast(em_lme_mu_data.mi_e, diff_new_vs_old))
confint(contrast(em_lme_mu_data.mi_c, diff_new_vs_old))

#NOTE: I have only shown how to implement MICE to a dataset WITHOUT getting
#bootstrap results. This is because combining MI with bootstrap is not simple
#and things can be derived not in an unique way. I think the code above is enough
#for one time.


