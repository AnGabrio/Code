#this file contains the code for all functions used in the main ".qmd" file 
#you can automatically load these functions in your R workspace by simply typing 
#source("nameofthisfile.R") in your R console

#function and needed packages to generate bootstrap estimates for mean QALY/TC and incremental
#differences when running OLS or SUR models

library(data.table) #package to handle datasets more efficiently
library(bootstrap) #package to use bootstrap procedure 
library(rlang)
boot_ec <- function(data, B, QALYreg, TCreg, method = "OLS",
                    profile_QALY="default", profile_TC="default", trt_pos = 2){
  #the following lines are needed to make sure proper inputs are given
  if(!is.data.frame(data)){stop("data needs to be a data frame object")}
  if(!is.numeric(B)){stop("please provide number of bootstrap iterations")}
  if(B<=0 | !B%%1==0){stop("please provide number of bootstrap iterations")}
  if(!is_formula(QALYreg)){stop("please provide formula for QALY model")}
  if(!is_formula(TCreg)){stop("please provide formula for TC model")}
  if(!method %in% c("OLS","SUR")){stop("please provide valid method name")}
  if(!is.numeric(trt_pos) | length(trt_pos)!=1 | trt_pos<=0){stop("please provide valid trt indicator position in regressions")}
  n <- dim(data)[1] #original sample size
  #n covariates 
  nX_e <- dim(model.matrix(QALYreg, data))[2]
  nX_c <- dim(model.matrix(TCreg, data))[2]
  #extract name of trt indicator and outcomes from provided formula
  trt_name_e <- all.vars(QALYreg)[trt_pos]
  trt_name_c <- all.vars(TCreg)[trt_pos]
  if(trt_name_e != trt_name_c){stop("please provide same trt variable name and position in QALY and TC formuale")}
  QALY_name <- all.vars(QALYreg)[1]
  TC_name <- all.vars(TCreg)[1]
  #check if trt indicator is factor and store its levels
  if(is.factor(data[,trt_name_e])){
    trt_fact <- TRUE
    trt_lev <- levels(data[,trt_name_e])} else {
      trt_fact <- FALSE
      trt_lev <- unique(data[,trt_name_e])}
  if(length(trt_lev)!=2){stop("The function only allows comparison between two trt groups")}  
  #check that correct profile provided or set default
  if(profile_QALY != "default"){
    if(!is.vector(profile_QALY) | length(profile_QALY)!=nX_e){stop("provide valid profile for QALYreg")}}
  if(profile_TC != "default"){
    if(!is.vector(profile_TC) | length(profile_TC)!=nX_c){stop("provide valid profile for TCreg")}}
  #prepare empty objects to contain bootstrapped estimates
  data_ec_b_list <- list()
  coeff_e <- c()
  coeff_c <- c()
  em_e_ctr <- em_e_int <- c()
  em_c_ctr <- em_c_int <- c()
  dataset.dt <- data.table(data) #convert data into data.table object
  for(i in 1:B){
    #sample with replacement
    data_ec_b_list[[i]] <- dataset.dt[sample(.N, n, replace = T)]
    #fit model
    model_ec <- systemfit(list(QALYreg = QALYreg, TCreg = TCreg), 
                          method=method, data=data_ec_b_list[[i]])
    #extract covariate values
    X_e <- model.matrix(model_ec$eq[[1]])
    X_c <- model.matrix(model_ec$eq[[2]])
    #define QALYreg profile
    if(profile_QALY == "default"){
      profile_b_QALY <- apply(X_e, 2, mean, na.rm=T)
    } else {profile_b_QALY <- profile_QALY}
    profile_b_QALY_ctr <- profile_b_QALY_int <- profile_b_QALY
    profile_b_QALY_ctr[trt_pos] <- 0 #set profile for comparator
    profile_b_QALY_int[trt_pos] <- 1 #set profile for reference
    #define TCreg profile
    if(profile_TC == "default"){
      profile_b_TC <- apply(X_c, 2, mean, na.rm=T)
    } else {profile_b_TC <- profile_TC}
    profile_b_TC_ctr <- profile_b_TC_int <- profile_b_TC
    profile_b_TC_ctr[trt_pos] <- 0 #set profile for comparator
    profile_b_TC_int[trt_pos] <- 1 #set profile for reference
    #extract coefficient estimates from each model
    coeff_e[i] <- summary(model_ec$eq[[1]])$coefficients[trt_pos,"Estimate"]
    coeff_c[i] <- summary(model_ec$eq[[2]])$coefficients[trt_pos,"Estimate"]
    #compute linear combination of parameters
    em_e_ctr[i] <- t(profile_b_QALY_ctr) %*% summary(model_ec$eq[[1]])$coefficients[,"Estimate"] 
    em_e_int[i] <- t(profile_b_QALY_int) %*% summary(model_ec$eq[[1]])$coefficients[,"Estimate"] 
    em_c_ctr[i] <- t(profile_b_TC_ctr) %*% summary(model_ec$eq[[2]])$coefficients[,"Estimate"] 
    em_c_int[i] <- t(profile_b_TC_int) %*% summary(model_ec$eq[[2]])$coefficients[,"Estimate"] 
  }
  #create list objects to store all results 
  res_e_b_list <-list("Delta_e"=coeff_e,"mu_e_ctr"=em_e_ctr,"mu_e_int"=em_e_int)
  res_c_b_list <-list("Delta_c"=coeff_c,"mu_c_ctr"=em_c_ctr,"mu_c_int"=em_c_int)
  input_list <- list("data"=data, "method"=method, "trt_pos"=trt_pos, "QALYreg"=QALYreg,
                     "TCreg"=TCreg,"profile_QALY_ctr"=profile_b_QALY_ctr,
                     "profile_QALY_int"=profile_b_QALY_int,"profile_TC_ctr"=profile_b_TC_ctr,
                     "profile_TC_int"=profile_b_TC_int)
  #compute overall list and return it as output from the function
  res_ec_b_list <- list("QALY_boot"=res_e_b_list,"TC_boot"=res_c_b_list,"inputs"=input_list)
  class(res_ec_b_list) <- "bootCE"
  return(res_ec_b_list)
}



#functions and needed packages to generate bootstrap confidence intervals based on
#bootstrapped estimates obtained through the function boot_ec

#jackknife sampling function (used to compute BCa interval inside main boot_ci function)
jk_ec <- function(data,QALYreg,TCreg,trt_pos,profile_QALY_ctr,profile_QALY_int,
                  profile_TC_ctr,profile_TC_int,or_method){
  n <- dim(data)[1]
  #prepare objects to store results
  jk_delta_c_i <- jk_delta_e_i <- c()
  jk_mu0_c_i <- jk_mu0_e_i <- c()
  jk_mu1_c_i <- jk_mu1_e_i <- c()
  for(i in 1:n){
    #apply jackknife re-sampling
    data_i <- data[-i,]
    #obtain estimates of interest
    model_ec_i <- systemfit(list(QALYreg = QALYreg, TCreg = TCreg), 
                            method=or_method, data=data_i)
    jk_delta_e_i[i] <- summary(model_ec_i$eq[[1]])$coefficients[trt_pos,"Estimate"]
    jk_delta_c_i[i] <- summary(model_ec_i$eq[[2]])$coefficients[trt_pos,"Estimate"]
    jk_mu0_e_i[i] <- as.numeric(t(profile_QALY_ctr) %*% summary(model_ec_i$eq[[1]])$coefficients[,"Estimate"])
    jk_mu1_e_i[i] <- as.numeric(t(profile_QALY_int) %*% summary(model_ec_i$eq[[1]])$coefficients[,"Estimate"])
    jk_mu0_c_i[i] <- as.numeric(t(profile_TC_ctr) %*% summary(model_ec_i$eq[[2]])$coefficients[,"Estimate"])
    jk_mu1_c_i[i] <- as.numeric(t(profile_TC_int) %*% summary(model_ec_i$eq[[2]])$coefficients[,"Estimate"])
  }
  jk_est_i <- cbind.data.frame(jk_delta_e_i,jk_delta_c_i,jk_mu0_e_i,jk_mu1_e_i,jk_mu0_c_i,jk_mu1_c_i)
  names(jk_est_i) <- c("jk_Delta_e","jk_Delta_c","jk_mu_e_ctr","jk_mu_e_int","jk_mu_c_ctr","jk_mu_c_int")
  return(jk_est_i)
}

boot_ci <- function(x, method = "perc", confidence = 0.95){
  #the following lines are needed to make sure proper inputs are given
  if(!inherits(x, c("bootCE","tsbootCE"))) {stop("Only objects of class 'bootCE' or can 'tsbootCE' be used")}
  if(!method %in% c("perc","BCa")){stop("please provide valid method name")}
  if(!is.numeric(confidence)){stop("please provide valid confidence level")}
  if(confidence<=0 | confidence>=1){stop("please provide valid confidence level")}
  #extract information from inputs
  B <- length(x$QALY_boot$Delta_e)
  mu_e_ctr <- x$QALY_boot$mu_e_ctr
  mu_e_int <- x$QALY_boot$mu_e_int
  Delta_e <- x$QALY_boot$Delta_e
  mu_c_ctr <- x$TC_boot$mu_c_ctr
  mu_c_int <- x$TC_boot$mu_c_int
  Delta_c <- x$TC_boot$Delta_c
  alpha <- 1 - confidence
  #compute CI bounds according to method chosen
  if(method == "perc"){
    ci_mu_e_ctr <- quantile(mu_e_ctr, probs = c(alpha/2,(1-alpha/2)))
    ci_mu_e_int <- quantile(mu_e_int, probs = c(alpha/2,(1-alpha/2)))
    ci_Delta_e <- quantile(Delta_e, probs = c(alpha/2,(1-alpha/2)))
    ci_mu_c_ctr <- quantile(mu_c_ctr, probs = c(alpha/2,(1-alpha/2)))
    ci_mu_c_int <- quantile(mu_c_int, probs = c(alpha/2,(1-alpha/2)))
    ci_Delta_c <- quantile(Delta_c, probs = c(alpha/2,(1-alpha/2)))
  }
  if(method == "BCa"){
    or_method <- x$inputs$method
    data <- x$inputs$data
    trt_pos <- x$inputs$trt_pos
    QALYreg <- x$inputs$QALYreg
    TCreg <- x$inputs$TCreg
    profile_QALY_ctr <- x$inputs$profile_QALY_ctr
    profile_QALY_int <- x$inputs$profile_QALY_int
    profile_TC_ctr <- x$inputs$profile_TC_ctr
    profile_TC_int <- x$inputs$profile_TC_int
    #obtain avg BCa estimates based on original sample
    model_ec <- systemfit(list(QALYreg = QALYreg, TCreg = TCreg), 
                          method=or_method, data=data)
    avg_BCa_Delta_e <- summary(model_ec$eq[[1]])$coefficients[trt_pos,"Estimate"]
    avg_BCa_Delta_c <- summary(model_ec$eq[[2]])$coefficients[trt_pos,"Estimate"]
    avg_BCa_mu_e_ctr <- as.numeric(t(profile_QALY_ctr) %*% summary(model_ec$eq[[1]])$coefficients[,"Estimate"])
    avg_BCa_mu_e_int <- as.numeric(t(profile_QALY_int) %*% summary(model_ec$eq[[1]])$coefficients[,"Estimate"])
    avg_BCa_mu_c_ctr <- as.numeric(t(profile_TC_ctr) %*% summary(model_ec$eq[[2]])$coefficients[,"Estimate"])
    avg_BCa_mu_c_int <- as.numeric(t(profile_TC_int) %*% summary(model_ec$eq[[2]])$coefficients[,"Estimate"])
    #compute proportion of samples below avg estimates
    plower_Delta_e <- length(Delta_e[Delta_e<avg_BCa_Delta_e])/B
    plower_Delta_c <- length(Delta_c[Delta_c<avg_BCa_Delta_c])/B
    plower_mu_e_ctr <- length(mu_e_ctr[mu_e_ctr<avg_BCa_mu_e_ctr])/B
    plower_mu_e_int <- length(mu_e_int[mu_e_int<avg_BCa_mu_e_int])/B
    plower_mu_c_ctr <- length(mu_c_ctr[mu_c_ctr<avg_BCa_mu_c_ctr])/B
    plower_mu_c_int <- length(mu_c_int[mu_c_int<avg_BCa_mu_c_int])/B
    #compute bias-correction term
    z0_Delta_e <- qnorm(plower_Delta_e, mean = 0, sd = 1, lower.tail = TRUE)
    z0_Delta_c <- qnorm(plower_Delta_c, mean = 0, sd = 1, lower.tail = TRUE)
    z0_mu_e_ctr <- qnorm(plower_mu_e_ctr, mean = 0, sd = 1, lower.tail = TRUE)
    z0_mu_e_int <- qnorm(plower_mu_e_int, mean = 0, sd = 1, lower.tail = TRUE)
    z0_mu_c_ctr <- qnorm(plower_mu_c_ctr, mean = 0, sd = 1, lower.tail = TRUE)
    z0_mu_c_int <- qnorm(plower_mu_c_int, mean = 0, sd = 1, lower.tail = TRUE)
    #apply jackknife sampling functions to get jeckknife estimates
    jk_res <- jk_ec(data = data, QALYreg=QALYreg,TCreg=TCreg,trt_pos=trt_pos,
                    profile_QALY_ctr=profile_QALY_ctr,profile_QALY_int=profile_QALY_int,
                    profile_TC_ctr=profile_TC_ctr,profile_TC_int=profile_TC_int,or_method=or_method)
    #compute avg of jk estimates
    jk_res_avg <- apply(jk_res, 2, mean, na.rm=T) 
    #compute skewness correction term
    a_Delta_e <- sum((jk_res_avg["jk_Delta_e"] - jk_res[,"jk_Delta_e"])^3) / (6*(sum((jk_res_avg["jk_Delta_e"] - jk_res[,"jk_Delta_e"])^2))^(3/2))
    a_Delta_c <- sum((jk_res_avg["jk_Delta_c"] - jk_res[,"jk_Delta_c"])^3) / (6*(sum((jk_res_avg["jk_Delta_c"] - jk_res[,"jk_Delta_c"])^2))^(3/2))
    a_mu_e_ctr <- sum((jk_res_avg["jk_mu_e_ctr"] - jk_res[,"jk_mu_e_ctr"])^3) / (6*(sum((jk_res_avg["jk_mu_e_ctr"] - jk_res[,"jk_mu_e_ctr"])^2))^(3/2))
    a_mu_e_int <- sum((jk_res_avg["jk_mu_e_int"] - jk_res[,"jk_mu_e_int"])^3) / (6*(sum((jk_res_avg["jk_mu_e_int"] - jk_res[,"jk_mu_e_int"])^2))^(3/2))
    a_mu_c_ctr <- sum((jk_res_avg["jk_mu_c_ctr"] - jk_res[,"jk_mu_c_ctr"])^3) / (6*(sum((jk_res_avg["jk_mu_c_ctr"] - jk_res[,"jk_mu_c_ctr"])^2))^(3/2))
    a_mu_c_int <- sum((jk_res_avg["jk_mu_c_int"] - jk_res[,"jk_mu_c_int"])^3) / (6*(sum((jk_res_avg["jk_mu_c_int"] - jk_res[,"jk_mu_c_int"])^2))^(3/2))    
    #compute adjusted probs for getting desired confidence level
    z_alpha1 <- qnorm(alpha/2, mean = 0, sd = 1, lower.tail = TRUE)
    z_alpha2 <- qnorm(1-alpha/2, mean = 0, sd = 1, lower.tail = TRUE)
    ci_l_Delta_e <- pnorm(z0_Delta_e + ((z0_Delta_e+z_alpha1)/(1-a_Delta_e*(z0_Delta_e+z_alpha1))))
    ci_u_Delta_e <- pnorm(z0_Delta_e + ((z0_Delta_e+z_alpha2)/(1-a_Delta_e*(z0_Delta_e+z_alpha2))))
    ci_l_Delta_c <- pnorm(z0_Delta_c + ((z0_Delta_c+z_alpha1)/(1-a_Delta_c*(z0_Delta_c+z_alpha1))))
    ci_u_Delta_c <- pnorm(z0_Delta_c + ((z0_Delta_c+z_alpha2)/(1-a_Delta_c*(z0_Delta_c+z_alpha2))))
    ci_l_mu_e_ctr <- pnorm(z0_mu_e_ctr + ((z0_mu_e_ctr+z_alpha1)/(1-a_mu_e_ctr*(z0_mu_e_ctr+z_alpha1))))
    ci_u_mu_e_ctr <- pnorm(z0_mu_e_ctr + ((z0_mu_e_ctr+z_alpha2)/(1-a_mu_e_ctr*(z0_mu_e_ctr+z_alpha2))))
    ci_l_mu_e_int <- pnorm(z0_mu_e_int + ((z0_mu_e_int+z_alpha1)/(1-a_mu_e_int*(z0_mu_e_int+z_alpha1))))
    ci_u_mu_e_int <- pnorm(z0_mu_e_int + ((z0_mu_e_int+z_alpha2)/(1-a_mu_e_int*(z0_mu_e_int+z_alpha2))))
    ci_l_mu_c_ctr <- pnorm(z0_mu_c_ctr + ((z0_mu_c_ctr+z_alpha1)/(1-a_mu_c_ctr*(z0_mu_c_ctr+z_alpha1))))
    ci_u_mu_c_ctr <- pnorm(z0_mu_c_ctr + ((z0_mu_c_ctr+z_alpha2)/(1-a_mu_c_ctr*(z0_mu_c_ctr+z_alpha2))))
    ci_l_mu_c_int <- pnorm(z0_mu_c_int + ((z0_mu_c_int+z_alpha1)/(1-a_mu_c_int*(z0_mu_c_int+z_alpha1))))
    ci_u_mu_c_int <- pnorm(z0_mu_c_int + ((z0_mu_c_int+z_alpha2)/(1-a_mu_c_int*(z0_mu_c_int+z_alpha2))))
    #obtain quantiles on original scale
    ci_mu_e_ctr <- quantile(mu_e_ctr, probs = c(ci_l_mu_e_ctr,ci_u_mu_e_ctr))
    ci_mu_e_int <- quantile(mu_e_int, probs = c(ci_l_mu_e_int,ci_u_mu_e_int))
    ci_Delta_e <- quantile(Delta_e, probs = c(ci_l_Delta_e,ci_u_Delta_e))
    ci_mu_c_ctr <- quantile(mu_c_ctr, probs = c(ci_l_mu_c_ctr,ci_u_mu_c_ctr))
    ci_mu_c_int <- quantile(mu_c_int, probs = c(ci_l_mu_c_int,ci_u_mu_c_int))
    ci_Delta_c <- quantile(Delta_c, probs = c(ci_l_Delta_c,ci_u_Delta_c))
  }
  #organise and return results
  res_list <- list("Delta_e"=ci_Delta_e,"Delta_c"=ci_Delta_c,"mu_e_ctr"=ci_mu_e_ctr,"mu_e_int"=ci_mu_e_int,"mu_c_ctr"=ci_mu_c_ctr,"mu_c_int"=ci_mu_c_int)
  return(res_list)
}



#function and needed packages to generate bootstrap estimates for mean QALY/TC and incremental
#differences when running GLMs
#note: that compared to the function developed before for OLS/SUR models, this one
#does not require to provide profiles for the computation of mean estimates but instead relies on the emmeans function to obtain these estimates without manual computation. Because of this the function is less flexible and can only evaluate mean estimates assuming a single profile for both QALY and TC models (the one used by emmeans to compute these quantities). However, the function needs the user to specify different distributions and link functions for the QALY and TC models. 
library(data.table)
library(rlang)
library(mfx) 
library(MASS)
library(bootstrap)
library(emmeans)
boot_ec_glm <- function(data, B, QALYreg, TCreg, QALY_dist, TC_dist, 
                        QALY_link, TC_link, trt_pos = 2){
  #the following lines are needed to make sure proper inputs are given
  if(!is.data.frame(data)){stop("data needs to be a data frame object")}
  if(!is.numeric(B)){stop("please provide number of bootstrap iterations")}
  if(B<=0 | !B%%1==0){stop("please provide number of bootstrap iterations")}
  if(!is_formula(QALYreg)){stop("please provide formula for QALY model")}
  if(!is_formula(TCreg)){stop("please provide formula for TC model")}
  if(!QALY_dist %in% c("Beta","Binomial","NegBinomial","Gamma","InvGaussian","Poisson","Gaussian")){stop("please provide valid distribution name")}
  if(!TC_dist %in% c("Beta","Binomial","NegBinomial","Gamma","InvGaussian","Poisson","Gaussian")){stop("please provide valid distribution name")}
  if(!QALY_link %in% c("logit","probit","cauchit", "cloglog", "identity", "log", "sqrt", "1/mu^2", "inverse")){stop("please provide valid link function name")}
  if(!TC_link %in% c("logit","probit","cauchit", "cloglog", "identity", "log", "sqrt", "1/mu^2", "inverse")){stop("please provide valid link function name")}
  if(!is.numeric(trt_pos) | length(trt_pos)!=1 | trt_pos<=0){stop("please provide valid trt indicator position in regressions")}
  n <- dim(data)[1] #original sample size
  #n covariates 
  nX_e <- dim(model.matrix(QALYreg, data))[2]
  nX_c <- dim(model.matrix(TCreg, data))[2]
  #extract name of trt indicator and outcomes from provided formula
  trt_name_e <- all.vars(QALYreg)[trt_pos]
  trt_name_c <- all.vars(TCreg)[trt_pos]
  if(trt_name_e != trt_name_c){stop("please provide same trt variable name and position in QALY and TC formuale")}
  QALY_name <- all.vars(QALYreg)[1]
  TC_name <- all.vars(TCreg)[1]
  #check if trt indicator is factor and store its levels
  if(is.factor(data[,trt_name_e])){
    trt_fact <- TRUE
    trt_lev <- levels(data[,trt_name_e])} else {
      trt_fact <- FALSE
      trt_lev <- unique(data[,trt_name_e])}
  if(length(trt_lev)!=2){stop("The function only allows comparison between two trt groups")}  
  #prepare empty objects to contain bootstrapped estimates
  data_ec_b_list <- list()
  coeff_e <- c()
  coeff_c <- c()
  em_e_ctr <- em_e_int <- c()
  em_c_ctr <- em_c_int <- c()
  dataset.dt <- data.table(data) #convert data into data.table object
  for(i in 1:B){
    #sample with replacement
    data_ec_b_list[[i]] <- dataset.dt[sample(.N, n, replace = T)]
    #select and fit GLM based on distribution and link function (QALY)
    if(QALY_dist=="Beta"){
      glm_e <- betareg(QALYreg, data = data_ec_b_list[[i]], link = QALY_link)}
    if(QALY_dist=="NegBinomial"){
      glm_e <- glm.nb(QALYreg, data = data_ec_b_list[[i]], link = QALY_link)}
    if(QALY_dist=="Binomial"){
      glm_e <- glm(QALYreg, data = data_ec_b_list[[i]], family = binomial(link = QALY_link))}
    if(QALY_dist=="Gamma"){
      glm_e <- glm(QALYreg, data = data_ec_b_list[[i]], family = Gamma(link = QALY_link))}
    if(QALY_dist=="InvGaussian"){
      glm_e <- glm(QALYreg, data = data_ec_b_list[[i]], family = inverse.gaussian(link = QALY_link))}
    if(QALY_dist=="Poisson"){
      glm_e <- glm(QALYreg, data = data_ec_b_list[[i]], family = poisson(link = QALY_link))} 
    if(QALY_dist=="Gaussian"){
      glm_e <- glm(QALYreg, data = data_ec_b_list[[i]], family = gaussian(link = QALY_link))}
    #select and fit GLM based on distribution and link function (TC)
    if(TC_dist=="Beta"){
      glm_c <- betareg(TCreg, data = data_ec_b_list[[i]], link = TC_link)}
    if(TC_dist=="NegBinomial"){
      glm_c <- glm.nb(TCreg, data = data_ec_b_list[[i]], link = TC_link)}
    if(TC_dist=="Binomial"){
      glm_c <- glm(TCreg, data = data_ec_b_list[[i]], family = binomial(link = TC_link))}
    if(TC_dist=="Gamma"){
      glm_c <- glm(TCreg, data = data_ec_b_list[[i]], family = Gamma(link = TC_link))}
    if(TC_dist=="InvGaussian"){
      glm_c <- glm(TCreg, data = data_ec_b_list[[i]], family = inverse.gaussian(link = TC_link))}
    if(TC_dist=="Poisson"){
      glm_c <- glm(TCreg, data = data_ec_b_list[[i]], family = poisson(link = TC_link))} 
    if(TC_dist=="Gaussian"){
      glm_c <- glm(TCreg, data = data_ec_b_list[[i]], family = gaussian(link = TC_link))}
    #use emmeans function to get mean outcomes for each arm
    glm_e.em <- emmeans(glm_e, trt_name_e, type = "response", data = data_ec_b_list[[i]])
    glm_c.em <- emmeans(glm_c, trt_name_c, type = "response", data = data_ec_b_list[[i]])
    em_e_ctr[i] <- summary(glm_e.em)[1,2]
    em_e_int[i] <- summary(glm_e.em)[2,2]
    em_c_ctr[i] <- summary(glm_c.em)[1,2]
    em_c_int[i] <- summary(glm_c.em)[2,2]
    #specify and compute mean differences between groups
    coeff_e[i] <- em_e_int[i] - em_e_ctr[i]
    coeff_c[i] <- em_c_int[i] - em_c_ctr[i]
  }
  #create list objects to store all results 
  res_e_b_list <-list("Delta_e"=coeff_e,"mu_e_ctr"=em_e_ctr,"mu_e_int"=em_e_int)
  res_c_b_list <-list("Delta_c"=coeff_c,"mu_c_ctr"=em_c_ctr,"mu_c_int"=em_c_int)
  input_list <- list("data"=data, "trt_pos"=trt_pos, "QALYreg"=QALYreg,
                     "TCreg"=TCreg,"QALY_link"=QALY_link,"QALY_dist"=QALY_dist,
                     "TC_dist"=TC_dist,"TC_link"=TC_link)
  #compute overall list and return it as output from the function
  res_ec_b_list <- list("QALY_boot"=res_e_b_list,"TC_boot"=res_c_b_list,"inputs"=input_list)
  class(res_ec_b_list) <- "bootCE_glm"
  return(res_ec_b_list)
}


#functions and needed packages to generate bootstrap confidence intervals based on
#bootstrapped estimates obtained through the function boot_ec_glm

#jackknife sampling function (used to compute BCa interval inside main boot_ci function)
jk_ec_glm <- function(data,QALYreg,TCreg,trt_pos,QALY_dist,TC_dist,
                      QALY_link,TC_link){
  n <- dim(data)[1]
  #extract name of trt indicator from provided formula
  trt_name_e <- all.vars(QALYreg)[trt_pos]
  trt_name_c <- all.vars(TCreg)[trt_pos]
  #prepare objects to store results
  jk_delta_c_i <- jk_delta_e_i <- c()
  jk_mu0_c_i <- jk_mu0_e_i <- c()
  jk_mu1_c_i <- jk_mu1_e_i <- c()
  for(i in 1:n){
    #apply jackknife re-sampling
    data_i <- data[-i,]
    #obtain estimates of interest
    if(QALY_dist=="Beta"){
      glm_e <- betareg(QALYreg, data = data_i, link = QALY_link)}
    if(QALY_dist=="NegBinomial"){
      glm_e <- glm.nb(QALYreg, data = data_i, link = QALY_link)}
    if(QALY_dist=="Binomial"){
      glm_e <- glm(QALYreg, data = data_i, family = binomial(link = QALY_link))}
    if(QALY_dist=="Gamma"){
      glm_e <- glm(QALYreg, data = data_i, family = Gamma(link = QALY_link))}
    if(QALY_dist=="InvGaussian"){
      glm_e <- glm(QALYreg, data = data_i, family = inverse.gaussian(link = QALY_link))}
    if(QALY_dist=="Poisson"){
      glm_e <- glm(QALYreg, data = data_i, family = poisson(link = QALY_link))} 
    if(QALY_dist=="Gaussian"){
      glm_e <- glm(QALYreg, data = data_i, family = gaussian(link = QALY_link))}
    #select and fit GLM based on distribution and link function (TC)
    if(TC_dist=="Beta"){
      glm_c <- betareg(TCreg, data = data_i, link = TC_link)}
    if(TC_dist=="NegBinomial"){
      glm_c <- glm.nb(TCreg, data = data_i, link = TC_link)}
    if(TC_dist=="Binomial"){
      glm_c <- glm(TCreg, data = data_i, family = binomial(link = TC_link))}
    if(TC_dist=="Gamma"){
      glm_c <- glm(TCreg, data = data_i, family = Gamma(link = TC_link))}
    if(TC_dist=="InvGaussian"){
      glm_c <- glm(TCreg, data = data_i, family = inverse.gaussian(link = TC_link))}
    if(TC_dist=="Poisson"){
      glm_c <- glm(TCreg, data = data_i, family = poisson(link = TC_link))} 
    if(TC_dist=="Gaussian"){
      glm_c <- glm(TCreg, data = data_i, family = gaussian(link = TC_link))}    
    #use emmeans function to get mean outcomes for each arm
    glm_e.em <- emmeans(glm_e, trt_name_e, type = "response", data = data_i)
    glm_c.em <- emmeans(glm_c, trt_name_c, type = "response", data = data_i)
    jk_mu0_e_i[i] <- summary(glm_e.em)[1,2]
    jk_mu1_e_i[i] <- summary(glm_e.em)[2,2]
    jk_mu0_c_i[i] <- summary(glm_c.em)[1,2]
    jk_mu1_c_i[i] <- summary(glm_c.em)[2,2]
    #specify and compute mean differences between groups
    jk_delta_e_i[i] <- jk_mu1_e_i[i] - jk_mu0_e_i[i]
    jk_delta_c_i[i] <- jk_mu1_c_i[i] - jk_mu0_c_i[i]   
  }
  jk_est_i <- cbind.data.frame(jk_delta_e_i,jk_delta_c_i,jk_mu0_e_i,jk_mu1_e_i,jk_mu0_c_i,jk_mu1_c_i)
  names(jk_est_i) <- c("jk_Delta_e","jk_Delta_c","jk_mu_e_ctr","jk_mu_e_int","jk_mu_c_ctr","jk_mu_c_int")
  return(jk_est_i)
}

boot_ci_glm <- function(x, method = "perc", confidence = 0.95){
  #the following lines are needed to make sure proper inputs are given
  if(!inherits(x, c("bootCE_glm","tsbootCE_glm"))) {stop("Only objects of class 'bootCE_glm' or 'tsbootCE_glm' can be used")}
  if(!method %in% c("perc","BCa")){stop("please provide valid method name")}
  if(!is.numeric(confidence)){stop("please provide valid confidence level")}
  if(confidence<=0 | confidence>=1){stop("please provide valid confidence level")}
  #extract information from inputs
  B <- length(x$QALY_boot$Delta_e)
  mu_e_ctr <- x$QALY_boot$mu_e_ctr
  mu_e_int <- x$QALY_boot$mu_e_int
  Delta_e <- x$QALY_boot$Delta_e
  mu_c_ctr <- x$TC_boot$mu_c_ctr
  mu_c_int <- x$TC_boot$mu_c_int
  Delta_c <- x$TC_boot$Delta_c
  alpha <- 1 - confidence
  #compute CI bounds according to method chosen
  if(method == "perc"){
    ci_mu_e_ctr <- quantile(mu_e_ctr, probs = c(alpha/2,(1-alpha/2)))
    ci_mu_e_int <- quantile(mu_e_int, probs = c(alpha/2,(1-alpha/2)))
    ci_Delta_e <- quantile(Delta_e, probs = c(alpha/2,(1-alpha/2)))
    ci_mu_c_ctr <- quantile(mu_c_ctr, probs = c(alpha/2,(1-alpha/2)))
    ci_mu_c_int <- quantile(mu_c_int, probs = c(alpha/2,(1-alpha/2)))
    ci_Delta_c <- quantile(Delta_c, probs = c(alpha/2,(1-alpha/2)))
  }
  if(method == "BCa"){
    QALY_dist <- x$inputs$QALY_dist
    TC_dist <- x$inputs$TC_dist
    QALY_link <- x$inputs$QALY_link
    TC_link <- x$inputs$TC_link
    data <- x$inputs$data
    trt_pos <- x$inputs$trt_pos
    QALYreg <- x$inputs$QALYreg
    TCreg <- x$inputs$TCreg
    #obtain avg BCa estimates based on original sample
    if(QALY_dist=="Beta"){
      glm_e <- betareg(QALYreg, data = data, link = QALY_link)}
    if(QALY_dist=="NegBinomial"){
      glm_e <- glm.nb(QALYreg, data = data, link = QALY_link)}
    if(QALY_dist=="Binomial"){
      glm_e <- glm(QALYreg, data = data, family = binomial(link = QALY_link))}
    if(QALY_dist=="Gamma"){
      glm_e <- glm(QALYreg, data = data, family = Gamma(link = QALY_link))}
    if(QALY_dist=="InvGaussian"){
      glm_e <- glm(QALYreg, data = data, family = inverse.gaussian(link = QALY_link))}
    if(QALY_dist=="Poisson"){
      glm_e <- glm(QALYreg, data = data, family = poisson(link = QALY_link))} 
    if(QALY_dist=="Gaussian"){
      glm_e <- glm(QALYreg, data = data, family = gaussian(link = QALY_link))}
    #select and fit GLM based on distribution and link function (TC)
    if(TC_dist=="Beta"){
      glm_c <- betareg(TCreg, data = data, link = TC_link)}
    if(TC_dist=="NegBinomial"){
      glm_c <- glm.nb(TCreg, data = data, link = TC_link)}
    if(TC_dist=="Binomial"){
      glm_c <- glm(TCreg, data = data, family = binomial(link = TC_link))}
    if(TC_dist=="Gamma"){
      glm_c <- glm(TCreg, data = data, family = Gamma(link = TC_link))}
    if(TC_dist=="InvGaussian"){
      glm_c <- glm(TCreg, data = data, family = inverse.gaussian(link = TC_link))}
    if(TC_dist=="Poisson"){
      glm_c <- glm(TCreg, data = data, family = poisson(link = TC_link))} 
    if(TC_dist=="Gaussian"){
      glm_c <- glm(TCreg, data = data, family = gaussian(link = TC_link))}    
    #extract name of trt indicator from provided formula
    trt_name_e <- all.vars(QALYreg)[trt_pos]
    trt_name_c <- all.vars(TCreg)[trt_pos]
    #use emmeans function to get mean outcomes for each arm
    glm_e.em <- emmeans(glm_e, trt_name_e, type = "response", data = data)
    glm_c.em <- emmeans(glm_c, trt_name_c, type = "response", data = data)
    avg_BCa_mu_e_ctr <- summary(glm_e.em)[1,2]
    avg_BCa_mu_e_int <- summary(glm_e.em)[2,2]
    avg_BCa_mu_c_ctr <- summary(glm_c.em)[1,2]
    avg_BCa_mu_c_int <- summary(glm_c.em)[2,2]
    #specify and compute mean differences between groups
    avg_BCa_Delta_e <- avg_BCa_mu_e_int - avg_BCa_mu_e_ctr
    avg_BCa_Delta_c <- avg_BCa_mu_c_int - avg_BCa_mu_c_ctr    
    #compute proportion of samples below avg estimates
    plower_Delta_e <- length(Delta_e[Delta_e<avg_BCa_Delta_e])/B
    plower_Delta_c <- length(Delta_c[Delta_c<avg_BCa_Delta_c])/B
    plower_mu_e_ctr <- length(mu_e_ctr[mu_e_ctr<avg_BCa_mu_e_ctr])/B
    plower_mu_e_int <- length(mu_e_int[mu_e_int<avg_BCa_mu_e_int])/B
    plower_mu_c_ctr <- length(mu_c_ctr[mu_c_ctr<avg_BCa_mu_c_ctr])/B
    plower_mu_c_int <- length(mu_c_int[mu_c_int<avg_BCa_mu_c_int])/B
    #compute bias-correction term
    z0_Delta_e <- qnorm(plower_Delta_e, mean = 0, sd = 1, lower.tail = TRUE)
    z0_Delta_c <- qnorm(plower_Delta_c, mean = 0, sd = 1, lower.tail = TRUE)
    z0_mu_e_ctr <- qnorm(plower_mu_e_ctr, mean = 0, sd = 1, lower.tail = TRUE)
    z0_mu_e_int <- qnorm(plower_mu_e_int, mean = 0, sd = 1, lower.tail = TRUE)
    z0_mu_c_ctr <- qnorm(plower_mu_c_ctr, mean = 0, sd = 1, lower.tail = TRUE)
    z0_mu_c_int <- qnorm(plower_mu_c_int, mean = 0, sd = 1, lower.tail = TRUE)
    #apply jackknife sampling functions to get jeckknife estimates
    jk_res_glm <- jk_ec_glm(data = data, QALYreg=QALYreg,TCreg=TCreg,trt_pos=trt_pos,
                            QALY_dist=QALY_dist,TC_dist=TC_dist,
                            QALY_link=QALY_link,TC_link=TC_link)
    #compute avg of jk estimates
    jk_res_avg <- apply(jk_res_glm, 2, mean, na.rm = T) 
    #compute skewness correction term
    a_Delta_e <- sum((jk_res_avg["jk_Delta_e"] - jk_res_glm[,"jk_Delta_e"])^3) / (6*(sum((jk_res_avg["jk_Delta_e"] - jk_res_glm[,"jk_Delta_e"])^2))^(3/2))
    a_Delta_c <- sum((jk_res_avg["jk_Delta_c"] - jk_res_glm[,"jk_Delta_c"])^3) / (6*(sum((jk_res_avg["jk_Delta_c"] - jk_res_glm[,"jk_Delta_c"])^2))^(3/2))
    a_mu_e_ctr <- sum((jk_res_avg["jk_mu_e_ctr"] - jk_res_glm[,"jk_mu_e_ctr"])^3) / (6*(sum((jk_res_avg["jk_mu_e_ctr"] - jk_res_glm[,"jk_mu_e_ctr"])^2))^(3/2))
    a_mu_e_int <- sum((jk_res_avg["jk_mu_e_int"] - jk_res_glm[,"jk_mu_e_int"])^3) / (6*(sum((jk_res_avg["jk_mu_e_int"] - jk_res_glm[,"jk_mu_e_int"])^2))^(3/2))
    a_mu_c_ctr <- sum((jk_res_avg["jk_mu_c_ctr"] - jk_res_glm[,"jk_mu_c_ctr"])^3) / (6*(sum((jk_res_avg["jk_mu_c_ctr"] - jk_res_glm[,"jk_mu_c_ctr"])^2))^(3/2))
    a_mu_c_int <- sum((jk_res_avg["jk_mu_c_int"] - jk_res_glm[,"jk_mu_c_int"])^3) / (6*(sum((jk_res_avg["jk_mu_c_int"] - jk_res_glm[,"jk_mu_c_int"])^2))^(3/2))    
    #compute adjusted probs for getting desired confidence level
    z_alpha1 <- qnorm(alpha/2, mean = 0, sd = 1, lower.tail = TRUE)
    z_alpha2 <- qnorm(1-alpha/2, mean = 0, sd = 1, lower.tail = TRUE)
    ci_l_Delta_e <- pnorm(z0_Delta_e + ((z0_Delta_e+z_alpha1)/(1-a_Delta_e*(z0_Delta_e+z_alpha1))))
    ci_u_Delta_e <- pnorm(z0_Delta_e + ((z0_Delta_e+z_alpha2)/(1-a_Delta_e*(z0_Delta_e+z_alpha2))))
    ci_l_Delta_c <- pnorm(z0_Delta_c + ((z0_Delta_c+z_alpha1)/(1-a_Delta_c*(z0_Delta_c+z_alpha1))))
    ci_u_Delta_c <- pnorm(z0_Delta_c + ((z0_Delta_c+z_alpha2)/(1-a_Delta_c*(z0_Delta_c+z_alpha2))))
    ci_l_mu_e_ctr <- pnorm(z0_mu_e_ctr + ((z0_mu_e_ctr+z_alpha1)/(1-a_mu_e_ctr*(z0_mu_e_ctr+z_alpha1))))
    ci_u_mu_e_ctr <- pnorm(z0_mu_e_ctr + ((z0_mu_e_ctr+z_alpha2)/(1-a_mu_e_ctr*(z0_mu_e_ctr+z_alpha2))))
    ci_l_mu_e_int <- pnorm(z0_mu_e_int + ((z0_mu_e_int+z_alpha1)/(1-a_mu_e_int*(z0_mu_e_int+z_alpha1))))
    ci_u_mu_e_int <- pnorm(z0_mu_e_int + ((z0_mu_e_int+z_alpha2)/(1-a_mu_e_int*(z0_mu_e_int+z_alpha2))))
    ci_l_mu_c_ctr <- pnorm(z0_mu_c_ctr + ((z0_mu_c_ctr+z_alpha1)/(1-a_mu_c_ctr*(z0_mu_c_ctr+z_alpha1))))
    ci_u_mu_c_ctr <- pnorm(z0_mu_c_ctr + ((z0_mu_c_ctr+z_alpha2)/(1-a_mu_c_ctr*(z0_mu_c_ctr+z_alpha2))))
    ci_l_mu_c_int <- pnorm(z0_mu_c_int + ((z0_mu_c_int+z_alpha1)/(1-a_mu_c_int*(z0_mu_c_int+z_alpha1))))
    ci_u_mu_c_int <- pnorm(z0_mu_c_int + ((z0_mu_c_int+z_alpha2)/(1-a_mu_c_int*(z0_mu_c_int+z_alpha2))))
    #obtain quantiles on original scale
    ci_mu_e_ctr <- quantile(mu_e_ctr, probs = c(ci_l_mu_e_ctr,ci_u_mu_e_ctr))
    ci_mu_e_int <- quantile(mu_e_int, probs = c(ci_l_mu_e_int,ci_u_mu_e_int))
    ci_Delta_e <- quantile(Delta_e, probs = c(ci_l_Delta_e,ci_u_Delta_e))
    ci_mu_c_ctr <- quantile(mu_c_ctr, probs = c(ci_l_mu_c_ctr,ci_u_mu_c_ctr))
    ci_mu_c_int <- quantile(mu_c_int, probs = c(ci_l_mu_c_int,ci_u_mu_c_int))
    ci_Delta_c <- quantile(Delta_c, probs = c(ci_l_Delta_c,ci_u_Delta_c))
  }
  #organise and return results
  res_list <- list("Delta_e"=ci_Delta_e,"Delta_c"=ci_Delta_c,"mu_e_ctr"=ci_mu_e_ctr,"mu_e_int"=ci_mu_e_int,"mu_c_ctr"=ci_mu_c_ctr,"mu_c_int"=ci_mu_c_int)
  return(res_list)
}


#function and needed packages to generate bootstrap estimates for mean QALY/TC and incremental
#differences when running MLMs

library(data.table) #package to handle datasets more efficiently
library(bootstrap) #package to use bootstrap procedure 
library(rlang)
library(lme4)
library(nlme)
library(emmeans)
boot_ec_mlm <- function(data, B, QALYreg, TCreg, QALYrandom, TCrandom, trt_pos = 2){
  #the following lines are needed to make sure proper inputs are given
  if(!is.data.frame(data)){stop("data needs to be a data frame object")}
  if(!is.numeric(B)){stop("please provide number of bootstrap iterations")}
  if(B<=0 | !B%%1==0){stop("please provide number of bootstrap iterations")}
  if(!is_formula(QALYreg)){stop("please provide formula for QALY model")}
  if(!is_formula(TCreg)){stop("please provide formula for TC model")}
  if(!is_formula(QALYrandom)){stop("please provide formula for QALY random effects")}
  if(!is_formula(TCrandom)){stop("please provide formula for TC random effects")}
  if(!is.numeric(trt_pos) | length(trt_pos)!=1 | trt_pos<=0){stop("please provide valid trt indicator position in regressions")}
  n <- dim(data)[1] #original sample size
  #n covariates 
  nX_e <- dim(model.matrix(QALYreg, data))[2]
  nX_c <- dim(model.matrix(TCreg, data))[2]
  #extract name of trt indicator and outcomes from provided formula
  trt_name_e <- all.vars(QALYreg)[trt_pos]
  trt_name_c <- all.vars(TCreg)[trt_pos]
  if(trt_name_e != trt_name_c){stop("please provide same trt variable name and position in QALY and TC formuale")}
  QALY_name <- all.vars(QALYreg)[1]
  TC_name <- all.vars(TCreg)[1]
  #check if trt indicator is factor and store its levels
  if(is.factor(data[,trt_name_e])){
    trt_fact <- TRUE
    trt_lev <- levels(data[,trt_name_e])} else {
      trt_fact <- FALSE
      trt_lev <- unique(data[,trt_name_e])}
  if(length(trt_lev)!=2){stop("The function only allows comparison between two trt groups")}  
  #prepare empty objects to contain bootstrapped estimates
  data_ec_b_list <- list()
  coeff_e <- c()
  coeff_c <- c()
  em_e_ctr <- em_e_int <- c()
  em_c_ctr <- em_c_int <- c()
  dataset.dt <- data.table(data) #convert data into data.table object
  for(i in 1:B){
    #sample with replacement
    data_ec_b_list[[i]] <- dataset.dt[sample(.N, n, replace = T)]
    #fit models for QALY and TC
    mlm_e <- lme(QALYreg, random = QALYrandom, data = data_ec_b_list[[i]])
    mlm_c <- lme(TCreg, random = TCrandom, data = data_ec_b_list[[i]])
    #use emmeans function to get mean outcomes for each arm
    mlm_e.em <- emmeans(mlm_e, trt_name_e, type = "response", data = data_ec_b_list[[i]])
    mlm_c.em <- emmeans(mlm_c, trt_name_c, type = "response", data = data_ec_b_list[[i]])
    em_e_ctr[i] <- summary(mlm_e.em)[1,2]
    em_e_int[i] <- summary(mlm_e.em)[2,2]
    em_c_ctr[i] <- summary(mlm_c.em)[1,2]
    em_c_int[i] <- summary(mlm_c.em)[2,2]
    #specify and compute mean differences between groups
    coeff_e[i] <- em_e_int[i] - em_e_ctr[i]
    coeff_c[i] <- em_c_int[i] - em_c_ctr[i]
  }
  #create list objects to store all results 
  res_e_b_list <-list("Delta_e"=coeff_e,"mu_e_ctr"=em_e_ctr,"mu_e_int"=em_e_int)
  res_c_b_list <-list("Delta_c"=coeff_c,"mu_c_ctr"=em_c_ctr,"mu_c_int"=em_c_int)
  input_list <- list("data"=data, "trt_pos"=trt_pos, "QALYreg"=QALYreg,
                     "TCreg"=TCreg,"QALYrandom"=QALYrandom,"TCrandom"=TCrandom)
  #compute overall list and return it as output from the function
  res_ec_b_list <- list("QALY_boot"=res_e_b_list,"TC_boot"=res_c_b_list,"inputs"=input_list)
  class(res_ec_b_list) <- "bootCE_mlm"
  return(res_ec_b_list)
}

#functions and needed packages to generate bootstrap confidence intervals based on
#bootstrapped estimates obtained through the function boot_ec_mlm

#jackknife sampling function (used to compute BCa interval inside main boot_ci function)
jk_ec_mlm <- function(data,QALYreg,TCreg,trt_pos,QALYrandom,TCrandom){
  n <- dim(data)[1]
  #extract name of trt indicator from provided formula
  trt_name_e <- all.vars(QALYreg)[trt_pos]
  trt_name_c <- all.vars(TCreg)[trt_pos]
  #prepare objects to store results
  jk_delta_c_i <- jk_delta_e_i <- c()
  jk_mu0_c_i <- jk_mu0_e_i <- c()
  jk_mu1_c_i <- jk_mu1_e_i <- c()
  for(i in 1:n){
    #apply jackknife re-sampling
    data_i <- data[-i,]
    #fit models
    mlm_e <- lme(QALYreg, random = QALYrandom, data = data_i)
    mlm_c <- lme(TCreg, random = TCrandom, data = data_i)
    #use emmeans function to get mean outcomes for each arm
    mlm_e.em <- emmeans(mlm_e, trt_name_e, type = "response", data = data_i)
    mlm_c.em <- emmeans(mlm_c, trt_name_c, type = "response", data = data_i)
    jk_mu0_e_i[i] <- summary(mlm_e.em)[1,2]
    jk_mu1_e_i[i] <- summary(mlm_e.em)[2,2]
    jk_mu0_c_i[i] <- summary(mlm_c.em)[1,2]
    jk_mu1_c_i[i] <- summary(mlm_c.em)[2,2]
    #specify and compute mean differences between groups
    jk_delta_e_i[i] <- jk_mu1_e_i[i] - jk_mu0_e_i[i]
    jk_delta_c_i[i] <- jk_mu1_c_i[i] - jk_mu0_c_i[i]    
  }
  jk_est_i <- cbind.data.frame(jk_delta_e_i,jk_delta_c_i,jk_mu0_e_i,jk_mu1_e_i,jk_mu0_c_i,jk_mu1_c_i)
  names(jk_est_i) <- c("jk_Delta_e","jk_Delta_c","jk_mu_e_ctr","jk_mu_e_int","jk_mu_c_ctr","jk_mu_c_int")
  return(jk_est_i)
}

boot_ci_mlm <- function(x, method = "perc", confidence = 0.95){
  #the following lines are needed to make sure proper inputs are given
  if(!inherits(x, c("bootCE_mlm"))) {stop("Only objects of class 'bootCE_mlm' can be used")}
  if(!method %in% c("perc","BCa")){stop("please provide valid method name")}
  if(!is.numeric(confidence)){stop("please provide valid confidence level")}
  if(confidence<=0 | confidence>=1){stop("please provide valid confidence level")}
  #extract information from inputs
  B <- length(x$QALY_boot$Delta_e)
  mu_e_ctr <- x$QALY_boot$mu_e_ctr
  mu_e_int <- x$QALY_boot$mu_e_int
  Delta_e <- x$QALY_boot$Delta_e
  mu_c_ctr <- x$TC_boot$mu_c_ctr
  mu_c_int <- x$TC_boot$mu_c_int
  Delta_c <- x$TC_boot$Delta_c
  alpha <- 1 - confidence
  #compute CI bounds according to method chosen
  if(method == "perc"){
    ci_mu_e_ctr <- quantile(mu_e_ctr, probs = c(alpha/2,(1-alpha/2)))
    ci_mu_e_int <- quantile(mu_e_int, probs = c(alpha/2,(1-alpha/2)))
    ci_Delta_e <- quantile(Delta_e, probs = c(alpha/2,(1-alpha/2)))
    ci_mu_c_ctr <- quantile(mu_c_ctr, probs = c(alpha/2,(1-alpha/2)))
    ci_mu_c_int <- quantile(mu_c_int, probs = c(alpha/2,(1-alpha/2)))
    ci_Delta_c <- quantile(Delta_c, probs = c(alpha/2,(1-alpha/2)))
  }
  if(method == "BCa"){
    data <- x$inputs$data
    trt_pos <- x$inputs$trt_pos
    QALYreg <- x$inputs$QALYreg
    TCreg <- x$inputs$TCreg
    QALYrandom <- x$inputs$QALYrandom
    TCrandom <- x$inputs$TCrandom
    #obtain avg BCa estimates based on original sample
    mlm_e <- lme(QALYreg, random = QALYrandom, data = data)
    mlm_c <- lme(TCreg, random = TCrandom, data = data)
    #extract name of trt indicator from provided formula
    trt_name_e <- all.vars(QALYreg)[trt_pos]
    trt_name_c <- all.vars(TCreg)[trt_pos]
    #use emmeans function to get mean outcomes for each arm
    mlm_e.em <- emmeans(mlm_e, trt_name_e, type = "response", data = data)
    mlm_c.em <- emmeans(mlm_c, trt_name_c, type = "response", data = data)
    avg_BCa_mu_e_ctr <- summary(mlm_e.em)[1,2]
    avg_BCa_mu_e_int <- summary(mlm_e.em)[2,2]
    avg_BCa_mu_c_ctr <- summary(mlm_c.em)[1,2]
    avg_BCa_mu_c_int <- summary(mlm_c.em)[2,2]
    #specify and compute mean differences between groups
    avg_BCa_Delta_e <- avg_BCa_mu_e_int - avg_BCa_mu_e_ctr
    avg_BCa_Delta_c <- avg_BCa_mu_c_int - avg_BCa_mu_c_ctr     
    #compute proportion of samples below avg estimates
    plower_Delta_e <- length(Delta_e[Delta_e<avg_BCa_Delta_e])/B
    plower_Delta_c <- length(Delta_c[Delta_c<avg_BCa_Delta_c])/B
    plower_mu_e_ctr <- length(mu_e_ctr[mu_e_ctr<avg_BCa_mu_e_ctr])/B
    plower_mu_e_int <- length(mu_e_int[mu_e_int<avg_BCa_mu_e_int])/B
    plower_mu_c_ctr <- length(mu_c_ctr[mu_c_ctr<avg_BCa_mu_c_ctr])/B
    plower_mu_c_int <- length(mu_c_int[mu_c_int<avg_BCa_mu_c_int])/B
    #compute bias-correction term
    z0_Delta_e <- qnorm(plower_Delta_e, mean = 0, sd = 1, lower.tail = TRUE)
    z0_Delta_c <- qnorm(plower_Delta_c, mean = 0, sd = 1, lower.tail = TRUE)
    z0_mu_e_ctr <- qnorm(plower_mu_e_ctr, mean = 0, sd = 1, lower.tail = TRUE)
    z0_mu_e_int <- qnorm(plower_mu_e_int, mean = 0, sd = 1, lower.tail = TRUE)
    z0_mu_c_ctr <- qnorm(plower_mu_c_ctr, mean = 0, sd = 1, lower.tail = TRUE)
    z0_mu_c_int <- qnorm(plower_mu_c_int, mean = 0, sd = 1, lower.tail = TRUE)
    #apply jackknife sampling functions to get jeckknife estimates
    jk_res_mlm <- jk_ec_mlm(data = data, QALYreg=QALYreg,TCreg=TCreg,trt_pos=trt_pos,
                            QALYrandom=QALYrandom,TCrandom=TCrandom)
    #compute avg of jk estimates
    jk_res_avg <- apply(jk_res_mlm, 2, mean, na.rm=T) 
    #compute skewness correction term
    a_Delta_e <- sum((jk_res_avg["jk_Delta_e"] - jk_res_mlm[,"jk_Delta_e"])^3) / (6*(sum((jk_res_avg["jk_Delta_e"] - jk_res_mlm[,"jk_Delta_e"])^2))^(3/2))
    a_Delta_c <- sum((jk_res_avg["jk_Delta_c"] - jk_res_mlm[,"jk_Delta_c"])^3) / (6*(sum((jk_res_avg["jk_Delta_c"] - jk_res_mlm[,"jk_Delta_c"])^2))^(3/2))
    a_mu_e_ctr <- sum((jk_res_avg["jk_mu_e_ctr"] - jk_res_mlm[,"jk_mu_e_ctr"])^3) / (6*(sum((jk_res_avg["jk_mu_e_ctr"] - jk_res_mlm[,"jk_mu_e_ctr"])^2))^(3/2))
    a_mu_e_int <- sum((jk_res_avg["jk_mu_e_int"] - jk_res_mlm[,"jk_mu_e_int"])^3) / (6*(sum((jk_res_avg["jk_mu_e_int"] - jk_res_mlm[,"jk_mu_e_int"])^2))^(3/2))
    a_mu_c_ctr <- sum((jk_res_avg["jk_mu_c_ctr"] - jk_res_mlm[,"jk_mu_c_ctr"])^3) / (6*(sum((jk_res_avg["jk_mu_c_ctr"] - jk_res_mlm[,"jk_mu_c_ctr"])^2))^(3/2))
    a_mu_c_int <- sum((jk_res_avg["jk_mu_c_int"] - jk_res_mlm[,"jk_mu_c_int"])^3) / (6*(sum((jk_res_avg["jk_mu_c_int"] - jk_res_mlm[,"jk_mu_c_int"])^2))^(3/2))    
    #compute adjusted probs for getting desired confidence level
    z_alpha1 <- qnorm(alpha/2, mean = 0, sd = 1, lower.tail = TRUE)
    z_alpha2 <- qnorm(1-alpha/2, mean = 0, sd = 1, lower.tail = TRUE)
    ci_l_Delta_e <- pnorm(z0_Delta_e + ((z0_Delta_e+z_alpha1)/(1-a_Delta_e*(z0_Delta_e+z_alpha1))))
    ci_u_Delta_e <- pnorm(z0_Delta_e + ((z0_Delta_e+z_alpha2)/(1-a_Delta_e*(z0_Delta_e+z_alpha2))))
    ci_l_Delta_c <- pnorm(z0_Delta_c + ((z0_Delta_c+z_alpha1)/(1-a_Delta_c*(z0_Delta_c+z_alpha1))))
    ci_u_Delta_c <- pnorm(z0_Delta_c + ((z0_Delta_c+z_alpha2)/(1-a_Delta_c*(z0_Delta_c+z_alpha2))))
    ci_l_mu_e_ctr <- pnorm(z0_mu_e_ctr + ((z0_mu_e_ctr+z_alpha1)/(1-a_mu_e_ctr*(z0_mu_e_ctr+z_alpha1))))
    ci_u_mu_e_ctr <- pnorm(z0_mu_e_ctr + ((z0_mu_e_ctr+z_alpha2)/(1-a_mu_e_ctr*(z0_mu_e_ctr+z_alpha2))))
    ci_l_mu_e_int <- pnorm(z0_mu_e_int + ((z0_mu_e_int+z_alpha1)/(1-a_mu_e_int*(z0_mu_e_int+z_alpha1))))
    ci_u_mu_e_int <- pnorm(z0_mu_e_int + ((z0_mu_e_int+z_alpha2)/(1-a_mu_e_int*(z0_mu_e_int+z_alpha2))))
    ci_l_mu_c_ctr <- pnorm(z0_mu_c_ctr + ((z0_mu_c_ctr+z_alpha1)/(1-a_mu_c_ctr*(z0_mu_c_ctr+z_alpha1))))
    ci_u_mu_c_ctr <- pnorm(z0_mu_c_ctr + ((z0_mu_c_ctr+z_alpha2)/(1-a_mu_c_ctr*(z0_mu_c_ctr+z_alpha2))))
    ci_l_mu_c_int <- pnorm(z0_mu_c_int + ((z0_mu_c_int+z_alpha1)/(1-a_mu_c_int*(z0_mu_c_int+z_alpha1))))
    ci_u_mu_c_int <- pnorm(z0_mu_c_int + ((z0_mu_c_int+z_alpha2)/(1-a_mu_c_int*(z0_mu_c_int+z_alpha2))))
    #obtain quantiles on original scale
    ci_mu_e_ctr <- quantile(mu_e_ctr, probs = c(ci_l_mu_e_ctr,ci_u_mu_e_ctr))
    ci_mu_e_int <- quantile(mu_e_int, probs = c(ci_l_mu_e_int,ci_u_mu_e_int))
    ci_Delta_e <- quantile(Delta_e, probs = c(ci_l_Delta_e,ci_u_Delta_e))
    ci_mu_c_ctr <- quantile(mu_c_ctr, probs = c(ci_l_mu_c_ctr,ci_u_mu_c_ctr))
    ci_mu_c_int <- quantile(mu_c_int, probs = c(ci_l_mu_c_int,ci_u_mu_c_int))
    ci_Delta_c <- quantile(Delta_c, probs = c(ci_l_Delta_c,ci_u_Delta_c))
  }
  #organise and return results
  res_list <- list("Delta_e"=ci_Delta_e,"Delta_c"=ci_Delta_c,"mu_e_ctr"=ci_mu_e_ctr,"mu_e_int"=ci_mu_e_int,"mu_c_ctr"=ci_mu_c_ctr,"mu_c_int"=ci_mu_c_int)
  return(res_list)
}

#function and needed packages to generate TS bootstrap estimates for mean QALY/TC and incremental
#differences when running OLS or SUR models

library(data.table) #package to handle datasets more efficiently
library(bootstrap) #package to use bootstrap procedure 
library(rlang)
tsboot_ec <- function(data, B, QALYreg, TCreg, method = "OLS", cluster, unbalclus="donner",
                      profile_QALY="default", profile_TC="default", trt_pos = 2){
  #the following lines are needed to make sure proper inputs are given
  if(!is.data.frame(data)){stop("data needs to be a data frame object")}
  if(!is.numeric(B)){stop("please provide number of bootstrap iterations")}
  if(B<=0 | !B%%1==0){stop("please provide number of bootstrap iterations")}
  if(!is_formula(QALYreg)){stop("please provide formula for QALY model")}
  if(!is_formula(TCreg)){stop("please provide formula for TC model")}
  if(!method %in% c("OLS","SUR")){stop("please provide valid method name")}
  if(!is.numeric(trt_pos) | length(trt_pos)!=1 | trt_pos<=0){stop("please provide valid trt indicator position in regressions")}
  if(!is.character(cluster)){stop("please provide valid cluster variable name")}
  if(!cluster %in% names(data)){stop("please provide valid cluster variable name")}
  if(!unbalclus %in% c("donner","median","mean")){stop("please provide valid method to compute avg cluster size when standardising")}
  #convert cluster as factor and then numeric 
  data[,cluster] <- as.numeric(as.factor(data[,cluster]))
  #check that cluster variable is integer 
  if(!all(data[,cluster] - floor(data[,cluster]) == 0)){stop("cluster values should be integers")}
  n_size <- dim(data)[1] #original sample size
  #n covariates 
  nX_e <- dim(model.matrix(QALYreg, data))[2]
  nX_c <- dim(model.matrix(TCreg, data))[2]
  #check that correct profile provided or set default
  if(profile_QALY != "default"){
    if(!is.vector(profile_QALY) | length(profile_QALY)!=nX_e){stop("provide valid profile for QALYreg")}}
  if(profile_TC != "default"){
    if(!is.vector(profile_TC) | length(profile_TC)!=nX_c){stop("provide valid profile for TCreg")}}
  #extract name of trt indicator and outcomes from provided formula
  trt_name_e <- all.vars(QALYreg)[trt_pos]
  trt_name_c <- all.vars(TCreg)[trt_pos]
  if(trt_name_e != trt_name_c){stop("please provide same trt variable name and position in QALY and TC formuale")}
  QALY_name <- all.vars(QALYreg)[1]
  TC_name <- all.vars(TCreg)[1]
  #check if trt indicator is factor and store its levels
  if(is.factor(data[,trt_name_e])){
    trt_fact <- TRUE
    trt_lev <- levels(data[,trt_name_e])} else {
      trt_fact <- FALSE
      trt_lev <- unique(data[,trt_name_e])}
  if(length(trt_lev)!=2){stop("The function only allows comparison between two trt groups")}
  #prepare empty objects to contain bootstrapped estimates
  coeff_e <- c()
  coeff_c <- c()
  em_e_ctr <- em_e_int <- c()
  em_c_ctr <- em_c_int <- c()
  for(i in 1:B){
    count <- 0 #set count for while loop across strata
    n.strata <- length(unique(data[,trt_name_e]))
    shrunk.data <- c()
    while (count<n.strata){
      count <- count+1
      data1 <- data.frame(data[data[,trt_name_e]==unique(data[,trt_name_e])[count],])
      clus.size <- table(data1[,cluster])
      cost.x <- tapply(data1[,TC_name],data1[,cluster],mean) # calc cluster means
      qaly.x <- tapply(data1[,QALY_name],data1[,cluster],mean) # calc cluster means
      # STANDARDIZE Z: calc b for standardiwing z
      a <- length(unique(data1[,cluster]))
      if (var(clus.size)==0){
        b <- unique(clus.size)
      } else {
        if (unbalclus=="donner"){
          ifelse(warning,print("'average' clus size = Donner"),NA)
          n <- sum(clus.size)
          b <- (n-(sum(clus.size^2)/n))/(a-1)
        } else if (unbalclus=="median"){
          ifelse(warning,print("'average' clus size = median"),NA)
          b <- median(clus.size)
        } else if (unbalclus=="mean"){
          ifelse(warning,print("'average' clus size = mean"),NA)
          b <- mean(clus.size)
        } else {}
      } # End of 'else'
      # standardise z using cluster means (dfm = deviation from cluster mean)
      cost.dfm <- data1[,TC_name]-rep(cost.x,times=clus.size)
      qaly.dfm <- data1[,QALY_name]-rep(qaly.x,times=clus.size)
      cost.z <- (cost.dfm)/sqrt(1-1/b)
      qaly.z <- (qaly.dfm)/sqrt(1-1/b)
      # SHRINKAGE: calc c for shrinking x
      cost.ssw <- sum(cost.dfm^2)
      qaly.ssw <- sum(qaly.dfm^2)
      cost.ssb <- sum((cost.x-mean(cost.x))^2)
      qaly.ssb <- sum((qaly.x-mean(qaly.x))^2)
      cost.rhs <- a/(a-1) - cost.ssw/(b*(b-1)*cost.ssb)
      qaly.rhs <- a/(a-1) - qaly.ssw/(b*(b-1)*qaly.ssb)
      ifelse(cost.rhs<0, cost.c<-1, cost.c<-1-sqrt(cost.rhs))
      ifelse(qaly.rhs<0, qaly.c<-1, qaly.c<-1-sqrt(qaly.rhs))
      ## re-calc x
      cost.x <- cost.c*mean(data1[,TC_name]) + (1-cost.c)*cost.x
      qaly.x <- qaly.c*mean(data1[,QALY_name]) + (1-qaly.c)*qaly.x
      # TWO-STAGE SAMPLING & RE-CONSTRUCT OBS WITH SHRUNKEN MEANS AND STANDARDISED RESIDUALS
      # gen random clus (order) id with replacement
      sampled.x.cid <- sample(1:length(unique(data1[,cluster])),replace=T)
      sampled.z.iid <- sample(1:length(cost.z),sum(clus.size[sampled.x.cid]),replace=T) # chosen ind ids     for varying stratum sizes
      sampled.cost <- rep(cost.x[sampled.x.cid],times=clus.size[sampled.x.cid])+cost.z[sampled.z.iid]
      sampled.qaly <- rep(qaly.x[sampled.x.cid],times=clus.size[sampled.x.cid])+qaly.z[sampled.z.iid]
      # bind data from multiple strata together
      shrunk.data <- as.data.frame(rbind(shrunk.data,cbind(sampled.cost,sampled.qaly,
                                                           rep(unique(data1[,cluster])[sampled.x.cid],times=clus.size[sampled.x.cid]),
                                                           rep(unique(data[,trt_name_e])[count],times=sum(clus.size[sampled.x.cid])))))
    } # end of while
    #rename variables
    colnames(shrunk.data) <- c(TC_name,QALY_name,cluster,trt_name_e)
    #copy trt levels if factor
    if(is_true(trt_fact)){ 
      shrunk.data[,trt_name_e] <- factor(shrunk.data[,trt_name_e], levels=sort(unique(shrunk.data[,trt_name_e])), labels = trt_lev)}
    #create a dt object
    dataset_tsb.dt <- data.table(shrunk.data)
    #fit model
    model_ec <- systemfit(list(QALYreg = QALYreg, TCreg = TCreg), 
                          method=method, data=dataset_tsb.dt)
    #extract covariate values
    X_e <- model.matrix(model_ec$eq[[1]])
    X_c <- model.matrix(model_ec$eq[[2]])
    #define QALYreg profile
    if(profile_QALY == "default"){
      profile_b_QALY <- apply(X_e, 2, mean, na.rm=T)
    } else {profile_b_QALY <- profile_QALY}
    profile_b_QALY_ctr <- profile_b_QALY_int <- profile_b_QALY
    profile_b_QALY_ctr[trt_pos] <- 0 #set profile for comparator
    profile_b_QALY_int[trt_pos] <- 1 #set profile for reference
    #define TCreg profile
    if(profile_TC == "default"){
      profile_b_TC <- apply(X_c, 2, mean, na.rm=T)
    } else {profile_b_TC <- profile_TC}
    profile_b_TC_ctr <- profile_b_TC_int <- profile_b_TC
    profile_b_TC_ctr[trt_pos] <- 0 #set profile for comparator
    profile_b_TC_int[trt_pos] <- 1 #set profile for reference
    #extract coefficient estimates from each model
    coeff_e[i] <- summary(model_ec$eq[[1]])$coefficients[trt_pos,"Estimate"]
    coeff_c[i] <- summary(model_ec$eq[[2]])$coefficients[trt_pos,"Estimate"]
    #compute linear combination of parameters
    em_e_ctr[i] <- t(profile_b_QALY_ctr) %*% summary(model_ec$eq[[1]])$coefficients[,"Estimate"] 
    em_e_int[i] <- t(profile_b_QALY_int) %*% summary(model_ec$eq[[1]])$coefficients[,"Estimate"] 
    em_c_ctr[i] <- t(profile_b_TC_ctr) %*% summary(model_ec$eq[[2]])$coefficients[,"Estimate"] 
    em_c_int[i] <- t(profile_b_TC_int) %*% summary(model_ec$eq[[2]])$coefficients[,"Estimate"] 
  }
  #create list objects to store all results 
  res_e_tsb_list <-list("Delta_e"=coeff_e,"mu_e_ctr"=em_e_ctr,"mu_e_int"=em_e_int)
  res_c_tsb_list <-list("Delta_c"=coeff_c,"mu_c_ctr"=em_c_ctr,"mu_c_int"=em_c_int)
  input_list <- list("data"=data, "method"=method, "trt_pos"=trt_pos, "QALYreg"=QALYreg,
                     "TCreg"=TCreg,"profile_QALY_ctr"=profile_b_QALY_ctr,
                     "profile_QALY_int"=profile_b_QALY_int,"profile_TC_ctr"=profile_b_TC_ctr,
                     "profile_TC_int"=profile_b_TC_int, "cluster"=cluster, "unbalclus"=unbalclus)
  #compute overall list and return it as output from the function
  res_ec_tsb_list <- list("QALY_boot"=res_e_tsb_list,"TC_boot"=res_c_tsb_list,"inputs"=input_list)
  class(res_ec_tsb_list) <- "tsbootCE"
  return(res_ec_tsb_list)
}

#function and needed packages to generate TS bootstrap estimates for mean QALY/TC and incremental
#differences when running GLMs
#note: that compared to the function developed before for OLS/SUR models, this one
#does not require to provide profiles for the computation of mean estimates but instead relies on the emmeans function to obtain these estimates without manual computation. Because of this the function is less flexible and can only evaluate mean estimates assuming a single profile for both QALY and TC models (the one used by emmeans to compute these quantities). However, the function needs the user to specify different distributions and link functions for the QALY and TC models. 
library(data.table)
library(rlang)
library(mfx) 
library(MASS)
library(bootstrap)
library(emmeans)
tsboot_ec_glm <- function(data, B, QALYreg, TCreg, QALY_dist, TC_dist, 
                          QALY_link, TC_link, cluster, unbalclus="donner", trt_pos = 2){
  #the following lines are needed to make sure proper inputs are given
  if(!is.data.frame(data)){stop("data needs to be a data frame object")}
  if(!is.numeric(B)){stop("please provide number of bootstrap iterations")}
  if(B<=0 | !B%%1==0){stop("please provide number of bootstrap iterations")}
  if(!is_formula(QALYreg)){stop("please provide formula for QALY model")}
  if(!is_formula(TCreg)){stop("please provide formula for TC model")}
  if(!QALY_dist %in% c("Beta","Binomial","NegBinomial","Gamma","InvGaussian","Poisson","Gaussian")){stop("please provide valid distribution name")}
  if(!TC_dist %in% c("Beta","Binomial","NegBinomial","Gamma","InvGaussian","Poisson","Gaussian")){stop("please provide valid distribution name")}
  if(!QALY_link %in% c("logit","probit","cauchit", "cloglog", "identity", "log", "sqrt", "1/mu^2", "inverse")){stop("please provide valid link function name")}
  if(!TC_link %in% c("logit","probit","cauchit", "cloglog", "identity", "log", "sqrt", "1/mu^2", "inverse")){stop("please provide valid link function name")}
  if(!is.numeric(trt_pos) | length(trt_pos)!=1 | trt_pos<=0){stop("please provide valid trt indicator position in regressions")}
  if(!is.character(cluster)){stop("please provide valid cluster variable name")}
  if(!cluster %in% names(data)){stop("please provide valid cluster variable name")}
  if(!unbalclus %in% c("donner","median","mean")){stop("please provide valid method to compute avg cluster size when standardising")}
  #convert cluster as factor and then numeric 
  data[,cluster] <- as.numeric(as.factor(data[,cluster]))
  #check that cluster variable is integer 
  if(!all(data[,cluster] - floor(data[,cluster]) == 0)){stop("cluster values should be integers")}
  n_size <- dim(data)[1] #original sample size
  #n covariates 
  nX_e <- dim(model.matrix(QALYreg, data))[2]
  nX_c <- dim(model.matrix(TCreg, data))[2]
  #extract name of trt indicator and outcomes from provided formula
  trt_name_e <- all.vars(QALYreg)[trt_pos]
  trt_name_c <- all.vars(TCreg)[trt_pos]
  if(trt_name_e != trt_name_c){stop("please provide same trt variable name and position in QALY and TC formuale")}
  QALY_name <- all.vars(QALYreg)[1]
  TC_name <- all.vars(TCreg)[1]
  #check if trt indicator is factor and store its levels
  if(is.factor(data[,trt_name_e])){
    trt_fact <- TRUE
    trt_lev <- levels(data[,trt_name_e])} else {
      trt_fact <- FALSE
      trt_lev <- unique(data[,trt_name_e])}
  if(length(trt_lev)!=2){stop("The function only allows comparison between two trt groups")}
  #prepare empty objects to contain bootstrapped estimates
  coeff_e <- c()
  coeff_c <- c()
  em_e_ctr <- em_e_int <- c()
  em_c_ctr <- em_c_int <- c()
  for(i in 1:B){
    count <- 0 #set count for while loop across strata
    n.strata <- length(unique(data[,trt_name_e]))
    shrunk.data <- c()
    while (count<n.strata){
      count <- count+1
      data1 <- data.frame(data[data[,trt_name_e]==unique(data[,trt_name_e])[count],])
      clus.size <- table(data1[,cluster])
      cost.x <- tapply(data1[,TC_name],data1[,cluster],mean) # calc cluster means
      qaly.x <- tapply(data1[,QALY_name],data1[,cluster],mean) # calc cluster means
      # STANDARDIZE Z: calc b for standardiwing z
      a <- length(unique(data1[,cluster]))
      if (var(clus.size)==0){
        b <- unique(clus.size)
      } else {
        if (unbalclus=="donner"){
          ifelse(warning,print("'average' clus size = Donner"),NA)
          n <- sum(clus.size)
          b <- (n-(sum(clus.size^2)/n))/(a-1)
        } else if (unbalclus=="median"){
          ifelse(warning,print("'average' clus size = median"),NA)
          b <- median(clus.size)
        } else if (unbalclus=="mean"){
          ifelse(warning,print("'average' clus size = mean"),NA)
          b <- mean(clus.size)
        } else {}
      } # End of 'else'
      # standardise z using cluster means (dfm = deviation from cluster mean)
      cost.dfm <- data1[,TC_name]-rep(cost.x,times=clus.size)
      qaly.dfm <- data1[,QALY_name]-rep(qaly.x,times=clus.size)
      cost.z <- (cost.dfm)/sqrt(1-1/b)
      qaly.z <- (qaly.dfm)/sqrt(1-1/b)
      # SHRINKAGE: calc c for shrinking x
      cost.ssw <- sum(cost.dfm^2)
      qaly.ssw <- sum(qaly.dfm^2)
      cost.ssb <- sum((cost.x-mean(cost.x))^2)
      qaly.ssb <- sum((qaly.x-mean(qaly.x))^2)
      cost.rhs <- a/(a-1) - cost.ssw/(b*(b-1)*cost.ssb)
      qaly.rhs <- a/(a-1) - qaly.ssw/(b*(b-1)*qaly.ssb)
      ifelse(cost.rhs<0, cost.c<-1, cost.c<-1-sqrt(cost.rhs))
      ifelse(qaly.rhs<0, qaly.c<-1, qaly.c<-1-sqrt(qaly.rhs))
      ## re-calc x
      cost.x <- cost.c*mean(data1[,TC_name]) + (1-cost.c)*cost.x
      qaly.x <- qaly.c*mean(data1[,QALY_name]) + (1-qaly.c)*qaly.x
      # TWO-STAGE SAMPLING & RE-CONSTRUCT OBS WITH SHRUNKEN MEANS AND STANDARDISED RESIDUALS
      # gen random clus (order) id with replacement
      sampled.x.cid <- sample(1:length(unique(data1[,cluster])),replace=T)
      sampled.z.iid <- sample(1:length(cost.z),sum(clus.size[sampled.x.cid]),replace=T) # chosen ind ids     for varying stratum sizes
      sampled.cost <- rep(cost.x[sampled.x.cid],times=clus.size[sampled.x.cid])+cost.z[sampled.z.iid]
      sampled.qaly <- rep(qaly.x[sampled.x.cid],times=clus.size[sampled.x.cid])+qaly.z[sampled.z.iid]
      # bind data from multiple strata together
      shrunk.data <- as.data.frame(rbind(shrunk.data,cbind(sampled.cost,sampled.qaly,
                                                           rep(unique(data1[,cluster])[sampled.x.cid],times=clus.size[sampled.x.cid]),
                                                           rep(unique(data[,trt_name_e])[count],times=sum(clus.size[sampled.x.cid])))))
    } # end of while
    #rename variables
    colnames(shrunk.data) <- c(TC_name,QALY_name,cluster,trt_name_e)
    #copy trt levels if factor
    if(is_true(trt_fact)){ 
      shrunk.data[,trt_name_e] <- factor(shrunk.data[,trt_name_e], levels=sort(unique(shrunk.data[,trt_name_e])), labels = trt_lev)}
    #create a dt object
    dataset_tsb.dt <- data.table(shrunk.data)
    #select and fit GLM based on distribution and link function (QALY)
    if(QALY_dist=="Beta"){
      glm_e <- betareg(QALYreg, data = dataset_tsb.dt, link = QALY_link)}
    if(QALY_dist=="NegBinomial"){
      glm_e <- glm.nb(QALYreg, data = dataset_tsb.dt, link = QALY_link)}
    if(QALY_dist=="Binomial"){
      glm_e <- glm(QALYreg, data = dataset_tsb.dt, family = binomial(link = QALY_link))}
    if(QALY_dist=="Gamma"){
      glm_e <- glm(QALYreg, data = dataset_tsb.dt, family = Gamma(link = QALY_link))}
    if(QALY_dist=="InvGaussian"){
      glm_e <- glm(QALYreg, data = dataset_tsb.dt, family = inverse.gaussian(link = QALY_link))}
    if(QALY_dist=="Poisson"){
      glm_e <- glm(QALYreg, data = dataset_tsb.dt, family = poisson(link = QALY_link))} 
    if(QALY_dist=="Gaussian"){
      glm_e <- glm(QALYreg, data = dataset_tsb.dt, family = gaussian(link = QALY_link))}
    #select and fit GLM based on distribution and link function (TC)
    if(TC_dist=="Beta"){
      glm_c <- betareg(TCreg, data = dataset_tsb.dt, link = TC_link)}
    if(TC_dist=="NegBinomial"){
      glm_c <- glm.nb(TCreg, data = dataset_tsb.dt, link = TC_link)}
    if(TC_dist=="Binomial"){
      glm_c <- glm(TCreg, data = dataset_tsb.dt, family = binomial(link = TC_link))}
    if(TC_dist=="Gamma"){
      glm_c <- glm(TCreg, data = dataset_tsb.dt, family = Gamma(link = TC_link))}
    if(TC_dist=="InvGaussian"){
      glm_c <- glm(TCreg, data = dataset_tsb.dt, family = inverse.gaussian(link = TC_link))}
    if(TC_dist=="Poisson"){
      glm_c <- glm(TCreg, data = dataset_tsb.dt, family = poisson(link = TC_link))} 
    if(TC_dist=="Gaussian"){
      glm_c <- glm(TCreg, data = dataset_tsb.dt, family = gaussian(link = TC_link))}
    #use emmeans function to get mean outcomes for each arm
    glm_e.em <- emmeans(glm_e, trt_name_e, type = "response", data = dataset_tsb.dt)
    glm_c.em <- emmeans(glm_c, trt_name_c, type = "response", data = dataset_tsb.dt)
    em_e_ctr[i] <- summary(glm_e.em)[1,2]
    em_e_int[i] <- summary(glm_e.em)[2,2]
    em_c_ctr[i] <- summary(glm_c.em)[1,2]
    em_c_int[i] <- summary(glm_c.em)[2,2]
    #specify and compute mean differences between groups
    coeff_e[i] <- em_e_int[i] - em_e_ctr[i]
    coeff_c[i] <- em_c_int[i] - em_c_ctr[i]
  }
  #create list objects to store all results 
  res_e_tsb_list <-list("Delta_e"=coeff_e,"mu_e_ctr"=em_e_ctr,"mu_e_int"=em_e_int)
  res_c_tsb_list <-list("Delta_c"=coeff_c,"mu_c_ctr"=em_c_ctr,"mu_c_int"=em_c_int)
  input_list <- list("data"=data, "trt_pos"=trt_pos, "QALYreg"=QALYreg,
                     "TCreg"=TCreg,"QALY_link"=QALY_link,"QALY_dist"=QALY_dist,
                     "TC_dist"=TC_dist,"TC_link"=TC_link,"cluster"=cluster, "unbalclus"=unbalclus)
  #compute overall list and return it as output from the function
  res_ec_tsb_list <- list("QALY_boot"=res_e_tsb_list,"TC_boot"=res_c_tsb_list,"inputs"=input_list)
  class(res_ec_tsb_list) <- "tsbootCE_glm"
  return(res_ec_tsb_list)
}

