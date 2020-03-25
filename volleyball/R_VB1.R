#R code for fitting the basic model (VB1) using STAN and summarise
#results using plots and tables

library(plyr)
library(dplyr)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#load data
data<-read.csv(file = "Volley_prepro_v1.csv",header = T)
data$home.team<-as.factor(data$home.team)
data$away.team<-as.factor(data$away.team)
data$bloeff1<-(data$bloper1-data$bloinv1)/(data$blo1+data$bloper1)
data$bloeff2<-(data$bloper2-data$bloinv2)/(data$blo2+data$bloper2)
data_stat1<-ddply(data, .(h), summarise, meanatteff1=mean(atteff1),meansereff1=mean(sereff1),
                  meandefeff1=mean(defeff1),meanbloeff1=mean(bloeff1))
data_stat2<-ddply(data, .(a), summarise, meanatteff2=mean(atteff2),meansereff2=mean(sereff2),
                  meandefeff2=mean(defeff2),meanbloeff2=mean(bloeff2))
#center cov grand mean
atteff1.cen<-data$atteff1-mean(data$atteff1)
sereff1.cen<-data$sereff1-mean(data$sereff1)
bloeff1.cen<-data$bloeff1-mean(data$bloeff1)
defeff1.cen<-data$defeff1-mean(data$defeff1)

atteff2.cen<-data$atteff2-mean(data$atteff2)
sereff2.cen<-data$sereff2-mean(data$sereff2)
bloeff2.cen<-data$bloeff2-mean(data$bloeff2)
defeff2.cen<-data$defeff2-mean(data$defeff2)

#prepare data
y1<-data$y1
y2<-data$y2
ngames<-max(data$Game)
nteams<-max(data$h)
home_team<-data$h
away_team<-data$a

att_eff1<-data$atteff1
att_eff2<-data$atteff2
ser_eff1<-data$sereff1
ser_eff2<-data$sereff2
blo_eff1<-data$bloeff1
blo_eff2<-data$bloeff2
def_eff1<-data$defeff1
def_eff2<-data$defeff2

att_eff1<-atteff1.cen
att_eff2<-atteff2.cen
ser_eff1<-sereff1.cen
ser_eff2<-sereff2.cen
blo_eff1<-bloeff1.cen
blo_eff2<-bloeff2.cen
def_eff1<-defeff1.cen
def_eff2<-defeff2.cen

ds<-ifelse(data$settot==5,1,0)
dg<-ifelse(data$set1>data$set2,1,0)

#pre-processing
datalist <- list(y1=y1,y2=y2,ngames=ngames,nteams=nteams,home_team=home_team,away_team=away_team,
                    att_eff1=att_eff1,att_eff2=att_eff2,def_eff1=def_eff1,def_eff2=def_eff2,
                    ser_eff1=ser_eff1,ser_eff2=ser_eff2,blo_eff1=blo_eff1,blo_eff2=blo_eff2,
                    ds=ds, dg=dg)
params <- c("mu","home","gamma","delta","beta0_att","beta0_def","beta1_att","beta1_def",
            "beta1_ser","beta1_blo","theta1","theta2","loglik_y1","loglik_y2","loglik_ds","loglik_dg")
burnInSteps = 5000
nChains = 2
numSavedSteps = 20000
thinSteps = 1
nIter = ceiling((numSavedSteps * thinSteps)/nChains)

set.seed(3456)
model_stan<- stan(data = datalist, file = "STAN_VB1.stan", 
                       chains = nChains, pars = params, iter = nIter, 
                       warmup = burnInSteps, thin = thinSteps)
print(model_stan, pars = c("mu","home","gamma","delta","beta0_att","beta0_def","beta1_att","beta1_def",
                           "beta1_ser","beta1_blo"))

model_stan_par<-extract(model_stan)



#plot att vs def effext by team
beta0.att<-apply(model_stan_par$beta0_att, 2, mean)
beta0.def<-apply(model_stan_par$beta0_def, 2, mean)
names<-unique(data.frame(data$home.team,data$h)[order(data$h),][,1])

plot(beta0.att,beta0.def, main = "", type = "n", xlab = "Mean attack effect", ylab = "Mean defence effect", 
     xlim=c(-0.13,0.08), ylim=c(-0.18,0.1))
points(beta0.att,beta0.def, pch=16, col="red",cex = 1.2)
abline(v=0)
abline(h=0)
text(beta0.att,beta0.def,names,cex = 0.8, adj = c(0.4,1.3))


#prediction joint
library(boot)
y1.pred<-y2.pred<-matrix(NA,length(model_stan_par$home),132)
ds.pred<-dg.pred<-matrix(NA,length(model_stan_par$home),132)
pi.s<-pi.g<-matrix(NA,length(model_stan_par$home),132)
y1.mat<-y2.mat<-ds.mat<-matrix(NA,length(model_stan_par$home),132)
for(i in 1:length(model_stan_par$home)){
 y1.mat[i,]<-y1
 y2.mat[i,]<-y2
 ds.mat[i,]<-ds
 pi.s[i,]<-inv.logit(model_stan_par$gamma[i,1] + model_stan_par$gamma[i,2]*y1.mat[i,] + model_stan_par$gamma[i,3]*y2.mat[i,])
 pi.g[i,]<-inv.logit(model_stan_par$delta[i,1] + model_stan_par$delta[i,2]*y1.mat[i,] + model_stan_par$delta[i,3]*y2.mat[i,] + model_stan_par$delta[i,4]*ds.mat[i,])
}


set.seed(3456)
for(i in 1:length(model_stan_par$home)){
  y1.pred[i,]<-rpois(n=132,lambda = model_stan_par$theta1[i,])
  y2.pred[i,]<-rpois(n=132,lambda = model_stan_par$theta2[i,])
  ds.pred[i,]<-rbinom(n=132, size = 1,prob = pi.s[i,])
  dg.pred[i,]<-rbinom(n=132, size = 1,prob = pi.g[i,])
}

#predictions joint
hist(apply(y1.pred, 1, mean))
abline(v=mean(data$y1), col="red", lwd=2)

hist(apply(y2.pred, 1, mean))
abline(v=mean(data$y2), col="red", lwd=2)

hist(apply(ds.pred, 1, sum))
abline(v=sum(ds), col="red", lwd=2)

hist(apply(dg.pred, 1, sum))
abline(v=sum(dg), col="red", lwd=2)

#prediction points
results<-list()
for(i in 1:1000){
  results[[i]]<-data.frame(y1.pred[i,],y2.pred[i,],ds.pred[i,],dg.pred[i,],data$h,data$a) 
  results[[i]]$points.home<-ifelse(results[[i]]$ds.pred.i...==0 & results[[i]]$dg.pred.i...==1,3,0)
  results[[i]]$points.home<-ifelse(results[[i]]$ds.pred.i...==1 & results[[i]]$dg.pred.i...==1,2,results[[i]]$points.home)
  results[[i]]$points.home<-ifelse(results[[i]]$ds.pred.i...==1 & results[[i]]$dg.pred.i...==0,1,results[[i]]$points.home)
  
  results[[i]]$points.away<-ifelse(results[[i]]$ds.pred.i...==0 & results[[i]]$dg.pred.i...==0,3,0)
  results[[i]]$points.away<-ifelse(results[[i]]$ds.pred.i...==1 & results[[i]]$dg.pred.i...==0,2,results[[i]]$points.away)
  results[[i]]$points.away<-ifelse(results[[i]]$ds.pred.i...==1 & results[[i]]$dg.pred.i...==1,1,results[[i]]$points.away)
}

#compare results for scores by team
tot.y1.list<-tot.y2.list<-list()
for(i in 1:1000){
  tot.y1.list[[i]]<-ddply(results[[i]], .(data.h), summarise, totscorehome=sum(y1.pred.i...))
  tot.y2.list[[i]]<-ddply(results[[i]], .(data.a), summarise, totscoreaway=sum(y2.pred.i...))
}
tot.y1.list.neg<-tot.y2.list.neg<-list()
for(i in 1:1000){
  tot.y1.list.neg[[i]]<-ddply(results[[i]], .(data.h), summarise, totscorehomeneg=sum(y2.pred.i...))
  tot.y2.list.neg[[i]]<-ddply(results[[i]], .(data.a), summarise, totscoreawayneg=sum(y1.pred.i...))
}


tot.y1.list.mat<-tot.y2.list.mat<-matrix(NA,1000,12)
tot.y1.list.neg.mat<-tot.y2.list.neg.mat<-matrix(NA,1000,12)
for(i in 1:1000){
  tot.y1.list.mat[i,]<-tot.y1.list[[i]][,2]
  tot.y2.list.mat[i,]<-tot.y2.list[[i]][,2]
  tot.y1.list.neg.mat[i,]<-tot.y1.list.neg[[i]][,2]
  tot.y2.list.neg.mat[i,]<-tot.y2.list.neg[[i]][,2]
}
tot.y1.obs<-ddply(data, .(h), summarise, totscorehomeobs=sum(y1))
tot.y2.obs<-ddply(data, .(a), summarise, totscorehomeobs=sum(y2))
tot.y1.obs.neg<-ddply(data, .(h), summarise, totscorehomeobs=sum(y2))
tot.y2.obs.neg<-ddply(data, .(a), summarise, totscorehomeobs=sum(y1))

cbind(tot.y1.obs,apply(tot.y1.list.mat,2,median))
cbind(tot.y2.obs,apply(tot.y2.list.mat,2,median))
cbind(tot.y1.obs.neg,apply(tot.y1.list.neg.mat,2,median))
cbind(tot.y2.obs.neg,apply(tot.y2.list.neg.mat,2,median))

#scored
tot.y.obs<-tot.y1.obs[,2]+tot.y2.obs[,2]
tot.y.pred<-apply(tot.y1.list.mat,2,median)+apply(tot.y2.list.mat,2,median)
res.y<-cbind(tot.y.obs,tot.y.pred)
rownames(res.y)<-names
res.y<-round(res.y,digits = 0)

#conceded
tot.y.obs.neg<-tot.y1.obs.neg[,2]+tot.y2.obs.neg[,2]
tot.y.pred.neg<-apply(tot.y1.list.neg.mat,2,median)+apply(tot.y2.list.neg.mat,2,median)
res.y.neg<-cbind(tot.y.obs.neg,tot.y.pred.neg)
rownames(res.y.neg)<-names
res.y.neg<-round(res.y.neg,digits = 0)

#compare results for points
data.points.list<-list()
for(i in 1:1000){
  data.points.list[[i]]<-data.frame(data$Game)
  data.points.list[[i]]$Game<-data$Game
  data.points.list[[i]]$h<-data$h
  data.points.list[[i]]$a<-data$a
  data.points.list[[i]]$points.home<-results[[i]]$points.home
  data.points.list[[i]]$points.away<-results[[i]]$points.away
}
tot.home.list<-tot.away.list<-tot.team.list<-list()
for(i in 1:1000){
  tot.home.list[[i]]<-ddply(data.points.list[[i]], .(h), summarise, totpointhome=sum(points.home))
  tot.away.list[[i]]<-ddply(data.points.list[[i]], .(a), summarise, totpointaway=sum(points.away))
}

for(i in 1:1000){
  tot.team.list[[i]]<-data.frame(levels(data$home.team), tot.home.list[[i]]$h)
  tot.team.list[[i]]$tot.team<-tot.home.list[[i]]$totpointhome + tot.away.list[[i]]$totpointaway
  tot.team.list[[i]]$true.team<-c(19, 39, 23, 50, 19, 11, 37, 51, 32, 33, 27, 50)
  tot.team.list[[i]]<-tot.team.list[[i]][order(tot.team.list[[i]]$tot.team, decreasing = TRUE),]
  tot.team.list[[i]]$rank<-rep(1:12)
}

#plot total scores by team
tot.scores<-matrix(NA,1000,12)
colnames(tot.scores)<-names
for(i in 1:1000){
  for(j in 1:12){
    tot.scores[i,j]<-tot.team.list[[i]]$tot.team[tot.team.list[[i]]$tot.home.list..i...h==j] 
  }
}
tot.scores.obs<-c(19, 39, 23, 50, 19, 11, 37, 51, 32, 33, 27, 50)
tot.scores.med<-apply(tot.scores, 2, median)
tot.scores.final<-cbind(tot.scores.obs,tot.scores.med)

#plot total wins 
tot.wins<-matrix(NA,1000,12)
colnames(tot.wins)<-names
for(i in 1:1000){
  for(j in 1:12){
    tot.wins[i,j]<-length(data.points.list[[i]]$points.home[data.points.list[[i]]$points.home>1 & 
                                                              data.points.list[[i]]$h==j]) + length(data.points.list[[i]]$points.away[data.points.list[[i]]$points.away>1 & 
                                                                                                                                        data.points.list[[i]]$a==j])
  }
}
tot.wins.obs<-c(7,12,6,17,7,5,13,17,10,12,8,18)
tot.wins.prop<-tot.wins/22
tot.wins.obs.prop<-tot.wins.obs/22
tot.wins.med<-apply(tot.wins, 2, median)
tot.wins.final<-cbind(tot.wins.obs,tot.wins.med)

#summarise pred results
res.final.obs<-cbind(res.y[,1],res.y.neg[,1],tot.wins.final[,1],tot.scores.final[,1])
res.final.pred<-cbind(res.y[,2],res.y.neg[,2],tot.wins.final[,2],tot.scores.final[,2])
res.final<-cbind(res.final.obs,res.final.pred)

library(xtable)
xtable(res.final, digits = 0)

#####################pic
#library(loo)
#log_lik<-model$BUGSoutput$sims.list$loglik_y1 + model$BUGSoutput$sims.list$loglik_y2 +
#  model$BUGSoutput$sims.list$loglik_ds + model$BUGSoutput$sims.list$loglik_dg
#waic<-waic(log_lik)
#looic<-loo(log_lik)

res.matrix<-matrix(NA,length(tot.team.list),12)
colnames(res.matrix)<-names
for(i in 1:1000){
  res.matrix[i,1]<-tot.team.list[[i]]$rank[tot.team.list[[i]]$tot.home.list..i...h==1]
  res.matrix[i,2]<-tot.team.list[[i]]$rank[tot.team.list[[i]]$tot.home.list..i...h==2]
  res.matrix[i,3]<-tot.team.list[[i]]$rank[tot.team.list[[i]]$tot.home.list..i...h==3]
  res.matrix[i,4]<-tot.team.list[[i]]$rank[tot.team.list[[i]]$tot.home.list..i...h==4]
  res.matrix[i,5]<-tot.team.list[[i]]$rank[tot.team.list[[i]]$tot.home.list..i...h==5]
  res.matrix[i,6]<-tot.team.list[[i]]$rank[tot.team.list[[i]]$tot.home.list..i...h==6]
  res.matrix[i,7]<-tot.team.list[[i]]$rank[tot.team.list[[i]]$tot.home.list..i...h==7]
  res.matrix[i,8]<-tot.team.list[[i]]$rank[tot.team.list[[i]]$tot.home.list..i...h==8]
  res.matrix[i,9]<-tot.team.list[[i]]$rank[tot.team.list[[i]]$tot.home.list..i...h==9]
  res.matrix[i,10]<-tot.team.list[[i]]$rank[tot.team.list[[i]]$tot.home.list..i...h==10]
  res.matrix[i,11]<-tot.team.list[[i]]$rank[tot.team.list[[i]]$tot.home.list..i...h==11]
  res.matrix[i,12]<-tot.team.list[[i]]$rank[tot.team.list[[i]]$tot.home.list..i...h==12]
}

#create stacked barplot of results
data.barplot<-data.frame(rep(1:c(1000*12)))
names(data.barplot)<-c("Game")
data.barplot$position<-as.factor(c(res.matrix[,1],res.matrix[,2],res.matrix[,3],res.matrix[,4],
                                   res.matrix[,5],res.matrix[,6],res.matrix[,7],res.matrix[,8],
                                   res.matrix[,9],res.matrix[,10],res.matrix[,11],res.matrix[,12]))
data.barplot$team<-rep(NA,1000*12)
data.barplot$team[1:1000]<-rep(paste(names[1]),1000)
data.barplot$team[1001:2000]<-rep(paste(names[2]),1000)
data.barplot$team[2001:3000]<-rep(paste(names[3]),1000)
data.barplot$team[3001:4000]<-rep(paste(names[4]),1000)
data.barplot$team[4001:5000]<-rep(paste(names[5]),1000)
data.barplot$team[5001:6000]<-rep(paste(names[6]),1000)
data.barplot$team[6001:7000]<-rep(paste(names[7]),1000)
data.barplot$team[7001:8000]<-rep(paste(names[8]),1000)
data.barplot$team[8001:9000]<-rep(paste(names[9]),1000)
data.barplot$team[9001:10000]<-rep(paste(names[10]),1000)
data.barplot$team[10001:11000]<-rep(paste(names[11]),1000)
data.barplot$team[11001:12000]<-rep(paste(names[12]),1000)
#data.barplot$team<-as.factor(data.barplot$team)
data.barplot$team <-factor(data.barplot$team, levels = c("Novara", "Scandicci","Conegliano", "Monza","Busto Arsizio",
                                                         "Pesaro","Piacenza", "San Casciano","Casalmaggiore", "Bergamo",
                                                         "Filottrano", "Legnano"))



data.barplot$match<-c(rep(1,1000),rep(2,1000),rep(3,1000),rep(4,1000),rep(5,1000),rep(6,1000),
                      rep(7,1000),rep(8,1000),rep(9,1000),rep(10,1000),rep(11,1000),rep(12,1000))
data.barplot$match<-as.factor(data.barplot$match)
data.barplot$Game<-rep(1,nrow(data.barplot))
data.barplot$area<-ifelse(data.barplot$position==1|data.barplot$position==2|data.barplot$position==3,"high","middle")
data.barplot$area<-ifelse(data.barplot$position==12|data.barplot$position==11|data.barplot$position==10,"low",data.barplot$area)
data.barplot$area<-as.factor(data.barplot$area)
data.barplot$area<-ordered(data.barplot$area,levels=c("low","middle","high"))
data.barplot$team<-factor(data.barplot$team,levels = rev(levels(data.barplot$team)))

df.summary1<-ddply(data.barplot,.(team,position),summarise,count=sum(Game), percent=sum(Game)/1000)
df.summary2<-ddply(data.barplot,.(team,area),summarise,count=sum(Game), percent=sum(Game)/1000)

library(ggplot2)
library(scales)
library(RColorBrewer)
library(viridis)

plot1<-ggplot(df.summary1, aes(x=team, y=percent, fill=position)) +
  geom_bar(stat="identity", width = 0.7, colour="black", lwd=0.1) +
  geom_text(aes(label=ifelse(percent >= 0.1, paste0(sprintf("%.0f", percent*100),"%"),"")),
            position=position_stack(vjust=0.5), colour="white") +
  coord_flip() + scale_y_continuous(labels = percent_format()) +
  labs(y="", x="") + scale_fill_viridis(discrete = T) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"))

#########plot for cumlative points over simulated season
points.list<-cum.points.list<-list()
for(i in 1:nrow(res.matrix)){
  points.list[[i]]<-matrix(NA,22,12)
  cum.points.list[[i]]<-matrix(NA,22,12)
  for(j in 1:12){
    points.list[[i]][,j]<-c(data.points.list[[i]]$points.home[data.points.list[[i]]$h==j],data.points.list[[i]]$points.away[data.points.list[[i]]$a==j])
  }
  colnames(points.list[[i]])<-unique(levels(data$home.team))
  rownames(points.list[[i]])<-rep(1:22)
  cum.points.list[[i]]<-apply(points.list[[i]], 2, cumsum)
}

#########plot for cumlative points over simulated season
points.list<-cum.points.list<-list()
for(i in 1:nrow(res.matrix)){
  points.list[[i]]<-matrix(NA,22,12)
  cum.points.list[[i]]<-matrix(NA,22,12)
  for(j in 1:12){
    points.list[[i]][,j]<-c(data.points.list[[i]]$points.home[data.points.list[[i]]$h==j],data.points.list[[i]]$points.away[data.points.list[[i]]$a==j])
  }
  colnames(points.list[[i]])<-unique(levels(data$home.team))
  rownames(points.list[[i]])<-rep(1:22)
  cum.points.list[[i]]<-apply(points.list[[i]], 2, cumsum)
}

#plot cumulative points obs vs pred
obs.cum.points<-matrix(NA,22,12)
colnames(obs.cum.points)<-unique(levels(data$home.team))
rownames(obs.cum.points)<-rep(1:22)
obs.cum.points[1,]<-c(0,1,0,3,3,3,0,2,0,0,3,3)
obs.cum.points[2,]<-c(0,3,0,3,0,0,3,3,0,3,0,3)
obs.cum.points[3,]<-c(0,2,3,3,0,0,0,3,3,1,0,3)
obs.cum.points[4,]<-c(0,3,0,3,0,0,3,3,1,2,0,3)
obs.cum.points[5,]<-c(0,3,0,3,0,3,0,3,3,0,0,3)
obs.cum.points[6,]<-c(0,3,1,3,0,0,0,3,1,3,2,2)
obs.cum.points[7,]<-c(2,1,3,2,0,1,1,3,2,2,1,0)
obs.cum.points[8,]<-c(0,2,1,1,0,0,3,2,3,2,1,3)
obs.cum.points[9,]<-c(3,3,1,3,0,2,2,1,0,0,3,0)
obs.cum.points[10,]<-c(1,2,1,2,1,2,1,1,2,2,1,2)
obs.cum.points[11,]<-c(0,1,0,3,0,0,3,3,0,3,3,2)
obs.cum.points[12,]<-c(3,0,0,3,0,0,3,3,3,0,0,3)
obs.cum.points[13,]<-c(0,1,3,3,2,0,3,2,1,3,0,0)
obs.cum.points[14,]<-c(2,3,1,3,0,0,0,3,3,0,0,3)
obs.cum.points[15,]<-c(2,1,0,3,2,0,3,1,0,3,0,3)
obs.cum.points[16,]<-c(0,3,0,3,0,0,0,3,3,0,3,3)
obs.cum.points[17,]<-c(2,0,3,0,0,3,3,1,0,3,0,3)
obs.cum.points[18,]<-c(0,0,0,3,3,1,2,3,1,2,3,0)
obs.cum.points[19,]<-c(3,3,3,2,0,0,3,0,0,1,0,3)
obs.cum.points[20,]<-c(0,0,2,0,3,1,0,3,3,0,3,3)
obs.cum.points[21,]<-c(0,3,0,1,3,0,2,2,0,1,3,3)
obs.cum.points[22,]<-c(1,1,1,0,2,0,2,3,3,2,1,2)

obs.cum.points<-apply(obs.cum.points, 2, cumsum)

par(mar=c(2.1, 3.1, 3.1, 3.1))
par(mfrow=c(4,3), mai = c(0.4, 0.4, 0.1, 0.2))
for(i in 1:12){
  plot(rep(1:22),obs.cum.points[,i], axes = F, type = "n", xlab = "games", ylab = "points",xlim = c(0,23),ylim = c(0,54))
  axis(1,at=c(0,5,10,15,20,25),labels = c(0,5,10,15,20,25))
  axis(2,at=c(0,10,20,30,40,50),labels = c(0,10,20,30,40,50))
  lines(rep(1:22),obs.cum.points[,i], lty=1,lwd=1,col="black")
  lines(rep(1:22), cum.points.list[[7]][,i], lty=1,lwd=1,col="red")
  text(5,50,unique(levels(names))[i],pos = 1, cex = 1)
}












