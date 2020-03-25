#R code used to generate marginal mean estimates after fitting the model
#in JAGS and saving the posterior values of the parameters

# N.sim=20,000 number of the iterations 
# N.rep=40,000 replicated samples at each iteration
# r1=completers pattern
# r2=non-completers pattern 
# r.con=9  observed patterns in the control
# r.int=6  observed patterns in the intervention

## compute Monte Carlo estimates for mean uj and cj at each time j=0,1,2

mu.c0.r1<-mu.c1.r1<-mu.c2.r1<-matrix(NA,N.sim,2)
mu.u0.r1<-mu.u1.r1<-mu.u2.r1<-matrix(NA,N.sim,2)
mu.c0.r2<-mu.c1.r2<-mu.c2.r2<-matrix(NA,N.sim,2)
mu.u0.r2<-mu.u1.r2<-mu.u2.r2<-matrix(NA,N.sim,2)

# control

for(k in 1:N.sim){
  mu.c0.r1[k,1]<-mean(rlnorm(N.rep,nu0.c.r1[k,1],sigma0.c.r1[k,1]))
  mu.c1.r1[k,1]<-mean(rlnorm(N.rep,nu1.c.r1[k,1],sigma1.c.r1[k,1])) 
  mu.c2.r1[k,1]<-mean(rlnorm(N.rep,nu2.c.r1[k,1],sigma2.c.r1[k,1])) 
  mu.u0.r1[k,1]<-mean(rbeta(N.rep,a0.r1[k,1],b0.r1[k,1]))
  mu.u1.r1[k,1]<-mean(rbeta(N.rep,a1.r1[k,1],b1.r1[k,1])) 
  mu.u2.r1[k,1]<-mean(rbeta(N.rep,a2.r1[k,1],b2.r1[k,1])) 
  
  mu.c0.r2[k,1]<-mean(rlnorm(N.rep,nu0.c.r2[k,1],sigma0.c.r2[k,1]))
  mu.c1.r2[k,1]<-mean(rlnorm(N.rep,nu1.c.r2[k,1],sigma1.c.r2[k,1])) 
  mu.c2.r2[k,1]<-mean(rlnorm(N.rep,nu2.c.r2[k,1],sigma2.c.r2[k,1])) 
  mu.u0.r2[k,1]<-mean(rbeta(N.rep,a0.r2[k,1],b0.r2[k,1]))
  mu.u1.r2[k,1]<-mean(rbeta(N.rep,a1.r2[k,1],b1.r2[k,1])) 
  mu.u2.r2[k,1]<-mean(rbeta(N.rep,a2.r2[k,1],b2.r2[k,1])) 
  
  
  # intervention
  
  mu.c0.r1[k,2]<-mean(rlnorm(N.rep,nu0.c.r1[k,2],sigma0.c.r1[k,2]))
  mu.c1.r1[k,2]<-mean(rlnorm(N.rep,nu1.c.r1[k,2],sigma1.c.r1[k,2])) 
  mu.c2.r1[k,2]<-mean(rlnorm(N.rep,nu2.c.r1[k,2],sigma2.c.r1[k,2])) 
  mu.u0.r1[k,2]<-mean(rbeta(N.rep,a0.r1[k,2],b0.r1[k,2]))
  mu.u1.r1[k,2]<-mean(rbeta(N.rep,a1.r1[k,2],b1.r1[k,2])) 
  mu.u2.r1[k,2]<-mean(rbeta(N.rep,a2.r1[k,2],b2.r1[k,2])) 
  
  mu.c0.r2[k,2]<-mean(rlnorm(N.rep,nu0.c.r2[k,2],sigma0.c.r2[k,2]))
  mu.c1.r2[k,2]<-mean(rlnorm(N.rep,nu1.c.r2[k,2],sigma1.c.r2[k,2])) 
  mu.c2.r2[k,2]<-mean(rlnorm(N.rep,nu2.c.r2[k,2],sigma2.c.r2[k,2])) 
  mu.u0.r2[k,2]<-mean(rbeta(N.rep,a0.r2[k,2],b0.r2[k,2]))
  mu.u1.r2[k,2]<-mean(rbeta(N.rep,a1.r2[k,2],b1.r2[k,2])) 
  mu.u2.r2[k,2]<-mean(rbeta(N.rep,a2.r2[k,2],b2.r2[k,2])) 
}


## priors on the sensitivity parameters 

# flat prior

Delta.c0.flat<-runif(N.sim,0,2*sd(c0))
Delta.u0.flat<-runif(N.sim,-2*sd(u0),0)
Delta.c1.flat<-runif(N.sim,0,2*sd(c1))
Delta.u1.flat<-runif(N.sim,-2*sd(u1),0)
Delta.c2.flat<-runif(N.sim,0,2*sd(c2))
Delta.u2.flat<-runif(N.sim,-2*sd(u2),0)

# skew0 prior

Delta.c0.skew0<- 2*sd(c0)*(1-sqrt(runif(N.sim,0,1)))
Delta.u0.skew0<- -2*sd(u0)*(1-sqrt(runif(N.sim,0,1)))
Delta.c1.skew0<- 2*sd(c1)*(1-sqrt(runif(N.sim,0,1)))
Delta.u1.skew0<- -2*sd(u1)*(1-sqrt(runif(N.sim,0,1)))
Delta.c2.skew0<- 2*sd(c2)*(1-sqrt(runif(N.sim,0,1)))
Delta.u2.skew0<- -2*sd(u2)*(1-sqrt(runif(N.sim,0,1)))

#skew1 prior
Delta.c0.skew1<- 2*sd(c0)*sqrt(runif(N.sim,0,1))
Delta.u0.skew1<- -2*sd(u0)*sqrt(runif(N.sim,0,1))
Delta.c1.skew1<- 2*sd(c1)*sqrt(runif(N.sim,0,1))
Delta.u1.skew1<- -2*sd(u1)*sqrt(runif(N.sim,0,1))
Delta.c2.skew1<- 2*sd(c2)*sqrt(runif(N.sim,0,1))
Delta.u2.skew1<- -2*sd(u2)*sqrt(runif(N.sim,0,1))
                                                                                                                                                                                                                                                                                                                                                                  
## stratify the marginal means by pattern and arm under Delta=0
# control

mu.c0.r.con<-mu.c1.r.con<-mu.c2.r.con<-matrix(NA,N.sim,9)
mu.c0.r.con[,1]<-mu.c0.r1[,1]
mu.c0.r.con[,2:9]<-mu.c0.r2[,1]
mu.c1.r.con[,1]<-mu.c1.r1[,1]
mu.c1.r.con[,2:9]<-mu.c1.r2[,1]
mu.c2.r.con[,1]<-mu.c2.r1[,1]
mu.c2.r.con[,2:9]<-mu.c2.r2[,1]

mu.u0.r.con<-mu.u1.r.con<-mu.u2.r.con<-matrix(NA,N.sim,9)
mu.u0.r.con[,1]<-mu.u0.r1[,1]
mu.u0.r.con[,2:9]<-mu.u0.r2[,1]
mu.u1.r.con[,1]<-mu.u1.r1[,1]
mu.u1.r.con[,2:9]<-mu.u1.r2[,1]
mu.u2.r.con[,1]<-mu.u2.r1[,1]
mu.u2.r.con[,2:9]<-mu.u2.r2[,1]

# intervention

mu.c0.r.int<-mu.c1.r.int<-mu.c2.r.int<-matrix(NA,N.sim,6)
mu.c0.r.int[,1]<-mu.c0.r1[,2]
mu.c0.r.int[,2:6]<-mu.c0.r2[,2]
mu.c1.r.int[,1]<-mu.c1.r1[,2]
mu.c1.r.int[,2:6]<-mu.c1.r2[,2]
mu.c2.r.int[,1]<-mu.c2.r1[,2]
mu.c2.r.int[,2:6]<-mu.c2.r2[,2]

mu.u0.r.int<-mu.u1.r.int<-mu.u2.r.int<-matrix(NA,N.sim,6)
mu.u0.r.int[,1]<-mu.u0.r1[,2]
mu.u0.r.int[,2:6]<-mu.u0.r2[,2]
mu.u1.r.int[,1]<-mu.u1.r1[,2]
mu.u1.r.int[,2:6]<-mu.u1.r2[,2]
mu.u2.r.int[,1]<-mu.u2.r1[,2]
mu.u2.r.int[,2:6]<-mu.u2.r2[,2]

## compute marginal means for missing data in each pattern and arm under each prior

mu.c1.r.con.flat<-mu.c1.r.con.skew0<-mu.c1.r.con.skew1<-mu.c1.r.con
mu.c2.r.con.flat<-mu.c2.r.con.skew0<-mu.c2.r.con.skew1<-mu.c2.r.con
mu.u0.r.con.flat<-mu.u0.r.con.skew0<-mu.u0.r.con.skew1<-mu.u0.r.con
mu.u1.r.con.flat<-mu.u1.r.con.skew0<-mu.u1.r.con.skew1<-mu.u1.r.con
mu.u2.r.con.flat<-mu.u2.r.con.skew0<-mu.u2.r.con.skew1<-mu.u2.r.con
mu.c1.r.int.flat<-mu.c1.r.int.skew0<-mu.c1.r.int.skew1<-mu.c1.r.int
mu.c2.r.int.flat<-mu.c2.r.int.skew0<-mu.c2.r.int.skew1<-mu.c2.r.int
mu.u0.r.int.flat<-mu.u0.r.int.skew0<-mu.u0.r.int.skew1<-mu.u0.r.int
mu.u1.r.int.flat<-mu.u1.r.int.skew0<-mu.u1.r.int.skew1<-mu.u1.r.int
mu.u2.r.int.flat<-mu.u2.r.int.skew0<-mu.u2.r.int.skew1<-mu.u2.r.int

# control

# j=1 (missing for r.con=8,9)

mu.c1.r.con.flat[,c(8,9)]<-mu.c1.r2[,1]+Delta.c1.flat
mu.c1.r.con.skew0[,c(8,9)]<-mu.c1.r2[,1]+Delta.c1.skew0
mu.c1.r.con.skew1[,c(8,9)]<-mu.c1.r2[,1]+Delta.c1.skew1
# j=2 (missing for r.con=6,9)
mu.c2.r.con.flat[,c(6,9)]<-mu.c2.r2[,1]+Delta.c2.flat
mu.c2.r.con.skew0[,c(6,9)]<-mu.c2.r2[,1]+Delta.c2.skew0
mu.c2.r.con.skew1[,c(6,9)]<-mu.c2.r2[,1]+Delta.c2.skew1

# j=0 (missing for r.con=2,5)

mu.u0.r.con.flat[,c(2,5)]<-mu.u0.r2[,1]+Delta.u0.flat
mu.u0.r.con.skew0[,c(2,5)]<-mu.u0.r2[,1]+Delta.u0.skew0
mu.u0.r.con.skew1[,c(2,5)]<-mu.u0.r2[,1]+Delta.u0.skew1

# j=1  (missing for r.con=3,5,7,8,9)

mu.u1.r.con.flat[,c(3,5,7,8,9)]<-mu.u1.r2[,1]+Delta.u1.flat
mu.u1.r.con.skew0[,c(3,5,7,8,9)]<-mu.u1.r2[,1]+Delta.u1.skew0
mu.u1.r.con.skew1[,c(3,5,7,8,9)]<-mu.u1.r2[,1]+Delta.u1.skew1

# j=2  (missing for r.con=4,6,7,9)

mu.u2.r.con.flat[,c(4,6,7,9)]<-mu.u2.r2[,1]+Delta.u2.flat
mu.u2.r.con.skew0[,c(4,6,7,9)]<-mu.u2.r2[,1]+Delta.u2.skew0
mu.u2.r.con.skew1[,c(4,6,7,9)]<-mu.u2.r2[,1]+Delta.u2.skew1

# intervention

# j=1 (missing for r.int=5,6)

mu.c1.r.int.flat[,c(5,6)]<-mu.c1.r2[,2]+Delta.c1.flat
mu.c1.r.int.skew0[,c(5,6)]<-mu.c1.r2[,2]+Delta.c1.skew0
mu.c1.r.int.skew1[,c(5,6)]<-mu.c1.r2[,2]+Delta.c1.skew1

# j=2 (missing for r.int=6)

mu.c2.r.int.flat[,6]<-mu.c2.r2[,2]+Delta.c2.flat
mu.c2.r.int.skew0[,6]<-mu.c2.r2[,2]+Delta.c2.skew0
mu.c2.r.int.skew1[,6]<-mu.c2.r2[,2]+Delta.c2.skew1

# j=0 (missing for r.int=2)

mu.u0.r.int.flat[,2]<-mu.u0.r2[,2]+Delta.u0.flat
mu.u0.r.int.skew0[,2]<-mu.u0.r2[,2]+Delta.u0.skew0
mu.u0.r.int.skew1[,2]<-mu.u0.r2[,2]+Delta.u0.skew1

# j=1  (missing for r.int=3,5,6)

mu.u1.r.int.flat[,c(3,5,6)]<-mu.u1.r2[,2]+Delta.u1.flat
mu.u1.r.int.skew0[,c(3,5,6)]<-mu.u1.r2[,2]+Delta.u1.skew0
mu.u1.r.int.skew1[,c(3,5,6)]<-mu.u1.r2[,2]+Delta.u1.skew1

# j=2  (missing for r.int=4,6,7,9)

mu.u2.r.int.flat[,c(4,6)]<-mu.u2.r2[,2]+Delta.u2.flat
mu.u2.r.int.skew0[,c(4,6)]<-mu.u2.r2[,2]+Delta.u2.skew0
mu.u2.r.int.skew1[,c(4,6)]<-mu.u2.r2[,2]+Delta.u2.skew1

## compute marginal means across patterns in each group under each prior

mu.c0.bench<-mu.c1.bench<-mu.c2.bench<-matrix(NA,N.sim,2)
mu.u0.bench<-mu.u1.bench<-mu.u2.bench<-matrix(NA,N.sim,2)
mu.c0.flat<-mu.c1.flat<-mu.c2.flat<-matrix(NA,N.sim,2)
mu.u0.flat<-mu.u1.flat<-mu.u2.flat<-matrix(NA,N.sim,2)
mu.c0.skew0<-mu.c1.skew0<-mu.c2.skew0<-matrix(NA,N.sim,2)
mu.u0.skew0<-mu.u1.skew0<-mu.u2.skew0<-matrix(NA,N.sim,2)
mu.c0.skew1<-mu.c1.skew1<-mu.c2.skew1<-matrix(NA,N.sim,2)
mu.u0.skew1<-mu.u1.skew1<-mu.u2.skew1<-matrix(NA,N.sim,2)

# control

mu.c0.bench[,1]<-rowSums(mu.c0.r.con*psi[,1])
mu.c1.bench[,1]<-rowSums(mu.c1.r.con*psi[,1])
mu.c2.bench[,1]<-rowSums(mu.c2.r.con*psi[,1])
mu.u0.bench[,1]<-rowSums(mu.u0.r.con*psi[,1])
mu.u1.bench[,1]<-rowSums(mu.u1.r.con*psi[,1])
mu.u2.bench[,1]<-rowSums(mu.u2.r.con*psi[,1])
mu.c1.flat[,1]<-rowSums(mu.c1.r.con.flat*psi[,1])
mu.c2.flat[,1]<-rowSums(mu.c2.r.con.flat*psi[,1])
mu.u0.flat[,1]<-rowSums(mu.u0.r.con.flat*psi[,1])
mu.u1.flat[,1]<-rowSums(mu.u1.r.con.flat*psi[,1])
mu.u2.flat[,1]<-rowSums(mu.u2.r.con.flat*psi[,1])
mu.c1.skew0[,1]<-rowSums(mu.c1.r.con.skew0*psi[,1])
mu.c2.skew0[,1]<-rowSums(mu.c2.r.con.skew0*psi[,1])
mu.u0.skew0[,1]<-rowSums(mu.u0.r.con.skew0*psi[,1])
mu.u1.skew0[,1]<-rowSums(mu.u1.r.con.skew0*psi[,1])
mu.u2.skew0[,1]<-rowSums(mu.u2.r.con.skew0*psi[,1])
mu.c1.skew1[,1]<-rowSums(mu.c1.r.con.skew1*psi[,1])
mu.c2.skew1[,1]<-rowSums(mu.c2.r.con.skew1*psi[,1])
mu.u0.skew1[,1]<-rowSums(mu.u0.r.con.skew1*psi[,1])
mu.u1.skew1[,1]<-rowSums(mu.u1.r.con.skew1*psi[,1])
mu.u2.skew1[,1]<-rowSums(mu.u2.r.con.skew1*psi[,1])

# intervention

mu.c0.bench[,2]<-rowSums(mu.c0.r.int*psi[,2])
mu.c1.bench[,2]<-rowSums(mu.c1.r.int*psi[,2])
mu.c2.bench[,2]<-rowSums(mu.c2.r.int*psi[,2])
mu.u0.bench[,2]<-rowSums(mu.u0.r.int*psi[,2])
mu.u1.bench[,2]<-rowSums(mu.u1.r.int*psi[,2])
mu.u2.bench[,2]<-rowSums(mu.u2.r.int*psi[,2])
mu.c1.flat[,2]<-rowSums(mu.c1.r.int.flat*psi[,2])
mu.c2.flat[,2]<-rowSums(mu.c2.r.int.flat*psi[,2])
mu.u0.flat[,2]<-rowSums(mu.u0.r.int.flat*psi[,2])
mu.u1.flat[,2]<-rowSums(mu.u1.r.int.flat*psi[,2])
mu.u2.flat[,2]<-rowSums(mu.u2.r.int.flat*psi[,2])
mu.c1.skew0[,2]<-rowSums(mu.c1.r.int.skew0*psi[,2])
mu.c2.skew0[,2]<-rowSums(mu.c2.r.int.skew0*psi[,2])
mu.u0.skew0[,2]<-rowSums(mu.u0.r.int.skew0*psi[,2])
mu.u1.skew0[,2]<-rowSums(mu.u1.r.int.skew0*psi[,2])
mu.u2.skew0[,2]<-rowSums(mu.u2.r.int.skew0*psi[,2])
mu.c1.skew1[,2]<-rowSums(mu.c1.r.int.skew1*psi[,2])
mu.c2.skew1[,2]<-rowSums(mu.c2.r.int.skew1*psi[,2])
mu.u0.skew1[,2]<-rowSums(mu.u0.r.int.skew1*psi[,2])
mu.u1.skew1[,2]<-rowSums(mu.u1.r.int.skew1*psi[,2])
mu.u2.skew1[,2]<-rowSums(mu.u2.r.int.skew1*psi[,2])                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                                                  
                                                                                                                                                                                                                                                                                                                                                                  
                                                                                                                                                                                                                                                                                                                                                                  