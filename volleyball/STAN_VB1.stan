
data{
int<lower=1> nteams; // number of teams
int<lower=1> ngames; // number of games
int<lower=1, upper=nteams> home_team[ngames]; // home team ID (1, ..., 12)
int<lower=1, upper=nteams> away_team[ngames]; // away team ID (1, ..., 12)
vector [ngames] att_eff1; // in game statistics for number of attacks (home)
vector [ngames] ser_eff1; // in game statistics for number of serves (home)
vector [ngames] def_eff2; // in game statistics for number of digs (away)
vector [ngames] blo_eff2; // in game statistics for number of blocks (away)
vector [ngames] att_eff2; // in game statistics for number of attacks (away)
vector [ngames] ser_eff2; // in game statistics for number of serves (away)
vector [ngames] def_eff1; // in game statistics for number of digs (home)
vector [ngames] blo_eff1; // in game statistics for number of blocks (home)
int<lower=0> y1[ngames]; // number of points scored by home team
int<lower=0> y2[ngames]; // number of points scored by away team
int<lower=0, upper=1> ds[ngames]; // indicator for number of sets played (3/4 or 5)
int<lower=0, upper=1> dg[ngames]; // indicator for winning of the match for home team
}
parameters{
real home;
real mu;
real mu0_att;
real mu0_def;
real<lower=0> tau0_att;
real<lower=0> tau0_def;
real mu1_att;
real mu1_def;
real<lower=0> tau1_att;
real<lower=0> tau1_def;
real mu1_ser;
real mu1_blo;
real<lower=0> tau1_ser;
real<lower=0> tau1_blo;
vector [nteams] beta0_att_star;
vector [nteams] beta1_att_star;
vector [nteams] beta1_ser_star;
vector [nteams] beta0_def_star;
vector [nteams] beta1_def_star;
vector [nteams] beta1_blo_star;
real gamma[3];
real delta[4];
}
transformed parameters{
//Trick to code the sum-to-zero constraint
vector [nteams] beta0_att;
vector [nteams] beta1_att;
vector [nteams] beta1_ser;
vector [nteams] beta0_def;
vector [nteams] beta1_def;
vector [nteams] beta1_blo;
vector<lower=0>[ngames] theta1;
vector<lower=0>[ngames] theta2;
beta0_att = beta0_att_star - mean(beta0_att_star);
beta1_att = beta1_att_star - mean(beta1_att_star);
beta1_ser = beta1_ser_star - mean(beta1_ser_star);
beta0_def = beta0_def_star - mean(beta0_def_star);
beta1_def = beta1_def_star - mean(beta1_def_star);
beta1_blo = beta1_blo_star - mean(beta1_blo_star);
 for (g in 1:ngames) {
    theta1[g] = exp(home + mu + beta0_att[home_team[g]] + beta1_att[home_team[g]]*att_eff1[g] + beta1_ser[home_team[g]]*ser_eff1[g] + 
             beta0_def[away_team[g]] + beta1_def[away_team[g]]*def_eff2[g] + beta1_blo[away_team[g]]*blo_eff2[g]); 
    theta2[g] = exp(home + mu + beta0_att[away_team[g]] + beta1_att[away_team[g]]*att_eff1[g] + beta1_ser[away_team[g]]*ser_eff1[g] + 
             beta0_def[home_team[g]] + beta1_def[home_team[g]]*def_eff2[g] + beta1_blo[home_team[g]]*blo_eff2[g]); 
 }
}
model {
//priors
home ~ normal(0, 100);
mu ~ normal(0, 100);
mu0_att ~ normal(0, 100);
mu0_def ~ normal(0, 100);
tau0_att ~ gamma(0.01, 0.01);
tau0_def ~ gamma(0.01, 0.01);
mu1_att ~ normal(0, 100);
mu1_def ~ normal(0, 100);
tau1_att ~ gamma(0.01, 0.01);
tau1_def ~ gamma(0.01, 0.01);
mu1_ser ~ normal(0, 100);
mu1_blo ~ normal(0, 100);
tau1_ser ~ gamma(0.01, 0.01);
tau1_blo ~ gamma(0.01, 0.01);
gamma ~ normal(0, 100);
delta ~ normal(0, 100);
beta0_att_star ~ normal(mu0_att, tau0_att);
beta0_def_star ~ normal(mu0_def, tau0_def);
beta1_att_star ~ normal(mu1_att, tau1_att);
beta1_def_star ~ normal(mu1_def, tau1_def);
beta1_ser_star ~ normal(mu1_ser, tau1_ser);
beta1_blo_star ~ normal(mu1_blo, tau1_blo);
// likelihood
 for (g in 1:ngames) {
  y1[g] ~ poisson(theta1[g]);
  y2[g] ~ poisson(theta2[g]);
  ds[g] ~ bernoulli_logit(gamma[1] + gamma[2]*y1[g] + gamma[3]*y2[g]);
  dg[g] ~ bernoulli_logit(delta[1] + delta[2]*y1[g] + delta[3]*y2[g] + delta[4]*ds[g]);
 }
}
generated quantities{
// loglikelihood 
vector[ngames] loglik_y1;
vector[ngames] loglik_y2;
vector[ngames] loglik_ds;
vector[ngames] loglik_dg;
 for (g in 1:ngames) {
  loglik_y1[g] = poisson_lpmf(y1[g]| theta1[g]);
  loglik_y2[g] = poisson_lpmf(y2[g]| theta2[g]);
  loglik_ds[g] = bernoulli_logit_lpmf(ds[g]| gamma[1] + gamma[2]*y1[g] + gamma[3]*y2[g]);
  loglik_dg[g] = bernoulli_logit_lpmf(dg[g]| delta[1] + delta[2]*y1[g] + delta[3]*y2[g] + delta[4]*ds[g]);
 }
}


