
data{
int<lower=1> nteams; // number of teams
int<lower=1> ngames; // number of games
int<lower=1, upper=nteams> home_team[ngames]; // home team ID (1, ..., 12)
int<lower=1, upper=nteams> away_team[ngames]; // away team ID (1, ..., 12)
int<lower=1> df; // degress of freedom for IW 
cov_matrix[3] IW_att; // scale matrix of IW fpr att parameters
cov_matrix[3] IW_def; // scale matrix of IW fpr def parameters
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
real xi_beta0_att;
real xi_beta0_def;
real xi_beta1_att;
real xi_beta1_def;
real xi_beta1_ser;
real xi_beta1_blo;
real mu_raw_beta0_att;
real mu_raw_beta1_att;
real mu_raw_beta1_ser;
real mu_raw_beta0_def;
real mu_raw_beta1_def;
real mu_raw_beta1_blo;
matrix [nteams,3] B_raw_att;
matrix [nteams,3] B_raw_def;
cov_matrix[3] Tau_B_raw_att;
cov_matrix[3] Tau_B_raw_def;
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
vector [nteams] beta0_att_star;
vector [nteams] beta1_att_star;
vector [nteams] beta1_ser_star;
vector [nteams] beta0_def_star;
vector [nteams] beta1_def_star;
vector [nteams] beta1_blo_star;
matrix [nteams,3] B_raw_hat_att;
matrix [nteams,3] B_raw_hat_def;
vector<lower=0>[ngames] theta1;
vector<lower=0>[ngames] theta2;
 for (t in 1:nteams) {
  B_raw_hat_att[t,1] = mu_raw_beta0_att;
  B_raw_hat_att[t,2] = mu_raw_beta1_att;
  B_raw_hat_att[t,3] = mu_raw_beta1_ser;
  B_raw_hat_def[t,1] = mu_raw_beta0_def;
  B_raw_hat_def[t,2] = mu_raw_beta1_def;
  B_raw_hat_def[t,3] = mu_raw_beta1_blo;
 }
beta0_att_star = xi_beta0_att*B_raw_att[,1];
beta1_att_star = xi_beta1_att*B_raw_att[,2];
beta1_ser_star = xi_beta1_ser*B_raw_att[,3];
beta0_def_star = xi_beta0_def*B_raw_def[,1];
beta1_def_star = xi_beta1_def*B_raw_def[,2];
beta1_blo_star = xi_beta1_blo*B_raw_def[,3];
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
Tau_B_raw_att ~ wishart(df, IW_att);
Tau_B_raw_def ~ wishart(df, IW_def);
mu_raw_beta0_att ~ normal(0, 100);
mu_raw_beta1_att ~ normal(0, 100);
mu_raw_beta1_ser ~ normal(0, 100);
mu_raw_beta0_def ~ normal(0, 100);
mu_raw_beta1_def ~ normal(0, 100);
mu_raw_beta1_blo ~ normal(0, 100);
  for (t in 1:nteams) {
   B_raw_att[t,] ~ multi_normal_prec(B_raw_hat_att[t,], Tau_B_raw_att);
   B_raw_def[t,] ~ multi_normal_prec(B_raw_hat_def[t,], Tau_B_raw_def);
  }
xi_beta0_att ~ uniform(0, 100);
xi_beta1_att ~ uniform(0, 100);
xi_beta1_ser ~ uniform(0, 100);
xi_beta0_def ~ uniform(0, 100);
xi_beta1_def ~ uniform(0, 100);
xi_beta1_blo ~ uniform(0, 100);
gamma ~ normal(0, 10);
delta ~ normal(0, 10);
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


