
data{
 // length of data and country indexes
 int nData;
 int nCountry;
 int nUK;
 int nFin;
 int nSpan;

 // outcome var in model
  vector[nData] ICOV ;

 // main x var
  vector[nData] SpeICOV;
  vector[nData] SynIndex;
  vector[nData] LgMeanCom;
  // intercept
 vector[nData] Intercept ;

 // random effects
  int CountryID[nData];
  int SiteID[nData];

////
 }

parameters {
  // overall sigma
  real<lower=0> sigma;

  // 2 is the number of the intercept + params
  matrix[4,nCountry] z_p;  // non centered version of above

  vector<lower=0>[4] sigma_p; //sigma intercept and slope(s)
  vector[4] hyperbeta; // hyper prior intercept and slope(s)

  cholesky_factor_corr[4] L_p; // cholesky version of above

  }

transformed parameters {
  // make model matrix
  matrix[nData, 4] M;
  matrix[4,nCountry] z; // non-centered version of beta_p

  M[,1] = Intercept;
  M[,2] = SpeICOV;
  M[,3] = SynIndex;
  M[,4] = LgMeanCom;

  z = (diag_pre_multiply(sigma_p, L_p) * z_p);

}

model{
 // mean
 vector[nData] mu;

 // hyperpriors
 hyperbeta ~ normal(0,1); // hprior on intercept and slope
 sigma_p ~ exponential(1); //  hprior on intercept and slope

 L_p ~ lkj_corr_cholesky(2);
 to_vector(z_p) ~ normal(0, 1);

 sigma ~ exponential(1);

 for(i in 1:nData) {
  if (CountryID[i] == 1) {
   mu[i] =  hyperbeta[1] + z[1,CountryID[i]] + (hyperbeta[2] + z[2,CountryID[i]]) * M[i,2] +  (hyperbeta[3] + z[3,CountryID[i]]) * M[i,3]+ (hyperbeta[4] + z[4,CountryID[i]]) * M[i,4] ;
  }
  if (CountryID[i] == 2) {
   mu[i] =  hyperbeta[1] + z[1,CountryID[i]] + (hyperbeta[2] + z[2,CountryID[i]]) * M[i,2] +  (hyperbeta[3] + z[3,CountryID[i]]) * M[i,3]+ (hyperbeta[4] + z[4,CountryID[i]]) * M[i,4] ;
  }
  if (CountryID[i] == 3) {
     mu[i] =  hyperbeta[1] + z[1,CountryID[i]] + (hyperbeta[2] + z[2,CountryID[i]]) * M[i,2] +  (hyperbeta[3] + z[3,CountryID[i]]) * M[i,3]+ (hyperbeta[4] + z[4,CountryID[i]]) * M[i,4] ;
  }
  }

   ICOV ~ normal(mu,sigma);


}
//
// calculate loglik for the WAIC
generated quantities{
   vector[nData] log_lik;
   vector[nData] mu;
   matrix[4, nCountry] beta_p;
   matrix[4, 4] Omega;
   for(i in 1:nData) {
    if (CountryID[i] == 1) {
       mu[i] =  hyperbeta[1] + z[1,CountryID[i]] + (hyperbeta[2] + z[2,CountryID[i]]) * M[i,2] +  (hyperbeta[3] + z[3,CountryID[i]]) * M[i,3]+ (hyperbeta[4] + z[4,CountryID[i]]) * M[i,4];
    }
    if (CountryID[i] == 2) {
       mu[i] =  hyperbeta[1] + z[1,CountryID[i]] + (hyperbeta[2] + z[2,CountryID[i]]) * M[i,2] +  (hyperbeta[3] + z[3,CountryID[i]]) * M[i,3]+ (hyperbeta[4] + z[4,CountryID[i]]) * M[i,4] ;
    }
    if (CountryID[i] == 3) {
       mu[i] =  hyperbeta[1] + z[1,CountryID[i]] + (hyperbeta[2] + z[2,CountryID[i]]) * M[i,2] +  (hyperbeta[3] + z[3,CountryID[i]]) * M[i,3]+ (hyperbeta[4] + z[4,CountryID[i]]) * M[i,4];
    }
    }

    // transform omega back
    Omega = multiply_lower_tri_self_transpose(L_p);

    // backtransform my beta_p for ease
    for(i in 1:4){
       for (j in 1:3) {
        beta_p[i,j] = hyperbeta[i] + z[i,j];
            }
        }

    for (i in 1:nData ) log_lik[i] = normal_lpdf( ICOV[i] | mu[i],sigma );
}
