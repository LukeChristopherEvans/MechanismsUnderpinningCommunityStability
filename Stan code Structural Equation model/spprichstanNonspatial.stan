
data{
 // length of data and country indexes
 int nData;
 int nCountry;
 int nUK;
 int nFin;
 int nSpan;

 // outcome var in first model
 vector[nData] SppRich ;

 // main x vars
 vector[nData] OverallNiDist ;

 // intercept
 vector[nData] Intercept ;

 // random effects
  int CountryID[nData];
  int SiteID[nData];


 }

parameters {
  // overall sigma
   real<lower=0> sigma;

  // 2 is the number of the intercept + params
//  vector[2] beta_p[nCountry];  // parameter matrix for each country [param x country]
    matrix[2,nCountry] z_p;  // non centered version of above

   vector<lower=0>[2] sigma_p; //sigma intercept and slope(s)
   vector[2] hyperbeta; // hyper prior intercept and slope(s)
// corr_matrix[2] Omega; // correlation matrix [p x p]
  cholesky_factor_corr[2] L_p; // cholesky version of above

  }

transformed parameters {
  // make model matrix
  matrix[nData, 2] M;
  matrix[2,nCountry] z; // non-centered version of beta_p

  M[,1] = Intercept;
  M[,2] = OverallNiDist;


  z = (diag_pre_multiply(sigma_p, L_p) * z_p);

}

model{
 // mean
 vector[nData] mu;

 // hyperpriors
 hyperbeta ~ normal(0,1); // hprior on intercept and slope
 sigma_p ~ exponential(1); //  hprior on intercept and slope
// Omega ~ lkj_corr(2);

 L_p ~ lkj_corr_cholesky(2);
 to_vector(z_p) ~ normal(0, 1);
 //beta_p ~ multi_normal(hyperbeta, quad_form_diag(Omega, sigma_p));

 sigma ~ exponential(1);

 for(i in 1:nData) {
  if (CountryID[i] == 1) {
  // mu[i] = M[i] * (beta_p[CountryID[i]]) + spK[SiteID[i]]; // row of M matrix multi x beta_p column
   mu[i] =  hyperbeta[1] + z[1,CountryID[i]] + (hyperbeta[2] + z[2,CountryID[i]]) * M[i,2];
  }
  if (CountryID[i] == 2) {
   //mu[i] = M[i] * (beta_p[CountryID[i]]) + finK[SiteID[i]];
   mu[i] =  hyperbeta[1] + z[1,CountryID[i]] + (hyperbeta[2] + z[2,CountryID[i]]) * M[i,2];
  }
  if (CountryID[i] == 3) {
  // mu[i] = M[i] * (beta_p[CountryID[i]]) + ukK[SiteID[i]];
   mu[i] =  hyperbeta[1] + z[1,CountryID[i]] + (hyperbeta[2] + z[2,CountryID[i]]) * M[i,2];
  }

  }

   SppRich ~ normal(mu,sigma);

}

// calculate loglik for the WAIC
generated quantities {

    vector[nData] log_lik;
    vector[nData] mu;
    matrix[2, nCountry] beta_p;
    matrix[2, 2] Omega;

    for ( i in 1:nData ) {
     if (CountryID[i] == 1) {
    //  mu[i] =  M[i] * (beta_p[CountryID[i]]) + spK[SiteID[i]];
      mu[i] =  hyperbeta[1] + z[1,CountryID[i]] + (hyperbeta[2] + z[2,CountryID[i]]) * M[i,2];
     }
     if (CountryID[i] == 2) {
    //  mu[i] = M[i] * (beta_p[CountryID[i]])  + finK[SiteID[i]];
       mu[i] =  hyperbeta[1] + z[1,CountryID[i]] + (hyperbeta[2] + z[2,CountryID[i]]) * M[i,2];
     }
     if (CountryID[i] == 3) {
    //  mu[i] =  M[i] * (beta_p[CountryID[i]]) +  ukK[SiteID[i]];
       mu[i] =  hyperbeta[1] + z[1,CountryID[i]] + (hyperbeta[2] + z[2,CountryID[i]]) * M[i,2];
     }
     }

    // transform omega back
    Omega = multiply_lower_tri_self_transpose(L_p);

    // backtransform my beta_p for ease
    for( i in 1:2){
       for (j in 1:3) {
        beta_p[i,j] = hyperbeta[i] + z[i,j];
       }
    }
    for ( i in 1:nData ) log_lik[i] = normal_lpdf( SppRich[i] | mu[i],sigma );
}
