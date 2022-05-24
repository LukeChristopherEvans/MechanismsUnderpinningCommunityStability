functions{
    matrix cov_GPL2(matrix x, real sq_alpha, real sq_rho, real delta) {
        int N = dims(x)[1];
        matrix[N, N] K;
        for (i in 1:(N-1)) {
          K[i, i] = sq_alpha + delta;
          for (j in (i + 1):N) {
            K[i, j] = sq_alpha * exp(-sq_rho * square(x[i,j]) );
            K[j, i] = K[i, j];
          }
        }
        K[N, N] = sq_alpha + delta;
        return K;
    }
}

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
 vector[nData] MeanVol ; // dsep

 // intercept
 vector[nData] Intercept ;

 // random effects
  int CountryID[nData];
  int SiteID[nData];

  matrix[nUK, nUK] UKMat;
  matrix[nFin, nFin] FinMat;
  matrix[nSpan, nSpan] SpainMat;

 }
  parameters {
  // overall sigma
   real<lower=0> sigma;

  // 2 is the number of the intercept + params
//  vector[3] beta_p[nCountry];  // parameter matrix for each country [param x country]
    matrix[3,nCountry] z_p;  // non centered version of above

   vector<lower=0>[3] sigma_p; //sigma intercept and slope(s)
   vector[3] hyperbeta; // hyper prior intercept and slope(s)
// corr_matrix[2] Omega; // correlation matrix [p x p]
  cholesky_factor_corr[3] L_p; // cholesky version of above

// gaussian params
  vector<lower=0>[nCountry] etasq;
  vector<lower=0>[nCountry] rhosq;

  vector[nUK] zUK;
  vector[nSpan] zSp;
  vector[nFin] zFin;
  }

transformed parameters {
// Gauss intercepts
  vector[nUK] ukK;
  vector[nSpan] spK;
  vector[nFin] finK;

  matrix[nSpan, nSpan] SigmaMSp;
  matrix[nFin, nFin] SigmaMFi;
  matrix[nUK, nUK] SigmaMUK;

  matrix[nSpan, nSpan] LSigmaMSp;
  matrix[nFin, nFin] LSigmaMFi;
  matrix[nUK, nUK] LSigmaMUK;

  // make model matrix
  matrix[nData, 3] M;
  matrix[3,nCountry] z; // non-centered version of beta_p

  M[,1] = Intercept;
  M[,2] = OverallNiDist;
  M[,3] = MeanVol;

  SigmaMSp = cov_GPL2(SpainMat, etasq[1], rhosq[1], 0.1);
  SigmaMFi = cov_GPL2(FinMat, etasq[2], rhosq[2], 0.1);
  SigmaMUK = cov_GPL2(UKMat, etasq[3], rhosq[3], 0.1);

  LSigmaMSp= cholesky_decompose( SigmaMSp);
  spK = LSigmaMSp* zSp;
  LSigmaMFi= cholesky_decompose( SigmaMFi);
  finK = LSigmaMFi* zFin;
  LSigmaMUK= cholesky_decompose( SigmaMUK);
  ukK = LSigmaMUK* zUK;


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

// Gaussian process matrices and params
  rhosq ~ exponential( 1 );
  etasq ~ exponential( 1);

  zSp ~ normal( 0 , 1);
  zFin~ normal( 0, 1 );
  zUK~  normal(0 , 1);

 for(i in 1:nData) {
  if (CountryID[i] == 1) {
  // mu[i] = M[i] * (beta_p[CountryID[i]]) + spK[SiteID[i]]; // row of M matrix multi x beta_p column
   mu[i] =  hyperbeta[1] + z[1,CountryID[i]] + (hyperbeta[2] + z[2,CountryID[i]]) * M[i,2] +(hyperbeta[3] + z[3,CountryID[i]]) * M[i,3] +  spK[SiteID[i]];
  }
  if (CountryID[i] == 2) {
   //mu[i] = M[i] * (beta_p[CountryID[i]]) + finK[SiteID[i]];
   mu[i] =  hyperbeta[1] + z[1,CountryID[i]] + (hyperbeta[2] + z[2,CountryID[i]]) * M[i,2] +(hyperbeta[3] + z[3,CountryID[i]]) * M[i,3] +  finK[SiteID[i]];
  }
  if (CountryID[i] == 3) {
  // mu[i] = M[i] * (beta_p[CountryID[i]]) + ukK[SiteID[i]];
   mu[i] =  hyperbeta[1] + z[1,CountryID[i]] + (hyperbeta[2] + z[2,CountryID[i]]) * M[i,2] +(hyperbeta[3] + z[3,CountryID[i]]) * M[i,3] +  ukK[SiteID[i]];
  }

  }

   SppRich ~ normal(mu,sigma);

}

// calculate loglik for the WAIC
generated quantities{
    vector[nData] log_lik;
    vector[nData] mu;
    matrix[3, nCountry] beta_p;
    matrix[3, 3] Omega;
    vector[4] ADSEP;

    for(i in 1:nData) {
     if (CountryID[i] == 1) {
     // mu[i] = M[i] * (beta_p[CountryID[i]]) + spK[SiteID[i]]; // row of M matrix multi x beta_p column
      mu[i] =  hyperbeta[1] + z[1,CountryID[i]] + (hyperbeta[2] + z[2,CountryID[i]]) * M[i,2] +(hyperbeta[3] + z[3,CountryID[i]]) * M[i,3] +  spK[SiteID[i]];
     }
     if (CountryID[i] == 2) {
      //mu[i] = M[i] * (beta_p[CountryID[i]]) + finK[SiteID[i]];
      mu[i] =  hyperbeta[1] + z[1,CountryID[i]] + (hyperbeta[2] + z[2,CountryID[i]]) * M[i,2] +(hyperbeta[3] + z[3,CountryID[i]]) * M[i,3] +  finK[SiteID[i]];
     }
     if (CountryID[i] == 3) {
     // mu[i] = M[i] * (beta_p[CountryID[i]]) + ukK[SiteID[i]];
      mu[i] =  hyperbeta[1] + z[1,CountryID[i]] + (hyperbeta[2] + z[2,CountryID[i]]) * M[i,2] +(hyperbeta[3] + z[3,CountryID[i]]) * M[i,3] +  ukK[SiteID[i]];
     }

     }

    // transform omega back
    Omega = multiply_lower_tri_self_transpose(L_p);

    // backtransform my beta_p for ease
    for( i in 1:3){
       for (j in 1:3) {
        beta_p[i,j] = hyperbeta[i] + z[i,j];
       }
    }
    // get the key parameters up top
    ADSEP[1] = hyperbeta[3];
    ADSEP[2] = beta_p[3,1];
    ADSEP[3] = beta_p[3,2];
    ADSEP[4] = beta_p[3,3];
    for ( i in 1:nData ) log_lik[i] = normal_lpdf( SppRich[i] | mu[i],sigma );
}
