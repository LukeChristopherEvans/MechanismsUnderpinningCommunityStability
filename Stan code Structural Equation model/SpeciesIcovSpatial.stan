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

///
data{
 // length of data and country indexes
 int nData;
 int nCountry;
 int nUK;
 int nFin;
 int nSpan;

 // outcome var in model
  vector[nData] SpeICOV ;

 // main x var
  vector[nData] Mismatch ;
  vector[nData] MeanVol ;
  vector[nData] LgMeanAb ;
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
  matrix[4,nCountry] z_p;  // parameter matrix for each country [param x country]

   vector<lower=0>[4] sigma_p; //sigma intercept and slope(s)
   vector[4] hyperbeta; // hyper prior intercept and slope(s)

   cholesky_factor_corr[4] L_p; // cholesky version of above

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
  matrix[nData, 4] M;
  matrix[4,nCountry] z; // non-centered version of beta_p

  M[,1] = Intercept;
  M[,2] = Mismatch;
  M[,3] = MeanVol;
  M[,4] = LgMeanAb;

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

 L_p ~ lkj_corr_cholesky(2);
 to_vector(z_p) ~ normal(0, 1);

 sigma ~ exponential(1);

// Gaussian process matrices and params
 rhosq ~ exponential(1);
 etasq ~ exponential(1);

 zSp ~ normal(0,1);
 zFin~ normal(0,1);
 zUK~  normal(0,1);


  for(i in 1:nData) {
   if (CountryID[i] == 1) {
    mu[i] =  hyperbeta[1] + z[1,CountryID[i]] + (hyperbeta[2] + z[2,CountryID[i]]) * M[i,2] +  (hyperbeta[3] + z[3,CountryID[i]]) * M[i,3]+ (hyperbeta[4] + z[4,CountryID[i]]) * M[i,4] +  spK[SiteID[i]];
   }
   if (CountryID[i] == 2) {
    mu[i] =  hyperbeta[1] + z[1,CountryID[i]] + (hyperbeta[2] + z[2,CountryID[i]]) * M[i,2] +  (hyperbeta[3] + z[3,CountryID[i]]) * M[i,3]+ (hyperbeta[4] + z[4,CountryID[i]]) * M[i,4] + finK[SiteID[i]];
   }
   if (CountryID[i] == 3) {
      mu[i] =  hyperbeta[1] + z[1,CountryID[i]] + (hyperbeta[2] + z[2,CountryID[i]]) * M[i,2] +  (hyperbeta[3] + z[3,CountryID[i]]) * M[i,3]+ (hyperbeta[4] + z[4,CountryID[i]]) * M[i,4] +  ukK[SiteID[i]];
   }
   }

   SpeICOV ~ normal(mu,sigma);

}

// calculate loglik for the WAIC
generated quantities{
    vector[nData] log_lik;
    vector[nData] mu;
    matrix[4, nCountry] beta_p;
    matrix[4, 4] Omega;
    for(i in 1:nData) {
     if (CountryID[i] == 1) {
        mu[i] =  hyperbeta[1] + z[1,CountryID[i]] + (hyperbeta[2] + z[2,CountryID[i]]) * M[i,2] +  (hyperbeta[3] + z[3,CountryID[i]]) * M[i,3]+ (hyperbeta[4] + z[4,CountryID[i]]) * M[i,4] +  spK[SiteID[i]];
     }
     if (CountryID[i] == 2) {
        mu[i] =  hyperbeta[1] + z[1,CountryID[i]] + (hyperbeta[2] + z[2,CountryID[i]]) * M[i,2] +  (hyperbeta[3] + z[3,CountryID[i]]) * M[i,3]+ (hyperbeta[4] + z[4,CountryID[i]]) * M[i,4] + finK[SiteID[i]];
     }
     if (CountryID[i] == 3) {
        mu[i] =  hyperbeta[1] + z[1,CountryID[i]] + (hyperbeta[2] + z[2,CountryID[i]]) * M[i,2] +  (hyperbeta[3] + z[3,CountryID[i]]) * M[i,3]+ (hyperbeta[4] + z[4,CountryID[i]]) * M[i,4] +  ukK[SiteID[i]];
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

    for (i in 1:nData ) log_lik[i] = normal_lpdf(SpeICOV[i] | mu[i],sigma );
}
