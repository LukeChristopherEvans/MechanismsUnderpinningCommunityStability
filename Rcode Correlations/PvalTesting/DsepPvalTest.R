

library(rstan)
library(rethinking)


bayesp =c()
freqp = c()

for (i in 1:25) {
  
y =  rnorm(100)
x = -1 + 0.5 * y + rnorm(100,0,2.5)


lmmod = 
"
data {
  
  int <lower=0> n;
  vector[n] y;
  vector[n] x;

}

parameters {

 real beta;
 real alpha;
 real<lower=0> sigma; 

}

model {
 sigma ~ exponential(1);
 beta~std_normal();
 alpha ~std_normal();
 y ~ normal(alpha + beta * x,sigma);

}

"

dtlist = list(
  y = y,
  x= x,
  n = length(x)
)


fit =stan(model_code = lmmod,data=dtlist,iter=1000)
samps = extract.samples(fit,pars="beta")

bp =2*(1- max(c(sum(samps$beta > 0) / 2000, sum(samps$beta  < 0) /2000)))
fp = summary(lm(y~x))$coefficients[,"Pr(>|t|)"][[2]]

bayesp =c(bayesp,bp)
freqp = c(freqp,fp)
}


plot(bayesp,freqp)
abline(0,1)
