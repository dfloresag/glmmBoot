// Regression model with random slope and intercept.
  #include <TMB.hpp>
  #include <cmath>
  template<class Type>
  Type objective_function<Type>::operator() ()
  {
  DATA_FACTOR(group);
  DATA_VECTOR(y);
  DATA_MATRIX(X);
  DATA_VECTOR(w);         
  DATA_VECTOR(t);      
  
  DATA_INTEGER(family);
  PARAMETER_MATRIX(u);     
  PARAMETER_VECTOR(betas);   
  PARAMETER_VECTOR(log_sigma);    
  PARAMETER(tan_rho);
  
  using namespace density;
  int nobs=y.size();
  int ngroups=NLEVELS(group);

  Type res(0.0);
  Type P_I(3.141592653589793238463);
  
  vector<Type> sigma(exp(log_sigma));
  vector<Type> sigma2(exp(2*log_sigma));
  Type rho((2/P_I)*atan(tan_rho));
  
  // Type rho()
  ADREPORT(sigma);
  ADREPORT(sigma2);
  ADREPORT(rho);
  
  matrix<Type> Sigma(2,2);
  Sigma.fill(rho);             
  Sigma.diagonal() *= 1.0/rho;

  /* Prior: intercept~N(mu,sd) */
  for(int i=0;i<ngroups;i++){
    res += w[i]*MVNORM(Sigma)(vector<Type>(u.row(i)));
  }

  vector<Type> xtb = X*betas;
  vector<Type> eta(nobs);
  
  for(int ij=0;ij<nobs;ij++){
    int i=group[ij];
    eta[ij] = xtb[ij]+exp(log_sigma[0])*u(i,0)+t[ij]*exp(log_sigma[1])*u(i,1);
    if(family==0){
      res-= w[i]*dbinom(y[ij],Type(1),Type(1/(1+exp(-eta[ij]))),1);
    } else if (family==1){
      res -= w[i]*dpois(y[ij], Type(exp(eta[ij])),1); 
    } 
    // else if (family==3){
    //   res -= w[i]*dnorm(y[ij], Type(eta[ij]), Type phi ,1); 
    // }
      
  }
  return res;
  }

