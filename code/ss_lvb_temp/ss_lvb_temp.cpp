#include <TMB.hpp>
#include <iostream>

#define see(object) std::cout << #object ":\n" << object << "\n";

//same as ss_temp_growth_v5 except tidied up.

template<class Type>
Type dbinorm(vector<Type> x, vector<Type> mu, matrix<Type> Sigma, int do_log)
{
  int dim = 2;
	Type cond_mu;
	Type cond_var;
	Type d = dnorm(x(0), mu(0), exp(0.5*log(Sigma(0,0))),1);
  cond_mu = mu(1) + (x(0) - mu(0))*Sigma(1,0)/Sigma(0,0);
  cond_var = Sigma(1,1) - Sigma(1,0)*Sigma(1,0)/Sigma(0,0);
  d += dnorm(x(1), cond_mu, exp(0.5*log(cond_var)),1);
	if(do_log == 1) return(d);
	else return(exp(d));
}

template<class Type>
vector<Type> rbinorm(vector<Type> mu, matrix<Type> Sigma)
{
  int dim = 2;
	Type cond_mu;
	Type cond_var;
	vector<Type> x(dim);
	x(0) = rnorm(mu(0), exp(0.5*log(Sigma(0,0))));
  cond_mu = mu(1) + (x(0) - mu(0))*Sigma(1,0)/Sigma(0,0);
  cond_var = Sigma(1,1) - Sigma(1,0)*Sigma(1,0)/Sigma(0,0);
  x(1) = rnorm(cond_mu, exp(0.5*log(cond_var)));
	return(x);
}


template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(weight); //number mature
  DATA_IVECTOR(iswt); // weight present
  DATA_VECTOR(len); //number of mature + not mature
  DATA_IVECTOR(islen); // length present
  DATA_IVECTOR(age_obs); //age for each observation
  DATA_IVECTOR(year_obs); //age for each observation
  DATA_IVECTOR(cy_obs); //cohort year for each observation (year of birth)
  DATA_IVECTOR(Ecov_maxages_b); //last year of Ecov to use (nobs)
  DATA_IVECTOR(Ecov_maxages_a); //last year of Ecov to use (nobs)
  DATA_IVECTOR(Ecov_maxages_Linf); //last year of Ecov to use (nobs)
  DATA_MATRIX(Ecov_obs);
  DATA_IMATRIX(use_Ecov_obs);
  DATA_MATRIX(Ecov_obs_sigma);
  DATA_INTEGER(fit_growth);
  DATA_INTEGER(fit_Ecov);
  DATA_INTEGER(fit_Ecov_obs);
  DATA_INTEGER(fit_k);
  DATA_INTEGER(laa_matrix_max_age); //max age for reporting predicted length and weight at age over time
  DATA_MATRIX(X_Linf) //design matrix for Linf (nobs x ncov)
  DATA_INTEGER(obs_for_Linf_y) // which row of X_Linf to use for yearly predicted Linf
    
  PARAMETER(beta_b);
  PARAMETER(beta_a);
  PARAMETER_VECTOR(beta_k);
  PARAMETER_VECTOR(beta_Linf);
  PARAMETER(t0);
  PARAMETER_MATRIX(beta_Ecov_b);
  PARAMETER_MATRIX(beta_Ecov_a);
  PARAMETER_MATRIX(beta_Ecov_k);
  PARAMETER_MATRIX(beta_Ecov_Linf);
  PARAMETER_VECTOR(Ecov_mu);
  PARAMETER_MATRIX(Ecov_AR_pars); //AR1 process (2)
  PARAMETER_MATRIX(Ecov_re);
  PARAMETER_VECTOR(k_AR_pars); //AR1 process (2)
  PARAMETER_VECTOR(k_re);
  PARAMETER_VECTOR(log_Ecov_obs_sig_scale);
  PARAMETER(beta_sig_Linf);
  PARAMETER(beta_sig_L_obs);
  PARAMETER(beta_sig_W_obs);

  using namespace density;

  int n_obs = weight.size();
  int n_years_Ecov = Ecov_obs.rows();
  int n_Ecov = Ecov_obs.cols();
  Type zero = Type(0);
  Type one = Type(1);
  Type nll= zero; //negative log-likelihood
  matrix<Type> Ecov_y(n_years_Ecov,n_Ecov);
  vector<Type> nll_Ecov(n_Ecov);
  nll_Ecov.setZero();
  vector<Type> Ecov_phi(n_Ecov), Ecov_sig(n_Ecov);
  for(int p = 0; p < n_Ecov; p++) 
  {
    Ecov_phi(p) = -1 + 2/(1 + exp(-Ecov_AR_pars(0,p)));
    Ecov_sig(p) = exp(Ecov_AR_pars(1,p));
    for(int y = 0; y < n_years_Ecov; y++) Ecov_y(y,p) = Ecov_mu(p) + Ecov_re(y,p);
    nll_Ecov(p) -= dnorm(Ecov_re(0,p), zero, Ecov_sig(p)*exp(-Type(0.5) * log(one - pow(Ecov_phi(p),Type(2)))), 1);
  //see(nll_Ecov);
    SIMULATE
    {
      Ecov_re(0,p) = rnorm(zero, Ecov_sig(p)*exp(-Type(0.5) * log(one - pow(Ecov_phi(p),Type(2)))));
      Ecov_y(0,p) = Ecov_mu(p) + Ecov_re(0,p);
    }
    for(int y = 1; y < n_years_Ecov; y++) 
    {
      nll_Ecov(p) -= dnorm(Ecov_re(y,p), Ecov_phi(p) * Ecov_re(y-1,p), Ecov_sig(p), 1);
      SIMULATE
      {
        Ecov_re(y,p) = rnorm(Ecov_phi(p) * Ecov_re(y-1,p), Ecov_sig(p));
        Ecov_y(y,p) = Ecov_mu(p) + Ecov_re(y,p);
      }
    }
  }  
  if(fit_Ecov == 1) nll += nll_Ecov.sum();
  SIMULATE REPORT(Ecov_re);
  
  vector<Type> nll_Ecov_obs(n_Ecov);
  nll_Ecov_obs.setZero();
  for(int p = 0; p < n_Ecov; p++) 
  {
    for(int y = 0; y < n_years_Ecov; y++)
    {
      if(use_Ecov_obs(y,p) == 1) 
      {
        nll_Ecov_obs(p) -= dnorm(Ecov_obs(y,p), Ecov_y(y,p), Ecov_obs_sigma(y,p)*exp(log_Ecov_obs_sig_scale(p)), 1);
        SIMULATE Ecov_obs(y,p) = rnorm(Ecov_y(y,p), Ecov_obs_sigma(y,p)*exp(log_Ecov_obs_sig_scale(p)));
      }
    }
  }
  if(fit_Ecov_obs == 1) nll += nll_Ecov_obs.sum();
  //see(nll_Ecov_obs);
  SIMULATE REPORT(Ecov_obs);
  //k AR(1) process
  Type k_phi;
  Type k_sig;
  Type nll_k = zero;
  k_phi = -one + Type(2)/(one + exp(-k_AR_pars(0)));
  k_sig = exp(k_AR_pars(1));
  nll_k -= dnorm(k_re(0), zero, k_sig*exp(-Type(0.5) * log(one - pow(k_phi,Type(2)))), 1);
  for(int y = 1; y < n_years_Ecov; y++) nll_k -= dnorm(k_re(y), k_phi * k_re(y-1), k_sig, 1);
  if(fit_k == 1) 
  {
    nll += nll_k;
    //see(nll_k);
    SIMULATE
    {
      k_re(0) = rnorm(zero, k_sig*exp(-Type(0.5) * log(one - pow(k_phi,Type(2)))));
      for(int y = 1; y < n_years_Ecov; y++) 
      {
        k_re(y) = rnorm(k_phi * k_re(y-1), k_sig);
      }
      REPORT(k_re);
    }
  }
  
  vector<Type> l_b(n_obs), l_a(n_obs), log_l_pred(n_obs), log_w_pred(n_obs), ll(n_obs), v_log_w(n_obs), k_sum(n_obs);
  k_sum.setZero(); ll.setZero();
  vector <Type> l_Linf = X_Linf * beta_Linf;
  Type v_log_l = exp(Type(2) * beta_sig_L_obs) + exp(Type(2) * beta_sig_Linf);
  for(int i = 0; i < n_obs; i++) 
  {
    l_b(i) = beta_b;
    l_a(i) = beta_a;
    //for(int x = 0; x < X_Linf.cols(); x++) l_Linf(i) += beta_Linf(x) * X_Linf(i,x);
    for(int p = 0; p < n_Ecov; p++)
    {
      for(int y = 0; y <= Ecov_maxages_b(i); y++)       l_b(i) += Ecov_y(year_obs(i) - 1 - age_obs(i) + y,p) * beta_Ecov_b(y,p);
      for(int y = 0; y <= Ecov_maxages_a(i); y++)       l_a(i) += Ecov_y(year_obs(i) - 1 - age_obs(i) + y,p) * beta_Ecov_a(y,p);
      for(int y = 0; y <= Ecov_maxages_Linf(i); y++) l_Linf(i) += Ecov_y(year_obs(i) - 1 - age_obs(i) + y,p) * beta_Ecov_Linf(y,p);
    }
    //do LVB k something like Millar et al. 1999. Need to specify length of beta_Ecov_k = maximum observed age and specify map to fix older ages at 0
    Type log_k0 = beta_k(0) + k_re(year_obs(i) -1 - age_obs(i) + 0);
    for(int p = 0; p < n_Ecov; p++) log_k0 += Ecov_y(year_obs(i) - 1 - age_obs(i) + 0,p) * beta_Ecov_k(0,p);
    //k0 += exp(beta_k(0) + k_re(year_obs(i) -1 - age_obs(i) + 0) + Ecov_y(year_obs(i) - 1 - age_obs(i) + 0,p) * beta_Ecov_k(0,p));
    Type k0 = exp(log_k0);
    k_sum(i) = k0;
    Type log_ky = 0.0;
    for(int y = 1; y < age_obs(i); y++) 
    {
      Type log_ky = beta_k(y) + k_re(year_obs(i) -1 - age_obs(i) + y);
      for(int p = 0; p < n_Ecov; p++) log_ky += Ecov_y(year_obs(i) - 1 - age_obs(i) + y,p) * beta_Ecov_k(y,p);
      k_sum(i) += exp(log_ky);
    }
    //k_sum(i) += exp(beta_k(y) + k_re(year_obs(i) -1 - age_obs(i) + y) + Ecov_y(year_obs(i) - 1 - age_obs(i) + y,p) * beta_Ecov_k(y,p));
    log_l_pred(i) = l_Linf(i) + log(Type(1) - exp(t0*k0 - k_sum(i)));
    log_w_pred(i) = l_a(i) + exp(l_b(i)) * log_l_pred(i);
    v_log_w(i) = exp(Type(2) * beta_sig_W_obs) + exp(Type(2) * (l_b(i) + beta_sig_Linf));
    matrix<Type> S(2,2);
    S(0,0) = v_log_l; S(1,1) = v_log_w(i);
    S(0,1) = S(1,0) = exp(l_b(i) + Type(2) * beta_sig_Linf);
    vector<Type> M(2); M(0) = log_l_pred(i); M(1) = log_w_pred(i); 
    vector<Type> obs(2); obs(0) = log(len(i)); obs(1) = log(weight(i));
    if(islen(i) == 1 & iswt(i) == 1) ll(i) += dbinorm(obs,M,S,1);
    else
    {
      if(islen(i) == 1) ll(i) += dnorm(obs(0), M(0), exp(Type(0.5)*log(S(0,0))), 1);
      if(iswt(i) == 1) ll(i) += dnorm(obs(1), M(1), exp(Type(0.5)*log(S(1,1))), 1);
    }
    SIMULATE
    {
      obs = rbinorm(M,S);
      len(i) = exp(obs(0)); weight(i) = exp(obs(1));
    }
  }
  //see(sum(ll));
  if(fit_growth == 1) nll -= sum(ll);
  SIMULATE 
  {
    REPORT(len);
    REPORT(weight);
  }
  matrix<Type> log_a(n_years_Ecov,beta_Ecov_a.rows()), log_b(n_years_Ecov,beta_Ecov_b.rows());
  log_a.setZero();
  for(int i = 0; i < beta_Ecov_a.rows(); i++) //age_obs
  {
    for(int y = i+1; y < n_years_Ecov; y++) //year_obs
    {
      log_a(y,i) = beta_a;
      for(int p = 0; p < n_Ecov; p++) for(int j = 0; j <= i; j++) log_a(y,i) += beta_Ecov_a(j,p) * Ecov_y(y - 1 - i + j,p);
    }
  }
  log_b.setZero();
  for(int i = 0; i < beta_Ecov_b.rows(); i++) //age_obs
  {
    for(int y = i+1; y < n_years_Ecov; y++) //year_obs
    {
      log_b(y,i) = beta_b;
      for(int p = 0; p < n_Ecov; p++) for(int j = 0; j <= i; j++) log_b(y,i) += beta_Ecov_b(j,p) * Ecov_y(y - 1 - i + j,p);
    }
  }
  matrix<Type> log_k_ya(n_years_Ecov,beta_Ecov_k.rows());
  log_k_ya.setZero();
  for(int y = 1; y < n_years_Ecov; y++) //start in second year because 
  {
    for(int a = 0; a < y; a++) if(a < beta_Ecov_k.rows()) //age
    {
      log_k_ya(y,a) += beta_k(a) + k_re(y - 1);
      for(int p = 0; p < n_Ecov; p++) log_k_ya(y,a) += beta_Ecov_k(a,p) * Ecov_y(y - 1,p);
    }
  }
  matrix <Type> log_Linf(n_years_Ecov,beta_Ecov_Linf.rows());
  log_Linf.setZero();
  for(int i = 0; i < beta_Ecov_Linf.rows(); i++) //age_obs
  {
    for(int y = i+1; y < n_years_Ecov; y++) //year_obs
    {
      for(int x = 0; x < X_Linf.cols(); x++) log_Linf(y,i) += beta_Linf(x) * X_Linf(obs_for_Linf_y-1,x);
      for(int p = 0; p < n_Ecov; p++) for(int j = 0; j <= i; j++) log_Linf(y,i) += beta_Ecov_Linf(j,p) * Ecov_y(y - 1 - i + j,p);
    }
  }
  matrix<Type> log_a_Ecov(n_years_Ecov,beta_Ecov_a.rows()), log_b_Ecov(n_years_Ecov,beta_Ecov_b.rows());
  matrix<Type> log_k_Ecov(n_years_Ecov,beta_Ecov_k.rows()), log_Linf_Ecov(n_years_Ecov,beta_Ecov_Linf.rows());
  log_b_Ecov.fill(beta_b);
  log_a_Ecov.fill(beta_a);
  log_k_Ecov.setZero();
  Type temp_Linf_pred = 0.0;
  for(int x = 0; x < X_Linf.cols(); x++) temp_Linf_pred += beta_Linf(x) * X_Linf(obs_for_Linf_y-1,x);
  log_Linf_Ecov.fill(temp_Linf_pred);
  
  for(int y = 0; y < n_years_Ecov; y++)  
  {
    for(int p = 0; p < n_Ecov; p++)
    {
      for(int i = 0; i < beta_Ecov_b.rows(); i++) for(int j = 0; j <= i; j++) log_b_Ecov(y,i) += beta_Ecov_b(j,p) * Ecov_y(y,p);
      for(int i = 0; i < beta_Ecov_a.rows(); i++) for(int j = 0; j <= i; j++) log_a_Ecov(y,i) += beta_Ecov_a(j,p) * Ecov_y(y,p);
      for(int i = 0; i < beta_Ecov_Linf.rows(); i++) for(int j = 0; j <= i; j++) log_Linf_Ecov(y,i) += beta_Ecov_Linf(j,p) * Ecov_y(y,p);
    }
    for(int i = 0; i < beta_Ecov_k.rows(); i++)
    {
      log_k_Ecov(y,i) += beta_k(i);
      for(int p = 0; p < n_Ecov; p++) log_k_Ecov(y,i) += beta_Ecov_k(i,p) * Ecov_y(y,p);
    }
  }
  
  //get weight at age
  int max_age = laa_matrix_max_age;
  matrix<Type> k_age(n_years_Ecov,max_age), log_Linf_age(n_years_Ecov,max_age), log_a_age(n_years_Ecov,max_age), log_b_age(n_years_Ecov,max_age);
  matrix<Type> log_laa(n_years_Ecov,max_age), log_waa(n_years_Ecov,max_age);
  k_age.setZero(); log_Linf_age.setZero(); log_a_age.setZero(); log_b_age.setZero(); log_laa.setZero(); log_waa.setZero();
  log_Linf_age.fill(temp_Linf_pred);
  for(int y = 1; y < n_years_Ecov; y++) for(int a = 0; a < max_age; a++) if(y-1>=a)
  {
    log_a_age(y,a) = beta_a;
    log_b_age(y,a) = beta_b;
    int maxage_a = a+1;
    if(a+1 > beta_Ecov_a.rows()) maxage_a = beta_Ecov_a.rows();
    for(int i = 0; i < maxage_a; i++) for(int p = 0; p < n_Ecov; p++) log_a_age(y,a) += beta_Ecov_a(i,p) * Ecov_y(y - 1 - a + i,p);
    int maxage_b = a+1;
    if(a+1 > beta_Ecov_b.rows()) maxage_b = beta_Ecov_b.rows();
    for(int i = 0; i < maxage_b; i++) for(int p = 0; p < n_Ecov; p++) log_b_age(y,a) += beta_Ecov_b(i,p) * Ecov_y(y - 1 - a + i,p);
    int maxage_Linf = a+1;
    if(a+1 > beta_Ecov_Linf.rows()) maxage_Linf = beta_Ecov_Linf.rows();
    for(int i = 0; i < maxage_Linf; i++) for(int p = 0; p < n_Ecov; p++) log_Linf_age(y,a) += beta_Ecov_Linf(i,p) * Ecov_y(y - 1 - a + i,p);
    Type log_k0 = beta_k(0) + k_re(y -1 - a + 0);
    for(int p = 0; p < n_Ecov; p++) log_k0 += Ecov_y(y - 1 - a + 0,p) * beta_Ecov_k(0,p);
    k_age(y,a) = exp(log_k0);
    for(int i = 1; i < a + 1; i++) 
    {
      Type log_kya = beta_k(i) + k_re(y -1 - a + i);
      for(int p = 0; p < n_Ecov; p++) log_kya += beta_Ecov_k(i,p) * Ecov_y(y - 1 - a + i,p);
      k_age(y,a) += exp(log_kya); 
    }
    log_laa(y,a) = log_Linf_age(y,a) + log(one - exp(t0 * exp(log_k0) - k_age(y,a)));
    log_waa(y,a) = log_a_age(y,a) + exp(log_b_age(y,a)) * log_laa(y,a);
  }

  matrix<Type> Ecov_dev = Ecov_y - Ecov_obs;
  ADREPORT(Ecov_dev);
  ADREPORT(Ecov_y);
  ADREPORT(Ecov_phi);
  ADREPORT(Ecov_sig);
  ADREPORT(k_phi);
  ADREPORT(k_sig);
  ADREPORT(log_a);
  ADREPORT(log_b);
  ADREPORT(log_k_ya);
  ADREPORT(log_Linf);
  ADREPORT(log_laa);
  ADREPORT(log_waa);
  REPORT(Ecov_y);
  REPORT(log_a);
  REPORT(log_b);
  REPORT(log_k_ya);
  REPORT(log_Linf);
  REPORT(log_a_Ecov);
  REPORT(log_b_Ecov);
  REPORT(log_k_Ecov);
  REPORT(log_Linf_Ecov);
  REPORT(log_l_pred);
  REPORT(log_w_pred);
  REPORT(nll_Ecov);
  REPORT(nll_Ecov_obs);
  REPORT(log_laa);
  REPORT(log_waa);
  REPORT(log_Linf_age);
  REPORT(k_age);
  REPORT(log_a_age);
  REPORT(log_b_age);
  REPORT(ll);
  return nll;
}

