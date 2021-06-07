// This was edited from ssm_env_ssb_v4.cpp fro GB cod paper. Fits code/ss_lvb_temp simultaneously with a state-space age-structured model.
// Uses predicted weight at age in estimating SSB and reference points, BUT NOT FOR CATCH AND SURVEYS.

#include <TMB.hpp>
#include <iostream>
#include "helper_functions.hpp"


template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(n_years);
  DATA_INTEGER(n_ages);
  DATA_INTEGER(n_fleets);
  DATA_INTEGER(n_indices);
  DATA_INTEGER(n_ages_pop); //greater than n_ages to better model plus group (n_ages) in observations
  DATA_INTEGER(n_selblocks);
  DATA_IVECTOR(selblock_models);
  DATA_IMATRIX(selblock_pointer_fleets);
  DATA_IMATRIX(selblock_pointer_indices);
  DATA_IVECTOR(age_comp_model_fleets);
  DATA_IVECTOR(n_age_comp_pars_fleets);
  DATA_IVECTOR(age_comp_model_indices);
  DATA_IVECTOR(n_age_comp_pars_indices);
  DATA_VECTOR(fracyr_SSB);
  DATA_MATRIX(mature);
  DATA_IVECTOR(waa_pointer_fleets);
  DATA_INTEGER(waa_pointer_totcatch);
  DATA_IVECTOR(waa_pointer_indices);
  DATA_INTEGER(waa_pointer_ssb);
  DATA_INTEGER(waa_pointer_jan1);
  DATA_ARRAY(waa);
  DATA_MATRIX(Maa);
  DATA_MATRIX(agg_catch);
  DATA_MATRIX(agg_catch_sigma);
  DATA_MATRIX(catch_paa);
  DATA_IMATRIX(use_catch_paa);
  DATA_MATRIX(catch_Neff);
  DATA_IMATRIX(catch_aref);
  DATA_IVECTOR(units_indices);
  DATA_MATRIX(fracyr_indices);
  DATA_MATRIX(agg_indices);
  DATA_IMATRIX(use_indices);
  DATA_MATRIX(agg_index_sigma);
  DATA_IVECTOR(units_index_paa);
  DATA_MATRIX(index_paa);
  DATA_IMATRIX(use_index_paa);
  DATA_MATRIX(index_Neff);
  DATA_IMATRIX(index_aref);
  DATA_VECTOR(q_lower);
  DATA_VECTOR(q_upper);
  DATA_INTEGER(n_estimated_selpars);
  DATA_INTEGER(n_other_selpars);
  DATA_VECTOR(other_selpars);
  DATA_IVECTOR(estimated_selpar_pointers);
  DATA_IVECTOR(other_selpar_pointers);
  DATA_VECTOR(selpars_lower);
  DATA_VECTOR(selpars_upper);
  DATA_INTEGER(n_NAA_sigma);
  DATA_IVECTOR(NAA_sigma_pointers);
  DATA_INTEGER(recruit_model);
  //DATA_INTEGER(use_mat_model);
  DATA_INTEGER(use_growth_model);
  DATA_SCALAR(percentSPR);

  //DATA_VECTOR(Y); //number mature
  //DATA_VECTOR(N); //number of mature + not mature
  //DATA_IVECTOR(age_obs); //age for each observation
  //DATA_IVECTOR(year_obs); //age for each observation
  //DATA_IVECTOR(Ecov_maxages_k); //last year of Ecov to use (nobs)
  //DATA_IVECTOR(Ecov_maxages_a50); //last year of Ecov to use (nobs)
  DATA_MATRIX(Ecov_obs);
  DATA_IMATRIX(use_Ecov_obs);
  DATA_MATRIX(Ecov_obs_sigma);
  //DATA_INTEGER(binomial);
  
  DATA_IVECTOR(model_years); //which of the years in the environmental time series (Which may be longer) are the years of the assessment model

  DATA_IVECTOR(age_obs_growth); //age for each observation
  DATA_IVECTOR(year_obs_growth); //age for each observation
  DATA_VECTOR(weight); //number mature
  DATA_IVECTOR(iswt); // weight present
  DATA_VECTOR(len); //number of mature + not mature
  DATA_IVECTOR(islen); // length present
  DATA_IVECTOR(Ecov_maxages_b); //last year of Ecov to use (nobs)
  DATA_IVECTOR(Ecov_maxages_a); //last year of Ecov to use (nobs)
  //DATA_IVECTOR(Ecov_maxages_k); //last year of Ecov to use (nobs)
  DATA_IVECTOR(Ecov_maxages_Linf); //last year of Ecov to use (nobs)
  DATA_INTEGER(fit_k_LVB);
  //DATA_INTEGER(laa_matrix_max_age); //max age for reporting predicted length and weight at age over time
  DATA_MATRIX(X_Linf); //design matrix for Linf (nobs x ncov)
  DATA_INTEGER(obs_for_Linf_y); // which row of X_Linf to use for yearly predicted Linf
  DATA_INTEGER(do_growth);

  PARAMETER_VECTOR(mean_rec_pars);
  PARAMETER_VECTOR(logit_q);
  PARAMETER_VECTOR(log_F1);
  PARAMETER_MATRIX(F_devs);
  PARAMETER_VECTOR(log_N1);
  PARAMETER_VECTOR(log_NAA_sigma);
  PARAMETER_VECTOR(estimated_selpars);
  PARAMETER_VECTOR(catch_paa_pars);
  PARAMETER_VECTOR(index_paa_pars);
  PARAMETER_MATRIX(log_NAA);
  PARAMETER_VECTOR(N1_re); //random walk re for initial NAA
  PARAMETER(log_sig_N1); //var par for random walk re for initial NAA

  //maturity parameters
 /* PARAMETER(beta_k);
  PARAMETER(beta_a50);
  PARAMETER_VECTOR(beta_Ecov_k);
  PARAMETER_VECTOR(beta_Ecov_a50);
  PARAMETER(beta_phi);
  PARAMETER_VECTOR(k_AR_pars); //AR1 process (2)
  PARAMETER_VECTOR(k_re);
  PARAMETER_VECTOR(a50_AR_pars); //AR1 process (2)
  PARAMETER_VECTOR(a50_re);
*/
  //growth parameters
  PARAMETER(beta_b);
  PARAMETER(beta_a);
  PARAMETER_VECTOR(beta_k_LVB);
  PARAMETER_VECTOR(beta_Linf);
  PARAMETER(t0);
  PARAMETER_MATRIX(beta_Ecov_b);
  PARAMETER_MATRIX(beta_Ecov_a);
  PARAMETER_MATRIX(beta_Ecov_k_LVB);
  PARAMETER_MATRIX(beta_Ecov_Linf);
  PARAMETER_VECTOR(k_LVB_AR_pars); //AR1 process (2)
  PARAMETER_VECTOR(k_LVB_re);
  PARAMETER(beta_sig_Linf);
  PARAMETER(beta_sig_L_obs);
  PARAMETER(beta_sig_W_obs);
  
  //Environmental covariate parameters
  PARAMETER_VECTOR(Ecov_mu);
  PARAMETER_MATRIX(Ecov_AR_pars); //AR1 process (2)
  PARAMETER_MATRIX(Ecov_re);
  PARAMETER_VECTOR(log_Ecov_obs_sig_scale);

  Type zero = Type(0);
  Type one = Type(1);
  Type half = Type(0.5);
  Type two = Type(2);
  vector<int> any_index_age_comp(n_indices);
  vector<int> any_fleet_age_comp(n_fleets);
  for(int i = 0; i < n_indices; i++)
  {
    any_index_age_comp(i) = 0;
    for(int y = 0; y < n_years; y++) if(use_index_paa(y,i) == 1) any_index_age_comp(i) = 1;
  }
  for(int i = 0; i < n_fleets; i++)
  {
    any_fleet_age_comp(i) = 0;
    for(int y = 0; y < n_years; y++) if(use_catch_paa(y,i) == 1) any_fleet_age_comp(i) = 1;
  }
  vector<Type> sigma2_log_NAA = exp(log_NAA_sigma*two);
  vector<Type> SSB(n_years), SSB_E(n_years);
  matrix<Type> F(n_years,n_fleets);
  matrix<Type> log_F(n_years,n_fleets);
  array<Type> pred_CAA(n_years,n_fleets,n_ages);
  array<Type> pred_catch_paa(n_years,n_fleets,n_ages);
  matrix<Type> pred_catch(n_years,n_fleets);
  array<Type> pred_IAA(n_years,n_indices,n_ages);
  array<Type> pred_index_paa(n_years,n_indices,n_ages);
  matrix<Type> pred_indices(n_years,n_indices);
  matrix<Type> log_pred_catch(n_years,n_fleets);
  matrix<Type> NAA(n_years,n_ages_pop);
  matrix<Type> NAA_obs(n_years,n_ages); //form to provide for predicting observations
  matrix<Type> pred_NAA(n_years,n_ages_pop);
  array<Type> FAA(n_years,n_fleets,n_ages);
  matrix<Type> FAA_tot(n_years,n_ages);
  matrix<Type> ZAA(n_years,n_ages);
  array<Type> QAA(n_years,n_indices,n_ages);
  matrix<Type> selblocks(n_selblocks,n_ages);
  vector<Type> q(n_indices);
  Type nll = zero; //negative log-likelihood
  vector<Type> t_paa(n_ages), t_pred_paa(n_ages);
  
  //do Ecov process and observations. Fill out Ecov_out to use later
  //int n_obs_mat = Y.size();
  int n_obs = weight.size();
  int n_years_Ecov = Ecov_obs.rows();// + max_ecov_ages;
  int n_Ecov = Ecov_obs.cols();
   
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
  nll += nll_Ecov.sum();
  //see(nll_Ecov);
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
  nll += nll_Ecov_obs.sum();
  //see(nll_Ecov_obs);
  SIMULATE REPORT(Ecov_obs);
  
  //do fit for growth data to get parameters to specify weight at age for SSB
  //k AR(1) process
  Type k_LVB_phi;
  Type k_LVB_sig;
  Type nll_k_LVB = zero;
  k_LVB_phi = -1.0 + 2.0/(1.0 + exp(-k_LVB_AR_pars(0)));
  k_LVB_sig = exp(k_LVB_AR_pars(1));
  nll_k_LVB -= dnorm(k_LVB_re(0), zero, k_LVB_sig*exp(-Type(0.5) * log(one - pow(k_LVB_phi,Type(2)))), 1);
  for(int y = 1; y < n_years_Ecov; y++) nll_k_LVB -= dnorm(k_LVB_re(y), k_LVB_phi * k_LVB_re(y-1), k_LVB_sig, 1);
  if(fit_k_LVB == 1) SIMULATE {
    k_LVB_re(0) = rnorm(zero, k_LVB_sig*exp(-Type(0.5) * log(one - pow(k_LVB_phi,Type(2)))));
    for(int y = 1; y < n_years_Ecov; y++) k_LVB_re(y) = rnorm(k_LVB_phi * k_LVB_re(y-1), k_LVB_sig);
    REPORT(k_LVB_re);
  }
  if(fit_k_LVB == 1) nll += nll_k_LVB;
  //see(nll_k_LVB);
  REPORT(nll_k_LVB);

  vector<Type> l_b(n_obs), l_a(n_obs), log_l_pred(n_obs), log_w_pred(n_obs), ll_growth(n_obs), v_log_w(n_obs);
  vector<Type> k_sum(n_obs);
  k_sum.setZero(); ll_growth.setZero();
  vector <Type> l_Linf = X_Linf * beta_Linf;
  Type v_log_l = exp(2.0 * beta_sig_L_obs) + exp(2.0 * beta_sig_Linf);
  for(int i = 0; i < n_obs; i++) 
  {
    l_b(i) = beta_b; 
    l_a(i) = beta_a;
    for(int p = 0; p < n_Ecov; p++)
    {
      for(int y = 0; y <= Ecov_maxages_b(i); y++)       l_b(i) += Ecov_y(year_obs_growth(i) - 1 - age_obs_growth(i) + y,p) * beta_Ecov_b(y,p);
      for(int y = 0; y <= Ecov_maxages_a(i); y++)       l_a(i) += Ecov_y(year_obs_growth(i) - 1 - age_obs_growth(i) + y,p) * beta_Ecov_a(y,p);
      for(int y = 0; y <= Ecov_maxages_Linf(i); y++) l_Linf(i) += Ecov_y(year_obs_growth(i) - 1 - age_obs_growth(i) + y,p) * beta_Ecov_Linf(y,p);
    }

    //do LVB k something like Millar et al. 1999. Need to specify length of beta_Ecov_k = maximum observed age and specify map to fix older ages at 0
    Type log_k0 = beta_k_LVB(0) + k_LVB_re(year_obs_growth(i) -1 - age_obs_growth(i) + 0);
    for(int p = 0; p < n_Ecov; p++) log_k0 += Ecov_y(year_obs_growth(i) - 1 - age_obs_growth(i) + 0,p) * beta_Ecov_k_LVB(0,p);
    //Type k0 = exp(beta_k_LVB(0) + k_LVB_re(year_obs_growth(i) - 1 - age_obs_growth(i) + 0) + Ecov_y(year_obs_growth(i) - 1 - age_obs_growth(i) + 0) * beta_Ecov_k_LVB(0));
    Type k0 = exp(log_k0);
    k_sum(i) = k0;
    Type log_ky = 0.0;
    for(int y = 1; y < age_obs_growth(i); y++) 
    {
      Type log_ky = beta_k_LVB(y) + k_LVB_re(year_obs_growth(i) -1 - age_obs_growth(i) + y);
      for(int p = 0; p < n_Ecov; p++) log_ky += Ecov_y(year_obs_growth(i) - 1 - age_obs_growth(i) + y,p) * beta_Ecov_k_LVB(y,p);
      k_sum(i) += exp(log_ky);
    }
    //k_sum(i) += exp(beta_k_LVB(y) + k_LVB_re(year_obs_growth(i) - 1 - age_obs_growth(i) + y) + Ecov_y(year_obs_growth(i) - 1 - age_obs_growth(i) + y) * beta_Ecov_k_LVB(y));
    log_l_pred(i) = l_Linf(i) + log(1.0 - exp(t0*k0 - k_sum(i)));
    log_w_pred(i) = l_a(i) + exp(l_b(i)) * log_l_pred(i);
    v_log_w(i) = exp(2.0 * beta_sig_W_obs) + exp(2.0 * (l_b(i) + beta_sig_Linf));
    matrix<Type> S(2,2);
    S(0,0) = v_log_l; S(1,1) = v_log_w(i);
    S(0,1) = S(1,0) = exp(l_b(i) + 2.0 * beta_sig_Linf);
    vector<Type> M(2); M(0) = log_l_pred(i); M(1) = log_w_pred(i); 
    vector<Type> obs(2); obs(0) = log(len(i)); obs(1) = log(weight(i));
    if(islen(i) == 1 & iswt(i) == 1) {
      ll_growth(i) += dbinorm(obs,M,S,1);
      SIMULATE obs = rbinorm(M,S);
    }
    else
    {
      if(islen(i) == 1) {
        ll_growth(i) += dnorm(obs(0), M(0), exp(0.5*log(S(0,0))), 1);
        SIMULATE obs(0) = rnorm(M(0), exp(0.5*log(S(0,0))));
      }
      if(iswt(i) == 1) {
        ll_growth(i) += dnorm(obs(1), M(1), exp(0.5*log(S(1,1))), 1);
        SIMULATE obs(1) = rnorm(M(1), exp(0.5*log(S(1,1))));
      }
    }
    SIMULATE{
      len(i) = exp(obs(0)); weight(i) = exp(obs(1));
    }
  }
  SIMULATE{
    REPORT(len);
    REPORT(weight);
  }
  //see(sum(ll_growth));
  if(do_growth == 1) nll -= sum(ll_growth);
  REPORT(ll_growth);
  //get weight at age
  //int n_ages = beta_Ecov_k_LVB.size();
  //int n_years_model = n_years_Ecov-n_ages;
  //int max_age = n_ages_pop;
  matrix<Type> k_age(n_years,n_ages_pop), log_Linf_age(n_years,n_ages_pop), log_a_age(n_years,n_ages_pop), log_b_age(n_years,n_ages_pop);
  matrix<Type> log_laa(n_years,n_ages_pop), log_waa(n_years,n_ages_pop);
  k_age.setZero(); log_Linf_age.setZero(); log_a_age.setZero(); log_b_age.setZero(); log_laa.setZero(); log_waa.setZero();
  vector<Type> fix(6);
  Type temp_Linf_pred = 0.0;
  for(int x = 0; x < X_Linf.cols(); x++) temp_Linf_pred += beta_Linf(x) * X_Linf(obs_for_Linf_y-1,x);
  log_Linf_age.fill(temp_Linf_pred);
  for(int y = 0; y < n_years; y++) 
  {
    for(int a = 0; a < n_ages_pop; a++) 
    {
      //log_Linf_age(y,a) = beta_Linf;
      log_a_age(y,a) = beta_a;
      log_b_age(y,a) = beta_b;
      int maxage_a = a + 1;
      if(a + 1 > beta_Ecov_a.rows()) maxage_a = beta_Ecov_a.rows();
      for(int i = 0; i < maxage_a; i++) for(int p = 0; p < n_Ecov; p++) log_a_age(y,a) += beta_Ecov_a(i,p) * Ecov_y(model_years(y)-1 - a + i - 1,p);
      int maxage_b = a + 1;
      if(a + 1 > beta_Ecov_b.rows()) maxage_b = beta_Ecov_b.rows();
      for(int i = 0; i < maxage_b; i++) for(int p = 0; p < n_Ecov; p++) log_b_age(y,a) += beta_Ecov_b(i,p) * Ecov_y(model_years(y) - 1 - a + i - 1,p);
      int maxage_Linf = a + 1;
      if(a + 1 > beta_Ecov_Linf.rows()) maxage_Linf = beta_Ecov_Linf.rows();
      for(int i = 0; i < maxage_Linf; i++) for(int p = 0; p < n_Ecov; p++) log_Linf_age(y,a) += beta_Ecov_Linf(i,p) * Ecov_y(model_years(y) - 1 - a + i - 1,p);
      Type log_k0 = beta_k_LVB(0) + k_LVB_re(model_years(y) - 1 - a + 0 - 1);
      for(int p = 0; p < n_Ecov; p++) log_k0 += beta_Ecov_k_LVB(0,p) * Ecov_y(model_years(y) - 1 - a + 0 - 1,p);
      //Type k0 = exp(beta_k_LVB(0) + k_LVB_re(model_years(y) - 1 - a + 0 - 1) + Ecov_y(model_years(y) - 1 - a + 0 - 1) * beta_Ecov_k_LVB(0));
      k_age(y,a) = exp(log_k0);
      for(int i = 1; i < a + 1; i++) {
        int ik = i; //size of beta_k_LVB is defined by the maximum age of the observations which may be less than n_ages_pop;
        if(ik >= beta_k_LVB.size()) ik = beta_k_LVB.size() - 1;
        //k_age(y,a) += exp(beta_k_LVB(i) + k_LVB_re(model_years(y) - 1 - a + i - 1) + beta_Ecov_k_LVB(i) * Ecov_y(model_years(y) - 1 - a + i - 1));
        Type log_kya = beta_k_LVB(ik) + k_LVB_re(model_years(y) - 1 - a + i - 1);
        for(int p = 0; p < n_Ecov; p++) log_kya += beta_Ecov_k_LVB(ik,p) * Ecov_y(model_years(y) - 1 - a + i - 1,p);
        k_age(y,a) += exp(log_kya); 
      }
      log_laa(y,a) = log_Linf_age(y,a) + log(1.0 - exp(t0 * exp(log_k0) - k_age(y,a)));
      log_waa(y,a) = log_a_age(y,a) + exp(log_b_age(y,a)) * log_laa(y,a);
      if(y==n_years-1 & a == 0)
      {
          fix(0) = k_LVB_re(model_years(y) - 1 - a + 0 - 1);
          fix(1) = Ecov_y(model_years(y) - 1 - a + 0 - 1,0);
          fix(2) = beta_k_LVB(0);
          fix(3) = beta_Ecov_k_LVB(0,0);
          fix(4) = exp(log_k0);
          fix(5) = k_age(y,a);
      }
    }
  }
  REPORT(fix);
  //REPORT(fix_pmat);
  ADREPORT(Ecov_y);
  ADREPORT(Ecov_phi);
  ADREPORT(Ecov_sig);
  ADREPORT(log_laa);
  ADREPORT(log_waa);
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

  selblocks = get_selblocks(n_ages, n_selblocks, n_estimated_selpars, n_other_selpars, selblock_models, estimated_selpar_pointers, other_selpar_pointers, estimated_selpars, other_selpars, selpars_lower, selpars_upper);
  for(int i = 0; i < n_indices; i++)
  {
    q(i) = q_lower(i) + (q_upper(i) - q_lower(i))/(1 + exp(-logit_q(i)));
    for(int y = 0; y < n_years; y++) 
    {
      for(int a = 0; a < n_ages; a++) QAA(y,i,a) = q(i) * selblocks(selblock_pointer_indices(y,i)-1,a);
    }
  }
 
  for(int y = 0; y < n_years; y++) for(int a = 0; a < n_ages; a++) FAA_tot(y,a) = zero;
  for(int f = 0; f < n_fleets; f++)
  {
    log_F(0,f) = log_F1(f);
    F(0,f) = exp(log_F(0,f));
    for(int a = 0; a < n_ages; a++) 
    {
      FAA(0,f,a) = F(0,f) * selblocks(selblock_pointer_fleets(0,f)-1,a);
      FAA_tot(0,a) = FAA_tot(0,a) + FAA(0,f,a);
    }
    for(int y = 1; y < n_years; y++) 
    {
      log_F(y,f) = log_F(y-1,f) + F_devs(y-1,f);
      F(y,f) = exp(log_F(y,f));
      for(int a = 0; a < n_ages; a++) 
      {
        FAA(y,f,a) = F(y,f) * selblocks(selblock_pointer_fleets(y,f)-1,a);
        FAA_tot(y,a) = FAA_tot(y,a) + FAA(y,f,a);
      }
    }
  }

  ZAA = FAA_tot + Maa;
  Type nll_N1 = 0;
  //Initial log_NAA is a random walk from initial log_recruitment
  nll_N1 -= dnorm(N1_re(0), log_N1(0), exp(log_sig_N1),1);
  SIMULATE N1_re(0) = rnorm(log_N1(0), exp(log_sig_N1));
  NAA(0,0) = exp(log_N1(0));
  NAA_obs(0,0) = NAA(0,0);
  for(int a = 1; a < n_ages_pop-1; a++) {
    nll_N1 -= dnorm(N1_re(a), N1_re(a-1), exp(log_sig_N1),1);
    SIMULATE N1_re(a) = rnorm(N1_re(a-1),exp(log_sig_N1));
  }
  for(int a = 1; a < n_ages_pop; a++) {
    NAA(0,a) = exp(N1_re(a-1));
    if(a < n_ages) NAA_obs(0,a) = NAA(0,a);
    else NAA_obs(0,n_ages-1) += NAA(0,a);
  }
  nll += nll_N1;
  //see(nll_N1);
  REPORT(nll_N1);
  
  SSB.setZero();
  SSB_E.setZero();
  for(int a = 0; a < n_ages; a++) SSB(0) += NAA_obs(0,a) * waa(waa_pointer_ssb-1,0,a) * mature(0,a) * exp(-ZAA(0,a) * fracyr_SSB(0));
  for(int a = 0; a < n_ages_pop; a++) 
  {
    int tage = a;
    if(a > n_ages-1) tage = n_ages-1;
    Type ssb_naa = NAA(0,a) * exp(-ZAA(0,tage)*fracyr_SSB(0)) * mature(0,tage);
    if(use_growth_model == 1) SSB_E(0) += ssb_naa * exp(log_waa(0,a));
    else SSB_E(0) += ssb_naa * waa(waa_pointer_ssb-1,0,tage);
  }
  for(int y = 1; y < n_years; y++) 
  {
    for(int a = 0; a < n_ages_pop; a++) {
      NAA(y,a) = exp(log_NAA(y-1,a)); //random effects NAA
      if(a < n_ages) NAA_obs(y,a) = NAA(y,a);
      else NAA_obs(y,n_ages-1) += NAA(y,a);
    }
    for(int a = 0; a < n_ages; a++) 
    {
      SSB(y) += NAA(y,a) * waa(waa_pointer_ssb-1,y,a) * mature(y,a) * exp(-ZAA(y,a)*fracyr_SSB(y));
    }
    for(int a = 0; a < n_ages_pop; a++) 
    {
      int tage = a;
      if(a > n_ages-1) tage = n_ages-1;
      Type ssb_naa = NAA(y,a) * exp(-ZAA(y,tage)*fracyr_SSB(y)) * mature(y,tage);
      if(use_growth_model == 1) SSB_E(y) += ssb_naa * exp(log_waa(y,a));
      else SSB_E(y) += ssb_naa * waa(waa_pointer_ssb-1,y,tage);
    }
  }
      
  for(int a = 0; a < n_ages_pop; a++) pred_NAA(0,a) = NAA(0,a);
  
  Type nll_NAA = 0;
  for(int y = 1; y < n_years; y++)
  {
    if(recruit_model == 1) pred_NAA(y,0) = NAA(y-1,0); //random walkNAA(y,1)
    else
    {
      if(recruit_model == 2) pred_NAA(y,0) = exp(mean_rec_pars(0)); //random about mean
      else //BH stock recruit
      {
        Type tSSB = SSB(y-1);
        if(use_growth_model == 1) tSSB = SSB_E(y-1);
        if(recruit_model == 3) //BH stock recruit
        {
          pred_NAA(y,0) = exp(mean_rec_pars(0) + log(tSSB) - log(one + exp(mean_rec_pars(1)) * tSSB));
        }
        else //Ricker stock recruit
        {
          pred_NAA(y,0) = exp(mean_rec_pars(0) + log(tSSB) - exp(mean_rec_pars(1)) * tSSB); 
        }
      }
    }
    for(int a = 1; a < n_ages_pop-1; a++) 
    {
      if(a < n_ages) pred_NAA(y,a) = NAA(y-1,a-1) * exp(-ZAA(y-1,a-1));
      else pred_NAA(y,a) = NAA(y-1,a-1) * exp(-ZAA(y-1,n_ages-1));
    }
    //note ZAA is the same in both components below because n_ages << n_ages_pop
    pred_NAA(y,n_ages_pop-1) = NAA(y-1,n_ages_pop-2) * exp(-ZAA(y-1,n_ages-1)) + NAA(y-1,n_ages_pop-1) * exp(-ZAA(y-1,n_ages-1));
  }
  
  for(int y = 1; y < n_years; y++)
  {
    for(int a = 0; a < n_ages_pop; a++)
      nll_NAA -= dnorm(log(NAA(y,a)), log(pred_NAA(y,a)), exp(log_NAA_sigma(NAA_sigma_pointers(a)-1)), 1);
  }
  SIMULATE {
    for(int y = 1; y < n_years; y++)
    {
      for(int a = 0; a < n_ages_pop; a++){
        log_NAA(y-1,a) = rnorm(log(pred_NAA(y,a)), exp(log_NAA_sigma(NAA_sigma_pointers(a)-1)));
        NAA(y,a) = exp(log_NAA(y-1,a));
        if(a < n_ages) NAA_obs(y,a) = NAA(y,a);
        else NAA_obs(y,n_ages-1) += NAA(y,a);
      }
    }
  }
  //see(nll_NAA);
  nll += nll_NAA;
  REPORT(nll_NAA);
  Type nll_agg_catch = zero;
  matrix<Type> nll_catch_acomp(n_years, n_fleets);
  nll_catch_acomp.setZero();
  for(int y = 0; y < n_years; y++)
  {
    int acomp_par_count = 0;
    for(int f = 0; f < n_fleets; f++)
    {
      pred_catch(y,f) = zero;
      Type tsum = zero;
      for(int a = 0; a < n_ages; a++) 
      {
        pred_CAA(y,f,a) =  NAA_obs(y,a) * FAA(y,f,a) * (1 - exp(-ZAA(y,a)))/ZAA(y,a);
        pred_catch(y,f) += waa(waa_pointer_fleets(f)-1,y,a) * pred_CAA(y,f,a);
        tsum += pred_CAA(y,f,a);
      }
      nll_agg_catch -= dnorm(log(agg_catch(y,f)), log(pred_catch(y,f)), agg_catch_sigma(y,f), 1);
      SIMULATE agg_catch(y,f) = exp(rnorm(log(pred_catch(y,f)), agg_catch_sigma(y,f)));
      if(any_fleet_age_comp(f) == 1)
      {
        vector<Type> acomp_pars(n_age_comp_pars_fleets(f));
        for(int j = 0; j < n_age_comp_pars_fleets(f); j++) 
        {
          acomp_pars(j) = catch_paa_pars(acomp_par_count);
          acomp_par_count++;
        }
        if(use_catch_paa(y,f) == 1) 
        {
          for(int a = 0; a < n_ages; a++)
          {
            pred_catch_paa(y,f,a) = pred_CAA(y,f,a)/tsum;
            t_pred_paa(a) = pred_catch_paa(y,f,a);
            t_paa(a) = catch_paa(f * n_years + y,a);
          }
/*          if(y == 0)
          {
            see(f);
            see(t_pred_paa);
            see(t_paa);
            see(catch_aref(y,f));
            see(catch_Neff(y,f));
            see(NAA_obs.row(y));
            see(ZAA.row(y));
            see(FAA_tot.row(y));
            see(pred_catch.row(y));
          }*/
          nll_catch_acomp(y,f) -= get_acomp_ll(y, n_ages, catch_Neff(y,f), age_comp_model_fleets(f), t_paa, t_pred_paa, acomp_pars, catch_aref(y,f));
          SIMULATE {
            t_paa = sim_acomp(y, n_ages, catch_Neff(y,f), age_comp_model_fleets(f), t_paa, t_pred_paa, acomp_pars, catch_aref(y,f));
            for(int a = 0; a < n_ages; a++) catch_paa(f * n_years + y, a) = t_paa(a);
          }

        }
      }
    }
  }
  //see(nll_agg_catch);
  nll += nll_agg_catch;
  //see(nll_catch_acomp);
  nll += nll_catch_acomp.sum();
  REPORT(nll_agg_catch);
  REPORT(nll_catch_acomp);
  REPORT(pred_catch);
  REPORT(pred_catch_paa);
  REPORT(selblocks);

  vector<Type> nll_agg_indices(n_indices), nll_index_acomp(n_indices);
  nll_agg_indices.setZero();
  nll_index_acomp.setZero();
  for(int y = 0; y < n_years; y++)
  {
    int acomp_par_count = 0;
    for(int i = 0; i < n_indices; i++) 
    {
      pred_indices(y,i) = zero;
      Type tsum = zero;
      for(int a = 0; a < n_ages; a++) 
      {
        pred_IAA(y,i,a) =  NAA_obs(y,a) * QAA(y,i,a) * exp(-ZAA(y,a) * fracyr_indices(y,i));
        if(units_indices(i) == 1) pred_indices(y,i) += waa(waa_pointer_indices(i)-1,y,a) * pred_IAA(y,i,a);
        else pred_indices(y,i) += pred_IAA(y,i,a);
      }
      for(int a = 0; a < n_ages; a++) 
      {
        if(units_index_paa(i) == 1) pred_IAA(y,i,a) = waa(waa_pointer_indices(i)-1,y,a) * pred_IAA(y,i,a);
        tsum += pred_IAA(y,i,a);
      }
      
      if(use_indices(y,i) == 1)
      {
        nll_agg_indices(i) -= dnorm(log(agg_indices(y,i)), log(pred_indices(y,i)), agg_index_sigma(y,i), 1);
        SIMULATE agg_indices(y,i) = exp(rnorm(log(pred_indices(y,i)), agg_index_sigma(y,i)));
      }
      if(any_index_age_comp(i) == 1)
      {
        vector<Type> acomp_pars(n_age_comp_pars_indices(i));
        for(int j = 0; j < n_age_comp_pars_indices(i); j++) 
        {
          acomp_pars(j) = index_paa_pars(acomp_par_count);
          acomp_par_count++;
        }
        if(use_index_paa(y,i) > 0)
        {
          for(int a = 0; a < n_ages; a++)
          {
            pred_index_paa(y,i,a) = pred_IAA(y,i,a)/tsum;
            t_pred_paa(a) = pred_index_paa(y,i,a);
            t_paa(a) = index_paa(i * n_years + y,a);
          }
          
          nll_index_acomp(i) -= get_acomp_ll(y, n_ages, index_Neff(y,i), age_comp_model_indices(i), t_paa, t_pred_paa, acomp_pars, index_aref(y,i));
          SIMULATE {
            t_paa = sim_acomp(y, n_ages, index_Neff(y,i), age_comp_model_indices(i), t_paa, t_pred_paa, acomp_pars, index_aref(y,i));
            for(int a = 0; a < n_ages; a++) index_paa(i * n_years + y, a) = t_paa(a);
          }
        }
      }
    }
  }
  //see(nll_agg_indices);
  nll += nll_agg_indices.sum();
  //see(nll_index_acomp);
  nll += nll_index_acomp.sum();
  REPORT(nll_agg_indices);
  REPORT(nll_index_acomp);
  REPORT(pred_indices);
  REPORT(pred_index_paa);
  ////////////////////////////////////////////////////////////////////
  
  
  //////////////////////////////////////////////////////////////////////
  vector<Type> log_SSB = log(SSB);
  vector<Type> log_SSB_E = log(SSB_E);
/*  vector<Type> SSB_E(n_years);
  SSB_E.setZero();
  for(int y = 0; y < n_years; y++) for(int a = 0; a < n_ages_pop; a++) 
  {
    Type SSB_a = 0;
    if(a < n_ages) SSB_a = NAA(y,a) * exp(-ZAA(y,a)*fracyr_SSB(y)) * mature(y,a);
    else SSB_a = NAA(y,a) * exp(-ZAA(y,n_ages-1)*fracyr_SSB(y)) * mature(y,n_ages-1);
    if(use_growth_model == 1) SSB_a *= exp(log_waa(y,a));
    else {
      if(a < n_ages) SSB_a *= waa(waa_pointer_ssb-1,y,a);
      else SSB_a *= waa(waa_pointer_ssb-1,y,n_ages-1);
    }
    SSB_E(y) += SSB_a;
  }
  vector<Type> log_SSB_E = log(SSB_E);
*/  
  vector<Type> wssb(n_ages), wcatch(n_ages);
  vector<Type> log_SPR_0(n_years), log_SPR_0_E(n_years), log_SPR40(n_years), log_SPR40_E(n_years), log_F40(n_years), log_F40_E(n_years);
  vector<Type> log_SSB40(n_years), log_SSB40_E(n_years), log_YPR40(n_years), log_YPR40_E(n_years), log_Y40(n_years), log_Y40_E(n_years);
  vector<Type> msy_proxy(3);
  for(int a = 0; a < n_ages; a++) 
  {
    wssb(a) = waa(waa_pointer_ssb-1, n_years-1, a);
    wcatch(a) = waa(waa_pointer_totcatch-1, n_years-1, a);
  }
  vector<Type> M2 = Maa.row(n_years-1);
  vector<Type> maturity = mature.row(n_years-1);
  vector<Type> sel_tot = FAA_tot.row(n_years-1)/FAA_tot(n_years-1,n_ages-1);
  matrix<Type> wssb_E(n_years,n_ages);
  for(int y = 0; y < n_years; y++)
  {
    if(use_growth_model == 1) {
      vector<Type> paaplus = NAA.row(y).tail(n_ages_pop-n_ages);
      paaplus = paaplus / paaplus.sum();
      vector<Type> waaplus = log_waa.row(y).tail(n_ages_pop-n_ages);
      waaplus = exp(waaplus);
      waaplus = waaplus * paaplus;
      wssb_E(y,n_ages-1) = waaplus.sum();
      for(int a= 0; a < (n_ages-1); a++) 
      {
        wssb_E(y,a) = exp(log_waa(y,a));
      }
    }
    else for(int a = 0; a < n_ages; a++) wssb_E(y,a) = waa(waa_pointer_ssb-1, y , a);
    vector<Type> wssb_E_y = wssb_E.row(y);
    log_SPR_0(y) = log(get_SPR_0(M2, maturity, wssb, fracyr_SSB(n_years-1)));
    msy_proxy = getFXSPR(exp(log_SPR_0(y)), percentSPR, M2, sel_tot, maturity, wssb, wcatch, fracyr_SSB(n_years-1), 10);
    log_F40(y) = log(msy_proxy(0));
    log_SPR40(y) = log(msy_proxy(1));
    log_YPR40(y) = log(msy_proxy(2));
    log_SSB40(y) = log(msy_proxy(1)*NAA.col(0).mean());
    log_Y40(y) = log(msy_proxy(2)*NAA.col(0).mean());
    log_SPR_0_E(y) = log(get_SPR_0(M2, maturity, wssb_E_y, fracyr_SSB(n_years-1)));
    Type SPR0 = exp(log_SPR_0_E(y));
    msy_proxy = getFXSPR(SPR0, percentSPR, M2, sel_tot, maturity, wssb_E_y, wcatch, fracyr_SSB(n_years-1), 10);
    log_SPR40_E(y) = log(msy_proxy(1));
    log_YPR40_E(y) = log(msy_proxy(2));
    log_F40_E(y) = log(msy_proxy(0));
    log_SSB40_E(y) = log(msy_proxy(1)*NAA.col(0).mean());
    log_Y40_E(y) = log(msy_proxy(2)*NAA.col(0).mean());
  }
  if(recruit_model>2){
    vector<Type> log_SPR_MSY(n_years), log_FMSY(n_years), log_YPR_MSY(n_years), log_R_MSY(n_years), log_SSB_MSY(n_years);
    Type F_init = 0.1;
    for(int y = 0; y < n_years; y++){
      vector<Type> wssb_y = wssb_E.row(y);
      log_FMSY(y) = get_FMSY(mean_rec_pars(0), mean_rec_pars(1), M2, sel_tot, wcatch, wssb_y,
        maturity, fracyr_SSB(n_years-1), log_SPR_0_E(y), recruit_model, F_init);
      log_SPR_MSY(y) = log(get_SPR(log_FMSY(y), M2, sel_tot, maturity, wssb_y, fracyr_SSB(n_years-1)));
      log_YPR_MSY(y) = log(get_YPR(log_FMSY(y), M2, sel_tot, wcatch));
      if(recruit_model == 3) log_R_MSY(y) = log((exp(mean_rec_pars(0)) - 1/exp(log_SPR_MSY(y))) / exp(mean_rec_pars(1))); //bh
      else log_R_MSY(y) = log(mean_rec_pars(0) + log_SPR_MSY(y)) - mean_rec_pars(1) - log_SPR_MSY(y); //ricker
      log_SSB_MSY(y) = log_SPR_MSY(y) + log_R_MSY(y);
    }
    REPORT(log_FMSY);
    REPORT(log_SPR_MSY);
    REPORT(log_YPR_MSY);
    REPORT(log_R_MSY);
    REPORT(log_SSB_MSY);
    ADREPORT(log_FMSY);
    ADREPORT(log_SPR_MSY);
    ADREPORT(log_YPR_MSY);
    ADREPORT(log_R_MSY);
    ADREPORT(log_SSB_MSY);
  }

  //if(reportMode==0){
    REPORT(NAA);
    REPORT(NAA_obs);
    REPORT(pred_NAA);
    REPORT(SSB);
    REPORT(selblocks);
    REPORT(q);
    REPORT(F);
    REPORT(log_waa);
    REPORT(log_SSB_E);
    REPORT(log_SSB);
    REPORT(log_SPR_0);
    REPORT(log_SPR40);
    REPORT(log_YPR40);
    REPORT(log_F40);
    REPORT(log_SSB40);
    REPORT(log_Y40);
    REPORT(log_SPR_0_E);
    REPORT(log_SPR40_E);
    REPORT(log_YPR40_E);
    REPORT(log_F40_E);
    REPORT(log_SSB40_E);
    REPORT(log_Y40_E);

    //REPORT(pmat);
    ADREPORT(log_SSB);
    ADREPORT(log_SSB_E);
    ADREPORT(log_F);
    ADREPORT(log_SPR_0);
    ADREPORT(log_SPR40);
    ADREPORT(log_YPR40);
    ADREPORT(log_F40);
    ADREPORT(log_SSB40);
    ADREPORT(log_Y40);
    ADREPORT(log_SPR_0_E);
    ADREPORT(log_SPR40_E);
    ADREPORT(log_YPR40_E);
    ADREPORT(log_F40_E);
    ADREPORT(log_SSB40_E);
    ADREPORT(log_Y40_E);
    //ADREPORT(pmat);
  //}
  REPORT(nll);
  //see(nll);
  return nll;
}

