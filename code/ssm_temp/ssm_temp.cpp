// This was edited from ssm_env_ssb_v4.cpp fro GB cod paper. Fits code/ss_lvb_temp simultaneously with a state-space age-structured model.
// Uses predicted weight at age in estimating SSB and reference points, BUT NOT FOR CATCH AND SURVEYS.

#include <TMB.hpp>
#include <iostream>

#define see(object) std::cout << #object ":\n" << object << "\n";

template <class Type> 
Type square(Type x){return x*x;}

template<class Type>
Type dbetabinom(Type x, Type n, Type mu, Type phi, int do_log)
{
  Type ll = lgamma(n + 1.0) - lgamma(x + 1.0) - lgamma(n - x + 1.0) + 
	  lgamma(x + mu*phi) + lgamma(n - x +(1-mu)*phi) - lgamma(n + phi) +
	  lgamma(phi) - lgamma(mu*phi) - lgamma((1-mu)*phi);
  if(do_log == 1) return(ll);
  else return(exp(ll));  
}

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

template <class Type>
Type get_SPR(Type log_F, vector<Type> M, vector<Type> sel, vector<Type> mat, vector<Type> waassb, Type fracyearSSB)
{
  int n_ages = M.size();
  Type zero = Type(0);
  Type one = Type(1);
  Type SPR = zero, ntemp = one;
  vector<Type> F(n_ages), Z(n_ages);
 
  F = exp(log_F) * sel;
  Z = F + M;
  for(int age=0; age<n_ages-1; age++)
  {
    SPR += ntemp * mat(age) * waassb(age) * exp(-fracyearSSB * Z(age));
    ntemp *= exp(-Z(age));
  }
  ntemp /= one-exp(-Z(n_ages-1));
  SPR += ntemp * mat(n_ages-1) * waassb(n_ages-1) * exp(-fracyearSSB*Z(n_ages-1));
  
  return SPR;
}
template <class Type>
Type get_YPR(Type log_F, vector<Type> M, vector<Type> sel, vector<Type> waacatch)
{
  int n_ages = M.size();
  Type zero = Type(0);
  Type one = Type(1);
  Type YPR = zero, ntemp = one;
  vector<Type> F(n_ages), Z(n_ages);
 
  F = exp(log_F) * sel;
  Z = F + M;
  for(int age=0; age<n_ages-1; age++)
  {
    YPR += ntemp * F(age) * waacatch(age) * (one - exp(-Z(age)))/Z(age);
    ntemp *= exp(-Z(age));
  }
  ntemp /= one - exp(-Z(n_ages-1));
  YPR += ntemp * F(n_ages-1) * waacatch(n_ages-1) * (one - exp(-Z(n_ages-1)))/Z(n_ages-1);
  
  return YPR;
}
template <class Type>
Type get_SPR_0(vector<Type> M, vector<Type> mat, vector<Type> waassb, Type fracyearSSB)
{
  int n_ages = M.size();
  Type SPR_0 = Type(0.0);
  Type ntemp0 = Type(1.0);
  for (int a = 0; a < n_ages - 1; a++)
  {
    SPR_0 += ntemp0 * mat(a) * waassb(a) * exp(-fracyearSSB * M(a));
    ntemp0 *= exp(-M(a));
  }
  ntemp0 /= Type(1.0)-exp(-M(n_ages-1));
  SPR_0 += ntemp0 * mat(n_ages-1) * waassb(n_ages-1) * exp(-fracyearSSB * M(n_ages-1));
  return SPR_0;
}

template <class Type>
Type get_dSPRdF(Type log_F, vector<Type> M, vector<Type> sel, vector<Type> mat, vector<Type> waassb, Type fracyearSSB)
{
  int n_ages = M.size();
  Type zero = Type(0);
  Type one = Type(1);
  vector<Type> S(n_ages), dSdF(n_ages), Sa(n_ages), dSadF(n_ages), F(n_ages), Z(n_ages);
  Type T = zero, dTdF = zero, dSPRdF = zero, cum_sel = zero;
  F = exp(log_F) * sel;
  Z = F + M;
  T = one/(one - exp(-Z(n_ages-1)));
  dTdF = - sel(n_ages-1) * exp(-Z(n_ages-1)) * T * T;
  for(int a = 0; a < n_ages; a++)
  {
    Sa(a) = exp(-fracyearSSB * Z(a));
    dSadF(a) = - fracyearSSB * sel(a) * Sa(a);
    if(a == 0) 
    {
      S(a) = one;
      dSdF(a) = zero;
    }
    else 
    {
      S(a) = S(a-1) * exp(-Z(a-1));
      dSdF(a) = - cum_sel * S(a);
    }
    cum_sel += sel(a);
    
    if(a < n_ages-1) dSPRdF += waassb(a) * mat(a) * (Sa(a)* dSdF(a) + S(a) * dSadF(a));
    else dSPRdF += waassb(a) * mat(a) * ((Sa(a)* dSdF(a) + S(a) * dSadF(a)) * T + S(a) * Sa(a) *  dTdF);
  }
  return dSPRdF;
}
template <class Type>
vector<Type> getFXSPR(Type SPR_0, Type percentSPR, vector<Type> M, vector<Type> sel, vector<Type> mat, vector<Type> waassb, vector<Type> waacatch, Type fracyearSSB, int n)
{
  vector<Type> log_FXSPR(n), FXSPR(n), SPR(n), YPR(n), dSPRdF(n), SPR_vals(3);
  log_FXSPR(0) = log(Type(0.2));
  FXSPR(0) = exp(log_FXSPR(0));
  SPR(0) = get_SPR(log_FXSPR(0), M, sel, mat, waassb, fracyearSSB);
  dSPRdF(0) = get_dSPRdF(log_FXSPR(0), M, sel, mat, waassb, fracyearSSB);
  for (int i=0; i<n-1; i++)
  {
    //need change of variable because of log scale
    Type dSPRdf = FXSPR(i) * dSPRdF(i);
    log_FXSPR(i+1) = log_FXSPR(i) - (SPR(i) - Type(0.01)*percentSPR*SPR_0)/dSPRdf; 
    FXSPR(i+1) = exp(log_FXSPR(i+1));
    SPR(i+1) = get_SPR(log_FXSPR(i+1), M, sel, mat, waassb, fracyearSSB);
    YPR(i+1) = get_YPR(log_FXSPR(i+1), M, sel, waacatch);
    dSPRdF(i+1) = get_dSPRdF(log_FXSPR(i+1), M, sel, mat, waassb, fracyearSSB);
  }
  SPR_vals(0) = FXSPR(n-1);
  SPR_vals(1) = SPR(n-1);
  SPR_vals(2) = YPR(n-1);
  return SPR_vals;
}

template <class Type>
matrix<Type> get_selblocks(int n_ages, int n_selblocks, int n_estimated_selpars, int n_other_selpars, vector<int> selblock_models, vector<int> estimated_pointers, vector<int> other_pointers, vector<Type> estimated_selpars, vector<Type> other_selpars, vector<Type> selpars_lower, vector<Type> selpars_upper)
{
  Type zero = Type(0);
  Type one = Type(1);
  Type half = Type(0.5);
  Type two = Type(2);
  int n_selpars = n_estimated_selpars + n_other_selpars;
  Type a50_1 = zero, k_1 = zero, a50_2 = zero, k_2 = zero;
  matrix<Type> selectivity_blocks(n_selblocks,n_ages);
  vector<Type> selpars(n_selpars);
  if(n_estimated_selpars>0) for(int i = 0; i < n_estimated_selpars; i++)
  {
    selpars(estimated_pointers(i)-1) = selpars_lower(i) + (selpars_upper(i) - selpars_lower(i))/(1 + exp(-estimated_selpars(i)));
  }
  if(n_other_selpars>0) for(int i = 0; i < n_other_selpars; i++) selpars(other_pointers(i)-1) = other_selpars(i);
  int count = 0;
  for(int i = 0; i < n_selblocks; i++) 
  {
    if (selblock_models(i)==1)
    { //proportions at age
      for(int a = 0; a < n_ages; a++) 
      {
        selectivity_blocks(i,a) = selpars(count);
        count++;
      }
    }
    else
    { //logistic or double-logistic
      a50_1 = selpars(count); // a50 parameter
      count++;
      k_1 = selpars(count); //  1/slope
      count++;
      if (selblock_models(i)==2) 
      { //increasing logistic
        Type age = zero;
        for (int a = 0; a < n_ages; a++) 
        {
          age += one;
          selectivity_blocks(i,a) = one/(one + exp(-(age - a50_1)/k_1));
        }
        //for (int a = 0; a < n_ages; a++) selectivity_blocks(i,a) = selectivity_blocks(i,a)/selectivity_blocks(i,n_ages-1);
      }
      else
      { //double logistic
        a50_2 = selpars(count);
        count++;
        k_2 = selpars(count);
        count++;
        Type age = zero;
        for (int a = 0; a < n_ages; a++)
        {
          age += one;
 	        selectivity_blocks(i,a) = one/(one + exp(-(age - a50_1)/k_1));
          selectivity_blocks(i,a) *= one/(one + exp((age - a50_2)/k_2)); //1-p
        }
      }
    }
  }
  return selectivity_blocks;
}

template<class Type>
Type get_acomp_ll(int year, int n_ages, Type Neff, int age_comp_model, vector<Type> paa_obs, vector<Type> paa_pred, vector<Type> age_comp_pars, int aref)
{
  Type zero = Type(0);
  Type one = Type(1);
  Type half = Type(0.5);
  Type two = Type(2);
  vector<Type> temp_n(n_ages); 
  Type temp_Neff = zero, ll = zero;
  if(age_comp_model == 1) //multinomial
  { 
    temp_Neff = Neff * exp(age_comp_pars(0));
    temp_n = temp_Neff * paa_obs;
    ll = lgamma(temp_Neff + one);
    for(int a = 0; a < n_ages; a++) ll += -lgamma(temp_n(a) + one) + temp_n(a) * log(paa_pred(a));
  }
  if(age_comp_model == 2) //dirichlet-multinomial
  { 
    temp_Neff = Neff;
    temp_n = temp_Neff * paa_obs;
    ll = lgamma(temp_Neff + one) + lgamma(exp(age_comp_pars(0))) - lgamma(temp_Neff + exp(age_comp_pars(0)));
    for(int a = 0; a < n_ages; a++) ll += -lgamma(temp_n(a) + one) + lgamma(temp_n(a) + exp(age_comp_pars(0)) * paa_pred(a)) -
      lgamma(exp(age_comp_pars(0)) * paa_pred(a));
  } 
  if(age_comp_model == 3) //dirichlet
  { 
    Type obs = zero, pred = zero, obs_2 = zero, pred_2 = zero;
    for(int a = aref-1; a < n_ages; a++) 
    {
      obs_2 += paa_obs(a);
      pred_2 += paa_pred(a);
    }
    ll = lgamma(exp(age_comp_pars(0)));
    for(int a = 0; a < aref-1; a++) 
    {
      pred += paa_pred(a);
      obs += paa_obs(a);
      if(paa_obs(a) > Type(1.0e-15))
      {
        ll +=  -lgamma(exp(age_comp_pars(0)) * pred) + (exp(age_comp_pars(0)) * pred - one) * log(obs);
        pred = zero;
        obs = zero;
      }
      //else pooling with next age
    }
    //add in the last age class(es).
    ll += -lgamma(exp(age_comp_pars(0)) * pred_2) + (exp(age_comp_pars(0)) * pred_2 - one) * log(obs_2); 
  }
  if(age_comp_model == 4) //zero-one inflated logistic normal. Inspired by zero-one inflated beta in Ospina and Ferrari (2012).
  {
    vector<Type> X(n_ages), p0(n_ages);
    Type mu = zero, sd = zero, pos_obs = zero, pos_pred = zero, pos_obs_l = zero, pos_pred_l = zero, pos_obs_sum = zero;
    Type pos_pred_sum = zero, y = zero;
    X = log(paa_pred + Type(1.0e-15)) - log(one - paa_pred + Type(1.0e-15));
    p0 = one/(one + exp(exp(age_comp_pars(1))*(X - age_comp_pars(0)))); //prob of zero declines with proportion caught
    sd = exp(age_comp_pars(2));
    int last_pos = 0;
    pos_obs_sum = sum(paa_obs); 
    for(int a = 0; a < n_ages; a++) if(paa_obs(a) > Type(1.0e-15))
    {
      pos_pred_sum += paa_pred(a);
      last_pos = a;
    }
    //logistic applies only to proportions of non-zero observations
    pos_obs_l = paa_obs(last_pos)/pos_obs_sum; 
    pos_pred_l = paa_pred(last_pos)/pos_pred_sum;
    for(int a = 0; a < n_ages; a++)
    {
      if(paa_obs(a) < Type(1.0e-15)) ll += log(p0(a) + Type(1.0e-15));
      else
      {
        ll += log(one - p0(a) + Type(1.0e-15));
        if(a < last_pos) //add in logistic-normal for positive observations less than last observed age class
        {
          pos_pred = paa_pred(a)/pos_pred_sum;
          pos_obs = paa_obs(a)/pos_obs_sum;
          y = log(pos_obs) - log(pos_obs_l);
          mu = log(pos_pred + Type(1.0e-15)) - log(pos_pred_l + Type(1.0e-15));
          ll += -half * (log(two * M_PI) + square((y - mu)/sd)) - log(sd) - log(pos_obs);
        }
      }
    }
    ll -= log(pos_obs_l); //add in the last observed age class(es).
  }
  if(age_comp_model == 5) //logistic normal. Pool zero observations with adjacent age classes.
  {
    Type mu = zero, sd = zero, obs = zero, pred = zero, obs_2 = zero, pred_2 = zero, y = zero;
    for(int a = aref-1; a < n_ages; a++) 
    {
      obs_2 += paa_obs(a);
      pred_2 += paa_pred(a);
    }
    for(int a = 0; a < aref-1; a++) 
    {
      pred += paa_pred(a);
      obs += paa_obs(a);
      if(paa_obs(a) > Type(1.0e-15))
      {
        sd = exp(age_comp_pars(0)-half*log(Neff));
        y = log(obs) - log(obs_2);
        mu = log(pred + Type(1.0e-15)) - log(pred_2 + Type(1.0e-15));
        ll += -half * (log(two * M_PI) + square((y - mu)/sd)) - log(sd) - log(obs);
        pred = zero;
        obs = zero;
      }
      //else pooling with next age
    }
    ll -= log(obs_2); //add in the last age class(es).
  }
  if(age_comp_model == 6) //zero-one inflated logistic normal where p0 is a function of binomial sample size. 2 parameters
  {
    vector<Type> p0(n_ages);
    Type n_e = zero, mu = zero, sd = zero, pos_obs = zero, pos_pred = zero, pos_obs_l = zero, pos_pred_l = zero, pos_obs_sum = zero; 
    Type pos_pred_sum = zero, y = zero;
    n_e = exp(age_comp_pars(0));
    p0 = exp(n_e * log(one-paa_pred + Type(1.0e-15))); //prob of zero declines with proportion caught
    sd = exp(age_comp_pars(1));
    int last_pos = 0;
    pos_obs_sum = sum(paa_obs); 
    for(int a = 0; a < n_ages; a++) if(paa_obs(a) > Type(1.0e-15))
    {
      pos_pred_sum += paa_pred(a);
      last_pos = a;
    }
    //logistic applies only to proportions of non-zero observations
    pos_obs_l = paa_obs(last_pos)/pos_obs_sum; 
    pos_pred_l = paa_pred(last_pos)/pos_pred_sum;
    for(int a = 0; a < n_ages; a++)
    {
      if(paa_obs(a) < Type(1.0e-15)) ll += log(p0(a) + Type(1.0e-15));
      else
      {
        ll += log(one - p0(a) + Type(1.0e-15));
        if(a < last_pos) //add in logistic-normal for positive observations less than last observed age class
        {
          pos_pred = paa_pred(a)/pos_pred_sum;
          pos_obs = paa_obs(a)/pos_obs_sum;
          y = log(pos_obs) - log(pos_obs_l);
          mu = log(pos_pred + Type(1.0e-15)) - log(pos_pred_l + Type(1.0e-15));
          ll += -half * (log(two * M_PI) + square((y - mu)/sd)) - log(sd) - log(pos_obs);
        }
      }
    }
    ll -= log(pos_obs_l); //add in the last observed age class(es).
  }
  return ll;
}

template<class Type>
vector<Type> rmultinom(Type N, vector<Type> p)
{
  //multinomial
  int dim = p.size();
  vector<Type> x(dim);
  int Nint = CppAD::Integer(N);
  x.setZero();
  for(int i = 0; i < Nint; i++)
  {
    Type y = runif(0.0,1.0);
    for(int a = 0; a < dim; a++) if(y < p.head(a+1).sum()) 
    {
      x(a) += 1.0;
      break;
    }
  }
  return x;
}
template<class Type>
vector<Type> rdirichlet(vector<Type> p, Type phi)
{
  vector<Type> alpha = p * phi;
  vector<Type> obs = rgamma(alpha,Type(1.0));
  obs = obs/obs.sum();
  return obs;
}

template<class Type>
vector<Type> rdirmultinom(Type N, vector<Type> p, Type phi) //dirichlet generated from iid gammas
{
  int Nint = CppAD::Integer(N);
  int dim = p.size();
  vector<Type> obs(dim);
  obs.setZero();
  for(int i = 0; i < Nint; i++)
  {
    vector<Type> dp = rdirichlet(p, phi);
    obs = obs + rmultinom(Type(1),dp);
  }
  return(obs);
}

template<class Type>
vector<Type> sim_acomp(int year, int n_ages, Type Neff, int age_comp_model, vector<Type> paa_obs, vector<Type> paa_pred, vector<Type> age_comp_pars, int aref)
{
  vector<Type> obs(n_ages);
  obs.setZero();
  if(age_comp_model == 1) 
  { 
    //int N = CppAD::Integer(Neff);
    obs = rmultinom(Neff, paa_pred);
    obs = obs/obs.sum();// proportions
  }
  if(age_comp_model == 2) //dirichlet-multinomial. dirichlet generated from iid gammas and multinomial from uniform
  { 
    //int N = CppAD::Integer(Neff);
    obs = rdirmultinom(Neff,paa_pred,exp(age_comp_pars(0)));
    obs = obs/obs.sum();// proportions
  } 
  if(age_comp_model == 3) //dirichlet generated from iid gammas
  { 
    Type obs_2 = 0.0;
    vector<Type> best_obs = rdirichlet(paa_pred, exp(age_comp_pars(0)));
    obs_2 = best_obs.tail(n_ages-aref+1).sum(); // .tail last n_ages-aref+1 components
    for(int a = aref-1; a < n_ages; a++) if(paa_obs(a) > Type(1.0e-15)) obs(a) = obs_2;
    obs_2 = 0.0;
    for(int a = 0; a < aref-1; a++) 
    {
      obs_2 += best_obs(a);
      if(paa_obs(a) > Type(1.0e-15))
      {
        obs(a) = obs_2;
        obs_2 = 0.0;
      }
      else obs(a) = 0.0;
      //else pooling with next age
    }
  }
  if(age_comp_model == 4) //zero-one inflated logistic normal. Inspired by zero-one inflated beta in Ospina and Ferrari (2012).
  {
    vector<Type> X = log(paa_pred + Type(1.0e-15)) - log(1.0 - paa_pred + Type(1.0e-15));
    vector<Type> p0 = 1.0/(1.0 + exp(exp(age_comp_pars(1))*(X - age_comp_pars(0)))); //prob of zero declines with proportion caught
    Type sd = exp(age_comp_pars(2));
    for(int a = 0; a < n_ages; a++) obs(a) = rbinom(Type(1.0), Type(1.0) - p0(a)); // generate instances of positive observations
    int n_pos = 0;
    for(int a = 0; a < n_ages; a++) if(obs(a) > 0.5) n_pos++;
    if(n_pos>0)
    {
      vector<Type> pos_pred(n_pos);
      int k = 0;
      for(int a = 0; a < n_ages; a++) if(obs(a) > 0.5)
      {
        pos_pred(k) = paa_pred(a);
        k++;
      }
      vector<Type> pos_obs(n_pos);
      pos_obs.setZero();
      for(int a = 0; a < n_pos-1; a++) 
      {
        pos_obs(a) = exp(rnorm(log(pos_pred(a)) - log(pos_pred(n_pos-1)), sd));
      }
      pos_obs = pos_obs/(1.0 + pos_obs.sum());
      pos_obs(n_pos-1) = 1.0 - pos_obs.sum();
      k = 0;
      for(int a = 0; a < n_ages; a++) if(obs(a) > 0.5)
      {
        obs(a) = pos_obs(k);
        k++;
      }
    }
  }
  if(age_comp_model == 5) //logistic normal. Pool zero observations with adjacent age classes.
  {
    vector<Type> best_obs(n_ages);
    best_obs.setZero();
    Type sd = exp(age_comp_pars(0)-0.5*log(Neff));
    for(int a = 0; a < n_ages-1; a++) best_obs(a) = exp(rnorm(log(paa_pred(a)) - log(paa_pred(n_ages-1)), sd));
    best_obs = best_obs/(1.0 + best_obs.sum());
    best_obs(n_ages-1) = 1.0 - best_obs.sum();
    
    Type obs_2 = best_obs.tail(n_ages-aref+1).sum(); // .tail last n_ages-aref+1 components
    for(int a = aref-1; a < n_ages; a++) if(paa_obs(a) > Type(1.0e-15)) obs(a) = obs_2;
    obs_2 = 0.0;
    for(int a = 0; a < aref-1; a++) 
    {
      obs_2 += best_obs(a);
      if(paa_obs(a) > Type(1.0e-15))
      {
        obs(a) = obs_2;
        obs_2 = 0.0;
      }
      else obs(a) = 0.0;
      //else pooling with next age
    }
  }
  if(age_comp_model == 6) //zero-one inflated logistic normal where p0 is a function of binomial sample size. 2 parameters
  {
    Type n_e = exp(age_comp_pars(0));
    vector<Type> p0 = exp(n_e * log(1.0-paa_pred + Type(1.0e-15))); //prob of zero declines with proportion caught
    Type sd = exp(age_comp_pars(1));
    
    for(int a = 0; a < n_ages; a++) obs(a) = rbinom(Type(1.0), Type(1.0) - p0(a)); // generate instances of positive observations
    int n_pos = 0;
    for(int a = 0; a < n_ages; a++) if(obs(a) > 0.5) n_pos++;
    if(n_pos>0)
    {
      vector<Type> pos_pred(n_pos);
      int k = 0;
      for(int a = 0; a < n_ages; a++) if(obs(a) > 0.5)
      {
        pos_pred(k) = paa_pred(a);
        k++;
      }
      vector<Type> pos_obs(n_pos);
      pos_obs.setZero();
      for(int a = 0; a < n_pos-1; a++) 
      {
        pos_obs(a) = exp(rnorm(log(pos_pred(a)) - log(pos_pred(n_pos-1)), sd));
      }
      pos_obs = pos_obs/(1.0 + pos_obs.sum());
      pos_obs(n_pos-1) = 1.0 - pos_obs.sum();
      k = 0;
      for(int a = 0; a < n_ages; a++) if(obs(a) > 0.5)
      {
        obs(a) = pos_obs(k);
        k++;
      }
    }
  }
  if(age_comp_model == 7) //logistic normal treating 0 observations as missing. One parameter.
  {
    Type sd = exp(age_comp_pars(0));
    
    int n_pos = 0;
    for(int a = 0; a < n_ages; a++) if(paa_obs(a) > 1.0e-15) n_pos++;
    if(n_pos>0)
    {
      vector<Type> pos_pred(n_pos);
      int k = 0;
      for(int a = 0; a < n_ages; a++) if(paa_obs(a) > 1.0e-15)
      {
        pos_pred(k) = paa_pred(a);
        k++;
      }
      vector<Type> pos_obs(n_pos);
      pos_obs.setZero();
      for(int a = 0; a < n_pos-1; a++) 
      {
        pos_obs(a) = exp(rnorm(log(pos_pred(a)) - log(pos_pred(n_pos-1)), sd));
      }
      pos_obs = pos_obs/(1.0 + pos_obs.sum());
      pos_obs(n_pos-1) = 1.0 - pos_obs.sum();
      k = 0;
      for(int a = 0; a < n_ages; a++) if(paa_obs(a) > 1.0e-15)
      {
        obs(a) = pos_obs(k);
        k++;
      }
    }
  }
  return obs;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(n_years);
  DATA_INTEGER(n_ages);
  DATA_INTEGER(n_fleets);
  DATA_INTEGER(n_indices);
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
  DATA_IVECTOR(Ecov_maxages_k); //last year of Ecov to use (nobs)
  //DATA_IVECTOR(Ecov_maxages_a50); //last year of Ecov to use (nobs)
  DATA_VECTOR(Ecov_obs);
  DATA_IVECTOR(use_Ecov_obs);
  DATA_VECTOR(Ecov_obs_sigma);
  DATA_INTEGER(fit_k);
  DATA_INTEGER(fit_a50);
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
  PARAMETER(beta_Linf);
  PARAMETER(t0);
  PARAMETER_VECTOR(beta_Ecov_b);
  PARAMETER_VECTOR(beta_Ecov_a);
  PARAMETER_VECTOR(beta_Ecov_k_LVB);
  PARAMETER_VECTOR(beta_Ecov_Linf);
  PARAMETER_VECTOR(k_LVB_AR_pars); //AR1 process (2)
  PARAMETER_VECTOR(k_LVB_re);
  PARAMETER(beta_sig_Linf);
  PARAMETER(beta_sig_L_obs);
  PARAMETER(beta_sig_W_obs);
  
  //Environmental covariate parameters
  PARAMETER(Ecov_mu);
  PARAMETER_VECTOR(Ecov_AR_pars); //AR1 process (2)
  PARAMETER_VECTOR(Ecov_re);

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
  vector<Type> SSB(n_years);
  matrix<Type> F(n_years,n_fleets);
  matrix<Type> log_F(n_years,n_fleets);
  array<Type> pred_CAA(n_years,n_fleets,n_ages);
  array<Type> pred_catch_paa(n_years,n_fleets,n_ages);
  matrix<Type> pred_catch(n_years,n_fleets);
  array<Type> pred_IAA(n_years,n_indices,n_ages);
  array<Type> pred_index_paa(n_years,n_indices,n_ages);
  matrix<Type> pred_indices(n_years,n_indices);
  matrix<Type> log_pred_catch(n_years,n_fleets);
  matrix<Type> NAA(n_years,n_ages);
  matrix<Type> pred_NAA(n_years,n_ages);
  array<Type> FAA(n_years,n_fleets,n_ages);
  matrix<Type> FAA_tot(n_years,n_ages);
  matrix<Type> ZAA(n_years,n_ages);
  array<Type> QAA(n_years,n_indices,n_ages);
  matrix<Type> selblocks(n_selblocks,n_ages);
  vector<Type> q(n_indices);
  Type nll = zero; //negative log-likelihood
  vector<Type> t_paa(n_ages), t_pred_paa(n_ages);
  
  //do Ecov process and observations. Fill out Ecov_out to use later
  int n_obs_mat = Y.size();
  int n_obs = weight.size();
  int n_years_Ecov = Ecov_obs.rows();// + max_ecov_ages;
    
  vector<Type> Ecov_y(n_years_Ecov);
  Type Ecov_phi;
  Type Ecov_sig;
  Type nll_Ecov = zero;
  Ecov_phi = -one + Type(2)/(one + exp(-Ecov_AR_pars(0)));
  Ecov_sig = exp(Ecov_AR_pars(1));
  for(int y = 0; y < n_years_Ecov; y++) Ecov_y(y) = Ecov_mu + Ecov_re(y);
  nll_Ecov -= dnorm(Ecov_re(0), zero, Ecov_sig*exp(-Type(0.5) * log(one - pow(Ecov_phi,Type(2)))), 1);
  for(int y = 1; y < n_years_Ecov; y++) nll -= dnorm(Ecov_re(y), Ecov_phi * Ecov_re(y-1), Ecov_sig, 1);
  nll += nll_Ecov;
  see(nll_Ecov);
  SIMULATE {
    Ecov_re(0) = rnorm(zero, Ecov_sig*exp(-Type(0.5) * log(one - pow(Ecov_phi,Type(2)))));
    for(int y = 1; y < n_years_Ecov; y++) Ecov_re(y) = rnorm(Ecov_phi * Ecov_re(y-1), Ecov_sig);
    for(int y = 0; y < n_years_Ecov; y++) Ecov_y(y) = Ecov_mu + Ecov_re(y);
  }
  Type nll_Ecov_obs = zero;
  for(int y = 0; y < Ecov_obs.size(); y++) if(use_Ecov_obs(y) == 1) nll_Ecov_obs -= dnorm(Ecov_obs(y), Ecov_y(y), Ecov_obs_sigma(y), 1);
  nll += nll_Ecov_obs;
  see(nll_Ecov_obs);
  SIMULATE {
    for(int y = 0; y < Ecov_obs.size(); y++) if(use_Ecov_obs(y) == 1) Ecov_obs(y) = rnorm(Ecov_y(y),Ecov_obs_sigma(y));
  }
 /*
  //do fit for maturity data to get parameters to specify maturity at age for SSB
  //k AR(1) process
  Type k_phi;
  Type k_sig;
  Type nll_k = zero;
  k_phi = -one + Type(2)/(one + exp(-k_AR_pars(0)));
  k_sig = exp(k_AR_pars(1));
  nll_k -= dnorm(k_re(0), zero, k_sig*exp(-Type(0.5) * log(one - pow(k_phi,Type(2)))), 1);
  for(int y = 1; y < n_years_Ecov; y++) nll_k -= dnorm(k_re(y), k_phi * k_re(y-1), k_sig, 1);
  if(fit_k == 1) nll += nll_k;

  //a50 AR(1) process
  Type a50_phi;
  Type a50_sig;
  Type nll_a50 = zero;
  a50_phi = -one + Type(2)/(one + exp(-a50_AR_pars(0)));
  a50_sig = exp(a50_AR_pars(1));
  nll_a50 -= dnorm(a50_re(0), zero, a50_sig*exp(-Type(0.5) * log(one - pow(a50_phi,Type(2)))), 1);
  for(int y = 1; y < n_years_Ecov; y++) nll_a50 -= dnorm(a50_re(y), a50_phi * a50_re(y-1), a50_sig, 1);
  if(fit_a50 == 1) nll += nll_a50;
  
  vector<Type> l_k(n_obs_mat);
  vector<Type> l_a50(n_obs_mat);
  vector<Type> logit_mat(n_obs_mat);
  for(int i = 0; i < n_obs_mat; i++) 
  {
    l_k(i) = beta_k; 
    for(int y = 0; y <= Ecov_maxages_k(i); y++) l_k(i) += Ecov_y(year_obs(i) - 1 - age_obs(i) + y) * beta_Ecov_k(y) + k_re(year_obs(i) -1 - age_obs(i) + y);
    l_a50(i) = beta_a50;
    for(int y = 0; y <= Ecov_maxages_a50(i); y++) l_a50(i) += Ecov_y(year_obs(i) - 1 - age_obs(i) + y) * beta_Ecov_a50(y) + a50_re(year_obs(i) -1 -age_obs(i) + y);
    logit_mat(i) = exp(l_k(i)) * (Type(age_obs(i)) - exp(l_a50(i)));
  }
  vector<Type> mat = one/(one + exp(-logit_mat));
  Type phi = exp(beta_phi);
  vector<Type> ll_mat(n_obs_mat); 
  ll_mat.setZero();
  if(binomial == 1) 
  {
    ll_mat = dbinom(Y, N, mat, 1); //binomial?
  }
  else 
  {
    for(int i = 0; i < n_obs_mat; i++) ll_mat(i) = dbetabinom(Y(i), N(i), mat(i), phi, 1);
  }
  //see(sum(ll_mat));
  nll -= sum(ll_mat);
  see(n_obs_mat);

  //get maturity at age
  matrix<Type> logit_pmat(n_years,n_ages), pmat(n_years,n_ages);
  vector<Type> fix_pmat(10);
  for(int y = 0; y < n_years; y++) for(int a = 0; a < n_ages; a++) 
  {
    Type temp = beta_k;
    int maxage = a+1;
    if(a+1 > beta_Ecov_k.size()) maxage = beta_Ecov_k.size();
    for(int i = 0; i < maxage; i++) temp += beta_Ecov_k(i) * Ecov_y(model_years(y)-1 - a + i - 1) + k_re(model_years(y) - 1 - a + i - 1);
    logit_pmat(y,a) = exp(temp);
    if(y == n_years-1 & a == 0) fix_pmat(8) = exp(temp);
    temp = beta_a50;
    maxage = a+1;
    if(a+1 > beta_Ecov_a50.size()) maxage = beta_Ecov_a50.size();
    for(int i = 0; i < maxage; i++) temp += beta_Ecov_a50(i) * Ecov_y(model_years(y)-1 - a + i - 1) + a50_re(model_years(y) - 1 - a + i - 1);
    if(y == n_years-1 & a == 0) fix_pmat(9) = exp(temp);
    logit_pmat(y,a) *= Type(a+1) - exp(temp);
    pmat(y,a) = one/(one + exp(-logit_pmat(y,a)));
    if(y == n_years-1 & a == 0) 
    {
      fix_pmat(0) = beta_k;
      fix_pmat(1) = beta_Ecov_k(0);
      fix_pmat(2) = Ecov_y(model_years(y)-1 - a + 0 - 1);
      fix_pmat(3) = k_re(model_years(y) - 1 - a + 0 - 1);
      fix_pmat(4) = beta_a50;
      fix_pmat(5) = beta_Ecov_a50(0);
      fix_pmat(6) = a50_re(model_years(y) - 1 - a + 0 - 1);
      fix_pmat(7) = logit_pmat(y,a);
    }
  }
*/
  //do fit for growth data to get parameters to specify weight at age for SSB
  //k AR(1) process
  Type k_LVB_phi;
  Type k_LVB_sig;
  Type nll_k_LVB = zero;
  k_LVB_phi = -one + Type(2)/(one + exp(-k_LVB_AR_pars(0)));
  k_LVB_sig = exp(k_LVB_AR_pars(1));
  nll_k_LVB -= dnorm(k_LVB_re(0), zero, k_LVB_sig*exp(-Type(0.5) * log(one - pow(k_LVB_phi,Type(2)))), 1);
  for(int y = 1; y < n_years_Ecov; y++) nll_k_LVB -= dnorm(k_LVB_re(y), k_LVB_phi * k_LVB_re(y-1), k_LVB_sig, 1);
  SIMULATE {
    k_LVB_re(0) = rnorm(zero, k_LVB_sig*exp(-Type(0.5) * log(one - pow(k_LVB_phi,Type(2)))));
    for(int y = 1; y < n_years_Ecov; y++) k_LVB_re(y) = rnorm(k_LVB_phi * k_LVB_re(y-1), k_LVB_sig);
  }
  if(fit_k_LVB == 1) nll += nll_k_LVB;
  
  vector<Type> l_b(n_obs), l_a(n_obs), l_Linf(n_obs), log_l_pred(n_obs), log_w_pred(n_obs), ll_growth(n_obs), v_log_w(n_obs);
  vector<Type> k_sum(n_obs);
  k_sum.setZero(); ll_growth.setZero();
  Type v_log_l = exp(Type(2) * beta_sig_L_obs) + exp(Type(2) * beta_sig_Linf);
  for(int i = 0; i < n_obs; i++) 
  {
    l_b(i) = beta_b; 
    for(int y = 0; y <= Ecov_maxages_b(i); y++) l_b(i) += Ecov_y(year_obs_growth(i) - 1 - age_obs_growth(i) + y) * beta_Ecov_b(y);
    l_a(i) = beta_a;// - Type(11);
    for(int y = 0; y <= Ecov_maxages_a(i); y++) l_a(i) += Ecov_y(year_obs_growth(i) - 1 - age_obs_growth(i) + y) * beta_Ecov_a(y);
    //do LVB k something like Millar et al. 1999. Need to specify length of beta_Ecov_k = maximum observed age and specify map to fix older ages at 0
    Type k0 = exp(beta_k_LVB(0) + k_LVB_re(year_obs_growth(i) - 1 - age_obs_growth(i) + 0) + Ecov_y(year_obs_growth(i) - 1 - age_obs_growth(i) + 0) * beta_Ecov_k_LVB(0));
    k_sum(i) = k0;
    for(int y = 1; y < age_obs_growth(i); y++) k_sum(i) += exp(beta_k_LVB(y) + k_LVB_re(year_obs_growth(i) - 1 - age_obs_growth(i) + y) + Ecov_y(year_obs_growth(i) - 1 - age_obs_growth(i) + y) * beta_Ecov_k_LVB(y));
    l_Linf(i) = beta_Linf;
    for(int y = 0; y <= Ecov_maxages_Linf(i); y++) l_Linf(i) += Ecov_y(year_obs_growth(i) - 1 - age_obs_growth(i) + y) * beta_Ecov_Linf(y);
    log_l_pred(i) = l_Linf(i) + log(Type(1) - exp(t0*k0 - k_sum(i)));
    log_w_pred(i) = l_a(i) + exp(l_b(i)) * log_l_pred(i);
    v_log_w(i) = exp(Type(2) * beta_sig_W_obs) + exp(Type(2) * (l_b(i) + beta_sig_Linf));
    matrix<Type> S(2,2);
    S(0,0) = v_log_l; S(1,1) = v_log_w(i);
    S(0,1) = S(1,0) = exp(l_b(i) + Type(2) * beta_sig_Linf);
    vector<Type> M(2); M(0) = log_l_pred(i); M(1) = log_w_pred(i); 
    vector<Type> obs(2); obs(0) = log(len(i)); obs(1) = log(weight(i));
    if(islen(i) == 1 & iswt(i) == 1) {
      ll_growth(i) += dbinorm(obs,M,S,1);
      SIMULATE obs = rbinorm(M,S);
    }
    else
    {
      if(islen(i) == 1) {
        ll_growth(i) += dnorm(obs(0), M(0), exp(Type(0.5)*log(S(0,0))), 1);
        SIMULATE obs(0) = rnorm(M(0), exp(Type(0.5)*log(S(0,0))));
      }
      if(iswt(i) == 1) {
        ll_growth(i) += dnorm(obs(1), M(1), exp(Type(0.5)*log(S(1,1))), 1);
        SIMULATE obs(1) = rnorm(M(1), exp(Type(0.5)*log(S(1,1))));
    }
  }
  //see(sum(ll_growth));
  nll -= sum(ll_growth);

  //get weight at age
  //int n_ages = beta_Ecov_k_LVB.size();
  //int n_years_model = n_years_Ecov-n_ages;
  matrix<Type> k_age(n_years,n_ages), log_Linf_age(n_years,n_ages), log_a_age(n_years,n_ages), log_b_age(n_years,n_ages);
  matrix<Type> log_laa(n_years,n_ages), log_waa(n_years,n_ages);
  k_age.setZero(); log_Linf_age.setZero(); log_a_age.setZero(); log_b_age.setZero(); log_laa.setZero(); log_waa.setZero();
  vector<Type> fix(6);
  for(int y = 0; y < n_years; y++) 
  {
    for(int a = 0; a < n_ages; a++) 
    {
      log_Linf_age(y,a) = beta_Linf;
      log_a_age(y,a) = beta_a;
      log_b_age(y,a) = beta_b;
      int maxage_a = a + 1;
      if(a + 1 > beta_Ecov_a.size()) maxage_a = beta_Ecov_a.size();
      for(int i = 0; i < maxage_a; i++) log_a_age(y,a) += beta_Ecov_a(i) * Ecov_y(model_years(y)-1 - a + i - 1);
      int maxage_b = a + 1;
      if(a + 1 > beta_Ecov_b.size()) maxage_b = beta_Ecov_b.size();
      for(int i = 0; i < maxage_b; i++) log_b_age(y,a) += beta_Ecov_b(i) * Ecov_y(model_years(y) - 1 - a + i - 1);
      int maxage_Linf = a + 1;
      if(a + 1 > beta_Ecov_Linf.size()) maxage_Linf = beta_Ecov_Linf.size();
      for(int i = 0; i < maxage_Linf; i++) log_Linf_age(y,a) += beta_Ecov_Linf(i) * Ecov_y(model_years(y) - 1 - a + i - 1);
      Type k0 = exp(beta_k_LVB(0) + k_LVB_re(model_years(y) - 1 - a + 0 - 1) + Ecov_y(model_years(y) - 1 - a + 0 - 1) * beta_Ecov_k_LVB(0));
      k_age(y,a) = k0;
      for(int i = 1; i < a + 1; i++) k_age(y,a) += exp(beta_k_LVB(i) + k_LVB_re(model_years(y) - 1 - a + i - 1) + beta_Ecov_k_LVB(i) * Ecov_y(model_years(y) - 1 - a + i - 1));
      log_laa(y,a) = log_Linf_age(y,a) + log(one - exp(t0 * k0 - k_age(y,a)));
      log_waa(y,a) = log_a_age(y,a) + exp(log_b_age(y,a)) * log_laa(y,a);
      if(y==n_years-1 & a == 0)
      {
          fix(0) = k_LVB_re(model_years(y) - 1 - a + 0 - 1);
          fix(1) = Ecov_y(model_years(y) - 1 - a + 0 - 1);
          fix(2) = beta_k_LVB(0);
          fix(3) = beta_Ecov_k_LVB(0);
          fix(4) = k0;
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
  for(int y = 0; y < n_years; y++) SSB(y) = zero;

  ZAA = FAA_tot + Maa;
  for(int a = 0; a < n_ages; a++) NAA(0,a) = exp(log_N1(a));
  
  for(int a = 0; a < n_ages; a++) SSB(0) += NAA(0,a) * waa(waa_pointer_ssb-1,0,a) * mature(0,a) * exp(-ZAA(0,a) * fracyr_SSB(0));

  for(int y = 1; y < n_years; y++) 
  {
    for(int a = 0; a < n_ages; a++) 
    {
      NAA(y,a) = exp(log_NAA(y-1,a)); //random effects NAA
      SSB(y) += NAA(y,a) * waa(waa_pointer_ssb-1,y,a) * mature(y,a) * exp(-ZAA(y,a)*fracyr_SSB(y));
    }
  }
      
  for(int a = 0; a < n_ages; a++) pred_NAA(0,a) = exp(log_N1(a));
  
  Type nll_NAA = zero;
  for(int y = 1; y < n_years; y++)
  {
    if(recruit_model == 1) pred_NAA(y,0) = NAA(y-1,0); //random walkNAA(y,1)
    else
    {
      if(recruit_model == 2) pred_NAA(y,0) = exp(mean_rec_pars(0)); //random about mean
      else //BH stock recruit
      {
        if(recruit_model == 3) //BH stock recruit
        {
          pred_NAA(y,0) = exp(mean_rec_pars(0) + log(SSB(y-1)) - log(one + exp(mean_rec_pars(1)) * SSB(y-1)));
        }
        else //Ricker stock recruit
        {
          pred_NAA(y,0) = exp(mean_rec_pars(0) + log(SSB(y-1)) - exp(mean_rec_pars(1)) * SSB(y-1)); 
        }
      }
    }
    for(int a = 1; a < n_ages-1; a++) 
    {
      pred_NAA(y,a) = NAA(y-1,a-1) * exp(-ZAA(y-1,a-1));
    }
    pred_NAA(y,n_ages-1) = NAA(y-1,n_ages-2) * exp(-ZAA(y-1,n_ages-2)) + NAA(y-1,n_ages-1) * exp(-ZAA(y-1,n_ages-1));
  }
  
  for(int y = 1; y < n_years; y++)
  {
    for(int a = 0; a < n_ages; a++)
      nll_NAA -= dnorm(log(NAA(y,a)), log(pred_NAA(y,a)), exp(log_NAA_sigma(NAA_sigma_pointers(a)-1)), 1);
  }
  SIMULATE {
    for(int y = 1; y < n_years; y++)
    {
      for(int a = 0; a < n_ages; a++){
        log_NAA(y-1,a) = rnorm(log(pred_NAA(y,a)), exp(log_NAA_sigma(NAA_sigma_pointers(a)-1)), 1);
        NAA(y,a) = exp(log_NAA(y-1,a));
      }
    }
  }
  //see(nll_NAA);
  nll += nll_NAA;

  Type nll_agg_catch = zero, nll_catch_acomp = zero;
  for(int y = 0; y < n_years; y++)
  {
    int acomp_par_count = 0;
    for(int f = 0; f < n_fleets; f++)
    {
      pred_catch(y,f) = zero;
      Type tsum = zero;
      for(int a = 0; a < n_ages; a++) 
      {
        pred_CAA(y,f,a) =  NAA(y,a) * FAA(y,f,a) * (1 - exp(-ZAA(y,a)))/ZAA(y,a);
        pred_catch(y,f) += waa(waa_pointer_fleets(f)-1,y,a) * pred_CAA(y,f,a);
        tsum += pred_CAA(y,f,a);
      }
      nll_agg_catch -= dnorm(log(agg_catch(y,f)), log(pred_catch(y,f)), agg_catch_sigma(y,f), 1);
      SIMULATE agg_catch(y_f) = exp(rnorm(log(pred_catch(y,f)), agg_catch_sigma(y,f)));
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
          nll_catch_acomp -= get_acomp_ll(y, n_ages, catch_Neff(y,f), age_comp_model_fleets(f), t_paa, t_pred_paa, acomp_pars, catch_aref(y,f));
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
  nll += nll_catch_acomp;
  
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
        pred_IAA(y,i,a) =  NAA(y,a) * QAA(y,i,a) * exp(-ZAA(y,a) * fracyr_indices(y,i));
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
  ////////////////////////////////////////////////////////////////////
  
  
  //////////////////////////////////////////////////////////////////////
  vector<Type> log_SSB = log(SSB);
  vector<Type> SSB_E(n_years);
  SSB_E.setZero();
  for(int y = 0; y < n_years; y++) for(int a = 0; a < n_ages; a++) 
  {
    Type SSB_a = NAA(y,a) * exp(-ZAA(y,a)*fracyr_SSB(y));
    //if(use_mat_model == 1) SSB_a *= pmat(y,a);
    //else 
    SSB_a *= mature(y,a);
    if(use_growth_model == 1) SSB_a *= exp(log_waa(y,a));
    else SSB_a *= waa(waa_pointer_ssb-1,y,a);
    SSB_E(y) += SSB_a;
  }
  vector<Type> log_SSB_E = log(SSB_E);
  
  vector<Type> wssb(n_ages), wcatch(n_ages);
  vector<Type> SPR_0(n_years), SPR_0_E(n_years), SPR40(n_years), SPR40_E(n_years), log_F40(n_years), log_F40_E(n_years);
  vector<Type> log_SSB40(n_years), log_SSB40_E(n_years), YPR40(n_years), YPR40_E(n_years), log_Y40(n_years), log_Y40_E(n_years);
  vector<Type> msy_proxy(3);
  for(int a = 0; a < n_ages; a++) 
  {
    wssb(a) = waa(waa_pointer_ssb-1, n_years-1, a);
    wcatch(a) = waa(waa_pointer_totcatch-1, n_years-1, a);
  }
  vector<Type> M2 = Maa.row(n_years-1);
  vector<Type> sel_tot = FAA_tot.row(n_years-1)/FAA_tot(n_years-1,n_ages-1);
  for(int y = 0; y < n_years; y++)
  {
    vector<Type> maturity_E(n_ages);
    //if(use_mat_model == 1) maturity_E = pmat.row(y);
    //else  
    maturity_E = mature.row(y);
    vector<Type> wssb_E(n_ages);
    if(use_growth_model == 1) for(int a= 0; a < n_ages; a++) wssb_E(a) = exp(log_waa(y,a));
    else for(int a = 0; a < n_ages; a++) wssb_E(a) = waa(waa_pointer_ssb-1, y , a);
    vector<Type> maturity = mature.row(y);
    SPR_0(y) = get_SPR_0(M2, maturity, wssb, fracyr_SSB(n_years-1));
    msy_proxy = getFXSPR(SPR_0(y), percentSPR, M2, sel_tot, maturity, wssb, wcatch, fracyr_SSB(n_years-1), 10);
    log_F40(y) = log(msy_proxy(0));
    SPR40(y) = msy_proxy(1);
    YPR40(y) = msy_proxy(2);
    log_SSB40(y) = log(msy_proxy(1)*NAA.col(0).mean());
    log_Y40(y) = log(msy_proxy(2)*NAA.col(0).mean());
    SPR_0_E(y) = get_SPR_0(M2, maturity_E, wssb_E, fracyr_SSB(n_years-1));
    msy_proxy = getFXSPR(SPR_0_E(y), percentSPR, M2, sel_tot, maturity_E, wssb_E, wcatch, fracyr_SSB(n_years-1), 10);
    SPR40_E(y) = msy_proxy(1);
    YPR40_E(y) = msy_proxy(2);
    log_F40_E(y) = log(msy_proxy(0));
    log_SSB40_E(y) = log(msy_proxy(1)*NAA.col(0).mean());
    log_Y40_E(y) = log(msy_proxy(2)*NAA.col(0).mean());
  }
  //if(reportMode==0){
    REPORT(NAA);
    REPORT(pred_NAA);
    REPORT(SSB);
    REPORT(selblocks);
    REPORT(q);
    REPORT(F);
    //REPORT(pmat);
    ADREPORT(log_SSB);
    ADREPORT(log_SSB_E);
    ADREPORT(log_F);
    ADREPORT(SPR_0);
    ADREPORT(SPR40);
    ADREPORT(YPR40);
    ADREPORT(log_F40);
    ADREPORT(log_SSB40);
    ADREPORT(log_Y40);
    ADREPORT(SPR_0_E);
    ADREPORT(SPR40_E);
    ADREPORT(YPR40_E);
    ADREPORT(log_F40_E);
    ADREPORT(log_SSB40_E);
    ADREPORT(log_Y40_E);
    //ADREPORT(pmat);
  //}
  
  return nll;
}

