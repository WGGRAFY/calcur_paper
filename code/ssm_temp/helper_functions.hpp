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
        for (int a = 0; a < n_ages; a++) selectivity_blocks(i,a) = selectivity_blocks(i,a)/selectivity_blocks(i,n_ages-1);
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
    for(int a = 0; a < n_ages; a++) ll += -lgamma(temp_n(a) + one) + temp_n(a) * log(paa_pred(a)+1e-10);
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

/* calculate beverton-holt or ricker equilibrium yield */
template<class Type>
struct sr_yield {
  /* Data and parameter objects for yield calculation: */
  Type SR_a;
  Type SR_b;
  vector<Type> M;
  vector<Type> sel;
  vector<Type> mat;
  vector<Type> waassb;
  vector<Type> waacatch;
  Type fracyearSSB;
  int sr_type;

  /* Constructor */
  sr_yield(Type SR_a_, Type SR_b_,
  vector<Type> M_,
  vector<Type> sel_,
  vector<Type> mat_,
  vector<Type> waassb_,
  vector<Type> waacatch_,
  Type fracyearSSB_, int sr_type_) :
    SR_a(SR_a_), SR_b(SR_b_), M(M_), sel(sel_), mat(mat_),
  waassb(waassb_), waacatch(waacatch_), fracyearSSB(fracyearSSB_), sr_type(sr_type_) {}

  template <typename T> //I think this allows you to differentiate the function wrt whatever is after operator() on line below
  T operator()(vector<T> log_F) { //find such that it maximizes yield
    int n_ages = M.size();
    T YPR = 0, SPR = 0, ntemp = 1, R;
    vector<T> F(n_ages), Z(n_ages);

    F = exp(log_F(0)) * sel.template cast<T>();
    Z = F + M.template cast<T>();
    for(int age=0; age<n_ages-1; age++)
    {
      YPR += ntemp * F(age) * T(waacatch(age)) * (1- exp(-Z(age)))/Z(age);
      SPR += ntemp * T(mat(age) * waassb(age)) * exp(-T(fracyearSSB) * Z(age));
      ntemp *= exp(-Z(age));
    }
    ntemp /= 1 - exp(-Z(n_ages-1));
    YPR += ntemp * F(n_ages-1) * T(waacatch(n_ages-1)) * (1 - exp(-Z(n_ages-1)))/Z(n_ages-1);
    SPR += ntemp * T(mat(n_ages-1) * waassb(n_ages-1)) * exp(-T(fracyearSSB)*Z(n_ages-1));

    //Type SPR = get_SPR(x, M, sel, mat, waassb, fracyearSSB);
    if(sr_type == 0) R = (T(SR_a) - 1/SPR) / T(SR_b); //beverton-holt
    if(sr_type == 1) R = log(T(SR_a) * SPR)/(T(SR_b) * SPR); //ricker
    T Y = YPR * R;
    return Y;
  }
};

template <class Type>
Type get_FMSY(Type log_a, Type log_b, vector<Type> M, vector<Type> sel, vector<Type> waacatch, vector<Type> waassb,
  vector<Type> mat, Type fracyr_SSB, Type log_SPR0, int recruit_model, Type F_init)
{    
  int n = 10;
  vector<Type> log_FMSY_i(1);
  vector<Type> log_FMSY_iter(n);
  log_FMSY_iter.fill(log(F_init)); //starting value
  Type a = exp(log_a);
  Type b = exp(log_b);
  int sr_type = 0; //recruit_model = 3, B-H
  if(recruit_model == 4) sr_type = 1; //recruit_model = 4, Ricker
  sr_yield<Type> sryield(a, b, M, sel, mat, waassb, waacatch, fracyr_SSB, sr_type);
  for (int i=0; i<n-1; i++)
  {
    log_FMSY_i(0) = log_FMSY_iter(i);
    vector<Type> grad_sr_yield = autodiff::gradient(sryield,log_FMSY_i);
    matrix<Type> hess_sr_yield = autodiff::hessian(sryield,log_FMSY_i);
    log_FMSY_iter(i+1) = log_FMSY_iter(i) - grad_sr_yield(0)/hess_sr_yield(0,0);
  }
  Type FMSY = exp(log_FMSY_iter(n-1));
  return FMSY;
}
