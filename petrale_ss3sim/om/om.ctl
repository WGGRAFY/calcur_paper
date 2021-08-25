#V3.30.13-safe;_2019_03_09;_Stock_Synthesis_by_Richard_Methot_(NOAA)_using_ADMB_12.0
#Stock Synthesis (SS) is a work of the U.S. Government and is not subject to copyright protection in the United States.
#Foreign copyrights may apply. See copyright.txt for more information.
#_user_support_available_at:NMFS.Stock.Synthesis@noaa.gov
#_user_info_available_at:https://vlab.ncep.noaa.gov/group/stock-synthesis
#_data_and_control_files: 2019_petrale.dat // 2019_petrale.ctl
0  # 0 means do not read wtatage.ss; 1 means read and use wtatage.ss and also read and use growth parameters
1  #_N_Growth_Patterns
1 #_N_platoons_Within_GrowthPattern 
#_Cond 1 #_Morph_between/within_stdev_ratio (no read if N_morphs=1)
#_Cond  1 #vector_Morphdist_(-1_in_first_val_gives_normal_approx)
#
4 # recr_dist_method for parameters:  2=main effects for GP, Area, Settle timing; 3=each Settle entity; 4=none (only when N_GP*Nsettle*pop==1)
1 # not yet implemented; Future usage: Spawner-Recruitment: 1=global; 2=by area
1 #  number of recruitment settlement assignments 
0 # unused option
#GPattern month  area  age (for each settlement assignment)
 1 1 1 0
#
#_Cond 0 # N_movement_definitions goes here if Nareas > 1
#_Cond 1.0 # first age that moves (real age at begin of season, not integer) also cond on do_migration>0
#_Cond 1 1 1 2 4 10 # example move definition for seas=1, morph=1, source=1 dest=2, age1=4, age2=10
#
5 #_Nblock_Patterns
 5 3 3 1 1 #_blocks_per_pattern 
# begin and end years of blocks
 1973 1982 1983 1992 1993 2002 2003 2010 2011 2018
 2003 2009 2010 2010 2011 2018
 2003 2008 2009 2010 2011 2018
 1875 1875
 2004 2009
#
# controls for all timevary parameters 
1 #_env/block/dev_adjust_method for all time-vary parms (1=warn relative to base parm bounds; 3=no bound check)
#
# AUTOGEN
1 1 1 1 1 # autogen: 1st element for biology, 2nd for SR, 3rd for Q, 4th reserved, 5th for selex
# where: 0 = autogen all time-varying parms; 1 = read each time-varying parm line; 2 = read then autogen if parm min==-12345
#
#_Available timevary codes
#_Block types: 0: P_block=P_base*exp(TVP); 1: P_block=P_base+TVP; 2: P_block=TVP; 3: P_block=P_block(-1) + TVP
#_Block_trends: -1: trend bounded by base parm min-max and parms in transformed units (beware); -2: endtrend and infl_year direct values; -3: end and infl as fraction of base range
#_EnvLinks:  1: P(y)=P_base*exp(TVP*env(y));  2: P(y)=P_base+TVP*env(y);  3: null;  4: P(y)=2.0/(1.0+exp(-TVP1*env(y) - TVP2))
#_DevLinks:  1: P(y)*=exp(dev(y)*dev_se;  2: P(y)+=dev(y)*dev_se;  3: random walk;  4: zero-reverting random walk with rho;  21-24 keep last dev for rest of years
#
#
#
# setup for M, growth, maturity, fecundity, recruitment distibution, movement 
#
0 #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate
  #_no additional input for selected M option; read 1P per morph
#
1 # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K_incr; 4=age_specific_K_decr; 5=age_specific_K_each; 6=NA; 7=NA; 8=growth cessation
2 #_Age(post-settlement)_for_L1;linear growth below this
17 #_Growth_Age_for_L2 (999 to use as Linf)
-999 #_exponential decay for growth above maxage (value should approx initial Z; -999 replicates 3.24; -998 to not allow growth above maxage)
0  #_placeholder for future growth feature
#
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
0 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
#
1 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
3 #_First_Mature_Age
1 #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0 #_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
1 #_parameter_offset_approach (1=none, 2= M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)
#
#_growth_parms
#_ LO HI INIT PRIOR PR_SD PR_type PHASE env_var&link dev_link dev_minyr dev_maxyr dev_PH Block Block_Fxn
# Sex: 1  BioPattern: 1  NatMort
 0.005 0.5 0.158704 -1.7793 0.438 3 2 0 0 0 0 0.5 0 0 # NatM_p_1_Fem_GP_1
# Sex: 1  BioPattern: 1  Growth
 10 45 15.6515 17.18 10 0 3 0 0 0 0 0.5 0 0 # L_at_Amin_Fem_GP_1
 35 80 53.1167 54.2 10 0 3 0 0 0 0 0.5 0 0 # L_at_Amax_Fem_GP_1
 0.04 0.5 0.141731 0.157 0.8 0 3 0 0 0 0 0.5 0 0 # VonBert_K_Fem_GP_1
 0.01 1 0.186051 3 0.8 0 3 0 0 0 0 0.5 0 0 # CV_young_Fem_GP_1
 0.01 1 0.0351949 0 1 0 4 0 0 0 0 0 0 0 # CV_old_Fem_GP_1
# Sex: 1  BioPattern: 1  WtLen
 -3 3 1.986e-06 1.99e-06 0.8 6 -3 0 0 0 0 0.5 0 0 # Wtlen_1_Fem_GP_1
 1 5 3.484 3.478 0.8 6 -3 0 0 0 0 0.5 0 0 # Wtlen_2_Fem_GP_1
# Sex: 1  BioPattern: 1  Maturity&Fecundity
 10 50 33.1 33.1 0.8 6 -3 0 0 0 0 0.5 0 0 # Mat50%_Fem_GP_1
 -3 3 -0.743 -0.743 0.8 6 -3 0 0 0 0 0.5 0 0 # Mat_slope_Fem_GP_1
 -3 3 1 1 1 6 -3 0 0 0 0 0.5 0 0 # Eggs/kg_inter_Fem_GP_1
 -3 3 0 0 1 6 -3 0 0 0 0 0.5 0 0 # Eggs/kg_slope_wt_Fem_GP_1
# Sex: 2  BioPattern: 1  NatMort
 0.005 0.6 0.164428 -1.6809 0.438 3 2 0 0 0 0 0.5 0 0 # NatM_p_1_Mal_GP_1
# Sex: 2  BioPattern: 1  Growth
 10 45 16.1562 17.18 10 0 3 0 0 0 0 0.5 0 0 # L_at_Amin_Mal_GP_1
 35 80 40.8281 41.1 10 0 3 0 0 0 0 0.5 0 0 # L_at_Amax_Mal_GP_1
 0.04 0.5 0.238148 0.247 0.8 0 3 0 0 0 0 0.5 0 0 # VonBert_K_Mal_GP_1
 0.01 1 0.136371 3 0.8 0 3 0 0 0 0 0.5 0 0 # CV_young_Mal_GP_1
 0.01 1 0.0597473 0 1 0 4 0 0 0 0 0 0 0 # CV_old_Mal_GP_1
# Sex: 2  BioPattern: 1  WtLen
 -3 3 2.983e-06 2.98e-06 0.8 6 -3 0 0 0 0 0.5 0 0 # Wtlen_1_Mal_GP_1
 -3 5 3.363 3.363 0.8 6 -3 0 0 0 0 0.5 0 0 # Wtlen_2_Mal_GP_1
# Hermaphroditism
#  Recruitment Distribution  
#  Cohort growth dev base
 0 1 1 1 0 0 -4 0 0 0 0 0 0 0 # CohortGrowDev
#  Movement
#  Age Error from parameters
#  catch multiplier
#  fraction female, by GP
 0.01 0.99 0.5 0.5 0.5 0 -99 0 0 0 0 0 0 0 # FracFemale_GP_1
#
#_no timevary MG parameters
#
#_seasonal_effects_on_biology_parms
 0 0 0 0 0 0 0 0 0 0 #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K
#_ LO HI INIT PRIOR PR_SD PR_type PHASE
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no seasonal MG parameters
#
3 #_Spawner-Recruitment; Options: 2=Ricker; 3=std_B-H; 4=SCAA; 5=Hockey; 6=B-H_flattop; 7=survival_3Parm; 8=Shepherd_3Parm; 9=RickerPower_3parm
0  # 0/1 to use steepness in initial equ recruitment calculation
0  #  future feature:  0/1 to make realized sigmaR a function of SR curvature
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn #  parm_name
             5            20       9.92138             9            10             0          1          0          0          0          0          0          0          0 # SR_LN(R0)
           0.2             1      0.841493           0.8          0.09             6          5          0          0          0          0          0          0          0 # SR_BH_steep
             0             2           0.4           0.9             5             6        -99          0          0          0          0          0          0          0 # SR_sigmaR
            -5             5             0             0           0.2             6         -2          0          0          0          0          0          0          0 # SR_regime
             0             0             0             0             0             0        -99          0          0          0          0          0          0          0 # SR_autocorr
1 #do_recdev:  0=none; 1=devvector (R=F(SSB)+dev); 2=deviations (R=F(SSB)+dev); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
1959 # first year of main recr_devs; early devs can preceed this era
2016 # last year of main recr_devs; forecast devs start in following year
1 #_recdev phase 
1 # (0/1) to read 13 advanced options
 0 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
 3 #_recdev_early_phase
 0 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
 1 #_lambda for Fcast_recr_like occurring before endyr+1
 1953 #_last_yr_nobias_adj_in_MPD; begin of ramp
 1964 #_first_yr_fullbias_adj_in_MPD; begin of plateau
 2015 #_last_yr_fullbias_adj_in_MPD
 2018 #_end_yr_for_ramp_in_MPD (can be in forecast to shape ramp, but SS sets bias_adj to 0.0 for fcast yrs)
 0.7647 #_max_bias_adj_in_MPD (-1 to override ramp and set biasadj=1.0 for all estimated recdevs)
 0 #_period of cycles in recruitment (N parms read below)
 -4 #min rec_dev
 4 #max rec_dev
 0 #_read_recdevs
#_end of advanced SR options
#
#_placeholder for full parameter lines for recruitment cycles
# read specified recr devs
#_Yr Input_value
#
# all recruitment deviations
#  1845E 1846E 1847E 1848E 1849E 1850E 1851E 1852E 1853E 1854E 1855E 1856E 1857E 1858E 1859E 1860E 1861E 1862E 1863E 1864E 1865E 1866E 1867E 1868E 1869E 1870E 1871E 1872E 1873E 1874E 1875E 1876E 1877E 1878E 1879E 1880E 1881E 1882E 1883E 1884E 1885E 1886E 1887E 1888E 1889E 1890E 1891E 1892E 1893E 1894E 1895E 1896E 1897E 1898E 1899E 1900E 1901E 1902E 1903E 1904E 1905E 1906E 1907E 1908E 1909E 1910E 1911E 1912E 1913E 1914E 1915E 1916E 1917E 1918E 1919E 1920E 1921E 1922E 1923E 1924E 1925E 1926E 1927E 1928E 1929E 1930E 1931E 1932E 1933E 1934E 1935E 1936E 1937E 1938E 1939E 1940E 1941E 1942E 1943E 1944E 1945E 1946E 1947E 1948E 1949E 1950E 1951E 1952E 1953E 1954E 1955E 1956E 1957E 1958E 1959R 1960R 1961R 1962R 1963R 1964R 1965R 1966R 1967R 1968R 1969R 1970R 1971R 1972R 1973R 1974R 1975R 1976R 1977R 1978R 1979R 1980R 1981R 1982R 1983R 1984R 1985R 1986R 1987R 1988R 1989R 1990R 1991R 1992R 1993R 1994R 1995R 1996R 1997R 1998R 1999R 2000R 2001R 2002R 2003R 2004R 2005R 2006R 2007R 2008R 2009R 2010R 2011R 2012R 2013R 2014R 2015R 2016R 2017F 2018F 2019F 2020F 2021F 2022F 2023F 2024F 2025F 2026F 2027F 2028F 2029F 2030F
#  1.94064e-07 2.2766e-07 2.6352e-07 3.11448e-07 3.63083e-07 4.20272e-07 4.95381e-07 5.76002e-07 6.73192e-07 7.87185e-07 9.19947e-07 1.07607e-06 1.25601e-06 1.46407e-06 1.70636e-06 1.98596e-06 2.30604e-06 2.68141e-06 3.10519e-06 3.59611e-06 4.15824e-06 4.7992e-06 5.5294e-06 6.35756e-06 7.28903e-06 8.33796e-06 9.51955e-06 1.08453e-05 1.23508e-05 1.40634e-05 1.60097e-05 1.82211e-05 2.07391e-05 2.35948e-05 2.68407e-05 3.05289e-05 3.47201e-05 3.94783e-05 4.48813e-05 5.10066e-05 5.79606e-05 6.5859e-05 7.48064e-05 8.49618e-05 9.64766e-05 0.000109523 0.000124318 0.000141083 0.000160081 0.000181601 0.00020597 0.000233564 0.000264801 0.000300157 0.000340159 0.00038541 0.000436599 0.000494498 0.000559974 0.00063401 0.000717725 0.000812361 0.000919314 0.0010402 0.00117592 0.00132873 0.00150129 0.00169564 0.00191435 0.00216051 0.00243607 0.00274094 0.00307603 0.00344411 0.00385133 0.00430066 0.00479361 0.00533061 0.00590979 0.00653366 0.00720382 0.00792259 0.00869292 0.00951447 0.0104021 0.0113927 0.0126171 0.0143244 0.0169103 0.0209393 0.0268821 0.0348427 0.0440279 0.0516114 0.0514958 0.0353084 -0.00141791 -0.0523266 -0.0977662 -0.110361 -0.0792079 -0.0395338 -0.0329372 -0.0391267 -0.0418615 -0.0384161 -0.0279741 -0.0148726 -0.0178005 -0.0431722 -0.0819993 -0.127811 -0.187962 -0.221907 -0.194692 -0.00332773 0.200232 -0.234196 -0.301866 0.170322 -0.0844305 0.69661 -0.131004 -0.0652386 0.0491235 0.175753 0.0632826 -0.162676 -0.2646 -0.0461246 -0.0145873 0.275821 0.361202 0.163352 -0.085461 -0.0206643 -0.163668 -0.175993 -0.0520474 0.301509 -0.178777 -0.603732 -0.444628 0.0125347 0.324476 0.291344 -0.158251 -0.574091 0.0585736 0.240043 -0.3434 -0.228037 -0.22206 0.654334 0.216223 -0.14471 -0.221353 -0.142628 -0.38465 -0.163691 -0.0859676 0.517299 0.722214 1.00819 0.123725 -0.134388 -0.00226759 0.338751 -0.238828 -0.260996 -0.32955 -0.102334 -0.121855 0 0 0 0 0 0 0 0 0 0 0 0 0
# implementation error by year in forecast:  0 0 0 0 0 0 0 0 0 0 0 0
#
#Fishing Mortality info 
0.3 # F ballpark
-2001 # F ballpark year (neg value to disable)
2 # F_Method:  1=Pope; 2=instan. F; 3=hybrid (hybrid is recommended)
4 # max F or harvest rate, depends on F_Method
# no additional F input needed for Fmethod 1
# if Fmethod=2; read overall start F value; overall phase; N detailed inputs to read
0 1 1
1 1 1 0 0.01 -1
# if Fmethod=3; read N iterations for tuning for Fmethod 3
# N iterations for tuning F in hybrid method (recommend 3 to 7)
#
#_initial_F_parms; count = 0
#_ LO HI INIT PRIOR PR_SD  PR_type  PHASE
#2030 2039
# F rates by fleet
# Yr:  1876 1877 1878 1879 1880 1881 1882 1883 1884 1885 1886 1887 1888 1889 1890 1891 1892 1893 1894 1895 1896 1897 1898 1899 1900 1901 1902 1903 1904 1905 1906 1907 1908 1909 1910 1911 1912 1913 1914 1915 1916 1917 1918 1919 1920 1921 1922 1923 1924 1925 1926 1927 1928 1929 1930 1931 1932 1933 1934 1935 1936 1937 1938 1939 1940 1941 1942 1943 1944 1945 1946 1947 1948 1949 1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960 1961 1962 1963 1964 1965 1966 1967 1968 1969 1970 1971 1972 1973 1974 1975 1976 1977 1978 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 2023 2024 2025 2026 2027 2028 2029 2030
# seas:  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
# WinterN 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5.49145e-05 0.00016578 0.000280765 0.000402834 0.000469745 0.000597679 0.000858911 0.00114584 0.00132966 0.00146259 0.00170029 0.00198329 0.00226609 0.00252683 0.00284836 0.00323548 0.00372797 0.00429842 0.010064 0.0084568 0.00993236 0.00358524 0.00210624 0.00463041 0.011155 0.0140348 0.0083634 0.0284518 0.021005 0.0199208 0.029362 0.0724671 0.0335015 0.0462066 0.0470142 0.0537254 0.0379522 0.0254142 0.055817 0.103112 0.104276 0.0938883 0.103164 0.12047 0.125116 0.180812 0.223631 0.277584 0.181955 0.427549 0.677685 0.249089 0.14687 0.132096 0.160024 0.316204 0.306835 0.306695 0.297601 0.411819 0.376121 0.486733 0.31687 0.287473 0.248165 0.252743 0.208669 0.172983 0.217645 0.245617 0.22999 0.195655 0.308825 0.276189 0.15084 0.259645 0.243178 0.254217 0.082945 0.0492573 0.0663447 0.0627274 0.0839112 0.0884789 0.067289 0.0852698 0.0711812 0.0783513 0.0806309 0.116077 0.115393 0.114866 0.114346 0.113682 0.113147 0.1126 0.111909 0.111353 0.110797
# SummerN 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5.51089e-06 4.52372e-06 3.53258e-06 3.44377e-06 3.35578e-06 3.26748e-06 3.17886e-06 3.08876e-06 2.99944e-06 2.90975e-06 2.81853e-06 2.72807e-06 2.63721e-06 2.5448e-06 2.45312e-06 2.36101e-06 2.26847e-06 2.17434e-06 2.08092e-06 1.98706e-06 1.89148e-06 1.79893e-06 1.70603e-06 1.60866e-06 1.50632e-06 1.40441e-06 1.30531e-06 1.2066e-06 1.10986e-06 1.01354e-06 9.14989e-07 8.579e-07 1.21527e-08 3.81032e-05 3.04216e-05 0.00203324 0.00630043 0.0103328 0.0145917 0.0170999 0.020645 0.0288024 0.0335016 0.0452832 0.053066 0.0570116 0.0961928 0.0996558 0.0741698 0.0693676 0.0967564 0.0730682 0.102618 0.0893564 0.1297 0.102423 0.0894884 0.048314 0.0600735 0.0584538 0.051525 0.075535 0.0808788 0.0674592 0.0996799 0.11832 0.125183 0.0940048 0.106277 0.0960148 0.0997222 0.0830641 0.0862264 0.0942274 0.114882 0.0940842 0.110807 0.162951 0.219222 0.236935 0.159566 0.134169 0.223865 0.272573 0.294122 0.156417 0.278466 0.282116 0.221631 0.160144 0.209061 0.183261 0.215006 0.215914 0.179264 0.174457 0.181272 0.165379 0.124603 0.130723 0.116235 0.115156 0.153516 0.125983 0.141185 0.143465 0.161922 0.160094 0.145027 0.200304 0.192503 0.104314 0.0717901 0.141252 0.0572129 0.0782819 0.0667132 0.108373 0.0737225 0.0819372 0.0830965 0.0879077 0.0876739 0.0925621 0.0875968 0.135825 0.13504 0.134435 0.133835 0.133064 0.13244 0.131804 0.130999 0.130352 0.129705
# WinterS 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.00144909 0.000836522 0.000892311 0.00327311 0.00373017 0.00233311 0.00185637 0.00122544 0.000817684 0.00445464 0.00314428 0.000716926 0.00221101 0.00273854 0.00333493 0.00245864 0.00555732 0.0106197 0.0266355 0.0204948 0.0109374 0.0177205 0.0293022 0.0445092 0.0297232 0.0287557 0.0300531 0.0331528 0.0224613 0.0313945 0.0345793 0.0355137 0.0400043 0.0299936 0.0301791 0.0264822 0.0444187 0.0305284 0.0352051 0.038514 0.0437929 0.058801 0.0342413 0.0578501 0.0540116 0.0761933 0.0563155 0.0542492 0.0531607 0.0860975 0.0619332 0.047024 0.0714264 0.0681174 0.0695355 0.10396 0.0872893 0.106241 0.12574 0.109038 0.172445 0.115527 0.155597 0.137477 0.105049 0.158226 0.17095 0.0830033 0.100704 0.125102 0.10336 0.1086 0.0536604 0.0324344 0.0611375 0.0226863 0.0748201 0.100117 0.0924273 0.0143949 0.00567389 0.013046 0.0105683 0.0188183 0.0135259 0.0141844 0.0118104 0.012991 0.0086508 0.0155882 0.0275251 0.0273671 0.0272454 0.0271239 0.0269675 0.0268408 0.0267116 0.0265487 0.0264178 0.0262871
# SummerS 2.02703e-05 2.02707e-05 2.0271e-05 2.02713e-05 0.000234161 0.000448189 0.000662477 0.000877137 0.00109227 0.00130795 0.00152426 0.00174127 0.00195903 0.00217759 0.00239698 0.00261724 0.00283839 0.00306068 0.0032837 0.00350768 0.00373265 0.00395862 0.00418558 0.00441355 0.00464254 0.00487255 0.0051036 0.00533569 0.00556882 0.005803 0.00603824 0.00627453 0.00651189 0.00675031 0.00698979 0.00723033 0.00747194 0.00771461 0.00795833 0.00820311 0.00834186 0.0114028 0.00921421 0.00725805 0.00501225 0.00637825 0.0092269 0.00929946 0.0116275 0.0115744 0.0114634 0.0139458 0.0137482 0.0157356 0.0147602 0.0119461 0.0117636 0.00893063 0.0207188 0.0183693 0.0103803 0.0181756 0.0224402 0.0268685 0.0161804 0.00925686 0.00624817 0.0104807 0.0141021 0.0141913 0.0375913 0.0389455 0.0682651 0.0771411 0.0780595 0.0514091 0.04727 0.0532247 0.0589278 0.0587372 0.0447354 0.0578246 0.0551842 0.0423223 0.0384075 0.0648066 0.058539 0.0723592 0.0727737 0.0680948 0.0745402 0.0721887 0.0729274 0.0701972 0.0864818 0.0815751 0.0828893 0.0604542 0.077119 0.0860442 0.078041 0.0571193 0.10313 0.142202 0.113455 0.154725 0.0903576 0.0674519 0.0625388 0.0886341 0.0683084 0.114282 0.0885254 0.0932744 0.0913002 0.0863415 0.090694 0.0865116 0.0919071 0.0740176 0.0950476 0.108272 0.0723946 0.0606791 0.0532695 0.0569644 0.0395408 0.0380145 0.0498398 0.0956759 0.0827328 0.0885967 0.0804372 0.0493466 0.0207758 0.00967936 0.00986253 0.0201281 0.0240086 0.020429 0.0130739 0.0216315 0.022531 0.0238832 0.0227054 0.0429088 0.0426656 0.042478 0.0422898 0.0420467 0.0418495 0.0416487 0.0413956 0.0411924 0.0409896
#
#_Q_setup for fleets with cpue or survey data
#_1:  fleet number
#_2:  link type: (1=simple q, 1 parm; 2=mirror simple q, 1 mirrored parm; 3=q and power, 2 parm; 4=mirror with offset, 2 parm)
#_3:  extra input for link, i.e. mirror fleet# or dev index number
#_4:  0/1 to select extra sd parameter
#_5:  0/1 for biasadj or not
#_6:  0/1 to float
#_   fleet      link link_info  extra_se   biasadj     float  #  fleetname
         1         3         0         0         1         0  #  WinterN
         3         3         0         0         1         0  #  WinterS
         5         1         0         1         0         1  #  TriEarly
         6         1         0         1         0         1  #  TriLate
         7         1         0         0         0         1  #  NWFSC
-9999 0 0 0 0 0
#
#_Q_parms(if_any);Qunits_are_ln(q)
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
           -20             5      -7.01915             0            99             0          1          0          0          0          0          0          5          1  #  LnQ_base_WinterN(1)
            -5             5     -0.120823             0            -1             0          3          0          0          0          0          0          0          0  #  Q_power_WinterN(1)
           -20             5       -1.3472             0            99             0          1          0          0          0          0          0          5          1  #  LnQ_base_WinterS(3)
            -5             5     -0.852654             0            -1             0          3          0          0          0          0          0          0          0  #  Q_power_WinterS(3)
           -15            15     -0.861191             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_TriEarly(5)
         0.001             2      0.218021          0.22            -1             0          5          0          0          0          0          0          0          0  #  Q_extraSD_TriEarly(5)
           -15            15     -0.430897             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_TriLate(6)
         0.001             2      0.312559          0.16            -1             0          4          0          0          0          0          0          0          0  #  Q_extraSD_TriLate(6)
           -15            15        1.0476             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_NWFSC(7)
# timevary Q parameters 
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type     PHASE  #  parm_name
         -0.99          0.99      0.490021             0           0.5             6      3  # LnQ_base_WinterN(1)_BLK5add_2004
         -0.99          0.99      0.619915             0           0.5             6      3  # LnQ_base_WinterS(3)_BLK5add_2004
# info on dev vectors created for Q parms are reported with other devs after tag parameter section 
#
#_size_selex_patterns
#Pattern:_0; parm=0; selex=1.0 for all sizes
#Pattern:_1; parm=2; logistic; with 95% width specification
#Pattern:_5; parm=2; mirror another size selex; PARMS pick the min-max bin to mirror
#Pattern:_15; parm=0; mirror another age or length selex
#Pattern:_6; parm=2+special; non-parm len selex
#Pattern:_43; parm=2+special+2;  like 6, with 2 additional param for scaling (average over bin range)
#Pattern:_8; parm=8; New doublelogistic with smooth transitions and constant above Linf option
#Pattern:_9; parm=6; simple 4-parm double logistic with starting length; parm 5 is first length; parm 6=1 does desc as offset
#Pattern:_21; parm=2+special; non-parm len selex, read as pairs of size, then selex
#Pattern:_22; parm=4; double_normal as in CASAL
#Pattern:_23; parm=6; double_normal where final value is directly equal to sp(6) so can be >1.0
#Pattern:_24; parm=6; double_normal with sel(minL) and sel(maxL), using joiners
#Pattern:_25; parm=3; exponential-logistic in size
#Pattern:_27; parm=3+special; cubic spline 
#Pattern:_42; parm=2+special+3; // like 27, with 2 additional param for scaling (average over bin range)
#_discard_options:_0=none;_1=define_retention;_2=retention&mortality;_3=all_discarded_dead;_4=define_dome-shaped_retention
#_Pattern Discard Male Special
 24 1 3 0 # 1 WinterN
 24 1 3 0 # 2 SummerN
 24 1 3 0 # 3 WinterS
 24 1 3 0 # 4 SummerS
 24 0 3 0 # 5 TriEarly
 24 0 3 0 # 6 TriLate
 24 0 3 0 # 7 NWFSC
#
#_age_selex_patterns
#Pattern:_0; parm=0; selex=1.0 for ages 0 to maxage
#Pattern:_10; parm=0; selex=1.0 for ages 1 to maxage
#Pattern:_11; parm=2; selex=1.0  for specified min-max age
#Pattern:_12; parm=2; age logistic
#Pattern:_13; parm=8; age double logistic
#Pattern:_14; parm=nages+1; age empirical
#Pattern:_15; parm=0; mirror another age or length selex
#Pattern:_16; parm=2; Coleraine - Gaussian
#Pattern:_17; parm=nages+1; empirical as random walk  N parameters to read can be overridden by setting special to non-zero
#Pattern:_41; parm=2+nages+1; // like 17, with 2 additional param for scaling (average over bin range)
#Pattern:_18; parm=8; double logistic - smooth transition
#Pattern:_19; parm=6; simple 4-parm double logistic with starting age
#Pattern:_20; parm=6; double_normal,using joiners
#Pattern:_26; parm=3; exponential-logistic in age
#Pattern:_27; parm=3+special; cubic spline in age
#Pattern:_42; parm=2+special+3; // cubic spline; with 2 additional param for scaling (average over bin range)
#_Pattern Discard Male Special
 10 0 0 0 # 1 WinterN
 10 0 0 0 # 2 SummerN
 10 0 0 0 # 3 WinterS
 10 0 0 0 # 4 SummerS
 10 0 0 0 # 5 TriEarly
 10 0 0 0 # 6 TriLate
 10 0 0 0 # 7 NWFSC
#
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
# 1   WinterN LenSelex
            15            75       48.6805          43.1             5             0          2          0          0          0          0        0.5          1          1  #  Size_DblN_peak_WinterN(1)
            -5             3             3           0.7             5             0         -3          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_WinterN(1)
            -4            12       4.30771          3.42             5             0          3          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_WinterN(1)
            -2            15            14          0.21             5             0         -3          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_WinterN(1)
           -15             5          -999          -8.9             5             0         -4          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_WinterN(1)
            -5             5          -999          0.15             5             0         -4          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_WinterN(1)
            10            40       28.0301            15             9             0          2          0          0          0          0          0          2          1  #  Retain_L_infl_WinterN(1)
           0.1            10        1.8503             3             9             0          4          0          0          0          0          0          2          1  #  Retain_L_width_WinterN(1)
           -10            10        8.3732            10             9             0          4          0          0          0          0          0          2          2  #  Retain_L_asymptote_logit_WinterN(1)
           -10            10             0             0             9             0         -2          0          0          0          0          0          0          0  #  Retain_L_maleoffset_WinterN(1)
           -15            15      -11.8861             0             5             0          4          0          0          0          0        0.5          0          0  #  SzSel_Male_Peak_WinterN(1)
           -15            15      -1.45306             0             5             0          4          0          0          0          0        0.5          0          0  #  SzSel_Male_Ascend_WinterN(1)
           -15            15             0             0             5             0         -4          0          0          0          0        0.5          0          0  #  SzSel_Male_Descend_WinterN(1)
           -15            15             0             0             5             0         -4          0          0          0          0        0.5          0          0  #  SzSel_Male_Final_WinterN(1)
           -15            15             1             0             5             0         -4          0          0          0          0        0.5          0          0  #  SzSel_Male_Scale_WinterN(1)
# 2   SummerN LenSelex
            15            75       48.4299          43.1             5             0          2          0          0          0          0        0.5          1          1  #  Size_DblN_peak_SummerN(2)
            -5             3             3           0.7             5             0         -3          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_SummerN(2)
            -4            12       5.29851          3.42             5             0          3          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_SummerN(2)
            -2            15            14          0.21             5             0         -3          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_SummerN(2)
           -15             5          -999          -8.9             5             0         -4          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_SummerN(2)
            -5             5          -999          0.15             5             0         -4          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_SummerN(2)
            10            40       30.6729            15             9             0          2          0          0          0          0          0          3          1  #  Retain_L_infl_SummerN(2)
           0.1            10       1.31436             3             9             0          4          0          0          0          0          0          3          1  #  Retain_L_width_SummerN(2)
           -10            10       9.37198            10             9             0          4          0          0          0          0          0          3          2  #  Retain_L_asymptote_logit_SummerN(2)
           -10            10             0             0             9             0         -2          0          0          0          0          0          0          0  #  Retain_L_maleoffset_SummerN(2)
           -20            15      -12.7368             0            -5             0          4          0          0          0          0        0.5          0          0  #  SzSel_Male_Peak_SummerN(2)
           -15            15      -1.89766             0            -5             0          4          0          0          0          0        0.5          0          0  #  SzSel_Male_Ascend_SummerN(2)
           -15            15             0             0             5             0         -4          0          0          0          0        0.5          0          0  #  SzSel_Male_Descend_SummerN(2)
           -15            15             0             0             5             0         -4          0          0          0          0        0.5          0          0  #  SzSel_Male_Final_SummerN(2)
           -15            15             1             0             5             0         -4          0          0          0          0        0.5          0          0  #  SzSel_Male_Scale_SummerN(2)
# 3   WinterS LenSelex
            15            75       38.4882          43.1             5             0          2          0          0          0          0        0.5          1          1  #  Size_DblN_peak_WinterS(3)
            -5             3             3           0.7             5             0         -3          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_WinterS(3)
            -4            12       4.41185          3.42             5             0          3          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_WinterS(3)
            -2            15            14          0.21             5             0         -3          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_WinterS(3)
           -15             5          -999          -8.9             5             0         -4          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_WinterS(3)
            -5             5          -999          0.15             5             0         -4          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_WinterS(3)
            10            40       28.8815            15             9             0          2          0          0          0          0          0          2          1  #  Retain_L_infl_WinterS(3)
           0.1            10       1.35726             3             9             0          3          0          0          0          0          0          2          1  #  Retain_L_width_WinterS(3)
           -10            10       3.97227            10             9             0          4          0          0          0          0          0          2          2  #  Retain_L_asymptote_logit_WinterS(3)
           -10            10             0             0             9             0         -2          0          0          0          0          0          0          0  #  Retain_L_maleoffset_WinterS(3)
           -15            15      -12.7221             0             5             0          4          0          0          0          0        0.5          0          0  #  SzSel_Male_Peak_WinterS(3)
           -15            15      -1.86133             0             5             0          4          0          0          0          0        0.5          0          0  #  SzSel_Male_Ascend_WinterS(3)
           -15            15             0             0             5             0         -4          0          0          0          0        0.5          0          0  #  SzSel_Male_Descend_WinterS(3)
           -15            15             0             0             5             0         -4          0          0          0          0        0.5          0          0  #  SzSel_Male_Final_WinterS(3)
           -15            15             1             0             5             0         -4          0          0          0          0        0.5          0          0  #  SzSel_Male_Scale_WinterS(3)
# 4   SummerS LenSelex
            15            75       40.6429          43.1             5             0          2          0          0          0          0        0.5          1          1  #  Size_DblN_peak_SummerS(4)
            -5             3             3           0.7             5             0         -3          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_SummerS(4)
            -4            12       4.89772          3.42             5             0          3          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_SummerS(4)
            -2            15            14          0.21             5             0         -3          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_SummerS(4)
           -15             5          -999          -8.9             5             0         -4          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_SummerS(4)
            -5             5          -999          0.15             5             0         -4          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_SummerS(4)
            10            40       28.8753            15             9             0          3          0          0          0          0          0          3          1  #  Retain_L_infl_SummerS(4)
           0.1            10       1.07128             3             9             0          3          0          0          0          0          0          3          1  #  Retain_L_width_SummerS(4)
           -10            10        9.5208            10             9             0          4          0          0          0          0          0          3          2  #  Retain_L_asymptote_logit_SummerS(4)
           -10            10             0             0             9             0         -2          0          0          0          0          0          0          0  #  Retain_L_maleoffset_SummerS(4)
           -15            15       -12.548             0             5             0          4          0          0          0          0        0.5          0          0  #  SzSel_Male_Peak_SummerS(4)
           -15            15      -1.89491             0             5             0          4          0          0          0          0        0.5          0          0  #  SzSel_Male_Ascend_SummerS(4)
           -15            15             0             0             5             0         -4          0          0          0          0        0.5          0          0  #  SzSel_Male_Descend_SummerS(4)
           -15            15             0             0             5             0         -4          0          0          0          0        0.5          0          0  #  SzSel_Male_Final_SummerS(4)
           -15            15             1             0             5             0         -4          0          0          0          0        0.5          0          0  #  SzSel_Male_Scale_SummerS(4)
# 5   TriEarly LenSelex
            15            61       35.3503          43.1             5             0          2          0          0          0          0        0.5          0          0  #  Size_DblN_peak_TriEarly(5)
            -5             3             3           0.7             5             0         -2          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_TriEarly(5)
            -4            12       4.21179          3.42             5             0          2          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_TriEarly(5)
            -2            15            14          0.21             5             0         -2          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_TriEarly(5)
           -15             5          -999          -8.9             5             0         -4          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_TriEarly(5)
            -5             5          -999          0.15             5             0         -4          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_TriEarly(5)
           -15            15      -3.88585             0             5             0          3          0          0          0          0        0.5          0          0  #  SzSel_Male_Peak_TriEarly(5)
           -15            15     -0.561008             0             5             0          3          0          0          0          0        0.5          0          0  #  SzSel_Male_Ascend_TriEarly(5)
           -15            15             0             0             5             0         -3          0          0          0          0        0.5          0          0  #  SzSel_Male_Descend_TriEarly(5)
           -15            15             0             0             5             0         -3          0          0          0          0        0.5          0          0  #  SzSel_Male_Final_TriEarly(5)
           -15            15             1             0             5             0         -4          0          0          0          0        0.5          0          0  #  SzSel_Male_Scale_TriEarly(5)
# 6   TriLate LenSelex
            15            61       36.5056          43.1             5             0          2          0          0          0          0        0.5          0          0  #  Size_DblN_peak_TriLate(6)
            -5             3             3           0.7             5             0         -2          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_TriLate(6)
            -4            12       4.64265          3.42             5             0          2          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_TriLate(6)
            -2            15            14          0.21             5             0         -2          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_TriLate(6)
           -15             5          -999          -8.9             5             0         -4          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_TriLate(6)
            -5             5          -999          0.15             5             0         -4          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_TriLate(6)
           -15            15      -2.23813             0             5             0          3          0          0          0          0        0.5          0          0  #  SzSel_Male_Peak_TriLate(6)
           -15            15    -0.0352576             0             5             0          3          0          0          0          0        0.5          0          0  #  SzSel_Male_Ascend_TriLate(6)
           -15            15             0             0             5             0         -3          0          0          0          0        0.5          0          0  #  SzSel_Male_Descend_TriLate(6)
           -15            15             0             0             5             0         -3          0          0          0          0        0.5          0          0  #  SzSel_Male_Final_TriLate(6)
           -15            15             1             0             5             0         -4          0          0          0          0        0.5          0          0  #  SzSel_Male_Scale_TriLate(6)
# 7   NWFSC LenSelex
            15            61       43.0085          43.1             5             0          2          0          0          0          0        0.5          0          0  #  Size_DblN_peak_NWFSC(7)
            -5             3             3           0.7             5             0         -2          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_NWFSC(7)
            -4            12       5.14971          3.42             5             0          2          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_NWFSC(7)
            -2            15            14          0.21             5             0         -2          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_NWFSC(7)
           -15             5          -999          -8.9             5             0         -4          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_NWFSC(7)
            -5             5          -999          0.15             5             0         -4          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_NWFSC(7)
           -15            15       -5.0654             0             5             0          3          0          0          0          0        0.5          0          0  #  SzSel_Male_Peak_NWFSC(7)
           -15            15     -0.410501             0             5             0          3          0          0          0          0        0.5          0          0  #  SzSel_Male_Ascend_NWFSC(7)
           -15            15             0             0             5             0         -3          0          0          0          0        0.5          0          0  #  SzSel_Male_Descend_NWFSC(7)
           -15            15             0             0             5             0         -3          0          0          0          0        0.5          0          0  #  SzSel_Male_Final_NWFSC(7)
           -15            15             1             0             5             0         -4          0          0          0          0        0.5          0          0  #  SzSel_Male_Scale_NWFSC(7)
# 1   WinterN AgeSelex
# 2   SummerN AgeSelex
# 3   WinterS AgeSelex
# 4   SummerS AgeSelex
# 5   TriEarly AgeSelex
# 6   TriLate AgeSelex
# 7   NWFSC AgeSelex
# timevary selex parameters 
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type    PHASE  #  parm_name
         -31.6          28.4       2.05549             0          14.2             6      5  # Size_DblN_peak_WinterN(1)_BLK1add_1973
         -31.6          28.4      -1.98219             0          14.2             6      5  # Size_DblN_peak_WinterN(1)_BLK1add_1983
         -31.6          28.4     -0.790166             0          14.2             6      5  # Size_DblN_peak_WinterN(1)_BLK1add_1993
         -31.6          28.4      0.383774             0          14.2             6      5  # Size_DblN_peak_WinterN(1)_BLK1add_2003
         -31.6          28.4      0.878534             0          14.2             6      5  # Size_DblN_peak_WinterN(1)_BLK1add_2011
        -16.19         13.81      -2.26483             0         6.905             6      5  # Retain_L_infl_WinterN(1)_BLK2add_2003
        -16.19         13.81       1.76566             0         6.905             6      5  # Retain_L_infl_WinterN(1)_BLK2add_2010
        -16.19         13.81      -3.23884             0         6.905             6      5  # Retain_L_infl_WinterN(1)_BLK2add_2011
        -1.601         8.299       0.12285             0        0.8005             6      5  # Retain_L_width_WinterN(1)_BLK2add_2003
        -1.601         8.299      0.393503             0        0.8005             6      5  # Retain_L_width_WinterN(1)_BLK2add_2010
        -1.601         8.299     -0.697801             0        0.8005             6      5  # Retain_L_width_WinterN(1)_BLK2add_2011
           -10            10       6.63678            10             9             0      5  # Retain_L_asymptote_logit_WinterN(1)_BLK2repl_2003
           -10            10       2.11235            10             9             0      5  # Retain_L_asymptote_logit_WinterN(1)_BLK2repl_2010
           -10            10        9.9881            10             9             0      5  # Retain_L_asymptote_logit_WinterN(1)_BLK2repl_2011
         -38.8          21.2        1.9558             0          10.6             6      5  # Size_DblN_peak_SummerN(2)_BLK1add_1973
         -38.8          21.2     -0.399466             0          10.6             6      5  # Size_DblN_peak_SummerN(2)_BLK1add_1983
         -38.8          21.2      -2.38392             0          10.6             6      5  # Size_DblN_peak_SummerN(2)_BLK1add_1993
         -38.8          21.2    -0.0936437             0          10.6             6      5  # Size_DblN_peak_SummerN(2)_BLK1add_2003
         -38.8          21.2       3.26642             0          10.6             6      5  # Size_DblN_peak_SummerN(2)_BLK1add_2011
       -20.679         9.321     -0.419922             0        4.6605             6      5  # Retain_L_infl_SummerN(2)_BLK3add_2003
       -20.679         9.321       1.29282             0        4.6605             6      5  # Retain_L_infl_SummerN(2)_BLK3add_2009
       -20.679         9.321      -1.92562             0        4.6605             6      5  # Retain_L_infl_SummerN(2)_BLK3add_2011
       -1.0278        8.8722      0.161789             0        0.5139             6      5  # Retain_L_width_SummerN(2)_BLK3add_2003
       -1.0278        8.8722      0.138985             0        0.5139             6      5  # Retain_L_width_SummerN(2)_BLK3add_2009
       -1.0278        8.8722       0.21153             0        0.5139             6      5  # Retain_L_width_SummerN(2)_BLK3add_2011
           -10            10       5.45846            10             9             0      5  # Retain_L_asymptote_logit_SummerN(2)_BLK3repl_2003
           -10            10       7.55883            10             9             0      5  # Retain_L_asymptote_logit_SummerN(2)_BLK3repl_2009
           -10            10       6.15877            10             9             0      5  # Retain_L_asymptote_logit_SummerN(2)_BLK3repl_2011
       -25.422        34.578      -15.0996             0        12.711             6      5  # Size_DblN_peak_WinterS(3)_BLK1add_1973
       -25.422        34.578       5.08526             0        12.711             6      5  # Size_DblN_peak_WinterS(3)_BLK1add_1983
       -25.422        34.578       9.18368             0        12.711             6      5  # Size_DblN_peak_WinterS(3)_BLK1add_1993
       -25.422        34.578       6.98535             0        12.711             6      5  # Size_DblN_peak_WinterS(3)_BLK1add_2003
       -25.422        34.578       8.30205             0        12.711             6      5  # Size_DblN_peak_WinterS(3)_BLK1add_2011
       -18.816        11.184       -2.0172             0         5.592             6      5  # Retain_L_infl_WinterS(3)_BLK2add_2003
       -18.816        11.184       1.52688             0         5.592             6      5  # Retain_L_infl_WinterS(3)_BLK2add_2010
       -18.816        11.184      -4.29967             0         5.592             6      5  # Retain_L_infl_WinterS(3)_BLK2add_2011
       -1.0443        8.8557      0.366784             0       0.52215             6      5  # Retain_L_width_WinterS(3)_BLK2add_2003
       -1.0443        8.8557       0.13891             0       0.52215             6      5  # Retain_L_width_WinterS(3)_BLK2add_2010
       -1.0443        8.8557    -0.0497998             0       0.52215             6      5  # Retain_L_width_WinterS(3)_BLK2add_2011
           -10            10       7.80958            10             9             0      5  # Retain_L_asymptote_logit_WinterS(3)_BLK2repl_2003
           -10            10       5.53053            10             9             0      5  # Retain_L_asymptote_logit_WinterS(3)_BLK2repl_2010
           -10            10       7.87413            10             9             0      5  # Retain_L_asymptote_logit_WinterS(3)_BLK2repl_2011
      -28.0793       31.9207      -5.16048             0       14.0397             6      5  # Size_DblN_peak_SummerS(4)_BLK1add_1973
      -28.0793       31.9207      -6.38506             0       14.0397             6      5  # Size_DblN_peak_SummerS(4)_BLK1add_1983
      -28.0793       31.9207       3.57295             0       14.0397             6      5  # Size_DblN_peak_SummerS(4)_BLK1add_1993
      -28.0793       31.9207        6.2503             0       14.0397             6      5  # Size_DblN_peak_SummerS(4)_BLK1add_2003
      -28.0793       31.9207       6.03167             0       14.0397             6      5  # Size_DblN_peak_SummerS(4)_BLK1add_2011
       -19.055        10.945      -1.40609             0        5.4725             6      5  # Retain_L_infl_SummerS(4)_BLK3add_2003
       -19.055        10.945      -1.68655             0        5.4725             6      5  # Retain_L_infl_SummerS(4)_BLK3add_2009
       -19.055        10.945       -2.0893             0        5.4725             6      5  # Retain_L_infl_SummerS(4)_BLK3add_2011
        -0.876         9.024      0.604487             0         0.438             6      5  # Retain_L_width_SummerS(4)_BLK3add_2003
        -0.876         9.024       0.47059             0         0.438             6      5  # Retain_L_width_SummerS(4)_BLK3add_2009
        -0.876         9.024      0.581508             0         0.438             6      5  # Retain_L_width_SummerS(4)_BLK3add_2011
           -10            10       7.55833            10             9             0      5  # Retain_L_asymptote_logit_SummerS(4)_BLK3repl_2003
           -10            10       8.89461            10             9             0      5  # Retain_L_asymptote_logit_SummerS(4)_BLK3repl_2009
           -10            10       7.68067            10             9             0      5  # Retain_L_asymptote_logit_SummerS(4)_BLK3repl_2011
# info on dev vectors created for selex parms are reported with other devs after tag parameter section 
#
0   #  use 2D_AR1 selectivity(0/1):  experimental feature
#_no 2D_AR1 selex offset used
#
# Tag loss and Tag reporting parameters go next
0  # TG_custom:  0=no read; 1=read if tags exist
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
#
# deviation vectors for timevary parameters
#  base   base first block   block  env  env   dev   dev   dev   dev   dev
#  type  index  parm trend pattern link  var  vectr link _mnyr  mxyr phase  dev_vector
#      3     1     1     5     1     0     0     0     0     0     0     0
#      3     3     2     5     1     0     0     0     0     0     0     0
#      5     1     3     1     1     0     0     0     0     0     0     0
#      5     7     8     2     1     0     0     0     0     0     0     0
#      5     8    11     2     1     0     0     0     0     0     0     0
#      5     9    14     2     2     0     0     0     0     0     0     0
#      5    16    17     1     1     0     0     0     0     0     0     0
#      5    22    22     3     1     0     0     0     0     0     0     0
#      5    23    25     3     1     0     0     0     0     0     0     0
#      5    24    28     3     2     0     0     0     0     0     0     0
#      5    31    31     1     1     0     0     0     0     0     0     0
#      5    37    36     2     1     0     0     0     0     0     0     0
#      5    38    39     2     1     0     0     0     0     0     0     0
#      5    39    42     2     2     0     0     0     0     0     0     0
#      5    46    45     1     1     0     0     0     0     0     0     0
#      5    52    50     3     1     0     0     0     0     0     0     0
#      5    53    53     3     1     0     0     0     0     0     0     0
#      5    54    56     3     2     0     0     0     0     0     0     0
     #
# Input variance adjustments factors: 
 #_1=add_to_survey_CV
 #_2=add_to_discard_stddev
 #_3=add_to_bodywt_CV
 #_4=mult_by_lencomp_N
 #_5=mult_by_agecomp_N
 #_6=mult_by_size-at-age_N
 #_7=mult_by_generalized_sizecomp
#_Factor  Fleet  Value
      2      1      0.02
      2      2      0.02
      2      3      0.02
      2      4      0.02
      4      1     1.366
      4      2     1.039
      4      3     1.017
      4      4     1.169
      4      5     1.807
      4      6     1.285
      4      7     0.579
      5      1     2.926
      5      2      2.45
      5      3     1.756
      5      4     1.601
      5      7     0.215
 -9999   1    0  # terminator
#
15 #_maxlambdaphase
1 #_sd_offset; must be 1 if any growthCV, sigmaR, or survey extraSD is an estimated parameter
# read 10 changes to default Lambdas (default value is 1.0)
# Like_comp codes:  1=surv; 2=disc; 3=mnwt; 4=length; 5=age; 6=SizeFreq; 7=sizeage; 8=catch; 9=init_equ_catch; 
# 10=recrdev; 11=parm_prior; 12=parm_dev; 13=CrashPen; 14=Morphcomp; 15=Tag-comp; 16=Tag-negbin; 17=F_ballpark; 18=initEQregime
#like_comp fleet  phase  value  sizefreq_method
 1 1 1 1 1
 1 3 1 1 1
 5 1 1 0.5 1
 5 2 1 0.5 1
 5 3 1 0.5 1
 5 4 1 0.5 1
 4 1 1 0.5 1
 4 2 1 0.5 1
 4 3 1 0.5 1
 4 4 1 0.5 1
-9999  1  1  1  1  #  terminator
#
# lambdas (for info only; columns are phases)
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_CPUE/survey:_1
#  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 #_CPUE/survey:_2
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_CPUE/survey:_3
#  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 #_CPUE/survey:_4
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_CPUE/survey:_5
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_CPUE/survey:_6
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_CPUE/survey:_7
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_discard:_1
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_discard:_2
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_discard:_3
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_discard:_4
#  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 #_discard:_5
#  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 #_discard:_6
#  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 #_discard:_7
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_meanbodywt:1
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_meanbodywt:2
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_meanbodywt:3
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_meanbodywt:4
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_meanbodywt:5
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_meanbodywt:6
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_meanbodywt:7
#  0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 #_lencomp:_1
#  0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 #_lencomp:_2
#  0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 #_lencomp:_3
#  0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 #_lencomp:_4
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_lencomp:_5
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_lencomp:_6
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_lencomp:_7
#  0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 #_agecomp:_1
#  0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 #_agecomp:_2
#  0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 #_agecomp:_3
#  0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 #_agecomp:_4
#  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 #_agecomp:_5
#  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 #_agecomp:_6
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_agecomp:_7
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_init_equ_catch
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_recruitments
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_parameter-priors
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_parameter-dev-vectors
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_crashPenLambda
#  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 # F_ballpark_lambda
0 # (0/1) read specs for more stddev reporting 
 # 0 0 0 0 0 0 0 0 0 # placeholder for # selex_fleet, 1=len/2=age/3=both, year, N selex bins, 0 or Growth pattern, N growth ages, 0 or NatAge_area(-1 for all), NatAge_yr, N Natages
 # placeholder for vector of selex bins to be reported
 # placeholder for vector of growth ages to be reported
 # placeholder for vector of NatAges ages to be reported
999

