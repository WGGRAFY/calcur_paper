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
#LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env_var&link	dev_link	dev_minyr	dev_maxyr	dev_PH	Block	Block_Fxn		
0.005	0.5	0.156244	-1.7793	0.438	3	2	0	0	0	0	0.5	0	0	#	NatM_p_1_Fem_GP_1
10	45	15.6596	17.18	10	0	3	0	0	0	0	0.5	0	0	#	L_at_Amin_Fem_GP_1
35	80	53.1137	54.2	10	0	3	0	0	0	0	0.5	0	0	#	L_at_Amax_Fem_GP_1
0.04	0.5	0.141811	0.157	0.8	0	3	0	0	0	0	0.5	0	0	#	VonBert_K_Fem_GP_1
0.01	1	0.186061	3	0.8	0	3	0	0	0	0	0.5	0	0	#	CV_young_Fem_GP_1
0.01	1	0.0351848	0	1	0	4	0	0	0	0	0	0	0	#	CV_old_Fem_GP_1
-3	3	1.99E-06	1.99E-06	0.8	6	-3	0	0	0	0	0.5	0	0	#	Wtlen_1_Fem_GP_1
1	5	3.484	3.478	0.8	6	-3	0	0	0	0	0.5	0	0	#	Wtlen_2_Fem_GP_1
10	50	33.1	33.1	0.8	6	-3	0	0	0	0	0.5	0	0	#	Mat50%_Fem_GP_1
-3	3	-0.743	-0.743	0.8	6	-3	0	0	0	0	0.5	0	0	#	Mat_slope_Fem_GP_1
-3	3	1	1	1	6	-3	0	0	0	0	0.5	0	0	#	Eggs/kg_inter_Fem_GP_1
-3	3	0	0	1	6	-3	0	0	0	0	0.5	0	0	#	Eggs/kg_slope_wt_Fem_GP_1
0.005	0.6	0.161585	-1.6809	0.438	3	2	0	0	0	0	0.5	0	0	#	NatM_p_1_Mal_GP_1
10	45	16.1562	17.18	10	0	3	0	0	0	0	0.5	0	0	#	L_at_Amin_Mal_GP_1
35	80	40.8118	41.1	10	0	3	0	0	0	0	0.5	0	0	#	L_at_Amax_Mal_GP_1
0.04	0.5	0.23867	0.247	0.8	0	3	0	0	0	0	0.5	0	0	#	VonBert_K_Mal_GP_1
0.01	1	0.136254	3	0.8	0	3	0	0	0	0	0.5	0	0	#	CV_young_Mal_GP_1
0.01	1	0.059832	0	1	0	4	0	0	0	0	0	0	0	#	CV_old_Mal_GP_1
-3	3	2.98E-06	2.98E-06	0.8	6	-3	0	0	0	0	0.5	0	0	#	Wtlen_1_Mal_GP_1
-3	5	3.363	3.363	0.8	6	-3	0	0	0	0	0.5	0	0	#	Wtlen_2_Mal_GP_1
0	1	1	1	0	0	-4	0	0	0	0	0	0	0	#	CohortGrowDev
0.01	0.99	0.5	0.5	0.5	0	-99	0	0	0	0	0	0	0	#	FracFemale_GP_1
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
             5            20       9.88647             9            10             0          1          0          0          0          0          0          0          0 # SR_LN(R0)
           0.2             1      0.844442           0.8          0.09             6          5          0          0          0          0          0          0          0 # SR_BH_steep
             0             2           0.4           0.9             5             6        -99          0          0          0          0          0          0          0 # SR_sigmaR
            -5             5             0             0           0.2             6         -2          0          0          0          0          0          0          0 # SR_regime
             0             0             0             0             0             0        -99          0          0          0          0          0          0          0 # SR_autocorr
1 #do_recdev:  0=none; 1=devvector (R=F(SSB)+dev); 2=deviations (R=F(SSB)+dev); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
1959 # first year of main recr_devs; early devs can preceed this era
2016 # last year of main recr_devs; forecast devs start in following year
1 #_recdev phase 
1 # (0/1) to read 13 advanced options
 1845 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
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
#  2.29437e-07 2.68076e-07 3.13194e-07 3.65836e-07 4.27217e-07 4.9881e-07 5.82199e-07 6.79313e-07 7.92313e-07 9.23762e-07 1.07647e-06 1.25371e-06 1.4592e-06 1.69721e-06 1.97244e-06 2.29027e-06 2.65672e-06 3.0783e-06 3.56241e-06 4.11708e-06 4.75108e-06 5.47389e-06 6.29558e-06 7.22631e-06 8.27567e-06 9.45285e-06 1.0772e-05 1.22583e-05 1.39435e-05 1.5857e-05 1.80294e-05 2.04951e-05 2.32936e-05 2.64697e-05 3.0074e-05 3.41635e-05 3.88022e-05 4.40629e-05 5.00276e-05 5.6789e-05 6.4452e-05 7.31351e-05 8.29723e-05 9.41144e-05 0.000106732 0.000121018 0.000137191 0.000155495 0.000176207 0.000199642 0.000226145 0.000256117 0.000290004 0.000328307 0.000371591 0.000420491 0.000475734 0.000538139 0.000608617 0.000688204 0.000778073 0.000879529 0.000994029 0.00112326 0.00126817 0.0014311 0.00161483 0.00182146 0.00205364 0.00231456 0.00260619 0.00292828 0.00328178 0.00366953 0.00409787 0.00456989 0.00508703 0.00564967 0.00625572 0.00690782 0.00760748 0.00835696 0.00915926 0.010014 0.0109359 0.0119629 0.0132273 0.0149784 0.0176117 0.0216917 0.027686 0.0356991 0.0449477 0.0526263 0.0526913 0.0368429 0.000643012 -0.0496157 -0.0944988 -0.106979 -0.076277 -0.0368617 -0.0300844 -0.0361086 -0.0387217 -0.035171 -0.0246007 -0.0112969 -0.0139742 -0.0390071 -0.077414 -0.12275 -0.182509 -0.216091 -0.187412 0.00484388 0.20947 -0.226096 -0.294029 0.179089 -0.0765225 0.706844 -0.122426 -0.0565472 0.0562605 0.177536 0.0534492 -0.187231 -0.296688 -0.0829469 -0.0514846 0.240056 0.319464 0.119279 -0.124719 -0.0452619 -0.171391 -0.173281 -0.04517 0.309539 -0.173765 -0.596159 -0.43574 0.0190181 0.332212 0.301322 -0.15017 -0.566764 0.0667041 0.24912 -0.334158 -0.217835 -0.213209 0.663645 0.225236 -0.136543 -0.213387 -0.135606 -0.37868 -0.158423 -0.0818691 0.520858 0.725727 1.01221 0.128247 -0.129426 0.00348581 0.345438 -0.231406 -0.25369 -0.322962 -0.098052 -0.120943 0 0 0 0 0 0 0 0 0 0 0 0 0
# implementation error by year in forecast:  0 0 0 0 0 0 0 0 0 0 0 0
#
#Fishing Mortality info 
0.3 # F ballpark
-2001 # F ballpark year (neg value to disable)
3 # F_Method:  1=Pope; 2=instan. F; 3=hybrid (hybrid is recommended)
4 # max F or harvest rate, depends on F_Method
# no additional F input needed for Fmethod 1
# if Fmethod=2; read overall start F value; overall phase; N detailed inputs to read
# if Fmethod=3; read N iterations for tuning for Fmethod 3
5  # N iterations for tuning F in hybrid method (recommend 3 to 7)
#
#_initial_F_parms; count = 0
#_ LO HI INIT PRIOR PR_SD  PR_type  PHASE
#2030 2039
# F rates by fleet
# Yr:  1876 1877 1878 1879 1880 1881 1882 1883 1884 1885 1886 1887 1888 1889 1890 1891 1892 1893 1894 1895 1896 1897 1898 1899 1900 1901 1902 1903 1904 1905 1906 1907 1908 1909 1910 1911 1912 1913 1914 1915 1916 1917 1918 1919 1920 1921 1922 1923 1924 1925 1926 1927 1928 1929 1930 1931 1932 1933 1934 1935 1936 1937 1938 1939 1940 1941 1942 1943 1944 1945 1946 1947 1948 1949 1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960 1961 1962 1963 1964 1965 1966 1967 1968 1969 1970 1971 1972 1973 1974 1975 1976 1977 1978 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 2023 2024 2025 2026 2027 2028 2029 2030
# seas:  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
# WinterN 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5.46664e-05 0.000165034 0.000279507 0.000401036 0.000467665 0.000595063 0.000855212 0.00114102 0.00132426 0.00145696 0.00169418 0.00197685 0.00225971 0.00252088 0.00284304 0.00323092 0.00372443 0.00429659 0.0100661 0.00846451 0.00994839 0.00359361 0.00211275 0.00464813 0.0112034 0.0141024 0.00840812 0.0286147 0.021132 0.0200521 0.0295771 0.0730589 0.0337993 0.0466339 0.0474541 0.054229 0.0383025 0.0256401 0.0563045 0.104018 0.105219 0.095207 0.104674 0.122399 0.127297 0.184088 0.227766 0.282851 0.185341 0.435103 0.68952 0.251245 0.147746 0.132434 0.160008 0.31582 0.306499 0.306522 0.29759 0.412023 0.376542 0.487398 0.317504 0.288131 0.248816 0.253546 0.209429 0.173627 0.218466 0.246602 0.230985 0.19656 0.310254 0.277621 0.15174 0.261406 0.245122 0.256684 0.0838213 0.0497287 0.0669269 0.063249 0.0845883 0.0891694 0.0677853 0.0858546 0.0716311 0.0787989 0.0810431 0.113137 0.112471 0.111955 0.111447 0.110799 0.110277 0.109746 0.109074 0.108533 0.107992
# SummerN 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5.49243e-06 4.50864e-06 3.52085e-06 3.43238e-06 3.34474e-06 3.25678e-06 3.16851e-06 3.07875e-06 2.98977e-06 2.90043e-06 2.80956e-06 2.71944e-06 2.62892e-06 2.53685e-06 2.4455e-06 2.35372e-06 2.26152e-06 2.16772e-06 2.07463e-06 1.98109e-06 1.88583e-06 1.79359e-06 1.701e-06 1.60395e-06 1.50194e-06 1.40036e-06 1.30157e-06 1.20315e-06 1.10669e-06 1.01066e-06 9.12392e-07 8.55471e-07 1.21184e-08 3.7996e-05 3.03365e-05 0.00202759 0.00628304 0.0103045 0.0145521 0.0170541 0.0205908 0.0287288 0.0334192 0.045178 0.0529524 0.0569032 0.0960394 0.099536 0.0741141 0.0693459 0.0967684 0.0731059 0.10271 0.0894792 0.129962 0.102706 0.0898012 0.0485165 0.060366 0.0587747 0.0518287 0.0760088 0.081423 0.0679336 0.100405 0.119234 0.126229 0.0948574 0.107302 0.0969607 0.100704 0.0838777 0.0870564 0.0951046 0.115935 0.0949498 0.111849 0.164755 0.221743 0.239883 0.161666 0.13591 0.226645 0.275835 0.297366 0.157956 0.281038 0.283406 0.22227 0.160269 0.2089 0.18303 0.214787 0.215783 0.179233 0.174511 0.181428 0.165919 0.12508 0.13124 0.116715 0.115681 0.154269 0.126603 0.141883 0.144209 0.162822 0.16082 0.145695 0.201345 0.193655 0.105024 0.072365 0.14263 0.0578138 0.0789822 0.0672647 0.109225 0.0742853 0.0825415 0.0836756 0.0884776 0.0881981 0.0930632 0.088026 0.132383 0.131617 0.131025 0.130438 0.129686 0.129078 0.128459 0.127676 0.127046 0.126417
# WinterS 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.00144861 0.000836267 0.000892065 0.00327235 0.0037295 0.00233286 0.00185636 0.0012256 0.000817938 0.00445709 0.00314693 0.000717786 0.00221465 0.00274434 0.00334331 0.00246569 0.00557497 0.0106568 0.0267414 0.0205913 0.0109983 0.0178328 0.0295071 0.0448455 0.0299625 0.0289957 0.0303129 0.0334521 0.022669 0.0316891 0.034916 0.0358779 0.0404374 0.0303309 0.0305235 0.0267853 0.0449208 0.030866 0.0355837 0.0389263 0.0442731 0.0594545 0.0345663 0.0584069 0.054548 0.0769491 0.0568194 0.0546706 0.0535306 0.0866413 0.0622667 0.0472446 0.0718613 0.068352 0.069591 0.103885 0.0872112 0.106196 0.125725 0.109037 0.17251 0.115681 0.154482 0.136588 0.104385 0.157243 0.169977 0.0826105 0.100298 0.124626 0.102979 0.108237 0.0536355 0.0324098 0.0611075 0.022697 0.074951 0.100441 0.0928902 0.0144752 0.0057039 0.0131074 0.0106116 0.0188848 0.0135693 0.0142297 0.0118493 0.0130324 0.00867561 0.0156264 0.0268275 0.0266733 0.0265543 0.0264355 0.026283 0.0261595 0.0260339 0.0258754 0.0257481 0.025621
# SummerS 2.02374e-05 2.02378e-05 2.02381e-05 2.02385e-05 0.000233781 0.000447462 0.000661403 0.000875714 0.00109049 0.00130583 0.0015218 0.00173846 0.00195588 0.0021741 0.00239315 0.00261308 0.00283391 0.00305588 0.00327859 0.00350227 0.00372694 0.00395262 0.0041793 0.004407 0.00463572 0.00486549 0.00509629 0.00532814 0.00556105 0.00579502 0.00603005 0.00626615 0.00650332 0.00674156 0.00698087 0.00722126 0.00746272 0.00770524 0.00794883 0.00819348 0.00833224 0.0113899 0.00920394 0.00725012 0.00500687 0.00637153 0.00921731 0.00928988 0.0116156 0.0115626 0.0114518 0.0139319 0.0137346 0.0157203 0.0147461 0.0119349 0.0117528 0.00892271 0.0207012 0.0183545 0.0103726 0.0181638 0.0224284 0.0268588 0.0161781 0.00925806 0.00625109 0.0104901 0.014121 0.0142158 0.0376696 0.0390385 0.068449 0.0773834 0.0783561 0.0516444 0.0475208 0.0535408 0.059311 0.0591478 0.0450612 0.0582628 0.0556232 0.0426684 0.0387276 0.06537 0.059078 0.0730649 0.0735088 0.0687847 0.0752924 0.0729139 0.0736448 0.0708634 0.0872826 0.0823377 0.083686 0.0610772 0.0779352 0.086993 0.0788989 0.0576885 0.104043 0.143362 0.114335 0.155832 0.0909499 0.06778 0.0627702 0.0888718 0.0684639 0.114556 0.0887654 0.0935629 0.091603 0.0866467 0.0910688 0.0866608 0.0921151 0.0742165 0.0953453 0.108656 0.0726709 0.060916 0.0534813 0.0572003 0.0397108 0.0381928 0.0500778 0.0961877 0.0832395 0.0892135 0.0810939 0.0498321 0.0209955 0.0097776 0.00995706 0.0203139 0.0242234 0.0206045 0.0131797 0.0217953 0.0226898 0.0240377 0.0228423 0.0418206 0.041583 0.0413995 0.0412154 0.0409783 0.0407862 0.0405908 0.0403447 0.040147 0.0399497
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
           -20             5      -7.06261             0            99             0          1          0          0          0          0          0          5          1  #  LnQ_base_WinterN(1)
            -5             5     -0.114972             0            -1             0          3          0          0          0          0          0          0          0  #  Q_power_WinterN(1)
           -20             5      -1.36181             0            99             0          1          0          0          0          0          0          5          1  #  LnQ_base_WinterS(3)
            -5             5     -0.850878             0            -1             0          3          0          0          0          0          0          0          0  #  Q_power_WinterS(3)
           -15            15     -0.857253             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_TriEarly(5)
         0.001             2      0.216409          0.22            -1             0          5          0          0          0          0          0          0          0  #  Q_extraSD_TriEarly(5)
           -15            15     -0.425227             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_TriLate(6)
         0.001             2      0.313165          0.16            -1             0          4          0          0          0          0          0          0          0  #  Q_extraSD_TriLate(6)
           -15            15       1.05559             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_NWFSC(7)
# timevary Q parameters 
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type     PHASE  #  parm_name
         -0.99          0.99      0.491681             0           0.5             6      3  # LnQ_base_WinterN(1)_BLK5add_2004
         -0.99          0.99      0.619542             0           0.5             6      3  # LnQ_base_WinterS(3)_BLK5add_2004
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
            15            75       48.6767          43.1             5             0          2          0          0          0          0        0.5          1          1  #  Size_DblN_peak_WinterN(1)
            -5             3             3           0.7             5             0         -3          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_WinterN(1)
            -4            12       4.30851          3.42             5             0          3          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_WinterN(1)
            -2            15            14          0.21             5             0         -3          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_WinterN(1)
           -15             5          -999          -8.9             5             0         -4          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_WinterN(1)
            -5             5          -999          0.15             5             0         -4          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_WinterN(1)
            10            40       28.0498            15             9             0          2          0          0          0          0          0          2          1  #  Retain_L_infl_WinterN(1)
           0.1            10       1.84942             3             9             0          4          0          0          0          0          0          2          1  #  Retain_L_width_WinterN(1)
           -10            10       8.45982            10             9             0          4          0          0          0          0          0          2          2  #  Retain_L_asymptote_logit_WinterN(1)
           -10            10             0             0             9             0         -2          0          0          0          0          0          0          0  #  Retain_L_maleoffset_WinterN(1)
           -15            15      -11.8875             0             5             0          4          0          0          0          0        0.5          0          0  #  SzSel_Male_Peak_WinterN(1)
           -15            15      -1.45432             0             5             0          4          0          0          0          0        0.5          0          0  #  SzSel_Male_Ascend_WinterN(1)
           -15            15             0             0             5             0         -4          0          0          0          0        0.5          0          0  #  SzSel_Male_Descend_WinterN(1)
           -15            15             0             0             5             0         -4          0          0          0          0        0.5          0          0  #  SzSel_Male_Final_WinterN(1)
           -15            15             1             0             5             0         -4          0          0          0          0        0.5          0          0  #  SzSel_Male_Scale_WinterN(1)
# 2   SummerN LenSelex
            15            75       48.4248          43.1             5             0          2          0          0          0          0        0.5          1          1  #  Size_DblN_peak_SummerN(2)
            -5             3             3           0.7             5             0         -3          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_SummerN(2)
            -4            12       5.30037          3.42             5             0          3          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_SummerN(2)
            -2            15            14          0.21             5             0         -3          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_SummerN(2)
           -15             5          -999          -8.9             5             0         -4          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_SummerN(2)
            -5             5          -999          0.15             5             0         -4          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_SummerN(2)
            10            40       30.6749            15             9             0          2          0          0          0          0          0          3          1  #  Retain_L_infl_SummerN(2)
           0.1            10       1.31459             3             9             0          4          0          0          0          0          0          3          1  #  Retain_L_width_SummerN(2)
           -10            10       9.36966            10             9             0          4          0          0          0          0          0          3          2  #  Retain_L_asymptote_logit_SummerN(2)
           -10            10             0             0             9             0         -2          0          0          0          0          0          0          0  #  Retain_L_maleoffset_SummerN(2)
           -20            15      -12.7281             0            -5             0          4          0          0          0          0        0.5          0          0  #  SzSel_Male_Peak_SummerN(2)
           -15            15      -1.90014             0            -5             0          4          0          0          0          0        0.5          0          0  #  SzSel_Male_Ascend_SummerN(2)
           -15            15             0             0             5             0         -4          0          0          0          0        0.5          0          0  #  SzSel_Male_Descend_SummerN(2)
           -15            15             0             0             5             0         -4          0          0          0          0        0.5          0          0  #  SzSel_Male_Final_SummerN(2)
           -15            15             1             0             5             0         -4          0          0          0          0        0.5          0          0  #  SzSel_Male_Scale_SummerN(2)
# 3   WinterS LenSelex
            15            75       38.5376          43.1             5             0          2          0          0          0          0        0.5          1          1  #  Size_DblN_peak_WinterS(3)
            -5             3             3           0.7             5             0         -3          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_WinterS(3)
            -4            12       4.38359          3.42             5             0          3          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_WinterS(3)
            -2            15            14          0.21             5             0         -3          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_WinterS(3)
           -15             5          -999          -8.9             5             0         -4          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_WinterS(3)
            -5             5          -999          0.15             5             0         -4          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_WinterS(3)
            10            40       28.8241            15             9             0          2          0          0          0          0          0          2          1  #  Retain_L_infl_WinterS(3)
           0.1            10       1.36925             3             9             0          3          0          0          0          0          0          2          1  #  Retain_L_width_WinterS(3)
           -10            10       3.96554            10             9             0          4          0          0          0          0          0          2          2  #  Retain_L_asymptote_logit_WinterS(3)
           -10            10             0             0             9             0         -2          0          0          0          0          0          0          0  #  Retain_L_maleoffset_WinterS(3)
           -15            15      -12.4859             0             5             0          4          0          0          0          0        0.5          0          0  #  SzSel_Male_Peak_WinterS(3)
           -15            15      -1.82575             0             5             0          4          0          0          0          0        0.5          0          0  #  SzSel_Male_Ascend_WinterS(3)
           -15            15             0             0             5             0         -4          0          0          0          0        0.5          0          0  #  SzSel_Male_Descend_WinterS(3)
           -15            15             0             0             5             0         -4          0          0          0          0        0.5          0          0  #  SzSel_Male_Final_WinterS(3)
           -15            15             1             0             5             0         -4          0          0          0          0        0.5          0          0  #  SzSel_Male_Scale_WinterS(3)
# 4   SummerS LenSelex
            15            75       40.6489          43.1             5             0          2          0          0          0          0        0.5          1          1  #  Size_DblN_peak_SummerS(4)
            -5             3             3           0.7             5             0         -3          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_SummerS(4)
            -4            12       4.89926          3.42             5             0          3          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_SummerS(4)
            -2            15            14          0.21             5             0         -3          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_SummerS(4)
           -15             5          -999          -8.9             5             0         -4          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_SummerS(4)
            -5             5          -999          0.15             5             0         -4          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_SummerS(4)
            10            40       28.8597            15             9             0          3          0          0          0          0          0          3          1  #  Retain_L_infl_SummerS(4)
           0.1            10       1.07215             3             9             0          3          0          0          0          0          0          3          1  #  Retain_L_width_SummerS(4)
           -10            10       9.51624            10             9             0          4          0          0          0          0          0          3          2  #  Retain_L_asymptote_logit_SummerS(4)
           -10            10             0             0             9             0         -2          0          0          0          0          0          0          0  #  Retain_L_maleoffset_SummerS(4)
           -15            15      -12.5358             0             5             0          4          0          0          0          0        0.5          0          0  #  SzSel_Male_Peak_SummerS(4)
           -15            15      -1.89572             0             5             0          4          0          0          0          0        0.5          0          0  #  SzSel_Male_Ascend_SummerS(4)
           -15            15             0             0             5             0         -4          0          0          0          0        0.5          0          0  #  SzSel_Male_Descend_SummerS(4)
           -15            15             0             0             5             0         -4          0          0          0          0        0.5          0          0  #  SzSel_Male_Final_SummerS(4)
           -15            15             1             0             5             0         -4          0          0          0          0        0.5          0          0  #  SzSel_Male_Scale_SummerS(4)
# 5   TriEarly LenSelex
            15            61       35.2825          43.1             5             0          2          0          0          0          0        0.5          0          0  #  Size_DblN_peak_TriEarly(5)
            -5             3             3           0.7             5             0         -2          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_TriEarly(5)
            -4            12       4.20438          3.42             5             0          2          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_TriEarly(5)
            -2            15            14          0.21             5             0         -2          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_TriEarly(5)
           -15             5          -999          -8.9             5             0         -4          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_TriEarly(5)
            -5             5          -999          0.15             5             0         -4          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_TriEarly(5)
           -15            15      -3.83912             0             5             0          3          0          0          0          0        0.5          0          0  #  SzSel_Male_Peak_TriEarly(5)
           -15            15     -0.554495             0             5             0          3          0          0          0          0        0.5          0          0  #  SzSel_Male_Ascend_TriEarly(5)
           -15            15             0             0             5             0         -3          0          0          0          0        0.5          0          0  #  SzSel_Male_Descend_TriEarly(5)
           -15            15             0             0             5             0         -3          0          0          0          0        0.5          0          0  #  SzSel_Male_Final_TriEarly(5)
           -15            15             1             0             5             0         -4          0          0          0          0        0.5          0          0  #  SzSel_Male_Scale_TriEarly(5)
# 6   TriLate LenSelex
            15            61       36.4612          43.1             5             0          2          0          0          0          0        0.5          0          0  #  Size_DblN_peak_TriLate(6)
            -5             3             3           0.7             5             0         -2          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_TriLate(6)
            -4            12        4.6407          3.42             5             0          2          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_TriLate(6)
            -2            15            14          0.21             5             0         -2          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_TriLate(6)
           -15             5          -999          -8.9             5             0         -4          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_TriLate(6)
            -5             5          -999          0.15             5             0         -4          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_TriLate(6)
           -15            15       -2.2152             0             5             0          3          0          0          0          0        0.5          0          0  #  SzSel_Male_Peak_TriLate(6)
           -15            15    -0.0318622             0             5             0          3          0          0          0          0        0.5          0          0  #  SzSel_Male_Ascend_TriLate(6)
           -15            15             0             0             5             0         -3          0          0          0          0        0.5          0          0  #  SzSel_Male_Descend_TriLate(6)
           -15            15             0             0             5             0         -3          0          0          0          0        0.5          0          0  #  SzSel_Male_Final_TriLate(6)
           -15            15             1             0             5             0         -4          0          0          0          0        0.5          0          0  #  SzSel_Male_Scale_TriLate(6)
# 7   NWFSC LenSelex
            15            61       42.9488          43.1             5             0          2          0          0          0          0        0.5          0          0  #  Size_DblN_peak_NWFSC(7)
            -5             3             3           0.7             5             0         -2          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_NWFSC(7)
            -4            12       5.14849          3.42             5             0          2          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_NWFSC(7)
            -2            15            14          0.21             5             0         -2          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_NWFSC(7)
           -15             5          -999          -8.9             5             0         -4          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_NWFSC(7)
            -5             5          -999          0.15             5             0         -4          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_NWFSC(7)
           -15            15      -5.02527             0             5             0          3          0          0          0          0        0.5          0          0  #  SzSel_Male_Peak_NWFSC(7)
           -15            15     -0.406796             0             5             0          3          0          0          0          0        0.5          0          0  #  SzSel_Male_Ascend_NWFSC(7)
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
         -31.6          28.4       2.09015             0          14.2             6      5  # Size_DblN_peak_WinterN(1)_BLK1add_1973
         -31.6          28.4      -1.99134             0          14.2             6      5  # Size_DblN_peak_WinterN(1)_BLK1add_1983
         -31.6          28.4     -0.797911             0          14.2             6      5  # Size_DblN_peak_WinterN(1)_BLK1add_1993
         -31.6          28.4      0.381086             0          14.2             6      5  # Size_DblN_peak_WinterN(1)_BLK1add_2003
         -31.6          28.4      0.875568             0          14.2             6      5  # Size_DblN_peak_WinterN(1)_BLK1add_2011
        -16.19         13.81      -2.28584             0         6.905             6      5  # Retain_L_infl_WinterN(1)_BLK2add_2003
        -16.19         13.81       1.74534             0         6.905             6      5  # Retain_L_infl_WinterN(1)_BLK2add_2010
        -16.19         13.81      -3.25655             0         6.905             6      5  # Retain_L_infl_WinterN(1)_BLK2add_2011
        -1.601         8.299      0.122615             0        0.8005             6      5  # Retain_L_width_WinterN(1)_BLK2add_2003
        -1.601         8.299      0.393058             0        0.8005             6      5  # Retain_L_width_WinterN(1)_BLK2add_2010
        -1.601         8.299     -0.697336             0        0.8005             6      5  # Retain_L_width_WinterN(1)_BLK2add_2011
           -10            10       6.63553            10             9             0      5  # Retain_L_asymptote_logit_WinterN(1)_BLK2repl_2003
           -10            10       2.11173            10             9             0      5  # Retain_L_asymptote_logit_WinterN(1)_BLK2repl_2010
           -10            10       9.98781            10             9             0      5  # Retain_L_asymptote_logit_WinterN(1)_BLK2repl_2011
         -38.8          21.2       1.97543             0          10.6             6      5  # Size_DblN_peak_SummerN(2)_BLK1add_1973
         -38.8          21.2     -0.426336             0          10.6             6      5  # Size_DblN_peak_SummerN(2)_BLK1add_1983
         -38.8          21.2      -2.39885             0          10.6             6      5  # Size_DblN_peak_SummerN(2)_BLK1add_1993
         -38.8          21.2     -0.112112             0          10.6             6      5  # Size_DblN_peak_SummerN(2)_BLK1add_2003
         -38.8          21.2       3.25066             0          10.6             6      5  # Size_DblN_peak_SummerN(2)_BLK1add_2011
       -20.679         9.321      -0.42375             0        4.6605             6      5  # Retain_L_infl_SummerN(2)_BLK3add_2003
       -20.679         9.321       1.29123             0        4.6605             6      5  # Retain_L_infl_SummerN(2)_BLK3add_2009
       -20.679         9.321      -1.92402             0        4.6605             6      5  # Retain_L_infl_SummerN(2)_BLK3add_2011
       -1.0278        8.8722      0.161921             0        0.5139             6      5  # Retain_L_width_SummerN(2)_BLK3add_2003
       -1.0278        8.8722      0.138787             0        0.5139             6      5  # Retain_L_width_SummerN(2)_BLK3add_2009
       -1.0278        8.8722      0.210417             0        0.5139             6      5  # Retain_L_width_SummerN(2)_BLK3add_2011
           -10            10       5.45837            10             9             0      5  # Retain_L_asymptote_logit_SummerN(2)_BLK3repl_2003
           -10            10       7.55493            10             9             0      5  # Retain_L_asymptote_logit_SummerN(2)_BLK3repl_2009
           -10            10       6.15652            10             9             0      5  # Retain_L_asymptote_logit_SummerN(2)_BLK3repl_2011
       -25.422        34.578       -9.1907             0        12.711             6      5  # Size_DblN_peak_WinterS(3)_BLK1add_1973
       -25.422        34.578       4.89511             0        12.711             6      5  # Size_DblN_peak_WinterS(3)_BLK1add_1983
       -25.422        34.578       8.89684             0        12.711             6      5  # Size_DblN_peak_WinterS(3)_BLK1add_1993
       -25.422        34.578       6.70388             0        12.711             6      5  # Size_DblN_peak_WinterS(3)_BLK1add_2003
       -25.422        34.578       8.02639             0        12.711             6      5  # Size_DblN_peak_WinterS(3)_BLK1add_2011
       -18.816        11.184      -1.97623             0         5.592             6      5  # Retain_L_infl_WinterS(3)_BLK2add_2003
       -18.816        11.184       1.57306             0         5.592             6      5  # Retain_L_infl_WinterS(3)_BLK2add_2010
       -18.816        11.184      -4.34263             0         5.592             6      5  # Retain_L_infl_WinterS(3)_BLK2add_2011
       -1.0443        8.8557      0.360241             0       0.52215             6      5  # Retain_L_width_WinterS(3)_BLK2add_2003
       -1.0443        8.8557      0.136063             0       0.52215             6      5  # Retain_L_width_WinterS(3)_BLK2add_2010
       -1.0443        8.8557     -0.057897             0       0.52215             6      5  # Retain_L_width_WinterS(3)_BLK2add_2011
           -10            10       7.90011            10             9             0      5  # Retain_L_asymptote_logit_WinterS(3)_BLK2repl_2003
           -10            10       5.51172            10             9             0      5  # Retain_L_asymptote_logit_WinterS(3)_BLK2repl_2010
           -10            10       7.95146            10             9             0      5  # Retain_L_asymptote_logit_WinterS(3)_BLK2repl_2011
      -28.0793       31.9207        -5.094             0       14.0397             6      5  # Size_DblN_peak_SummerS(4)_BLK1add_1973
      -28.0793       31.9207      -6.33587             0       14.0397             6      5  # Size_DblN_peak_SummerS(4)_BLK1add_1983
      -28.0793       31.9207         3.536             0       14.0397             6      5  # Size_DblN_peak_SummerS(4)_BLK1add_1993
      -28.0793       31.9207       6.22288             0       14.0397             6      5  # Size_DblN_peak_SummerS(4)_BLK1add_2003
      -28.0793       31.9207       6.00477             0       14.0397             6      5  # Size_DblN_peak_SummerS(4)_BLK1add_2011
       -19.055        10.945      -1.38928             0        5.4725             6      5  # Retain_L_infl_SummerS(4)_BLK3add_2003
       -19.055        10.945      -1.66543             0        5.4725             6      5  # Retain_L_infl_SummerS(4)_BLK3add_2009
       -19.055        10.945      -2.06749             0        5.4725             6      5  # Retain_L_infl_SummerS(4)_BLK3add_2011
        -0.876         9.024      0.604245             0         0.438             6      5  # Retain_L_width_SummerS(4)_BLK3add_2003
        -0.876         9.024      0.470079             0         0.438             6      5  # Retain_L_width_SummerS(4)_BLK3add_2009
        -0.876         9.024      0.580645             0         0.438             6      5  # Retain_L_width_SummerS(4)_BLK3add_2011
           -10            10       7.55709            10             9             0      5  # Retain_L_asymptote_logit_SummerS(4)_BLK3repl_2003
           -10            10       8.90125            10             9             0      5  # Retain_L_asymptote_logit_SummerS(4)_BLK3repl_2009
           -10            10       7.67804            10             9             0      5  # Retain_L_asymptote_logit_SummerS(4)_BLK3repl_2011
# info on dev vectors created for selex parms are reported with other devs after tag parameter section 
#
0   #  use 2D_AR1 selectivity(0/1):  experimental feature
#_no 2D_AR1 selex offset used
#
# Tag loss and Tag reporting parameters go next
0  # TG_custom:  0=no read; 1=read if tags exist
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
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
0 # (0/1) read specs for more stddev reporting 
 # 0 0 0 0 0 0 0 0 0 # placeholder for # selex_fleet, 1=len/2=age/3=both, year, N selex bins, 0 or Growth pattern, N growth ages, 0 or NatAge_area(-1 for all), NatAge_yr, N Natages
 # placeholder for vector of selex bins to be reported
 # placeholder for vector of growth ages to be reported
 # placeholder for vector of NatAges ages to be reported
999

