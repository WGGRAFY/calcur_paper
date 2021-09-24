#V3.30
#C file created using the SS_writectl function in the R package r4ss
#C file write time: 2021-09-23 13:34:28
#
0 # 0 means do not read wtatage.ss; 1 means read and usewtatage.ss and also read and use growth parameters
1 #_N_Growth_Patterns
1 #_N_platoons_Within_GrowthPattern
2 # recr_dist_method for parameters
1 # not yet implemented; Future usage:Spawner-Recruitment; 1=global; 2=by area
1 # number of recruitment settlement assignments 
0 # unused option
# for each settlement assignment:
#_GPattern	month	area	age
1	1	1	0	#_recr_dist_pattern1
#
#_Cond 0 # N_movement_definitions goes here if N_areas > 1
#_Cond 1.0 # first age that moves (real age at begin of season, not integer) also cond on do_migration>0
#_Cond 1 1 1 2 4 10 # example move definition for seas=1, morph=1, source=1 dest=2, age1=4, age2=10
#
0 #_Nblock_Patterns
#_Cond 0 #_blocks_per_pattern
# begin and end years of blocks
#
# controls for all timevary parameters 
1 #_env/block/dev_adjust_method for all time-vary parms (1=warn relative to base parm bounds; 3=no bound check)
#
# AUTOGEN
0 0 0 0 0 # autogen: 1st element for biology, 2nd for SR, 3rd for Q, 4th reserved, 5th for selex
# where: 0 = autogen all time-varying parms; 1 = read each time-varying parm line; 2 = read then autogen if parm min==-12345
#
# setup for M, growth, maturity, fecundity, recruitment distibution, movement
#
0 #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate
#_no additional input for selected M option; read 1P per morph
1 # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K_incr; 4=age_specific_K_decr;5=age_specific_K_each; 6=NA; 7=NA; 8=growth cessation
1 #_Age(post-settlement)_for_L1;linear growth below this
999 #_Growth_Age_for_L2 (999 to use as Linf)
-999 #_exponential decay for growth above maxage (value should approx initial Z; -999 replicates 3.24; -998 to not allow growth above maxage)
0 #_placeholder for future growth feature
#
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
0 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
1 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
0 #_First_Mature_Age
1 #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0 #_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
1 #_parameter_offset_approach (1=none, 2= M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)
#
#_growth_parms
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env_var&link	dev_link	dev_minyr	dev_maxyr	dev_PH	Block	Block_Fxn
 1.0e-02	  1.800000	 2.000e-01	1.000e-01	0.8	0	 -3	0	0	0	0	0	0	0	#_NatM_p_1_Fem_GP_1        
 1.0e+00	100.000000	 2.000e+01	3.080e+01	0.2	0	  4	0	0	0	0	0	0	0	#_L_at_Amin_Fem_GP_1       
 6.6e+00	660.000000	 1.320e+02	1.201e+02	0.2	0	  4	0	0	0	0	0	0	0	#_L_at_Amax_Fem_GP_1       
 1.0e-03	  1.000000	 2.000e-01	2.500e-01	0.8	0	  4	0	0	0	0	0	0	0	#_VonBert_K_Fem_GP_1       
 1.0e-03	  0.500000	 1.000e-01	1.000e-01	0.8	0	  5	0	0	0	0	0	0	0	#_CV_young_Fem_GP_1        
 1.0e-03	  0.500000	 1.000e-01	1.000e-01	0.8	0	  5	0	0	0	0	0	0	0	#_CV_old_Fem_GP_1          
 0.0e+00	  3.000000	 6.800e-06	6.800e-06	0.0	0	 -3	0	0	0	0	0	0	0	#_Wtlen_1_Fem_GP_1         
 2.5e+00	  3.500000	 3.101e+00	3.101e+00	0.2	0	 -3	0	0	0	0	0	0	0	#_Wtlen_2_Fem_GP_1         
 1.0e+01	 50.000000	 3.818e+01	0.000e+00	0.0	0	-99	0	0	0	0	0	0	0	#_Mat50%_Fem_GP_1          
-2.0e+00	  2.000000	-2.760e-01	0.000e+00	0.0	0	-99	0	0	0	0	0	0	0	#_Mat_slope_Fem_GP_1       
-3.0e+00	  3.000000	 1.000e+00	0.000e+00	0.0	0	-99	0	0	0	0	0	0	0	#_Eggs/kg_inter_Fem_GP_1   
-3.0e+00	  4.000000	 0.000e+00	0.000e+00	0.0	0	-99	0	0	0	0	0	0	0	#_Eggs/kg_slope_wt_Fem_GP_1
-4.0e+00	  4.000000	 0.000e+00	0.000e+00	0.0	0	-99	0	0	0	0	0	0	0	#_RecrDist_GP_1            
-4.0e+00	  4.000000	 0.000e+00	0.000e+00	0.0	0	-99	0	0	0	0	0	0	0	#_RecrDist_Area_1          
-4.0e+00	  4.000000	 0.000e+00	0.000e+00	0.0	0	-99	0	0	0	0	0	0	0	#_RecrDist_month_1         
-4.0e+00	  4.000000	 1.000e+00	0.000e+00	0.0	0	-99	0	0	0	0	0	0	0	#_CohortGrowDev            
 1.0e-06	  0.999999	 5.000e-01	5.000e-01	0.5	0	-99	0	0	0	0	0	0	0	#_FracFemale_GP_1          
#_no timevary MG parameters
#
#_seasonal_effects_on_biology_parms
0 0 0 0 0 0 0 0 0 0 #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K
#_ LO HI INIT PRIOR PR_SD PR_type PHASE
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no seasonal MG parameters
#
3 #_Spawner-Recruitment; 2=Ricker; 3=std_B-H; 4=SCAA;5=Hockey; 6=B-H_flattop; 7=survival_3Parm;8=Shepard_3Parm
0 # 0/1 to use steepness in initial equ recruitment calculation
0 # future feature: 0/1 to make realized sigmaR a function of SR curvature
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn # parm_name
 4.0	20	18.70	10.3	10.00	0	  1	0	0	0	0	0	0	0	#_SR_LN(R0)  
 0.2	 1	 0.65	 0.7	 0.05	0	 -3	0	0	0	0	0	0	0	#_SR_BH_steep
 0.0	 2	 0.40	 0.8	 0.80	0	-99	0	0	0	0	0	0	0	#_SR_sigmaR  
-5.0	 5	 0.00	 0.0	 1.00	0	-99	0	0	0	0	0	0	0	#_SR_regime  
 0.0	 0	 0.00	 0.0	 0.00	0	 -6	0	0	0	0	0	0	0	#_SR_autocorr
#_no timevary SR parameters
1 #do_recdev:  0=none; 1=devvector (R=F(SSB)+dev); 2=deviations (R=F(SSB)+dev); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
1 # first year of main recr_devs; early devs can preceed this era
100 # last year of main recr_devs; forecast devs start in following year
-2 #_recdev phase
1 # (0/1) to read 13 advanced options
0 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
-4 #_recdev_early_phase
0 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
1 #_lambda for Fcast_recr_like occurring before endyr+1
1 #_last_yr_nobias_adj_in_MPD; begin of ramp
1 #_first_yr_fullbias_adj_in_MPD; begin of plateau
100 #_last_yr_fullbias_adj_in_MPD
100 #_end_yr_for_ramp_in_MPD (can be in forecast to shape ramp, but SS sets bias_adj to 0.0 for fcast yrs)
0 #_max_bias_adj_in_MPD (-1 to override ramp and set biasadj=1.0 for all estimated recdevs)
0 #_period of cycles in recruitment (N parms read below)
-10 #min rec_dev
10 #max rec_dev
101 #_read_recdevs
#_end of advanced SR options
#
#_placeholder for full parameter lines for recruitment cycles
#_Year	recdev
  1	-0.69543204	#_recdev_input1  
  2	 0.40155322	#_recdev_input2  
  3	-0.11989007	#_recdev_input3  
  4	 0.44542403	#_recdev_input4  
  5	 0.09875037	#_recdev_input5  
  6	 0.14913202	#_recdev_input6  
  7	-0.52834561	#_recdev_input7  
  8	-0.09010645	#_recdev_input8  
  9	 0.30076704	#_recdev_input9  
 10	-0.23368550	#_recdev_input10 
 11	 0.09008947	#_recdev_input11 
 12	-0.85262873	#_recdev_input12 
 13	 0.71960186	#_recdev_input13 
 14	-0.04725786	#_recdev_input14 
 15	-0.58730660	#_recdev_input15 
 16	-0.14559294	#_recdev_input16 
 17	 0.34612772	#_recdev_input17 
 18	-0.95711702	#_recdev_input18 
 19	-0.14140997	#_recdev_input19 
 20	-0.54420166	#_recdev_input20 
 21	-0.25381445	#_recdev_input21 
 22	 0.03755734	#_recdev_input22 
 23	 0.09947142	#_recdev_input23 
 24	-0.24477016	#_recdev_input24 
 25	 0.06190721	#_recdev_input25 
 26	 0.07043711	#_recdev_input26 
 27	 0.76999696	#_recdev_input27 
 28	 0.10140199	#_recdev_input28 
 29	 0.75154202	#_recdev_input29 
 30	 0.39113510	#_recdev_input30 
 31	-0.21230768	#_recdev_input31 
 32	 0.45259482	#_recdev_input32 
 33	-0.42707020	#_recdev_input33 
 34	-1.17946958	#_recdev_input34 
 35	-0.77180126	#_recdev_input35 
 36	 0.08221242	#_recdev_input36 
 37	 0.11506277	#_recdev_input37 
 38	 0.16482754	#_recdev_input38 
 39	-0.21525646	#_recdev_input39 
 40	-0.00168606	#_recdev_input40 
 41	-0.10350247	#_recdev_input41 
 42	 0.37938474	#_recdev_input42 
 43	 0.30446763	#_recdev_input43 
 44	-0.15882585	#_recdev_input44 
 45	 0.00788690	#_recdev_input45 
 46	 0.53776510	#_recdev_input46 
 47	-0.14926907	#_recdev_input47 
 48	 0.09569090	#_recdev_input48 
 49	-0.30489262	#_recdev_input49 
 50	-0.03287142	#_recdev_input50 
 51	-0.65363353	#_recdev_input51 
 52	-0.82200094	#_recdev_input52 
 53	-0.28992383	#_recdev_input53 
 54	-0.23754570	#_recdev_input54 
 55	 0.05153852	#_recdev_input55 
 56	 0.05008056	#_recdev_input56 
 57	 0.01637693	#_recdev_input57 
 58	 0.59864222	#_recdev_input58 
 59	-0.09306596	#_recdev_input59 
 60	-0.92008426	#_recdev_input60 
 61	-0.28953630	#_recdev_input61 
 62	-0.75012287	#_recdev_input62 
 63	 0.08061874	#_recdev_input63 
 64	-0.18164443	#_recdev_input64 
 65	 0.42223309	#_recdev_input65 
 66	-0.34078488	#_recdev_input66 
 67	-0.15975218	#_recdev_input67 
 68	 0.20999046	#_recdev_input68 
 69	 0.28494684	#_recdev_input69 
 70	-0.00424632	#_recdev_input70 
 71	 0.25446390	#_recdev_input71 
 72	-0.39439575	#_recdev_input72 
 73	-0.22346661	#_recdev_input73 
 74	-0.72201965	#_recdev_input74 
 75	 0.17477281	#_recdev_input75 
 76	 0.55191783	#_recdev_input76 
 77	-0.12310962	#_recdev_input77 
 78	-0.18685552	#_recdev_input78 
 79	 0.34363791	#_recdev_input79 
 80	 0.18795634	#_recdev_input80 
 81	-0.58879414	#_recdev_input81 
 82	-0.20234175	#_recdev_input82 
 83	-0.33299625	#_recdev_input83 
 84	-0.75299177	#_recdev_input84 
 85	 0.11234514	#_recdev_input85 
 86	 0.17053662	#_recdev_input86 
 87	-0.50985624	#_recdev_input87 
 88	-0.03706445	#_recdev_input88 
 89	 0.00403095	#_recdev_input89 
 90	-0.33218636	#_recdev_input90 
 91	-0.11547221	#_recdev_input91 
 92	 0.13270286	#_recdev_input92 
 93	-0.29405885	#_recdev_input93 
 94	-0.55842264	#_recdev_input94 
 95	 0.11047923	#_recdev_input95 
 96	-0.29014330	#_recdev_input96 
 97	 0.27905948	#_recdev_input97 
 98	 0.64669549	#_recdev_input98 
 99	 0.25063171	#_recdev_input99 
100	-0.35081508	#_recdev_input100
101	-0.21360362	#_recdev_input101
#
#Fishing Mortality info
0.3 # F ballpark
-2001 # F ballpark year (neg value to disable)
2 # F_Method:  1=Pope; 2=instan. F; 3=hybrid (hybrid is recommended)
4 # max F or harvest rate, depends on F_Method
#_overall start F value; overall phase; N detailed inputs to read
0 1 61 #_F_setup
#_fleet	yr	seas	Fvalue	se	phase
1	1960	1	0.1	0.005	1	#_F_setup21 
1	1961	1	0.1	0.005	1	#_F_setup22 
1	1962	1	0.1	0.005	1	#_F_setup23 
1	1963	1	0.1	0.005	1	#_F_setup24 
1	1964	1	0.1	0.005	1	#_F_setup25 
1	1965	1	0.1	0.005	1	#_F_setup26 
1	1966	1	0.1	0.005	1	#_F_setup27 
1	1967	1	0.1	0.005	1	#_F_setup28 
1	1968	1	0.1	0.005	1	#_F_setup29 
1	1969	1	0.1	0.005	1	#_F_setup210
1	1970	1	0.1	0.005	1	#_F_setup211
1	1971	1	0.1	0.005	1	#_F_setup212
1	1972	1	0.1	0.005	1	#_F_setup213
1	1973	1	0.1	0.005	1	#_F_setup214
1	1974	1	0.1	0.005	1	#_F_setup215
1	1975	1	0.1	0.005	1	#_F_setup216
1	1976	1	0.1	0.005	1	#_F_setup217
1	1977	1	0.1	0.005	1	#_F_setup218
1	1978	1	0.1	0.005	1	#_F_setup219
1	1979	1	0.1	0.005	1	#_F_setup220
1	1980	1	0.1	0.005	1	#_F_setup221
1	1981	1	0.1	0.005	1	#_F_setup222
1	1982	1	0.1	0.005	1	#_F_setup223
1	1983	1	0.1	0.005	1	#_F_setup224
1	1984	1	0.1	0.005	1	#_F_setup225
1	1985	1	0.1	0.005	1	#_F_setup226
1	1986	1	0.1	0.005	1	#_F_setup227
1	1987	1	0.1	0.005	1	#_F_setup228
1	1988	1	0.1	0.005	1	#_F_setup229
1	1989	1	0.1	0.005	1	#_F_setup230
1	1990	1	0.1	0.005	1	#_F_setup231
1	1991	1	0.1	0.005	1	#_F_setup232
1	1992	1	0.1	0.005	1	#_F_setup233
1	1993	1	0.1	0.005	1	#_F_setup234
1	1994	1	0.1	0.005	1	#_F_setup235
1	1995	1	0.1	0.005	1	#_F_setup236
1	1996	1	0.1	0.005	1	#_F_setup237
1	1997	1	0.1	0.005	1	#_F_setup238
1	1998	1	0.1	0.005	1	#_F_setup239
1	1999	1	0.1	0.005	1	#_F_setup240
1	2000	1	0.1	0.005	1	#_F_setup241
1	2001	1	0.1	0.005	1	#_F_setup242
1	2002	1	0.1	0.005	1	#_F_setup243
1	2003	1	0.1	0.005	1	#_F_setup244
1	2004	1	0.1	0.005	1	#_F_setup245
1	2005	1	0.1	0.005	1	#_F_setup246
1	2006	1	0.1	0.005	1	#_F_setup247
1	2007	1	0.1	0.005	1	#_F_setup248
1	2008	1	0.1	0.005	1	#_F_setup249
1	2009	1	0.1	0.005	1	#_F_setup250
1	2010	1	0.1	0.005	1	#_F_setup251
1	2011	1	0.1	0.005	1	#_F_setup252
1	2012	1	0.1	0.005	1	#_F_setup253
1	2013	1	0.1	0.005	1	#_F_setup254
1	2014	1	0.1	0.005	1	#_F_setup255
1	2015	1	0.1	0.005	1	#_F_setup256
1	2016	1	0.1	0.005	1	#_F_setup257
1	2017	1	0.1	0.005	1	#_F_setup258
1	2018	1	0.1	0.005	1	#_F_setup259
1	2019	1	0.1	0.005	1	#_F_setup260
1	2020	1	0.1	0.005	1	#_F_setup261
#
#_initial_F_parms; count = 0
#
#_Q_setup for fleets with cpue or survey data
#_fleet	link	link_info	extra_se	biasadj	float  #  fleetname
    2	1	0	0	0	0	#_Survey    
-9999	0	0	0	0	0	#_terminator
#_Q_parms(if_any);Qunits_are_ln(q)
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn  #  parm_name
-20	20	0	0	99	0	-5	0	0	0	0	0	0	0	#_LnQ_base_Survey(2)
#_no timevary Q parameters
#
#_size_selex_patterns
#_Pattern	Discard	Male	Special
24	0	0	0	#_1 Fishery
24	0	0	0	#_2 Survey 
#
#_age_selex_patterns
#_Pattern	Discard	Male	Special
11	0	0	0	#_1 Fishery
11	0	0	0	#_2 Survey 
#
#_SizeSelex
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn  #  parm_name
   5	 102	  50.8	  50.8	0.05	0	  2	0	0	0	0	0.5	0	0	#_SizeSel_P_1_Fishery(1)
  -5	   3	  -3.0	  -3.0	0.05	0	-99	0	0	0	0	0.5	0	0	#_SizeSel_P_2_Fishery(1)
   0	  26	   5.1	   5.1	0.05	0	  3	0	0	0	0	0.5	0	0	#_SizeSel_P_3_Fishery(1)
  -2	  16	  15.0	  15.0	0.05	0	-99	0	0	0	0	0.5	0	0	#_SizeSel_P_4_Fishery(1)
 -15	   5	-999.0	-999.0	0.05	0	-99	0	0	0	0	0.5	0	0	#_SizeSel_P_5_Fishery(1)
  -5	1000	 999.0	 999.0	0.05	0	-99	0	0	0	0	0.5	0	0	#_SizeSel_P_6_Fishery(1)
  20	 102	  41.8	  41.8	0.05	0	  2	0	0	0	0	0.5	0	0	#_SizeSel_P_1_Survey(2) 
  -5	   3	  -4.0	  -4.0	0.05	0	-99	0	0	0	0	0.5	0	0	#_SizeSel_P_2_Survey(2) 
  -4	  26	   5.2	   5.2	0.05	0	  3	0	0	0	0	0.5	0	0	#_SizeSel_P_3_Survey(2) 
  -2	  16	  14.0	  14.0	0.05	0	-99	0	0	0	0	0.5	0	0	#_SizeSel_P_4_Survey(2) 
-100	 100	 -99.0	 -99.0	0.05	0	-99	0	0	0	0	0.5	0	0	#_SizeSel_P_5_Survey(2) 
-100	 100	  99.0	  99.0	0.05	0	-99	0	0	0	0	0.5	0	0	#_SizeSel_P_6_Survey(2) 
#_AgeSelex
0	  1	  0.1	  0.1	99	0	-3	0	0	0	0	0.5	0	0	#_AgeSel_P_1_Fishery(1)
0	101	100.0	100.0	99	0	-3	0	0	0	0	0.5	0	0	#_AgeSel_P_2_Fishery(1)
0	  1	  0.1	  0.1	99	0	-3	0	0	0	0	0.5	0	0	#_AgeSel_P_1_Survey(2) 
0	101	100.0	100.0	99	0	-3	0	0	0	0	0.5	0	0	#_AgeSel_P_2_Survey(2) 
#_no timevary selex parameters
#
0 #  use 2D_AR1 selectivity(0/1):  experimental feature
#_no 2D_AR1 selex offset used
# Tag loss and Tag reporting parameters go next
0 # TG_custom:  0=no read; 1=read if tags exist
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
#
# Input variance adjustments factors: 
#_Factor Fleet Value
-9999 1 0 # terminator
#
4 #_maxlambdaphase
1 #_sd_offset; must be 1 if any growthCV, sigmaR, or survey extraSD is an estimated parameter
# read 0 changes to default Lambdas (default value is 1.0)
-9999 0 0 0 0 # terminator
#
0 # 0/1 read specs for more stddev reporting
#
999
