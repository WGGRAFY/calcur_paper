#install
install.packages("here")
remotes::install_github("r4ss/r4ss", ref="development")
require(r4ss)
remotes::install_github("ss3sim/ss3sim",
                        ref = "main", build_vignettes = TRUE, dependencies = TRUE)
require(ss3sim)
require(purrr)
require(here)
require(dplyr)
require(tidyr)
require(ggplot2)


#look for ss.exe
system.file("bin", package = "ss3sim")
ss3sim_dir <- paste0(substr(here(),1,32),"ss3sim")
devtools::load_all(ss3sim_dir)
### Creating scenarios
# Scenario are based on:
# 1)	Age/size selective sampling
# 2)	Cohort and initial size effects
# 3)	Time-varying maturation
# 4)	Density-dependent growth effects
# 5)	Fishing rate

### What do we need? length at age information



## Now read in the dat file to extract some important info
# Directory needs to be named "om" to work with ss3sim
  Datfile = readLines(here("om/om.dat"))
  Nage = as.numeric(substr(Datfile[grep("#_Nages", Datfile, fixed=TRUE)], 1, 2))
  # if length bin method is 2 (like in this case) do below:
  Nlength =  (as.numeric(substr(Datfile[grep("# maximum size ", Datfile, fixed=TRUE)], 1, 3)) - as.numeric(substr(Datfile[grep("# minimum size ", Datfile, fixed=TRUE)], 1, 2)))/as.numeric(substr(Datfile[grep("# binwidth for population", Datfile, fixed=TRUE)], 1, 2))+1
  StartYr = as.numeric(substr(Datfile[grep("#_StartYr", Datfile, fixed=TRUE)], 1, 1))
  EndYr = as.numeric(substr(Datfile[grep("#_EndYr", Datfile, fixed=TRUE)], 1, 3))
  
  source(here("helper_functions.R"))
  
  ## now read in the Report.sso file
  rawrep <- read.table(file=here("om/Report.sso"),col.names=1:200,fill=TRUE,quote="", colClasses="character",nrows=-1,comment.char="")
  
  # Extract the general place where the AGE_LENGTH_KEY is written 
  LAA <- matchfun2(string1 ="MEAN_SIZE_TIMESERIES", adjust1=1, string2="mean_size_Jan_1_for_sex", adjust2=-2, header = TRUE)  %>% as_tibble()
  
  #Split by morph and subseas
  LAA_list <- LAA %>% filter(Yr>=as.numeric(StartYr), Yr<=as.numeric(EndYr)) %>% split(list(.$Morph,.$SubSeas))
  save(LAA_list,file="petraleLAA.Rds")
  
  # We then need to add some observation error on top of them to create some data to feed in the estimation model
  CV <- 0.1
  SD = sqrt(log(CV^2+1)) # for lognromal distribution
  LAA_dat <- t(apply(LAA_list[[1]][,-c(1:4)],1,function(x) as.numeric(as.character(x))))*rlnorm(prod(dim(LAA_list[[1]][,-c(1:4)])), 0, SD)
  
    #Read in the values for L at a max
    
    input <- lapply(here(paste0("om/G",c("M","F"),"3-pet.csv")), read.csv, header=FALSE)
    input[[2]]$Fval <-scale(input[[2]]$V1) 
    setwd(here("om"))
    
    # 1. Change in VBk (monotonic trend, decadal) (cases 1-2)
    # 2. Change in VB L inf  (monotonic trend, decadal)(cases 3-4)
    # 3. Change in maturation parameters (A50, L50) (TODO)
    # 4. CV of size at age  (TODO)
    # 5. Sample size of survey per year (TODO)
    # 6. Sample size of fishery data (TODO) 
    # 7. Selectivity  (TODO)
    # 8. Level of recruitment variation (CV of recruitment) (TODO)
    # 9. Level of fishing pressure (TODO)
    # 10. CV of survey index (TODO)
    df <-make_data_in(tv_Linf = "1")
   
    
    iterations <- 1
    
    
    arg_list <- setup_scenarios(df)
    # arg_list <- setup_scenarios(setup_scenarios_defaults())
    
    # Below, I am using the "new" ss3sim_base function where I have added the "skim_em" argument to skip the EM part. 
    # As we do not need the EM part, this saves us some time. 

    source(here("ss3sim_base_new.R"))
    
    out <- lapply(arg_list, function(x) {
      do.call("ss3sim_base_new", c(x, list(iterations = iterations), skip_em=TRUE))
    })
    unlink(df[,"scenarios"], recursive=TRUE)
    
    output_path <- file.path("om/D1-L1-A0-K0-M0-R0-F0-C0-pet", "1")
    
    #Compare output OM and EM
    r_om <- r4ss::SS_output(here(output_path,"om"),
                            verbose = FALSE, printstats = FALSE, covar = FALSE)
    #Compare output EM
    r_em <- r4ss::SS_readdat(file.path(output_path,"em","ss3.dat"),
                            verbose = FALSE)
    r4ss::SS_plots(r_om)
    
    #EM data file has the simulated comps we need
    em_dat<- r4ss::SS_readdat_3.30(here(output_path,"em/ss3.dat"))
    
    #Function to turn length comps and CAAL comps to samples
    Length_sims_1 <- turn_ss_out_to_data(em_dat, flt=1)

