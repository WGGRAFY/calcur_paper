#install
install.packages("here")
remotes::install_github("r4ss/r4ss", ref="development")
require(r4ss)
remotes::install_github("ss3sim/ss3sim",
                        ref = "main", build_vignettes = TRUE, dependencies = TRUE)
require(ss3sim)
require(purrr)
require(here)


#look for ss.exe
system.file("bin", package = "ss3sim")


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
  StartYr = as.numeric(substr(Datfile[grep("#_StartYr", Datfile, fixed=TRUE)], 1, 4))
  EndYr = as.numeric(substr(Datfile[grep("#_EndYr", Datfile, fixed=TRUE)], 1, 4))
  
## now read in the Report.sso file
  rawrep <- read.table(file=paste0(getwd(), "/calcur_paper-master1/petrale_ss3sim/om_test/Report.sso"),col.names=1:200,fill=TRUE,quote="", colClasses="character",nrows=-1,comment.char="")

# Extract the general place where the AGE_LENGTH_KEY is written 
  LAA <- matchfun2(string1 ="MEAN_SIZE_TIMESERIES", adjust1=1, string2="mean_size_Jan_1_for_sex", adjust2=-2, header = TRUE)

  LAA %>% purrr::map(c(1,2), filter(Morph==.x))
# size from subseason 1 and morph 1
  LAA_sea1_morph1 <- LAA %>% filter(Morph==1, SubSeas==1, Yr>=StartYr, Yr <= EndYr) 
 
# size from subseason 1 and morph 2
  LAA_sea1_morph2 <- LAA %>% filter(Morph==1, SubSeas==2, Yr>=StartYr, Yr <= EndYr)
  
# size from subseason 2 and morph 1
  LAA_sea2_morph1 <- LAA %>% filter(Morph==2, SubSeas==1, Yr>=StartYr, Yr <= EndYr)
  
# size from subseason 2 and morph 2
  LAA_sea2_morph2 <- LAA %>% filter(Morph==2, SubSeas==2, Yr>=StartYr, Yr <= EndYr)
  

# We then need to add some observation error on top of them to create some data to feed in the estimation model
  
  CV <- 0.1
  SD = sqrt(log(CV^2+1)) # for lognromal distribution
  LAA_dat <- t(apply(LAA_sea1_morph1[,-c(1:4)],1,function(x) as.numeric(as.character(x))))*rlnorm(prod(dim(LAA_sea1_morph1[,-c(1:4)])), 0, SD)
  
  
  

df <- setup_scenarios_defaults()
df[, "si.years.2"] <- "seq(76, 100, by = 2)"
df <- rbind(df, df)
df[, "si.sds_obs.2"] <- c(0.1, 0.4)
df[, "bias_adjust"] <- FALSE
df[, "hess_always"] <- TRUE