#install
remotes::install_github("r4ss/r4ss", ref="development")
require(r4ss)
remotes::install_github("ss3sim/ss3sim",
                        ref = "main", build_vignettes = TRUE, dependencies = TRUE)
require(ss3sim)


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


## Now some function to extract the LAA info from the SS3 runs 
## (the 2 functions below are taken from https://github.com/r4ss/testing/blob/master/r4ss%20tables.r)
  emptytest <- function(x){ sum(!is.na(x) & x=="")/length(x) }
  
  matchfun2 <- function(string1,adjust1,string2,adjust2,cols="nonblank",matchcol1=1,matchcol2=1,
                        objmatch=rawrep,objsubset=rawrep,substr1=TRUE,substr2=TRUE,header=FALSE)
  {
    # return a subset of values from the report file (or other file)
    # subset is defined by character strings at the start and end, with integer
    # adjustments of the number of lines to above/below the two strings
    line1 <- match(string1,if(substr1){substring(objmatch[,matchcol1],1,nchar(string1))}else{objmatch[,matchcol1]})
    line2 <- match(string2,if(substr2){substring(objmatch[,matchcol2],1,nchar(string2))}else{objmatch[,matchcol2]})
    if(is.na(line1) | is.na(line2)) return("absent")
    
    if(is.numeric(cols))    out <- objsubset[(line1+adjust1):(line2+adjust2),cols]
    if(cols[1]=="all")      out <- objsubset[(line1+adjust1):(line2+adjust2),]
    if(cols[1]=="nonblank"){
      # returns only columns that contain at least one non-empty value
      out <- objsubset[(line1+adjust1):(line2+adjust2),]
      out <- out[,apply(out,2,emptytest) < 1]
    }
    if(header && nrow(out)>0){
      out[1,out[1,]==""] <- "NoName"
      names(out) <- out[1,]
      out <- out[-1,]
    }
    return(out)
  }

## Now read in the dat file to extract some important info
  Datfile = readLines(con=paste0(getwd(), "/calcur_paper-master1/petrale_ss3sim/om_test/om.dat", sep=""))
  Nage = as.numeric(substr(Datfile[grep("#_Nages", Datfile, fixed=TRUE)], 1, 2))
  # if length bin method is 2 (like in this case) do below:
  Nlength =  (as.numeric(substr(Datfile[grep("# maximum size ", Datfile, fixed=TRUE)], 1, 3)) - as.numeric(substr(Datfile[grep("# minimum size ", Datfile, fixed=TRUE)], 1, 2)))/as.numeric(substr(Datfile[grep("# binwidth for population", Datfile, fixed=TRUE)], 1, 2))+1

## now read in the Report.sso file
  rawrep <- read.table(file=paste0(getwd(), "/calcur_paper-master1/petrale_ss3sim/om_test/Report.sso"),col.names=1:200,fill=TRUE,quote="", colClasses="character",nrows=-1,comment.char="")

# Extract the general place where the AGE_LENGTH_KEY is written 
  LAA <- matchfun2(string1 ="MEAN_SIZE_TIMESERIES", adjust1=1, string2="mean_size_Jan_1_for_sex", adjust2=-2, header = TRUE)

# size from subseason 1 and morph 1
  LAA_sea1_morph1 <- LAA %>% filter(Morph==1, SubSeas==1)
 
# size from subseason 1 and morph 2
  LAA_sea1_morph2 <- LAA %>% filter(Morph==1, SubSeas==2)
  
# size from subseason 2 and morph 1
  LAA_sea2_morph1 <- LAA %>% filter(Morph==2, SubSeas==1)
  
# size from subseason 2 and morph 2
  LAA_sea2_morph2 <- LAA %>% filter(Morph==2, SubSeas==2)
  

# We then need to add some observation error on top of them to create some data to feed in the estimation model
  
  
  
  
  

df <- setup_scenarios_defaults()
df[, "si.years.2"] <- "seq(76, 100, by = 2)"
df <- rbind(df, df)
df[, "si.sds_obs.2"] <- c(0.1, 0.4)
df[, "bias_adjust"] <- FALSE
df[, "hess_always"] <- TRUE