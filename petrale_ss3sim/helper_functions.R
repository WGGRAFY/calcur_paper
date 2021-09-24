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


make_data_in <- function(tv_params_arg = 1){
  df <- setup_scenarios_defaults()
  
  #f value years
  df[,"cf.years.1"] <- "20:100"
  #F values - should match dimensions of above
  df[,"cf.fvals.1"] <- "rep(0.1, 81)"
  #years to sample for index, length, age samples
  df[, "si.years.2"] <- 
    df[,"sl.years.1"] <-
    df[,"sl.years.2"] <-
    df[,"sa.years.1"] <-
    df[,"sl.years.2"] <-
    df[,"sa.years.2"] <- "50:100"
  
  df <- rbind(df, df)
  df[, "si.sds_obs.2"] <- c(0.1, 0.4)
  
  #Use change_tv to make VB L inf a decadal trend
  df[,"tv_params"] <- switch(tv_params_arg,
                             #regime for L inf
                             "1" = "list(L_at_Amax_Fem_GP_1 = c(rep(0,50),input[[2]]$V1))",
                             #monotonic increase l2
                             "2" = "list(L_at_Amax_Fem_GP_1 = seq(0,4,length.out=100)",
                             #Regime for k
                             "3" = "list(VonBert_K_Fem_GP_1 = rep(c(rep(0,25),rep(0.1,25)),2)",
                             #monotonic increase k
                             "4" = "list(VonBert_K_Fem_GP_1 = seq(0,0.1,length.out=100)")
  scname <- c(paste0("D1-E", tv_params_arg, "-F0-pet"),
              paste0("D2-E", tv_params_arg, "-F0-pet"))
  df[,"scenarios"] <- scname
  df[, "bias_adjust"] <- FALSE
  df[, "hess_always"] <- FALSE
  return(df)
}

