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


make_data_in <- function(tv_k = "0",
                         tv_Linf = "0",
                         tv_mat = "0", sigma_R = "0",
                         f_scen = "0", nsamp_ages = "0",
                         sel = "0", cvL= "0"){
  df <- setup_scenarios_defaults()
  
  #f value years
  df[,"cf.years.1"] <- "20:100"
  #F values - should match dimensions of above
  df[,"cf.fvals.1"] <- switch(f_scen, 
                              "0"="rep(0.1, 81)",
                              #regime of F matching Linf
                              "1"= "c(rep(0,31),rep(0.4,13),rep(0,9)), rep(0.4,16),rep(0,12))", #monotonic increase in F
                              "2" = "seq(0,0.4,length.out=81)")
  #years to sample for index, length, age samples
  df[, "si.years.2"] <- 
    df[,"sl.years.1"] <-
    df[,"sl.years.2"] <-
    df[,"sl.years.2"] <- 
    df[,"sa.years.2"] <-
    df[,"sa.years.1"] <- "50:100"
  
  #0 = default - scenario is 50 age samples per year and 100 survey age samples per year
  df[,"sa.Nsamp.1"] <- switch(nsamp_ages,"0" = 50,
                              "1" = 10)
  #1 = 10 age samples per year from fishery and 30 from survey
  df[,"sa.Nsamp.2"] <- switch(nsamp_ages,"0" = 100,
                              "1" = 30)
  
  # Now adding some conditional age at length information
  df[,"sc.years.1"] <- df[, "sa.years.1"]   # same sampling year as the age comps 1
  df[,"sc.years.2"] <- df[, "sa.years.2"]   # same sampling year as the age comps 2
  df[,"sc.Nsamp_lengths.1"] <- df[,"sl.Nsamp.1"]    # same sample size as the length comps 1
  df[,"sc.Nsamp_lengths.2"] <- df[,"sl.Nsamp.2"]    # same sample size as the length comps 2
  df[,"sc.Nsamp_ages.1"] <- df[,"sa.Nsamp.1"]    # same sample size as the age comps 1
  df[,"sc.Nsamp_ages.2"] <- df[,"sa.Nsamp.2"]    # same sample size as the age comps 2
    
  df <- rbind(df, df)
  df[, "si.sds_obs.2"] <- c(0.1, 0.4)
  
  df[,"ct.Size_DblN_peak_WinterN(1)"] <- switch(sel,
                                         "0" = NULL,
                                         #time varying sel
                                         "1" = "Size_DblN_peak_WinterN(1) = c(rep(0,50),input[[2]]$V1)")
  
  #Use change_tv to make VB L inf a decadal trend
  df[,"ct.L_at_Amax_Fem_GP_1"] <- switch(tv_Linf,
                            "0" = NULL,
                             #regime for L inf
                             "1" = "L_at_Amax_Fem_GP_1 = c(rep(0,50),rep(2.8924,13),
                            rep(3.0652,10),rep(1.1188,16),rep(-0.5686,11))",
                             #monotonic increase l2
                             "2" = "L_at_Amax_Fem_GP_1 = seq(0,4,length.out=100)")
  
  df[,"ct.VonBert_K_Fem_GP_1"] <- switch(tv_k,
                                         "0" = NULL, #Regime for k
  "1" = "VonBert_K_Fem_GP_1 = rep(c(rep(0,25),rep(0.1,25)),2)",
  #monotonic increase k
  "2" = "VonBert_K_Fem_GP_1 = seq(0,0.1,length.out=100)")

  #This doesn't work I think because of the %
  df[,"ct.Mat50%_Fem_GP_1"] <- switch(tv_mat,
                                      "0"=NULL,
                                      #Declining length at maturity
                                      "1"="Mat50%_Fem_GP_1=seq(0,-5,length.out=100)")
  
  df[,"co.par_name"] <- "SR_sigmaR"
  df[,"co.par_int"] <- switch(sigma_R, 
                              "0"=0.1,
                              "1"=0.4) 
  df[,"co.par_name"] <- "CV_old_Fem_GP_1"
  df[,"co.par_int"] <- switch(cvL,
                              "0" = .035,
                              "1" = .085)
  
  
  scname <- c(paste0("D1-L", tv_Linf, 
                     "-A",nsamp_ages,
                     "-K",tv_k,
                     "-M",tv_mat,
                     "-R",sigma_R,
                     "-F", f_scen,
                     "-C",cvL,"-pet"),
              paste0("D2-L", tv_Linf, 
                     "-A",nsamp_ages,
                     "-K",tv_k,
                     "-M",tv_mat,
                     "-R",sigma_R,
                     "-F", f_scen, 
                     "-C",cvL,"-pet"))
  df[,"scenarios"] <- scname
  df[, "bias_adjust"] <- FALSE
  df[, "hess_always"] <- FALSE
  return(df)
}

turn_ss_out_to_data <- function(em_dat, flt = 1){
  len_matrix <- em_dat$lencomp 
  #Multiply length comp by Nsamp and split into fleets
  len_mult_flt <- len_matrix %>%
    mutate(across(.cols = l20:l152,~.x*Nsamp)) %>%
    split(len_matrix$FltSvy)
  
  #Remove marginal age comps
  caal_matrix <- em_dat$agecomp %>% 
    filter(Lbin_lo != (-1))
  
  #Get CAAL by age
  caal_matrix_flt <- split(caal_matrix, caal_matrix$FltSvy)
  out_tibble <- NULL
   for(i in len_mult[[flt]]$Yr){
     #there are two conflicting rows being written, bug in ss3sim? the first one seems right
     lenrow <- filter(len_mult[[flt]], Yr==i)[1,]
     caalrows <- filter(caal_matrix_flt[[flt]], Yr==i)

     #len bin column regex
    regmatches<- grep("^[l0-9]",names(lenrow))
     
     #fill missing bins with 0 

       #get all bins
       lbins <- gsub("l","",names(lenrow)[regmatches])
       full_bins <- expand_grid(lbins = as.integer(lbins))
       full_matches<- left_join(full_bins, caalrows, by=c("lbins"="Lbin_lo"), keep = TRUE) %>%
         select(starts_with("a", ignore.case=FALSE)) %>%
         mutate(across(everything(),replace_na,0)) 
       full_matches$lencomp <- t(lenrow[regmatches])
       
       #This matrix is p age length with age as columns and length bins as rows
      p_age_length <- full_matches %>% mutate(across(starts_with("a")),.*lencomp)
      p_age_length$length <- full_bins 
      p_age_length$Yr <- i
      p_age_length <- p_age_length[,c("Yr","length","lencomp",paste0("a",1:25))]
      out_tibble <- bind_rows(out_tibble, p_age_length)
   }
     
   
  return(out_tibble)
}
