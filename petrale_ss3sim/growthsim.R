#install packages
remotes::install_github("r4ss/r4ss", ref = "development")
remotes::install_github("ss3sim/ss3sim", ref = "development")
remotes::install_github("nmfs-fish-tools/SSMSE", ref="main")


require(SSMSE)
require(here)

# create paths
run_SSMSE_dir <- file.path("run_SSMSE-ex")
dir.create(run_SSMSE_dir)
petrale_mod_path <- here("om_SSMSE")


# develop_OMs will save a model called "cod_SR_BH_steep_1" in the out_dir
# specified
develop_OMs(OM_in_dir = petrale_mod_path, out_dir = run_SSMSE_dir, par_name = "L_at_Amax_Fem_GP_1",
            par_vals = 80, refit_OMs = FALSE, hess = FALSE)
# OM model for scenario 2
#cod_1_path <- file.path(run_SSMSE_dir, "om_L_at_Amax_Fem_GP_1")

#find dat file
datfile <- file.path(petrale_mod_path, "om.dat")

#creates a sample structure for the datfile
sample_struct <- create_sample_struct(dat = datfile, nyrs = 6) # note warning

#For years that are in the original model
sample_struct_hist_list
sample_struct_list <- list("base" = sample_struct)

#Add in values for NA quantities: base$CPUE$SE, base$lencomp, base$agecomp, remove mean size at age
sample_struct_list$base$meanbodywt$SE <- rep(0.1, 24)
sample_struct_list$base$CPUE$SE <- rep(0.01, 22)
sample_struct_list$base$lencomp$Sex[is.na(sample_struct_list$base$lencomp$Sex)] <-  rep(3,3)
sample_struct_list$base$lencomp$Part[is.na(sample_struct_list$base$lencomp$Part)] <-  rep(2,2)
sample_struct_list$base$lencomp$Nsamp <- rep(100,22)
sample_struct_list$base$lencomp <- sample_struct_list$base$lencomp[7:22,]

#Remove fishery age comp
sample_struct_list$base$agecomp <- sample_struct_list$base$agecomp[-c(1:4),]
sample_struct_list$base$agecomp$Nsamp <- rep(100,12)
sample_struct_list$base$agecomp$Lbin_lo <- sample_struct_list$base$agecomp$Lbin_hi <-  rep(-1,12)
sample_struct_list$base$agecomp$Sex <-  rep(3,12)
sample_struct_list$base$agecomp$Ageerr <-  rep(2,12)
sample_struct_list$base$MeanSize_at_Age_obs <- NULL
sample_struct_list$base$meanbodywt <- NULL


run_res_path <- file.path(run_SSMSE_dir, "results")
dir.create(run_res_path)
run <- run_SSMSE(scen_name_vec = c("base"),# name of the scenario
          out_dir_scen_vec = run_res_path, # directory in which to run the scenario
          iter_vec = c(5), # run with 5 iterations each
          OM_name_vec = NULL, # specify directories instead
          OM_in_dir_vec = c(petrale_mod_path), # OM files
          EM_in_dir_vec = c(petrale_mod_path),
          EM_name_vec = c("petrale"), # cod is included in package data
          MS_vec = c("EM"),       # The management strategy is specified in the EM
          use_SS_boot_vec = c(TRUE), # use the SS bootstrap module for sampling
          nyrs_vec = c(6),        # Years to project OM forward
          nyrs_assess_vec = c(3), # Years between assessments
          run_EM_last_yr = FALSE, # Run the EM in 106
          run_parallel = FALSE, # Run iterations in parallel
          sample_struct_list = sample_struct_list, # How to sample data for running the EM.
          seed = 12345) #Set a fixed integer seed that allows replication



######################
######################Try inbuilt model
######################
tv_mod_path <- system.file("extdata", "models", "growth_timevary", package = "SSMSE")
develop_OMs(OM_name = "growth_timevary", out_dir = run_SSMSE_dir, par_name = "L_at_Amax_Fem_GP_1",
            par_vals = 60, refit_OMs = FALSE, hess = FALSE)
