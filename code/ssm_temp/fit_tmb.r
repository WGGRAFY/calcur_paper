fit_tmb = function(model, n.newton=3, do.sdrep=TRUE, do.check=FALSE, save.sdrep=FALSE)
{
  model$opt <- stats::nlminb(model$par, model$fn, model$gr, control = list(iter.max = 1000, eval.max = 1000))
  if(n.newton){ # Take a few extra newton steps
    # print("is n.newton")
    tryCatch(for(i in 1:n.newton) { 
      g <- as.numeric(model$gr(model$opt$par))
      h <- stats::optimHess(model$opt$par, model$fn, model$gr)
      model$opt$par <- model$opt$par - solve(h, g)
      model$opt$objective <- model$fn(model$opt$par)
    }, error = function(e) {model$err <<- conditionMessage(e)}) # still want fit_tmb to return model if newton steps error out
  }
  #assigning model$err already does the below if statement.
  #if(exists("err", inherits = FALSE)){
  #  model$err <- err # store error message to print out in fit_wham
  #  rm("err")
  #}  
  #model$env$parList() gives error when there are no random effects
  is.re = length(model$env$random)>0
  fe = model$opt$par
  if(is.re) model$env$last.par.best[-c(model$env$random)] = fe
  else  model$env$last.par.best = fe
  #fe = model$env$last.par.best
  #if(is.re) fe = fe[-c(model$env$random)]
  
  Gr = model$gr(fe)
  if(do.check){
    if(any(Gr > 0.01)){
      df <- data.frame(param = names(fe),
                       MLE = fe,
                       gr.at.MLE = Gr)
      ind.hi <- which(Gr > 0.01)
      model$badpar <- df[ind.hi,]
      warning(paste("","Some parameter(s) have high gradients at the MLE:","",
        paste(capture.output(print(model$badpar)), collapse = "\n"), sep="\n"))
    } else {
      test <- TMBhelper::check_estimability(model)
      if(length(test$WhichBad) > 0){
        bad.par <- as.character(test$BadParams$Param[test$BadParams$Param_check=='Bad'])
        bad.par.grep <- grep(bad.par, test$BadParams$Param)
        model$badpar <- test$BadParams[bad.par.grep,]
        warning(paste("","Some fixed effect parameter(s) are not identifiable.",
          "Consider 1) removing them from the model by fixing input$par and input$map = NA, or",
          "2) changing your model configuration.","",
          paste(capture.output(print(test$BadParams[bad.par.grep,])), collapse = "\n"), sep="\n"))    
      }
    }
  }

  model$date = Sys.time()
  model$dir = getwd()
  model$rep <- model$report()
  # model$TMB_version = packageVersion("TMB")
  #ver <- sessioninfo::package_info() %>% as.data.frame %>% dplyr::filter(package=="TMB") %>% dplyr::select(loadedversion, source) %>% unname
  #model$TMB_version <- paste0(ver, collapse=" / ")
  model$parList = model$env$parList(x = fe)
  model$final_gradient = Gr

  # if(do.sdrep & !exists("err")) # only do sdrep if no error
  if(do.sdrep) # only do sdrep if no error
  {
    model$sdrep <- try(TMB::sdreport(model))
    model$is_sdrep = !is.character(model$sdrep)
    if(model$is_sdrep) model$na_sdrep = any(is.na(summary(model$sdrep,"fixed")[,2])) else model$na_sdrep = NA
    if(!save.sdrep) model$sdrep <- summary(model$sdrep) # only save summary to reduce model object size
  } else {
    model$is_sdrep = FALSE
    model$na_sdrep = NA
  }

  return(model)
}
