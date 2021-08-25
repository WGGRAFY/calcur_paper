#install
remotes::install_github("r4ss/r4ss", ref="development")
require(r4ss)
remotes::install_github("ss3sim/ss3sim",
                        ref = "main", build_vignettes = TRUE, dependencies = TRUE)
require(ss3sim)


#look for ss.exe
system.file("bin", package = "ss3sim")
