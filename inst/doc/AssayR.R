## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----usage, eval=FALSE---------------------------------------------------
#  # change directory to where your raw files are
#  setwd(path.to.raw.files)
#  # you can also do this in the RGui or RStudio menu.
#  
#  # convert raw files to mzML if they are not already
#  msconvert_all()
#  # this extract pos and neg and puts them in mzMLpos and mzMLneg
#  # directories
#  
#  # generate TICs/EICs/XICs from config file
#  path.to.tics <- run.config.tics(path.to.config, path.to.mzML)
#  # more about the config file below
#  
#  # analyse the peaks (interactively)
#  results <- run.config.peaks(path.to.config, path.to.tics)
#  # this is interactive, more below.
#  # you could write the result to csv
#  # write.csv(results, "results.csv")
#  
#  # rearrange and output the results, and output plots of the
#  # isotopes:
#  assay.plotter(results)
#  # check the current working directory for outputs
#  

