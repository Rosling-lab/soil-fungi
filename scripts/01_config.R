# read configuration files for use in the pipeline
# Brendan Furneaux
# April 2021

configdir <- "config"
regions_meta <- readr::read_csv(file.path(configdir, "default.regions.csv"))
