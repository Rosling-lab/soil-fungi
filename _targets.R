# Master script for running targets pipeline
# normally this should be run from Snakemake using `R -e "targets::tar_make()"`
# Brendan Furneaux
# April 2021

library(targets)
library(tarchetypes)
library(magrittr)

# options for targets
tar_option_set(
    format = "qs", # default format qs: smaller and faster than default rds
    error = "continue" # if there is an error, keep trying other targets
)

# define necessary functions
source("scripts/00_functions.R")
# load configuration files
source("scripts/01_config.R")
# define the pipeline targets
source("scripts/11_LSUx.R")
source("scripts/12_dada2.R")
source("scripts/13_tzara.R")

c(
    config_targets,
    tar_map(
        values = list(seqrun = seqruns),
        lsux_targets,
        dada2_targets,
        tzara_targets
    )
)
