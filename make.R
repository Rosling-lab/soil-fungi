# make the targets pipeline
#
# Brendan Furneaux
# May 2021

# intended to be sourced in a fresh R session, bypasses callr so that a calling
# snakemake object may be available

library(targets)
tar_make(callr_function = NULL)

# error if tar_make encountered an error
stopifnot(all(is.na(tar_meta()$error)))
