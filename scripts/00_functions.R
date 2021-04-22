# functions for use in the pipeline
# Brendan Furneaux
# April 2021


#### utility functions ####
# are we running slurm?
is_slurm <- function() nchar(Sys.which("sbatch")) > 0
is_local <- function() !is_slurm()

# are we running snakemake?
is_snakemake <- function() !interactive() && exists("snakemake")

# how many cpus do we have on the local machine?
# if we're not running on the cluster, leave one cpu free.
local_cpus <- function() {
    if (is_snakemake()) {
        snakemake@threads
    } else if (is_slurm()) {
        out <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))
        assertthat::assert_that(assertthat::is.count(out))
        out
    } else {
        max(parallel::detectCores() - 1, 1)
    }
}
