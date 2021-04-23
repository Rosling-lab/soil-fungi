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

#### dada2 functions ####

# call dada2::derepFastq() on an R object containing sequences instead of a file.
derep <- function(reads, n, verbose, qualityType) {
    UseMethod("derep")
}

derep.ShortReadQ <- function(reads, n = 1e+06, verbose = FALSE, qualityType = "Auto") {
    if (length(reads) == 0) return(NULL)
    fname <- tempfile("reads", fileext = ".fastq.gz")
    on.exit(unlink(fname))
    ShortRead::writeFastq(reads, fname, compress = TRUE, qualityType = qualityType)
    dada2::derepFastq(fname, n = n, verbose = verbose, qualityType = qualityType)
}

derep.QualityScaledXStringSet <- function(reads, n = 1e+06, verbose = FALSE, qualityType = "Auto") {
    if (length(reads) == 0) return(NULL)
    fname <- tempfile("reads", fileext = ".fastq.gz")
    on.exit(unlink(fname))
    Biostrings::writeQualityScaledXStringSet(reads, fname, compress = TRUE)
    dada2::derepFastq(fname, n = n, verbose = verbose, qualityType = qualityType)
}
