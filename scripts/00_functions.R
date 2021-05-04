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

# wrapper for dada2::dada
# always returns a list of dada objects,
# handles NULL inputs,
# doesn't re-fit if err is already a dada object from this region
do_dada <- function(derep, err, region, err_region, derepname = NA, ...) {
    # if our current region is the region that the error model was fit from,
    # then check if we have a valid dada object already, and don't refit.
    if (region == err_region) {
        if (methods::is(err, "dada")) {
            out <- list(err)
        } else if (is.list(err) &&
                   all(vapply(err, methods::is, TRUE, "dada") |
                       vapply(err, is.null, TRUE))) {
            return(err)
        }
    }
    # if we have only a single derep object, then wrap it in a list
    # and optionally name it.
    if (methods::is(derep, "derep")) {
        derep <- list(derep)
        if (!is.na(derepname)) names(derep) <- derepname
    }
    # find and remove any NULL derep objects.
    oldnames <- names(derep)
    nullderep <- vapply(derep, is.null, TRUE)
    derep <- derep[!nullderep]
    # run dada2
    out <- dada2::dada(derep, err, ...)
    # put back NULLs
    out <- out[oldnames]
    names(out) <- oldnames
    out
}
