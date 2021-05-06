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
        get("snakemake")@threads
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
do_dada <- function(derep, err, region, err_region, derepname = NULL, ...) {
    # if our current region is the region that the error model was fit from,
    # then check if we have a valid dada object already, and don't refit.
    if (region == err_region) {
        if (methods::is(err, "dada")) {
            return(list(err))
        } else if (is.list(err) &&
                   all(vapply(err, methods::is, TRUE, "dada") |
                       vapply(err, is.null, TRUE))) {
            return(err)
        }
    }
    # if we have only a single derep object, then wrap it in a list
    if (methods::is(derep, "derep")) {
        derep <- list(derep)
    }
    # name the dereps if necessary
    if (!is.null(derepname) && is.null(names(derep))) names(derep) <- derepname
    # make an all-NULL output
    out <- vector("list", length(derep))
    names(out) <- names(derep)
    # run dada2 on non-null entries only
    nullderep <- vapply(derep, is.null, TRUE)
    dadaout <- dada2::dada(derep[!nullderep], err, ...)
    if (methods::is(dadaout, "dada")) dadaout <- list(dadaout)
    out[!nullderep] <- dadaout
    out
}

join_dadamaps <- function(..., key_region = "ITS2") {
    key_region_q <- parse(text = key_region)
    dplyr::bind_rows(...) %>%
        dplyr::select(name, seq_id, region, derep_seq, dada_seq) %>%
        dplyr::filter(region != key_region | !is.na(dada_seq)) %>%
        dplyr::mutate(dada_seq = dplyr::coalesce(dada_seq, derep_seq)) %>%
        dplyr::select(-derep_seq) %>%
        tidyr::pivot_wider(names_from = region, values_from = dada_seq) %>%
        dplyr::filter(!is.na(!!key_region_q)) %>%
        dplyr::select(-seq_id) %>%
        dplyr::group_by_all() %>%
        dplyr::summarise(n_read = dplyr::n())

}
