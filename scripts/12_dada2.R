# Denoise homologous regions using DADA2
# Brendan Furneaux
# April 2021

dada2_targets <- c(
    tar_map(
        values = regions_meta,
        names = region_name,
        tar_file(
            extracted,
            {
                outdir <- file.path("process", seqrun, "regions", region_name)
                if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
                outfile <- file.path(outdir, basename(trimmed_files))
                tzara::extract_region(
                    seq = trimmed_files,
                    positions = positions,
                    region = start_region,
                    region2 = end_region,
                    outfile = outfile
                )
                outfile
            },
            pattern = map(trimmed_files, positions)
        ),
        tar_files(
            filtered,
            {
                outdir <- file.path("process", seqrun, "filtered", region_name)
                if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
                outfile <- file.path(outdir, basename(extracted))
                dada2::filterAndTrim(
                    extracted,
                    filt = outdir,
                    truncQ = 0,
                    maxLen = max_length,
                    minLen = min_length,
                    maxEE = 3,
                    multithread = local_cpus()
                )
                outfile
            }
        ),
        tar_target(
            dereplicated,
            dada2::derepFastq(
                filtered,
                qualityType = "FastqQuality",
                verbose = TRUE
            ),
            pattern = map(filtered),
            iteration = "list"
        ),
        tar_target(
            dadaobj,
            do_dada(
                dereplicated,
                err,
                region = region_name,
                err_region = "full",
                errorEstimationFunction = dada2::PacBioErrfun,
                multithread = local_cpus(),
                verbose = TRUE,
                DETECT_SINGLETONS = singletons,
                # use alignment scores as used internally in PacBio software
                MATCH = 1,
                MISMATCH = -2,
                GAP_PENALTY = -2,
                # reduce the homopolymer gap penalty
                HOMOPOLYMER_GAP_PENALTY = -1,
                # increase the band size
                BAND_SIZE = 64,
                pool = TRUE
            )
        ),
        tar_fst_tbl(
            dadamapping,
            tzara::dadamap(dereplicated, dadaobj, trimmed_files)
        )
    ),
    list(
        err = tar_target(
            err,
            dada2::dada(
                dereplicated_full,
                err = NULL,
                selfConsist = TRUE,
                errorEstimationFunction = dada2::PacBioErrfun,
                multithread = local_cpus(),
                verbose = TRUE,
                MAX_CONSIST = 50,
                DETECT_SINGLETONS = FALSE,
                # use alignment scores as used internally in PacBio software
                MATCH = 1,
                MISMATCH = -2,
                GAP_PENALTY = -2,
                # reduce the homopolymer gap penalty
                HOMOPOLYMER_GAP_PENALTY = -1,
                # increase the band size
                BAND_SIZE = 64
            )
        )
    )
)
