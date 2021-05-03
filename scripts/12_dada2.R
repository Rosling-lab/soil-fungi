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
                outfile <- file.path("process", seqrun, "regions", region_name, basename(trimmed_files))
                outdir <- dirname(outfile)
                if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
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
        tar_target(
            dereplicated,
            dada2::derepFastq(extracted, qualityType = "FastqQuality", verbose = TRUE),
            pattern = map(extracted),
            iteration = "list"
        ),
        tar_target(
            dadaobj,
            dada2::dada(
                dereplicated,
                err,
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
                BAND_SIZE = 64
            ),
            pattern = map(dereplicated, err),
            iteration = "list"
        ),
        tar_fst_tbl(
            dadamapping,
            tzara::dadamap(dereplicated, dadaobj, trimmed_files),
            pattern = map(dereplicated, dadaobj, trimmed_files),
            iteration = "list"
        )
    ),
    list(
        err = tar_target(
            err,
            dada2::learnErrors(
                dereplicated_full,
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
            ),
            pattern = map(dereplicated_full),
            iteration = "list"
        )
    )
)
