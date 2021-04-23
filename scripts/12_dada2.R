# Denoise homologous regions using DADA2
# Brendan Furneaux
# April 2021

regions_meta <- tibble::tribble(
    ~region_name,    ~start_region,     ~end_region,
    "ITS1",     "ITS1",            "ITS1",
    "5_8S",     "5_8S",            "5_8S",
    "ITS2",     "ITS2",            "ITS2",
    "LSU1",     "LSU1",            "LSU1",
    "D1",       "V2",              "V2",
    "LSU2",     "LSU2",            "LSU2",
    "D2",       "V3",              "V3",
    "LSU3",     "LSU3",            "LSU3",
    "D3",       "V4",              "V4",
    "LSU4",     "LSU4",            "LSU4"
)

dada2_targets <- c(
    tar_map(
        values = regions_meta,
        names = region_name,
        tar_target(
            extracted,
            tzara::extract_region(
                seq = trimmed_files,
                positions = positions,
                region = start_region,
                region2 = end_region
            ),
            pattern = map(trimmed_files, positions),
            iteration = "list"
        ),
        tar_target(
            dereplicated,
            derep(extracted, qualityType = "FastqQuality", verbose = TRUE),
            pattern = map(extracted),
            iteration = "list"
        ),
        tar_target(
            dadaobj,
            dada2::dada(
                dereplicated,
                err,
                errorEstimationFunction = dada2::PacBioErrfun(),
                multithread = local_cpus(),
                verbose = TRUE,
                DETECT_SINGLETONS = TRUE,
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
                dereplicated_5_8S,
                errorEstimationFunction = dada2::PacBioErrfun(),
                multithread = local_cpus(),
                verbose = TRUE,
                DETECT_SINGLETONS = TRUE,
                # use alignment scores as used internally in PacBio software
                MATCH = 1,
                MISMATCH = -2,
                GAP_PENALTY = -2,
                # reduce the homopolymer gap penalty
                HOMOPOLYMER_GAP_PENALTY = -1,
                # increase the band size
                BAND_SIZE = 64
            ),
            pattern = map(dereplicated_5_8S),
            iteration = "list"
        )
    )
)
