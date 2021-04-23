# Denoise homologous regions using DADA2
# Brendan Furneaux
# April 2021

regions_meta <- tibble::tribble(
    ~region,    ~start_region,     ~end_region,
    "ITS1",     "ITS1",            "ITS1",
    "5_8S",     "5_8S",            "5_8S",
    "ITS2",     "ITS2",            "ITS2",
    "LSU1",     "LSU1",            "LSU1",
    "D1",       "V2",              "V2",
    "LSU2",     "LSU2",            "LSU2",
    "D2",       "V3",              "V3",
    "LSU3",     "LSU3",            "LSU3",
    "D3",       "D3",              "D3",
    "LSU4",     "LSU4",            "LSU4"
)

dada2_targets <- c(
    tar_map(
        values = regions_meta,
        names = region,
        tar_fst_tbl(
            extracted,
            tzara::extract_region(
                seq = trimmed_files,
                positions = positions,
                region = start_region,
                region2 = end_region
            ),
            pattern = map(trimmed_files, positions)
        ),
        tar_target(
            derep,
            derep(extracted, qualityType = "FastqQuality", verbose = TRUE),
            pattern = map(extracted)
        ),
        tar_target(
            dada,
            dada2::dada(
                derep,
                err,
                errorEstimationFunction = dada2::PacBioErrfun(),
                multithread = local_cpus(),
                verbose = TRUE
            ),
            pattern = map(derep, err)
        ),
        tar_fst_tbl(
            dadamap,
            tzara::dadamap(derep, dada, trimmed_files),
            pattern = map(derep, dada, trimmed_files)
        )
    ),
    list(
        err = tar_target(
            err,
            dada2::learnErrors(
                derep_5_8S,
                errorEstimationFunction = dada2::PacBioErrfun(),
                multithread = local_cpus(),
                verbose = TRUE
            ),
            pattern = map(derep_5_8S)
        )
    )
)
