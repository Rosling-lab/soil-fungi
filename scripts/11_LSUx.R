# Cut raw reads into homologous regions using LSUx
# Brendan Furneaux
# April 2021

lsux_targets <- list(
    # find all the trimmed fastq files
    trimmed_files = tar_files(
            trimmed_files,
            list.files(
                path = "process",
                pattern = ".*trimmed.fastq.gz",
                full.names = TRUE
            )
        ),

    # get the path for the CM which is truncated at the LR5 primer site
    # (included in LSUx)
    cm_32S_trunc = tar_file(
        cm_32S_trunc,
        system.file(
            file.path("extdata", "fungi_32S_LR5.cm"),
            package = "LSUx"
        )
    ),

    # run LSUx on each input file
    positions = tar_fst_tbl(
        positions,
        LSUx::lsux(
            trimmed_files,
            # use the truncated cm which stops at LR5
            cm_32 = cm_32S_trunc,
            # include ITS1
            ITS1 = TRUE,
            cpu = local_cpus(),
            # allow 2 Gb ram (per process)
            mxsize = 2048
        ) %>%
            tibble::add_column(file = trimmed_files),
        pattern = map(trimmed_files),
        iteration = "list"
    )
)
