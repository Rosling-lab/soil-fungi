# read configuration files for use in the pipeline
# Brendan Furneaux
# April 2021

configdir <- "config"
regions_meta <- readr::read_csv(file.path(configdir, "default.regions.csv"))

seqruns <- list.dirs("process", recursive = FALSE, full.names = FALSE)

config_targets <- list{
    # get the path for the CM which is truncated at the LR5 primer site
    # (included in LSUx)
    cm_32S_trunc = tar_file(
        cm_32S_trunc,
        system.file(
            file.path("extdata", "fungi_32S_LR5.cm"),
            package = "LSUx"
        )
    )
}
