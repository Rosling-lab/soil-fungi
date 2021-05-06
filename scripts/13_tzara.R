# Recombine ASVs and generate consensus sequences using Tzara
# Brendan Furneaux
# April 2021

dadamap_exprs <- purrr::map(dada2_targets$dadamapping, ~parse(text = .$settings$name))

tzara_targets = list(
    combomap = tar_combine_raw(
        "combomap",
        dada2_targets$dadamapping,
        command = quote(join_dadamaps(!!!.x)),
        pattern = bquote(map(..(dadamap_exprs)), splice = TRUE)
    )
)
