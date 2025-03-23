

data <- data.frame(
    smina_scores_fname = list.files(
        "intermediate_data/",
        pattern = "smina_scores.tsv",
        recursive = TRUE)) |>
    dplyr::do({
        smina_scores_fname <- .$smina_scores_fname[1]
        smina_scores <- readr::read_tsv(
            smina_scores_fname,
            show_col_types = FALSE)
        smina_scores |>
            dplyr::mutate(
                fname = smina_scores_fname,
                .before = 1)
    })
