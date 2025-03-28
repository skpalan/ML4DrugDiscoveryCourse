


data <- readr::read_tsv(
    file = "intermediate_data/approved_template_1HXW/rmsd.tsv",
    show_col_types = FALSE) |>
    dplyr::mutate(
        ligand_id = output_path |> stringr::str_replace("intermediate_data/approved_template_1HXW/", ""))

plot <- ggplot2::ggplot() +
    ggplot2::geom_point(
        data = data,
        mapping = ggplot2::aes(
            x = rmsd,
            y = affinity - affinity_best)) +
    ggplot2::facet_wrap(
        facets = dplyr::vars(ligand_id))

ggplot2::ggsave(
    filename = "product/score_vs_rmsd_20250328.pdf",
    width = 5,
    height = 5,
    useDingbats = FALSE)


