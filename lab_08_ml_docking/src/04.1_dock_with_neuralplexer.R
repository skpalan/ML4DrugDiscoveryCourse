

# APPROVED should be 6 ligands
data <- readr::read_tsv(
    "data/1HXW/CHEMBL243_approved.tsv",
    show_col_types = FALSE) |>
    dplyr::distinct(Smiles, .keep_all = TRUE)


use_template <- FALSE
data |>
    dplyr::rowwise() |>
    dplyr::do({
        ligand <- .
        substance_id <- ligand$`Molecule ChEMBL ID`[1]
        smiles <- ligand$Smiles[1]
        model_path <- "${HOME}/turbo_bioinf595/opt/NeuralPLexer/neuralplexermodels_downstream_datasets_predictions/models"
        data_path <- "data/1HXW"
        if (use_template) {
            output_path <- paste0("intermediate_data/1HXW_approved_template/", substance_id)
        } else {
            output_path <- paste0("intermediate_data/1HXW_approved_no_template/", substance_id)
        }
        if (!dir.exists(output_path)) {
            dir.create(output_path, recursive = TRUE)
            cat("docking ", substance_id, "\n", sep = "")
            cmd <- paste0(
                "neuralplexer-inference ",
                "--task=batched_structure_sampling ",
                "--input-receptor ", data_path, "/protein.pdb ",
                "--input-ligand '", smiles, "' ",
                ifelse(
                    use_template,
                    paste0("--use-template  --input-template ", data_path, "/protein.pdb "),
                    ""),
                "--out-path ", output_path, " ",
                "--model-checkpoint ", model_path, "/complex_structure_prediction.ckpt ",
                "--n-samples 16 ",
                "--chunk-size 4 ",
                "--num-steps=40 ",
                "--cuda ",
                "--discard-sdf-coords ",
                "--sampler=langevin_simulated_annealing ",
                "--separate-pdbs")
            cat(cmd, "\n", sep = "")
            system(cmd)
            data.frame()
        } else {
            cat(
                "Skipping docking ", substance_id, " ",
                "because output path ", output_path, " exists\n")
            data.frame()
        }
    })


# inactives
data <- readr::read_tsv(
    "data/1HXW/CHEMBL243_inactive.tsv",
    show_col_types = FALSE) |>
    dplyr::distinct(Smiles, .keep_all = TRUE) |>
    dplyr::filter(
        `Standard Relation` == "'>'",
        `#RO5 Violations` %in% c("None", "1"))


use_template <- TRUE
data |>
    dplyr::rowwise() |>
    dplyr::do({
        ligand <- .
        substance_id <- ligand$`Molecule ChEMBL ID`[1]
        smiles <- ligand$Smiles[1]
        model_path <- "${HOME}/turbo_bioinf595/opt/NeuralPLexer/neuralplexermodels_downstream_datasets_predictions/models"
        data_path <- "data/1HXW"
        if (use_template) {
            output_path <- paste0("intermediate_data/1HXW_inactive_template/", substance_id)
        } else {
            output_path <- paste0("intermediate_data/1HXW_inactive_no_template/", substance_id)
        }
        if (!dir.exists(output_path)) {
            dir.create(output_path, recursive = TRUE)
            cat("docking ", substance_id, "\n", sep = "")
            cmd <- paste0(
                "neuralplexer-inference ",
                "--task=batched_structure_sampling ",
                "--input-receptor ", data_path, "/protein.pdb ",
                "--input-ligand '", smiles, "' ",
                ifelse(
                    use_template,
                    paste0("--use-template  --input-template ", data_path, "/protein.pdb "),
                    ""),
                "--out-path ", output_path, " ",
                "--model-checkpoint ", model_path, "/complex_structure_prediction.ckpt ",
                "--n-samples 16 ",
                "--chunk-size 4 ",
                "--num-steps=40 ",
                "--cuda ",
                "--discard-sdf-coords ",
                "--sampler=langevin_simulated_annealing ",
                "--separate-pdbs")
            cat(cmd, "\n", sep = "")
            system(cmd)
            data.frame()
        } else {
            cat(
                "Skipping docking ", substance_id, " ",
                "because output path ", output_path, " exists\n")
            data.frame()
        }
    })
