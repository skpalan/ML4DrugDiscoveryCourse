


model_checkpoint <- "${HOME}/turbo_bioinf595/opt/NeuralPLexer/neuralplexermodels_downstream_datasets_predictions/models/complex_structure_prediction.ckpt"

model_checkpoint <- "${HOME}/turbo/maom/opt/NeuralPLexer/neuralplexermodels_downstream_datasets_predictions/models/complex_structure_prediction.ckpt"





dock_neuralplexer <- function(
    input_ligand,
    input_receptor,
    output_path,
    model_checkpoint,
    input_template = NULL,
    application_path = "neuralplexer-inference --task=batched_structure_sampling ",
    n_samples = 16,
    chunk_size = 4,
    num_steps = 40,
    discard_sdf_coords = TRUE,
    overwrite = FALSE,
    verbose = FALSE) {


    if (!file.exists(input_receptor)) {
        stop("Input receptor path '", input_receptor, "' does not exist\n")
    }

    if (dir.exists(output_path)) {
        if (overwrite) {
            if (verbose) {
                cat(
                    "Outputpath '", output_path, "' exists and ovewrite is requested \n",
                    sep = "")
            }
        } else {
            if (verbose) {
                cat(
                    "Outputpath '", output_path, "' exists and overwrite not requested, ",
                    "skipping ...\n", sep = "")
            }
            return(1)
        } # overwrite
    } else {
        if (verbose) {
            cat ("Creating output path '", output_path, "'\n", sep = "")
        }
        dir.create(output_path, recursive = TRUE)
    }
    cmd <- paste0(
        application_path, " ",
        "--input-receptor '", input_receptor, "' ",
        "--input-ligand '", input_ligand, "' ",
        ifelse(
            !is.null(input_template),
            paste0("--use-template  --input-template ", input_template, " "),
            ""),
        "--out-path ", output_path, " ",
        "--model-checkpoint ", model_checkpoint, " ",
        "--n-samples ", n_samples, " ",
        "--chunk-size ", chunk_size, " ",
        "--num-steps ", num_steps, " ",
        "--cuda ",
        ifelse(discard_sdf_coords, "--discard-sdf-coords ", ""),
        "--sampler=langevin_simulated_annealing ",
        "--separate-pdb")
    if (verbose) {
        cat(cmd, "\n", sep = "")
    }
    system(cmd)
    return(1)
}

# APPROVED, should be 6 ligands
data <- readr::read_tsv(
    "data/1HXW/CHEMBL243_approved.tsv",
    show_col_types = FALSE) |>
    dplyr::distinct(Smiles, .keep_all = TRUE)

# dock in to holo receptor
data |>
    dplyr::rowwise() |>
    dplyr::do({
        ligand <- .
        substance_id <- ligand$`Molecule ChEMBL ID`[1]
        dock_neuralplexer(
            input_ligand = ligand$Smiles[1],
            input_receptor = "data/1HXW/protein.pdb",
            output_path = paste0("intermediate_data/approved_template_1HXW/", substance_id),
            model_checkpoint = model_checkpoint,
            input_template = "data/1HXW/protein.pdb",
            verbose = TRUE)
        data.frame()
    })

# dock into apo receptor
data |>
    dplyr::rowwise() |>
    dplyr::do({
        ligand <- .
        substance_id <- ligand$`Molecule ChEMBL ID`[1]
        dock_neuralplexer(
            input_ligand = ligand$Smiles[1],
            input_receptor = "data/1HHP/protein.pdb",
            output_path = paste0("intermediate_data/approved_template_1HHP/", substance_id),
            model_checkpoint = model_checkpoint,
            input_template = "data/1HHP/protein.pdb",            
            verbose = TRUE)
        data.frame()
    })

# dock into no template
data |>
    dplyr::rowwise() |>
    dplyr::do({
        ligand <- .
        substance_id <- ligand$`Molecule ChEMBL ID`[1]
        dock_neuralplexer(
            input_ligand = ligand$Smiles[1],
            input_receptor = "data/1HXW/protein.pdb",
            output_path = paste0("intermediate_data/approved_no_template/", substance_id),
            model_checkpoint = model_checkpoint,
            verbose = TRUE)
        data.frame()
    })



# inactives
data <- readr::read_tsv(
    "data/1HXW/CHEMBL243_inactive.tsv",
    show_col_types = FALSE) |>
    dplyr::distinct(Smiles, .keep_all = TRUE) |>
    dplyr::filter(
        `Standard Relation` == "'>'",
        `#RO5 Violations` %in% c("None", "1"))


# dock in to holo receptor
data |>
    dplyr::rowwise() |>
    dplyr::do({
        ligand <- .
        substance_id <- ligand$`Molecule ChEMBL ID`[1]
        dock_neuralplexer(
            input_ligand = ligand$Smiles[1],
            input_receptor = "data/1HXW/protein.pdb",
            output_path = paste0("intermediate_data/inactive_template_1HXW/", substance_id),
            model_checkpoint = model_checkpoint,
            input_template = "data/1HXW/protein.pdb",
            verbose = TRUE)
        data.frame()
    })

# dock into apo receptor
data |>
    dplyr::rowwise() |>
    dplyr::do({
        ligand <- .
        substance_id <- ligand$`Molecule ChEMBL ID`[1]
        dock_neuralplexer(
            input_ligand = ligand$Smiles[1],
            input_receptor = "data/1HHP/protein.pdb",
            output_path = paste0("intermediate_data/inactive_template_1HHP/", substance_id),
            model_checkpoint = model_checkpoint,
            input_template = "data/1HHP/protein.pdb",            
            verbose = TRUE)
        data.frame()
    })

# dock into no template
data |>
    dplyr::rowwise() |>
    dplyr::do({
        ligand <- .
        substance_id <- ligand$`Molecule ChEMBL ID`[1]
        dock_neuralplexer(
            input_ligand = ligand$Smiles[1],
            input_receptor = "data/1HXW/protein.pdb",
            output_path = paste0("intermediate_data/inactive_no_template/", substance_id),
            model_checkpoint = model_checkpoint,
            verbose = TRUE)
        data.frame()
    })
