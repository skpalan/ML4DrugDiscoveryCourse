
# Machine Learning Docking

The goal of this lab is to use machine learning docking to dock
known active and inactive compounds into a target. Here we will use
NeuralPlexer model,


    Qiao Z, Nie W, Vahdat A, Miller III TF, Anandkumar A.
    State-specific protein-ligand complex structure prediction with a multi-scale deep
    generative model. Nature Machine Intelligence, 2024. https://doi.org/10.1038/s42256-024-00792-z.
    https://github.com/zrqiao/NeuralPLexer


We will revisit the HIV-1 Protease as the target and dock approved drugs
and a selection of inactive ligands curated from ChEMBL. NeuralPlexer
does "co-folding" in that it predicts the protein structure along with
the docked ligand conformation.

To score bound ligand conformations, we will use gnina, a GNN based ligand docking scoring model

    GNINA 1.0: Molecular docking with deep learning
    A McNutt, P Francoeur, R Aggarwal, T Masuda, R Meli, M Ragoza, J Sunseri, DR Koes.
    J. Cheminformatics, 2021 https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00522-2



## Lab steps

1) Install NeuralPlexer and dependencies
      This is a little tricky so follow the instructions carefully!
      src/01_install_packages.sh

2) Download and install targets and ligands and the models
      02.1_download_target_and_ligands.sh

   Using the instructions from Lab 1 download approved ligands
   and inactive ligands from ChEMBL for target CHEMBL243.
   Pre-downloaded versions are here, which you can use
      data/1HXW/CHEMBL243_approved.tsv
      data/1HXW/CHEMBL243_inactives.tsv

3) Test that NeuralPlexer and smina are working on your system
   03.1_test_neuralplex_docking.sh
   03.2_test_scoring.sh

4) Dock the approve and inactive compounds into HIV-1 protease
   04.1_dock_with_neuralplexer.R
     Consider docking
       * with 1HXW where only the sequence is used from the pdb
         (exclude the --use-template --input-template <path>) flags. 
       * with 1HXW that is used as a template
       * with 1HHP that is used as a template
     => so 8 different docking runs:
           (approved, inactive) * (use-template, no template) * (1HXW, 1HHP)

5) Score docked poses using the SMINA molecular mechanics focefield
     04.2_score_predictions.R

6) Analyze the resulting scores to see if any of template/no-template
     NeuralPlexer docking runs are able to discriminate the actives
     from inactives. Plot the results and use a T-test to test
     for significance.

Note: that some of the results have been computed as an example
but the rest will need to be computed.


## Writeup

   Summarize the docking that you did and test if the NeuralPlexer
   docking poses combined with SMINA docking scores are able to
   explain the experimentally characterized activities.
