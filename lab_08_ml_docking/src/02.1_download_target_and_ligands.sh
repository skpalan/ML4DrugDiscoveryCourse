

# HIV Protease apo structure
# THE THREE-DIMENSIONAL STRUCTURE OF THE ASPARTYL PROTEASE FROM THE HIV-1 ISOLATE BRU
mkdir data/1HHP
cd data/1HHP
wget https://files.rcsb.org/download/1HHP.pdb
cd ../..

# HIV Protease holo structure
# HIV-1 PROTEASE DIMER COMPLEXED WITH A-84538
mkdir data/1HXW
cd data/1HXW
wget https://files.rcsb.org/download/1HXW.pdb
grep 'ATOM' 1HXW.pdb > protein.pdb
grep 'HETATM' 1HXW.pdb > ligand.pdb
#obabel -i pdb ligand.pdb -o sdf > ligand.sdf
cd ../..


