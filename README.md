# ANNalog with docking
This project aims to bias the generation of the [ANNalog](https://github.com/DVNecromancer/ANNalog/tree/main) Seq2Seq model towards a given protein target by means of selecting a subset of generated molecules with good docking scores, as docked with Autodock Vina. This subset of molecules will be used for a subsequent round of generation with ANNalog, and the process will repeat for as many times as the user specifies.

# Installation
Using the [pixi package manager](https://github.com/prefix-dev/pixi), simply type the following to install the project environment:
```
pixi install
```
Ensure that this command is run in the same directory as the `pixi.lock` and `pixi.toml` files.

# Repository structure
- The structures tested and used in this work are in the `Structure` directory. The structure used in the MSc thesis is `refined_cnr2_human_5ZTY_inactive.pdb`. 
- The `redocking_active.ipynb` and `redocking_inactive.ipynb` files are the notebooks used for redocking cognate ligands for the refined structures from the GPCRdb of an inactive CB2 structure (PDB ID: 5ZTY) and active structure (PDB ID: 8GUR).
- The notebook used to develop and run the generation-docking workflow is `workflow_development.ipynb`.
- The `known_cb2_ligands` directory contains .`csv` files that contains information on known CB2 ligands from the ChEMBL database. The files used to build the known ligand set used for analysis in the MSc thesis are `CHEMBL253_CB2_Homo_sapiens_IC50_170726.csv` and `CHEMBL253_CB2_Homo_sapiens_Ki_140726.csv`, so as to focus on antagonists like the 5ZTY cognate ligand, AM10257, which is a potent CB2 antagonist.
- The results used in the MSc thesis are in the `results` directory. Inside, there is a `figures` directory which contains figures used in the work, as well as a `workflow` directory which contains the results files. There is a `protein` directory which contains the prepped protein files as well as `beam` and `sampling` directories containing five rounds of generation and docking results for both the random and docking score candidate selection methods.