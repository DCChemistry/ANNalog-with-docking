# ANNalog with docking
This project aims to bias the generation of the [ANNalog](https://github.com/DVNecromancer/ANNalog/tree/main) Seq2Seq model towards a given protein target by means of selecting a subset of generated molecules with good docking scores, as docked with Autodock Vina. This subset of molecules will be used for a subsequent round of generation with ANNalog, and the process will repeat for as many times as the user specifies.

# Installation
Using the [pixi package manager](https://github.com/prefix-dev/pixi), simply type the following to install the project environment to run the workflow:
```
pixi install
```
Ensure that this command is run in the same directory as the `pixi.lock` and `pixi.toml` files.

To install the separate conda environment needed for plotting, use the following command after creating a new conda environment (tested with Python 3.14.7):
```
conda install -c conda-forge rdkit umap-learn matplotlib seaborn pandas numpy scipy mdanalysis && pip install prolif starbars statannotations
```

# Repository structure
- The structures tested and used in this work are in the `Structure` directory. The structure used in the MSc thesis is `refined_cnr2_human_5ZTY_inactive.pdb`. 
- The `redocking_active.ipynb` and `redocking_inactive.ipynb` files are the notebooks used for redocking cognate ligands for the refined structures from the GPCRdb of an inactive CB2 structure (PDB ID: 5ZTY) and active structure (PDB ID: 8GUR).
- The notebook used to develop and run the generation-docking workflow is `workflow_development.ipynb`.
- The `known_cb2_ligands` directory contains `.csv` files that contains information on known CB2 ligands from the ChEMBL database. The files used to build the known ligand set used for analysis in the MSc thesis are `CHEMBL253_CB2_Homo_sapiens_IC50_170726.csv` and `CHEMBL253_CB2_Homo_sapiens_Ki_140726.csv`, so as to focus on antagonists like the 5ZTY cognate ligand, AM10257, which is a potent CB2 antagonist. As the filenames suggest, this ChEMBL data combines Ki and IC50 data.
- `redocking` contains the results from the `redocking_active.ipynb` and `redocking_inactive.ipynb` files. The inactive cognate ligand pose is used for some plotting in `plotting_results.ipynb`.
- `PLIP_results` contains .pse files obtained from [PLIP](https://plip-tool.biotec.tu-dresden.de/plip-web/plip/index) of the experimental cognate ligand of PDB 5ZTY, and the docked poses of the best scoring beam search-generated ligand and the best scoring multinomial sampling-generated ligand. These show the protein-ligand interactions as analysed by PLIP.
- `plotting_results.ipynb` contains the code used for making the majority of figures used in the MSc thesis associated with this project. This requires the conda environment mentioned above (as opposed to the pixi environment used for the workflow).