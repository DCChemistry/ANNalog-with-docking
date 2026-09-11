#!/usr/bin/env python3

import pandas as pd
import os
import shutil
from vina import Vina
from annalog_with_docking_functions import prepare_protein, vina_setup_receptor, prepare_ligands, docking, generate_and_deduplicate, cluster_and_select_candidates, get_generation_candidates_smi, _select_candidates_random

# workflow settings
workflow_rounds = 5 # run the workflow workflow_rounds times
gen_candidates_to_save = 5 # maximum number of candidates to save for the next generation-docking round
cluster_distance_threshold = 0.45 # distance threshold for clustering generated candidates

## preparation settings (protein+ligand)
input_smiles = "Cc1c(n(nc1C(=O)NC23CC4CC(C2)CC(C4)C3)CCCCCO)c5ccccc5" # PDB ligand ID: 9JU, bound to PDB ID 5ZTY
pH = 7.4 # physiological pH - specified for ligand and protein protonation
pH_deviation = 2 # deviation of +/- 2 log units for ligand protonation
enumerate_only_unassigned = True # only enumerate undefined stereocentres.
max_opt_iter = 1000 # when generating 3D molecules, this is the maximum number of optimisation steps for the 3D geometry optimisation

## docking settings
exhaustiveness = 24 # increased from the default of 8
n_poses = 1 # only save the best pose per ligand.
grid_buffer = 5 # buffer in Angstroms to add to the receptor grid dimensions around the co-crystal input structure cognate ligand.

## generation settings
num_of_mols_to_gen = 100
internal_multiplier = 2 # internally, generate 2 times the number of molecules to account for duplicates, which will be removed later

# input and output directories for the workflow
workflow_results_dir = os.path.join("results", "workflow") # working directory for the workflow results, which will contain subdirectories for each round of generation and docking.
input_protein_structure_path = os.path.join("Structures", "refined_cnr2_human_5ZTY_inactive.pdb") # protein structure to dock to (refined version of PDB ID: 5ZTY from GPCRdb - CB2 receptor)

# protein preparation
print("Preparing protein structure for docking...")
if(not os.path.exists("results/workflow/protein/_refined_cnr2_human_5ZTY_inactive.pdbqt")):
    protein_pdbqt_path = prepare_protein(input_protein_structure_path, workflow_results_dir, pH=pH)
else:
    protein_pdbqt_path = "results/workflow/protein/_refined_cnr2_human_5ZTY_inactive.pdbqt"

v=Vina(sf_name="vina", verbosity=0) # using the default vina scoring function, and setting the verbosity to 0 to silence Vina output
v=vina_setup_receptor(v, input_protein_structure_path, protein_pdbqt_path, grid_buffer=grid_buffer)


# generation method for ANNalog. Options are: "beam" (default), "BF-beam", and "sampling". See ANNalog documentation for more details.
gen_algo = "beam"
run_directory = os.path.join(workflow_results_dir, gen_algo) # naming the directory to save results to after the generation algorithm used


selected_candidates_filepath = None # initialise the selected candidates filepath to None, which will be updated after each round of generation and docking.
selection_method = "docking_score" # selection method for generation candidates

for round_num in range(1, workflow_rounds+1): # 1-indexing to make printing nicer
    print(f"Starting round {round_num} of {workflow_rounds}...")
    print("-"*60)

    generation_round_dir = os.path.join(run_directory, f"round_{round_num}")
    os.makedirs(generation_round_dir)

    if(round_num > 1): # for subsequent rounds, the input SMILES is the generation candidates from the previous round
        input_smiles = get_generation_candidates_smi(selected_candidates_filepath, generation_round_dir)

    print(f"Generating {num_of_mols_to_gen} molecules...")
    generated_smiles = generate_and_deduplicate(num_of_mols_to_gen, input_smiles, generation_round_dir, method=gen_algo, internal_multiplier=internal_multiplier)
    print(f"{len(generated_smiles)}/{num_of_mols_to_gen} molecules remain after deduplication.")
    print("-"*60)

    print("Preparing ligands for docking...")
    ligand_pdbqt_filepaths, enumerated_smiles_path = prepare_ligands(generated_smiles,
                                                                    generation_round_dir,
                                                                    pH=pH,
                                                                    pH_deviation=pH_deviation,
                                                                    enumerate_only_unassigned=enumerate_only_unassigned,
                                                                    max_iter=max_opt_iter)

    print("Docking...")
    docking_results_filepath = docking(ligand_pdbqt_filepaths, v, n_poses, exhaustiveness, enumerated_smiles_path, generation_round_dir)

    selected_candidates_filepath = cluster_and_select_candidates(docking_results_filepath, generation_round_dir, selection_method=selection_method, candidates_to_save=gen_candidates_to_save, clustering_threshold=cluster_distance_threshold)
print("-"*60+"\n"+"-"*60)
print("Workflow completed.")
# moving results to a docking score selection directory to make room for the random selection run
final_results_dir = os.path.join(run_directory, selection_method)
print(f"Moving results to {final_results_dir}")
os.makedirs(final_results_dir)
for round_num in range(1, workflow_rounds+1):
    generation_round_dir = os.path.join(run_directory, f"round_{round_num}")
    shutil.move(generation_round_dir, final_results_dir)

# loading in results from docking score round 1
round_1_gen_dock_dir = os.path.join(final_results_dir, "round_1") # the last generation-docking round directory will be used as the input for the next run
round_1_gen_dock_clustered_results_filepath = os.path.join(round_1_gen_dock_dir, "clustering", "clustered_results.csv")
results_df_with_clusters = pd.read_csv(round_1_gen_dock_clustered_results_filepath)
# selected candidates for the next generation-docking round randomly, where each candidate comes from a different cluster.
selected_candidates_df = _select_candidates_random(results_df_with_clusters, gen_candidates_to_save)
clustering_results_directory = os.path.join(run_directory, "round_1")
os.makedirs(clustering_results_directory)
selected_candidates_filepath = os.path.join(clustering_results_directory, "generation_candidates.csv")
selected_candidates_df.to_csv(selected_candidates_filepath, index=False)

selection_method = "random" # selection method for generation candidates

for round_num in range(2, workflow_rounds+1): # starting from round 2, since round 1 has already been completed with docking score selection
    print(f"Starting round {round_num} of {workflow_rounds}...")
    print("-"*60)

    generation_round_dir = os.path.join(run_directory, f"round_{round_num}")
    os.makedirs(generation_round_dir)

    if(round_num > 1): # for subsequent rounds, the input SMILES is the generation candidates from the previous round
        input_smiles = get_generation_candidates_smi(selected_candidates_filepath, generation_round_dir)

    print(f"Generating {num_of_mols_to_gen} molecules...")
    generated_smiles = generate_and_deduplicate(num_of_mols_to_gen, input_smiles, generation_round_dir, method=gen_algo, internal_multiplier=internal_multiplier)
    print(f"{len(generated_smiles)}/{num_of_mols_to_gen} molecules remain after deduplication.")
    print("-"*60)

    print("Preparing ligands for docking...")
    ligand_pdbqt_filepaths, enumerated_smiles_path = prepare_ligands(generated_smiles,
                                                                    generation_round_dir,
                                                                    pH=pH,
                                                                    pH_deviation=pH_deviation,
                                                                    enumerate_only_unassigned=enumerate_only_unassigned,
                                                                    max_iter=max_opt_iter)

    print("Docking...")
    docking_results_filepath = docking(ligand_pdbqt_filepaths, v, n_poses, exhaustiveness, enumerated_smiles_path, generation_round_dir)

    selected_candidates_filepath = cluster_and_select_candidates(docking_results_filepath, generation_round_dir, selection_method=selection_method, clustering_threshold=cluster_distance_threshold)
print("-"*60+"\n"+"-"*60)
print("Workflow completed.")
# moving results to a docking score selection directory to make room for the random selection run
final_results_dir = os.path.join(run_directory, selection_method)
print(f"Moving results to {final_results_dir}")
os.makedirs(final_results_dir)
for round_num in range(1, workflow_rounds+1):
    generation_round_dir = os.path.join(run_directory, f"round_{round_num}")
    shutil.move(generation_round_dir, final_results_dir)



### SAMPLING RUN
# redefining the input SMILES for the next run of the workflow
input_smiles = "Cc1c(n(nc1C(=O)NC23CC4CC(C2)CC(C4)C3)CCCCCO)c5ccccc5" # PDB ligand ID: 9JU, bound to PDB ID 5ZTY

# generation method for ANNalog. Options are: "beam" (default), "BF-beam", and "sampling". See ANNalog documentation for more details.
gen_algo = "sampling"
run_directory = os.path.join(workflow_results_dir, gen_algo) # naming the directory to save results to after the generation algorithm used


selected_candidates_filepath = None # initialise the selected candidates filepath to None, which will be updated after each round of generation and docking.
selection_method = "docking_score" # selection method for generation candidates

for round_num in range(1, workflow_rounds+1): # 1-indexing to make printing nicer
    print(f"Starting round {round_num} of {workflow_rounds}...")
    print("-"*60)

    generation_round_dir = os.path.join(run_directory, f"round_{round_num}")
    os.makedirs(generation_round_dir)

    if(round_num > 1): # for subsequent rounds, the input SMILES is the generation candidates from the previous round
        input_smiles = get_generation_candidates_smi(selected_candidates_filepath, generation_round_dir)

    print(f"Generating {num_of_mols_to_gen} molecules...")
    generated_smiles = generate_and_deduplicate(num_of_mols_to_gen, input_smiles, generation_round_dir, method=gen_algo, internal_multiplier=internal_multiplier)
    print(f"{len(generated_smiles)}/{num_of_mols_to_gen} molecules remain after deduplication.")
    print("-"*60)

    print("Preparing ligands for docking...")
    ligand_pdbqt_filepaths, enumerated_smiles_path = prepare_ligands(generated_smiles,
                                                                    generation_round_dir,
                                                                    pH=pH,
                                                                    pH_deviation=pH_deviation,
                                                                    enumerate_only_unassigned=enumerate_only_unassigned,
                                                                    max_iter=max_opt_iter)

    print("Docking...")
    docking_results_filepath = docking(ligand_pdbqt_filepaths, v, n_poses, exhaustiveness, enumerated_smiles_path, generation_round_dir)

    selected_candidates_filepath = cluster_and_select_candidates(docking_results_filepath, generation_round_dir, selection_method=selection_method, candidates_to_save=gen_candidates_to_save, clustering_threshold=cluster_distance_threshold)
print("-"*60+"\n"+"-"*60)
print("Workflow completed.")
# moving results to a docking score selection directory to make room for the random selection run
final_results_dir = os.path.join(run_directory, selection_method)
print(f"Moving results to {final_results_dir}")
os.makedirs(final_results_dir)
for round_num in range(1, workflow_rounds+1):
    generation_round_dir = os.path.join(run_directory, f"round_{round_num}")
    shutil.move(generation_round_dir, final_results_dir)


# loading in results from docking score round 1
round_1_gen_dock_dir = os.path.join(final_results_dir, "round_1") # the last generation-docking round directory will be used as the input for the next run
round_1_gen_dock_clustered_results_filepath = os.path.join(round_1_gen_dock_dir, "clustering", "clustered_results.csv")
results_df_with_clusters = pd.read_csv(round_1_gen_dock_clustered_results_filepath)
# selected candidates for the next generation-docking round randomly, where each candidate comes from a different cluster.
selected_candidates_df = _select_candidates_random(results_df_with_clusters, gen_candidates_to_save)
clustering_results_directory = os.path.join(run_directory, "round_1")
os.makedirs(clustering_results_directory)
selected_candidates_filepath = os.path.join(clustering_results_directory, "generation_candidates.csv")
selected_candidates_df.to_csv(selected_candidates_filepath, index=False)

selection_method = "random" # selection method for generation candidates

for round_num in range(2, workflow_rounds+1): # starting from round 2, since round 1 has already been completed with docking score selection
    print(f"Starting round {round_num} of {workflow_rounds}...")
    print("-"*60)

    generation_round_dir = os.path.join(run_directory, f"round_{round_num}")
    os.makedirs(generation_round_dir)

    if(round_num > 1): # for subsequent rounds, the input SMILES is the generation candidates from the previous round
        input_smiles = get_generation_candidates_smi(selected_candidates_filepath, generation_round_dir)

    print(f"Generating {num_of_mols_to_gen} molecules...")
    generated_smiles = generate_and_deduplicate(num_of_mols_to_gen, input_smiles, generation_round_dir, method=gen_algo, internal_multiplier=internal_multiplier)
    print(f"{len(generated_smiles)}/{num_of_mols_to_gen} molecules remain after deduplication.")
    print("-"*60)

    print("Preparing ligands for docking...")
    ligand_pdbqt_filepaths, enumerated_smiles_path = prepare_ligands(generated_smiles,
                                                                    generation_round_dir,
                                                                    pH=pH,
                                                                    pH_deviation=pH_deviation,
                                                                    enumerate_only_unassigned=enumerate_only_unassigned,
                                                                    max_iter=max_opt_iter)

    print("Docking...")
    docking_results_filepath = docking(ligand_pdbqt_filepaths, v, n_poses, exhaustiveness, enumerated_smiles_path, generation_round_dir)

    selected_candidates_filepath = cluster_and_select_candidates(docking_results_filepath, generation_round_dir, selection_method=selection_method, clustering_threshold=cluster_distance_threshold)
print("-"*60+"\n"+"-"*60)
print("Workflow completed.")
# moving results to a docking score selection directory to make room for the random selection run
final_results_dir = os.path.join(run_directory, selection_method)
print(f"Moving results to {final_results_dir}")
os.makedirs(final_results_dir)
for round_num in range(1, workflow_rounds+1):
    generation_round_dir = os.path.join(run_directory, f"round_{round_num}")
    shutil.move(generation_round_dir, final_results_dir)