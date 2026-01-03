#pip install pymican
#--> MIT License https://github.com/ShintaroMinami/mican

# cif conversion related modules
from Bio.PDB import *
import numpy as np


# import superimposer related modules
from pymican import mican
import os
from collections import namedtuple
import multiprocessing as mp
from functools import partial
from multiprocessing import Pool, Manager
import os

import json



def cifToPdbBiopy(cif_file, pdb_file):
    """Converts a CIF file to PDB format using BioPython's parser
    Args:
        cif_file (str): path to input CIF file
        pdb_file (str): path where output PDB file will be saved
    Returns:
        None: saves converted file to specified path
    """
    parser = MMCIFParser()
    structure = parser.get_structure('structure', cif_file)
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_file)


def superimposeMican(pdb1, pdb2, export=False, export_path=None):
    """Superimposes two protein structures using MICAN
    Args:
        pdb1 (str): path to first PDB file
        pdb2 (str): path to second PDB file
        export (bool, optional): save superimposed structure. Defaults to False
    Returns:
        rmsd (float): Root Mean Square Deviation
        tmscore (float): TM-score of the alignment
    """
    # create object
    m = mican()
    #print(pdb1)
    
    if export:
        # align with export
        if export_path is not None:
            aln = m.align(pdb1=pdb1, pdb2=pdb2, options=f"-f -o {export_path}")
        else:
            aln = m.align(pdb1=pdb1, pdb2=pdb2, options='-f -o super.pdb')
    else:
        #align only
        aln = m.align(pdb1=pdb1, pdb2=pdb2, options='-s') 
    
    # get TM-score, RMSD, etc.
    tmscore = aln.TMscore
    rmsd  = aln.rmsd
    
    # get TM-score, RMSD, etc.
    return aln, rmsd, tmscore


def batchSuperimposeMican(folder, export=False):
    """Performs batch superimposition of protein structures
    Processes pairs of structures from predictions and groundtruth folders
    Args:
        folder (str): path containing predictions and groundtruth folders
        export (bool, optional): save superimposed structures. Defaults to False
    Returns:
        results (list): list of tuples containing:
            - id (str): structure identifier
            - rmsd (float): Root Mean Square Deviation
            - tm_score (float): TM-score of alignment
    """
    # Dictionary to hold the pairs of files
    file_pairs = {}
    
    # Check and get files from predictions subfolder
    pred_folder = os.path.join(folder, "predictions")
    gt_folder = os.path.join(folder, "groundtruth")
    
    # Check if both directories exist
    if not os.path.exists(pred_folder) or not os.path.exists(gt_folder):
        raise ValueError(f"One or both directories missing:\nPredictions: {pred_folder}\nGround truth: {gt_folder}")
    
    # Get files from predictions folder
    pred_files = os.listdir(pred_folder)
    for file in pred_files:
        if file.endswith('.pdb'):
            base_name = file[:-4]  # Remove the .pdb extension
            file_pairs[base_name] = {'pred': file}
    
    # Get files from ground truth folder and match with predictions
    gt_files = os.listdir(gt_folder)
    for file in gt_files:
        if file.endswith('.pdb'):
            base_name = file[:-4]  # Remove the .pdb extension
            if base_name.endswith('_gt'):
                id = base_name[:-3]  # Remove the _gt suffix
                if id in file_pairs:
                    file_pairs[id]['gt'] = file
    
    # Process each pair of files
    results = []
    missing_files = []
    
    for id, files in file_pairs.items():
        # Check if both prediction and ground truth files exist
        if 'pred' not in files:
            missing_files.append(f"Missing prediction file for {id}")
            continue
        if 'gt' not in files:
            missing_files.append(f"Missing ground truth file for {id}")
            continue
            
        pred_path = os.path.join(pred_folder, files['pred'])
        gt_path = os.path.join(gt_folder, files['gt'])
        
        # Additional check for file existence
        if not os.path.isfile(pred_path):
            missing_files.append(f"Prediction file not found: {pred_path}")
            continue
        if not os.path.isfile(gt_path):
            missing_files.append(f"Ground truth file not found: {gt_path}")
            continue
        
        # Calculate scores
        _, rmsd, tm_score = superimposeMican(pred_path, gt_path, export)
        
        # Store results
        results.append((id, rmsd, tm_score))
    
    # Report any missing files
    if missing_files:
        print("Warning: The following files were missing:")
        for msg in missing_files:
            print(f"- {msg}")
    
    return results

######### Parallel version with Mican #########
def process_pair(args):
    """Helper function to process a single pair of structures
    Args:
        args (tuple): contains (id, files, pred_folder, gt_folder, export)
    Returns:
        tuple: (id, result_tuple) or (id, error_message)
    """
    id, files, pred_folder, gt_folder, export, method = args
    
    if 'pred' not in files:
        return (id, f"Missing prediction file for {id}")
    if 'gt' not in files:
        return (id, f"Missing ground truth file for {id}")
        
    pred_path = os.path.join(pred_folder, files['pred'])
    gt_path = os.path.join(gt_folder, files['gt'])
    
    if not os.path.isfile(pred_path):
        return (id, f"Prediction file not found: {pred_path}")
    if not os.path.isfile(gt_path):
        return (id, f"Ground truth file not found: {gt_path}")
    
    try:
        if method == 'mican':
            _, rmsd, tm_score = superimposeMican(pred_path, gt_path, export)
            return (id, (rmsd, tm_score))

        if method == 'svd':
            _, _, rmsd, tm_score = SuperimposeSVD(pred_path, gt_path, export)
            return (id, (rmsd, tm_score))

    except Exception as e:
        return (id, f"Error processing {id}: {str(e)}")


def batchSuperimposeParallel(folder, export=False, n_processes=None, method='mican'):
    """Performs parallel batch superimposition of protein structures
    Processes pairs of structures from predictions and groundtruth folders

    Args:
        folder (str): path containing predictions and groundtruth folders
        export (bool, optional): save superimposed structures. Defaults to False
        n_processes (int, optional): number of processes to use. Defaults to CPU count
        method (str, optional): superimposition method. 'mican' or 'svd'. Defaults to 'mican'

    Returns:
        results (list): list of tuples containing:
            - id (str): structure identifier
            - rmsd (float): Root Mean Square Deviation
            - tm_score (float): TM-score of alignment
    """
    # Dictionary to hold the pairs of files
    file_pairs = {}
    
    # Check and get files from predictions subfolder
    pred_folder = os.path.join(folder, "predictions")
    gt_folder = os.path.join(folder, "groundtruth")
    
    # Check if both directories exist
    if not os.path.exists(pred_folder) or not os.path.exists(gt_folder):
        raise ValueError(f"One or both directories missing:\nPredictions: {pred_folder}\nGround truth: {gt_folder}")
    
    # Get files from predictions folder
    pred_files = os.listdir(pred_folder)
    for file in pred_files:
        if file.endswith('.pdb'):
            base_name = file[:-4]  # Remove the .pdb extension
            file_pairs[base_name] = {'pred': file}
    
    # Get files from ground truth folder and match with predictions
    gt_files = os.listdir(gt_folder)
    for file in gt_files:
        if file.endswith('.pdb'):
            base_name = file[:-4]  # Remove the .pdb extension
            if base_name.endswith('_gt'):
                id = base_name[:-3]  # Remove the _gt suffix
                if id in file_pairs:
                    file_pairs[id]['gt'] = file
    
    # Prepare arguments for parallel processing
    process_args = [
        (id, files, pred_folder, gt_folder, export, method)
        for id, files in file_pairs.items()
    ]
    
    # Initialize results and missing files lists
    results = []
    missing_files = []
    
    # Process files in parallel
    with Pool(processes=n_processes) as pool:
        processed_results = pool.map(process_pair, process_args)
        
        # Sort results and errors
        for id, result in processed_results:
            if isinstance(result, tuple):  # Success case
                results.append((id,) + result)
            else:  # Error case
                missing_files.append(result)
    
    # Report any missing files
    if missing_files:
        print("Warning: The following files were missing or had errors:")
        for msg in missing_files:
            print(f"- {msg}")
    
    return results

def batchSuperimposeParallel2(pred_folder, gt_folder, export=False, n_processes=None, method='mican'):
    """Performs parallel batch superimposition of protein structures
    Processes pairs of structures from predictions and groundtruth folders

    Args:
        pred_folder (str) path containing groundtruth 
        gt_folder (str): path containing predictions
        export (bool, optional): save superimposed structures. Defaults to False
        n_processes (int, optional): number of processes to use. Defaults to CPU count
        method (str, optional): superimposition method. 'mican' or 'svd'. Defaults to 'mican'

    Returns:
        results (list): list of tuples containing:
            - id (str): structure identifier
            - rmsd (float): Root Mean Square Deviation
            - tm_score (float): TM-score of alignment
    """
    # Dictionary to hold the pairs of files
    file_pairs = {}
    
    # Check if both directories exist
    if not os.path.exists(pred_folder) or not os.path.exists(gt_folder):
        raise ValueError(f"One or both directories missing:\nPredictions: {pred_folder}\nGround truth: {gt_folder}")
    
    # Get files from predictions folder
    pred_files = os.listdir(pred_folder)
    for file in pred_files:
        if file.endswith('.pdb'):
            base_name = file[:-4]  # Remove the .pdb extension
            file_pairs[base_name] = {'pred': file}
    
    # Get files from ground truth folder and match with predictions
    gt_files = os.listdir(gt_folder)
    for file in gt_files:
        if file.endswith('.pdb'):
            base_name = file[:-4]  # Remove the .pdb extension
            file_pairs[base_name]['gt'] = file
    
    # Prepare arguments for parallel processing
    process_args = [
        (id, files, pred_folder, gt_folder, export, method)
        for id, files in file_pairs.items()
    ]
    
    # Initialize results and missing files lists
    results = []
    missing_files = []
    
    # Process files in parallel
    with Pool(processes=n_processes) as pool:
        processed_results = pool.map(process_pair, process_args)
        
        # Sort results and errors
        for id, result in processed_results:
            if isinstance(result, tuple):  # Success case
                results.append((id,) + result)
            else:  # Error case
                missing_files.append(result)
    
    # Report any missing files
    if missing_files:
        print("Warning: The following files were missing or had errors:")
        for msg in missing_files:
            print(f"- {msg}")
    
    return results

#########
def CalculateGdtFromSuperimposedFile(super_pdb_path):
    """reads a superimposed pdb file 
    and calculates GDT_TS, GTD-scores, RMSD, and TM-score

    Args:
        super_pdb_path (str): filepath to the superimposed pdb file

    Returns:
        gdt_values (list): list of GDT scores for cutoffs [1.0, 2.0, 4.0, 8.0]
        gdt_ts (float): GDT total score (mean of GDT scores)
        rmsd (float): Root Mean Square Deviation
        tm_score (float): TM-score of alignment
    """
    # Parse the PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('pdb', super_pdb_path)

    # Initialize lists to hold alpha carbon atoms for each model
    atoms_model1 = []
    atoms_model2 = []

    # Iterate through the models in the structure
    for model_index, model in enumerate(structure):
        print(f"\nStructure {model_index + 1} chains:")
        for chain in model:
            n_residues = len([res for res in chain if res.has_id('CA')])
            print(f"Chain {chain.id}: {n_residues} CA atoms")
            
            # Extract alpha carbon atoms for the current model
            atoms = [atom for residue in chain if residue.has_id('CA') for atom in residue if atom.get_id() == 'CA']
            
            # Append atoms to the corresponding model list
            if model_index == 0:
                atoms_model1.extend(atoms)
            elif model_index == 1:
                atoms_model2.extend(atoms)


    # Get coordinates of superimposed atoms
    coords1 = np.array([atom.get_coord() for atom in atoms_model1])
    coords2 = np.array([atom.get_coord() for atom in atoms_model2])

    # Calculate distances
    distances = np.linalg.norm(coords1 - coords2, axis=1)

    # Define cutoffs for GDT
    cutoffs = [1.0, 2.0, 4.0, 8.0]

    # Calculate GDT for each cutoff
    gdt_values = []
    for cutoff in cutoffs:
        gdt = np.mean(distances <= cutoff)
        gdt_values.append(gdt)

    # Calculate GDT_TS (total score)
    gdt_ts = np.mean(gdt_values)

    # Calculate RMSD
    rmsd = np.sqrt(np.mean((coords1 - coords2) ** 2))

    # Calculate TM-score
    L = len(atoms_model1)
    d0 = 1.24 * (L - 15) ** (1/3) - 1.8
    if d0 < 0.5:
        d0 = 0.5

    tm_score = np.sum(1.0 / (1.0 + (distances / d0) ** 2)) / L

    return gdt_values, gdt_ts, rmsd, tm_score


#biopython method for superimpose
def SuperimposeSVD(pdb1_path, pdb2_path, verbose=False):
    """superimposes two protein structures using SVD method
    and calculates GDT_TS, GTD-scores, RMSD, and TM-score

    Args:
        pdb1_path (str): filepath to first PDB file
        pdb2_path (str): filepath to second PDB file
        verbose (bool, optional): print chain information. Defaults to False

    Returns:
        gdt_values (list): list of GDT scores for cutoffs [1.0, 2.0, 4.0, 8.0]
        gdt_ts (float): GDT total score (mean of GDT scores)
        rmsd (float): Root Mean Square Deviation
        tm_score (float): TM-score of alignment
    """
    # Parse PDB files
    parser = PDBParser(QUIET=True)
    structure1 = parser.get_structure('pdb1', pdb1_path)
    structure2 = parser.get_structure('pdb2', pdb2_path)

    # Print chain information
    if verbose:
        print("Structure 1 chains:")
        for model in structure1:
            for chain in model:
                n_residues = len([res for res in chain if res.has_id('CA')])
                print(f"Chain {chain.id}: {n_residues} CA atoms")
                
        print("\nStructure 2 chains:")
        for model in structure2:
            for chain in model:
                n_residues = len([res for res in chain if res.has_id('CA')])
                print(f"Chain {chain.id}: {n_residues} CA atoms")
                
    
    # Extract alpha carbon atoms
    atoms1 = [atom for model in structure1 for chain in model for residue in chain if residue.has_id('CA') for atom in residue if atom.get_id() == 'CA']
    atoms2 = [atom for model in structure2 for chain in model for residue in chain if residue.has_id('CA') for atom in residue if atom.get_id() == 'CA']

    # Superimpose structures
    superimposer = Superimposer()
    superimposer.set_atoms(atoms1, atoms2)
    superimposer.apply(structure2.get_atoms())

    # Get coordinates of superimposed atoms
    coords1 = np.array([atom.get_coord() for atom in atoms1])
    coords2 = np.array([atom.get_coord() for atom in atoms2])

    # Calculate distances
    distances = np.linalg.norm(coords1 - coords2, axis=1)

    # Define cutoffs for GDT
    cutoffs = [1.0, 2.0, 4.0, 8.0]

    # Calculate GDT for each cutoff
    gdt_values = []
    for cutoff in cutoffs:
        gdt = np.mean(distances <= cutoff)
        gdt_values.append(round(float(gdt), 4))

    # Calculate GDT_TS (total score)
    gdt_ts = np.mean(gdt_values)
    gdt_ts = round(float(gdt_ts), 4)

    # Calculate RMSD
    rmsd = np.sqrt(np.mean((coords1 - coords2) ** 2))
    rmsd = round(float(rmsd), 4)
    
    # Calculate TM-score
    L = len(atoms1)
    d0 = 1.24 * (L - 15) ** (1/3) - 1.8
    if d0 < 0.5:
        d0 = 0.5

    tm_score = np.sum(1.0 / (1.0 + (distances / d0) ** 2)) / L
    tm_score = round(float(tm_score), 4)
    
    return gdt_values, gdt_ts, rmsd, tm_score


def extract_af_metrics(json_file, print_results=False):
    """Extracts metrics from AlphaFold JSON file
    Args:
        json_file (str): path to AlphaFold JSON file
    Returns:
        dict: dictionary containing the extracted metrics
    """
    
    with open(json_file, 'r') as file:
        json_data = json.load(file)  # Parse JSON file into a dictionary

    results = {
        'chain_iptm': [],       # List to store individual chain iPTM values
        'iptm': None,           # Single value for global iPTM
        'chain_ptm': [],        # List to store individual chain PTM values
        'ptm': None,            # Single value for global PTM
        'has_clash': None,      # Single value for has_clash
        'num_recycles' : None,  # Number of recycles'
        'ranking_score': None   # Ranking score
    }

    # Extract chain iPTM values (filtering out None values)
    chain_iptm = json_data.get('chain_iptm', [])
    results['chain_iptm'] = [x for x in chain_iptm if x is not None]

    # Extract global iPTM
    results['iptm'] = json_data.get('iptm')
    
    # Extract chain PTM values (filtering out None values)
    chain_ptm = json_data.get('chain_ptm', [])
    results['chain_ptm'] = [x for x in chain_ptm if x is not None]

    # Extract global PTM
    results['ptm'] = json_data.get('ptm')

    # Extract pLDDT
    results['has_clash'] = json_data.get('has_clash')

    # Extract ranking score
    results['num_recycles'] = json_data.get('num_recycles')

    # Extract ranking score
    results['ranking_score'] = json_data.get('ranking_score')
    
    # Print results if requested
    if print_results:
        if results:
            for key, value in results.items():
                print(f"{key}: {value}")
    
    return results




# helper function for development
def calculate_chain_scores(pdb1_path, pdb2_path):
    """calculates RMSD scores for each matching chain
    between two protein structures

    Args:
        pdb1_path (str): filepath to first PDB file
        pdb2_path (str): filepath to second PDB file

    Returns:
        results (dict): dictionary with chain IDs as keys and RMSD scores as values
    """
    parser = PDBParser(QUIET=True)
    structure1 = parser.get_structure('pdb1', pdb1_path)
    structure2 = parser.get_structure('pdb2', pdb2_path)
    
    results = {}
    
    # Compare each chain separately
    for chain1 in structure1[0]:  # First model
        chain_id = chain1.id
        if chain_id in [chain.id for chain in structure2[0]]:
            # Get atoms for this chain only
            atoms1 = [atom for residue in chain1 if residue.has_id('CA') 
                     for atom in residue if atom.get_id() == 'CA']
            atoms2 = [atom for residue in structure2[0][chain_id] if residue.has_id('CA') 
                     for atom in residue if atom.get_id() == 'CA']
            
            sup = Superimposer()
            sup.set_atoms(atoms1, atoms2)
            rmsd = sup.rms
            
            results[chain_id] = rmsd
            
    return results