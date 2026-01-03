import os
import glob
import pandas as pd
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def combine_protein_ligand_structures(df, 
                                     protein_dir, 
                                     ligand_base_dir, 
                                     output_dir,
                                     protein_col='Protein',
                                     ligand_col='Ligand'):
    """
    Combine protein PDB files with ligand SDF docking poses.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing protein and ligand names
    protein_dir : str
        Directory containing protein PDB files
    ligand_base_dir : str
        Base directory containing ligand pose SDF files in subdirectories
    output_dir : str
        Base directory for output combined PDB files
    protein_col : str, optional
        Name of the column in the dataframe containing protein names (default: 'Protein')
    ligand_col : str, optional
        Name of the column in the dataframe containing ligand names (default: 'Ligand')
    
    Returns:
    --------
    dict
        A summary of the processing results
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Validate input
    if not os.path.exists(protein_dir):
        logger.error(f"Protein directory does not exist: {protein_dir}")
        return {"status": "error", "message": f"Protein directory does not exist: {protein_dir}"}
    
    if not os.path.exists(ligand_base_dir):
        logger.error(f"Ligand base directory does not exist: {ligand_base_dir}")
        return {"status": "error", "message": f"Ligand base directory does not exist: {ligand_base_dir}"}
    
    # Initialize results summary
    results_summary = {
        "total_combinations": 0,
        "successful": 0,
        "failed": 0,
        "details": []
    }
    
    # Process each row in the dataframe
    for idx, row in df.iterrows():
        protein_name = row[protein_col]
        ligand_name = row[ligand_col]
        
        protein_file = os.path.join(protein_dir, f"{protein_name}.pdb")
        ligand_dir = os.path.join(ligand_base_dir, ligand_name)
        
        # Check if protein file exists
        if not os.path.exists(protein_file):
            error_msg = f"Protein file not found: {protein_file}"
            logger.error(error_msg)
            results_summary["failed"] += 1
            results_summary["details"].append({
                "protein": protein_name,
                "ligand": ligand_name,
                "status": "error",
                "message": error_msg
            })
            continue
        
        # Check if ligand directory exists
        if not os.path.exists(ligand_dir):
            error_msg = f"Ligand directory not found: {ligand_dir}"
            logger.error(error_msg)
            results_summary["failed"] += 1
            results_summary["details"].append({
                "protein": protein_name,
                "ligand": ligand_name,
                "status": "error",
                "message": error_msg
            })
            continue
        
        # Find all pose files for this ligand
        pose_files = glob.glob(os.path.join(ligand_dir, "pose_*.sdf"))
        
        if not pose_files:
            error_msg = f"No pose files found in: {ligand_dir}"
            logger.warning(error_msg)
            results_summary["failed"] += 1
            results_summary["details"].append({
                "protein": protein_name,
                "ligand": ligand_name,
                "status": "error",
                "message": error_msg
            })
            continue
        
        # Create output directory for this ligand
        ligand_output_dir = os.path.join(output_dir, ligand_name)
        os.makedirs(ligand_output_dir, exist_ok=True)
        
        # Read protein structure
        try:
            # Set proximityBonding=False to avoid creating bonds automatically
            protein = Chem.MolFromPDBFile(protein_file, sanitize=False, proximityBonding=False)
            if protein is None:
                error_msg = f"Failed to read protein file: {protein_file}"
                logger.error(error_msg)
                results_summary["failed"] += 1
                results_summary["details"].append({
                    "protein": protein_name,
                    "ligand": ligand_name,
                    "status": "error",
                    "message": error_msg
                })
                continue
        except Exception as e:
            error_msg = f"Error reading protein file {protein_file}: {str(e)}"
            logger.error(error_msg)
            results_summary["failed"] += 1
            results_summary["details"].append({
                "protein": protein_name,
                "ligand": ligand_name,
                "status": "error",
                "message": error_msg
            })
            continue
        
        # Process each pose file
        successful_poses = 0
        for pose_file in pose_files:
            pose_filename = os.path.basename(pose_file)
            output_filename = pose_filename.replace('.sdf', '.pdb')
            output_file = os.path.join(ligand_output_dir, output_filename)
            
            try:
                # Read ligand structure
                ligand = Chem.SDMolSupplier(pose_file, sanitize=False)[0]
                if ligand is None:
                    error_msg = f"Failed to read ligand file: {pose_file}"
                    logger.error(error_msg)
                    continue
                
                # Combine structures by creating a combined file
                # We'll write both to separate PDB blocks and then combine them in a single file
                # This preserves connectivity information better than creating a combined molecule
                
                # Write protein to PDB block
                protein_block = Chem.MolToPDBBlock(protein)
                
                # Add chain ID to ligand atoms (use 'L' for ligand)
                for atom in ligand.GetAtoms():
                    atom.SetProp("_CIPCode", "L")
                
                # Convert ligand to PDB format
                ligand_block = Chem.MolToPDBBlock(ligand)
                
                # Combine blocks (remove END from protein block)
                combined_block = protein_block.replace("END", "") + ligand_block
                
                # Write to output file
                with open(output_file, 'w') as f:
                    f.write(combined_block)
                
                successful_poses += 1
                logger.info(f"Created combined structure: {output_file}")
                
            except Exception as e:
                error_msg = f"Error processing {pose_file}: {str(e)}"
                logger.error(error_msg)
                results_summary["failed"] += 1
                results_summary["details"].append({
                    "protein": protein_name,
                    "ligand": ligand_name,
                    "pose": pose_filename,
                    "status": "error",
                    "message": error_msg
                })
                continue
        
        results_summary["total_combinations"] += len(pose_files)
        results_summary["successful"] += successful_poses
        results_summary["failed"] += (len(pose_files) - successful_poses)
        
        results_summary["details"].append({
            "protein": protein_name,
            "ligand": ligand_name,
            "total_poses": len(pose_files),
            "successful_poses": successful_poses,
            "status": "success" if successful_poses > 0 else "error",
            "output_dir": ligand_output_dir
        })
    
    # Save summary report
    summary_path = os.path.join(output_dir, "processing_summary.json")
    import json
    with open(summary_path, "w") as f:
        json.dump(results_summary, f, indent=2)
    
    logger.info(f"Processing complete. Successful: {results_summary['successful']}, Failed: {results_summary['failed']}")
    logger.info(f"Summary saved to {summary_path}")
    
    return results_summary

# FG: How to use in a notebook
"""
import pandas as pd
from combine_structures import combine_protein_ligand_structures

# Load your dataframe with protein and ligand information
df = pd.read_csv('your_combinations.csv')

# Define directories
protein_dir = "data/raw/docking/prep_Flare"
ligand_base_dir = "data/interim/docking/Res_Diffdock_Flare_R1"
output_dir = "data/interim/docking/Combined_Structures"

# Run the combination process
results = combine_protein_ligand_structures(
    df=df,
    protein_dir=protein_dir,
    ligand_base_dir=ligand_base_dir,
    output_dir=output_dir,
    protein_col='Protein',  # Customize if your column names are different
    ligand_col='Ligand'
)

# Display summary
print(f"Successfully combined {results['successful']} out of {results['total_combinations']} structures")
"""

# Alternative approach using a more advanced merge method
def combine_protein_ligand_advanced(df, 
                                   protein_dir, 
                                   ligand_base_dir, 
                                   output_dir,
                                   protein_col='Protein',
                                   ligand_col='Ligand'):
    """
    An alternative approach using RDKit's CombineMols function to merge structures.
    Note: This may not preserve all connectivity information as well as the block-based approach.
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize results summary
    results_summary = {
        "total_combinations": 0,
        "successful": 0,
        "failed": 0,
        "details": []
    }
    
    # Process each row in the dataframe
    for idx, row in df.iterrows():
        protein_name = row[protein_col]
        ligand_name = row[ligand_col]
        
        protein_file = os.path.join(protein_dir, f"{protein_name}.pdb")
        ligand_dir = os.path.join(ligand_base_dir, ligand_name)
        
        # Skip if files don't exist (with error logging)
        if not os.path.exists(protein_file) or not os.path.exists(ligand_dir):
            continue
        
        # Find all pose files for this ligand
        pose_files = glob.glob(os.path.join(ligand_dir, "pose_*.sdf"))
        
        if not pose_files:
            continue
        
        # Create output directory for this ligand
        ligand_output_dir = os.path.join(output_dir, ligand_name)
        os.makedirs(ligand_output_dir, exist_ok=True)
        
        # Read protein structure
        try:
            protein = Chem.MolFromPDBFile(protein_file, sanitize=False, proximityBonding=False)
            if protein is None:
                continue
                
            # Process each pose file
            for pose_file in pose_files:
                pose_filename = os.path.basename(pose_file)
                output_filename = pose_filename.replace('.sdf', '.pdb')
                output_file = os.path.join(ligand_output_dir, output_filename)
                
                try:
                    # Read ligand structure
                    ligand = Chem.SDMolSupplier(pose_file, sanitize=False)[0]
                    if ligand is None:
                        continue
                    
                    # Combine molecules
                    combined = Chem.CombineMols(protein, ligand)
                    
                    # Write to output file
                    Chem.MolToPDBFile(combined, output_file)
                    
                    results_summary["successful"] += 1
                    
                except Exception as e:
                    results_summary["failed"] += 1
                    continue
                    
        except Exception as e:
            results_summary["failed"] += len(pose_files)
            continue
    
    return results_summary