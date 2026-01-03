import time
import os
import shutil
import json
import requests
import subprocess
import time
import logging
import pandas as pd
from pathlib import Path
from loguru import logger
from datetime import datetime


# preparing and formatting the input arguments for DiffDock inference request, and saving the results in output.json file
def run_diffdock(query_url, protein_file_path, ligand_file_path, num_poses=10, time_divisions=20, steps=18, save_trajectory=False, is_staged=False):
    """
    Main function to run the molecular docking
    :param query_url: str, the url to send the request to
    :param protein_file_path: str, path to the protein file
    :param ligand_file_path: str, path to the ligand file
    return JSON response
    """

    protein_bytes = file_to_json_compatible_string(protein_file_path)
    ligand_bytes = file_to_json_compatible_string(ligand_file_path)

    data = {
        "ligand": ligand_bytes,
        "ligand_file_type": "sdf",
        "protein": protein_bytes,
        "num_poses": num_poses,
        "time_divisions": time_divisions,
        "steps": steps,
        "save_trajectory": save_trajectory,  # diffusion trajectory
        "is_staged": is_staged
    }

    headers = {"Content-Type": "application/json"}

    response = requests.post(query_url, headers=headers, json=data)

    if response.status_code == 200:
        print("Request successful, output saved to output.json")
    else:
        print(f"Request failed with status code {response.status_code}")
        print("Response:", response.text)

    return response.json()

# reading in the input PDB/SDF/SMILES files as a string to be used for JSON request
def file_to_json_compatible_string(file_path):
    """
    Convert PDB file and sdf file to JSON
    """
    with open(file_path, 'r') as file:
        content_str = file.read()
    return content_str

############################################################# BATCH Processing



def run_diffdock_batch(csv_path, 
                       ligand_dir, 
                       protein_dir, 
                       output_base_dir, 
                       query_url,
                       ligand_col='ligand_name',
                       protein_col='protein_name',
                       num_poses=10, 
                       time_divisions=20, 
                       steps=20, 
                       save_trajectory=False,
                       is_staged=False):
    """
    Run DiffDock for each ligand-protein combination specified in a CSV file.
    
    Parameters:
    -----------
    csv_path : str
        Path to the CSV file containing ligand and protein combinations
    ligand_dir : str
        Path to the directory containing ligand SDF files
    protein_dir : str
        Path to the directory containing protein PDB files
    output_base_dir : str
        Base directory for output results
    query_url : str
        URL for the DiffDock inference API
    ligand_col : str, optional
        Name of the column in the CSV file that contains ligand names (default: 'ligand_name')
    protein_col : str, optional
        Name of the column in the CSV file that contains protein names (default: 'protein_name')
    
    Returns:
    --------
    dict
        A summary of the processing results
    """

    # Get the current date and time
    now = datetime.now()
    # Format the timestamp
    timestamp = now.strftime("%Y%m%d_%H%M%S")
    # Define the filename
    filename = f"docking_experiment_{timestamp}.log"
    # Configure the logger
    logger.add(filename, rotation="1 MB", level="DEBUG", format="{time} - {level} - {message}")

    
    # Ensure directories exist
    os.makedirs(output_base_dir, exist_ok=True)
    
    # Load combinations from CSV
    try:
        df = pd.read_csv(csv_path)
        logger.info(f"Loaded {len(df)} ligand-protein combinations from {csv_path}")
    except Exception as e:
        logger.error(f"Failed to load CSV file: {e}")
        return {"status": "error", "message": f"Failed to load CSV file: {e}"}
    
    # Validate required columns
    required_columns = [ligand_col, protein_col]
    if not all(col in df.columns for col in required_columns):
        logger.error(f"CSV missing required columns. Expected: {required_columns}")
        return {"status": "error", "message": f"CSV missing required columns. Expected: {required_columns}"}
    
    # Process each combination
    results_summary = {
        "total": len(df),
        "successful": 0,
        "failed": 0,
        "details": []
    }
    
    for index, row in df.iterrows():
        ligand_name = row[ligand_col]
        protein_name = row[protein_col]
        
        # Define file paths
        ligand_file = f"{ligand_name}.sdf"
        protein_file = f"{protein_name}.pdb"
        
        ligand_filepath = os.path.join(ligand_dir, ligand_file)
        protein_filepath = os.path.join(protein_dir, protein_file)
        
        # Organize output by protein
        #protein_output_dir = os.path.join(output_base_dir, protein_name)
        ligand_output_dir = os.path.join(output_base_dir, ligand_name)
        
        # Validate input files exist
        if not os.path.exists(ligand_filepath):
            error_msg = f"Ligand file not found: {ligand_filepath}"
            logger.error(error_msg)
            results_summary["failed"] += 1
            results_summary["details"].append({
                "ligand": ligand_name,
                "protein": protein_name,
                "status": "error",
                "message": error_msg
            })
            continue
            
        if not os.path.exists(protein_filepath):
            error_msg = f"Protein file not found: {protein_filepath}"
            logger.error(error_msg)
            results_summary["failed"] += 1
            results_summary["details"].append({
                "ligand": ligand_name,
                "protein": protein_name,
                "status": "error",
                "message": error_msg
            })
            continue
        
        # Create output directories
        os.makedirs(ligand_output_dir, exist_ok=True)
        
        logger.info(f"Processing {index+1}/{len(df)}: {ligand_name} with {protein_name}")
        
        try:
            # Start measuring time
            start = time.time()
            
            # Run DiffDock
            result = run_diffdock(
                query_url=query_url,
                protein_file_path=protein_filepath,
                ligand_file_path=ligand_filepath,
                num_poses=num_poses,
                time_divisions=time_divisions,
                steps=steps,
                save_trajectory=save_trajectory,
                is_staged=is_staged
            )
            
            # End measuring time
            end = time.time()
            processing_time = end - start
            logger.debug(f"{ligand_name} with {protein_name} took {processing_time:.2f} seconds")
            
            # Save result to output.json
            with open(f"{ligand_output_dir}/output.json", "w") as f:
                json.dump(result, f)
            
            # Save ligand positions
            for i, ligand_geometry in enumerate(result["ligand_positions"]):
                with open(f"{ligand_output_dir}/pose_{i}.sdf", "w") as f:
                    f.write(ligand_geometry)
            
            # Update results summary
            results_summary["successful"] += 1
            results_summary["details"].append({
                "ligand": ligand_name,
                "protein": protein_name,
                "status": "success",
                "processing_time": processing_time,
                "output_dir": ligand_output_dir
            })
            
        except Exception as e:
            logger.error(f"Error processing {ligand_name} with {protein_name}: {e}")
            results_summary["failed"] += 1
            results_summary["details"].append({
                "ligand": ligand_name,
                "protein": protein_name,
                "status": "error",
                "message": str(e)
            })
    
    # Save summary report
    summary_path = os.path.join(output_base_dir, "processing_summary.json")
    with open(summary_path, "w") as f:
        json.dump(results_summary, f, indent=2)
    
    logger.info(f"Processing complete. Successful: {results_summary['successful']}, Failed: {results_summary['failed']}")
    logger.info(f"Summary saved to {summary_path}")
    
    return results_summary

# Example usage
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Run DiffDock for multiple ligand-protein combinations')
    parser.add_argument('--csv', required=True, help='Path to CSV file with ligand-protein combinations')
    parser.add_argument('--ligand_dir', required=True, help='Directory containing ligand SDF files')
    parser.add_argument('--protein_dir', required=True, help='Directory containing protein PDB files')
    parser.add_argument('--output_dir', required=True, help='Base directory for outputs')
    parser.add_argument('--query_url', required=True, help='DiffDock API query URL')
    parser.add_argument('--ligand_col', default='ligand_name', help='Column name for ligand in CSV (default: ligand_name)')
    parser.add_argument('--protein_col', default='protein_name', help='Column name for protein in CSV (default: protein_name)')
    
    args = parser.parse_args()
    
    # Run the batch processing
    run_diffdock_batch(
        csv_path=args.csv,
        ligand_dir=args.ligand_dir,
        protein_dir=args.protein_dir,
        output_base_dir=args.output_dir,
        query_url=args.query_url,
        ligand_col=args.ligand_col,
        protein_col=args.protein_col
    )