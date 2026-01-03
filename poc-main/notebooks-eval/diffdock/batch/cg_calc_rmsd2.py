import os
import glob
import re
import warnings
import pandas as pd
import numpy as np
from tqdm import tqdm  # Optional, for progress bar
from rdkit import Chem
from rdkit.Chem import AllChem
from Bio.PDB import PDBParser, Superimposer, Selection
from Bio.PDB.PDBExceptions import PDBConstructionWarning

# Suppress BioPython warnings
warnings.simplefilter('ignore', PDBConstructionWarning)

###################################
def calculate_docking_rmsd(pdb_result1, ligand_result1, pdb_result2, ligand_resname=None, model_id=0, output_pdb=None):
    """
    Calculate RMSD between two docking poses after aligning the protein structures.
    
    Parameters:
    -----------
    pdb_result1 : str
        Path to the PDB file of the protein from docking result 1
    ligand_result1 : str
        Path to the ligand file (SDF or PDB) from docking result 1
    pdb_result2 : str
        Path to the PDB file containing both protein and ligand from docking result 2
    ligand_resname : str, optional
        Residue name of the ligand in the PDB files. If None, will try to auto-detect
        a numeric 3-character residue name from HETATM records.
    model_id : int
        Model ID to use from the PDB file (default 0)
    output_pdb : str, optional
        Path to save the superimposed structures as a PDB file. If None, no file is saved.
        
    Returns:
    --------
    float
        RMSD value between the two ligand poses after protein alignment
    """
    # 1. Parse protein structures using BioPython
    parser = PDBParser()
    structure1 = parser.get_structure("result1", pdb_result1)
    structure2 = parser.get_structure("result2", pdb_result2)
    
    # Use specified model
    model1 = structure1[model_id]
    model2 = structure2[model_id]
    
    # 2. Extract protein atoms for alignment (exclude HETATM)
    protein_atoms1 = []
    protein_atoms2 = []
    
    # Get CA atoms from each protein for alignment
    # Only include ATOM records (not HETATM)
    for chain in model1:
        for residue in chain:
            if residue.id[0] == " ":  # Regular amino acid (not HETATM, water, etc.)
                if "CA" in residue:
                    protein_atoms1.append(residue["CA"])
                        
    for chain in model2:
        for residue in chain:
            if residue.id[0] == " ":  # Regular amino acid (not HETATM, water, etc.)
                if "CA" in residue:
                    protein_atoms2.append(residue["CA"])
    
    # 3. Align proteins using BioPython Superimposer
    super_imposer = Superimposer()
    
    # Make sure we have the same number of atoms for alignment
    min_atoms = min(len(protein_atoms1), len(protein_atoms2))
    super_imposer.set_atoms(protein_atoms1[:min_atoms], protein_atoms2[:min_atoms])
    
    # Apply rotation/translation to entire second structure
    super_imposer.apply(model2.get_atoms())
    
    # 4. Extract ligand coordinates from first structure
    # Check file extension to determine how to handle ligand_result1
    ligand_file_ext = os.path.splitext(ligand_result1.lower())[1]
    
    if ligand_file_ext == '.sdf':
        # For SDF files, read using RDKit
        ligand1 = Chem.SDMolSupplier(ligand_result1, removeHs=False)[0]
        if ligand1 is None:
            raise ValueError(f"Could not read molecule from {ligand_result1}")
        
        # Get coordinates from first ligand
        ligand1_conf = ligand1.GetConformer()
        ligand1_coords = np.array([ligand1_conf.GetAtomPosition(i) for i in range(ligand1.GetNumAtoms())])
        
        # Get element types for atom matching
        sdf_elements = [atom.GetSymbol() for atom in ligand1.GetAtoms()]
        
    elif ligand_file_ext == '.pdb':
        # For PDB files, read using BioPython
        ligand1_structure = parser.get_structure("ligand1", ligand_result1)
        ligand1_model = ligand1_structure[model_id]
        
        # Extract all HETATM atoms from the PDB file
        ligand1_atoms = []
        sdf_elements = []
        
        for chain in ligand1_model:
            for residue in chain:
                if residue.id[0] != " ":  # HETATM records
                    for atom in residue:
                        ligand1_atoms.append(atom)
                        sdf_elements.append(atom.element)
        
        if not ligand1_atoms:
            raise ValueError(f"No ligand atoms (HETATM) found in {ligand_result1}")
        
        ligand1_coords = np.array([atom.coord for atom in ligand1_atoms])
    else:
        raise ValueError(f"Unsupported ligand file format: {ligand_file_ext}. Must be .sdf or .pdb")
    
    # 5. Extract ligand from second structure (already aligned)
    # If ligand_resname is not provided, try to auto-detect it
    ligand2_coords = []
    ligand2_atoms = []
    detected_ligand_resname = None
    
    if ligand_resname is None:
        # First, scan for a numeric 3-character residue name in HETATM records
        potential_ligands = {}
        for chain in model2:
            for residue in chain:
                # Check if it's a HETATM
                if residue.id[0] != " ":
                    # Try to check if it's potentially a numeric resname - SAFELY
                    try:
                        # Only consider it if all characters are digits
                        if all(c.isdigit() for c in residue.resname) and len(residue.resname) == 3:
                            if residue.resname not in potential_ligands:
                                potential_ligands[residue.resname] = 0
                            potential_ligands[residue.resname] += len(residue)
                    except:
                        # Skip any residue that causes an error
                        pass
        
        # Choose the ligand with the most atoms (most likely the actual ligand, not water, etc.)
        if potential_ligands:
            detected_ligand_resname = max(potential_ligands.items(), key=lambda x: x[1])[0]
            print(f"Auto-detected ligand residue name: {detected_ligand_resname}")
            ligand_resname = detected_ligand_resname
        else:
            # If no numeric resname found, try to find largest HETATM residue
            hetatm_residues = {}
            for chain in model2:
                for residue in chain:
                    if residue.id[0] != " ":  # HETATM record
                        if residue.resname not in hetatm_residues:
                            hetatm_residues[residue.resname] = 0
                        hetatm_residues[residue.resname] += len(residue)
            
            if hetatm_residues:
                # Exclude common waters and ions
                common_waters_ions = {"HOH", "WAT", "H2O", "NA", "CL", "MG", "CA", "ZN", "K"}
                filtered_hetatms = {resname: count for resname, count in hetatm_residues.items() 
                                   if resname not in common_waters_ions}
                
                if filtered_hetatms:
                    detected_ligand_resname = max(filtered_hetatms.items(), key=lambda x: x[1])[0]
                    print(f"Auto-detected non-numeric ligand residue name: {detected_ligand_resname}")
                    ligand_resname = detected_ligand_resname
    
    if ligand_resname is None:
        raise ValueError("Could not auto-detect ligand residue name. Please specify it manually.")
    
    # Now extract the ligand with the determined resname
    for chain in model2:
        for residue in chain:
            # Check if it's a HETATM and matches the ligand residue name
            if residue.id[0] != " " and residue.resname == ligand_resname:
                for atom in residue:
                    ligand2_atoms.append(atom)
                    ligand2_coords.append(atom.coord)
    
    if not ligand2_coords:
        if detected_ligand_resname:
            raise ValueError(f"Error accessing auto-detected ligand with residue name {ligand_resname}")
        else:
            raise ValueError(f"No ligand with residue name {ligand_resname} found in {pdb_result2}")
    
    ligand2_coords = np.array(ligand2_coords)
    
    # 6. Extract element types from PDB ligand
    pdb_elements = [atom.element for atom in ligand2_atoms]
    
    # 7. Create mapping between atoms based on element
    # This is a simple approach - for complex molecules you'd need more sophisticated matching
    atom_mapping = []
    used_pdb_indices = set()
    
    for i, sdf_elem in enumerate(sdf_elements):
        for j, pdb_elem in enumerate(pdb_elements):
            if j not in used_pdb_indices and sdf_elem == pdb_elem:
                atom_mapping.append((i, j))
                used_pdb_indices.add(j)
                break
    
    if not atom_mapping:
        print("Warning: Could not match atoms between structures")
        return None
    
    # 8. Calculate RMSD using matched atoms
    squared_diff = 0
    for sdf_idx, pdb_idx in atom_mapping:
        squared_diff += np.sum((ligand1_coords[sdf_idx] - ligand2_coords[pdb_idx])**2)
    
    rmsd = np.sqrt(squared_diff / len(atom_mapping))
    
    # 9. Export superimposed structures if output_pdb is provided
    if output_pdb:
        # Create a new structure to hold the superimposed complexes
        from Bio.PDB.Structure import Structure
        from Bio.PDB.Model import Model
        
        # Create new structure and model
        combined = Structure("combined")
        combined_model = Model(0)
        combined.add(combined_model)
        
        # Add protein from structure1 (renamed chains to avoid conflicts)
        # Use single character chain IDs and cycle through letters if needed
        import string
        available_chain_ids = list(string.ascii_uppercase + string.ascii_lowercase + string.digits)
        used_chain_ids = set()
        
        # Add protein from structure1
        for i, chain in enumerate(model1):
            # Generate a new unique single-character chain ID
            new_chain_id = available_chain_ids[i % len(available_chain_ids)]
            while new_chain_id in used_chain_ids:
                i += 1
                new_chain_id = available_chain_ids[i % len(available_chain_ids)]
            
            used_chain_ids.add(new_chain_id)
            
            # Make a deep copy to avoid modifying the original
            import copy
            new_chain = copy.deepcopy(chain)
            new_chain.id = new_chain_id
            combined_model.add(new_chain)
        
        # Add protein from structure2 with different chain IDs
        for i, chain in enumerate(model2):
            # Generate a new unique single-character chain ID
            offset = len(used_chain_ids)  # Start from where structure1 ended
            new_chain_id = available_chain_ids[(i + offset) % len(available_chain_ids)]
            while new_chain_id in used_chain_ids:
                i += 1
                new_chain_id = available_chain_ids[(i + offset) % len(available_chain_ids)]
            
            used_chain_ids.add(new_chain_id)
            
            # Make a deep copy to avoid modifying the original
            import copy
            new_chain = copy.deepcopy(chain)
            new_chain.id = new_chain_id
            combined_model.add(new_chain)
        
        # Use PDBIO to write out the combined structure
        from Bio.PDB import PDBIO
        io = PDBIO()
        io.set_structure(combined)
        io.save(output_pdb)
        print(f"Superimposed structures saved to {output_pdb}")
    
    return rmsd


#################################################
def calculate_rmsd_for_directory(
    pdb_result,
    ligand_directory,
    pdb_complex,
    ligand_resname=None,
    file_pattern=None,
    output_dir=None,
    save_superimposed=False
    ):
    """
    Calculate RMSD for multiple ligand files in a directory against the same reference PDB complex.
    
    Parameters:
    -----------
    pdb_result : str
        Path to the PDB file of the protein from docking result
    ligand_directory : str
        Path to directory containing ligand files (SDF or PDB) to evaluate
    pdb_complex : str
        Path to the PDB file containing both protein and ligand complex (reference)
    ligand_resname : str
        Residue name of the ligand in the PDB files (default None = Auto determined)
    file_pattern : str or list, optional
        Pattern(s) to match ligand files. If None, will match "*.sdf" and "*.pdb"
    output_dir : str, optional
        Directory to save superimposed PDB files. If None and save_superimposed is True,
        they will be saved in a subfolder named "superimposed" within ligand_directory.
    save_superimposed : bool
        Whether to save superimposed structures as PDB files (default False)
    
    Returns:
    --------
    pandas.DataFrame
        DataFrame containing file names and their corresponding RMSD values
    """
    # If file_pattern is not specified, use both SDF and PDB patterns
    if file_pattern is None:
        file_patterns = ["*.sdf", "*.pdb"]
    elif isinstance(file_pattern, str):
        file_patterns = [file_pattern]
    else:
        file_patterns = file_pattern
    
    # Find all matching files in the directory, ensuring no duplicates
    ligand_files = set()  # Use a set to automatically remove duplicates
    for pattern in file_patterns:
        ligand_file_pattern = os.path.join(ligand_directory, pattern)
        ligand_files.update(glob.glob(ligand_file_pattern))
    
    # Convert back to a sorted list for deterministic processing order
    ligand_files = sorted(list(ligand_files))
    
    if not ligand_files:
        patterns_str = ", ".join(file_patterns)
        print(f"No files matching '{patterns_str}' found in {ligand_directory}")
        return None
    
    # Create output directory for superimposed structures if needed
    if save_superimposed:
        if output_dir is None:
            output_dir = os.path.join(ligand_directory, "superimposed")
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    
    # Calculate RMSD for each ligand file
    results = []
    
    for i, ligand_file in enumerate(ligand_files):
        file_basename = os.path.basename(ligand_file)
        print(f"Processing {i+1}/{len(ligand_files)}: {file_basename}")
        
        # Determine output PDB path if saving superimposed structures
        output_pdb = None
        if save_superimposed:
            output_pdb_name = os.path.splitext(file_basename)[0] + "_superimposed.pdb"
            output_pdb = os.path.join(output_dir, output_pdb_name)
        
        try:
            # Calculate RMSD
            rmsd = calculate_docking_rmsd(
                pdb_result,
                ligand_file,
                pdb_complex,
                ligand_resname=ligand_resname,
                output_pdb=output_pdb
            )
            
            results.append({
                'file': file_basename,
                'full_path': ligand_file,
                'rmsd': rmsd
            })
            
            print(f"  RMSD: {rmsd:.3f} Å")
            
        except Exception as e:
            print(f"  Error processing {file_basename}: {str(e)}")
            results.append({
                'file': file_basename,
                'full_path': ligand_file,
                'rmsd': None,
                'error': str(e)
            })
    
    # Convert results to DataFrame
    results_df = pd.DataFrame(results)
    
    # Sort by RMSD if available
    if 'rmsd' in results_df.columns and not results_df['rmsd'].isnull().all():
        results_df = results_df.sort_values('rmsd')
    
    return results_df


#################################################
def batch_calculate_rmsd(
    pairs_df,
    prep_pdb_dir,
    ligand_base_dir,
    gt_pdb_dir,
    ligand_resname=None,
    file_pattern=None,
    output_base_dir=None,
    save_superimposed=False,
    save_individual_results=True
    ):
    """
    Calculate RMSD for multiple protein-ligand pairs specified in a DataFrame.
    
    Parameters:
    -----------
    pairs_df : pandas.DataFrame
        DataFrame with columns 'Protein' (PDB filename without extension) and 
        'Ligand' (ligand subfolder name)
    prep_pdb_dir : str
        Directory containing prepared PDB files for docking
    ligand_base_dir : str
        Base directory containing ligand result subfolders
    gt_pdb_dir : str
        Directory containing ground truth protein-ligand complex PDB files
    ligand_resname : str
        Residue name of the ligand in the PDB files (default None = autodetermined)
    file_pattern : str or list, optional
        Pattern(s) to match ligand files. If None, will use both "*.sdf" and "*.pdb"
    output_base_dir : str, optional
        Base directory for saving results. If None, results will be saved in ligand_base_dir
    save_superimposed : bool
        Whether to save superimposed structures as PDB files (default False)
    save_individual_results : bool
        Whether to save individual results for each protein-ligand pair (default True)
    
    Returns:
    --------
    pandas.DataFrame
        Combined DataFrame with RMSD results for all protein-ligand pairs
    """
    # Validate input DataFrame
    required_columns = ['Protein', 'Ligand']
    for col in required_columns:
        if col not in pairs_df.columns:
            raise ValueError(f"Input DataFrame must contain a '{col}' column")
    
    # Set up output directory
    if output_base_dir is None:
        output_base_dir = ligand_base_dir
    
    if not os.path.exists(output_base_dir):
        os.makedirs(output_base_dir)
    
    # List to store all results
    all_results = []
    
    # Process each protein-ligand pair
    for idx, row in tqdm(pairs_df.iterrows(), total=len(pairs_df), desc="Processing protein-ligand pairs"):
        protein_name = row['Protein']
        ligand_folder = row['Ligand']
        
        print(f"\nProcessing {protein_name} with ligand {ligand_folder} ({idx+1}/{len(pairs_df)})")
        
        # Construct file paths
        pdb_result = os.path.join(prep_pdb_dir, f"{protein_name}.pdb")
        ligand_directory = os.path.join(ligand_base_dir, ligand_folder)
        pdb_complex = os.path.join(gt_pdb_dir, f"{protein_name}.pdb")
        
        # Check if files exist
        if not os.path.exists(pdb_result):
            print(f"Warning: Prepared PDB file not found: {pdb_result}")
            continue
            
        if not os.path.exists(pdb_complex):
            print(f"Warning: Ground truth PDB complex not found: {pdb_complex}")
            continue
            
        if not os.path.exists(ligand_directory):
            print(f"Warning: Ligand directory not found: {ligand_directory}")
            continue
        
        # Set up file pattern based on folder name
        # If file_pattern is not provided, try to detect based on naming conventions
        local_file_pattern = file_pattern
        if local_file_pattern is None:
            # Try both old and new naming conventions
            # Old: pose_*.sdf, New: BCS-XXX_*_*.pdb
            bcs_id = None
            match = re.search(r'(BCS-[\w\d]+)', ligand_folder)
            if match:
                bcs_id = match.group(1)
                # Be more specific with patterns to reduce duplicates
                local_file_pattern = [f"{bcs_id}*.pdb", "*.sdf"] #[f"{bcs_id}_*_*.pdb", "pose_*.sdf"]
            else:
                local_file_pattern = ["*.sdf", "*.pdb"] #["pose_*.sdf", "*.pdb"]
        
        # Determine output directory for superimposed structures
        output_dir = None
        if save_superimposed:
            output_dir = os.path.join(output_base_dir, 'superimposed', protein_name, ligand_folder)
        
        try:
            # Calculate RMSD for this protein-ligand pair
            results_df = calculate_rmsd_for_directory(
                pdb_result,
                ligand_directory,
                pdb_complex,
                ligand_resname=ligand_resname,
                file_pattern=local_file_pattern,
                output_dir=output_dir,
                save_superimposed=save_superimposed
            )
            
            if results_df is not None and not results_df.empty:
                # Add protein and ligand information to results
                results_df['protein'] = protein_name
                results_df['ligand_folder'] = ligand_folder
                
                # Save individual results if requested
                if save_individual_results:
                    results_dir = os.path.join(output_base_dir, 'results')
                    if not os.path.exists(results_dir):
                        os.makedirs(results_dir)
                    results_csv = os.path.join(results_dir, f"{protein_name}_{ligand_folder}_rmsd.csv")
                    results_df.to_csv(results_csv, index=False)
                    print(f"Results saved to {results_csv}")
                
                # Add to combined results
                all_results.append(results_df)
            
        except Exception as e:
            print(f"Error processing {protein_name} with {ligand_folder}: {str(e)}")
            # You could add an error entry to all_results here if desired
    
    # Combine all results into a single DataFrame
    if all_results:
        combined_results = pd.concat(all_results, ignore_index=True)
        
        # Save combined results
        combined_csv = os.path.join(output_base_dir, 'combined_rmsd_results.csv')
        combined_results.to_csv(combined_csv, index=False)
        print(f"\nCombined results saved to {combined_csv}")
        
        return combined_results
    else:
        print("No results were generated")
        return pd.DataFrame()


#################################################
def print_combined_results(df=None):
    combined_results = df.copy()
    # Display summary of results
    if not combined_results.empty:
        # First, drop rows with NaN RMSD values
        combined_results_clean = combined_results.dropna(subset=['rmsd'])
        
        if combined_results_clean.empty:
            print("\nNo valid RMSD values found in the results. All calculations resulted in errors or NaN values.")
        else:
            # Group by protein and ligand, get mean and min RMSD
            summary = combined_results_clean.groupby(['protein', 'ligand_folder']).agg({
                'rmsd': ['mean', 'min', 'count']
            }).reset_index()
            
            print("\nSummary of Results:")
            print(summary.to_string())
            
            # Get the best pose for each protein-ligand pair - safely handling groups with all NaN values
            best_poses_indices = combined_results_clean.groupby(['protein', 'ligand_folder'])['rmsd'].idxmin()
            best_poses = combined_results_clean.loc[best_poses_indices]
            
            print("\nBest Poses for Each Protein-Ligand Pair:")
            print(best_poses[['protein', 'ligand_folder', 'file', 'rmsd']].to_string(index=False))
            
        # Also display information about failed calculations
        failed_calculations = combined_results[combined_results['rmsd'].isna()]
        if not failed_calculations.empty:
            print(f"\nNumber of failed RMSD calculations: {len(failed_calculations)}")
            print("\nFailed Calculations by Protein-Ligand Pair:")
            failed_summary = failed_calculations.groupby(['protein', 'ligand_folder']).size().reset_index(name='count')
            print(failed_summary.to_string(index=False))
            
            # If there are error messages in the DataFrame, display them
            if 'error' in failed_calculations.columns:
                print("\nSample of error messages:")
                for error in failed_calculations['error'].dropna().unique()[:5]:  # Show up to 5 unique error messages
                    print(f"- {error}")


def create_wide_format_rmsd_results(combined_results):
    """
    Convert combined RMSD results to a wide format with one column per pose.
    
    Parameters:
    -----------
    combined_results : pandas.DataFrame
        The combined RMSD results DataFrame with protein, ligand_folder, file, and rmsd columns
        
    Returns:
    --------
    pandas.DataFrame
        Wide-format DataFrame with Protein, Ligand and one column per pose RMSD
    """
    # Check if the combined_results DataFrame is empty
    if combined_results.empty:
        print("Combined results DataFrame is empty. No wide format can be created.")
        return pd.DataFrame()
    
    # Make a copy to avoid modifying the original DataFrame
    df = combined_results.copy()
    
    # Extract pose number from filename, handling both old and new patterns
    def extract_pose_number(filename):
        # First try old pattern: "pose_X.sdf"
        match = re.search(r'pose_(\d+)', filename)
        if match:
            return int(match.group(1))
        
        # Then try new pattern: "BCS-XXX_Y_Z.pdb"
        match = re.search(r'_(\d+)\.pdb$', filename)
        if match:
            return int(match.group(1))
            
        # If no match, try to extract any number from the filename
        match = re.search(r'(\d+)', filename)
        if match:
            return int(match.group(1))
        else:
            # If still no match, use the filename as is
            return filename
    
    # Add pose number column
    df['pose_number'] = df['file'].apply(extract_pose_number)
    
    # Create column names for each pose
    df['pose_column'] = df['pose_number'].apply(lambda x: f'RMSD_pose{x}')
    
    # Pivot the DataFrame
    # First check if we have the necessary columns
    required_cols = ['protein', 'ligand_folder', 'pose_column', 'rmsd']
    if not all(col in df.columns for col in required_cols):
        missing = [col for col in required_cols if col not in df.columns]
        print(f"Missing required columns: {missing}")
        print("Available columns:", df.columns.tolist())
        print("Creating a default wide format with available columns.")
        
        # Create a simple pivot with whatever columns we have
        if 'protein' in df.columns and 'ligand_folder' in df.columns and 'pose_column' in df.columns and 'rmsd' in df.columns:
            pivot_df = df.pivot(index=['protein', 'ligand_folder'], columns='pose_column', values='rmsd')
            pivot_df = pivot_df.reset_index()
        else:
            # Create a minimal pivot with whatever we have
            index_cols = [col for col in ['protein', 'ligand_folder'] if col in df.columns]
            if not index_cols:
                print("Cannot create pivot: missing both 'protein' and 'ligand_folder' columns")
                return df
            
            value_col = 'rmsd' if 'rmsd' in df.columns else df.columns[-1]
            col_col = 'pose_column' if 'pose_column' in df.columns else 'file'
            
            pivot_df = df.pivot(index=index_cols, columns=col_col, values=value_col)
            pivot_df = pivot_df.reset_index()
    else:
        # Create the pivot table with protein and ligand as indices, and pose columns as columns
        pivot_df = df.pivot(index=['protein', 'ligand_folder'], columns='pose_column', values='rmsd')
        pivot_df = pivot_df.reset_index()
    
    # Rename columns for clarity
    pivot_df = pivot_df.rename(columns={'protein': 'Protein', 'ligand_folder': 'Ligand'})
    
    # Fill NaN values with actual NaN rather than empty string
    pivot_df = pivot_df.replace('', np.nan)
    
    return pivot_df

def summarize_wide_format_results(wide_format_df):
    """
    Add summary statistics to the wide format RMSD results.
    
    Parameters:
    -----------
    wide_format_df : pandas.DataFrame
        The wide format RMSD results DataFrame
        
    Returns:
    --------
    pandas.DataFrame
        Wide-format DataFrame with added summary columns
    """
    if wide_format_df.empty:
        return wide_format_df
        
    # Make a copy to avoid modifying the original
    df = wide_format_df.copy()
    
    # Find RMSD columns
    rmsd_cols = [col for col in df.columns if col.startswith('RMSD_pose')]
    
    if rmsd_cols:
        # Add min, max, mean, and median RMSD columns
        df['Min_RMSD'] = df[rmsd_cols].min(axis=1)
        df['Max_RMSD'] = df[rmsd_cols].max(axis=1)
        df['Mean_RMSD'] = df[rmsd_cols].mean(axis=1)
        df['Median_RMSD'] = df[rmsd_cols].median(axis=1)
        
        # Find the best pose (minimum RMSD) for each protein-ligand pair
        def get_best_pose(row):
            rmsd_values = [row[col] for col in rmsd_cols]
            if all(pd.isna(val) for val in rmsd_values):
                return np.nan
            
            min_rmsd = min((val for val in rmsd_values if not pd.isna(val)), default=np.nan)
            for col in rmsd_cols:
                if row[col] == min_rmsd:
                    return col.replace('RMSD_', '')
            return np.nan
        
        df['Best_Pose'] = df.apply(get_best_pose, axis=1)
    
    return df




######################################### Example usage:
if __name__ == "__main__":
    # Example 1: Process a single ligand file
    print("Example 1: Single ligand file processing")
    
    # Input files
    pdb_result = "../../../data/interim/docking/prep_Flare/bcsin_235.pdb" 
    ligand_file = "../../../data/interim/docking/Res_Diffdock_Flare_R2/BCS-CN76566/pose_0.sdf"
    # For new format:
    # ligand_file = "../../../data/interim/docking/Res_Diffdock_Flare_R2/BCS-CN76566/BCS-CN76566_S_0.pdb"
    pdb_complex = "../../../data/raw/docking_old/GT/bcsin_235.pdb"
    
    # Calculate RMSD
    try:
        rmsd = calculate_docking_rmsd(
            pdb_result,
            ligand_file,
            pdb_complex,
            #ligand_resname="566",
            output_pdb=None
        )
        print(f"RMSD: {rmsd:.3f} Å")
    except Exception as e:
        print(f"Error: {str(e)}")
    
    # Example 2: Process a directory of ligand files
    print("\nExample 2: Directory processing")
    
    # Input paths
    pdb_result = "../../../data/interim/docking/prep_Flare/bcsin_235.pdb" 
    ligand_directory = "../../../data/interim/docking/Res_Diffdock_Flare_R2/BCS-CN76566"
    pdb_complex = "../../../data/raw/docking_old/GT/bcsin_235.pdb"
    
    # Calculate RMSD for all ligand files in the directory
    # Can handle both old "pose_*.sdf" and new "BCS-CN76566_*_*.pdb" patterns
    results_df = calculate_rmsd_for_directory(
        pdb_result,
        ligand_directory,
        pdb_complex,
        #ligand_resname="566",
        file_pattern=None,  # Will match both *.sdf and *.pdb
        save_superimposed=True
    )
    
    # Display results
    if results_df is not None and not results_df.empty:
        print("\nRMSD Results:")
        print(results_df[['file', 'rmsd']].to_string(index=False))
        
        # Optionally save results to CSV
        results_csv = os.path.join(ligand_directory, "rmsd_results.csv")
        results_df.to_csv(results_csv, index=False)
        print(f"\nResults saved to {results_csv}")
    
    # Example 3: Batch processing with a DataFrame
    print("\nExample 3: Batch processing with DataFrame")
    
    # Create a DataFrame with protein-ligand pairs
    pairs_data = [
        {'Protein': 'bcsin_235', 'Ligand': 'BCS-CN76566'},
        {'Protein': 'bcsin_236', 'Ligand': 'BCS-CN76567'},
        # Add more pairs as needed
    ]
    pairs_df = pd.DataFrame(pairs_data)
    
    # Directories
    prep_pdb_dir = "../../../data/interim/docking/prep_Flare"
    ligand_base_dir = "../../../data/interim/docking/Res_Diffdock_Flare_R2"
    gt_pdb_dir = "../../../data/raw/docking_old/GT"
    
    # Calculate RMSD for all pairs
    # Will automatically detect file patterns based on folder names
    combined_results = batch_calculate_rmsd(
        pairs_df,
        prep_pdb_dir,
        ligand_base_dir,
        gt_pdb_dir,
        #ligand_resname="566",
        file_pattern=None,  # Auto-detect both sdf and new pdb naming patterns
        save_superimposed=False
    )
    
    # Print summary of results
    if combined_results is not None and not combined_results.empty:
        print_combined_results(combined_results)
        
        # Create wide format results
        wide_format = create_wide_format_rmsd_results(combined_results)
        
        # Add summary statistics
        wide_format_with_summary = summarize_wide_format_results(wide_format)
        
        # Save wide format results
        wide_format_csv = os.path.join(ligand_base_dir, "wide_format_rmsd_results.csv")
        wide_format_with_summary.to_csv(wide_format_csv, index=False)
        print(f"\nWide format results saved to {wide_format_csv}")