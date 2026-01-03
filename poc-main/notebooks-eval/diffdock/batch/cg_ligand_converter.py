#!/usr/bin/env python3
"""
PDB to SDF Converter

This script converts PDB files to SDF format, preserving the connectivity information
from CONECT records. It can process a single file or all PDB files in a directory.

"""

import os
import sys
import argparse
from pathlib import Path
import pandas as pd
import warnings
warnings.filterwarnings("ignore")


# Try to import the required libraries
try:
    import rdkit
    from rdkit import Chem
    from rdkit.Chem import AllChem
    HAVE_RDKIT = True
except ImportError:
    HAVE_RDKIT = False
    print("Warning: RDKit not found. RDKit-based conversion will not be available.")

try:
    from Bio.PDB import PDBParser
    from Bio.PDB.PDBIO import PDBIO
    import tempfile
    HAVE_BIOPYTHON = True
except ImportError:
    HAVE_BIOPYTHON = False
    print("Warning: Biopython not found. Biopython-based conversion will not be available.")


def convert_pdb_to_sdf_rdkit(pdb_file, sdf_file, ligand_only=False, verbose=False):
    """
    Convert PDB file to SDF using RDKit.
    
    Args:
        pdb_file (str): Path to the input PDB file
        sdf_file (str): Path to the output SDF file
        ligand_only (bool): If True, extract only ligands (HETATM)
        verbose (bool): If True, print detailed progress info
    
    Returns:
        bool: True if successful, False otherwise
    """
    if not HAVE_RDKIT:
        print("Error: RDKit is required for this conversion method.")
        return False
    
    if verbose:
        print(f"Reading PDB file: {pdb_file}")
    
    if ligand_only:
        # Read only HETATM entries (ligands)
        mol = Chem.MolFromPDBFile(pdb_file, removeHs=False, sanitize=False, proximityBonding=False)
        
        # If conversion failed and returned None, try alternative approach
        if mol is None:
            # Manual approach: read file, extract HETATM lines, write temp file, then read with RDKit
            try:
                with open(pdb_file, 'r') as f:
                    lines = f.readlines()
                
                # Keep only HETATM lines and CONECT records
                hetatm_lines = [line for line in lines if line.startswith('HETATM') or line.startswith('CONECT')]
                
                if not hetatm_lines:
                    if verbose:
                        print(f"No HETATM entries found in {pdb_file}")
                    return False
                
                # Write to temporary file
                import tempfile
                with tempfile.NamedTemporaryFile(mode='w+', suffix='.pdb', delete=False) as temp_file:
                    temp_file_path = temp_file.name
                    temp_file.writelines(hetatm_lines)
                    temp_file.write("END\n")  # Add END record
                
                # Try to read the temporary file
                mol = Chem.MolFromPDBFile(temp_file_path, removeHs=False, sanitize=False, proximityBonding=False)
                
                # Clean up
                try:
                    os.unlink(temp_file_path)
                except:
                    pass
            except Exception as e:
                if verbose:
                    print(f"Error in alternative ligand extraction: {e}")
                return False
    else:
        # Read all atoms including protein
        reader = Chem.PDBFileReader()
        reader.SetReadAllModels(False)  # Only read the first model
        reader.SetSanitize(False)  # Don't sanitize to preserve original data
        reader.SetRemoveHs(False)  # Don't remove hydrogens
        reader.SetProximityBonding(False)  # Don't add bonds based on proximity
        mol = reader.ReadFile(pdb_file)
    
    if mol is None:
        print(f"Error: Failed to read PDB file: {pdb_file}")
        return False
    
    # Check if we have a valid molecule
    if mol is None:
        print(f"Error: Failed to create a valid molecule from {pdb_file}")
        return False
        
    # Try to use the connectivity information from the PDB
    try:
        mol.GetConformer()
    except:
        # If there's an issue with conformers, try to compute 3D coordinates
        try:
            AllChem.EmbedMolecule(mol)
        except:
            if verbose:
                print("Warning: Could not generate 3D coordinates. Some information might be lost.")
    
    # Make sure the molecule has a name
    try:
        if mol.GetProp("_Name") == "":
            mol_name = os.path.splitext(os.path.basename(pdb_file))[0]
            mol.SetProp("_Name", mol_name)
    except KeyError:
        # If the _Name property doesn't exist, set it
        mol_name = os.path.splitext(os.path.basename(pdb_file))[0]
        mol.SetProp("_Name", mol_name)
    
    # Write to SDF
    if verbose:
        print(f"Writing SDF file: {sdf_file}")
    
    writer = Chem.SDWriter(sdf_file)
    writer.write(mol)
    writer.close()
    
    if verbose:
        print(f"Conversion completed successfully.")
    
    return True


def convert_pdb_to_sdf_biopython(pdb_file, sdf_file, ligand_only=False, verbose=False):
    """
    Convert PDB file to SDF using Biopython and Open Babel.
    
    Args:
        pdb_file (str): Path to the input PDB file
        sdf_file (str): Path to the output SDF file
        ligand_only (bool): If True, extract only ligands (HETATM)
        verbose (bool): If True, print detailed progress info
    
    Returns:
        bool: True if successful, False otherwise
    """
    if not HAVE_BIOPYTHON:
        print("Error: Biopython is required for this conversion method.")
        return False
    
    try:
        from openbabel import openbabel
        HAVE_OPENBABEL = True
    except ImportError:
        print("Error: OpenBabel Python bindings are required for Biopython conversion.")
        return False
    
    if verbose:
        print(f"Reading PDB file: {pdb_file}")
    
    # Parse the PDB file
    parser = PDBParser(QUIET=not verbose)
    try:
        structure = parser.get_structure("structure", pdb_file)
    except Exception as e:
        print(f"Error parsing PDB file: {e}")
        return False
    
    # If ligand_only, create a new structure with only HETATM
    if ligand_only:
        from Bio.PDB.Structure import Structure
        from Bio.PDB.Model import Model
        
        # Create a new structure to hold only ligands
        ligand_structure = Structure("ligands")
        ligand_model = Model(0)
        ligand_structure.add(ligand_model)
        
        # Add all hetero residues to the new structure
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.id[0] != ' ':  # Non-blank insertion code means hetero
                        ligand_model.add(chain)
                        break
        
        structure = ligand_structure
    
    # Create a temporary PDB file
    with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
        tmp_pdb = tmp.name
    
    # Write the structure to the temporary PDB file
    io = PDBIO()
    io.set_structure(structure)
    io.save(tmp_pdb)
    
    # Convert to SDF using OpenBabel
    if verbose:
        print(f"Converting to SDF using OpenBabel")
    
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdb", "sdf")
    
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, tmp_pdb)
    
    # Write to SDF
    obConversion.WriteFile(mol, sdf_file)
    
    # Clean up the temporary file
    try:
        os.unlink(tmp_pdb)
    except:
        pass
    
    if verbose:
        print(f"Conversion completed successfully.")
    
    return True


def process_file(input_file, output_file, use_rdkit=True, ligand_only=False, verbose=False):
    """
    Convert a single PDB file to SDF.
    
    Args:
        input_file (str): Path to the input PDB file
        output_file (str): Path to the output SDF file
        use_rdkit (bool): If True, use RDKit for conversion, otherwise use Biopython
        ligand_only (bool): If True, extract only ligands (HETATM)
        verbose (bool): If True, print detailed progress info
    
    Returns:
        bool: True if successful, False otherwise
    """
    if use_rdkit:
        if not HAVE_RDKIT:
            print("Error: RDKit is required but not installed.")
            return False
        return convert_pdb_to_sdf_rdkit(input_file, output_file, ligand_only, verbose)
    else:
        if not HAVE_BIOPYTHON:
            print("Error: Biopython is required but not installed.")
            return False
        return convert_pdb_to_sdf_biopython(input_file, output_file, ligand_only, verbose)


def converter_process_directory(input_dir, output_dir, use_rdkit=True, ligand_only=False, verbose=False):
    """
    Convert all PDB files in a directory to SDF.
    
    Args:
        input_dir (str): Path to the input directory
        output_dir (str): Path to the output directory
        use_rdkit (bool): If True, use RDKit for conversion, otherwise use Biopython
        ligand_only (bool): If True, extract only ligands (HETATM)
        verbose (bool): If True, print detailed progress info
    """
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    
    # Create output directory if it doesn't exist
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Get all PDB files in the input directory
    pdb_files = list(input_path.glob('*.pdb'))
    
    if not pdb_files:
        print(f"No PDB files found in {input_dir}")
        return
    
    total_files = len(pdb_files)
    successful = 0
    
    for i, pdb_file in enumerate(pdb_files, 1):
        output_file = output_path / f"{pdb_file.stem}.sdf"
        
        if verbose:
            print(f"Processing file {i}/{total_files}: {pdb_file.name}")
        
        if process_file(str(pdb_file), str(output_file), use_rdkit, ligand_only, verbose):
            successful += 1
    
    print(f"Processed {total_files} PDB files. {successful} conversions successful.")




def rename_sdf_files(directory, dataframe, protein_column, ligand_column, suffix="_ligands.sdf"):
    """
    Rename SDF files based on a mapping between protein and ligand names in a dataframe.
    
    Parameters:
    -----------
    directory : str
        Path to the directory containing SDF files
    dataframe : pandas.DataFrame
        DataFrame containing the mapping information
    protein_column : str
        Name of the column containing protein identifiers
    ligand_column : str
        Name of the column containing ligand identifiers
    suffix : str, optional
        Suffix to remove from filenames, default is "_ligand.sdf"
    """
    # Create a direct mapping from Protein names to Ligand names
    protein_to_ligand = {
        row[protein_column]: row[ligand_column] 
        for _, row in dataframe.iterrows()
    }
    
    # Get list of all files to check existence before renaming
    all_files = os.listdir(directory)
    renamed_count = 0
    skipped_count = 0
    
    # Iterate through files in the directory
    for filename in all_files:
        if filename.endswith('.sdf'):
            # Extract the Protein name part (without the suffix)
            file_base = filename.replace(suffix, '')
            print(f"Processing file: {filename}, extracted base: {file_base}")
            
            # Check if this protein is in our mapping
            if file_base in protein_to_ligand:
                # Get the new name
                new_name = protein_to_ligand[file_base] + '.sdf'
                
                # Full paths
                old_path = os.path.join(directory, filename)
                new_path = os.path.join(directory, new_name)
                
                # Check if destination file already exists
                if os.path.exists(new_path) and old_path != new_path:
                    print(f"Skipped: {filename} → {new_name} (destination file already exists)")
                    skipped_count += 1
                else:
                    # Rename the file
                    os.rename(old_path, new_path)
                    print(f"Renamed: {filename} → {new_name}")
                    renamed_count += 1
    
    print(f"Summary: Renamed {renamed_count} files, skipped {skipped_count} files")
    
    return renamed_count, skipped_count


####################### MAIN
def main():
    parser = argparse.ArgumentParser(description="Convert PDB files to SDF format")
    parser.add_argument("input_path", help="Input PDB file or directory")
    parser.add_argument("output_path", help="Output SDF file or directory")
    parser.add_argument("--batch", action="store_true", help="Process all PDB files in input directory")
    parser.add_argument("--ligand-only", action="store_true", help="Extract and convert only ligands (HETATM)")
    
    # Conversion method options (mutually exclusive)
    method_group = parser.add_mutually_exclusive_group()
    method_group.add_argument("--rdkit", action="store_true", help="Use RDKit for conversion (default)")
    method_group.add_argument("--biopython", action="store_true", help="Use Biopython for conversion")
    
    parser.add_argument("--verbose", action="store_true", help="Print detailed information")
    
    args = parser.parse_args()
    
    # Determine which method to use
    use_rdkit = not args.biopython  # Default to RDKit unless Biopython is explicitly requested
    
    # Check if required libraries are available
    if use_rdkit and not HAVE_RDKIT:
        print("Error: RDKit is required but not installed.")
        print("Please install RDKit or use --biopython option.")
        return 1
    
    if not use_rdkit and not HAVE_BIOPYTHON:
        print("Error: Biopython is required but not installed.")
        print("Please install Biopython or use --rdkit option.")
        return 1
    
    if args.batch:
        # Process all PDB files in directory
        converter_process_directory(
            args.input_path, 
            args.output_path, 
            use_rdkit, 
            args.ligand_only, 
            args.verbose
        )
    else:
        # Process a single file
        try:
            success = process_file(
                args.input_path, 
                args.output_path, 
                use_rdkit, 
                args.ligand_only, 
                args.verbose
            )
            
            if success:
                print(f"Conversion successful: {args.output_path}")
            else:
                print("Conversion failed.")
                return 1
                
        except Exception as e:
            print(f"Error processing file: {e}")
            return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())