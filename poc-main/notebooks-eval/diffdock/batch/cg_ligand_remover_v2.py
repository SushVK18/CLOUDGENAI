#!/usr/bin/env python3
"""
PDB Ligand Remover

"""

import os
import sys
import argparse
import re
from pathlib import Path


def is_numeric_resname(resname):
    """
    Check if a residue name is numeric (contains only digits).
    
    Args:
        resname (str): Residue name to check
    
    Returns:
        bool: True if the residue name is numeric, False otherwise
    """
    return resname.isdigit()


def remove_numeric_ligands_from_pdb(input_file, output_file, save_ligands=False, ligands_file=None, 
                                    verbose=False):
    """
    Remove HETATM entries with numeric residue names and their CONNECT entries from a PDB file.
    Optionally save removed ligands to a separate file.
    
    Args:
        input_file (str): Path to the input PDB file
        output_file (str): Path to the output PDB file
        save_ligands (bool): Whether to save removed ligands to a separate file
        ligands_file (str): Path to save the removed ligands (if None, generate from output_file)
        verbose (bool): Whether to print detailed information
    
    Returns:
        tuple: (num_hetatm_removed, num_connect_removed, removed_ligands)
    """
    # Store atom serial numbers of the ligands to remove
    ligand_serial_numbers = set()
    
    # Read the input file
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # Generate ligands file path if not provided
    if save_ligands and ligands_file is None:
        output_path = Path(output_file)
        ligands_file = str(output_path.parent / f"{output_path.stem}_ligands{output_path.suffix}")
    
    # First pass: identify the ligand atoms with numeric residue names
    removed_ligands = set()
    removed_lines = []  # Store lines that will be removed for saving later
    
    for line in lines:
        if line.startswith('HETATM'):
            resname = line[17:20].strip()
            
            if is_numeric_resname(resname):
                # Extract the atom serial number (columns 7-11)
                try:
                    serial_number = int(line[6:11].strip())
                    ligand_serial_numbers.add(serial_number)
                    removed_ligands.add(resname)
                    if save_ligands:
                        removed_lines.append(line)
                except ValueError:
                    if verbose:
                        print(f"Warning: Could not parse atom serial number in line: {line.strip()}")
    
    if verbose:
        print(f"Found {len(ligand_serial_numbers)} atoms for numeric ligands: {', '.join(removed_ligands)}")
    
    # Filter out ligand HETATM and CONNECT entries
    output_lines = []
    hetatm_removed = 0
    connect_removed = 0
    
    for line in lines:
        if line.startswith('HETATM'):
            resname = line[17:20].strip()
            
            if is_numeric_resname(resname):
                hetatm_removed += 1
                continue
        
        if line.startswith('CONECT'):
            # Parse the CONECT line
            fields = line.split()
            if len(fields) >= 2:
                try:
                    # First number in CONECT is the central atom
                    central_atom = int(fields[1])
                    # If the central atom is part of a ligand, skip this CONECT record
                    if central_atom in ligand_serial_numbers:
                        connect_removed += 1
                        if save_ligands:
                            removed_lines.append(line)
                        continue
                except ValueError:
                    pass
        
        output_lines.append(line)
    
    # Write the output file (with ligands removed)
    with open(output_file, 'w') as f:
        f.writelines(output_lines)
    
    # Write the ligands file if requested
    if save_ligands and removed_lines:
        # Add header and END record
        ligand_file_content = []
        
        # Add REMARK with information about the source file
        ligand_file_content.append(f"REMARK   1 LIGAND FILE EXTRACTED FROM {os.path.basename(input_file)}\n")
        ligand_file_content.append(f"REMARK   2 CONTAINS NUMERIC LIGANDS: {', '.join(removed_ligands)}\n")
        
        # Add all removed HETATM and CONECT lines
        ligand_file_content.extend(removed_lines)
        
        # Add END record
        ligand_file_content.append("END\n")
        
        with open(ligands_file, 'w') as f:
            f.writelines(ligand_file_content)
        
        if verbose:
            print(f"Saved {len(removed_lines)} removed lines to {ligands_file}")
    
    return hetatm_removed, connect_removed, removed_ligands


def remover_process_directory(input_dir, output_dir, save_ligands=False, ligands_dir=None,
                              verbose=False):
    """
    Process all PDB files in a directory.
    
    Args:
        input_dir (str): Path to the input directory
        output_dir (str): Path to the output directory
        save_ligands (bool): Whether to save removed ligands to separate files
        ligands_dir (str): Directory to save ligand files (if None, use output_dir)
        verbose (bool): Whether to print detailed information
    """
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    
    # Create output directory if it doesn't exist
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Create ligands directory if needed and doesn't exist
    if save_ligands:
        if ligands_dir is None:
            ligands_path = output_path
        else:
            ligands_path = Path(ligands_dir)
            ligands_path.mkdir(parents=True, exist_ok=True)
    
    # Get all PDB files in the input directory
    pdb_files = list(input_path.glob('*.pdb'))
    
    if not pdb_files:
        print(f"No PDB files found in {input_dir}")
        return
    
    total_files = len(pdb_files)
    total_hetatm_removed = 0
    total_connect_removed = 0
    all_removed_ligands = set()
    
    for i, pdb_file in enumerate(pdb_files, 1):
        output_file = output_path / pdb_file.name
        
        # Generate ligands file path if saving ligands
        if save_ligands:
            if ligands_dir is None:
                ligands_file = output_path / f"{pdb_file.stem}_ligands{pdb_file.suffix}"
            else:
                ligands_file = ligands_path / f"{pdb_file.stem}_ligands{pdb_file.suffix}"
        else:
            ligands_file = None
        
        if verbose:
            print(f"Processing file {i}/{total_files}: {pdb_file.name}")
        
        hetatm_removed, connect_removed, removed_ligands = remove_numeric_ligands_from_pdb(
            str(pdb_file), str(output_file), save_ligands, str(ligands_file) if ligands_file else None, verbose
        )
        
        total_hetatm_removed += hetatm_removed
        total_connect_removed += connect_removed
        all_removed_ligands.update(removed_ligands)
        
        if verbose:
            print(f"  Removed {hetatm_removed} HETATM and {connect_removed} CONNECT entries")
            if removed_ligands:
                print(f"  Ligands removed: {', '.join(removed_ligands)}")
    
    print(f"Processed {total_files} PDB files")
    print(f"Total: Removed {total_hetatm_removed} HETATM and {total_connect_removed} CONNECT entries")
    print(f"All ligands removed: {', '.join(all_removed_ligands)}")


def remove_specific_ligands_from_pdb(input_file, output_file, residue_names=None, numeric_only=True, 
                                     save_ligands=False, ligands_file=None, verbose=False):
    """
    Remove HETATM entries with specified residue names (or all numeric residue names if not specified)
    and their CONNECT entries from a PDB file. Optionally save removed ligands to a separate file.
    
    Args:
        input_file (str): Path to the input PDB file
        output_file (str): Path to the output PDB file
        residue_names (list or set, optional): Specific residue names to remove. If None and numeric_only=True,
                                              all numeric residue names will be removed.
        numeric_only (bool): Whether to only remove numeric residue names (when residue_names=None)
        save_ligands (bool): Whether to save removed ligands to a separate file
        ligands_file (str): Path to save the removed ligands (if None, generate from output_file)
        verbose (bool): Whether to print detailed information
    
    Returns:
        tuple: (num_hetatm_removed, num_connect_removed, removed_ligands)
    """
    from pathlib import Path
    import os
    
    # Helper function to check if a residue name is numeric
    def is_numeric_resname(resname):
        # Try to convert to integer or float
        try:
            int(resname)
            return True
        except ValueError:
            try:
                float(resname)
                return True
            except ValueError:
                return False
    
    # Store atom serial numbers of the ligands to remove
    ligand_serial_numbers = set()
    
    # Read the input file
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # Generate ligands file path if not provided
    if save_ligands and ligands_file is None:
        output_path = Path(output_file)
        ligands_file = str(output_path.parent / f"{output_path.stem}_ligands{output_path.suffix}")
    
    # First pass: identify the ligand atoms to remove
    removed_ligands = set()
    removed_lines = []  # Store lines that will be removed for saving later
    
    for line in lines:
        if line.startswith('HETATM'):
            resname = line[17:20].strip()
            
            # Check if the residue should be removed
            should_remove = False
            
            if residue_names is not None:
                # Remove specified residue names
                should_remove = resname in residue_names
            elif numeric_only:
                # Remove all numeric residue names
                should_remove = is_numeric_resname(resname)
            
            if should_remove:
                # Extract the atom serial number (columns 7-11)
                try:
                    serial_number = int(line[6:11].strip())
                    ligand_serial_numbers.add(serial_number)
                    removed_ligands.add(resname)
                    if save_ligands:
                        removed_lines.append(line)
                except ValueError:
                    if verbose:
                        print(f"Warning: Could not parse atom serial number in line: {line.strip()}")
    
    if verbose:
        print(f"Found {len(ligand_serial_numbers)} atoms for ligands: {', '.join(removed_ligands)}")
    
    # Filter out ligand HETATM and CONNECT entries
    output_lines = []
    hetatm_removed = 0
    connect_removed = 0
    
    for line in lines:
        if line.startswith('HETATM'):
            resname = line[17:20].strip()
            
            # Check if the residue should be removed
            should_remove = False
            
            if residue_names is not None:
                # Remove specified residue names
                should_remove = resname in residue_names
            elif numeric_only:
                # Remove all numeric residue names
                should_remove = is_numeric_resname(resname)
            
            if should_remove:
                hetatm_removed += 1
                continue
        
        if line.startswith('CONECT'):
            # Parse the CONECT line
            fields = line.split()
            if len(fields) >= 2:
                try:
                    # First number in CONECT is the central atom
                    central_atom = int(fields[1])
                    # If the central atom is part of a ligand, skip this CONECT record
                    if central_atom in ligand_serial_numbers:
                        connect_removed += 1
                        if save_ligands:
                            removed_lines.append(line)
                        continue
                except ValueError:
                    pass
        
        output_lines.append(line)
    
    # Write the output file (with ligands removed)
    with open(output_file, 'w') as f:
        f.writelines(output_lines)
    
    # Write the ligands file if requested
    if save_ligands and removed_lines:
        # Add header and END record
        ligand_file_content = []
        
        # Add REMARK with information about the source file
        ligand_file_content.append(f"REMARK   1 LIGAND FILE EXTRACTED FROM {os.path.basename(input_file)}\n")
        ligand_file_content.append(f"REMARK   2 CONTAINS REMOVED LIGANDS: {', '.join(removed_ligands)}\n")
        
        # Add all removed HETATM and CONECT lines
        ligand_file_content.extend(removed_lines)
        
        # Add END record
        ligand_file_content.append("END\n")
        
        with open(ligands_file, 'w') as f:
            f.writelines(ligand_file_content)
        
        if verbose:
            print(f"Saved {len(removed_lines)} removed lines to {ligands_file}")
    
    return hetatm_removed, connect_removed, removed_ligands

def main():
    parser = argparse.ArgumentParser(description="Remove ligands with numeric residue names from PDB files")
    parser.add_argument("input_path", help="Input PDB file or directory")
    parser.add_argument("output_path", help="Output PDB file or directory")
    parser.add_argument("--batch", action="store_true", help="Process all PDB files in input directory")
    parser.add_argument("--save-ligands", action="store_true", help="Save removed ligands to separate PDB files")
    parser.add_argument("--ligands-dir", help="Directory to save ligand files (default: same as output files)")
    parser.add_argument("--verbose", action="store_true", help="Print detailed information")
    
    args = parser.parse_args()
    
    if args.batch:
        # Process all PDB files in directory
        remover_process_directory(
            args.input_path, 
            args.output_path, 
            args.save_ligands, 
            args.ligands_dir, 
            args.verbose
        )
    else:
        # Process a single file
        try:
            ligands_file = None
            if args.save_ligands and args.ligands_dir:
                # Create ligands directory if it doesn't exist
                ligands_path = Path(args.ligands_dir)
                ligands_path.mkdir(parents=True, exist_ok=True)
                
                # Generate ligands file path
                input_filename = os.path.basename(args.input_path)
                ligands_filename = f"{os.path.splitext(input_filename)[0]}_ligands.pdb"
                ligands_file = str(ligands_path / ligands_filename)
            
            hetatm_removed, connect_removed, removed_ligands = remove_numeric_ligands_from_pdb(
                args.input_path, 
                args.output_path, 
                args.save_ligands,
                ligands_file,
                args.verbose
            )
            
            print(f"Removed {hetatm_removed} HETATM and {connect_removed} CONNECT entries")
            if removed_ligands:
                print(f"Ligands removed: {', '.join(removed_ligands)}")
                
            if args.save_ligands:
                if ligands_file:
                    print(f"Saved removed ligands to: {ligands_file}")
                else:
                    output_dir = os.path.dirname(args.output_path) or '.'
                    output_name = os.path.basename(args.output_path)
                    ligands_filename = f"{os.path.splitext(output_name)[0]}_ligands.pdb"
                    print(f"Saved removed ligands to: {os.path.join(output_dir, ligands_filename)}")
                    
        except Exception as e:
            print(f"Error processing file: {e}")
            sys.exit(1)
    
    

if __name__ == "__main__":
    main()