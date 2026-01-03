from Bio import PDB
import numpy as np

def analyze_pdb_structure(pdb_file, verbose=True):
    """
    Analyze a PDB file to identify proteins and ligands.
    
    Parameters:
    pdb_file (str): Path to the PDB file
    verbose (bool): Whether to print detailed information
    
    Returns:
    dict: A dictionary containing protein and ligand information and the structure object
    """
    # Parse PDB
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('PDB_structure', pdb_file)
    
    # Initialize dictionaries to store protein and ligand information
    proteins = {}
    ligands = {}
    waters = []
    other_molecules = []
    
    # Get header information if available
    header = structure.header if hasattr(structure, 'header') else {}
    
    # Iterate through the model, chains, and residues
    for model in structure:
        for chain in model:
            chain_id = chain.id
            if chain_id not in proteins:
                proteins[chain_id] = []
            
            for residue in chain:
                res_id = residue.id
                res_name = residue.get_resname()
                
                # Check if the residue is a standard amino acid (protein)
                if PDB.is_aa(residue):
                    proteins[chain_id].append(residue)
                
                # Check for water molecules
                elif res_name in ['HOH', 'WAT', 'H2O']:
                    waters.append(residue)
                
                # Check if the residue is a ligand (hetero atom that's not water)
                elif res_id[0] != ' ':  # Hetero atoms have a non-space character in the first position
                    # Get center of mass for the ligand
                    try:
                        coords = np.array([atom.get_coord() for atom in residue])
                        center_of_mass = coords.mean(axis=0)
                    except:
                        center_of_mass = None
                    
                    # Count atoms in the ligand
                    atom_count = len(list(residue.get_atoms()))
                    
                    # Store ligand information
                    if res_name not in ligands:
                        ligands[res_name] = []
                    
                    ligands[res_name].append({
                        'chain': chain_id,
                        'residue_id': res_id,
                        'atom_count': atom_count,
                        'center_of_mass': center_of_mass
                    })
                
                # Other molecules (e.g., modified residues)
                else:
                    other_molecules.append((chain_id, res_id, res_name))
    
    # Calculate total counts
    total_protein_residues = sum(len(chain_residues) for chain_residues in proteins.values())
    total_ligands = sum(len(ligand_instances) for ligand_instances in ligands.values())
    
    # Output the results if verbose
    if verbose:
        print(f"Structure: {pdb_file}")
        if header:
            try:
                print(f"Title: {header.get('name', 'N/A')}")
                print(f"Classification: {header.get('head', 'N/A')}")
                print(f"Deposition date: {header.get('deposition_date', 'N/A')}")
            except:
                pass
        
        print(f"\nProtein chains: {len(proteins)}")
        for chain_id, residues in proteins.items():
            if residues:  # Only show non-empty chains
                print(f" - Chain {chain_id}: {len(residues)} residues")
        
        print(f"\nTotal protein residues: {total_protein_residues}")
        print(f"Total ligands: {total_ligands}")
        print(f"Water molecules: {len(waters)}")
        
        if ligands:
            print("\nLigands found:")
            for ligand_name, instances in ligands.items():
                print(f" - {ligand_name} ({len(instances)} instances)")
                for i, instance in enumerate(instances):
                    print(f"   Instance {i+1}: Chain {instance['chain']}, Residue ID {instance['residue_id']}, Atoms: {instance['atom_count']}")
                    if instance['center_of_mass'] is not None:
                        com = instance['center_of_mass']
                        print(f"   Position: X={com[0]:.2f}, Y={com[1]:.2f}, Z={com[2]:.2f}")
        else:
            print("\nNo ligands found.")
        
        if other_molecules:
            print("\nOther molecules (modified residues):")
            for chain_id, res_id, res_name in other_molecules:
                print(f" - {res_name} (Chain: {chain_id}, Residue ID: {res_id})")
    
    # Return the results as a dictionary
    return {
        'structure': structure,  # Return the structure object so it can be used by other functions
        'proteins': proteins,
        'ligands': ligands,
        'waters': waters,
        'other_molecules': other_molecules,
        'total_protein_residues': total_protein_residues,
        'total_ligands': total_ligands,
        'header': header
    }

def identify_binding_sites(structure_info, distance_cutoff=5.0):
    """
    Identify potential binding sites based on proximity of protein residues to ligands.
    
    Parameters:
    structure_info (dict): Output from analyze_pdb_structure
    distance_cutoff (float): Maximum distance in Angstroms to consider a residue part of binding site
    
    Returns:
    dict: Dictionary of binding site residues for each ligand
    """
    binding_sites = {}
    
    structure = structure_info['structure']  # Get the structure object from the input dictionary
    proteins = structure_info['proteins']
    ligands = structure_info['ligands']
    
    for ligand_name, instances in ligands.items():
        binding_sites[ligand_name] = []
        
        for instance in instances:
            ligand_chain = instance['chain']
            ligand_res_id = instance['residue_id']
            
            # Get the actual ligand residue
            for model in structure:
                if ligand_chain in model:
                    chain = model[ligand_chain]
                    if ligand_res_id in chain:
                        ligand_residue = chain[ligand_res_id]
                        
                        # Find nearby protein residues
                        nearby_residues = []
                        
                        for chain_id, chain_residues in proteins.items():
                            for protein_residue in chain_residues:
                                # Check distance between ligand and protein residue
                                for ligand_atom in ligand_residue:
                                    for protein_atom in protein_residue:
                                        distance = ligand_atom - protein_atom
                                        if distance < distance_cutoff:
                                            nearby_residues.append({
                                                'chain': chain_id,
                                                'residue': protein_residue.get_resname(),
                                                'residue_id': protein_residue.id[1],
                                                'distance': distance
                                            })
                                            break  # One close atom is enough to consider the residue part of binding site
                        
                        # Add binding site info
                        binding_sites[ligand_name].append({
                            'ligand_instance': instance,
                            'binding_site_residues': nearby_residues
                        })
    
    # Print binding site information
    print("\nBinding Site Analysis:")
    for ligand_name, instances in binding_sites.items():
        print(f"\nBinding sites for {ligand_name}:")
        
        for i, instance in enumerate(instances):
            site_residues = instance['binding_site_residues']
            print(f" - Instance {i+1}: {len(site_residues)} residues in binding site")
            
            if site_residues:
                # Group by chain
                chains = {}
                for res in site_residues:
                    chain_id = res['chain']
                    if chain_id not in chains:
                        chains[chain_id] = []
                    chains[chain_id].append(res)
                
                # Print by chain
                for chain_id, residues in chains.items():
                    print(f"   Chain {chain_id}: {len(residues)} residues")
                    # Optionally print detailed residue information
                    # for res in residues:
                    #     print(f"     {res['residue']}{res['residue_id']} (distance: {res['distance']:.2f}Å)")
    
    return binding_sites

# Example usage
# result = analyze_pdb_structure('your_pdb_file.pdb')
# binding_sites = identify_binding_sites(result)

def integrate_with_py3dmol(pdb_file, highlight_binding_sites=True, binding_site_distance=5.0, print_analysis=True):
    """
    Analyze a PDB file and create a py3Dmol visualization with proteins as cartoons
    and ligands as sticks. Optionally highlight binding sites.
    
    Parameters:
    pdb_file (str): Path to the PDB file
    highlight_binding_sites (bool): Whether to highlight binding site residues
    binding_site_distance (float): Distance cutoff for binding site residues
    print_analysis (bool): Whether to print PDB analysis information
    
    Returns:
    py3Dmol.view: The visualization view object
    """
    import py3Dmol
    
    # Analyze the PDB structure
    result = analyze_pdb_structure(pdb_file, verbose=print_analysis)
    
    # Check if we have a valid structure
    if not result['total_protein_residues'] and not result['total_ligands']:
        print("WARNING: Neither protein residues nor ligands found in the PDB file!")
    
    # Create py3Dmol view
    view = py3Dmol.view(width=800, height=600)
    
    # Load the PDB file
    with open(pdb_file, 'r') as f:
        pdb_data = f.read()
    view.addModel(pdb_data, 'pdb')
    
    # First apply a base style to everything (thin lines)
    view.setStyle({}, {'line': {'colorscheme': 'whiteCarbon', 'opacity': 0.1}})
    
    # Apply protein cartoon style if protein residues exist
    if result['total_protein_residues'] > 0:
        view.setStyle({'polymer': True, 'hetflag': False}, 
                     {'cartoon': {'color': 'spectrum', 'opacity': 0.9}})
        
        # If there are multiple chains, color them differently
        if len(result['proteins']) > 1:
            for chain_id in result['proteins'].keys():
                if len(result['proteins'][chain_id]) > 0:  # Only if chain has residues
                    view.setStyle({'chain': chain_id}, 
                                 {'cartoon': {'color': chain_id}})
    
    # Apply ligand stick style if ligands exist
    if result['total_ligands'] > 0:
        # First set all ligands to sticks with default coloring
        view.setStyle({'hetflag': True, 'not': {'resn': ['HOH', 'WAT']}}, 
                     {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.3}})
        
        # If there are multiple ligand types, color them differently
        ligand_colors = {
            0: 'greenCarbon', 
            1: 'cyanCarbon', 
            2: 'magentaCarbon', 
            3: 'yellowCarbon',
            4: 'whiteCarbon',
            5: 'orangeCarbon'
        }
        
        # Apply different colors to different ligand types
        for i, ligand_name in enumerate(result['ligands'].keys()):
            color_scheme = ligand_colors.get(i % len(ligand_colors), 'greenCarbon')
            view.setStyle({'resn': ligand_name}, 
                         {'stick': {'colorscheme': color_scheme, 'radius': 0.3}})
    
    # Highlight binding sites if requested
    if highlight_binding_sites and result['total_ligands'] > 0:
        binding_sites = identify_binding_sites(result, distance_cutoff=binding_site_distance)
        
        # Get all binding site residue IDs
        binding_residues = set()
        
        for ligand_name, instances in binding_sites.items():
            for instance in instances:
                for res in instance['binding_site_residues']:
                    # Format: chain:residue_number
                    res_id = f"{res['chain']}:{res['residue_id']}"
                    binding_residues.add(res_id)
        
        # Highlight binding site residues
        if binding_residues:
            for res_id in binding_residues:
                chain, resnum = res_id.split(':')
                view.addStyle({'chain': chain, 'resi': resnum},
                             {'stick': {'colorscheme': 'yellowCarbon', 'radius': 0.2}})
    
    # Set view settings
    view.zoomTo()
    
    return view