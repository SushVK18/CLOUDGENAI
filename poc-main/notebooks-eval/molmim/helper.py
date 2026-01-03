from rdkit import Chem
from rdkit.Chem import AllChem, Draw, QED, Crippen
import pandas as pd
import requests
import json



# Function to calculate LogP
def calculate_logp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Crippen.MolLogP(mol)
    else:
        return None  # Return None if the SMILES is invalid
        
# Function to calculate QED
def calculate_qed(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return QED.qed(mol)
    else:
        return None  # Return None if the SMILES is invalid

        
# Function to calculate QED
def calculate_qeh(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return QED.qeh(mol)
    else:
        return None  # Return None if the SMILES is invalid



## molmim specific 
def visualize_molecules(molmim_response):
    """ Visualizing 2-D structures of MolMIM-generated compounds
    # Example usage:
    # result_image = visualize_molecules(molmim_response)
    """
    mols_from_gen_smis = []
    gen_legends = []
    
    for index, i in enumerate(molmim_response['generated']):
        # Convert SMILES to molecule and ensure uniqueness
        mol = Chem.MolFromSmiles(i['smiles'])
        if mol and mol not in mols_from_gen_smis:
            mols_from_gen_smis.append(mol)
            gen_legends.append(f"Generated Molecule {index + 1} \n Score: {i['score']:.3f}")

    #uniques
    unique_subset = list(set(mols_from_gen_smis))
    
    # Output the number of unique molecules
    print("Total unique molecule designs obtained:", len(unique_subset))

    # Draw the molecules in a grid
    return Draw.MolsToGridImage(unique_subset, legends=gen_legends, molsPerRow=5, subImgSize=(350, 350))


def df_from_molmim_response(json_data):
    """
    Create a DataFrame from the given JSON data.

    Parameters:
    json_data (dict): The JSON data containing the 'generated' key.

    Returns:
    pd.DataFrame: A DataFrame containing the 'smiles' and 'score' data.
    """
    # Check if 'generated' key exists in the JSON data
    if 'generated' in json_data:
        # Create a DataFrame from the 'generated' key
        df = pd.DataFrame(json_data['generated'])
        return df
    else:
        raise ValueError("The provided JSON data does not contain the 'generated' key.")

    