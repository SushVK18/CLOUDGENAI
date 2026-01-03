""" Capgemini Bioengineering Toolkit
#  
#   @@ All Rights Reserved @@
#   This file is part of the Capgemini Engineerings Bioengineering Toolkit
#   
#   Please read the accompanying LICENSE file
#
#   This module provides the core Molecule class that wraps RDKit functionality
#   with an easier-to-use interface and additional methods specific to the toolkit,
#   including fingerprint generation and similarity calculations.
#
#
#   not included in the toolkit are external dependencies e.g.:
#   pymican - MIT License https://github.com/ShintaroMinami/mican
#   rdkit - 3-Clause BSD License
#   etc.
#   Author: fabigr8
#
"""

from typing import Optional, Union, List, Dict, Any, Tuple
import os
import math
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski, MACCSkeys
from rdkit import DataStructs

class Molecule:
    """
    A class representing a small molecule (ligand, drug, etc.).
    
    This class wraps an RDKit molecule object and provides additional
    functionality for drug design applications.
    
    Attributes
    ----------
    name : str
        Name of the molecule
    """
    
    def __init__(
        self,
        input_data: Optional[Union[str, 'Chem.Mol']] = None,
        name: Optional[str] = None,
        sanitize: bool = True
    ):
        """
        Initialize a Molecule object.
        
        Parameters
        ----------
        input_data : str or rdkit.Chem.Mol, optional
            Input molecule data, which can be:
            - SMILES string
            - Path to a molecule file (.mol, .mol2, .sdf, .pdb)
            - RDKit Mol object
        name : str, optional
            Name of the molecule
        sanitize : bool, default=True
            Whether to sanitize the molecule during creation
        
        Examples
        --------
        >>> # From SMILES
        >>> mol = Molecule("CC(=O)OC1=CC=CC=C1C(=O)O")
        >>> # From file
        >>> mol = Molecule("path/to/molecule.sdf")
        >>> # From RDKit Mol
        >>> rdkit_mol = Chem.MolFromSmiles("CCO")
        >>> mol = Molecule(rdkit_mol)
        """
        self._rdkit_mol = None
        self._properties = {}
        self.name = name
        
        if input_data is not None:
            self._rdkit_mol = self._initialize_rdkit_mol(input_data, sanitize)
            
            # Set name if not provided but available in the RDKit molecule
            if name is None and self._rdkit_mol.HasProp("_Name"):
                self.name = self._rdkit_mol.GetProp("_Name")
        
        # If name is still None, set a default
        if self.name is None:
            self.name = "Unnamed Molecule"
    
    def _initialize_rdkit_mol(
        self, 
        input_data: Union[str, 'Chem.Mol'], 
        sanitize: bool
    ) -> 'Chem.Mol':
        """
        Initialize an RDKit Mol object from various input types.
        
        Parameters
        ----------
        input_data : str or rdkit.Chem.Mol
            Input data for the molecule
        sanitize : bool
            Whether to sanitize the molecule
            
        Returns
        -------
        rdkit.Chem.Mol
            RDKit Mol object
            
        Raises
        ------
        ValueError
            If the molecule cannot be created from the input data
        TypeError
            If the input data type is not supported
        """
        # If input is already an RDKit Mol
        if isinstance(input_data, Chem.Mol):
            return input_data
        
        # If input is a string, it could be SMILES or a file path
        elif isinstance(input_data, str):
            # Check if it's a file path
            if os.path.exists(input_data):
                # Determine file type by extension
                _, ext = os.path.splitext(input_data.lower())
                
                if ext in ['.mol', '.sdf', '.mdl']:
                    mol = Chem.MolFromMolFile(input_data, sanitize=sanitize)
                elif ext == '.mol2':
                    mol = Chem.MolFromMol2File(input_data, sanitize=sanitize)
                elif ext == '.pdb':
                    mol = Chem.MolFromPDBFile(input_data, sanitize=sanitize)
                else:
                    raise ValueError(f"Unsupported file format: {ext}")
            else:
                # Assume it's a SMILES string
                mol = Chem.MolFromSmiles(input_data, sanitize=sanitize)
            
            if mol is None:
                raise ValueError(f"Failed to create molecule from: {input_data}")
            
            return mol
        
        else:
            raise TypeError(f"Unsupported input type: {type(input_data)}")
    
    @property
    def mol(self) -> 'Chem.Mol':
        """Get the underlying RDKit Mol object."""
        return self._rdkit_mol
    
    @property
    def rdkit_mol(self) -> 'Chem.Mol':
        """Get the underlying RDKit Mol object."""
        return self._rdkit_mol
    
    @property
    def smiles(self) -> str:
        """Get the canonical SMILES representation of the molecule."""
        if self._rdkit_mol is None:
            return ""
        return Chem.MolToSmiles(self._rdkit_mol)
    
    def to_smiles(self) -> str:
        """Alias for the smiles property for compatibility."""
        return self.smiles
    
    @property
    def formula(self) -> str:
        """Get the molecular formula."""
        if self._rdkit_mol is None:
            return ""
        return Chem.rdMolDescriptors.CalcMolFormula(self._rdkit_mol)
    
    @property
    def molecular_weight(self) -> float:
        """Get the molecular weight in g/mol."""
        if self._rdkit_mol is None:
            return 0.0
        return Descriptors.MolWt(self._rdkit_mol)
    
    @property
    def num_atoms(self) -> int:
        """Get the number of atoms."""
        if self._rdkit_mol is None:
            return 0
        return self._rdkit_mol.GetNumAtoms()
    
    @property
    def num_bonds(self) -> int:
        """Get the number of bonds."""
        if self._rdkit_mol is None:
            return 0
        return self._rdkit_mol.GetNumBonds()
    
    @property
    def coordinates(self) -> np.ndarray:
        """
        Get the 3D coordinates of all atoms in the molecule.
        
        Returns
        -------
        numpy.ndarray
            Array of shape (n_atoms, 3) with the x, y, z coordinates
            
        Raises
        ------
        ValueError
            If the molecule has no 3D coordinates
        """
        if self._rdkit_mol is None:
            return np.array([])
        
        # Check if molecule has 3D coordinates
        if not self._rdkit_mol.GetNumConformers():
            raise ValueError("Molecule has no 3D coordinates")
        
        # Get the first conformer
        conf = self._rdkit_mol.GetConformer()
        
        # Extract coordinates
        coords = []
        for i in range(self._rdkit_mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            coords.append([pos.x, pos.y, pos.z])
        
        return np.array(coords)
    
    def generate_3d_coordinates(self, num_conformers: int = 1) -> bool:
        """
        Generate 3D coordinates for the molecule.
        
        Parameters
        ----------
        num_conformers : int, default=1
            Number of conformers to generate
            
        Returns
        -------
        bool
            True if coordinates were successfully generated
            
        Examples
        --------
        >>> mol = Molecule("CC(=O)OC1=CC=CC=C1C(=O)O")
        >>> success = mol.generate_3d_coordinates()
        """
        if self._rdkit_mol is None:
            return False
        
        try:
            # Add hydrogens
            mol_h = Chem.AddHs(self._rdkit_mol)
            
            # Generate conformers
            result = AllChem.EmbedMultipleConfs(
                mol_h, 
                numConfs=num_conformers,
                randomSeed=42,
                useExpTorsionAnglePrefs=True,
                useBasicKnowledge=True
            )
            
            # Optimize conformers
            for conf_id in range(mol_h.GetNumConformers()):
                AllChem.MMFFOptimizeMolecule(mol_h, confId=conf_id)
            
            # Remove hydrogens
            self._rdkit_mol = Chem.RemoveHs(mol_h)
            
            return self._rdkit_mol.GetNumConformers() > 0
        
        except Exception as e:
            print(f"Error generating 3D coordinates: {e}")
            return False
    
    def calculate_descriptors(self) -> Dict[str, float]:
        """
        Calculate common molecular descriptors.
        
        Returns
        -------
        dict
            Dictionary of descriptor name: value pairs
        
        Examples
        --------
        >>> mol = Molecule("CC(=O)OC1=CC=CC=C1C(=O)O")
        >>> descriptors = mol.calculate_descriptors()
        >>> print(f"LogP: {descriptors['LogP']}")
        """
        if self._rdkit_mol is None:
            return {}
        
        descriptors = {
            "MolecularWeight": Descriptors.MolWt(self._rdkit_mol),
            "LogP": Descriptors.MolLogP(self._rdkit_mol),
            "NumHDonors": Lipinski.NumHDonors(self._rdkit_mol),
            "NumHAcceptors": Lipinski.NumHAcceptors(self._rdkit_mol),
            "NumRotatableBonds": Descriptors.NumRotatableBonds(self._rdkit_mol),
            "TPSA": Descriptors.TPSA(self._rdkit_mol),
            "NumRings": Descriptors.RingCount(self._rdkit_mol),
            "NumAromaticRings": Descriptors.NumAromaticRings(self._rdkit_mol),
        }
        
        return descriptors
    
    def get_fingerprint(self, fp_type: str = 'morgan', radius: int = 2, nBits: int = 2048, as_numpy: bool = False) -> Union[np.ndarray, 'DataStructs.ExplicitBitVect', 'DataStructs.SparseBitVect']:
        """
        Generate molecular fingerprint.
        
        Parameters
        ----------
        fp_type : str, default='morgan'
            Type of fingerprint ('morgan', 'ecfp', 'fcfp', 'morgan_feat', 'maccs', 'topological')
        radius : int, default=2
            Radius for Morgan fingerprint variants
        nBits : int, default=2048
            Number of bits for fingerprint
        as_numpy : bool, default=False
            Whether to convert the fingerprint to a numpy array
        
        Returns
        -------
        Union[numpy.ndarray, rdkit.DataStructs object]
            Fingerprint in the requested format
        """
        if self._rdkit_mol is None:
            return np.array([]) if as_numpy else None
        
        # Generate fingerprint using your existing implementations
        if fp_type == 'morgan-old':
            fp = AllChem.GetMorganFingerprintAsBitVect(self._rdkit_mol, radius, nBits=nBits)
        elif fp_type == 'morgan':
            fpgen = AllChem.GetMorganGenerator(radius=radius)
            fp = fpgen.GetFingerprint(self._rdkit_mol)
        elif fp_type == 'morgan_feat':
            fpgen = AllChem.GetMorganGenerator(radius=radius)
            fp = fpgen.GetSparseCountFingerprint(self._rdkit_mol)
        elif fp_type == 'ecfp':
            invgen = AllChem.GetMorganAtomInvGen() 
            fpgen = AllChem.GetMorganGenerator(radius=radius, atomInvariantsGenerator=invgen)
            fp = fpgen.GetSparseCountFingerprint(self._rdkit_mol)
        elif fp_type == 'fcfp':
            invgen = AllChem.GetMorganFeatureAtomInvGen() 
            fpgen = AllChem.GetMorganGenerator(radius=radius, atomInvariantsGenerator=invgen)
            fp = fpgen.GetSparseCountFingerprint(self._rdkit_mol)
        elif fp_type == 'maccs':
            fp = MACCSkeys.GenMACCSKeys(self._rdkit_mol)
        elif fp_type == 'topological':
            fpgen = AllChem.GetRDKitFPGenerator()
            fp = Chem.GetFingerprint(self._rdkit_mol)
        else:
            raise ValueError(f"Unsupported fingerprint type: {fp_type}")
        
        # Handle output format
        if as_numpy:
            # Try to convert to numpy array based on the type of fingerprint
            try:
                # For bit vectors and explicit bit vectors
                arr = np.zeros((1,))
                DataStructs.ConvertToNumpyArray(fp, arr)
                return arr
            except Exception as e:
                # For sparse count fingerprints or other types
                if hasattr(fp, "GetNonzeroElements"):
                    # Convert sparse fingerprint to numpy array
                    sparse_dict = fp.GetNonzeroElements()
                    max_idx = max(sparse_dict.keys()) if sparse_dict else 0
                    arr = np.zeros(max_idx + 1)
                    for idx, count in sparse_dict.items():
                        arr[idx] = count
                    return arr
                else:
                    raise ValueError(f"Cannot convert fingerprint of type {type(fp)} to numpy array: {e}")
        
        # Return the RDKit fingerprint object if numpy conversion not requested
        return fp
    
    @staticmethod
    def tanimoto_similarity(fp1: np.ndarray, fp2: np.ndarray) -> float: # FG Not needed? 
        """
        Calculate Tanimoto similarity between two fingerprints.
        
        Parameters
        ----------
        fp1 : numpy.ndarray
            First fingerprint as numpy array
        fp2 : numpy.ndarray
            Second fingerprint as numpy array
        
        Returns
        -------
        float
            Tanimoto similarity score (0-1)
        """
        
        #intersection = np.sum(np.logical_and(fp1, fp2))
        #union = np.sum(np.logical_or(fp1, fp2))
        #return float(intersection) / union if union > 0 else 0.0
        return DataStructs.TanimotoSimilarity(fp1,fp2)
    
    def calculate_similarity(self, other_mol: 'Molecule', fp_type: str = 'morgan', sim_type: str = 'tanimoto') -> float:
        """
        Calculate similarity between this molecule and another.
        
        Parameters
        ----------
        other_mol : Molecule
            Another Molecule object
        fp_type : str, default='morgan'
            Type of fingerprint to use
        
        Returns
        -------
        float
            Similarity score (0-1)
            
        Examples
        --------
        >>> mol1 = Molecule("CC(=O)OC1=CC=CC=C1C(=O)O")
        >>> mol2 = Molecule("CC(=O)OCCC(=O)O")
        >>> similarity = mol1.calculate_similarity(mol2)
        """
        
        fp1 = self.get_fingerprint(fp_type=fp_type)
        fp2 = other_mol.get_fingerprint(fp_type=fp_type)
        
        if sim_type == 'tanimoto':
            return DataStructs.TanimotoSimilarity(fp1,fp2)
        elif sim_type == 'dice':
            return DataStructs.DiceSimilarity(fp1, fp2)
        elif sim_type == 'cosine':
            return DataStructs.CosineSimilarity(fp1, fp2)
        elif sim_type == 'sokal':
            return DataStructs.SokalSimilarity(fp1, fp2)
        elif sim_type == 'russel':  
            return DataStructs.RusselSimilarity(fp1, fp2)
        else:
            raise ValueError(f"Unsupported similarity calc type: {sim_type}")
        
    
    def __str__(self) -> str:
        """String representation of the molecule."""
        if self._rdkit_mol is None:
            return f"{self.name} (Empty)"
        
        return f"{self.name}: {self.smiles} ({self.formula}, MW={self.molecular_weight:.2f})"
    
    def __repr__(self) -> str:
        """Detailed string representation of the molecule."""
        return f"Molecule(name='{self.name}', smiles='{self.smiles}')"

# FG: Added MoleculeSet class
class MoleculeSet:
    """
    Container class for working with sets of molecules.
    Provides batch operations and dataset-level metrics.
    
    Attributes
    ----------
    molecules : list
        List of Molecule objects
    """
    
    def __init__(self, molecules: Optional[List[Molecule]] = None):
        """
        Initialize a set of molecules.
        
        Parameters
        ----------
        molecules : list of Molecule, optional
            List of Molecule objects
            
        Examples
        --------
        >>> mol_set = MoleculeSet()
        >>> mol_set.add_molecule(Molecule("CC(=O)OC1=CC=CC=C1C(=O)O", name="Aspirin"))
        >>> mol_set.add_molecule(Molecule("CC1=CC=C(C=C1)C(=O)OC", name="Methyl p-toluate"))
        """
        self.molecules = molecules or []
    
    def add_molecule(self, molecule: Molecule) -> None:
        """
        Add a molecule to the set.
        
        Parameters
        ----------
        molecule : Molecule
            Molecule to add to the set
        """
        self.molecules.append(molecule)
    
    def add_from_smiles_list(self, smiles_list: List[str], names: Optional[List[str]] = None) -> None:
        """
        Add multiple molecules from a list of SMILES strings.
        
        Parameters
        ----------
        smiles_list : list of str
            List of SMILES strings
        names : list of str, optional
            List of names for the molecules (must match length of smiles_list)
            
        Examples
        --------
        >>> mol_set = MoleculeSet()
        >>> smiles_list = ["CC(=O)OC1=CC=CC=C1C(=O)O", "CC1=CC=C(C=C1)C(=O)OC"]
        >>> names = ["Aspirin", "Methyl p-toluate"]
        >>> mol_set.add_from_smiles_list(smiles_list, names)
        """
        if names and len(names) != len(smiles_list):
            raise ValueError("If provided, names list must have same length as smiles_list")
            
        for i, smiles in enumerate(smiles_list):
            try:
                name = names[i] if names else None
                mol = Molecule(smiles, name=name)
                self.molecules.append(mol)
            except ValueError:
                print(f"Warning: Failed to create molecule from SMILES: {smiles}")
    
    def get_fingerprints(self, fp_type: str = 'morgan', **fp_kwargs) -> np.ndarray:
        """
        Get fingerprints for all molecules in the set.
        
        Parameters
        ----------
        fp_type : str, default='morgan'
            Type of fingerprint to generate
        **fp_kwargs : 
            Additional keyword arguments to pass to get_fingerprint method
        
        Returns
        -------
        list
            List of fingerprints in the appropriate format
        """
        if not self.molecules:
            return []
        
        # Get fingerprints for all molecules
        fps = []
        for mol in self.molecules:
            fps.append(mol.get_fingerprint(fp_type=fp_type, **fp_kwargs))
        
        return fps

    def calculate_pairwise_similarity(self, fp_type: str = 'morgan', sim_type: str = 'tanimoto') -> np.ndarray:
        """
        Calculate pairwise similarity matrix for all molecules in the set.
        
        Parameters
        ----------
        fp_type : str, default='morgan'
            Type of fingerprint to use
        sim_type : str, default='tanimoto'
            Type of similarity metric to use ('tanimoto', 'dice')
        
        Returns
        -------
        numpy.ndarray
            Square matrix of pairwise similarities
        """
        n = len(self.molecules)
        similarity_matrix = np.zeros((n, n))
        
        # Get fingerprints - keeping them in native format for faster similarity calculations
        fingerprints = self.get_fingerprints(fp_type=fp_type)
        
        for i in range(n):
            for j in range(i, n):
                if sim_type == 'tanimoto':
                    sim = DataStructs.TanimotoSimilarity(fingerprints[i], fingerprints[j])
                elif sim_type == 'dice':
                    sim = DataStructs.DiceSimilarity(fingerprints[i], fingerprints[j])
                else:
                    raise ValueError(f"Unsupported similarity type: {sim_type}")
                
                similarity_matrix[i, j] = sim
                similarity_matrix[j, i] = sim
        
        return similarity_matrix
    
    def calculate_shannon_entropy(self, fp_type: str = 'morgan') -> float:
        """
        Calculate Shannon entropy as a measure of molecular diversity.
        
        Parameters
        ----------
        fp_type : str, default='morgan'
            Type of fingerprint to use
        
        Returns
        -------
        float
            Shannon entropy value
        """
        if len(self.molecules) <= 1:
            return 0.0
        
        # Get pairwise similarities
        sim_matrix = self.calculate_pairwise_similarity(fp_type=fp_type)
        
        # Convert similarities to probabilities (normalize)
        sim_sum = np.sum(sim_matrix)
        if sim_sum == 0:
            return 0.0
        
        prob_matrix = sim_matrix / sim_sum
        
        # Calculate entropy
        entropy = 0.0
        for prob in prob_matrix.flatten():
            if prob > 0:
                entropy -= prob * math.log2(prob)
        
        return entropy
    
    def calculate_mean_diversity(self, fp_type: str = 'morgan') -> float:
        """
        Calculate mean diversity (1 - mean similarity) of the molecule set.
        
        Parameters
        ----------
        fp_type : str, default='morgan'
            Type of fingerprint to use
        
        Returns
        -------
        float
            Mean diversity score (0-1)
            
        Examples
        --------
        >>> mol_set = MoleculeSet()
        >>> mol_set.add_from_smiles_list(["CC(=O)OC1=CC=CC=C1C(=O)O", "CC1=CC=C(C=C1)C(=O)OC", "CCO"])
        >>> diversity = mol_set.calculate_mean_diversity()
        """
        if len(self.molecules) <= 1:
            return 0.0
        
        sim_matrix = self.calculate_pairwise_similarity(fp_type=fp_type)
        
        # Exclude self-similarities (diagonal)
        n = len(self.molecules)
        total_sim = sim_matrix.sum() - np.trace(sim_matrix)
        num_pairs = n * (n - 1)
        
        mean_sim = total_sim / num_pairs if num_pairs > 0 else 0.0
        return 1.0 - mean_sim
    
    def __len__(self) -> int:
        """Return the number of molecules in the set."""
        return len(self.molecules)
    
    def __getitem__(self, index: int) -> Molecule:
        """Get a molecule by index."""
        return self.molecules[index]
    
    def __iter__(self):
        """Iterator over molecules in the set."""
        return iter(self.molecules)
    
    def __str__(self) -> str:
        """String representation of the molecule set."""
        return f"MoleculeSet with {len(self.molecules)} molecules"