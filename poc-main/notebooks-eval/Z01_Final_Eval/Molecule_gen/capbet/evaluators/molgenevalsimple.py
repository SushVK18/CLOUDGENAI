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

import json
import os
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from typing import List, Dict, Tuple, Optional, Union, Set
from rdkit import Chem, DataStructs
from datetime import datetime
from typing import Optional, Union, List, Dict, Any, Tuple
from rdkit.Chem import AllChem, Descriptors, Lipinski, MACCSkeys

# import Moleculeset and Molecule classes from core.wrappers.molecule
from capbet.core.wrappers.molecule import Molecule, MoleculeSet

# You already have Molecule and MoleculeSet classes, so we'll reference those

class MoleculeSetEvaluatorSimple:
    """
    A class to evaluate metrics for molecular datasets.
    
    This evaluator calculates various metrics including:
    - Validity: Percentage of valid SMILES strings
    - Uniqueness: Percentage of unique valid molecules
    - Novelty: Percentage of molecules not in a reference set
    - Non-Identicality: Percentage of molecules not in the input/seed set
    - Effective Novelty: Percentage of molecules that are both novel and non-identical
    
    Attributes
    ----------
    generated_set : MoleculeSet
        Set of generated molecules to evaluate
    reference_set : MoleculeSet
        Reference set for novelty comparison
    input_set : MoleculeSet, optional
        Input/seed set for non-identicality comparison
    comparison_method : str
        Method for molecule comparison ('smiles' or 'fingerprint')
    fingerprint_type : str
        Type of fingerprint to use if comparison_method is 'fingerprint'
    similarity_threshold : float
        Threshold for considering molecules similar with fingerprint comparison
    """
    
    def __init__(
        self,
        generated_set: 'MoleculeSet',
        reference_set: 'MoleculeSet',
        input_set: Optional['MoleculeSet'] = None,
        comparison_method: str = 'smiles',
        fingerprint_type: str = 'morgan',
        similarity_threshold: float = 0.85
    ):
        """
        Initialize the MoleculeSetEvaluator.
        
        Parameters
        ----------
        generated_set : MoleculeSet
            Set of generated molecules to evaluate
        reference_set : MoleculeSet
            Reference set for novelty comparison
        input_set : MoleculeSet, optional
            Input/seed set for non-identicality comparison
        comparison_method : str, default='smiles'
            Method for molecule comparison ('smiles' or 'fingerprint')
        fingerprint_type : str, default='morgan'
            Type of fingerprint to use if comparison_method is 'fingerprint'
        similarity_threshold : float, default=0.85
            Threshold for considering molecules similar with fingerprint comparison
        """
        self.generated_set = generated_set
        self.reference_set = reference_set
        self.input_set = input_set if input_set is not None else None
        self.comparison_method = comparison_method
        self.fingerprint_type = fingerprint_type
        self.similarity_threshold = similarity_threshold
        
        # Initialize results dictionary to store all metrics
        self.results = {}
        
        # For storing multiple run results
        self.all_runs = []
    
    def _get_valid_molecules(self, molecule_set: 'MoleculeSet') -> List['Molecule']:
        """
        Get all valid molecules from a molecule set.
        
        Parameters
        ----------
        molecule_set : MoleculeSet
            Set of molecules to check
            
        Returns
        -------
        list
            List of valid Molecule objects
        """
        valid_molecules = []
        for mol in molecule_set.molecules:
            # Check if the molecule has a valid RDKit representation
            if mol.rdkit_mol is not None:
                valid_molecules.append(mol)
        return valid_molecules
    
    def _get_unique_molecules(self, molecules: List['Molecule']) -> List['Molecule']:
        """
        Get unique molecules from a list of molecules.
        
        Parameters
        ----------
        molecules : list
            List of Molecule objects
            
        Returns
        -------
        list
            List of unique Molecule objects
        """
        if self.comparison_method == 'smiles':
            # Use SMILES string for comparison
            unique_smiles = set()
            unique_molecules = []
            
            for mol in molecules:
                smiles = mol.smiles
                if smiles not in unique_smiles:
                    unique_smiles.add(smiles)
                    unique_molecules.append(mol)
            
            return unique_molecules
        
        elif self.comparison_method == 'fingerprint':
            # Use fingerprints for comparison
            unique_molecules = []
            fingerprints = []
            
            for mol in molecules:
                fp = mol.get_fingerprint(fp_type=self.fingerprint_type)
                
                # Check if this molecule is similar to any already added
                is_unique = True
                for existing_fp in fingerprints:
                    similarity = DataStructs.TanimotoSimilarity(fp, existing_fp)
                    if similarity >= self.similarity_threshold:
                        is_unique = False
                        break
                
                if is_unique:
                    unique_molecules.append(mol)
                    fingerprints.append(fp)
            
            return unique_molecules
        
        else:
            raise ValueError(f"Unsupported comparison method: {self.comparison_method}")
    
    def _is_molecule_in_set(self, molecule: 'Molecule', molecule_set: 'MoleculeSet') -> bool:
        """
        Check if a molecule is present in a molecule set.
        
        Parameters
        ----------
        molecule : Molecule
            Molecule to check
        molecule_set : MoleculeSet
            Set of molecules to search in
            
        Returns
        -------
        bool
            True if the molecule is in the set, False otherwise
        """
        if self.comparison_method == 'smiles':
            # Use SMILES string for comparison
            molecule_smiles = molecule.smiles
            for mol in molecule_set.molecules:
                if mol.smiles == molecule_smiles:
                    return True
            return False
        
        elif self.comparison_method == 'fingerprint':
            # Use fingerprints for comparison
            molecule_fp = molecule.get_fingerprint(fp_type=self.fingerprint_type)
            
            for mol in molecule_set.molecules:
                other_fp = mol.get_fingerprint(fp_type=self.fingerprint_type)
                similarity = DataStructs.TanimotoSimilarity(molecule_fp, other_fp)
                if similarity >= self.similarity_threshold:
                    return True
            
            return False
        
        else:
            raise ValueError(f"Unsupported comparison method: {self.comparison_method}")
    
    def calculate_validity(self) -> float:
        """
        Calculate the validity metric.
        
        Returns
        -------
        float
            Percentage of valid molecules in the generated set
        """
        total_molecules = len(self.generated_set.molecules)
        if total_molecules == 0:
            return 0.0
        
        valid_molecules = self._get_valid_molecules(self.generated_set)
        validity = len(valid_molecules) / total_molecules if total_molecules > 0 else 0.0
        
        # Store in results
        self.results['validity'] = validity
        self.results['valid_molecules'] = valid_molecules
        self.results['total_molecules'] = total_molecules
        
        return validity
    
    def calculate_uniqueness(self) -> float:
        """
        Calculate the uniqueness metric.
        
        Returns
        -------
        float
            Percentage of unique molecules among valid molecules
        """
        # Ensure validity has been calculated
        if 'valid_molecules' not in self.results:
            self.calculate_validity()
        
        valid_molecules = self.results['valid_molecules']
        if len(valid_molecules) == 0:
            return 0.0
        
        unique_molecules = self._get_unique_molecules(valid_molecules)
        uniqueness = len(unique_molecules) / len(valid_molecules) if len(valid_molecules) > 0 else 0.0
        
        # Store in results
        self.results['uniqueness'] = uniqueness
        self.results['unique_molecules'] = unique_molecules
        
        return uniqueness
    
    def calculate_novelty(self) -> float:
        """
        Calculate the novelty metric.
        
        Returns
        -------
        float
            Percentage of unique valid molecules not present in the reference set
        """
        # Ensure uniqueness has been calculated
        if 'unique_molecules' not in self.results:
            self.calculate_uniqueness()
        
        unique_molecules = self.results['unique_molecules']
        if len(unique_molecules) == 0:
            return 0.0
        
        novel_molecules = []
        for mol in unique_molecules:
            if not self._is_molecule_in_set(mol, self.reference_set):
                novel_molecules.append(mol)
        
        novelty = len(novel_molecules) / len(unique_molecules) if len(unique_molecules) > 0 else 0.0
        
        # Store in results
        self.results['novelty'] = novelty
        self.results['novel_molecules'] = novel_molecules
        
        return novelty
    
    def calculate_non_identicality(self) -> float:
        """
        Calculate the non-identicality metric.
        
        Returns
        -------
        float
            Percentage of valid molecules not present in the input set
        """
        # Ensure validity has been calculated
        if 'valid_molecules' not in self.results:
            self.calculate_validity()
        
        valid_molecules = self.results['valid_molecules']
        if len(valid_molecules) == 0 or self.input_set is None:
            return 0.0
        
        non_identical_molecules = []
        for mol in valid_molecules:
            if not self._is_molecule_in_set(mol, self.input_set):
                non_identical_molecules.append(mol)
        
        non_identicality = len(non_identical_molecules) / len(valid_molecules) if len(valid_molecules) > 0 else 0.0
        
        # Store in results
        self.results['non_identicality'] = non_identicality
        self.results['non_identical_molecules'] = non_identical_molecules
        
        return non_identicality
    
    def calculate_effective_novelty(self) -> float:
        """
        Calculate the effective novelty metric.
        
        Returns
        -------
        float
            Percentage of molecules that are both novel and non-identical
        """
        # Ensure novelty and non-identicality have been calculated
        if 'novel_molecules' not in self.results:
            self.calculate_novelty()
        
        if 'non_identical_molecules' not in self.results and self.input_set is not None:
            self.calculate_non_identicality()
        
        total_molecules = self.results['total_molecules']
        if total_molecules == 0:
            return 0.0
        
        # If input set is None, use only novelty for effective novelty
        if self.input_set is None:
            novel_molecules = self.results['novel_molecules']
            effective_novelty = len(novel_molecules) / total_molecules
        else:
            # Get molecules that are both novel and non-identical
            novel_smiles = {mol.smiles for mol in self.results['novel_molecules']}
            non_identical_smiles = {mol.smiles for mol in self.results['non_identical_molecules']}
            
            effective_smiles = novel_smiles.intersection(non_identical_smiles)
            effective_novelty = len(effective_smiles) / total_molecules
        
        # Store in results
        self.results['effective_novelty'] = effective_novelty
        
        return effective_novelty
    
    def evaluate_all_metrics(self) -> Dict[str, float]:
        """
        Calculate all metrics and return the results.
        
        Returns
        -------
        dict
            Dictionary containing all calculated metrics
        """
        self.calculate_validity()
        self.calculate_uniqueness()
        self.calculate_novelty()
        
        if self.input_set is not None:
            self.calculate_non_identicality()
        
        self.calculate_effective_novelty()
        
        # Extract just the numerical metrics
        metrics = {
            'validity': self.results['validity'],
            'uniqueness': self.results['uniqueness'],
            'novelty': self.results['novelty'],
            'effective_novelty': self.results['effective_novelty']
        }
        
        if self.input_set is not None:
            metrics['non_identicality'] = self.results['non_identicality']
        
        return metrics
    
    def add_run(self, run_name: str = None) -> None:
        """
        Add the current results as a new run for statistical analysis.
        
        Parameters
        ----------
        run_name : str, optional
            Name to identify this run
        """
        # Ensure all metrics are calculated
        metrics = self.evaluate_all_metrics()
        
        # Add run name if provided
        if run_name is not None:
            metrics['run_name'] = run_name
        else:
            metrics['run_name'] = f"Run {len(self.all_runs) + 1}"
        
        # Add timestamp
        metrics['timestamp'] = datetime.now().isoformat()
        
        # Store the run
        self.all_runs.append(metrics)
    
    def get_stats(self) -> Dict[str, Dict[str, float]]:
        """
        Calculate statistics across all runs.
        
        Returns
        -------
        dict
            Dictionary containing mean, median, min, max, and std for each metric
        """
        if not self.all_runs:
            return {}
        
        # Convert to DataFrame for easy statistics
        df = pd.DataFrame(self.all_runs)
        
        # Remove non-metric columns
        non_metric_cols = ['run_name', 'timestamp']
        metric_cols = [col for col in df.columns if col not in non_metric_cols]
        
        stats = {}
        for metric in metric_cols:
            stats[metric] = {
                'mean': df[metric].mean(),
                'median': df[metric].median(),
                'min': df[metric].min(),
                'max': df[metric].max(),
                'std': df[metric].std()
            }
        
        return stats
    
    def plot_metrics(self, title: str = "Molecule Set Metrics", save_path: Optional[str] = None) -> None:
        """
        Create a spider chart of the calculated metrics.
        
        Parameters
        ----------
        title : str, default="Molecule Set Metrics"
            Title for the plot
        save_path : str, optional
            Path to save the plot to a file
        """
        # Ensure all metrics are calculated
        metrics = self.evaluate_all_metrics()
        
        # Set up the spider chart
        categories = list(metrics.keys())
        values = list(metrics.values())
        
        # Number of variables
        N = len(categories)
        
        # Create angles for each metric (divide the plot into equal parts)
        angles = [n / N * 2 * np.pi for n in range(N)]
        angles += angles[:1]  # Close the loop
        
        # Add values to complete the loop
        values += values[:1]
        
        # Create the plot
        fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))
        
        # Draw the outline of the chart
        ax.plot(angles, values, linewidth=2, linestyle='solid')
        
        # Fill the area
        ax.fill(angles, values, alpha=0.25)
        
        # Set the labels for each axis
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(categories)
        
        # Add chart title
        ax.set_title(title, size=15, pad=20)
        
        # Set y-axis limits
        ax.set_ylim(0, 1)
        
        # Add gridlines
        ax.grid(True)
        
        # Save if requested
        if save_path:
            plt.savefig(save_path, bbox_inches='tight')
        
        plt.show()
    
    def plot_multiple_runs(self, title: str = "Comparison of Runs", save_path: Optional[str] = None, fontsize = 14) -> None:
        """
        Create a spider chart comparing multiple runs.
        
        Parameters
        ----------
        title : str, default="Comparison of Runs"
            Title for the plot
        save_path : str, optional
            Path to save the plot to a file
        """
        if not self.all_runs:
            print("No runs to plot. Add runs using add_run() method first.")
            return
        
        # Create DataFrame from all runs
        df = pd.DataFrame(self.all_runs)
        
        # Remove non-metric columns
        non_metric_cols = ['timestamp']
        metric_cols = [col for col in df.columns if col not in non_metric_cols and col != 'run_name']
        
        # Number of variables
        N = len(metric_cols)
        
        # Create angles for each metric
        angles = [n / N * 2 * np.pi for n in range(N)]
        angles += angles[:1]  # Close the loop
        
        # Create the plot
        fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(polar=True))
        
        # Define a color cycle
        colors = plt.cm.tab10(np.linspace(0, 1, len(self.all_runs)))
        
        # Plot each run
        for i, run in enumerate(self.all_runs):
            values = [run[metric] for metric in metric_cols]
            values += values[:1]  # Close the loop
            
            ax.plot(angles, values, linewidth=2, linestyle='solid', color=colors[i], 
                    label=run['run_name'])
            ax.fill(angles, values, alpha=0.1, color=colors[i])
        
        # Set the labels for each axis
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(metric_cols, fontsize=fontsize)
        
        # Add chart title and legend
        ax.set_title(title, size=15, pad=20)
        ax.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))
        
        # Set y-axis limits
        ax.set_ylim(0, 1)
        
        # Add gridlines
        ax.grid(True)
        
        # Save if requested
        if save_path:
            plt.savefig(save_path, bbox_inches='tight')
        
        plt.show()
    
    def save_results(self, filepath: str) -> None:
        """
        Save results to a JSON file.
        
        Parameters
        ----------
        filepath : str
            Path to save the results to
        """
        # Prepare results for saving
        save_data = {
            'metrics': self.evaluate_all_metrics(),
            'comparison_method': self.comparison_method,
            'fingerprint_type': self.fingerprint_type,
            'similarity_threshold': self.similarity_threshold,
            'timestamp': datetime.now().isoformat()
        }
        
        # Add all runs if available
        if self.all_runs:
            save_data['runs'] = self.all_runs
            save_data['statistics'] = self.get_stats()
        
        # Save to file
        with open(filepath, 'w') as f:
            json.dump(save_data, f, indent=4)
    
    def print_report(self) -> None:
        """
        Print a formatted report of all metrics.
        """
        # Ensure all metrics are calculated
        metrics = self.evaluate_all_metrics()
        
        print("=" * 60)
        print(f"MOLECULE SET EVALUATION REPORT")
        print("=" * 60)
        print(f"Comparison Method: {self.comparison_method}")
        if self.comparison_method == 'fingerprint':
            print(f"Fingerprint Type: {self.fingerprint_type}")
            print(f"Similarity Threshold: {self.similarity_threshold}")
        print("-" * 60)
        
        # Print each metric
        for metric, value in metrics.items():
            print(f"{metric.capitalize()}: {value:.4f} ({value*100:.2f}%)")
        
        print("-" * 60)
        print(f"Total molecules: {self.results['total_molecules']}")
        print(f"Valid molecules: {len(self.results['valid_molecules'])}")
        print(f"Unique molecules: {len(self.results['unique_molecules'])}")
        print(f"Novel molecules: {len(self.results['novel_molecules'])}")
        
        if self.input_set is not None:
            print(f"Non-identical molecules: {len(self.results['non_identical_molecules'])}")
        
        print("=" * 60)
        
        # Print statistics if multiple runs available
        if self.all_runs:
            stats = self.get_stats()
            print("\nSTATISTICS ACROSS RUNS")
            print("-" * 60)
            
            for metric, stat_values in stats.items():
                print(f"{metric.capitalize()}:")
                for stat_name, stat_value in stat_values.items():
                    print(f"  {stat_name}: {stat_value:.4f}")
                print()
            
            print("=" * 60)


# Example usage
if __name__ == "__main__":
    # This is just a placeholder for demonstration 
    from rdkit import Chem
    
    # Create sample molecule sets
    generated_smiles = ["CC(=O)OC1=CC=CC=C1C(=O)O", "CCO", "CCCC", "CCC(=O)O", "CC(C)CC"]
    reference_smiles = ["CC(=O)OC1=CC=CC=C1C(=O)O", "CCO", "CCCCCO"]
    input_smiles = ["CC(=O)OC1=CC=CC=C1C(=O)O"]
    
    # Create MoleculeSet objects
    generated_set = MoleculeSet()
    generated_set.add_from_smiles_list(generated_smiles)
    
    reference_set = MoleculeSet()
    reference_set.add_from_smiles_list(reference_smiles)
    
    input_set = MoleculeSet()
    input_set.add_from_smiles_list(input_smiles)
    
    # Create evaluator
    evaluator = MoleculeSetEvaluator(
        generated_set=generated_set,
        reference_set=reference_set,
        input_set=input_set,
        comparison_method='smiles'
    )
    
    # Run evaluation
    evaluator.print_report()
    
    # Add multiple runs
    evaluator.add_run("Run 1")
    
    # Change comparison method and add another run
    evaluator.comparison_method = 'fingerprint'
    evaluator.add_run("Run 2")
    
    # Plot comparison
    evaluator.plot_multiple_runs()
    
    # Save results
    evaluator.save_results("molecule_evaluation_results.json")