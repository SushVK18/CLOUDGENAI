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

from typing import List, Dict, Optional, Union, Tuple, Any
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold
from collections import Counter
import math
from scipy import stats


class MoleculeSetEvaluator:
    """
    A class for evaluating and comparing molecule sets.
    
    This class provides methods for calculating diversity metrics,
    comparing sets of molecules, and analyzing their properties.
    
    Attributes
    ----------
    target_set : MoleculeSet
        The primary set of molecules to evaluate
    reference_set : MoleculeSet, optional
        A secondary set to compare against
    seed_molecules : List[Molecule]
        List of reference/seed molecules for novelty calculations
    """
    
    def __init__(
        self,
        target_set,
        reference_set=None,
        seed_molecules=None
    ):
        """
        Initialize a MoleculeSetEvaluator.
        
        Parameters
        ----------
        target_set : MoleculeSet
            The primary set of molecules to evaluate
        reference_set : MoleculeSet, optional
            A secondary set to compare against
        seed_molecules : List[Union[Molecule, str]], optional
            List of reference/seed molecules or SMILES strings for novelty calculations
            
        Examples
        --------
        >>> target = MoleculeSet([Molecule("CC(=O)OC1=CC=CC=C1C(=O)O", name="Aspirin")])
        >>> reference = MoleculeSet([Molecule("CC1=CC=C(C=C1)C(=O)OC", name="Methyl p-toluate")])
        >>> # Using a molecule object as seed
        >>> seed_mol = Molecule("CCO", name="Ethanol")
        >>> evaluator1 = MoleculeSetEvaluator(target, reference, seed_molecules=[seed_mol])
        >>> # Using SMILES strings as seeds
        >>> evaluator2 = MoleculeSetEvaluator(target, reference, seed_molecules=["CCO", "CCCO"])
        """
        self.target_set = target_set
        self.reference_set = reference_set
        
        # Handle seed molecules, converting SMILES strings to Molecule objects if needed
        self.seed_molecules = []
        if seed_molecules:
            for seed in seed_molecules:
                if isinstance(seed, str):
                    # Convert SMILES string to Molecule object
                    from rdkit import Chem
                    # Import Molecule class dynamically to avoid circular imports
                    # This assumes Molecule class is available in the current context
                    try:
                        # Try to get the Molecule class from the same module as target_set
                        if hasattr(target_set, 'molecules') and len(target_set.molecules) > 0:
                            Molecule = target_set.molecules[0].__class__
                            self.seed_molecules.append(Molecule(seed, name=f"Seed_{len(self.seed_molecules)+1}"))
                        else:
                            # Fallback - get the first module-level reference to Molecule class
                            from inspect import getmodule
                            from sys import modules
                            module = getmodule(target_set.__class__)
                            if hasattr(module, 'Molecule'):
                                Molecule = module.Molecule
                                self.seed_molecules.append(Molecule(seed, name=f"Seed_{len(self.seed_molecules)+1}"))
                            else:
                                # Last resort - try to use a direct import
                                try:
                                    from molecule import Molecule
                                    self.seed_molecules.append(Molecule(seed, name=f"Seed_{len(self.seed_molecules)+1}"))
                                except ImportError:
                                    raise ValueError(f"Cannot convert SMILES to Molecule: {seed}. Please provide Molecule objects.")
                    except:
                        raise ValueError(f"Failed to convert SMILES to Molecule: {seed}. Please provide Molecule objects.")
                else:
                    # Assume it's already a Molecule object
                    self.seed_molecules.append(seed)
        
        # Initialize cache for expensive computations
        self._cache = {}
    
    def calculate_diversity_metrics(self, fp_type: str = 'morgan', sim_type: str = 'tanimoto') -> Dict[str, float]:
        """
        Calculate diversity metrics for the target set.
        
        Parameters
        ----------
        fp_type : str, default='morgan'
            Type of fingerprint to use for calculations
        sim_type : str, default='tanimoto'
            Type of similarity metric to use
            
        Returns
        -------
        Dict[str, float]
            Dictionary of diversity metrics
            
        Examples
        --------
        >>> evaluator = MoleculeSetEvaluator(molecule_set)
        >>> diversity_metrics = evaluator.calculate_diversity_metrics()
        >>> print(f"Internal diversity: {diversity_metrics['mean_diversity']}")
        """
        # Check if we have enough molecules
        if len(self.target_set) <= 1:
            return {
                "mean_diversity": 0.0,
                "median_diversity": 0.0,
                "min_diversity": 0.0,
                "max_diversity": 0.0,
                "shannon_entropy": 0.0,
                "simpson_diversity": 0.0,
                "unique_scaffolds_fraction": 0.0,
                "num_molecules": len(self.target_set)
            }
        
        # Check cache to avoid recalculation
        cache_key = f"div_metrics_{fp_type}_{sim_type}"
        if cache_key in self._cache:
            return self._cache[cache_key]
        
        # Get pairwise similarity matrix
        sim_matrix = self.target_set.calculate_pairwise_similarity(fp_type=fp_type)
        
        # Convert to distance matrix (1 - similarity)
        dist_matrix = 1.0 - sim_matrix
        
        # Exclude self-comparisons (diagonal)
        n = len(self.target_set)
        mask = ~np.eye(n, dtype=bool)  # mask to exclude diagonal
        pairwise_distances = dist_matrix[mask].reshape(n, n-1)
        
        # Calculate basic statistics on pairwise distances
        mean_div = np.mean(pairwise_distances)
        median_div = np.median(pairwise_distances)
        min_div = np.min(pairwise_distances)
        max_div = np.max(pairwise_distances)
        
        # Calculate Shannon entropy
        # Normalize similarity matrix
        sim_sum = np.sum(sim_matrix)
        prob_matrix = sim_matrix / sim_sum if sim_sum > 0 else np.zeros_like(sim_matrix)
        
        # Entropy calculation
        entropy = 0.0
        for prob in prob_matrix.flatten():
            if prob > 0:
                entropy -= prob * math.log2(prob)
        
        # Calculate Simpson's Diversity Index (1 - sum(p_i^2))
        # Where p_i is the proportional abundance of the ith type
        # Use similarity as a proxy for similarity of types
        simpson = 1.0 - np.sum(prob_matrix**2)
        
        # Calculate scaffold diversity
        scaffolds = []
        for mol in self.target_set:
            try:
                # Get Murcko scaffold as SMILES
                scaffold = MurckoScaffold.MurckoScaffoldSmiles(
                    mol=mol.rdkit_mol, 
                    includeChirality=False
                )
                scaffolds.append(scaffold)
            except:
                # In case of failure (e.g., no scaffold), use empty string
                scaffolds.append("")
        
        # Count unique scaffolds
        unique_scaffolds = len(set(scaffolds))
        scaffold_fraction = unique_scaffolds / len(self.target_set)
        
        # Prepare results
        result = {
            "mean_diversity": mean_div,
            "median_diversity": median_div,
            "min_diversity": min_div,
            "max_diversity": max_div,
            "shannon_entropy": entropy,
            "simpson_diversity": simpson,
            "unique_scaffolds_fraction": scaffold_fraction,
            "unique_scaffolds_count": unique_scaffolds,
            "num_molecules": len(self.target_set)
        }
        
        # Cache results
        self._cache[cache_key] = result
        
        return result
    
    def _get_scaffolds(self, molecule_set) -> List[str]:
        """
        Extract Murcko scaffolds from a molecule set.
        
        Parameters
        ----------
        molecule_set : MoleculeSet
            Set of molecules to extract scaffolds from
            
        Returns
        -------
        List[str]
            List of scaffold SMILES strings
        """
        scaffolds = []
        for mol in molecule_set:
            try:
                scaffold = MurckoScaffold.MurckoScaffoldSmiles(
                    mol=mol.rdkit_mol, 
                    includeChirality=False
                )
                scaffolds.append(scaffold)
            except:
                scaffolds.append("")
        return scaffolds
    
    def calculate_novelty_from_seeds(self, fp_type: str = 'morgan', sim_type: str = 'tanimoto') -> Dict[str, float]:
        """
        Calculate novelty metrics comparing target set to seed molecules.
        
        Parameters
        ----------
        fp_type : str, default='morgan'
            Type of fingerprint to use for calculations
        sim_type : str, default='tanimoto'
            Type of similarity metric to use
            
        Returns
        -------
        Dict[str, float]
            Dictionary of novelty metrics
            
        Examples
        --------
        >>> seed = Molecule("CC(=O)OC1=CC=CC=C1C(=O)O", name="Aspirin")
        >>> evaluator = MoleculeSetEvaluator(molecule_set, seed_molecules=[seed])
        >>> novelty_metrics = evaluator.calculate_novelty_from_seeds()
        >>> print(f"Average distance from seed: {novelty_metrics['mean_distance']}")
        """
        # Check if we have seed molecules
        if not self.seed_molecules:
            return {
                "error": "No seed molecules provided",
                "mean_distance": None,
                "median_distance": None,
                "min_distance": None,
                "max_distance": None,
                "novel_scaffold_fraction": None
            }
            
        # Check cache to avoid recalculation
        cache_key = f"novelty_metrics_{fp_type}_{sim_type}"
        if cache_key in self._cache:
            return self._cache[cache_key]
            
        # Calculate nearest distance from each target molecule to any seed
        distances = []
        nearest_seeds = []
        
        for target_mol in self.target_set:
            # Find most similar seed molecule
            max_sim = -1
            nearest_seed = None
            
            for seed_mol in self.seed_molecules:
                sim = target_mol.calculate_similarity(seed_mol, fp_type=fp_type, sim_type=sim_type)
                if sim > max_sim:
                    max_sim = sim
                    nearest_seed = seed_mol
            
            # Convert similarity to distance
            distance = 1.0 - max_sim
            distances.append(distance)
            nearest_seeds.append(nearest_seed)
        
        # Calculate basic statistics
        mean_dist = np.mean(distances) if distances else 0.0
        median_dist = np.median(distances) if distances else 0.0
        min_dist = np.min(distances) if distances else 0.0
        max_dist = np.max(distances) if distances else 0.0
        
        # Calculate scaffold novelty
        # Get scaffolds of target molecules
        target_scaffolds = self._get_scaffolds(self.target_set)
        
        # Get scaffolds of seed molecules
        seed_scaffolds = []
        for seed in self.seed_molecules:
            try:
                scaffold = MurckoScaffold.MurckoScaffoldSmiles(
                    mol=seed.rdkit_mol, 
                    includeChirality=False
                )
                seed_scaffolds.append(scaffold)
            except:
                seed_scaffolds.append("")
        
        # Count novel scaffolds (not in seeds)
        novel_scaffolds = 0
        seed_scaffold_set = set(seed_scaffolds)
        
        for scaffold in target_scaffolds:
            if scaffold and scaffold not in seed_scaffold_set:
                novel_scaffolds += 1
        
        # Calculate novel scaffold fraction
        novel_scaffold_fraction = novel_scaffolds / len(target_scaffolds) if target_scaffolds else 0.0
        
        # Prepare distance distribution metrics
        result = {
            "mean_distance": mean_dist,
            "median_distance": median_dist,
            "min_distance": min_dist,
            "max_distance": max_dist,
            "novel_scaffold_fraction": novel_scaffold_fraction,
            "novel_scaffold_count": novel_scaffolds,
            "total_scaffolds": len(target_scaffolds)
        }
        
        # Add percentage of molecules with distance > thresholds
        thresholds = [0.3, 0.4, 0.5, 0.6, 0.7]
        for threshold in thresholds:
            result[f"fraction_distance_gt_{threshold}"] = np.mean([d > threshold for d in distances]) if distances else 0.0
        
        # Cache results
        self._cache[cache_key] = result
        
        return result
    
    def visualize_fp_comparison_old(self, 
                               fp_type1: str = 'ecfp', 
                               fp_type2: str = 'maccs',
                               title: str = 'Fingerprint Dissimilarity Comparison',
                               figsize: tuple = (12, 10),
                               hist_bins: int = 40,
                               scatter_bins: int = 100,
                               cmap: str = 'viridis',
                               point_size: int = 5,
                               alpha: float = 0.7):
        """
        Visualize the comparison of dissimilarity values between two fingerprint types,
        similar to the visualization in the reference image.
        
        Parameters
        ----------
        fp_type1 : str, default='ecfp'
            First fingerprint type
        fp_type2 : str, default='maccs'
            Second fingerprint type
        title : str, default='Fingerprint Dissimilarity Comparison'
            Title of the plot
        figsize : tuple, default=(12, 10)
            Figure size (width, height)
        hist_bins : int, default=40
            Number of bins for the histograms
        scatter_bins : int, default=100
            Number of bins for the scatter plot density calculation
        cmap : str, default='viridis'
            Colormap for the density plot
        point_size : int, default=5
            Size of scatter points
        alpha : float, default=0.7
            Transparency of the scatter points
            
        Returns
        -------
        tuple
            (fig, axs) - matplotlib figure and axes objects
            
        Examples
        --------
        >>> evaluator = MoleculeSetEvaluator(molecule_set)
        >>> fig, axs = evaluator.visualize_fp_comparison('ecfp', 'fcfp')
        >>> plt.show()
        """
        # Try to import visualization libraries
        try:
            import matplotlib.pyplot as plt
            from matplotlib.colors import LinearSegmentedColormap, LogNorm
            import matplotlib.gridspec as gridspec
            from scipy.stats import gaussian_kde
        except ImportError as e:
            raise ImportError(f"Required visualization libraries not available: {e}. Install with 'pip install matplotlib scipy'")
            
        # Get comparison data
        data = self.calculate_pairwise_fp_comparison(fp_types=[fp_type1, fp_type2])
        
        # Extract data
        key = f"{fp_type1}_vs_{fp_type2}"
        
        # Ensure the key exists in the data
        if key not in data['comparison_data']:
            key = f"{fp_type2}_vs_{fp_type1}"
            if key not in data['comparison_data']:
                raise ValueError(f"No comparison data found for {fp_type1} vs {fp_type2}")
            # Swap x and y if needed
            x = data['comparison_data'][key]['dissimilarity2']
            y = data['comparison_data'][key]['dissimilarity1']
        else:
            x = data['comparison_data'][key]['dissimilarity1']
            y = data['comparison_data'][key]['dissimilarity2']
            
        # Create figure with custom grid
        fig = plt.figure(figsize=figsize)
        gs = gridspec.GridSpec(2, 2, width_ratios=[3, 1], height_ratios=[1, 3])
        
        # Main scatter plot
        ax_scatter = fig.add_subplot(gs[1, 0])
        
        # Top histogram
        ax_hist_top = fig.add_subplot(gs[0, 0], sharex=ax_scatter)
        
        # Right histogram/colorbar
        ax_hist_right = fig.add_subplot(gs[1, 1], sharey=ax_scatter)
        
        # Calculate point density for coloring
        # Use a 2D histogram approach for efficiency
        hist, xedges, yedges = np.histogram2d(x, y, bins=scatter_bins)
        
        # Find the bin for each point
        x_indices = np.digitize(x, xedges) - 1
        y_indices = np.digitize(y, yedges) - 1
        
        # Clip indices to valid range
        x_indices = np.clip(x_indices, 0, len(xedges) - 2)
        y_indices = np.clip(y_indices, 0, len(yedges) - 2)
        
        # Get density for each point
        point_densities = hist[x_indices, y_indices]
        
        # Create custom colormap similar to the image (purple to yellow-green)
        colors = [(0.2, 0.1, 0.3),    # Dark purple
                 (0.35, 0.1, 0.6),   # Purple
                 (0.0, 0.6, 0.6),    # Teal
                 (0.0, 0.8, 0.3),    # Green
                 (0.9, 0.9, 0.0)]    # Yellow
        
        custom_cmap = LinearSegmentedColormap.from_list('purple_to_yellow', colors)
        
        # Create scatter plot with density coloring
        scatter = ax_scatter.scatter(x, y, c=point_densities, 
                                   cmap=custom_cmap,
                                   s=point_size, 
                                   alpha=alpha,
                                   norm=LogNorm())
        
        # Calculate linear regression line
        slope, intercept = np.polyfit(x, y, 1)
        x_line = np.array([0, 1])
        y_line = slope * x_line + intercept
        
        # Add regression line
        ax_scatter.plot(x_line, y_line, 'r-', linewidth=1.5)
        
        # Add diagonal identity line
        ax_scatter.plot([0, 1], [0, 1], 'k--', alpha=0.5, linewidth=1)
        
        # Set axis labels
        ax_scatter.set_xlabel(f'Tanimoto dissimilarity ({fp_type1})', fontsize=12)
        ax_scatter.set_ylabel(f'Tanimoto dissimilarity ({fp_type2})', fontsize=12)
        
        # Set limits
        ax_scatter.set_xlim(0, 1)
        ax_scatter.set_ylim(0, 1)
        
        # Create top histogram
        counts, edges, _ = ax_hist_top.hist(x, bins=hist_bins, color='navy', 
                                         alpha=0.7, density=True)
        
        # Label the top histogram count
        ax_hist_top.text(0.98, 0.98, f"{len(x)}", 
                      transform=ax_hist_top.transAxes,
                      horizontalalignment='right',
                      verticalalignment='top',
                      fontsize=12)
        
        # Add vertical line at median
        median_x = np.median(x)
        ax_hist_top.axvline(median_x, color='red', linestyle='--', linewidth=1.5)
        
        # Replace y-axis on right with a colorbar
        cbar = fig.colorbar(scatter, cax=ax_hist_right)
        cbar.set_label('Count', fontsize=12)
        
        # Add correlation coefficient
        corr = data['comparison_data'][key]['correlation']
        ax_scatter.text(0.05, 0.95, f'r = {corr:.3f}', 
                      transform=ax_scatter.transAxes,
                      horizontalalignment='left',
                      verticalalignment='top',
                      fontsize=12,
                      bbox=dict(facecolor='white', alpha=0.8))
        
        # Hide tick labels on histograms
        ax_hist_top.set_yticks([])
        
        # Add grid to main plot
        ax_scatter.grid(True, linestyle='--', alpha=0.3)
        
        # Set title
        fig.suptitle(title, fontsize=14)
        
        # Add subplot label
        ax_scatter.text(-0.1, -0.1, '(a)', transform=ax_scatter.transAxes,
                      fontsize=14, fontweight='bold')
        
        # Tight layout
        plt.tight_layout()
        
        return fig, (ax_scatter, ax_hist_top, ax_hist_right)
    
    def visualize_fp_comparison(self, 
                               fp_type1: str = 'ecfp', 
                               fp_type2: str = 'maccs',
                               title: str = 'Fingerprint Dissimilarity Comparison',
                               figsize: tuple = (12, 10),
                               hist_bins: int = 40,
                               scatter_bins: int = 100,
                               cmap: str = 'viridis',
                               point_size: int = 5,
                               alpha: float = 0.7):
        """
        Visualize the comparison of dissimilarity values between two fingerprint types,
        similar to the visualization in the reference image.
        
        Parameters
        ----------
        fp_type1 : str, default='ecfp'
            First fingerprint type
        fp_type2 : str, default='maccs'
            Second fingerprint type
        title : str, default='Fingerprint Dissimilarity Comparison'
            Title of the plot
        figsize : tuple, default=(12, 10)
            Figure size (width, height)
        hist_bins : int, default=40
            Number of bins for the histograms
        scatter_bins : int, default=100
            Number of bins for the scatter plot density calculation
        cmap : str, default='viridis'
            Colormap for the density plot
        point_size : int, default=5
            Size of scatter points
        alpha : float, default=0.7
            Transparency of the scatter points
            
        Returns
        -------
        tuple
            (fig, axs) - matplotlib figure and axes objects
            
        Examples
        --------
        >>> evaluator = MoleculeSetEvaluator(molecule_set)
        >>> fig, axs = evaluator.visualize_fp_comparison('ecfp', 'fcfp')
        >>> plt.show()
        """
        # Try to import visualization libraries
        try:
            import matplotlib.pyplot as plt
            from matplotlib.colors import LinearSegmentedColormap, LogNorm
            import matplotlib.gridspec as gridspec
            from scipy.stats import gaussian_kde
            import numpy as np
        except ImportError as e:
            raise ImportError(f"Required visualization libraries not available: {e}. Install with 'pip install matplotlib scipy numpy'")
            
        # Get comparison data
        data = self.calculate_pairwise_fp_comparison(fp_types=[fp_type1, fp_type2])
        
        # Extract data
        key = f"{fp_type1}_vs_{fp_type2}"
        
        # Ensure the key exists in the data
        if key not in data['comparison_data']:
            key = f"{fp_type2}_vs_{fp_type1}"
            if key not in data['comparison_data']:
                raise ValueError(f"No comparison data found for {fp_type1} vs {fp_type2}")
            # Swap x and y if needed
            x = data['comparison_data'][key]['dissimilarity2']
            y = data['comparison_data'][key]['dissimilarity1']
            # Also swap labels
            fp_type1, fp_type2 = fp_type2, fp_type1
        else:
            x = data['comparison_data'][key]['dissimilarity1']
            y = data['comparison_data'][key]['dissimilarity2']
        
        # Create figure with custom grid
        fig = plt.figure(figsize=figsize)
        
        # Create a gridspec with 2x2 grid
        gs = gridspec.GridSpec(2, 2, width_ratios=[3, 1], height_ratios=[1, 3],
                              wspace=0.05, hspace=0.05)
        
        # Main scatter plot
        ax_scatter = fig.add_subplot(gs[1, 0])
        
        # Top histogram
        ax_hist_top = fig.add_subplot(gs[0, 0], sharex=ax_scatter)
        
        # Right histogram 
        ax_hist_right = fig.add_subplot(gs[1, 1], sharey=ax_scatter)
        
        # Calculate point density
        try:
            # Try to use KDE for smooth density estimation
            xy = np.vstack([x, y])
            z = gaussian_kde(xy)(xy)
            
            # Sorting the points by density allows higher density points to appear on top
            idx = z.argsort()
            x_sorted, y_sorted, z_sorted = x[idx], y[idx], z[idx]
            
            # Create custom colormap (purple to yellow-green)
            colors = [(0.2, 0.1, 0.3),    # Dark purple
                     (0.35, 0.1, 0.6),   # Purple
                     (0.0, 0.6, 0.6),    # Teal
                     (0.0, 0.8, 0.3),    # Green
                     (0.9, 0.9, 0.0)]    # Yellow
            
            custom_cmap = LinearSegmentedColormap.from_list('purple_to_yellow', colors)
            
            # Create scatter plot with density coloring
            scatter = ax_scatter.scatter(x_sorted, y_sorted, c=z_sorted,
                                       cmap=custom_cmap, 
                                       s=point_size,
                                       alpha=alpha)
        except Exception as e:
            # Fallback to simple scatter without density coloring
            print(f"Warning: Could not calculate density: {e}. Using simple scatter.")
            scatter = ax_scatter.scatter(x, y, s=point_size, alpha=alpha)
        
        # Create top histogram for fp_type1 (x-axis)
        ax_hist_top.hist(x, bins=hist_bins, color='#4472C4', alpha=0.7, density=True)
        
        # Create right histogram for fp_type2 (y-axis)
        ax_hist_right.hist(y, bins=hist_bins, color='#4472C4', alpha=0.7, 
                         density=True, orientation='horizontal')
        
        # Add identity line (diagonal)
        ax_scatter.plot([0, 1], [0, 1], 'r--', linewidth=1.5, label="Identity")
        
        # Calculate linear regression and add regression line
        slope, intercept = np.polyfit(x, y, 1)
        x_line = np.array([0, 1])
        y_line = slope * x_line + intercept
        ax_scatter.plot(x_line, y_line, 'k-', linewidth=1.0)
        
        # Set axis labels
        ax_scatter.set_xlabel(f'Tanimoto dissimilarity ({fp_type1})', fontsize=12)
        ax_scatter.set_ylabel(f'Tanimoto dissimilarity ({fp_type2})', fontsize=12)
        
        # Set limits
        ax_scatter.set_xlim(0, 1)
        ax_scatter.set_ylim(0, 1)
        
        # Add correlation coefficient
        corr = data['comparison_data'][key]['correlation']
        ax_scatter.text(0.05, 0.95, f'r = {corr:.3f}',
                      transform=ax_scatter.transAxes,
                      horizontalalignment='left',
                      verticalalignment='top',
                      fontsize=12,
                      bbox=dict(facecolor='white', alpha=0.8))
        
        # Show number of points on top histogram
        point_count = len(x)
        ax_hist_top.text(0.95, 0.95, f"{point_count:,}",
                       transform=ax_hist_top.transAxes,
                       horizontalalignment='right',
                       verticalalignment='top',
                       fontsize=12,
                       bbox=dict(facecolor='white', alpha=0.8))
        
        # Add colorbar for density
        if 'z_sorted' in locals():
            cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
            cbar = fig.colorbar(scatter, cax=cbar_ax)
            cbar.set_label('Density', fontsize=12)
        
        # Hide tick labels on histograms
        ax_hist_top.tick_params(axis='x', labelbottom=False)
        ax_hist_right.tick_params(axis='y', labelleft=False)
        
        # Set background grid
        ax_scatter.grid(True, linestyle='--', alpha=0.3)
        
        # Set title
        fig.suptitle(title, fontsize=14)
        
        # Adjust spacing
        plt.tight_layout()
        fig.subplots_adjust(top=0.92)  # Make room for the title
        
        return fig, (ax_scatter, ax_hist_top, ax_hist_right)
        
    def calculate_pairwise_fp_comparison(self, fp_types: List[str] = ['morgan', 'maccs'], sim_type: str = 'tanimoto') -> Dict:
        """
        Compare dissimilarity values between pairs of molecules using different fingerprint types.
        
        Parameters
        ----------
        fp_types : List[str], default=['morgan', 'maccs']
            List of fingerprint types to compare
        sim_type : str, default='tanimoto'
            Similarity metric to use (tanimoto, dice)
            
        Returns
        -------
        Dict
            Dictionary containing dissimilarity matrices and comparison data
            
        Examples
        --------
        >>> evaluator = MoleculeSetEvaluator(molecule_set)
        >>> comparison = evaluator.calculate_pairwise_fp_comparison(['ecfp', 'maccs'])
        """
        if len(fp_types) < 2:
            raise ValueError("At least two fingerprint types must be provided for comparison")
            
        # Cache key for this calculation
        cache_key = f"fp_comparison_{'-'.join(fp_types)}_{sim_type}"
        if cache_key in self._cache:
            return self._cache[cache_key]
            
        # Dictionary to store similarity matrices for each fingerprint type
        similarity_matrices = {}
        dissimilarity_matrices = {}
        
        # Calculate similarity matrix for each fingerprint type
        for fp_type in fp_types:
            sim_matrix = self.target_set.calculate_pairwise_similarity(fp_type=fp_type)
            similarity_matrices[fp_type] = sim_matrix
            # Convert to dissimilarity
            dissimilarity_matrices[fp_type] = 1.0 - sim_matrix
        
        # Prepare comparison data for visualization
        # For each pair of fingerprint types, create a flattened array of dissimilarity values
        comparison_data = {}
        for i, fp1 in enumerate(fp_types):
            for j, fp2 in enumerate(fp_types):
                if i < j:  # Only compare unique pairs
                    # Get the upper triangle indices (excluding diagonal)
                    n = len(self.target_set)
                    triu_indices = np.triu_indices(n, k=1)
                    
                    # Extract dissimilarity values
                    diss1 = dissimilarity_matrices[fp1][triu_indices]
                    diss2 = dissimilarity_matrices[fp2][triu_indices]
                    
                    # Create comparison data
                    comparison_data[f"{fp1}_vs_{fp2}"] = {
                        "dissimilarity1": diss1,
                        "dissimilarity2": diss2,
                        "correlation": np.corrcoef(diss1, diss2)[0, 1],
                        "mean_diff": np.mean(np.abs(diss1 - diss2)),
                        "count": len(diss1)
                    }
        
        # Create a compact version of the data for plotting
        plot_data = {
            "fingerprint_types": fp_types,
            "dissimilarity_matrices": dissimilarity_matrices,
            "comparison_data": comparison_data
        }
        
        # Cache the results
        self._cache[cache_key] = plot_data
        
        return plot_data
    
    def calculate_scaffold_diversity(self) -> Dict[str, Union[float, Dict]]:
        """
        Analyze scaffold diversity in the target set.
        
        Returns
        -------
        Dict
            Dictionary with scaffold diversity metrics and distributions
            
        Examples
        --------
        >>> evaluator = MoleculeSetEvaluator(molecule_set)
        >>> scaffold_metrics = evaluator.calculate_scaffold_diversity()
        >>> print(f"Unique scaffolds: {scaffold_metrics['unique_count']}")
        """
        scaffolds = self._get_scaffolds(self.target_set)
        
        # Count scaffolds
        scaffold_counter = Counter(scaffolds)
        total_count = len(scaffolds)
        unique_count = len(scaffold_counter)
        
        # Calculate statistics
        if total_count > 0:
            scaffold_diversity = unique_count / total_count
            
            # Shannon entropy of scaffold distribution
            entropy = 0.0
            for count in scaffold_counter.values():
                p = count / total_count
                entropy -= p * math.log2(p)
                
            # Gini coefficient (measure of inequality in distribution)
            # Higher values indicate more inequality (less diversity)
            sorted_counts = sorted(scaffold_counter.values())
            cumsum = np.cumsum(sorted_counts)
            gini = (np.sum((2 * np.arange(1, len(sorted_counts) + 1) - len(sorted_counts) - 1) 
                   * sorted_counts)) / (len(sorted_counts) * np.sum(sorted_counts))
        else:
            scaffold_diversity = 0.0
            entropy = 0.0
            gini = 0.0
        
        # Get most common scaffolds
        most_common = scaffold_counter.most_common(5)
        
        # Calculate percentage for each scaffold
        scaffold_percentages = {
            scaffold: (count / total_count * 100) for scaffold, count in scaffold_counter.most_common()
        }
        
        return {
            "unique_count": unique_count,
            "total_count": total_count,
            "diversity_ratio": scaffold_diversity,
            "shannon_entropy": entropy,
            "gini_coefficient": gini,
            "most_common": most_common,
            "scaffold_percentages": scaffold_percentages
        }
    
    def calculate_scaffold_novelty(self) -> Dict[str, Union[float, int, List]]:
        """
        Calculate scaffold-based novelty metrics comparing target set to seed molecules.
        
        Returns
        -------
        Dict
            Dictionary of scaffold novelty metrics
            
        Examples
        --------
        >>> seed = Molecule("CC(=O)OC1=CC=CC=C1C(=O)O", name="Aspirin")
        >>> evaluator = MoleculeSetEvaluator(molecule_set, seed_molecules=[seed])
        >>> novelty_metrics = evaluator.calculate_scaffold_novelty()
        """
        # Check if we have seed molecules
        if not self.seed_molecules:
            return {
                "error": "No seed molecules provided"
            }
        
        # Get scaffolds
        target_scaffolds = self._get_scaffolds(self.target_set)
        seed_scaffolds = self._get_scaffolds(self.seed_molecules)
        
        # Count scaffolds
        target_scaffold_counter = Counter(target_scaffolds)
        seed_scaffold_counter = Counter(seed_scaffolds)
        
        # Find novel scaffolds
        novel_scaffolds = []
        for scaffold in target_scaffold_counter:
            if scaffold and scaffold not in seed_scaffold_counter:
                novel_scaffolds.append(scaffold)
        
        # Calculate novelty metrics
        total_target_scaffolds = len(target_scaffolds)
        unique_target_scaffolds = len(target_scaffold_counter)
        unique_seed_scaffolds = len(seed_scaffold_counter)
        novel_scaffold_count = len(novel_scaffolds)
        
        # Calculate fractions
        if total_target_scaffolds > 0:
            novel_scaffold_fraction = novel_scaffold_count / total_target_scaffolds
            unique_novel_scaffold_fraction = novel_scaffold_count / unique_target_scaffolds if unique_target_scaffolds > 0 else 0.0
        else:
            novel_scaffold_fraction = 0.0
            unique_novel_scaffold_fraction = 0.0
        
        # Count molecules with novel scaffolds
        molecules_with_novel_scaffolds = sum(1 for s in target_scaffolds if s in novel_scaffolds)
        
        # Calculate Tanimoto scaffold novelty for more nuanced comparison
        # Convert scaffolds to fingerprints for comparison
        scaffold_fps = {}
        
        for scaffold in set(target_scaffolds + seed_scaffolds):
            if scaffold:  # Skip empty scaffolds
                mol = Chem.MolFromSmiles(scaffold)
                if mol:
                    # uses morgan fingerprint with radius 2
                    fpgen = AllChem.GetMorganGenerator(radius=2)
                    scaffold_fps[scaffold] = fpgen.GetFingerprint(mol)
                    # scaffold_fps[scaffold] = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
        
        # Calculate average scaffold similarity
        scaffold_similarities = []
        
        for t_scaffold in set(target_scaffolds):
            if t_scaffold and t_scaffold in scaffold_fps:
                max_sim = 0.0
                for s_scaffold in set(seed_scaffolds):
                    if s_scaffold and s_scaffold in scaffold_fps:
                        sim = DataStructs.TanimotoSimilarity(
                            scaffold_fps[t_scaffold], 
                            scaffold_fps[s_scaffold]
                        )
                        max_sim = max(max_sim, sim)
                scaffold_similarities.append(max_sim)
        
        avg_scaffold_similarity = np.mean(scaffold_similarities) if scaffold_similarities else 0.0
        scaffold_novelty_score = 1.0 - avg_scaffold_similarity
        
        return {
            "novel_scaffold_count": novel_scaffold_count,
            "novel_scaffold_fraction": novel_scaffold_fraction,
            "unique_novel_scaffold_fraction": unique_novel_scaffold_fraction,
            "molecules_with_novel_scaffolds": molecules_with_novel_scaffolds,
            "molecules_with_novel_scaffolds_fraction": molecules_with_novel_scaffolds / len(target_scaffolds) if target_scaffolds else 0.0,
            "scaffold_novelty_score": scaffold_novelty_score,
            "avg_scaffold_similarity": avg_scaffold_similarity,
            "novel_scaffolds": novel_scaffolds[:10] if len(novel_scaffolds) > 10 else novel_scaffolds  # Return top 10 novel scaffolds
        }
    
    def calculate_property_distributions(self, properties: List[str] = None) -> pd.DataFrame:
        """
        Calculate property distributions for the target set of molecules.
        
        Parameters
        ----------
        properties : List[str], optional
            List of specific properties to calculate. If None, uses a default set.
            
        Returns
        -------
        pandas.DataFrame
            DataFrame containing property values for all molecules
            
        Examples
        --------
        >>> evaluator = MoleculeSetEvaluator(molecule_set)
        >>> property_df = evaluator.calculate_property_distributions()
        >>> print(property_df.describe())
        """
        # Check cache
        cache_key = f"property_distributions_{'-'.join(properties) if properties else 'default'}"
        if cache_key in self._cache:
            return self._cache[cache_key]
        
        # Default set of properties if none specified
        if properties is None:
            properties = [
                "MolecularWeight", "LogP", "NumHDonors", "NumHAcceptors", 
                "NumRotatableBonds", "TPSA", "NumRings", "NumAromaticRings"
            ]
        
        # Collect property values for all molecules
        property_data = []
        for i, mol in enumerate(self.target_set):
            # Get molecule name
            mol_name = getattr(mol, 'name', f"Molecule_{i+1}")
            
            # Calculate descriptors using the Molecule class's method
            descriptors = mol.calculate_descriptors()
            
            # Create a dictionary with molecule name and properties
            mol_props = {'MoleculeID': mol_name, 'SMILES': mol.smiles}
            
            # Add requested properties
            for prop in properties:
                if prop in descriptors:
                    mol_props[prop] = descriptors[prop]
                else:
                    mol_props[prop] = None  # Property not available
            
            property_data.append(mol_props)
        
        # Convert to DataFrame
        df = pd.DataFrame(property_data)
        
        # Cache the result
        self._cache[cache_key] = df
        
        return df
    
    def calculate_property_statistics(self, properties: List[str] = None) -> Dict[str, Dict[str, float]]:
        """
        Calculate statistical measures for molecular properties in the target set.
        
        Parameters
        ----------
        properties : List[str], optional
            List of specific properties to calculate statistics for.
            
        Returns
        -------
        Dict[str, Dict[str, float]]
            Dictionary of property statistics
            
        Examples
        --------
        >>> evaluator = MoleculeSetEvaluator(molecule_set)
        >>> stats = evaluator.calculate_property_statistics()
        >>> print(f"Average molecular weight: {stats['MolecularWeight']['mean']:.2f}")
        """
        # Get property distributions
        df = self.calculate_property_distributions(properties)
        
        # Exclude non-numeric columns
        numeric_cols = df.select_dtypes(include=np.number).columns.tolist()
        
        # Calculate statistics for each property
        statistics = {}
        for prop in numeric_cols:
            if prop != 'MoleculeID':  # Skip identifier columns
                prop_stats = {
                    'mean': df[prop].mean(),
                    'median': df[prop].median(),
                    'std': df[prop].std(),
                    'min': df[prop].min(),
                    'max': df[prop].max(),
                    'q1': df[prop].quantile(0.25),
                    'q3': df[prop].quantile(0.75)
                }
                statistics[prop] = prop_stats
        
        return statistics
    
    def compare_property_distributions(self, properties: List[str] = None) -> Dict[str, Dict[str, Any]]:
        """
        Compare property distributions between target and reference sets.
        
        Parameters
        ----------
        properties : List[str], optional
            List of specific properties to compare.
            
        Returns
        -------
        Dict[str, Dict[str, Any]]
            Dictionary containing comparison results for each property
            
        Examples
        --------
        >>> target_set = MoleculeSet([...])
        >>> reference_set = MoleculeSet([...])
        >>> evaluator = MoleculeSetEvaluator(target_set, reference_set)
        >>> comparison = evaluator.compare_property_distributions()
        >>> for prop, results in comparison.items():
        >>>     print(f"{prop}: p-value = {results['ks_pvalue']:.4f}")
        """
        # Check if we have a reference set
        if self.reference_set is None:
            return {"error": "No reference set provided for comparison"}
        
        # Check cache
        cache_key = f"property_comparison_{'-'.join(properties) if properties else 'default'}"
        if cache_key in self._cache:
            return self._cache[cache_key]
        
        # Get property distributions for target set
        target_df = self.calculate_property_distributions(properties)
        
        # Create a temporary evaluator for the reference set to get its properties
        temp_evaluator = MoleculeSetEvaluator(self.reference_set)
        reference_df = temp_evaluator.calculate_property_distributions(properties)
        
        # Identify common numeric properties
        numeric_cols = set(target_df.select_dtypes(include=np.number).columns)
        numeric_cols = numeric_cols.intersection(set(reference_df.select_dtypes(include=np.number).columns))
        
        # Exclude identifier columns
        numeric_cols = [col for col in numeric_cols if col not in ['MoleculeID']]
        
        # Compare distributions for each property
        comparison_results = {}
        for prop in numeric_cols:
            # Get property values for both sets
            target_values = target_df[prop].dropna().values
            reference_values = reference_df[prop].dropna().values
            
            if len(target_values) > 0 and len(reference_values) > 0:
                # Perform Kolmogorov-Smirnov test
                ks_stat, ks_pvalue = stats.ks_2samp(target_values, reference_values)
                
                # Calculate effect size (Cohen's d)
                target_mean, target_std = np.mean(target_values), np.std(target_values)
                reference_mean, reference_std = np.mean(reference_values), np.std(reference_values)
                
                # Pooled standard deviation
                pooled_std = np.sqrt(((len(target_values) - 1) * target_std**2 + 
                                    (len(reference_values) - 1) * reference_std**2) / 
                                    (len(target_values) + len(reference_values) - 2))
                
                # Cohen's d
                cohens_d = (target_mean - reference_mean) / pooled_std if pooled_std > 0 else 0
                
                # Calculate overlap coefficient (overlapping area under the PDFs)
                # Approximate using normal distributions
                if target_std > 0 and reference_std > 0:
                    mean_diff = abs(target_mean - reference_mean)
                    overlap = 2 * stats.norm.cdf(-mean_diff / (2 * np.sqrt(target_std**2 + reference_std**2)))
                else:
                    overlap = 0.0 if target_mean != reference_mean else 1.0
                
                # Store results
                comparison_results[prop] = {
                    'target_mean': target_mean,
                    'reference_mean': reference_mean,
                    'target_std': target_std,
                    'reference_std': reference_std,
                    'ks_statistic': ks_stat,
                    'ks_pvalue': ks_pvalue,
                    'cohens_d': cohens_d,
                    'overlap_coefficient': overlap,
                    'target_count': len(target_values),
                    'reference_count': len(reference_values)
                }
            else:
                comparison_results[prop] = {'error': 'Insufficient data for comparison'}
        
        # Cache results
        self._cache[cache_key] = comparison_results
        
        return comparison_results
    
    def analyze_druglikeness(self) -> Dict[str, Dict[str, Union[float, int]]]:
        """
        Analyze druglikeness of molecules in the target set based on various rules.
        
        Returns
        -------
        Dict[str, Dict[str, Union[float, int]]]
            Dictionary containing druglikeness metrics
            
        Examples
        --------
        >>> evaluator = MoleculeSetEvaluator(molecule_set)
        >>> druglikeness = evaluator.analyze_druglikeness()
        >>> print(f"Lipinski compliant: {druglikeness['lipinski']['compliant_fraction']*100:.1f}%")
        """
        # Check cache
        cache_key = "druglikeness_analysis"
        if cache_key in self._cache:
            return self._cache[cache_key]
        
        # Get property distributions
        df = self.calculate_property_distributions()
        
        # Ensure required properties are available
        required_props = ["MolecularWeight", "LogP", "NumHDonors", "NumHAcceptors", "NumRotatableBonds", "TPSA"]
        missing_props = [prop for prop in required_props if prop not in df.columns]
        
        if missing_props:
            return {"error": f"Missing required properties: {', '.join(missing_props)}"}
        
        # Lipinski's Rule of Five analysis
        lipinski_rules = {
            'MW': df['MolecularWeight'] <= 500,
            'LogP': df['LogP'] <= 5,
            'HDonors': df['NumHDonors'] <= 5,
            'HAcceptors': df['NumHAcceptors'] <= 10
        }
        
        # Count violations for each molecule
        lipinski_violations = sum(~rule for rule in lipinski_rules.values())
        
        # Calculate compliance percentages
        lipinski_results = {
            'compliant_count': sum(lipinski_violations <= 1),
            'compliant_fraction': sum(lipinski_violations <= 1) / len(df),
            'violation_counts': {i: sum(lipinski_violations == i) for i in range(5)},
            'violation_fractions': {i: sum(lipinski_violations == i) / len(df) for i in range(5)}
        }
        
        # Veber rules analysis (Rotatable bonds <= 10, TPSA <= 140)
        veber_rules = {
            'RotBonds': df['NumRotatableBonds'] <= 10,
            'TPSA': df['TPSA'] <= 140
        }
        
        veber_violations = sum(~rule for rule in veber_rules.values())
        
        veber_results = {
            'compliant_count': sum(veber_violations == 0),
            'compliant_fraction': sum(veber_violations == 0) / len(df),
            'violation_counts': {i: sum(veber_violations == i) for i in range(3)},
            'violation_fractions': {i: sum(veber_violations == i) / len(df) for i in range(3)}
        }
        
        # Combined analysis
        combined_violations = lipinski_violations + veber_violations
        
        combined_results = {
            'compliant_count': sum(combined_violations <= 1),
            'compliant_fraction': sum(combined_violations <= 1) / len(df),
            'violation_histogram': {i: sum(combined_violations == i) for i in range(1, 8)},
            'violation_fractions': {i: sum(combined_violations == i) / len(df) for i in range(1, 8)}
        }
        
        results = {
            'lipinski': lipinski_results,
            'veber': veber_results,
            'combined': combined_results,
            'total_molecules': len(df)
        }
        
        # Cache results
        self._cache[cache_key] = results
        
        return results
    
    def visualize_property_distributions(self, 
                                        properties: List[str] = None, 
                                        figsize: tuple = (12, 10),
                                        bins: int = 30,
                                        kde: bool = True,
                                        include_reference: bool = True) -> tuple:
        """
        Visualize distributions of molecular properties.
        
        Parameters
        ----------
        properties : List[str], optional
            List of specific properties to visualize.
        figsize : tuple, default=(12, 10)
            Figure size (width, height)
        bins : int, default=30
            Number of bins for histograms
        kde : bool, default=True
            Whether to include kernel density estimation curves
        include_reference : bool, default=True
            Whether to include the reference set in the visualization
            
        Returns
        -------
        tuple
            (fig, axs) - matplotlib figure and axes objects
            
        Examples
        --------
        >>> evaluator = MoleculeSetEvaluator(molecule_set, reference_set)
        >>> fig, axs = evaluator.visualize_property_distributions(['MolecularWeight', 'LogP'])
        >>> plt.show()
        """
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
        except ImportError:
            raise ImportError("Matplotlib and seaborn are required for visualization. Install with 'pip install matplotlib seaborn'")
        
        # Default properties to visualize
        if properties is None:
            properties = ["MolecularWeight", "LogP", "NumHDonors", "NumHAcceptors", 
                         "NumRotatableBonds", "TPSA"]
        
        # Get property data for target set
        target_df = self.calculate_property_distributions(properties)
        
        # Check if we should include the reference set
        has_reference = include_reference and self.reference_set is not None
        
        if has_reference:
            # Create a temporary evaluator for the reference set to get its properties
            temp_evaluator = MoleculeSetEvaluator(self.reference_set)
            reference_df = temp_evaluator.calculate_property_distributions(properties)
        
        # Determine grid size for subplots
        n_props = len(properties)
        n_cols = min(3, n_props)
        n_rows = (n_props + n_cols - 1) // n_cols
        
        # Create figure and axes
        fig, axs = plt.subplots(n_rows, n_cols, figsize=figsize)
        
        # Flatten axes array for easier indexing
        if n_props > 1:
            axs = axs.flatten()
        else:
            axs = [axs]
        
        # Plot each property
        for i, prop in enumerate(properties):
            if prop in target_df.columns:
                # Target set histogram
                sns.histplot(target_df[prop].dropna(), ax=axs[i], bins=bins, 
                           kde=kde, color='blue', alpha=0.6, label='Target')
                
                # Reference set histogram (if available)
                if has_reference and prop in reference_df.columns:
                    sns.histplot(reference_df[prop].dropna(), ax=axs[i], bins=bins,
                               kde=kde, color='red', alpha=0.4, label='Reference')
                
                # Add mean vertical lines
                axs[i].axvline(target_df[prop].mean(), color='blue', linestyle='--', 
                             linewidth=1.5, label='Target Mean')
                
                if has_reference and prop in reference_df.columns:
                    axs[i].axvline(reference_df[prop].mean(), color='red', linestyle='--',
                                 linewidth=1.5, label='Reference Mean')
                
                # Add labels and legend
                axs[i].set_xlabel(prop)
                axs[i].set_ylabel('Frequency')
                axs[i].set_title(f'{prop} Distribution')
                axs[i].legend()
        
        # Hide unused subplots
        for i in range(n_props, len(axs)):
            axs[i].set_visible(False)
        
        # Add overall title
        title = 'Property Distributions'
        if has_reference:
            title += ' (Target vs Reference)'
        fig.suptitle(title, fontsize=16)
        
        # Adjust layout
        plt.tight_layout()
        fig.subplots_adjust(top=0.92)
        
        return fig, axs
    
    def visualize_property_radar(self, 
                               properties: List[str] = None,
                               figsize: tuple = (10, 8)) -> tuple:
        """
        Create a radar chart (spider plot) comparing property distributions
        between target and reference sets.
        
        Parameters
        ----------
        properties : List[str], optional
            List of specific properties to include in the radar chart.
        figsize : tuple, default=(10, 8)
            Figure size (width, height)
            
        Returns
        -------
        tuple
            (fig, ax) - matplotlib figure and axis objects
            
        Examples
        --------
        >>> evaluator = MoleculeSetEvaluator(target_set, reference_set)
        >>> fig, ax = evaluator.visualize_property_radar()
        >>> plt.show()
        """
        try:
            import matplotlib.pyplot as plt
            from matplotlib.path import Path
            from matplotlib.spines import Spine
            from matplotlib.transforms import Affine2D
        except ImportError:
            raise ImportError("Matplotlib is required for visualization. Install with 'pip install matplotlib'")
            
        # Check if we have a reference set
        if self.reference_set is None:
            raise ValueError("Reference set is required for radar chart visualization")
            
        # Default properties
        if properties is None:
            properties = ["MolecularWeight", "LogP", "NumHDonors", "NumHAcceptors", 
                         "NumRotatableBonds", "TPSA"]
            
        # Get property comparison data
        comparison = self.compare_property_distributions(properties)
        
        # Filter properties that are available in the comparison
        valid_props = [prop for prop in properties if prop in comparison 
                     and 'error' not in comparison[prop]]
        
        if not valid_props:
            raise ValueError("No valid properties found for radar chart visualization")
            
        # Number of variables
        N = len(valid_props)
        
        # Create angles for radar chart
        angles = np.linspace(0, 2*np.pi, N, endpoint=False).tolist()
        angles += angles[:1]  # Close the loop
        
        # Extend properties list to close the loop
        props = valid_props + [valid_props[0]]
        
        # Calculate normalized means for target and reference sets
        # Normalize means to [0, 1] range for radar chart
        target_means = []
        reference_means = []
        
        for prop in valid_props:
            target_mean = comparison[prop]['target_mean']
            reference_mean = comparison[prop]['reference_mean']
            
            # Determine min and max for normalization
            min_val = min(target_mean, reference_mean) * 0.8
            max_val = max(target_mean, reference_mean) * 1.2
            
            # Normalize
            if max_val > min_val:
                target_means.append((target_mean - min_val) / (max_val - min_val))
                reference_means.append((reference_mean - min_val) / (max_val - min_val))
            else:
                target_means.append(0.5)
                reference_means.append(0.5)
                
        # Close the loop
        target_means += target_means[:1]
        reference_means += reference_means[:1]
        
        # Create figure and polar axis
        fig, ax = plt.subplots(figsize=figsize, subplot_kw=dict(polar=True))
        
        # Plot data
        ax.plot(angles, target_means, 'o-', linewidth=2, label='Target')
        ax.fill(angles, target_means, alpha=0.25)
        ax.plot(angles, reference_means, 'o-', linewidth=2, label='Reference')
        ax.fill(angles, reference_means, alpha=0.25)
        
        # Add property labels
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(valid_props)
        
        # Add legend and title
        ax.legend(loc='upper right')
        plt.title('Property Comparison: Target vs Reference', size=15)
        
        return fig, ax
    
    
    def calculate_chemical_space_pca(self, 
                               fp_type: str = 'morgan', 
                               n_components: int = 2,
                               include_reference: bool = True) -> Dict[str, Any]:
        """
        Calculate PCA projection of chemical space based on molecular fingerprints.
        
        Parameters
        ----------
        fp_type : str, default='morgan'
            Type of fingerprint to use for the chemical space analysis
        n_components : int, default=2
            Number of principal components to calculate
        include_reference : bool, default=True
            Whether to include reference set in the analysis
            
        Returns
        -------
        Dict[str, Any]
            Dictionary containing PCA results and projections
            
        Examples
        --------
        >>> evaluator = MoleculeSetEvaluator(target_set, reference_set)
        >>> pca_results = evaluator.calculate_chemical_space_pca()
        >>> print(f"Explained variance: {pca_results['explained_variance_ratio']}")
        """
        try:
            from sklearn.decomposition import PCA
            from sklearn.preprocessing import StandardScaler
        except ImportError:
            raise ImportError("Scikit-learn is required for PCA analysis. Install with 'pip install scikit-learn'")
            
        # Check cache
        cache_key = f"chemical_space_pca_{fp_type}_{n_components}_{include_reference}"
        if cache_key in self._cache:
            return self._cache[cache_key]
            
        # Get fingerprints for target set
        target_fps = np.array([mol.get_fingerprint(fp_type=fp_type, as_numpy=True) 
                             for mol in self.target_set])
        
        # Initialize labels and data containers
        all_fps = []
        set_labels = []
        molecule_ids = []
        
        # Add target set data
        all_fps.append(target_fps)
        set_labels.extend(['target'] * len(target_fps))
        molecule_ids.extend([getattr(mol, 'name', f"Target_{i+1}") 
                           for i, mol in enumerate(self.target_set)])
        
        # Add reference set if available and requested
        if include_reference and self.reference_set is not None:
            reference_fps = np.array([mol.get_fingerprint(fp_type=fp_type, as_numpy=True) 
                                    for mol in self.reference_set])
            all_fps.append(reference_fps)
            set_labels.extend(['reference'] * len(reference_fps))
            molecule_ids.extend([getattr(mol, 'name', f"Reference_{i+1}")
                               for i, mol in enumerate(self.reference_set)])
        
        # Combine fingerprints
        combined_fps = np.vstack(all_fps)
        
        # Apply standardization (important for PCA)
        scaler = StandardScaler()
        scaled_fps = scaler.fit_transform(combined_fps)
        
        # Perform PCA
        pca = PCA(n_components=n_components)
        pca_result = pca.fit_transform(scaled_fps)
        
        # Create a DataFrame with the results
        pca_df = pd.DataFrame(data=pca_result, 
                             columns=[f'PC{i+1}' for i in range(n_components)])
        pca_df['set'] = set_labels
        pca_df['molecule_id'] = molecule_ids
        
        # Calculate explained variance
        explained_variance = pca.explained_variance_ratio_ * 100
        
        # Calculate chemical space overlap
        if include_reference and self.reference_set is not None:
            # Split PCA results by set
            target_pca = pca_df[pca_df['set'] == 'target'].iloc[:, :n_components].values
            reference_pca = pca_df[pca_df['set'] == 'reference'].iloc[:, :n_components].values
            
            # Calculate convex hull area (for 2D only)
            if n_components == 2:
                try:
                    from scipy.spatial import ConvexHull
                    
                    # Calculate convex hulls
                    if len(target_pca) >= 3:
                        target_hull = ConvexHull(target_pca)
                        target_area = target_hull.volume
                    else:
                        target_area = 0
                        
                    if len(reference_pca) >= 3:
                        reference_hull = ConvexHull(reference_pca)
                        reference_area = reference_hull.volume
                    else:
                        reference_area = 0
                        
                    # Calculate overlap (approximation)
                    # This is a crude approximation by combining points and comparing areas
                    if len(target_pca) >= 3 and len(reference_pca) >= 3:
                        combined_hull = ConvexHull(np.vstack([target_pca, reference_pca]))
                        combined_area = combined_hull.volume
                        
                        # Estimate overlap
                        if combined_area > 0:
                            overlap_area = target_area + reference_area - combined_area
                            overlap_ratio = max(0, overlap_area) / min(target_area, reference_area) if min(target_area, reference_area) > 0 else 0
                        else:
                            overlap_area = 0
                            overlap_ratio = 0
                    else:
                        overlap_area = 0
                        overlap_ratio = 0
                        
                    # Add to results
                    hull_info = {
                        'target_area': target_area,
                        'reference_area': reference_area,
                        'combined_area': combined_area if 'combined_area' in locals() else None,
                        'overlap_area': overlap_area,
                        'overlap_ratio': overlap_ratio
                    }
                except Exception as e:
                    hull_info = {'error': str(e)}
            else:
                hull_info = {'error': 'Convex hull analysis only available for 2D PCA'}
                
            # Calculate nearest neighbor statistics
            try:
                from sklearn.neighbors import NearestNeighbors
                
                # Target -> Reference
                nn_t2r = NearestNeighbors(n_neighbors=1)
                nn_t2r.fit(reference_pca)
                distances_t2r, _ = nn_t2r.kneighbors(target_pca)
                
                # Reference -> Target
                nn_r2t = NearestNeighbors(n_neighbors=1)
                nn_r2t.fit(target_pca)
                distances_r2t, _ = nn_r2t.kneighbors(reference_pca)
                
                # Calculate statistics
                nn_stats = {
                    'target_to_reference_mean_distance': np.mean(distances_t2r),
                    'target_to_reference_median_distance': np.median(distances_t2r),
                    'reference_to_target_mean_distance': np.mean(distances_r2t),
                    'reference_to_target_median_distance': np.median(distances_r2t),
                    'bidirectional_mean_distance': np.mean(np.concatenate([distances_t2r, distances_r2t]))
                }
            except Exception as e:
                nn_stats = {'error': str(e)}
        else:
            hull_info = {}
            nn_stats = {}
        
        # Prepare results
        results = {
            'pca_model': pca,
            'scaler': scaler,
            'explained_variance_ratio': pca.explained_variance_ratio_,
            'explained_variance_percentage': explained_variance,
            'cumulative_variance': np.cumsum(explained_variance),
            'projections': pca_df,
            'hull_analysis': hull_info if include_reference and self.reference_set is not None else {},
            'nearest_neighbor_stats': nn_stats if include_reference and self.reference_set is not None else {}
        }
        
        # Cache results
        self._cache[cache_key] = results
        
        return results
        
    def visualize_chemical_space(self,
                                fp_type: str = 'morgan',
                                color_by: str = 'set',
                                properties: List[str] = None,
                                figsize: tuple = (12, 10),
                                marker_size: int = 60,
                                alpha: float = 0.7,
                                include_reference: bool = True,
                                add_hulls: bool = True,
                                add_density: bool = True) -> tuple:
        """
        Visualize chemical space using PCA projection of molecular fingerprints.
        
        Parameters
        ----------
        fp_type : str, default='morgan'
            Type of fingerprint to use for the chemical space analysis
        color_by : str, default='set'
            How to color points: 'set' (target/reference) or property name
        properties : List[str], optional
            List of properties to include as point labels (hover text)
        figsize : tuple, default=(12, 10)
            Figure size (width, height)
        marker_size : int, default=60
            Size of scatter plot markers
        alpha : float, default=0.7
            Transparency of points
        include_reference : bool, default=True
            Whether to include reference set in the visualization
        add_hulls : bool, default=True
            Whether to add convex hulls around each set
        add_density : bool, default=True
            Whether to add density contours
            
        Returns
        -------
        tuple
            (fig, ax) - matplotlib figure and axis objects
            
        Examples
        --------
        >>> evaluator = MoleculeSetEvaluator(target_set, reference_set)
        >>> fig, ax = evaluator.visualize_chemical_space(color_by='LogP')
        >>> plt.show()
        """
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
            from matplotlib.colors import LinearSegmentedColormap
            from scipy.spatial import ConvexHull
            from scipy.stats import gaussian_kde
        except ImportError:
            raise ImportError("Required libraries not available. Install with 'pip install matplotlib seaborn scipy'")
            
        # Get PCA results
        pca_results = self.calculate_chemical_space_pca(
            fp_type=fp_type, 
            n_components=2, 
            include_reference=include_reference
        )
        
        # Get projections
        projections = pca_results['projections']
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Color by set or property
        if color_by == 'set':
            # Set color palette
            palette = {'target': 'blue', 'reference': 'red'}
            
            # Plot each set
            for set_name, set_data in projections.groupby('set'):
                ax.scatter(
                    set_data['PC1'], 
                    set_data['PC2'],
                    s=marker_size,
                    alpha=alpha,
                    label=set_name,
                    c=palette.get(set_name, 'gray')
                )
                
            # Add convex hulls if requested
            if add_hulls:
                for set_name, set_data in projections.groupby('set'):
                    points = set_data[['PC1', 'PC2']].values
                    if len(points) >= 3:  # Need at least 3 points for a hull
                        hull = ConvexHull(points)
                        hull_color = palette.get(set_name, 'gray')
                        for simplex in hull.simplices:
                            ax.plot(points[simplex, 0], points[simplex, 1], 
                                  '-', color=hull_color, alpha=0.5)
                            
            # Add density contours if requested
            if add_density:
                for set_name, set_data in projections.groupby('set'):
                    points = set_data[['PC1', 'PC2']].values
                    if len(points) >= 10:  # Need enough points for density estimation
                        x, y = points[:, 0], points[:, 1]
                        xy = np.vstack([x, y])
                        
                        # Calculate kernel density estimate
                        try:
                            kde = gaussian_kde(xy)
                            
                            # Create grid for contour plot
                            x_grid = np.linspace(x.min(), x.max(), 100)
                            y_grid = np.linspace(y.min(), y.max(), 100)
                            X, Y = np.meshgrid(x_grid, y_grid)
                            positions = np.vstack([X.ravel(), Y.ravel()])
                            
                            # Evaluate KDE at grid positions
                            Z = kde(positions).reshape(X.shape)
                            
                            # Plot contours
                            contour_color = palette.get(set_name, 'gray')
                            ax.contour(X, Y, Z, levels=5, colors=contour_color, alpha=0.5)
                        except:
                            # Skip density contours if estimation fails
                            pass
            
        else:
            # Color by property - first get property values
            property_df = self.calculate_property_distributions([color_by])
            
            if include_reference and self.reference_set is not None:
                # Get reference set properties
                temp_evaluator = MoleculeSetEvaluator(self.reference_set)
                ref_property_df = temp_evaluator.calculate_property_distributions([color_by])
                
                # Map molecule IDs to properties
                property_map = {**{row['MoleculeID']: row[color_by] 
                                 for _, row in property_df.iterrows() if color_by in row},
                               **{row['MoleculeID']: row[color_by]
                                 for _, row in ref_property_df.iterrows() if color_by in row}}
            else:
                # Map molecule IDs to properties (target set only)
                property_map = {row['MoleculeID']: row[color_by] 
                              for _, row in property_df.iterrows() if color_by in row}
            
            # Add property values to projections
            projections['property_value'] = projections['molecule_id'].map(property_map)
            
            # Create scatter plot colored by property
            scatter = ax.scatter(
                projections['PC1'],
                projections['PC2'],
                c=projections['property_value'],
                s=marker_size,
                alpha=alpha,
                cmap='viridis'
            )
            
            # Add colorbar
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.set_label(color_by)
        
        # Set labels and title
        explained_var = pca_results['explained_variance_percentage']
        ax.set_xlabel(f'PC1 ({explained_var[0]:.1f}%)', fontsize=12)
        ax.set_ylabel(f'PC2 ({explained_var[1]:.1f}%)', fontsize=12)
        
        title = f'Chemical Space Analysis: {fp_type} Fingerprints'
        if color_by != 'set':
            title += f' (colored by {color_by})'
        ax.set_title(title, fontsize=14)
        
        # Add legend if coloring by set
        if color_by == 'set':
            ax.legend()
            
        # Add grid
        ax.grid(True, linestyle='--', alpha=0.3)
        
        # Add annotation with variance explanation
        var_text = f'Total explained variance: {sum(explained_var[:2]):.1f}%'
        ax.text(0.02, 0.98, var_text, transform=ax.transAxes,
              verticalalignment='top', horizontalalignment='left',
              bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Add hull analysis stats if available
        if include_reference and add_hulls and 'hull_analysis' in pca_results and pca_results['hull_analysis']:
            hull_info = pca_results['hull_analysis']
            if 'overlap_ratio' in hull_info:
                overlap_text = f"Overlap ratio: {hull_info['overlap_ratio']:.2f}"
                ax.text(0.02, 0.93, overlap_text, transform=ax.transAxes,
                      verticalalignment='top', horizontalalignment='left',
                      bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Add nearest neighbor stats if available
        if include_reference and 'nearest_neighbor_stats' in pca_results and pca_results['nearest_neighbor_stats']:
            nn_stats = pca_results['nearest_neighbor_stats']
            if 'bidirectional_mean_distance' in nn_stats:
                nn_text = f"Mean NN distance: {nn_stats['bidirectional_mean_distance']:.3f}"
                ax.text(0.02, 0.88, nn_text, transform=ax.transAxes,
                      verticalalignment='top', horizontalalignment='left',
                      bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                
        # Tight layout
        plt.tight_layout()
        
        return fig, ax
    
    def calculate_chemical_space_tsne(self, 
                                fp_type: str = 'morgan', 
                                n_components: int = 2,
                                perplexity: float = 30.0,
                                early_exaggeration: float = 12.0,
                                learning_rate: float = 200.0,
                                n_iter: int = 1000,
                                include_reference: bool = True) -> Dict[str, Any]:
        """
        Calculate t-SNE projection of chemical space based on molecular fingerprints.
        
        Parameters
        ----------
        fp_type : str, default='morgan'
            Type of fingerprint to use for the chemical space analysis
        n_components : int, default=2
            Number of dimensions in the embedded space
        perplexity : float, default=30.0
            Related to the number of nearest neighbors used in t-SNE
        early_exaggeration : float, default=12.0
            Controls how tight clusters are in the embedding space
        learning_rate : float, default=200.0
            The learning rate for t-SNE
        n_iter : int, default=1000
            Maximum number of iterations for the optimization
        include_reference : bool, default=True
            Whether to include reference set in the analysis
            
        Returns
        -------
        Dict[str, Any]
            Dictionary containing t-SNE results and projections
            
        Examples
        --------
        >>> evaluator = MoleculeSetEvaluator(target_set, reference_set)
        >>> tsne_results = evaluator.calculate_chemical_space_tsne()
        >>> print(f"KL divergence: {tsne_results['kl_divergence']}")
        """
        try:
            from sklearn.manifold import TSNE
            from sklearn.preprocessing import StandardScaler
        except ImportError:
            raise ImportError("Scikit-learn is required for t-SNE analysis. Install with 'pip install scikit-learn'")
            
        # Check cache
        cache_key = f"chemical_space_tsne_{fp_type}_{n_components}_{perplexity}_{n_iter}_{include_reference}"
        if cache_key in self._cache:
            return self._cache[cache_key]
            
        # Get fingerprints for target set
        target_fps = np.array([mol.get_fingerprint(fp_type=fp_type, as_numpy=True) 
                             for mol in self.target_set])
        
        # Initialize labels and data containers
        all_fps = []
        set_labels = []
        molecule_ids = []
        
        # Add target set data
        all_fps.append(target_fps)
        set_labels.extend(['target'] * len(target_fps))
        molecule_ids.extend([getattr(mol, 'name', f"Target_{i+1}") 
                           for i, mol in enumerate(self.target_set)])
        
        # Add reference set if available and requested
        if include_reference and self.reference_set is not None:
            reference_fps = np.array([mol.get_fingerprint(fp_type=fp_type, as_numpy=True) 
                                    for mol in self.reference_set])
            all_fps.append(reference_fps)
            set_labels.extend(['reference'] * len(reference_fps))
            molecule_ids.extend([getattr(mol, 'name', f"Reference_{i+1}")
                               for i, mol in enumerate(self.reference_set)])
        
        # Combine fingerprints
        combined_fps = np.vstack(all_fps)
        
        # Apply standardization (can help t-SNE performance)
        scaler = StandardScaler()
        scaled_fps = scaler.fit_transform(combined_fps)
        
        # Perform t-SNE
        tsne = TSNE(
            n_components=n_components,
            perplexity=perplexity,
            early_exaggeration=early_exaggeration,
            learning_rate=learning_rate,
            max_iter=n_iter,
            metric='cosine',  # Cosine is often good for fingerprints
            init='pca',       # Initialize with PCA (more stable)
            random_state=42   # For reproducibility
        )
        
        tsne_result = tsne.fit_transform(scaled_fps)
        
        # Create a DataFrame with the results
        tsne_df = pd.DataFrame(data=tsne_result, 
                              columns=[f'TSNE{i+1}' for i in range(n_components)])
        tsne_df['set'] = set_labels
        tsne_df['molecule_id'] = molecule_ids
        
        # Calculate nearest neighbor statistics if reference set is included
        if include_reference and self.reference_set is not None:
            # Split t-SNE results by set
            target_tsne = tsne_df[tsne_df['set'] == 'target'].iloc[:, :n_components].values
            reference_tsne = tsne_df[tsne_df['set'] == 'reference'].iloc[:, :n_components].values
            
            # Calculate nearest neighbor statistics
            try:
                from sklearn.neighbors import NearestNeighbors
                
                # Target -> Reference
                nn_t2r = NearestNeighbors(n_neighbors=1)
                nn_t2r.fit(reference_tsne)
                distances_t2r, indices_t2r = nn_t2r.kneighbors(target_tsne)
                
                # Reference -> Target
                nn_r2t = NearestNeighbors(n_neighbors=1)
                nn_r2t.fit(target_tsne)
                distances_r2t, indices_r2t = nn_r2t.kneighbors(reference_tsne)
                
                # Calculate statistics
                nn_stats = {
                    'target_to_reference_mean_distance': np.mean(distances_t2r),
                    'target_to_reference_median_distance': np.median(distances_t2r),
                    'reference_to_target_mean_distance': np.mean(distances_r2t),
                    'reference_to_target_median_distance': np.median(distances_r2t),
                    'bidirectional_mean_distance': np.mean(np.concatenate([distances_t2r, distances_r2t]))
                }
                
                # Map each target molecule to its nearest reference
                nearest_references = []
                for i, idx in enumerate(indices_t2r.flatten()):
                    target_id = tsne_df[tsne_df['set'] == 'target']['molecule_id'].iloc[i]
                    reference_id = tsne_df[tsne_df['set'] == 'reference']['molecule_id'].iloc[idx]
                    distance = distances_t2r[i][0]
                    nearest_references.append({
                        'target_id': target_id,
                        'reference_id': reference_id,
                        'distance': distance
                    })
                
                # Map each reference molecule to its nearest target
                nearest_targets = []
                for i, idx in enumerate(indices_r2t.flatten()):
                    reference_id = tsne_df[tsne_df['set'] == 'reference']['molecule_id'].iloc[i]
                    target_id = tsne_df[tsne_df['set'] == 'target']['molecule_id'].iloc[idx]
                    distance = distances_r2t[i][0]
                    nearest_targets.append({
                        'reference_id': reference_id,
                        'target_id': target_id,
                        'distance': distance
                    })
                
                # Add nearest neighbor mappings to stats
                nn_stats['nearest_references'] = nearest_references
                nn_stats['nearest_targets'] = nearest_targets
                
            except Exception as e:
                nn_stats = {'error': str(e)}
        else:
            nn_stats = {}
        
        # Prepare results
        results = {
            'tsne_model': tsne,
            'scaler': scaler,
            'kl_divergence': tsne.kl_divergence_,
            'n_iter': tsne.n_iter_,
            'projections': tsne_df,
            'nearest_neighbor_stats': nn_stats if include_reference and self.reference_set is not None else {},
            'parameters': {
                'perplexity': perplexity,
                'early_exaggeration': early_exaggeration,
                'learning_rate': learning_rate,
                'n_iter': n_iter
            }
        }
        
        # Cache results
        self._cache[cache_key] = results
        
        return results
        
    def visualize_chemical_space_tsne(self,
                                     fp_type: str = 'morgan',
                                     color_by: str = 'set',
                                     properties: List[str] = None,
                                     figsize: tuple = (12, 10),
                                     marker_size: int = 60,
                                     alpha: float = 0.7,
                                     include_reference: bool = True,
                                     perplexity: float = 30.0,
                                     add_density: bool = True) -> tuple:
        """
        Visualize chemical space using t-SNE projection of molecular fingerprints.
        
        Parameters
        ----------
        fp_type : str, default='morgan'
            Type of fingerprint to use for the chemical space analysis
        color_by : str, default='set'
            How to color points: 'set' (target/reference) or property name
        properties : List[str], optional
            List of properties to include as point labels (hover text)
        figsize : tuple, default=(12, 10)
            Figure size (width, height)
        marker_size : int, default=60
            Size of scatter plot markers
        alpha : float, default=0.7
            Transparency of points
        include_reference : bool, default=True
            Whether to include reference set in the visualization
        perplexity : float, default=30.0
            Perplexity parameter for t-SNE
        add_density : bool, default=True
            Whether to add density contours
            
        Returns
        -------
        tuple
            (fig, ax) - matplotlib figure and axis objects
            
        Examples
        --------
        >>> evaluator = MoleculeSetEvaluator(target_set, reference_set)
        >>> fig, ax = evaluator.visualize_chemical_space_tsne(color_by='LogP')
        >>> plt.show()
        """
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
            from matplotlib.colors import LinearSegmentedColormap
            from scipy.stats import gaussian_kde
        except ImportError:
            raise ImportError("Required libraries not available. Install with 'pip install matplotlib seaborn scipy'")
            
        # Get t-SNE results
        tsne_results = self.calculate_chemical_space_tsne(
            fp_type=fp_type, 
            n_components=2, 
            perplexity=perplexity,
            include_reference=include_reference
        )
        
        # Get projections
        projections = tsne_results['projections']
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Color by set or property
        if color_by == 'set':
            # Set color palette
            palette = {'target': 'blue', 'reference': 'red'}
            
            # Plot each set
            for set_name, set_data in projections.groupby('set'):
                ax.scatter(
                    set_data['TSNE1'], 
                    set_data['TSNE2'],
                    s=marker_size,
                    alpha=alpha,
                    label=set_name,
                    c=palette.get(set_name, 'gray')
                )
                            
            # Add density contours if requested
            if add_density:
                for set_name, set_data in projections.groupby('set'):
                    points = set_data[['TSNE1', 'TSNE2']].values
                    if len(points) >= 10:  # Need enough points for density estimation
                        x, y = points[:, 0], points[:, 1]
                        xy = np.vstack([x, y])
                        
                        # Calculate kernel density estimate
                        try:
                            kde = gaussian_kde(xy)
                            
                            # Create grid for contour plot
                            x_grid = np.linspace(x.min(), x.max(), 100)
                            y_grid = np.linspace(y.min(), y.max(), 100)
                            X, Y = np.meshgrid(x_grid, y_grid)
                            positions = np.vstack([X.ravel(), Y.ravel()])
                            
                            # Evaluate KDE at grid positions
                            Z = kde(positions).reshape(X.shape)
                            
                            # Plot contours
                            contour_color = palette.get(set_name, 'gray')
                            ax.contour(X, Y, Z, levels=5, colors=contour_color, alpha=0.5)
                        except:
                            # Skip density contours if estimation fails
                            pass
                
        else:
            # Color by property - first get property values
            property_df = self.calculate_property_distributions([color_by])
            
            if include_reference and self.reference_set is not None:
                # Get reference set properties
                temp_evaluator = MoleculeSetEvaluator(self.reference_set)
                ref_property_df = temp_evaluator.calculate_property_distributions([color_by])
                
                # Map molecule IDs to properties
                property_map = {**{row['MoleculeID']: row[color_by] 
                                 for _, row in property_df.iterrows() if color_by in row},
                               **{row['MoleculeID']: row[color_by]
                                 for _, row in ref_property_df.iterrows() if color_by in row}}
            else:
                # Map molecule IDs to properties (target set only)
                property_map = {row['MoleculeID']: row[color_by] 
                              for _, row in property_df.iterrows() if color_by in row}
            
            # Add property values to projections
            projections['property_value'] = projections['molecule_id'].map(property_map)
            
            # Create scatter plot colored by property
            scatter = ax.scatter(
                projections['TSNE1'],
                projections['TSNE2'],
                c=projections['property_value'],
                s=marker_size,
                alpha=alpha,
                cmap='viridis'
            )
            
            # Add colorbar
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.set_label(color_by)
        
        # Set labels and title
        ax.set_xlabel('t-SNE dimension 1', fontsize=12)
        ax.set_ylabel('t-SNE dimension 2', fontsize=12)
        
        title = f't-SNE Chemical Space Analysis: {fp_type} Fingerprints'
        if color_by != 'set':
            title += f' (colored by {color_by})'
        ax.set_title(title, fontsize=14)
        
        # Add legend if coloring by set
        if color_by == 'set':
            ax.legend()
            
        # Add grid
        ax.grid(True, linestyle='--', alpha=0.3)
        
        # Add annotation with t-SNE parameters
        params_text = f'Perplexity: {tsne_results["parameters"]["perplexity"]}, '
        params_text += f'KL divergence: {tsne_results["kl_divergence"]:.2f}'
        ax.text(0.02, 0.98, params_text, transform=ax.transAxes,
              verticalalignment='top', horizontalalignment='left',
              bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Add nearest neighbor stats if available
        if include_reference and 'nearest_neighbor_stats' in tsne_results and tsne_results['nearest_neighbor_stats']:
            nn_stats = tsne_results['nearest_neighbor_stats']
            if 'bidirectional_mean_distance' in nn_stats:
                nn_text = f"Mean NN distance: {nn_stats['bidirectional_mean_distance']:.3f}"
                ax.text(0.02, 0.93, nn_text, transform=ax.transAxes,
                      verticalalignment='top', horizontalalignment='left',
                      bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                
        # Tight layout
        plt.tight_layout()
        
        return fig, ax
        
    def visualize_combined_chemical_space(self,
                                         fp_type: str = 'morgan',
                                         color_by: str = 'set',
                                         figsize: tuple = (16, 8),
                                         marker_size: int = 60,
                                         alpha: float = 0.7,
                                         perplexity: float = 30.0,
                                         include_reference: bool = True) -> tuple:
        """
        Create a side-by-side comparison of PCA and t-SNE chemical space visualizations.
        
        Parameters
        ----------
        fp_type : str, default='morgan'
            Type of fingerprint to use for the analysis
        color_by : str, default='set'
            How to color points: 'set' (target/reference) or property name
        figsize : tuple, default=(16, 8)
            Figure size (width, height)
        marker_size : int, default=60
            Size of scatter plot markers
        alpha : float, default=0.7
            Transparency of points
        include_reference : bool, default=True
            Whether to include reference set in the visualization
            
        Returns
        -------
        tuple
            (fig, (ax1, ax2)) - matplotlib figure and axes objects
            
        Examples
        --------
        >>> evaluator = MoleculeSetEvaluator(target_set, reference_set)
        >>> fig, (ax1, ax2) = evaluator.visualize_combined_chemical_space()
        >>> plt.show()
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("Matplotlib is required for visualization. Install with 'pip install matplotlib'")
            
        # Create figure with two subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
        
        # Get PCA and t-SNE results
        pca_results = self.calculate_chemical_space_pca(
            fp_type=fp_type, 
            include_reference=include_reference
        )
        
        tsne_results = self.calculate_chemical_space_tsne(
            fp_type=fp_type,
            include_reference=include_reference,
            perplexity=perplexity
        )
        
        # Get projections
        pca_projections = pca_results['projections']
        tsne_projections = tsne_results['projections']
        
        # Color by set or property
        if color_by == 'set':
            # Set color palette
            palette = {'target': 'blue', 'reference': 'red'}
            
            # Plot PCA
            for set_name, set_data in pca_projections.groupby('set'):
                ax1.scatter(
                    set_data['PC1'], 
                    set_data['PC2'],
                    s=marker_size,
                    alpha=alpha,
                    label=set_name,
                    c=palette.get(set_name, 'gray')
                )
                
            # Plot t-SNE
            for set_name, set_data in tsne_projections.groupby('set'):
                ax2.scatter(
                    set_data['TSNE1'], 
                    set_data['TSNE2'],
                    s=marker_size,
                    alpha=alpha,
                    label=set_name,
                    c=palette.get(set_name, 'gray')
                )
                
            # Add legends
            ax1.legend()
            ax2.legend()
            
        else:
            # Color by property - first get property values
            property_df = self.calculate_property_distributions([color_by])
            
            if include_reference and self.reference_set is not None:
                # Get reference set properties
                temp_evaluator = MoleculeSetEvaluator(self.reference_set)
                ref_property_df = temp_evaluator.calculate_property_distributions([color_by])
                
                # Map molecule IDs to properties
                property_map = {**{row['MoleculeID']: row[color_by] 
                                 for _, row in property_df.iterrows() if color_by in row},
                               **{row['MoleculeID']: row[color_by]
                                 for _, row in ref_property_df.iterrows() if color_by in row}}
            else:
                # Map molecule IDs to properties (target set only)
                property_map = {row['MoleculeID']: row[color_by] 
                              for _, row in property_df.iterrows() if color_by in row}
            
            # Add property values to projections
            pca_projections['property_value'] = pca_projections['molecule_id'].map(property_map)
            tsne_projections['property_value'] = tsne_projections['molecule_id'].map(property_map)
            
            # Find global min/max for consistent color scale
            vmin = min(pca_projections['property_value'].min(), tsne_projections['property_value'].min())
            vmax = max(pca_projections['property_value'].max(), tsne_projections['property_value'].max())
            
            # Plot PCA
            scatter1 = ax1.scatter(
                pca_projections['PC1'],
                pca_projections['PC2'],
                c=pca_projections['property_value'],
                s=marker_size,
                alpha=alpha,
                cmap='viridis',
                vmin=vmin,
                vmax=vmax
            )
            
            # Plot t-SNE
            scatter2 = ax2.scatter(
                tsne_projections['TSNE1'],
                tsne_projections['TSNE2'],
                c=tsne_projections['property_value'],
                s=marker_size,
                alpha=alpha,
                cmap='viridis',
                vmin=vmin,
                vmax=vmax
            )
            
            # Add colorbar
            cbar = fig.colorbar(scatter1, ax=[ax1, ax2])
            cbar.set_label(color_by)
        
        # Set labels and titles
        explained_var = pca_results['explained_variance_percentage']
        ax1.set_xlabel(f'PC1 ({explained_var[0]:.1f}%)', fontsize=12)
        ax1.set_ylabel(f'PC2 ({explained_var[1]:.1f}%)', fontsize=12)
        ax1.set_title('PCA Chemical Space', fontsize=14)
        
        ax2.set_xlabel('t-SNE dimension 1', fontsize=12)
        ax2.set_ylabel('t-SNE dimension 2', fontsize=12)
        ax2.set_title('t-SNE Chemical Space', fontsize=14)
        
        # Add grids
        ax1.grid(True, linestyle='--', alpha=0.3)
        ax2.grid(True, linestyle='--', alpha=0.3)
        
        # Set overall title
        overall_title = f'Chemical Space Analysis: {fp_type} Fingerprints'
        if color_by != 'set':
            overall_title += f' (colored by {color_by})'
        fig.suptitle(overall_title, fontsize=16)
        
        # Adjust layout
        #plt.tight_layout()
        fig.subplots_adjust(top=0.9)
        
        return fig, (ax1, ax2)