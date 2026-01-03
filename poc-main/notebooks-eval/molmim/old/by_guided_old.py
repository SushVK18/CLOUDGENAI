# qsar_scoring.py
import numpy as np
import boto3
from botocore.exceptions import ClientError
import json
import os
import atexit
from msal import ConfidentialClientApplication, SerializableTokenCache
import requests
from typing import List, Dict, Union, Tuple, Optional
import time
from rdkit import Chem
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity
from rdkit.Chem import QED
from rdkit.Chem.Crippen import MolLogP

class QSARScorer:
    """
    A class to handle QSAR-based scoring for molecule optimization.
    Combines QSAR predictions with other molecular properties.
    """
    
    TOKEN_FILE = os.path.join(os.path.expanduser('~'), ".insilicoclitoken")
    
    def __init__(self, 
                 secret_name: str = "PROPCALC_API_CREDS",
                 tenant_id: str = "fcb2b37b-5da0-466b-9b83-0014b67a7c78",
                 qsar_weight: float = 1.0, 
                 logp_weight: float = 0.5,
                 similarity_weight: float = 0.5,
                 qed_weight: float = 0.5,
                 qsar_threshold: float = 0.7,
                 logp_threshold: float = 0.9,
                 similarity_threshold: float = 0.4,
                 qed_threshold: float = 1.0,
                 qsar_endpoint: str = "https://chemprop.cropkey-np.ag/activity/ica/invocations_hsec",
                 reference_smiles: Optional[str] = None):
        """
        Initialize the QSAR scorer with configurable weights and thresholds.
        
        Args:
            secret_name: AWS secret name containing API credentials
            tenant_id: Microsoft tenant ID
            qsar_weight: Weight for QSAR predictions in the overall score
            logp_weight: Weight for logP contribution in the overall score
            similarity_weight: Weight for similarity contribution in the overall score
            qed_weight: Weight for QED (drug-likeness) contribution in the overall score
            qsar_threshold: Threshold for normalizing QSAR predictions
            logp_threshold: Threshold for normalizing logP values
            similarity_threshold: Threshold for normalizing similarity values
            qed_threshold: Threshold for normalizing QED values (typically 1.0 as QED is 0-1)
            qsar_endpoint: API endpoint for QSAR predictions
            reference_smiles: Reference SMILES for calculating similarity (if None, similarity must be pre-calculated)
        """
        self.secret_name = secret_name
        self.tenant_id = tenant_id
        self.qsar_weight = qsar_weight
        self.logp_weight = logp_weight
        self.similarity_weight = similarity_weight
        self.qed_weight = qed_weight
        self.qsar_threshold = qsar_threshold
        self.logp_threshold = logp_threshold
        self.similarity_threshold = similarity_threshold
        self.qed_threshold = qed_threshold
        self.qsar_endpoint = qsar_endpoint
        self.reference_smiles = reference_smiles
        
        # Initialize reference molecule for similarity calculations if provided
        self.reference_mol = None
        self.reference_fingerprint = None
        if reference_smiles:
            self._initialize_reference_fingerprint(reference_smiles)
            
        # Get API credentials
        self._initialize_api_connection()
    
    def _initialize_reference_fingerprint(self, reference_smiles: str):
        """Initialize reference molecule fingerprint for similarity calculations."""
        try:
            # Morgan fingerprint parameters
            self.fingerprint_radius = 2
            self.fingerprint_nbits = 2048
            
            # Generate reference fingerprint
            self.reference_mol = Chem.MolFromSmiles(reference_smiles)
            if self.reference_mol is None:
                print(f"Warning: Could not parse reference SMILES: {reference_smiles}")
                return
                
            self.reference_fingerprint = GetMorganFingerprintAsBitVect(
                self.reference_mol, 
                radius=self.fingerprint_radius, 
                nBits=self.fingerprint_nbits
            )
        except Exception as e:
            print(f"Error initializing reference fingerprint: {e}")
            self.reference_mol = None
            self.reference_fingerprint = None
            
    def _initialize_api_connection(self):
        """Initialize the API connection and authentication."""
        # Get API credentials from AWS Secrets Manager
        secret = self._get_secret(self.secret_name)
        self.client_id = secret.get("CLIENT_ID")
        self.client_secret = secret.get("CLIENT_SECRET")
        
        # Initialize token cache
        self.cache = SerializableTokenCache()
        if os.path.exists(self.TOKEN_FILE):
            with open(self.TOKEN_FILE, 'r') as f:
                self.cache.deserialize(f.read())
        
        # Register cache write on exit
        atexit.register(lambda: self._write_cache(self.cache, self.TOKEN_FILE))
        
        # Initialize MSAL application
        self.app = ConfidentialClientApplication(
            client_id=self.client_id,
            client_credential=self.client_secret,
            authority=f"https://login.microsoftonline.com/{self.tenant_id}",
            token_cache=self.cache
        )
    
    def _get_secret(self, secret_name: str) -> Dict:
        """Retrieve secret from AWS Secrets Manager."""
        region_name = "us-east-1"
        client = boto3.client(service_name="secretsmanager", region_name=region_name)
        try:
            get_secret_value_response = client.get_secret_value(SecretId=secret_name)
            secret = json.loads(get_secret_value_response["SecretString"])
            if secret is None:
                raise Exception("Secret is None")
            return secret
        except ClientError as e:
            print(f"Unable to retrieve secret: {e}")
            raise e
        except Exception as e:
            print(f"Unknown exception: Unable to retrieve secret: {e}")
            raise e
    
    def _write_cache(self, cache: SerializableTokenCache, file: str):
        """Write token cache to file if changed."""
        if cache.has_state_changed:
            with open(file, "w") as f:
                f.write(cache.serialize())
    
    def _post_data(self, bearer_token: str, url: str, data: Dict) -> Dict:
        """Post data to API endpoint."""
        headers = {
            'accept': 'application/json',
            'Content-Type': 'application/json',
            'Authorization': f'Bearer {bearer_token}'
        }
        response = requests.post(url, headers=headers, data=json.dumps(data))
        return response.json()
    
    def _get_access_token(self) -> str:
        """Get access token for API authentication."""
        result = self.app.acquire_token_for_client(scopes=[f"{self.client_id}/.default"])
        access_token = result.get("access_token")
        if not access_token:
            raise Exception("Failed to acquire access token")
        return access_token
    
    def calculate_tanimoto_similarity(self, smiles: str) -> float:
        """
        Calculate Tanimoto similarity between a SMILES string and the reference molecule.
        
        Args:
            smiles: SMILES string to compare with reference
            
        Returns:
            Tanimoto similarity score (0-1)
        """
        if self.reference_fingerprint is None:
            raise ValueError("Reference fingerprint not initialized. Provide a reference_smiles during initialization.")
            
        # Parse the molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 0.0
            
        # Generate fingerprint
        fingerprint = GetMorganFingerprintAsBitVect(
            mol, 
            radius=self.fingerprint_radius, 
            nBits=self.fingerprint_nbits
        )
        
        # Calculate similarity
        return TanimotoSimilarity(fingerprint, self.reference_fingerprint)
    
    def calculate_qed(self, smiles: str) -> float:
        """
        Calculate Quantitative Estimate of Drug-likeness (QED) for a molecule.
        
        Args:
            smiles: SMILES string of the molecule
            
        Returns:
            QED score (0-1, higher is more drug-like)
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 0.0
            
        try:
            return QED.qed(mol)
        except:
            return 0.0
    
    def get_qsar_predictions(self, smiles_list: List[str]) -> List[float]:
        """
        Get QSAR predictions for a list of SMILES strings.
        
        Args:
            smiles_list: List of SMILES strings representing molecules
            
        Returns:
            List of QSAR prediction scores
        """
        # Format input for API
        request_input = {
            "input": [
                {
                    "id": str(i+1),
                    "structure": smiles
                } for i, smiles in enumerate(smiles_list)
            ]
        }
        
        # Get access token
        access_token = self._get_access_token()
        
        # Make API request
        response = self._post_data(access_token, self.qsar_endpoint, request_input)
        
        # Process response
        predictions = []
        for result in response.get("result", []):
            # Handle potential errors in predictions
            if "errmsg" in result:
                predictions.append(0.0)  # Default low score for invalid structures
            else:
                # Get prediction value, defaulting to 0 if not found
                pred = result.get("prediction", {}).get("value", 0.0)
                predictions.append(float(pred))
        
        # Ensure we have a prediction for each input SMILES
        if len(predictions) != len(smiles_list):
            # Fill in missing predictions with zeros
            predictions.extend([0.0] * (len(smiles_list) - len(predictions)))
        
        return predictions
    
    def compute_scores(self, 
                      smiles_list: List[str], 
                      logp_scores: List[float], 
                      similarity_scores: Optional[List[float]] = None,
                      batch_size: int = 50) -> np.ndarray:
        """
        Compute combined scores using QSAR predictions, logP, similarity, and QED scores.
        
        Args:
            smiles_list: List of SMILES strings
            logp_scores: List of logP scores
            similarity_scores: List of similarity scores (optional if reference_smiles provided)
            batch_size: Batch size for API calls to avoid timeouts
            
        Returns:
            Array of combined scores
        """
        # Process QSAR predictions in batches to avoid timeouts
        all_qsar_scores = []
        for i in range(0, len(smiles_list), batch_size):
            batch_smiles = smiles_list[i:i+batch_size]
            batch_scores = self.get_qsar_predictions(batch_smiles)
            all_qsar_scores.extend(batch_scores)
            # Add a small delay between batches
            if i + batch_size < len(smiles_list):
                time.sleep(0.5)
        
        # Calculate similarity scores if not provided but reference SMILES is available
        if similarity_scores is None:
            if self.reference_mol is None:
                raise ValueError("No similarity scores provided and no reference SMILES initialized")
            similarity_scores = [self.calculate_tanimoto_similarity(smiles) for smiles in smiles_list]
        
        # Calculate QED scores
        qed_scores = [self.calculate_qed(smiles) for smiles in smiles_list]
        
        # Ensure all arrays have the same length
        n = len(smiles_list)
        qsar_array = np.array(all_qsar_scores[:n])
        logp_array = np.array(logp_scores[:n])
        similarity_array = np.array(similarity_scores[:n])
        qed_array = np.array(qed_scores[:n])
        
        # Normalize scores
        norm_qsar = np.clip(qsar_array / self.qsar_threshold, a_min=0.0, a_max=1.0)
        norm_logp = np.clip(logp_array / self.logp_threshold, a_min=0.0, a_max=1.0)
        norm_similarity = np.clip(similarity_array / self.similarity_threshold, a_min=0.0, a_max=1.0)
        norm_qed = np.clip(qed_array / self.qed_threshold, a_min=0.0, a_max=1.0)
        
        # Combine scores with weights (negative for minimization)
        combined_scores = -1.0 * (
            self.qsar_weight * norm_qsar + 
            self.logp_weight * norm_logp + 
            self.similarity_weight * norm_similarity +
            self.qed_weight * norm_qed
        )
        
        return combined_scores

    def score(self, 
             smiles_list: List[str], 
             logp_scores: List[float], 
             similarity_scores: Optional[List[float]] = None) -> np.ndarray:
        """
        Simplified scoring function interface compatible with optimization oracles.
        
        Args:
            smiles_list: List of SMILES strings
            logp_scores: List of logP scores (can be dummy values if logp_weight=0)
            similarity_scores: List of similarity scores (optional if reference_smiles provided)
            
        Returns:
            Array of scores (lower is better for minimization)
        """
        return self.compute_scores(smiles_list, logp_scores, similarity_scores)
        
    def score_without_logp(self,
                         smiles_list: List[str],
                         similarity_scores: Optional[List[float]] = None) -> np.ndarray:
        """
        Convenience method that scores without needing logP values.
        
        Args:
            smiles_list: List of SMILES strings
            similarity_scores: List of similarity scores (optional if reference_smiles provided)
            
        Returns:
            Array of scores (lower is better for minimization)
        """
        dummy_logp_scores = [0.0] * len(smiles_list)
        return self.compute_scores(smiles_list, dummy_logp_scores, similarity_scores)
        
    def score_qsar_similarity_qed(self,
                                smiles_list: List[str],
                                similarity_scores: Optional[List[float]] = None) -> np.ndarray:
        """
        Score using QSAR, similarity, and QED (no logP).
        
        Args:
            smiles_list: List of SMILES strings
            similarity_scores: List of similarity scores (optional if reference_smiles provided)
            
        Returns:
            Array of scores (lower is better for minimization)
        """
        # Create a temporary scorer with logP weight=0
        temp_scorer = QSARScorer(
            qsar_weight=self.qsar_weight,
            logp_weight=0.0,
            similarity_weight=self.similarity_weight,
            qed_weight=self.qed_weight,
            qsar_threshold=self.qsar_threshold,
            similarity_threshold=self.similarity_threshold,
            qed_threshold=self.qed_threshold,
            qsar_endpoint=self.qsar_endpoint,
            reference_smiles=self.reference_smiles
        )
        
        # Initialize API connection for the temporary scorer
        temp_scorer.client_id = self.client_id
        temp_scorer.client_secret = self.client_secret
        temp_scorer.app = self.app
        temp_scorer.cache = self.cache
        
        # Use dummy logP scores
        dummy_logp_scores = [0.0] * len(smiles_list)
        return temp_scorer.compute_scores(smiles_list, dummy_logp_scores, similarity_scores)


# Utility functions for direct import
def create_qsar_scorer(
    qsar_weight: float = 1.0,
    logp_weight: float = 0.5,
    similarity_weight: float = 0.5,
    qed_weight: float = 0.5,
    qsar_threshold: float = 0.7,
    logp_threshold: float = 0.9,
    similarity_threshold: float = 0.4,
    qed_threshold: float = 1.0,
    reference_smiles: Optional[str] = None
) -> QSARScorer:
    """
    Create a QSAR scorer with custom weights and thresholds.
    
    Args:
        qsar_weight: Weight for QSAR predictions
        logp_weight: Weight for logP contribution
        similarity_weight: Weight for similarity contribution
        qed_weight: Weight for QED (drug-likeness) contribution
        qsar_threshold: Threshold for normalizing QSAR predictions
        logp_threshold: Threshold for normalizing logP values
        similarity_threshold: Threshold for normalizing similarity values
        qed_threshold: Threshold for normalizing QED values
        reference_smiles: Reference molecule SMILES for similarity calculation
        
    Returns:
        QSARScorer instance
    """
    return QSARScorer(
        qsar_weight=qsar_weight,
        logp_weight=logp_weight,
        similarity_weight=similarity_weight,
        qed_weight=qed_weight,
        qsar_threshold=qsar_threshold,
        logp_threshold=logp_threshold,
        similarity_threshold=similarity_threshold,
        qed_threshold=qed_threshold,
        reference_smiles=reference_smiles
    )


def qsar_scoring_function(smiles_list: List[str], logp_scores: List[float], similarity_scores: List[float]) -> np.ndarray:
    """
    Direct drop-in replacement for the original scoring function that includes QSAR predictions.
    Uses the same input parameters as the original function for compatibility.
    
    Args:
        smiles_list: List of SMILES strings
        logp_scores: List of logP scores
        similarity_scores: List of similarity scores
        
    Returns:
        Array of scores (lower is better for minimization)
    """
    scorer = create_qsar_scorer()
    return scorer.score(smiles_list, logp_scores, similarity_scores)



def enhanced_scoring_function(smiles_list: List[str], 
                          reference_smiles: str,
                          include_qsar: bool = True,
                          include_logp: bool = True, 
                          include_similarity: bool = True,
                          include_qed: bool = True,
                          qsar_weight: float = 1.0,
                          logp_weight: float = 0.5,
                          similarity_weight: float = 0.5,
                          qed_weight: float = 0.5) -> np.ndarray:
    """
    Enhanced scoring function that calculates all needed properties from SMILES strings.
    
    Args:
        smiles_list: List of SMILES strings
        reference_smiles: Reference SMILES for similarity calculations
        include_qsar: Whether to include QSAR predictions
        include_logp: Whether to include LogP
        include_similarity: Whether to include similarity
        include_qed: Whether to include QED
        qsar_weight: Weight for QSAR component
        logp_weight: Weight for LogP component
        similarity_weight: Weight for similarity component
        qed_weight: Weight for QED component
        
    Returns:
        Array of scores (lower is better for minimization)
    """
    # Calculate logP if needed (otherwise use dummy values)
    if include_logp:
        logp_scores = [calculate_logp(smiles) for smiles in smiles_list]
        real_logp_weight = logp_weight
    else:
        logp_scores = [0.0] * len(smiles_list)
        real_logp_weight = 0.0
    
    # Calculate similarities if needed (otherwise use dummy values)
    if include_similarity:
        similarity_scores = [tanimoto_similarity(smiles, reference_smiles) for smiles in smiles_list]
        real_similarity_weight = similarity_weight
    else:
        similarity_scores = [0.0] * len(smiles_list)
        real_similarity_weight = 0.0
    
    # Set QSAR weight appropriately
    real_qsar_weight = qsar_weight if include_qsar else 0.0
    
    # Set QED weight appropriately
    real_qed_weight = qed_weight if include_qed else 0.0
    
    # Create scorer with appropriate weights
    scorer = create_qsar_scorer(
        qsar_weight=real_qsar_weight,
        logp_weight=real_logp_weight,
        similarity_weight=real_similarity_weight,
        qed_weight=real_qed_weight,
        reference_smiles=reference_smiles if include_similarity else None
    )
    
    # Calculate scores
    return scorer.score(smiles_list, logp_scores, similarity_scores)


def calculate_logp(smiles: str) -> float:
    """
    Calculate the LogP value for a given SMILES string using RDKit.
    
    Args:
        smiles: SMILES string of the molecule
        
    Returns:
        LogP value or 0.0 if molecule is invalid
    """
    mol = Chem.MolFromSmiles(smiles)
    # check if valid
    if mol is None:
        return 0.0
    
    #calculate LogP
    try:
        return MolLogP(mol)
    except:
        return 0.0


def tanimoto_similarity(smiles: str, reference: str) -> float:
    """
    Calculate Tanimoto similarity between a molecule and a reference molecule.
    Implementation matching the provided function for compatibility.
    
    Args:
        smiles: SMILES string of the molecule
        reference: Reference SMILES string
        
    Returns:
        Tanimoto similarity (0-1)
    """
    # Get fingerprint params
    fingerprint_radius_param = 2
    fingerprint_nbits = 2048
    
    # Handle the reference molecule
    reference_mol = Chem.MolFromSmiles(reference)
    if reference_mol is None:
        return 0.0
        
    reference_fingerprint = GetMorganFingerprintAsBitVect(
        reference_mol, radius=fingerprint_radius_param, nBits=fingerprint_nbits
    )
    
    # Validate the other molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0.0
        
    fingerprint = GetMorganFingerprintAsBitVect(
        mol, radius=fingerprint_radius_param, nBits=fingerprint_nbits
    )
    
    # Calculate and return the Tanimoto similarity
    return TanimotoSimilarity(fingerprint, reference_fingerprint)


# Legacy-compatible scoring function that behaves like the original
def original_scoring_function(logps, similarities) -> np.ndarray:
    """
    Original scoring function implementation for backward compatibility.
    
    Args:
        logps: List of logP values
        similarities: List of similarity values
        
    Returns:
        Array of scores
    """
    return -1. * (
        np.clip(np.array(logps) / 0.9, a_min=0.0, a_max=1.0) + 
        np.clip(np.array(similarities) / 0.4, a_min=0.0, a_max=1.0)
    )



# Example usage in a Jupyter notebook:
"""
from qsar_scoring import (
    tanimoto_similarity, 
    calculate_logp, 
    scoring_function,
    enhanced_scoring_function,
    qsar_scoring_function
)
from rdkit import Chem
from rdkit.Chem.Crippen import MolLogP

# Traditional usage pattern (drop-in replacement)
tanimoto_scores = []
logp_scores = []
for smile in molecules:
    mol = Chem.MolFromSmiles(smile)
    if mol:
        tanimoto = tanimoto_similarity(smile, m1)
        logp_score = MolLogP(mol)
        tanimoto_scores.append(tanimoto)
        logp_scores.append(logp_score)

# Original scoring function
original_scores = scoring_function(logp_scores, tanimoto_scores)

# QSAR-enhanced scoring function (drop-in replacement)
qsar_scores = qsar_scoring_function(molecules, logp_scores, tanimoto_scores)

# Enhanced scoring with QSAR, QED, similarity and logP
# (all calculated internally)
enhanced_scores = enhanced_scoring_function(
    molecules,
    reference_smiles=m1,
    include_qsar=True,
    include_qed=True
)

# Enhanced scoring without logP
without_logp_scores = enhanced_scoring_function(
    molecules,
    reference_smiles=m1,
    include_logp=False
)
"""