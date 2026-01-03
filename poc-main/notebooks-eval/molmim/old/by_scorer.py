# qsar_scoring.py
import numpy as np
import boto3
from botocore.exceptions import ClientError
import json
import os
import atexit
from msal import ConfidentialClientApplication, SerializableTokenCache
import requests
from typing import List, Dict, Union, Tuple, Optional, NamedTuple
import time
from rdkit import Chem
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity
from rdkit.Chem import QED
from rdkit.Chem.Crippen import MolLogP


class ScoringResults(NamedTuple):
    """Container for all scoring results"""
    combined_scores: np.ndarray     # Final combined scores
    logp_scores: List[float]        # Individual LogP scores
    similarity_scores: List[float]  # Individual similarity scores
    qsar_hsec_scores: List[float]   # Individual QSAR-HSEC prediction scores
    qsar_sip_scores: List[float]    # Individual QSAR-SIP prediction scores
    qed_scores: List[float]         # Individual QED scores


class MoleculeScorer:
    """
    A configurable scorer for molecular optimization that can combine different metrics:
    - QSAR predictions
    - Structural similarity
    - LogP
    - QED (drug-likeness)
    """
    
    TOKEN_FILE = os.path.join(os.path.expanduser('~'), ".insilicoclitoken")
    
    def __init__(self, 
                 reference_smiles: Optional[str] = None,
                 use_qsar_hsec: bool = True,
                 use_qsar_sip: bool = False, 
                 use_logp: bool = False,
                 use_similarity: bool = True,
                 use_qed: bool = False,
                 qsar_hsec_weight: float = 1.0,
                 qsar_sip_weight: float = 1.0,
                 logp_weight: float =  1.0,
                 similarity_weight: float = 1.0,
                 qed_weight: float = 1.0,
                 qsar_hsec_threshold: float = 0.9,
                 qsar_sip_threshold: float = 0.9,
                 logp_threshold: float = 0.9,
                 similarity_threshold: float = 0.4,
                 qed_threshold: float = 0.9,
                 secret_name: str = "PROPCALC_API_CREDS",
                 tenant_id: str = "fcb2b37b-5da0-466b-9b83-0014b67a7c78",
                 qsar_hsec_endpoint: str = "https://chemprop.cropkey-np.ag/activity/ica/invocations_hsec",
                 qsar_sip_endpoint: str = "https://chemprop.cropkey-np.ag/activity/ica/invocations"):
        """
        Initialize a configurable molecular scorer.
        
        Args:
            reference_smiles: Reference SMILES for similarity calculation
            use_qsar_hsec: Whether to include QSAR-HSEC predictions in scoring
            use_qsar_sip: Whether to include QSAR-SIP predictions in scoring
            use_logp: Whether to include LogP in scoring
            use_similarity: Whether to include structural similarity in scoring
            use_qed: Whether to include drug-likeness score in scoring
            qsar_hsec_weight: Weight for QSAR-HSEC predictions in the combined score
            qsar_sip_weight: Weight for QSAR-SIP predictions in the combined score
            logp_weight: Weight for LogP in the combined score
            similarity_weight: Weight for similarity in the combined score
            qed_weight: Weight for QED in the combined score
            qsar_hsec_threshold: Threshold for normalizing QSAR-HSEC predictions
            qsar_sip_threshold: Threshold for normalizing QSAR-SIP predictions
            logp_threshold: Threshold for normalizing LogP values
            similarity_threshold: Threshold for normalizing similarity values
            qed_threshold: Threshold for normalizing QED values (typically 1.0)
            secret_name: AWS secret name containing API credentials
            tenant_id: Microsoft tenant ID
            qsar_hsec_endpoint: API endpoint for QSAR-HSEC predictions
            qsar_sip_endpoint: API endpoint for QSAR-SIP predictions
        """
        # Store configuration
        self.reference_smiles = reference_smiles
        self.use_qsar_hsec = use_qsar_hsec
        self.use_qsar_sip = use_qsar_sip
        self.use_logp = use_logp
        self.use_similarity = use_similarity
        self.use_qed = use_qed
        
        # Store weights
        self.qsar_hsec_weight = qsar_hsec_weight if use_qsar_hsec else 0.0
        self.qsar_sip_weight = qsar_sip_weight if use_qsar_sip else 0.0
        self.logp_weight = logp_weight if use_logp else 0.0
        self.similarity_weight = similarity_weight if use_similarity else 0.0
        self.qed_weight = qed_weight if use_qed else 0.0
        
        # Store thresholds
        self.qsar_hsec_threshold = qsar_hsec_threshold
        self.qsar_sip_threshold = qsar_sip_threshold
        self.logp_threshold = logp_threshold
        self.similarity_threshold = similarity_threshold
        self.qed_threshold = qed_threshold
        
        # API credentials
        self.secret_name = secret_name
        self.tenant_id = tenant_id
        self.qsar_hsec_endpoint = qsar_hsec_endpoint
        self.qsar_sip_endpoint = qsar_sip_endpoint
        
        # Reference fingerprint cache
        self.reference_fingerprint_cache = {}
        
        # Initialize reference molecule for similarity calculations if provided
        if reference_smiles is not None:
            self._initialize_reference_fingerprint(reference_smiles)
        
        # Initialize API connection if using any QSAR endpoint
        if self.use_qsar_hsec or self.use_qsar_sip:
            self._initialize_api_connection()
            
    def _initialize_reference_fingerprint(self, reference_smiles: str):
        """Initialize reference molecule fingerprint for similarity calculations."""
        if not self.use_similarity:
            return None, None
            
        # Check if we already have this fingerprint cached
        if reference_smiles in self.reference_fingerprint_cache:
            return self.reference_fingerprint_cache[reference_smiles]
            
        try:
            # Morgan fingerprint parameters
            self.fingerprint_radius = 2
            self.fingerprint_nbits = 2048
            
            # Generate reference fingerprint
            reference_mol = Chem.MolFromSmiles(reference_smiles)
            if reference_mol is None:
                print(f"Warning: Could not parse reference SMILES: {reference_smiles}")
                self.reference_fingerprint_cache[reference_smiles] = (None, None)
                return None, None
                
            reference_fingerprint = GetMorganFingerprintAsBitVect(
                reference_mol, 
                radius=self.fingerprint_radius, 
                nBits=self.fingerprint_nbits
            )
            
            # Cache the result
            self.reference_fingerprint_cache[reference_smiles] = (reference_mol, reference_fingerprint)
            return reference_mol, reference_fingerprint
            
        except Exception as e:
            print(f"Error initializing reference fingerprint: {e}")
            self.reference_fingerprint_cache[reference_smiles] = (None, None)
            return None, None
            
    def _initialize_api_connection(self):
        """Initialize the API connection and authentication for QSAR predictions."""
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
    
    def calculate_tanimoto_similarity(self, smiles: str, reference_smiles: str) -> float:
        """
        Calculate Tanimoto similarity between a SMILES string and a reference molecule.
        
        Args:
            smiles: SMILES string to compare with reference
            reference_smiles: Reference SMILES for comparison
            
        Returns:
            Tanimoto similarity score (0-1)
        """
        if not self.use_similarity:
            return 0.0
            
        # Get reference fingerprint (from cache or generate new)
        _, reference_fingerprint = self._get_or_create_reference_fingerprint(reference_smiles)
        if reference_fingerprint is None:
            return 0.0
            
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
        return TanimotoSimilarity(fingerprint, reference_fingerprint)
        
    def _get_or_create_reference_fingerprint(self, reference_smiles: str):
        """Get reference fingerprint from cache or create new one."""
        if reference_smiles in self.reference_fingerprint_cache:
            return self.reference_fingerprint_cache[reference_smiles]
        else:
            return self._initialize_reference_fingerprint(reference_smiles)
    
    def calculate_logp(self, smiles: str) -> float:
        """
        Calculate LogP for a molecule.
        
        Args:
            smiles: SMILES string of the molecule
            
        Returns:
            LogP value or 0.0 if molecule is invalid
        """
        if not self.use_logp:
            return 0.0
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 0.0
            
        try:
            return MolLogP(mol)
        except:
            return 0.0
    
    def calculate_qed(self, smiles: str) -> float:
        """
        Calculate Quantitative Estimate of Drug-likeness (QED) for a molecule.
        
        Args:
            smiles: SMILES string of the molecule
            
        Returns:
            QED score (0-1, higher is more drug-like)
        """
        if not self.use_qed:
            return 0.0
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 0.0
            
        try:
            return QED.qed(mol)
        except:
            return 0.0

    
    def get_qsar_predictions(self, smiles_list: List[str], endpoint: str, batch_size: int = 50) -> List[float]:
        """
        Get QSAR predictions for a list of SMILES strings using a specific endpoint.
        """
        #print(f"Getting QSAR predictions for {len(smiles_list)} molecules using endpoint: {endpoint}")
        
        # Process in batches to avoid timeouts
        all_predictions = []
        num_batches = (len(smiles_list) + batch_size - 1) // batch_size  # Ceiling division
        
        for i in range(0, len(smiles_list), batch_size):
            batch_smiles = smiles_list[i:i+batch_size]
            #print(f"Processing batch {(i//batch_size)+1}/{num_batches} with {len(batch_smiles)} molecules")
            
            batch_predictions = self._get_batch_qsar_predictions(batch_smiles, endpoint)
            all_predictions.extend(batch_predictions)
            
            # Add a small delay between batches
            if i + batch_size < len(smiles_list):
                time.sleep(0.5)
        
        #print(f"Completed QSAR predictions, got {len(all_predictions)} results")
        return all_predictions

    
    def get_qsar_hsec_predictions(self, smiles_list: List[str], batch_size: int = 50) -> List[float]:
        """
        Get QSAR-HSEC predictions for a list of SMILES strings.
        
        Args:
            smiles_list: List of SMILES strings
            batch_size: Batch size for API calls
            
        Returns:
            List of QSAR-HSEC prediction scores
        """
        if not self.use_qsar_hsec:
            return [0.0] * len(smiles_list)
        return self.get_qsar_predictions(smiles_list, self.qsar_hsec_endpoint, batch_size)
    
    def get_qsar_sip_predictions(self, smiles_list: List[str], batch_size: int = 50) -> List[float]:
        """
        Get QSAR-SIP predictions for a list of SMILES strings.
        
        Args:
            smiles_list: List of SMILES strings
            batch_size: Batch size for API calls
            
        Returns:
            List of QSAR-SIP prediction scores
        """
        if not self.use_qsar_sip:
            return [0.0] * len(smiles_list)
        return self.get_qsar_predictions(smiles_list, self.qsar_sip_endpoint, batch_size)
    
    def _get_batch_qsar_predictions(self, smiles_list: List[str], endpoint: str) -> List[float]:
        """Get QSAR predictions for a batch of SMILES strings using specific endpoint."""
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
        response = self._post_data(access_token, endpoint, request_input)
        
        # Print raw response for debugging
        #print(f"Raw API response: {json.dumps(response, indent=2)}")
        
        # Process response
        predictions = []
        
        if not isinstance(response, dict):
            print(f"Unexpected response type: {type(response)}")
            return [0.0] * len(smiles_list)
        
        # Handle the new response format which has 'data' instead of 'result'
        data_items = response.get("data", [])
        
        if not data_items:
            print(f"No data items found in response")
            return [0.0] * len(smiles_list)
        
        # Create a mapping from id to prediction
        id_to_prediction = {}
        
        for item in data_items:
            # Extract the ID from input
            item_id = item.get("input", {}).get("id")
            
            # Check if predictions exist and status is successful
            if item.get("status", {}).get("success", False):
                predictions_dict = item.get("predictions", {})
                
                # Look for pred_hsec or other prediction fields
                if "pred_hsec" in predictions_dict and "hsec" in endpoint:
                    value = predictions_dict.get("pred_hsec")
                    id_to_prediction[item_id] = float(value) if value is not None else 0.0
                elif "pred_sip" in predictions_dict and "sip" in endpoint:
                    value = predictions_dict.get("pred_sip")
                    id_to_prediction[item_id] = float(value) if value is not None else 0.0
                else:
                    # Try to find any prediction field if specific ones not found
                    for key, value in predictions_dict.items():
                        if key.startswith("pred_"):
                            id_to_prediction[item_id] = float(value) if value is not None else 0.0
                            break
            else:
                print(f"Prediction failed for ID {item_id}: {item.get('status', {}).get('error', 'Unknown error')}")
        
        # Build the final predictions list in the same order as input SMILES
        for i in range(len(smiles_list)):
            item_id = str(i+1)
            if item_id in id_to_prediction:
                predictions.append(id_to_prediction[item_id])
            else:
                predictions.append(0.0)
                print(f"No prediction found for ID {item_id}")
        
        return predictions

    
    def score(self, smiles_list: List[str], reference_smiles: Optional[str] = None) -> ScoringResults:
        """
        Score a list of molecules using configured metrics.
        
        Args:
            smiles_list: List of SMILES strings to score
            reference_smiles: Reference SMILES for similarity calculation (overrides the one provided at initialization)
            
        Returns:
            ScoringResults containing combined scores and individual property scores
        """
        
        # Determine which reference SMILES to use
        ref_smiles = reference_smiles if reference_smiles is not None else self.reference_smiles
        
        if self.use_similarity and ref_smiles is None:
            raise ValueError("Reference SMILES must be provided either during initialization or when calling score()")
        
        # status
        print("#~now scoring~#")
        
        # Calculate individual property scores
        similarity_scores = [self.calculate_tanimoto_similarity(smiles, ref_smiles) for smiles in smiles_list] if self.use_similarity else [0.0] * len(smiles_list)
        logp_scores = [self.calculate_logp(smiles) for smiles in smiles_list] if self.use_logp else [0.0] * len(smiles_list)
        qed_scores = [self.calculate_qed(smiles) for smiles in smiles_list] if self.use_qed else [0.0] * len(smiles_list)
        qsar_hsec_scores = self.get_qsar_hsec_predictions(smiles_list) if self.use_qsar_hsec else [0.0] * len(smiles_list)
        qsar_sip_scores = self.get_qsar_sip_predictions(smiles_list) if self.use_qsar_sip else [0.0] * len(smiles_list)
        
        # Convert to numpy arrays for calculations
        similarity_array = np.array(similarity_scores)
        logp_array = np.array(logp_scores)
        qed_array = np.array(qed_scores)
        qsar_hsec_array = np.array(qsar_hsec_scores)
        qsar_sip_array = np.array(qsar_sip_scores)
        
        # Normalize scores based on thresholds
        norm_similarity = np.clip(similarity_array / self.similarity_threshold, a_min=0.0, a_max=1.0)
        norm_logp = np.clip(logp_array / self.logp_threshold, a_min=0.0, a_max=1.0)
        norm_qed = np.clip(qed_array / self.qed_threshold, a_min=0.0, a_max=1.0)
        norm_qsar_hsec = np.clip(qsar_hsec_array / self.qsar_hsec_threshold, a_min=0.0, a_max=1.0)
        norm_qsar_sip = np.clip(qsar_sip_array / self.qsar_sip_threshold, a_min=0.0, a_max=1.0)
        
        # Combine scores with weights (negative for minimization)
        combined_scores = -1.0 * (
            self.similarity_weight * norm_similarity +
            self.logp_weight * norm_logp +
            self.qed_weight * norm_qed +
            self.qsar_hsec_weight * norm_qsar_hsec +
            self.qsar_sip_weight * norm_qsar_sip
        )
        
        # status
        print("#~Done scoring~#")
        
        # Return structured results
        return ScoringResults(
            combined_scores=combined_scores,
            similarity_scores=similarity_scores,
            logp_scores=logp_scores,
            qed_scores=qed_scores,
            qsar_hsec_scores=qsar_hsec_scores,
            qsar_sip_scores=qsar_sip_scores
        )

    
    def score_from_properties(self, 
                    logp_values=None, 
                    similarity_values=None, 
                    qsar_hsec_values=None,
                    qsar_sip_values=None,
                    qed_values=None) -> np.ndarray:
        """
        Score directly from property values instead of calculating from SMILES.
        This is useful for visualization and analysis of the scoring landscape.
        
        Args:
            logp_values: LogP values (numpy array or scalar)
            similarity_values: Tanimoto similarity values (numpy array or scalar)
            qsar_hsec_values: QSAR-HSEC prediction values (numpy array or scalar)
            qsar_sip_values: QSAR-SIP prediction values (numpy array or scalar)
            qed_values: QED drug-likeness values (numpy array or scalar)
            
        Returns:
            Array of combined scores (negative values, lower is better for minimization)
        """
        # Initialize arrays to zeros with the shape of the first non-None input
        shape = None
        for values in [logp_values, similarity_values, qsar_hsec_values, qsar_sip_values, qed_values]:
            if values is not None:
                if hasattr(values, 'shape'):
                    shape = values.shape
                    break
                elif not np.isscalar(values):
                    shape = np.array(values).shape
                    break
        
        if shape is None:
            # No valid inputs, return a single zero
            return np.array([0.0])
        
        # Initialize arrays to zeros with the determined shape
        logp_array = np.zeros(shape)
        similarity_array = np.zeros(shape)
        qsar_hsec_array = np.zeros(shape)
        qsar_sip_array = np.zeros(shape)
        qed_array = np.zeros(shape)
        
        # Update arrays with provided values (if any)
        if logp_values is not None:
            logp_array = np.array(logp_values)
        if similarity_values is not None:
            similarity_array = np.array(similarity_values)
        if qsar_hsec_values is not None:
            qsar_hsec_array = np.array(qsar_hsec_values)
        if qsar_sip_values is not None:
            qsar_sip_array = np.array(qsar_sip_values)
        if qed_values is not None:
            qed_array = np.array(qed_values)
        
        # Normalize scores based on thresholds
        norm_similarity = np.clip(similarity_array / self.similarity_threshold, a_min=0.0, a_max=1.0)
        norm_logp = np.clip(logp_array / self.logp_threshold, a_min=0.0, a_max=1.0)
        norm_qed = np.clip(qed_array / self.qed_threshold, a_min=0.0, a_max=1.0)
        norm_qsar_hsec = np.clip(qsar_hsec_array / self.qsar_hsec_threshold, a_min=0.0, a_max=1.0)
        norm_qsar_sip = np.clip(qsar_sip_array / self.qsar_sip_threshold, a_min=0.0, a_max=1.0)
        
        # Combine scores with weights (negative for minimization)
        combined_scores = -1.0 * (
            self.similarity_weight * norm_similarity +
            self.logp_weight * norm_logp +
            self.qed_weight * norm_qed +
            self.qsar_hsec_weight * norm_qsar_hsec +
            self.qsar_sip_weight * norm_qsar_sip
        )
        
        return combined_scores


def create_scorer(
    reference_smiles: Optional[str] = None,
    use_qsar_hsec: bool = False,
    use_qsar_sip: bool = True,
    use_logp: bool = False,
    use_similarity: bool = True,
    use_qed: bool = False,
    qsar_hsec_weight: float = 1.0,
    qsar_sip_weight: float = 1.0,
    logp_weight: float = 1.0,
    similarity_weight: float = 1.0,
    qed_weight: float = 1.0,
    qsar_hsec_threshold: float = 0.9,
    qsar_sip_threshold: float = 0.9,
    logp_threshold: float = 0.9,
    similarity_threshold: float = 0.4,
    qed_threshold: float = 0.9,
) -> MoleculeScorer:
    """
    Create a configured MoleculeScorer.
    
    Args:
        reference_smiles: Reference SMILES for similarity calculation (optional, can be provided later)
        use_qsar_hsec: Whether to include QSAR-HSEC predictions
        use_qsar_sip: Whether to include QSAR-SIP predictions
        use_logp: Whether to include LogP
        use_similarity: Whether to include similarity
        use_qed: Whether to include QED
        qsar_hsec_weight: Weight for QSAR-HSEC component
        qsar_sip_weight: Weight for QSAR-SIP component
        logp_weight: Weight for LogP component
        similarity_weight: Weight for similarity component
        qed_weight: Weight for QED component
        
    Returns:
        Configured MoleculeScorer
    """
    return MoleculeScorer(
        reference_smiles=reference_smiles,
        use_qsar_hsec=use_qsar_hsec,
        use_qsar_sip=use_qsar_sip,
        use_logp=use_logp,
        use_similarity=use_similarity,
        use_qed=use_qed,
        qsar_hsec_weight=qsar_hsec_weight,
        qsar_sip_weight=qsar_sip_weight,
        logp_weight=logp_weight,
        similarity_weight=similarity_weight,
        qed_weight=qed_weight,
        qsar_hsec_threshold = qsar_hsec_threshold,
        qsar_sip_threshold = qsar_sip_threshold,
        logp_threshold = logp_threshold,
        similarity_threshold = similarity_threshold,
        qed_threshold = qed_threshold,
    )


def simple_scoring_function(
    smiles_list: List[str],
    reference_smiles: str,
    use_qsar_hsec: bool = False,
    use_qsar_sip: bool = True,
    use_logp: bool = False,
    use_similarity: bool = True,
    use_qed: bool = False,
    qsar_hsec_weight: float = 1.0,
    qsar_sip_weight: float = 1.0,
    logp_weight: float = 1.0,
    similarity_weight: float = 1.0,
    qed_weight: float = 1.0
) -> np.ndarray:
    """
    Simplified function for scoring molecules with common configurations.
    Returns only the combined scores.
    
    Args:
        smiles_list: List of SMILES strings
        reference_smiles: Reference SMILES for similarity calculations
        use_qsar_hsec: Whether to include QSAR-HSEC predictions
        use_qsar_sip: Whether to include QSAR-SIP predictions
        use_logp: Whether to include LogP
        use_similarity: Whether to include similarity
        use_qed: Whether to include QED
        qsar_hsec_weight: Weight for QSAR-HSEC component
        qsar_sip_weight: Weight for QSAR-SIP component
        logp_weight: Weight for LogP component
        similarity_weight: Weight for similarity component
        qed_weight: Weight for QED component
        
    Returns:
        Array of combined scores (lower is better for minimization)
    """
    scorer = create_scorer(
        use_qsar_hsec=use_qsar_hsec,
        use_qsar_sip=use_qsar_sip,
        use_logp=use_logp,
        use_similarity=use_similarity,
        use_qed=use_qed,
        qsar_hsec_weight=qsar_hsec_weight,
        qsar_sip_weight=qsar_sip_weight,
        logp_weight=logp_weight,
        similarity_weight=similarity_weight,
        qed_weight=qed_weight
    )
    
    results = scorer.score(smiles_list, reference_smiles=reference_smiles)
    return results.combined_scores

# Simple test function
def test_qsar_prediction(smiles: str, endpoint: str):
    scorer = MoleculeScorer(use_qsar_hsec=True, use_qsar_sip=True)
    if "hsec" in endpoint.lower():
        scores = scorer.get_qsar_hsec_predictions([smiles])
    else:
        scores = scorer.get_qsar_sip_predictions([smiles])
    print(f"SMILES: {smiles}")
    print(f"Endpoint: {endpoint}")
    print(f"Score: {scores[0]}")
    return scores[0]


# Example usage in a Jupyter notebook:
"""
from qsar_scoring import create_scorer, enhanced_scoring_function

# Option 1: Use the enhanced_scoring_function for a quick result
scores = enhanced_scoring_function(
    molecules,                # List of SMILES strings
    reference_smiles=m1,      # Reference molecule
    use_qsar_sip=True,        # Include QSAR-SIP model
    use_qsar_hsec=False,      # Exclude QSAR-HSEC model
    use_similarity=True,      # Include similarity
    use_logp=False,           # Exclude LogP
    use_qed=False             # Exclude QED
)

# Option 2: Create a scorer object for more detailed results
scorer = create_scorer(
    use_qsar_sip=True,
    use_similarity=True
)

# Get detailed scoring results with both QSAR models
results = scorer.score(molecules, reference_smiles=m1)

# Option 3: Reuse scorer with different reference molecules
results_alt = scorer.score(molecules, reference_smiles=m2)

# Access specific components
combined_scores = results.combined_scores
qsar_sip_scores = results.qsar_sip_scores
qsar_hsec_scores = results.qsar_hsec_scores
similarity_scores = results.similarity_scores

# Print high-scoring molecules
for i, (smiles, score, qsar_sip, sim) in enumerate(zip(
    molecules, 
    results.combined_scores, 
    results.qsar_sip_scores, 
    results.similarity_scores
)):
    # Note: scores are negative, so higher (less negative) is better
    if score > -1.5:
        print(f"Molecule {i}: {smiles}")
        print(f"  Combined score: {score:.2f}")
        print(f"  QSAR-SIP score: {qsar_sip:.2f}")
        print(f"  Similarity: {sim:.2f}")

# Use multiple metrics simultaneously
multi_metric_scorer = create_scorer(
    use_qsar_sip=True,        # Include QSAR-SIP predictions  
    use_qsar_hsec=True,       # Include QSAR-HSEC predictions
    use_similarity=True,      # Include structural similarity
    use_qed=True,             # Include drug-likeness
    qsar_sip_weight=1.0,      # Weight for QSAR-SIP
    qsar_hsec_weight=0.8,     # Weight for QSAR-HSEC
    similarity_weight=0.6,    # Weight for similarity
    qed_weight=0.4            # Weight for drug-likeness
)

multi_metric_results = multi_metric_scorer.score(molecules, reference_smiles=m1)
"""


def visualize_normalized_landscape(scorer, x_scores, y_scores, grid_size=100):
    """
    Create a visualization with normalized QSAR values (0-1 range).
    
    Args:
        scorer: The MoleculeScorer instance
        x_scores: QSAR values (will be normalized to 0-1)
        y_scores: Similarity values
        grid_size: Number of points in each dimension
        
    Returns:
        Xgrid, Ygrid, score_grid, normalized_x_scores (for plotting)
    """
    # Get data ranges for normalization
    x_min, x_max = min(x_scores), max(x_scores)
    y_min, y_max = min(y_scores), max(y_scores)
    
    # Normalize x_scores to 0-1 range
    normalized_x_scores = [(x - x_min) / (x_max - x_min) for x in x_scores]
    
    # Create grid in normalized space (0-1 for x)
    x_vals = np.linspace(0, 1, grid_size)
    y_vals = np.linspace(y_min, y_max, grid_size)
    Xgrid, Ygrid = np.meshgrid(x_vals, y_vals)
    
    # Initialize score grid
    score_grid = np.zeros_like(Xgrid)
    
    # Fill the grid with smooth transition across the full range
    for i in range(grid_size):
        for j in range(grid_size):
            x_val = Xgrid[i, j]  # Already normalized to 0-1
            y_val = Ygrid[i, j]
            
            # Normalize similarity (standard approach)
            norm_similarity = min(1.0, y_val / scorer.similarity_threshold)
            
            # Calculate combined score - smooth gradient across entire range
            weighted_x = scorer.qsar_hsec_weight * x_val
            weighted_y = scorer.similarity_weight * norm_similarity
            
            # Combine weighted components (use positive values for visualization)
            score_grid[i, j] = weighted_x + weighted_y
    
    return Xgrid, Ygrid, score_grid, normalized_x_scores

# Usage example:
"""
# Create the normalized visualization
Xgrid, Ygrid, score_grid, normalized_x_scores = visualize_normalized_landscape(
    scorer, x_scores, y_scores
)

# Plot with improved color mapping
plt.figure(figsize=(12, 10))
contour = plt.contourf(Xgrid, Ygrid, score_grid, levels=100, cmap='viridis')
plt.colorbar(contour, label='Score')

# Convert the reference point to normalized scale
normalized_reference = (imatinib_logp - min(x_scores)) / (max(x_scores) - min(x_scores))
plt.plot([normalized_reference, normalized_reference], [min(y_scores), max(y_scores)], 
         linestyle=':', c='k', label="Original seeds value")

# Plot data points using normalized x-coordinates
plt.scatter(normalized_x_scores, y_scores, c=results.combined_scores, 
            edgecolor='white', linewidths=0.5, cmap='viridis')

plt.xlabel("Normalized QSAR Score (0-1)")
plt.ylabel("Similarity")
plt.title("Molecular Scoring Landscape - Normalized View")
plt.legend()
plt.show()
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

def plot_molecular_landscape(results, scorer, seed_score, metric='qsar_hsec_scores'):
    """
    Create a molecular scoring landscape plot that works with different metrics and scales.
    
    Parameters:
    -----------
    results : DataFrame
        Contains the scoring data
    scorer : object
        Object with score_from_properties method for calculating background
    seed_score : float
        Reference value for original seeds
    metric : str
        Column name in results to use for x-axis (default: 'qsar_hsec_scores')
    """
    # Extract scores based on the chosen metric and convert to numpy arrays
    x_scores_raw = np.array(getattr(results, metric))
    y_scores = np.array(results.similarity_scores)
    
    # Get combined scores for point coloring
    color_values = -1 * np.array(results.combined_scores)
    
    # Normalize x_scores to 0-1 range for contour background
    x_min, x_max = np.min(x_scores_raw), np.max(x_scores_raw)
    x_scores_norm = (x_scores_raw - x_min) / (x_max - x_min) if x_max > x_min else x_scores_raw
    
    # Normalize seed score to the same 0-1 range, but handle out-of-range values
    # Clamp the seed score to the plot boundaries if it's outside the range
    clamped_seed_score = max(min(seed_score, x_max), x_min)
    seed_score_norm = (clamped_seed_score - x_min) / (x_max - x_min) if x_max > x_min else clamped_seed_score
    
    # Track if the seed score is out of range for annotation
    seed_score_out_of_range = seed_score != clamped_seed_score
    
    # Create a uniform grid from 0 to 1 for the normalized space
    deltas = np.linspace(0, 1, 100)
    Xgrid, Ygrid = np.meshgrid(deltas, deltas)
    
    # Calculate scores across the normalized grid
    score_background = -1 * scorer.score_from_properties(similarity_values=Ygrid, qed_values=Xgrid)
    
    # Create a colormap normalizer that works for both the contour and scatter plot
    norm = Normalize(vmin=min(color_values), vmax=max(color_values))
    
    # Create plot
    plt.figure(figsize=(10, 8))
    
    # Plot contour background (0-1 normalized range)
    contour = plt.contourf(Xgrid, Ygrid, score_background, levels=50, cmap='viridis')
    plt.colorbar(label='Combined Score')
    
    # Add two x-axes - one normalized (0-1) and one with the original scale
    ax1 = plt.gca()
    ax2 = ax1.twiny()  # Create a twin x-axis
    
    # Plot reference line for seed value in normalized space
    line_label = "Original seeds value"
    if seed_score_out_of_range:
        if seed_score > x_max:
            line_label = f"Original seeds value (off chart: {seed_score:.1f} > {x_max:.1f})"
            ax1.plot([1, 1], [0, 1], linestyle=':', c='r', label=line_label)
            # Add an arrow pointing to the right
            ax1.annotate('', xy=(1.05, 0.5), xytext=(0.95, 0.5),
                        arrowprops=dict(facecolor='red', shrink=0.05, width=1.5, headwidth=8))
        else:  # seed_score < x_min
            line_label = f"Original seeds value (off chart: {seed_score:.1f} < {x_min:.1f})"
            ax1.plot([0, 0], [0, 1], linestyle=':', c='r', label=line_label)
            # Add an arrow pointing to the left
            ax1.annotate('', xy=(-0.05, 0.5), xytext=(0.05, 0.5),
                        arrowprops=dict(facecolor='red', shrink=0.05, width=1.5, headwidth=8))
    else:
        ax1.plot([seed_score_norm, seed_score_norm], [0, 1], linestyle=':', c='k', label=line_label)
    
    # Plot the data points using normalized x-coordinates but with color from combined score
    scatter = ax1.scatter(x_scores_norm, y_scores, c=color_values, edgecolor='white', 
                         linewidths=0.5, cmap='viridis', norm=norm)
    
    # Set up the primary axis (normalized)
    ax1.set_xlabel("Normalized Score")
    ax1.set_ylabel("Similarity")
    ax1.set_title(f"Molecular Scoring Landscape ({metric})")
    ax1.legend()
    ax1.grid(alpha=0.3, linestyle='--')
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    
    # Set up the secondary axis (original scale)
    # Create tick positions that correspond to nice values in the original scale
    tick_positions = np.linspace(0, 1, 5)
    original_tick_values = x_min + tick_positions * (x_max - x_min)
    original_tick_labels = [f"{val:.1f}" for val in original_tick_values]
    
    ax2.set_xlim(0, 1)  # Same as the normalized axis
    ax2.set_xticks(tick_positions)
    ax2.set_xticklabels(original_tick_labels)
    ax2.set_xlabel(f"Original {metric} Scale")
    
    # Add annotations showing ranges and reference value
    plt.figtext(0.05, 0.07, f"{metric} range: [{np.min(x_scores_raw):.1f}, {np.max(x_scores_raw):.1f}]")
    
    # Modify seed value annotation if out of range
    seed_annotation = f"Reference seed value: {seed_score:.1f}"
    if seed_score_out_of_range:
        if seed_score > x_max:
            seed_annotation += f" (above max)"
        else:
            seed_annotation += f" (below min)"
    plt.figtext(0.05, 0.01, seed_annotation, color='red' if seed_score_out_of_range else 'black')
    
    plt.figtext(0.05, 0.04, f"Similarity range: [{np.min(y_scores):.2f}, {np.max(y_scores):.2f}]")
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.85)  # Make room for the second x-axis label
    
    return plt

# Example usage:
# plot = plot_molecular_landscape(results, scorer, seed_score, metric='qsar_hsec_scores')
# plot.show()

# Alternative usage with a different metric:
# plot = plot_molecular_landscape(results, scorer, seed_score, metric='alternative_metric')
# plot.show()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

def plot_molecular_landscape_old(results, scorer, seed_score, metric='qsar_hsec_scores'):
    """
    Create a molecular scoring landscape plot that works with different metrics and scales.
    
    Parameters:
    -----------
    results : DataFrame
        Contains the scoring data
    scorer : object
        Object with score_from_properties method for calculating background
    seed_score : float
        Reference value for original seeds
    metric : str
        Column name in results to use for x-axis (default: 'qsar_hsec_scores')
    """
    # Extract scores based on the chosen metric and convert to numpy arrays
    x_scores_raw = np.array(getattr(results, metric))
    y_scores = np.array(results.similarity_scores)
    
    # Get combined scores for point coloring
    color_values = -1 * np.array(results.combined_scores)
    
    # Normalize x_scores to 0-1 range for contour background
    x_min, x_max = np.min(x_scores_raw), np.max(x_scores_raw)
    x_scores_norm = (x_scores_raw - x_min) / (x_max - x_min) if x_max > x_min else x_scores_raw
    
    # Normalize seed score to the same 0-1 range
    seed_score_norm = (seed_score - x_min) / (x_max - x_min) if x_max > x_min else seed_score
    
    # Create a uniform grid from 0 to 1 for the normalized space
    deltas = np.linspace(0, 1, 100)
    Xgrid, Ygrid = np.meshgrid(deltas, deltas)
    
    # Calculate scores across the normalized grid
    score_background = -1 * scorer.score_from_properties(similarity_values=Ygrid, qed_values=Xgrid)
    
    # Create a colormap normalizer that works for both the contour and scatter plot
    norm = Normalize(vmin=min(color_values), vmax=max(color_values))
    
    # Create plot
    plt.figure(figsize=(10, 8))
    
    # Plot contour background (0-1 normalized range)
    contour = plt.contourf(Xgrid, Ygrid, score_background, levels=50, cmap='viridis')
    plt.colorbar(label='Combined Score')
    
    # Add two x-axes - one normalized (0-1) and one with the original scale
    ax1 = plt.gca()
    ax2 = ax1.twiny()  # Create a twin x-axis
    
    # Plot reference line for seed value in normalized space
    ax1.plot([seed_score_norm, seed_score_norm], [0, 1], linestyle=':', c='k', label="Original seeds value")
    
    # Plot the data points using normalized x-coordinates but with color from combined score
    scatter = ax1.scatter(x_scores_norm, y_scores, c=color_values, edgecolor='white', 
                         linewidths=0.5, cmap='viridis', norm=norm)
    
    # Set up the primary axis (normalized)
    ax1.set_xlabel("Normalized Score")
    ax1.set_ylabel("Similarity")
    ax1.set_title(f"Molecular Scoring Landscape ({metric})")
    ax1.legend()
    ax1.grid(alpha=0.3, linestyle='--')
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    
    # Set up the secondary axis (original scale)
    # Create tick positions that correspond to nice values in the original scale
    tick_positions = np.linspace(0, 1, 5)
    original_tick_values = x_min + tick_positions * (x_max - x_min)
    original_tick_labels = [f"{val:.1f}" for val in original_tick_values]
    
    ax2.set_xlim(0, 1)  # Same as the normalized axis
    ax2.set_xticks(tick_positions)
    ax2.set_xticklabels(original_tick_labels)
    ax2.set_xlabel(f"Original {metric} Scale")
    
    # Add annotations showing ranges and reference value
    plt.figtext(0.05, 0.07, f"{metric} range: [{min(x_scores_raw):.1f}, {max(x_scores_raw):.1f}]")
    plt.figtext(0.05, 0.01, f"Reference seed value: {seed_score:.1f}")
    plt.figtext(0.05, 0.04, f"Similarity range: [{min(y_scores):.2f}, {max(y_scores):.2f}]")
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.85)  # Make room for the second x-axis label
    
    return plt

# Example usage:
# plot = plot_molecular_landscape(results, scorer, seed_score, metric='qsar_hsec_scores')
# plot.show()

# Alternative usage with a different metric:
# plot = plot_molecular_landscape(results, scorer, seed_score, metric='alternative_metric')
# plot.show()


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

def plot_molecular_landscape_v1(results, scorer, seed_score, metric='qsar_hsec_scores'):
    """
    Create a molecular scoring landscape plot that works with different metrics and scales.
    
    Parameters:
    -----------
    results : DataFrame
        Contains the scoring data
    scorer : object
        Object with score_from_properties method for calculating background
    seed_score : float
        Reference value for original seeds
    metric : str
        Column name in results to use for x-axis (default: 'qsar_hsec_scores')
    """
    # Extract scores based on the chosen metric and convert to numpy arrays
    x_scores_raw = np.array(getattr(results, metric))
    y_scores = np.array(results.similarity_scores)
    
    # Get combined scores for point coloring
    color_values = -1 * np.array(results.combined_scores)
    
    # Determine min/max including the seed score to ensure it's in range
    # Use a margin to avoid values exactly at the edge of the plot
    margin = 0.05  # 5% margin
    
    # Calculate the data range
    data_min = np.min(x_scores_raw)
    data_max = np.max(x_scores_raw)
    
    # Expand range to include seed value plus margins
    plot_min = min(data_min, seed_score) - margin * abs(data_max - data_min)
    plot_max = max(data_max, seed_score) + margin * abs(data_max - data_min)
    
    # Normalize x_scores to 0-1 range for contour background, using the expanded range
    x_scores_norm = (x_scores_raw - plot_min) / (plot_max - plot_min) if plot_max > plot_min else x_scores_raw
    
    # Normalize seed score to the same 0-1 range
    seed_score_norm = (seed_score - plot_min) / (plot_max - plot_min) if plot_max > plot_min else seed_score
    
    # Create a uniform grid from 0 to 1 for the normalized space
    deltas = np.linspace(0, 1, 100)
    Xgrid, Ygrid = np.meshgrid(deltas, deltas)
    
    # Calculate scores across the normalized grid
    score_background = -1 * scorer.score_from_properties(similarity_values=Ygrid, qed_values=Xgrid)
    
    # Create a colormap normalizer that works for both the contour and scatter plot
    norm = Normalize(vmin=min(color_values), vmax=max(color_values))
    
    # Create plot
    plt.figure(figsize=(10, 8))
    
    # Plot contour background (0-1 normalized range)
    contour = plt.contourf(Xgrid, Ygrid, score_background, levels=50, cmap='viridis')
    plt.colorbar(label='Combined Score')
    
    # Add two x-axes - one normalized (0-1) and one with the original scale
    ax1 = plt.gca()
    ax2 = ax1.twiny()  # Create a twin x-axis
    
    # Plot reference line for seed value in normalized space
    ax1.plot([seed_score_norm, seed_score_norm], [0, 1], linestyle=':', c='k', label="Original seeds value")
    
    # Plot the data points using normalized x-coordinates but with color from combined score
    scatter = ax1.scatter(x_scores_norm, y_scores, c=color_values, edgecolor='white', 
                         linewidths=0.5, cmap='viridis', norm=norm)
    
    # Set up the primary axis (normalized)
    ax1.set_xlabel("Normalized Score")
    ax1.set_ylabel("Similarity")
    ax1.set_title(f"Molecular Scoring Landscape ({metric})")
    ax1.legend()
    ax1.grid(alpha=0.3, linestyle='--')
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    
    # Set up the secondary axis (original scale)
    # Create tick positions that correspond to nice values in the original scale
    tick_positions = np.linspace(0, 1, 5)
    original_tick_values = plot_min + tick_positions * (plot_max - plot_min)
    original_tick_labels = [f"{val:.1f}" for val in original_tick_values]
    
    ax2.set_xlim(0, 1)  # Same as the normalized axis
    ax2.set_xticks(tick_positions)
    ax2.set_xticklabels(original_tick_labels)
    ax2.set_xlabel(f"Original {metric} Scale")
    
    # Add annotations showing ranges and reference value
    plt.figtext(0.05, 0.07, f"{metric} range: [{np.min(x_scores_raw):.1f}, {np.max(x_scores_raw):.1f}]")
    plt.figtext(0.05, 0.01, f"Reference seed value: {seed_score:.1f}")
    plt.figtext(0.05, 0.04, f"Similarity range: [{np.min(y_scores):.2f}, {np.max(y_scores):.2f}]")
    
    #plt.tight_layout()
    plt.subplots_adjust(top=0.85)  # Make room for the second x-axis label
    
    return plt

# Example usage:
# plot = plot_molecular_landscape(results, scorer, seed_score, metric='qsar_hsec_scores')
# plot.show()

# Alternative usage with a different metric:
# plot = plot_molecular_landscape(results, scorer, seed_score, metric='alternative_metric')
# plot.show()


import numpy as np
import matplotlib.pyplot as plt

def plot_molecular_landscape_man(results, scorer, seed_score, metric='qsar_hsec_scores', x_min=0, x_max=1):
    """
    Create a molecular scoring landscape plot with manually specified x-axis limits.
    
    Parameters:
    -----------
    results : DataFrame
        Contains the scoring data
    scorer : object
        Object with score_from_properties method for calculating background
    seed_score : float
        Reference value for original seeds
    metric : str
        Column name in results to use for x-axis (default: 'qsar_hsec_scores')
    x_min : float
        Minimum value for x-axis (default: 0)
    x_max : float
        Maximum value for x-axis (default: 1)
    """
    # Extract scores based on the chosen metric and convert to numpy arrays if needed
    x_scores_raw = np.array(getattr(results, metric))
    y_scores = np.array(results.similarity_scores)
    
    # Get combined scores for point coloring
    color_values = -1 * np.array(results.combined_scores)
    
    # Create a uniform grid from 0 to 1 for the background
    deltas = np.linspace(0, 1, 100)
    Xgrid, Ygrid = np.meshgrid(deltas, deltas)
    
    # Calculate scores across the grid for background
    score_background = -1 * scorer.score_from_properties(similarity_values=Ygrid, qed_values=Xgrid)
    
    # Create plot
    plt.figure(figsize=(7, 6))
    
    # Plot contour background (0-1 range)
    contour = plt.contourf(Xgrid, Ygrid, score_background, levels=50, cmap='viridis', vmin=0, vmax=2.0)
    plt.colorbar(label='score', ticks=np.linspace(0, 2, 11))
    
    # Normalize seed to the 0-1 range for plotting
    seed_norm = (seed_score - x_min) / (x_max - x_min)
    plt.plot([seed_norm, seed_norm], [0, 1], linestyle=':', c='k', label="Original seeds value")
    
    # Normalize x_scores to 0-1 range for plotting on the background
    x_scores_norm = (x_scores_raw - x_min) / (x_max - x_min)
    
    # Plot data points, but filter out points outside our x range
    in_range = (x_scores_raw >= x_min) & (x_scores_raw <= x_max)
    plt.scatter(
        x_scores_norm[in_range], 
        y_scores[in_range], 
        c=color_values[in_range], 
        edgecolor='white', 
        linewidths=0.5, 
        cmap='viridis'
    )
    
    # Count how many points are out of range
    n_below = np.sum(x_scores_raw < x_min)
    n_above = np.sum(x_scores_raw > x_max)
    
    # Setup axis labels and title
    plt.xlabel(f"{metric.replace('_', ' ').title()}")
    plt.ylabel("Similarity")
    plt.title("Molecular Scoring Landscape")
    plt.legend()
    plt.grid(alpha=0.3, linestyle='--')
    
    # Set limits to properly display the 0-1 normalized plot
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    
    # Create custom x-tick labels to show the real scale
    tick_positions = np.linspace(0, 1, 5)
    tick_values = x_min + tick_positions * (x_max - x_min)
    plt.xticks(tick_positions, [f"{v:.1f}" for v in tick_values])
    
    # Add annotations showing ranges and reference value
    plt.annotate(f"{metric.replace('_', ' ').title()} range: [{np.min(x_scores_raw):.1f}, {np.max(x_scores_raw):.1f}]",
                xy=(0.05, 0.07), xycoords='figure fraction')
    
    # Add annotation for seed value
    plt.annotate(f"Reference seed value: {seed_score:.1f}", 
                xy=(0.05, 0.01), xycoords='figure fraction')
    
    plt.annotate(f"Similarity range: [{np.min(y_scores):.2f}, {np.max(y_scores):.2f}]",
                xy=(0.05, 0.04), xycoords='figure fraction')
    
    # Add note if points are out of range
    if n_below > 0 or n_above > 0:
        out_msg = []
        if n_below > 0:
            out_msg.append(f"{n_below} point{'s' if n_below > 1 else ''} below x_min")
        if n_above > 0:
            out_msg.append(f"{n_above} point{'s' if n_above > 1 else ''} above x_max")
        plt.annotate(f"Note: {' and '.join(out_msg)} not shown", 
                    xy=(0.05, 0.10), xycoords='figure fraction', color='blue', style='italic')
    
    plt.tight_layout()
    return plt

# Example usage:
# For the typical 0-1 range:
# plot_molecular_landscape(results, scorer, seed_score, metric='qsar_hsec_scores', x_min=0, x_max=1)

# For a different range (e.g., -40 to 1200):
# plot_molecular_landscape(results, scorer, seed_score, metric='different_metric', x_min=-40, x_max=1200)