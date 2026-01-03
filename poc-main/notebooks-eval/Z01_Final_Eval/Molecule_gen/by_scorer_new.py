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
    qsar_aih_scores: List[float]    # Individual QSAR-AIH prediction scores
    qed_scores: List[float]         # Individual QED scores
    
    # Normalized scores (0-1 scale)
    norm_logp_scores: List[float]        # Normalized LogP scores
    norm_similarity_scores: List[float]  # Normalized similarity scores  
    norm_qsar_hsec_scores: List[float]   # Normalized QSAR-HSEC scores
    norm_qsar_aih_scores: List[float]    # Normalized QSAR-AIH scores
    norm_qed_scores: List[float]         # Normalized QED scores


class MoleculeScorer:
    """
    A configurable scorer for molecular optimization that can combine different metrics:
    - QSAR predictions
    - Structural similarity
    - LogP
    - QED (drug-likeness)
    """
    
    def create_scorer(
        reference_smiles: Optional[str] = None,
        use_qsar_hsec: bool = False,
        use_qsar_aih: bool = True,
        use_logp: bool = False,
        use_similarity: bool = True,
        use_qed: bool = False,
        qsar_hsec_weight: float = 1.0,
        qsar_aih_weight: float = 1.0,
        logp_weight: float = 1.0,
        similarity_weight: float = 1.0,
        qed_weight: float = 1.0,
        qsar_hsec_threshold: float = 0.9,
        qsar_aih_threshold: float = 0.9,
        logp_threshold: float = 0.9,
        similarity_threshold: float = 0.9,
        qed_threshold: float = 0.9,
        # Keep domain parameters for HSEC and AIH
        hsec_min: float = -50,
        hsec_max: float = 500,
        aih_min: float = -50,
        aih_max: float = 500,
        # Remove use_domain_normalization parameter
        logp_min: float = -3,  # Still keep these as they might be used elsewhere
        logp_max: float = 7
    ) -> 'MoleculeScorer':
        """
        Create a configured MoleculeScorer.
        
        Args:
            reference_smiles: Reference SMILES for similarity calculation (optional, can be provided later)
            use_qsar_hsec: Whether to include QSAR-HSEC predictions
            use_qsar_aih: Whether to include QSAR-AIH predictions
            use_logp: Whether to include LogP
            use_similarity: Whether to include similarity
            use_qed: Whether to include QED
            qsar_hsec_weight: Weight for QSAR-HSEC component
            qsar_aih_weight: Weight for QSAR-AIH component
            logp_weight: Weight for LogP component
            similarity_weight: Weight for similarity component
            qed_weight: Weight for QED component
            qsar_hsec_threshold: Threshold for historical compatibility (not used for HSEC)
            qsar_aih_threshold: Threshold for historical compatibility (not used for AIH)
            logp_threshold: Threshold for normalizing LogP
            similarity_threshold: Threshold for normalizing similarity
            qed_threshold: Threshold for normalizing QED
            hsec_min: Minimum expected HSEC value for domain normalization
            hsec_max: Maximum expected HSEC value for domain normalization
            aih_min: Minimum expected AIH value for domain normalization
            aih_max: Maximum expected AIH value for domain normalization
            logp_min: Minimum expected LogP (for historical compatibility)
            logp_max: Maximum expected LogP (for historical compatibility)
            
        Returns:
            Configured MoleculeScorer
        """
        return MoleculeScorer(
            reference_smiles=reference_smiles,
            use_qsar_hsec=use_qsar_hsec,
            use_qsar_aih=use_qsar_aih,
            use_logp=use_logp,
            use_similarity=use_similarity,
            use_qed=use_qed,
            qsar_hsec_weight=qsar_hsec_weight,
            qsar_aih_weight=qsar_aih_weight,
            logp_weight=logp_weight,
            similarity_weight=similarity_weight,
            qed_weight=qed_weight,
            qsar_hsec_threshold=qsar_hsec_threshold,
            qsar_aih_threshold=qsar_aih_threshold,
            logp_threshold=logp_threshold,
            similarity_threshold=similarity_threshold,
            qed_threshold=qed_threshold,
            hsec_min=hsec_min,
            hsec_max=hsec_max,
            aih_min=aih_min,
            aih_max=aih_max,
            logp_min=logp_min,
            logp_max=logp_max
            )
    
    def simple_scoring_function(
        smiles_list: List[str],
        reference_smiles: str,
        use_qsar_hsec: bool = False,
        use_qsar_aih: bool = True,
        use_logp: bool = False,
        use_similarity: bool = True,
        use_qed: bool = False,
        qsar_hsec_weight: float = 1.0,
        qsar_aih_weight: float = 1.0,
        logp_weight: float = 1.0,
        similarity_weight: float = 0.4,
        qed_weight: float = 1.0,
        # Remove use_domain_normalization parameter
        hsec_min: float = -50,
        hsec_max: float = 500,
        aih_min: float = -50,
        aih_max: float = 500,
        logp_min: float = -3,
        logp_max: float = 7,
        logp_threshold: float = 0.9,
        similarity_threshold: float = 0.4
    ) -> np.ndarray:
        """
        Simplified function for scoring molecules with common configurations.
        Returns only the combined scores.
        
        Args:
            smiles_list: List of SMILES strings
            reference_smiles: Reference SMILES for similarity calculations
            use_qsar_hsec: Whether to include QSAR-HSEC predictions
            use_qsar_aih: Whether to include QSAR-AIH predictions
            use_logp: Whether to include LogP
            use_similarity: Whether to include similarity
            use_qed: Whether to include QED
            qsar_hsec_weight: Weight for QSAR-HSEC component
            qsar_aih_weight: Weight for QSAR-AIH component
            logp_weight: Weight for LogP component
            similarity_weight: Weight for similarity component
            qed_weight: Weight for QED component
            hsec_min: Minimum expected HSEC value for domain normalization
            hsec_max: Maximum expected HSEC value for domain normalization
            aih_min: Minimum expected AIH value for domain normalization
            aih_max: Maximum expected AIH value for domain normalization
            logp_min: Minimum expected LogP (for historical compatibility)
            logp_max: Maximum expected LogP (for historical compatibility)
            logp_threshold: Threshold for normalizing LogP
            similarity_threshold: Threshold for normalizing similarity
            
        Returns:
            Array of combined scores (lower is better for minimization)
        """
        scorer = MoleculeScorer.create_scorer(
            use_qsar_hsec=use_qsar_hsec,
            use_qsar_aih=use_qsar_aih,
            use_logp=use_logp,
            use_similarity=use_similarity,
            use_qed=use_qed,
            qsar_hsec_weight=qsar_hsec_weight,
            qsar_aih_weight=qsar_aih_weight,
            logp_weight=logp_weight,
            similarity_weight=similarity_weight,
            qed_weight=qed_weight,
            hsec_min=hsec_min,
            hsec_max=hsec_max,
            aih_min=aih_min,
            aih_max=aih_max,
            logp_min=logp_min,
            logp_max=logp_max,
            logp_threshold=logp_threshold,
            similarity_threshold=similarity_threshold
        )
        
        results = scorer.score(smiles_list, reference_smiles=reference_smiles)
        return results.combined_scores
    
    TOKEN_FILE = os.path.join(os.path.expanduser('~'), ".insilicoclitoken")
    
    def __init__(self, 
                 reference_smiles: Optional[str] = None,
                 use_qsar_hsec: bool = True,
                 use_qsar_aih: bool = False, 
                 use_logp: bool = False,
                 use_similarity: bool = True,
                 use_qed: bool = False,
                 qsar_hsec_weight: float = 1.0,
                 qsar_aih_weight: float = 1.0,
                 logp_weight: float = 1.0,
                 similarity_weight: float = 0.4,
                 qed_weight: float = 1.0,
                 qsar_hsec_threshold: float = 0.9,
                 qsar_aih_threshold: float = 0.9,
                 logp_threshold: float = 0.9,
                 similarity_threshold: float = 0.8,
                 qed_threshold: float = 0.9,
                 # Add domain knowledge min/max parameters with defaults
                 hsec_min: float = -50,
                 hsec_max: float = 500,
                 aih_min: float = -50,
                 aih_max: float = 500,
                 logp_min: float = -3,
                 logp_max: float = 7,
                 # Remove use_domain_normalization parameter
                 secret_name: str = "PROPCALC_API_CREDS",
                 tenant_id: str = "fcb2b37b-5da0-466b-9b83-0014b67a7c78",
                 qsar_hsec_endpoint: str = "https://chemprop.cropkey-np.ag/activity/ica/invocations_hsec",
                 qsar_aih_endpoint: str = "https://chemprop.cropkey-np.ag/activity/ica/invocations"):
        """
        Initialize a configurable molecular scorer.
        
        Args:
            reference_smiles: Reference SMILES for similarity calculation
            use_qsar_hsec: Whether to include QSAR-HSEC predictions in scoring
            use_qsar_aih: Whether to include QSAR-AIH predictions in scoring
            use_logp: Whether to include LogP in scoring
            use_similarity: Whether to include structural similarity in scoring
            use_qed: Whether to include drug-likeness score in scoring
            qsar_hsec_weight: Weight for QSAR-HSEC predictions in the combined score
            qsar_aih_weight: Weight for QSAR-AIH predictions in the combined score
            logp_weight: Weight for LogP in the combined score
            similarity_weight: Weight for similarity in the combined score
            qed_weight: Weight for QED in the combined score
            qsar_hsec_threshold: Threshold for historical compatibility (not used for HSEC)
            qsar_aih_threshold: Threshold for historical compatibility (not used for AIH)
            logp_threshold: Threshold for normalizing LogP values
            similarity_threshold: Threshold for normalizing similarity values
            qed_threshold: Threshold for normalizing QED values
            hsec_min: Minimum expected value for HSEC scores for domain normalization
            hsec_max: Maximum expected value for HSEC scores for domain normalization
            aih_min: Minimum expected value for AIH scores for domain normalization
            aih_max: Maximum expected value for AIH scores for domain normalization
            logp_min: Minimum expected value for LogP (for historical compatibility)
            logp_max: Maximum expected value for LogP (for historical compatibility)
            secret_name: AWS secret name containing API credentials
            tenant_id: Microsoft tenant ID
            qsar_hsec_endpoint: API endpoint for QSAR-HSEC predictions
            qsar_aih_endpoint: API endpoint for QSAR-AIH predictions
        """
        # Store configuration
        self.reference_smiles = reference_smiles
        self.use_qsar_hsec = use_qsar_hsec
        self.use_qsar_aih = use_qsar_aih
        self.use_logp = use_logp
        self.use_similarity = use_similarity
        self.use_qed = use_qed
        
        # Store weights
        self.qsar_hsec_weight = qsar_hsec_weight if use_qsar_hsec else 0.0
        self.qsar_aih_weight = qsar_aih_weight if use_qsar_aih else 0.0
        self.logp_weight = logp_weight if use_logp else 0.0
        self.similarity_weight = similarity_weight if use_similarity else 0.0
        self.qed_weight = qed_weight if use_qed else 0.0
        
        # Store thresholds
        self.qsar_hsec_threshold = qsar_hsec_threshold
        self.qsar_aih_threshold = qsar_aih_threshold
        self.logp_threshold = logp_threshold
        self.similarity_threshold = similarity_threshold
        self.qed_threshold = qed_threshold
        
        # Store domain knowledge normalization parameters
        self.hsec_min = hsec_min
        self.hsec_max = hsec_max
        self.aih_min = aih_min
        self.aih_max = aih_max
        self.logp_min = logp_min
        self.logp_max = logp_max
        
        # API credentials
        self.secret_name = secret_name
        self.tenant_id = tenant_id
        self.qsar_hsec_endpoint = qsar_hsec_endpoint
        self.qsar_aih_endpoint = qsar_aih_endpoint
        
        # Reference fingerprint cache
        self.reference_fingerprint_cache = {}
        
        # Initialize reference molecule for similarity calculations if provided
        if reference_smiles is not None:
            self._initialize_reference_fingerprint(reference_smiles)
        
        # Initialize API connection if using any QSAR endpoint
        if self.use_qsar_hsec or self.use_qsar_aih:
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
    
    def get_qsar_aih_predictions(self, smiles_list: List[str], batch_size: int = 50) -> List[float]:
        """
        Get QSAR-AIH predictions for a list of SMILES strings.
        
        Args:
            smiles_list: List of SMILES strings
            batch_size: Batch size for API calls
            
        Returns:
            List of QSAR-AIH prediction scores
        """
        if not self.use_qsar_aih:
            return [0.0] * len(smiles_list)
        return self.get_qsar_predictions(smiles_list, self.qsar_aih_endpoint, batch_size)
    
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
                elif "pred_aih" in predictions_dict and "aih" in endpoint:
                    value = predictions_dict.get("pred_aih")
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

    
    def score_from_properties(self, 
                logp_values=None, 
                similarity_values=None, 
                qsar_hsec_values=None,
                qsar_aih_values=None,
                qed_values=None) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
        """
        Score directly from property values instead of calculating from SMILES.
        This is useful for visualization and analysis of the scoring landscape.
        
        Args:
            logp_values: LogP values (numpy array or scalar)
            similarity_values: Tanimoto similarity values (numpy array or scalar)
            qsar_hsec_values: QSAR-HSEC prediction values (numpy array or scalar)
            qsar_aih_values: QSAR-AIH prediction values (numpy array or scalar)
            qed_values: QED drug-likeness values (numpy array or scalar)
            
        Returns:
            Tuple containing:
            - Array of combined scores (negative values, lower is better for minimization)
            - Dictionary of normalized property arrays
        """
        # Initialize arrays to zeros with the shape of the first non-None input
        shape = None
        for values in [logp_values, similarity_values, qsar_hsec_values, qsar_aih_values, qed_values]:
            if values is not None:
                if hasattr(values, 'shape'):
                    shape = values.shape
                    break
                elif not np.isscalar(values):
                    shape = np.array(values).shape
                    break
        
        if shape is None:
            # No valid inputs, return a single zero
            return np.array([0.0]), {
                'norm_logp': np.array([0.0]),
                'norm_similarity': np.array([0.0]),
                'norm_qsar_hsec': np.array([0.0]),
                'norm_qsar_aih': np.array([0.0]),
                'norm_qed': np.array([0.0])
            }
        
        # Initialize arrays to zeros with the determined shape
        logp_array = np.zeros(shape)
        similarity_array = np.zeros(shape)
        qsar_hsec_array = np.zeros(shape)
        qsar_aih_array = np.zeros(shape)
        qed_array = np.zeros(shape)
        
        # Update arrays with provided values (if any)
        if logp_values is not None:
            logp_array = np.array(logp_values)
        if similarity_values is not None:
            similarity_array = np.array(similarity_values)
        if qsar_hsec_values is not None:
            qsar_hsec_array = np.array(qsar_hsec_values)
        if qsar_aih_values is not None:
            qsar_aih_array = np.array(qsar_aih_values)
        if qed_values is not None:
            qed_array = np.array(qed_values)
        
        # Normalize scores using different methods for different metrics:
        # - LogP, Similarity, QED: Use threshold-based normalization
        # - HSEC, AIH: Use domain-based min-max normalization
        
        # Threshold-based normalization for LogP, Similarity, QED
        norm_similarity = np.clip(similarity_array / self.similarity_threshold, a_min=0.0, a_max=1.0)
        norm_logp = np.clip(logp_array / self.logp_threshold, a_min=0.0, a_max=1.0)
        norm_qed = np.clip(qed_array / self.qed_threshold, a_min=0.0, a_max=1.0)
        
        # Domain-based min-max normalization for HSEC and AIH
        norm_qsar_hsec = np.clip((qsar_hsec_array - self.hsec_min) / (self.hsec_max - self.hsec_min), a_min=0.0, a_max=1.0)
        norm_qsar_aih = np.clip((qsar_aih_array - self.aih_min) / (self.aih_max - self.aih_min), a_min=0.0, a_max=1.0)
        
        # Combine scores with weights (negative for minimization)
        combined_scores = -1.0 * (
            self.similarity_weight * norm_similarity +
            self.logp_weight * norm_logp +
            self.qed_weight * norm_qed +
            self.qsar_hsec_weight * norm_qsar_hsec +
            self.qsar_aih_weight * norm_qsar_aih
        )
        
        # Return combined scores and normalized property arrays
        normalized_properties = {
            'norm_logp': norm_logp,
            'norm_similarity': norm_similarity,
            'norm_qed': norm_qed,
            'norm_qsar_hsec': norm_qsar_hsec,
            'norm_qsar_aih': norm_qsar_aih
        }
        
        return combined_scores, normalized_properties

    def score(self, smiles_list: List[str], reference_smiles: Optional[str] = None) -> ScoringResults:
        """
        Score a list of molecules using configured metrics.
        
        Args:
            smiles_list: List of SMILES strings to score
            reference_smiles: Reference SMILES for similarity calculation (overrides the one provided at initialization)
            
        Returns:
            ScoringResults containing combined scores and individual property scores (both raw and normalized)
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
        qsar_aih_scores = self.get_qsar_aih_predictions(smiles_list) if self.use_qsar_aih else [0.0] * len(smiles_list)
        
        # Convert to numpy arrays for calculations
        similarity_array = np.array(similarity_scores)
        logp_array = np.array(logp_scores)
        qed_array = np.array(qed_scores)
        qsar_hsec_array = np.array(qsar_hsec_scores)
        qsar_aih_array = np.array(qsar_aih_scores)
        
        # Use score_from_properties to get combined scores and normalized properties
        combined_scores, norm_props = self.score_from_properties(
            logp_values=logp_array,
            similarity_values=similarity_array,
            qsar_hsec_values=qsar_hsec_array,
            qsar_aih_values=qsar_aih_array,
            qed_values=qed_array
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
            qsar_aih_scores=qsar_aih_scores,
            norm_similarity_scores=norm_props['norm_similarity'].tolist(),
            norm_logp_scores=norm_props['norm_logp'].tolist(),
            norm_qed_scores=norm_props['norm_qed'].tolist(),
            norm_qsar_hsec_scores=norm_props['norm_qsar_hsec'].tolist(),
            norm_qsar_aih_scores=norm_props['norm_qsar_aih'].tolist()
        )

#### scorer end

def test_qsar_prediction(smiles: str, endpoint: str):
    scorer = MoleculeScorer(use_qsar_hsec=True, use_qsar_aih=True)
    if "hsec" in endpoint.lower():
        scores = scorer.get_qsar_hsec_predictions([smiles])
    else:
        scores = scorer.get_qsar_aih_predictions([smiles])
    print(f"SMILES: {smiles}")
    print(f"Endpoint: {endpoint}")
    print(f"Score: {scores[0]}")
    return scores[0]