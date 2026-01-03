import numpy as np
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt
from rdkit.Chem.Crippen import MolLogP
import warnings
from typing import Dict, List, Callable, Union, Tuple, Optional
import time
import boto3
from botocore.exceptions import ClientError
import json
import os
import atexit
from msal import ConfidentialClientApplication, SerializableTokenCache
import requests

class MultiObjectiveOracle:
    """
    Oracle for multi-objective optimization with customizable scoring approaches.
    Supports 
    * weighted product, 
    * Harrington's desirability function, 
    * and enhanced property-specific **desirability functions** with dynamic weighting.
    
    Includes direct integration with QSAR APIs for AIH and HSEC models.
    """
    
    SCORING_METHODS = {
        "weighted_product": 1,
        "harrington": 2,
        "enhanced": 3
    }
    
    TOKEN_FILE = os.path.join(os.path.expanduser('~'), ".insilicoclitoken")
    
    def __init__(
        self,
        scoring_method: str = "weighted_product",
        # Default ranges for common properties
        logp_min: float = 1.45, # 1x STD - MEAN von eval datensatz
        logp_max: float = 3.0, # 1x STD + MEAN von eval datensatz
        mw_min: float = 335.0, # 1x STD - MEAN von eval datensatz 285
        mw_max: float = 375.0, # 1x STD + MEAN von eval datensatz 440
        qsar_aih_min: float = -50.0,
        qsar_aih_max: float = 500.0,
        qsar_hsec_min: float = -50.0,
        qsar_hsec_max: float = 450.0,
        # Weights
        logp_weight: float = 1.0,
        mw_weight: float = 1.0,
        qsar_weight: float = 1.0,
        similarity_weight: float = 0.0,  # By default, similarity is not used (weight 0)
        # Desirability configuration
        logp_optimal_min: float = 1.55,
        logp_optimal_max: float = 3.0,
        mw_optimal_min: float = 335.0, # 1x STD - MEAN von eval datensatz
        mw_optimal_max: float = 375.0, # 1x STD + MEAN von eval datensatz
        # Property maximization flags (True = higher is better)
        maximize_logp: bool = True,
        maximize_mw: bool = False,  # Usually we want lower MW
        maximize_qsar: bool = True,  # Higher QSAR scores are typically better
        # Similarity parameters
        use_similarity: bool = False,
        reference_smiles: str = None,
        similarity_min: float = 0.2,
        similarity_max: float = 0.9,
        fingerprint_radius: int = 2,
        fingerprint_size: int = 2048,
        # QSAR specific parameters
        use_qsar: bool = False,
        qsar_model: str = "AIH",  # "AIH" or "HSEC"
        secret_name: str = "PROPCALC_API_CREDS",
        tenant_id: str = "fcb2b37b-5da0-466b-9b83-0014b67a7c78",
        qsar_hsec_endpoint: str = "https://chemprop.cropkey-np.ag/activity/ica/invocations_hsec",
        qsar_aih_endpoint: str = "https://chemprop.cropkey-np.ag/activity/ica/invocations",
        # Verbose output
        verbose: bool = False
    ):
        """
        Initialize the multi-objective oracle with simplified parameters.
        """
        # Validate scoring method
        if scoring_method not in self.SCORING_METHODS:
            raise ValueError(f"Scoring method must be one of {list(self.SCORING_METHODS.keys())}")
        self.scoring_method = scoring_method
        self.verbose = verbose
        
        # Store property ranges
        self.logp_min = logp_min
        self.logp_max = logp_max
        self.mw_min = mw_min
        self.mw_max = mw_max
        self.qsar_aih_min = qsar_aih_min
        self.qsar_aih_max = qsar_aih_max
        self.qsar_hsec_min = qsar_hsec_min
        self.qsar_hsec_max = qsar_hsec_max
        
        # Store property weights
        self.logp_weight = logp_weight
        self.mw_weight = mw_weight
        self.qsar_weight = qsar_weight
        self.similarity_weight = similarity_weight
        
        # Store desirability parameters
        self.logp_optimal_min = logp_optimal_min
        self.logp_optimal_max = logp_optimal_max
        self.mw_optimal_min = mw_optimal_min
        self.mw_optimal_max = mw_optimal_max
        
        # Store maximization flags
        self.maximize_logp = maximize_logp
        self.maximize_mw = maximize_mw
        self.maximize_qsar = maximize_qsar
        
        # Setup fingerprint similarity
        self.use_similarity = use_similarity
        self.reference_smiles = reference_smiles
        self.similarity_min = similarity_min
        self.similarity_max = similarity_max
        self.fingerprint_radius = fingerprint_radius
        self.fingerprint_size = fingerprint_size
        self.reference_fingerprint = None
        
        # Initialize reference fingerprint if provided
        if self.use_similarity and self.reference_smiles:
            self._initialize_reference_fingerprint()
        
        # QSAR configuration
        self.use_qsar = use_qsar
        self.qsar_model = qsar_model.upper()
        if self.qsar_model not in ["AIH", "HSEC"]:
            raise ValueError("QSAR model must be either 'AIH' or 'HSEC'")
        
        # Consolidated ranges dictionary for internal use
        self.ranges = {
            "LogP": (self.logp_min, self.logp_max),
            "MW": (self.mw_min, self.mw_max),
            "QSAR_AIH": (self.qsar_aih_min, self.qsar_aih_max),
            "QSAR_HSEC": (self.qsar_hsec_min, self.qsar_hsec_max),
            "Similarity": (self.similarity_min, self.similarity_max)
        }
        
        # Consolidated weights dictionary for internal use
        self.weights = {
            "LogP": self.logp_weight,
            "MW": self.mw_weight,
            "QSAR_AIH": self.qsar_weight if self.qsar_model == "AIH" else 0.0,
            "QSAR_HSEC": self.qsar_weight if self.qsar_model == "HSEC" else 0.0,
            "Similarity": self.similarity_weight if self.use_similarity else 0.0
        }
        
        # Normalize weights to sum to 1
        weight_sum = sum(weight for weight in self.weights.values() if weight > 0)
        if weight_sum > 0:
            for key in self.weights:
                self.weights[key] /= weight_sum
        
        # Consolidated maximize dictionary for internal use
        self.maximize = {
            "LogP": self.maximize_logp,
            "MW": self.maximize_mw,
            "QSAR_AIH": self.maximize_qsar,
            "QSAR_HSEC": self.maximize_qsar,
            "Similarity": True  # Higher similarity is always better
        }
        
        # Consolidated desirability parameters for internal use
        self.desirability_params = {
            "LogP": {"optimal_min": self.logp_optimal_min, "optimal_max": self.logp_optimal_max, "steepness": 1.0},
            "MW": {"optimal_min": self.mw_optimal_min, "optimal_max": self.mw_optimal_max, "steepness": 1.0},
            "QSAR_AIH": {"scale_factor": 100, "steepness": 1.0},
            "QSAR_HSEC": {"scale_factor": 100, "steepness": 1.0},
            "Similarity": {"scale_factor": 0.5, "steepness": 1.0}  # Scale factor for similarity curve
        }
        
        # Initialize QSAR API if needed
        if self.use_qsar:
            # Initialize QSAR API connection
            self.secret_name = secret_name
            self.tenant_id = tenant_id
            self.qsar_hsec_endpoint = qsar_hsec_endpoint
            self.qsar_aih_endpoint = qsar_aih_endpoint
            self._initialize_api_connection()
            
            # Initialize QSAR cache
            self.qsar_cache = {}
    
    def _calculate_logp(self, smiles: str) -> float:
        """Calculate LogP for a molecule using RDKit."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 0.0
        try:
            return MolLogP(mol)
        except:
            return 0.0
    
    def _calculate_mw(self, smiles: str) -> float:
        """Calculate molecular weight for a molecule using RDKit."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 0.0
        try:
            return MolWt(mol)
        except:
            return 0.0
    
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
            if self.verbose:
                print(f"Unable to retrieve secret: {e}")
            raise e
        except Exception as e:
            if self.verbose:
                print(f"Unknown exception: Unable to retrieve secret: {e}")
            raise e
            
    def _write_cache(self, cache: SerializableTokenCache, file: str):
        """Write token cache to file if changed."""
        if cache.has_state_changed:
            with open(file, "w") as f:
                f.write(cache.serialize())
                
    def _initialize_api_connection(self):
        """Initialize the API connection and authentication for QSAR predictions."""
        try:
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
            if self.verbose:
                print(f"QSAR API connection initialized for {self.qsar_model} model")
        except Exception as e:
            if self.verbose:
                print(f"Error initializing QSAR API connection: {e}")
            self.use_qsar = False
            
    def _get_access_token(self) -> str:
        """Get access token for API authentication."""
        result = self.app.acquire_token_for_client(scopes=[f"{self.client_id}/.default"])
        access_token = result.get("access_token")
        if not access_token:
            raise Exception("Failed to acquire access token")
        return access_token
    
    def _post_data(self, bearer_token: str, url: str, data: Dict) -> Dict:
        """Post data to API endpoint."""
        headers = {
            'accept': 'application/json',
            'Content-Type': 'application/json',
            'Authorization': f'Bearer {bearer_token}'
        }
        response = requests.post(url, headers=headers, data=json.dumps(data))
        return response.json()
    
    def _get_batch_qsar_predictions(self, smiles_list: List[str], endpoint: str) -> List[float]:
        """Get QSAR predictions for a batch of SMILES strings using specific endpoint."""
        if not smiles_list:
            return []
            
        try:
            # Format input for API
            request_input = {
                "input": [
                    {
                        "id": str(i+1),
                        "structure": smiles
                    } for i, smiles in enumerate(smiles_list)
                ]
            }
            
            # Get access token and make API request
            access_token = self._get_access_token()
            response = self._post_data(access_token, endpoint, request_input)
            
            # Initialize predictions with defaults
            predictions = [0.0] * len(smiles_list)
            
            # Process response
            if isinstance(response, dict):
                data_items = response.get("data", [])
                
                # Map predictions to their indices
                for item in data_items:
                    item_id = item.get("input", {}).get("id")
                    if not item_id:
                        continue
                        
                    if item.get("status", {}).get("success", False):
                        predictions_dict = item.get("predictions", {})
                        
                        # Extract the right prediction value
                        if "pred_hsec" in predictions_dict and "hsec" in endpoint:
                            value = predictions_dict.get("pred_hsec")
                        elif "pred_aih" in predictions_dict and "aih" in endpoint:
                            value = predictions_dict.get("pred_aih")
                        else:
                            # Try to find any prediction field
                            value = None
                            for key, val in predictions_dict.items():
                                if key.startswith("pred_"):
                                    value = val
                                    break
                                    
                        # Store prediction in the right position
                        if value is not None:
                            idx = int(item_id) - 1
                            if 0 <= idx < len(predictions):
                                predictions[idx] = float(value)
            
            return predictions
            
        except Exception as e:
            if self.verbose:
                print(f"Error in batch QSAR prediction: {e}")
            return [0.0] * len(smiles_list)
    
    def get_qsar_predictions(self, smiles_list: List[str], batch_size: int = 50) -> List[float]:
        """
        Get QSAR predictions for a list of SMILES strings using the appropriate endpoint.
        """
        if not self.use_qsar or not smiles_list:
            return [0.0] * len(smiles_list)
            
        # Determine which endpoint to use based on model
        is_hsec = (self.qsar_model == "HSEC")
        endpoint = self.qsar_hsec_endpoint if is_hsec else self.qsar_aih_endpoint
        
        # Process in batches to avoid timeouts
        all_predictions = []
        
        for i in range(0, len(smiles_list), batch_size):
            batch_smiles = smiles_list[i:i+batch_size]
            batch_predictions = self._get_batch_qsar_predictions(batch_smiles, endpoint)
            all_predictions.extend(batch_predictions)
            
            # Small delay between batches
            if i + batch_size < len(smiles_list):
                time.sleep(0.2)
        
        # Ensure we have the right number of predictions
        if len(all_predictions) != len(smiles_list):
            if len(all_predictions) < len(smiles_list):
                all_predictions.extend([0.0] * (len(smiles_list) - len(all_predictions)))
            else:
                all_predictions = all_predictions[:len(smiles_list)]
                
        return all_predictions
    
    def _get_cached_qsar_prediction(self, smiles: str) -> float:
        """Get QSAR prediction for a single SMILES string (cached)."""
        if not self.use_qsar:
            return 0.0
            
        # Check cache first
        if smiles in self.qsar_cache:
            return self.qsar_cache[smiles]
            
        try:
            # Make API call for single SMILES
            predictions = self.get_qsar_predictions([smiles])
            result = predictions[0] if predictions else 0.0
            
            # Cache the result
            self.qsar_cache[smiles] = result
            return result
        except:
            self.qsar_cache[smiles] = 0.0
            return 0.0
        
    def normalize_property(self, name: str, value: float) -> float:
        """Normalize a property value to [0,1] range using min-max scaling."""
        if value is None:
            return 0.0
            
        min_val, max_val = self.ranges[name]
        normalized = (value - min_val) / (max_val - min_val)
        
        # Clip to [0,1] to handle values outside the expected range
        return max(0.0, min(1.0, normalized))
    
    def desirability_function(self, name: str, normalized_value: float) -> float:
        """
        Convert a normalized property value to a desirability score based on
        whether we want to maximize or minimize the property.
        """
        if normalized_value is None:
            return 0.0
            
        if self.maximize[name]:
            return normalized_value
        else:
            return 1.0 - normalized_value
            
    def property_specific_desirability(self, name: str, value: float) -> float:
        """
        Custom desirability functions tailored to each property type.
        For the enhanced scoring method.
        """
        if value is None:
            return 0.0
            
        if name == "LogP":
            params = self.desirability_params.get(name)
            opt_min = params["optimal_min"]
            opt_max = params["optimal_max"]
            steepness = params["steepness"]
            
            if value < opt_min:
                return np.exp(-steepness * (value - opt_min)**2 / 2)
            elif value <= opt_max:
                return 1.0
            else:
                return np.exp(-steepness * (value - opt_max)**2 / 2)
                
        elif name == "MW":
            params = self.desirability_params.get(name)
            opt_min = params["optimal_min"]
            opt_max = params["optimal_max"]
            steepness = params["steepness"]
            
            if value < opt_min:
                return np.exp(-steepness * (value - opt_min)**2 / 2)
            elif value <= opt_max:
                return 1.0
            else:
                return np.exp(-steepness * (value - opt_max)**2 / 2)
                
        elif name.startswith("QSAR_"):
            params = self.desirability_params.get(name)
            scale_factor = params["scale_factor"]
            
            # For QSAR, assuming higher is better
            return 1.0 - np.exp(-value / scale_factor)
            
        elif name == "Similarity":
            params = self.desirability_params.get(name)
            scale_factor = params["scale_factor"]
            
            # For similarity, we want to prioritize highly similar compounds
            # This function gives a steeper curve toward higher similarity values
            return 1.0 - np.exp(-value / scale_factor)
            
        else:
            # Generic desirability function for other properties
            normalized = self.normalize_property(name, value)
            return self.desirability_function(name, normalized)
    
    def compute_weighted_product(self, desirability_scores: Dict[str, float]) -> float:
        """Compute the weighted product of desirability scores."""
        product = 1.0
        for name, score in desirability_scores.items():
            if score <= 0:
                return 0.0  # If any property completely fails, the whole molecule fails
            product *= score ** self.weights.get(name, 1.0)
        return product
        
    def compute_harrington(self, desirability_scores: Dict[str, float]) -> float:
        """Compute Harrington's desirability function."""
        total_weight = sum(self.weights.values())
        product = 1.0
        for name, score in desirability_scores.items():
            if score <= 0:
                return 0.0
            product *= score ** self.weights.get(name, 1.0)
        return product ** (1.0 / total_weight) if total_weight > 0 else 0.0
        
    def compute_enhanced_old(self, property_values: Dict[str, float]) -> float:
        """
        Compute the enhanced scoring function using property-specific 
        desirability functions and dynamic weighting.
        """
        # Calculate property-specific desirability scores
        desirability_scores = {}
        for name, value in property_values.items():
            if value is None:
                return 0.0
            desirability_scores[name] = self.property_specific_desirability(name, value)
        
        # Apply dynamic weighting
        product = 1.0
        for name, score in desirability_scores.items():
            if score <= 0:
                return 0.0
            product *= (self.weights.get(name, 1.0) * score)
            
        return product

    def compute_enhanced_mean(self, property_values: Dict[str, float]) -> float:
        """
        Compute the enhanced scoring function using property-specific 
        desirability functions and arithmetic mean.
        """
        # Calculate property-specific desirability scores
        weighted_sum = 0.0
        total_weight = 0.0
        
        for name, value in property_values.items():
            if value is None:
                return 0.0
            
            desirability = self.property_specific_desirability(name, value)
            weight = self.weights.get(name, 0.0)
            
            weighted_sum += weight * desirability
            total_weight += weight
        
        if total_weight <= 0:
            return 0.0
            
        return weighted_sum / total_weight

    
    def compute_enhanced(self, property_values: Dict[str, float]) -> float:
        """
        Compute the enhanced scoring function using property-specific 
        desirability functions and dynamic weighting.
        """
        # Calculate property-specific desirability scores
        desirability_scores = {}
        for name, value in property_values.items():
            if value is None:
                return 0.0
            desirability_scores[name] = self.property_specific_desirability(name, value)
        
        # Apply weights using weighted geometric mean
        product = 1.0
        total_weight = sum(self.weights.get(name, 0.0) for name in desirability_scores)
        
        if total_weight <= 0:
            return 0.0
            
        for name, score in desirability_scores.items():
            weight = self.weights.get(name, 0.0)
            if weight > 0:
                if score <= 0:
                    return 0.0
                # Use weighted exponentiation (geometric mean)
                product *= score ** (weight / total_weight)
                
        return product
    
    def calculate_score(self, property_values: Dict[str, float]) -> float:
        """Calculate the final score based on the selected scoring method."""
        # Check for None values
        if any(value is None for value in property_values.values()):
            return 0.0
            
        if self.scoring_method == "weighted_product":
            # Normalize and get desirability scores
            desirability_scores = {
                name: self.desirability_function(name, self.normalize_property(name, value))
                for name, value in property_values.items()
            }
            return self.compute_weighted_product(desirability_scores)
            
        elif self.scoring_method == "harrington":
            # Normalize and get desirability scores
            desirability_scores = {
                name: self.desirability_function(name, self.normalize_property(name, value))
                for name, value in property_values.items()
            }
            return self.compute_harrington(desirability_scores)
            
        elif self.scoring_method == "enhanced":
            return self.compute_enhanced(property_values)
            
        else:
            raise ValueError(f"Unknown scoring method: {self.scoring_method}")
    
    def _initialize_reference_fingerprint(self):
        """Initialize the reference fingerprint for similarity calculations."""
        if not self.use_similarity or not self.reference_smiles:
            return False
            
        try:
            ref_mol = Chem.MolFromSmiles(self.reference_smiles)
            if ref_mol is None:
                if self.verbose:
                    print(f"Warning: Invalid reference SMILES: {self.reference_smiles}")
                return False
                
            # Generate fingerprint using the specified parameters
            mfpgen = Chem.rdFingerprintGenerator.GetMorganGenerator(
                radius=self.fingerprint_radius, 
                fpSize=self.fingerprint_size
            )
            self.reference_fingerprint = mfpgen.GetFingerprint(ref_mol)
            
            if self.verbose:
                print(f"Initialized reference fingerprint for similarity calculations")
            return True
            
        except Exception as e:
            if self.verbose:
                print(f"Error initializing reference fingerprint: {e}")
            return False
    
    def calculate_similarity(self, smiles: str) -> float:
        """Calculate molecular fingerprint similarity to the reference molecule."""
        if not self.use_similarity or self.reference_fingerprint is None:
            return 0.0
            
        try:
            # Parse molecule
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return 0.0
                
            # Generate fingerprint
            mfpgen = Chem.rdFingerprintGenerator.GetMorganGenerator(
                radius=self.fingerprint_radius, 
                fpSize=self.fingerprint_size
            )
            fp = mfpgen.GetFingerprint(mol)
            
            # Calculate Tanimoto similarity
            from rdkit import DataStructs
            similarity = DataStructs.TanimotoSimilarity(fp, self.reference_fingerprint)
            return similarity
            
        except:
            return 0.0
    
    def score_molecule(self, smiles: str) -> float:
        """Score a single molecule based on its SMILES string."""
        try:
            # Check if the molecule is valid
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return 0.0
                
            # Calculate all properties
            property_values = {
                "LogP": self._calculate_logp(smiles),
                "MW": self._calculate_mw(smiles)
            }
            
            # Add QSAR prediction if enabled
            if self.use_qsar:
                qsar_prop_name = f"QSAR_{self.qsar_model}"
                property_values[qsar_prop_name] = self._get_cached_qsar_prediction(smiles)
            
            # Add similarity if enabled
            if self.use_similarity and self.reference_fingerprint is not None:
                property_values["Similarity"] = self.calculate_similarity(smiles)
                    
            # Calculate final score
            return self.calculate_score(property_values)
            
        except:
            return 0.0
        
    def evaluate(self, molecules: List[str]) -> Dict[str, float]:
        """Evaluate a list of molecules and return their scores."""
        if self.verbose:
            print(f"Evaluating {len(molecules)} molecules using {self.scoring_method} method")
            if self.use_qsar:
                print(f"Including QSAR {self.qsar_model} predictions in scoring")
        
        # Handle duplicates properly
        result_list = []
        for idx, smiles in enumerate(molecules):
            score = self.score_molecule(smiles)
            result_list.append((smiles, score))
        
        # Convert to dictionary while preserving duplicates with unique keys
        results = {}
        seen_smiles = {}
        
        for idx, (smiles, score) in enumerate(result_list):
            if smiles in seen_smiles:
                # This is a duplicate SMILES
                seen_smiles[smiles] += 1
                # Add a suffix to make the key unique
                unique_key = f"{smiles}__DUP_{seen_smiles[smiles]}"
                results[unique_key] = score
            else:
                # First time seeing this SMILES
                seen_smiles[smiles] = 0
                results[smiles] = score
        
        if self.verbose and len(results) != len(molecules):
            print(f"Note: Input contained {len(molecules) - len(set(molecules))} duplicate SMILES strings")
            
        return results

    def debug_scoring(self, smiles: str) -> Dict:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"score": 0.0, "error": "Invalid SMILES"}
        
        # Calculate individual properties
        logp = self._calculate_logp(smiles)
        mw = self._calculate_mw(smiles)
        
        property_values = {
            "LogP": logp,
            "MW": mw
        }
        
        # Add similarity if enabled
        if self.use_similarity and self.reference_fingerprint is not None:
            property_values["Similarity"] = self.calculate_similarity(smiles)
        
        # Normalize properties and get desirability
        normalized = {
            name: self.normalize_property(name, value)
            for name, value in property_values.items()
        }
        
        desirability = {
            name: self.desirability_function(name, norm_value)
            for name, norm_value in normalized.items()
        }
        
        # Calculate final score
        score = self.calculate_score(property_values)

        # Enhancd final score 
        enh_score = self.compute_enhanced(property_values)
        return {
            "score": score,
            "enhanced score": enh_score,
            "properties": property_values,
            "normalized": normalized,
            "desirability": desirability,
            "weights": {k: v for k, v in self.weights.items() if v > 0}
        }

# Example usage
if __name__ == "__main__":
    # Example SMILES strings
    test_molecules = [
        "CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine
        "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"  # Ibuprofen
    ]
    
    # Basic usage
    oracle = MultiObjectiveOracle(
        scoring_method="weighted_product",
        logp_weight=1.0,
        mw_weight=1.0,
        use_qsar=True,
        qsar_model="AIH",
        verbose=True  # Set to False to suppress output
    )
    
    results = oracle.evaluate(test_molecules)
    print(results)