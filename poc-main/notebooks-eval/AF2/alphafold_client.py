import os
import requests
import argparse
from Bio import SeqIO
import logging
from tqdm import tqdm
import time

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("alphafold_batch.log"),
        logging.StreamHandler()
    ]
)

# Default AlphaFold2 server address
DEFAULT_AF2_HOST = 'http://localhost:8000'

# Function to read a FASTA file and return the protein sequence
def fasta_to_str(file_path):
    # Use Biopython to read the FASTA file
    with open(file_path, 'r') as file:
        for record in SeqIO.parse(file, 'fasta'):
            return str(record.seq)  # Return the sequence as a string

# Function to call AlphaFold2 for a single FASTA file
def call_alphafold(fasta_file, 
                   af2_nim_host, 
                   output_dir=None, 
                   databases=['uniref90', 'mgnify', 'small_bfd'], 
                   msa_algorithm='jackhmmer', 
                   e_value=0.0001, 
                   msa_iterations=1, 
                   relax_prediction=False):
    """
    Process a single FASTA file with AlphaFold2
    
    Args:
        fasta_file (str): Path to the FASTA file
        af2_nim_host (str): AlphaFold2 server address
        output_dir (str, optional): Directory to save output PDB files. If None, won't save to file.
        databases (list): Databases to use for MSA
        msa_algorithm (str): Algorithm for MSA generation
        e_value (float): E-value threshold
        msa_iterations (int): Number of MSA iterations
        relax_prediction (bool): Whether to relax the prediction
        
    Returns:
        tuple: (protein_id, success_status, message, structure_data)
    """
    # Extract protein ID from filename
    protein_id = os.path.splitext(os.path.basename(fasta_file))[0]
    
    try:
        # Read the protein sequence from the FASTA file
        protein_sequence = fasta_to_str(fasta_file)
        
        if not protein_sequence:
            return protein_id, False, "Empty or invalid FASTA file", None
        
        # Set default databases if none provided
        if databases is None:
            databases = ['uniref90', 'mgnify', 'small_bfd']
        
        # Query AlphaFold2 with the protein sequence
        logging.info(f"Processing {protein_id} ({len(protein_sequence)} residues)")
        
        af2_response = requests.post(
            f'{af2_nim_host}/protein-structure/alphafold2/predict-structure-from-sequence',
            json={
                'sequence': protein_sequence,
                'databases': databases,
                'msa_algorithm': msa_algorithm,
                'e_value': e_value,
                'bit_score': -1,  # -1 means to fallback to the e-value
                'msa_iterations': msa_iterations,
                'relax_prediction': relax_prediction,
            }
        )
        
        # Check for a successful response
        if af2_response.status_code == 200:
            result = af2_response.json()
            structure_data = result[0]
            
            # Save the result to a PDB file if output_dir is provided
            if output_dir:
                os.makedirs(output_dir, exist_ok=True)
                output_file = os.path.join(output_dir, f"{protein_id}_relaxed.pdb")
                with open(output_file, 'w') as file:
                    file.write(structure_data)
                
            return protein_id, True, "Success", structure_data
        else:
            error_msg = f"HTTP Error: {af2_response.status_code} - {af2_response.text}"
            logging.error(f"Failed processing {protein_id}: {error_msg}")
            return protein_id, False, error_msg, None
            
    except Exception as e:
        error_msg = str(e)
        logging.error(f"Error processing {protein_id}: {error_msg}")
        return protein_id, False, error_msg, None

class AlphaFoldBatchProcessor:
    """
    Batch processor for running multiple protein sequences through AlphaFold2
    """
    
    def __init__(self, af2_host=DEFAULT_AF2_HOST, databases=None, 
                 msa_algorithm='jackhmmer', e_value=0.0001, 
                 msa_iterations=1, relax_prediction=False):
        """
        Initialize the batch processor
        
        Args:
            af2_host (str): AlphaFold2 server address
            databases (list): Databases to use for MSA
            msa_algorithm (str): Algorithm for MSA generation
            e_value (float): E-value threshold
            msa_iterations (int): Number of MSA iterations
            relax_prediction (bool): Whether to relax the prediction
        """
        self.af2_host = af2_host
        self.databases = databases if databases else ['uniref90', 'mgnify', 'small_bfd']
        self.msa_algorithm = msa_algorithm
        self.e_value = e_value
        self.msa_iterations = msa_iterations
        self.relax_prediction = relax_prediction
        
        logging.info("AlphaFold2 Batch Processor initialized with:")
        logging.info(f"  AlphaFold2 host: {self.af2_host}")
        logging.info(f"  Databases: {self.databases}")
        logging.info(f"  MSA algorithm: {self.msa_algorithm}")
        logging.info(f"  E-value: {self.e_value}")
        logging.info(f"  MSA iterations: {self.msa_iterations}")
        logging.info(f"  Relax prediction: {self.relax_prediction}")
    
    def process_file(self, fasta_file, output_dir=None):
        """
        Process a single FASTA file
        
        Args:
            fasta_file (str): Path to the FASTA file
            output_dir (str, optional): Directory to save output PDB files
            
        Returns:
            tuple: (protein_id, success_status, message, structure_data)
        """
        return call_alphafold(
            fasta_file,
            self.af2_host,
            output_dir,
            self.databases,
            self.msa_algorithm,
            self.e_value,
            self.msa_iterations,
            self.relax_prediction
        )
    
    def process_directory(self, input_dir, output_dir=None, show_progress=True):
        """
        Process all FASTA files in a directory sequentially
        
        Args:
            input_dir (str): Directory containing FASTA files
            output_dir (str, optional): Directory to save output PDB files
            show_progress (bool): Whether to show a progress bar
            
        Returns:
            dict: Results with success and failure counts, and processed structures
        """
        # Find all FASTA files in the input directory
        fasta_files = []
        for file in os.listdir(input_dir):
            if file.endswith(('.fasta', '.fa')):
                fasta_files.append(os.path.join(input_dir, file))
        
        if not fasta_files:
            logging.warning(f"No FASTA files found in {input_dir}")
            return {"success": 0, "failure": 0, "total": 0, "structures": {}}
        
        logging.info(f"Found {len(fasta_files)} FASTA files to process")
        
        # Use the process_files method to handle the files
        return self.process_files(fasta_files, output_dir, show_progress)
        
    def process_files(self, fasta_files, output_dir=None, show_progress=True):
        """
        Process a list of FASTA files sequentially
        
        Args:
            fasta_files (list): List of paths to FASTA files
            output_dir (str, optional): Directory to save output PDB files
            show_progress (bool): Whether to show a progress bar
            
        Returns:
            dict: Results with success and failure counts, and processed structures
        """
        if not fasta_files:
            logging.warning("No FASTA files provided")
            return {"success": 0, "failure": 0, "total": 0, "structures": {}}
        
        logging.info(f"Processing {len(fasta_files)} FASTA files sequentially")
        
        # Create output directory if it doesn't exist and was specified
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
        
        results = {
            "success": 0, 
            "failure": 0, 
            "total": len(fasta_files),
            "structures": {},
            "failed_proteins": []
        }
        
        # Process files sequentially with an optional progress bar
        iterator = tqdm(fasta_files, desc="Processing proteins") if show_progress else fasta_files
        for fasta_file in iterator:
            protein_id, success, message, structure_data = self.process_file(fasta_file, output_dir)
            
            if success:
                results["success"] += 1
                results["structures"][protein_id] = structure_data
            else:
                results["failure"] += 1
                results["failed_proteins"].append((protein_id, message))
        
        # Log the summary
        logging.info(f"Processed {results['total']} proteins: "
                    f"{results['success']} successful, {results['failure']} failed")
        
        if results["failed_proteins"]:
            logging.warning("Failed proteins:")
            for protein_id, message in results["failed_proteins"]:
                logging.warning(f"  - {protein_id}: {message}")
        
        return results

def main():
    """Command line entry point"""
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Process multiple FASTA files with AlphaFold2')
    parser.add_argument('--input_dir', required=True, help='Directory containing FASTA files')
    parser.add_argument('--output_dir', required=True, help='Directory to save output PDB files')
    parser.add_argument('--af2_host', default=DEFAULT_AF2_HOST, help='AlphaFold2 server address')
    parser.add_argument('--databases', nargs='+', default=['uniref90', 'mgnify', 'small_bfd'], 
                        help='Databases to use for MSA')
    parser.add_argument('--msa_algorithm', default='jackhmmer', 
                        choices=['jackhmmer', 'hhblits'], help='Algorithm for MSA generation')
    parser.add_argument('--e_value', type=float, default=0.0001, help='E-value threshold')
    parser.add_argument('--msa_iterations', type=int, default=1, help='Number of MSA iterations')
    parser.add_argument('--relax', action='store_true', help='Relax the prediction')
    
    args = parser.parse_args()
    
    # Start timer
    start_time = time.time()
    
    # Initialize the processor
    processor = AlphaFoldBatchProcessor(
        af2_host=args.af2_host,
        databases=args.databases,
        msa_algorithm=args.msa_algorithm,
        e_value=args.e_value,
        msa_iterations=args.msa_iterations,
        relax_prediction=args.relax
    )
    
    # Process the directory
    results = processor.process_directory(
        args.input_dir,
        args.output_dir
    )
    
    # Log elapsed time
    elapsed_time = time.time() - start_time
    logging.info(f"Total processing time: {elapsed_time:.2f} seconds")
    
    return results

if __name__ == "__main__":
    main()