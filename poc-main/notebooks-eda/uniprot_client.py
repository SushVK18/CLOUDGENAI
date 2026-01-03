#   """ Capgeminis Bioengineering Toolkit """
#   @@ All Rights Reserved @@
#   Author: ~fabigr8~ 
#   This file is part of the Capgeminis Bioengineering Toolkit
#   Please read the accompanying LICENSE file
#
#   not included in the toolkit are external dependencies
#
#

import os
import re
import requests
from requests.adapters import HTTPAdapter, Retry
from datetime import datetime
from typing import Dict, List, Optional, Union


class UniProtQuerier:
    """Class to query UniProt database and download results."""
    
    def __init__(self, output_dir: str = "uniprot_results"):
        """
        Initialize UniProt querier.
        
        Args:
            output_dir (str): Directory to store downloaded results
        """
        self.base_url = "https://rest.uniprot.org/uniprotkb/search"
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        # Setup session with retry logic
        self.session = requests.Session()
        retries = Retry(total=5, backoff_factor=0.25, 
                       status_forcelist=[500, 502, 503, 504])
        self.session.mount("https://", HTTPAdapter(max_retries=retries))
        self.re_next_link = re.compile(r'<(.+)>; rel="next"')
        
    def _get_next_link(self, headers: Dict) -> Optional[str]:
        """Extract next link from response headers."""
        if "Link" in headers:
            match = self.re_next_link.match(headers["Link"])
            if match:
                return match.group(1)
        return None
        
    def _get_batch(self, batch_url: str) -> tuple:
        """Get a batch of results and yield response with total count."""
        while batch_url:
            response = self.session.get(batch_url)
            response.raise_for_status()
            total = response.headers.get("x-total-results", 0)
            yield response, total
            batch_url = self._get_next_link(response.headers)
        
    def build_query(self, 
                    gene: str = None, 
                    taxonomy: str = None, 
                    keyword: str = None,
                    entry: Union[str, List[str]] = None) -> str:
        
        """Build query string from parameters."""
        query_parts = []
        
        if entry:
            if isinstance(entry, list):
                # Handle list of entries using OR operator
                entry_parts = [f'accession:{e}' for e in entry]
                query_parts.append(f"({' OR '.join(entry_parts)})")
            else:
                # Single entry
                query_parts.append(f'accession:{entry}')
        if gene:
            query_parts.append(f'gene:"{gene}"')
        if taxonomy:
            # Extract taxonomy ID if provided in format "Name [ID]"
            if "[" in taxonomy and "]" in taxonomy:
                tax_id = taxonomy.split("[")[1].strip("]")
                query_parts.append(f'taxonomy_id:{tax_id}')
            else:
                query_parts.append(f'organism:"{taxonomy}"')
        if keyword:
            # Extract keyword ID if provided in format "Name [ID]"
            if "[" in keyword and "]" in keyword:
                kw_id = keyword.split("[")[1].strip("]")
                query_parts.append(f'keyword:{kw_id}')
            else:
                query_parts.append(f'keyword:"{keyword}"')
                
        return " AND ".join(query_parts)

    def query(self, 
        gene: Optional[str] = None, 
        taxonomy: Optional[str] = None, 
        keyword: Optional[str] = None,
        entry: Optional[Union[str, List[str]]] = None,
        format: str = "json",
        size: int = 500,
        max_entries: Optional[int] = None) -> Dict:
        """
        Query UniProt database.
        
        Args:
            gene (str): Gene name
            taxonomy (str): Taxonomy name or ID
            keyword (str): Keyword or keyword ID
            entry (str or list): UniProt accession ID or list of IDs
            format (str): Result format (json, tsv, fasta)
            size (int): Number of results per page (max 500)
            max_entries (int, optional): Maximum number of entries to fetch. If None, fetch all.
            
        Returns:
            Dict: Query results
        """
        query_string = self.build_query(gene, taxonomy, keyword, entry)
    
        params = {
            "query": query_string,
            "format": format,
            "size": size
        }
        
        url = f"{self.base_url}?{'&'.join(f'{k}={v}' for k, v in params.items())}"
        print(f"Requesting URL: {url}")
        all_results = []
        
        for batch, total in self._get_batch(url):
            if format == "json":
                batch_data = batch.json()
                if "results" in batch_data:
                    all_results.extend(batch_data["results"])
            else:
                all_results.extend(batch.text.splitlines()[1:])  # Skip header
                
            print(f"Fetched {len(all_results)} / {total} results...")
            
            # Check if we've reached the max_entries limit
            if max_entries and len(all_results) >= max_entries:
                all_results = all_results[:max_entries]  # Trim to exact limit
                print(f"Reached maximum entries limit of {max_entries}")
                break

        return {
            "query": query_string,
            "results": all_results
        }
    
    def _create_subfolder_name(self,
                       query_string: str, 
                       timestamp: str) -> str:
    
        """Create a subfolder name from query string and timestamp."""
        # Check if this is an accession query with multiple entries
        if query_string.startswith('(accession:') and ' OR ' in query_string:
            # For multiple accessions, use a generic name with count
            accession_count = query_string.count('accession:')
            clean_query = f"uniprot_entries_{accession_count}"
        else:
            # For other queries, clean up as before but be more aggressive with shortening
            clean_query = query_string.replace('"', '').replace(':', '_')
            clean_query = clean_query.replace(' AND ', '_').replace(' ', '_')
            clean_query = clean_query.replace('taxonomy', 'tax').replace('keyword', 'kw')
            clean_query = clean_query.replace('accession', 'acc')
            
            # Limit length of folder name
            if len(clean_query) > 30:
                # Take first 27 chars and add ellipsis
                clean_query = clean_query[:27] + "..."
        
        folder_name = f"{timestamp}_{clean_query}"
        print(f"Created folder: {folder_name}")
        return folder_name

    def _save_query_log(self, query_dir: str, query_string: str, results_count: int, timestamp: str) -> str:
        """
        Save query details to a log file.
        
        Args:
            query_dir: Directory to save the log file
            query_string: The query string used
            results_count: Number of results found
            timestamp: Timestamp of the query
            
        Returns:
            str: Path to the log file
        """
        log_file = os.path.join(query_dir, f"uniprot_results_query_log_{timestamp}.txt")
        
        with open(log_file, "w", encoding="utf-8") as f:
            f.write(f"Query timestamp: {timestamp}\n")
            f.write(f"Query string: {query_string}\n")
            f.write(f"Results count: {results_count}\n")
            
            # If it's a list of accessions, extract and list them
            if query_string.startswith('(accession:') and ' OR ' in query_string:
                accessions = []
                parts = query_string.strip('()').split(' OR ')
                for part in parts:
                    if part.startswith('accession:'):
                        acc = part.replace('accession:', '')
                        accessions.append(acc)
                
                f.write(f"\nQueried accessions ({len(accessions)}):\n")
                for i, acc in enumerate(accessions, 1):
                    f.write(f"{i}. {acc}\n")

            # add last line
            f.write("\n")
            f.write(f"generated via by Bioengineering Toolkit @fabigr8")
        
        return log_file

    def download_results(self, 
                   query_results: Dict, 
                   formats: List[str] = ["json", "fasta", "tsv"]) -> Dict[str, List[str]]:
        """
        Download query results in specified formats.
        
        Args:
            query_results (Dict): Query results from query method
            formats (List[str]): List of formats to download
            
        Returns:
            Dict[str, List[str]]: Dictionary mapping formats to lists of saved file paths
        """
        downloaded_files = {fmt: [] for fmt in formats}
        downloaded_files['log'] = []  # Add log file tracking
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        if not query_results.get("results"):
            return downloaded_files

        # Create query-specific subfolder
        subfolder_name = self._create_subfolder_name(query_results["query"], timestamp)
        query_dir = os.path.join(self.output_dir, subfolder_name)
        os.makedirs(query_dir, exist_ok=True)
        
        # Save query details to log file
        log_file = self._save_query_log(
            query_dir, 
            query_results["query"], 
            len(query_results["results"]), 
            timestamp
        )
        downloaded_files['log'].append(log_file)

        # For JSON and TSV formats, save all results in one file
        for fmt in [f for f in formats if f in ["json", "tsv"]]:
            filename = f"uniprot_results_{timestamp}.{fmt}"
            filepath = os.path.join(query_dir, filename)
            
            params = {
                "query": query_results["query"],
                "format": fmt
            }
            response = requests.get(self.base_url, params=params)
            response.raise_for_status()
            
            with open(filepath, "w", encoding="utf-8") as f:
                f.write(response.text)
            
            downloaded_files[fmt].append(filepath)

        # For FASTA format, create separate files for each entry
        if "fasta" in formats and query_results["results"]:
            for entry in query_results["results"]:
                if isinstance(entry, dict):
                    uniprot_id = entry.get("primaryAccession")
                    if not uniprot_id:
                        continue
                        
                    # Use only the UniProt ID as the filename
                    filename = f"{uniprot_id}.fasta"
                    filepath = os.path.join(query_dir, filename)
                    
                    params = {
                        "query": f"accession:{uniprot_id}",
                        "format": "fasta"
                    }
                    response = requests.get(self.base_url, params=params)
                    response.raise_for_status()
                    
                    with open(filepath, "w", encoding="utf-8") as f:
                        f.write(response.text)
                    
                    downloaded_files["fasta"].append(filepath)
                
        return downloaded_files

  
if __name__ == '__main__':
    # Example usage with pagination and max entries
    querier = UniProtQuerier(output_dir="data/raw/uniprot_results")

    # Query example with limit of 1000 results
    # edit example by edding gene="P5CS1",gene="OsDREB1F",
    #results = querier.query(
    #    taxonomy="Oryza sativa (Rice) [4530]",
    #    keyword="Stress response [KW-0346]",
    #    max_entries=10  # limit to 1000 entries
    #)
    # or directly by entry
    #results = querier.query(
    #    entry=["O49686", "Q9CAJ0"]  # list of IDs 
    #)

    # Download results
    print(f"Total results found: {len(results['results'])}")
    downloaded_files = querier.download_results(results, formats=["json", "fasta", "tsv"])

    # Print file locations
    for fmt, filepath in downloaded_files.items():
        print(f"Downloaded {fmt} file to: {filepath}")