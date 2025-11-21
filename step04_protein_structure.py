"""
Druggability Assessment Pipeline - Step 4: Protein Structure Retrieval
=======================================================================

This script takes Step3.csv and retrieves protein structures:
- Search PDB (Protein Data Bank) for experimental structures
- Query AlphaFold Database for predicted structures
- Download structure files (.pdb format)
- Analyze structure quality and coverage

OUTPUT: Step4.csv with all Step3 fields + protein structure data
"""

import requests
import pandas as pd
import time
from typing import List, Dict, Optional, Tuple
from pathlib import Path
import json
import os


class ProteinStructureRetriever:
    """Retrieve and analyze protein structures from PDB and AlphaFold."""
    
    def __init__(self, structures_dir: str = "protein_structures"):
        self.pdb_search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
        self.pdb_data_url = "https://data.rcsb.org/rest/v1/core/entry"
        self.pdb_download_url = "https://files.rcsb.org/download"
        self.alphafold_url = "https://alphafold.ebi.ac.uk/api"
        self.alphafold_files_url = "https://alphafold.ebi.ac.uk/files"
        
        # Create directories for structures
        self.structures_dir = Path(structures_dir)
        self.pdb_dir = self.structures_dir / "pdb_experimental"
        self.alphafold_dir = self.structures_dir / "alphafold_predicted"
        
        self.pdb_dir.mkdir(parents=True, exist_ok=True)
        self.alphafold_dir.mkdir(parents=True, exist_ok=True)
    
    def search_pdb_by_uniprot(self, uniprot_id: str) -> List[Dict]:
        """
        Search PDB for structures using UniProt ID.
        
        Args:
            uniprot_id: UniProt accession
            
        Returns:
            List of PDB entries with metadata
        """
        if not uniprot_id or pd.isna(uniprot_id):
            return []
        
        # PDB search query
        query = {
            "query": {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                    "operator": "exact_match",
                    "value": uniprot_id
                }
            },
            "return_type": "entry",
            "request_options": {
                "return_all_hits": True
            }
        }
        
        response = requests.post(
            self.pdb_search_url,
            json=query,
            headers={"Content-Type": "application/json"},
            timeout=15
        )
        
        if response.status_code == 200:
            data = response.json()
            if 'result_set' in data:
                return [hit['identifier'] for hit in data['result_set']]
        
        return []
    
    def get_pdb_metadata(self, pdb_id: str) -> Dict:
        """
        Get detailed metadata for a PDB entry.
        
        Args:
            pdb_id: PDB identifier
            
        Returns:
            Dictionary with PDB metadata
        """
        result = {
            'pdb_id': pdb_id,
            'title': None,
            'resolution': None,
            'experimental_method': None,
            'release_date': None,
            'organism': None
        }
        
        url = f"{self.pdb_data_url}/{pdb_id}"
        
        response = requests.get(url, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            
            # Get title
            if 'struct' in data:
                result['title'] = data['struct'].get('title')
            
            # Get experimental method and resolution
            if 'exptl' in data and len(data['exptl']) > 0:
                result['experimental_method'] = data['exptl'][0].get('method')
            
            if 'rcsb_entry_info' in data:
                entry_info = data['rcsb_entry_info']
                result['resolution'] = entry_info.get('resolution_combined', [None])[0]
                result['release_date'] = entry_info.get('deposit_date')
            
            # Get organism
            if 'rcsb_entity_source_organism' in data and len(data['rcsb_entity_source_organism']) > 0:
                result['organism'] = data['rcsb_entity_source_organism'][0].get('ncbi_scientific_name')
        
        return result
    
    def download_pdb_structure(self, pdb_id: str, output_dir: Path) -> Optional[str]:
        """
        Download PDB structure file.
        
        Args:
            pdb_id: PDB identifier
            output_dir: Directory to save file
            
        Returns:
            Path to downloaded file or None
        """
        url = f"{self.pdb_download_url}/{pdb_id}.pdb"
        output_file = output_dir / f"{pdb_id}.pdb"
        
        # Skip if already downloaded
        if output_file.exists():
            return str(output_file)
        
        response = requests.get(url, timeout=30)
        
        if response.status_code == 200:
            with open(output_file, 'w') as f:
                f.write(response.text)
            return str(output_file)
        
        return None
    
    def query_alphafold(self, uniprot_id: str) -> Dict:
        """
        Query AlphaFold database for predicted structure.
        
        Args:
            uniprot_id: UniProt accession
            
        Returns:
            Dictionary with AlphaFold data
        """
        result = {
            'found': False,
            'uniprot_id': uniprot_id,
            'alphafold_version': None,
            'model_url': None,
            'confidence_url': None,
            'pae_url': None,
            'avg_plddt': None,
            'structure_quality': None
        }
        
        if not uniprot_id or pd.isna(uniprot_id):
            return result
        
        # Query AlphaFold API
        url = f"{self.alphafold_url}/prediction/{uniprot_id}"
        
        response = requests.get(url, timeout=15)
        
        if response.status_code == 200:
            data = response.json()
            
            if len(data) > 0:
                entry = data[0]
                result['found'] = True
                result['alphafold_version'] = entry.get('modelCreatedDate')
                
                # Build file URLs
                result['model_url'] = f"{self.alphafold_files_url}/AF-{uniprot_id}-F1-model_v4.pdb"
                result['confidence_url'] = f"{self.alphafold_files_url}/AF-{uniprot_id}-F1-confidence_v4.json"
                result['pae_url'] = f"{self.alphafold_files_url}/AF-{uniprot_id}-F1-predicted_aligned_error_v4.json"
                
                # Get confidence scores if available
                if 'paeDocUrl' in entry:
                    try:
                        pae_response = requests.get(entry['paeDocUrl'], timeout=10)
                        if pae_response.status_code == 200:
                            pae_data = pae_response.json()
                            if 'avgPlddt' in pae_data:
                                result['avg_plddt'] = pae_data['avgPlddt']
                    except:
                        pass
        
        return result
    
    def download_alphafold_structure(self, uniprot_id: str, output_dir: Path) -> Optional[str]:
        """
        Download AlphaFold predicted structure.
        
        Args:
            uniprot_id: UniProt accession
            output_dir: Directory to save file
            
        Returns:
            Path to downloaded file or None
        """
        output_file = output_dir / f"AF-{uniprot_id}-F1-model_v4.pdb"
        
        # Skip if already downloaded
        if output_file.exists():
            return str(output_file)
        
        url = f"{self.alphafold_files_url}/AF-{uniprot_id}-F1-model_v4.pdb"
        
        response = requests.get(url, timeout=30)
        
        if response.status_code == 200:
            with open(output_file, 'w') as f:
                f.write(response.text)
            return str(output_file)
        
        return None
    
    def assess_structure_quality(self, avg_plddt: Optional[float]) -> str:
        """
        Assess AlphaFold structure quality based on pLDDT score.
        
        pLDDT (predicted Local Distance Difference Test) ranges 0-100:
        - >90: Very high confidence
        - 70-90: Confident
        - 50-70: Low confidence
        - <50: Very low confidence
        
        Args:
            avg_plddt: Average pLDDT score
            
        Returns:
            Quality assessment string
        """
        if avg_plddt is None:
            return 'Unknown'
        
        if avg_plddt >= 90:
            return 'Very High (>90)'
        elif avg_plddt >= 70:
            return 'High (70-90)'
        elif avg_plddt >= 50:
            return 'Low (50-70)'
        else:
            return 'Very Low (<50)'
    
    def process_gene_row(self, row: pd.Series) -> Dict:
        """
        Process a single gene and retrieve protein structures.
        
        Args:
            row: Row from Step3 DataFrame
            
        Returns:
            Dictionary with structure data
        """
        gene_symbol = row['gene_symbol']
        uniprot_id = row.get('uniprot_swissprot')
        
        print(f"  {gene_symbol}...", end=' ')
        
        result = {
            'gene_symbol': gene_symbol,
            
            # PDB experimental structures
            'has_pdb_structure': False,
            'num_pdb_structures': 0,
            'pdb_ids': None,
            'best_pdb_id': None,
            'best_pdb_resolution': None,
            'best_pdb_method': None,
            'pdb_structure_file': None,
            
            # AlphaFold predicted structure
            'has_alphafold_structure': False,
            'alphafold_version': None,
            'alphafold_avg_plddt': None,
            'alphafold_quality': None,
            'alphafold_structure_file': None,
            
            # Summary
            'structure_source': None,  # PDB, AlphaFold, or None
            'structure_file': None,
            'structure_quality_score': None
        }
        
        if not uniprot_id or pd.isna(uniprot_id):
            print("✗ No UniProt ID")
            return result
        
        # Search PDB
        pdb_ids = self.search_pdb_by_uniprot(uniprot_id)
        
        if pdb_ids:
            result['has_pdb_structure'] = True
            result['num_pdb_structures'] = len(pdb_ids)
            result['pdb_ids'] = '|'.join(pdb_ids[:10])  # Limit to 10
            
            # Get metadata for best structure (first one, usually highest quality)
            best_pdb_metadata = self.get_pdb_metadata(pdb_ids[0])
            result['best_pdb_id'] = best_pdb_metadata['pdb_id']
            result['best_pdb_resolution'] = best_pdb_metadata['resolution']
            result['best_pdb_method'] = best_pdb_metadata['experimental_method']
            
            # Download structure
            pdb_file = self.download_pdb_structure(pdb_ids[0], self.pdb_dir)
            if pdb_file:
                result['pdb_structure_file'] = pdb_file
                result['structure_source'] = 'PDB'
                result['structure_file'] = pdb_file
                
                # Quality score based on resolution
                if result['best_pdb_resolution']:
                    if result['best_pdb_resolution'] <= 2.0:
                        result['structure_quality_score'] = 'Excellent (<2.0Å)'
                    elif result['best_pdb_resolution'] <= 2.5:
                        result['structure_quality_score'] = 'Good (2.0-2.5Å)'
                    elif result['best_pdb_resolution'] <= 3.0:
                        result['structure_quality_score'] = 'Moderate (2.5-3.0Å)'
                    else:
                        result['structure_quality_score'] = f'Low (>{result["best_pdb_resolution"]:.1f}Å)'
            
            print(f"PDB:{len(pdb_ids)} Best:{pdb_ids[0]} Res:{result['best_pdb_resolution']}Å", end=' ')
        else:
            print("PDB:0", end=' ')
        
        time.sleep(0.3)
        
        # Query AlphaFold
        alphafold_data = self.query_alphafold(uniprot_id)
        
        if alphafold_data['found']:
            result['has_alphafold_structure'] = True
            result['alphafold_version'] = alphafold_data['alphafold_version']
            result['alphafold_avg_plddt'] = alphafold_data['avg_plddt']
            result['alphafold_quality'] = self.assess_structure_quality(alphafold_data['avg_plddt'])
            
            # Download structure if no PDB structure
            if not result['has_pdb_structure']:
                af_file = self.download_alphafold_structure(uniprot_id, self.alphafold_dir)
                if af_file:
                    result['alphafold_structure_file'] = af_file
                    result['structure_source'] = 'AlphaFold'
                    result['structure_file'] = af_file
                    result['structure_quality_score'] = result['alphafold_quality']
            else:
                # Still download AlphaFold for comparison
                af_file = self.download_alphafold_structure(uniprot_id, self.alphafold_dir)
                if af_file:
                    result['alphafold_structure_file'] = af_file
            
            print(f"AF:✓ pLDDT:{alphafold_data['avg_plddt']:.1f}" if alphafold_data['avg_plddt'] else "AF:✓", end=' ')
        else:
            print("AF:✗", end=' ')
        
        time.sleep(0.3)
        
        # Summary
        if result['structure_source']:
            print(f"✓ {result['structure_source']}")
        else:
            print("✗ No structure")
        
        return result
    
    def process_step3_file(self, input_file: str = "Step3.csv") -> pd.DataFrame:
        """
        Process Step3.csv and add protein structure data.
        
        Args:
            input_file: Path to Step3.csv
            
        Returns:
            DataFrame with all Step3 columns + structure data
        """
        print("="*70)
        print("PROTEIN STRUCTURE RETRIEVAL - STEP 4")
        print("="*70)
        
        # Read Step3.csv
        print(f"\nReading {input_file}...")
        step3_df = pd.read_csv(input_file)
        print(f"  ✓ Loaded {len(step3_df)} genes")
        
        print(f"\nStructure directories:")
        print(f"  PDB structures:       {self.pdb_dir}")
        print(f"  AlphaFold structures: {self.alphafold_dir}")
        
        # Process each gene
        print(f"\nRetrieving structures for {len(step3_df)} genes...")
        print("Format: Gene... PDB:X Best:ID Res:YÅ AF:✓ pLDDT:Z ✓ Source\n")
        
        structure_results = []
        
        for idx, row in step3_df.iterrows():
            print(f"[{idx+1}/{len(step3_df)}]", end=' ')
            result = self.process_gene_row(row)
            structure_results.append(result)
        
        # Create structure data DataFrame
        structure_df = pd.DataFrame(structure_results)
        
        # Merge with Step3 data
        structure_df = structure_df.drop(columns=['gene_symbol'])
        step4_df = pd.concat([step3_df, structure_df], axis=1)
        
        # Print summary
        self.print_summary(step4_df)
        
        return step4_df
    
    def print_summary(self, df: pd.DataFrame):
        """Print summary of protein structure retrieval."""
        print("\n" + "="*70)
        print("PROTEIN STRUCTURE RETRIEVAL SUMMARY")
        print("="*70)
        
        total_genes = len(df)
        
        # Count structures
        has_pdb = (df['has_pdb_structure'] == True).sum()
        has_alphafold = (df['has_alphafold_structure'] == True).sum()
        has_any_structure = (df['structure_source'].notna()).sum()
        no_structure = total_genes - has_any_structure
        
        print(f"\nTotal genes analyzed: {total_genes}")
        print(f"\nStructure Availability:")
        print(f"  ✓ PDB experimental:       {has_pdb} ({has_pdb/total_genes*100:.1f}%)")
        print(f"  ✓ AlphaFold predicted:    {has_alphafold} ({has_alphafold/total_genes*100:.1f}%)")
        print(f"  ✓ Any structure:          {has_any_structure} ({has_any_structure/total_genes*100:.1f}%)")
        print(f"  ✗ No structure:           {no_structure} ({no_structure/total_genes*100:.1f}%)")
        
        # PDB statistics
        if has_pdb > 0:
            total_pdb_structures = df['num_pdb_structures'].sum()
            avg_structures = df[df['has_pdb_structure'] == True]['num_pdb_structures'].mean()
            
            print(f"\n📊 PDB Statistics:")
            print(f"  Total PDB structures:     {int(total_pdb_structures)}")
            print(f"  Avg structures/gene:      {avg_structures:.1f}")
            
            # Resolution distribution
            resolutions = df[df['best_pdb_resolution'].notna()]['best_pdb_resolution']
            if len(resolutions) > 0:
                excellent = (resolutions <= 2.0).sum()
                good = ((resolutions > 2.0) & (resolutions <= 2.5)).sum()
                moderate = ((resolutions > 2.5) & (resolutions <= 3.0)).sum()
                low = (resolutions > 3.0).sum()
                
                print(f"\n  Resolution Distribution:")
                print(f"    Excellent (<2.0Å):      {excellent}")
                print(f"    Good (2.0-2.5Å):        {good}")
                print(f"    Moderate (2.5-3.0Å):    {moderate}")
                print(f"    Low (>3.0Å):            {low}")
                print(f"    Average resolution:     {resolutions.mean():.2f}Å")
            
            # Experimental methods
            methods = df[df['best_pdb_method'].notna()]['best_pdb_method'].value_counts()
            if len(methods) > 0:
                print(f"\n  Experimental Methods:")
                for method, count in methods.items():
                    print(f"    {method:30s}: {count}")
        
        # AlphaFold statistics
        if has_alphafold > 0:
            plddts = df[df['alphafold_avg_plddt'].notna()]['alphafold_avg_plddt']
            
            if len(plddts) > 0:
                very_high = (plddts >= 90).sum()
                high = ((plddts >= 70) & (plddts < 90)).sum()
                low = ((plddts >= 50) & (plddts < 70)).sum()
                very_low = (plddts < 50).sum()
                
                print(f"\n📊 AlphaFold Quality (pLDDT scores):")
                print(f"  Very High (>90):          {very_high} ({very_high/has_alphafold*100:.1f}%)")
                print(f"  High (70-90):             {high} ({high/has_alphafold*100:.1f}%)")
                print(f"  Low (50-70):              {low} ({low/has_alphafold*100:.1f}%)")
                print(f"  Very Low (<50):           {very_low} ({very_low/has_alphafold*100:.1f}%)")
                print(f"  Average pLDDT:            {plddts.mean():.1f}")
        
        # Structure source priority
        pdb_primary = (df['structure_source'] == 'PDB').sum()
        af_primary = (df['structure_source'] == 'AlphaFold').sum()
        
        print(f"\n📁 Structure Files Downloaded:")
        print(f"  Primary source PDB:       {pdb_primary}")
        print(f"  Primary source AlphaFold: {af_primary}")
        print(f"  Total files downloaded:   {pdb_primary + af_primary}")
        
        # Top genes with structures
        print(f"\n🎯 Genes with Best PDB Structures (by resolution):")
        best_structures = df[df['has_pdb_structure'] == True].sort_values(
            'best_pdb_resolution'
        ).head(10)
        
        if len(best_structures) > 0:
            for idx, row in best_structures.iterrows():
                pdb_id = row.get('best_pdb_id', 'N/A')
                resolution = row.get('best_pdb_resolution', 0)
                method = row.get('best_pdb_method', 'N/A')
                print(f"  • {row['gene_symbol']:12s} - {pdb_id} ({resolution:.2f}Å, {method})")
        else:
            print("  No PDB structures found.")
        
        print("="*70)


# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    
    # Check if Step3.csv exists
    if not Path("Step3.csv").exists():
        print("❌ Error: Step3.csv not found!")
        print("Please run Step 3 first to generate Step3.csv")
        exit(1)
    
    # Create protein structure retriever
    retriever = ProteinStructureRetriever(structures_dir="protein_structures")
    
    # Process Step3.csv
    step4_df = retriever.process_step3_file("Step3.csv")
    
    # Save to Step4.csv
    output_path = Path("Step4.csv")
    step4_df.to_csv(output_path, index=False)
    
    print(f"\n✅ Results saved to: {output_path.absolute()}")
    print(f"\nTotal columns: {len(step4_df.columns)}")
    print(f"Total rows: {len(step4_df)}")
    
    # Display sample
    print("\n📋 Sample Structure Data (first 5 genes):")
    structure_columns = ['gene_symbol', 'has_pdb_structure', 'num_pdb_structures',
                        'best_pdb_id', 'best_pdb_resolution', 'has_alphafold_structure',
                        'structure_source']
    if all(col in step4_df.columns for col in structure_columns):
        print(step4_df[structure_columns].head(5).to_string(index=False))
    
    print("\n✅ Step 4 Complete!")
    print("\n📁 Structure files saved in:")
    print(f"   • protein_structures/pdb_experimental/")
    print(f"   • protein_structures/alphafold_predicted/")
    print("\n💡 You can now use these structure files for:")
    print("   • Binding pocket prediction (Step 5)")
    print("   • Molecular docking studies")
    print("   • Structure-based drug design")
    print("   • Visualization in PyMOL, ChimeraX, etc.")
    print("\nNext: Proceed to Step 5 (Binding Pocket Prediction)")