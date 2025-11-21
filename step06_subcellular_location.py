"""
Step 6: Subcellular Location Assessment
Reads Step5.csv, counts pockets from binding_pockets/ directory, 
adds complete location data to ALL pocket rows, saves to Step6.csv
"""

import pandas as pd
import requests
import time
import re
from typing import List, Dict, Optional
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('step6_location_analysis.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


class PocketCounter:
    """Count actual pockets from binding_pockets directory structure"""
    
    def __init__(self, base_dir: str = 'binding_pockets'):
        self.base_dir = Path(base_dir)
    
    def count_pockets_for_gene(self, gene: str) -> Dict[str, int]:
        """
        Count actual pocket files for a specific gene
        
        Args:
            gene: Gene symbol (e.g., 'ACE', 'CALCA')
        
        Returns:
            Dict with pocket count information
        """
        gene_dir = self.base_dir / gene
        
        if not gene_dir.exists():
            logger.warning(f"Directory not found for gene: {gene}")
            return {
                'Total_Pockets': 0,
                'PDB_Structures': 0,
                'Pockets_List': []
            }
        
        total_pockets = 0
        pdb_structures = 0
        all_pockets = []
        
        # Iterate through PDB output directories (e.g., 1O86_out)
        pdb_dirs = [d for d in gene_dir.iterdir() if d.is_dir() and d.name.endswith('_out')]
        pdb_structures = len(pdb_dirs)
        
        for pdb_dir in pdb_dirs:
            pdb_name = pdb_dir.name.replace('_out', '')
            pockets_dir = pdb_dir / 'pockets'
            
            if pockets_dir.exists():
                # Count pocket*_atm.pdb files
                pocket_files = list(pockets_dir.glob('pocket*_atm.pdb'))
                num_pockets = len(pocket_files)
                total_pockets += num_pockets
                
                # Extract pocket numbers
                for pf in pocket_files:
                    match = re.search(r'pocket(\d+)_atm\.pdb', pf.name)
                    if match:
                        pocket_num = int(match.group(1))
                        all_pockets.append({
                            'PDB': pdb_name,
                            'Pocket_Number': pocket_num
                        })
                
                logger.info(f"  {gene}/{pdb_name}: {num_pockets} pockets")
        
        return {
            'Total_Pockets': total_pockets,
            'PDB_Structures': pdb_structures,
            'Pockets_List': all_pockets
        }
    
    def count_all_genes(self, gene_list: List[str]) -> pd.DataFrame:
        """
        Count pockets for all genes
        
        Args:
            gene_list: List of gene symbols
        
        Returns:
            DataFrame with pocket counts per gene
        """
        logger.info("\nCounting pockets from binding_pockets/ directory...")
        logger.info("="*80)
        
        results = []
        
        for gene in sorted(gene_list):
            logger.info(f"\nCounting pockets for {gene}:")
            pocket_info = self.count_pockets_for_gene(gene)
            
            results.append({
                'Gene': gene,
                'Step6_Total_Pockets': pocket_info['Total_Pockets'],
                'Step6_PDB_Structures': pocket_info['PDB_Structures'],
                'Step6_Avg_Pockets_Per_Structure': (
                    round(pocket_info['Total_Pockets'] / pocket_info['PDB_Structures'], 2)
                    if pocket_info['PDB_Structures'] > 0 else 0
                )
            })
        
        logger.info("\n" + "="*80)
        logger.info("POCKET COUNT SUMMARY")
        logger.info("="*80)
        
        for result in results:
            logger.info(f"{result['Gene']}: {result['Step6_Total_Pockets']} pockets across {result['Step6_PDB_Structures']} structures")
        
        return pd.DataFrame(results)


class SubcellularLocationAnalyzer:
    """Analyze subcellular locations for druggability assessment"""
    
    def __init__(self):
        self.uniprot_base_url = "https://rest.uniprot.org/uniprotkb/search"
        
    def load_step5_data(self, filepath: str = 'Step5.csv') -> pd.DataFrame:
        """Load Step 5 pocket data from CSV"""
        logger.info(f"Loading Step 5 data from {filepath}")
        
        try:
            df = pd.read_csv(filepath)
            logger.info(f"Successfully loaded {len(df)} pocket entries")
            logger.info(f"Columns: {df.columns.tolist()}")
            return df
            
        except FileNotFoundError:
            logger.error(f"File not found: {filepath}")
            raise
        except Exception as e:
            logger.error(f"Error loading {filepath}: {e}")
            raise
    
    def _find_gene_column(self, df: pd.DataFrame) -> str:
        """Find the gene identifier column"""
        possible_names = ['Gene', 'gene', 'Gene_Symbol', 'gene_symbol', 
                         'Gene_Name', 'gene_name', 'Protein', 'protein']
        
        for col in possible_names:
            if col in df.columns:
                return col
        
        raise ValueError(f"Could not find gene column. Available: {df.columns.tolist()}")
    
    def query_uniprot_batch(self, gene_list: List[str]) -> Dict[str, Dict]:
        """Query UniProt for multiple genes with rate limiting"""
        results = {}
        total = len(gene_list)
        
        logger.info(f"\nQuerying UniProt for {total} genes...")
        logger.info("="*80)
        
        for idx, gene in enumerate(gene_list, 1):
            logger.info(f"[{idx}/{total}] Querying {gene}...")
            
            params = {
                'query': f'(gene:{gene}) AND (reviewed:true) AND (organism_id:9606)',
                'fields': 'accession,gene_primary,gene_names,cc_subcellular_location,ft_topo_dom,cc_function,protein_name',
                'format': 'json',
                'size': 1
            }
            
            try:
                response = requests.get(self.uniprot_base_url, params=params, timeout=30)
                response.raise_for_status()
                data = response.json()
                
                if data['results']:
                    entry = data['results'][0]
                    location_info = self._extract_location_data(entry, gene)
                    results[gene] = location_info
                    logger.info(f"  ? Location: {location_info['Step6_Subcellular_Location']}")
                else:
                    logger.warning(f"  ? No UniProt entry found")
                    results[gene] = self._create_empty_result(gene)
                
                time.sleep(0.5)  # Rate limiting
                
            except Exception as e:
                logger.error(f"  ? Error: {e}")
                results[gene] = self._create_empty_result(gene, error=str(e))
        
        return results
    
    def _extract_location_data(self, entry: Dict, gene: str) -> Dict:
        """Extract complete subcellular location information from UniProt entry"""
        
        location_data = {
            'Gene': gene,
            'Step6_UniProt_ID': entry.get('primaryAccession', 'N/A'),
            'Step6_Protein_Name': '',
            'Step6_Subcellular_Location': '',
            'Step6_Topology': '',
            'Step6_Function': '',
            'Step6_Data_Source': 'UniProt'
        }
        
        # Extract protein name
        if 'proteinDescription' in entry:
            rec_name = entry['proteinDescription'].get('recommendedName', {})
            if 'fullName' in rec_name:
                location_data['Step6_Protein_Name'] = rec_name['fullName'].get('value', '')
        
        # Extract subcellular location
        locations = []
        topologies = []
        
        if 'comments' in entry:
            for comment in entry['comments']:
                if comment['commentType'] == 'SUBCELLULAR LOCATION':
                    for loc in comment.get('subcellularLocations', []):
                        if 'location' in loc:
                            locations.append(loc['location']['value'])
                        if 'topology' in loc:
                            topologies.append(loc['topology']['value'])
                
                elif comment['commentType'] == 'FUNCTION':
                    texts = comment.get('texts', [])
                    if texts:
                        function_text = texts[0].get('value', '')
                        # Truncate if too long
                        location_data['Step6_Function'] = function_text[:500] + ('...' if len(function_text) > 500 else '')
        
        # Handle topology from features
        if 'features' in entry:
            for feature in entry['features']:
                if feature['type'] == 'Topological domain':
                    topo_desc = feature.get('description', '')
                    if topo_desc and topo_desc not in topologies:
                        topologies.append(topo_desc)
        
        location_data['Step6_Subcellular_Location'] = '; '.join(locations) if locations else 'Unknown'
        location_data['Step6_Topology'] = '; '.join(topologies) if topologies else 'N/A'
        
        return location_data
    
    def _create_empty_result(self, gene: str, error: Optional[str] = None) -> Dict:
        """Create empty result for genes not found"""
        return {
            'Gene': gene,
            'Step6_UniProt_ID': 'Not Found' if not error else 'Error',
            'Step6_Protein_Name': 'N/A',
            'Step6_Subcellular_Location': 'Not Found' if not error else f'Error: {error}',
            'Step6_Topology': 'N/A',
            'Step6_Function': 'N/A',
            'Step6_Data_Source': 'UniProt'
        }
    
    def assign_accessibility_score(self, location_str: str, topology_str: str) -> Dict:
        """Assign comprehensive accessibility and druggability scores"""
        
        location_lower = location_str.lower()
        topology_lower = topology_str.lower()
        
        # Priority 1: Secreted/Extracellular (Highest druggability)
        if any(term in location_lower for term in [
            'secreted', 'extracellular space', 'extracellular region',
            'cell surface', 'extracellular matrix', 'extracellular'
        ]):
            return {
                'Step6_Accessibility_Score': 4,
                'Step6_Location_Category': 'Extracellular/Secreted',
                'Step6_Druggability_Tier': 'Tier 1 - High',
                'Step6_Accessibility_Level': 'Excellent',
                'Step6_Accessibility_Description': 'Direct access via bloodstream - no barrier crossing required',
                'Step6_Drug_Modality_Options': 'Small molecule; Antibody; Peptide; Protein therapeutic',
                'Step6_Targeting_Difficulty': 'Easy',
                'Step6_Clinical_Precedent': 'High - many approved drugs',
                'Step6_Druggability_Rationale': 'Extracellular targets are highly accessible and have numerous successful examples'
            }
        
        # Priority 2: Membrane proteins (High druggability)
        elif any(term in location_lower for term in [
            'membrane', 'cell membrane', 'plasma membrane', 'receptor', 'transmembrane'
        ]):
            # Check if extracellular domain exists
            has_extracellular = any(term in topology_lower for term in [
                'extracellular', 'outside', 'external'
            ])
            
            if has_extracellular:
                return {
                    'Step6_Accessibility_Score': 3,
                    'Step6_Location_Category': 'Membrane (Extracellular Domain)',
                    'Step6_Druggability_Tier': 'Tier 1 - High',
                    'Step6_Accessibility_Level': 'Very Good',
                    'Step6_Accessibility_Description': 'Extracellular domain accessible from bloodstream',
                    'Step6_Drug_Modality_Options': 'Small molecule; Antibody (for extracellular epitopes)',
                    'Step6_Targeting_Difficulty': 'Easy-Moderate',
                    'Step6_Clinical_Precedent': 'High - GPCRs, RTKs, ion channels',
                    'Step6_Druggability_Rationale': 'Membrane proteins with extracellular domains are prime drug targets'
                }
            else:
                return {
                    'Step6_Accessibility_Score': 3,
                    'Step6_Location_Category': 'Membrane (Intramembrane/Intracellular)',
                    'Step6_Druggability_Tier': 'Tier 2 - Medium-High',
                    'Step6_Accessibility_Level': 'Good',
                    'Step6_Accessibility_Description': 'Embedded in membrane - requires lipophilic drugs',
                    'Step6_Drug_Modality_Options': 'Small molecule (membrane permeable)',
                    'Step6_Targeting_Difficulty': 'Moderate',
                    'Step6_Clinical_Precedent': 'Moderate - some ion channels, transporters',
                    'Step6_Druggability_Rationale': 'Intramembrane targets accessible but require specific drug properties'
                }
        
        # Priority 3: Cytoplasmic (Moderate druggability)
        elif any(term in location_lower for term in [
            'cytoplasm', 'cytosol', 'cytoplasmic'
        ]):
            return {
                'Step6_Accessibility_Score': 2,
                'Step6_Location_Category': 'Cytoplasmic',
                'Step6_Druggability_Tier': 'Tier 3 - Medium',
                'Step6_Accessibility_Level': 'Moderate',
                'Step6_Accessibility_Description': 'Intracellular - drug must cross plasma membrane',
                'Step6_Drug_Modality_Options': 'Small molecule (cell permeable); RNA therapeutics',
                'Step6_Targeting_Difficulty': 'Moderate-Difficult',
                'Step6_Clinical_Precedent': 'Moderate - kinases, proteases',
                'Step6_Druggability_Rationale': 'Requires cell-permeable compounds; selectivity can be challenging',
            }
        
        # Priority 4: Nuclear (Lower druggability)
        elif any(term in location_lower for term in [
            'nucleus', 'nuclear', 'chromatin'
        ]):
            return {
                'Step6_Accessibility_Score': 1,
                'Step6_Location_Category': 'Nuclear',
                'Step6_Druggability_Tier': 'Tier 4 - Medium-Low',
                'Step6_Accessibility_Level': 'Poor',
                'Step6_Accessibility_Description': 'Nuclear - drug must cross plasma and nuclear membranes',
                'Step6_Drug_Modality_Options': 'Small molecule (nuclear permeable); Oligonucleotides',
                'Step6_Targeting_Difficulty': 'Difficult',
                'Step6_Clinical_Precedent': 'Low-Moderate - transcription factors, chromatin modifiers',
                'Step6_Druggability_Rationale': 'Requires crossing two membranes; historically challenging targets',
            }
        
        # Other organelles
        elif any(term in location_lower for term in [
            'mitochondrion', 'mitochondrial', 'endoplasmic reticulum', 'golgi',
            'lysosome', 'peroxisome', 'endosome'
        ]):
            organelle = 'Unknown'
            if 'mitochondr' in location_lower:
                organelle = 'Mitochondrial'
            elif 'endoplasmic' in location_lower:
                organelle = 'Endoplasmic Reticulum'
            elif 'golgi' in location_lower:
                organelle = 'Golgi Apparatus'
            elif 'lysosome' in location_lower or 'lysosomal' in location_lower:
                organelle = 'Lysosomal'
            
            return {
                'Step6_Accessibility_Score': 1,
                'Step6_Location_Category': f'Organellar ({organelle})',
                'Step6_Druggability_Tier': 'Tier 4 - Medium-Low',
                'Step6_Accessibility_Level': 'Poor',
                'Step6_Accessibility_Description': f'Organelle-specific - requires targeted delivery to {organelle}',
                'Step6_Drug_Modality_Options': 'Small molecule (organelle-targeted); Peptide conjugates',
                'Step6_Targeting_Difficulty': 'Difficult',
                'Step6_Clinical_Precedent': 'Low - limited approved examples',
                'Step6_Druggability_Rationale': 'Organellar targeting remains a significant challenge',
            }
        
        # Unknown
        else:
            return {
                'Step6_Accessibility_Score': 0,
                'Step6_Location_Category': 'Unknown',
                'Step6_Druggability_Tier': 'Unassessed',
                'Step6_Accessibility_Level': 'Unknown',
                'Step6_Accessibility_Description': 'Location not determined - requires experimental validation',
                'Step6_Drug_Modality_Options': 'To be determined',
                'Step6_Targeting_Difficulty': 'Unknown',
                'Step6_Clinical_Precedent': 'Unknown',
                'Step6_Druggability_Rationale': 'Cannot assess without location information',
            }
    
    def create_complete_step6(self, step5_df: pd.DataFrame, 
                            pocket_counts_df: pd.DataFrame) -> pd.DataFrame:
        """
        Create complete Step 6 CSV with all fields
        
        Args:
            step5_df: DataFrame from Step5.csv
            pocket_counts_df: DataFrame with actual pocket counts
        
        Returns:
            Complete Step6 DataFrame
        """
        # Find gene column
        gene_col = self._find_gene_column(step5_df)
        
        # Get unique genes
        unique_genes = step5_df[gene_col].unique().tolist()
        logger.info(f"\nProcessing {len(unique_genes)} unique genes")
        
        # Query UniProt for all genes
        location_results = self.query_uniprot_batch(unique_genes)
        location_df = pd.DataFrame([location_results[gene] for gene in unique_genes])
        
        # Assign accessibility scores
        logger.info("\n" + "="*80)
        logger.info("ASSIGNING ACCESSIBILITY AND DRUGGABILITY SCORES")
        logger.info("="*80)
        
        accessibility_data = location_df.apply(
            lambda row: self.assign_accessibility_score(
                row['Step6_Subcellular_Location'],
                row['Step6_Topology']
            ),
            axis=1
        )
        
        accessibility_df = pd.DataFrame(accessibility_data.tolist())
        
        # Combine all location data
        location_complete = pd.concat([location_df, accessibility_df], axis=1)
        location_complete = location_complete.rename(columns={'Gene': gene_col})
        
        # Merge with pocket counts
        location_complete = location_complete.merge(
            pocket_counts_df, 
            left_on=gene_col, 
            right_on='Gene', 
            how='left'
        )
        
        # Drop duplicate Gene column if it exists
        if 'Gene' in location_complete.columns and gene_col != 'Gene':
            location_complete = location_complete.drop('Gene', axis=1)
        
        # Merge with step5_df
        logger.info("\nMerging all data...")
        combined_df = step5_df.merge(location_complete, on=gene_col, how='left')
        
        logger.info(f"? Complete Step 6 data: {len(combined_df)} rows")
        
        return combined_df


def main():
    """Main execution function"""
    
    # ===== CONFIGURATION =====
    STEP5_INPUT = 'Step5.csv'
    STEP6_OUTPUT = 'Step6.csv'
    BINDING_POCKETS_DIR = 'binding_pockets'
    # =========================
    
    try:
        logger.info("="*80)
        logger.info("STARTING STEP 6: SUBCELLULAR LOCATION ANALYSIS")
        logger.info("="*80)
        
        # Initialize components
        analyzer = SubcellularLocationAnalyzer()
        pocket_counter = PocketCounter(BINDING_POCKETS_DIR)
        
        # Load Step 5 data
        step5_df = analyzer.load_step5_data(STEP5_INPUT)
        gene_col = analyzer._find_gene_column(step5_df)
        unique_genes = step5_df[gene_col].unique().tolist()
        
        # Count actual pockets from directory structure
        pocket_counts_df = pocket_counter.count_all_genes(unique_genes)
        
        # Create complete Step 6 data
        logger.info("\n" + "="*80)
        logger.info("CREATING COMPLETE STEP 6 DATASET")
        logger.info("="*80)
        
        step6_df = analyzer.create_complete_step6(step5_df, pocket_counts_df)
        
        # Save complete results
        logger.info(f"\nSaving complete results to {STEP6_OUTPUT}...")
        step6_df.to_csv(STEP6_OUTPUT, index=False)
        logger.info(f"? Saved {len(step6_df)} rows")
        
        # Generate protein-level summary
        logger.info("\n" + "="*80)
        logger.info("GENERATING SUMMARIES")
        logger.info("="*80)
        
        summary_cols = [
            gene_col, 'Step6_UniProt_ID', 'Step6_Protein_Name',
            'Step6_Total_Pockets', 'Step6_PDB_Structures',
            'Step6_Avg_Pockets_Per_Structure',
            'Step6_Subcellular_Location', 'Step6_Topology',
            'Step6_Location_Category', 'Step6_Accessibility_Score',
            'Step6_Accessibility_Level', 'Step6_Druggability_Tier',
            'Step6_Drug_Modality_Options', 'Step6_Targeting_Difficulty',
            'Step6_Clinical_Precedent', 'Step6_Function'
        ]
        
        available_cols = [col for col in summary_cols if col in step6_df.columns]
        protein_summary = step6_df[available_cols].drop_duplicates(subset=[gene_col])
        protein_summary = protein_summary.sort_values('Step6_Accessibility_Score', ascending=False)
        
        # Save protein summary
        summary_file = 'Step6_protein_summary.csv'
        protein_summary.to_csv(summary_file, index=False)
        logger.info(f"? Protein summary saved: {summary_file}")
        
        # Save priority targets
        priority = protein_summary[protein_summary['Step6_Accessibility_Score'] >= 3]
        if len(priority) > 0:
            priority_file = 'Step6_priority_targets.csv'
            priority.to_csv(priority_file, index=False)
            logger.info(f"? Priority targets saved: {priority_file} ({len(priority)} targets)")
        
        # Display summary statistics
        logger.info("\n" + "="*80)
        logger.info("ANALYSIS SUMMARY")
        logger.info("="*80)
        logger.info(f"\nTotal genes analyzed: {len(unique_genes)}")
        logger.info(f"Total pocket entries: {len(step6_df)}")
        
        logger.info(f"\nLocation Distribution:")
        loc_counts = protein_summary['Step6_Location_Category'].value_counts()
        for loc, count in loc_counts.items():
            logger.info(f"  {loc}: {count}")
        
        logger.info(f"\nAccessibility Scores:")
        score_counts = protein_summary['Step6_Accessibility_Score'].value_counts().sort_index(ascending=False)
        for score, count in score_counts.items():
            logger.info(f"  Score {score}: {count} genes")
        
        # Display top 10 targets
        logger.info("\n" + "="*80)
        logger.info("TOP 10 DRUGGABLE TARGETS")
        logger.info("="*80)
        
        for idx, (_, row) in enumerate(protein_summary.head(10).iterrows(), 1):
            logger.info(f"\n{idx}. {row[gene_col]} ({row.get('Step6_UniProt_ID', 'N/A')})")
            logger.info(f"   Protein: {row.get('Step6_Protein_Name', 'N/A')}")
            logger.info(f"   Total Pockets: {row.get('Step6_Total_Pockets', 0)}")
            logger.info(f"   Location: {row.get('Step6_Location_Category', 'Unknown')}")
            logger.info(f"   Accessibility: {row.get('Step6_Accessibility_Score', 0)}/4 ({row.get('Step6_Accessibility_Level', 'Unknown')})")
            logger.info(f"   Druggability: {row.get('Step6_Druggability_Tier', 'Unknown')}")
            logger.info(f"   Drug Modalities: {row.get('Step6_Drug_Modality_Options', 'N/A')}")
        
        # Final summary
        logger.info("\n" + "="*80)
        logger.info("STEP 6 COMPLETE!")
        logger.info("="*80)
        logger.info(f"\nGenerated files:")
        logger.info(f"  ? {STEP6_OUTPUT}")
        logger.info(f"  ? {summary_file}")
        if len(priority) > 0:
            logger.info(f"  ? {priority_file}")
        logger.info(f"  ? step6_location_analysis.log")
        

        
    except Exception as e:
        logger.error(f"Error during analysis: {e}", exc_info=True)
        raise


if __name__ == "__main__":
    main()