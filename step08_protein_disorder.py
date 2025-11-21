"""
Step 8: Protein Disorder Analysis
Uses multiple methods: MobiDB API, AlphaFold pLDDT scores, and composition-based prediction
"""

import pandas as pd
import requests
import time
from typing import Dict, Optional, List
from pathlib import Path
import logging
import re

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('step8_disorder_analysis.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


class ProteinDisorderAnalyzer:
    """Analyze protein disorder using multiple prediction methods"""
    
    def __init__(self):
        self.uniprot_base_url = "https://rest.uniprot.org/uniprotkb"
        self.mobidb_base_url = "https://mobidb.bio.unipd.it/api"
        self.alphafold_base_url = "https://alphafold.ebi.ac.uk/api"
        self.cache = {}
        
    def load_step7_data(self, filepath: str = 'Step7.csv') -> pd.DataFrame:
        """Load Step 7 data"""
        logger.info(f"Loading Step 7 data from {filepath}")
        df = pd.read_csv(filepath)
        logger.info(f"Successfully loaded {len(df)} entries")
        return df
    
    def _find_gene_column(self, df: pd.DataFrame) -> str:
        """Find gene identifier column"""
        for col in ['Gene', 'gene', 'Gene_Symbol', 'gene_symbol', 'Gene_Name', 'gene_name']:
            if col in df.columns:
                return col
        raise ValueError(f"Could not find gene column. Available: {df.columns.tolist()}")
    
    def get_protein_sequence(self, gene_symbol: str, uniprot_id: str = None) -> Optional[str]:
        """Get protein sequence from UniProt"""
        
        logger.info(f"  Getting protein sequence...")
        
        try:
            if uniprot_id and uniprot_id not in ['N/A', 'Not Found', 'Error', 'Unknown']:
                url = f"{self.uniprot_base_url}/{uniprot_id}.fasta"
                response = requests.get(url, timeout=10)
            else:
                url = f"{self.uniprot_base_url}/search"
                params = {
                    'query': f'(gene:{gene_symbol}) AND (reviewed:true) AND (organism_id:9606)',
                    'format': 'fasta',
                    'size': 1
                }
                response = requests.get(url, params=params, timeout=10)
            
            if response.status_code == 200 and response.text:
                lines = response.text.strip().split('\n')
                if len(lines) > 1:
                    sequence = ''.join(lines[1:])
                    logger.info(f"    ? Retrieved sequence ({len(sequence)} aa)")
                    return sequence
            
            return None
            
        except Exception as e:
            logger.warning(f"    ? Error: {e}")
            return None
    
    def predict_disorder_mobidb(self, uniprot_id: str) -> Optional[Dict]:
        """
        Predict disorder using MobiDB (Database of protein disorder and mobility annotations)
        MobiDB aggregates predictions from multiple sources
        """
        
        logger.info(f"  Querying MobiDB...")
        
        try:
            # MobiDB API endpoint
            url = f"{self.mobidb_base_url}/download/{uniprot_id}"
            
            response = requests.get(url, timeout=15)
            
            if response.status_code == 200:
                data = response.json()
                
                # MobiDB provides consensus disorder predictions
                if 'mobidb_consensus' in data:
                    consensus = data['mobidb_consensus']
                    
                    # Get disorder regions
                    disorder_regions = []
                    if 'disorder' in consensus and 'regions' in consensus['disorder']:
                        disorder_regions = consensus['disorder']['regions']
                    
                    # Calculate disorder
                    sequence_length = data.get('length', 0)
                    
                    if sequence_length > 0:
                        disordered_residues = 0
                        for region in disorder_regions:
                            start = region.get('start', 0)
                            end = region.get('end', 0)
                            disordered_residues += (end - start + 1)
                        
                        disorder_percentage = (disordered_residues / sequence_length) * 100
                        
                        logger.info(f"    ? MobiDB: {disorder_percentage:.1f}% disorder")
                        
                        return {
                            'disorder_percentage': disorder_percentage,
                            'disordered_residues': disordered_residues,
                            'total_residues': sequence_length,
                            'disordered_regions': [(r['start'], r['end']) for r in disorder_regions],
                            'num_disordered_regions': len(disorder_regions),
                            'method': 'MobiDB Consensus'
                        }
                
                logger.info(f"    ? No consensus disorder data in MobiDB")
                return None
                
            else:
                logger.info(f"    ? MobiDB returned {response.status_code}")
                return None
                
        except Exception as e:
            logger.warning(f"    ? MobiDB error: {e}")
            return None
    
    def predict_disorder_alphafold(self, uniprot_id: str) -> Optional[Dict]:
        """
        Use AlphaFold pLDDT (predicted Local Distance Difference Test) scores
        Low pLDDT (<50) often indicates disorder
        """
        
        logger.info(f"  Checking AlphaFold pLDDT scores...")
        
        try:
            # AlphaFold API
            url = f"{self.alphafold_base_url}/prediction/{uniprot_id}"
            
            response = requests.get(url, timeout=15)
            
            if response.status_code == 200:
                data = response.json()
                
                if data and len(data) > 0:
                    # Get pLDDT scores
                    entry = data[0]
                    
                    # Check if we can access the structure file
                    pdb_url = entry.get('pdbUrl')
                    
                    if pdb_url:
                        # Download PDB to get pLDDT scores from B-factor column
                        pdb_response = requests.get(pdb_url, timeout=15)
                        
                        if pdb_response.status_code == 200:
                            pdb_text = pdb_response.text
                            
                            # Extract pLDDT scores from B-factor column (appears as temperature factor)
                            plddt_scores = []
                            for line in pdb_text.split('\n'):
                                if line.startswith('ATOM'):
                                    # B-factor is in columns 61-66
                                    try:
                                        plddt = float(line[60:66].strip())
                                        plddt_scores.append(plddt)
                                    except:
                                        pass
                            
                            if plddt_scores:
                                # pLDDT < 50 suggests disorder
                                disordered_residues = sum(1 for score in plddt_scores if score < 50)
                                total_residues = len(plddt_scores)
                                disorder_percentage = (disordered_residues / total_residues) * 100
                                
                                logger.info(f"    ? AlphaFold: {disorder_percentage:.1f}% disorder (pLDDT<50)")
                                
                                return {
                                    'disorder_percentage': disorder_percentage,
                                    'disordered_residues': disordered_residues,
                                    'total_residues': total_residues,
                                    'disordered_regions': [],
                                    'num_disordered_regions': 0,
                                    'method': 'AlphaFold pLDDT'
                                }
                
                logger.info(f"    ? Could not extract pLDDT scores")
                return None
                
            else:
                logger.info(f"    ? AlphaFold returned {response.status_code}")
                return None
                
        except Exception as e:
            logger.warning(f"    ? AlphaFold error: {e}")
            return None
    
    def predict_disorder_composition(self, sequence: str) -> Dict:
        """
        Composition-based disorder prediction
        Uses known disorder-promoting amino acid properties
        """
        
        logger.info(f"  Using composition-based prediction...")
        
        # Disorder-promoting amino acids (charge, proline, small)
        disorder_promoting = {
            'R': 1.0, 'K': 1.0,  # Positive charge
            'D': 1.0, 'E': 1.0,  # Negative charge
            'Q': 0.8, 'S': 0.8, 'P': 1.0, 'G': 0.9,  # Flexible
            'A': 0.5  # Small
        }
        
        # Order-promoting amino acids (hydrophobic, aromatic)
        order_promoting = {
            'W': -1.0, 'F': -0.9, 'Y': -0.9,  # Aromatic
            'I': -0.8, 'L': -0.8, 'V': -0.8, 'M': -0.7,  # Hydrophobic
            'C': -0.6  # Disulfide bridges
        }
        
        # Calculate disorder propensity
        disorder_score = 0
        for aa in sequence:
            disorder_score += disorder_promoting.get(aa, 0)
            disorder_score += order_promoting.get(aa, 0)
        
        # Normalize to percentage (rough estimate)
        # Average disorder propensity per residue, scaled to percentage
        avg_propensity = disorder_score / len(sequence)
        
        # Map propensity to disorder percentage (calibrated estimate)
        if avg_propensity > 0.3:
            disorder_percentage = min(70 + (avg_propensity - 0.3) * 100, 90)
        elif avg_propensity > 0:
            disorder_percentage = 30 + (avg_propensity * 133)
        else:
            disorder_percentage = max(5, 30 + (avg_propensity * 60))
        
        # Estimate disordered residues
        disordered_residues = int(len(sequence) * disorder_percentage / 100)
        
        logger.info(f"    ! Estimated disorder: ~{disorder_percentage:.1f}%")
        
        return {
            'disorder_percentage': disorder_percentage,
            'disordered_residues': disordered_residues,
            'total_residues': len(sequence),
            'disordered_regions': [],
            'num_disordered_regions': 0,
            'method': 'Composition-Based'
        }
    
    def analyze_protein_disorder(self, gene_symbol: str, uniprot_id: str = None) -> Dict:
        """Complete disorder analysis using multiple methods"""
        
        cache_key = gene_symbol.upper()
        if cache_key in self.cache:
            return self.cache[cache_key]
        
        logger.info(f"\n[Disorder Analysis] {gene_symbol}")
        
        disorder_data = {
            'Gene': gene_symbol,
            'Sequence_Available': False,
            'Disorder_Analysis_Available': False,
            'Disorder_Percentage': 0.0,
            'Disordered_Residues': 0,
            'Total_Residues': 0,
            'Disorder_Category': 'Unknown',
            'Num_Disordered_Regions': 0,
            'Longest_Disordered_Region': 0,
            'Disorder_Method': 'None',
            'Druggability_Impact': 'Unknown'
        }
        
        # Get sequence
        sequence = self.get_protein_sequence(gene_symbol, uniprot_id)
        
        if not sequence:
            self.cache[cache_key] = disorder_data
            return disorder_data
        
        disorder_data['Sequence_Available'] = True
        disorder_data['Total_Residues'] = len(sequence)
        
        # Try methods in order of preference
        disorder_result = None
        
        # Method 1: MobiDB (if we have UniProt ID)
        if uniprot_id and uniprot_id not in ['N/A', 'Not Found', 'Error', 'Unknown']:
            disorder_result = self.predict_disorder_mobidb(uniprot_id)
        
        # Method 2: AlphaFold pLDDT (if we have UniProt ID)
        if not disorder_result and uniprot_id and uniprot_id not in ['N/A', 'Not Found', 'Error', 'Unknown']:
            disorder_result = self.predict_disorder_alphafold(uniprot_id)
        
        # Method 3: Composition-based (fallback)
        if not disorder_result:
            disorder_result = self.predict_disorder_composition(sequence)
        
        if disorder_result:
            disorder_data['Disorder_Analysis_Available'] = True
            disorder_data['Disorder_Percentage'] = round(disorder_result['disorder_percentage'], 2)
            disorder_data['Disordered_Residues'] = disorder_result['disordered_residues']
            disorder_data['Num_Disordered_Regions'] = disorder_result['num_disordered_regions']
            disorder_data['Disorder_Method'] = disorder_result['method']
            
            # Calculate longest region
            if disorder_result['disordered_regions']:
                region_lengths = [end - start + 1 for start, end in disorder_result['disordered_regions']]
                disorder_data['Longest_Disordered_Region'] = max(region_lengths)
            
            # Categorize
            disorder_pct = disorder_data['Disorder_Percentage']
            if disorder_pct < 10:
                disorder_data['Disorder_Category'] = 'Highly Ordered'
            elif disorder_pct < 30:
                disorder_data['Disorder_Category'] = 'Mostly Ordered'
            elif disorder_pct < 50:
                disorder_data['Disorder_Category'] = 'Moderately Disordered'
            elif disorder_pct < 70:
                disorder_data['Disorder_Category'] = 'Highly Disordered'
            else:
                disorder_data['Disorder_Category'] = 'Intrinsically Disordered'
            
            # Druggability impact
            if disorder_pct < 20:
                disorder_data['Druggability_Impact'] = 'Excellent - Well-structured'
            elif disorder_pct < 50:
                disorder_data['Druggability_Impact'] = 'Good - Druggable (<50% disorder)'
            elif disorder_pct < 70:
                disorder_data['Druggability_Impact'] = 'Challenging - High disorder (>50%)'
            else:
                disorder_data['Druggability_Impact'] = 'Poor - Mostly disordered (>70%)'
        
        self.cache[cache_key] = disorder_data
        time.sleep(1.0)
        
        return disorder_data
    
    def calculate_disorder_scores(self, disorder_data: Dict) -> Dict:
        """Calculate disorder-based druggability scores"""
        
        scores = {
            'Step8_Disorder_Score': 0,
            'Step8_Structural_Quality_Score': 0,
            'Step8_Passes_50_Percent_Rule': False,
            'Step8_Druggability_Confidence': 'Unknown'
        }
        
        if not disorder_data['Disorder_Analysis_Available']:
            scores['Step8_Druggability_Confidence'] = 'Low - No Data'
            return scores
        
        disorder_pct = disorder_data['Disorder_Percentage']
        
        # Apply 50% rule
        scores['Step8_Passes_50_Percent_Rule'] = (disorder_pct < 50)
        
        # Disorder score (0-5)
        if disorder_pct < 10:
            scores['Step8_Disorder_Score'] = 5
        elif disorder_pct < 20:
            scores['Step8_Disorder_Score'] = 4
        elif disorder_pct < 50:
            scores['Step8_Disorder_Score'] = 3
        elif disorder_pct < 70:
            scores['Step8_Disorder_Score'] = 2
        else:
            scores['Step8_Disorder_Score'] = 1
        
        # Structural quality
        num_regions = disorder_data['Num_Disordered_Regions']
        longest_region = disorder_data['Longest_Disordered_Region']
        total_residues = disorder_data['Total_Residues']
        
        if total_residues > 0 and longest_region > 0:
            longest_region_pct = (longest_region / total_residues) * 100
            
            if longest_region_pct < 10 and num_regions <= 2:
                scores['Step8_Structural_Quality_Score'] = 5
            elif longest_region_pct < 30 and num_regions <= 5:
                scores['Step8_Structural_Quality_Score'] = 4
            elif longest_region_pct < 50:
                scores['Step8_Structural_Quality_Score'] = 3
            else:
                scores['Step8_Structural_Quality_Score'] = 2
        else:
            scores['Step8_Structural_Quality_Score'] = 5
        
        # Confidence
        method = disorder_data['Disorder_Method']
        if 'MobiDB' in method:
            scores['Step8_Druggability_Confidence'] = 'High - MobiDB Consensus'
        elif 'AlphaFold' in method:
            scores['Step8_Druggability_Confidence'] = 'High - AlphaFold pLDDT'
        elif 'Composition' in method:
            scores['Step8_Druggability_Confidence'] = 'Medium - Estimated'
        
        return scores
    
    def add_disorder_data(self, step7_df: pd.DataFrame) -> pd.DataFrame:
        """Add disorder analysis to Step 7"""
        
        gene_col = self._find_gene_column(step7_df)
        unique_genes = step7_df[gene_col].unique().tolist()
        
        # Get UniProt IDs
        uniprot_col = None
        for col in ['Step6_UniProt_ID', 'UniProt_ID', 'uniprot_id']:
            if col in step7_df.columns:
                uniprot_col = col
                break
        
        logger.info(f"\nAnalyzing protein disorder for {len(unique_genes)} genes...")
        logger.info("="*80)
        
        all_disorder_data = []
        
        for idx, gene in enumerate(unique_genes, 1):
            uniprot_id = None
            if uniprot_col:
                gene_rows = step7_df[step7_df[gene_col] == gene]
                if len(gene_rows) > 0:
                    uniprot_id = gene_rows.iloc[0][uniprot_col]
            
            disorder_data = self.analyze_protein_disorder(gene, uniprot_id)
            scores = self.calculate_disorder_scores(disorder_data)
            
            combined = {
                gene_col: gene,
                'Step8_Sequence_Available': disorder_data['Sequence_Available'],
                'Step8_Disorder_Analysis_Available': disorder_data['Disorder_Analysis_Available'],
                'Step8_Disorder_Percentage': disorder_data['Disorder_Percentage'],
                'Step8_Disordered_Residues': disorder_data['Disordered_Residues'],
                'Step8_Total_Residues': disorder_data['Total_Residues'],
                'Step8_Disorder_Category': disorder_data['Disorder_Category'],
                'Step8_Num_Disordered_Regions': disorder_data['Num_Disordered_Regions'],
                'Step8_Longest_Disordered_Region': disorder_data['Longest_Disordered_Region'],
                'Step8_Disorder_Method': disorder_data['Disorder_Method'],
                'Step8_Druggability_Impact': disorder_data['Druggability_Impact'],
                **scores
            }
            
            all_disorder_data.append(combined)
            
            logger.info(f"  [{idx}/{len(unique_genes)}] {gene}: {disorder_data['Disorder_Percentage']:.1f}% disorder, " +
                       f"Passes <50% rule: {scores['Step8_Passes_50_Percent_Rule']}")
        
        disorder_df = pd.DataFrame(all_disorder_data)
        combined_df = step7_df.merge(disorder_df, on=gene_col, how='left')
        
        logger.info(f"\n? Added disorder data to {len(combined_df)} rows")
        return combined_df
    
    def generate_summary_report(self, df: pd.DataFrame):
        """Generate summary statistics"""
        
        gene_col = self._find_gene_column(df)
        protein_level = df.drop_duplicates(subset=[gene_col])
        
        logger.info("\n" + "="*80)
        logger.info("STEP 8 DISORDER ANALYSIS SUMMARY")
        logger.info("="*80)
        
        logger.info(f"\nData Availability:")
        available = protein_level['Step8_Disorder_Analysis_Available'].sum()
        logger.info(f"  Disorder data: {available}/{len(protein_level)}")
        
        logger.info(f"\n50% Disorder Rule:")
        passes_rule = protein_level['Step8_Passes_50_Percent_Rule'].sum()
        logger.info(f"  Passes (<50% disorder): {passes_rule}/{available}")
        logger.info(f"  Fails (=50% disorder): {available - passes_rule}/{available}")
        
        logger.info(f"\nDisorder Categories:")
        for cat, count in protein_level['Step8_Disorder_Category'].value_counts().items():
            logger.info(f"  {cat}: {count}")
        
        logger.info(f"\nMethods Used:")
        for method, count in protein_level['Step8_Disorder_Method'].value_counts().items():
            logger.info(f"  {method}: {count}")
        
        return protein_level


def main():
    STEP7_INPUT = 'Step7.csv'
    STEP8_OUTPUT = 'Step8.csv'
    
    try:
        logger.info("="*80)
        logger.info("STEP 8: PROTEIN DISORDER ANALYSIS")
        logger.info("="*80)
        
        analyzer = ProteinDisorderAnalyzer()
        step7_df = analyzer.load_step7_data(STEP7_INPUT)
        step8_df = analyzer.add_disorder_data(step7_df)
        
        logger.info(f"\nSaving to {STEP8_OUTPUT}...")
        step8_df.to_csv(STEP8_OUTPUT, index=False)
        logger.info(f"? Saved {len(step8_df)} rows")
        
        protein_level = analyzer.generate_summary_report(step8_df)
        
        gene_col = analyzer._find_gene_column(step8_df)
        summary = protein_level[[
            gene_col, 'Step8_Disorder_Percentage', 'Step8_Passes_50_Percent_Rule',
            'Step8_Disorder_Category', 'Step8_Druggability_Impact'
        ]].sort_values('Step8_Disorder_Percentage')
        
        summary.to_csv('Step8_protein_summary.csv', index=False)
        logger.info(f"\n? Summary saved")
        
        logger.info("\n" + "="*80)
        logger.info("STEP 8 COMPLETE!")
        logger.info("="*80)
        
    except Exception as e:
        logger.error(f"Error: {e}", exc_info=True)
        raise


if __name__ == "__main__":
    main()